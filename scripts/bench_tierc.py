"""A/B benchmarks for the Tier C (pure-Python) optimizations.

Each benchmark holds a **local copy of the 0.3.3 implementation** and runs it
against the current one, interleaved in a single process (see ``bench_ab``).
That keeps the comparison honest without needing two source trees — the trap
documented in the optimization plan's §6, where the editable install's meta-path
finder silently serves the working tree no matter what ``sys.path`` says.

Run::

    python scripts/bench_tierc.py
"""
from __future__ import annotations

import io
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from bench_ab import ab, header  # noqa: E402

from zna.cli import (  # noqa: E402
    _fragment_units, _pair_interleaved, _read_key, parse_fastq_keyed,
    parse_fastq_with_names, get_base_name, get_read_suffix_number,
)
from zna.core import (  # noqa: E402
    COMPRESSION_ZSTD, ZnaHeader, ZnaWriter, _RC_FULL_BY_ENDS, _flags_from_ends,
)

RNG = random.Random(20260811)
N_READS = 100_000
READ_LEN = 150


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

def make_reads(n, with_n_fraction=0.02):
    reads = []
    for i in range(n):
        seq = "".join(RNG.choices("ACGT", k=READ_LEN))
        if RNG.random() < with_n_fraction:
            pos = RNG.randrange(READ_LEN)
            seq = seq[:pos] + "N" + seq[pos + 1:]
        reads.append(seq)
    return reads


def make_interleaved_fastq(n_pairs):
    """Interleaved FASTQ bytes with realistic Illumina-style names."""
    out = []
    for i in range(n_pairs):
        base = f"INSTR:42:FLOWCELL:1:1101:{1000 + i}:{2000 + i}"
        for mate in (1, 2):
            seq = "".join(RNG.choices("ACGT", k=READ_LEN))
            out.append(f"@{base}/{mate} some comment here\n{seq}\n+\n{'I' * READ_LEN}\n")
    return "".join(out).encode("ascii")


SEQS = make_reads(N_READS)
#: The write-path benchmarks encode for real, and encoding an N without an
#: N-policy is an error, so those use a clean set.
SEQS_CLEAN = make_reads(N_READS, with_n_fraction=0.0)
UNITS_SE = [[(s, False, False, False)] for s in SEQS]
FASTQ_IL = make_interleaved_fastq(N_READS // 2)


# ---------------------------------------------------------------------------
# 1. N-drop filter
# ---------------------------------------------------------------------------

def bench_ndrop():
    units = UNITS_SE

    def old():
        kept = 0
        for unit in units:
            if any('N' in rec[0].upper() for rec in unit):
                continue
            kept += 1
        return kept

    def new():
        kept = 0
        for unit in units:
            if any(('N' in rec[0] or 'n' in rec[0]) for rec in unit):
                continue
            kept += 1
        return kept

    assert old() == new(), "N-drop filter changed which reads are dropped"
    return ab("N-drop filter (no .upper() copy)", old, new)


# ---------------------------------------------------------------------------
# 2. Single-end fragment grouping, measured through a real writer
# ---------------------------------------------------------------------------

def _hdr():
    return ZnaHeader(read_group="bench", seq_len_bytes=2,
                     compression_method=COMPRESSION_ZSTD, compression_level=3)


def bench_single_end_loop():
    seqs = SEQS_CLEAN

    def old():
        # Uses the *new* N-drop test, so the only difference from `new` is the
        # fragment grouping itself.  The `.upper()` copy is benched separately.
        buf = io.BytesIO()
        with ZnaWriter(buf, _hdr()) as w:
            count = 0
            stream = ((s, False, False, False) for s in seqs)
            for unit in _fragment_units(stream):
                if any(('N' in rec[0] or 'n' in rec[0]) for rec in unit):
                    continue
                for rec in unit:
                    w.write_record(rec[0], rec[1], rec[2], rec[3],
                                   is_full_fragment=False)
                count += len(unit)
        return buf.getbuffer().nbytes

    def new():
        buf = io.BytesIO()
        with ZnaWriter(buf, _hdr()) as w:
            count = 0
            write_record = w.write_record
            stream = ((s, False, False, False) for s in seqs)
            for rec in stream:
                seq = rec[0]
                if 'N' in seq or 'n' in seq:
                    continue
                write_record(seq, rec[1], rec[2], rec[3], is_full_fragment=False)
                count += 1
        return buf.getbuffer().nbytes

    assert old() == new(), "single-end fast path produced a different file"
    return ab("single-end encode loop (no fragment grouping)", old, new)


# ---------------------------------------------------------------------------
# 3. _flags_from_ends per record
# ---------------------------------------------------------------------------

def bench_flags_from_ends():
    # A realistic mix of the three reachable end states.
    ends = [(True, False), (False, True), (True, True)] * (N_READS // 3)

    def old():
        acc = 0
        for has_start, has_end in ends:
            is_rc, is_full = _flags_from_ends(has_start, has_end)
            acc += (8 if is_rc else 0) | (16 if is_full else 0)
        return acc

    def new():
        acc = 0
        table = _RC_FULL_BY_ENDS
        for has_start, has_end in ends:
            is_rc, is_full = table[(2 if has_start else 0) | (1 if has_end else 0)]
            acc += (8 if is_rc else 0) | (16 if is_full else 0)
        return acc

    assert old() == new(), "_flags_from_ends table disagrees with the function"
    return ab("_flags_from_ends -> table lookup (isolated)", old, new)


def bench_flags_from_ends_in_context():
    """The same change measured where it actually runs: a pass-through write
    loop, which is what ``shuffle`` and ``encode`` of a normalized file do."""
    recs = [(s, True, i % 2 == 0, i % 2 == 1, i % 3 != 1, i % 3 == 1)
            for i, s in enumerate(SEQS_CLEAN)]

    def hdr_norm():
        return ZnaHeader(read_group="bench", seq_len_bytes=2,
                         strand_normalized=True,
                         compression_method=COMPRESSION_ZSTD, compression_level=3)

    def old():
        buf = io.BytesIO()
        with ZnaWriter(buf, hdr_norm(), preserve_normalization=True) as w:
            write_record = w.write_record
            for rec in recs:
                is_rc, is_full = _flags_from_ends(rec[4], rec[5])
                write_record(rec[0], rec[1], rec[2], rec[3],
                             is_rc=is_rc, is_full_fragment=is_full)
        return buf.getbuffer().nbytes

    def new():
        buf = io.BytesIO()
        table = _RC_FULL_BY_ENDS
        with ZnaWriter(buf, hdr_norm(), preserve_normalization=True) as w:
            write_record = w.write_record
            for rec in recs:
                is_rc, is_full = table[(2 if rec[4] else 0) | (1 if rec[5] else 0)]
                write_record(rec[0], rec[1], rec[2], rec[3],
                             is_rc=is_rc, is_full_fragment=is_full)
        return buf.getbuffer().nbytes

    assert old() == new(), "table lookup changed the written file"
    return ab("_flags_from_ends in a pass-through write loop", old, new)


# ---------------------------------------------------------------------------
# 4. Interleaved pairing
# ---------------------------------------------------------------------------

def _old_parse_with_names(fh):
    """0.3.3 parse_fastq_with_names — decodes the whole header to str."""
    return parse_fastq_with_names(fh)


def _old_pair_interleaved(records):
    """0.3.3 _pair_interleaved — re-derives the base name per comparison."""
    prev = None
    for name, payload in records:
        if prev is None:
            prev = (name, payload)
            continue
        prev_name, prev_payload = prev
        if get_base_name(prev_name) == get_base_name(name):
            get_read_suffix_number(prev_name)
            get_read_suffix_number(name)
            yield prev_payload, True, True, False
            yield payload, True, False, True
            prev = None
        else:
            yield prev_payload, False, False, False
            prev = (name, payload)
    if prev is not None:
        yield prev[1], False, False, False


def bench_interleaved_pairing():
    data = FASTQ_IL

    def old():
        return sum(1 for _ in _old_pair_interleaved(
            _old_parse_with_names(io.BytesIO(data))))

    def new():
        return sum(1 for _ in _pair_interleaved(parse_fastq_keyed(io.BytesIO(data))))

    # Same pairing decisions, not just the same count.
    a = list(_old_pair_interleaved(_old_parse_with_names(io.BytesIO(data))))
    b = list(_pair_interleaved(parse_fastq_keyed(io.BytesIO(data))))
    assert a == b, "bytes-native pairing changed the pairing decisions"
    return ab("interleaved pairing (bytes-native)", old, new)


def bench_read_key():
    names = [f"INSTR:42:FLOWCELL:1:1101:{1000 + i}:{2000 + i}/1 some comment here"
             for i in range(20_000)]
    raw = [n.encode("ascii") for n in names]

    def old():
        return sum(len(get_base_name(n)) + get_read_suffix_number(n) for n in names)

    def new():
        acc = 0
        for r in raw:
            base, suffix = _read_key(r)
            acc += len(base) + suffix
        return acc

    assert old() == new(), "_read_key disagrees with get_base_name/suffix"
    return ab("_read_key vs get_base_name+get_read_suffix_number", old, new)


def main():
    header("Tier C — Python-only optimizations (0.3.3 baseline vs current)")
    results = [
        bench_ndrop(),
        bench_single_end_loop(),
        bench_flags_from_ends(),
        bench_flags_from_ends_in_context(),
        bench_read_key(),
        bench_interleaved_pairing(),
    ]
    print("  " + "-" * 92)
    sig = [r for r in results if r["significant"]]
    print(f"\n  {len(sig)}/{len(results)} cleared the {5.0}% noise floor.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
