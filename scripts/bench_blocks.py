"""What is ``blocks()`` actually worth, on top of an already-fixed decoder?

The optimization plan put ``blocks()`` last and heavily caveated: its headline
+274% was measured against the *old* decoder and a loader that did nothing but
append strings.  This measures it against the current one, with a per-record
step, and in the sharded shape khorana's loader actually uses.

Run::

    python scripts/bench_blocks.py
"""
from __future__ import annotations

import io
import itertools
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from bench_ab import ab, header  # noqa: E402

from zna.core import (  # noqa: E402
    COMPRESSION_ZSTD, FLAG_FIELDS, ZnaHeader, ZnaReader, ZnaWriter,
)

N_READS = 200_000


def make_file(read_len=150, n=N_READS, seed=7, block_size=None):
    rng = random.Random(seed)
    hdr = ZnaHeader(read_group="bench", seq_len_bytes=4,
                    compression_method=COMPRESSION_ZSTD, compression_level=3)
    buf = io.BytesIO()
    kw = {} if block_size is None else {"block_size": block_size}
    with ZnaWriter(buf, hdr, **kw) as w:
        for i in range(0, n - 1, 2):
            w.write_record("".join(rng.choices("ACGT", k=read_len)), True, True, False)
            w.write_record("".join(rng.choices("ACGT", k=read_len)), True, False, True)
    return buf.getvalue()


# --- consumers -------------------------------------------------------------

def collect_records(data):
    """Optimistic: a loader that only appends sequences."""
    return sum(len(r[0]) for r in ZnaReader(io.BytesIO(data)).records())


def collect_blocks(data):
    total = 0
    for seqs, _flags in ZnaReader(io.BytesIO(data)).blocks():
        for s in seqs:
            total += len(s)
    return total


MIN_LEN = 32


def realistic_records(data):
    """A per-record step in the shape of khorana's: length filter, flag-driven
    branch, and a slice of the sequence."""
    kept = 0
    for seq, is_paired, is_read1, _is_read2 in ZnaReader(io.BytesIO(data)).records():
        if len(seq) < MIN_LEN:
            continue
        has_sos = (not is_paired) or is_read1
        kept += len(seq[:300]) + has_sos
    return kept


def realistic_blocks(data):
    kept = 0
    fields = FLAG_FIELDS
    for seqs, flags in ZnaReader(io.BytesIO(data)).blocks():
        for seq, fl in zip(seqs, flags):
            if len(seq) < MIN_LEN:
                continue
            is_paired, is_read1, _is_read2 = fields[fl]
            has_sos = (not is_paired) or is_read1
            kept += len(seq[:300]) + has_sos
    return kept


def sharded_records(data, shard, n_shards):
    """khorana's pattern: stride over records. Decodes everything, keeps 1/N."""
    kept = 0
    stream = itertools.islice(
        ZnaReader(io.BytesIO(data)).records(), shard, None, n_shards)
    for seq, is_paired, is_read1, _is_read2 in stream:
        if len(seq) < MIN_LEN:
            continue
        kept += len(seq[:300]) + ((not is_paired) or is_read1)
    return kept


def sharded_blocks(data, shard, n_shards):
    """Stride over blocks: non-selected blocks are never decompressed."""
    kept = 0
    fields = FLAG_FIELDS
    for seqs, flags in ZnaReader(io.BytesIO(data)).blocks(
            stride=n_shards, offset=shard):
        for seq, fl in zip(seqs, flags):
            if len(seq) < MIN_LEN:
                continue
            is_paired, is_read1, _is_read2 = fields[fl]
            kept += len(seq[:300]) + ((not is_paired) or is_read1)
    return kept


def main():
    data = make_file()
    n_blocks = sum(1 for _ in ZnaReader(io.BytesIO(data)).blocks())
    print(f"\n  {N_READS} x 150 bp, {len(data) / 1e6:.1f} MB, {n_blocks} blocks")

    assert collect_records(data) == collect_blocks(data)
    assert realistic_records(data) == realistic_blocks(data)

    header("blocks() vs records(), same file, same process")
    ab("loader that only appends sequences", lambda: collect_records(data),
       lambda: collect_blocks(data))
    ab("loader with a realistic per-record step", lambda: realistic_records(data),
       lambda: realistic_blocks(data))

    # Block-level sharding needs many more blocks than shards, or workers get
    # wildly unequal shares (and some get none).  A real training file has
    # hundreds of blocks; this one is small, so shrink the block to match the
    # ratio rather than measure a degenerate case.
    sharded_data = make_file(block_size=192 * 1024)
    n_sharded_blocks = sum(1 for _ in ZnaReader(io.BytesIO(sharded_data)).blocks())
    header(f"sharded: stride over records vs blocks ({n_sharded_blocks} blocks)")
    data = sharded_data
    for n_shards in (2, 4, 8, 16):
        # One worker's share of the work, which is what wall-clock depends on.
        total_r = sum(sharded_records(data, s, n_shards) for s in range(n_shards))
        total_b = sum(sharded_blocks(data, s, n_shards) for s in range(n_shards))
        assert total_r == total_b == realistic_records(data), "shards lost records"
        ab(f"one worker of {n_shards} (whole file / this shard's blocks)",
           lambda ns=n_shards: sharded_records(data, 0, ns),
           lambda ns=n_shards: sharded_blocks(data, 0, ns))

    header("read-length decay (the plan flagged this)")
    for read_len, n in ((150, N_READS), (1000, 40_000), (10_000, 4_000)):
        d = make_file(read_len=read_len, n=n)
        assert realistic_records(d) == realistic_blocks(d)
        ab(f"{read_len} bp reads, realistic step",
           lambda dd=d: realistic_records(dd), lambda dd=d: realistic_blocks(dd))
    print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
