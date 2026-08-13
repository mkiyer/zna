"""Round-trip fuzz harness: arbitrary sequence must survive arbitrary codec
configuration bit-exactly.

Why this file exists
--------------------
The rest of the suite round-trips *fixed* inputs and checks API contracts.  It
does not assert that arbitrary sequence survives arbitrary codec configuration
bit-exactly — which is precisely the property a codec optimization threatens.
During the 0.3.4 performance sweep two proposed changes (a NEON packing kernel
and an encode/shuffle fusion) silently corrupted data **while all 282 tests
passed**.  A green suite was not evidence of a correct codec patch.

So: no codec change lands without this file green.

What it asserts
---------------
The matrix is sequence content × length × ``seq_len_bytes`` × N-policy × strand
configuration × compression × record layout × backend, and the assertions are:

1.  **Stored form is exact.**  Read back with ``with_rc=True`` and check the
    stored sequence equals the reference model's transform of the original,
    given the ``IS_RC`` the encoder actually chose.  The model is written here
    from the spec, independently of ``_pycodec``, so a bug in the Python backend
    cannot propagate into the expectation.
2.  **``restore_strand`` inverts exactly.**
3.  **Flags survive**, including the ``with_ends`` derivation.
4.  **Backend parity**: the Python and C++ backends produce byte-identical
    columnar streams from identical input, and identical records from identical
    streams.  This is the tightest check available for a codec rewrite — it
    needs no model at all.
5.  **Re-encode is a faithful copy** under ``preserve_normalization=True``.
6.  **Labels round-trip** per dtype.

Running deeper
--------------
Default settings keep this in the low seconds so it can gate every commit::

    python -m pytest tests/test_fuzz_roundtrip.py

For a long adversarial run (new seed each iteration)::

    python tests/test_fuzz_roundtrip.py --iters 200 --seed 12345

Environment overrides: ``ZNA_FUZZ_SEED``, ``ZNA_FUZZ_ITERS``.
"""
from __future__ import annotations

import contextlib
import gc
import io
import os
import random
import struct
import sys
import unittest

from zna.core import (
    COMPRESSION_NONE,
    COMPRESSION_ZSTD,
    FLAG_FIELDS,
    _BLOCK_HEADER_FMT,
    _BLOCK_HEADER_SIZE,
    ZnaHeader,
    ZnaReader,
    ZnaWriter,
    _ends_from_flags,
)
from zna.codec import available_backends, get_backend
from zna.dtypes import DTYPE_BY_CODE, LabelDef

# ---------------------------------------------------------------------------
# Knobs
# ---------------------------------------------------------------------------

SEED = int(os.environ.get("ZNA_FUZZ_SEED", "20260811"))
ITERS = int(os.environ.get("ZNA_FUZZ_ITERS", "1"))

#: Small enough that every fuzz file spans several blocks, so block boundaries,
#: the paired-R1 flush deferral, and multi-block decode are all on the hot path.
BLOCK_SIZE = 512

BACKENDS = tuple(available_backends())

#: ``(strand_specific, read1_antisense, read2_antisense, strand_normalized)``
STRAND_CONFIGS = (
    (False, False, False, False),  # no normalization at all
    (True, True, False, True),     # stranded, R1 antisense
    (True, False, True, True),     # stranded, R2 antisense
    (True, True, True, True),      # stranded, both antisense
    (False, False, False, True),   # unstranded -> random normalization
)

NPOLICIES = ("", "random")

LAYOUTS = ("se", "pe", "merged", "mixed")

#: Lengths for the main matrix: zero, sub-byte, every residue class mod 4, the
#: packing boundary, and the 1-byte width's own maximum.  The multi-byte width
#: maxima are large enough that generating them 24-per-config would dominate the
#: runtime, so they live in :meth:`TestRoundTripMatrix.test_large_lengths`
#: instead, where a handful of records covers them just as well.
_LENGTHS = {
    1: (0, 1, 2, 3, 4, 5, 7, 8, 13, 150, 151, 254, 255),
    2: (0, 1, 3, 4, 150, 255, 256, 1000),
    4: (0, 1, 3, 4, 150, 1000),
}

#: Lengths that exercise each width's boundary, checked separately.
_LARGE_LENGTHS = {
    1: (253, 254, 255),
    2: (65533, 65534, 65535),
    4: (65535, 65536, 70001),
}


# ---------------------------------------------------------------------------
# Reference model — written from the spec, deliberately not reusing _pycodec
# ---------------------------------------------------------------------------

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def model_revcomp(seq: str) -> str:
    """Reverse-complement.  Characters outside ACGTacgt reverse uncomplemented
    (both backends agree on this: N stays N, and position is mirrored)."""
    return seq.translate(_COMPLEMENT)[::-1]


def model_stored(orig: str, is_rc: bool, npolicy: str) -> tuple[str, list[int]]:
    """The sequence the file should hold, given the ``IS_RC`` the encoder chose.

    Returns ``(expected, wildcard_positions)``.  Under ``npolicy="random"`` the
    substituted base comes from a seeded, position-derived stream the model does not
    reproduce, so those positions are returned as wildcards and checked structurally
    (must be one of ACGT) instead.  The *positions* are fully determined, which is what
    this model is for.

    Order matters and mirrors both backends: reverse-complement first, N-policy
    second.  Decoding always emits uppercase, so the model uppercases last.
    """
    s = model_revcomp(orig) if is_rc else orig
    if not npolicy or ("N" not in s and "n" not in s):
        return s.upper(), []
    # Random replacement: the bases are unpredictable, so hand back the N
    # positions as wildcards for the caller to check structurally.
    wildcards = [i for i, c in enumerate(s) if c in ("N", "n")]
    return s.replace("N", "A").replace("n", "A").upper(), wildcards


# ---------------------------------------------------------------------------
# Generators
# ---------------------------------------------------------------------------

def gen_seq(rng: random.Random, length: int, *, allow_n: bool, allow_lower: bool) -> str:
    alphabet = "ACGT"
    if allow_lower:
        alphabet += "acgt"
    if allow_n:
        # Enough N to matter without swamping the sequence.
        alphabet += "N" * 2 + ("n" if allow_lower else "")
    if not length:
        return ""
    return "".join(rng.choices(alphabet, k=length))


def gen_records(rng: random.Random, layout: str, n: int, seq_len_bytes: int,
                *, allow_n: bool, allow_lower: bool) -> list[tuple]:
    """Build ``(seq, is_paired, is_read1, is_read2, is_full_fragment)`` records.

    Paired records are emitted strictly as adjacent R1/R2, which is what the
    unstranded normalizer requires to pick exactly one mate per fragment.
    """
    lengths = _LENGTHS[seq_len_bytes]
    recs: list[tuple] = []

    def one(is_paired, is_read1, is_read2, is_full=False):
        seq = gen_seq(rng, rng.choice(lengths), allow_n=allow_n, allow_lower=allow_lower)
        return (seq, is_paired, is_read1, is_read2, is_full)

    while len(recs) < n:
        if layout == "se":
            recs.append(one(False, True, False))
        elif layout == "pe":
            recs.append(one(True, True, False))
            recs.append(one(True, False, True))
        elif layout == "merged":
            # A merged read carries neither R1 nor R2 and spans its fragment.
            recs.append(one(False, False, False, is_full=True))
        else:  # mixed
            pick = rng.randrange(3)
            if pick == 0:
                recs.append(one(False, True, False))
            elif pick == 1:
                recs.append(one(True, True, False))
                recs.append(one(True, False, True))
            else:
                recs.append(one(False, False, False, is_full=rng.random() < 0.5))
    return recs[:n]


_LABEL_RANGES = {
    "A": (0, 255), "C": (0, 255), "c": (-128, 127),
    "S": (0, 65535), "s": (-32768, 32767),
    "I": (0, 2**32 - 1), "i": (-(2**31), 2**31 - 1),
    "Q": (0, 2**64 - 1), "q": (-(2**63), 2**63 - 1),
}


def gen_label_value(rng: random.Random, code: str):
    if code in _LABEL_RANGES:
        lo, hi = _LABEL_RANGES[code]
        return rng.randint(lo, hi)
    if code == "f":
        # Round-trip through float32 so the model matches what disk can hold.
        v = rng.uniform(-1e6, 1e6)
        return struct.unpack("<f", struct.pack("<f", v))[0]
    if code == "d":
        return rng.uniform(-1e12, 1e12)
    raise AssertionError(f"unhandled dtype code {code!r}")


# ---------------------------------------------------------------------------
# Harness
# ---------------------------------------------------------------------------

@contextlib.contextmanager
def force_backend(name: str):
    """Force the codec backend.

    ``zna.core`` resolves its backend at *import* time, so patching ``zna.codec``
    afterwards has no effect — the module globals are what must be swapped.
    """
    import zna.core as core
    saved_codec, saved_accel = core._codec, core._accel_mod
    try:
        core._codec = get_backend(name)
        if name == "accel":
            from zna import _accel as accel_mod
            core._accel_mod = accel_mod
        else:
            core._accel_mod = None
        yield
    finally:
        core._codec, core._accel_mod = saved_codec, saved_accel


def make_header(seq_len_bytes, strand_cfg, compression, labels=()):
    ss, r1a, r2a, sn = strand_cfg
    return ZnaHeader(
        read_group="fuzz",
        seq_len_bytes=seq_len_bytes,
        strand_specific=ss,
        read1_antisense=r1a,
        read2_antisense=r2a,
        strand_normalized=sn,
        compression_method=compression,
        compression_level=3,
        labels=labels,
    )


def write_file(header, recs, npolicy, labels_per_rec=None) -> bytes:
    buf = io.BytesIO()
    with ZnaWriter(buf, header, block_size=BLOCK_SIZE, npolicy=npolicy) as w:
        for i, (seq, is_paired, is_read1, is_read2, is_full) in enumerate(recs):
            w.write_record(
                seq, is_paired, is_read1, is_read2,
                labels=None if labels_per_rec is None else labels_per_rec[i],
                is_full_fragment=is_full,
            )
    return buf.getvalue()


class FuzzCase(unittest.TestCase):
    """Shared assertions.  Failures print the config and the offending record so
    a fuzz failure is directly reproducible rather than just "something broke"."""

    def assert_stored(self, ctx, idx, orig, got, is_rc, npolicy):
        expected, wildcards = model_stored(orig, is_rc, npolicy)
        self.assertEqual(
            len(got), len(expected),
            f"{ctx} rec{idx}: length {len(got)} != {len(expected)}",
        )
        if not wildcards:
            self.assertEqual(
                got, expected,
                f"{ctx} rec{idx}: stored sequence mismatch\n"
                f"  orig={orig!r}\n  is_rc={is_rc}\n  got={got!r}\n  exp={expected!r}",
            )
            return
        wild = set(wildcards)
        for i, (g, e) in enumerate(zip(got, expected)):
            if i in wild:
                self.assertIn(
                    g, "ACGT",
                    f"{ctx} rec{idx} pos{i}: random N-policy emitted {g!r}",
                )
            else:
                self.assertEqual(
                    g, e,
                    f"{ctx} rec{idx} pos{i}: {g!r} != {e!r} (orig={orig!r}, is_rc={is_rc})",
                )


class TestRoundTripMatrix(FuzzCase):
    """The main sweep: every configuration combination, random data in each."""

    def _run_config(self, rng, backend, seq_len_bytes, npolicy, strand_cfg,
                    compression, layout, n_records=24):
        ctx = (f"backend={backend} slb={seq_len_bytes} npolicy={npolicy!r} "
               f"strand={strand_cfg} comp={compression} layout={layout}")
        allow_n = bool(npolicy)
        recs = gen_records(rng, layout, n_records, seq_len_bytes,
                           allow_n=allow_n, allow_lower=True)
        header = make_header(seq_len_bytes, strand_cfg, compression)

        with force_backend(backend):
            data = write_file(header, recs, npolicy)

            # -- leg 1: stored form, with the encoder's own IS_RC ------------
            reader = ZnaReader(io.BytesIO(data))
            out = list(reader.records(with_rc=True))
            self.assertEqual(len(out), len(recs), f"{ctx}: record count")

            rc_flags = []
            for i, (rec, exp) in enumerate(zip(out, recs)):
                seq, is_paired, is_read1, is_read2, is_rc = rec
                self.assert_stored(ctx, i, exp[0], seq, is_rc, npolicy)
                self.assertEqual((is_paired, is_read1, is_read2), exp[1:4],
                                 f"{ctx} rec{i}: flags")
                rc_flags.append(is_rc)

            # An un-normalized file must never set IS_RC.
            if not strand_cfg[3]:
                self.assertFalse(any(rc_flags), f"{ctx}: IS_RC set without normalization")

            # -- leg 2: restore_strand inverts exactly -----------------------
            if strand_cfg[3]:
                reader = ZnaReader(io.BytesIO(data))
                restored = list(reader.records(restore_strand=True))
                self.assertEqual(len(restored), len(recs), f"{ctx}: restored count")
                for i, (rec, is_rc) in enumerate(zip(restored, rc_flags)):
                    stored, _ = model_stored(recs[i][0], is_rc, npolicy)
                    expected = model_revcomp(stored) if is_rc else stored
                    if npolicy == "random":
                        self.assertEqual(len(rec[0]), len(expected), f"{ctx} rec{i}")
                    else:
                        self.assertEqual(
                            rec[0], expected,
                            f"{ctx} rec{i}: restore_strand did not invert\n"
                            f"  orig={recs[i][0]!r} is_rc={is_rc}",
                        )

            # -- leg 3: with_ends derivation ---------------------------------
            reader = ZnaReader(io.BytesIO(data))
            ends_out = list(reader.records(with_ends=True))
            for i, (rec, is_rc) in enumerate(zip(ends_out, rc_flags)):
                has_start, has_end = rec[4], rec[5]
                self.assertEqual(
                    (has_start, has_end), _ends_from_flags(is_rc, recs[i][4]),
                    f"{ctx} rec{i}: with_ends mismatch (is_rc={is_rc}, "
                    f"is_full={recs[i][4]})",
                )

            # -- leg 4: default 4-tuple width is unchanged -------------------
            reader = ZnaReader(io.BytesIO(data))
            plain = list(reader.records())
            self.assertTrue(all(len(r) == 4 for r in plain),
                            f"{ctx}: default records() width must stay 4")

    def test_full_matrix(self):
        """Exhaustive over every config axis; random sequence content within."""
        rng = random.Random(SEED)
        n = 0
        for it in range(ITERS):
            for backend in BACKENDS:
                for seq_len_bytes in (1, 2, 4):
                    for npolicy in NPOLICIES:
                        for strand_cfg in STRAND_CONFIGS:
                            for compression in (COMPRESSION_NONE, COMPRESSION_ZSTD):
                                for layout in LAYOUTS:
                                    self._run_config(
                                        rng, backend, seq_len_bytes, npolicy,
                                        strand_cfg, compression, layout,
                                    )
                                    n += 1
        self.assertGreater(n, 0)

    def test_edge_lengths_exhaustive(self):
        """Every length 0..64 plus the width boundaries, one record per file so
        no neighbour can mask a tail bug."""
        rng = random.Random(SEED + 1)
        for backend in BACKENDS:
            for seq_len_bytes in (1, 2, 4):
                lengths = list(range(0, 65)) + list(_LARGE_LENGTHS[seq_len_bytes])
                for length in lengths:
                    seq = gen_seq(rng, length, allow_n=False, allow_lower=False)
                    header = make_header(seq_len_bytes, STRAND_CONFIGS[0],
                                         COMPRESSION_ZSTD)
                    with force_backend(backend):
                        data = write_file(header, [(seq, False, True, False, False)], "")
                        got = list(ZnaReader(io.BytesIO(data)).records())
                    self.assertEqual(len(got), 1)
                    self.assertEqual(
                        got[0][0], seq,
                        f"backend={backend} slb={seq_len_bytes} len={length}: "
                        f"tail corruption",
                    )

    def test_large_lengths(self):
        """The multi-byte length widths at their boundaries, under strand
        normalization and an N-policy so the reverse-complement and the packing
        tail are both exercised on a long record."""
        rng = random.Random(SEED + 9)
        for backend in BACKENDS:
            for seq_len_bytes in (1, 2, 4):
                for strand_cfg in (STRAND_CONFIGS[1], STRAND_CONFIGS[4]):
                    for npolicy in ("", "random"):
                        lengths = _LARGE_LENGTHS[seq_len_bytes]
                        recs = [
                            (gen_seq(rng, ln, allow_n=bool(npolicy), allow_lower=True),
                             False, True, False, False)
                            for ln in lengths
                        ]
                        header = make_header(seq_len_bytes, strand_cfg,
                                             COMPRESSION_ZSTD)
                        ctx = (f"backend={backend} slb={seq_len_bytes} "
                               f"strand={strand_cfg} npolicy={npolicy!r}")
                        with force_backend(backend):
                            data = write_file(header, recs, npolicy)
                            got = list(
                                ZnaReader(io.BytesIO(data)).records(with_rc=True)
                            )
                        self.assertEqual(len(got), len(recs), f"{ctx}: count")
                        for i, rec in enumerate(got):
                            self.assert_stored(ctx, i, recs[i][0], rec[0], rec[4],
                                               npolicy)


class TestBackendParity(FuzzCase):
    """Python and C++ must be bit-identical on deterministic configurations.

    This needs no reference model, so it stays valid however the codec is
    rewritten — and it is the check most likely to catch a fast-path kernel that
    is subtly wrong (the NEON packing bug would have failed here immediately).
    """

    def setUp(self):
        if "accel" not in BACKENDS:
            self.skipTest("C++ accel backend not available")

    def test_encode_streams_identical(self):
        from zna import _accel, _pycodec
        rng = random.Random(SEED + 2)
        for seq_len_bytes in (1, 2, 4):
            for npolicy in ("", "random"):   # `random` is seeded, so comparable
                for do_rc_r1, do_rc_r2 in ((False, False), (True, False),
                                           (False, True), (True, True)):
                    recs = gen_records(rng, "mixed", 40, seq_len_bytes,
                                       allow_n=bool(npolicy), allow_lower=True)
                    seqs = [r[0] for r in recs]
                    flags = [
                        (1 if r[2] else 0) | (2 if r[3] else 0) | (4 if r[1] else 0)
                        | (16 if r[4] else 0)
                        for r in recs
                    ]
                    py = _pycodec.encode_block(seqs, list(flags), seq_len_bytes,
                                               npolicy, do_rc_r1, do_rc_r2, False)
                    cc = _accel.encode_block(seqs, list(flags), seq_len_bytes,
                                             npolicy, do_rc_r1, do_rc_r2, False)
                    ctx = (f"slb={seq_len_bytes} npolicy={npolicy!r} "
                           f"rc1={do_rc_r1} rc2={do_rc_r2}")
                    self.assertEqual(py[0], cc[0], f"{ctx}: flags stream differs")
                    self.assertEqual(py[1], cc[1], f"{ctx}: lengths stream differs")
                    self.assertEqual(py[2], cc[2], f"{ctx}: sequence stream differs")

    def test_decode_records_identical(self):
        from zna import _accel, _pycodec
        rng = random.Random(SEED + 3)
        for seq_len_bytes in (1, 2, 4):
            for do_rc_r1 in (False, True):
                recs = gen_records(rng, "mixed", 60, seq_len_bytes,
                                   allow_n=False, allow_lower=False)
                seqs = [r[0] for r in recs]
                flags = [
                    (1 if r[2] else 0) | (2 if r[3] else 0) | (4 if r[1] else 0)
                    for r in recs
                ]
                fl, ln, sq = _pycodec.encode_block(seqs, list(flags), seq_len_bytes,
                                                   "", do_rc_r1, False, False)
                py = _pycodec.decode_block(fl, ln, sq, seq_len_bytes, len(seqs))
                cc = _accel.decode_block(fl, ln, sq, seq_len_bytes, len(seqs))
                ctx = f"slb={seq_len_bytes} rc1={do_rc_r1}"
                self.assertEqual(len(py), len(cc), f"{ctx}: record count differs")
                for i, (a, b) in enumerate(zip(py, cc)):
                    self.assertEqual(
                        tuple(a), tuple(b),
                        f"{ctx} rec{i}: decoded record differs\n  py={a!r}\n  cc={b!r}",
                    )

    def test_reverse_complement_identical(self):
        """Bases and their order must match exactly.

        Case is compared separately: the two backends genuinely differ on it
        (see :meth:`test_reverse_complement_case_is_backend_specific`), and that
        difference is invisible through a file because decoding always emits
        uppercase.  Comparing case-folded here keeps this test fully sensitive
        to the thing that would be corruption — a wrong base or a wrong
        position — without tripping over the known divergence.
        """
        from zna import _accel, _pycodec
        rng = random.Random(SEED + 4)
        for length in list(range(0, 40)) + [149, 150, 151, 1000]:
            for allow_n in (False, True):
                seq = gen_seq(rng, length, allow_n=allow_n, allow_lower=True)
                py = _pycodec.reverse_complement(seq)
                cc = _accel.reverse_complement(seq)
                self.assertEqual(
                    py.upper(), cc.upper(),
                    f"reverse_complement differs at len={length} seq={seq!r}",
                )
                self.assertEqual(
                    cc.upper(), model_revcomp(seq).upper(),
                    f"reverse_complement disagrees with model at seq={seq!r}",
                )

    def test_reverse_complement_case_is_backend_specific(self):
        """Pin the known case divergence so neither backend drifts unnoticed.

        ``_pycodec`` preserves case (it is a ``str.translate``); ``_accel``
        uppercases (it round-trips through the 2-bit encode table).  Nothing in
        the file format depends on this — decode always emits uppercase — but
        ``reverse_complement`` is re-exported from ``zna.core``, so it is
        observable to callers.
        """
        from zna import _accel, _pycodec
        self.assertEqual(_pycodec.reverse_complement("acgt"), "acgt")
        self.assertEqual(_pycodec.reverse_complement("aCgT"), "AcGt")
        self.assertEqual(_accel.reverse_complement("acgt"), "ACGT")
        self.assertEqual(_accel.reverse_complement("aCgT"), "ACGT")
        # Both leave non-ACGT alone, complementing nothing and only mirroring.
        self.assertEqual(_pycodec.reverse_complement("acgtN"), "Nacgt")
        self.assertEqual(_accel.reverse_complement("acgtN"), "NACGT")

    def test_encode_sequence_identical(self):
        from zna import _accel, _pycodec
        rng = random.Random(SEED + 5)
        for length in list(range(0, 40)) + [149, 150, 151, 1000]:
            seq = gen_seq(rng, length, allow_n=False, allow_lower=True)
            self.assertEqual(
                _pycodec.encode_sequence(seq), _accel.encode_sequence(seq),
                f"encode_sequence differs at len={length} seq={seq!r}",
            )


class TestReencode(FuzzCase):
    """A ZNA -> ZNA copy under ``preserve_normalization=True`` must be a faithful
    copy, and must stay faithful when repeated.  Applying orientation twice was a
    real defect: it un-normalized the data while the header still claimed
    otherwise, and nothing downstream could detect it."""

    def test_reencode_is_idempotent(self):
        rng = random.Random(SEED + 6)
        for backend in BACKENDS:
            for strand_cfg in STRAND_CONFIGS:
                for layout in ("pe", "mixed"):
                    recs = gen_records(rng, layout, 30, 2,
                                       allow_n=False, allow_lower=False)
                    header = make_header(2, strand_cfg, COMPRESSION_ZSTD)
                    ctx = f"backend={backend} strand={strand_cfg} layout={layout}"
                    with force_backend(backend):
                        data = write_file(header, recs, "")

                        # Copy twice, carrying the FLAG BYTE.
                        #
                        # This used to copy through `records(with_ends=True)` and compare
                        # the same view, which made it blind to the defect it exists to
                        # catch: (has_start, has_end) has three states where
                        # (IS_RC, IS_FULL_FRAGMENT) has four, so byte 24 decayed to 16
                        # on generation 2 — and the comparison, being in the same lossy
                        # view, saw two records that agreed. Comparing ZnaRecord
                        # compares the whole byte, so a dropped bit fails here now.
                        prev = list(ZnaReader(io.BytesIO(data)).copy_records())
                        for gen_i in (2, 3):
                            buf = io.BytesIO()
                            with ZnaWriter(buf, header, block_size=BLOCK_SIZE,
                                           preserve_normalization=True) as w:
                                for rec in prev:
                                    w.write_copy(rec)
                            nxt = list(
                                ZnaReader(io.BytesIO(buf.getvalue())).copy_records()
                            )
                            self.assertEqual(
                                nxt, prev,
                                f"{ctx}: re-encode generation {gen_i} diverged",
                            )
                            prev = nxt


class TestLabels(FuzzCase):
    """Label columns share the block payload with the sequence stream, so a
    decoder rewrite can shift them.  Every dtype must round-trip exactly."""

    def test_all_dtypes_roundtrip(self):
        rng = random.Random(SEED + 7)
        codes = sorted(DTYPE_BY_CODE)
        for backend in BACKENDS:
            for compression in (COMPRESSION_NONE, COMPRESSION_ZSTD):
                for seq_len_bytes in (1, 2):
                    labels = tuple(
                        LabelDef(label_id=i, name=f"L{c}", description=f"dtype {c}",
                                 dtype=DTYPE_BY_CODE[c])
                        for i, c in enumerate(codes)
                    )
                    recs = gen_records(rng, "mixed", 25, seq_len_bytes,
                                       allow_n=False, allow_lower=False)
                    values = [
                        tuple(gen_label_value(rng, c) for c in codes)
                        for _ in recs
                    ]
                    header = make_header(seq_len_bytes, STRAND_CONFIGS[0],
                                         compression, labels=labels)
                    ctx = (f"backend={backend} comp={compression} "
                           f"slb={seq_len_bytes}")
                    with force_backend(backend):
                        data = write_file(header, recs, "", labels_per_rec=values)
                        got = list(ZnaReader(io.BytesIO(data)).records())
                    self.assertEqual(len(got), len(recs), f"{ctx}: count")
                    for i, rec in enumerate(got):
                        self.assertEqual(len(rec), 5,
                                         f"{ctx}: labeled records() width must stay 5")
                        self.assertEqual(rec[0], recs[i][0], f"{ctx} rec{i}: sequence")
                        self.assertEqual(
                            tuple(rec[4]), values[i],
                            f"{ctx} rec{i}: labels\n  got={rec[4]!r}\n"
                            f"  exp={values[i]!r}",
                        )


class TestNPolicyErrors(FuzzCase):
    """With no N-policy, an N is an error on both backends — an optimization
    that silently packs it as A would be data corruption."""

    def test_n_without_policy_raises(self):
        rng = random.Random(SEED + 8)
        for backend in BACKENDS:
            for pos in (0, 1, 2, 3, 4, 7, 149):
                length = max(pos + 1, 150)
                seq = list(gen_seq(rng, length, allow_n=False, allow_lower=False))
                seq[pos] = "N"
                seq = "".join(seq)
                header = make_header(1, STRAND_CONFIGS[0], COMPRESSION_NONE)
                with force_backend(backend):
                    with self.assertRaises(
                        ValueError,
                        msg=f"backend={backend} pos={pos}: N accepted without policy",
                    ):
                        write_file(header, [(seq, False, True, False, False)], "")


class TestInvalidBaseHandling(FuzzCase):
    """Regressions for the first bug this harness found.

    ``_pycodec.encode_sequence`` validated a 4-base group by testing
    ``val > 255`` *after* OR-ing the four 2-bit values together.  An invalid
    base in the group's last slot ORs to exactly 255 and never exceeds it, so
    the group was packed as 0xFF and decoded back as "TTTT" — silently losing
    three valid bases with it.  Every index congruent to 3 mod 4 was affected.
    """

    def test_invalid_base_raises_at_every_offset(self):
        from zna import _accel, _pycodec
        base = "ACGTACGTACGTACGT"
        for pos in range(len(base)):
            for ch in "NRYSWKMBDHVXn":
                seq = base[:pos] + ch + base[pos + 1:]
                for name, mod in (("python", _pycodec), ("accel", _accel)):
                    with self.assertRaises(
                        ValueError,
                        msg=f"{name}: {ch!r} at index {pos} was accepted, "
                            f"not rejected",
                    ):
                        mod.encode_sequence(seq)

    def test_invalid_base_never_decodes_as_bases(self):
        """The specific corruption: "ACGN" must not become "TTTT"."""
        from zna import _pycodec
        with self.assertRaises(ValueError):
            _pycodec.encode_sequence("ACGN")

    def test_both_backends_reject_non_n_ambiguity_codes(self):
        """`--npolicy` governs `N` and nothing else, in BOTH backends.

        This used to document a live divergence: with a policy set, ``_accel``
        substituted the policy base for *any* unencodable character while ``_pycodec``
        substituted only ``N``/``n`` and raised on the rest — so an IUPAC code encoded on
        the compiled backend and raised on the reference one. Both now raise.
        """
        from zna import _pycodec
        mods = [_pycodec]
        try:
            from zna import _accel
            mods.append(_accel)
        except ImportError:
            pass
        for mod in mods:
            for policy in ("", "random"):
                for ch in "RYSWKMBDHV.-":
                    with self.subTest(mod=mod.__name__, policy=policy, ch=ch):
                        with self.assertRaises(ValueError):
                            mod.encode_block([f"AC{ch}GT"], [0], 2, policy,
                                             False, False, False)

    def test_blocks_match_records(self):
        rng = random.Random(SEED + 12)
        for backend in BACKENDS:
            for seq_len_bytes in (1, 2, 4):
                for strand_cfg in (STRAND_CONFIGS[0], STRAND_CONFIGS[1],
                                   STRAND_CONFIGS[4]):
                    for compression in (COMPRESSION_NONE, COMPRESSION_ZSTD):
                        recs = gen_records(rng, "mixed", 60, seq_len_bytes,
                                           allow_n=False, allow_lower=False)
                        header = make_header(seq_len_bytes, strand_cfg, compression)
                        ctx = (f"backend={backend} slb={seq_len_bytes} "
                               f"strand={strand_cfg} comp={compression}")
                        with force_backend(backend):
                            data = write_file(header, recs, "")
                            want = list(ZnaReader(io.BytesIO(data)).records())
                            got = []
                            for seqs, flags in ZnaReader(io.BytesIO(data)).blocks():
                                self.assertEqual(
                                    len(seqs), len(flags),
                                    f"{ctx}: sequences and flags out of step")
                                for seq, fl in zip(seqs, flags):
                                    got.append((seq, *FLAG_FIELDS[fl]))
                        self.assertEqual(got, want, f"{ctx}: blocks() != records()")

    def test_blocks_restore_strand_matches_records(self):
        rng = random.Random(SEED + 13)
        for backend in BACKENDS:
            for strand_cfg in (STRAND_CONFIGS[1], STRAND_CONFIGS[2],
                               STRAND_CONFIGS[4]):
                recs = gen_records(rng, "mixed", 50, 2,
                                   allow_n=False, allow_lower=False)
                header = make_header(2, strand_cfg, COMPRESSION_ZSTD)
                ctx = f"backend={backend} strand={strand_cfg}"
                with force_backend(backend):
                    data = write_file(header, recs, "")
                    want = [r[0] for r in
                            ZnaReader(io.BytesIO(data)).records(restore_strand=True)]
                    got = []
                    for seqs, _flags in ZnaReader(io.BytesIO(data)).blocks(
                            restore_strand=True):
                        got.extend(seqs)
                self.assertEqual(got, want, f"{ctx}: restore_strand differs")

    def test_sharding_partitions_exactly(self):
        """Every record lands in exactly one shard, for every stride."""
        rng = random.Random(SEED + 14)
        # A small block size so the file has many blocks to spread over shards.
        recs = gen_records(rng, "mixed", 400, 2, allow_n=False, allow_lower=False)
        header = make_header(2, STRAND_CONFIGS[0], COMPRESSION_ZSTD)
        for backend in BACKENDS:
            with force_backend(backend):
                data = write_file(header, recs, "")
                want = [r[0] for r in ZnaReader(io.BytesIO(data)).records()]
                n_blocks = sum(1 for _ in ZnaReader(io.BytesIO(data)).blocks())
                for stride in (1, 2, 3, 5, 8):
                    shards = []
                    for offset in range(stride):
                        got = []
                        for seqs, _f in ZnaReader(io.BytesIO(data)).blocks(
                                stride=stride, offset=offset):
                            got.extend(seqs)
                        shards.append(got)
                    merged = [s for shard in shards for s in shard]
                    self.assertEqual(
                        sorted(merged), sorted(want),
                        f"backend={backend} stride={stride}: shards do not "
                        f"partition the file ({n_blocks} blocks)",
                    )
                    self.assertEqual(
                        len(merged), len(want),
                        f"backend={backend} stride={stride}: record count changed",
                    )

    def test_sharding_works_on_a_non_seekable_stream(self):
        """Block skipping seeks when it can and reads through when it cannot."""
        class Unseekable(io.RawIOBase):
            def __init__(self, data):
                self._b = io.BytesIO(data)
            def readable(self):
                return True
            def read(self, n=-1):
                return self._b.read(n)
            def readinto(self, b):
                return self._b.readinto(b)
            def seek(self, *a):
                raise OSError("not seekable")
            def seekable(self):
                return False

        rng = random.Random(SEED + 15)
        recs = gen_records(rng, "mixed", 300, 2, allow_n=False, allow_lower=False)
        header = make_header(2, STRAND_CONFIGS[0], COMPRESSION_ZSTD)
        data = write_file(header, recs, "")
        want = [r[0] for r in ZnaReader(io.BytesIO(data)).records()]
        merged = []
        for offset in range(3):
            for seqs, _f in ZnaReader(Unseekable(data)).blocks(stride=3, offset=offset):
                merged.extend(seqs)
        self.assertEqual(sorted(merged), sorted(want),
                         "sharding over a non-seekable stream lost records")

    def test_rejects_labeled_files_and_bad_shards(self):
        labels = (LabelDef(label_id=0, name="NH", description="",
                           dtype=DTYPE_BY_CODE["C"]),)
        header = make_header(2, STRAND_CONFIGS[0], COMPRESSION_NONE, labels=labels)
        data = write_file(header, [("ACGT", False, True, False, False)], "",
                          labels_per_rec=[(1,)])
        with self.assertRaises(NotImplementedError):
            ZnaReader(io.BytesIO(data)).blocks()

        plain = write_file(make_header(2, STRAND_CONFIGS[0], COMPRESSION_NONE),
                           [("ACGT", False, True, False, False)], "")
        for stride, offset in ((0, 0), (-1, 0), (2, 2), (2, -1), (2, 5)):
            with self.assertRaises(ValueError, msg=f"stride={stride} offset={offset}"):
                ZnaReader(io.BytesIO(plain)).blocks(stride=stride, offset=offset)


class TestBlocksWithLabels(FuzzCase):
    """``blocks()`` on a labeled file.

    The risk here is not the API, it is the payload split.  A block is
    ``flags | labels | lengths | sequences``, and ``blocks()`` previously refused
    labeled files, so its splitting code had never had to account for the label
    bytes.  Getting that offset wrong decodes sequences from the wrong position
    and produces plausible garbage rather than an error — so every one of these
    checks it against ``records()``, which has always split correctly.
    """

    CODES = ("i", "i", "C")  # hulkrna's schema: ZI, ZJ, ZF

    def _labeled_file(self, rng, n=200, seq_len_bytes=2, strand_cfg=None,
                      compression=COMPRESSION_ZSTD, codes=None,
                      block_size=BLOCK_SIZE):
        codes = codes or self.CODES
        labels = tuple(
            LabelDef(label_id=i, name=f"L{i}", description="",
                     dtype=DTYPE_BY_CODE[c])
            for i, c in enumerate(codes)
        )
        recs = gen_records(rng, "mixed", n, seq_len_bytes,
                           allow_n=False, allow_lower=False)
        values = [tuple(gen_label_value(rng, c) for c in codes) for _ in recs]
        header = make_header(seq_len_bytes, strand_cfg or STRAND_CONFIGS[0],
                             compression, labels=labels)
        buf = io.BytesIO()
        with ZnaWriter(buf, header, block_size=block_size) as w:
            for (seq, is_paired, is_r1, is_r2, is_full), vals in zip(recs, values):
                w.write_record(seq, is_paired, is_r1, is_r2, labels=vals,
                               is_full_fragment=is_full)
        return buf.getvalue(), recs, values

    def test_default_still_refuses_labeled_files(self):
        rng = random.Random(SEED + 20)
        data, _, _ = self._labeled_file(rng)
        with self.assertRaises(NotImplementedError) as cm:
            ZnaReader(io.BytesIO(data)).blocks()
        self.assertIn("labels=", str(cm.exception),
                      "the refusal should name the way out")

    def test_labels_false_matches_records(self):
        """The split must be right even when the values are not wanted."""
        for backend in BACKENDS:
            for seq_len_bytes in (1, 2, 4):
                for compression in (COMPRESSION_NONE, COMPRESSION_ZSTD):
                    for codes in (("i", "i", "C"), ("C",), ("d", "q", "f", "A")):
                        rng = random.Random(SEED + 21)
                        data, recs, _ = self._labeled_file(
                            rng, seq_len_bytes=seq_len_bytes,
                            compression=compression, codes=codes)
                        ctx = (f"backend={backend} slb={seq_len_bytes} "
                               f"comp={compression} codes={codes}")
                        with force_backend(backend):
                            want = [(r[0], *r[1:4]) for r in
                                    ZnaReader(io.BytesIO(data)).records()]
                            got = []
                            for seqs, flags in ZnaReader(io.BytesIO(data)).blocks(
                                    labels=False):
                                self.assertEqual(len(seqs), len(flags), ctx)
                                for seq, fl in zip(seqs, flags):
                                    got.append((seq, *FLAG_FIELDS[fl]))
                        self.assertEqual(got, want, f"{ctx}: labels=False diverged")

    def test_labels_true_returns_columns(self):
        for backend in BACKENDS:
            for compression in (COMPRESSION_NONE, COMPRESSION_ZSTD):
                rng = random.Random(SEED + 22)
                data, recs, values = self._labeled_file(
                    rng, compression=compression)
                ctx = f"backend={backend} comp={compression}"
                got_seqs, got_vals = [], []
                with force_backend(backend):
                    for seqs, flags, columns in ZnaReader(io.BytesIO(data)).blocks(
                            labels=True):
                        self.assertEqual(len(columns), len(self.CODES),
                                         f"{ctx}: wrong column count")
                        for col in columns:
                            self.assertEqual(len(col), len(seqs),
                                             f"{ctx}: column length != record count")
                        got_seqs.extend(seqs)
                        got_vals.extend(zip(*columns))
                self.assertEqual(got_seqs, [r[0] for r in recs], f"{ctx}: sequences")
                self.assertEqual(got_vals, values, f"{ctx}: label values")

    def test_labels_true_agrees_with_records(self):
        """Columnar labels must equal what records() yields per record."""
        rng = random.Random(SEED + 23)
        data, _, _ = self._labeled_file(rng, n=300)
        per_record = [tuple(r[4]) for r in ZnaReader(io.BytesIO(data)).records()]
        columnar = []
        for _seqs, _flags, columns in ZnaReader(io.BytesIO(data)).blocks(labels=True):
            columnar.extend(zip(*columns))
        self.assertEqual(columnar, per_record)

    def test_labels_true_on_unlabeled_file_yields_no_columns(self):
        """len(label_columns) == len(header.labels) is the invariant, so an
        unlabeled file gives an empty tuple rather than an error."""
        rng = random.Random(SEED + 24)
        recs = gen_records(rng, "mixed", 50, 2, allow_n=False, allow_lower=False)
        header = make_header(2, STRAND_CONFIGS[0], COMPRESSION_ZSTD)
        data = write_file(header, recs, "")
        n = 0
        for seqs, flags, columns in ZnaReader(io.BytesIO(data)).blocks(labels=True):
            self.assertEqual(columns, ())
            n += len(seqs)
        self.assertEqual(n, len(recs))

    def test_labeled_sharding_and_restore_strand(self):
        """The label offset has to hold under every other blocks() feature."""
        rng = random.Random(SEED + 25)
        data, recs, values = self._labeled_file(
            rng, n=400, strand_cfg=STRAND_CONFIGS[1], block_size=256)
        index = ZnaReader(io.BytesIO(data)).block_index()
        self.assertGreater(len(index), 3, "need several blocks")

        # Sharded, labels skipped: shards must partition the file exactly.
        want = [r[0] for r in ZnaReader(io.BytesIO(data)).records()]
        merged = []
        for off in range(3):
            for seqs, _f in ZnaReader(io.BytesIO(data)).blocks(
                    stride=3, offset=off, labels=False):
                merged.extend(seqs)
        self.assertEqual(sorted(merged), sorted(want), "sharding lost records")

        # indices=, with labels: values must stay aligned to their sequences.
        want_pairs = [(r[0], tuple(r[4]))
                      for r in ZnaReader(io.BytesIO(data)).records()]
        by_block = []
        for seqs, _f, cols in ZnaReader(io.BytesIO(data)).blocks(labels=True):
            by_block.append(list(zip(seqs, zip(*cols))))
        picked = [0, 2]
        got = [p for i in picked for p in by_block[i]]
        flat = []
        for seqs, _f, cols in ZnaReader(io.BytesIO(data)).blocks(
                indices=picked, labels=True):
            flat.extend(zip(seqs, zip(*cols)))
        self.assertEqual(flat, got, "indices= misaligned labels")

        # restore_strand on a labeled file.
        want_rs = [r[0] for r in
                   ZnaReader(io.BytesIO(data)).records(restore_strand=True)]
        got_rs = []
        for seqs, _f in ZnaReader(io.BytesIO(data)).blocks(
                restore_strand=True, labels=False):
            got_rs.extend(seqs)
        self.assertEqual(got_rs, want_rs, "restore_strand diverged on a labeled file")


class TestBlockIndex(FuzzCase):
    """``block_index()`` is what lets a consumer size a subsample before reading
    anything, so its counts have to be exact — an estimate would silently skew a
    balanced corpus sample."""

    def _file(self, n=400, block_size=BLOCK_SIZE, seq_len_bytes=2):
        rng = random.Random(SEED + 16)
        recs = gen_records(rng, "mixed", n, seq_len_bytes,
                           allow_n=False, allow_lower=False)
        header = make_header(seq_len_bytes, STRAND_CONFIGS[0], COMPRESSION_ZSTD)
        buf = io.BytesIO()
        with ZnaWriter(buf, header, block_size=block_size) as w:
            for seq, is_paired, is_read1, is_read2, is_full in recs:
                w.write_record(seq, is_paired, is_read1, is_read2,
                               is_full_fragment=is_full)
        return buf.getvalue(), recs

    def test_counts_are_exact(self):
        data, recs = self._file()
        reader = ZnaReader(io.BytesIO(data))
        index = reader.block_index()
        self.assertEqual(sum(b.n_records for b in index), len(recs),
                         "block_index() record total disagrees with the file")
        # Agrees block-for-block with what blocks() actually yields.
        got = [len(seqs) for seqs, _f in ZnaReader(io.BytesIO(data)).blocks()]
        self.assertEqual([b.n_records for b in index], got)
        self.assertEqual([b.index for b in index], list(range(len(index))))

    def test_offsets_point_at_block_headers(self):
        data, _ = self._file()
        index = ZnaReader(io.BytesIO(data)).block_index()
        for b in index:
            fh = io.BytesIO(data)
            fh.seek(b.offset)
            raw = fh.read(_BLOCK_HEADER_SIZE)
            comp_size, uncomp_size, count, _fs, _ls = struct.unpack(
                _BLOCK_HEADER_FMT, raw)
            self.assertEqual(
                (comp_size, uncomp_size, count),
                (b.comp_size, b.uncomp_size, b.n_records),
                f"block {b.index}: offset does not point at its header",
            )

    def test_does_not_disturb_read_position(self):
        """The cursor is restored, so an in-flight iteration is unaffected.

        (Note this is about the *cursor*, not about ``records()`` rewinding — a
        second ``records()`` call on the same reader has always continued from
        wherever the first one stopped, independently of this method.)
        """
        data, recs = self._file()
        reader = ZnaReader(io.BytesIO(data))
        stream = reader.records()
        first = [next(stream) for _ in range(5)]

        pos_before = reader._fh.tell()
        index = reader.block_index()
        self.assertEqual(reader._fh.tell(), pos_before,
                         "block_index() left the cursor somewhere else")
        self.assertGreater(len(index), 0)

        rest = list(stream)          # the same generator carries on
        self.assertEqual(first + rest,
                         list(ZnaReader(io.BytesIO(data)).records()),
                         "block_index() perturbed an in-flight iteration")

    def test_works_before_any_read_and_after_full_read(self):
        data, recs = self._file()
        reader = ZnaReader(io.BytesIO(data))
        before = reader.block_index()
        self.assertEqual(sum(b.n_records for b in before), len(recs))
        list(reader.records())       # drain to EOF
        after = reader.block_index()
        self.assertEqual(before, after,
                         "block_index() differs once the stream is exhausted")

    def test_requires_a_seekable_stream(self):
        class Unseekable(io.RawIOBase):
            def __init__(self, data):
                self._b = io.BytesIO(data)
            def readable(self):
                return True
            def read(self, n=-1):
                return self._b.read(n)
            def readinto(self, b):
                return self._b.readinto(b)
            def tell(self):
                raise OSError("not seekable")
            def seekable(self):
                return False

        data, _ = self._file()
        reader = ZnaReader(Unseekable(data))
        with self.assertRaises(TypeError):
            reader.block_index()

    def test_indices_selects_exactly_those_blocks(self):
        data, _ = self._file(n=600, block_size=256)
        reader = ZnaReader(io.BytesIO(data))
        index = reader.block_index()
        self.assertGreater(len(index), 4, "need several blocks to select among")

        want = [b.index for b in index][::3]
        got_counts = [
            len(seqs) for seqs, _f in
            ZnaReader(io.BytesIO(data)).blocks(indices=want)
        ]
        self.assertEqual(got_counts, [index[i].n_records for i in want])

        # And the sequences match those blocks read individually.
        by_block = [seqs for seqs, _f in ZnaReader(io.BytesIO(data)).blocks()]
        selected = [s for i in want for s in by_block[i]]
        flat = [s for seqs, _f in
                ZnaReader(io.BytesIO(data)).blocks(indices=want) for s in seqs]
        self.assertEqual(flat, selected)

    def test_indices_edge_cases(self):
        data, _ = self._file(n=600, block_size=256)
        n_blocks = len(ZnaReader(io.BytesIO(data)).block_index())

        # Empty selection yields nothing (and warns, like an empty shard).
        with self.assertWarns(RuntimeWarning):
            self.assertEqual(
                list(ZnaReader(io.BytesIO(data)).blocks(indices=[])), [])

        # Out-of-range indices are simply never matched.
        with self.assertWarns(RuntimeWarning):
            self.assertEqual(
                list(ZnaReader(io.BytesIO(data)).blocks(indices=[n_blocks + 5])), [])

        # Duplicates collapse; order of the argument does not matter.
        a = [len(s) for s, _ in
             ZnaReader(io.BytesIO(data)).blocks(indices=[2, 0, 2])]
        b = [len(s) for s, _ in ZnaReader(io.BytesIO(data)).blocks(indices=[0, 2])]
        self.assertEqual(a, b)

        for bad in ([-1], [1.5], [True]):
            with self.assertRaises(ValueError, msg=f"indices={bad}"):
                ZnaReader(io.BytesIO(data)).blocks(indices=bad)
        with self.assertRaises(ValueError):
            ZnaReader(io.BytesIO(data)).blocks(stride=2, indices=[0])

    def test_subsampling_is_representative(self):
        """The property the corpus balancing depends on: taking k of n blocks
        yields about k/n of the records."""
        data, recs = self._file(n=4000, block_size=1024)
        index = ZnaReader(io.BytesIO(data)).block_index()
        total = sum(b.n_records for b in index)
        self.assertGreater(len(index), 10)

        rng = random.Random(SEED + 17)
        k = max(1, len(index) // 4)
        want = rng.sample([b.index for b in index], k)
        got = sum(len(s) for s, _f in
                  ZnaReader(io.BytesIO(data)).blocks(indices=want))
        expected = total * k / len(index)
        # Generous: the final block is a partial, so exactness is not the claim.
        self.assertLess(
            abs(got - expected) / expected, 0.35,
            f"sampling {k}/{len(index)} blocks gave {got} records, "
            f"expected about {expected:.0f}",
        )


class TestDecoderMemory(FuzzCase):
    """The C++ decoder builds its list, tuples and strings with the raw C API
    and GC-untracks the record tuples.  That is exactly where a rewrite leaks or
    over-frees without any test noticing, so check both directly."""

    def setUp(self):
        if "accel" not in BACKENDS:
            self.skipTest("C++ accel backend not available")

    def _block(self, n=400, seq_len_bytes=2):
        from zna import _accel
        rng = random.Random(SEED + 10)
        recs = gen_records(rng, "mixed", n, seq_len_bytes,
                           allow_n=False, allow_lower=False)
        seqs = [r[0] for r in recs]
        flags = [(1 if r[2] else 0) | (2 if r[3] else 0) | (4 if r[1] else 0)
                 for r in recs]
        fl, ln, sq = _accel.encode_block(seqs, list(flags), seq_len_bytes,
                                         "", True, False, False)
        return seqs, fl, ln, sq, len(seqs)

    def test_repeated_decode_does_not_leak(self):
        from zna import _accel
        seqs, fl, ln, sq, n = self._block()

        def once(**kw):
            recs = _accel.decode_block(fl, ln, sq, 2, n, **kw)
            self.assertEqual(len(recs), n)
            del recs

        for kw in ({}, {"with_rc": False}, {"with_rc": False, "restore_strand": True}):
            for _ in range(5):          # settle allocator and caches
                once(**kw)
            gc.collect()
            before = sys.getallocatedblocks()
            for _ in range(50):
                once(**kw)
            gc.collect()
            after = sys.getallocatedblocks()
            # A per-record leak over 50 x 400 records would be 20k blocks; a
            # small drift from caches and free lists is expected and fine.
            self.assertLess(
                after - before, 1000,
                f"decode_block({kw}) leaked {after - before} allocated blocks "
                f"over 50 iterations of {n} records",
            )

    def test_records_survive_collection_while_held(self):
        """Untracked tuples must still be kept alive by the list that holds
        them: an over-eager untrack shows up as freed memory under GC pressure."""
        from zna import _accel
        seqs, fl, ln, sq, n = self._block()
        recs = _accel.decode_block(fl, ln, sq, 2, n)
        for _ in range(3):
            gc.collect()
            # Churn the heap so a freed object would likely be reused.
            _ = [bytearray(4096) for _ in range(200)]
            gc.collect()
        for i, rec in enumerate(recs):
            self.assertEqual(len(rec), 5, f"rec{i}: width changed after GC")
            self.assertIsInstance(rec[0], str, f"rec{i}: sequence is not a str")
            self.assertEqual(len(rec[0]), len(seqs[i]), f"rec{i}: length changed")

    def test_labeled_decode_does_not_leak(self):
        from zna import _accel
        rng = random.Random(SEED + 11)
        codes = sorted(DTYPE_BY_CODE)
        n = 300
        recs = gen_records(rng, "mixed", n, 2, allow_n=False, allow_lower=False)
        seqs = [r[0] for r in recs]
        flags = [(1 if r[2] else 0) | (4 if r[1] else 0) for r in recs]
        col_data, col_sizes = [], []
        for c in codes:
            dt = DTYPE_BY_CODE[c]
            vals = [gen_label_value(rng, c) for _ in range(n)]
            col_data.append(struct.pack(f"<{n}{dt.struct_ch}", *vals))
            col_sizes.append(dt.size)
        fl, lb, ln, sq = _accel.encode_block_labeled(
            seqs, list(flags), 2, "", False, False, False, col_data, col_sizes)
        dtype_codes = "".join(codes)

        def once():
            out = _accel.decode_block_labeled(fl, ln, sq, 2, n, col_data,
                                              col_sizes, dtype_codes)
            self.assertEqual(len(out), n)
            self.assertEqual(len(out[0][5]), len(codes))
            del out

        for _ in range(5):
            once()
        gc.collect()
        before = sys.getallocatedblocks()
        for _ in range(30):
            once()
        gc.collect()
        after = sys.getallocatedblocks()
        self.assertLess(
            after - before, 1000,
            f"decode_block_labeled leaked {after - before} allocated blocks",
        )


def _main(argv):
    import argparse
    ap = argparse.ArgumentParser(description="Deep round-trip fuzz for the ZNA codec.")
    ap.add_argument("--iters", type=int, default=25,
                    help="matrix sweeps, each with fresh random data")
    ap.add_argument("--seed", type=int, default=SEED)
    args = ap.parse_args(argv)

    mod = sys.modules[__name__]
    for i in range(args.iters):
        # The test bodies read these as module globals, so rebind them there.
        mod.SEED, mod.ITERS = args.seed + i * 1009, 1
        suite = unittest.defaultTestLoader.loadTestsFromModule(mod)
        result = unittest.TextTestRunner(verbosity=1).run(suite)
        if not result.wasSuccessful():
            print(f"!! seed {mod.SEED} FAILED", file=sys.stderr)
            return 1
        print(f"-- seed {mod.SEED} ok ({i + 1}/{args.iters})")
    return 0


if __name__ == "__main__":
    sys.exit(_main(sys.argv[1:]))
