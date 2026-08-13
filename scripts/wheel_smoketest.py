"""Smoke-test an installed ZNA wheel.

Run by cibuildwheel against every wheel before it is published, on every
interpreter and platform.  Wheels were previously uploaded to PyPI without ever
being imported; this is the cheapest check that the compiled extension actually
loads and works on the interpreter it claims.

Kept as a file rather than an inline ``python -c`` so the quoting is the same on
cmd.exe and sh.

    python scripts/wheel_smoketest.py
"""
from __future__ import annotations

import io
import random
import sys


def main() -> int:
    import zna
    from zna import ZnaHeader, ZnaReader, ZnaWriter, FLAG_FIELDS, ENDS_BY_FLAG
    from zna.codec import available_backends, get_backend_name
    from zna.core import COMPRESSION_ZSTD
    from zna.dtypes import DTYPE_BY_CODE, LabelDef

    print(f"zna {zna.__version__} on {sys.implementation.name} "
          f"{'.'.join(map(str, sys.version_info[:3]))} ({sys.platform})")
    print(f"  backends: {available_backends()}  active: {get_backend_name()}")

    # The wheel exists to ship the compiled extension.  A wheel that silently
    # falls back to the pure-Python backend is a broken wheel.
    if get_backend_name() != "accel":
        print("FAIL: the C++ backend did not load from this wheel")
        return 1

    # There are TWO extensions, and both are imported as `_accel` within their own
    # package (`zna._accel`, `zna.merge._accel`).  Asserting only the codec catches the
    # merge target overwriting it, but not the reverse and not the merge target simply
    # failing to build on a platform -- in which case `zna merge` ships refusing to run.
    try:
        from zna.merge.backend import available_merge_backends
    except ImportError as exc:
        print(f"FAIL: zna.merge is not importable from this wheel: {exc}")
        return 1
    if "accel" not in available_merge_backends():
        print("FAIL: the compiled merge kernel (zna.merge._accel) is not in this wheel")
        return 1
    import zna.merge._accel as _merge_accel
    if not hasattr(_merge_accel, "scan") or hasattr(_merge_accel, "encode_block"):
        print("FAIL: zna.merge._accel is not the merge extension "
              "(the two _accel targets collided)")
        return 1

    rng = random.Random(20260812)
    n = 5000
    seqs = ["".join(rng.choices("ACGT", k=rng.choice([0, 1, 3, 150, 151])))
            for _ in range(n)]
    labels = (LabelDef(0, "ZI", "", DTYPE_BY_CODE["i"], missing=-1),
              LabelDef(1, "ZF", "", DTYPE_BY_CODE["C"], missing=0))
    values = [(i, i % 256) for i in range(n)]

    header = ZnaHeader(
        read_group="smoke", seq_len_bytes=2,
        strand_specific=True, read1_antisense=True, strand_normalized=True,
        compression_method=COMPRESSION_ZSTD, compression_level=3, labels=labels,
    )
    buf = io.BytesIO()
    with ZnaWriter(buf, header, block_size=64 * 1024) as w:
        for i, seq in enumerate(seqs):
            w.write_record(seq, True, i % 2 == 0, i % 2 == 1, labels=values[i])
    data = buf.getvalue()

    # records(): the default 5-tuple width for a labeled file
    recs = list(ZnaReader(io.BytesIO(data)).records())
    assert len(recs) == n, f"records(): {len(recs)} != {n}"
    assert all(len(r) == 5 for r in recs), "labeled records() width changed"

    # restore_strand must invert the encoder's reverse-complement
    restored = [r[0] for r in
                ZnaReader(io.BytesIO(data)).records(restore_strand=True)]
    assert restored == seqs, "restore_strand did not recover the input"

    # blocks(): columnar batch path, with and without labels
    got, cols = [], []
    for s, f, c in ZnaReader(io.BytesIO(data)).blocks(labels=True):
        assert len(f) == len(s)
        got.extend(s)
        cols.extend(zip(*c))
    assert got == [r[0] for r in recs], "blocks(labels=True) sequences differ"
    assert cols == values, "blocks(labels=True) label values differ"
    skipped = sum(len(s) for s, _f in
                  ZnaReader(io.BytesIO(data)).blocks(labels=False))
    assert skipped == n, "blocks(labels=False) record count differs"

    # block_index(): counts without decompressing
    index = ZnaReader(io.BytesIO(data)).block_index()
    assert sum(b.n_records for b in index) == n, "block_index() count differs"
    assert len(index) > 1, "expected several blocks at this block_size"

    # flag tables, and per-record fragment geometry
    flags0 = next(iter(ZnaReader(io.BytesIO(data)).blocks(labels=False)))[1]
    assert FLAG_FIELDS[flags0[0]][0] is True, "FLAG_FIELDS: expected is_paired"
    assert len(ENDS_BY_FLAG) == 256

    # An N with no policy must still be an error, not silent corruption.
    try:
        bad = io.BytesIO()
        with ZnaWriter(bad, ZnaHeader(read_group="x", seq_len_bytes=1)) as w:
            w.write_record("ACGN", False, True, False)
        print("FAIL: an invalid base was accepted without an N-policy")
        return 1
    except ValueError:
        pass

    print(f"  OK: {n} records through records(), blocks(), block_index(), "
          f"{len(index)} blocks, restore_strand exact")
    return 0


if __name__ == "__main__":
    sys.exit(main())
