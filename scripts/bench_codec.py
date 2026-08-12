"""End-to-end benchmarks for the C++ codec rewrite (Tier A/B).

Unlike the Tier C benchmarks, these cannot hold a local copy of the baseline —
the baseline is a different compiled extension.  So this script measures the
current build in absolute terms and, when pointed at a reference build with
``--baseline``, compares the two.

Building a reference:

    git stash && pip install -e . --no-build-isolation -q
    python scripts/bench_codec.py --save /tmp/base.json
    git stash pop && pip install -e . --no-build-isolation -q
    python scripts/bench_codec.py --baseline /tmp/base.json

Both halves run min-of-N in a fresh process on the same data, which is weaker
than the in-process interleaving used elsewhere; the numbers are reported with
that caveat attached.
"""
from __future__ import annotations

import argparse
import gc
import io
import json
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from zna.core import (  # noqa: E402
    COMPRESSION_ZSTD, ZnaHeader, ZnaReader, ZnaWriter,
)

READ_LEN = 150
N_READS = 200_000
ROUNDS = 9


def build_corpus(seed=20260811):
    rng = random.Random(seed)
    return ["".join(rng.choices("ACGT", k=READ_LEN)) for _ in range(N_READS)]


def make_file(seqs, *, strand_normalized=False, paired=True):
    hdr = ZnaHeader(
        read_group="bench", seq_len_bytes=2,
        strand_normalized=strand_normalized,
        compression_method=COMPRESSION_ZSTD, compression_level=3,
    )
    buf = io.BytesIO()
    with ZnaWriter(buf, hdr) as w:
        if paired:
            for i in range(0, len(seqs) - 1, 2):
                w.write_record(seqs[i], True, True, False)
                w.write_record(seqs[i + 1], True, False, True)
        else:
            for s in seqs:
                w.write_record(s, False, False, False)
    return buf.getvalue()


def timed(fn, rounds=ROUNDS):
    """Min-of-N wall time for *fn*."""
    best = float("inf")
    for _ in range(rounds):
        gc.collect()
        gc.disable()
        try:
            t0 = time.perf_counter()
            fn()
            best = min(best, time.perf_counter() - t0)
        finally:
            gc.enable()
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--save", help="write results as JSON")
    ap.add_argument("--baseline", help="compare against a saved JSON")
    args = ap.parse_args()

    seqs = build_corpus()
    plain = make_file(seqs)
    normalized = make_file(seqs, strand_normalized=True)

    results = {}

    def bench(name, fn):
        results[name] = timed(fn)

    bench("decode records() 4-tuple",
          lambda: sum(1 for _ in ZnaReader(io.BytesIO(plain)).records()))
    bench("decode records(with_rc=True)",
          lambda: sum(1 for _ in ZnaReader(io.BytesIO(plain)).records(with_rc=True)))
    bench("decode records(with_ends=True)",
          lambda: sum(1 for _ in ZnaReader(io.BytesIO(plain)).records(with_ends=True)))
    bench("decode records(restore_strand=True)",
          lambda: sum(1 for _ in
                      ZnaReader(io.BytesIO(normalized)).records(restore_strand=True)))
    bench("encode 200k paired reads", lambda: make_file(seqs))
    bench("encode 200k, strand-normalized",
          lambda: make_file(seqs, strand_normalized=True))

    from zna import _accel
    big = "".join(random.Random(3).choices("ACGT", k=1_000_000))
    bench("reverse_complement 1 Mbp x20",
          lambda: [_accel.reverse_complement(big) for _ in range(20)])

    # Raw backend calls.  The end-to-end figures above carry zstd, block framing
    # and generator overhead in the denominator; the optimization plan quotes
    # numbers "on the path changed", which is this.  Called positionally so the
    # pre-rewrite build, which has no with_rc/restore_strand parameters, runs the
    # identical benchmark.
    flags = [(1 if i % 2 == 0 else 2) | 4 for i in range(len(seqs))]
    fl, ln, sq = _accel.encode_block(seqs, list(flags), 2, "", False, False, False)
    n = len(seqs)
    bench("raw decode_block (200k records, one block)",
          lambda: _accel.decode_block(fl, ln, sq, 2, n))
    bench("raw encode_block (200k records)",
          lambda: _accel.encode_block(seqs, list(flags), 2, "", False, False, False))
    bench("raw encode_block, reverse-complementing",
          lambda: _accel.encode_block(seqs, list(flags), 2, "", True, False, False))

    base = json.load(open(args.baseline)) if args.baseline else None
    width = max(len(k) for k in results)
    print(f"\n  {'benchmark':<{width}}  {'time':>10s}" +
          ("   baseline      ratio  change" if base else ""))
    print("  " + "-" * (width + (44 if base else 12)))
    for k, v in results.items():
        line = f"  {k:<{width}}  {v * 1e3:8.2f} ms"
        if base and k in base:
            b = base[k]
            ratio = b / v
            pct = (ratio - 1) * 100
            tag = f"{pct:+.1f}%" if abs(pct) >= 5 else f"{pct:+.1f}% (NOISE)"
            line += f"  {b * 1e3:8.2f} ms  {ratio:5.2f}x  {tag}"
        print(line)
    print()

    if args.save:
        json.dump(results, open(args.save, "w"), indent=2)
        print(f"  saved -> {args.save}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
