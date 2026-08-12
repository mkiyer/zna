"""Argument definitions for ``zna merge`` — deliberately free of heavy imports.

Split out of :mod:`zna.merge.cli` so that *registering* the subcommand costs nothing.
``cli.py`` pulls in ``overlap.py``, which imports numba: measured on this repo,
``import zna.cli`` is 40 ms but ``import zna.cli, zna.merge.cli`` is 210 ms, and
``zna inspect --json`` advertises itself as fast enough to catalogue a whole corpus.
So ``zna/cli.py`` imports *this* module to build the parser and only reaches for
``cli.run`` once the user has actually asked for ``merge``.

This module imports nothing but ``argparse``. Keep it that way.
"""
from __future__ import annotations

import argparse
import os

#: More threads than this cannot help -- see --threads.
_DEFAULT_THREADS = min(4, os.cpu_count() or 1)

_DESCRIPTION = ("Overlap-merge paired-end reads (and trim residual overlap on unmerged "
                "pairs) into one mixed interleaved FASTQ for ZNA encoding.")

_EPILOG = """\
Every pair is scored ONCE: R1 is slid against revcomp(R2) over the single axis of
candidate fragment lengths, and each shift gets a log-likelihood ratio in BITS --
+2 bits per matching base (log2 4), -6.2 bits per mismatch at a 1% error rate. The
best-scoring shift (argmax, not first-accept) is then read at two thresholds:

  * score >= --threshold-merge                -> MERGE into one read
      (R1 wins the overlap; R2's non-overlapping tail is appended, reverse-
       complemented). Emitted as a single record with the /1,/2 suffix stripped.
  * --threshold-trim <= score < merge         -> KEEP BOTH, but trim the redundant
      bases off R2's 3' end so the overlap is not counted twice.
  * score < --threshold-trim                  -> KEEP BOTH, unchanged.

Both thresholds are on one calibrated scale: T bits tolerates a spurious rate of
about N * 2**-T over the N ~ 2*readlen candidate shifts, so the default 28 is one
spurious merge in 1e6 pairs AGAINST CHANCE ALIGNMENT (measured: 0 in 40,000
uniform-random pairs, at every read length from 50 to 300). It is not a bound
against real sequence, where reads share genuine homology and repeat content --
raising T there buys far less than the formula suggests. Trim keeps a much lower
threshold only because a wrong trim costs a few real bases while a wrong merge is
a chimera.

Output is ONE mixed interleaved FASTQ: merged reads are singles, unmerged pairs are
adjacent /1,/2 records. Feed it to `zna encode --interleaved`. Where the mates
overlap and DISAGREE, the consensus takes the better-supported call by posterior
from the two Phred scores (and derates its quality, because a contested base is
less certain) -- no cutoffs, nothing to tune. Defaults suit 2x150 bp data; you
normally set nothing.

Example:
  zna merge --in1 R1.fq.gz --in2 R2.fq.gz --out merged.fq.gz --json merge.json

PERFORMANCE: the merge kernel is compiled C++ and releases the GIL, so --threads
are real threads. It is not usually the bottleneck -- gzip decompression is -- so
2-3 threads saturate and more does nothing. Output is byte-identical at any thread
count.
"""


class _Fmt(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    """Show defaults on each option AND preserve the epilog's formatting."""


def add_merge_arguments(p):
    """Add every ``zna merge`` flag to an existing parser. Returns it."""
    p.add_argument("--in1", required=True, help="R1 FASTQ (optionally .gz)")
    p.add_argument("--in2", required=True, help="R2 FASTQ (optionally .gz)")
    p.add_argument("--out", required=True, help="output mixed interleaved FASTQ (.gz to gzip)")
    p.add_argument("--json", help="write run statistics (length histogram, counts) as JSON")
    p.add_argument("--threads", type=int, default=_DEFAULT_THREADS,
                   help="merge worker threads. The merge kernel releases the GIL, so "
                        "these are real. Compute is ~1 us/pair against a ~0.8 us/pair "
                        "gzip decompression floor, so 2-3 saturates and more does "
                        "nothing.")
    p.add_argument("--io-threads", type=int, default=4,
                   help="pigz threads for the gzipped OUTPUT. The reader always uses "
                        "one, because pigz cannot parallelise inflate and extra reader "
                        "threads only contend with the workers.")
    p.add_argument("--chunk-size", type=int, default=2000,
                   help="read pairs per work unit. Bounds memory and sets the "
                        "parallel granularity.")
    p.add_argument("--compress-level", type=int, default=1,
                   help="pigz level for --out. Default 1 (fast): the output is an "
                        "intermediate consumed by `zna encode`, so speed beats ratio. "
                        "Raise for archival/standalone use.")
    p.add_argument("--threshold-merge", type=float, default=28.0, dest="t_merge",
                   help="overlap score (bits) >= this -> merge the pair into one "
                        "full-fragment read. 28 bits = one spurious merge per 1e6 "
                        "pairs at 2x150; each extra bit halves that.")
    p.add_argument("--threshold-trim", type=float, default=8.0, dest="t_trim",
                   help="overlap score (bits) >= this (but below --threshold-merge) "
                        "-> keep both reads and trim the redundant overlap off R2's "
                        "3' end. Low on purpose: a wrong trim costs a few bases.")
    p.add_argument("--min-read-length", type=int, default=40, dest="min_read_length",
                   help="drop emitted reads shorter than this (after merge/trim; a "
                        "trimmed read can fall below it). MUST match the pipeline-wide "
                        "floor used by any earlier quality-trimming step.")
    p.add_argument("--no-sync-check", action="store_true",
                   help="skip the per-pair R1/R2 read-name consistency check")
    p.add_argument("--allow-empty", action="store_true",
                   help="exit 0 on an input with no read pairs. Off by default: an "
                        "empty input otherwise succeeds silently all the way to a "
                        "0-record .zna, and the library vanishes from the corpus with "
                        "every stage green.")
    p.add_argument("--backend", choices=("auto", "accel", "python"), default="auto",
                   help="merge kernel. `auto` uses the compiled backend and fails if it "
                        "is missing; `python` selects the reference implementation, "
                        "which is correct but ~50x slower and exists to be an oracle, "
                        "not a fallback. A silently-correct 50x slowdown on a cluster "
                        "is indistinguishable from a slow node, so it is never chosen "
                        "for you.")
    p.add_argument("-q", "--quiet", action="store_true", help="suppress progress logging")
    return p


def add_merge_parser(subparsers):
    """Register ``merge`` on zna's top-level subparsers. Returns the subparser."""
    p = subparsers.add_parser(
        "merge",
        help="Overlap-merge paired-end FASTQ into one interleaved FASTQ",
        formatter_class=_Fmt,
        description=_DESCRIPTION,
        epilog=_EPILOG,
    )
    return add_merge_arguments(p)


def build_parser() -> argparse.ArgumentParser:
    """Stand-alone parser, for ``python -m zna.merge`` and for the test suite.

    Kept as a thin wrapper over the same flag definitions the subcommand uses, so the
    two can never drift.
    """
    p = argparse.ArgumentParser(
        prog="zna merge",
        formatter_class=_Fmt,
        description=_DESCRIPTION,
        epilog=_EPILOG,
    )
    return add_merge_arguments(p)
