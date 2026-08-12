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

Example (production, multicore node):
  zna merge --in1 R1.fq.gz --in2 R2.fq.gz --out merged.fq.gz \\
            --json merge.json --threads 8 --processes 8

PERFORMANCE: the overlap scan is JIT-compiled with numba, which is an OPTIONAL zna
dependency -- install it with `pip install zna[merge]`. Without it the identical
scan runs as pure Python, correct but ~50x slower, so this command refuses to start
rather than turn a fast job into a silently slow one; pass --allow-slow if you mean
it. Set --processes to the number of allocated cores.
"""


class _Fmt(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    """Show defaults on each option AND preserve the epilog's formatting."""


def add_merge_arguments(p):
    """Add every ``zna merge`` flag to an existing parser. Returns it."""
    p.add_argument("--in1", required=True, help="R1 FASTQ (optionally .gz)")
    p.add_argument("--in2", required=True, help="R2 FASTQ (optionally .gz)")
    p.add_argument("--out", required=True, help="output mixed interleaved FASTQ (.gz to gzip)")
    p.add_argument("--json", help="write run statistics (length histogram, counts) as JSON")
    p.add_argument("--threads", type=int, default=4, help="pigz threads for gzip I/O")
    p.add_argument("--processes", type=int, default=1,
                   help="worker processes for the merge scan. 1 = single process. "
                        "Set to the allocated core count on a cluster node for a large "
                        "speedup (output order is irrelevant, so this is safe).")
    p.add_argument("--chunk-size", type=int, default=2000,
                   help="read pairs per work unit. Smaller chunks keep the workers fed "
                        "and bound the parent's memory; 50000 measured 1.24x slower "
                        "end to end because the pool drains between refills.")
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
    p.add_argument("--min-read-length", "--length-required", type=int, default=40,
                   dest="min_read_length",
                   help="drop emitted reads shorter than this (after merge/trim; a "
                        "trimmed read can fall below it). MUST match the pipeline-wide "
                        "floor used by the initial fastp run. `--length-required` is an alias.")
    p.add_argument("--no-sync-check", action="store_true",
                   help="skip the per-pair R1/R2 read-name consistency check")
    p.add_argument("--allow-empty", action="store_true",
                   help="exit 0 on an input with no read pairs. Off by default: an "
                        "empty input otherwise succeeds silently all the way to a "
                        "0-record .zna, and the library vanishes from the corpus with "
                        "every stage green.")
    p.add_argument("--allow-slow", action="store_true",
                   help="run without numba. Off by default: the pure-Python scan is "
                        "correct but ~50x slower, and a silently-correct 50x slowdown "
                        "on a cluster is indistinguishable from a slow node.")
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
