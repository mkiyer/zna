"""Argument definitions for ``zna merge`` — deliberately free of heavy imports.

Split out of :mod:`zna.merge.cli` so that *registering* the subcommand costs nothing:
``cli.py`` reaches the backend, the extension module and the 64 KiB consensus table,
none of which belong in the startup of ``zna inspect --json``, which advertises itself
as fast enough to catalogue a whole corpus. So ``zna/cli.py`` imports *this* module to
build the parser and only reaches for ``cli.run`` once the user has asked for
``merge``.

This module imports nothing but ``argparse``. Keep it that way.
"""
from __future__ import annotations

import argparse
import os

#: More threads than this cannot help -- see --threads.
_DEFAULT_THREADS = min(4, os.cpu_count() or 1)

_DESCRIPTION = ("Overlap-merge paired-end reads into one mixed interleaved FASTQ for ZNA "
                "encoding: each pair becomes one full-fragment read or stays two mates.")

_EPILOG = """\
Every pair is scored ONCE: R1 is slid against revcomp(R2) over the single axis of
candidate fragment lengths, and each shift gets a log-likelihood ratio in BITS --
~+2 bits per matching base (log2 4), ~-6.2 bits per mismatch at a 1% error rate. The
best-scoring shift (argmax, not first-accept) is MERGED when both hold:

  * its score reaches T = log2(N / alpha) bits, N = len1 + len2 - 1 candidate shifts
    (28.2 bits at 2x150, 26.6 at 2x50, 29.2 at 2x300): at most alpha chance merges of
    UNRELATED sequence per pair;
  * its mismatches are plausible as sequencing error: a true overlap of n bases shows
    more than dfit[n] of them with probability below the same alpha (dfit = 9 of 122
    at a 1% error rate). A long divergent repeat that outscores a short true overlap
    fails this, and the pair is KEPT -- nothing is searched for in its place.

Anything else keeps both reads unchanged. A merged record is R1 then R2's
non-overlapping tail, reverse-complemented; where the mates disagree the consensus
takes the better-supported call by posterior from the two Phred scores (and derates its
quality) -- no cutoffs, nothing to tune.

WHAT --alpha BOUNDS, AND WHAT IT DOES NOT. alpha caps chance merges of unrelated
sequence (measured: 0 in 40,000 random pairs at the default) and, through the
plausibility test, the rate at which true overlaps are refused. It does NOT bound merges
of genuinely homologous sequence -- a near-identical repeat (a perfect 15 bp match, or
7 mismatches in 64, plausible even at the default error rate) is plausible under any
error model and still merges; no threshold reaches zero there. Divergent repeats are
what the plausibility test catches, at the same alpha. Each factor of 10 in alpha moves T by 3.3 bits, i.e.
~1.7 matching bases.

--error-rate E (default 0.01) is the expected fraction of positions at which the two
mates DISAGREE where they truly overlap -- about twice the per-base sequencing error,
since either read can be wrong. It is not estimated; it is a setting, and it does two
different jobs:

  * in the SCORE it sets strictness: a larger E makes a mismatch cost less, so both
    sequencing errors and repeat divergence are tolerated more. It does NOT affect the
    alpha bound on chance merges, which holds at any E.
  * in the PLAUSIBILITY TEST it is a promise: a true overlap is refused with
    probability at most alpha only if E is at least the library's real disagreement
    rate. Set too low, true overlaps are refused and kept whole (never merged wrongly):
    a library that truly disagrees at 3% loses 0.07-0.6% of them at the default, at 5%
    1-13% (overlaps of 50-150 bases).

zna checks the setting against the data: --json reports
detected_overlap_mismatch_rate, the disagreement over every overlap the scan detected
(before the plausibility test), and expected_refused_true_overlap_fraction, the share of
those overlaps the test is expected to refuse at E were they all true and disagreeing at
that rate. The run warns when that share exceeds 0.1%, with the number of pairs the test
refused and a suggested value. It is a check, not an estimate: repeats inflate the rate
on clean libraries, and on degraded ones it can read below the truth. Production RNA
libraries measured ~0.9% on 0.5.3's comparable statistic, just under the default.

Raising E is a TRADE, not a free fix. It stops the plausibility test refusing true
overlaps, but it also makes every mismatch cost less in the score, so more short or
divergent overlaps merge -- false ones included. On a simulated 3'-degraded library
(true disagreement ~3.9%), E = 0.01 gave 240 wrong merges and E = 0.036 gave 745. Raise
it when the refused pairs matter more than that: typically poor or 3'-degraded
libraries (2-5%) whose run warns and refuses many pairs.

--adapter-trimmed DECLARES that no read extends past its molecule (adapters removed,
or reads clipped to the fragment). Read-through alignments then become impossible and
are never merged. zna checks the declaration on the first 100,000 pairs and warns when
more than 1% still look like strong read-through (honest declarations measure ~0.1%,
false ones 2-23%). WITHOUT the flag, read-through is inferred as before: that IS overlap-based
adapter removal -- a read contains adapter only when its insert is shorter than the
read, so the mates overlap fully and the merged record excludes the adapter. Declaring
it on raw reads is catastrophic (measured: 231,000 correct merges lost per million
raw hg38 pairs); leaving it off on trimmed reads only forgoes the gain.

Output is ONE mixed interleaved FASTQ: merged reads are singles, unmerged pairs are
adjacent /1,/2 records. Feed it to `zna encode --interleaved`, or merge in process with
`zna encode --merge-pairs`, which also records the policy in the file's prologue.

Example:
  zna merge --in1 R1.fq.gz --in2 R2.fq.gz --out merged.fq.gz --json merge.json

PERFORMANCE: the merge kernel is compiled C++ and releases the GIL, so --threads
are real threads. It is not usually the bottleneck -- gzip decompression is -- so
2-3 threads saturate and more does nothing. Output is byte-identical at any thread
count.
"""


class _Fmt(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    """Show defaults on each option AND preserve the epilog's formatting."""


def add_merge_algorithm_arguments(p):
    """The scoring knobs shared by ``zna merge`` and ``zna encode --merge-pairs``.

    One definition, two parsers: factoring these out is what keeps the two
    from drifting (MERGE_PAIRS_PLAN.md §5).  Only the knobs that parameterize
    the DECISION live here; I/O, threading and policy flags stay with their
    owning command (`zna encode` already has its own ``--npolicy``, ``--seed``
    and ``-q``, with identical semantics).
    """
    p.add_argument("--alpha", default="1e-6", metavar="ALPHA",
                   help="the policy's one statistical tolerance: at most ALPHA chance "
                        "merges of UNRELATED sequence per pair (the merge threshold is "
                        "log2((len1+len2-1)/ALPHA) bits, per pair), and at most ALPHA "
                        "true overlaps refused as implausible. It does not bound merges "
                        "of near-identical repeats, which no threshold can. Read as an "
                        "exact decimal; must be in [1e-300, 1).")
    p.add_argument("--error-rate", default="0.01", dest="error_rate", metavar="E",
                   help="expected fraction of positions at which the two mates disagree "
                        "in a TRUE overlap (~2x the per-base sequencing error). In the "
                        "score it sets how much a mismatch costs (strictness, not the "
                        "alpha guarantee); in the plausibility test it must be >= the "
                        "library's real rate, or true overlaps are refused and kept "
                        "whole. The run warns when, at the disagreement its detected "
                        "overlaps show, over 0.1%% of true overlaps would be refused; "
                        "raising it then recovers refused pairs but also merges more "
                        "short or divergent overlaps, false ones included (see below). "
                        "Read as an exact decimal; must be in [1e-9, 0.75).")
    p.add_argument("--adapter-trimmed", action="store_true", dest="adapter_trimmed",
                   help="declare that no read extends past its molecule (adapters "
                        "already removed, or reads clipped to the fragment): read-through "
                        "alignments become impossible. Checked on the first 100,000 "
                        "pairs, with a warning if the reads still look like they "
                        "contain adapter. Do NOT pass it for raw reads.")
    p.add_argument("--min-read-length", type=int, default=40, dest="min_read_length",
                   help="drop emitted reads shorter than this (a merged read is its "
                        "fragment; an unmerged pair is dropped whole if either mate is "
                        "short). MUST match the pipeline-wide floor used by any earlier "
                        "quality-trimming step.")
    p.add_argument("--no-sync-check", action="store_true",
                   help="skip the per-pair R1/R2 read-name consistency check. Only "
                        "for input whose mate names genuinely differ by design.")
    return p


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
    add_merge_algorithm_arguments(p)
    p.add_argument("--npolicy", choices=("trim3", "random"), default="trim3",
                   help="what to do with a no-call (N) that the overlap could not "
                        "rescue from the mate. trim3: cut the read at it, keeping "
                        "[0, first N) -- 3' only, so base 0 stays a true fragment "
                        "boundary however short the read gets, and the length filter "
                        "below discards what is left of a read that is mostly N. random: "
                        "substitute a base from a seeded stream (--seed), which "
                        "never costs a merge because it does not change a length. "
                        "Rescue from the mate happens first either way and costs "
                        "nothing. Same flag and same values as `zna encode`.")
    p.add_argument("--seed", type=int, default=42,
                   help="seed for --npolicy random. Substitution is derived from it and "
                        "the read's position, never from a running stream, so the output "
                        "does not depend on --chunk-size or --threads.")
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
