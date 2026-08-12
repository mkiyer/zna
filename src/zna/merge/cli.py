"""CLI for zna read-merge (the ``zna merge`` subcommand).

Reads two positionally-synced FASTQ files, merges/trims/keeps each pair, and writes one
mixed interleaved FASTQ stream for ``zna encode --interleaved`` (see
``docs/READ_MERGE_REDESIGN.md``). Also invocable as ``python -m zna.merge``.

Single-process by default. With ``--processes N`` (N>1) the CPU-bound merge work is
fanned out to a worker pool; output order is irrelevant (ZNA shuffles), so workers
process chunks independently and the main process interleaves their output blobs.
"""
from __future__ import annotations

import json
import logging
import platform
import shutil
import sys
import time

from .args import add_merge_arguments, add_merge_parser, build_parser  # noqa: F401
from .fastqio import FastqWriter, InputError, read_pairs
from .overlap import HAVE_NUMBA, score_weights, threshold_bits
from .pairs import MergeParams, PairOutcome, base_name, process_pair

logger = logging.getLogger("zna.merge")

_HIST_MAX = 1024   # length/overlap histograms are clamped to this many bins


# --------------------------------------------------------------------------- #
# accounting
#
# There is exactly ONE implementation of the per-pair bookkeeping (`_process_chunk`),
# used by both the single-process and the worker paths. It used to be duplicated —
# ~35 lines, once per path — which meant every statistic added had to be added twice
# and the parallel-vs-single test could only compare a hardcoded list of keys. The two
# copies were verified equivalent, but the duplication was the standing hazard, not the
# divergence.
# --------------------------------------------------------------------------- #

def _new_acc():
    # [n_pairs, merged, trimmed, kept, n_emitted, n_dropped, bases_trimmed,
    #  fragments_dropped_short_mate, length_hist, bases_consensus, overlap_hist,
    #  trim_guard_kept, sum_olen, sum_diff, insert_hist]
    return [0, 0, 0, 0, 0, 0, 0, 0, [0] * (_HIST_MAX + 1), 0,
            [0] * (_HIST_MAX + 1), 0, 0, 0, [0] * (_HIST_MAX + 1)]


def _fold(result, acc):
    """Fold one chunk's counters into the accumulator (in-place). Returns the blob."""
    (blob, merged, trimmed, kept, n_emitted, n_dropped, bases_trimmed,
     frags_short, hist, bases_consensus, ohist, trim_guard,
     sum_olen, sum_diff, ihist) = result
    acc[1] += merged; acc[2] += trimmed; acc[3] += kept
    acc[4] += n_emitted; acc[5] += n_dropped; acc[6] += bases_trimmed
    acc[7] += frags_short; acc[9] += bases_consensus; acc[11] += trim_guard
    acc[12] += sum_olen; acc[13] += sum_diff
    for src, dst in ((hist, acc[8]), (ohist, acc[10]), (ihist, acc[14])):
        for i, c in src.items():
            dst[i] += c
    return blob


# --------------------------------------------------------------------------- #
# the merge work itself. Runs in the main process (--processes 1) or in a forked
# worker; workers are set up via fork so they inherit the parent's already-JIT-
# compiled overlap kernel (no per-worker recompile).
# --------------------------------------------------------------------------- #

_W = {}  # worker-local state (set by _init_worker)


def _init_worker(params, check):
    _W["params"] = params
    _W["check"] = check


def _process_chunk(work):
    """Merge one chunk of pairs; return ``(output_blob, counters...)``.

    ``work`` is ``(base_index, chunk)``; ``base_index`` is the index of the chunk's
    first pair in the input, carried only so a desync can be reported by absolute pair
    number rather than "somewhere in this chunk".
    """
    base_index, chunk = work
    params = _W["params"]
    check = _W["check"]
    merged = trimmed = kept = n_emitted = n_dropped = bases_trimmed = frags_short = 0
    counters = [0, 0]                   # [bases_consensus_changed, trim_guard_kept]
    sum_olen = sum_diff = 0
    hist = [0] * (_HIST_MAX + 1)        # every emitted record's length
    ohist = [0] * (_HIST_MAX + 1)       # detected overlap lengths
    ihist = [0] * (_HIST_MAX + 1)       # inferred fragment length, merged pairs only
    parts = []
    ap = parts.append
    for i, (h1, s1, q1, h2, s2, q2) in enumerate(chunk):
        if check and base_name(h1) != base_name(h2):
            raise InputError(f"R1/R2 out of sync at pair {base_index + i + 1}: "
                             f"{base_name(h1)!r} != {base_name(h2)!r}")
        records, outcome, dropped, score, olen, diff = process_pair(
            h1, s1, q1, h2, s2, q2, params, counters)
        n_dropped += dropped
        if outcome == PairOutcome.MERGED:
            merged += 1
        elif outcome == PairOutcome.TRIMMED:
            trimmed += 1
            if records:
                # records is exactly [R1, R2] on this branch (pairs.py), so index it.
                # The old `next(r for r in records if r[0] != h1)` silently returned
                # None when the two mates carried identical headers — which happens
                # for any input not passed through `samtools fastq -N` — and then
                # charged the whole of R2 as trimmed, overstating by ~12x.
                bases_trimmed += len(s2) - len(records[1][1])
        else:
            kept += 1
        if not records and outcome != PairOutcome.MERGED:
            frags_short += 1
        if olen:
            ohist[olen if olen <= _HIST_MAX else _HIST_MAX] += 1
            sum_olen += olen
            sum_diff += diff
        for header, seq, qual in records:
            ap(b"@%b\n%b\n+\n%b\n" % (header, seq, qual))
            n_emitted += 1
            L = len(seq)
            hist[L if L <= _HIST_MAX else _HIST_MAX] += 1
            if outcome == PairOutcome.MERGED:
                # A merged record IS the fragment, so its length is the insert size.
                ihist[L if L <= _HIST_MAX else _HIST_MAX] += 1
    # Sparse: a chunk touches ~150 of 1025 bins, and this tuple is pickled per chunk.
    sparse = lambda h: {i: c for i, c in enumerate(h) if c}
    return (b"".join(parts), merged, trimmed, kept, n_emitted, n_dropped,
            bases_trimmed, frags_short, sparse(hist), counters[0], sparse(ohist),
            counters[1], sum_olen, sum_diff, sparse(ihist))


def _iter_chunks(pair_iter, size):
    """Yield ``(base_index, chunk)``; splits only on whole pairs."""
    chunk = []
    ap = chunk.append
    base = 0
    for r1, r2 in pair_iter:
        ap((r1[0], r1[1], r1[2], r2[0], r2[1], r2[2]))
        if len(chunk) >= size:
            yield base, chunk
            base += len(chunk)
            chunk = []
            ap = chunk.append
    if chunk:
        yield base, chunk


# --------------------------------------------------------------------------- #
# single-process path (default; exactly what the unit tests exercise)
# --------------------------------------------------------------------------- #

def _run_single(args, params):
    _init_worker(params, not args.no_sync_check)
    acc = _new_acc()
    with FastqWriter(args.out, threads=args.threads, level=args.compress_level) as w:
        try:
            for work in _iter_chunks(read_pairs(args.in1, args.in2, args.threads),
                                     args.chunk_size):
                w.write_raw(_fold(_process_chunk(work), acc))
                acc[0] = work[0] + len(work[1])
                if not args.quiet and acc[0] % 5_000_000 < args.chunk_size:
                    logger.info("processed %d pairs", acc[0])
        except InputError as e:                       # desync, or a malformed record
            raise SystemExit(str(e))
    acc[0] = acc[1] + acc[2] + acc[3]
    return acc


# --------------------------------------------------------------------------- #
# parallel path (--processes > 1)
# --------------------------------------------------------------------------- #

def _run_parallel(args, params):
    """Fan the chunks out to a worker pool.

    Uses ``ProcessPoolExecutor`` rather than ``multiprocessing.Pool`` because the
    latter cannot detect abrupt worker death: a SIGKILLed worker (cgroup/node OOM, a
    native crash in JIT code) leaves ``IMapUnorderedIterator.next()`` blocked forever
    while ``_maintain_pool`` quietly respawns the process. Measured: the parent sat
    alive at 0% CPU for 150 s with no output, no JSON and no error — at cluster scale
    that is indistinguishable from a slow node, and it burns the whole allocation.
    ``BrokenProcessPool`` surfaces the same failure in about a second.

    Submission is windowed at ``2 * processes`` chunks so the input is streamed rather
    than read into memory: peak parent RSS is bounded by ``chunk_size``, not by the
    library.
    """
    import concurrent.futures as cf
    import multiprocessing as mp
    # Warm the JIT in the parent so forked workers inherit the compiled kernel.
    process_pair(b"w", b"ACGT" * 10, b"I" * 40, b"w", b"ACGT" * 10, b"I" * 40, params)
    acc = _new_acc()
    # fork is deliberate: the workers inherit the parent's already-JIT-compiled kernel
    # (warmed just above) instead of each paying the numba compile. Windows has no fork
    # at all, so fall back rather than crash on --processes > 1 -- the kernels pass
    # cache=True, so after the first run each worker loads the compiled artifact from
    # the numba cache instead of rebuilding it.
    try:
        ctx = mp.get_context("fork")
    except ValueError:
        ctx = mp.get_context()
        logger.info("no fork on this platform: using %r, so each worker compiles the "
                    "overlap kernel once (set NUMBA_CACHE_DIR to a writable path if "
                    "that repeats every run)", ctx.get_start_method())
    # One decompression thread per reader: with N workers already resident, extra
    # pigz threads only contend for the same allocation (measured 1.13x).
    chunks = _iter_chunks(read_pairs(args.in1, args.in2, 1), args.chunk_size)
    window = max(2, 2 * args.processes)
    with FastqWriter(args.out, threads=args.threads, level=args.compress_level) as w, \
            cf.ProcessPoolExecutor(max_workers=args.processes, mp_context=ctx,
                                   initializer=_init_worker,
                                   initargs=(params, not args.no_sync_check)) as pool:
        pending = set()
        try:
            for work in chunks:
                pending.add(pool.submit(_process_chunk, work))
                if len(pending) >= window:
                    done, pending = cf.wait(pending, return_when=cf.FIRST_COMPLETED)
                    for fut in done:
                        w.write_raw(_fold(fut.result(), acc))
            while pending:
                done, pending = cf.wait(pending, return_when=cf.FIRST_COMPLETED)
                for fut in done:
                    w.write_raw(_fold(fut.result(), acc))
        except InputError as e:                       # desync raised in a worker
            raise SystemExit(str(e))
        except cf.process.BrokenProcessPool:
            raise SystemExit(
                "a merge worker died abruptly (killed by the OS, or a native crash). "
                "The output is incomplete; re-run, and if this repeats check the "
                "memory limit against --processes/--chunk-size."
            )
    acc[0] = acc[1] + acc[2] + acc[3]
    return acc


# --------------------------------------------------------------------------- #

def _assemble_stats(acc, params, elapsed=None):
    (n_pairs, merged, trimmed, kept, n_emitted, n_dropped, bases_trimmed,
     frags_short, hist, bases_consensus, ohist, trim_guard,
     sum_olen, sum_diff, ihist) = acc
    total_bases = sum(i * c for i, c in enumerate(hist))
    pct = (lambda n: round(100.0 * n / n_pairs, 3) if n_pairs else 0.0)
    match_w, mismatch_w = score_weights(params.err_rate)
    import zna
    stats = {
        # Provenance: the question every future corpus defect opens with is "which
        # build produced this file?". Config values are already cohort-queryable via
        # gather's pipeline tool; only code identity was missing.
        "tool": "zna-merge",
        "tool_version": zna.__version__,
        "numba": HAVE_NUMBA,          # a numba-less run is ~50x slower and SILENTLY correct
        "python": platform.python_version(),
        "input_pairs": n_pairs,
        "merged": merged,
        "trimmed_pairs": trimmed,
        "kept_pairs": kept,
        "merged_pct": pct(merged),
        "trimmed_pct": pct(trimmed),
        "kept_pct": pct(kept),
        "emitted_records": n_emitted,
        "dropped_below_min_length": n_dropped,
        "fragments_dropped_short_mate": frags_short,
        "bases_trimmed": bases_trimmed,
        "bases_consensus_changed": bases_consensus,
        # A trim that would have left R2 below min_read_length; both reads kept whole
        # instead of discarding the fragment. Expected to be ~0 (redesign §8).
        "trim_guard_kept_untrimmed": trim_guard,
        "mean_emitted_length": round(total_bases / n_emitted, 1) if n_emitted else 0.0,
        # Mismatches per aligned overlap base, over every detected overlap. This is a
        # direct calibration check on err_rate: it should sit near it (~0.009 measured
        # against the assumed 0.01). It is also the SENSITIVE degradation channel —
        # 5% per-base degradation moves it 5x while merged_pct moves 1% — whereas
        # merged_pct is the catastrophic one. Note it cannot approach the 0.75 of
        # chance alignment: the threshold caps the observable rate at ~0.22.
        "overlap_mismatch_rate": round(sum_diff / sum_olen, 6) if sum_olen else 0.0,
        "overlap_bases_compared": sum_olen,
        "params": {
            "threshold_merge_bits": params.t_merge,
            "threshold_trim_bits": params.t_trim,
            "err_rate": params.err_rate,
            "match_bits": round(match_w, 4),
            "mismatch_bits": round(-mismatch_w, 4),
            "min_read_length": params.min_read_length,
        },
        "length_histogram": {str(i): c for i, c in enumerate(hist) if c},
        # Detected overlap length per pair. This replaced a score histogram whose
        # integer bins aliased against the 1.9855-bits-per-base quantum (a 12.6:1
        # odd/even comb). Overlap length is the natural quantum, and the position of
        # its short-end cliff reads --threshold-merge directly: the shortest clean
        # overlap that can merge is ceil(t_merge / match_bits).
        "overlap_length_histogram": {str(i): c for i, c in enumerate(ohist) if c},
        # Inferred fragment length, merged pairs only (a merged record IS the
        # fragment). CENSORED AT BOTH ENDS, so do not read it as the library's insert
        # distribution without accounting for that: hard-floored at min_read_length
        # (shorter fragments are dropped, not observed) and hard-capped per pair at
        # len1 + len2 - ceil(t_merge / match_bits), beyond which the mates no longer
        # overlap enough to merge.
        "insert_size_histogram": {str(i): c for i, c in enumerate(ihist) if c},
        "insert_size_censoring": {
            "floor": params.min_read_length,
            "min_mergeable_overlap": int(-(-params.t_merge // match_w)),
        },
    }
    if elapsed is not None:
        stats["elapsed_s"] = round(elapsed, 1)
        # A cohort field that detects, for free, nodes where numba failed to load or
        # pigz was missing — both of which are silently correct and 10-50x slower.
        stats["pairs_per_second"] = round(n_pairs / elapsed) if elapsed > 0 else 0
    return stats


def _validate(args, params) -> None:
    """Reject argument combinations that would silently produce a wrong corpus.

    The realistic failure here is not a hand-typed flag, it is a config typo
    propagating into a whole-cohort run with nothing but ``merged_pct`` to signal it —
    ``--threshold-merge 0.5`` is accepted today and merges 94% of pairs instead of 82%.
    """
    if params.t_trim > params.t_merge:
        raise SystemExit("--threshold-trim must be <= --threshold-merge")
    if params.t_trim <= 0:
        raise SystemExit("--threshold-trim must be > 0 (it is evidence in bits; at 0 "
                         "every pair 'overlaps')")
    # Floor the merge threshold at the evidence needed to beat chance over the shift
    # space at all. threshold_bits() is the redesign's own derivation, so the bound is
    # executable rather than folklore; 50 bp is the shortest read worth merging.
    floor = threshold_bits(50, 0.05)
    if params.t_merge < floor:
        raise SystemExit(
            f"--threshold-merge {params.t_merge:g} is below {floor:.1f} bits, the point "
            f"at which chance alignments start passing even for 2x50 reads. The default "
            f"is 28. If you really mean it, you want a different tool."
        )
    if params.min_read_length < 1:
        raise SystemExit("--min-read-length must be >= 1")
    if args.processes < 1:
        raise SystemExit("--processes must be >= 1")
    if args.chunk_size < 1:
        raise SystemExit("--chunk-size must be >= 1")
    if args.threads < 1:
        raise SystemExit("--threads must be >= 1")


def run(args) -> dict:
    """Execute the merge and return the statistics dict.

    This is the in-process entry point: it does *not* refuse a numba-less run (see
    :func:`_require_numba`, which the command-line entry points call first), because
    a caller that reached this function chose the pure-Python scorer knowingly.
    """
    if not HAVE_NUMBA:
        logger.warning(
            "numba NOT found -> SLOW pure-Python overlap scan (~50x). "
            "Install it with `pip install zna[merge]` for production speed."
        )
    logger.info("numba: %s | gzip: %s | processes: %d",
                "active" if HAVE_NUMBA else "MISSING (slow)",
                "pigz" if shutil.which("pigz") else "stdlib gzip",
                max(1, args.processes))
    params = MergeParams(
        t_merge=args.t_merge,
        t_trim=args.t_trim,
        min_read_length=args.min_read_length,
    )
    _validate(args, params)

    t0 = time.perf_counter()
    acc = _run_parallel(args, params) if args.processes > 1 else _run_single(args, params)
    elapsed = time.perf_counter() - t0
    # An empty input is otherwise a silent success all the way down: rc=0 here, then a
    # 22-byte 0-record .zna, and a library disappears from the corpus with every stage
    # green. Cheaper to fail here than to find the hole in a trained model.
    if acc[0] == 0 and not args.allow_empty:
        raise SystemExit(
            f"no read pairs in {args.in1} / {args.in2}. If that is expected, pass "
            f"--allow-empty; otherwise the input is truncated or the wrong file."
        )
    stats = _assemble_stats(acc, params, elapsed)

    if args.json:
        with open(args.json, "w") as fh:
            json.dump(stats, fh, indent=2)
    if not args.quiet:
        logger.info(
            "done: %d pairs -> %d merged, %d trimmed, %d kept (%d records, %d dropped)",
            stats["input_pairs"], stats["merged"], stats["trimmed_pairs"],
            stats["kept_pairs"], stats["emitted_records"], stats["dropped_below_min_length"],
        )
    return stats


def _require_numba(args) -> None:
    """Refuse to start the pure-Python scan unless the user asked for it.

    The fallback in ``overlap.py`` is correct, and it is what keeps ``import
    zna.merge`` working on a plain ``pip install zna`` — but it is ~50x slower, and a
    *silently correct* 50x slowdown is the worst shape a failure can take at cluster
    scale: the job does not fail, it just looks like a slow node, and it burns the
    whole allocation before anyone looks. So the command line refuses; the library
    entry point (:func:`run`) does not.
    """
    if HAVE_NUMBA or args.allow_slow:
        return
    raise SystemExit(
        "numba is not installed, so the overlap scan would run as pure Python: "
        "correct, but about 50x slower. At cluster scale that is indistinguishable "
        "from a slow node. Install it with `pip install zna[merge]` (or add "
        "`numba >=0.61` to the conda environment). To run the slow path "
        "deliberately, pass --allow-slow."
    )


def run_command(args) -> int:
    """Entry point for ``zna merge`` and for ``python -m zna.merge``."""
    logging.basicConfig(
        level=logging.WARNING if args.quiet else logging.INFO,
        format="[zna merge] %(message)s",
        stream=sys.stderr,
    )
    _require_numba(args)
    run(args)
    return 0


def main(argv=None) -> int:
    return run_command(build_parser().parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
