"""CLI for zna read-merge (the ``zna merge`` subcommand).

Reads two positionally-synced FASTQ files, merges/trims/keeps each pair, and writes one
mixed interleaved FASTQ stream for ``zna encode --interleaved`` (see
``docs/METHODS.md``). Also invocable as ``python -m zna.merge``.

All per-pair work happens inside one backend call per chunk, which releases the GIL,
so ``--threads`` are real worker threads. Output is written in submission order and is
therefore byte-identical at any thread count.
"""
from __future__ import annotations

import json
import logging
import platform
import sys
import time

from .. import _gzip
from . import backend as _backend
from .args import add_merge_arguments, add_merge_parser, build_parser  # noqa: F401
from .fastqio import FastqWriter, InputError, _open_binary_read
from .params import DISAGREE_Q, SCALE, score_weights, threshold_bits
from .pairs import MergeParams

logger = logging.getLogger("zna.merge")

# --------------------------------------------------------------------------- #
# the merge loop
#
# All per-pair work -- parsing, scanning, consensus, construction, formatting and the
# histograms -- happens inside one backend call per chunk, which releases the GIL for
# its whole duration. That is what makes worker THREADS worth having where the previous
# design needed worker processes, and it deletes the fork context, the pickling of
# chunks and blobs, the per-worker globals, and the sparse-histogram workaround that
# only existed because dense ones were being pickled 
# --------------------------------------------------------------------------- #

#: Counter fields, in the order every backend returns them.
_N_COUNTERS = 15
(_N_PAIRS, _MERGED, _TRIMMED, _KEPT, _EMITTED, _DROPPED, _BASES_TRIMMED,
 _FRAGS_SHORT, _BASES_CONSENSUS, _TRIM_GUARD, _SUM_OLEN, _SUM_DIFF,
 _MAX_READ_LEN, _NPOLICY_BASES, _N_RESCUED) = range(_N_COUNTERS)

#: Reads longer than this get one informational line. The scan is O(L^2), so a
#: long-read FASTQ fed here by accident is slow rather than wrong, and saying so once
#: is the difference between "diagnosable" and "looks hung".
_LONG_READ_NOTICE = 1024


def _new_acc():
    return ([0] * _N_COUNTERS, [], [], [])


def _fold(counters, len_hist, olen_hist, insert_hist, acc):
    """Fold one chunk's statistics into the accumulator, in place.

    The histograms are uncapped and both backends return them with trailing zero bins
    dropped, so a chunk carrying a longer read than anything seen before extends the
    accumulator. Bin *i* is the count of value *i* in every one of them.
    """
    ac, al, ao, ai = acc
    for i in range(_MAX_READ_LEN):
        ac[i] += counters[i]
    if counters[_MAX_READ_LEN] > ac[_MAX_READ_LEN]:      # a maximum, not a sum
        ac[_MAX_READ_LEN] = counters[_MAX_READ_LEN]
    # Everything past the maximum is a plain sum again. Adding a counter without
    # extending this loop silently reports zero -- which is how the first version of the
    # N-policy counters read 0 on input that definitely had no-calls.
    for i in range(_MAX_READ_LEN + 1, _N_COUNTERS):
        ac[i] += counters[i]
    for src, dst in ((len_hist, al), (olen_hist, ao), (insert_hist, ai)):
        if len(src) > len(dst):
            dst.extend([0] * (len(src) - len(dst)))
        for i, c in enumerate(src):
            if c:
                dst[i] += c


#: Bytes per read from the decompressor. Modest on purpose: each read is overlapped
#: with merge work by a prefetch thread, so what matters is keeping the pipeline full
#: rather than amortising syscalls. Measured, 4 MiB blocks cost 45% against 256 KiB
#: because a large blocking read starves the workers while it completes.
_READ_BLOCK = 256 << 10


class _RawStream:
    """Raw byte reader over one (optionally gzipped) FASTQ.

    Holds an immutable ``bytes`` buffer and a read offset rather than slicing: chunks
    are handed to the backend as ``(buf, start, end)``, so a chunk costs no copy and
    worker threads can share one buffer safely (it is immutable, and rebinding
    ``self.buf`` on a refill leaves the old object alive for whoever still holds it).
    Only a refill copies, and then only the unconsumed tail.
    """

    def __init__(self, path, threads):
        import concurrent.futures as cf
        # The prefetch thread below is what makes ISA-L the right choice here; see
        # `_open_binary_read` and `zna._gzip.prefer_isal`.
        self._stream, self._proc = _open_binary_read(str(path), threads,
                                                     own_read_thread=True)
        self._path = str(path)
        self.buf = b""
        self.pos = 0
        self.eof = False
        # One prefetch thread per stream, so a blocking read overlaps with the merge
        # work instead of stalling it. `read` releases the GIL, so this is real overlap.
        # Without it the main thread alternates read-then-merge and the workers starve:
        # measured 1.71 -> 1.32 us/pair at --threads 2.
        self._pool = cf.ThreadPoolExecutor(max_workers=1,
                                           thread_name_prefix="zna-merge-read")
        self._pending = self._pool.submit(self._stream.read, _READ_BLOCK)

    @property
    def avail(self):
        return len(self.buf) - self.pos

    def fill(self, target):
        """Ensure at least *target* unconsumed bytes are buffered, if input remains."""
        if self.avail >= target or self.eof:
            return self.avail > 0
        blocks = [self.buf[self.pos:]] if self.pos else [self.buf]
        got = len(blocks[0])
        self.pos = 0
        while got < target and not self.eof:
            block = self._pending.result()
            if not block:
                self.eof = True
                self._pending = None
                break
            # Start the next read before touching this one, so it runs under the merge.
            self._pending = self._pool.submit(self._stream.read, _READ_BLOCK)
            blocks.append(block)
            got += len(block)
        self.buf = blocks[0] if len(blocks) == 1 else b"".join(blocks)
        return self.avail > 0

    def close(self):
        if self._pending is not None:
            try:
                self._pending.result()
            except Exception:
                pass
            self._pending = None
        self._pool.shutdown(wait=True)
        try:
            self._stream.close()
        except Exception:
            pass
        if self._proc is not None:
            self._proc.wait()
            # Positive exit = a real pigz error. Negative = killed by a signal (e.g.
            # -13 SIGPIPE when the consumer stopped early), which is not an error.
            if self._proc.returncode and self._proc.returncode > 0:
                raise IOError(
                    f"pigz failed reading {self._path} (exit {self._proc.returncode})")


def _run_merge(args, params):
    """Read, merge and write the whole input. Returns the accumulator."""
    backend = _backend.active()
    acc = _new_acc()
    kwargs = (params.match_q, params.step_q, params.t_merge_q, params.t_trim_q,
              params.min_read_length, DISAGREE_Q, not args.no_sync_check)
    # `base_index` is appended per call by the drivers; the policy code and seed
    # follow it.

    # Bytes to keep buffered per stream. Chunks are cut to whole records, so this only
    # has to be comfortably larger than one chunk's worth of them.
    target = max(_READ_BLOCK, args.chunk_size * 1024)

    # Both streams constructed inside the try, closed in a nested finally:
    # each owns a pigz child and a prefetch thread, and an unreadable in2 (or
    # an in1 whose close() raises on a corrupt .gz) used to leak in the other.
    r1 = r2 = None
    try:
        r1 = _RawStream(args.in1, 1)
        r2 = _RawStream(args.in2, 1)
        with FastqWriter(args.out, threads=args.io_threads,
                         level=args.compress_level) as w:
            if args.threads > 1:
                _drive_threaded(backend, r1, r2, w, acc, kwargs, target,
                                args.threads, args.chunk_size, args.quiet,
                                params.npolicy_code, params.rng_seed)
            else:
                _drive_serial(backend, r1, r2, w, acc, kwargs, target, args.quiet,
                              params.npolicy_code, params.rng_seed)
        # Both streams must run out together. A non-empty leftover here is the failure
        # the audit's prototype for this shipped silently: R1 ending first left the
        # trailing R2 records unread, and the desync check cannot see records that were
        # never read.
        _check_drained(backend, r1, "R1", "R2")
        _check_drained(backend, r2, "R2", "R1")
    finally:
        try:
            if r1 is not None:
                r1.close()
        finally:
            if r2 is not None:
                r2.close()
    return acc


def _check_drained(backend, stream, which, other):
    """Raise if *stream* has anything left after the merge.

    Distinguishes the two ways that happens, because they mean different things and a
    wrong message sends the next person to the wrong file: leftover bytes that do not
    form a complete record are a TRUNCATED input, while leftover whole records mean the
    other stream ran out first.
    """
    if not stream.buf[stream.pos:].strip():
        return
    _off, count = backend.split_records(stream.buf, stream.pos, 1)
    if count == 0:
        raise InputError(f"truncated FASTQ record at the end of {which}")
    raise InputError(f"{other} exhausted before {which} (unequal read counts)")


def _drive_serial(backend, r1, r2, w, acc, kwargs, target, quiet, npolicy=1,
                  rng_seed=42):
    while True:
        r1.fill(target)
        r2.fill(target)
        if not r1.avail or not r2.avail:
            break
        blob, c1, c2, counters, lh, oh, ih = backend.merge_chunk(
            r1.buf, r1.pos, len(r1.buf), r2.buf, r2.pos, len(r2.buf),
            *kwargs, acc[0][_N_PAIRS], npolicy, rng_seed)
        if not c1 and not c2:
            break                      # neither stream holds a complete record
        w.write_raw(blob)
        _fold(counters, lh, oh, ih, acc)
        r1.pos += c1
        r2.pos += c2
        if not quiet and acc[0][_N_PAIRS] % 5_000_000 < 100_000:
            logger.info("processed %d pairs", acc[0][_N_PAIRS])


def _drive_threaded(backend, r1, r2, w, acc, kwargs, target, n_threads, chunk_size,
                    quiet, npolicy=1, rng_seed=42):
    """Fan chunks out to worker threads, writing results in SUBMISSION order.

    Ordered output makes the file a pure function of the input and the parameters, so
    `zna merge` produces the same bytes at any thread count -- which is what a corpus
    tool should do, and lets the tests compare whole files across `--threads`. The cost
    is head-of-line blocking bounded by the variance in per-chunk compute time, which for
    fixed-size chunks of Illumina reads is a few percent.

    Submission is windowed so the input is streamed rather than read into memory.
    """
    import concurrent.futures as cf
    from collections import deque

    pending = deque()
    window = max(2, 2 * n_threads)

    def drain(upto):
        while len(pending) > upto:
            blob, c1, c2, counters, lh, oh, ih = pending.popleft().result()
            w.write_raw(blob)
            _fold(counters, lh, oh, ih, acc)

    with cf.ThreadPoolExecutor(max_workers=n_threads) as pool:
        submitted_pairs = 0
        while True:
            r1.fill(target)
            r2.fill(target)
            if not r1.avail or not r2.avail:
                break
            o1, k1 = backend.split_records(r1.buf, r1.pos, chunk_size)
            o2, k2 = backend.split_records(r2.buf, r2.pos, chunk_size)
            k = k1 if k1 < k2 else k2
            if k == 0:
                break                  # neither stream holds a complete record
            if k1 != k:
                o1, _ = backend.split_records(r1.buf, r1.pos, k)
            if k2 != k:
                o2, _ = backend.split_records(r2.buf, r2.pos, k)
            # No slicing: the workers read [pos, o) of a buffer they share.
            pending.append(pool.submit(backend.merge_chunk,
                                       r1.buf, r1.pos, o1, r2.buf, r2.pos, o2,
                                       *kwargs, submitted_pairs, npolicy, rng_seed))
            r1.pos, r2.pos = o1, o2
            submitted_pairs += k
            drain(window)
            if not quiet and submitted_pairs % 5_000_000 < chunk_size:
                logger.info("processed %d pairs", submitted_pairs)
        drain(0)


def _assemble_stats(acc, params, elapsed=None, inflate="unknown"):
    counters, hist, ohist, ihist = acc
    n_pairs = counters[_N_PAIRS]
    merged, trimmed, kept = counters[_MERGED], counters[_TRIMMED], counters[_KEPT]
    n_emitted, n_dropped = counters[_EMITTED], counters[_DROPPED]
    bases_trimmed, frags_short = counters[_BASES_TRIMMED], counters[_FRAGS_SHORT]
    bases_consensus, trim_guard = counters[_BASES_CONSENSUS], counters[_TRIM_GUARD]
    sum_olen, sum_diff = counters[_SUM_OLEN], counters[_SUM_DIFF]
    max_read_len = counters[_MAX_READ_LEN]
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
        # Which kernel ran. The reference backend is ~50x slower and silently correct,
        # so a run that quietly fell back to it looks like a slow node, not a mistake.
        "backend": _backend.active_name(),
        # Which decompressor fed the run. Inflate is the largest single cost of a merge,
        # so a wall-clock number in this dict is not comparable across runs without it.
        # Passed in rather than re-derived: this function is also called by `zna encode
        # --merge-pairs`, which reaches the reader by a different route.
        "inflate": inflate,
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
        # What the N policy did. Reported unconditionally, because the failure this
        # guards against is silent: one dark cycle can make a policy eat most of a
        # library while the run still finishes with "Done".
        "npolicy": params.npolicy,
        "n_rescued_from_mate": counters[_N_RESCUED],
        "npolicy_bases": counters[_NPOLICY_BASES],
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
        # Longest input read. There is no read-length limit -- buffers size themselves
        # -- but the scan is O(L^2), so this is what explains an unexpectedly slow run.
        "max_read_length": max_read_len,
        "params": {
            "threshold_merge_bits": params.t_merge,
            "threshold_trim_bits": params.t_trim,
            "err_rate": params.err_rate,
            "match_bits": round(match_w, 4),
            "mismatch_bits": round(-mismatch_w, 4),
            "min_read_length": params.min_read_length,
            # The exact integers the scan actually used. Recorded so a corpus can be
            # audited against them rather than against a float that was re-derived
            # somewhere else with a different libm. See zna/merge/params.py.
            "score_scale": SCALE,
            "match_q": params.match_q,
            "step_q": params.step_q,
            "threshold_merge_q": params.t_merge_q,
            "threshold_trim_q": params.t_trim_q,
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
        # Only the upper bound needs publishing: the lower one is `min_read_length`,
        # already in `params` above, and a second copy of it here was one more thing to
        # keep in step for no information.
        "insert_size_censoring": {
            "min_mergeable_overlap": -(-params.t_merge_q // params.match_q),
        },
    }
    if elapsed is not None:
        stats["elapsed_s"] = round(elapsed, 1)
        # A cohort field that detects, for free, nodes where the compiled backend or
        # pigz was missing — both of which are silently correct and much slower.
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
    if args.chunk_size < 1:
        raise SystemExit("--chunk-size must be >= 1")
    if args.threads < 1:
        raise SystemExit("--threads must be >= 1")
    if args.io_threads < 1:
        raise SystemExit("--io-threads must be >= 1")


def inflate_backend_for(path1, path2) -> str:
    """Which decompressor will inflate this run's two inputs.

    Both `zna merge` and `zna encode --merge-pairs` read through :class:`_RawStream`,
    which prefetches on its own thread, so both pass ``own_read_thread=True``; see
    :mod:`zna._gzip`. Returns one name when the two inputs agree and ``"a/b"`` when they
    do not -- one gzipped input and one plain is unusual enough to be worth seeing rather
    than hiding behind a single label.
    """
    n1 = _gzip.inflate_backend_name(str(path1), own_read_thread=True)
    n2 = _gzip.inflate_backend_name(str(path2), own_read_thread=True)
    return n1 if n1 == n2 else f"{n1}/{n2}"


def run(args) -> dict:
    """Execute the merge and return the statistics dict."""
    _backend.use(getattr(args, "backend", "auto"))
    # Name the decompressor that will actually run, not the one that is installed. This
    # line said "pigz" whenever the binary existed, which stopped being true when the
    # reader gained an ISA-L path -- and inflate is the largest single cost of a merge,
    # so a throughput number is not interpretable without it.
    inflate = inflate_backend_for(args.in1, args.in2)
    logger.info("backend: %s | inflate: %s | threads: %d",
                _backend.active_name(), inflate, max(1, args.threads))
    params = MergeParams(
        t_merge=args.t_merge,
        t_trim=args.t_trim,
        min_read_length=args.min_read_length,
        npolicy=getattr(args, "npolicy", "trim3"),
        rng_seed=getattr(args, "seed", 42),
    )
    _validate(args, params)

    t0 = time.perf_counter()
    try:
        acc = _run_merge(args, params)
    except InputError as e:                       # desync, or a malformed record
        raise SystemExit(str(e))
    elapsed = time.perf_counter() - t0

    # An empty input is otherwise a silent success all the way down: rc=0 here, then a
    # 22-byte 0-record .zna, and a library disappears from the corpus with every stage
    # green. Cheaper to fail here than to find the hole in a trained model.
    if acc[0][_N_PAIRS] == 0 and not args.allow_empty:
        raise SystemExit(
            f"no read pairs in {args.in1} / {args.in2}. If that is expected, pass "
            f"--allow-empty; otherwise the input is truncated or the wrong file."
        )
    if acc[0][_MAX_READ_LEN] > _LONG_READ_NOTICE:
        logger.info(
            "longest read %d bp: the overlap scan is O(L^2), so expect it to be slow "
            "in proportion (no limit is imposed)", acc[0][_MAX_READ_LEN])
    stats = _assemble_stats(acc, params, elapsed, inflate=inflate)

    if args.json:
        with open(args.json, "w") as fh:
            json.dump(stats, fh, indent=2)
    if not args.quiet:
        logger.info(
            "done: %d pairs -> %d merged, %d trimmed, %d kept (%d records, %d dropped)",
            stats["input_pairs"], stats["merged"], stats["trimmed_pairs"],
            stats["kept_pairs"], stats["emitted_records"], stats["dropped_below_min_length"],
        )
        # Always say what the N policy did. The failure this guards against is silent:
        # a single dark cycle can make a policy consume most of a library while the run
        # still ends with "done".
        n_bases = stats["npolicy_bases"]
        logger.info(
            "no-calls: %d rescued from the mate; --npolicy %s then %s %d base%s",
            stats["n_rescued_from_mate"], stats["npolicy"],
            "removed" if stats["npolicy"] == "trim3" else "substituted",
            n_bases, "" if n_bases == 1 else "s",
        )
        total_bases = sum(i * c for i, c in enumerate(acc[1]))
        if total_bases and n_bases / total_bases > 0.01:
            logger.warning(
                "--npolicy %s affected %.1f%% of emitted bases. That is high enough to "
                "be a run problem (a dark cycle, a failed tile) rather than ordinary "
                "no-calls -- check the input before using this library.",
                stats["npolicy"], 100.0 * n_bases / total_bases,
            )
    return stats


def _require_backend(args) -> None:
    """Refuse to start on the reference kernel unless it was asked for by name.

    The Python backend is correct and ~50x slower. It exists to be an oracle, not a
    fallback, and a *silently correct* 50x slowdown is the worst shape a failure can
    take at cluster scale: the job does not fail, it just looks like a slow node, and it
    burns the whole allocation before anyone looks. So the command line refuses; the
    library entry point (:func:`run`) does not.
    """
    if getattr(args, "backend", "auto") != "auto":
        return
    from .backend import available_merge_backends
    if "accel" in available_merge_backends():
        return
    raise SystemExit(
        "the compiled merge backend is not available, so the scan would run as pure "
        "Python: correct, but about 50x slower. At cluster scale that is "
        "indistinguishable from a slow node. Reinstall zna with a working C++ "
        "toolchain, or pass --backend python if you really mean it."
    )


def run_command(args) -> int:
    """Entry point for ``zna merge`` and for ``python -m zna.merge``."""
    logging.basicConfig(
        level=logging.WARNING if args.quiet else logging.INFO,
        format="[zna merge] %(message)s",
        stream=sys.stderr,
    )
    _require_backend(args)
    run(args)
    return 0


def main(argv=None) -> int:
    return run_command(build_parser().parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
