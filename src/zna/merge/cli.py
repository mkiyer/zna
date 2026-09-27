"""CLI for zna read-merge (the ``zna merge`` subcommand).

Reads two positionally-synced FASTQ files, merges or keeps each pair, and writes one
mixed interleaved FASTQ stream for ``zna encode --interleaved`` (see
``docs/METHODS.md``). Also invocable as ``python -m zna.merge``.

Every parameter a decision uses is fixed before the first pair is read
(:func:`params_from_args`): ``--alpha``, ``--error-rate`` and ``--adapter-trimmed``, and
what :mod:`zna.merge.params` derives from them. The input is only ever *checked* against
them, by two run-level diagnostics that ride along in the chunk loop and never change a
decision (``docs/archive/MERGE_ACCURACY_PLAN.md`` §4): the disagreement rate of every overlap
the scan detects, and from it the share of true overlaps the plausibility test is
expected to refuse at ``--error-rate``, and the read-through share of the first
:data:`~zna.merge.params.READTHROUGH_CHECK_PAIRS` pairs, against ``--adapter-trimmed``.
``zna encode --merge-pairs`` reports and warns through the same functions
(:func:`log_policy`, :func:`run_warnings`, :func:`_assemble_stats`), so the two commands
cannot disagree about a library.

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
from decimal import ROUND_CEILING, Decimal, localcontext
from fractions import Fraction

from .. import _gzip
from . import backend as _backend
from . import params as _params
from .args import add_merge_arguments, add_merge_parser, build_parser  # noqa: F401
from .fastqio import FastqWriter, InputError, _open_binary_read
from .params import DISAGREE_Q, SCALE, decimal_str, score_weights
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

#: Counter fields, in the order every backend returns them (``_pymerge._Tally``).
_N_COUNTERS = 16
(N_PAIRS, MERGED, KEPT, EMITTED, DROPPED, FRAGS_SHORT, BASES_CONSENSUS,
 IMPLAUSIBLE, SUM_OLEN, SUM_DIFF, MAX_READ_LEN, NPOLICY_BASES,
 N_RESCUED, DET_BASES, DET_MISMATCHES, RT_STRONG) = range(_N_COUNTERS)

#: Reads longer than this get one informational line. The scan is O(L^2), so a
#: long-read FASTQ fed here by accident is slow rather than wrong, and saying so once
#: is the difference between "diagnosable" and "looks hung".
_LONG_READ_NOTICE = 1024


def _new_acc():
    """``(counters, len_hist, olen_hist, insert_hist, det_olen_hist)``."""
    return ([0] * _N_COUNTERS, [], [], [], [])


def _fold(counters, len_hist, olen_hist, insert_hist, det_olen_hist, acc):
    """Fold one chunk's statistics into the accumulator, in place.

    The histograms are uncapped and both backends return them with trailing zero bins
    dropped, so a chunk carrying a longer read than anything seen before extends the
    accumulator. Bin *i* is the count of value *i* in every one of them.
    """
    ac, al, ao, ai, ad = acc
    for i in range(MAX_READ_LEN):
        ac[i] += counters[i]
    if counters[MAX_READ_LEN] > ac[MAX_READ_LEN]:      # a maximum, not a sum
        ac[MAX_READ_LEN] = counters[MAX_READ_LEN]
    # Everything past the maximum is a plain sum again. Adding a counter without
    # extending this loop silently reports zero -- which is how the first version of the
    # N-policy counters read 0 on input that definitely had no-calls.
    for i in range(MAX_READ_LEN + 1, _N_COUNTERS):
        ac[i] += counters[i]
    for src, dst in ((len_hist, al), (olen_hist, ao), (insert_hist, ai),
                     (det_olen_hist, ad)):
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


# --------------------------------------------------------------------------- #
# the diagnostics: what the run says about its own parameters
# --------------------------------------------------------------------------- #

#: Measured share of pairs whose unrestricted best alignment is a strong read-through:
#: ~0.1% where ``--adapter-trimmed`` is honest, 2-23% where it is false. 1% separates
#: them. Exact, like the other comparisons against a threshold.
_WARN_READTHROUGH = Fraction(1, 100)


#: The ``--error-rate`` warning fires when the plausibility test is expected to refuse
#: more than this share of the run's true overlaps: one in a thousand. A documented
#: definition, not a fit. For scale, at the default ``e``: the test's own promise is
#: ``alpha`` (1e-6) at the declared rate; a library at the rate the default was measured
#: on (~0.009 detected) expects 1.3e-7; on the dev panel's 46 benches only the extreme
#: 5% 3'-ramp set crosses it (0.0355 detected, 0.40% expected), and the next highest
#: expects 0.017% (R2 at 10x the error, 0.0228 detected). The trigger it replaced,
#: "detected rate above ``e``", fired on 12 of the 46, ten of them expecting under
#: 0.02%. Exact, like the other thresholds.
_WARN_REFUSED = Fraction(1, 1000)


def refusal_probability(n: int, dfit_n: int, rate: Fraction) -> Decimal:
    """``P(Binom(n, rate) > dfit_n)``: the chance that a TRUE overlap of *n* bases
    whose mates disagree at *rate* per base carries more mismatches than the
    plausibility test allows, and is refused.

    The upper tail, summed directly rather than as one minus the lower tail, so a tiny
    probability keeps its digits instead of cancelling against 1. Terms are built by the
    pmf recurrence ``p(k+1) = p(k) (n-k)/(k+1) r/(1-r)`` from ``p(0) = (1-r)^n``, in
    :mod:`decimal` at the 50 digits the policy tables use (software, the same result on
    every platform), and the sum stops past the mode once a term no longer changes it.
    A diagnostic, never a decision.
    """
    if dfit_n >= n or rate <= 0:
        return Decimal(0)
    if rate >= 1:
        return Decimal(1)
    with localcontext(_params._CTX):
        r = Decimal(rate.numerator) / Decimal(rate.denominator)
        q = 1 - r
        ratio = r / q
        mode = rate * n
        pmf = q ** n
        tail = Decimal(0)
        for k in range(n + 1):
            if k > dfit_n:
                grown = tail + pmf
                if grown == tail and k > mode:
                    break
                tail = grown
            pmf = pmf * (n - k) / (k + 1) * ratio
        return tail


def expected_refused_fraction(rate: Fraction, det_olen_hist, dfit_table) -> Decimal:
    """The share of the run's detected overlaps the plausibility test is expected to
    refuse if every one of them were TRUE and disagreed at *rate*::

        sum_n  count[n] * P(Binom(n, rate) > dfit[n])  /  sum_n count[n]

    over the histogram of detected overlap lengths (bin *n* = overlaps of *n* bases).
    *rate* is the run's detected disagreement, before the gate; *dfit_table* must cover
    the longest detected overlap. Counts may be weights (floats): the panel evaluation
    projects a weighted sample through this same function. 0 when nothing was detected.

    What it answers is the question ``--error-rate`` poses -- how many true overlaps a
    too-low setting costs -- in the currency of the plan's §8 closed form, rather than
    "is the detected rate above the setting", which fires on every clean library whose
    repeats nudge the rate past it (hg38-noisy at full scale: 0.0051 detected at 0.0041 true) and says
    nothing about how much is lost. It inherits the detected rate's biases in both
    directions (see ``detected_overlap_mismatch_rate`` in :func:`_assemble_stats`), and
    it treats every detected overlap as true, which repeats are not.
    """
    total = sum(det_olen_hist)
    if not total:
        return Decimal(0)
    with localcontext(_params._CTX):
        acc = Decimal(0)
        for n, c in enumerate(det_olen_hist):
            if c:
                acc += Decimal(c) * refusal_probability(n, dfit_table[n], rate)
        return acc / Decimal(total)


def run_refused_fraction(acc, params: MergeParams) -> Decimal:
    """:func:`expected_refused_fraction` for a finished run, at its exact detected
    rate (the reported ``detected_overlap_mismatch_rate`` is this rounded)."""
    c, det_hist = acc[0], acc[4]
    if not c[DET_BASES] or not det_hist:
        return Decimal(0)
    params.ensure(len(det_hist) - 1)
    return expected_refused_fraction(Fraction(c[DET_MISMATCHES], c[DET_BASES]),
                                     det_hist, params.dfit_table)


def readthrough_check_pairs(acc) -> int:
    """How many pairs the ``--adapter-trimmed`` check covered: the input's first
    :data:`~zna.merge.params.READTHROUGH_CHECK_PAIRS`, or all of them."""
    return min(acc[0][N_PAIRS], _params.READTHROUGH_CHECK_PAIRS)


def _suggest_error_rate(rate: Fraction):
    """*rate* as a value to pass back: two significant figures, rounded UP, so the
    suggestion is never below what the data showed ("0.0137..." -> "0.014"). None when
    that is not a value ``--error-rate`` accepts (``>= 0.75``: overlaps that disagree
    that much agree no better than chance, and no setting describes them)."""
    d = _params._CTX.divide(Decimal(rate.numerator), Decimal(rate.denominator))
    q = d.quantize(Decimal(1).scaleb(d.adjusted() - 1), rounding=ROUND_CEILING)
    if q >= Decimal("0.75"):
        return None
    return format(q.normalize(), "f")


def min_mergeable_overlap(acc, params: MergeParams) -> int:
    """The shortest CLEAN overlap that reaches the floor for a pair of the run's longest
    reads, ``ceil(T(2 * max_read_len - 1) / match_bits)`` (0 on an empty run). Shorter
    reads need slightly less, as the floor rises with the number of shifts."""
    longest = acc[0][MAX_READ_LEN]
    if not longest:
        return 0
    return -(-params.t_q(2 * longest - 1) // params.match_q)


def run_warnings(acc, params: MergeParams) -> list[str]:
    """What the finished run says about its own parameters, one message each.

    Diagnostics, never decisions: all are computed alongside the merge and none changed a
    pair. Called at the END of a run by both ``zna merge`` and ``zna encode
    --merge-pairs``, because the disagreement rate is over the whole input.

    **The ``--error-rate`` warning** fires when the plausibility test is expected to
    refuse more than :data:`_WARN_REFUSED` of the run's true overlaps
    (:func:`run_refused_fraction`), and it states the trade rather than just the fix.
    Raising the rate is what stops the gate refusing true overlaps, but the same number
    sets the score's mismatch cost, so it also merges more short or divergent overlaps,
    false ones included -- and on the dev panel's 3'-ramp sets that cost outweighed the
    gain (MERGE_ACCURACY_PLAN.md §3). The refused count is in the message because it
    bounds what raising the rate can recover.
    """
    c = acc[0]
    out = []
    det_d, det_n = c[DET_MISMATCHES], c[DET_BASES]
    refused = run_refused_fraction(acc, params)
    if Fraction(refused) > _WARN_REFUSED:
        rate = Fraction(det_d, det_n)
        suggestion = _suggest_error_rate(rate)
        if suggestion is None:
            advice = ("No --error-rate describes overlaps that disagree this much (it "
                      "must be below 0.75): these reads do not behave like overlapping "
                      "mates.")
        elif rate > params.e:
            advice = (
                f"Rerun with --error-rate {suggestion} (the detected rate, rounded up) "
                f"to recover them -- but the rate also sets how little a mismatch costs "
                f"in the score, so more short or divergent overlaps merge too, false "
                f"ones included (a simulated 3'-degraded library: 240 wrong merges at "
                f"0.01, 745 at 0.036). Raise it only if the refused pairs matter more "
                f"than that.")
        else:
            # Only reachable with a large --alpha: at a rate within the setting the
            # test refuses a true overlap with probability below alpha, by construction.
            advice = (
                f"The detected rate is within --error-rate, so this is --alpha "
                f"{decimal_str(params.alpha_exact)} itself: the test refuses up to that "
                f"share of true overlaps by design. A smaller --alpha refuses fewer.")
        out.append(
            f"the overlaps this run detected disagree at {float(rate):.4g} (informative "
            f"mismatches per compared base, {det_d} in {det_n}, before the plausibility "
            f"test), so at --error-rate {decimal_str(params.e)} the test is expected to "
            f"refuse {100 * float(refused):.2g}% of true overlaps (warned above "
            f"{100 * float(_WARN_REFUSED):g}%). A refused pair is kept whole, never "
            f"merged wrongly; this run refused {c[IMPLAUSIBLE]} as implausible, true "
            f"overlaps and repeats together. {advice}")
    checked = readthrough_check_pairs(acc)
    if (params.adapter_trimmed and checked
            and Fraction(c[RT_STRONG], checked) > _WARN_READTHROUGH):
        out.append(
            f"--adapter-trimmed was declared, but {c[RT_STRONG] / checked:.1%} of the "
            f"first {checked} pairs align best as a strong read-through (an honest "
            f"declaration measures ~0.1%): the reads appear to contain adapter. The "
            f"declaration disables read-through merges; if these reads are raw, drop the "
            f"flag.")
    longest = c[MAX_READ_LEN]
    need = min_mergeable_overlap(acc, params)
    if need > longest:
        out.append(
            f"no pair in this run could merge: at --alpha "
            f"{decimal_str(params.alpha_exact)} and --error-rate {decimal_str(params.e)} "
            f"even a perfect overlap must be {need} bases long, and the longest read is "
            f"{longest}. Every pair was kept whole.")
    return out


def log_policy(params: MergeParams, emit) -> None:
    """Report the run's policy through *emit(level, message)*, before the first pair."""
    emit(logging.INFO,
         f"alpha {decimal_str(params.alpha_exact)}, error rate "
         f"{decimal_str(params.e)} (--error-rate), adapter-trimmed "
         f"{'declared' if params.adapter_trimmed else 'not declared'}")


# --------------------------------------------------------------------------- #
# the chunk loop
#
# **Table capacity.** A chunk stops in front of a pair whose read is longer than the
# policy tables cover and returns `need`; the driver grows the tables (doubling) and
# resumes from where the chunk stopped. The tables' prefix does not change when they
# grow, so the output is identical whether or not, and wherever, that happens. They
# start at 256 bases, so a 2x150 library never regrows.
# --------------------------------------------------------------------------- #

def _run_merge(args, params):
    """Read, merge and write the whole input. Returns the accumulator."""
    backend = _backend.active()
    acc = _new_acc()
    check_sync = not args.no_sync_check

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
                _drive_threaded(backend, r1, r2, w, acc, params, check_sync,
                                target, args.threads, args.chunk_size, args.quiet)
            else:
                _drive_serial(backend, r1, r2, w, acc, params, check_sync,
                              target, args.quiet)
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


def _merge_args(params, check_sync, buf1, s1, e1, buf2, s2, e2, base):
    """``merge_chunk``'s arguments, with the tables as they stand NOW.

    Resolved on the main thread, at submission: the tables are only ever grown there, so
    a worker never builds or reads a table while another thread replaces it. ``base``
    numbers the chunk's first pair in the whole input, which is what places the
    read-through check on the same pairs at any chunking.
    """
    return (buf1, s1, e1, buf2, s2, e2, *params.kernel_args(), params.min_read_length,
            DISAGREE_Q, check_sync, base, params.npolicy_code, params.rng_seed,
            _params.READTHROUGH_CHECK_PAIRS)


def _drive_serial(backend, r1, r2, w, acc, params, check_sync, target, quiet):
    while True:
        r1.fill(target)
        r2.fill(target)
        if not r1.avail or not r2.avail:
            break
        blob, c1, c2, counters, lh, oh, ih, dh, need = backend.merge_chunk(
            *_merge_args(params, check_sync, r1.buf, r1.pos, len(r1.buf), r2.buf,
                         r2.pos, len(r2.buf), acc[0][N_PAIRS]))
        w.write_raw(blob)
        _fold(counters, lh, oh, ih, dh, acc)
        r1.pos += c1
        r2.pos += c2
        if need:
            params.ensure(need)        # and resume in front of the pair that needed it
            continue
        if not c1 and not c2:
            break                      # neither stream holds a complete record
        if not quiet and acc[0][N_PAIRS] % 5_000_000 < 100_000:
            logger.info("processed %d pairs", acc[0][N_PAIRS])


def _merge_job(merge_chunk, args, job):
    """One threaded chunk: ``(merge_chunk's result, job if it stopped for capacity)``.

    *args* are resolved by the caller, on the main thread (see `_merge_args`). *job*
    comes back only when the driver must resume the chunk, so a finished chunk pins
    no input buffer.
    """
    res = merge_chunk(*args)
    return res, (job if res[-1] else None)


def _drive_threaded(backend, r1, r2, w, acc, params, check_sync, target, n_threads,
                    chunk_size, quiet):
    """Fan chunks out to worker threads, writing results in SUBMISSION order.

    Ordered output makes the file a pure function of the input and the parameters, so
    `zna merge` produces the same bytes at any thread count -- which is what a corpus
    tool should do, and lets the tests compare whole files across `--threads`. The cost
    is head-of-line blocking bounded by the variance in per-chunk compute time, which for
    fixed-size chunks of Illumina reads is a few percent.

    Submission is windowed so the input is streamed rather than read into memory. A
    chunk that stopped short for table capacity is finished here, in order, on the
    main thread -- the only thread that grows the tables.

    **Only a stopped chunk keeps its input.** ``pending`` holds futures and nothing
    else, and a worker hands its buffers back only when it stopped for capacity (see
    `_merge_job`). `_RawStream.fill` makes a new ~2 MB buffer about once per chunk, so a
    window that held each chunk's buffers pinned up to 2x``window`` of them: peak
    footprint on 1M chr22 pairs went 30 / 69 / 194 MB at 1 / 4 / 16 threads, against
    0.5.3's flat 30-36 MB, whose finished futures had already dropped their arguments.
    """
    import concurrent.futures as cf
    from collections import deque

    pending = deque()
    window = max(2, 2 * n_threads)

    def drain(upto):
        while len(pending) > upto:
            (blob, c1, c2, counters, lh, oh, ih, dh, need), job = \
                pending.popleft().result()
            if need:
                buf1, s1, e1, buf2, s2, e2, base = job
            while True:
                w.write_raw(blob)
                _fold(counters, lh, oh, ih, dh, acc)
                if not need:
                    break
                params.ensure(need)
                s1, s2, base = s1 + c1, s2 + c2, base + counters[N_PAIRS]
                blob, c1, c2, counters, lh, oh, ih, dh, need = backend.merge_chunk(
                    *_merge_args(params, check_sync, buf1, s1, e1, buf2, s2, e2, base))

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
            job = (r1.buf, r1.pos, o1, r2.buf, r2.pos, o2, submitted_pairs)
            pending.append(pool.submit(_merge_job, backend.merge_chunk,
                                       _merge_args(params, check_sync, *job), job))
            r1.pos, r2.pos = o1, o2
            submitted_pairs += k
            drain(window)
            if not quiet and submitted_pairs % 5_000_000 < chunk_size:
                logger.info("processed %d pairs", submitted_pairs)
        drain(0)


# --------------------------------------------------------------------------- #
# statistics
# --------------------------------------------------------------------------- #

def _assemble_stats(acc, params: MergeParams, elapsed=None, inflate="unknown"):
    """The run's statistics. Every value is finite and of a fixed type whatever the
    input -- hulkrna's cohort gather rejects ``Infinity`` and type changes -- so an empty
    input reports zeros, never ``NaN``, and ratios are floats even when integral."""
    counters, hist, ohist, ihist, dhist = acc
    n_pairs = counters[N_PAIRS]
    merged, kept = counters[MERGED], counters[KEPT]
    n_emitted, n_dropped = counters[EMITTED], counters[DROPPED]
    sum_olen, sum_diff = counters[SUM_OLEN], counters[SUM_DIFF]
    max_read_len = counters[MAX_READ_LEN]
    det_n, det_d = counters[DET_BASES], counters[DET_MISMATCHES]
    rt_pairs = readthrough_check_pairs(acc)
    total_bases = sum(i * c for i, c in enumerate(hist))
    pct = (lambda n: round(100.0 * n / n_pairs, 3) if n_pairs else 0.0)
    match_w, mismatch_w = score_weights(params.e)
    min_overlap = min_mergeable_overlap(acc, params)
    refused = run_refused_fraction(acc, params)
    import zna
    stats = {
        # Provenance: the question every future corpus defect opens with is "which
        # build produced this file?". Config values are already cohort-queryable via
        # gather's pipeline tool; only code identity was missing.
        "tool": "zna-merge",
        "tool_version": zna.__version__,
        "policy": _params.POLICY,
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
        "kept_pairs": kept,
        "merged_pct": pct(merged),
        "kept_pct": pct(kept),
        "emitted_records": n_emitted,
        "dropped_below_min_length": n_dropped,
        "fragments_dropped_short_mate": counters[FRAGS_SHORT],
        "bases_consensus_changed": counters[BASES_CONSENSUS],
        # Pairs whose best alignment reached T but whose mismatches are implausible as
        # sequencing error at alpha: repeats, kept whole instead of merged.
        "implausible_refused": counters[IMPLAUSIBLE],
        # What the N policy did. Reported unconditionally, because the failure this
        # guards against is silent: one dark cycle can make a policy eat most of a
        # library while the run still finishes with "Done".
        "npolicy": params.npolicy,
        "n_rescued_from_mate": counters[N_RESCUED],
        "npolicy_bases": counters[NPOLICY_BASES],
        "mean_emitted_length": round(total_bases / n_emitted, 1) if n_emitted else 0.0,
        # The error rate every decision used: --error-rate, exactly (the prologue's
        # merge record carries it as a decimal string).
        "error_rate": float(params.e),
        # ...and what the data say about it: informative mismatches per informative
        # compared base over EVERY overlap the scan detected (score >= T), before the
        # plausibility gate -- merged and refused alike. It is a check, not an
        # estimate, and it errs both ways: repeats that reached T inflate it on clean
        # libraries (hg38-noisy at full scale: 0.0051 detected against 0.0041 true), while on a degraded one
        # the most divergent true overlaps never reach T, so it can read below the truth
        # (a 3'-ramp set: 0.0355 against 0.0385).
        "detected_overlap_mismatch_rate": round(det_d / det_n, 6) if det_n else 0.0,
        "detected_overlap_bases": det_n,
        # What that rate costs at `error_rate`: the share of the detected overlaps the
        # gate is expected to refuse were they all true, sum_n count[n] * P(Binom(n,
        # rate) > dfit[n]) / sum_n count[n] over the histogram below
        # (expected_refused_fraction). run_warnings warns above 0.001. 4 significant
        # figures: at the default it is ~1e-7 on a library that fits it.
        "expected_refused_true_overlap_fraction": float(format(refused, ".4g")),
        # Every detected overlap's length, merged and refused alike -- the n the gate
        # looked dfit[n] up at. With the rate above it reproduces the fraction.
        "detected_overlap_length_histogram": {str(i): c for i, c in enumerate(dhist)
                                               if c},
        "adapter_trimmed": bool(params.adapter_trimmed),
        # Share of the first `readthrough_check_pairs` pairs whose unrestricted best
        # alignment is a read-through reaching T: ~0.1% on honestly trimmed input,
        # 2-23% on raw reads. Computed with or without the declaration; only a declared
        # run warns on it.
        "readthrough_check_pairs": rt_pairs,
        "readthrough_check_strong_fraction":
            round(counters[RT_STRONG] / rt_pairs, 6) if rt_pairs else 0.0,
        # Mismatches per aligned base over ADMITTED overlaps -- the alignments pairs
        # were merged from, after the plausibility gate -- counting a one-sided N as a
        # mismatch, which the detected rate above does not. Still the sensitive
        # degradation channel: per-base degradation moves it long before merged_pct
        # moves.
        "overlap_mismatch_rate": round(sum_diff / sum_olen, 6) if sum_olen else 0.0,
        "overlap_bases_compared": sum_olen,
        # Longest input read. There is no read-length limit -- buffers and tables size
        # themselves -- but the scan is O(L^2), so this explains an unexpectedly slow run.
        "max_read_length": max_read_len,
        "params": {
            "alpha": float(params.alpha_exact),
            "match_bits": round(match_w, 4),
            "mismatch_bits": round(-mismatch_w, 4),
            "min_read_length": params.min_read_length,
            # The exact integers the scan actually used. Recorded so a corpus can be
            # audited against them rather than against a float that was re-derived
            # somewhere else. The per-pair floor is T_q[N] = to_q(log2(N / alpha)),
            # reproducible from `alpha` alone; see zna/merge/params.py.
            "score_scale": SCALE,
            "match_q": params.match_q,
            "step_q": params.step_q,
        },
        "length_histogram": {str(i): c for i, c in enumerate(hist) if c},
        # Admitted overlap length per merged pair, in its natural quantum (bases). Its
        # short-end cliff reads the floor directly: the shortest clean overlap that can
        # merge is ceil(T / match_bits).
        "overlap_length_histogram": {str(i): c for i, c in enumerate(ohist) if c},
        # Inferred fragment length, merged pairs only (a merged record IS the
        # fragment). CENSORED AT BOTH ENDS, so do not read it as the library's insert
        # distribution without accounting for that: hard-floored at min_read_length
        # (shorter fragments are dropped, not observed) and hard-capped per pair at
        # len1 + len2 - ceil(T(N) / match_bits), beyond which the mates no longer
        # overlap enough to merge.
        "insert_size_histogram": {str(i): c for i, c in enumerate(ihist) if c},
        # The cap for a pair of the run's longest reads (the floor rises with N, so
        # shorter reads need slightly less). The lower bound is `min_read_length`.
        "insert_size_censoring": {
            "min_mergeable_overlap": min_overlap,
            "at_read_length": max_read_len,
        },
    }
    if elapsed is not None:
        stats["elapsed_s"] = round(elapsed, 1)
        # A cohort field that detects, for free, nodes where the compiled backend or
        # pigz was missing — both of which are silently correct and much slower.
        stats["pairs_per_second"] = round(n_pairs / elapsed) if elapsed > 0 else 0
    return stats


def params_from_args(args) -> MergeParams:
    """The :class:`MergeParams` a parsed command line asks for, or a clean exit.

    Shared by ``zna merge`` and ``zna encode --merge-pairs``, whose algorithm flags are
    one definition (:func:`zna.merge.args.add_merge_algorithm_arguments`). The realistic
    failure is not a hand-typed flag, it is a config typo propagating into a
    whole-cohort run, so a value the policy cannot mean is refused rather than clamped.
    """
    try:
        params = MergeParams(
            alpha=getattr(args, "alpha", _params.DEFAULT_ALPHA),
            error_rate=getattr(args, "error_rate", _params.DEFAULT_ERROR_RATE),
            adapter_trimmed=bool(getattr(args, "adapter_trimmed", False)),
            min_read_length=getattr(args, "min_read_length", 40),
            npolicy=getattr(args, "npolicy", None) or "trim3",
            rng_seed=42 if getattr(args, "seed", None) is None else args.seed,
        )
    except ValueError as e:
        # The message names the parameter ("alpha ..." / "error rate ...").
        raise SystemExit(f"invalid merge parameter: {e}") from None
    if params.min_read_length < 1:
        raise SystemExit("--min-read-length must be >= 1")
    return params


def _validate(args) -> None:
    """Reject I/O arguments that cannot mean anything."""
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
    params = params_from_args(args)
    _validate(args)
    log_policy(params, logger.log)

    t0 = time.perf_counter()
    try:
        acc = _run_merge(args, params)
    except InputError as e:                       # desync, or a malformed record
        raise SystemExit(str(e))
    elapsed = time.perf_counter() - t0

    # An empty input is otherwise a silent success all the way down: rc=0 here, then a
    # 22-byte 0-record .zna, and a library disappears from the corpus with every stage
    # green. Cheaper to fail here than to find the hole in a trained model.
    if acc[0][N_PAIRS] == 0 and not args.allow_empty:
        raise SystemExit(
            f"no read pairs in {args.in1} / {args.in2}. If that is expected, pass "
            f"--allow-empty; otherwise the input is truncated or the wrong file."
        )
    if acc[0][MAX_READ_LEN] > _LONG_READ_NOTICE:
        logger.info(
            "longest read %d bp: the overlap scan is O(L^2), so expect it to be slow "
            "in proportion (no limit is imposed)", acc[0][MAX_READ_LEN])
    stats = _assemble_stats(acc, params, elapsed, inflate=inflate)
    # Warnings go out whatever -q says: each one means a parameter may be wrong for
    # this library, which is exactly what a quiet cluster run must not swallow.
    for msg in run_warnings(acc, params):
        logger.warning(msg)

    if args.json:
        with open(args.json, "w") as fh:
            json.dump(stats, fh, indent=2)
    if not args.quiet:
        logger.info(
            "done: %d pairs -> %d merged, %d kept (%d refused as implausible); "
            "%d records, %d dropped",
            stats["input_pairs"], stats["merged"], stats["kept_pairs"],
            stats["implausible_refused"], stats["emitted_records"],
            stats["dropped_below_min_length"],
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
