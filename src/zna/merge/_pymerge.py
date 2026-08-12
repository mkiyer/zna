"""The reference merge backend: the scan, in readable Python.

This is the **oracle** the accelerated backend is defined to agree with, not a fallback
for when the fast one is missing. It is never deleted and never optimised at the cost of
clarity — see :mod:`zna.merge.backend`.

It is JIT-compiled by numba when numba is installed, and runs as identical pure Python
when it is not (correct, ~50x slower). Both are the same source; the ``njit`` decorator
below degrades to a no-op. When the C++ backend lands, numba goes away entirely and this
becomes plain Python, which is what a reference implementation should have been all
along.

Scores are integers throughout — see :mod:`zna.merge.params` for why, and
``docs/MERGE_CPP_DESIGN.md`` §5 for the argmax total order the visiting order realises.
"""
from __future__ import annotations

from .fastqio import InputError
from .names import base_name, strip_pair_suffix

try:  # numba JIT-compiles the scan; without it the same code runs as pure Python.
    from numba import njit
    HAVE_NUMBA = True
except ImportError:  # pragma: no cover - exercised only when numba is missing
    HAVE_NUMBA = False

    def njit(*args, **kwargs):
        """No-op stand-in for numba.njit supporting both @njit and @njit(...)."""
        if len(args) == 1 and callable(args[0]) and not kwargs:
            return args[0]

        def _decorate(func):
            return func

        return _decorate


# Sentinel for "this shift cannot beat the incumbent". A finite int64 keeps numba's
# typing simple and stays far below any reachable score.
_REJECT = -(1 << 62)

# Complement table. A/C/G/T/N in both cases; everything else passes through
# UNCOMPLEMENTED -- deliberate, so an IUPAC ambiguity code survives as itself and the
# kernel's N-vs-N semantics are unchanged: rc(b"RYKMSWBDHVN") == b"NVHDBWSMKYR".
_COMPLEMENT = bytes.maketrans(b"ACGTNacgtn", b"TGCANtgcan")


def reverse_complement(seq: bytes) -> bytes:
    """Reverse-complement a nucleotide sequence (bytes in, bytes out)."""
    return seq.translate(_COMPLEMENT)[::-1]


@njit(cache=True)
def _shift_score(s1, s2rc, s, n, match_q, step_q, best):
    """Score one candidate shift. Returns ``(score_q, mismatches)``.

    ``step_q = match_q + mismatch_q`` is the score given up per mismatch relative to an
    all-match overlap, so ``score_q = n * match_q - d * step_q``. That is monotone in
    ``d`` alone, so the largest mismatch count that could still beat ``best`` is known
    up front and the loop bails the moment it is exceeded. Returns ``_REJECT`` when the
    shift cannot beat ``best``.

    ``dmax`` is the largest ``d`` that can still *strictly* beat ``best``, and in
    integers it is exact: ``score_q > best`` iff ``d * step_q < ceiling - best`` iff
    ``d <= (ceiling - best - 1) // step_q``. The float version of this line truncated,
    which let a shift that could only tie survive to be rejected one comparison later —
    same answer, more work, and the one place a float decided control flow.

    Mismatches are accumulated **branchlessly in blocks of 8**, with the bail tested
    once per block. On a wrong shift the comparison is a coin flip, so a per-position
    ``if`` mispredicts constantly; hoisting the branch out of the block is worth 3.7x
    on the whole scan (measured) and costs at most 7 extra comparisons per rejected
    shift. The result is bit-for-bit identical either way: overshooting ``dmax`` inside
    a block still rejects, and a surviving shift has its exact ``d``.
    """
    ceiling = n * match_q
    if ceiling <= best:
        return _REJECT, 0
    dmax = (ceiling - best - 1) // step_q
    i1 = s if s > 0 else 0          # overlap start in R1
    i2 = -s if s < 0 else 0         # ...and in revcomp(R2)
    d = 0
    k = 0
    lim = n - 7
    while k < lim:
        d += ((s1[i1 + k] != s2rc[i2 + k])
              + (s1[i1 + k + 1] != s2rc[i2 + k + 1])
              + (s1[i1 + k + 2] != s2rc[i2 + k + 2])
              + (s1[i1 + k + 3] != s2rc[i2 + k + 3])
              + (s1[i1 + k + 4] != s2rc[i2 + k + 4])
              + (s1[i1 + k + 5] != s2rc[i2 + k + 5])
              + (s1[i1 + k + 6] != s2rc[i2 + k + 6])
              + (s1[i1 + k + 7] != s2rc[i2 + k + 7]))
        if d > dmax:
            return _REJECT, d
        k += 8
    while k < n:
        d += s1[i1 + k] != s2rc[i2 + k]
        k += 1
    if d > dmax:
        return _REJECT, d
    return ceiling - d * step_q, d


@njit(cache=True)
def scan(s1, s2rc, len1, len2, match_q, step_q, floor_q):
    """Best-scoring shift over ``s in [-(len2-1), len1-1]``.

    Returns ``(shift, score_q, overlap_len, mismatches)`` on the signed single axis
    (``shift < 0`` is read-through). ``overlap_len == 0`` means no shift reached
    ``floor_q``. Shifts are visited in decreasing overlap length, so the scan can stop
    outright once the remaining ceiling cannot beat the incumbent.

    The visiting order — plateau first at maximal ``n`` and ascending ``s``, then the
    flanks at decreasing ``n``, read-through side (the smaller ``s``) before the normal
    side — combined with strict ``>`` is what realises the specified argmax order
    (maximise score, then minimise ``s``). Do not reorder these loops.
    """
    best = floor_q - 1             # a score exactly equal to `floor_q` must win
    best_s = 0
    best_n = 0
    best_d = 0
    nmax = len1 if len1 < len2 else len2
    if nmax <= 0:
        return 0, 0, 0, 0

    # Shifts achieving the maximal overlap: a plateau of width |len1 - len2| + 1.
    plo = 0 if len1 >= len2 else len1 - len2
    phi = len1 - len2 if len1 >= len2 else 0
    s = plo
    while s <= phi:
        sc, d = _shift_score(s1, s2rc, s, nmax, match_q, step_q, best)
        if sc > best:
            best = sc
            best_s = s
            best_n = nmax
            best_d = d
        s += 1

    # Then both flanks, in lockstep, at decreasing overlap length.
    n = nmax - 1
    while n > 0:
        if n * match_q <= best:
            break
        s = n - len2                       # read-through flank (s < plo)
        sc, d = _shift_score(s1, s2rc, s, n, match_q, step_q, best)
        if sc > best:
            best = sc
            best_s = s
            best_n = n
            best_d = d
        s = len1 - n                       # normal-overlap flank (s > phi)
        sc, d = _shift_score(s1, s2rc, s, n, match_q, step_q, best)
        if sc > best:
            best = sc
            best_s = s
            best_n = n
            best_d = d
        n -= 1

    if best_n == 0:
        return 0, 0, 0, 0
    return best_s, best, best_n, best_d


# =========================================================================== #
# Level 2: one pair -- consensus, decision, record construction.
#
# The accelerated backend mirrors this exactly (src/zna/merge/merge_core.hpp), and
# tests/test_merge.py compares them record by record.
# =========================================================================== #

MERGED, TRIMMED, KEPT = 0, 1, 2


def _consensus_r1_overlap(s1, q1, s2rc, q2r, s, olen, disagree_q):
    """Resolve overlap disagreements into R1, by posterior. Returns ``(s1, q1, n)``.

    Only R1 is written because the output's overlap always comes from R1 (merge
    concatenates R1 through the overlap; trim keeps R1 and discards R2's copy), so
    editing R2's discarded copy could never reach the corpus. New ``bytes`` if anything
    changed, else the originals unchanged.
    """
    a0 = s if s > 0 else 0        # mirrors the scan's overlap alignment
    b0 = -s if s < 0 else 0
    s1b = q1b = None
    n = 0
    for i in range(olen):
        a = a0 + i
        b = b0 + i
        qa = q1[a]
        qb = q2r[b]
        if s1[a] != s2rc[b]:
            if qb > qa:                       # R2 is the better-supported call
                if s1b is None:
                    s1b = bytearray(s1)
                    q1b = bytearray(q1)
                s1b[a] = s2rc[b]
                q1b[a] = disagree_q[qb * 256 + qa]
                n += 1
            else:                             # R1 wins, but it is contested: derate it
                if s1b is None:
                    s1b = bytearray(s1)
                    q1b = bytearray(q1)
                q1b[a] = disagree_q[qa * 256 + qb]
    if s1b is None:
        return s1, q1, 0
    return bytes(s1b), bytes(q1b), n


def _build_merged(s, s1, q1, s2rc, q2, len1, len2):
    """Build the merged sequence/quality from the **fragment span** (R1-wins overlap).

    The scan infers exactly one quantity -- the signed offset ``s`` of revcomp(R2) on the
    shared axis -- and the fragment is therefore ``[0, L)`` with ``L = s + len2``,
    uniformly, for every geometry. So build from ``L`` directly rather than
    case-analysing the direction: take R1 from its 5' end as far as it reaches, then let
    revcomp(R2) supply whatever of the fragment R1 does not cover.

    A per-direction construction carried an implicit ``len1 >= len2`` assumption and
    truncated 374 of 137,796 merged records (0.271%) on a production library, each one
    stamped IS_FULL_FRAGMENT while missing bases. ``len(seq) == L`` identically here, by
    construction. Truncation required ``len1 < len2`` strictly, which is why every
    equal-length fixture in the suite was blind to it.

    Returns ``(seq, qual, n1, n2)``, the bases contributed by R1 and R2.
    """
    L = s + len2                                    # fragment length
    take1 = len1 if len1 < L else L                 # R1 covers [0, take1)
    take2 = L - take1                               # R2rc covers [take1, L)
    seq = s1[:take1]
    qual = q1[:take1]
    if take2:
        b = take1 - s                               # ...at this index into revcomp(R2)
        seq = seq + s2rc[b:b + take2]
        qual = qual + q2[::-1][b:b + take2]         # reversed only when actually needed
    return seq, qual, take1, take2


def process_pair(h1, s1, q1, h2, s2, q2, match_q, step_q, t_merge_q, t_trim_q,
                 min_read_length, disagree_q):
    """Classify one pair and build its output records.

    Returns ``(records, outcome, n_dropped, score_q, overlap_len, mismatches,
    bases_consensus_changed, trim_guard_fired)``.

    **The decision** is a single ``argmax`` shift read at two thresholds::

        score >= t_merge_q            -> merge (one full-fragment record)
        t_trim_q <= score < t_merge_q -> keep both, trim the redundant overlap off R2's 3'
        score <  t_trim_q             -> keep both unchanged

    Trimming applies only to the normal (``s >= 0``) geometry: in a read-through the
    redundant bases are R2's *5'* fragment copy and its 3' end is adapter, so there is
    nothing sensible to cut -- such a pair is either merged or kept whole.

    **Trim guard:** a trim that would leave R2 below ``min_read_length`` keeps both reads
    *untrimmed* instead, turning a would-be whole-fragment discard into a no-op.

    **Pair integrity:** an unmerged pair is emitted all-or-nothing. A lone surviving mate
    would be encoded as a spurious "single" -- a full molecule with both endpoints --
    corrupting the fragment-end supervision. A merged read is a genuine full molecule and
    is filtered on its own.
    """
    len1, len2 = len(s1), len(s2)
    s2rc = reverse_complement(s2)
    shift, score, olen, diff = scan(s1, s2rc, len1, len2, match_q, step_q, t_trim_q)

    n_consensus = trim_guard = 0
    if diff > 0:
        s1, q1, n_consensus = _consensus_r1_overlap(
            s1, q1, s2rc, q2[::-1], shift, olen, disagree_q)

    lr = min_read_length
    if olen > 0 and score >= t_merge_q:
        seq, qual, n1, n2 = _build_merged(shift, s1, q1, s2rc, q2, len1, len2)
        name = strip_pair_suffix(h1) + b" merged_%d_%d" % (n1, n2)
        cand = [(name, seq, qual)]
        paired, outcome = False, MERGED
    elif olen > 0 and shift >= 0 and score >= t_trim_q and len2 - olen >= lr:
        keep2 = len2 - olen
        cand = [(h1, s1, q1), (h2, s2[:keep2], q2[:keep2])]
        paired, outcome = True, TRIMMED
    else:
        if olen > 0 and shift >= 0 and score >= t_trim_q:
            trim_guard = 1                                  # guard fired
        cand = [(h1, s1, q1), (h2, s2, q2)]
        paired, outcome = True, KEPT

    if paired:
        kept = cand if (len(cand[0][1]) >= lr and len(cand[1][1]) >= lr) else []
    else:
        kept = [r for r in cand if len(r[1]) >= lr]
    return (kept, outcome, len(cand) - len(kept), score, olen, diff,
            n_consensus, trim_guard)


# =========================================================================== #
# Level 3: a slab of raw FASTQ text in, formatted FASTQ text out.
#
# The production path in the accelerated backend; here, the oracle its blob and its
# counters are compared against byte for byte.
#
# `merge_chunk` consumes only WHOLE PAIRS and reports how many bytes it took from each
# stream separately, so the two buffers may carry different leftovers and the caller
# never scans for record boundaries. A partial record at the end of a buffer is not an
# error -- it is simply not consumed, and at EOF the caller checks that both buffers
# came out empty.
# =========================================================================== #

HIST_MAX = 1024


def _next_record(buf, pos, limit, which):
    """Parse one record at *pos*. Returns ``(header, seq, qual, new_pos)`` or None when
    no COMPLETE record remains (not an error: the caller refills and retries)."""
    if pos >= limit:
        return None
    e1 = buf.find(b"\n", pos, limit)
    if e1 < 0:
        return None
    e2 = buf.find(b"\n", e1 + 1, limit)
    if e2 < 0:
        return None
    e3 = buf.find(b"\n", e2 + 1, limit)
    if e3 < 0:
        return None
    e4 = buf.find(b"\n", e3 + 1, limit)
    if e4 < 0:
        return None
    if buf[pos] != 0x40:                       # b'@'
        raise InputError(f"malformed FASTQ header in {which}")
    h = buf[pos + 1:e1].rstrip(b"\r")
    s = buf[e1 + 1:e2].rstrip(b"\r").upper()
    q = buf[e3 + 1:e4].rstrip(b"\r")
    if len(s) != len(q):
        # A file truncated inside its LAST quality line otherwise looks like a complete
        # record and is emitted malformed with a zero exit status.
        raise InputError(f"FASTQ record in {which} has {len(s)} bases but {len(q)} "
                         f"quality scores (truncated or malformed)")
    return h, s, q, e4 + 1


def merge_chunk(buf1, start1, end1, buf2, start2, end2, match_q, step_q, t_merge_q,
                t_trim_q, min_read_length, disagree_q, check_sync, base_index):
    """Merge every whole pair available in both buffers.

    Returns ``(blob, consumed1, consumed2, counters, len_hist, olen_hist,
    insert_hist)``.
    """
    parts = []
    n_pairs = merged = trimmed = kept = emitted = dropped = 0
    bases_trimmed = frags_short = bases_consensus = trim_guard = 0
    sum_olen = sum_diff = 0
    len_hist = [0] * (HIST_MAX + 1)
    olen_hist = [0] * (HIST_MAX + 1)
    insert_hist = [0] * (HIST_MAX + 1)
    pos1, pos2 = start1, start2

    while True:
        a = _next_record(buf1, pos1, end1, "R1")
        if a is None:
            break
        b = _next_record(buf2, pos2, end2, "R2")
        if b is None:
            break
        h1, s1, q1, try1 = a
        h2, s2, q2, try2 = b

        if check_sync and base_name(h1) != base_name(h2):
            raise InputError(
                f"R1/R2 out of sync at pair {base_index + n_pairs + 1}: "
                f"'{base_name(h1).decode('latin-1')}' != "
                f"'{base_name(h2).decode('latin-1')}'")

        (records, outcome, n_dropped, score, olen, diff,
         n_consensus, guard) = process_pair(
            h1, s1, q1, h2, s2, q2, match_q, step_q, t_merge_q, t_trim_q,
            min_read_length, disagree_q)

        n_pairs += 1
        dropped += n_dropped
        bases_consensus += n_consensus
        trim_guard += guard
        if outcome == MERGED:
            merged += 1
        elif outcome == TRIMMED:
            trimmed += 1
            if records:
                bases_trimmed += len(s2) - len(records[1][1])
        else:
            kept += 1
        if not records and outcome != MERGED:
            frags_short += 1
        if olen:
            olen_hist[olen if olen <= HIST_MAX else HIST_MAX] += 1
            sum_olen += olen
            sum_diff += diff
        for header, seq, qual in records:
            parts.append(b"@%b\n%b\n+\n%b\n" % (header, seq, qual))
            emitted += 1
            L = len(seq)
            len_hist[L if L <= HIST_MAX else HIST_MAX] += 1
            if outcome == MERGED:
                insert_hist[L if L <= HIST_MAX else HIST_MAX] += 1
        pos1, pos2 = try1, try2

    counters = (n_pairs, merged, trimmed, kept, emitted, dropped, bases_trimmed,
                frags_short, bases_consensus, trim_guard, sum_olen, sum_diff)
    return (b"".join(parts), pos1 - start1, pos2 - start2, counters,
            len_hist, olen_hist, insert_hist)


def split_records(buf, start, max_records):
    """Byte offset just past *max_records* complete records, and how many were found.

    Lets the caller cut a buffer into whole-record chunks for parallel workers without
    parsing anything itself. A trailing partial record is not counted and its bytes are
    left for the next chunk.
    """
    pos = start
    found = 0
    n = len(buf)
    while max_records <= 0 or found < max_records:
        p = pos
        complete = True
        for _ in range(4):
            nl = buf.find(b"\n", p)
            if nl < 0:
                complete = False
                break
            p = nl + 1
        if not complete:
            break
        pos = p
        found += 1
        if pos >= n:
            break
    return pos, found
