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
