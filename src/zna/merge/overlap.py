"""Overlap detection kernel — the performance-critical core.

One axis, one scan, one score. Write R1 and revcomp(R2) on a common coordinate axis:
revcomp(R2)'s fragment portion always ends at its own end, so its offset relative to
R1's base 0 is ``s = L - len2`` for a fragment of length ``L``. ``s >= 0`` is a normal
overlap, ``s < 0`` is read-through, ``s == 0`` is exact full overlap. There is exactly
one unknown (``s``), so there is exactly one scan.

Each candidate ``s`` is scored as a log-likelihood ratio in **bits**::

    score(s) = matches * log2((1 - e) / 0.25) + mismatches * log2(e / 0.75)

i.e. a matching base is worth ~2 bits (log2 4: the information in agreeing on one of
four bases) and a mismatch costs ~6.2 bits at ``e = 1%``. Both weights fall out of the
error rate; neither is tuned. The decision is ``argmax`` over ``s`` — not fastp's
first-accept — which is what stops a spurious short hit from preempting the real
offset. See ``docs/READ_MERGE_REDESIGN.md``.

**Pruning.** Because ``score = n * match_q - d * step_q`` depends only on the
overlap length ``n`` and the mismatch count ``d``, the best score still reachable inside
a shift depends only on ``d`` — so the per-shift mismatch budget can be computed *once*
from the incumbent best (``_shift_score``), and the inner loop is a plain compare-and-
count with an early bail, exactly as before. Shifts are visited in decreasing ``n``, so
once ``n * match_q`` cannot beat the incumbent the whole scan terminates.

**The scan is exactly reproducible.** Scores are integers in the fixed-point scale of
:mod:`zna.merge.params` — no float takes part in a comparison, a bail bound or the
argmax — and the argmax is a specified total order, not an artifact of iteration order:

    maximise ``score``; among ties, minimise ``s``.

Shifts are visited in decreasing overlap length and, within that, ascending ``s``, and
improvement is strict ``>``, which realises exactly that order. Ties can only arise
between shifts of *equal* overlap length: a tie across different ``n`` would need
``dn * match_q == dd * step_q``, whose minimal solution is ``step_q / gcd(match_q,
step_q)`` and that gcd is 1, so ``dn`` would have to exceed 1.3e8. See
``docs/MERGE_CPP_DESIGN.md`` §5.

A single ``@njit`` kernel operates directly on ``bytes`` (indexing yields ints). When
numba is installed it JIT-compiles to C speed; when absent the identical function runs
as pure Python (correct, ~50x slower — install numba for production). Operating on
``bytes`` avoids per-pair ``np.frombuffer`` and is ~2.5x faster than the ndarray path.
"""
from __future__ import annotations

from .params import (  # noqa: F401  (score_weights/threshold_bits are re-exported API)
    MergeParams, P_NULL, SCALE, score_weights, threshold_bits, to_bits, to_q,
)

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


# Complement table for A/C/G/T/N (both cases). Bytes outside that set are reversed but
# NOT complemented — `maketrans` passes anything unlisted through — so an IUPAC ambiguity
# code survives as itself: rc(b"RYKMSWBDHVN") == b"NVHDBWSMKYR". Deliberately left alone:
# mapping them to N would change the kernel's N-vs-N scoring semantics (an N pair
# currently earns a full +match_q), and the exposure is nil — of 167,784 real pairs from
# a production BAM, every non-ACGT byte was already N.
_COMPLEMENT = bytes.maketrans(b"ACGTNacgtn", b"TGCANtgcan")

# Direction codes returned by the kernel.
NO_OVERLAP = 0
FORWARD = 1
REVERSE = -1

# Sentinel for "this shift cannot beat the incumbent". A finite int64 keeps numba's
# typing simple and stays far below any reachable score.
_REJECT = -(1 << 62)

#: Default parameters, so ``find_overlap(s1, s2rc)`` needs no ceremony.
_DEFAULTS = MergeParams()


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
def _scan(s1, s2rc, len1, len2, match_q, step_q, floor_q):
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


def find_overlap(seq1: bytes, seq2rc: bytes, params: MergeParams = _DEFAULTS):
    """Detect the overlap between R1 (``seq1``) and revcomp(R2) (``seq2rc``).

    Returns ``(direction, shift, overlap_len, diff, score_q)``:

    * ``FORWARD`` (1): normal overlap; ``shift`` is the R1 offset where R2rc begins;
      the fragment length is ``shift + len(R2)``.
    * ``REVERSE`` (-1): read-through; the fragment is ``seq1[:overlap_len]``.
    * ``NO_OVERLAP`` (0): nothing reached ``params.t_trim``; ``score_q`` is 0.

    ``score_q`` is in the fixed-point scale of :mod:`zna.merge.params` — an integer,
    ``SCALE`` units per bit. Use :func:`~zna.merge.params.to_bits` to report it; never
    convert it to make a decision.

    ``params.t_trim`` is the lower of the two decision thresholds: overlaps that cannot
    reach it are of no interest to any caller, so it doubles as the pruning floor.
    """
    s, score_q, olen, diff = _scan(seq1, seq2rc, len(seq1), len(seq2rc),
                                   params.match_q, params.step_q, params.t_trim_q)
    if olen == 0:
        return NO_OVERLAP, 0, 0, 0, 0
    if s >= 0:
        return FORWARD, s, olen, diff, score_q
    return REVERSE, -s, olen, diff, score_q
