"""Overlap detection: the public interface over a selectable kernel backend.

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

**This module is the public interface, not the kernel.** The scan itself lives in a
backend — :mod:`zna.merge._pymerge` (the reference oracle) or the accelerated
extension — selected by :mod:`zna.merge.backend` and resolved on first use, so importing
this module costs nothing. Backends operate directly on ``bytes`` (indexing yields
ints), which avoids a per-pair ``np.frombuffer`` and is ~2.5x faster than an ndarray
path.
"""


from .backend import get_merge_backend, get_merge_backend_name
from .params import (  # noqa: F401  (score_weights/threshold_bits are re-exported API)
    MergeParams, P_NULL, SCALE, score_weights, threshold_bits, to_bits, to_q,
)

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


def reverse_complement(seq: bytes) -> bytes:
    """Reverse-complement a nucleotide sequence (bytes in, bytes out)."""
    return seq.translate(_COMPLEMENT)[::-1]


#: Default parameters, so ``find_overlap(s1, s2rc)`` needs no ceremony.
_DEFAULTS = MergeParams()

#: The selected backend's scan, resolved on first use so that importing this module
#: costs nothing (the Python backend pulls in numba; the accel one, an extension).
_SCAN = None
_BACKEND_NAME = None


def use_backend(name=None) -> str:
    """Select the scan backend (``"accel"``, ``"python"``, or ``None``/``"auto"``).

    Returns its canonical name. Raises ``ImportError`` if it cannot be loaded.
    """
    global _SCAN, _BACKEND_NAME
    _SCAN = get_merge_backend(name).scan
    _BACKEND_NAME = get_merge_backend_name(name)
    return _BACKEND_NAME


def backend_name() -> str:
    """Canonical name of the backend in use, selecting the default if none is yet."""
    if _SCAN is None:
        use_backend()
    return _BACKEND_NAME


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
    if _SCAN is None:
        use_backend()
    s, score_q, olen, diff = _SCAN(seq1, seq2rc, len(seq1), len(seq2rc),
                                   params.match_q, params.step_q, params.t_trim_q)
    if olen == 0:
        return NO_OVERLAP, 0, 0, 0, 0
    if s >= 0:
        return FORWARD, s, olen, diff, score_q
    return REVERSE, -s, olen, diff, score_q
