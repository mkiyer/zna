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
error rate ``e`` (``--error-rate``, default 1%; :mod:`zna.merge.params`); neither is
tuned. The candidate is ``argmax`` over ``s`` — not fastp's first-accept — which is
what stops a spurious short hit from preempting the real offset. See
``docs/METHODS.md``.

**The decision** (``docs/archive/MERGE_ACCURACY_PLAN.md`` §2) reads that argmax ``W`` against
two tests at one tolerance ``alpha``: ``W`` must reach the pair's floor ``T =
log2((len1 + len2 - 1) / alpha)`` -- at most ``alpha`` chance merges of unrelated
sequence per pair -- and its mismatches, not counting positions where exactly one mate
is ``N``, must be plausible as sequencing error, ``d <= dfit[n]`` -- at most ``alpha``
true overlaps refused. An implausible ``W`` is a repeat outscoring the truth; the pair
then has no overlap and nothing is searched for in its place. Under ``adapter_trimmed``
only shifts with ``L >= max(len1, len2)`` are eligible. :func:`find_overlap` returns
that decision; :func:`scan_unrestricted` is the bare argmax, for diagnostics only.

**Pruning.** Because ``score = n * match_q - d * step_q`` depends only on the
overlap length ``n`` and the mismatch count ``d``, the best score still reachable inside
a shift depends only on ``d`` — so the per-shift mismatch budget can be computed *once*
from the incumbent best (``_shift_score``), and the inner loop is a plain compare-and-
count with an early bail, exactly as before. Shifts are visited in decreasing ``n``, so
once ``n * match_q`` cannot beat the incumbent the whole scan terminates. The floor
``T`` is the incumbent the scan starts from, so at ~28 bits it prunes from the first
shift.

**The scan is exactly reproducible.** Scores are integers in the fixed-point scale of
:mod:`zna.merge.params` — no float takes part in a comparison, a bail bound or the
argmax — and the argmax is a specified total order, not an artifact of iteration order:

    maximise ``score``; among ties, minimise ``s``.

Shifts are visited in decreasing overlap length and, within that, ascending ``s``, and
improvement is strict ``>``, which realises exactly that order. Ties can only arise
between shifts of *equal* overlap length: a tie across different ``n`` would need
``dn * match_q == dd * step_q``, whose minimal solution is ``step_q / gcd(match_q,
step_q)``, and for every error rate checked that exceeds 1e7 bases. See
``docs/METHODS.md``.

**This module is the public interface, not the kernel.** The scan itself lives in a
backend — :mod:`zna.merge._pymerge` (the reference oracle) or the accelerated
extension — selected by :mod:`zna.merge.backend` and resolved on first use, so importing
this module costs nothing. Backends operate directly on ``bytes`` (indexing yields
ints), which avoids a per-pair ``np.frombuffer`` and is ~2.5x faster than an ndarray
path.
"""
from typing import NamedTuple

from . import backend as _backend
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

#: Verdicts, as :func:`find_overlap` reports them. The backends use integer codes
#: (``VERDICT_*`` in ``_pymerge``); index order must match.
MERGE = "merge"
NONE = "none"
IMPLAUSIBLE = "implausible"
_VERDICTS = (NONE, MERGE, IMPLAUSIBLE)


class Overlap(NamedTuple):
    """One pair's overlap decision.

    ``verdict`` is ``"merge"`` (the pair is merged from this alignment), ``"none"``
    (no eligible shift reached the pair's floor ``T``) or ``"implausible"`` (the best
    one did, but its mismatches are implausible as sequencing error: a repeat, so the
    pair is kept whole). The alignment fields describe the alignment the verdict is
    about -- the one used for ``"merge"``, the one refused for ``"implausible"`` -- and
    are all zero for ``"none"``.

    ``shift`` is signed on the single axis (``< 0`` is read-through), so the inferred
    fragment is ``fragment_length = shift + len(R2)``. ``score_q`` is fixed point
    (:func:`~zna.merge.params.to_bits` to report it; never convert it to decide).
    ``informative_mismatches`` excludes positions where exactly one base is ``N``.
    """
    verdict: str
    shift: int
    overlap_len: int
    mismatches: int
    informative_mismatches: int
    score_q: int
    fragment_length: int


class Alignment(NamedTuple):
    """A bare argmax: ``overlap_len == 0`` means nothing reached the floor."""
    shift: int
    overlap_len: int
    mismatches: int
    score_q: int


def reverse_complement(seq: bytes) -> bytes:
    """Reverse-complement a nucleotide sequence (bytes in, bytes out)."""
    return seq.translate(_COMPLEMENT)[::-1]


#: Default parameters, so ``find_overlap(s1, s2rc)`` needs no ceremony: ``--alpha`` and
#: ``--error-rate``'s defaults, no declaration.
_DEFAULTS = MergeParams()


def use_backend(name=None) -> str:
    """Select the merge backend (``"accel"``, ``"python"``, or ``None``/``"auto"``).

    Returns its canonical name. Raises ``ImportError`` if it cannot be loaded.
    """
    return _backend.use(name)


def backend_name() -> str:
    """Canonical name of the backend in use, selecting the default if none is yet."""
    return _backend.active_name()


def find_overlap(seq1: bytes, seq2rc: bytes, params: MergeParams = _DEFAULTS) -> Overlap:
    """The authoritative overlap decision for R1 (``seq1``) and revcomp(R2) (``seq2rc``).

    Exactly the decision :func:`zna.merge.pairs.process_pair` acts on -- the same
    backend call -- under ``params``' ``alpha``, error rate and ``adapter_trimmed``.
    """
    len1, len2 = len(seq1), len(seq2rc)
    params.ensure(max(len1, len2))
    verdict, shift, score_q, olen, diff, informative = _backend.active().overlap(
        seq1, seq2rc, len1, len2, *params.kernel_args())
    return Overlap(_VERDICTS[verdict], shift, olen, diff, informative, score_q,
                   shift + len2 if olen else 0)


def scan_unrestricted(seq1: bytes, seq2rc: bytes,
                      params: MergeParams = _DEFAULTS) -> Alignment:
    """DIAGNOSTIC: the best shift over every ``s``, at the pair's floor ``T``.

    No contract, no plausibility gate -- 0.5.3's detection, with 0.6's per-pair floor.
    This is what the ``--adapter-trimmed`` check counts (a read-through here is a read
    that appears to run into adapter), and what an analysis wants when it asks "what
    would the scan have picked". It is not a decision; do not merge from it.
    """
    len1, len2 = len(seq1), len(seq2rc)
    if len1 <= 0 or len2 <= 0:
        return Alignment(0, 0, 0, 0)
    s, score_q, olen, diff = _backend.active().scan(
        seq1, seq2rc, len1, len2, params.match_q, params.step_q,
        params.t_q(len1 + len2 - 1), 0)
    return Alignment(s, olen, diff, score_q)
