"""Per-pair classification and merged-sequence construction.

This module is the public interface; the work happens in the selected backend —
:mod:`zna.merge._pymerge` (the reference oracle) or the accelerated extension. All
sequences are ``bytes``; headers carry no leading ``@`` and no newline.

A pair has two outcomes (``docs/archive/MERGE_ACCURACY_PLAN.md`` §2): **merged** into one
full-fragment record when :func:`zna.merge.overlap.find_overlap`'s verdict is
``merge``, otherwise **kept** as its two mates, unchanged apart from the N policy. 0.5.x
had a third, a trim band that cut the redundant overlap off both mates' 3' ends; it
removed ~0.55 duplicated bases per pair while causing every wrong trim, and is gone.

The emitted overlap comes from R1, but where the two mates *disagree* the base is
resolved by posterior from the two Phred scores — the quality-aware consensus that
removed the last of fastp's tuning knobs. Its table lives in :mod:`params`, built once
in Python and handed to whichever backend runs.
"""
from __future__ import annotations

from typing import NamedTuple

from . import backend as _backend
from .names import base_name  # noqa: F401  (re-exported: callers look for it here)
from .params import DISAGREE_Q, MergeParams  # noqa: F401  (MergeParams re-exported)


class PairOutcome:
    """Decision categories (for statistics)."""
    MERGED = "merged"
    KEPT = "kept"


#: Backends return an integer outcome; this maps it back. Index order must match the
#: MERGED/KEPT constants in _pymerge and merge_core.hpp.
_OUTCOMES = (PairOutcome.MERGED, PairOutcome.KEPT)


class PairResult(NamedTuple):
    """What :func:`process_pair` did with one pair.

    * ``records`` — ``(header, seq, qual)`` bytes triples to emit (after the
      minimum-read-length filter). A MERGED result is one single; a KEPT result is a
      two-record mate pair, emitted all-or-nothing.
    * ``outcome`` — a :class:`PairOutcome` value.
    * ``n_dropped`` — reads removed by the length filter.
    * ``score``, ``olen``, ``diff``, ``shift`` — the alignment the pair was merged
      from: fixed-point score (:mod:`zna.merge.params`), overlap length, mismatches and
      signed shift (fragment length ``shift + len(R2)``). All 0 when no overlap was
      admitted -- including when one was refused as implausible.
    * ``implausible`` — the best alignment reached the floor but failed the
      plausibility gate, so the pair was kept.
    * ``detected_bases``, ``detected_mismatches`` — informative positions and
      informative mismatches of the best alignment whenever it reached the floor,
      merged or refused (0 when nothing did): the overlap as detected, BEFORE the gate.
    * ``detected_overlap_len`` — that alignment's length (0 when nothing reached the
      floor): the ``n`` the gate looked ``dfit[n]`` up at.

    Over a library ``sum(diff) / sum(olen)`` is the post-admission disagreement rate,
    reported as ``overlap_mismatch_rate``, and ``sum(detected_mismatches) /
    sum(detected_bases)`` the pre-gate one ``--error-rate`` is checked against,
    ``detected_overlap_mismatch_rate``; the histogram of ``detected_overlap_len``,
    with that rate, gives ``expected_refused_true_overlap_fraction``.
    """
    records: list
    outcome: str
    n_dropped: int
    score: int
    olen: int
    diff: int
    shift: int
    implausible: bool
    detected_bases: int
    detected_mismatches: int
    detected_overlap_len: int


def process_pair(h1, s1, q1, h2, s2, q2, p: MergeParams, counters=None) -> PairResult:
    """Classify one pair and build output records under ``p``.

    ``counters`` (optional list of two ints) accumulates ``[bases_consensus_changed,
    implausible_refused]``; leave ``None`` to skip counting.
    """
    p.ensure(max(len(s1), len(s2)))
    (records, outcome, n_dropped, shift, score, olen, diff, n_consensus, implausible,
     _npolicy_bases, _n_rescued, det_bases, det_mismatches,
     _rt, det_len) = _backend.active().process_pair(
        h1, s1, q1, h2, s2, q2, *p.kernel_args(), p.min_read_length, DISAGREE_Q,
        p.npolicy_code, p.rng_seed)
    if counters is not None:
        counters[0] += n_consensus
        counters[1] += implausible
    return PairResult(records, _OUTCOMES[outcome], n_dropped, score, olen, diff, shift,
                      bool(implausible), det_bases, det_mismatches, det_len)
