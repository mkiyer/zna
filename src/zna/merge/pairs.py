"""Per-pair classification and merged/trimmed sequence construction.

This module is the public interface; the work happens in the selected backend —
:mod:`zna.merge._pymerge` (the reference oracle) or the accelerated extension. All
sequences are ``bytes``; headers carry no leading ``@`` and no newline.

The decision reads **one** likelihood-ratio score (from :mod:`overlap`, in the
fixed-point scale of :mod:`params`) at **two** thresholds, because a wrong merge and a
wrong trim cost different amounts: a wrong merge produces a chimera — actively false
sequence — while a wrong trim removes a few real bases from a read tail. See
``docs/READ_MERGE_REDESIGN.md`` §4.

The emitted overlap comes from R1, but where the two mates *disagree* the base is
resolved by posterior from the two Phred scores — the §7 quality-aware consensus the
redesign wanted, which removed the last of fastp's tuning knobs. Its table lives in
:mod:`params`, built once in Python and handed to whichever backend runs.
"""
from __future__ import annotations

from . import backend as _backend
from .names import base_name  # noqa: F401  (re-exported: callers look for it here)
from .params import DISAGREE_Q, MergeParams  # noqa: F401  (MergeParams re-exported)


class PairOutcome:
    """Decision categories (for statistics)."""
    MERGED = "merged"
    TRIMMED = "trimmed"
    KEPT = "kept"


#: Backends return an integer outcome; this maps it back. Index order must match the
#: MERGED/TRIMMED/KEPT constants in _pymerge and merge_core.hpp.
_OUTCOMES = (PairOutcome.MERGED, PairOutcome.TRIMMED, PairOutcome.KEPT)


def process_pair(h1, s1, q1, h2, s2, q2, p: MergeParams, counters=None):
    """Classify one pair and build output records.

    ``counters`` (optional list of two ints) accumulates ``[bases_consensus_changed,
    trim_guard_kept]``; leave ``None`` to skip counting.

    Returns ``(records, outcome, n_dropped, score, olen, diff)``:

    * ``records`` — list of ``(header, seq, qual)`` bytes triples to emit (after the
      minimum-read-length filter). A MERGED result is one single; a KEPT/TRIMMED result
      is a two-record mate pair, emitted all-or-nothing.
    * ``outcome`` — a :class:`PairOutcome` value.
    * ``n_dropped`` — reads removed by the length filter.
    * ``score`` — the winning shift's log-likelihood ratio as a fixed-point integer
      (0 if nothing reached ``t_trim``). See :mod:`zna.merge.params`.
    * ``olen``/``diff`` — the winning overlap's length and mismatch count (0 if none).
      Their ratio, accumulated over a library, is a direct calibration check on
      ``err_rate``: it should sit near it, and a chance-alignment regime drives it up.
    """
    (records, outcome, n_dropped, score, olen, diff,
     n_consensus, trim_guard) = _backend.active().process_pair(
        h1, s1, q1, h2, s2, q2,
        p.match_q, p.step_q, p.t_merge_q, p.t_trim_q, p.min_read_length, DISAGREE_Q)
    if counters is not None:
        counters[0] += n_consensus
        counters[1] += trim_guard
    return records, _OUTCOMES[outcome], n_dropped, score, olen, diff
