"""Per-pair classification and merged/trimmed sequence construction.

Given a paired-end read, decide merge / trim / keep and build the output record(s).
Pure Python (bytes slicing) on top of the overlap kernel in :mod:`overlap`. All
sequences are handled as ``bytes``; headers carry no leading ``@`` and no newline.

The decision reads **one** likelihood-ratio score (in bits, from :mod:`overlap`) at
**two** thresholds, because a wrong merge and a wrong trim cost different amounts: a
wrong merge produces a chimera — actively false sequence — while a wrong trim removes a
few real bases from a read tail. See ``docs/READ_MERGE_REDESIGN.md`` §4.

The emitted overlap comes from R1, but where the two mates *disagree* the base is
resolved by posterior from the two Phred scores (see ``_consensus_r1_overlap``) — this
is the §7 quality-aware consensus the redesign wanted, and it removed the last of
fastp's tuning knobs.
"""
from __future__ import annotations

from .overlap import FORWARD, NO_OVERLAP, REVERSE, find_overlap, reverse_complement
from .params import MergeParams  # noqa: F401  (re-exported: this is where callers look)


# --------------------------------------------------------------------------- #
# overlap consensus
# --------------------------------------------------------------------------- #
#
# Where the two mates overlap they are two independent reads of the same physical base,
# so a disagreement means exactly one of them is wrong. Which one is not a judgement
# call — the sequencer already told us, in the two Phred scores. With p_i = 10**(-Q_i/10)
# and "exactly one is wrong":
#
#     P(R1 is the wrong one) = p1(1-p2) / (p1(1-p2) + p2(1-p1))
#
# so the consensus base is simply the higher-Q call, and the posterior error of that
# call is the expression above — which is *worse* than the winner's own Q, because a
# contested base is less certain than an uncontested one. Both come out of numbers we
# already have; there is nothing to tune.
#
# This replaces fastp's ``BaseCorrector::correctByOverlapAnalysis``, a pair of hard
# cutoffs (act only when R1 <= Q14 and R2 >= Q30, otherwise R1 always wins). Those
# cutoffs leave the whole middle of the (Q1,Q2) plane defaulting to R1-wins, including
# bands where R2 is overwhelmingly the better call: Q11-vs-Q25 and Q25-vs-Q37 are ~5% of
# real overlap mismatches, R1-wins is ~95% wrong on them, and neither trips the gate.
# Measured residual errors per 100k merged bases: 196 (R1-wins), 77.5 (Q14/Q30), 51.1
# (this rule). It also deletes three knobs -- ``correction``, GOOD_QUAL, BAD_QUAL -- which
# is the point of docs/READ_MERGE_REDESIGN.md §7 and §9.
#
# **Quality is recomputed only where the base changed.** The agreement case would also
# justify a *higher* Q (two independent observations), but the merged FASTQ is a temp
# file that only ``zna encode`` reads, and ZNA stores no quality at all — so a full
# per-position recompute would cost throughput to produce a number with no consumer. If
# that ever changes, the agreement case is Q_out = Q1 + Q2 (probabilities multiply) and
# belongs here.

_MAXQ = 94                       # Phred+33 covers chr(33)..chr(126)


def _build_disagree_table():
    """``_DISAGREE_Q[q_win][q_lose]`` = Phred+33 quality of the winning call.

    Indexed by *raw byte value* so the inner loop is two subscripts and no arithmetic.
    """
    import math
    tbl = [bytearray(_MAXQ + 33) for _ in range(_MAXQ + 33)]
    for qw in range(33, _MAXQ + 33):
        for ql in range(33, _MAXQ + 33):
            pw = 10.0 ** (-(qw - 33) / 10.0)
            pl = 10.0 ** (-(ql - 33) / 10.0)
            num = pw * (1.0 - pl)
            den = num + pl * (1.0 - pw)
            post = 0.5 if den <= 0.0 else num / den
            post = min(max(post, 1e-10), 0.9999)
            q = int(round(-10.0 * math.log10(post)))
            tbl[qw][ql] = min(max(q + 33, 33), _MAXQ + 32)
    return tbl


_DISAGREE_Q = _build_disagree_table()


def _consensus_r1_overlap(s1, q1, s2rc, q2r, direction, shift, olen):
    """Resolve overlap disagreements into R1, by posterior. Returns ``(s1, q1, n)``.

    Only R1 is written because the output's overlap always comes from R1 (merge
    concatenates R1 through the overlap; trim keeps R1 and discards R2's copy), so
    editing R2's discarded copy could never reach the corpus. New ``bytes`` if anything
    changed, else the originals unchanged.
    """
    # Overlap alignment mirrors the scan in overlap._scan:
    #   FORWARD: R1[shift+i] vs R2rc[i];  REVERSE: R1[i] vs R2rc[shift+i]
    a0, b0 = (shift, 0) if direction == FORWARD else (0, shift)
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
                q1b[a] = _DISAGREE_Q[qb][qa]
                n += 1
            else:                             # R1 wins, but it is contested: derate it
                if s1b is None:
                    s1b = bytearray(s1)
                    q1b = bytearray(q1)
                q1b[a] = _DISAGREE_Q[qa][qb]
    if s1b is None:
        return s1, q1, 0
    return bytes(s1b), bytes(q1b), n


class PairOutcome:
    """Decision categories (for statistics)."""
    MERGED = "merged"
    TRIMMED = "trimmed"
    KEPT = "kept"


def _id_end(header: bytes) -> int:
    """Index of the first whitespace (space or tab), or len(header) if none.

    Uses C-level ``bytes.find`` rather than a Python per-char loop (this runs on
    every read).
    """
    sp = header.find(b" ")
    tab = header.find(b"\t")
    if sp == -1:
        return tab if tab != -1 else len(header)
    if tab == -1:
        return sp
    return sp if sp < tab else tab


def base_name(header: bytes) -> bytes:
    """Read ID without the ``/1``/``/2`` pair suffix and without comments/tags.

    Mirrors ZNA's ``get_base_name`` so the R1/R2 sync check matches ZNA's pairing.
    """
    idtok = header[:_id_end(header)]
    slash = idtok.rfind(b"/")
    return idtok[:slash] if slash != -1 else idtok


def _strip_pair_suffix(header: bytes) -> bytes:
    """Drop a trailing ``/1`` or ``/2`` from the ID token, preserving tags/separator."""
    cut = _id_end(header)
    idtok, rest = header[:cut], header[cut:]
    if idtok.endswith(b"/1") or idtok.endswith(b"/2"):
        idtok = idtok[:-2]
    return idtok + rest


def _build_merged(direction, shift, s1, q1, s2rc, q2, len1, len2):
    """Build the merged sequence/quality from the **fragment span** (R1-wins overlap).

    The scan infers exactly one quantity — the offset ``s`` of revcomp(R2) on the shared
    axis — and the fragment is therefore ``[0, L)`` with ``L = s + len2``, uniformly, for
    every geometry. So build from ``L`` directly rather than case-analysing the
    direction: take R1 from its 5' end as far as it reaches, then let revcomp(R2) supply
    whatever of the fragment R1 does not cover.

    This replaces a per-direction construction (a port of fastp's
    ``OverlapAnalysis::merge``) that only ever emitted R2's tail on the forward path with
    ``shift > 0``. That carried an implicit ``len1 >= len2`` assumption, so whenever R1
    was quality-trimmed shorter than the fragment — forward with ``shift == 0`` and
    ``len2 > len1``, or **any** read-through with ``len1 < L`` — the merged read was
    truncated to R1 while ``--treat-unpaired-as-merged`` still declared both of its edges
    true fragment boundaries. Measured on 167,784 real pairs: 374 of 137,796 merged
    records (0.271%), 307 of them via read-through. Truncation required ``len1 < len2``
    strictly, which is why every equal-length fixture in the suite was blind to it.

    Returns ``(seq, qual, n1, n2)`` where ``n1``/``n2`` are the bases contributed by R1
    and R2 respectively (used to build fastp's ``merged_<n1>_<n2>`` name suffix).
    ``len(seq) == L`` identically, by construction.
    """
    s = shift if direction == FORWARD else -shift   # signed offset of R2rc on the axis
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


def process_pair(h1, s1, q1, h2, s2, q2, p: MergeParams, counters=None):
    """Classify one pair and build output records.

    ``counters`` (optional list of two ints) accumulates ``[bases_corrected,
    trim_guard_kept]``; leave ``None`` to skip counting.

    Returns ``(records, outcome, n_dropped, score, olen, diff)``:

    * ``records`` — list of ``(header, seq, qual)`` bytes triples to emit (after the
      minimum-read-length filter). A MERGED result is one single; a KEPT/TRIMMED result
      is a two-record mate pair.
    * ``outcome`` — a :class:`PairOutcome` value (the classification).
    * ``n_dropped`` — reads removed by the length filter.
    * ``score`` — the winning shift's log-likelihood ratio, as a fixed-point integer
      (0 if nothing reached ``t_trim``). See :mod:`zna.merge.params`.
    * ``olen``/``diff`` — the winning overlap's length and mismatch count (0 if none).
      Their ratio, accumulated over a library, is a direct calibration check on
      ``err_rate``: it should sit near it, and a chance-alignment regime drives it up.

    **The decision** is a single ``argmax`` shift read at two thresholds::

        score >= t_merge_q            -> merge (one full-fragment record)
        t_trim_q <= score < t_merge_q -> keep both, trim the redundant overlap off R2's 3'
        score <  t_trim_q             -> keep both unchanged

    Trimming applies only to the normal (``FORWARD``) geometry: in a read-through the
    redundant bases are R2's *5'* fragment copy and its 3' end is adapter, so there is
    nothing sensible to cut — such a pair is either merged or kept whole.

    **Trim guard** (redesign §8): a trim that would leave R2 below ``min_read_length``
    keeps both reads *untrimmed* instead, turning a would-be whole-fragment discard into
    a no-op. This also subsumes the old ``keep2 <= 0`` case (R2 entirely inside the
    overlap), which now only merges — via ``score >= t_merge`` — or is left alone.

    **Pair integrity:** a merged read is a genuine full molecule and is filtered on its
    own. For an *unmerged* pair, the two reads must stay together — if we dropped only
    one short mate, the survivor would become a lone read that ZNA encodes as a spurious
    "single" (a full molecule with both endpoints), corrupting the pretraining labels.
    So the pair is emitted **all-or-nothing**: if either mate is below ``min_read_length``
    the whole fragment is discarded (matching the upstream fastp fragment-discard policy).
    """
    len1, len2 = len(s1), len(s2)
    s2rc = reverse_complement(s2)
    direction, shift, olen, diff, score = find_overlap(s1, s2rc, p)
    # Resolve overlap disagreements by posterior. Nothing to do when the overlap is
    # clean, so this costs nothing on the majority of pairs. Applies in the trim band
    # too: there R2's overlap copy is discarded, so folding its evidence into R1 is the
    # only way that evidence reaches the corpus at all.
    if diff > 0:
        s1, q1, _nc = _consensus_r1_overlap(s1, q1, s2rc, q2[::-1], direction, shift, olen)
        if counters is not None:
            counters[0] += _nc
    lr = p.min_read_length

    if direction != NO_OVERLAP and score >= p.t_merge_q:
        # Enough evidence to assert the fragment: collapse to one read.
        seq, qual, n1, n2 = _build_merged(direction, shift, s1, q1, s2rc, q2, len1, len2)
        # fastp-style merged name: "<id> merged_<n1>_<n2>" (pair suffix stripped, tags
        # preserved). ZNA ignores it; a merged read is a single with both endpoints.
        name = _strip_pair_suffix(h1) + b" merged_%d_%d" % (n1, n2)
        cand = [(name, seq, qual)]
        paired, outcome = False, PairOutcome.MERGED
    elif direction == FORWARD and score >= p.t_trim_q and len2 - olen >= lr:
        # Real overlap, but not enough evidence to risk a chimera: keep both reads and
        # trim the redundant ``olen`` bp off R2's 3' end so it is not counted twice.
        # R1 (full) + R2[:keep2] tile the fragment exactly once. Only R2's 3' end is
        # ever cut, so base 0 of each mate stays a true fragment boundary.
        keep2 = len2 - olen
        cand = [(h1, s1, q1), (h2, s2[:keep2], q2[:keep2])]
        paired, outcome = True, PairOutcome.TRIMMED
    else:
        # No detectable overlap, an unmergeable read-through, or a trim blocked by the
        # guard: keep both reads exactly as they are.
        if counters is not None and direction == FORWARD and score >= p.t_trim_q:
            counters[1] += 1                                # trim guard fired
        cand = [(h1, s1, q1), (h2, s2, q2)]
        paired, outcome = True, PairOutcome.KEPT

    if paired:
        # All-or-nothing: keep the mate pair only if BOTH reads pass (no lone reads).
        kept = cand if (len(cand[0][1]) >= lr and len(cand[1][1]) >= lr) else []
    else:
        # Merged single: filter on its own (drops tiny read-through inserts).
        kept = [r for r in cand if len(r[1]) >= lr]
    return kept, outcome, len(cand) - len(kept), score, olen, diff
