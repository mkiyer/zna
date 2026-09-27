"""The reference merge backend: the scan, the decision, in readable Python.

This is the **oracle** the accelerated backend is defined to agree with, not a fallback
for when the fast one is missing. It is never deleted and never optimised at the cost of
clarity — see :mod:`zna.merge.backend`.

Plain Python: no numba, no JIT, nothing between the reader and the algorithm. It is
~50x slower than the compiled backend, and that is the correct trade for an oracle --
speed here would only buy the ability to be wrong in the same way as the thing it
checks.

Scores are integers throughout — see :mod:`zna.merge.params` for why, and
``docs/METHODS.md`` for the argmax total order the visiting order realises. The policy
the functions below implement is ``docs/archive/MERGE_ACCURACY_PLAN.md`` §2; every parameter
arrives already derived, as integers and two ``int64`` tables (``T_q`` indexed by the
pair's shift count, ``dfit`` by overlap length).
"""
from __future__ import annotations

from .fastqio import InputError
from .names import base_name, strip_pair_suffix

#: The backend contract this module implements. :mod:`zna.merge.backend` refuses a
#: backend whose marker differs, which is what stops a compiled extension built for an
#: older policy (0.5.x has no marker at all) from running with the new arguments.
POLICY_ABI = 2

#: Which mate of the pair an emitted record came from -- the whole geometry
#: transfer of ``zna encode --merge-pairs``.  Mirrors ``Slot`` in
#: ``fastq_chunk.hpp``.
SLOT_MERGED = 0
SLOT_MATE1 = 1
SLOT_MATE2 = 2

# Sentinel for "this shift cannot beat the incumbent"; far below any reachable score.
_REJECT = -(1 << 62)

#: The input bounds of every entry point, mirrored by ``MAX_WEIGHT_Q``/``MAX_FLOOR_Q``
#: in ``_accel.cpp``. Python's integers are unbounded and the kernel's are ``int64``, so
#: past these the two backends would disagree (and the kernel's arithmetic would be
#: undefined): with reads under 2^31 bases, a weight ``<= 2^30`` keeps every score under
#: 2^61, and a floor within ``+-2^61`` keeps ``floor - 1`` above ``_REJECT`` and
#: ``ceiling - best - 1`` under 2^63. Real values are far inside: ``step_q`` is 1.4e8
#: at ``e = 0.01`` and 5.3e8 at the smallest accepted, 1e-9; ``T`` is ~2^29.
MAX_WEIGHT_Q = 1 << 30
MAX_FLOOR_Q = 1 << 61


def _check_weights(match_q, step_q, fn):
    if not (0 < match_q <= MAX_WEIGHT_Q and 0 < step_q <= MAX_WEIGHT_Q):
        raise ValueError(f"{fn}: match_q and step_q must be in [1, 2^30]")


def _check_tables(t_table, dfit_table, fn):
    """``T_q`` needs entries 0 and 1 (one base per mate), ``dfit`` entry 0, and every
    ``T_q`` entry is a floor within ``+-MAX_FLOOR_Q``. Checked once per call, over the
    whole table, as the compiled backend does when it borrows it."""
    if len(t_table) < 2 or len(dfit_table) < 1:
        raise ValueError(f"{fn}: the T and dfit tables are empty")
    if min(t_table) < -MAX_FLOOR_Q or max(t_table) > MAX_FLOOR_Q:
        raise ValueError(f"{fn}: a T table entry is out of range")


def _check_npolicy(npolicy, fn):
    if npolicy not in (0, 1, 2):          # NPOLICY_KEEP, _TRIM3, _RANDOM
        raise ValueError(f"{fn}: npolicy must be 0, 1 or 2")

# Complement table. A/C/G/T/N in both cases; everything else passes through
# UNCOMPLEMENTED -- deliberate, so an IUPAC ambiguity code survives as itself and the
# kernel's N-vs-N semantics are unchanged: rc(b"RYKMSWBDHVN") == b"NVHDBWSMKYR".
_COMPLEMENT = bytes.maketrans(b"ACGTNacgtn", b"TGCANtgcan")

_N = 0x4E                          # ord('N'); the parser upper-cases every sequence


def reverse_complement(seq: bytes) -> bytes:
    """Reverse-complement a nucleotide sequence (bytes in, bytes out)."""
    return seq.translate(_COMPLEMENT)[::-1]


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


def scan(s1, s2rc, len1, len2, match_q, step_q, floor_q, adapter_trimmed=0):
    """Best-scoring eligible shift, with ``floor_q`` as the least score that counts.

    Returns ``(shift, score_q, overlap_len, mismatches)`` on the signed single axis
    (``shift < 0`` is read-through). ``overlap_len == 0`` means no eligible shift reached
    ``floor_q``. Shifts are visited in decreasing overlap length, so the scan can stop
    outright once the remaining ceiling cannot beat the incumbent.

    **Eligible shifts.** All of ``s in [-(len2-1), len1-1]`` by default. Under
    ``adapter_trimmed`` -- the declaration that no read extends past its molecule --
    only ``s >= max(0, len1 - len2)``, i.e. inferred fragment ``L = s + len2 >=
    max(len1, len2)``: a shorter ``L`` would put a read past the fragment's end, which
    the declaration makes impossible. That is not merely ``s >= 0``: with ``len1 >
    len2``, ``0 <= s < len1 - len2`` puts R1 past the end. In the loops below the
    restriction is exactly two things -- the plateau shrinks to its last shift, and the
    read-through flank is never visited (every shift on it is ``< plo <= 0``).

    The visiting order — plateau first at maximal ``n`` and ascending ``s``, then the
    flanks at decreasing ``n``, read-through side (the smaller ``s``) before the normal
    side — combined with strict ``>`` is what realises the specified argmax order
    (maximise score, then minimise ``s``). Restricting eligibility removes shifts from
    that order without reordering the rest, so the same holds over the eligible set. Do
    not reorder these loops.
    """
    _check_weights(match_q, step_q, "scan()")
    if not -MAX_FLOOR_Q <= floor_q <= MAX_FLOOR_Q:
        raise ValueError("scan(): floor_q is out of range")
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
    if adapter_trimmed:
        plo = phi                  # the one plateau shift with L >= max(len1, len2)
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
        if not adapter_trimmed:
            s = n - len2                   # read-through flank (s < plo)
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


def _n_positions(s1, s2rc, s, n):
    """``(one_sided, both)``: positions of the overlap at shift *s* where exactly one
    mate reads ``N``, and where both do.

    Neither says anything about whether the mates agree. A one-sided ``N`` is a
    mismatch in the scan (a no-call never equals a call), so the plausibility gate
    discounts it from the mismatches; ``N`` against ``N`` is a match in the scan, so the
    detected-overlap disagreement rate discounts it from the compared bases (plan §3:
    "informative positions only"). The gate's ``dfit`` stays indexed by the overlap
    length, as §2 specifies. Byte ``N`` only: lower case never reaches here from the
    parser, and IUPAC codes compare as themselves (they carry partial information).
    """
    if _N not in s1 and _N not in s2rc:
        return 0, 0
    i1 = s if s > 0 else 0
    i2 = -s if s < 0 else 0
    one = both = 0
    for k in range(n):
        a = s1[i1 + k] == _N
        b = s2rc[i2 + k] == _N
        one += a != b
        both += a and b
    return one, both


def table_capacity(t_table, dfit_table):
    """The longest read the two tables cover: ``T_q`` needs ``N = len1 + len2 - 1 <
    len(t_table)``, ``dfit`` needs ``n <= len(dfit_table) - 1``."""
    cap = len(t_table) // 2
    dcap = len(dfit_table) - 1
    return cap if cap < dcap else dcap


#: `overlap` verdicts. Mirror the VERDICT_* constants in merge_core.hpp.
VERDICT_NONE, VERDICT_MERGE, VERDICT_IMPLAUSIBLE = 0, 1, 2


def overlap(s1, s2rc, len1, len2, match_q, step_q, t_table, dfit_table,
            adapter_trimmed):
    """The authoritative overlap decision for one pair.

    Returns ``(verdict, shift, score_q, overlap_len, mismatches, informative)``:

    * ``W`` = the best eligible shift (:func:`scan`) with the pair's own floor
      ``T_q[len1 + len2 - 1]``: nothing reaching it is ``VERDICT_NONE``.
    * ``informative`` = ``W``'s mismatches minus the positions where exactly one base is
      ``N`` (:func:`_n_positions`).
    * ``informative > dfit[n_W]`` is ``VERDICT_IMPLAUSIBLE``: that many disagreements
      would happen in a true overlap of this length with probability below ``alpha``,
      so ``W`` is a repeat, not the fragment. The pair has no overlap. Nothing is
      searched for in its place -- a runner-up in a repetitive context re-placed 37-38%
      of caught wrong merges onto *another* wrong shift on the gene-disjoint holdout.
    * otherwise ``VERDICT_MERGE``: ``W`` is the alignment the pair is built from.

    The alignment fields describe ``W`` for both MERGE and IMPLAUSIBLE -- for the latter
    it is the refused alignment, reported for diagnostics; nothing is built from it --
    and are all zero for NONE. The tables must cover ``max(len1, len2)``.
    """
    _check_weights(match_q, step_q, "overlap()")
    _check_tables(t_table, dfit_table, "overlap()")
    return _decide(s1, s2rc, len1, len2, match_q, step_q, t_table, dfit_table,
                   adapter_trimmed)[:6]


def _decide(s1, s2rc, len1, len2, match_q, step_q, t_table, dfit_table,
            adapter_trimmed):
    """:func:`overlap`, plus a seventh field: ``W``'s positions where both mates read
    ``N`` (:func:`_n_positions`), which the detected-overlap diagnostic discounts."""
    if len1 <= 0 or len2 <= 0:
        return VERDICT_NONE, 0, 0, 0, 0, 0, 0
    longest = len1 if len1 > len2 else len2
    if longest > table_capacity(t_table, dfit_table):
        raise ValueError(f"read of {longest} bases exceeds the policy tables "
                         f"(capacity {table_capacity(t_table, dfit_table)})")
    shift, score, olen, diff = scan(s1, s2rc, len1, len2, match_q, step_q,
                                    t_table[len1 + len2 - 1], adapter_trimmed)
    if olen == 0:
        return VERDICT_NONE, 0, 0, 0, 0, 0, 0
    one_sided, both_n = _n_positions(s1, s2rc, shift, olen)
    informative = diff - one_sided
    verdict = VERDICT_IMPLAUSIBLE if informative > dfit_table[olen] else VERDICT_MERGE
    return verdict, shift, score, olen, diff, informative, both_n


# =========================================================================== #
# Level 2: one pair -- consensus, decision, record construction.
#
# The accelerated backend mirrors this exactly (src/zna/merge/merge_core.hpp), and
# tests/test_merge.py compares them record by record.
# =========================================================================== #

#: Pair outcomes: a pair is merged into one record or kept as two. There is no third
#: outcome -- 0.6 removed the trim band (MERGE_ACCURACY_PLAN.md §2).
MERGED, KEPT = 0, 1

#: What to do with a no-call the overlap could not rescue. Same vocabulary as
#: ``zna encode --npolicy``, deliberately: one flag, one meaning, both tools.
NPOLICY_KEEP, NPOLICY_TRIM3, NPOLICY_RANDOM = 0, 1, 2

#: Per-record provenance bits, emitted as the ``ZN:i:<bits>`` header tag. Mirrors the
#: ``PROV_*`` constants in ``merge_core.hpp``; see there for why the byte exists and why
#: there is deliberately no "merged" bit. Bit 1 was ``PROV_TRIMMED`` until 0.6 removed
#: the trim band; it is retired rather than reused, so a set bit in any corpus still
#: means one thing.
PROV_RESCUED, PROV_NTRIMMED, PROV_NSUBBED = 2, 4, 8

_M64 = 0xFFFFFFFFFFFFFFFF
_SUB = b"ACGT"


def _merge_mix64(x):
    """splitmix64's finalizer — the same function as ``merge_mix64`` in merge_core.hpp.

    Substitution is position-derived rather than drawn from a running stream, so it
    cannot depend on how pairs were batched into chunks.
    """
    x = (x + 0x9E3779B97F4A7C15) & _M64
    x = ((x ^ (x >> 30)) * 0xBF58476D1CE4E5B9) & _M64
    x = ((x ^ (x >> 27)) * 0x94D049BB133111EB) & _M64
    return x ^ (x >> 31)


def _sub_n(seq, seed, rec):
    """Replace every ``N`` in *seq* with a base from the seeded, position-derived stream."""
    if b"N" not in seq:
        return seq, 0
    out = bytearray(seq)
    n = 0
    for i, c in enumerate(out):
        if c == 0x4E:
            out[i] = _SUB[_merge_mix64(
                (seed + 0xBF58476D1CE4E5B9 * (rec + 1)
                      + 0x94D049BB133111EB * (i + 1)) & _M64) & 3]
            n += 1
    return bytes(out), n


def _consensus_r1_overlap(s1, q1, s2rc, q2r, s, olen, disagree_q):
    """Resolve overlap disagreements by posterior, into R1's copy of the overlap.

    Returns ``(s1, q1, n, rescued)``: ``n`` counts R1 bases changed and ``rescued`` the
    no-calls among them recovered from R2. New ``bytes`` if anything changed, else the
    originals.

    R1 alone, because the merged record is the only record built from the overlap: it
    takes the overlap from R1 and R2 contributes only outside it, so R2's copy is
    discarded. (0.5.x also wrote R2 on its trim path, where each mate kept part of the
    overlap; the trim path is gone.) A KEPT pair gets no consensus at all -- nothing
    about it depends on the alignment being right.

    The decision is symmetric — the better-supported base by posterior from the two
    Phred scores, with the winner's quality derated because a contested base is less
    certain. On equal quality R1 stands (derated).

    **N rescue.**  An ``N`` carries no base information, so a real call on the other
    mate beats it whatever the two qualities say, and the rescued base keeps the
    surviving mate's own quality rather than a contested-base derating — there was no
    contest.  Without this the rescue happened only by luck, because an instrument
    usually assigns an N a low quality; a *high*-quality N beat a real base and survived
    into the corpus.  Only ``N`` is rescued, not the IUPAC codes, which do carry partial
    information.

    Rescue does not touch the *scan*: an N still counts as a mismatch when the shift is
    scored, so which shift wins is unchanged. It is discounted only where the plausibility
    gate asks whether the mismatches look like sequencing error (:func:`_n_positions`).
    """
    a0 = s if s > 0 else 0        # mirrors the scan's overlap alignment
    b0 = -s if s < 0 else 0
    s1b = q1b = None
    n = rescued = 0
    for i in range(olen):
        a = a0 + i
        b = b0 + i
        if s1[a] != s2rc[b]:
            if s1b is None:
                s1b, q1b = bytearray(s1), bytearray(q1)
            a_is_n = s1[a] == _N
            b_is_n = s2rc[b] == _N
            if a_is_n != b_is_n:              # rescue: a real call beats an N
                if a_is_n:                    # R2 rescues R1
                    s1b[a] = s2rc[b]
                    q1b[a] = q2r[b]
                    n += 1
                    rescued += 1
                # R1 rescuing R2 would write a copy nothing emits: skip it.
            elif a_is_n:                      # both are N: nothing to rescue from
                pass
            elif q2r[b] > q1[a]:              # R2 is the better-supported call
                s1b[a] = s2rc[b]
                q1b[a] = disagree_q[q2r[b] * 256 + q1[a]]
                n += 1
            else:                             # R1 wins, but it is contested: derate it
                q1b[a] = disagree_q[q1[a] * 256 + q2r[b]]
    if s1b is None:
        return s1, q1, 0, 0
    return bytes(s1b), bytes(q1b), n, rescued


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


def _prov_name(header, bits, trim3_n, subn_n, rescued_n):
    """*header* with this record's provenance tokens appended. Mirrors ``build_name``.

    **Tags pass through untouched.** The header is copied verbatim and the tokens are
    *appended*; nothing is removed or rewritten, so ``zna encode --label`` reads the same
    ``KEY:T:VALUE`` tags off an emitted record that it would have read off the input.

    Returns *header* itself when there is nothing to say, which is the common case — and
    on the accelerated side that is what keeps an untouched record zero-copy.

    The colon-less tokens are skipped by ZNA's tag parser, which requires ``KEY:T:VALUE``.
    ``ZN:i:<bits>`` is the one meant to be read, and it is absent when no bit is set, so
    it resolves to 0 through the label machinery's own missing-value path.
    """
    if not bits:
        return header
    out = header + b" ZN:i:%d" % bits
    if trim3_n:
        out += b" trim3_%d" % trim3_n
    if subn_n:
        out += b" subn_%d" % subn_n
    if rescued_n:
        out += b" rescued_%d" % rescued_n
    return out


def _trim3(seq, qual):
    """Cut a read at its first ``N``, keeping ``[0, first_N)``.

    3' only. Base 0 is a true fragment boundary — the read starts at a fragment end and
    runs inward — so cutting from the far end never disturbs it, and the emitted read
    stays honestly anchored however short it gets.

    Returns the original objects untouched when there is no ``N``, so the common case
    allocates nothing.
    """
    k = seq.find(b"N")
    if k < 0:
        return seq, qual, len(seq)
    return seq[:k], qual[:k], k


def _npolicy_mate(seq, qual, npolicy, seed, rec):
    """The N policy on one mate. Returns ``(seq, qual, k)``, ``k`` the bases it touched:
    substituted under ``random`` (``rec`` = ``2 * pair_index + mate``, which is what
    makes substitution position-derived), cut off under ``trim3``, none under ``keep``.
    """
    if npolicy == NPOLICY_RANDOM:
        seq, k = _sub_n(seq, seed, rec)
        return seq, qual, k
    if npolicy == NPOLICY_TRIM3:
        t_seq, t_qual, k = _trim3(seq, qual)
        return t_seq, t_qual, len(seq) - k
    return seq, qual, 0


def process_pair(h1, s1, q1, h2, s2, q2, match_q, step_q, t_table, dfit_table,
                 adapter_trimmed, min_read_length, disagree_q, npolicy=NPOLICY_TRIM3,
                 rng_seed=0, pair_index=0, rt_check=0):
    """Classify one pair and build its output records.

    Returns ``(records, outcome, n_dropped, shift, score_q, overlap_len, mismatches,
    bases_consensus_changed, implausible, npolicy_bases, n_rescued, detected_bases,
    detected_mismatches, readthrough_strong, detected_overlap_len)``, with each record a
    ``(header, seq, qual)`` tuple. ``shift``/``score_q``/``overlap_len``/``mismatches``
    are the alignment the pair was merged from (``overlap_len == 0``: none); a refused
    implausible alignment reports zeros there and ``implausible == 1``. The last four
    are the run's
    diagnostics (see :func:`_process_pair_ex`); none of them affects a decision. The thin
    public shim over :func:`_process_pair_ex`, which additionally carries each record's
    PROV_* byte -- the record adapter reads the bits there directly, with no ``ZN:i:``
    tag round-trip, mirroring ``PairResult::prov`` in the C++ core.
    """
    _check_weights(match_q, step_q, "process_pair()")
    _check_npolicy(npolicy, "process_pair()")
    _check_tables(t_table, dfit_table, "process_pair()")
    records, *rest = _process_pair_ex(
        h1, s1, q1, h2, s2, q2, match_q, step_q, t_table, dfit_table, adapter_trimmed,
        min_read_length, disagree_q, npolicy, rng_seed, pair_index, rt_check)
    return ([r[:3] for r in records], *rest)


def _process_pair_ex(h1, s1, q1, h2, s2, q2, match_q, step_q, t_table, dfit_table,
                     adapter_trimmed, min_read_length, disagree_q,
                     npolicy=NPOLICY_TRIM3, rng_seed=0, pair_index=0, rt_check=0):
    """:func:`process_pair` with records as ``(header, seq, qual, prov)``.

    **The decision** is :func:`overlap`'s verdict, and there are two outcomes::

        verdict MERGE           -> one full-fragment record (R1 wins ties in the
                                   posterior consensus), unless trim3 has cut the mates
                                   so far that they no longer tile the fragment, and
                                   then as below
        verdict NONE/IMPLAUSIBLE -> both mates, unchanged apart from the N policy --
                                   never the consensus, which only a merged record uses

    A merged record shorter than ``min_read_length`` is dropped (it is its fragment, so
    nothing is lost that was not already too short).

    **Two diagnostics ride along; neither changes anything above** (plan §4):

    * ``detected_bases``/``detected_mismatches`` -- informative positions and
      informative mismatches of the best alignment (every position where either mate
      reads ``N`` left out, :func:`_n_positions`) whenever it reached ``T``, i.e. for the
      MERGE *and* the IMPLAUSIBLE verdict: every overlap the scan detected, before the
      gate. Summed over a run they are the disagreement
      rate ``--error-rate`` is checked against. Before the gate, because the gate is
      what a too-low ``e`` makes wrong: a rate measured on its survivors could never
      exceed the ``e`` that filtered them.
    * ``detected_overlap_len`` -- that alignment's length, the ``n`` the gate looked
      ``dfit[n]`` up at (0 when nothing reached ``T``). Its histogram over a run is what
      the expected share of refused true overlaps is summed over
      (:func:`zna.merge.cli.expected_refused_fraction`).
    * ``readthrough_strong`` -- with *rt_check* set, 1 when the UNRESTRICTED best shift
      reaches ``T`` as a read-through (``L = s + len2 < max(len1, len2)``): the
      ``--adapter-trimmed`` check. Without the declaration every shift is eligible, so
      that shift is the scan's own winner and costs nothing; under it, the shifts it
      forbids were never visited and one extra unrestricted scan runs. The caller sets
      *rt_check* for the input's first pairs only.

    **Pair integrity:** an unmerged pair is emitted all-or-nothing. A lone surviving mate
    would be encoded as a spurious "single" -- a full molecule with both endpoints --
    corrupting the fragment-end supervision. A merged read is a genuine full molecule and
    is filtered on its own.
    """
    len1, len2 = len(s1), len(s2)
    s2rc = reverse_complement(s2)
    verdict, shift, score, olen, diff, informative, both_n = _decide(
        s1, s2rc, len1, len2, match_q, step_q, t_table, dfit_table, adapter_trimmed)
    implausible = 1 if verdict == VERDICT_IMPLAUSIBLE else 0
    # Compared positions minus the uninformative ones: the one-sided N positions (all
    # mismatches, so `diff - informative` counts them) and the N-against-N ones.
    det_bases = olen - (diff - informative) - both_n
    det_mismatches = informative
    det_len = olen
    rt_strong = 0
    if rt_check and len1 and len2:
        if adapter_trimmed:
            rs, _rsc, rn, _rd = scan(s1, s2rc, len1, len2, match_q, step_q,
                                     t_table[len1 + len2 - 1], 0)
        else:
            rs, rn = shift, olen
        rt_strong = 1 if rn and rs + len2 < (len1 if len1 > len2 else len2) else 0
    if verdict != VERDICT_MERGE:
        shift = score = olen = diff = 0    # nothing is built from a refused alignment

    lr = min_read_length
    L = shift + len2                       # the inferred fragment length, if merging

    # The consensus is written only into R1, and only on the merge verdict: the merged
    # record takes the overlap from R1. A pair with no admitted overlap is emitted
    # untouched -- an alignment too suspect to merge on is too suspect to rewrite bases
    # on (measured under 0.5.x: of 3,068 kept pairs with a detected overlap, zero had
    # found the true shift, and writing R1 there turned 1,379 correct bases wrong to fix
    # 78).
    s1_in, q1_in = s1, q1
    n_consensus = n_rescued = 0
    if diff > 0:
        s1, q1, n_consensus, n_rescued = _consensus_r1_overlap(
            s1, q1, s2rc, q2[::-1], shift, olen, disagree_q)

    # ---- the N policy, after the rescue, so a no-call the mate could answer costs
    #      nothing. trim3 is 3' only, so both 5' anchors -- the two fragment termini --
    #      are untouched however short the reads get.
    s1p, q1p, npolicy_1 = _npolicy_mate(s1, q1, npolicy, rng_seed, pair_index * 2)
    s2p, q2p, npolicy_2 = _npolicy_mate(s2, q2, npolicy, rng_seed, pair_index * 2 + 1)

    # ---- merge on GEOMETRY, reusing the evidence ------------------------------
    #
    # The pair still tiles the fragment iff len1 + len2 >= L. When it does, the
    # reconstruction IS the fragment, exactly and N-free. Nothing is re-scored: trimming
    # cuts 3' ends, which is where a normal overlap lives, so a re-scan would refuse
    # merges it had ample evidence for a moment earlier. Only trim3 changes a length, so
    # only trim3 can turn a merge verdict into a kept pair here.
    will_merge = olen > 0 and len(s1p) + len(s2p) >= L

    if not will_merge and s1 is not s1_in:
        # A merge verdict that trim3 cut below tiling: the pair is KEPT, and a kept mate
        # is the input with the N policy applied and nothing else (plan §2, and §8's
        # "kept-mate substitutions are zero by construction"). The consensus -- its
        # substitutions, derated qualities and rescues -- existed only to build the
        # merged record, so R1 is re-derived from the input and none of it is counted.
        # Under 0.5.x the rewritten R1 was emitted here (167 kept pairs on the dev
        # panel's N benches, 10 of their substitutions to a wrong base).
        s1p, q1p, npolicy_1 = _npolicy_mate(s1_in, q1_in, npolicy, rng_seed,
                                            pair_index * 2)
        n_consensus = n_rescued = 0

    # The run-level counters are the per-mate ones summed.
    npolicy_bases = npolicy_1 + npolicy_2
    # Which policy bit a touched record earns, and which token carries its count.
    rnd = npolicy == NPOLICY_RANDOM
    npolicy_bit = PROV_NSUBBED if rnd else PROV_NTRIMMED

    if will_merge:
        # `shift` is the offset of revcomp(R2) on the shared axis, so it is tied to R2's
        # length. R2 keeps its 5' anchor at fragment position L-1, so a trimmed mate
        # covers [L - len2', L) and the offset becomes L - len2'. L itself is unchanged
        # -- that is the whole point.
        s2prc = s2rc if s2p is s2 else reverse_complement(s2p)
        seq, qual, n1, n2 = _build_merged(L - len(s2p), s1p, q1p, s2prc, q2p,
                                          len(s1p), len(s2p))
        # A merged record is built from BOTH mates, so its provenance is the pair's: the
        # policy counts are the two summed, and the rescues are R1's, the only ones that
        # reached the emitted bases.
        bits = ((PROV_RESCUED if n_rescued else 0)
                | (npolicy_bit if npolicy_bases else 0))
        # fastp-style merged name, pair suffix stripped, tags preserved. Keeping
        # `merged_<n1>_<n2>` LAST is fastp's convention and costs nothing, so the
        # provenance tokens go before it.
        name = _prov_name(strip_pair_suffix(h1), bits,
                          0 if rnd else npolicy_bases,
                          npolicy_bases if rnd else 0,
                          n_rescued) + b" merged_%d_%d" % (n1, n2)
        cand = [(name, seq, qual, bits)]
        paired, outcome = False, MERGED
    else:
        # A kept record never carries PROV_RESCUED -- only the N policy has touched it.
        b1 = npolicy_bit if npolicy_1 else 0
        b2 = npolicy_bit if npolicy_2 else 0
        cand = [(_prov_name(h1, b1, 0 if rnd else npolicy_1,
                            npolicy_1 if rnd else 0, 0), s1p, q1p, b1),
                (_prov_name(h2, b2, 0 if rnd else npolicy_2,
                            npolicy_2 if rnd else 0, 0), s2p, q2p, b2)]
        paired, outcome = True, KEPT

    if paired:
        kept = cand if (len(cand[0][1]) >= lr and len(cand[1][1]) >= lr) else []
    else:
        kept = [r for r in cand if len(r[1]) >= lr]
    return (kept, outcome, len(cand) - len(kept), shift, score, olen, diff,
            n_consensus, implausible, npolicy_bases, n_rescued, det_bases,
            det_mismatches, rt_strong, det_len)


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
#
# **The read-through check** runs on the pairs numbered below `rt_check_pairs` -- the
# input's first ones, counted from `base_index`, not from the chunk -- so it covers the
# same pairs at any chunk size or thread count (params.READTHROUGH_CHECK_PAIRS).
#
# **Table capacity.** The policy tables cover reads up to a capacity (params.py). A pair
# with a longer read is not consumed: the chunk stops in front of it and returns
# `need`, the read length the tables must cover, and the caller grows them and calls
# again from where this one stopped. `need == 0` means the chunk ran to its end. The
# table prefix never changes when it grows, so where a chunk happened to stop leaves no
# trace in the output.
# =========================================================================== #

#: Counter fields, in the order every chunk function returns them. Mirrors
#: `ChunkCounters` in fastq_chunk.hpp and `_N_COUNTERS` in merge/cli.py.
N_COUNTERS = 16


def _bump(hist, i):
    """Count one observation of value *i*, growing *hist* to fit.

    The histograms are uncapped (the compiled backend sizes its dense arrays from the
    scratch arena; see ``fastq_chunk.hpp``), so both backends return a list whose last
    element is non-zero. Growing to exactly ``i + 1`` and then incrementing ``i`` keeps
    that true here by construction.
    """
    if i >= len(hist):
        hist.extend([0] * (i + 1 - len(hist)))
    hist[i] += 1


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


def _check_sync(h1, h2, index):
    if base_name(h1) != base_name(h2):
        raise InputError(
            f"R1/R2 out of sync at pair {index + 1}: "
            f"'{base_name(h1).decode('latin-1')}' != "
            f"'{base_name(h2).decode('latin-1')}'")


class _Tally:
    """The per-chunk counters and histograms both chunk adapters accumulate."""

    __slots__ = ("n_pairs", "merged", "kept", "emitted", "dropped", "frags_short",
                 "bases_consensus", "implausible", "sum_olen", "sum_diff",
                 "max_read_len", "npolicy_bases", "n_rescued", "det_bases",
                 "det_mismatches", "rt_strong",
                 "len_hist", "olen_hist", "insert_hist", "det_olen_hist")

    def __init__(self):
        for name in self.__slots__[:N_COUNTERS]:
            setattr(self, name, 0)
        self.len_hist, self.olen_hist, self.insert_hist = [], [], []
        self.det_olen_hist = []

    def pair(self, outcome, records, n_dropped, olen, diff, n_consensus, implausible,
             npol_bases, rescued, det_bases, det_mismatches, rt_strong, det_len):
        self.n_pairs += 1
        self.dropped += n_dropped
        self.bases_consensus += n_consensus
        self.implausible += implausible
        self.npolicy_bases += npol_bases
        self.n_rescued += rescued
        self.det_bases += det_bases
        self.det_mismatches += det_mismatches
        self.rt_strong += rt_strong
        # Every overlap that reached T, merged and refused alike: the lengths the gate
        # looked dfit up at, for the expected share of true overlaps it refused.
        if det_len:
            _bump(self.det_olen_hist, det_len)
        if outcome == MERGED:
            self.merged += 1
        else:
            self.kept += 1
            if not records:
                self.frags_short += 1
        # The overlap statistics are over ADMITTED overlaps -- the alignments pairs were
        # merged from, after the plausibility gate -- so `sum_diff / sum_olen` is the
        # post-admission disagreement rate. The pre-gate one is det_* above.
        if olen:
            _bump(self.olen_hist, olen)
            self.sum_olen += olen
            self.sum_diff += diff

    def record(self, outcome, length):
        self.emitted += 1
        _bump(self.len_hist, length)
        if outcome == MERGED:
            _bump(self.insert_hist, length)

    def counters(self):
        return (self.n_pairs, self.merged, self.kept, self.emitted, self.dropped,
                self.frags_short, self.bases_consensus, self.implausible,
                self.sum_olen, self.sum_diff, self.max_read_len, self.npolicy_bases,
                self.n_rescued, self.det_bases, self.det_mismatches, self.rt_strong)


def merge_chunk(buf1, start1, end1, buf2, start2, end2, match_q, step_q, t_table,
                dfit_table, adapter_trimmed, min_read_length, disagree_q, check_sync,
                base_index, npolicy=NPOLICY_TRIM3, rng_seed=0, rt_check_pairs=0):
    """Merge every whole pair available in both buffers (see "Table capacity" and "The
    read-through check" above).

    Returns ``(blob, consumed1, consumed2, counters, len_hist, olen_hist, insert_hist,
    det_olen_hist, need)``: ``olen_hist`` bins the ADMITTED overlaps (the merged
    pairs'), ``det_olen_hist`` every DETECTED one, before the gate.
    """
    _check_weights(match_q, step_q, "merge_chunk()")
    _check_npolicy(npolicy, "merge_chunk()")
    _check_tables(t_table, dfit_table, "merge_chunk()")
    parts = []
    tally = _Tally()
    cap = table_capacity(t_table, dfit_table)
    need = 0
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
        longest = len(s1) if len(s1) > len(s2) else len(s2)
        if longest > cap:
            need = longest
            break
        if longest > tally.max_read_len:
            tally.max_read_len = longest
        if check_sync:
            _check_sync(h1, h2, base_index + tally.n_pairs)

        index = base_index + tally.n_pairs
        (records, outcome, n_dropped, _shift, _score, olen, diff, n_consensus,
         implausible, npol_bases, rescued, det_b, det_d, rt,
         det_len) = _process_pair_ex(
            h1, s1, q1, h2, s2, q2, match_q, step_q, t_table, dfit_table,
            adapter_trimmed, min_read_length, disagree_q, npolicy, rng_seed, index,
            1 if index < rt_check_pairs else 0)
        tally.pair(outcome, records, n_dropped, olen, diff, n_consensus, implausible,
                   npol_bases, rescued, det_b, det_d, rt, det_len)
        for header, seq, qual, _prov in records:
            parts.append(b"@%b\n%b\n+\n%b\n" % (header, seq, qual))
            tally.record(outcome, len(seq))
        pos1, pos2 = try1, try2

    return (b"".join(parts), pos1 - start1, pos2 - start2, tally.counters(),
            tally.len_hist, tally.olen_hist, tally.insert_hist, tally.det_olen_hist,
            need)


def merge_chunk_records(buf1, start1, end1, buf2, start2, end2, match_q, step_q,
                        t_table, dfit_table, adapter_trimmed, min_read_length,
                        disagree_q, check_sync, base_index, want_headers,
                        npolicy=NPOLICY_TRIM3, rng_seed=0, rt_check_pairs=0):
    """Merge every whole pair available, emitting RECORDS instead of FASTQ text.

    The reference half of the ``zna encode --merge-pairs`` adapter; the
    specification the C++ ``merge_chunk_records`` must match element for
    element.  Returns ``(seqs, ends, consumed1, consumed2, counters, len_hist,
    olen_hist, insert_hist, det_olen_hist, need)`` where *seqs* is one bytes blob and each end
    is ``(seq_off, seq_len, hdr_off, hdr_len, slot, prov)``.

    Conventions mirror the text adapter exactly, and the two differ on purpose:
    *consumed* counts are RELATIVE to ``start`` (the caller does ``pos += c``),
    while ``hdr_off`` is ABSOLUTE into the caller's buffer, like
    :func:`split_records`'s return -- buf1 for MERGED and MATE1 records, buf2
    for MATE2, selected by the slot.  ``hdr_len`` is 0 when *want_headers* is
    false.  ``prov`` is the record's PROV_* byte taken directly from the pair
    result -- no ``ZN:i:`` tag round-trip.
    """
    _check_weights(match_q, step_q, "merge_chunk_records()")
    _check_npolicy(npolicy, "merge_chunk_records()")
    _check_tables(t_table, dfit_table, "merge_chunk_records()")
    seq_parts, ends = [], []
    seq_off = 0
    tally = _Tally()
    cap = table_capacity(t_table, dfit_table)
    need = 0
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
        longest = len(s1) if len(s1) > len(s2) else len(s2)
        if longest > cap:
            need = longest
            break
        if longest > tally.max_read_len:
            tally.max_read_len = longest
        if check_sync:
            _check_sync(h1, h2, base_index + tally.n_pairs)

        index = base_index + tally.n_pairs
        (records, outcome, n_dropped, _shift, _score, olen, diff, n_consensus,
         implausible, npol_bases, rescued, det_b, det_d, rt,
         det_len) = _process_pair_ex(
            h1, s1, q1, h2, s2, q2, match_q, step_q, t_table, dfit_table,
            adapter_trimmed, min_read_length, disagree_q, npolicy, rng_seed, index,
            1 if index < rt_check_pairs else 0)
        tally.pair(outcome, records, n_dropped, olen, diff, n_consensus, implausible,
                   npol_bases, rescued, det_b, det_d, rt, det_len)

        for i, (_header, seq, _qual, prov) in enumerate(records):
            slot = (SLOT_MERGED if outcome == MERGED
                    else (SLOT_MATE1 if i == 0 else SLOT_MATE2))
            hdr_off = hdr_len = 0
            if want_headers:
                # The record's OWN source header (MERGED reads R1's): the
                # located record starts at pos with '@', so the header is the
                # next byte.  Length excludes any trailing CR, matching
                # _next_record's rstrip.
                if slot == SLOT_MATE2:
                    hdr_off, hdr_len = pos2 + 1, len(h2)
                else:
                    hdr_off, hdr_len = pos1 + 1, len(h1)
            ends.append((seq_off, len(seq), hdr_off, hdr_len, slot, prov))
            seq_parts.append(seq)
            seq_off += len(seq)
            tally.record(outcome, len(seq))
        pos1, pos2 = try1, try2

    return (b"".join(seq_parts), ends, pos1 - start1, pos2 - start2, tally.counters(),
            tally.len_hist, tally.olen_hist, tally.insert_hist, tally.det_olen_hist,
            need)


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
