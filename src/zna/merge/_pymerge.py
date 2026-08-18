"""The reference merge backend: the scan, in readable Python.

This is the **oracle** the accelerated backend is defined to agree with, not a fallback
for when the fast one is missing. It is never deleted and never optimised at the cost of
clarity — see :mod:`zna.merge.backend`.

Plain Python: no numba, no JIT, nothing between the reader and the algorithm. It is
~50x slower than the compiled backend, and that is the correct trade for an oracle --
speed here would only buy the ability to be wrong in the same way as the thing it
checks.

Scores are integers throughout — see :mod:`zna.merge.params` for why, and
``docs/METHODS.md`` for the argmax total order the visiting order realises.
"""
from __future__ import annotations

from .fastqio import InputError
from .names import base_name, strip_pair_suffix

#: Which mate of the pair an emitted record came from -- the whole geometry
#: transfer of ``zna encode --merge-pairs``.  Mirrors ``Slot`` in
#: ``fastq_chunk.hpp``.
SLOT_MERGED = 0
SLOT_MATE1 = 1
SLOT_MATE2 = 2

# Sentinel for "this shift cannot beat the incumbent"; far below any reachable score.
_REJECT = -(1 << 62)

# Complement table. A/C/G/T/N in both cases; everything else passes through
# UNCOMPLEMENTED -- deliberate, so an IUPAC ambiguity code survives as itself and the
# kernel's N-vs-N semantics are unchanged: rc(b"RYKMSWBDHVN") == b"NVHDBWSMKYR".
_COMPLEMENT = bytes.maketrans(b"ACGTNacgtn", b"TGCANtgcan")


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

#: What to do with a no-call the overlap could not rescue. Same vocabulary as
#: ``zna encode --npolicy``, deliberately: one flag, one meaning, both tools.
NPOLICY_KEEP, NPOLICY_TRIM3, NPOLICY_RANDOM = 0, 1, 2

#: Per-record provenance bits, emitted as the ``ZN:i:<bits>`` header tag. Mirrors the
#: ``PROV_*`` constants in ``merge_core.hpp``; see there for why the byte exists and why
#: there is deliberately no "merged" bit.
PROV_TRIMMED, PROV_RESCUED, PROV_NTRIMMED, PROV_NSUBBED = 1, 2, 4, 8

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


def _consensus_pair_overlap(s1, q1, s2, q2, s2rc, q2r, s, olen, disagree_q,
                            write_r2):
    """Resolve overlap disagreements by posterior, into every emitted copy.

    Returns ``(s1, q1, s2, q2, s2rc, n, rescued1, rescued2)``; ``n`` counts bases changed
    in R1, and the two rescue counts are charged to the mate that was *repaired*, so each
    emitted record's ``rescued_<n>`` token counts only its own recovered no-calls. Their
    sum is the run-level counter, unchanged.

    The decision is symmetric — the better-supported base by posterior from the two
    Phred scores, with the winner's quality derated because a contested base is less
    certain.  The only question is which emitted copies of the overlap receive it, and
    the rule is *every copy that reaches the corpus*.  :func:`process_pair` decides and
    passes ``write_r2``; its comment gives the four cases.

    **N rescue.**  An ``N`` carries no base information, so a real call on the other
    mate beats it whatever the two qualities say, and the rescued base keeps the
    surviving mate's own quality rather than a contested-base derating — there was no
    contest.  Without this the rescue happened only by luck, because an instrument
    usually assigns an N a low quality; a *high*-quality N beat a real base and survived
    into the corpus.  Only ``N`` is rescued, not the IUPAC codes, which do carry partial
    information.

    Rescue does not touch the *scan*: an N still counts as a mismatch when the shift is
    scored, so which shift wins — and therefore whether the pair merges — is unchanged.

    ``s2``/``q2`` are R2 in its own orientation, which is how it is emitted: overlap slot
    ``b`` on the reverse-complemented axis is R2 index ``len2 - 1 - b``, holding the
    complement of the resolved call.  New ``bytes`` if anything changed, else the
    originals unchanged.
    """
    a0 = s if s > 0 else 0        # mirrors the scan's overlap alignment
    b0 = -s if s < 0 else 0
    len2 = len(s2)
    s1b = q1b = s2b = q2b = s2rcb = None
    n = rescued1 = rescued2 = 0
    for i in range(olen):
        a = a0 + i
        b = b0 + i
        qa = q1[a]
        qb = q2r[b]
        if s1[a] != s2rc[b]:
            if s1b is None:
                s1b, q1b = bytearray(s1), bytearray(q1)
                if write_r2:
                    s2b, q2b = bytearray(s2), bytearray(q2)
                    s2rcb = bytearray(s2rc)
            j2 = len2 - 1 - b                 # the same base, in R2's own frame
            a_is_n = s1[a] == 0x4E            # ord('N'); sequences are upper-cased
            b_is_n = s2rc[b] == 0x4E
            if a_is_n != b_is_n:              # rescue: a real call beats an N
                if a_is_n:                    # R2 rescues R1
                    s1b[a] = s2rc[b]
                    q1b[a] = qb
                    n += 1                    # `n` counts R1 only, by contract
                    rescued1 += 1
                elif write_r2:                # R1 rescues R2
                    s2rcb[b] = s1b[a]
                    s2b[j2] = _COMPLEMENT[s1b[a]]
                    q2b[j2] = qa
                    rescued2 += 1             # only where it is actually WRITTEN
            elif a_is_n:                      # both are N: nothing to rescue from
                pass
            elif qb > qa:                     # R2 is the better-supported call
                nq = disagree_q[qb * 256 + qa]
                s1b[a] = s2rc[b]              # R2's base already stands in s2rc/s2
                q1b[a] = nq
                if write_r2:
                    q2b[j2] = nq
                n += 1
            else:                             # R1 wins, but it is contested: derate it
                nq = disagree_q[qa * 256 + qb]
                q1b[a] = nq
                if write_r2:
                    s2rcb[b] = s1b[a]         # keep s2rc the exact revcomp of s2
                    s2b[j2] = _COMPLEMENT[s1b[a]]
                    q2b[j2] = nq
    if s1b is None:
        return s1, q1, s2, q2, s2rc, 0, 0, 0
    if not write_r2:
        return bytes(s1b), bytes(q1b), s2, q2, s2rc, n, rescued1, rescued2
    return (bytes(s1b), bytes(q1b), bytes(s2b), bytes(q2b), bytes(s2rcb), n,
            rescued1, rescued2)


def _balanced_split(L, len1, len2):
    """Emitted lengths for a trimmed pair: as close to equal as the geometry allows.

    The pair must tile the fragment exactly once, so ``keep1 + keep2 == L`` is forced and
    the only freedom is where the cut falls. Splitting the overlap down the middle rather
    than taking all of it off R2 keeps the two emitted reads the same length -- what
    downstream aligners and models expect -- and discards the *last* cycles of both
    reads, the lowest-quality bases in the pair, instead of one read's entire copy.

    ``keep1`` is clamped into ``[L - len2, len1]``, the range in which both reads can
    supply their share; for equal-length mates the clamp never binds and this is exactly
    "cut ``olen / 2`` from each".
    """
    k = (L + 1) // 2                  # an odd overlap leaves the extra base on R1
    k = max(k, L - len2)
    k = min(k, len1)
    return k, L - k


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


def _trim_is_allowed(L, len1, len2, lr):
    """May this pair be trimmed? Each mate must reach at least ``lr`` past the other's
    3' end — which both keeps every emitted read above the length filter and caps the
    overlap the trim band may act on."""
    return (L - len1) >= lr and (L - len2) >= lr


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


def process_pair(h1, s1, q1, h2, s2, q2, match_q, step_q, t_merge_q, t_trim_q,
                 min_read_length, disagree_q, npolicy=NPOLICY_TRIM3, rng_seed=0,
                 pair_index=0):
    """Classify one pair and build its output records.

    Returns ``(records, outcome, n_dropped, score_q, overlap_len, mismatches,
    bases_consensus_changed, trim_guard_fired, npolicy_bases, n_rescued)``,
    with each record a ``(header, seq, qual)`` tuple.  The thin public shim over
    :func:`_process_pair_ex`, which additionally carries each record's PROV_*
    byte -- the record adapter reads the bits there directly, with no
    ``ZN:i:`` tag round-trip, mirroring ``PairResult::prov`` in the C++ core.
    """
    records, *rest = _process_pair_ex(
        h1, s1, q1, h2, s2, q2, match_q, step_q, t_merge_q, t_trim_q,
        min_read_length, disagree_q, npolicy, rng_seed, pair_index)
    return ([r[:3] for r in records], *rest)


def _process_pair_ex(h1, s1, q1, h2, s2, q2, match_q, step_q, t_merge_q, t_trim_q,
                     min_read_length, disagree_q, npolicy=NPOLICY_TRIM3, rng_seed=0,
                     pair_index=0):
    """:func:`process_pair` with records as ``(header, seq, qual, prov)``.

    The last two are pair totals. Their per-mate splits stay local: they exist only to
    build each record's provenance tokens, and summing them here keeps the run-level
    counter tuple at its existing width — see :func:`_prov_name` and trap 6 in
    ``docs/HANDOFF_0.4.0.md``.

    **The decision** is a single ``argmax`` shift read at two thresholds::

        score >= t_merge_q            -> merge (one full-fragment record)
        t_trim_q <= score < t_merge_q -> keep both, split the redundant overlap between
                                         their 3' ends so the fragment is tiled once
        score <  t_trim_q             -> keep both unchanged

    **The trim is symmetric.** The overlap sits at the 3' end of *both* mates (each read
    starts at a fragment end and reads inward), so cutting half from each tiles the
    fragment exactly once just as taking it all off R2 did -- and leaves the two emitted
    reads the same length, which is what downstream aligners and models expect, while
    discarding the last cycles of both reads rather than one read's whole copy. Because
    both mates now keep part of the overlap, the consensus is written into both.

    Trimming applies only to the normal (``s >= 0``) geometry: in a read-through the
    redundant bases are R2's *5'* fragment copy and its 3' end is adapter, so there is
    nothing sensible to cut -- such a pair is either merged or kept whole.

    **Trim guard:** a trim that would leave *either* read below ``min_read_length`` keeps
    both *untrimmed* instead, turning a would-be whole-fragment discard into a no-op. It
    binds far less often than it did, because balancing puts both reads near ``L / 2``.

    **Pair integrity:** an unmerged pair is emitted all-or-nothing. A lone surviving mate
    would be encoded as a spurious "single" -- a full molecule with both endpoints --
    corrupting the fragment-end supervision. A merged read is a genuine full molecule and
    is filtered on its own.
    """
    len1, len2 = len(s1), len(s2)
    s2rc = reverse_complement(s2)
    shift, score, olen, diff = scan(s1, s2rc, len1, len2, match_q, step_q, t_trim_q)

    lr = min_read_length
    L = shift + len2                       # the inferred fragment length; 0 if no overlap

    # PROVISIONAL decision, on the full reads. It settles where the consensus is written
    # (below) and carries the evidence forward: the score is a statement about this
    # pair's fragment length, and trimming interior bases cannot change a fragment's
    # length. The final decision is re-taken after trimming, in `_decide` -- but on
    # GEOMETRY, never by re-scoring. See docs/NPOLICY_PLAN.md D4a.
    prov_merge = olen > 0 and score >= t_merge_q
    prov_band = olen > 0 and shift >= 0 and t_trim_q <= score < t_merge_q
    prov_trim = prov_band and _trim_is_allowed(L, len1, len2, lr)

    # WHERE the consensus is written: into the records whose CONSTRUCTION depends on
    # the overlap being real, and nowhere else.
    #
    #   merged  -> R1 alone. R1's overlap region becomes the merged record; R2
    #              contributes only outside it, so its copy is discarded.
    #   trimmed -> both. Each mate keeps part of the overlap, so both copies are emitted
    #              and both must carry the same call.
    #   kept    -> neither. Nothing emitted depends on the alignment being right.
    #
    # The kept case is what the measurements say, not a symmetry nicety. A detection
    # that lands in KEPT is spurious almost by construction: at `shift >= 0` it is here
    # only because `trim_is_allowed` refused it, which needs an inferred overlap over
    # ~110 bases -- and a genuine one of that length scores ~218 bits and would have
    # merged; at `shift < 0` a genuine read-through overlap equals the fragment length,
    # so scoring in [8, 28) needs a 5-14 bp fragment. An overlap too suspect to CUT on is
    # too suspect to REWRITE BASES on.
    #
    # Measured on 1M ground-truth pairs: of 3,068 kept pairs with a detected overlap,
    # ZERO found the true shift and 97.3% had no true overlap at all. Wrong emitted bases
    # in the overlap window -- correct-neither 208, correct-R1-only (the old behaviour)
    # 1,509, correct-both 17,870. The old R1 write turned 1,379 correct bases wrong to
    # fix 78.
    write_r1 = prov_merge or prov_trim
    write_r2 = prov_trim

    n_consensus = trim_guard = 0
    npolicy_1 = npolicy_2 = rescued_1 = rescued_2 = 0
    if diff > 0 and write_r1:
        s1, q1, s2, q2, s2rc, n_consensus, rescued_1, rescued_2 = \
            _consensus_pair_overlap(
                s1, q1, s2, q2, s2rc, q2[::-1], shift, olen, disagree_q, write_r2)

    # ---- trim3: cut each read at its first SURVIVING N -------------------------
    #
    # After the rescue, so a no-call the mate could answer costs nothing. 3' only, so
    # both 5' anchors -- the two fragment termini -- are untouched however short the
    # reads get.
    if npolicy == NPOLICY_RANDOM:
        # Substitution does not change a length, so the coverage test below is
        # unaffected and `random` never costs a merge -- unlike trim3.
        s1n, k1 = _sub_n(s1, rng_seed, pair_index * 2)
        s2n, k2 = _sub_n(s2, rng_seed, pair_index * 2 + 1)
        npolicy_1, npolicy_2 = k1, k2
        if k1 or k2:
            s1, s2 = s1n, s2n
            s2rc = reverse_complement(s2)
    elif npolicy == NPOLICY_TRIM3:
        s1t, q1t, k1 = _trim3(s1, q1)
        s2t, q2t, k2 = _trim3(s2, q2)
        if k1 != len1 or k2 != len2:
            npolicy_1, npolicy_2 = len1 - k1, len2 - k2
            s1, q1, s2, q2 = s1t, q1t, s2t, q2t
            len1, len2 = k1, k2
            s2rc = reverse_complement(s2)
            # `shift` is the offset of revcomp(R2) on the shared axis, so it is tied to
            # len2. R2 keeps its 5' anchor at fragment position L-1, so the trimmed mate
            # covers [L - len2, L) and the offset becomes L - len2. L itself is unchanged
            # -- that is the whole point.
            shift = L - len2

    # ---- the final decision, on GEOMETRY, reusing the original evidence --------
    #
    # The pair still tiles the fragment iff len1 + len2 >= L. When it does, the
    # reconstruction IS the fragment, exactly and N-free. Nothing is re-scored: trimming
    # cuts 3' ends, which is where a normal overlap lives, so a re-scan would refuse
    # merges it had ample evidence for a moment earlier.
    covers = olen > 0 and (len1 + len2) >= L
    will_merge = prov_merge and covers
    will_trim = (not will_merge) and prov_band and (len1 + len2) > L \
        and _trim_is_allowed(L, len1, len2, lr)
    if will_trim:
        keep1, keep2 = _balanced_split(L, len1, len2)

    # The run-level counters are the per-mate ones summed, so they keep their exact
    # previous values and the 15-field counter tuple does not grow -- `_fold` in
    # merge/cli.py sums a fixed prefix, and a counter added past it reports zero.
    npolicy_bases = npolicy_1 + npolicy_2
    n_rescued = rescued_1 + rescued_2
    # Which policy bit a touched record earns, and which token carries its count.
    rnd = npolicy == NPOLICY_RANDOM
    npolicy_bit = PROV_NSUBBED if rnd else PROV_NTRIMMED

    if will_merge:
        seq, qual, n1, n2 = _build_merged(shift, s1, q1, s2rc, q2, len1, len2)
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
    elif will_trim:
        b1 = (PROV_TRIMMED | (PROV_RESCUED if rescued_1 else 0)
              | (npolicy_bit if npolicy_1 else 0))
        b2 = (PROV_TRIMMED | (PROV_RESCUED if rescued_2 else 0)
              | (npolicy_bit if npolicy_2 else 0))
        cand = [(_prov_name(h1, b1, 0 if rnd else npolicy_1,
                            npolicy_1 if rnd else 0, rescued_1),
                 s1[:keep1], q1[:keep1], b1),
                (_prov_name(h2, b2, 0 if rnd else npolicy_2,
                            npolicy_2 if rnd else 0, rescued_2),
                 s2[:keep2], q2[:keep2], b2)]
        paired, outcome = True, TRIMMED
    else:
        if prov_band and not prov_trim:
            trim_guard = 1                                  # guard fired
        # No consensus is written on this path, so a kept record can never carry
        # PROV_RESCUED -- only the N policy can have touched it.
        b1 = npolicy_bit if npolicy_1 else 0
        b2 = npolicy_bit if npolicy_2 else 0
        cand = [(_prov_name(h1, b1, 0 if rnd else npolicy_1,
                            npolicy_1 if rnd else 0, 0), s1, q1, b1),
                (_prov_name(h2, b2, 0 if rnd else npolicy_2,
                            npolicy_2 if rnd else 0, 0), s2, q2, b2)]
        paired, outcome = True, KEPT

    if paired:
        kept = cand if (len(cand[0][1]) >= lr and len(cand[1][1]) >= lr) else []
    else:
        kept = [r for r in cand if len(r[1]) >= lr]
    return (kept, outcome, len(cand) - len(kept), score, olen, diff,
            n_consensus, trim_guard, npolicy_bases, n_rescued)


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


def merge_chunk(buf1, start1, end1, buf2, start2, end2, match_q, step_q, t_merge_q,
                t_trim_q, min_read_length, disagree_q, check_sync, base_index,
                npolicy=NPOLICY_TRIM3, rng_seed=0):
    """Merge every whole pair available in both buffers.

    Returns ``(blob, consumed1, consumed2, counters, len_hist, olen_hist,
    insert_hist)``.
    """
    parts = []
    n_pairs = merged = trimmed = kept = emitted = dropped = 0
    bases_trimmed = frags_short = bases_consensus = trim_guard = 0
    sum_olen = sum_diff = max_read_len = 0
    npolicy_bases = n_rescued_tot = 0
    len_hist, olen_hist, insert_hist = [], [], []
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
        if longest > max_read_len:
            max_read_len = longest

        if check_sync and base_name(h1) != base_name(h2):
            raise InputError(
                f"R1/R2 out of sync at pair {base_index + n_pairs + 1}: "
                f"'{base_name(h1).decode('latin-1')}' != "
                f"'{base_name(h2).decode('latin-1')}'")

        (records, outcome, n_dropped, score, olen, diff,
         n_consensus, guard, npol_bases, rescued) = process_pair(
            h1, s1, q1, h2, s2, q2, match_q, step_q, t_merge_q, t_trim_q,
            min_read_length, disagree_q, npolicy, rng_seed, base_index + n_pairs)

        n_pairs += 1
        dropped += n_dropped
        bases_consensus += n_consensus
        trim_guard += guard
        npolicy_bases += npol_bases
        n_rescued_tot += rescued
        if outcome == MERGED:
            merged += 1
        elif outcome == TRIMMED:
            trimmed += 1
            # Both mates are cut now, so charge both: the counter is "redundant bases
            # removed", which is the overlap length either way.
            if records:
                bases_trimmed += ((len(s1) - len(records[0][1]))
                                  + (len(s2) - len(records[1][1])))
        else:
            kept += 1
        if not records and outcome != MERGED:
            frags_short += 1
        if olen:
            _bump(olen_hist, olen)
            sum_olen += olen
            sum_diff += diff
        for header, seq, qual in records:
            parts.append(b"@%b\n%b\n+\n%b\n" % (header, seq, qual))
            emitted += 1
            L = len(seq)
            _bump(len_hist, L)
            if outcome == MERGED:
                _bump(insert_hist, L)
        pos1, pos2 = try1, try2

    counters = (n_pairs, merged, trimmed, kept, emitted, dropped, bases_trimmed,
                frags_short, bases_consensus, trim_guard, sum_olen, sum_diff,
                max_read_len, npolicy_bases, n_rescued_tot)
    return (b"".join(parts), pos1 - start1, pos2 - start2, counters,
            len_hist, olen_hist, insert_hist)


def merge_chunk_records(buf1, start1, end1, buf2, start2, end2, match_q, step_q,
                        t_merge_q, t_trim_q, min_read_length, disagree_q,
                        check_sync, base_index, want_headers,
                        npolicy=NPOLICY_TRIM3, rng_seed=0):
    """Merge every whole pair available, emitting RECORDS instead of FASTQ text.

    The reference half of the ``zna encode --merge-pairs`` adapter; the
    specification the C++ ``merge_chunk_records`` must match element for
    element.  Returns ``(seqs, ends, consumed1, consumed2, counters, len_hist,
    olen_hist, insert_hist)`` where *seqs* is one bytes blob and each end is
    ``(seq_off, seq_len, hdr_off, hdr_len, slot, prov)``.

    Conventions mirror the text adapter exactly, and the two differ on purpose:
    *consumed* counts are RELATIVE to ``start`` (the caller does ``pos += c``),
    while ``hdr_off`` is ABSOLUTE into the caller's buffer, like
    :func:`split_records`'s return -- buf1 for MERGED and MATE1 records, buf2
    for MATE2, selected by the slot.  ``hdr_len`` is 0 when *want_headers* is
    false.  ``prov`` is the record's PROV_* byte taken directly from the pair
    result -- no ``ZN:i:`` tag round-trip.
    """
    seq_parts, ends = [], []
    seq_off = 0
    n_pairs = merged = trimmed = kept = emitted = dropped = 0
    bases_trimmed = frags_short = bases_consensus = trim_guard = 0
    sum_olen = sum_diff = max_read_len = 0
    npolicy_bases = n_rescued_tot = 0
    len_hist, olen_hist, insert_hist = [], [], []
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
        if longest > max_read_len:
            max_read_len = longest

        if check_sync and base_name(h1) != base_name(h2):
            raise InputError(
                f"R1/R2 out of sync at pair {base_index + n_pairs + 1}: "
                f"'{base_name(h1).decode('latin-1')}' != "
                f"'{base_name(h2).decode('latin-1')}'")

        (records, outcome, n_dropped, score, olen, diff,
         n_consensus, guard, npol_bases, rescued) = _process_pair_ex(
            h1, s1, q1, h2, s2, q2, match_q, step_q, t_merge_q, t_trim_q,
            min_read_length, disagree_q, npolicy, rng_seed, base_index + n_pairs)

        n_pairs += 1
        dropped += n_dropped
        bases_consensus += n_consensus
        trim_guard += guard
        npolicy_bases += npol_bases
        n_rescued_tot += rescued
        if outcome == MERGED:
            merged += 1
        elif outcome == TRIMMED:
            trimmed += 1
            if records:
                bases_trimmed += ((len(s1) - len(records[0][1]))
                                  + (len(s2) - len(records[1][1])))
        else:
            kept += 1
        if not records and outcome != MERGED:
            frags_short += 1
        if olen:
            _bump(olen_hist, olen)
            sum_olen += olen
            sum_diff += diff

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
            emitted += 1
            L = len(seq)
            _bump(len_hist, L)
            if outcome == MERGED:
                _bump(insert_hist, L)
        pos1, pos2 = try1, try2

    counters = (n_pairs, merged, trimmed, kept, emitted, dropped, bases_trimmed,
                frags_short, bases_consensus, trim_guard, sum_olen, sum_diff,
                max_read_len, npolicy_bases, n_rescued_tot)
    return (b"".join(seq_parts), ends, pos1 - start1, pos2 - start2, counters,
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
