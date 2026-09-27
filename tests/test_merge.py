"""Tests for src/zna/merge: single-axis LR overlap scoring and the 0.6 merge policy.

Runs with or without the compiled backend. The reference kernel in ``_pymerge`` is the
oracle the compiled one is defined to agree with, so most of what is checked here is
checked against both; the cross-backend classes skip when no compatible extension is
built (a 0.5.x build is refused by its missing ``POLICY_ABI``).

The policy is ``docs/archive/MERGE_ACCURACY_PLAN.md`` §2: merge the best eligible shift when it
reaches the pair's floor ``T = log2((len1 + len2 - 1) / alpha)`` and its informative
mismatches are plausible at the same ``alpha`` (``d <= dfit[n]``); keep the pair whole
otherwise. The suite covers, in order:

  1. exact derivations (T, dfit, weights)     7. find_overlap and the contract range
  2. the run's diagnostics                    8. the plausibility gate
  3. backend selection and the ABI marker     9. detection, read-through, boundaries
  4. cross-backend equivalence               10. process_pair
  5. the fixed-point scale                   11. the CLI, its warnings, table growth
  6. the argmax total order                  12. the merge review's 19 cases

0.5.3's trim band is gone, and with it every test of the balanced split, the trim guard,
the trim-path consensus, PROV_TRIMMED and the ``--threshold-*`` flags.
"""
import gzip
import json
import math
import random
import sys
from fractions import Fraction
from pathlib import Path

import pytest

from zna.merge import cli
from zna.merge import params as zparams
from zna.merge.cli import (
    BASES_CONSENSUS, DET_BASES, DET_MISMATCHES, IMPLAUSIBLE, MAX_READ_LEN, MERGED,
    N_PAIRS, N_RESCUED, NPOLICY_BASES, RT_STRONG,
)
from zna.merge.overlap import (
    IMPLAUSIBLE as V_IMPLAUSIBLE, MERGE as V_MERGE, NONE as V_NONE,
    find_overlap, reverse_complement, scan_unrestricted,
)
from zna.merge.params import (
    DISAGREE_Q, SCALE, MergeParams, binom_cap, decimal_str, exact, score_weights,
    threshold_bits, threshold_q, to_q, weights_q,
)
from zna.merge.pairs import PairOutcome, base_name, process_pair


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #

ADAPTER1 = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER2 = b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

# Score weights at e = 0.01 -- the default, and the rate every fixture below pins -- used
# to build exact expectations.
MATCH_W, MISMATCH_W = score_weights(0.01)


def _fast_backend():
    """The reference backend is ~50x slower than the compiled one, so the statistical
    sweeps below size themselves to whichever is running."""
    from zna.merge.backend import available_merge_backends
    return "accel" in available_merge_backends()


# The kernel scores in fixed point (zna/merge/params.py), so every expectation here is
# an exact integer -- no pytest.approx, and no float anywhere near a decision.
_P = MergeParams(error_rate="0.01")

# e = 0.01, min_read_length=1 so tiny test reads survive the QC filter.
P = MergeParams(error_rate="0.01", min_read_length=1)


def t_q(len1, len2, p=_P):
    """The pair's merge floor in fixed point."""
    return p.t_q(len1 + len2 - 1)


def rc(seq: bytes) -> bytes:
    return reverse_complement(seq)


def qual(seq: bytes) -> bytes:
    return b"I" * len(seq)


def make_pair(fragment: bytes, read_len: int, name=b"frag"):
    """Simulate a PE read of `fragment`: R1 = 5' end, R2 = revcomp of 3' end."""
    r1 = fragment[:read_len]
    r2 = rc(fragment[-read_len:])
    return (name + b"/1", r1, qual(r1)), (name + b"/2", r2, qual(r2))


def rand_seq(n: int, seed: int) -> bytes:
    return bytes("".join(random.Random(seed).choices("ACGT", k=n)), "ascii")


def draw(rng, n) -> bytes:
    return bytes("".join(rng.choices("ACGT", k=n)), "ascii")


def cycle_pair(fragment: bytes, read_len: int, rng, name=b"frag"):
    """Full-cycle reads: fragment, then adapter, then random filler to the cycle length.

    A shorter-than-cycle read makes the true shift ``L - len(read2)`` rather than
    ``L - readlen``; building reads at full length keeps the geometry unambiguous.
    """
    r1 = (fragment + ADAPTER1 + draw(rng, read_len))[:read_len]
    r2 = (rc(fragment) + ADAPTER2 + draw(rng, read_len))[:read_len]
    return (name + b"/1", r1, qual(r1)), (name + b"/2", r2, qual(r2))


def mutate(seq: bytes, rng, err: float) -> bytes:
    b = bytearray(seq)
    for i in range(len(b)):
        if rng.random() < err:
            b[i] = ord(rng.choice("ACGT"))
    return bytes(b)


def flip(b: int) -> int:
    """A base guaranteed to differ from *b*."""
    return ord("A") if b != ord("A") else ord("C")


def score_of(matches: int, mismatches: int = 0) -> int:
    """The exact fixed-point score of an overlap with these match/mismatch counts."""
    return (matches + mismatches) * _P.match_q - mismatches * _P.step_q


def score_bits(matches: int, mismatches: int = 0) -> float:
    """The same score in bits, for the arithmetic tests that pin the weights."""
    return matches * MATCH_W - mismatches * MISMATCH_W


def min_matches(threshold_q: int, mismatches: int) -> int:
    """Fewest matching bases that reach `threshold_q` with `mismatches` mismatches."""
    m = 0
    while score_of(m, mismatches) < threshold_q:
        m += 1
    return m


def one_n(s1, s2rc, s, n):
    """Positions of the overlap at shift s where exactly one base is N (brute force)."""
    i1, i2 = max(s, 0), max(-s, 0)
    return sum((s1[i1 + k] == 78) != (s2rc[i2 + k] == 78) for k in range(n))


# --------------------------------------------------------------------------- #
# the pre-redesign rule, kept ONLY to pin parity where parity is expected
# --------------------------------------------------------------------------- #

def legacy_scan(s1, s2rc, len1, len2, require, diff_limit, diff_pct):
    """fastp's first-accept scan, as this tool shipped it before the LR redesign."""
    off = 0
    while off < len1 - require:
        olen = len1 - off
        if olen > len2:
            olen = len2
        dl = int(olen * diff_pct)
        if dl > diff_limit:
            dl = diff_limit
        diff = 0
        k = 0
        while k < olen:
            if s1[off + k] != s2rc[k]:
                diff += 1
                if diff > dl:
                    break
            k += 1
        if diff <= dl:
            return 1, off, olen, diff
        off += 1
    sh = 1
    while sh < len2 - require:
        olen = len1 if len1 <= len2 - sh else len2 - sh
        dl = int(olen * diff_pct)
        if dl > diff_limit:
            dl = diff_limit
        diff = 0
        k = 0
        while k < olen:
            if s1[k] != s2rc[sh + k]:
                diff += 1
                if diff > dl:
                    break
            k += 1
        if diff <= dl:
            return -1, sh, olen, diff
        sh += 1
    return 0, 0, 0, 0


# --------------------------------------------------------------------------- #
# reverse_complement
# --------------------------------------------------------------------------- #

class TestReverseComplement:
    def test_basic(self):
        assert rc(b"ACGT") == b"ACGT"
        assert rc(b"AAAA") == b"TTTT"
        assert rc(b"ACGTACGTAAAA") == b"TTTTACGTACGT"

    def test_n_is_preserved(self):
        assert rc(b"ACGTN") == b"NACGT"

    def test_roundtrip(self):
        s = rand_seq(100, 7)
        assert rc(rc(s)) == s


# --------------------------------------------------------------------------- #
# 1. exact derivations -- the policy's integers must not silently drift
# --------------------------------------------------------------------------- #

def _binom_tail(n, e, d):
    """P(Binom(n, e) >= d), exactly, the slow way."""
    return sum(math.comb(n, j) * e ** j * (1 - e) ** (n - j) for j in range(d, n + 1))


class TestDerivations:
    """Every number the kernel compares against is derived once, exactly, in params.py.

    These goldens are the integers a corpus was made with. A failure here is not a
    rounding nuisance: it means the same FASTQ now produces a different corpus.
    """

    def test_weights_at_one_percent_are_0_5_3s_integers(self):
        """0.5.3 derived these with libm's log2 and pinned them per platform; exact
        decimal arithmetic lands on the same two integers."""
        assert weights_q(Fraction(1, 100)) == (33_311_170, 137_813_407)
        assert (_P.match_q, _P.step_q) == (33_311_170, 137_813_407)

    def test_float_weights_for_humans(self):
        assert round(MATCH_W, 4) == 1.9855
        assert round(MISMATCH_W, 4) == 6.2288

    @pytest.mark.parametrize("n_shifts,golden", [
        (99, 445_618_379),     # 2x50: 26.56 bits
        (199, 462_517_532),    # 2x100: 27.57
        (299, 472_372_084),    # 2x150: 28.16 -- 0.5.3's fixed 28 was this, rounded
        (599, 489_189_741),    # 2x300: 29.16
    ])
    def test_the_floor_goldens(self, n_shifts, golden):
        assert threshold_q(n_shifts, "1e-6") == golden
        assert _P.t_q(n_shifts) == golden

    def test_the_floor_is_log2_n_over_alpha(self):
        assert round(threshold_bits(299, "1e-6"), 3) == 28.156
        assert round(threshold_bits(99, "1e-6"), 2) == 26.56
        assert round(threshold_bits(599, "1e-6"), 2) == 29.16
        # each factor of 10 in alpha is log2(10) = 3.32 bits
        assert round(threshold_bits(299, "1e-7") - threshold_bits(299, "1e-6"), 4) \
            == round(math.log2(10), 4)

    def test_the_floor_table_is_monotone_and_starts_empty(self):
        tab = _P.t_table
        assert tab[0] == 0                                  # N = 0 is never looked up
        assert all(tab[i] < tab[i + 1] for i in range(1, len(tab) - 1))

    @pytest.mark.parametrize("n,golden", [(20, 5), (64, 7), (122, 9), (150, 10)])
    def test_dfit_goldens_at_one_percent(self, n, golden):
        """The review's C01 (19 mismatches in 122 bases) is 10 past dfit[122] = 9."""
        assert _P.dfit(n) == golden

    def test_dfit_head(self):
        assert list(_P.dfit_table[:10]) == [0, 1, 2, 3, 3, 3, 3, 3, 3, 4]

    @pytest.mark.parametrize("e,alpha", [("0.01", "1e-6"), ("0.003", "1e-6"),
                                         ("0.05", "1e-3"), ("0.00016", "1e-6"),
                                         ("0.2", "0.01")])
    def test_dfit_is_the_definition(self, e, alpha):
        """max d with P(Binom(n, e) >= d) >= alpha, checked against brute force."""
        e, alpha = Fraction(e), Fraction(alpha)
        tab = binom_cap(60, e, alpha)
        for n in range(61):
            d = tab[n]
            assert _binom_tail(n, e, d) >= alpha, (n, d)
            if d < n:
                assert _binom_tail(n, e, d + 1) < alpha, (n, d)
        assert all(tab[i] <= tab[i + 1] for i in range(60))     # monotone in n

    def test_tables_grow_by_doubling_and_keep_their_prefix(self):
        p = MergeParams(error_rate="0.0123", alpha="1e-5")
        t0, d0 = list(p.t_table), list(p.dfit_table)
        assert p.capacity == 256 and len(t0) == 512 and len(d0) == 257
        p.ensure(300)
        assert p.capacity == 512
        assert list(p.t_table[:512]) == t0 and list(p.dfit_table[:257]) == d0
        p.ensure(1500)
        assert p.capacity == 2048 and len(p.dfit_table) == 2049
        # a fresh table built straight to that size is the same table
        assert list(p.dfit_table) == binom_cap(2048, Fraction("0.0123"),
                                               Fraction("1e-5"))

    def test_dfit_moves_one_step_per_n_and_is_the_definition_at_long_overlaps(self):
        """`binom_cap` tests one candidate per n (a trial adds at most one success, so
        dfit[n+1] is dfit[n] or dfit[n] + 1). Checked against the definition, in exact
        integers, far past the brute-force range above, and resumed from prefixes."""
        e, alpha = Fraction("0.0087"), Fraction("1e-6")
        a, b, A, B = e.numerator, e.denominator, alpha.numerator, alpha.denominator
        tab = binom_cap(1200, e, alpha)
        assert all(tab[i + 1] - tab[i] in (0, 1) for i in range(1200))

        def reaches(n, d):             # P(Binom(n, e) >= d) >= alpha
            tail = sum(math.comb(n, j) * a ** j * (b - a) ** (n - j)
                       for j in range(d, n + 1))
            return tail * B >= b ** n * A
        for n in (257, 700, 1200):
            assert reaches(n, tab[n]) and not reaches(n, tab[n] + 1), n
        for cut in (0, 1, 2, 99, 600):
            assert binom_cap(1200, e, alpha, start=tab[:cut + 1]) == tab

    def test_the_gate_refuses_a_true_overlap_with_probability_at_most_alpha(self):
        """The closed form of the gate's sensitivity (plan §8, layer 1): a true overlap
        of n bases at true rate e' is refused with P(Binom(n, e') > dfit[n]). At e' = e
        that is <= alpha by construction; setting e below the library's true rate
        costs in proportion -- which is what the detected-rate warning is for."""
        e, alpha = Fraction(1, 100), Fraction(1, 10 ** 6)
        for n in (20, 50, 100, 150):
            d = _P.dfit(n)
            assert _binom_tail(n, e, d + 1) <= alpha
        # the plan's numbers: at 3%, 0.07-0.6%; at 5%, 1-13% (n = 50..150)
        at3 = [float(_binom_tail(n, Fraction(3, 100), _P.dfit(n) + 1)) for n in (50, 150)]
        at5 = [float(_binom_tail(n, Fraction(5, 100), _P.dfit(n) + 1)) for n in (50, 150)]
        assert 5e-4 < at3[0] < 1e-3 and 5e-3 < at3[1] < 7e-3, at3
        assert 0.01 < at5[0] < 0.02 and 0.12 < at5[1] < 0.14, at5

    def test_inputs_are_exact_rationals(self):
        assert exact("1e-6") == Fraction(1, 10 ** 6)
        assert exact(0.01) == Fraction(1, 100)              # through repr, not binary
        assert exact("0.01") == Fraction(1, 100)
        assert decimal_str(Fraction(1, 100)) == "0.01"
        assert decimal_str(Fraction(1, 10 ** 6)) == "0.000001"
        assert decimal_str(Fraction(87404, 10 ** 7)) == "0.0087404"

    def test_the_error_rate_defaults_to_one_percent(self):
        """0.5.3's hidden constant, now a documented default: the goldens above are
        the default's."""
        assert MergeParams().e == Fraction(1, 100)
        assert MergeParams() == MergeParams(error_rate="0.01")

    def test_the_error_rate_is_used_exactly_as_given(self):
        """Nothing is rounded: the value typed is the value derived from and recorded."""
        a = MergeParams(error_rate="0.0087404")
        b = MergeParams(error_rate="0.00874")
        assert a.e == Fraction(87404, 10 ** 7) and b.e == Fraction(874, 10 ** 5)
        assert a.match_q != b.match_q or a.step_q != b.step_q
        assert a.merge_record()["error_rate"] == "0.0087404"
        # a float means its shortest repr, not the nearest double
        assert MergeParams(error_rate=0.03).e == Fraction(3, 100)
        # a tiny rate is a legitimate setting, not a rounding casualty
        assert MergeParams(error_rate="4e-7").e == Fraction(4, 10 ** 7)

    @pytest.mark.parametrize("kw,match", [
        (dict(alpha="0"), "alpha"), (dict(alpha="1"), "alpha"),
        (dict(alpha="-1e-6"), "alpha"), (dict(alpha="abc"), "alpha"),
        (dict(error_rate="0"), "must be > 0"), (dict(error_rate="-0.01"), "must be > 0"),
        (dict(error_rate="0.75"), "must be < 0.75"),
        (dict(error_rate="0.9"), "must be < 0.75"), (dict(error_rate="x"), "not a number"),
        (dict(error_rate=None), "not a number"),
        (dict(error_rate=Fraction(1, 3)), "no exact decimal form"),
        (dict(alpha=Fraction(1, 3 * 10 ** 6)), "no exact decimal form"),
        # a decimal string, and only that: Fraction() alone would take all three
        (dict(error_rate="1/100"), "not a number"),
        (dict(error_rate="1_0e-3"), "not a number"),
        (dict(error_rate="\u0661\u0660"), "not a number"),
        # a mistyped exponent is refused, not left building tables
        (dict(alpha="1e-301"), "alpha must be >= 1e-300"),
        (dict(error_rate="9e-10"), "error rate must be >= 1e-9"),
        (dict(alpha="1e-10000000"), "out of range"),
        (dict(error_rate="1e+10000000"), "out of range"),
    ])
    def test_values_the_policy_cannot_mean_are_refused(self, kw, match):
        with pytest.raises(ValueError, match=match):
            MergeParams(**kw)

    def test_the_bounds_themselves_are_accepted(self):
        assert MergeParams(alpha="1e-300").alpha_exact == Fraction(1, 10 ** 300)
        assert MergeParams(error_rate="1e-9").e == Fraction(1, 10 ** 9)
        assert MergeParams(error_rate="0.7499").e == Fraction(7499, 10 ** 4)

    def test_the_merge_record_is_exact_and_complete(self):
        p = MergeParams(alpha="1e-7", error_rate="0.0042", adapter_trimmed=True,
                        min_read_length=35, npolicy="random")
        rec = p.merge_record()
        assert rec == {
            "policy": "zna-merge-0.6", "zna_version": rec["zna_version"],
            "alpha": "0.0000001", "error_rate": "0.0042",
            "adapter_trimmed": True, "min_read_length": 35, "npolicy": "random",
        }
        import zna
        assert rec["zna_version"] == zna.__version__

    def test_the_version_that_writes_the_record_is_the_policys_release(self):
        """The 0.6 policy is a breaking release: a merge record, a stats JSON or a
        prologue written by this tree must not claim a 0.5.x version, or a consumer that
        pins versions would read 0.6 output as 0.5.3's."""
        import zna
        major, minor, patch = (int(x) for x in zna.__version__.split(".")[:3])
        assert (major, minor, patch) >= (0, 6, 0), zna.__version__

    def test_the_conda_recipe_carries_the_package_version(self):
        """Two places spell the version (pyproject reads ``__init__``); the recipe's
        sha256 is set at release, the version is not. Skipped where the recipe is not
        shipped (an sdist)."""
        import re
        import zna
        meta = Path(__file__).resolve().parents[1] / "conda" / "meta.yaml"
        if not meta.exists():
            pytest.skip("no conda/meta.yaml beside the tests")
        m = re.search(r'\{% set version = "([^"]+)" %\}', meta.read_text())
        assert m and m.group(1) == zna.__version__


# --------------------------------------------------------------------------- #
# 2. the run's diagnostics: the detected-overlap rate and the read-through check
# --------------------------------------------------------------------------- #

def _bufs(pairs, tag=b""):
    """Two FASTQ buffers from ``(s1, s2)`` pairs, in input order."""
    b1 = b"".join(b"@p%d/1%b\n%b\n+\n%b\n" % (i, tag, s1, qual(s1))
                  for i, (s1, _s2) in enumerate(pairs))
    b2 = b"".join(b"@p%d/2%b\n%b\n+\n%b\n" % (i, tag, s2, qual(s2))
                  for i, (_s1, s2) in enumerate(pairs))
    return b1, b2


def _chunk_out(pairs, p=_P, base=0, rt_check_pairs=0, backend="python"):
    """One reference-backend chunk over *pairs*; returns merge_chunk's whole result."""
    from zna.merge.backend import get_merge_backend
    b1, b2 = _bufs(pairs)
    out = get_merge_backend(backend).merge_chunk(
        b1, 0, len(b1), b2, 0, len(b2), *_chunk_args(p, lr=1, base=base), 1, 0,
        rt_check_pairs)
    assert out[8] == 0
    return out


def _counters(pairs, p=_P, base=0, rt_check_pairs=0, backend="python"):
    """One reference-backend chunk over *pairs*; returns its counters."""
    return _chunk_out(pairs, p, base, rt_check_pairs, backend)[3]


def _readthrough_pair(rng, frag_len=60):
    """A raw 2x150 pair whose fragment is shorter than the reads: a strong read-through,
    which is also the fragment -- it merges when read-through is allowed."""
    frag = draw(rng, frag_len)
    return ((frag + ADAPTER1 + draw(rng, 150))[:150],
            (rc(frag) + ADAPTER2 + draw(rng, 150))[:150])


class TestTheRunDiagnostics:
    """Two diagnostics ride along with the merge and never change a decision
    (plan §4). They are counted in the kernel, per pair, and summed like the other
    counters, so they need no buffer and cannot depend on how the input was chunked."""

    def test_the_detected_rate_counts_informative_positions_only(self):
        """An N against a call is neither a mismatch nor a compared base: it says
        nothing about whether the mates agree."""
        frag = rand_seq(200, 3)
        r1 = bytearray(frag[:150])
        r1[120] = ord("N")                          # inside the 100-base overlap
        r1[130] = flip(r1[130])                     # one real disagreement
        c = _counters([(bytes(r1), rc(frag[50:]))])
        assert (c[MERGED], c[DET_BASES], c[DET_MISMATCHES]) == (1, 99, 1)

    def test_n_against_n_is_not_a_compared_base_either(self):
        """Both mates reading N is a match in the scan, but it says no more about
        agreement than a one-sided N does: the detected rate leaves it out of the
        denominator. (The gate's dfit stays indexed by the overlap length, plan §2.)"""
        frag = rand_seq(200, 3)
        r1 = bytearray(frag[:150])
        r2 = bytearray(rc(frag[50:]))
        r1[120] = ord("N")
        r2[len(r2) - 1 - (120 - 50)] = ord("N")      # the same fragment position
        r1[130] = flip(r1[130])
        c = _counters([(bytes(r1), bytes(r2))])
        assert (c[DET_BASES], c[DET_MISMATCHES]) == (99, 1)
        # an all-N pair "overlaps" perfectly in the scan and adds nothing here
        c = _counters([(b"N" * 150, b"N" * 150)])
        assert c[DET_BASES] == c[DET_MISMATCHES] == 0

    def test_the_detected_rate_is_taken_before_the_gate(self):
        """A refused alignment is still a DETECTED overlap, and it counts: the rate is
        what `--error-rate` is checked against, and a rate measured only on what the
        gate let through could never exceed the rate the gate was built from."""
        frag = TestPlausibilityGate._repeat_beats_truth()
        L = len(frag)
        r1, r2 = frag[:150], rc(frag[L - 150:])
        o = find_overlap(r1, rc(r2), P)
        assert o.verdict == V_IMPLAUSIBLE
        out = _chunk_out([(r1, r2)])
        c, olen_hist, det_olen_hist = out[3], out[5], out[7]
        assert (c[MERGED], c[IMPLAUSIBLE]) == (0, 1)
        assert (c[DET_BASES], c[DET_MISMATCHES]) == (o.overlap_len,
                                                     o.informative_mismatches)
        # its length is binned among the detected overlaps, where dfit was looked up...
        assert det_olen_hist == [0] * o.overlap_len + [1]
        # ...while the post-admission statistics see nothing
        assert c[cli.SUM_OLEN] == c[cli.SUM_DIFF] == 0 and olen_hist == []

    def test_nothing_detected_counts_nothing(self):
        out = _chunk_out([(rand_seq(100, 4), rand_seq(100, 5)), (b"", b"ACGT")])
        c = out[3]
        assert c[N_PAIRS] == 2 and c[DET_BASES] == c[DET_MISMATCHES] == 0
        assert out[7] == []

    def test_the_readthrough_check_counts_the_unrestricted_winner(self):
        """Undeclared, the scan's own winner IS the unrestricted one; declared, the
        read-through side was never visited and the check scans it separately -- and
        counts the same pair either way."""
        rng = random.Random(8)
        rt_pair = _readthrough_pair(rng)
        normal = make_pair(draw(rng, 250), 150)
        pairs = [rt_pair, (normal[0][1], normal[1][1])]
        free = _counters(pairs, rt_check_pairs=10)
        declared = _counters(pairs, MergeParams(adapter_trimmed=True), rt_check_pairs=10)
        assert free[RT_STRONG] == declared[RT_STRONG] == 1
        assert (free[MERGED], declared[MERGED]) == (2, 1)
        # the check never touches a decision: with it off, the same outcomes
        off = _counters(pairs, MergeParams(adapter_trimmed=True), rt_check_pairs=0)
        assert off[RT_STRONG] == 0 and off[:RT_STRONG] == declared[:RT_STRONG]

    def test_the_readthrough_check_covers_the_first_pairs_of_the_INPUT(self):
        """Pairs are numbered from base_index, not from the chunk: a chunk that starts
        at pair 3 of a 4-pair check window checks exactly its first pair."""
        rng = random.Random(9)
        pairs = [_readthrough_pair(rng) for _ in range(5)]
        declared = MergeParams(adapter_trimmed=True)
        assert _counters(pairs, declared, base=0, rt_check_pairs=4)[RT_STRONG] == 4
        assert _counters(pairs, declared, base=3, rt_check_pairs=4)[RT_STRONG] == 1
        assert _counters(pairs, declared, base=4, rt_check_pairs=4)[RT_STRONG] == 0


def _brute_refusal(n, dfit_n, rate):
    """``P(Binom(n, rate) > dfit_n)`` by exact rational enumeration of every term --
    no recurrence, no rounding: what :func:`cli.refusal_probability` must equal."""
    return sum(Fraction(math.comb(n, k)) * rate ** k * (1 - rate) ** (n - k)
               for k in range(dfit_n + 1, n + 1))


class TestTheExpectedRefusedFraction:
    """The ``--error-rate`` check (plan §3, §8): the share of the run's detected
    overlaps the plausibility gate is expected to refuse were they all true and
    disagreeing at the detected rate, ``sum_n count[n] P(Binom(n, rate) > dfit[n]) /
    sum_n count[n]``. A diagnostic, computed once at the end of a run."""

    @pytest.mark.parametrize("rate", ["0.009", "0.01", "0.0123457", "0.03", "0.0355",
                                      "0.05", "0.25", "0.7"])
    def test_the_tail_equals_exact_enumeration(self, rate):
        """50 digits, like the tables: agreement to 1e-45, relative, over every n a
        2x150 run can detect, at the default dfit."""
        r = Fraction(rate)
        dfit = _P.dfit_table
        for n in range(1, 151):
            got = Fraction(cli.refusal_probability(n, dfit[n], r))
            want = _brute_refusal(n, dfit[n], r)
            assert abs(got - want) <= want * Fraction(1, 10 ** 45), (rate, n)

    def test_the_plans_closed_form_numbers(self):
        """§8 layer 1, and the --error-rate help text: a rule built for 1% refuses a
        true overlap ~1e-6 of the time at 1%, 0.07-0.6% at 3% and 1-13% at 5%, over
        overlaps of 50-150 bases."""
        dfit = _P.dfit_table

        def at(n, rate):
            return float(cli.refusal_probability(n, dfit[n], Fraction(rate)))
        assert all(1e-7 < at(n, "0.01") < 1e-6 for n in (50, 100, 150))
        assert round(100 * at(50, "0.03"), 2) == 0.07
        assert round(100 * at(150, "0.03"), 1) == 0.6
        assert round(100 * at(50, "0.05")) == 1 and round(100 * at(150, "0.05")) == 13

    def test_the_edges(self):
        assert cli.refusal_probability(10, 10, Fraction(1, 2)) == 0   # d > n impossible
        assert cli.refusal_probability(10, 3, Fraction(0)) == 0
        assert cli.refusal_probability(10, 3, Fraction(1)) == 1
        assert cli.refusal_probability(0, 0, Fraction(1, 2)) == 0
        assert cli.expected_refused_fraction(Fraction(1, 2), [], [0]) == 0
        assert cli.expected_refused_fraction(Fraction(1, 2), [0, 0], [0, 0]) == 0

    def test_the_fraction_is_the_histogram_weighted_mean(self):
        """Integer counts (a run) and float weights (the panel's projection) alike."""
        dfit = _P.dfit_table
        rate = Fraction(355, 10000)
        hist = [0] * 151
        hist[40], hist[97], hist[150] = 3, 11, 2
        want = sum(c * _brute_refusal(n, dfit[n], rate)
                   for n, c in enumerate(hist) if c) / sum(hist)
        got = cli.expected_refused_fraction(rate, hist, dfit)
        assert abs(Fraction(got) - want) <= want * Fraction(1, 10 ** 45)
        weighted = [c * 2.5 for c in hist]                   # a uniform weight cancels
        assert abs(Fraction(cli.expected_refused_fraction(rate, weighted, dfit)) - want) \
            <= want * Fraction(1, 10 ** 45)

    def test_a_uniform_run_reduces_to_one_tail(self):
        """Every detected overlap the same length: the fraction IS that length's
        refusal probability, whatever the count."""
        acc = cli._new_acc()
        acc[0][DET_BASES], acc[0][DET_MISMATCHES] = 60 * 50, 240          # 8%
        acc[4].extend([0] * 50 + [60])
        assert cli.run_refused_fraction(acc, MergeParams()) == \
            cli.refusal_probability(50, 6, Fraction(8, 100))
        assert cli.run_refused_fraction(cli._new_acc(), MergeParams()) == 0


# --------------------------------------------------------------------------- #
# 3. backend selection
# --------------------------------------------------------------------------- #

class TestExtensionsAreDistinct:
    """`zna._accel` and `zna.merge._accel` are two different extensions.

    Both are imported as `_accel` within their package, so both CMake targets emit a
    file called `_accel.cpython-*.so`. Without separate build output directories they
    collide and one overwrites the other — which happened during development: the merge
    scan ended up installed as `zna._accel`, `zna.is_accelerated()` returned False, and
    the codec silently fell back to pure Python. Nothing in the merge suite noticed,
    because the merge suite was perfectly happy.
    """

    # `exc_type=ImportError` because "the extension is not there" and "the extension is
    # there but will not load" are the same fact for these tests, and only the first is
    # a ModuleNotFoundError. A half-built environment -- a stale .so against a newer
    # interpreter, a missing runtime dependency -- raises the second, and without this
    # the tests error out instead of skipping. It is also required from pytest 9.1.
    def test_the_codec_extension_is_the_codec(self):
        accel = pytest.importorskip("zna._accel", exc_type=ImportError)
        assert hasattr(accel, "encode_block"), "zna._accel is not the codec extension"
        assert not hasattr(accel, "scan")

    def test_the_merge_extension_is_the_scan(self):
        accel = pytest.importorskip("zna.merge._accel", exc_type=ImportError)
        assert hasattr(accel, "scan"), "zna.merge._accel is not the scan extension"
        assert not hasattr(accel, "encode_block")

    def test_the_codec_is_still_accelerated(self):
        """The whole point of the package. A build that quietly loses this passes every
        functional test it has, which is why the conda recipe asserts it too."""
        pytest.importorskip("zna._accel", exc_type=ImportError)
        import zna
        assert zna.is_accelerated()


class TestBackendSelection:
    """The kernel is selectable, the same way the codec's is.

    The Python backend is the reference oracle rather than a fallback, so it must be
    available unconditionally — an environment where only the accelerated one loads
    would have nothing to check the accelerated one against.
    """

    def test_the_reference_backend_is_always_available(self):
        from zna.merge.backend import available_merge_backends
        assert "python" in available_merge_backends()

    def test_auto_prefers_accel_when_it_is_built(self, merge_backend_option):
        from zna.merge.backend import available_merge_backends, get_merge_backend_name
        if merge_backend_option == "python":
            # tests/conftest.py narrowed the preference on purpose: configuration 2.
            assert get_merge_backend_name() == "python"
            return
        expected = "accel" if "accel" in available_merge_backends() else "python"
        assert get_merge_backend_name() == expected

    def test_the_session_runs_the_backend_it_was_asked_for(self, merge_backend_option):
        """What ``--merge-backend`` (tests/conftest.py) promises: under ``python``,
        ``zna merge --backend auto`` resolves to the reference kernel on every run, and
        an explicit ``accel`` still loads the compiled one where it is built."""
        from zna.merge import backend
        if merge_backend_option != "python":
            pytest.skip("configuration 2 only (--merge-backend=python)")
        assert backend._PREFERENCE == ("python",)
        assert backend.use("auto") == "python" and backend.active_name() == "python"
        if "accel" in backend.available_merge_backends():
            assert backend.get_merge_backend_name("accel") == "accel"

    def test_an_unknown_backend_is_a_loud_error(self):
        from zna.merge.backend import get_merge_backend
        with pytest.raises(ImportError, match="unknown merge backend"):
            get_merge_backend("hopeful")

    def test_selection_round_trips_and_restores(self):
        from zna.merge import overlap
        original = overlap.backend_name()
        try:
            assert overlap.use_backend("python") == "python"
            frag = rand_seq(40, 3)
            (_, s1, _), (_, s2, _) = make_pair(frag, 30)
            assert find_overlap(s1, rc(s2)).verdict == V_MERGE
        finally:
            overlap.use_backend(original)
        assert overlap.backend_name() == original

    def test_a_backend_built_for_another_policy_is_refused(self, monkeypatch):
        """A 0.5.x extension takes different arguments under the same names. Calling it
        with 0.6's would be a crash at best and a silently different corpus at worst,
        so it is refused like a missing build -- by its ABI marker, or its absence."""
        import types
        from zna.merge import backend
        from zna.merge import _pymerge
        stale = types.ModuleType("zna_merge_stale_test")
        for name in backend._REQUIRED_FUNCTIONS:
            setattr(stale, name, getattr(_pymerge, name))       # everything but the ABI
        monkeypatch.setitem(sys.modules, "zna_merge_stale_test", stale)
        monkeypatch.setitem(backend._BACKEND_MODULES, "stale", "zna_merge_stale_test")
        with pytest.raises(ImportError, match="ABI None"):
            backend.get_merge_backend("stale")
        assert "stale" not in backend.available_merge_backends()
        stale.POLICY_ABI = backend.POLICY_ABI + 1
        with pytest.raises(ImportError, match="ABI"):
            backend.get_merge_backend("stale")

    def test_the_installed_extension_is_used_only_if_it_speaks_this_abi(self):
        from zna.merge.backend import POLICY_ABI, available_merge_backends
        try:
            import zna.merge._accel as accel
        except ImportError:
            pytest.skip("no compiled merge extension")
        speaks = getattr(accel, "POLICY_ABI", None) == POLICY_ABI
        assert ("accel" in available_merge_backends()) is speaks


# --------------------------------------------------------------------------- #
# 4. cross-backend equivalence: the accelerated kernel must agree EXACTLY
# --------------------------------------------------------------------------- #

def _backend_pair(fn):
    """``fn`` from the reference and the compiled backend, or skip."""
    from zna.merge.backend import available_merge_backends, get_merge_backend
    if "accel" not in available_merge_backends():
        pytest.skip("no C++ merge backend built for this policy ABI")
    return getattr(get_merge_backend("python"), fn), getattr(get_merge_backend("accel"), fn)


def _backends():
    return _backend_pair("scan")


def _chunk_backends():
    """The two `merge_chunk` implementations, for the level-3 differential."""
    return _backend_pair("merge_chunk")


@pytest.fixture(params=["python", "accel"])
def any_backend(request):
    """Run a test once per available backend.

    Without this, anything that goes through ``find_overlap`` silently tests only
    whichever backend is *selected* — which is ``accel`` wherever it is built, leaving
    the reference oracle unchecked in exactly the environment that ships. The oracle is
    what the accelerated kernel is defined to agree with, so an unchecked oracle makes
    the whole cross-backend suite circular.
    """
    from zna.merge import overlap
    from zna.merge.backend import available_merge_backends
    if request.param not in available_merge_backends():
        pytest.skip(f"{request.param} merge backend not available")
    original = overlap.backend_name()
    overlap.use_backend(request.param)
    try:
        yield request.param
    finally:
        overlap.use_backend(original)


#: A deliberately LOW floor for the scan differentials: 8 bits, so that far more
#: shifts survive to be compared than the policy's ~28 would leave.
LOW_FLOOR = to_q(8)


def _chunk_args(p=_P, lr=40, check_sync=True, base=0):
    return (*p.kernel_args(), lr, DISAGREE_Q, check_sync, base)


def exhaustive_scan(s1, s2rc, floor_q, adapter_trimmed=False):
    """Every eligible shift, no pruning — the slow truth to check a scan against."""
    out = []
    len1, len2 = len(s1), len(s2rc)
    for s in range(-(len2 - 1), len1):
        if adapter_trimmed and s + len2 < max(len1, len2):
            continue
        lo, hi = max(s, 0), min(len1, s + len2)
        n = hi - lo
        if n <= 0:
            continue
        off = lo - s
        d = sum(s1[lo + k] != s2rc[off + k] for k in range(n))
        sc = score_of(n - d, d)
        if sc >= floor_q:
            out.append((s, n, d, sc))
    return out


def argmax_by_rule(scored):
    """The specified order: maximise score, then minimise s. Returns (winner, n_ties)."""
    top = max(t[3] for t in scored)
    ties = [t for t in scored if t[3] == top]
    return min(ties, key=lambda t: t[0]), len(ties)


#: Periodic content on unequal-length mates: the plateau then holds several shifts of
#: equal overlap and equal mismatch count, which is the only way ties arise in practice.
#: Random sequence essentially never ties -- a sweep over 7,000 random and adversarial
#: pairs produced exactly zero -- so any tie-break assertion must build them like this.
def tie_fixtures():
    # (a) PLATEAU ties: unequal-length mates, so several maximal-overlap shifts exist,
    #     and periodic content makes some of them score identically.
    for period in (b"CA", b"ACG", b"AT", b"A", b"CAG", b"ACGT"):
        for l1 in range(20, 46, 3):
            for l2 in range(20, 46, 3):
                seq = period * 50
                yield seq[:l1], seq[:l2], f"plateau-{period.decode()}-{l1}x{l2}"

    # (b) FLANK ties: equal-length mates, periodic, mutually OUT of phase. The plateau
    #     (s=0) then mismatches everywhere and is rejected, while the two flanks at the
    #     same overlap length both come into phase — s = -k and s = +k, tied at the top,
    #     one read-through and one normal. This is the only construction that
    #     distinguishes the two flanks' visiting order, and without it swapping them
    #     passes every other test in this file while changing the winner on every tied
    #     pair. Verified to produce shifts {-1, +1} for CA/AC and {-2, +2} for ACGT/GTAC.
    for period, rot in ((b"CA", b"AC"), (b"ACGT", b"GTAC"), (b"AATT", b"TTAA")):
        for L in (24, 30, 40, 41, 50):
            yield (period * 40)[:L], (rot * 40)[:L], f"flank-{period.decode()}-L{L}"


class TestPopcount:
    """The 16-bit popcount under the SSE2 kernel, and why it is tested from everywhere.

    `neq16`'s x86 path used to count equal lanes with a popcount. That was
    `__builtin_popcount`, which is a GCC/Clang extension — MSVC does not have it, and the
    first Windows build of this extension failed on exactly that line with C3861. It went
    unnoticed because `zna merge` is new in 0.4.0, so no MSVC had ever compiled the file,
    and because an arm64 developer machine takes the NEON path, where the popcount does
    not appear at all.

    **Since 0.5.2 the scan has no popcount at all**: the x86 path reduces with `psadbw`,
    which is baseline SSE2, because the same `__builtin_popcount` turned out to compile to
    `callq __popcountdi2@plt` on a baseline build and cost the kernel 1.88x. The portable
    fold is still compiled on *every* platform and still checked here — it is the
    primitive a future kernel should reach for, and the ways of getting a popcount that
    *look* free are exactly the ones that are not. MSVC's `__popcnt16` remains deliberately
    unused: it emits the POPCNT instruction, which is not baseline x86-64, and would turn a
    build error into an illegal-instruction fault on an older CPU.
    """

    def _fn(self):
        accel = pytest.importorskip("zna.merge._accel", exc_type=ImportError)
        if not hasattr(accel, "_popcount16_portable"):
            pytest.skip("extension predates the portable popcount")
        return accel._popcount16_portable

    def test_exhaustive_over_every_16_bit_input(self):
        """65,536 inputs is small enough to check all of them, so check all of them."""
        fn = self._fn()
        bad = [x for x in range(1 << 16) if fn(x) != bin(x).count("1")]
        assert not bad, f"{len(bad)} mismatches, first at {bad[:5]}"

    def test_the_boundaries_it_would_plausibly_get_wrong(self):
        """A SWAR fold fails at the carry boundaries or nowhere; name them anyway."""
        fn = self._fn()
        assert fn(0x0000) == 0 and fn(0xFFFF) == 16
        assert fn(0x00FF) == 8 and fn(0xFF00) == 8      # the byte-sum step
        assert fn(0x5555) == 8 and fn(0xAAAA) == 8      # the pair step
        assert fn(0x8000) == 1 and fn(0x0001) == 1

    def test_bits_above_16_are_ignored(self):
        """It is documented as a 16-bit count and `neq16` relies on that: the movemask
        result is 16 bits, but nothing stops a caller passing a wider value."""
        fn = self._fn()
        assert fn(0xFFFF0000) == 0
        assert fn(0xDEAD_FFFF) == 16


def _library(rng, n, lmin=30, lmax=151, err=0.01, with_n=0.0, tags=b""):
    """A realistic-ish chunk: fragments 40-320, independent mate lengths, varied Q."""
    r1s, r2s = [], []
    for i in range(n):
        frag = draw(rng, rng.randrange(40, 320))
        l1, l2 = rng.randrange(lmin, lmax), rng.randrange(lmin, lmax)
        s1 = bytearray(mutate((frag + ADAPTER1 + draw(rng, 160))[:l1], rng, err))
        s2 = bytearray(mutate((rc(frag) + ADAPTER2 + draw(rng, 160))[:l2], rng, err))
        for sb in (s1, s2):
            if sb and rng.random() < with_n:
                sb[rng.randrange(len(sb))] = ord("N")
        q1 = bytes(rng.choice((70, 58, 44, 35)) for _ in range(len(s1)))
        q2 = bytes(rng.choice((70, 58, 44, 35)) for _ in range(len(s2)))
        r1s.append(b"@f%d/1%b\n%b\n+\n%b\n" % (i, tags, bytes(s1), q1))
        r2s.append(b"@f%d/2%b\n%b\n+\n%b\n" % (i, tags, bytes(s2), q2))
    return b"".join(r1s), b"".join(r2s)


class TestCrossBackend:
    """The oracle and the accelerated kernel are one algorithm with two implementations.

    Equality here is exact, not approximate, and that is only possible because the score
    is an integer (params.py) and the argmax is a specified total order rather than an
    artifact of iteration order (docs/METHODS.md). Asserting the weaker
    "returned *an* argmax" would let a tie-break divergence through, and a tie-break
    divergence changes which bases a merged read is built from.
    """

    def _agree(self, s1, s2rc, label, floor=LOW_FLOOR):
        py, cc = _backends()
        out = None
        for at in (0, 1):
            a = py(s1, s2rc, len(s1), len(s2rc), _P.match_q, _P.step_q, floor, at)
            b = cc(s1, s2rc, len(s1), len(s2rc), _P.match_q, _P.step_q, floor, at)
            assert a == b, (label, at, a, b)
            out = a if at == 0 else out
        return out

    def test_overlapping_pairs(self):
        rng = random.Random(11)
        for i in range(400):
            frag = draw(rng, rng.randrange(40, 320))
            l1 = rng.randrange(20, 151)
            l2 = rng.randrange(20, 151)
            r1 = mutate((frag + ADAPTER1 + draw(rng, 160))[:l1], rng, 0.01)
            r2 = mutate((rc(frag) + ADAPTER2 + draw(rng, 160))[:l2], rng, 0.01)
            self._agree(r1, rc(r2), f"ovl{i}")
            self._agree(r1, rc(r2), f"ovl{i}@T", floor=t_q(l1, l2))

    def test_unrelated_pairs_exercise_the_rejection_path(self):
        """Where the scan spends nearly all its time: every shift bails early."""
        rng = random.Random(12)
        for i in range(400):
            self._agree(draw(rng, rng.randrange(1, 200)),
                        draw(rng, rng.randrange(1, 200)), f"unrel{i}")

    @pytest.mark.parametrize("s1,s2rc", [
        (b"", b""), (b"", b"ACGT"), (b"ACGT", b""), (b"A", b"A"), (b"A", b"T"),
        (b"ACGT" * 4, b"ACGT" * 4),                       # exactly one vector wide
        (b"ACGT" * 8, b"ACGT" * 8),                       # exactly one bail block
        (b"ACGT" * 8 + b"A", b"ACGT" * 8 + b"A"),         # one past a bail block
        (b"ACGT" * 4 + b"A", b"ACGT" * 4 + b"A"),         # one past a vector
        (b"N" * 40, b"N" * 40),                           # N vs N is a match
        (b"N" * 40, b"A" * 40),                           # N vs A is not
        (b"RYKMSWBDHVN" * 4, b"RYKMSWBDHVN" * 4),         # IUPAC compares as itself
        (b"acgt" * 10, b"ACGT" * 10),                     # case is significant
    ])
    def test_edges_and_non_acgt(self, s1, s2rc):
        """Byte comparison IS the reference semantics, so none of this needs a special
        path — but that claim is exactly what a packed kernel would break, so pin it."""
        self._agree(s1, s2rc, repr((s1[:12], s2rc[:12])))

    def test_a_mismatch_at_every_position_of_a_bail_block(self):
        """The vector loop's own block-loop test.

        `test_block_loop_sees_every_position` guards the reference's 8-wide unrolled
        loop; this guards the accelerated 16-byte/32-base one. A stride bug, a wrong
        vector boundary, or a tail that starts one byte late mis-scores only overlaps
        whose mismatch sits at particular offsets — so sweep it across a vector
        boundary, a bail-block boundary, and into the scalar tail.
        """
        rng = random.Random(13)
        for n in (16, 31, 32, 33, 40, 47, 48, 49, 64, 65, 150):
            frag = draw(rng, n)
            for i in range(n):
                r1 = bytearray(frag)
                r1[i] = flip(r1[i])
                got = self._agree(bytes(r1), frag, f"n{n}pos{i}")
                assert got[3] == 1, (n, i, got)      # exactly one mismatch, found

    def test_lengths_around_the_vector_and_block_boundaries(self):
        """Reads whose length lands on, either side of, and between the 16- and 32-byte
        boundaries — where an off-by-one in the tail would hide."""
        rng = random.Random(14)
        for l1 in range(1, 70):
            for l2 in (l1, max(1, l1 - 1), l1 + 1, 16, 32, 33):
                frag = draw(rng, l1 + l2)
                self._agree(frag[:l1], rc(rc(frag[-l2:])), f"{l1}x{l2}")

    def test_deliberate_ties_agree_and_follow_the_specified_order(self):
        """The test this class was missing, and the one that matters most.

        Everything above is built from random or adversarial sequence, which essentially
        never ties — so *reversing the flank visiting order in the C++ kernel passes
        every other test in this class*, while silently changing which shift wins on
        every tied pair, and with it which bases a merged read is built from.

        Ties need periodic content on unequal-length mates. Here both backends must not
        only agree with each other but land on the specified winner: maximise score,
        then minimise s -- over every shift, and over the contract's eligible ones.
        """
        py, cc = _backends()
        n_tied = 0
        for at in (0, 1):
            for s1, s2rc, label in tie_fixtures():
                args = (len(s1), len(s2rc), _P.match_q, _P.step_q, LOW_FLOOR, at)
                a = py(s1, s2rc, *args)
                b = cc(s1, s2rc, *args)
                assert a == b, (label, at, a, b)
                scored = exhaustive_scan(s1, s2rc, LOW_FLOOR, bool(at))
                if not scored:
                    assert a[2] == 0, label
                    continue
                (want_s, want_n, want_d, want_sc), ties = argmax_by_rule(scored)
                assert a == (want_s, want_sc, want_n, want_d), (label, at, a, ties)
                n_tied += ties - 1
        assert n_tied >= 100, f"only {n_tied} ties exercised; the tie-break is untested"

    def test_random_bytes_not_just_nucleotides(self):
        """The kernel promises raw byte semantics; hold it to that on arbitrary input."""
        rng = random.Random(15)
        for i in range(300):
            n1, n2 = rng.randrange(0, 80), rng.randrange(0, 80)
            self._agree(bytes(rng.randrange(256) for _ in range(n1)),
                        bytes(rng.randrange(256) for _ in range(n2)), f"bytes{i}")

    def test_the_overlap_decision_agrees(self):
        """`overlap` -- scan, floor lookup, informative count, gate -- on pairs built to
        land in all three verdicts, under both contracts and several error rates."""
        py, cc = _backend_pair("overlap")
        rng = random.Random(21)
        seen = set()
        for e in ("0.01", "0.0003", "0.03"):
            for at in (False, True):
                p = MergeParams(error_rate=e, adapter_trimmed=at)
                for i in range(300):
                    unit = draw(rng, rng.randrange(8, 40))
                    frag = (unit * 20)[:rng.randrange(40, 320)]      # repeat-rich
                    frag = mutate(frag, rng, 0.08)
                    l1, l2 = rng.randrange(20, 151), rng.randrange(20, 151)
                    s1 = bytearray(mutate((frag + ADAPTER1 + draw(rng, 160))[:l1],
                                          rng, 0.02))
                    s2 = (rc(frag) + ADAPTER2 + draw(rng, 160))[:l2]
                    if rng.random() < 0.3 and l1 > 10:
                        k = rng.randrange(l1 - 5)
                        s1[k:k + 5] = b"NNNNN"
                    s1, s2rc = bytes(s1), rc(s2)
                    args = (len(s1), len(s2rc), *p.kernel_args())
                    a, b = py(s1, s2rc, *args), cc(s1, s2rc, *args)
                    assert a == b, (e, at, i, a, b)
                    seen.add(a[0])
        assert seen == {0, 1, 2}, f"the fixture only reached verdicts {seen}"

    def test_the_diagnostics_agree(self):
        """The detected-overlap counters and the read-through check, with the check
        window ending mid-chunk and mid-input, under both contracts, on input with N
        runs (the informative count) and read-through (the check)."""
        py, cc = _chunk_backends()
        rng = random.Random(22)
        buf1, buf2 = _library(rng, 300, lmin=20, with_n=0.3)
        for at in (False, True):
            for base, window in ((0, 250), (100, 250), (0, 0)):
                args = (*_chunk_args(MergeParams(adapter_trimmed=at), base=base), 1, 0,
                        window)
                a = py(buf1, 0, len(buf1), buf2, 0, len(buf2), *args)
                b = cc(buf1, 0, len(buf1), buf2, 0, len(buf2), *args)
                assert a[3] == b[3], (at, base, window)
                assert a[0] == b[0]
                if window:
                    assert a[3][RT_STRONG] > 0, "no read-through: not exercised"
                assert a[3][DET_BASES] > 0 and a[3][DET_MISMATCHES] > 0

    def test_chunks_agree_blob_for_blob(self):
        """Level 3: the production path. Same bytes out, same counters, same histograms.

        This is the strongest check available — it needs no model of what the answer
        should be, only that two independent implementations of one specification agree
        completely on a realistic input.
        """
        py, cc = _chunk_backends()
        buf1, buf2 = _library(random.Random(17), 400)
        for p in (_P, MergeParams(error_rate="0.03", adapter_trimmed=True)):
            args = (*_chunk_args(p), 1, 0, 300)      # read-through check on 300 of 400
            a = py(buf1, 0, len(buf1), buf2, 0, len(buf2), *args)
            b = cc(buf1, 0, len(buf1), buf2, 0, len(buf2), *args)
            assert a[0] == b[0], "blobs differ"
            assert a[1:3] == b[1:3], "consumed byte counts differ"
            assert a[3] == b[3], "counters differ"
            assert [list(x) for x in a[4:8]] == [list(x) for x in b[4:8]], \
                "histograms differ"
            assert a[8] == b[8] == 0
            # a detected overlap is merged or refused (no N here, so trim3 demotes none)
            assert sum(a[7]) == a[3][MERGED] + a[3][IMPLAUSIBLE]
            # This library carries raw adapter read-through, which the declaration makes
            # unmergeable: 169 merged undeclared, 79 declared.
            assert a[3][N_PAIRS] == 400, a[3]
            assert a[3][MERGED] > (50 if p.adapter_trimmed else 100), a[3]
            assert a[3][BASES_CONSENSUS] > 0, "no consensus changes: not exercised"

    def test_record_adapter_agrees_with_merge_chunk(self):
        """The record adapter and the FASTQ adapter share one inner loop; this
        holds them to it: same sequences record for record, same slot the
        FASTQ names imply, same consumed counts, all counters, all four
        histograms.  (MERGE_PAIRS_PLAN.md §4 step 1.)"""
        from zna.merge.backend import available_merge_backends, get_merge_backend
        if "accel" not in available_merge_backends():
            pytest.skip("no C++ merge backend built for this policy ABI")
        for name in ("python", "accel"):
            be = get_merge_backend(name)
            rng = random.Random(23)
            r1s, r2s = [], []
            for i in range(300):
                frag = draw(rng, rng.randrange(40, 320))
                l1, l2 = rng.randrange(30, 151), rng.randrange(30, 151)
                s1 = mutate((frag + ADAPTER1 + draw(rng, 160))[:l1], rng, 0.01)
                s2 = mutate((rc(frag) + ADAPTER2 + draw(rng, 160))[:l2], rng, 0.01)
                q1 = bytes(rng.choice((70, 58, 44, 35)) for _ in range(len(s1)))
                q2 = bytes(rng.choice((70, 58, 44, 35)) for _ in range(len(s2)))
                r1s.append(b"@f%d/1 XA:i:%d\n%b\n+\n%b\n" % (i, i, s1, q1))
                r2s.append(b"@f%d/2 XB:i:%d\n%b\n+\n%b\n" % (i, i, s2, q2))
            buf1, buf2 = b"".join(r1s), b"".join(r2s)
            args = _chunk_args()
            blob, c1, c2, counters, lh, oh, ih, dh, need = be.merge_chunk(
                buf1, 0, len(buf1), buf2, 0, len(buf2), *args, 1, 0, 200)
            seqs, ends, rc1, rc2, rcounters, rlh, roh, rih, rdh, rneed = \
                be.merge_chunk_records(buf1, 0, len(buf1), buf2, 0, len(buf2), *args,
                                       True, 1, 0, 200)
            assert (c1, c2, need) == (rc1, rc2, rneed)
            assert counters == rcounters
            assert (list(lh), list(oh), list(ih), list(dh)) == \
                (list(rlh), list(roh), list(rih), list(rdh))
            fastq = [ln for ln in blob.split(b"\n")[1::4] if ln]
            recs = [seqs[o:o + l] for (o, l, _ho, _hl, _slot, _p) in ends]
            assert fastq == recs
            names = [ln[1:] for ln in blob.split(b"\n")[0::4] if ln]
            for name_, (o, l, ho, hl, slot, prov) in zip(names, ends):
                implied = (1 if b"/1" in name_.split(b" ")[0]
                           else 2 if b"/2" in name_.split(b" ")[0] else 0)
                assert slot == implied, (name_, slot)
                src = buf2 if slot == 2 else buf1
                hdr = src[ho:ho + hl]
                assert hdr.startswith(b"f") and (b"XB:" in hdr if slot == 2
                                                 else b"XA:" in hdr)

    def test_record_chunks_agree_across_backends(self):
        """Cross-backend differential for the record adapter: seqs blob, ends
        (offsets, slots, prov bytes), consumed counts, counters, histograms --
        element for element.  (MERGE_PAIRS_PLAN.md §4 step 2.)"""
        py, cc = _backend_pair("merge_chunk_records")
        buf1, buf2 = _library(random.Random(29), 250, lmin=20, with_n=0.3, tags=b" t")
        for npolicy in (1, 0, 2):
            args = (*_chunk_args(), True, npolicy, 11, 200)
            a = py(buf1, 0, len(buf1), buf2, 0, len(buf2), *args)
            b = cc(buf1, 0, len(buf1), buf2, 0, len(buf2), *args)
            assert a[0] == b[0], f"npolicy={npolicy}: seq blobs differ"
            assert [tuple(e) for e in a[1]] == [tuple(e) for e in b[1]], \
                f"npolicy={npolicy}: ends differ"
            assert a[2:5] == b[2:5], f"npolicy={npolicy}: consumed/counters differ"
            assert [list(x) for x in a[5:9]] == [list(x) for x in b[5:9]]
            assert a[9] == b[9]
            assert any(e[5] for e in a[1]), f"npolicy={npolicy}: no prov bits set"

    @pytest.mark.parametrize("npolicy", [1, 0, 2])   # trim3, keep, random
    def test_chunks_with_no_calls_agree_blob_for_blob(self, npolicy):
        """The differential, on input that actually contains `N`.

        The suite had no such fixture, and it cost twice: an N-rescue that was wrong in
        the compiled backend only, and a rescue counter that double-counted. Both passed
        every cross-backend test that existed. N runs are also what the gate's
        informative count exists for.
        """
        py, cc = _chunk_backends()
        rng = random.Random(19)
        r1s, r2s = [], []
        for i in range(400):
            frag = draw(rng, rng.randrange(40, 320))
            l1, l2 = rng.randrange(30, 151), rng.randrange(30, 151)
            a = bytearray((frag + ADAPTER1 + draw(rng, 160))[:l1])
            b = bytearray((rc(frag) + ADAPTER2 + draw(rng, 160))[:l2])
            for _ in range(rng.randrange(0, 4)):
                if a: a[rng.randrange(len(a))] = ord("N")
                if b: b[rng.randrange(len(b))] = ord("N")
            if rng.random() < 0.2 and len(a) > 20:             # an N RUN
                k = rng.randrange(len(a) - 12)
                a[k:k + 12] = b"N" * 12
            q1 = bytes(rng.choice((70, 58, 44, 35)) for _ in range(len(a)))
            q2 = bytes(rng.choice((70, 58, 44, 35)) for _ in range(len(b)))
            r1s.append(b"@f%d/1\n%b\n+\n%b\n" % (i, bytes(a), q1))
            r2s.append(b"@f%d/2\n%b\n+\n%b\n" % (i, bytes(b), q2))
        buf1, buf2 = b"".join(r1s), b"".join(r2s)
        args = (*_chunk_args(), npolicy, 42)
        a = py(buf1, 0, len(buf1), buf2, 0, len(buf2), *args)
        b = cc(buf1, 0, len(buf1), buf2, 0, len(buf2), *args)
        assert a[0] == b[0], "blobs differ on input containing N"
        assert a[1:3] == b[1:3] and a[3] == b[3], f"counters differ: {a[3]} vs {b[3]}"
        assert [list(x) for x in a[4:8]] == [list(x) for x in b[4:8]]
        emitted = b"".join(a[0].split(b"\n")[1::4])
        assert (b"N" not in emitted) == bool(npolicy), "trim3 must leave no N"

    @pytest.mark.parametrize("npolicy", [1, 0, 2])   # trim3, keep, random
    def test_provenance_tokens_agree_across_backends(self, npolicy):
        """The header tokens are built in two places and must agree byte for byte.

        `merge_core.hpp::build_name` and `_pymerge._prov_name` are one specification with
        two implementations, exactly like the scan. This fixture drives every token —
        rescued, trim3/subn — through both, and compares the HEADERS specifically so a
        failure names the token rather than "blobs differ".
        """
        py, cc = _chunk_backends()
        rng = random.Random(23)
        r1s, r2s = [], []
        for i in range(400):
            frag = draw(rng, rng.randrange(40, 320))
            l1, l2 = rng.randrange(60, 151), rng.randrange(60, 151)
            a = bytearray(mutate((frag + ADAPTER1 + draw(rng, 160))[:l1], rng, 0.01))
            b = bytearray(mutate((rc(frag) + ADAPTER2 + draw(rng, 160))[:l2], rng, 0.01))
            # N in the overlap gets rescued; N past it survives to meet the policy.
            for _ in range(rng.randrange(0, 3)):
                if a: a[rng.randrange(len(a))] = ord("N")
                if b: b[rng.randrange(len(b))] = ord("N")
            q1 = bytes(rng.choice((70, 58, 44, 35)) for _ in range(len(a)))
            q2 = bytes(rng.choice((70, 58, 44, 35)) for _ in range(len(b)))
            r1s.append(b"@f%d/1\tZI:i:%d\n%b\n+\n%b\n" % (i, i, bytes(a), q1))
            r2s.append(b"@f%d/2\tZI:i:%d\n%b\n+\n%b\n" % (i, i, bytes(b), q2))
        buf1, buf2 = b"".join(r1s), b"".join(r2s)
        args = (*_chunk_args(), npolicy, 42)
        a = py(buf1, 0, len(buf1), buf2, 0, len(buf2), *args)
        b = cc(buf1, 0, len(buf1), buf2, 0, len(buf2), *args)

        ha = a[0].split(b"\n")[0::4]
        hb = b[0].split(b"\n")[0::4]
        assert ha == hb, "provenance tokens differ between backends"
        assert a[0] == b[0] and a[3] == b[3]

        seen = set()
        for h in ha:
            for tok in h.split()[1:]:
                if b"_" in tok and not tok.startswith(b"merged_"):
                    seen.add(tok.split(b"_")[0])
        want = {b"rescued"} | ({b"subn"} if npolicy == 2 else
                               {b"trim3"} if npolicy == 1 else set())
        assert want <= seen, f"fixture produced only {seen}"
        assert any(t.startswith(b"ZN:i:") for h in ha for t in h.split()), "no ZN tag"

    @pytest.mark.parametrize("payload", [
        b"",                                             # empty buffer
        b"@r1\nACGT\n+\nIIII\n",                         # a single complete record
        b"@r1\nACGT\n+\nIIII\n@r2\nAC",                  # trailing partial record
        b"@r1\r\nACGT\r\n+\r\nIIII\r\n",                 # CRLF must not reach the output
        b"@r1\nacgt\n+\nIIII\n",                         # lower case is upper-cased
    ])
    def test_parsing_edges_agree(self, payload):
        """The parser is where the audit's prototype had its four defects; CRLF
        surviving into the sequence was one of them."""
        py, cc = _chunk_backends()
        args = _chunk_args(lr=1, check_sync=False)
        a = py(payload, 0, len(payload), payload, 0, len(payload), *args)
        b = cc(payload, 0, len(payload), payload, 0, len(payload), *args)
        assert a[0] == b[0] and a[1:4] == b[1:4], (a, b)
        assert b"\r" not in a[0], "CRLF leaked into the output"

    @pytest.mark.parametrize("readlen", [150, 1023, 1024, 1025, 1600, 4100])
    def test_reads_longer_than_the_arena(self, readlen):
        """There is no read-length limit, and the arena has to grow to whatever turns up.

        It did not: the compiled backend sized its scratch to 1024 bases *before*
        parsing and threw "read longer than the scratch buffer" at 1025, while the
        reference merged the same pair happily — the two backends silently disagreed on
        an entire class of input. Sweep the boundary in both directions (with policy
        tables sized for it, which is the driver's job -- see TestTableGrowth).
        """
        py, cc = _chunk_backends()
        rng = random.Random(1000 + readlen)
        frag = draw(rng, readlen * 3 // 2)
        s1, s2 = frag[:readlen], rc(frag[-readlen:])
        b1 = b"@x/1\n%b\n+\n%b\n" % (s1, b"I" * len(s1))
        b2 = b"@x/2\n%b\n+\n%b\n" % (s2, b"I" * len(s2))
        p = MergeParams(error_rate="0.01")
        p.ensure(readlen)
        args = _chunk_args(p)
        a = py(b1, 0, len(b1), b2, 0, len(b2), *args)
        b = cc(b1, 0, len(b1), b2, 0, len(b2), *args)
        assert a[0] == b[0] and a[3] == b[3], (readlen, a[3], b[3])
        assert a[3][MERGED] == 1, f"the fixture stopped merging at {readlen}"
        assert a[3][MAX_READ_LEN] == readlen, "max_read_length is not being reported"

        # The histograms are uncapped too, and both backends bin identically. They were
        # `uint32_t[1025]` with every index clamped to the last bin, so a 2400 bp merged
        # record was counted as 1024 and the distributions silently aggregated at
        # exactly the length where the arena fix had just made long reads work.
        assert [list(x) for x in a[4:8]] == [list(x) for x in b[4:8]], "histograms"
        L = len(frag)                                  # the merged record IS the fragment
        assert a[4][L] == 1 and len(a[4]) == L + 1, "length histogram is clamped"
        assert a[6][L] == 1, "insert histogram is clamped"
        assert a[5][2 * readlen - L] == 1, "overlap histogram is clamped"
        assert a[7] == a[5], "the detected-overlap histogram is clamped"

    @pytest.mark.parametrize("hdrlen", [64, 1024, 1088, 2000, 16000, 70000])
    def test_headers_longer_than_the_read_arena(self, hdrlen):
        """A merged record's name is built from the HEADER; its buffer was sized from the
        READ arena.

        `Scratch::name` was resized inside `ensure()` as `cap + 64`, where `cap` is the
        read-length arena (minimum 1024). But the name is R1's header, pair suffix
        stripped, plus " merged_<n1>_<n2>" -- so any FASTQ whose headers outrun its reads
        wrote past the end of the heap block. Measured before the fix: a 16 KB header
        against 51 bp reads aborted under malloc's heap check on one run and returned
        well-formed output on an identical rerun, which is the signature of an overflow
        that is usually silently corrupting whatever follows it.
        """
        py, cc = _chunk_backends()
        rng = random.Random(4242 + hdrlen)
        frag = draw(rng, 78)                       # 52 bp reads: arena stays at 1024
        s1, s2 = frag[:52], rc(frag[-52:])
        h = b"x" * hdrlen
        b1 = b"@%b/1\tZI:i:7\n%b\n+\n%b\n" % (h, s1, b"I" * len(s1))
        b2 = b"@%b/2\tZI:i:7\n%b\n+\n%b\n" % (h, s2, b"I" * len(s2))
        args = _chunk_args()
        a = py(b1, 0, len(b1), b2, 0, len(b2), *args)
        b = cc(b1, 0, len(b1), b2, 0, len(b2), *args)
        assert a[3][MERGED] == 1, "the fixture stopped merging"
        assert a[0] == b[0], f"blobs differ at header length {hdrlen}"
        assert a[1:4] == b[1:4]
        take1 = min(52, len(frag))                 # R1's share, then R2 supplies the rest
        expected = b"@" + h + b"\tZI:i:7 merged_%d_%d\n" % (take1, len(frag) - take1)
        assert a[0].startswith(expected), a[0][:120]

    def test_a_growing_arena_does_not_corrupt_the_pairs_after_it(self):
        """The arena grows mid-chunk, so everything already merged in that chunk must
        survive it. Feed a short pair, then a long one, then short ones again."""
        py, cc = _chunk_backends()
        rng = random.Random(77)
        recs1, recs2 = [], []
        for i, readlen in enumerate((80, 80, 2200, 80, 90, 3000, 100)):
            frag = draw(rng, readlen * 3 // 2)
            s1, s2 = frag[:readlen], rc(frag[-readlen:])
            recs1.append(b"@r%d/1\n%b\n+\n%b\n" % (i, s1, b"I" * len(s1)))
            recs2.append(b"@r%d/2\n%b\n+\n%b\n" % (i, s2, b"I" * len(s2)))
        b1, b2 = b"".join(recs1), b"".join(recs2)
        p = MergeParams(error_rate="0.01")
        p.ensure(3000)
        args = _chunk_args(p)
        a = py(b1, 0, len(b1), b2, 0, len(b2), *args)
        b = cc(b1, 0, len(b1), b2, 0, len(b2), *args)
        assert a[0] == b[0], "blobs differ once the arena has grown"
        assert a[3] == b[3]
        assert a[3][N_PAIRS] == 7 and a[3][MERGED] == 7, a[3]

    def test_both_backends_stop_at_the_same_pair_for_table_capacity(self):
        py, cc = _chunk_backends()
        rng = random.Random(78)
        recs1, recs2 = [], []
        for i, readlen in enumerate((80, 100, 300, 90)):
            frag = draw(rng, readlen * 3 // 2)
            s1, s2 = frag[:readlen], rc(frag[-readlen:])
            recs1.append(b"@r%d/1\n%b\n+\n%b\n" % (i, s1, b"I" * len(s1)))
            recs2.append(b"@r%d/2\n%b\n+\n%b\n" % (i, s2, b"I" * len(s2)))
        b1, b2 = b"".join(recs1), b"".join(recs2)
        args = _chunk_args(MergeParams(error_rate="0.01"))       # capacity 256
        a = py(b1, 0, len(b1), b2, 0, len(b2), *args)
        b = cc(b1, 0, len(b1), b2, 0, len(b2), *args)
        assert a == b
        assert a[3][N_PAIRS] == 2 and a[8] == 300

    def test_split_records_agrees(self):
        py, cc = _backend_pair("split_records")
        buf = b"".join(b"@r%d\nACGTAC\n+\nIIIIII\n" % i for i in range(7)) + b"@part\nAC"
        for start in (0, 22, 44):
            for n in (0, 1, 3, 7, 99):
                assert py(buf, start, n) == cc(buf, start, n), (start, n)

    def test_whole_pairs_agree_through_process_pair(self):
        """One level up: the same decisions, records and counters from either backend."""
        from zna.merge import overlap
        _backends()          # skips when the extension is not built, like the rest
        rng = random.Random(16)
        pairs = []
        for _ in range(300):
            frag = draw(rng, rng.randrange(40, 320))
            l1, l2 = rng.randrange(30, 151), rng.randrange(30, 151)
            r1 = mutate((frag + ADAPTER1 + draw(rng, 160))[:l1], rng, 0.01)
            r2 = mutate((rc(frag) + ADAPTER2 + draw(rng, 160))[:l2], rng, 0.01)
            pairs.append((r1, qual(r1), r2, qual(r2)))

        def run(name):
            original = overlap.backend_name()
            try:
                overlap.use_backend(name)
                return [process_pair(b"x/1", s1, q1, b"x/2", s2, q2,
                                     MergeParams(error_rate="0.01", min_read_length=40))
                        for s1, q1, s2, q2 in pairs]
            finally:
                overlap.use_backend(original)

        assert run("python") == run("accel")


class TestInputValidationAgrees:
    """Inputs no driver produces, refused the same way by both backends.

    The drivers pass tables from params.py, weights ~2^25-2^29 and an npolicy code from
    MergeParams, so none of this is reachable from `zna merge`; but the backends are
    compared call for call, and past these bounds they disagreed (an unknown npolicy was
    'random' in the kernel and 'keep' in the reference; int64 overflowed where Python
    did not)."""

    S = b"ACGTTGCAAC" * 5

    def _both(self, fn, *args):
        out = []
        for f in _backend_pair(fn):
            try:
                out.append(f(*args))
            except Exception as e:                  # noqa: BLE001 -- compared by type
                out.append(type(e).__name__)
        assert out[0] == out[1], out
        return out[0]

    def test_an_unknown_npolicy_is_refused(self):
        assert self._both("process_pair", b"a", b"ACGTN", b"IIIII", b"b", b"TTTTT",
                          b"IIIII", *_P.kernel_args(), 0, DISAGREE_Q, 3, 42) \
            == "ValueError"
        assert self._both("merge_chunk", b"", 0, 0, b"", 0, 0, *_chunk_args(), -1, 0) \
            == "ValueError"

    @pytest.mark.parametrize("mq,sq,floor", [
        (1 << 58, (1 << 58) + 5, 0),                # n * match_q overflowed int64
        (1 << 31, 1 << 32, 0),
        (_P.match_q, _P.step_q, -(1 << 63) + 1),    # floor - 1 overflowed
        (_P.match_q, _P.step_q, 1 << 62),
    ])
    def test_weights_and_floors_past_int64_safety_are_refused(self, mq, sq, floor):
        assert self._both("scan", self.S, self.S, 50, 50, mq, sq, floor, 0) \
            == "ValueError"

    def test_the_largest_accepted_weights_still_agree(self):
        r = self._both("scan", self.S, self.S, 50, 50, 1 << 30, 1 << 30, -(1 << 61), 0)
        assert r == (0, 50 << 30, 50, 0)

    def test_a_bad_table_is_refused(self):
        from array import array
        t, d = array("q", _P.t_table), _P.dfit_table
        t[99] = -(1 << 63)
        assert self._both("overlap", self.S, self.S, 50, 50, _P.match_q, _P.step_q,
                          t, d, 0) == "ValueError"
        short = array("q", [0])
        assert self._both("merge_chunk", b"", 0, 0, b"", 0, 0, _P.match_q, _P.step_q,
                          short, d, 0, 40, DISAGREE_Q, True, 0) == "ValueError"

    def test_every_trailing_cr_is_stripped(self):
        b1 = b"@r/1\r\r\nACGT\r\r\n+\r\r\nIIII\r\r\n"
        b2 = b"@r/2\r\r\nACGT\r\r\n+\r\r\nIIII\r\r\n"
        r = self._both("merge_chunk", b1, 0, len(b1), b2, 0, len(b2),
                       *_chunk_args(lr=1))
        assert b"\r" not in r[0] and r[3][MAX_READ_LEN] == 4

    def test_a_desync_on_a_non_utf8_name_is_an_input_error_in_both(self):
        from zna.merge.fastqio import InputError
        c1, c2 = b"@r\xff/1\nACGT\n+\nIIII\n", b"@q/2\nACGT\n+\nIIII\n"
        msgs = []
        for f in _backend_pair("merge_chunk"):
            with pytest.raises(InputError) as ei:
                f(c1, 0, len(c1), c2, 0, len(c2), *_chunk_args())
            msgs.append(str(ei.value))
        assert msgs[0] == msgs[1] == "R1/R2 out of sync at pair 1: 'r\xff' != 'q'"

    def test_a_refused_table_is_released(self):
        """The compiled backend borrows a table through the buffer protocol. A table it
        refuses AFTER borrowing (wrong item type) must still be released: it used to
        leak the export, and the array could never be resized again."""
        from array import array
        _py, accel = _backend_pair("overlap")
        wrong = array("i", range(600))
        before = sys.getrefcount(wrong)
        for _ in range(100):
            with pytest.raises(TypeError):
                accel(self.S, self.S, 50, 50, _P.match_q, _P.step_q, wrong,
                      _P.dfit_table, 0)
        assert sys.getrefcount(wrong) == before
        wrong.append(1)                             # BufferError if still exported


# --------------------------------------------------------------------------- #
# 5. the fixed-point scale
# --------------------------------------------------------------------------- #

class TestFixedPointScale:
    """The score is computed in integers so the argmax is reproducible everywhere.

    Two things have to hold for that to be worth anything: the integers must be the
    same integers on every platform, and quantising must not move a decision.
    """

    def test_the_scale(self):
        assert SCALE == 1 << 24

    def test_step_is_the_sum_of_the_two_quantised_weights(self):
        """`score = n*match - d*step` and `score = (n-d)*match - d*mismatch` must agree,
        which they only do if `step` is quantised as the sum rather than separately."""
        from zna.merge.params import log2_exact, P_NULL
        for e in (Fraction(1, 100), Fraction(874, 10 ** 5), Fraction(3, 100)):
            match_q, step_q = weights_q(e)
            assert step_q == match_q + to_q(log2_exact((1 - P_NULL) / e))
        for n, d in ((40, 0), (40, 1), (150, 7), (19, 2)):
            assert n * _P.match_q - d * _P.step_q == score_of(n - d, d)

    @pytest.mark.parametrize("e", ["0.01", "0.0087", "0.001"])
    @pytest.mark.parametrize("n_shifts", [99, 299, 2047])
    def test_quantisation_flips_no_merge_decision_over_the_reachable_domain(
            self, e, n_shifts):
        """The enumeration params.py's docstring reports, as a test.

        A decision flips only where the exact score sits within the quantisation error
        of the pair's floor ``log2(N / alpha)``, and over integer (n, d) that is
        exhaustively checkable. Checked here to an overlap of 4,000 bases for nine
        (e, N) settings; the smallest disagreement over the 30 settings the docstring
        reports is at n = 10,951.
        """
        from decimal import Decimal
        from zna.merge.params import _CTX, P_NULL, log2_exact
        ef = Fraction(e)
        lm, lmm = log2_exact((1 - ef) / P_NULL), log2_exact((1 - P_NULL) / ef)
        lt = log2_exact(Fraction(n_shifts) / Fraction(1, 10 ** 6))
        mq, sq = weights_q(ef)
        tq = threshold_q(n_shifts, "1e-6")
        fm, fs, ft = float(lm), float(lm + lmm), float(lt)
        for n in range(1, 4001):
            # only d values that put the score anywhere near the floor matter
            dc = (n * fm - ft) / fs
            for d in range(max(0, int(dc) - 1), min(n, int(dc) + 2) + 1):
                exact_score = _CTX.subtract(_CTX.multiply(Decimal(n - d), lm),
                                            _CTX.multiply(Decimal(d), lmm))
                assert (exact_score >= lt) == (n * mq - d * sq >= tq), (n, d)


# --------------------------------------------------------------------------- #
# 6. the argmax total order
# --------------------------------------------------------------------------- #

class TestArgmaxTotalOrder:
    """`maximise score, then minimise s` — a specification, not an iteration artifact.

    This is what lets a rewritten kernel be tested for byte-exact equality instead of
    the weaker "returned *an* argmax". Random sequence essentially never ties, so the
    ties here are built deliberately.
    """

    def _check(self, s1, s2rc, label, adapter_trimmed=False):
        from zna.merge.overlap import _backend
        got = _backend.active().scan(s1, s2rc, len(s1), len(s2rc), _P.match_q,
                                     _P.step_q, LOW_FLOOR, int(adapter_trimmed))
        allsc = exhaustive_scan(s1, s2rc, LOW_FLOOR, adapter_trimmed)
        if not allsc:
            assert got[2] == 0, label
            return 0
        want, n_ties = argmax_by_rule(allsc)
        assert got == (want[0], want[3], want[1], want[2]), (label, got, want, n_ties)
        return n_ties - 1

    @pytest.mark.parametrize("adapter_trimmed", [False, True])
    def test_matches_an_unpruned_scan_on_real_reads(self, any_backend, adapter_trimmed):
        rng = random.Random(9)
        for i in range(300):
            frag = draw(rng, rng.randrange(40, 110))
            l1, l2 = rng.randrange(20, 60), rng.randrange(20, 60)
            r1 = (frag + ADAPTER1 + draw(rng, 60))[:l1]
            r2 = (rc(frag) + ADAPTER2 + draw(rng, 60))[:l2]
            self._check(r1, rc(r2), f"real{i}", adapter_trimmed)

    @pytest.mark.parametrize("adapter_trimmed", [False, True])
    def test_ties_are_broken_towards_the_smallest_shift(self, any_backend,
                                                        adapter_trimmed):
        """Periodic content on unequal-length mates makes the plateau tie exactly.

        Without this the tie-break is untested: an earlier sweep over 7,000 random and
        adversarial pairs produced *zero* ties and proved nothing about it.
        """
        tied = sum(self._check(s1, s2rc, label, adapter_trimmed)
                   for s1, s2rc, label in tie_fixtures())
        if adapter_trimmed:
            # Under the declaration every overlap length has exactly ONE eligible shift
            # (the plateau shrinks to its last shift, the read-through flank is gone),
            # and ties across lengths are unreachable -- so the argmax is unique.
            assert tied == 0
        else:
            assert tied >= 100, f"only {tied} ties exercised; the tie-break is untested"

    def test_a_tie_across_different_overlap_lengths_is_unreachable(self):
        """Why the rule needs no `n` key.

        Two shifts tie iff `dn * match_q == dd * step_q`, whose minimal solution is
        `dn = step_q / gcd(match_q, step_q)`. For every error rate checked that is far
        larger than any conceivable read -- ties can only ever occur at equal `n`.
        """
        for err in ("0.00016", "0.001", "0.005", "0.0087", "0.01", "0.02", "0.05",
                    "0.1", "0.3"):
            p = MergeParams(error_rate=err)
            dn = p.step_q // math.gcd(p.match_q, p.step_q)
            assert dn > 10_000_000, (err, dn)


# --------------------------------------------------------------------------- #
# 7. find_overlap, and the contract range
# --------------------------------------------------------------------------- #

class TestFindOverlap:
    def test_forward_normal_overlap(self):
        frag = rand_seq(40, 1)          # insert 40, read 30 -> overlap 20 at offset 10
        (_, r1, _), (_, r2, _) = make_pair(frag, 30)
        o = find_overlap(r1, rc(r2), P)
        assert (o.verdict, o.shift, o.overlap_len, o.mismatches) == (V_MERGE, 10, 20, 0)
        assert o.score_q == score_of(20) and o.fragment_length == 40

    def test_full_overlap(self):
        frag = rand_seq(30, 2)          # insert == read len -> full overlap, shift 0
        (_, r1, _), (_, r2, _) = make_pair(frag, 30)
        o = find_overlap(r1, rc(r2), P)
        assert (o.verdict, o.shift, o.overlap_len, o.mismatches) == (V_MERGE, 0, 30, 0)

    def test_no_overlap(self):
        o = find_overlap(rand_seq(50, 3), rc(rand_seq(50, 4)), P)
        assert o == (V_NONE, 0, 0, 0, 0, 0, 0)

    def test_mismatch_within_budget_is_accepted(self):
        frag = rand_seq(40, 5)
        (_, r1, _), (_, r2, _) = make_pair(frag, 30)
        r2 = bytearray(r2)
        # The overlap is R2's 3' end (its last bases map to the start of revcomp(R2)).
        r2[-1] = flip(r2[-1])                                    # 1 error in overlap
        o = find_overlap(r1, rc(bytes(r2)), P)
        assert (o.verdict, o.overlap_len, o.mismatches) == (V_MERGE, 20, 1)

    def test_noise_below_the_floor_is_no_overlap(self):
        frag = rand_seq(40, 6)
        (_, r1, _), (_, r2, _) = make_pair(frag, 30)
        r2 = bytearray(r2)
        for i in range(1, 12):                 # 11 errors in a 20 bp overlap: 9*1.99
            r2[-i] = flip(r2[-i])              # - 11*6.23 < 0
        assert find_overlap(r1, rc(bytes(r2)), P).verdict == V_NONE

    def test_read_through_is_a_negative_shift(self):
        # insert (20) shorter than read (30): both reads run past into adapter.
        insert = rand_seq(20, 8)
        r1 = (insert + ADAPTER1)[:30]
        r2 = (rc(insert) + ADAPTER2)[:30]
        o = find_overlap(r1, rc(r2), P)
        assert (o.verdict, o.shift, o.overlap_len, o.fragment_length) == \
            (V_MERGE, -10, 20, 20)

    def test_argmax_beats_a_spurious_short_hit(self):
        """A real 40 bp overlap wins outright over any chance 4-mer earlier in the scan.

        First-accept could be captured by the short hit; argmax cannot.
        """
        rng = random.Random(99)
        for _ in range(200):
            frag = draw(rng, 160)                       # insert 160, L 100 -> overlap 40
            (_, r1, _), (_, r2, _) = make_pair(frag, 100)
            o = find_overlap(r1, rc(r2), P)
            assert (o.verdict, o.shift, o.overlap_len) == (V_MERGE, 60, 40)
            assert o.score_q == score_of(40)

    def test_block_loop_sees_every_position(self, any_backend):
        """Sweep a single mismatch across every position of a 40 bp overlap.

        ``_shift_score`` accumulates mismatches branchlessly in blocks of 8 and only
        tests the bail once per block. A stride bug there (``k += 7``, an off-by-one in
        the 8-term unroll, a wrong ``lim``) silently mis-scores any overlap containing an
        error — it changes 6.34% of scores and 0.26% of merge/trim/keep decisions on real
        data, and every other test in this suite passes with it in place, because they
        use clean overlaps or put the mismatch in one fixed position.
        """
        rng = random.Random(24680)
        frag = draw(rng, 40)                       # insert == readlen -> full overlap
        r2rc = rc(rc(frag))
        for i in range(40):
            r1 = bytearray(frag)
            r1[i] = flip(r1[i])
            o = find_overlap(bytes(r1), r2rc, P)
            assert (o.verdict, o.shift, o.overlap_len, o.mismatches) == \
                (V_MERGE, 0, 40, 1), i
            assert o.score_q == score_of(39, 1), i

    def test_unequal_read_lengths_need_no_special_case(self):
        """s is defined by the offset; the compared region is just the intersection."""
        rng = random.Random(7)
        frag = draw(rng, 180)
        r1 = frag[:120]                     # R1 quality-trimmed to 120
        r2 = rc(frag[-90:])                 # R2 quality-trimmed to 90
        o = find_overlap(r1, rc(r2), P)
        assert (o.verdict, o.shift, o.overlap_len, o.mismatches) == (V_MERGE, 90, 30, 0)

    def test_the_diagnostic_scan_is_the_bare_argmax(self):
        """`scan_unrestricted` ignores the contract and the gate, and says so in its
        name: it is what the read-through check counts, never what a pair is built
        from."""
        insert = rand_seq(40, 9)
        r1 = (insert + ADAPTER1 + rand_seq(80, 10))[:100]
        r2 = (rc(insert) + ADAPTER2 + rand_seq(80, 11))[:100]
        p = MergeParams(error_rate="0.01", adapter_trimmed=True)
        a = scan_unrestricted(r1, rc(r2), p)
        assert (a.shift, a.overlap_len) == (-60, 40)
        assert find_overlap(r1, rc(r2), p).verdict == V_NONE


class TestContractRange:
    """``--adapter-trimmed``: only ``L = s + len2 >= max(len1, len2)`` is eligible.

    That is not just ``s >= 0``. With ``len1 > len2`` the plateau shifts ``0 <= s <
    len1 - len2`` put R1 past the fragment's end -- as impossible under the declaration
    as a read-through -- so eligibility starts at ``len1 - len2``. With ``len1 < len2``
    it starts at 0, and the plateau's negative shifts are read-through for R2.
    """

    @staticmethod
    def _eligible(len1, len2):
        return [s for s in range(-(len2 - 1), len1) if s + len2 >= max(len1, len2)]

    @pytest.mark.parametrize("len1,len2", [(100, 60), (60, 100), (80, 80), (31, 7),
                                           (7, 31)])
    def test_the_eligible_range(self, len1, len2):
        assert min(self._eligible(len1, len2)) == max(0, len1 - len2)

    def test_len1_longer_the_plateau_start_is_ineligible(self):
        """R1 = 100, R2 = 60, and the true fragment is 80: R1 runs 20 bases past it.

        Unrestricted, the scan finds s = 20 (L = 80, a plateau shift with 0 <= s <
        len1 - len2); declared trimmed, that is impossible, and nothing else aligns."""
        rng = random.Random(31)
        frag = draw(rng, 80)
        r1 = (frag + ADAPTER1)[:100]
        r2 = rc(frag)[:60]                                  # R2 = frag[20:80], revcomp
        free = find_overlap(r1, rc(r2), P)
        assert (free.verdict, free.shift, free.fragment_length) == (V_MERGE, 20, 80)
        assert 0 <= free.shift < len(r1) - len(r2)
        declared = MergeParams(error_rate="0.01", min_read_length=1,
                               adapter_trimmed=True)
        assert find_overlap(r1, rc(r2), declared).verdict == V_NONE

    def test_len1_shorter_a_plateau_read_through_is_ineligible(self):
        """R1 = 60, R2 = 100, fragment 80: s = -20 lies on the plateau (n = 60) but puts
        R2 past the fragment."""
        rng = random.Random(32)
        frag = draw(rng, 80)
        r1 = frag[:60]
        r2 = (rc(frag) + ADAPTER2)[:100]
        free = find_overlap(r1, rc(r2), P)
        assert (free.verdict, free.shift, free.overlap_len) == (V_MERGE, -20, 60)
        declared = MergeParams(error_rate="0.01", adapter_trimmed=True)
        assert find_overlap(r1, rc(r2), declared).verdict == V_NONE

    @pytest.mark.parametrize("len1,len2", [(100, 60), (60, 100), (90, 90)])
    def test_eligible_geometries_still_merge(self, len1, len2):
        """Every L >= max(len1, len2) with enough overlap merges under the declaration,
        at exactly the shift it merges at without it."""
        rng = random.Random(33 + len1)
        declared = MergeParams(error_rate="0.01", min_read_length=1,
                               adapter_trimmed=True)
        for L in range(max(len1, len2), len1 + len2 - 16):
            frag = draw(rng, L)
            r1, r2 = frag[:len1], rc(frag[L - len2:])
            a, b = find_overlap(r1, rc(r2), P), find_overlap(r1, rc(r2), declared)
            assert a == b and a.verdict == V_MERGE and a.fragment_length == L, (L, a, b)

    def test_the_contract_is_the_argmax_over_eligible_shifts_not_a_veto(self):
        """A read-through winner does not block a lower-scoring eligible alignment: the
        declared scan is the argmax over the eligible set, which here is the true L."""
        rng = random.Random(34)
        frag = draw(rng, 250)
        r1, r2 = bytearray(frag[:150]), bytearray(rc(frag[100:]))
        # plant a strong read-through: R1's first 60 bases equal R2rc's last 60
        r2rc = bytearray(rc(bytes(r2)))
        r2rc[-60:] = r1[:60]
        r2 = rc(bytes(r2rc))
        free = find_overlap(bytes(r1), rc(r2), P)
        declared = find_overlap(bytes(r1), rc(r2),
                                MergeParams(error_rate="0.01", adapter_trimmed=True))
        assert free.shift < 0 and free.verdict == V_MERGE
        assert (declared.verdict, declared.fragment_length) == (V_MERGE, 250)


# --------------------------------------------------------------------------- #
# 8. the plausibility gate
# --------------------------------------------------------------------------- #

class TestPlausibilityGate:
    """``d_inf > dfit[n]`` -> the pair has no overlap. Abstain; never re-place."""

    @staticmethod
    def _repeat_beats_truth(seed=41):
        """The review's C01, rebuilt: a 262-base fragment read 2x150 (true overlap 38)
        through a diverged tandem repeat of period 84, so R1[28:150] and revcomp(R2)'s
        first 122 bases are two copies of it. Every 7th base differs between copies: 18
        mismatches in 122, 94.4 bits at 1% -- beating the true, perfect 38-base overlap
        (75.4 bits) -- and 18 > dfit[122] = 9."""
        rng = random.Random(seed)
        u = draw(rng, 84)
        u1 = bytearray(u)
        for i in range(0, 84, 7):
            u1[i] = flip(u1[i])
        u2 = bytearray(u1[:38])
        for i in range(0, 38, 7):
            u2[i] = flip(u2[i])
        frag = draw(rng, 28) + u + bytes(u1) + bytes(u2) + draw(rng, 28)
        assert len(frag) == 262
        return frag

    def test_a_divergent_repeat_is_refused_not_merged(self):
        frag = self._repeat_beats_truth()
        L = len(frag)
        r1, r2 = frag[:150], rc(frag[L - 150:])
        o = find_overlap(r1, rc(r2), P)
        assert o.verdict == V_IMPLAUSIBLE and o.fragment_length != L
        assert o.informative_mismatches > P.dfit(o.overlap_len)
        # ...and nothing is searched for in its place: process_pair keeps both mates,
        # untouched, and reports the refusal.
        counters = [0, 0]
        res = process_pair(b"c/1", r1, qual(r1), b"c/2", r2, qual(r2), P, counters)
        assert res.outcome == PairOutcome.KEPT and res.implausible
        assert [r[1] for r in res.records] == [r1, r2]
        assert (res.olen, res.diff, res.score) == (0, 0, 0)
        assert counters == [0, 1]

    def test_the_gate_is_at_dfit_exactly(self):
        """A 64-base full overlap with k scattered mismatches: dfit[64] = 7 at 1%."""
        rng = random.Random(42)
        frag = draw(rng, 64)
        for k, want in ((7, V_MERGE), (8, V_IMPLAUSIBLE)):
            r1 = bytearray(frag)
            for i in range(k):
                r1[3 + 8 * i] = flip(r1[3 + 8 * i])
            o = find_overlap(bytes(r1), frag, P)
            assert (o.overlap_len, o.mismatches, o.verdict) == (64, k, want), k

    def test_an_n_run_is_not_evidence_of_disagreement(self):
        """Ten N in a 100-base true overlap: raw d = 10 > dfit[100] = 8, informative
        d = 0. Refusing it cost ~1,500 correct merges per 200k pairs on reads with N
        runs (the robustness review's one correctness regression)."""
        frag = rand_seq(200, 43)
        r1 = bytearray(frag[:150])
        r1[80:90] = b"N" * 10
        r2 = rc(frag[50:])
        o = find_overlap(bytes(r1), rc(r2), P)
        assert (o.overlap_len, o.mismatches, o.informative_mismatches) == (100, 10, 0)
        assert o.mismatches > P.dfit(100) and o.verdict == V_MERGE
        res = process_pair(b"n/1", bytes(r1), qual(r1), b"n/2", r2, qual(r2), P)
        assert res.outcome == PairOutcome.MERGED and res.records[0][1] == frag

    def test_the_same_run_in_r2_is_discounted_too(self):
        frag = rand_seq(200, 44)
        r1 = frag[:150]
        r2 = bytearray(rc(frag[50:]))
        r2[70:80] = b"N" * 10                       # R2's own frame, inside the overlap
        o = find_overlap(r1, rc(bytes(r2)), P)
        assert (o.mismatches, o.informative_mismatches, o.verdict) == (10, 0, V_MERGE)

    def test_n_against_n_is_a_match_and_not_discounted(self):
        """Both mates N at the same positions: the scan scores N==N as agreement, so
        there is nothing to discount; the count is of positions where EXACTLY one is N."""
        frag = bytearray(rand_seq(200, 45))
        frag[100:110] = b"N" * 10
        frag = bytes(frag)
        o = find_overlap(frag[:150], frag[50:], P)
        assert (o.mismatches, o.informative_mismatches, o.verdict) == (0, 0, V_MERGE)

    def test_n_is_discounted_only_where_it_meets_a_call(self):
        rng = random.Random(46)
        frag = draw(rng, 120)
        s1 = bytearray(frag)
        s1[10:14] = b"NNNN"                         # 4 one-sided N
        for i in (40, 50, 60, 70, 80, 90, 100, 110, 115):
            s1[i] = flip(s1[i])                     # 9 real mismatches > dfit[120] = 9?
        o = find_overlap(bytes(s1), frag, P)
        assert (o.mismatches, o.informative_mismatches) == (13, 9)
        assert one_n(bytes(s1), frag, 0, 120) == 4
        assert o.verdict == (V_MERGE if 9 <= P.dfit(120) else V_IMPLAUSIBLE)

    def test_a_refused_pair_gets_no_consensus(self):
        """Nothing is rewritten on a pair whose alignment was refused: its apparent
        disagreements are the repeat's, not sequencing errors."""
        frag = self._repeat_beats_truth(47)
        L = len(frag)
        r1, r2 = frag[:150], rc(frag[L - 150:])
        q1, q2 = b"#" * 150, b"~" * 150              # R2 would win every contest
        res = process_pair(b"c/1", r1, q1, b"c/2", r2, q2, P)
        assert res.implausible
        assert res.records[0][1:] == (r1, q1) and res.records[1][1:] == (r2, q2)

    def test_the_error_rate_moves_the_gate_and_the_argmax(self):
        """The review's observation, in both directions: at the library's own
        near-zero rate a mismatch costs ~12 bits, so the divergent repeat no longer
        even outscores the true perfect overlap; at 1% it wins the argmax and the gate
        refuses it."""
        frag = self._repeat_beats_truth(48)
        L = len(frag)
        r1, r2 = frag[:150], rc(frag[L - 150:])
        low = MergeParams(error_rate="0.00016", min_read_length=1)
        o = find_overlap(r1, rc(r2), low)
        assert o.verdict == V_MERGE and o.fragment_length == L and o.mismatches == 0
        assert find_overlap(r1, rc(r2), P).verdict == V_IMPLAUSIBLE


class TestDegenerateInputs:
    """Zero-length and 1-base reads survived every zna audit; keep it that way."""

    @pytest.mark.parametrize("s1,s2", [
        (b"", b"ACGTACGTAC"), (b"ACGTACGTAC", b""), (b"", b""), (b"A", b"T"),
    ])
    def test_no_overlap_and_no_crash(self, s1, s2):
        assert find_overlap(s1, rc(s2), P) == (V_NONE, 0, 0, 0, 0, 0, 0)

    def test_empty_pair_is_kept_not_merged(self):
        res = process_pair(b"z/1", b"", b"", b"z/2", b"", b"",
                           MergeParams(min_read_length=1))
        assert res.outcome == PairOutcome.KEPT and res.score == 0 and res.records == []

    def test_n_is_scored_as_an_ordinary_base(self):
        """N vs N counts as a 2-bit match in the SCORE -- inherited from the old kernel,
        and left alone because the scan is byte comparison. Only the plausibility gate
        and the detected-overlap rate discount a one-sided N."""
        o = find_overlap(b"N" * 20 + b"ACGT" * 20, rc(b"N" * 20 + b"ACGT" * 20), P)
        assert o.mismatches == 0 and o.score_q == score_of(o.overlap_len)


# --------------------------------------------------------------------------- #
# 9. detection, read-through, the boundary invariant
# --------------------------------------------------------------------------- #

class TestSpuriousDetection:
    """alpha bounds chance merges of unrelated sequence: none may appear."""

    def test_unrelated_pairs_do_not_merge(self):
        n = 20000 if _fast_backend() else 2000   # the reference kernel is ~50x slower
        rng = random.Random(20260811)
        merged = sum(find_overlap(draw(rng, 150), rc(draw(rng, 150)), P).verdict
                     != V_NONE for _ in range(n))
        # 1e-6 per pair: none in 20k. (0.5.3's trim band detected 0.2% of these.)
        assert merged == 0, f"{merged} spurious alignments reached T"

    def test_polya_does_not_merge(self):
        """Low-complexity tails must not produce a confident merge on their own: the
        score has no low-complexity correction, and a 20-base polyA scores 39.7 bits
        clean -- above T -- only where it is the whole overlap."""
        rng = random.Random(5150)
        merged = 0
        for _ in range(200):
            core1, core2 = draw(rng, 130), draw(rng, 130)
            r1 = core1 + b"A" * 20                    # polyA tail on both mates
            r2 = core2 + b"A" * 20
            merged += find_overlap(r1, rc(r2), P).verdict == V_MERGE
        assert merged == 0


class TestDetection:
    def test_known_overlaps_recover_the_true_shift(self):
        """Overlaps 4..40 at 0.5% per-base error, 2x100.

        Sensitivity is set by the arithmetic, not by tuning: at 2x100 the floor is
        27.57 bits, so 14 clean bases (27.80) reach it and 13 (25.81) cannot; one
        mismatch costs 8.2 bits, so it takes 19 bases (29.51) for an error to fit.
        """
        L = 100
        rng = random.Random(31337)
        wrong = trials = 0
        detected_by_olen = {}
        for olen in range(4, 41):
            hits = 0
            reps = 100
            for _ in range(reps):
                frag = draw(rng, 2 * L - olen)
                (_, r1, _), (_, r2, _) = make_pair(frag, L)
                r1 = mutate(r1, rng, 0.005)
                r2 = mutate(r2, rng, 0.005)
                o = find_overlap(r1, rc(r2), P)
                trials += 1
                if o.verdict != V_MERGE:
                    continue
                if o.shift == L - olen:
                    hits += 1
                else:
                    wrong += 1                      # a chance shift outscored the truth
            detected_by_olen[olen] = hits / reps
        assert wrong / trials < 0.002, f"{wrong}/{trials} pairs merged at a wrong shift"
        assert min_matches(t_q(L, L), 0) == 14
        for olen in range(4, 14):
            assert detected_by_olen[olen] == 0.0, olen
        assert min_matches(t_q(L, L), 1) == 18   # 18 matches + 1 mismatch = 19 bases
        for olen in range(14, 19):               # zero-mismatch band: P(no error)
            assert detected_by_olen[olen] >= 0.80, (olen, detected_by_olen[olen])
        for olen in range(19, 41):               # one mismatch now fits
            assert detected_by_olen[olen] >= 0.95, (olen, detected_by_olen[olen])

    @pytest.mark.parametrize("read_len,shortest", [(50, 14), (100, 14), (150, 15),
                                                   (300, 15)])
    def test_the_shortest_mergeable_clean_overlap(self, read_len, shortest):
        """ceil(T(N) / match_bits): 14-15 bases from 2x50 to 2x300 (plan §2)."""
        assert min_matches(t_q(read_len, read_len), 0) == shortest
        rng = random.Random(6161 + read_len)
        for olen, want in ((shortest - 1, False), (shortest, True), (30, True)):
            frag = draw(rng, 2 * read_len - olen)
            (h1, s1, q1), (h2, s2, q2) = make_pair(frag, read_len)
            res = process_pair(h1, s1, q1, h2, s2, q2, P)
            assert (res.outcome == PairOutcome.MERGED) is want, (olen, res)
            if want:
                assert res.score == score_of(olen)


class TestReadThrough:
    @pytest.mark.parametrize("insert", list(range(40, 100, 7)))
    def test_insert_shorter_than_read_merges_to_the_fragment(self, insert):
        """Without the declaration read-through IS adapter removal: the merged record
        is [0, L) and excludes the adapter."""
        rng = random.Random(1000 + insert)
        frag = draw(rng, insert)
        (h1, s1, q1), (h2, s2, q2) = cycle_pair(frag, 100, rng)
        p = MergeParams(error_rate="0.01", min_read_length=40)
        res = process_pair(h1, s1, q1, h2, s2, q2, p)
        assert res.outcome == PairOutcome.MERGED and res.n_dropped == 0
        assert [r[1] for r in res.records] == [frag]      # adapter and filler both gone
        assert len(res.records[0][2]) == len(frag)

    @pytest.mark.parametrize("insert", list(range(40, 100, 13)))
    def test_the_declaration_keeps_the_same_pair_whole(self, insert):
        """Declared trimmed on reads that are NOT: the read-through is impossible, so the
        pair is kept, adapter and all. This is why the declaration is checked."""
        rng = random.Random(1000 + insert)
        frag = draw(rng, insert)
        (h1, s1, q1), (h2, s2, q2) = cycle_pair(frag, 100, rng)
        p = MergeParams(error_rate="0.01", min_read_length=40, adapter_trimmed=True)
        res = process_pair(h1, s1, q1, h2, s2, q2, p)
        assert res.outcome == PairOutcome.KEPT
        assert [r[1] for r in res.records] == [s1, s2]


class TestBoundaryInvariant:
    @pytest.mark.parametrize("insert", list(range(45, 320, 9)))
    def test_base_zero_is_always_a_true_fragment_boundary(self, insert):
        """Nothing is ever removed from a read, and a merged read is the fragment
        exactly (both of its edges are true boundaries). With the trim band gone an
        unmerged mate is emitted whole."""
        rng = random.Random(2000 + insert)
        frag = draw(rng, insert)
        (h1, s1, q1), (h2, s2, q2) = cycle_pair(frag, 150, rng)
        p = MergeParams(error_rate="0.01", min_read_length=40)
        res = process_pair(h1, s1, q1, h2, s2, q2, p)
        if res.outcome == PairOutcome.MERGED:
            assert [r[1] for r in res.records] == [frag]
        else:
            assert [r[1] for r in res.records] == [s1, s2]
            assert s1[0] == frag[0] and s2[0] == rc(frag)[0]

    def test_merged_read_is_in_r1s_frame(self):
        """zna's single/merged normalization assumes merged reads are R1-framed."""
        rng = random.Random(808)
        frag = draw(rng, 160)
        (h1, s1, q1), (h2, s2, q2) = cycle_pair(frag, 100, rng)
        res = process_pair(h1, s1, q1, h2, s2, q2, P)
        assert res.outcome == PairOutcome.MERGED
        assert res.records[0][1].startswith(s1)          # R1 first, not revcomp(R2)


class TestParityWithLegacyRule:
    def test_clean_long_overlaps_merge_at_the_same_shift(self):
        """Overlap >= 30 and clean: old and new must agree on the shift."""
        L = 100
        rng = random.Random(4711)
        for olen in range(30, 81, 5):
            for _ in range(20):
                frag = draw(rng, 2 * L - olen)
                (_, r1, _), (_, r2, _) = make_pair(frag, L)
                r2rc = rc(r2)
                old_dir, old_shift, _o, _d = legacy_scan(r1, r2rc, len(r1), len(r2rc),
                                                         3, 3, 0.20)
                o = find_overlap(r1, r2rc, P)
                assert (old_dir, old_shift) == (1, L - olen)
                assert (o.verdict, o.shift) == (V_MERGE, L - olen)

    def test_legacy_rule_accepts_chance_four_mers_and_the_new_one_does_not(self):
        """Pins WHY the rules differ: same unrelated pairs, 5.17% vs none."""
        n = 5000 if _fast_backend() else 400
        rng = random.Random(2468)
        old_hits = new_hits = 0
        for _ in range(n):
            r1 = draw(rng, 150)
            r2rc = rc(draw(rng, 150))
            old_hits += legacy_scan(r1, r2rc, 150, 150, 3, 3, 0.20)[0] != 0
            new_hits += find_overlap(r1, r2rc, P).verdict != V_NONE
        assert old_hits / n > 0.03            # the defect the LR redesign removed
        assert new_hits == 0


# --------------------------------------------------------------------------- #
# 10. process_pair
# --------------------------------------------------------------------------- #

class TestProcessPair:
    def test_full_overlap_merges_to_r1(self):
        frag = rand_seq(30, 10)
        (h1, s1, q1), (h2, s2, q2) = make_pair(frag, 30)
        res = process_pair(h1, s1, q1, h2, s2, q2, P)
        assert res.outcome == PairOutcome.MERGED and res.n_dropped == 0
        assert [r[1] for r in res.records] == [frag]        # merged == fragment
        assert res.records[0][0] == b"frag merged_30_0"     # /1 stripped + fastp token
        assert (res.shift, res.olen, res.implausible) == (0, 30, False)

    def test_partial_overlap_reconstructs_fragment(self):
        frag = rand_seq(40, 11)
        (h1, s1, q1), (h2, s2, q2) = make_pair(frag, 30)
        res = process_pair(h1, s1, q1, h2, s2, q2, P)
        assert res.outcome == PairOutcome.MERGED
        assert res.records[0][1] == frag and len(res.records[0][2]) == len(frag)
        assert res.shift + len(s2) == 40

    def test_r1_wins_keeps_r1_base_on_r2_error(self):
        frag = rand_seq(40, 12)
        (h1, s1, q1), (h2, s2, q2) = make_pair(frag, 30)
        s2 = bytearray(s2)
        s2[-1] = flip(s2[-1])                           # error in R2's overlap end
        res = process_pair(h1, s1, q1, h2, bytes(s2), q2, P)
        assert res.outcome == PairOutcome.MERGED
        assert res.records[0][1] == frag                # R1's (correct) base wins

    def _mismatch_pair(self, seed, q_r1, q_r2, pos=15):
        """Insert 40, read 30 -> overlap R1[10:30]. R1 carries an error at `pos` called
        at `q_r1`; R2 has the correct base at `q_r2`. Returns (args, fragment)."""
        frag = rand_seq(40, seed)
        r1 = bytearray(frag[:30]); q1 = bytearray(bytes([q_r1 + 33]) * 30)
        r1[pos] = flip(r1[pos])
        r2 = rc(frag[10:40]); q2 = bytes([q_r2 + 33]) * 30
        return (bytes(r1), bytes(q1), r2, q2), frag

    def test_an_n_is_rescued_regardless_of_its_quality(self):
        """An N carries no base information, so a real call beats it — whatever the
        quality scores say. Sweep the N's quality across the whole range."""
        frag = rand_seq(200, 7)
        for qN in (0, 2, 20, 40, 60, 93):
            r1 = bytearray(frag[:150])
            pos, true_base = 120, frag[120:121]
            r1[pos] = ord("N")
            q1 = bytearray(bytes([30 + 33]) * 150)
            q1[pos] = qN + 33
            r2 = rc(frag[50:200])
            res = process_pair(b"n/1", bytes(r1), bytes(q1), b"n/2", r2,
                               bytes([30 + 33]) * 150, P)
            assert res.outcome == PairOutcome.MERGED
            assert res.records[0][1][pos:pos + 1] == true_base, (
                f"an N at Q{qN} was not rescued from the mate")
            assert b"N" not in res.records[0][1]
            # the merged record carries the rescue in its provenance
            assert res.records[0][0] == b"n ZN:i:2 rescued_1 merged_150_50"

    def test_both_mates_are_n_so_there_is_nothing_to_rescue_from(self):
        """With no real call opposite it, an N cannot be rescued — so trim3 removes it."""
        frag = rand_seq(200, 8)
        r1 = bytearray(frag[:150]); r2 = bytearray(rc(frag[50:200]))
        r1[120] = ord("N")
        r2[len(r2) - 1 - (120 - 50)] = ord("N")
        args = (b"n/1", bytes(r1), b"I" * 150, b"n/2", bytes(r2), b"I" * 150)

        res = process_pair(*args, MergeParams(error_rate="0.01", min_read_length=1,
                                              npolicy="keep"))
        assert res.records[0][1][120:121] == b"N", "an N with no real call must stay"

        res = process_pair(*args, P)
        assert all(b"N" not in r[1] for r in res.records)

    def test_trim3_cuts_at_the_first_surviving_n_and_keeps_the_5_anchor(self):
        """3' only. Base 0 is a fragment terminus and must survive any trim."""
        frag = rand_seq(400, 13)
        r1 = bytearray(frag[:150]); r1[37] = ord("N")
        r2 = rc(frag[250:400])                      # no overlap: the pair is kept
        res = process_pair(b"t/1", bytes(r1), b"I" * 150, b"t/2", r2, b"I" * 150, P)
        assert res.outcome == PairOutcome.KEPT
        assert res.records[0][1] == frag[:37], "trim3 must cut exactly at the first N"
        assert res.records[1][1] == r2, "the mate with no N must be untouched"
        assert res.records[0][0] == b"t/1 ZN:i:4 trim3_113"

    @pytest.mark.parametrize("npos,merges", [(40, False), (80, True), (120, True)])
    def test_a_trimmed_pair_can_still_merge_and_is_still_the_whole_fragment(
            self, npos, merges):
        """After trim3, merge on GEOMETRY, reusing the original evidence.

        trim3 removes interior bases and leaves both 5' anchors, so R1' covers [0, k1)
        and R2' covers [L-k2, L). The pair still tiles the fragment iff k1 + k2 >= L —
        and when it does the reconstruction is the fragment exactly, N-free. No re-scan:
        trim3 cuts 3' ends, which is where a normal overlap lives.
        """
        L, RL = 220, 150
        frag = rand_seq(L, 17)
        r1 = bytearray((frag + ADAPTER1 + rand_seq(RL, 18))[:RL])
        r2 = (rc(frag) + ADAPTER2 + rand_seq(RL, 19))[:RL]
        r1[npos] = ord("N")
        res = process_pair(b"m/1", bytes(r1), b"I" * RL, b"m/2", r2, b"I" * RL, P)
        assert all(b"N" not in r[1] for r in res.records)
        if merges:
            assert res.outcome == PairOutcome.MERGED
            assert res.records[0][1] == frag
        else:
            assert res.outcome == PairOutcome.KEPT, "coverage failed: keep the pair"
            assert len(res.records[0][1]) == npos and res.records[1][1] == r2

    @staticmethod
    def _merge_verdict_then_trim3_breaks_tiling(r1_edit):
        """2x150, fragment 200 (overlap R1[50:150]). R2 has an N at its own index 45 --
        fragment position 154, outside the overlap -- so trim3 leaves it 45 bases and
        150 + 45 < 200: the verdict is MERGE but the mates no longer tile."""
        frag = rand_seq(200, 44)
        r1, q1 = bytearray(frag[:150]), bytearray(b"I" * 150)
        r1_edit(r1, q1, frag)
        r2 = bytearray(rc(frag[50:]))
        r2[45] = ord("N")
        r1, q1, r2 = bytes(r1), bytes(q1), bytes(r2)
        assert find_overlap(r1, rc(r2), _P).verdict == V_MERGE
        counters = [0, 0]
        res = process_pair(b"p/1", r1, q1, b"p/2", r2, b"I" * 150, _P, counters)
        assert res.outcome == PairOutcome.KEPT and len(res.records) == 2
        assert res.records[1] == (b"p/2 ZN:i:4 trim3_105", r2[:45], b"I" * 45)
        return frag, r1, q1, res, counters

    def test_a_merge_verdict_kept_by_trim3_emits_r1_as_it_was_read(self):
        """A kept mate is the input with the N policy applied and nothing else (plan §8:
        kept-mate substitutions are zero by construction). The consensus had already
        rewritten R1's low-quality error from R2 -- correctly here, but only because
        the alignment is right, and a kept pair is not one the merger vouches for.
        0.5.x emitted the rewritten R1."""
        def error_at_120(r1, q1, frag):
            r1[120] = flip(r1[120])
            q1[120] = ord("&")                      # Q5: R2's Q40 call would win
        _frag, r1, q1, res, counters = self._merge_verdict_then_trim3_breaks_tiling(
            error_at_120)
        assert res.records[0] == (b"p/1", r1, q1)
        assert counters == [0, 0], "no consensus reached an emitted base"

    def test_a_merge_verdict_kept_by_trim3_does_not_keep_a_rescue(self):
        """The same for an N rescued from the mate: the kept R1 has its N back, so trim3
        cuts it there, and the record carries no rescue."""
        def n_at_120(r1, q1, frag):
            r1[120] = ord("N")
        frag, _r1, _q1, res, _c = self._merge_verdict_then_trim3_breaks_tiling(n_at_120)
        assert res.records[0] == (b"p/1 ZN:i:4 trim3_30", frag[:120], b"I" * 120)

    def test_the_n_policy_counters_are_reported(self):
        """A policy that quietly eats a library is the failure mode being guarded.

        `_fold` sums a fixed prefix of the counter tuple and special-cases the maximum;
        a counter added past that point reports zero on every input. Assert they move."""
        from zna.merge.backend import get_merge_backend
        be = get_merge_backend("python")
        rng = random.Random(5)
        r1s, r2s = [], []
        for i in range(200):
            frag = draw(rng, rng.randrange(60, 300))
            a = bytearray((frag + ADAPTER1 + draw(rng, 160))[:150])
            b = bytearray((rc(frag) + ADAPTER2 + draw(rng, 160))[:150])
            a[rng.randrange(40, 150)] = ord("N")
            b[rng.randrange(40, 150)] = ord("N")
            r1s.append(b"@f%d/1\n%b\n+\n%b\n" % (i, bytes(a), b"I" * 150))
            r2s.append(b"@f%d/2\n%b\n+\n%b\n" % (i, bytes(b), b"I" * 150))
        buf1, buf2 = b"".join(r1s), b"".join(r2s)
        seen = {}
        for policy in (1, 2):                       # trim3, random
            a = be.merge_chunk(buf1, 0, len(buf1), buf2, 0, len(buf2), *_chunk_args(),
                               policy, 42)
            assert a[3][NPOLICY_BASES] > 0, "the N-policy counter never moved"
            assert a[3][N_RESCUED] > 0, "no rescue: the fixture is not exercising it"
            seen[policy] = a[3][NPOLICY_BASES]
        assert seen[1] > seen[2], (seen, "trim3 should cost more bases than random")

    def test_a_read_through_kept_pair_is_emitted_untouched(self):
        """A 20-base fragment inside 150 bp reads with two mismatches scores 23.3 bits,
        under the 2x150 floor of 28.2: no overlap, so neither mate is rewritten."""
        rng = random.Random(31)
        frag = draw(rng, 20)
        r1 = bytearray((frag + ADAPTER1 + draw(rng, 150))[:150])
        r2 = bytearray((rc(frag) + ADAPTER2 + draw(rng, 150))[:150])
        for i in (5, 12):
            r1[i] = flip(r1[i])
        r1, r2 = bytes(r1), bytes(r2)
        q1, q2 = b"!" * 150, b"~" * 150         # R2 far higher quality: it would win
        res = process_pair(b"r/1", r1, q1, b"r/2", r2, q2,
                           MergeParams(error_rate="0.01", min_read_length=40))
        assert res.outcome == PairOutcome.KEPT and res.olen == 0
        assert res.records[0][1:] == (r1, q1) and res.records[1][1:] == (r2, q2)

    def test_consensus_takes_the_higher_quality_call(self):
        """The better-supported base wins, wherever it sits in the (Q1,Q2) plane."""
        (r1, q1, r2, q2), frag = self._mismatch_pair(90, q_r1=10, q_r2=40)
        res = process_pair(b"c/1", r1, q1, b"c/2", r2, q2, P)
        assert res.outcome == PairOutcome.MERGED
        assert res.records[0][1] == frag                 # R1's error resolved from R2

    def test_consensus_acts_in_the_band_fastps_cutoffs_never_touched(self):
        """Q11 vs Q25: R2 is ~25x better supported, but fastp's gate (R1<=Q14 AND
        R2>=Q30) does not fire, so the old rule silently kept R1's error."""
        (r1, q1, r2, q2), frag = self._mismatch_pair(93, q_r1=11, q_r2=25)
        assert process_pair(b"c/1", r1, q1, b"c/2", r2, q2, P).records[0][1] == frag
        (r1, q1, r2, q2), frag = self._mismatch_pair(94, q_r1=37, q_r2=25)
        out = process_pair(b"c/1", r1, q1, b"c/2", r2, q2, P).records[0][1]
        assert out != frag and out[15] == r1[15]

    def test_equal_quality_disagreement_keeps_r1(self):
        """A tie carries no information, so nothing is rewritten (R1 is the frame)."""
        (r1, q1, r2, q2), _frag = self._mismatch_pair(91, q_r1=40, q_r2=40)
        assert process_pair(b"c/1", r1, q1, b"c/2", r2, q2, P).records[0][1][15] == r1[15]

    def test_a_contested_base_is_derated_either_way(self):
        """The output quality is the POSTERIOR of the winning call, which is always
        worse than the winner's own Q — a disputed base is less certain."""
        (r1, q1, r2, q2), _frag = self._mismatch_pair(95, q_r1=10, q_r2=40)
        rec = process_pair(b"c/1", r1, q1, b"c/2", r2, q2, P).records[0]
        assert 33 < rec[2][15] < 40 + 33
        (r1, q1, r2, q2), _frag = self._mismatch_pair(96, q_r1=37, q_r2=30)
        rec = process_pair(b"c/1", r1, q1, b"c/2", r2, q2, P).records[0]
        assert rec[2][15] < 37 + 33
        assert rec[2][0] == 37 + 33                      # uncontested: untouched

    def test_consensus_counts_via_out_param(self):
        (r1, q1, r2, q2), _frag = self._mismatch_pair(92, q_r1=10, q_r2=40)
        counters = [0, 0]
        process_pair(b"c/1", r1, q1, b"c/2", r2, q2, P, counters)
        assert counters == [1, 0]

    def test_consensus_posterior_table_is_symmetric_and_monotone(self):
        q = lambda w, l: DISAGREE_Q[(w + 33) * 256 + (l + 33)] - 33
        assert q(40, 10) > q(40, 30) > q(40, 39)         # wider gap -> higher confidence
        assert q(20, 20) <= 4                            # a tie is ~50/50, i.e. ~Q3
        for w in (20, 30, 40):
            for l in (5, 15, 25):
                assert q(w, l) <= w                      # never more certain than the call

    def test_disjoint_keeps_both(self):
        s1, s2 = rand_seq(50, 14), rand_seq(50, 15)
        res = process_pair(b"x/1", s1, qual(s1), b"x/2", s2, qual(s2), P)
        assert res.outcome == PairOutcome.KEPT and res.score == 0
        assert [r[1] for r in res.records] == [s1, s2]
        assert [r[0] for r in res.records] == [b"x/1", b"x/2"]      # /1,/2 preserved

    def test_a_short_overlap_is_kept_whole_not_trimmed(self):
        """0.5.3 trimmed a 12-base overlap (23.8 bits) off both mates. 0.6 keeps both
        whole: khorana trains on one mate of an unmerged pair and wants it entire, and
        every wrong trim came from that band."""
        frag = rand_seq(48, 13)
        (h1, s1, q1), (h2, s2, q2) = make_pair(frag, 30)
        res = process_pair(h1, s1, q1, h2, s2, q2, P)
        assert res.outcome == PairOutcome.KEPT
        assert [r[:2] for r in res.records] == [(h1, s1), (h2, s2)]

    def test_read_through_collapses_to_insert(self):
        insert = rand_seq(20, 16)
        s1 = (insert + ADAPTER1)[:30]
        s2 = (rc(insert) + ADAPTER2)[:30]
        res = process_pair(b"rt/1", s1, qual(s1), b"rt/2", s2, qual(s2), P)
        assert res.outcome == PairOutcome.MERGED
        assert res.records[0][1] == insert                         # adapter gone

    def test_length_filter_drops_short_merged(self):
        insert = rand_seq(20, 17)
        s1 = (insert + ADAPTER1)[:30]
        s2 = (rc(insert) + ADAPTER2)[:30]
        p = MergeParams(error_rate="0.01", min_read_length=40)
        res = process_pair(b"rt/1", s1, qual(s1), b"rt/2", s2, qual(s2), p)
        assert res.outcome == PairOutcome.MERGED and res.records == []
        assert res.n_dropped == 1

    def test_a_pair_with_a_short_mate_is_dropped_whole(self):
        """All-or-nothing: no lone mate is ever emitted as a spurious single."""
        s1, s2 = rand_seq(100, 40), rand_seq(45, 41)
        res = process_pair(b"k/1", s1, qual(s1), b"k/2", s2, qual(s2),
                           MergeParams(error_rate="0.01", min_read_length=50))
        assert res.outcome == PairOutcome.KEPT
        assert res.records == [] and res.n_dropped == 2

    def test_cli_defaults(self):
        # 40 = the pipeline-wide floor (must match the initial fastp run).
        args = cli.build_parser().parse_args(["--in1", "a", "--in2", "b", "--out", "c"])
        assert args.min_read_length == 40
        assert (args.alpha, args.error_rate, args.adapter_trimmed) == \
            ("1e-6", "0.01", False)
        assert not hasattr(args, "t_merge") and not hasattr(args, "t_trim")

    def test_fully_redundant_r2_collapses_to_merged_insert(self):
        # R2 (15 bp) fully inside R1 -> R1 spans the short insert -> merged single.
        frag = rand_seq(40, 18)
        r1 = frag[:35]
        r2 = rc(frag[5:20])          # s2rc == frag[5:20] aligns at R1 offset 5, olen 15
        res = process_pair(b"e/1", r1, qual(r1), b"e/2", r2, qual(r2), P)
        assert res.outcome == PairOutcome.MERGED
        assert res.records[0][1] == frag[:20]          # clean insert (5 + 15)
        assert base_name(res.records[0][0]) == b"e"


class TestUnequalReadLengths:
    """Truncation requires ``len1 < len2`` strictly, so equal-length fixtures
    structurally cannot express it. Both mechanisms are pinned here."""

    def test_forward_shift_zero_with_r2_longer(self):
        rng = random.Random(1234)
        frag = draw(rng, 150)
        r1, r2 = frag[:100], rc(frag)
        res = process_pair(b"u/1", r1, qual(r1), b"u/2", r2, qual(r2), P)
        assert res.outcome == PairOutcome.MERGED
        assert res.records[0][1] == frag                   # not truncated to R1
        assert res.records[0][0].endswith(b"merged_100_50")

    def test_read_through_with_r1_shorter_than_the_fragment(self):
        rng = random.Random(5678)
        frag = draw(rng, 120)
        r1 = frag[:100]                                     # R1 stops 20 bases short
        r2 = (rc(frag) + ADAPTER2 + draw(rng, 150))[:150]
        res = process_pair(b"v/1", r1, qual(r1), b"v/2", r2, qual(r2), P)
        assert res.outcome == PairOutcome.MERGED
        assert res.records[0][1] == frag                    # R2 supplies frag[100:120]

    @pytest.mark.parametrize("insert", list(range(60, 300, 11)))
    def test_span_invariant_over_independent_read_lengths(self, insert):
        rng = random.Random(9000 + insert)
        frag = draw(rng, insert)
        for l1, l2 in ((150, 150), (100, 150), (150, 100), (90, 140), (140, 90)):
            r1 = (frag + ADAPTER1 + draw(rng, 200))[:l1]
            r2 = (rc(frag) + ADAPTER2 + draw(rng, 200))[:l2]
            res = process_pair(b"w/1", r1, qual(r1), b"w/2", r2, qual(r2),
                               MergeParams(error_rate="0.01", min_read_length=40))
            if res.outcome != PairOutcome.MERGED or not res.records:
                continue
            span = res.shift + len(r2)
            assert len(res.records[0][1]) == len(res.records[0][2]) == span
            assert res.records[0][1] == frag and span == insert


class TestNames:
    def test_base_name(self):
        assert base_name(b"SRR123.5/1") == b"SRR123.5"
        assert base_name(b"SRR123.5/2") == b"SRR123.5"
        assert base_name(b"SRR123.5") == b"SRR123.5"
        assert base_name(b"SRR123.5/1\tRX:Z:ACGT") == b"SRR123.5"
        assert base_name(b"SRR123.5 merged_150_87") == b"SRR123.5"

    def test_merged_name_is_fastp_style(self):
        frag = rand_seq(40, 22)                       # insert 40, L 30 -> n1=30, n2=10
        (h1, s1, q1), (h2, s2, q2) = make_pair(frag, 30, name=b"R.7")
        res = process_pair(h1, s1, q1, h2, s2, q2, P)
        name = res.records[0][0]
        assert name == b"R.7 merged_30_10"
        # khorana seq.py parse_merged_fastq requires the last token to start with 'merged'
        assert name.rsplit(None, 1)[-1].startswith(b"merged")
        assert base_name(name) == b"R.7"

    def test_merged_strips_suffix_preserves_tags(self):
        frag = rand_seq(30, 20)
        res = process_pair(b"SRR.1/1\tRX:Z:ACGT", frag, qual(frag),
                           b"SRR.1/2\tRX:Z:ACGT", rc(frag), qual(frag), P)
        assert res.records[0][0] == b"SRR.1\tRX:Z:ACGT merged_30_0"


class TestProperty:
    @pytest.mark.parametrize("insert", list(range(30, 61)))
    def test_merge_or_keep_never_double_counts(self, insert):
        """2x30 over every insert: merged exactly when the clean overlap reaches the
        floor (T(59) = 25.814 bits -> 14 bases; 13 fall 0.003 bits short), otherwise
        both reads verbatim."""
        L = 30
        frag = rand_seq(60, 999)[:insert]
        (h1, s1, q1), (h2, s2, q2) = make_pair(frag, L)
        res = process_pair(h1, s1, q1, h2, s2, q2, P)
        if 2 * L - insert >= min_matches(t_q(L, L), 0):
            assert res.outcome == PairOutcome.MERGED
            assert res.records[0][1] == frag
        else:
            assert res.outcome == PairOutcome.KEPT
            assert [r[1] for r in res.records] == [s1, s2]


# --------------------------------------------------------------------------- #
# 11. end-to-end CLI: the stats, the warnings, table growth
# --------------------------------------------------------------------------- #

def _write_fastq_gz(path, records):
    with gzip.open(path, "wb") as fh:
        for name, seq in records:
            fh.write(b"@%b\n%b\n+\n%b\n" % (name, seq, b"I" * len(seq)))


def _run(tmp_path, in1, in2, *flags, out="o.fastq"):
    args = cli.build_parser().parse_args([
        "--in1", str(in1), "--in2", str(in2), "--out", str(tmp_path / out), "-q",
        *flags])
    return cli.run(args)


def _library_files(tmp_path, n, seed, read_len=100, with_n=False, name="lib"):
    rng = random.Random(seed)
    r1s, r2s = [], []
    for i in range(n):
        frag = draw(rng, rng.randrange(60, 2 * read_len + 40))
        (h1, s1, _), (h2, s2, _) = cycle_pair(frag, read_len, rng,
                                              name=b"P%d" % i)
        s1 = mutate(s1, rng, 0.005)
        if with_n and rng.random() < 0.3:
            s1 = s1[:70] + b"N" + s1[71:]
        r1s.append((h1, s1)); r2s.append((h2, s2))
    in1, in2 = tmp_path / f"{name}_1.fastq.gz", tmp_path / f"{name}_2.fastq.gz"
    _write_fastq_gz(in1, r1s); _write_fastq_gz(in2, r2s)
    return in1, in2


class TestCLI:
    def test_end_to_end_mixed_stream(self, tmp_path):
        merge_frag = rand_seq(40, 30)                    # overlap 20 -> merge
        short_frag = rand_seq(48, 31)                    # overlap 12 -> keep (< 14)
        disj1, disj2 = rand_seq(50, 32), rand_seq(50, 33)  # -> keep both

        r1s, r2s = [], []
        (h1, s1, _), (h2, s2, _) = make_pair(merge_frag, 30, name=b"M")
        r1s.append((h1, s1)); r2s.append((h2, s2))
        (h1, s1, _), (h2, s2, _) = make_pair(short_frag, 30, name=b"T")
        r1s.append((h1, s1)); r2s.append((h2, s2))
        r1s.append((b"D/1", disj1)); r2s.append((b"D/2", disj2))
        in1, in2 = tmp_path / "r1.fastq.gz", tmp_path / "r2.fastq.gz"
        _write_fastq_gz(in1, r1s); _write_fastq_gz(in2, r2s)
        js = tmp_path / "stats.json"

        stats = _run(tmp_path, in1, in2, "--json", str(js), "--threads", "1",
                     "--min-read-length", "10", "--error-rate", "0.01")
        assert stats["input_pairs"] == 3
        assert stats["merged"] == 1 and stats["kept_pairs"] == 2
        assert stats["emitted_records"] == 1 + 2 + 2      # merged(1) + kept(2) + kept(2)
        assert stats["error_rate"] == 0.01
        assert stats["params"]["alpha"] == 1e-6
        assert stats["policy"] == "zna-merge-0.6"
        # only the admitted overlap is binned; the merged record's length is the insert
        assert stats["overlap_length_histogram"] == {"20": 1}
        assert stats["insert_size_histogram"] == {"40": 1}
        assert stats["overlap_mismatch_rate"] == 0.0      # clean synthetic overlaps

        lines = (tmp_path / "o.fastq").read_bytes().splitlines()
        headers = [lines[i][1:] for i in range(0, len(lines), 4)]
        ids = [h.split()[0] for h in headers]
        assert b"M" in ids                                # merged single, suffix stripped
        assert {b"T/1", b"T/2", b"D/1", b"D/2"} <= set(ids)
        seqs = {ids[k]: lines[k * 4 + 1] for k in range(len(ids))}
        assert seqs[b"M"] == merge_frag
        assert seqs[b"T/1"] == make_pair(short_frag, 30)[0][1]     # whole, not trimmed
        assert json.loads(js.read_text())["merged"] == 1

    def test_the_stats_are_finite_and_type_stable(self, tmp_path):
        """hulkrna's cohort gather rejects Infinity/NaN and type changes, so the stats
        of an empty input and of a real one have the same keys and value types."""
        in1, in2 = _library_files(tmp_path, 30, 1)
        full = _run(tmp_path, in1, in2)
        e1, e2 = tmp_path / "e1.fastq.gz", tmp_path / "e2.fastq.gz"
        _write_fastq_gz(e1, []); _write_fastq_gz(e2, [])
        empty = _run(tmp_path, e1, e2, "--allow-empty", out="e.fastq")
        json.dumps(full, allow_nan=False)
        json.dumps(empty, allow_nan=False)

        def types(d):
            return {k: (types(v) if isinstance(v, dict) and not k.endswith("histogram")
                        else type(v).__name__) for k, v in d.items()}
        assert types(full) == types(empty)
        for k in ("error_rate", "detected_overlap_mismatch_rate",
                  "expected_refused_true_overlap_fraction",
                  "readthrough_check_strong_fraction", "overlap_mismatch_rate",
                  "merged_pct"):
            assert isinstance(full[k], float) and isinstance(empty[k], float), k
        for gone in ("trimmed_pairs", "trimmed_pct", "bases_trimmed",
                     "trim_guard_kept_untrimmed", "error_rate_source",
                     "error_sample_pairs", "error_prior_share"):
            assert gone not in full
        assert empty["error_rate"] == 0.01
        assert (empty["detected_overlap_bases"], empty["readthrough_check_pairs"],
                empty["detected_overlap_mismatch_rate"],
                empty["expected_refused_true_overlap_fraction"],
                empty["detected_overlap_length_histogram"]) == (0, 0, 0.0, 0.0, {})
        assert full["readthrough_check_pairs"] == 30 and full["detected_overlap_bases"]

    def test_the_detected_rate_is_every_detected_overlap_before_the_gate(self,
                                                                          tmp_path):
        """``detected_overlap_mismatch_rate`` over the whole run equals the sum, pair by
        pair, of what ``find_overlap`` detects -- merged or refused -- at the run's own
        parameters, informative positions only."""
        in1, in2 = _library_files(tmp_path, 60, 2, with_n=True)
        stats = _run(tmp_path, in1, in2)
        assert stats["error_rate"] == 0.01                 # the default reached the run
        d = n = 0
        with gzip.open(in1) as f1, gzip.open(in2) as f2:
            l1, l2 = f1.read().splitlines(), f2.read().splitlines()
        for s1, s2 in zip(l1[1::4], l2[1::4]):
            o = find_overlap(s1, rc(s2), _P)
            if o.verdict != V_NONE:
                d += o.informative_mismatches
                n += o.overlap_len - (o.mismatches - o.informative_mismatches)
        assert d > 0 and stats["detected_overlap_bases"] == n
        assert stats["detected_overlap_mismatch_rate"] == round(d / n, 6)

    def test_the_expected_refusals_are_every_detected_overlap_by_brute_force(
            self, tmp_path):
        """``expected_refused_true_overlap_fraction`` over a whole run equals, by exact
        enumeration, the mean over every pair ``find_overlap`` detects of ``P(Binom(n,
        rate) > dfit[n])`` at the run's detected rate -- and the histogram it is summed
        over is those pairs' overlap lengths, merged or refused. Noisy enough (4%, a
        3'-degraded library) that the value is not a rounding-level zero."""
        rng = random.Random(41)
        r1s, r2s = [], []
        for i in range(80):
            frag = draw(rng, rng.randrange(110, 190))
            (h1, s1, _), (h2, s2, _) = make_pair(frag, 100, name=b"B%d" % i)
            r1s.append((h1, mutate(s1, rng, 0.04))); r2s.append((h2, s2))
        r1s.append((b"R/1", TestPlausibilityGate._repeat_beats_truth()[:150]))
        r2s.append((b"R/2", rc(TestPlausibilityGate._repeat_beats_truth()[-150:])))
        in1, in2 = tmp_path / "b1.fastq.gz", tmp_path / "b2.fastq.gz"
        _write_fastq_gz(in1, r1s); _write_fastq_gz(in2, r2s)
        stats = _run(tmp_path, in1, in2)
        det_d = det_n = 0
        lengths = []
        for (_h1, s1), (_h2, s2) in zip(r1s, r2s):
            o = find_overlap(s1, rc(s2), _P)
            if o.verdict != V_NONE:
                det_d += o.informative_mismatches
                det_n += o.overlap_len - (o.mismatches - o.informative_mismatches)
                lengths.append(o.overlap_len)
        assert stats["implausible_refused"] >= 1          # the repeat is in there
        assert stats["detected_overlap_length_histogram"] == {
            str(n): lengths.count(n) for n in sorted(set(lengths))}
        rate = Fraction(det_d, det_n)
        want = sum(_brute_refusal(n, _P.dfit(n), rate) for n in lengths) / len(lengths)
        assert 1e-4 < float(want) < 1e-1
        assert stats["expected_refused_true_overlap_fraction"] == \
            float(format(float(want), ".4g"))

    def test_the_error_rate_reaches_the_kernel_exactly(self, tmp_path):
        in1, in2 = _library_files(tmp_path, 60, 2)
        base = _run(tmp_path, in1, in2)
        user = _run(tmp_path, in1, in2, "--error-rate", "0.0123")
        assert user["error_rate"] == 0.0123
        assert (user["params"]["match_q"], user["params"]["step_q"]) == \
            weights_q(Fraction("0.0123"))
        assert (base["params"]["match_q"], base["params"]["step_q"]) == \
            weights_q(Fraction("0.01"))

    @pytest.mark.parametrize("threads,chunk", [(1, 7), (3, 7), (4, 13), (2, 1000),
                                               (3, 1)])
    @pytest.mark.parametrize("declared", [False, True])
    def test_the_output_and_the_diagnostics_ignore_threads_and_chunking(
            self, tmp_path, monkeypatch, threads, chunk, declared):
        """The read-through check covers the input's first READTHROUGH_CHECK_PAIRS
        pairs -- here 17 of 60, so the window ends inside a chunk at most chunk sizes --
        and the detected-overlap counters are plain sums, so neither can move with the
        chunking or the thread count. Neither can the output."""
        monkeypatch.setattr(zparams, "READTHROUGH_CHECK_PAIRS", 17)
        in1, in2 = _library_files(tmp_path, 60, 5, with_n=True)
        flags = ("--adapter-trimmed",) if declared else ()
        base = _run(tmp_path, in1, in2, "--threads", "1", "--chunk-size", "50000",
                    *flags, out="base.fastq")
        other = _run(tmp_path, in1, in2, "--threads", str(threads), "--chunk-size",
                     str(chunk), *flags, out="other.fastq")
        for s in (base, other):
            for wallclock in ("elapsed_s", "pairs_per_second"):
                s.pop(wallclock, None)
        assert other == base
        assert (tmp_path / "other.fastq").read_bytes() == \
            (tmp_path / "base.fastq").read_bytes()
        # ...and the check counted exactly the first 17, by brute force
        assert base["readthrough_check_pairs"] == 17
        with gzip.open(in1) as f1, gzip.open(in2) as f2:
            l1, l2 = f1.read().splitlines(), f2.read().splitlines()
        strong = 0
        for s1, s2 in list(zip(l1[1::4], l2[1::4]))[:17]:
            a = scan_unrestricted(s1, rc(s2), _P)
            strong += bool(a.overlap_len) and a.shift + len(s2) < max(len(s1), len(s2))
        assert 0 < strong < 17
        assert base["readthrough_check_strong_fraction"] == round(strong / 17, 6)

    def test_alpha_reaches_the_kernel(self, tmp_path):
        """A 20-base clean overlap (39.7 bits) merges at the default alpha and not at
        alpha = 1e-12 (T = 45.7 bits at 2x30)."""
        frag = rand_seq(40, 34)
        (h1, s1, _), (h2, s2, _) = make_pair(frag, 30, name=b"M")
        in1, in2 = tmp_path / "r1.fastq.gz", tmp_path / "r2.fastq.gz"
        _write_fastq_gz(in1, [(h1, s1)]); _write_fastq_gz(in2, [(h2, s2)])
        common = ("--min-read-length", "10", "--error-rate", "0.01")
        assert _run(tmp_path, in1, in2, *common)["merged"] == 1
        tight = _run(tmp_path, in1, in2, *common, "--alpha", "1e-12")
        assert tight["merged"] == 0 and tight["params"]["alpha"] == 1e-12

    def test_the_adapter_trimmed_declaration_reaches_the_kernel_and_is_checked(
            self, tmp_path, caplog):
        """Declared on raw reads full of read-through: those pairs stop merging, and
        the sample check says why, loudly."""
        rng = random.Random(6)
        r1s, r2s = [], []
        for i in range(40):
            frag = draw(rng, rng.randrange(50, 90))          # all read-through at 100
            (h1, s1, _), (h2, s2, _) = cycle_pair(frag, 100, rng, name=b"R%d" % i)
            r1s.append((h1, s1)); r2s.append((h2, s2))
        in1, in2 = tmp_path / "r1.fastq.gz", tmp_path / "r2.fastq.gz"
        _write_fastq_gz(in1, r1s); _write_fastq_gz(in2, r2s)
        free = _run(tmp_path, in1, in2, "--error-rate", "0.01")
        assert free["merged"] == 40 and not free["adapter_trimmed"]
        assert free["readthrough_check_strong_fraction"] == 1.0
        assert not [r for r in caplog.records if r.levelname == "WARNING"]
        with caplog.at_level("WARNING", logger="zna.merge"):
            declared = _run(tmp_path, in1, in2, "--error-rate", "0.01",
                            "--adapter-trimmed")
        assert declared["merged"] == 0 and declared["adapter_trimmed"]
        assert any("--adapter-trimmed was declared" in r.getMessage()
                   for r in caplog.records)

    @staticmethod
    def _noisy_files(tmp_path, n=60, seed=7):
        """2x100 pairs with a 50-base true overlap carrying 4 disagreements: 8%."""
        rng = random.Random(seed)
        r1s, r2s = [], []
        for i in range(n):
            frag = draw(rng, 150)
            s1 = bytearray(frag[:100])
            for k in range(50, 100, 16):
                s1[k] = flip(s1[k])
            r1s.append((b"N%d/1" % i, bytes(s1)))
            r2s.append((b"N%d/2" % i, rc(frag[50:])))
        in1, in2 = tmp_path / "n1.fastq.gz", tmp_path / "n2.fastq.gz"
        _write_fastq_gz(in1, r1s); _write_fastq_gz(in2, r2s)
        return in1, in2

    def test_a_library_noisier_than_the_error_rate_is_warned_about(self, tmp_path,
                                                                    caplog):
        """8% disagreement against the default 1%: every detected overlap is 50 bases,
        so the gate is expected to refuse P(Binom(50, 0.08) > dfit[50] = 6) = 10% of
        true overlaps. The run says so, names the rate it saw, and suggests a value --
        and at that value it is quiet."""
        in1, in2 = self._noisy_files(tmp_path)
        with caplog.at_level("WARNING", logger="zna.merge"):
            stats = _run(tmp_path, in1, in2)
        assert stats["detected_overlap_mismatch_rate"] == 0.08
        assert stats["detected_overlap_bases"] == 60 * 50
        assert stats["detected_overlap_length_histogram"] == {"50": 60}
        expected = cli.refusal_probability(50, 6, Fraction(8, 100))
        assert stats["expected_refused_true_overlap_fraction"] == \
            float(format(expected, ".4g")) == 0.1019
        msgs = [r.getMessage() for r in caplog.records if r.levelname == "WARNING"]
        assert len(msgs) == 1
        assert "at --error-rate 0.01 the test is expected to refuse 10% of true " \
            "overlaps (warned above 0.1%)" in msgs[0]
        assert "Rerun with --error-rate 0.08 " in msgs[0]
        # ...with what the gate refused, and the price of following the advice
        assert f"this run refused {stats['implausible_refused']} as implausible" in msgs[0]
        assert "false ones included" in msgs[0] and "745 at 0.036" in msgs[0]

        caplog.clear()
        with caplog.at_level("WARNING", logger="zna.merge"):
            again = _run(tmp_path, in1, in2, "--error-rate", "0.08")
        assert not [r for r in caplog.records if r.levelname == "WARNING"]
        # A diagnostic, not a decision: every pair here is 4 in 50 and dfit[50] is 6
        # at 1%, so both runs merged every pair -- the warning fires on the rate alone.
        assert stats["merged"] == again["merged"] == 60

    def test_a_rate_above_the_setting_is_not_by_itself_a_warning(self, tmp_path,
                                                                  caplog):
        """The warning is about what the rate COSTS, not whether it is above the
        setting. 1.2% against the default 1% (36 of 60 overlaps with one mismatch in 50)
        is above it, and the gate is expected to refuse 2e-6 of true overlaps: silent.
        So is a clean library, and 8% against a setting of 8.01%."""
        rng = random.Random(7)
        r1s, r2s = [], []
        for i in range(60):
            frag = draw(rng, 150)
            s1 = bytearray(frag[:100])
            if i < 36:
                s1[70] = flip(s1[70])
            r1s.append((b"S%d/1" % i, bytes(s1)))
            r2s.append((b"S%d/2" % i, rc(frag[50:])))
        in1, in2 = tmp_path / "s1.fastq.gz", tmp_path / "s2.fastq.gz"
        _write_fastq_gz(in1, r1s); _write_fastq_gz(in2, r2s)
        noisy1, noisy2 = self._noisy_files(tmp_path)
        clean1, clean2 = _library_files(tmp_path, 40, 11, name="clean")
        with caplog.at_level("WARNING", logger="zna.merge"):
            slight = _run(tmp_path, in1, in2)
            _run(tmp_path, noisy1, noisy2, "--error-rate", "0.0801", out="n.fastq")
            clean = _run(tmp_path, clean1, clean2, out="c.fastq")
        assert not [r for r in caplog.records if r.levelname == "WARNING"]
        assert slight["detected_overlap_mismatch_rate"] == 0.012
        assert slight["expected_refused_true_overlap_fraction"] == float(format(
            cli.refusal_probability(50, 6, Fraction(12, 1000)), ".4g")) == 2.277e-06
        assert 0 < clean["detected_overlap_mismatch_rate"] < 0.01

    def test_the_threshold_is_strict_and_exact(self, monkeypatch):
        """Warned ABOVE one in a thousand: a run whose fraction equals the threshold is
        silent, and one a hair over it is not. Compared as exact rationals."""
        acc = cli._new_acc()
        acc[0][DET_BASES], acc[0][DET_MISMATCHES], acc[0][MAX_READ_LEN] = 3000, 240, 100
        acc[4].extend([0] * 50 + [60])
        at = Fraction(cli.run_refused_fraction(acc, _P))
        monkeypatch.setattr(cli, "_WARN_REFUSED", at)
        assert cli.run_warnings(acc, _P) == []
        monkeypatch.setattr(cli, "_WARN_REFUSED", at - Fraction(1, 10 ** 60))
        (msg,) = cli.run_warnings(acc, _P)
        assert "expected to refuse" in msg

    def test_a_large_alpha_is_named_rather_than_the_error_rate(self):
        """At a rate within the setting the gate refuses a true overlap with
        probability below alpha, by construction -- so only an --alpha above 1e-3 can
        warn there, and raising --error-rate is not the advice. At --alpha 0.01 and a
        detected 1% on 100-base overlaps it is 0.34%."""
        acc = cli._new_acc()
        acc[0][DET_BASES], acc[0][DET_MISMATCHES], acc[0][MAX_READ_LEN] = 10_000, 100, 150
        acc[4].extend([0] * 100 + [100])
        loose = MergeParams(alpha="0.01")
        assert round(float(cli.run_refused_fraction(acc, loose)), 4) == 0.0034
        (msg,) = cli.run_warnings(acc, loose)
        assert "this is --alpha 0.01 itself" in msg and "Rerun" not in msg
        assert cli.run_warnings(acc, _P) == []
        # ...and the rate is what decides it, not its rounded suggestion: 0.0122 under a
        # setting of 0.0123 rounds UP to 0.013, above the setting, and is still not a
        # reason to raise it.
        acc[0][DET_MISMATCHES] = 122
        (msg,) = cli.run_warnings(acc, MergeParams(alpha="0.01", error_rate="0.0123"))
        assert "itself" in msg and "Rerun" not in msg

    def test_a_suggestion_is_only_ever_a_value_the_flag_accepts(self):
        assert cli._suggest_error_rate(Fraction(137, 10000)) == "0.014"
        assert cli._suggest_error_rate(Fraction(7, 10)) == "0.7"
        assert cli._suggest_error_rate(Fraction(746, 1000)) is None   # would be 0.75
        acc = cli._new_acc()
        acc[0][DET_BASES], acc[0][DET_MISMATCHES], acc[0][MAX_READ_LEN] = 100, 76, 150
        acc[4].extend([0] * 100 + [1])
        (msg,) = cli.run_warnings(acc, _P)
        assert "No --error-rate describes" in msg and "Rerun" not in msg

    def test_a_policy_under_which_nothing_can_merge_says_so(self):
        """Not a failure -- every pair is kept whole, correctly -- but almost certainly
        not what was meant, and otherwise silent."""
        acc = cli._new_acc()
        acc[0][N_PAIRS], acc[0][MAX_READ_LEN] = 10, 150
        assert cli.run_warnings(acc, _P) == []
        (msg,) = cli.run_warnings(acc, MergeParams(alpha="1e-300"))
        assert "no pair in this run could merge" in msg and "longest read is 150" in msg
        assert cli.run_warnings(cli._new_acc(), MergeParams(alpha="1e-300")) == []

    def test_the_warnings_survive_quiet(self, tmp_path, capsys):
        """A warning means a parameter may be wrong for this library, so -q does not
        silence it -- through the real entry point, whose logging -q configures."""
        in1, in2 = self._noisy_files(tmp_path)
        backend = "accel" if _fast_backend() else "python"
        argv = ["--in1", str(in1), "--in2", str(in2), "--out", str(tmp_path / "o.fq"),
                "-q", "--backend", backend]
        import logging
        root = logging.getLogger()
        saved_handlers, saved_level = root.handlers[:], root.level
        root.handlers[:] = []                  # so run_command's basicConfig applies
        try:
            assert cli.main(argv) == 0
        finally:
            root.handlers[:] = saved_handlers
            root.setLevel(saved_level)
        assert "at --error-rate 0.01 the test is expected to refuse" in \
            capsys.readouterr().err

    def test_histograms_are_not_capped_at_1024(self, tmp_path):
        """Every histogram bins the real value, however long the reads are -- and the
        policy tables grow to match (700 bp reads need capacity 1,024)."""
        readlen = 700
        fragments = [rand_seq(L, 900 + L) for L in (900, 1050, 1200, 1350)]
        r1s, r2s = [], []
        for frag in fragments:
            (h1, s1, _), (h2, s2, _) = make_pair(frag, readlen, name=b"L%d" % len(frag))
            r1s.append((h1, s1)); r2s.append((h2, s2))
        in1, in2 = tmp_path / "r1.fastq.gz", tmp_path / "r2.fastq.gz"
        _write_fastq_gz(in1, r1s); _write_fastq_gz(in2, r2s)
        stats = _run(tmp_path, in1, in2, "--threads", "2", "--chunk-size", "1")
        assert stats["merged"] == 4, stats
        lengths = {str(len(f)): 1 for f in fragments}
        assert stats["length_histogram"] == lengths
        assert stats["insert_size_histogram"] == lengths
        assert stats["overlap_length_histogram"] == \
            {str(2 * readlen - len(f)): 1 for f in fragments}
        assert stats["max_read_length"] == readlen
        assert stats["insert_size_censoring"]["at_read_length"] == readlen

    def test_insert_size_censoring_reports_the_cap_at_the_longest_reads(self, tmp_path):
        frag = rand_seq(40, 30)
        (h1, s1, _), (h2, s2, _) = make_pair(frag, 30)
        in1, in2 = tmp_path / "r1.fastq.gz", tmp_path / "r2.fastq.gz"
        _write_fastq_gz(in1, [(h1, s1)]); _write_fastq_gz(in2, [(h2, s2)])
        stats = _run(tmp_path, in1, in2, "--min-read-length", "25", "--error-rate",
                     "0.01")
        assert stats["insert_size_censoring"] == {
            "min_mergeable_overlap": min_matches(t_q(30, 30), 0), "at_read_length": 30}
        assert stats["params"]["min_read_length"] == 25

    def test_short_mate_fragment_dropped_and_logged(self, tmp_path):
        good1, good2 = rand_seq(100, 60), rand_seq(100, 61)
        smate1, smate2 = rand_seq(100, 62), rand_seq(30, 63)   # R2 30 bp (< 50)
        in1, in2 = tmp_path / "r1.fastq.gz", tmp_path / "r2.fastq.gz"
        _write_fastq_gz(in1, [(b"G/1", good1), (b"S/1", smate1)])
        _write_fastq_gz(in2, [(b"G/2", good2), (b"S/2", smate2)])
        stats = _run(tmp_path, in1, in2, "--min-read-length", "50")
        assert stats["fragments_dropped_short_mate"] == 1
        assert stats["dropped_below_min_length"] == 2
        assert stats["emitted_records"] == 2
        ids = [l[1:].split()[0] for l in (tmp_path / "o.fastq").read_bytes()
               .splitlines()[0::4]]
        assert set(ids) == {b"G/1", b"G/2"}

    def test_empty_input_fails_loudly(self, tmp_path):
        in1, in2 = tmp_path / "r1.fastq.gz", tmp_path / "r2.fastq.gz"
        _write_fastq_gz(in1, []); _write_fastq_gz(in2, [])
        with pytest.raises(SystemExit):
            _run(tmp_path, in1, in2)
        assert _run(tmp_path, in1, in2, "--allow-empty")["input_pairs"] == 0

    def _run_on(self, tmp_path, r1_bytes, r2_bytes):
        in1, in2 = tmp_path / "r1.fastq", tmp_path / "r2.fastq"
        in1.write_bytes(r1_bytes)
        in2.write_bytes(r2_bytes)
        return _run(tmp_path, in1, in2)

    @staticmethod
    def _records(n, suffix, seed0):
        return b"".join(b"@r%d/%b\n%b\n+\n%b\n" % (i, suffix, rand_seq(40, seed0 + i),
                                                    b"I" * 40) for i in range(n))

    def test_a_short_final_quality_line_is_rejected(self, tmp_path):
        good = self._records(3, b"1", 0)
        broken = good[:-13] + b"\n"
        with pytest.raises(SystemExit, match="quality"):
            self._run_on(tmp_path, broken, self._records(3, b"2", 100))

    def test_a_truncated_final_record_is_rejected(self, tmp_path):
        good = self._records(3, b"1", 0)
        with pytest.raises(SystemExit, match="truncated"):
            self._run_on(tmp_path, good[:-12], self._records(3, b"2", 100))

    def test_unequal_read_counts_are_rejected(self, tmp_path):
        with pytest.raises(SystemExit, match="unequal read counts"):
            self._run_on(tmp_path, self._records(4, b"1", 0), self._records(3, b"2", 100))
        with pytest.raises(SystemExit, match="unequal read counts"):
            self._run_on(tmp_path, self._records(3, b"1", 0), self._records(4, b"2", 100))

    @pytest.mark.parametrize("flags", [
        ["--alpha", "0"], ["--alpha", "1"], ["--alpha", "lots"],
        ["--error-rate", "0"], ["--error-rate", "0.8"], ["--error-rate", "-0.01"],
        ["--error-rate", "one percent"],
        ["--min-read-length", "-5"],
        ["--threads", "0"],
        ["--chunk-size", "0"],
    ])
    def test_nonsense_arguments_are_rejected(self, tmp_path, flags):
        in1, in2 = tmp_path / "r1.fastq.gz", tmp_path / "r2.fastq.gz"
        _write_fastq_gz(in1, [(b"A/1", rand_seq(40, 1))])
        _write_fastq_gz(in2, [(b"A/2", rand_seq(40, 2))])
        with pytest.raises(SystemExit):
            _run(tmp_path, in1, in2, *flags)

    @pytest.mark.parametrize("flag", ["--threshold-merge", "--threshold-trim"])
    def test_the_0_5_thresholds_are_gone(self, flag):
        with pytest.raises(SystemExit):
            cli.build_parser().parse_args(["--in1", "a", "--in2", "b", "--out", "c",
                                           flag, "28"])

    def test_sync_check_raises_on_desync(self, tmp_path):
        in1, in2 = tmp_path / "r1.fastq.gz", tmp_path / "r2.fastq.gz"
        _write_fastq_gz(in1, [(b"A/1", rand_seq(40, 1))])
        _write_fastq_gz(in2, [(b"B/2", rand_seq(40, 2))])
        with pytest.raises(SystemExit, match="out of sync"):
            _run(tmp_path, in1, in2)


class TestTableGrowth:
    """The policy tables cover a capacity and grow by doubling; a chunk that meets a
    longer read stops in front of it and the driver resumes after growing them. Where
    that happens must leave no trace in the output."""

    def test_a_chunk_stops_in_front_of_the_pair_that_needs_more(self):
        from zna.merge.backend import get_merge_backend
        be = get_merge_backend("python")
        rng = random.Random(81)
        recs1, recs2 = [], []
        for i, readlen in enumerate((80, 100, 300, 90)):
            frag = draw(rng, readlen * 3 // 2)
            s1, s2 = frag[:readlen], rc(frag[-readlen:])
            recs1.append(b"@r%d/1\n%b\n+\n%b\n" % (i, s1, b"I" * len(s1)))
            recs2.append(b"@r%d/2\n%b\n+\n%b\n" % (i, s2, b"I" * len(s2)))
        b1, b2 = b"".join(recs1), b"".join(recs2)
        p = MergeParams(error_rate="0.01")
        a = be.merge_chunk(b1, 0, len(b1), b2, 0, len(b2), *_chunk_args(p))
        assert a[3][N_PAIRS] == 2 and a[8] == 300
        assert b1[a[1]:].startswith(b"@r2/1")
        p.ensure(a[8])
        rest = be.merge_chunk(b1, a[1], len(b1), b2, a[2], len(b2), *_chunk_args(p, base=2))
        assert rest[3][N_PAIRS] == 2 and rest[8] == 0
        whole = be.merge_chunk(b1, 0, len(b1), b2, 0, len(b2), *_chunk_args(p))
        assert a[0] + rest[0] == whole[0]

    @pytest.mark.parametrize("threads,chunk", [(1, 2000), (1, 3), (3, 2), (2, 1)])
    def test_a_long_read_regrows_the_tables_mid_run(self, tmp_path, monkeypatch,
                                                    threads, chunk):
        """The tables start at 256 bases; a longer read forces a regrow inside the chunk
        loop, serial or threaded. Output equals the run whose tables were big enough
        from the start."""
        rng = random.Random(82)
        r1s, r2s = [], []
        for i, readlen in enumerate([100] * 6 + [600] + [100] * 5 + [1100, 90]):
            frag = draw(rng, readlen * 3 // 2)
            (h1, s1, _), (h2, s2, _) = make_pair(frag, readlen, name=b"G%d" % i)
            r1s.append((h1, s1)); r2s.append((h2, s2))
        in1, in2 = tmp_path / "g1.fastq.gz", tmp_path / "g2.fastq.gz"
        _write_fastq_gz(in1, r1s); _write_fastq_gz(in2, r2s)
        flags = ("--error-rate", "0.01", "--threads", str(threads), "--chunk-size",
                 str(chunk))
        late = _run(tmp_path, in1, in2, *flags, out="late.fastq")
        monkeypatch.setattr(zparams, "_MIN_CAPACITY", 2048)     # never regrows
        early = _run(tmp_path, in1, in2, *flags, out="early.fastq")
        assert late["merged"] == 14 and late["max_read_length"] == 1100
        for s in (late, early):
            for k in ("elapsed_s", "pairs_per_second"):
                s.pop(k, None)
        assert late == early
        assert (tmp_path / "late.fastq").read_bytes() == \
            (tmp_path / "early.fastq").read_bytes()


# --------------------------------------------------------------------------- #
# the compiled backend, and the CLI's refusal to quietly do without it
# --------------------------------------------------------------------------- #

class TestTheCompiledBackendIsRequiredByTheCLI:
    """The reference kernel is correct, so a run on it does not fail — it just takes
    ~50x longer, which at cluster scale is indistinguishable from a slow node and burns
    the whole allocation before anyone looks. So the command line refuses to choose it
    for you, while `run()` — the in-process API this suite uses — does not.
    """

    def _args(self, tmp_path, extra=()):
        in1, in2 = tmp_path / "r1.fastq.gz", tmp_path / "r2.fastq.gz"
        _write_fastq_gz(in1, [(b"A/1", rand_seq(60, 1))])
        _write_fastq_gz(in2, [(b"A/2", rand_seq(60, 2))])
        return cli.build_parser().parse_args(
            ["--in1", str(in1), "--in2", str(in2), "--out", str(tmp_path / "o.fastq"),
             "--min-read-length", "10", "-q", *extra])

    def test_the_cli_refuses_when_the_compiled_backend_is_missing(self, tmp_path,
                                                                  monkeypatch):
        monkeypatch.setattr("zna.merge.backend.available_merge_backends",
                            lambda: ["python"])
        with pytest.raises(SystemExit, match="backend python"):
            cli.run_command(self._args(tmp_path))

    def test_asking_for_the_reference_backend_by_name_is_allowed(self, tmp_path,
                                                                 monkeypatch):
        monkeypatch.setattr("zna.merge.backend.available_merge_backends",
                            lambda: ["python"])
        assert cli.run_command(self._args(tmp_path, ["--backend", "python"])) == 0

    def test_the_library_entry_point_never_refuses(self, tmp_path, monkeypatch):
        """`run()` is what every other test here calls. The guard belongs to the CLI."""
        monkeypatch.setattr("zna.merge.backend.available_merge_backends",
                            lambda: ["python"])
        assert cli.run(self._args(tmp_path))["input_pairs"] == 1

    def test_it_is_registered_as_a_zna_subcommand(self, tmp_path):
        """`zna merge` must reach this tool through zna's own top-level parser.

        Run out of process against the real entry point: an in-process check would
        miss the dispatch in zna/cli.py entirely.
        """
        import subprocess
        in1, in2 = tmp_path / "r1.fastq.gz", tmp_path / "r2.fastq.gz"
        out = tmp_path / "o.fastq"
        frag = rand_seq(40, 55)                       # overlap 20 -> merges
        (h1, s1, _), (h2, s2, _) = make_pair(frag, 30, name=b"S")
        _write_fastq_gz(in1, [(h1, s1)]); _write_fastq_gz(in2, [(h2, s2)])
        from zna.merge.backend import available_merge_backends
        # Without the compiled backend the CLI refuses by design, so name the one that
        # is actually there -- the point of this test is the dispatch, not the kernel.
        backend = "accel" if "accel" in available_merge_backends() else "python"
        proc = subprocess.run(
            [sys.executable, "-m", "zna.cli", "merge", "--in1", str(in1),
             "--in2", str(in2), "--out", str(out), "--min-read-length", "10",
             "--error-rate", "0.01", "--backend", backend, "-q"],
            capture_output=True, text=True, timeout=300)
        assert proc.returncode == 0, proc.stderr
        assert out.read_bytes().splitlines()[1] == frag


# --------------------------------------------------------------------------- #
# 12. the 19 cases of khorana's chr22 merge review (its §5-6), as fixtures
# --------------------------------------------------------------------------- #

_CASES = Path(__file__).parent / "data" / "report_cases"


def _report_cases():
    truth = json.loads((_CASES / "truth.json").read_text())["cases"]

    def fq(path):
        lines = path.read_bytes().splitlines()
        return [(lines[i][1:], lines[i + 1], lines[i + 3])
                for i in range(0, len(lines), 4)]
    return truth, fq(_CASES / "R1.fastq"), fq(_CASES / "R2.fastq")


def _classify_case(res, frag, s2):
    """M+ / M- / LOST / DROP-ok / K / LOSTp against the true molecule."""
    L = len(frag)
    if res.outcome == PairOutcome.MERGED:
        if res.shift + len(s2) == L:
            if not res.records:
                return "DROP-ok"
            # the merged record IS the molecule (in R1's orientation)
            assert res.records[0][1] in (frag, rc(frag))
            return "M+"
        return "M-" if res.records else "LOST"
    return "K" if res.records else ("LOSTp" if L >= 40 else "DROP-ok")


#: Expected outcomes under --adapter-trimmed (khorana's simulator clips reads to the
#: molecule, so the declaration is true), --min-read-length 40, alpha 1e-6.
#:
#: AT THE LIBRARY'S OWN ERROR RATE -- 1.6e-4, the disagreement of the error-free chr22
#: simulation's true overlaps (0.0003 gives the same outcomes). A mismatch then costs
#: ~12.2 bits, so a divergent repeat can no longer outscore a clean true overlap:
#: C01-C06 merge CORRECTLY, including C01's
#: true 38-base overlap that its 122/19 repeat used to beat (-> 206 - 232 < 0 bits).
#: C14's 7 mismatches in 64 exceed dfit[64] = 2: refused, kept. C08-C19 have no
#: mergeable true overlap (C19's is 11 bases, under the 15 T needs) and are kept whole.
#:
#: KNOWN RESIDUAL: C07, a PERFECT 15-base repeat (29.8 bits, over T = 28.2). It is
#: plausible under every error model and remains a wrong merge (plan §9).
EXPECTED_AT_LIBRARY_RATE = {
    "C01": "M+", "C02": "M+", "C03": "M+", "C04": "M+", "C05": "M+", "C06": "M+",
    "C07": "M-",                                                    # known residual
    "C08": "K", "C09": "K", "C10": "K", "C11": "K", "C12": "K", "C13": "K",
    "C14": "K", "C15": "K", "C16": "K", "C17": "K", "C18": "K", "C19": "K",
}
#: AT THE DEFAULT 1% -- what a run gets without --error-rate. The divergent repeats win
#: the argmax, and the gate refuses them (C01: 19 > dfit[122] = 9; C02 9 > 8; C05,
#: C09, C12, C15 11-14 > 7; C18 8 > 7): kept whole, never re-placed, so C01 and C02's
#: true overlaps are forgone rather than merged. The declaration removes the false
#: read-throughs (C03, C04, C06, C11, C17), and C03, C04, C06 then merge correctly.
#:
#: KNOWN RESIDUALS: C07 as above, and C14 -- 7 mismatches in 64 is exactly dfit[64] = 7
#: at 1%, plausible under that error model (plan §9).
EXPECTED_AT_DEFAULT = {
    "C01": "K", "C02": "K", "C03": "M+", "C04": "M+", "C05": "K", "C06": "M+",
    "C07": "M-",                                                    # known residual
    "C08": "K", "C09": "K", "C10": "K", "C11": "K", "C12": "K", "C13": "K",
    "C14": "M-",                                                    # known residual
    "C15": "K", "C16": "K", "C17": "K", "C18": "K", "C19": "K",
}
IMPLAUSIBLE_AT_DEFAULT = {"C01", "C02", "C05", "C09", "C12", "C15", "C18"}


class TestReportCases:
    """The 19 pairs of khorana's merge review, replayed under the 0.6 policy.

    Run at two error rates: the library's own and the default. The expectations are
    derived in the comments above. These 19 were SELECTED because 0.5.3 got them wrong,
    so their detected overlaps are mostly repeats and disagree far more than any real
    library -- which the run's own check reports, a statement about the selection, not
    the library (see the CLI test below).
    """

    @pytest.mark.parametrize("e,expected", [("0.00016", EXPECTED_AT_LIBRARY_RATE),
                                            ("0.0003", EXPECTED_AT_LIBRARY_RATE),
                                            ("0.01", EXPECTED_AT_DEFAULT)])
    def test_each_case(self, e, expected):
        truth, r1, r2 = _report_cases()
        assert [c["case_id"] for c in truth] == sorted(expected)
        p = MergeParams(error_rate=e, adapter_trimmed=True, min_read_length=40)
        got, implausible = {}, set()
        for c, (h1, s1, q1), (h2, s2, q2) in zip(truth, r1, r2):
            assert h1.startswith(c["read1_name"][:-2].encode())
            res = process_pair(h1, s1, q1, h2, s2, q2, p)
            got[c["case_id"]] = _classify_case(res, c["fragment_mrna"].encode(), s2)
            if res.implausible:
                implausible.add(c["case_id"])
        assert got == expected
        if e == "0.01":
            assert implausible == IMPLAUSIBLE_AT_DEFAULT
        else:
            assert implausible == {"C14"}

    def test_no_case_is_lost_and_no_mate_is_rewritten(self):
        """0.5.3 lost C03, C04 and C17 to false read-throughs and rewrote bases in the
        kept mates of C08, C10, C13 and C16 through its trim path. Under the policy,
        every kept pair is emitted exactly as read."""
        truth, r1, r2 = _report_cases()
        for e in ("0.00016", "0.01"):
            p = MergeParams(error_rate=e, adapter_trimmed=True, min_read_length=40)
            for c, (h1, s1, q1), (h2, s2, q2) in zip(truth, r1, r2):
                res = process_pair(h1, s1, q1, h2, s2, q2, p)
                assert res.records, c["case_id"]
                if res.outcome == PairOutcome.KEPT:
                    assert [r[1:] for r in res.records] == [(s1, q1), (s2, q2)]

    @pytest.mark.parametrize("e,strong,outcomes,detected", [
        # C03, C04, C06, C11, C17: every selected false read-through reaches T
        ("0.01", 5, (5, 14, 7), 0.137931),
        # at ~12 bits a mismatch, only C17's perfect 18-base one still does
        ("0.00016", 1, (7, 12, 1), 0.033654),
    ])
    def test_through_the_cli_the_fixture_trips_both_checks(self, tmp_path, caplog, e,
                                                            strong, outcomes, detected):
        """Five of the 19 (C03, C04, C06, C11, C17) are selected FALSE read-throughs, and
        most of the rest are repeats, so a declared run over this file warns twice --
        each check doing its job on an input where it cannot know the selection. The
        read-through check scores with the run's own weights, so how many of the five
        still reach T depends on --error-rate; both counts are far past 1%."""
        with caplog.at_level("WARNING", logger="zna.merge"):
            stats = _run(tmp_path, _CASES / "R1.fastq", _CASES / "R2.fastq",
                         "--adapter-trimmed", "--error-rate", e)
        assert stats["readthrough_check_pairs"] == 19
        assert stats["readthrough_check_strong_fraction"] == round(strong / 19, 6)
        assert stats["detected_overlap_mismatch_rate"] == detected
        msgs = [r.getMessage() for r in caplog.records if r.levelname == "WARNING"]
        assert any("--adapter-trimmed was declared" in m for m in msgs)
        assert any(f"at --error-rate {e} the test is expected to refuse" in m
                   for m in msgs)
        assert (stats["merged"], stats["kept_pairs"], stats["implausible_refused"]) \
            == outcomes
