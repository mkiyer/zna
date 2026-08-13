"""Tests for src/zna/merge (single-axis LR overlap scoring + merge/trim/keep).

Runs with or without the compiled backend. The reference kernel in ``_pymerge`` is the
oracle the compiled one is defined to agree with, so most of what is checked here is
checked against both.

The suite mirrors docs/READ_MERGE_REDESIGN.md §8c:

  1. threshold arithmetic      5. boundary invariant
  2. spurious detection rate   6. trim guard
  3. detection / shift recovery 7. parity with the old rule where it should hold
  4. read-through
"""
import gzip
import json
import random
import sys

import pytest

from zna.merge.overlap import (
    FORWARD, NO_OVERLAP, REVERSE, find_overlap, reverse_complement,
)
from zna.merge.params import (
    DISAGREE_Q, SCALE, MergeParams, score_weights, threshold_bits, to_bits, to_q,
)
from zna.merge.pairs import PairOutcome, base_name, process_pair
from zna.merge import cli


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #

ADAPTER1 = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER2 = b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

# Score weights at the default error rate, used to build exact expectations.
MATCH_W, MISMATCH_W = score_weights(0.01)

#: The reference backend is ~50x slower than the compiled one, so the statistical
#: sweeps below size themselves to whichever is running.
def _fast_backend():
    from zna.merge.backend import available_merge_backends
    return "accel" in available_merge_backends()

# The kernel scores in fixed point (zna/merge/params.py), so every expectation here is
# an exact integer -- no pytest.approx, and no float anywhere near a decision.
_P = MergeParams()
T_MERGE_Q, T_TRIM_Q = _P.t_merge_q, _P.t_trim_q


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
    ``L - readlen``; building reads at full length keeps the geometry unambiguous
    (redesign §8c.5).
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


# Default thresholds, but min_read_length=1 so tiny test reads survive the QC filter.
P = MergeParams(min_read_length=1)


# --------------------------------------------------------------------------- #
# the pre-redesign rule, kept ONLY to pin parity where parity is expected (§8c.7)
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
# §8c.1 threshold arithmetic — the score must not silently drift
# --------------------------------------------------------------------------- #

class TestScoreArithmetic:
    def test_weights(self):
        """+2 bits per match (log2 4), -6.23 per mismatch at e = 1%."""
        assert round(MATCH_W, 4) == 1.9855
        assert round(MISMATCH_W, 4) == 6.2288

    @pytest.mark.parametrize("threshold,expected", [
        (18.0, [10, 13, 16, 19]),      # redesign §5b table, mismatches 0..3
        (21.0, [11, 14, 17, 20]),
        (28.0, [15, 18, 21, 24]),
    ])
    def test_matches_needed(self, threshold, expected):
        assert [min_matches(to_q(threshold), d) for d in range(4)] == expected

    def test_threshold_is_derived_not_chosen(self):
        """T = log2(N / alpha) over N ~ 2*readlen shifts: 2x150 at alpha=1e-6 -> 28."""
        assert round(threshold_bits(150, 1e-6)) == 28
        assert round(threshold_bits(150, 1e-3)) == 18

    def test_kernel_score_matches_the_arithmetic(self):
        """The kernel's reported score is exactly matches*w+ - mismatches*w-."""
        rng = random.Random(4242)
        frag = draw(rng, 60)
        (_, r1, _), (_, r2, _) = make_pair(frag, 40)      # insert 60, L 40 -> overlap 20
        r2 = bytearray(r2)
        r2[-1] = ord("A") if r2[-1] != ord("A") else ord("C")   # 1 mismatch in overlap
        direction, shift, olen, diff, score = find_overlap(r1, rc(bytes(r2)))
        assert (direction, shift, olen, diff) == (FORWARD, 20, 20, 1)
        assert score == score_of(19, 1)

    def test_four_perfect_bases_falls_just_under_the_trim_threshold(self):
        """4 matches = 7.94 bits, just below the 8-bit default; 5 = 9.93, above.

        The redesign calls T_trim=8 "4 perfect matches"; with the exact match weight
        (1.9855, the value its own §5b arithmetic table requires) 4 bases fall 0.06
        bits short. Pinned so the boundary is a decision, not an accident.
        """
        assert score_bits(4) < 8.0 < score_bits(5)
        assert score_of(4) < T_TRIM_Q < score_of(5)      # ...and in fixed point too


# --------------------------------------------------------------------------- #
# backend selection
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

    def test_the_codec_extension_is_the_codec(self):
        accel = pytest.importorskip("zna._accel")
        assert hasattr(accel, "encode_block"), "zna._accel is not the codec extension"
        assert not hasattr(accel, "scan")

    def test_the_merge_extension_is_the_scan(self):
        accel = pytest.importorskip("zna.merge._accel")
        assert hasattr(accel, "scan"), "zna.merge._accel is not the scan extension"
        assert not hasattr(accel, "encode_block")

    def test_the_codec_is_still_accelerated(self):
        """The whole point of the package. A build that quietly loses this passes every
        functional test it has, which is why the conda recipe asserts it too."""
        pytest.importorskip("zna._accel")
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

    def test_auto_prefers_accel_when_it_is_built(self):
        from zna.merge.backend import available_merge_backends, get_merge_backend_name
        expected = "accel" if "accel" in available_merge_backends() else "python"
        assert get_merge_backend_name() == expected

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
            assert find_overlap(s1, rc(s2))[0] == FORWARD
        finally:
            overlap.use_backend(original)
        assert overlap.backend_name() == original


# --------------------------------------------------------------------------- #
# cross-backend equivalence: the accelerated scan must agree EXACTLY
# --------------------------------------------------------------------------- #

def _backends():
    from zna.merge.backend import available_merge_backends, get_merge_backend
    if "accel" not in available_merge_backends():
        pytest.skip("C++ merge backend not built")
    return get_merge_backend("python").scan, get_merge_backend("accel").scan


def _chunk_backends():
    """The two `merge_chunk` implementations, for the level-3 differential."""
    from zna.merge.backend import available_merge_backends, get_merge_backend
    if "accel" not in available_merge_backends():
        pytest.skip("C++ merge backend not built")
    return (get_merge_backend("python").merge_chunk,
            get_merge_backend("accel").merge_chunk)


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


def exhaustive_scan(s1, s2rc, floor_q):
    """Every shift, no pruning, no early exit — the slow truth to check a scan against."""
    out = []
    for s in range(-(len(s2rc) - 1), len(s1)):
        lo, hi = max(s, 0), min(len(s1), s + len(s2rc))
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


class TestCrossBackend:
    """The oracle and the accelerated kernel are one algorithm with two implementations.

    Equality here is exact, not approximate, and that is only possible because the score
    is an integer (params.py) and the argmax is a specified total order rather than an
    artifact of iteration order (docs/MERGE_CPP_DESIGN.md §5). Asserting the weaker
    "returned *an* argmax" would let a tie-break divergence through, and a tie-break
    divergence changes which bases a merged read is built from.
    """

    ARGS = (_P.match_q, _P.step_q, _P.t_trim_q)

    def _agree(self, s1, s2rc, label):
        py, cc = _backends()
        a = py(s1, s2rc, len(s1), len(s2rc), *self.ARGS)
        b = cc(s1, s2rc, len(s1), len(s2rc), *self.ARGS)
        assert a == b, (label, a, b)
        return a

    def test_overlapping_pairs(self):
        rng = random.Random(11)
        for i in range(400):
            frag = draw(rng, rng.randrange(40, 320))
            l1 = rng.randrange(20, 151)
            l2 = rng.randrange(20, 151)
            r1 = mutate((frag + ADAPTER1 + draw(rng, 160))[:l1], rng, 0.01)
            r2 = mutate((rc(frag) + ADAPTER2 + draw(rng, 160))[:l2], rng, 0.01)
            self._agree(r1, rc(r2), f"ovl{i}")

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
                r1[i] = ord("A") if r1[i] != ord("A") else ord("C")
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
        then minimise s.
        """
        py, cc = _backends()
        n_tied = 0
        for s1, s2rc, label in tie_fixtures():
            a = py(s1, s2rc, len(s1), len(s2rc), *self.ARGS)
            b = cc(s1, s2rc, len(s1), len(s2rc), *self.ARGS)
            assert a == b, (label, a, b)
            scored = exhaustive_scan(s1, s2rc, _P.t_trim_q)
            if not scored:
                assert a[2] == 0, label
                continue
            (want_s, want_n, want_d, want_sc), ties = argmax_by_rule(scored)
            assert a == (want_s, want_sc, want_n, want_d), (label, a, want_s, ties)
            n_tied += ties - 1
        assert n_tied >= 100, f"only {n_tied} ties exercised; the tie-break is untested"

    def test_random_bytes_not_just_nucleotides(self):
        """The kernel promises raw byte semantics; hold it to that on arbitrary input."""
        rng = random.Random(15)
        for i in range(300):
            n1, n2 = rng.randrange(0, 80), rng.randrange(0, 80)
            self._agree(bytes(rng.randrange(256) for _ in range(n1)),
                        bytes(rng.randrange(256) for _ in range(n2)), f"bytes{i}")

    def test_chunks_agree_blob_for_blob(self):
        """Level 3: the production path. Same bytes out, same counters, same histograms.

        This is the strongest check available — it needs no model of what the answer
        should be, only that two independent implementations of one specification agree
        completely on a realistic input.
        """
        py, cc = _chunk_backends()
        rng = random.Random(17)
        r1s, r2s = [], []
        for i in range(400):
            frag = draw(rng, rng.randrange(40, 320))
            l1, l2 = rng.randrange(30, 151), rng.randrange(30, 151)
            s1 = mutate((frag + ADAPTER1 + draw(rng, 160))[:l1], rng, 0.01)
            s2 = mutate((rc(frag) + ADAPTER2 + draw(rng, 160))[:l2], rng, 0.01)
            # quality varies so the posterior consensus actually has work to do
            q1 = bytes(rng.choice((70, 58, 44, 35)) for _ in range(len(s1)))
            q2 = bytes(rng.choice((70, 58, 44, 35)) for _ in range(len(s2)))
            r1s.append(b"@f%d/1 tag\n%b\n+\n%b\n" % (i, s1, q1))
            r2s.append(b"@f%d/2 tag\n%b\n+\n%b\n" % (i, s2, q2))
        buf1, buf2 = b"".join(r1s), b"".join(r2s)
        args = (_P.match_q, _P.step_q, _P.t_merge_q, _P.t_trim_q, 40, DISAGREE_Q, True, 0)

        a = py(buf1, 0, len(buf1), buf2, 0, len(buf2), *args)
        b = cc(buf1, 0, len(buf1), buf2, 0, len(buf2), *args)
        assert a[0] == b[0], "blobs differ"
        assert a[1:3] == b[1:3], "consumed byte counts differ"
        assert a[3] == b[3], "counters differ"
        assert list(a[4]) == list(b[4]) and list(a[5]) == list(b[5]) \
            and list(a[6]) == list(b[6]), "histograms differ"
        assert a[3][0] == 400 and a[3][1] > 100, a[3]      # the fixture proves nothing
        assert a[3][8] > 0, "no consensus changes: the fixture is not exercising it"

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
        args = (_P.match_q, _P.step_q, _P.t_merge_q, _P.t_trim_q, 1, DISAGREE_Q, False, 0)
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
        an entire class of input. Nothing here fed merge_chunk a long read, so nothing
        caught it. Sweep the boundary in both directions.
        """
        py, cc = _chunk_backends()
        rng = random.Random(1000 + readlen)
        frag = draw(rng, readlen * 3 // 2)
        s1, s2 = frag[:readlen], rc(frag[-readlen:])
        b1 = b"@x/1\n%b\n+\n%b\n" % (s1, b"I" * len(s1))
        b2 = b"@x/2\n%b\n+\n%b\n" % (s2, b"I" * len(s2))
        args = (_P.match_q, _P.step_q, _P.t_merge_q, _P.t_trim_q, 40, DISAGREE_Q, True, 0)
        a = py(b1, 0, len(b1), b2, 0, len(b2), *args)
        b = cc(b1, 0, len(b1), b2, 0, len(b2), *args)
        assert a[0] == b[0] and a[3] == b[3], (readlen, a[3], b[3])
        assert a[3][1] == 1, f"the fixture stopped merging at {readlen}"
        assert a[3][12] == readlen, "max_read_length is not being reported"

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
        args = (_P.match_q, _P.step_q, _P.t_merge_q, _P.t_trim_q, 40, DISAGREE_Q, True, 0)
        a = py(b1, 0, len(b1), b2, 0, len(b2), *args)
        b = cc(b1, 0, len(b1), b2, 0, len(b2), *args)
        assert a[0] == b[0], "blobs differ once the arena has grown"
        assert a[3] == b[3]
        assert a[3][0] == 7 and a[3][1] == 7, a[3]      # all seven merged

    def test_split_records_agrees(self):
        py = __import__("zna.merge._pymerge", fromlist=["x"]).split_records
        from zna.merge.backend import available_merge_backends, get_merge_backend
        if "accel" not in available_merge_backends():
            pytest.skip("C++ merge backend not built")
        cc = get_merge_backend("accel").split_records
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
                                     MergeParams(min_read_length=40))
                        for s1, q1, s2, q2 in pairs]
            finally:
                overlap.use_backend(original)

        assert run("python") == run("accel")


# --------------------------------------------------------------------------- #
# the fixed-point scale, and the argmax total order
# --------------------------------------------------------------------------- #

class TestFixedPointScale:
    """The score is computed in integers so the argmax is reproducible everywhere.

    Two things have to hold for that to be worth anything: the integers must be the
    same integers on every platform, and quantising must not move a decision.
    """

    def test_the_default_weights_are_exactly_these_integers(self):
        """Pin them. They are derived from `log2`, which is not correctly rounded and
        differs between libm implementations — so a platform where this fails would
        silently produce a different corpus from the same FASTQ. Fail loudly instead."""
        p = MergeParams()
        assert (SCALE, p.match_q, p.step_q) == (1 << 24, 33_311_170, 137_813_407)
        assert (p.t_merge_q, p.t_trim_q) == (469_762_048, 134_217_728)

    def test_step_is_the_sum_of_the_two_quantised_weights(self):
        """`score = n*match - d*step` and `score = (n-d)*match - d*mismatch` must agree,
        which they only do if `step` is quantised as the sum rather than separately."""
        p = MergeParams()
        assert p.step_q == to_q(MATCH_W) + to_q(MISMATCH_W)
        for n, d in ((40, 0), (40, 1), (150, 7), (19, 2)):
            assert n * p.match_q - d * p.step_q == score_of(n - d, d)

    def test_quantisation_flips_no_decision_over_the_reachable_domain(self):
        """The §4 enumeration, as a test.

        A decision flips only where the true float score sits within the quantisation
        error of a threshold, and over integer (n, d) that is exhaustively checkable.
        At SCALE = 2**24 the first disagreement is at an *overlap* of 32,830 bases; this
        sweeps everything up to 4,000, which is an order of magnitude past any read the
        tool will see. At 2**20 this test fails at n = 2,575, which is why the scale is
        not 2**20.
        """
        p = MergeParams()
        for threshold, tq in ((8.0, p.t_trim_q), (28.0, p.t_merge_q)):
            for n in range(1, 4001):
                # only d values that put the score anywhere near the threshold matter
                lo = max(0, int((n * MATCH_W - threshold - 5) / (MATCH_W + MISMATCH_W)) - 1)
                hi = min(n, int((n * MATCH_W - threshold + 5) / (MATCH_W + MISMATCH_W)) + 1)
                for d in range(lo, hi + 1):
                    exact = (n - d) * MATCH_W - d * MISMATCH_W
                    quant = n * p.match_q - d * p.step_q
                    assert (exact >= threshold) == (quant >= tq), (n, d, threshold)


class TestArgmaxTotalOrder:
    """`maximise score, then minimise s` — a specification, not an iteration artifact.

    This is what lets a rewritten kernel be tested for byte-exact equality instead of
    the weaker "returned *an* argmax". Random sequence essentially never ties, so the
    ties here are built deliberately.
    """

    def _check(self, s1, s2rc, label):
        p = MergeParams()
        direction, shift, olen, diff, score = find_overlap(s1, s2rc, p)
        got = None if direction == NO_OVERLAP else (
            (shift if direction == FORWARD else -shift), olen, diff)
        allsc = exhaustive_scan(s1, s2rc, p.t_trim_q)
        if not allsc:
            assert got is None, label
            return 0
        want, n_ties = argmax_by_rule(allsc)
        assert got == (want[0], want[1], want[2]), (label, got, want, n_ties)
        assert score == want[3], label
        return n_ties - 1

    def test_matches_an_unpruned_scan_on_real_reads(self, any_backend):
        rng = random.Random(9)
        for i in range(300):
            frag = draw(rng, rng.randrange(40, 110))
            l1, l2 = rng.randrange(20, 60), rng.randrange(20, 60)
            r1 = (frag + ADAPTER1 + draw(rng, 60))[:l1]
            r2 = (rc(frag) + ADAPTER2 + draw(rng, 60))[:l2]
            self._check(r1, rc(r2), f"real{i}")

    def test_ties_are_broken_towards_the_smallest_shift(self, any_backend):
        """Periodic content on unequal-length mates makes the plateau tie exactly.

        Without this the tie-break is untested: an earlier sweep over 7,000 random and
        adversarial pairs produced *zero* ties and proved nothing about it.
        """
        tied = sum(self._check(s1, s2rc, label) for s1, s2rc, label in tie_fixtures())
        assert tied >= 100, f"only {tied} ties exercised; the tie-break is untested"

    def test_a_tie_across_different_overlap_lengths_is_unreachable(self):
        """Why the rule needs no `n` key.

        Two shifts tie iff `dn * match_q == dd * step_q`, whose minimal solution is
        `dn = step_q / gcd(match_q, step_q)`. That gcd is 1, so `dn` would have to be
        larger than any conceivable read — ties can only ever occur at equal `n`.
        """
        from math import gcd
        for err in (0.001, 0.005, 0.01, 0.02, 0.05):
            p = MergeParams(err_rate=err)
            dn = p.step_q // gcd(p.match_q, p.step_q)
            assert dn > 10_000_000, (err, dn)


# --------------------------------------------------------------------------- #
# find_overlap
# --------------------------------------------------------------------------- #

class TestFindOverlap:
    def test_forward_normal_overlap(self):
        frag = rand_seq(40, 1)          # insert 40, read 30 -> overlap 20 at offset 10
        (_, r1, _), (_, r2, _) = make_pair(frag, 30)
        direction, shift, olen, diff, score = find_overlap(r1, rc(r2))
        assert direction == FORWARD
        assert shift == 10 and olen == 20 and diff == 0
        assert score == score_of(20)

    def test_full_overlap(self):
        frag = rand_seq(30, 2)          # insert == read len -> full overlap, shift 0
        (_, r1, _), (_, r2, _) = make_pair(frag, 30)
        direction, shift, olen, diff, _ = find_overlap(r1, rc(r2))
        assert direction == FORWARD and shift == 0 and olen == 30 and diff == 0

    def test_no_overlap(self):
        r1 = rand_seq(50, 3)
        r2 = rand_seq(50, 4)
        direction, _, _, _, score = find_overlap(r1, rc(r2))
        assert direction == NO_OVERLAP and score == 0

    def test_mismatch_within_budget_is_accepted(self):
        frag = rand_seq(40, 5)
        (_, r1, _), (_, r2, _) = make_pair(frag, 30)
        r2 = bytearray(r2)
        # The overlap is R2's 3' end (its last bases map to the start of revcomp(R2)).
        r2[-1] = ord("A") if r2[-1] != ord("A") else ord("C")   # 1 error in overlap
        direction, _, olen, diff, _ = find_overlap(r1, rc(bytes(r2)))
        assert direction == FORWARD and olen == 20 and diff == 1

    def test_noise_below_threshold_is_rejected(self):
        frag = rand_seq(40, 6)
        (_, r1, _), (_, r2, _) = make_pair(frag, 30)
        r2 = bytearray(r2)
        for i in range(1, 12):                 # 11 errors in a 20 bp overlap: 9*1.99
            r2[-i] = ord("A") if r2[-i] != ord("A") else ord("C")   # - 11*6.23 < 0
        direction, _, _, _, _ = find_overlap(r1, rc(bytes(r2)))
        assert direction == NO_OVERLAP

    def test_read_through_reverse(self):
        # insert (20) shorter than read (30): both reads run past into adapter.
        insert = rand_seq(20, 8)
        r1 = (insert + ADAPTER1)[:30]
        r2 = (rc(insert) + ADAPTER2)[:30]
        direction, shift, olen, diff, _ = find_overlap(r1, rc(r2))
        assert direction == REVERSE and olen == 20 and shift == 10

    def test_argmax_beats_a_spurious_short_hit(self):
        """A real 40 bp overlap wins outright over any chance 4-mer earlier in the scan.

        First-accept could be captured by the short hit; argmax cannot. This is the
        structural fix that lets T_trim stay low without suppressing merges.
        """
        rng = random.Random(99)
        for _ in range(200):
            frag = draw(rng, 160)                       # insert 160, L 100 -> overlap 40
            (_, r1, _), (_, r2, _) = make_pair(frag, 100)
            direction, shift, olen, _, score = find_overlap(r1, rc(r2))
            assert (direction, shift, olen) == (FORWARD, 60, 40)
            assert score == score_of(40)

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
            r1[i] = ord("A") if r1[i] != ord("A") else ord("C")
            direction, shift, olen, diff, score = find_overlap(bytes(r1), r2rc)
            assert (direction, shift, olen, diff) == (FORWARD, 0, 40, 1), i
            assert score == score_of(39, 1), i

    def test_unequal_read_lengths_need_no_special_case(self):
        """s is defined by the offset; the compared region is just the intersection."""
        rng = random.Random(7)
        frag = draw(rng, 180)
        r1 = frag[:120]                     # R1 quality-trimmed to 120
        r2 = rc(frag[-90:])                 # R2 quality-trimmed to 90
        direction, shift, olen, diff, _ = find_overlap(r1, rc(r2))
        assert (direction, shift, olen, diff) == (FORWARD, 90, 30, 0)   # s = 180 - 90


class TestDegenerateInputs:
    """Zero-length and 1-base reads survived every zna audit; keep it that way."""

    @pytest.mark.parametrize("s1,s2", [
        (b"", b"ACGTACGTAC"), (b"ACGTACGTAC", b""), (b"", b""), (b"A", b"T"),
    ])
    def test_no_overlap_and_no_crash(self, s1, s2):
        assert find_overlap(s1, rc(s2)) == (NO_OVERLAP, 0, 0, 0, 0)

    def test_empty_pair_is_kept_not_merged(self):
        recs, outcome, _d, score, _olen, _diff = process_pair(
            b"z/1", b"", b"", b"z/2", b"", b"", MergeParams(min_read_length=1))
        assert outcome == PairOutcome.KEPT and score == 0 and recs == []

    def test_n_is_scored_as_an_ordinary_base(self):
        """N vs N counts as a 2-bit match — inherited from the old kernel and left
        alone here, because `zna encode --npolicy drop` discards any fragment
        containing an N before it can reach the corpus. If that ever changes, N
        should stop earning evidence."""
        _d, _s, olen, diff, score = find_overlap(b"N" * 20 + b"ACGT" * 20,
                                                 rc(b"N" * 20 + b"ACGT" * 20))
        assert diff == 0 and score == score_of(olen)


# --------------------------------------------------------------------------- #
# §8c.2 spurious detection on genuinely unrelated pairs
# --------------------------------------------------------------------------- #

class TestSpuriousDetection:
    """The regression guard for the 5.17% false-positive rate of the old rule.

    2x150 unrelated reads. The old rule accepted a chance 4-mer (budget
    floor(0.20*4) = 0) at any of ~146 offsets; the LR score prices that at 7.9 bits
    and declines it.
    """

    def test_unrelated_pairs_rate(self):
        n = 20000 if _fast_backend() else 2000   # the reference kernel is ~50x slower
        rng = random.Random(20260811)
        detected = merged = 0
        for _ in range(n):
            r1 = draw(rng, 150)
            r2 = draw(rng, 150)
            direction, _, _, _, score = find_overlap(r1, rc(r2))
            if direction != NO_OVERLAP:
                detected += 1
                if score >= T_MERGE_Q:
                    merged += 1
        assert detected / n < 0.005, f"spurious detection {detected / n:.4%} (was 5.17%)"
        # A spurious *merge* is the expensive error — a chimera — and 28 bits prices
        # it at ~1e-6 per pair, so none may appear in 20k.
        assert merged == 0, f"{merged} spurious merges at 28 bits"

    def test_polya_does_not_merge(self):
        """Low-complexity tails must not produce a confident merge.

        The score has no explicit low-complexity guard (redesign §7 adds one via a
        composition-aware null), so this pins the current behaviour: a polyA run does
        pass the trim threshold but must never reach the merge threshold on its own.
        """
        rng = random.Random(5150)
        merged = 0
        for _ in range(200):
            core1, core2 = draw(rng, 130), draw(rng, 130)
            r1 = core1 + b"A" * 20                    # polyA tail on both mates
            r2 = core2 + b"A" * 20
            _d, _s, _o, _df, score = find_overlap(r1, rc(r2))
            merged += score >= T_MERGE_Q
        assert merged == 0


# --------------------------------------------------------------------------- #
# §8c.3 detection: the recovered shift must be the true one
# --------------------------------------------------------------------------- #

class TestDetection:
    def test_known_overlaps_recover_the_true_shift(self):
        """Overlaps 4..40 at 0.5% per-base error, 2x100.

        Sensitivity is set by the arithmetic, not by tuning: an overlap of 8 or less
        cannot afford a single mismatch (7*1.9855 - 6.2288 = 7.67 < 8), so detection
        there is just P(no error in the overlap) ~ 95%. From 9 bases up one mismatch
        fits inside the threshold and detection is essentially total.
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
                direction, shift, _o, _d, score = find_overlap(r1, rc(r2))
                trials += 1
                if direction == NO_OVERLAP:
                    continue
                if direction == FORWARD and shift == L - olen:
                    hits += 1
                else:
                    wrong += 1                      # a chance shift outscored the truth
            detected_by_olen[olen] = hits / reps
        # Chance wins are bounded by the same spurious rate as §8c.2.
        assert wrong / trials < 0.01, f"{wrong}/{trials} pairs detected at a wrong shift"
        assert detected_by_olen[4] == 0.0        # 4 perfect bases = 7.94 bits < 8
        for olen in range(5, 9):                 # zero-mismatch band
            assert detected_by_olen[olen] >= 0.85, (olen, detected_by_olen[olen])
        for olen in range(9, 41):                # one mismatch now fits
            assert detected_by_olen[olen] >= 0.95, (olen, detected_by_olen[olen])

    def test_merge_band_is_reached_at_15_clean_bases(self):
        rng = random.Random(6161)
        for olen, want_merge in ((14, False), (15, True), (30, True)):
            frag = draw(rng, 200 - olen)
            (h1, s1, q1), (h2, s2, q2) = make_pair(frag, 100)
            recs, outcome, _dropped, score, _olen, _diff = process_pair(h1, s1, q1, h2, s2, q2, P)
            assert score == score_of(olen)
            assert (outcome == PairOutcome.MERGED) is want_merge
            if not want_merge:
                assert outcome == PairOutcome.TRIMMED


# --------------------------------------------------------------------------- #
# §8c.4 read-through
# --------------------------------------------------------------------------- #

class TestReadThrough:
    @pytest.mark.parametrize("insert", list(range(40, 100, 7)))
    def test_insert_shorter_than_read_merges_to_the_fragment(self, insert):
        rng = random.Random(1000 + insert)
        frag = draw(rng, insert)
        (h1, s1, q1), (h2, s2, q2) = cycle_pair(frag, 100, rng)
        p = MergeParams(min_read_length=40)
        recs, outcome, dropped, _score, _olen, _diff = process_pair(h1, s1, q1, h2, s2, q2, p)
        assert outcome == PairOutcome.MERGED and dropped == 0
        assert len(recs) == 1
        assert recs[0][1] == frag                     # adapter and filler both gone
        assert len(recs[0][2]) == len(frag)


# --------------------------------------------------------------------------- #
# §8c.5 boundary invariant — what ZNA's fragment-end supervision depends on
# --------------------------------------------------------------------------- #

class TestBoundaryInvariant:
    @pytest.mark.parametrize("insert", list(range(45, 320, 9)))
    def test_base_zero_is_always_a_true_fragment_boundary(self, insert):
        """Nothing is ever removed from a read's 5' end, and a merged read is the
        fragment exactly (both of its edges are true boundaries)."""
        rng = random.Random(2000 + insert)
        frag = draw(rng, insert)
        (h1, s1, q1), (h2, s2, q2) = cycle_pair(frag, 150, rng)
        p = MergeParams(min_read_length=40)
        recs, outcome, _dropped, _score, _olen, _diff = process_pair(h1, s1, q1, h2, s2, q2, p)
        if outcome == PairOutcome.MERGED:
            assert len(recs) == 1
            assert recs[0][1] == frag                 # exactly the fragment
        else:
            assert len(recs) == 2
            # Only 3' bases may be removed, so each mate is a prefix of its read...
            assert s1.startswith(recs[0][1]) and s2.startswith(recs[1][1])
            # ...and base 0 of each mate is a true fragment end.
            assert recs[0][1][0] == frag[0]
            assert recs[1][1][0] == rc(frag)[0]

    # 2x150 with insert 286..292 puts the overlap at 8..14 bases: real, scoring
    # 15.9-27.8 bits, i.e. squarely inside the trim band and clear of chance.
    @pytest.mark.parametrize("insert", list(range(286, 293)))
    def test_trim_band_tiles_the_fragment_exactly_once(self, insert):
        """A trimmed pair covers the fragment once — no duplicated span, nothing lost
        from either 5' end."""
        rng = random.Random(3000 + insert)
        frag = draw(rng, insert)
        (h1, s1, q1), (h2, s2, q2) = cycle_pair(frag, 150, rng)
        recs, outcome, _d, score, _olen, _diff = process_pair(h1, s1, q1, h2, s2, q2,
                                                MergeParams(min_read_length=40))
        assert outcome == PairOutcome.TRIMMED and T_TRIM_Q <= score < T_MERGE_Q
        assert recs[0][1] + rc(recs[1][1]) == frag
        assert len(recs[1][1]) == len(recs[1][2])     # quality trimmed alongside

    def test_merged_read_is_in_r1s_frame(self):
        """zna's single/merged normalization assumes merged reads are R1-framed."""
        rng = random.Random(808)
        frag = draw(rng, 160)
        (h1, s1, q1), (h2, s2, q2) = cycle_pair(frag, 100, rng)
        recs, outcome, _d, _s, _olen, _diff = process_pair(h1, s1, q1, h2, s2, q2, P)
        assert outcome == PairOutcome.MERGED
        assert recs[0][1].startswith(s1)              # starts with R1, not revcomp(R2)


# --------------------------------------------------------------------------- #
# §8c.6 trim guard
# --------------------------------------------------------------------------- #

class TestTrimGuard:
    def _pair(self):
        # R1 100, R2 60, clean forward overlap 12 -> score 23.8 bits (trim band).
        rng = random.Random(515)
        frag = draw(rng, 148)                 # s = 88, L = s + 60 = 148
        r1 = frag[:100]
        r2 = rc(frag[88:])
        return r1, r2

    def test_trim_that_would_shorten_r2_below_min_keeps_both_untrimmed(self):
        r1, r2 = self._pair()
        counters = [0, 0]
        # keep2 would be 60 - 12 = 48, below min_read_length 50 -> guard fires.
        p = MergeParams(min_read_length=50)
        recs, outcome, dropped, score, _olen, _diff = process_pair(
            b"g/1", r1, qual(r1), b"g/2", r2, qual(r2), p, counters)
        assert T_TRIM_Q <= score < T_MERGE_Q
        assert outcome == PairOutcome.KEPT
        assert [r[1] for r in recs] == [r1, r2]       # both intact, fragment not lost
        assert dropped == 0
        assert counters[1] == 1                       # guard hit is counted

    def test_below_the_guard_the_trim_happens_normally(self):
        r1, r2 = self._pair()
        counters = [0, 0]
        p = MergeParams(min_read_length=40)           # 48 >= 40 -> trim as usual
        recs, outcome, dropped, _score, _olen, _diff = process_pair(
            b"g/1", r1, qual(r1), b"g/2", r2, qual(r2), p, counters)
        assert outcome == PairOutcome.TRIMMED
        assert len(recs[1][1]) == 48 and dropped == 0
        assert counters[1] == 0

    def test_a_pair_short_before_trimming_is_still_dropped_whole(self):
        """The guard rescues a trim, not a genuinely short mate: all-or-nothing holds."""
        s1 = rand_seq(100, 40)
        s2 = rand_seq(45, 41)                         # disjoint, and below the floor
        recs, outcome, dropped, _score, _olen, _diff = process_pair(
            b"k/1", s1, qual(s1), b"k/2", s2, qual(s2), MergeParams(min_read_length=50))
        assert outcome == PairOutcome.KEPT
        assert recs == [] and dropped == 2            # no lone read emitted


# --------------------------------------------------------------------------- #
# §8c.7 parity with the old rule where parity is expected
# --------------------------------------------------------------------------- #

class TestParityWithLegacyRule:
    def test_clean_long_overlaps_merge_at_the_same_shift(self):
        """Overlap >= 30 and clean: old and new must agree on the shift.

        Divergence is expected only in the short/noisy band, which is exactly where
        the old rule was measured to be wrong.
        """
        L = 100
        rng = random.Random(4711)
        for olen in range(30, 81, 5):
            for _ in range(20):
                frag = draw(rng, 2 * L - olen)
                (_, r1, _), (_, r2, _) = make_pair(frag, L)
                r2rc = rc(r2)
                old_dir, old_shift, _o, _d = legacy_scan(r1, r2rc, len(r1), len(r2rc),
                                                         3, 3, 0.20)
                new_dir, new_shift, _ol, _df, score = find_overlap(r1, r2rc)
                assert (old_dir, old_shift) == (new_dir, new_shift) == (FORWARD, L - olen)
                assert score >= T_MERGE_Q

    def test_legacy_rule_accepts_chance_four_mers_and_the_new_one_does_not(self):
        """Pins WHY the rules differ: same 20k unrelated pairs, 5.17% vs ~0.2%."""
        n = 5000 if _fast_backend() else 400
        rng = random.Random(2468)
        old_hits = new_hits = 0
        for _ in range(n):
            r1 = draw(rng, 150)
            r2rc = rc(draw(rng, 150))
            old_hits += legacy_scan(r1, r2rc, 150, 150, 3, 3, 0.20)[0] != 0
            new_hits += find_overlap(r1, r2rc)[0] != NO_OVERLAP
        assert old_hits / n > 0.03            # the defect this redesign removes
        assert new_hits / n < 0.005


# --------------------------------------------------------------------------- #
# process_pair
# --------------------------------------------------------------------------- #

class TestProcessPair:
    def test_full_overlap_merges_to_r1(self):
        frag = rand_seq(30, 10)
        (h1, s1, q1), (h2, s2, q2) = make_pair(frag, 30)
        recs, outcome, dropped, _score, _olen, _diff = process_pair(h1, s1, q1, h2, s2, q2, P)
        assert outcome == PairOutcome.MERGED and dropped == 0
        assert len(recs) == 1
        assert recs[0][1] == frag                       # merged == fragment
        assert recs[0][0] == b"frag merged_30_0"        # /1 stripped + fastp merged token

    def test_partial_overlap_reconstructs_fragment(self):
        frag = rand_seq(40, 11)
        (h1, s1, q1), (h2, s2, q2) = make_pair(frag, 30)
        recs, outcome, _d, _s, _olen, _diff = process_pair(h1, s1, q1, h2, s2, q2, P)
        assert outcome == PairOutcome.MERGED
        assert recs[0][1] == frag
        assert len(recs[0][2]) == len(frag)             # quality length matches

    def test_r1_wins_keeps_r1_base_on_r2_error(self):
        frag = rand_seq(40, 12)
        (h1, s1, q1), (h2, s2, q2) = make_pair(frag, 30)
        s2 = bytearray(s2)
        s2[-1] = ord("A") if s2[-1] != ord("A") else ord("C")  # error in R2's overlap end
        recs, outcome, _d, _s, _olen, _diff = process_pair(h1, s1, q1, h2, bytes(s2), q2, P)
        assert outcome == PairOutcome.MERGED
        assert recs[0][1] == frag                       # R1's (correct) base wins

    def _mismatch_pair(self, seed, q_r1, q_r2, pos=15):
        """Insert 40, read 30 -> overlap R1[10:30]. R1 carries an error at `pos` called
        at `q_r1`; R2 has the correct base at `q_r2`. Returns (args, fragment)."""
        frag = rand_seq(40, seed)
        r1 = bytearray(frag[:30]); q1 = bytearray(bytes([q_r1 + 33]) * 30)
        r1[pos] = ord("A") if r1[pos] != ord("A") else ord("C")
        r2 = rc(frag[10:40]); q2 = bytes([q_r2 + 33]) * 30
        return (bytes(r1), bytes(q1), r2, q2), frag

    def test_consensus_takes_the_higher_quality_call(self):
        """The whole point: the better-supported base wins, wherever it sits in the
        (Q1,Q2) plane."""
        (r1, q1, r2, q2), frag = self._mismatch_pair(90, q_r1=10, q_r2=40)
        recs, outcome, _d, _s, _olen, _diff = process_pair(b"c/1", r1, q1, b"c/2", r2, q2,
                                             MergeParams(min_read_length=1))
        assert outcome == PairOutcome.MERGED
        assert recs[0][1] == frag                        # R1's error resolved from R2

    def test_consensus_acts_in_the_band_fastps_cutoffs_never_touched(self):
        """Q11 vs Q25: R2 is ~25x better supported, but fastp's gate (R1<=Q14 AND
        R2>=Q30) does not fire, so the old rule silently kept R1's error. This band is
        ~5% of real overlap mismatches and R1-wins is ~95% wrong on it."""
        (r1, q1, r2, q2), frag = self._mismatch_pair(93, q_r1=11, q_r2=25)
        recs, _o, _d, _s, _olen, _diff = process_pair(b"c/1", r1, q1, b"c/2", r2, q2,
                                        MergeParams(min_read_length=1))
        assert recs[0][1] == frag
        # ...and the reverse band: R1 better than R2 -> R1 stands.
        (r1, q1, r2, q2), frag = self._mismatch_pair(94, q_r1=37, q_r2=25)
        recs, _o, _d, _s, _olen, _diff = process_pair(b"c/1", r1, q1, b"c/2", r2, q2,
                                        MergeParams(min_read_length=1))
        assert recs[0][1] != frag and recs[0][1][15] == r1[15]

    def test_equal_quality_disagreement_keeps_r1(self):
        """A tie carries no information, so nothing is rewritten (R1 is the frame)."""
        (r1, q1, r2, q2), frag = self._mismatch_pair(91, q_r1=40, q_r2=40)
        recs, _o, _d, _s, _olen, _diff = process_pair(b"c/1", r1, q1, b"c/2", r2, q2,
                                        MergeParams(min_read_length=1))
        assert recs[0][1][15] == r1[15]

    def test_a_contested_base_is_derated_either_way(self):
        """The output quality is the POSTERIOR of the winning call, which is always
        worse than the winner's own Q — a disputed base is less certain than an
        uncontested one. This is what fastp's copy-the-mate's-Q rule got wrong."""
        # R2 wins (Q40 vs Q10): posterior error ~1e-4/(1e-4+0.1) -> ~Q30, not Q40.
        (r1, q1, r2, q2), _frag = self._mismatch_pair(95, q_r1=10, q_r2=40)
        recs, _o, _d, _s, _olen, _diff = process_pair(b"c/1", r1, q1, b"c/2", r2, q2,
                                        MergeParams(min_read_length=1))
        assert 33 < recs[0][2][15] < 40 + 33
        # R1 wins (Q37 vs Q30) but is still contested, so its Q comes down too.
        (r1, q1, r2, q2), _frag = self._mismatch_pair(96, q_r1=37, q_r2=30)
        recs, _o, _d, _s, _olen, _diff = process_pair(b"c/1", r1, q1, b"c/2", r2, q2,
                                        MergeParams(min_read_length=1))
        assert recs[0][2][15] < 37 + 33
        # An uncontested position keeps its original quality untouched.
        assert recs[0][2][0] == 37 + 33

    def test_consensus_counts_via_out_param(self):
        """The counter accumulates bases the consensus actually changed."""
        (r1, q1, r2, q2), _frag = self._mismatch_pair(92, q_r1=10, q_r2=40)
        counters = [0, 0]
        process_pair(b"c/1", r1, q1, b"c/2", r2, q2,
                     MergeParams(min_read_length=1), counters)
        assert counters[0] == 1

    def test_consensus_posterior_table_is_symmetric_and_monotone(self):
        """Pin the table itself: a bigger quality gap means more confidence in the
        winner, and the roles are symmetric."""
        from zna.merge.params import DISAGREE_Q as T
        q = lambda w, l: T[(w + 33) * 256 + (l + 33)] - 33
        assert q(40, 10) > q(40, 30) > q(40, 39)         # wider gap -> higher confidence
        assert q(30, 10) == q(30, 10)
        assert q(20, 20) <= 4                            # a tie is ~50/50, i.e. ~Q3
        for w in (20, 30, 40):
            for l in (5, 15, 25):
                assert q(w, l) <= w                      # never more certain than the call

    def test_short_overlap_trims_r2(self):
        # insert 48, read 30 -> overlap 12 -> 23.8 bits: real, but under 28 -> trim.
        frag = rand_seq(48, 13)
        (h1, s1, q1), (h2, s2, q2) = make_pair(frag, 30)
        recs, outcome, _d, score, _olen, _diff = process_pair(h1, s1, q1, h2, s2, q2, P)
        assert outcome == PairOutcome.TRIMMED and T_TRIM_Q <= score < T_MERGE_Q
        assert len(recs) == 2
        r1_out, r2_out = recs[0], recs[1]
        assert r1_out[1] == s1 and r1_out[0] == h1                 # R1 full, /1 kept
        assert r2_out[0] == h2                                     # /2 kept
        assert len(r2_out[1]) == 30 - 12                           # trimmed by overlap
        # R1 + trimmed-R2 tile the fragment exactly once, no duplicated span:
        assert s1 + rc(r2_out[1]) == frag
        assert len(r2_out[1]) == len(r2_out[2])                    # qual trimmed too

    def test_trim_removes_the_full_overlap_including_mismatches(self):
        """Trimming removes the FULL detected overlap, not up to the first mismatch.
        19 bp overlap with 2 mismatches scores 21.3 bits: real, but below 28."""
        L = 100
        frag = rand_seq(181, 77)                # insert 181, L 100 -> 19 bp overlap
        r1 = frag[:L]
        r2 = bytearray(rc(frag[-L:]))           # s2rc[:19] == frag[81:100] (the overlap)
        # inject 2 mismatches into the overlap: s2rc[0],s2rc[1] map to r2[-1],r2[-2]
        r2[-1] = ord("A") if r2[-1] != ord("A") else ord("C")
        r2[-2] = ord("A") if r2[-2] != ord("A") else ord("C")
        r2 = bytes(r2)
        recs, outcome, _d, score, _olen, _diff = process_pair(b"m/1", r1, qual(r1), b"m/2", r2, qual(r2),
                                                MergeParams(min_read_length=1))
        assert outcome == PairOutcome.TRIMMED
        assert score == score_of(17, 2)
        assert len(recs[1][1]) == L - 19                    # ALL 19 overlap bp trimmed off R2
        assert r1 + rc(recs[1][1]) == frag                  # clean tail retained; tiles once

    def test_overlap_below_the_trim_threshold_is_not_trimmed(self):
        """An overlap whose evidence is under T_trim leaves both reads untouched."""
        L = 100
        frag = rand_seq(181, 78)
        r1 = frag[:L]
        r2 = bytearray(rc(frag[-L:]))
        for i in range(1, 6):                   # 5 mismatches in 19 bp: 27.8 - 31.1 < 0
            r2[-i] = ord("A") if r2[-i] != ord("A") else ord("C")
        r2 = bytes(r2)
        recs, outcome, _d, _s, _olen, _diff = process_pair(b"n/1", r1, qual(r1), b"n/2", r2, qual(r2),
                                             MergeParams(min_read_length=1))
        assert outcome == PairOutcome.KEPT                  # overlap not accepted
        assert [r[1] for r in recs] == [r1, r2]             # both kept full, untrimmed

    def test_disjoint_keeps_both(self):
        s1 = rand_seq(50, 14)
        s2 = rand_seq(50, 15)
        h1, h2 = b"x/1", b"x/2"
        recs, outcome, _d, score, _olen, _diff = process_pair(h1, s1, qual(s1), h2, s2, qual(s2), P)
        assert outcome == PairOutcome.KEPT and score == 0
        assert [r[1] for r in recs] == [s1, s2]
        assert [r[0] for r in recs] == [h1, h2]                    # /1,/2 preserved

    def test_read_through_collapses_to_insert(self):
        insert = rand_seq(20, 16)
        s1 = (insert + ADAPTER1)[:30]
        s2 = (rc(insert) + ADAPTER2)[:30]
        recs, outcome, _d, _s, _olen, _diff = process_pair(b"rt/1", s1, qual(s1), b"rt/2", s2, qual(s2), P)
        assert outcome == PairOutcome.MERGED
        assert recs[0][1] == insert                                # adapter gone

    def test_length_filter_drops_short_merged(self):
        insert = rand_seq(20, 17)
        s1 = (insert + ADAPTER1)[:30]
        s2 = (rc(insert) + ADAPTER2)[:30]
        p = MergeParams(min_read_length=40)
        recs, outcome, dropped, _s, _olen, _diff = process_pair(b"rt/1", s1, qual(s1), b"rt/2", s2, qual(s2), p)
        assert outcome == PairOutcome.MERGED and recs == [] and dropped == 1

    def test_min_read_length_default_is_40(self):
        # 40 = the pipeline-wide floor (must match the initial fastp run).
        from zna.merge.cli import build_parser
        args = build_parser().parse_args(["--in1", "a", "--in2", "b", "--out", "c"])
        assert args.min_read_length == 40
        assert (args.t_merge, args.t_trim) == (28.0, 8.0)

    def test_fully_redundant_r2_collapses_to_merged_insert(self):
        # R2 (15 bp) fully inside the overlap -> R1 spans the short insert -> merged single.
        # The old `keep2 <= 0` special case, now reached through `score >= t_merge`.
        frag = rand_seq(40, 18)
        r1 = frag[:35]
        r2 = rc(frag[5:20])          # s2rc == frag[5:20] aligns at R1 offset 5, olen 15
        recs, outcome, _d, _s, _olen, _diff = process_pair(b"e/1", r1, qual(r1), b"e/2", r2, qual(r2), P)
        assert outcome == PairOutcome.MERGED
        assert len(recs) == 1
        assert recs[0][1] == frag[:20]                 # clean insert (shift 5 + len2 15)
        assert base_name(recs[0][0]) == b"e"           # single; merged name


# --------------------------------------------------------------------------- #
# unequal read lengths: the regime every equal-length fixture is blind to
# --------------------------------------------------------------------------- #

class TestUnequalReadLengths:
    """Truncation requires ``len1 < len2`` strictly, so `make_pair`/`cycle_pair` — which
    build both mates at one length — structurally cannot express it. That is why 135
    tests were green over a defect affecting 0.271% of merged records in production.
    Both mechanisms are pinned here.
    """

    def test_forward_shift_zero_with_r2_longer(self):
        """R1 quality-trimmed; R2 spans the whole fragment (s == 0, len2 > len1)."""
        rng = random.Random(1234)
        frag = draw(rng, 150)
        r1 = frag[:100]
        r2 = rc(frag)
        recs, outcome, _d, _s, _olen, _diff = process_pair(b"u/1", r1, qual(r1), b"u/2", r2, qual(r2), P)
        assert outcome == PairOutcome.MERGED
        assert recs[0][1] == frag                       # not truncated to R1
        assert recs[0][0].endswith(b"merged_100_50")    # 100 from R1, 50 from R2

    def test_read_through_with_r1_shorter_than_the_fragment(self):
        """The majority mechanism: read-through (s < 0) where R1 does not reach L."""
        rng = random.Random(5678)
        frag = draw(rng, 120)
        r1 = frag[:100]                                  # R1 stops 20 bases short of L
        r2 = (rc(frag) + ADAPTER2 + draw(rng, 150))[:150]
        recs, outcome, _d, _s, _olen, _diff = process_pair(b"v/1", r1, qual(r1), b"v/2", r2, qual(r2), P)
        assert outcome == PairOutcome.MERGED
        assert recs[0][1] == frag                        # R2 supplies frag[100:120]

    @pytest.mark.parametrize("insert", list(range(60, 300, 11)))
    def test_span_invariant_over_independent_read_lengths(self, insert):
        """For every merged record, len(seq) == the fragment span the scan inferred,
        and the record equals the fragment exactly."""
        rng = random.Random(9000 + insert)
        frag = draw(rng, insert)
        for l1, l2 in ((150, 150), (100, 150), (150, 100), (90, 140), (140, 90)):
            r1 = (frag + ADAPTER1 + draw(rng, 200))[:l1]
            r2 = (rc(frag) + ADAPTER2 + draw(rng, 200))[:l2]
            direction, shift, _olen, _diff, score = find_overlap(r1, rc(r2))
            recs, outcome, _d, _s, _olen, _diff = process_pair(
                b"w/1", r1, qual(r1), b"w/2", r2, qual(r2), MergeParams(min_read_length=40))
            if outcome != PairOutcome.MERGED or not recs:
                continue
            span = (shift if direction == FORWARD else -shift) + len(r2)
            assert len(recs[0][1]) == span, (insert, l1, l2)
            assert len(recs[0][2]) == span                       # quality tracks sequence
            assert recs[0][1] == frag[:span] and span == insert  # exactly the fragment


# --------------------------------------------------------------------------- #
# name helpers
# --------------------------------------------------------------------------- #

class TestNames:
    def test_base_name(self):
        assert base_name(b"SRR123.5/1") == b"SRR123.5"
        assert base_name(b"SRR123.5/2") == b"SRR123.5"
        assert base_name(b"SRR123.5") == b"SRR123.5"
        assert base_name(b"SRR123.5/1\tRX:Z:ACGT") == b"SRR123.5"
        assert base_name(b"SRR123.5 merged_150_87") == b"SRR123.5"

    def test_merged_name_is_fastp_style(self):
        """Merged read carries fastp's `<id> merged_<n1>_<n2>` token; ZNA base name intact."""
        frag = rand_seq(40, 22)                       # insert 40, L 30 -> n1=30, n2=10
        (h1, s1, q1), (h2, s2, q2) = make_pair(frag, 30, name=b"R.7")
        recs, outcome, _d, _s, _olen, _diff = process_pair(h1, s1, q1, h2, s2, q2, P)
        assert outcome == PairOutcome.MERGED
        name = recs[0][0]
        assert name == b"R.7 merged_30_10"
        # khorana seq.py parse_merged_fastq requires the last token to start with 'merged'
        assert name.rsplit(None, 1)[-1].startswith(b"merged")
        assert base_name(name) == b"R.7"              # ZNA pairing still sees "R.7"

    def test_merged_strips_suffix_preserves_tags(self):
        frag = rand_seq(30, 20)
        h1 = b"SRR.1/1\tRX:Z:ACGT"
        recs, outcome, _d, _s, _olen, _diff = process_pair(
            h1, frag, qual(frag), b"SRR.1/2\tRX:Z:ACGT", rc(frag), qual(frag), P
        )
        assert outcome == PairOutcome.MERGED
        # /1 gone, tag preserved, fastp merged token appended (base name still SRR.1 for ZNA)
        assert recs[0][0] == b"SRR.1\tRX:Z:ACGT merged_30_0"
        assert base_name(recs[0][0]) == b"SRR.1"


# --------------------------------------------------------------------------- #
# property test: tile-or-merge over the whole insert-size range
# --------------------------------------------------------------------------- #

class TestProperty:
    @pytest.mark.parametrize("insert", list(range(30, 61)))
    def test_merge_or_trim_never_double_counts(self, insert):
        L = 30
        frag = rand_seq(60, 999)[:insert]
        (h1, s1, q1), (h2, s2, q2) = make_pair(frag, L)
        recs, outcome, _d, _s, _olen, _diff = process_pair(h1, s1, q1, h2, s2, q2, P)
        overlap = 2 * L - insert
        if overlap >= min_matches(T_MERGE_Q, 0):  # 15 clean bases reach the merge threshold
            assert outcome == PairOutcome.MERGED
            assert recs[0][1] == frag           # reconstructs the molecule
        elif overlap >= min_matches(T_TRIM_Q, 0):  # 5..14 -> trim band
            assert outcome == PairOutcome.TRIMMED
            r1_out, r2_out = recs[0], recs[1]
            assert r1_out[1] + rc(r2_out[1]) == frag       # tile exactly once
        else:                                   # <= 4 bases: not enough evidence
            assert outcome == PairOutcome.KEPT
            assert [r[1] for r in recs] == [s1, s2]


# --------------------------------------------------------------------------- #
# end-to-end CLI (mixed interleaved output + stats JSON)
# --------------------------------------------------------------------------- #

def _write_fastq_gz(path, records):
    with gzip.open(path, "wb") as fh:
        for name, seq in records:
            fh.write(b"@%b\n%b\n+\n%b\n" % (name, seq, b"I" * len(seq)))


class TestCLI:
    def test_end_to_end_mixed_stream(self, tmp_path):
        merge_frag = rand_seq(40, 30)                    # overlap 20 -> merge
        trim_frag = rand_seq(48, 31)                     # overlap 12 -> trim
        disj1, disj2 = rand_seq(50, 32), rand_seq(50, 33)  # -> keep both

        r1s, r2s = [], []
        (h1, s1, _), (h2, s2, _) = make_pair(merge_frag, 30, name=b"M")
        r1s.append((h1, s1)); r2s.append((h2, s2))
        (h1, s1, _), (h2, s2, _) = make_pair(trim_frag, 30, name=b"T")
        r1s.append((h1, s1)); r2s.append((h2, s2))
        r1s.append((b"D/1", disj1)); r2s.append((b"D/2", disj2))

        in1 = tmp_path / "r1.fastq.gz"
        in2 = tmp_path / "r2.fastq.gz"
        out = tmp_path / "merged.fastq"       # plain (no pigz dependency in test)
        js = tmp_path / "stats.json"
        _write_fastq_gz(in1, r1s)
        _write_fastq_gz(in2, r2s)

        args = cli.build_parser().parse_args([
            "--in1", str(in1), "--in2", str(in2), "--out", str(out),
            "--json", str(js), "--threads", "1", "--min-read-length", "10", "-q",
        ])
        stats = cli.run(args)

        assert stats["input_pairs"] == 3
        assert stats["merged"] == 1
        assert stats["trimmed_pairs"] == 1
        assert stats["kept_pairs"] == 1
        assert stats["emitted_records"] == 1 + 2 + 2      # merged(1) + trim(2) + keep(2)
        assert stats["params"]["threshold_merge_bits"] == 28.0
        assert stats["params"]["threshold_trim_bits"] == 8.0

        # the overlap-length histogram counts only pairs where an overlap was found,
        # in the natural quantum (bases), and the merged record's length is the insert
        ohist = stats["overlap_length_histogram"]
        assert ohist == {"20": 1, "12": 1}                # the disjoint pair is absent
        assert stats["insert_size_histogram"] == {"40": 1}
        assert stats["overlap_mismatch_rate"] == 0.0      # clean synthetic overlaps

        # Parse the output stream and check names/pairing structure.
        lines = out.read_bytes().splitlines()
        headers = [lines[i][1:] for i in range(0, len(lines), 4)]
        ids = [h.split()[0] for h in headers]             # base id (drop " merged_.." token)
        assert b"M" in ids                                # merged single, suffix stripped
        assert b"T/1" in ids and b"T/2" in ids            # trimmed pair, adjacent
        assert b"D/1" in ids and b"D/2" in ids            # kept pair
        # the merged read carries fastp's "merged_<n1>_<n2>" name token
        m_hdr = next(h for h in headers if h.split()[0] == b"M")
        assert m_hdr.split(b" ", 1)[1].startswith(b"merged_")
        # merged read reconstructs the molecule
        seqs = {ids[k]: lines[k * 4 + 1] for k in range(len(ids))}
        assert seqs[b"M"] == merge_frag

        # stats JSON round-trips
        assert json.loads(js.read_text())["merged"] == 1

    def test_thresholds_are_settable(self, tmp_path):
        """Raising --threshold-merge turns a merge into a trim; the bands are live."""
        frag = rand_seq(40, 34)                           # overlap 20 -> 39.7 bits
        (h1, s1, _), (h2, s2, _) = make_pair(frag, 30, name=b"M")
        in1, in2 = tmp_path / "r1.fastq.gz", tmp_path / "r2.fastq.gz"
        _write_fastq_gz(in1, [(h1, s1)]); _write_fastq_gz(in2, [(h2, s2)])

        def run(t_merge):
            args = cli.build_parser().parse_args([
                "--in1", str(in1), "--in2", str(in2), "--out", str(tmp_path / "o.fastq"),
                "--threshold-merge", str(t_merge), "--min-read-length", "10", "-q",
            ])
            return cli.run(args)

        assert run(28.0)["merged"] == 1
        assert run(50.0)["trimmed_pairs"] == 1            # 39.7 bits now below merge

    def test_trim_threshold_above_merge_is_rejected(self, tmp_path):
        args = cli.build_parser().parse_args([
            "--in1", "a", "--in2", "b", "--out", "c",
            "--threshold-merge", "10", "--threshold-trim", "20", "-q",
        ])
        with pytest.raises(SystemExit):
            cli.run(args)

    def test_threads_match_single(self, tmp_path):
        """--threads N must emit the same records as one thread -- in fact the same
        BYTES, because chunks are written in submission order."""
        r1s, r2s = [], []
        for i in range(40):
            insert = 30 + (i % 31)            # spans merge / trim / keep bands (L=100)
            frag = rand_seq(260, 500 + i)[:100 + insert]   # varied inserts, read len 100
            (h1, s1, _), (h2, s2, _) = make_pair(frag, 100, name=f"F{i}".encode())
            r1s.append((h1, s1)); r2s.append((h2, s2))
        in1, in2 = tmp_path / "r1.fastq.gz", tmp_path / "r2.fastq.gz"
        _write_fastq_gz(in1, r1s); _write_fastq_gz(in2, r2s)

        def record_set(procs):
            out = tmp_path / f"o{procs}.fastq"
            args = cli.build_parser().parse_args([
                "--in1", str(in1), "--in2", str(in2), "--out", str(out),
                "--threads", str(procs), "--chunk-size", "7", "-q",
            ])
            stats = cli.run(args)
            lines = out.read_bytes().splitlines()
            recs = frozenset(
                (lines[i], lines[i + 1], lines[i + 3]) for i in range(0, len(lines), 4)
            )
            return recs, stats

        single, s_stats = record_set(1)
        par3, p_stats = record_set(3)
        assert single == par3                       # identical record set
        assert s_stats["input_pairs"] == p_stats["input_pairs"] == 40
        for k in ("merged", "trimmed_pairs", "kept_pairs", "emitted_records",
                  "bases_trimmed", "dropped_below_min_length",
                  "fragments_dropped_short_mate", "trim_guard_kept_untrimmed"):
            assert s_stats[k] == p_stats[k], k       # identical aggregate stats
        assert s_stats["overlap_length_histogram"] == p_stats["overlap_length_histogram"]

    def test_short_mate_fragment_dropped_and_logged(self, tmp_path):
        """A pair with one below-min mate is dropped whole, logged, and emits no lone read."""
        good1, good2 = rand_seq(100, 60), rand_seq(100, 61)    # disjoint pair, both long
        smate1, smate2 = rand_seq(100, 62), rand_seq(30, 63)   # R2 30 bp (< 50): drop frag
        in1, in2 = tmp_path / "r1.fastq.gz", tmp_path / "r2.fastq.gz"
        out, js = tmp_path / "o.fastq", tmp_path / "s.json"
        _write_fastq_gz(in1, [(b"G/1", good1), (b"S/1", smate1)])
        _write_fastq_gz(in2, [(b"G/2", good2), (b"S/2", smate2)])
        args = cli.build_parser().parse_args([
            "--in1", str(in1), "--in2", str(in2), "--out", str(out),
            "--json", str(js), "--min-read-length", "50", "-q",
        ])
        stats = cli.run(args)
        assert stats["fragments_dropped_short_mate"] == 1
        assert stats["dropped_below_min_length"] == 2          # both mates of that pair
        assert stats["emitted_records"] == 2                   # only the good pair survives
        ids = [l[1:].split()[0] for l in out.read_bytes().splitlines()[0::4]]
        assert set(ids) == {b"G/1", b"G/2"}                    # no lone S read emitted

    def test_stats_are_identical_across_thread_counts(self, tmp_path):
        """TOTAL equality, not a hardcoded key list.

        All counting happens in one place (the backend's merge_chunk), so this is what
        pins it. An older version enumerated 8 keys and silently stopped covering a
        histogram the moment it was added.
        """
        r1s, r2s = [], []
        for i in range(60):
            frag = rand_seq(300, 700 + i)[:120 + (i % 40)]
            (h1, s1, _), (h2, s2, _) = make_pair(frag, 100, name=f"P{i}".encode())
            r1s.append((h1, s1)); r2s.append((h2, s2))
        in1, in2 = tmp_path / "r1.fastq.gz", tmp_path / "r2.fastq.gz"
        _write_fastq_gz(in1, r1s); _write_fastq_gz(in2, r2s)

        def run(procs, chunk):
            args = cli.build_parser().parse_args([
                "--in1", str(in1), "--in2", str(in2),
                "--out", str(tmp_path / f"o{procs}_{chunk}.fastq"),
                "--threads", str(procs), "--chunk-size", str(chunk),
                "-q"])
            stats = cli.run(args)
            for wallclock in ("elapsed_s", "pairs_per_second"):
                stats.pop(wallclock, None)          # timing, not a merge result
            return stats

        base = run(1, 50000)
        for procs, chunk in ((1, 7), (3, 7), (4, 13), (2, 1000)):
            assert run(procs, chunk) == base, (procs, chunk)

    def test_bases_trimmed_is_right_when_both_mates_share_a_header(self, tmp_path):
        """Without `samtools fastq -N` the two mates carry identical headers. The old
        `next(r for r in records if r[0] != h1)` then returned None and charged the
        whole of R2 as trimmed."""
        frag = rand_seq(48, 4242)                     # overlap 12 -> trim band
        (_, s1, _), (_, s2, _) = make_pair(frag, 30)
        in1, in2 = tmp_path / "r1.fastq.gz", tmp_path / "r2.fastq.gz"
        _write_fastq_gz(in1, [(b"SAME", s1)])         # no /1,/2 suffix: identical names
        _write_fastq_gz(in2, [(b"SAME", s2)])
        args = cli.build_parser().parse_args([
            "--in1", str(in1), "--in2", str(in2), "--out", str(tmp_path / "o.fastq"),
            "--min-read-length", "10", "-q"])
        stats = cli.run(args)
        assert stats["trimmed_pairs"] == 1
        assert stats["bases_trimmed"] == 12           # the overlap, not all 30 of R2

    def test_empty_input_fails_loudly(self, tmp_path):
        """An empty library must not sail through to a 0-record .zna."""
        in1, in2 = tmp_path / "r1.fastq.gz", tmp_path / "r2.fastq.gz"
        _write_fastq_gz(in1, []); _write_fastq_gz(in2, [])
        argv = ["--in1", str(in1), "--in2", str(in2),
                "--out", str(tmp_path / "o.fastq"), "-q"]
        with pytest.raises(SystemExit):
            cli.run(cli.build_parser().parse_args(argv))
        # ...unless it is expected.
        stats = cli.run(cli.build_parser().parse_args(argv + ["--allow-empty"]))
        assert stats["input_pairs"] == 0

    def _run_on(self, tmp_path, r1_bytes, r2_bytes):
        in1, in2 = tmp_path / "r1.fastq", tmp_path / "r2.fastq"
        in1.write_bytes(r1_bytes)
        in2.write_bytes(r2_bytes)
        return cli.run(cli.build_parser().parse_args(
            ["--in1", str(in1), "--in2", str(in2),
             "--out", str(tmp_path / "o.fastq"), "-q"]))

    @staticmethod
    def _records(n, suffix, seed0):
        return b"".join(b"@r%d/%b\n%b\n+\n%b\n" % (i, suffix, rand_seq(40, seed0 + i),
                                                    b"I" * 40) for i in range(n))

    def test_a_short_final_quality_line_is_rejected(self, tmp_path):
        """The truncation the old reader accepted: a short-but-nonempty quality line
        satisfied `if not q`, so the record was emitted malformed with rc=0."""
        good = self._records(3, b"1", 0)
        # drop 12 bases from the last quality line but KEEP its newline, so the record
        # still looks complete and only the length check can catch it
        broken = good[:-13] + b"\n"
        with pytest.raises(SystemExit, match="quality"):
            self._run_on(tmp_path, broken, self._records(3, b"2", 100))

    def test_a_truncated_final_record_is_rejected(self, tmp_path):
        """Cut so the last record loses its final newline: it is no longer a complete
        record, and must be reported as truncation rather than as a count mismatch."""
        good = self._records(3, b"1", 0)
        with pytest.raises(SystemExit, match="truncated"):
            self._run_on(tmp_path, good[:-12], self._records(3, b"2", 100))

    def test_unequal_read_counts_are_rejected(self, tmp_path):
        """Whole extra records on one side -- a different fault, and it must say so.
        This is the one the audit's raw-blob prototype returned success for."""
        with pytest.raises(SystemExit, match="unequal read counts"):
            self._run_on(tmp_path, self._records(4, b"1", 0), self._records(3, b"2", 100))
        with pytest.raises(SystemExit, match="unequal read counts"):
            self._run_on(tmp_path, self._records(3, b"1", 0), self._records(4, b"2", 100))

    @pytest.mark.parametrize("flags", [
        ["--threshold-merge", "0.5", "--threshold-trim", "0.5"],   # merges 94% silently
        ["--min-read-length", "-5"],
        ["--threads", "0"],
        ["--chunk-size", "0"],
        ["--threshold-trim", "40"],                                # trim > merge
    ])
    def test_nonsense_arguments_are_rejected(self, tmp_path, flags):
        in1, in2 = tmp_path / "r1.fastq.gz", tmp_path / "r2.fastq.gz"
        _write_fastq_gz(in1, [(b"A/1", rand_seq(40, 1))])
        _write_fastq_gz(in2, [(b"A/2", rand_seq(40, 2))])
        with pytest.raises(SystemExit):
            cli.run(cli.build_parser().parse_args(
                ["--in1", str(in1), "--in2", str(in2),
                 "--out", str(tmp_path / "o.fastq"), "-q"] + flags))

    def test_sync_check_raises_on_desync(self, tmp_path):
        in1 = tmp_path / "r1.fastq.gz"
        in2 = tmp_path / "r2.fastq.gz"
        out = tmp_path / "o.fastq"
        _write_fastq_gz(in1, [(b"A/1", rand_seq(40, 1))])
        _write_fastq_gz(in2, [(b"B/2", rand_seq(40, 2))])   # mismatched base name
        args = cli.build_parser().parse_args(
            ["--in1", str(in1), "--in2", str(in2), "--out", str(out), "-q"]
        )
        with pytest.raises(SystemExit):
            cli.run(args)


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
        """`run()` is what `zna encode --merge-pairs` will call in-process, and what
        every other test here calls. The guard belongs to the CLI, not to it."""
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
        backend = "auto" if "accel" in available_merge_backends() else "python"
        proc = subprocess.run(
            [sys.executable, "-m", "zna.cli", "merge", "--in1", str(in1),
             "--in2", str(in2), "--out", str(out), "--min-read-length", "10",
             "--backend", backend, "-q"],
            capture_output=True, text=True, timeout=300)
        assert proc.returncode == 0, proc.stderr
        assert out.read_bytes().splitlines()[1] == frag


# --------------------------------------------------------------------------- #
# worker death: must fail, not hang
# --------------------------------------------------------------------------- #

# The process pool is gone: the merge kernel releases the GIL, so workers are threads.
# Two tests went with it, and it is worth recording what they covered rather than
# quietly dropping them:
#
#   * `test_a_killed_worker_fails_the_run_instead_of_hanging` -- mp.Pool could not
#     detect abrupt worker death and blocked forever. There are no worker processes to
#     kill now; a fatal fault takes the whole process down, loudly.
#   * `test_the_parallel_path_survives_a_platform_without_fork` -- there is no fork
#     context to be missing.
#
# What replaced them is stronger than either: output is written in submission order, so
# `test_threads_match_single` compares whole files across thread counts.
