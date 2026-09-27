"""Merge parameters: what is chosen, and what is derived from it exactly.

The 0.6 policy (``docs/archive/MERGE_ACCURACY_PLAN.md`` §2) has one statistical tolerance,
``alpha``, one model parameter, the error rate ``e`` (``--error-rate``, default 0.01),
and one declaration, ``adapter_trimmed``. Everything the kernel compares against is
*derived* from those, here, once per run:

========================  ==========================================================
``match_q``, ``step_q``   the score weights, ``log2((1-e)/0.25)`` and
                          ``log2(0.75/e)``, in fixed point (below)
``T_q[N]``                the per-pair merge floor ``log2(N / alpha)``, ``N = len1 +
                          len2 - 1`` candidate shifts: a union bound over the shifts
                          caps a chance merge of unrelated sequence at ``alpha`` per pair
``dfit[n]``               the most mismatches a TRUE overlap of ``n`` compared bases
                          shows with probability ``>= alpha``: ``max d`` with
                          ``P(Binom(n, e) >= d) >= alpha`` -- the plausibility gate
========================  ==========================================================

**The error rate is a parameter, not an estimate** (plan §3). It enters twice, in
different roles. In the *score* it sets strictness only: under unrelated sequence the
likelihood ratio averages exactly 1 whatever ``e`` the overlap hypothesis assumes, so the
``alpha`` bound on chance merges holds at any ``e``. In the *gate* it is a promise: a
true overlap is refused with probability ``<= alpha`` only if ``e`` is at least the
library's real disagreement rate. A per-library estimate was built and measured: fed
into the score, a noisy library's own higher rate softened the mismatch penalty and let
MORE divergent repeats and short false read-throughs through (3'-ramp sets: 560 -> 741
wrong merges at an estimated 0.035, against 560 -> 240 at a fixed 0.01), and it needed a
buffered sample, a prior and a caveat about sorted input. So ``e`` is fixed and
documented, and the run *checks* it instead: :mod:`zna.merge.cli` reports the rate the
detected overlaps actually show, the share of true overlaps the gate is expected to
refuse at that rate, and warns when that share exceeds 0.1%.

**No libm anywhere in a decision.** ``e`` and ``alpha`` are exact rationals, each parsed
from the decimal string the user typed, so the weights and ``T`` are computed with
:mod:`decimal` (software, correctly rounded, 50 significant digits) and ``dfit`` with pure
integers. What crosses into either backend is integers: two weights and two ``int64``
tables. 0.5.3 computed its weights with ``math.log2`` and had to pin the resulting
integers per platform; exact arithmetic lands on the same two at ``e = 0.01``
(33,311,170 and 137,813,407), and the tests pin them.

**The tables are grown, never capped.** Both are indexed by read geometry, so they are
built for a *capacity* -- the longest read they cover -- and doubled whenever a longer
read turns up (:meth:`MergeParams.ensure`). A backend given too small a table stops
before the pair that needs more and says so; the driver grows the tables and resumes
(:mod:`zna.merge.cli`). A table's prefix never changes when it grows, so the output does
not depend on when growth happened. Measured cost: ``dfit`` to 512 at ``e = 0.01`` is
0.5 ms, to 2,048 at ``e = 0.0087`` 13 ms, to 16,384 0.41 s; ``T`` to ``N = 1,024`` is
37 ms and to ``N = 32,768`` about 1.2 s, so one 10 kb read costs ~1.6 s of table
building, once. ``dfit``'s integers grow with the digits of ``e``'s denominator, so a
many-digit ``e`` costs proportionally more, once.

The score and the fixed-point scale
-----------------------------------

The score is a log-likelihood ratio in **bits** (see ``docs/METHODS.md``), but it is
*computed* in **integers**::

    score_q = n * match_q - d * step_q          # int64, SCALE units per bit

No float appears anywhere in the kernel or in the merge decision, so the argmax is
bit-identical across compilers, optimisation levels and platforms. That is not a
micro-optimisation: this is training data, and a given FASTQ must produce the same
output everywhere.

**Why SCALE = 2**24.** Quantising the weights and the floor perturbs the score by at
most ``(n + 2d + 1)/2 * 2**-24`` bits relative to the floor, so the merge decision
``score_q >= T_q[N]`` could in principle differ from the exact real-valued one where the
true score sits within that of ``log2(N / alpha)``. Whether it ever does is exhaustively
enumerable over integer ``(n, d)`` for a given ``(e, alpha, N)``, and
``tests/test_merge.py`` carries the enumeration. With 0.5.3's two fixed thresholds (8
and 28 bits at ``e = 0.01``) the first disagreement was at an overlap of 32,830 bases.
Re-measured for the per-pair floor -- ``alpha = 1e-6``, ``N`` = 99, 199, 299, 599 and
2,047, ``e`` = 0.00016, 0.001, 0.0087, 0.01, 0.03 and 0.1 -- no decision differs at any
overlap up to 10,000 bases. The smallest overlap at which any of those 30 settings
disagrees is n = 10,951 (``e`` = 0.0087, ``N`` = 2,047), and 22 of the 30 show no
disagreement at all up to n = 20,000. That is still two orders of magnitude past any
Illumina read, and it is a statement about the settings enumerated,
not about every rational ``e``: a new ``(e, N)`` could put a near-threshold ``(n, d)``
closer to a rounding boundary. What holds for *every* setting is the bound itself -- a
flip needs the exact score within ``(n + 2d + 1)/2`` units of ``2**-24`` bits of the
floor -- and the plausibility gate is an integer comparison (``d_inf <= dfit[n]``) that
no scale affects. ``int64`` cannot overflow on any input that fits in memory
(``n * match_q`` overflows at n = 2.8e11).

(0.5.3 measured the scale choice at its fixed thresholds: 2**20, the obvious "millionths
of a bit", first disagreed at an overlap of 2,575 bases, which is not comfortably out of
reach; 2**24 at 32,830. The ordering across scales is not monotonic -- it depends on
where the near-threshold ``(n, d)`` fall relative to each rounding.)

**Everything is derived here and only here.** One call site; integers cross the
boundary; the values are echoed into the JSON stats and the ``.zna`` prologue so any
output can be audited against the exact numbers that produced it.
"""
from __future__ import annotations

import re
from array import array
from dataclasses import dataclass, field
from decimal import ROUND_HALF_EVEN, Context, Decimal
from fractions import Fraction

#: Null hypothesis: two unrelated bases agree 1 time in 4.
P_NULL = Fraction(1, 4)

#: Fixed-point units per bit of log-likelihood ratio. See the module docstring.
SCALE_BITS = 24
SCALE = 1 << SCALE_BITS

#: The policy these parameters implement, as written into the ``.zna`` prologue's merge
#: record. Consumers compare it for equality; it changes whenever the decision does.
POLICY = "zna-merge-0.6"

#: ``--alpha``'s default: at most one chance merge of unrelated sequence, and at most one
#: true overlap refused by the plausibility gate, per million pairs.
DEFAULT_ALPHA = "1e-6"

#: ``--error-rate``'s default: 0.5.3's hidden ``err_rate``, now documented. Production
#: RNA libraries measured ~0.009 on 0.5.3's overlap disagreement statistic (hulkrna),
#: just under it; a library whose detected overlaps disagree enough more that the gate
#: is expected to refuse over 0.1% of its true overlaps is warned (plan §3).
DEFAULT_ERROR_RATE = "0.01"

#: The ``--adapter-trimmed`` check scans the input's first this-many pairs, in input
#: order (plan §4). A module constant, not a flag, and indexed by pair number rather than
#: by chunk, so the check is the same at any thread count and chunk size. Tests override
#: it.
READTHROUGH_CHECK_PAIRS = 100_000

#: The smallest values accepted, far past anything meaningful, so that a mistyped
#: exponent is refused rather than left building tables: ``dfit``'s integers grow with
#: the digits of ``e`` and ``alpha`` (``--alpha 1e-5000`` took 5.7 s, ``1e-100000`` did
#: not finish in 40 s). At ``alpha = 1e-300`` a merge already needs a ~500-base perfect
#: overlap; at ``e = 1e-9`` one mismatch costs 29.5 bits.
MIN_ALPHA = Fraction(1, 10 ** 300)
MIN_ERROR_RATE = Fraction(1, 10 ** 9)

#: The smallest table capacity (read length) ever built. 2x150 fits without a regrow.
_MIN_CAPACITY = 256

# 50 significant digits, round-half-even: every `decimal` operation below is correctly
# rounded in this context, so the derived integers do not depend on the platform.
_CTX = Context(prec=50, rounding=ROUND_HALF_EVEN)
_LN2 = _CTX.ln(Decimal(2))


# --------------------------------------------------------------------------- #
# exact arithmetic
# --------------------------------------------------------------------------- #

#: What a parameter typed as text may look like: a plain or scientific ASCII decimal.
#: ``Fraction`` alone would also take ``"1/100"``, ``"1_0e-3"`` and non-ASCII digits,
#: none of which is what the help text, the prologue or a config file means by a number.
_DECIMAL = re.compile(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?",
                      re.ASCII)


def exact(value, name: str = "value") -> Fraction:
    """Parse a user-facing number into an exact :class:`~fractions.Fraction`.

    A string is read as the decimal it spells (``"0.01"`` is exactly 1/100), which is why
    the CLI passes strings; a float is read through its shortest ``repr`` for the same
    reason (``0.01`` means 1/100, not the nearest double). Raises ``ValueError`` naming
    *name* on anything else, including a string that is not a plain or scientific
    decimal (:data:`_DECIMAL`).
    """
    if isinstance(value, Fraction):
        return value
    try:
        if isinstance(value, float):
            return Fraction(repr(value))
        if isinstance(value, (int, Decimal)):
            return Fraction(value)
        text = str(value).strip()
        if not _DECIMAL.fullmatch(text):
            raise ValueError
        d = Decimal(text)
    except (ValueError, ZeroDivisionError, TypeError, OverflowError):
        raise ValueError(f"{name} {value!r} is not a number") from None
    # Checked before the exact conversion, which for 1e-10000000 alone takes 8 s.
    if d and abs(d.adjusted()) > 1000:
        raise ValueError(f"{name} {value!r} is out of range")
    return Fraction(d)


def decimal_str(x: Fraction) -> str:
    """A terminating rational as its plain decimal string (``Fraction(1, 100)`` -> "0.01").

    The prologue and the stats carry ``alpha`` and ``e`` this way so the exact value
    survives a JSON round trip. Raises ``ValueError`` for a non-terminating decimal (or
    one past 50 significant digits); :class:`MergeParams` refuses such a parameter up
    front, so a parameter parsed from a decimal string never reaches that.
    """
    d = _CTX.divide(Decimal(x.numerator), Decimal(x.denominator))
    if Fraction(d) != x:
        raise ValueError(f"{x} has no finite decimal expansion")
    s = format(d.normalize(), "f")
    return s


def log2_exact(x: Fraction) -> Decimal:
    """``log2(x)`` for a positive rational, correctly rounded to 50 digits."""
    ln = _CTX.subtract(_CTX.ln(Decimal(x.numerator)), _CTX.ln(Decimal(x.denominator)))
    return _CTX.divide(ln, _LN2)


def to_q(bits) -> int:
    """Bits -> fixed point, round-half-even. The single rounding rule; do not inline it.

    Takes a :class:`~decimal.Decimal` from :func:`log2_exact` (the derivations) or any
    number a test wants to express in bits; a float is converted exactly.
    """
    if not isinstance(bits, Decimal):
        bits = Decimal(bits) if not isinstance(bits, Fraction) else \
            _CTX.divide(Decimal(bits.numerator), Decimal(bits.denominator))
    return int(_CTX.multiply(bits, Decimal(SCALE)).to_integral_value(
        rounding=ROUND_HALF_EVEN))


def to_bits(q: int) -> float:
    """Fixed point -> bits, for humans, logs and the JSON stats. Never for a decision."""
    return q / SCALE


def score_weights(err_rate):
    """Return ``(match_w, mismatch_w)`` in bits, as floats, for humans.

    ``err_rate`` is the probability that an aligned pair of *real* overlap bases
    disagrees -- roughly twice the per-base sequencing error, since either read can be
    wrong. Both weights are log-likelihood ratios against the chance-alignment null (bases
    agree with probability ``P_NULL``); ``mismatch_w`` is a positive magnitude that the
    kernel *subtracts*. At ``err_rate = 0.01``: ``match_w = 1.9855``, ``mismatch_w =
    6.2288``. The kernel never sees these floats -- see :func:`weights_q`.
    """
    e = exact(err_rate, "error rate")
    return (float(log2_exact((1 - e) / P_NULL)),
            float(log2_exact((1 - P_NULL) / e)))


def weights_q(e: Fraction) -> tuple[int, int]:
    """``(match_q, step_q)``: the score weights in fixed point, exactly.

    The two weights are quantised independently and then added: ``step`` is the score
    given up per mismatch relative to an all-match overlap, and it must be exactly the
    sum of the two quantised weights or ``score = (n-d)*match - d*mismatch`` and ``score
    = n*match - d*step`` stop agreeing.
    """
    match_q = to_q(log2_exact((1 - e) / P_NULL))
    return match_q, match_q + to_q(log2_exact((1 - P_NULL) / e))


def threshold_q(n_shifts: int, alpha) -> int:
    """``T_q = to_q(log2(N / alpha))``: the merge floor for a pair with ``N`` shifts.

    Under unrelated sequence a shift reaches ``T`` bits with probability ``<= 2**-T``
    (``E[LR] = 1``), so a union over the ``N = len1 + len2 - 1`` shifts caps a chance
    merge at ``N * 2**-T = alpha`` per pair. 0.5.3's fixed 28 bits was this at 2x150.
    """
    return to_q(log2_exact(Fraction(n_shifts) / exact(alpha, "alpha")))


def threshold_bits(n_shifts: int, alpha) -> float:
    """:func:`threshold_q` in bits, for humans: 2x50 -> 26.6, 2x150 -> 28.2, 2x300 -> 29.2."""
    return float(log2_exact(Fraction(n_shifts) / exact(alpha, "alpha")))


def binom_cap(nmax: int, e: Fraction, alpha: Fraction, start=None) -> list[int]:
    """``dfit[n]`` for ``n = 0..nmax``: the largest ``d`` with ``P(Binom(n, e) >= d) >=
    alpha``, in pure integers.

    With ``e = a/b``, ``c = b - a`` and ``alpha = A/B``, write the upper tail and the
    point mass over the common denominator ``b^n``::

        U_n(k) = sum_{j >= k} t_n(j),    t_n(j) = C(n, j) a^j c^(n-j)
        P(X_n >= k) >= alpha   <=>   U_n(k) * B >= b^n * A     (no rounding anywhere)

    **One step per n, not one per mismatch.** Adding a trial adds at most one success, so
    ``P(X_{n+1} >= d+2) <= P(X_n >= d+1) < alpha`` and ``P(X_{n+1} >= d) >= P(X_n >= d)
    >= alpha``: ``dfit[n+1]`` is ``d`` or ``d + 1`` where ``d = dfit[n]``, and only
    ``U_{n+1}(d+1)`` needs testing. It follows from the state ``(U_n(d+1), t_n(d))`` by
    ``U_{n+1}(k) = b U_n(k) + a t_n(k-1)`` (condition on the last trial), and the point
    masses move by exact divisions -- both sides of each are the integer ``t`` named::

        t_{n+1}(d)   = t_n(d) * c (n+1) / (n+1-d)
        t_{n+1}(d+1) = t_{n+1}(d) * (n+1-d) a / ((d+1) c)

    The previous form rebuilt each ``n``'s tail from ``b^n`` down, ``O(dfit[n])``
    big-integer steps per ``n``, and cost 6.7x per doubling: one 10 kb read (capacity
    16,384) spent 44 s here. This one is O(1) steps per ``n`` and gives the same table
    (checked against it for 23 ``(e, alpha)`` to capacity 2,048).

    ``start`` is an existing prefix to extend: growth recomputes nothing already known,
    but the state at its last ``n`` is rebuilt once, in ``O(dfit[n])``.
    Ported from the policy study's ``bench/kernel_final/exact_tables.py``, generalised
    from ``alpha = 10**-k`` to a rational ``alpha``.
    """
    if not (0 < e < 1) or not (0 < alpha <= 1):
        raise ValueError("need 0 < e < 1 and 0 < alpha <= 1")
    a, b = e.numerator, e.denominator
    c = b - a
    A, B = alpha.numerator, alpha.denominator
    out = list(start) if start else [0]
    n = len(out) - 1
    d = out[-1]
    # The state at n: t = t_n(d) and u = U_n(d+1) = b^n - sum_{j <= d} t_n(j).
    bn = b ** n
    t = c ** n                                          # t_n(0)
    u = bn - t
    for j in range(d):
        t = t * (n - j) * a // ((j + 1) * c)            # t_n(j+1), exact
        u -= t
    while n < nmax:
        u = b * u + a * t                               # U_{n+1}(d+1)
        t = t * c * (n + 1) // (n + 1 - d)              # t_{n+1}(d)
        n += 1
        bn *= b
        if u * B >= bn * A:                             # dfit[n] = d + 1
            t = t * (n - d) * a // ((d + 1) * c)        # t_n(d+1)
            d += 1
            u -= t                                      # U_n(d+1)
        out.append(d)
    return out


# --------------------------------------------------------------------------- #
# the two tables, shared across MergeParams and grown by doubling
#
# Keyed by the exact rationals they depend on and the capacity. A larger table is built
# by extending a smaller one (its prefix is final), never by resizing it in place: a
# worker thread may be reading the old one through the buffer protocol, and an `array`
# exporting a buffer refuses to resize anyway. Capacities are powers of two from 256, so
# a process holds a handful of tables per (e, alpha).
# --------------------------------------------------------------------------- #

_T_TABLES: dict[tuple[Fraction, int], array] = {}
_DFIT_TABLES: dict[tuple[Fraction, Fraction, int], array] = {}


def _capacity_for(read_len: int) -> int:
    cap = _MIN_CAPACITY
    while cap < read_len:
        cap *= 2
    return cap


def _largest(cache: dict, key: tuple):
    """The largest cached table for *key* (a prefix of every larger one), or None."""
    have = [k for k in cache if k[:-1] == key]
    return cache[max(have, key=lambda k: k[-1])] if have else None


def t_table(alpha: Fraction, cap: int) -> array:
    """``T_q[N]`` for ``N = 0 .. 2*cap - 1`` (entry 0 unused: an empty mate has no shift).

    Exactly ``2 * cap`` long, whatever else is cached, so a table's length IS its
    capacity and a backend's behaviour cannot depend on what other runs in the process
    happened to grow.
    """
    tab = _T_TABLES.get((alpha, cap))
    if tab is None:
        prefix = _largest(_T_TABLES, (alpha,))
        tab = array("q", prefix[:2 * cap] if prefix is not None else [0])
        for n in range(max(1, len(tab)), 2 * cap):
            tab.append(threshold_q(n, alpha))
        _T_TABLES[(alpha, cap)] = tab
    return tab


def dfit_table(e: Fraction, alpha: Fraction, cap: int) -> array:
    """``dfit[n]`` for ``n = 0 .. cap``, exactly ``cap + 1`` long."""
    tab = _DFIT_TABLES.get((e, alpha, cap))
    if tab is None:
        prefix = _largest(_DFIT_TABLES, (e, alpha))
        start = None if prefix is None else list(prefix[:cap + 1])
        tab = array("q", binom_cap(cap, e, alpha, start=start))
        _DFIT_TABLES[(e, alpha, cap)] = tab
    return tab


# --------------------------------------------------------------------------- #
# MergeParams
# --------------------------------------------------------------------------- #

@dataclass
class MergeParams:
    """The merge policy's inputs, and everything derived from them.

    ``alpha`` is the one statistical tolerance: at most ``alpha`` chance merges of
    unrelated sequence per pair, and -- through the plausibility gate at the same
    ``alpha`` -- at most ``alpha`` true overlaps refused per pair. ``error_rate`` is the
    expected fraction of disagreeing positions between the two mates where they truly
    overlap (~2x the per-base sequencing error, since either read can be wrong); see the
    module docstring for its two roles. ``adapter_trimmed`` declares that no read extends
    past its molecule, which makes read-through alignments impossible.

    ``alpha`` and ``error_rate`` accept decimal strings (exact), Fractions, ints and
    floats (read through ``repr``), and are used exactly as given -- nothing is rounded.
    A value with no finite decimal expansion (``Fraction(1, 3)``) is refused, because
    the ``.zna`` prologue records both as decimal strings and must record what was used.

    Derived fields are not constructor arguments and are excluded from equality, so two
    MergeParams with the same inputs compare equal.
    """
    alpha: object = DEFAULT_ALPHA
    error_rate: object = DEFAULT_ERROR_RATE
    adapter_trimmed: bool = False
    min_read_length: int = 40    # drop emitted reads shorter than this
    npolicy: str = "trim3"       # 'trim3' | 'random' | 'keep' -- a surviving N
    rng_seed: int = 42           # seeds --npolicy random; see merge_core.hpp

    alpha_exact: Fraction = field(init=False, repr=False, compare=False, default=None)
    e: Fraction = field(init=False, repr=False, compare=False, default=None)
    match_q: int = field(init=False, repr=False, compare=False, default=0)
    step_q: int = field(init=False, repr=False, compare=False, default=0)
    capacity: int = field(init=False, repr=False, compare=False, default=_MIN_CAPACITY)

    #: ``npolicy`` name -> the integer code both backends take.
    _NPOLICY_CODES = {"keep": 0, "trim3": 1, "random": 2}

    def __post_init__(self) -> None:
        a = exact(self.alpha, "alpha")
        if not 0 < a < 1:
            raise ValueError(f"alpha must be in (0, 1); got {self.alpha!r}")
        if a < MIN_ALPHA:
            raise ValueError(
                f"alpha must be >= 1e-300; got {self.alpha!r}. Even 1e-300 needs a "
                f"~500-base perfect overlap to merge, so a smaller value is a typo")
        e = exact(self.error_rate, "error rate")
        if e <= 0:
            # A mismatch would cost log2(0.75/0) bits: one sequencing error would veto
            # any overlap, and dfit would be 0 everywhere.
            raise ValueError(
                f"error rate must be > 0; got {self.error_rate!r}. It is the fraction of "
                f"positions at which two mates disagree in a TRUE overlap, and no real "
                f"library has none")
        if e < MIN_ERROR_RATE:
            raise ValueError(
                f"error rate must be >= 1e-9; got {self.error_rate!r}. It is the fraction "
                f"of positions at which two mates disagree in a TRUE overlap; no "
                f"sequencer comes within orders of magnitude of that")
        if e >= Fraction(3, 4):
            # The mismatch weight log2(0.75/e) would not be positive: at e >= 0.75 the
            # "overlap" hypothesis agrees no better than chance, so the reads are not
            # behaving like overlapping mates at all.
            raise ValueError(
                f"error rate must be < 0.75; got {self.error_rate!r}. At 0.75 two "
                f"overlapping mates agree no better than unrelated sequence")
        for name, raw, x in (("alpha", self.alpha, a),
                             ("error rate", self.error_rate, e)):
            try:
                decimal_str(x)
            except ValueError:
                raise ValueError(
                    f"{name} {raw!r} has no exact decimal form of at most 50 digits, so "
                    f"the prologue could not record the value used") from None
        self.alpha_exact = a
        self.e = e
        self.match_q, self.step_q = weights_q(e)
        self.capacity = _MIN_CAPACITY

    # ------------------------------------------------------------------ #

    @property
    def npolicy_code(self) -> int:
        try:
            return self._NPOLICY_CODES[self.npolicy]
        except KeyError:
            raise ValueError(
                f"unknown npolicy {self.npolicy!r}; expected 'trim3' or 'random'"
            ) from None

    def ensure(self, read_len: int) -> None:
        """Grow both tables to cover reads up to *read_len* bases (doubling)."""
        if read_len > self.capacity:
            self.capacity = _capacity_for(read_len)

    @property
    def t_table(self) -> array:
        return t_table(self.alpha_exact, self.capacity)

    @property
    def dfit_table(self) -> array:
        return dfit_table(self.e, self.alpha_exact, self.capacity)

    def t_q(self, n_shifts: int) -> int:
        """``T_q[N]``, growing the table if needed."""
        self.ensure((n_shifts + 2) // 2)
        return self.t_table[n_shifts]

    def dfit(self, n: int) -> int:
        """``dfit[n]``, growing the table if needed."""
        self.ensure(n)
        return self.dfit_table[n]

    def kernel_args(self) -> tuple:
        """``(match_q, step_q, t_table, dfit_table, adapter_trimmed)``: the policy as the
        backends take it, integers and two ``int64`` buffers."""
        return (self.match_q, self.step_q, self.t_table, self.dfit_table,
                1 if self.adapter_trimmed else 0)

    def merge_record(self) -> dict:
        """The ``merge`` object of a ``.zna`` prologue (``MERGE_ACCURACY_PLAN.md`` §5).

        ``alpha`` and ``error_rate`` are decimal strings, so the exact values survive
        JSON; everything else is what a consumer needs to decide whether a file was made
        by the policy it expects. It is a function of the parameters alone -- nothing
        here is read from the input -- so the writer can emit it before the first record.
        """
        from .. import __version__
        return {
            "policy": POLICY,
            "zna_version": __version__,
            "alpha": decimal_str(self.alpha_exact),
            "error_rate": decimal_str(self.e),
            "adapter_trimmed": bool(self.adapter_trimmed),
            "min_read_length": int(self.min_read_length),
            "npolicy": self.npolicy,
        }


# --------------------------------------------------------------------------- #
# the overlap-consensus posterior table
# --------------------------------------------------------------------------- #
#
# Where the two mates overlap they are two independent reads of the same physical base,
# so a disagreement means exactly one of them is wrong. Which one is not a judgement
# call -- the sequencer already said, in the two Phred scores. With p_i = 10**(-Q_i/10)
# and "exactly one is wrong":
#
#     P(R1 is the wrong one) = p1(1-p2) / (p1(1-p2) + p2(1-p1))
#
# so the consensus base is the higher-Q call, and the posterior error of that call is
# the expression above -- which is *worse* than the winner's own Q, because a contested
# base is less certain than an uncontested one. Nothing to tune.
#
# **This table is built here, in Python, and passed to whichever backend needs it.** It
# is computed with `pow` and `log10`, and libm differs between platforms, so deriving it
# independently on a C++ side would be a licence for the two implementations to disagree
# by a quality unit on some cell. One source of truth, 64 KiB, once per process.

DISAGREE_Q_DIM = 256


def _build_disagree_table() -> bytes:
    """``DISAGREE_Q[q_win * 256 + q_lose]`` = Phred+33 quality of the winning call.

    Indexed by *raw byte value* so the inner loop is one subscript and no arithmetic.

    The table covers all 256 byte values rather than just the legal Phred+33 range
    (33..126). Quality bytes come from a FASTQ file and nothing upstream guarantees they
    are in range; a 127x127 table would make a malformed byte an out-of-bounds read in
    C++ and an ``IndexError`` (or a silent NUL, for bytes under 33) in Python. Covering
    the whole byte space costs 64 KiB once and makes every input defined and identical
    across backends.
    """
    import math
    tbl = bytearray(DISAGREE_Q_DIM * DISAGREE_Q_DIM)
    for qw in range(DISAGREE_Q_DIM):
        pw = 10.0 ** (-(qw - 33) / 10.0)
        base = qw * DISAGREE_Q_DIM
        for ql in range(DISAGREE_Q_DIM):
            pl = 10.0 ** (-(ql - 33) / 10.0)
            num = pw * (1.0 - pl)
            den = num + pl * (1.0 - pw)
            post = 0.5 if den <= 0.0 else num / den
            post = min(max(post, 1e-10), 0.9999)
            q = int(round(-10.0 * math.log10(post)))
            tbl[base + ql] = min(max(q + 33, 33), 126)
    return bytes(tbl)


#: Built at import. ~5 ms, paid only by code that actually merges -- `zna/cli.py`
#: reaches this package through `args.py`, which imports nothing.
DISAGREE_Q = _build_disagree_table()
