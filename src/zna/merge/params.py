"""Merge parameters, and the fixed-point scale the overlap score is computed in.

The score is a log-likelihood ratio in **bits** (see ``docs/READ_MERGE_REDESIGN.md``),
but it is *computed* in **integers**::

    score_q = n * match_q - d * step_q          # int64, SCALE units per bit

No float appears anywhere in the kernel or in the merge/trim decision, so the argmax is
bit-identical across compilers, optimisation levels and platforms. That is not a
micro-optimisation: this is training data, and a given FASTQ must produce the same
output everywhere. It also removes the one place a float used to decide control flow —
the per-shift mismatch budget, which was ``int((ceiling - best) / step)`` and is now the
exact integer ``(ceiling - best - 1) // step_q``.

**Why SCALE = 2**24.** Quantising the weights perturbs the score by at most
``(n + 2d)/2 * 2**-24`` bits, so a decision could in principle flip where the true score
sits within that of a threshold. Whether it ever does is exhaustively enumerable over
integer ``(n, d)`` -- ``tests/test_merge.py`` carries the enumeration. Measured:

smallest overlap ``n`` at which the fixed-point and float decisions differ at all:

=========  =====================  =====================
SCALE      at ``t_trim`` = 8      at ``t_merge`` = 28
=========  =====================  =====================
2**16      n = 2,718              n = 3,299
2**20      n = 4,377              n = 2,575
**2**24**  **n = 34,632**         **n = 32,830**
2**30      n = 176,375            n = 174,573
=========  =====================  =====================

(The ordering is not monotonic in SCALE because it depends on where the near-threshold
``(n, d)`` pairs happen to fall relative to each rounding; that is not a typo.)

2**20 -- the obvious "millionths of a bit" choice -- disagrees at an overlap of 2,575
bases, which is not comfortably out of reach. 2**24 needs an overlap of **32,830**,
which is two orders of magnitude past any Illumina read and sits on a scan whose O(L^2)
cost makes such input impractical anyway. So no read-length cap is needed: ``int64``
cannot overflow on any input that fits in memory (``n * match_q`` overflows at
n = 2.8e11), and the fixed-point score is exact and deterministic at every length.

**The weights are derived here and only here.** ``log2`` is not correctly rounded and
libm differs between platforms, so computing ``match_q`` on both sides of a Python/C++
boundary would invite a 1-ULP disagreement that silently changes a corpus. One call
site; integers cross the boundary; the values are echoed into the JSON stats so any
output can be audited against the exact integers that produced it.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from math import log2

#: Null hypothesis: two unrelated bases agree 1 time in 4.
P_NULL = 0.25

#: Fixed-point units per bit of log-likelihood ratio. See the module docstring.
SCALE_BITS = 24
SCALE = 1 << SCALE_BITS


def score_weights(err_rate: float):
    """Return ``(match_w, mismatch_w)`` in bits for a per-position error rate.

    ``err_rate`` is the probability that an aligned pair of *real* overlap bases
    disagrees — roughly twice the per-base sequencing error, since either read can be
    wrong. Both weights are log-likelihood ratios against the chance-alignment null
    (bases agree with probability ``P_NULL``); ``mismatch_w`` is returned as a positive
    magnitude and is *subtracted* by the kernel.

    At ``err_rate = 0.01``: ``match_w = 1.9855``, ``mismatch_w = 6.2288``.
    """
    return log2((1.0 - err_rate) / P_NULL), log2((1.0 - P_NULL) / err_rate)


def threshold_bits(read_len: int, alpha: float) -> float:
    """Threshold in bits for a target spurious-detection rate ``alpha``.

    ``T = log2(N / alpha)`` over ``N ~ 2 * read_len`` candidate shifts (union bound;
    ``E[LR] = 1`` under the null gives ``P(score > T) <= 2**-T`` per shift). Not used by
    the kernel — it is how the defaults were derived, kept here so the derivation is
    executable rather than folklore.
    """
    return log2(2.0 * read_len / alpha)


def to_q(bits: float) -> int:
    """Bits -> fixed point. The single rounding rule; do not inline it elsewhere."""
    return round(bits * SCALE)


def to_bits(q: int) -> float:
    """Fixed point -> bits, for humans, logs and the JSON stats. Never for a decision."""
    return q / SCALE


@dataclass
class MergeParams:
    """Decision thresholds, in **bits** of log-likelihood ratio against chance.

    ``t_merge`` and ``t_trim`` are two readings of one calibrated scale: a threshold
    ``T`` tolerates a spurious rate of about ``alpha = N * 2**-T`` over the ``N ~ 2 *
    readlen`` candidate shifts, so ``T = 28`` is ``alpha = 1e-6`` at 2x150. ``err_rate``
    is the probability that two aligned *real* overlap bases disagree (~2x the per-base
    sequencing error); it sets the two score weights and is not a tuning knob.

    The ``*_q`` fields are those same numbers in fixed point, derived once here. The
    kernel and every threshold comparison use them; nothing downstream re-derives them.
    """
    t_merge: float = 28.0        # score >= this -> merge into one full-fragment read
    t_trim: float = 8.0          # [t_trim, t_merge) -> keep both, trim R2's redundant 3'
    err_rate: float = 0.01       # per-position mismatch rate under "real overlap"
    min_read_length: int = 40    # drop emitted reads shorter than this (post merge/trim)

    # Derived, in fixed point. Not constructor arguments, and excluded from equality so
    # two MergeParams with the same bits compare equal.
    match_q: int = field(init=False, repr=False, compare=False, default=0)
    step_q: int = field(init=False, repr=False, compare=False, default=0)
    t_merge_q: int = field(init=False, repr=False, compare=False, default=0)
    t_trim_q: int = field(init=False, repr=False, compare=False, default=0)

    def __post_init__(self) -> None:
        match_w, mismatch_w = score_weights(self.err_rate)
        # Quantise the two weights independently, then add: `step` is the score given up
        # per mismatch relative to an all-match overlap, and it must be exactly the sum
        # of the two quantised weights or `score = (n-d)*match - d*mismatch` and
        # `score = n*match - d*step` stop agreeing.
        self.match_q = to_q(match_w)
        self.step_q = self.match_q + to_q(mismatch_w)
        self.t_merge_q = to_q(self.t_merge)
        self.t_trim_q = to_q(self.t_trim)


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
