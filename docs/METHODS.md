# zna: Methods

This document describes the algorithms zna implements and the reasoning behind the constants they use: how paired reads are merged, how that scan is made bit-reproducible, what per-record fragment geometry the format carries, and how records are packed, normalized, labeled and shuffled. It is written for two readers — a bioinformatician deciding whether the output can be trusted, and a maintainer who needs to know why a threshold is 28 bits rather than 34. It is not a user guide (the README covers the CLI and the API), not a wire-format specification beyond what the algorithms depend on, and not a record of how the code arrived here: no history, no status, no roadmap. Where prose and code disagree the code is the authority, and load-bearing claims below carry file and line citations so each can be checked against it.

**Measurements.** Unless a passage says otherwise, every measured merge number comes from one ground-truth benchmark (`docs/MERGE_BENCHMARK_RESULTS.md`, run at `710957d`): 1,000,000 read pairs simulated from hg38 at 2×150, fragment lengths uniform on [60, 450], 0.2% per-base sequencing error, a position-degrading NovaSeq quality profile.

The counts that recur below, and how they reconcile: of 1,000,000 input pairs, **582,866**
merged and 417,134 were kept or trimmed. 504 of the merges produced a record below
`--min-read-length 40` and were dropped entirely — every one of them a wrong merge — so
**582,362** merged records and **834,268** mate records are emitted, **1,416,630** in
total. Kernel timings and compression ratios come from separate measurements, described
where they appear.

## Contents

- [1. The overlap score and the merge/trim decision](#1-the-overlap-score-and-the-mergetrim-decision)
  - [1.1 One unknown](#11-one-unknown) · [1.2 The score](#12-the-score) · [1.3 Two thresholds on one calibrated scale](#13-two-thresholds-on-one-calibrated-scale) · [1.4 The decision](#14-the-decision) · [1.5 Merge](#15-merge) · [1.6 Trim, symmetric](#16-trim-symmetric) · [1.7 The overlap consensus](#17-the-overlap-consensus) · [1.8 Pair integrity](#18-pair-integrity)
- [2. The kernel: determinism, fixed point, and the search](#2-the-kernel-determinism-fixed-point-and-the-search)
  - [2.1 Scores are integers](#21-scores-are-integers) · [2.2 The argmax is a total order](#22-the-argmax-is-a-total-order) · [2.3 Pruning is exact in integers](#23-pruning-is-exact-in-integers) · [2.4 The inner loop compares raw bytes](#24-the-inner-loop-compares-raw-bytes) · [2.5 Two backends, one specification](#25-two-backends-one-specification) · [2.6 Layering and the consumption protocol](#26-layering-and-the-consumption-protocol) · [2.7 Deterministic output under threading](#27-deterministic-output-under-threading)
- [3. Fragment geometry: what the flags mean](#3-fragment-geometry-what-the-flags-mean)
  - [3.1 The invariant](#31-the-invariant) · [3.2 The flag byte](#32-the-flag-byte) · [3.3 What `IS_RC` records](#33-what-is_rc-records) · [3.4 `IS_FULL_FRAGMENT`](#34-is_full_fragment) · [3.5 The projection: read views and copy paths](#35-the-projection-read-views-and-copy-paths) · [3.6 Where records land](#36-where-records-land)
- [4. The container and the codec](#4-the-container-and-the-codec)
  - [4.1 2-bit packing](#41-2-bit-packing) · [4.2 Blocks and the columnar payload](#42-blocks-and-the-columnar-payload) · [4.3 Strand normalization](#43-strand-normalization) · [4.4 Labels](#44-labels) · [4.5 The shuffle](#45-the-shuffle)

---

## 1. The overlap score and the merge/trim decision

### 1.1 One unknown

Write R1 and revcomp(R2) on a single axis. The only quantity that is not observed is the fragment length `L`, and it is equivalent to the signed offset `s` of revcomp(R2) relative to R1:

```
s = L - len2
```

uniformly, for every geometry: `s > 0` is a normal partial overlap, `s = 0` is exact full overlap, and `s < 0` is read-through, where each mate runs past the far end of the fragment into adapter. There is one axis and one scan, over `s ∈ [-(len2-1), len1-1]`; each `s` fixes the compared region as the intersection of the two reads, so mates of unequal length (after quality trimming) need no special case.

That is the whole geometry. Everything downstream is a function of `s` alone — the merged record is the span `[0, L)` with `L = s + len2`, which mate supplies which base follows from where `len1` falls inside that span, and the redundant overlap is `len1 + len2 - L`. No per-direction case analysis is needed anywhere, and none is used.

### 1.2 The score

For a candidate `s` the scan sees `n` aligned base comparisons, `m` of them matching and `d = n - m` mismatching, and weighs two hypotheses:

- **H₁, real overlap.** A disagreement is a sequencing error, with probability `e` — `err_rate`, default 0.01. It is the probability that two aligned *real* overlap bases disagree, roughly twice the per-base error rate, since either read can be wrong.
- **H₀, chance alignment.** Unrelated bases agree one time in four (`P_NULL = 0.25`).

The score is the log-likelihood ratio in bits:

```
score(s) = m · log₂((1-e)/0.25)  -  d · log₂(0.75/e)
         = m · match_w  -  d · mismatch_w
```

At `e = 0.01`, `match_w = 1.9855` and `mismatch_w = 6.2288` bits (`src/zna/merge/params.py:60-71`). Neither is a tuning knob: both fall out of `e`, which is itself auditable against the data — the JSON stats report mismatches per aligned overlap base over every detected overlap, ~0.009 on production data against the assumed 0.01.

A matching base is worth just under `log₂ 4 = 2` bits, the information in two reads agreeing on one of four bases, and the shortfall is load-bearing rather than cosmetic: at `match_w = 1.9855`, four perfect matches score 7.94 bits and five score 9.93, which is what puts the shortest trimmable overlap at 5 bases rather than 4, and the shortest mergeable overlap at `ceil(28 / 1.9855) = 15` rather than 14. Both minima are pinned by `tests/test_merge.py:191-223`, and the 15-base cliff is visible in the benchmark (§1.3).

Comparison is raw byte equality on upper-cased sequence: `N` agrees with `N`, and IUPAC codes agree with themselves. There is no separate ambiguity case in the kernel, so no fast path can disagree with the reference implementation on any input; §2.4 explains why byte comparison is the specification rather than an approximation of it.

The score is computed in integers, never in floating point: `int64` fixed point at 2⁻²⁴ bits per unit, with both weights derived once in Python and crossed into C++ as integers. No float takes part in a comparison, a bail bound or the argmax, so a given FASTQ produces the same output on every compiler, optimization level and platform — this is training data, and the corpus must not depend on the machine. §2.1 gives the arithmetic and the exhaustive check that quantization never changes a decision at any reachable read length.

**Quality scores do not enter the score.** The scan takes two scalar integer weights derived from a single `err_rate` (`params.py:60-71`; scan signature `merge_core.hpp:125-126`) against a fixed null, `P_NULL = 0.25` (`params.py:53`). There is no per-position (Q1, Q2) score table and no composition-aware null. Quality is used only in the consensus (§1.7), through a precomputed posterior table — which is what removes the need for quality cutoffs anywhere in the decision. One consequence of the fixed null is worth naming: there is no low-complexity correction, so a polyA tail clears the trim threshold but does not on its own reach the merge threshold. That is pinned behaviour, not an accident.

### 1.3 Two thresholds on one calibrated scale

Under H₀ any likelihood ratio has `E[LR] = 1`, so `P(score > T) ≤ 2⁻ᵀ` per shift. Over `N ≈ 2·readlen` candidate shifts a union bound gives a spurious detection rate of about `α = N · 2⁻ᵀ`, i.e.

```
T = log₂(N / α)
```

so `T` is not chosen, it is read off a tolerated error rate: at 2×150, `T = 28` is `α ≈ 1×10⁻⁶`, and each additional bit halves it. Evidence accumulates linearly while the false-positive rate falls exponentially, so being conservative is cheap — tightening `α` from 10⁻³ to 10⁻⁹ costs ten extra matching bases.

Concretely, the matches required to clear a threshold (pinned by the test suite so the arithmetic cannot drift):

| mismatches in the overlap | matches needed at T = 28 |
|---|---|
| 0 | 15 |
| 1 | 18 |
| 2 | 21 |
| 3 | 24 |

Each mismatch costs `mismatch_w / match_w = 3.14` matching bases.

The resulting sensitivity cliff is exactly where the arithmetic puts it. On the benchmark, merges go from ~1% below 15 bases of true overlap to 92.55% in the 15–19 bin and 99.76% at 20–29. The shortfall in the first bin is arithmetic rather than slack: a 15-base overlap scores 29.8 bits and a single mismatch costs 8.21, so any error in a minimal overlap puts it under threshold. Sensitivity over all true overlaps ≥ 15 is 99.83%.

**Read the guarantee for what it is.** `α` bounds detection against *chance alignment of unrelated sequence*, and there it holds with room to spare: 0 spurious merges in 40,000 uniform-random pairs at `T = 28`, at every read length from 50 to 300, and 0 in the 20,000-pair 2×150 case the test suite pins. It is not a bound against real sequence, where two fragment ends can be genuinely homologous. On the benchmark, 1.23% of pairs with no true overlap merge at `T = 28`, and 1,403 wrong merges per million pairs remain at `T = 100` — every one a fragment whose two ends really are homologous, with the scan choosing the correct argmax over a region that genuinely matches. Raising `T` buys far less than the formula suggests, because past a point the score is measuring real similarity rather than chance. §3.4 gives the same population from the consumer's side: what `IS_FULL_FRAGMENT` does and does not assert.

**The two defaults.**

`T_merge = 28` minimizes false positives plus false negatives on the benchmark: 6,603 errors per million pairs, rising monotonically to 13,584 at 34, 45,616 at 60 and 96,512 at 100. The curve is steeply asymmetric because sensitivity falls much faster than the chimera rate does — going from 28 to 34 prevents 708 wrong merges and causes 7,689 missed ones, an exchange rate of about 11:1, reaching 22.5:1 at 100 bits. Raising the threshold therefore pays only if a chimera costs more than eleven missed merges, and a missed merge is not lost data: the pair is still emitted, correctly bounded, with its redundant overlap trimmed. The scale also prices a familiar knob — `--threshold-merge 60` is `ceil(60/1.9855) = 31` clean bases, and lands within a rounding error of fastp's default operating point (0.597% chimera rate at 92.63% sensitivity, against 0.621% and 92.98%), so fastp's 30-base floor is worth about 60 bits of evidence.

`T_trim = 8` sits far lower because the errors cost different amounts. A wrong merge invents sequence — a chimera with a false 3' boundary. A wrong trim deletes real bases from a read tail, which loses signal but invents none. Measured over the same 1M pairs, the trim band removes 255,764 bases of genuinely duplicated sequence and deletes 57,537 real ones, a benefit ratio of 4.45:1; a wrong trim removes a median of 9 bases, mean 20, up to 110. Sweeping 8/12/16/20 raises the ratio and shrinks the benefit, so the decision needs the two harms added — duplicated-left plus deleted is 85,611 / 97,053 / 127,398 / 167,445 — and **8 minimizes the sum**. It stops being optimal only if a deleted base is judged more than 1.63× as harmful as a duplicated one. The band's edges show in the exact-trim curve: nothing below 5 bases of true overlap can clear 8 bits, 96.0% of kept pairs at 5–9 and 99.8% at 10–14 are trimmed to exactly their true overlap, and above ~15 the pair merges instead.

Both thresholds are validated at startup: `0 < T_trim ≤ T_merge`, and `T_merge` must be at least `log₂(2·50 / 0.05) ≈ 11.0` bits, the point at which chance alignments start passing even for 2×50 reads.

### 1.4 The decision

```
for every shift s:  compute score(s)
s* = argmax score(s)

score(s*) ≥ T_merge           ->  merge at s*, one full-fragment record
T_trim ≤ score(s*) < T_merge  ->  keep both reads, split the redundant overlap
score(s*) < T_trim            ->  keep both reads unchanged
```

`argmax`, not first-accept: the scan cannot be captured by a spurious short hit before it reaches the real offset. A chance 4-base hit scores 7.9 bits; a real 40-base overlap scores 79 and wins outright. This is what lets `T_trim` be permissive without costing merges — aggressive trimming and reliable merging do not trade off against each other.

The argmax is a specified total order — maximize score, then minimize `s` — rather than whatever the loops happen to produce; §2.2 gives the order, the visiting sequence that realizes it, and why periodic sequence makes it necessary. Ties are settled by that order and not reported: no margin is computed. The JSON carries an overlap-length histogram rather than a score histogram, because integer bins alias against the 1.9855-bit quantum (`src/zna/merge/cli.py:338-344`).

`T_trim` doubles as the scan's floor. The incumbent starts at `T_trim - 1` (`merge_core.hpp:390-391`, `_pymerge.py:301`), so a best shift below the trim threshold is reported as no overlap at all (`overlap_len == 0`) rather than as a sub-threshold argmax — the same decision, reached without the scan ever returning a shift it does not believe in. The score is monotone in `d` alone, which additionally gives every shift an exact integer mismatch budget to bail on (§2.3).

### 1.5 Merge

`L = s + len2` is the fragment. R1 supplies `[0, min(len1, L))` from its 5' end — its copy of the overlap, carrying the consensus — and revcomp(R2) supplies whatever remains, which is by construction outside the overlap. The record has length `L` identically for every geometry, including read-through, where the adapter falls outside `[0, L)` and is discarded without ever being matched against an adapter sequence. A merged read is a genuine full molecule, so it is length-filtered on its own (§1.8) and carries both fragment ends (§3.4).

### 1.6 Trim, symmetric

The overlap sits at the 3' end of *both* mates — each read starts at a fragment end and reads inward, so they converge in the middle. The pair must tile the fragment exactly once, which forces `keep1 + keep2 = L`; the only freedom is where the cut falls, and `balanced_split` puts it in the middle (`merge_core.hpp:344-351`, `_pymerge.py:217-233`):

```
keep1 = clamp(ceil(L/2), L - len2, len1)
keep2 = L - keep1
```

For equal-length mates the clamp never binds and this is exactly "cut `olen/2` from each", with the odd base staying on R1. Splitting rather than taking the whole overlap off R2 leaves the two emitted reads the same length, which is what downstream aligners and models expect, and discards the *last* cycles of both reads — the lowest-quality bases in the pair — instead of one read's entire copy. Every decision and every accuracy metric is identical under the two rules: same merges, same trims, same duplicated and deleted bases, 0 boundary violations. What changes is the shape of the output — mean `|len(R1) - len(R2)|` on trimmed pairs falls from the whole overlap (10.8) to **0.51**, maximum 1.

Trimming requires `s ≥ 0`. In a read-through geometry the redundant bases are R2's *5'* copy of the fragment and its 3' end is adapter, so there is nothing sensible to cut; such a pair is either merged or kept whole.

**The guard** (`trim_is_allowed`, `merge_core.hpp:372-374`, mirrored at `_pymerge.py:314`): each mate must reach at least `min_read_length` past the other's 3' end,

```
(L - len1) ≥ min_read_length   and   (L - len2) ≥ min_read_length
```

which does two things at once.

1. It bounds the emitted lengths. The split clamps `keep1` into `[L - len2, len1]`, so `keep1 ≥ L - len2` and `keep2 ≥ L - len1`; requiring both of those to clear the length filter puts both emitted reads above it by construction, and a trimmed pair can never be dropped for being short. A trim that would fail the test keeps both reads *untrimmed* instead of turning a usable pair into a discarded fragment, and the JSON counts how often that happens.
2. It caps the overlap the trim band may act on. For equal-length mates the condition is exactly `len - olen ≥ min_read_length`, i.e. an inferred overlap covering nearly the whole read is refused. Such an overlap is spurious almost by definition: 145 clean bases score 288 bits and would have merged outright, so reaching the 8–28 bit band at that length takes a pile of mismatches. Dropping the rule admits exactly those pairs, at a measured cost of 133 extra false trims, 17,214 deleted bases and 9 corrupted 5' ends per million pairs — all of them fragments with no true overlap at all.

Trimming only ever removes bases from a 3' end, so base 0 of every emitted read stays a true fragment boundary. That is the invariant the whole format rests on, and §3.1 gives it and its measurement against simulated truth.

### 1.7 The overlap consensus

Where the mates overlap they are two independent reads of the same physical base, so a disagreement means exactly one of them is wrong. Which one is not a judgement call — the sequencer already said, in the two Phred scores. With `pᵢ = 10^(-Qᵢ/10)` and "exactly one is wrong",

```
P(R1 is the wrong one) = p₁(1-p₂) / (p₁(1-p₂) + p₂(1-p₁))
```

so the consensus base is the higher-Q call, and the posterior error of *that* call is the expression above — which is worse than the winner's own Q, because a contested base is less certain than an uncontested one. For well-separated qualities the derated value is close to `Q_win - Q_lose`: a Q40 base contested by a Q20 base is emitted at ~Q20, not Q40. Nothing to tune, and no quality cutoffs.

The posterior is precomputed as a **256×256 byte table indexed by raw quality byte** (`params.py:152-185`), built once in Python and passed to whichever backend runs, so the two implementations cannot disagree by a quality unit on some cell. The table is 256×256 rather than 94×94 (the legal Phred+33 range) because quality bytes come from a FASTQ file and nothing upstream guarantees they are legal: a 94×94 table would make a malformed byte an out-of-bounds read in C++ and an `IndexError` in Python, where a full table makes it defined. Both backends reject a table that is not exactly 65,536 bytes (`src/zna/merge/_accel.cpp:62`).

Measured against ground truth on the merged records, this recovers **90.35%** of the overlap errors that were recoverable at all (a position where both mates are wrong is not), reducing 182,184 errors under an R1-wins rule to 84,469 against an oracle floor of 74,037, and leaving 895 of 577,275 records worse than doing nothing — a ratio of 109:1.

**N rescue.** An `N` is a no-call: it carries no base information, so a real call on the
other mate beats it whatever the two quality scores say, and the rescued base keeps the
surviving mate's own quality rather than a contested-base derating — there was no
contest. Only `N` is rescued; the IUPAC ambiguity codes do carry partial information and
are left alone. Rescue does not touch the *scan* — an N still counts as a mismatch when a
shift is scored — so which shift wins, and therefore whether the pair merges, is
unchanged.

**Where the consensus is written is part of the decision, not an optimization.** It goes
into the records whose *construction depends on the overlap being real*, and nowhere
else:

| outcome | copies of the overlap emitted | corrected |
|---|---|---|
| merged | R1's only — R2 contributes outside it | R1 |
| trimmed | both, each mate keeping about half | both |
| kept | both, in full | **neither** |

The first two follow from what is emitted. On a merge, R2's copy is never emitted, so
writing it would be dead work. On a trim, R2 keeps roughly half the overlap, so writing
only R1 would hand back part of the consensus exactly where the mates disagreed.

The third is measured rather than deduced. A detection that lands in *kept* is spurious
almost by construction, for two independent reasons:

- at `s ≥ 0` the pair is kept only because `trim_is_allowed` refused it, which requires an
  inferred overlap of more than about 110 bases — and a genuine overlap that long scores
  ~218 bits and would have merged at 28. This is the trim guard's own argument (§1.6),
  and it applies verbatim: an overlap too suspect to *cut* on is too suspect to *rewrite
  bases* on;
- at `s < 0` a genuine read-through overlap equals the fragment length, so scoring inside
  [8, 28) needs a 5–14 bp fragment.

Over the 1M-pair benchmark, of the 3,068 kept pairs carrying a detected overlap **none
found the true shift** and 97.3% had no true overlap at all. Wrong emitted bases in the
overlap window under the three candidate rules:

| rule | wrong bases |
|---|---:|
| correct neither *(shipped)* | **208** |
| correct R1 only | 1,509 |
| correct both | 17,870 |

Correcting R1 alone turned 1,379 correct bases wrong to fix 78.

**A caution for anyone re-deriving this:** `s ≥ 0` is *not* a safe proxy for "the write
stays clear of R2's 5' end". The real bound is

```
min written R2 index = max(0, len2 - len1 + s)
```

which for unequal-length mates reaches index 0 at non-negative shifts. Of the 237
corrupted 5' ends measured when R2 was written unconditionally, **9 are at `s ≥ 0`** —
they are not all read-through, and it is the trim guard rather than the sign of `s` that
keeps trims away from that boundary today.

### 1.8 Pair integrity

`min_read_length` is an orthogonal QC filter applied to the emitted records, after any merge or trim. It is applied **all-or-nothing to an unmerged pair**: if either mate falls below the floor, both are dropped. A lone surviving mate would be encoded as a spurious "single" — a full molecule with two true endpoints — which is precisely the supervision signal the corpus exists to carry (§3), so a half-pair emitted alone is not merely a lost read but a fabricated fragment. A merged record is a genuine full molecule and is therefore filtered on its own.

---

## 2. The kernel: determinism, fixed point, and the search

The scan answers exactly one question — where does revcomp(R2) sit on R1's axis — by scoring every candidate offset `s` and taking the best. Two properties make that answer *reproducible* rather than merely repeatable: the score is an integer, and the argmax is a specified total order. Most of what follows falls out of those two.

### 2.1 Scores are integers

A shift's score is a log-likelihood ratio in bits (§1.2), but it is computed as

```
score_q = n * match_q - d * step_q          # int64, SCALE units per bit
```

with `n` the overlap length, `d` the mismatch count, and `SCALE = 2**24`. No float takes part in a comparison, in the per-shift bail bound, or in the merge/trim decision, so the argmax cannot depend on the compiler, the optimization level or `-ffast-math`. That is what turns "the two implementations agree" from an approximate claim into an exact one: they can be asserted equal shift for shift and byte for byte, because neither of them rounds anything.

The weights are derived **once, in Python** (`src/zna/merge/params.py`):

```
match_q = round(log2((1 - e) / 0.25) * SCALE)      #  33_311_170 at e = 0.01
step_q  = match_q + round(log2(0.75 / e) * SCALE)  # 137_813_407
```

and crossed as integers. `log2` is not correctly rounded and libm differs between platforms, so deriving `match_q` on the C++ side as well would invite a 1-ULP disagreement that silently changes a corpus. One call site; a test pins both constants, so a platform where they differ fails loudly rather than producing a different corpus from the same FASTQ; and the exact integers used are echoed into the JSON stats (`score_scale`, `match_q`, `step_q`, `threshold_merge_q`, `threshold_trim_q`) so any output can be audited against them. The 256×256 consensus posterior table (§1.7) is passed in for the same reason — it is built with `pow`/`log10`, and both backends reject a table of the wrong size.

`step_q` is quantized as the *sum* of the two quantized weights, not independently, so that `n*match_q - d*step_q` and `(n-d)*match_q - d*mismatch_q` agree exactly. `int64` is required rather than stylistic: `match_q ≈ 2**25`, so `n * match_q` leaves `int32` at n = 64, and a 150 bp overlap needs 33 bits.

**Why `SCALE = 2**24`.** Quantizing perturbs the score by at most `(n + 2d)/2 * 2**-24` bits, so a decision could in principle flip where the true score sits within that of a threshold. Whether it ever does is exhaustively enumerable over integer `(n, d)`, and `tests/test_merge.py` carries the enumeration. Smallest overlap `n` at which the fixed-point and float decisions differ at all:

| SCALE | at `t_trim` = 8 | at `t_merge` = 28 |
|---|---:|---:|
| 2^16 | 2,718 | 3,299 |
| 2^20 | 4,377 | 2,575 |
| **2^24** | **34,632** | **32,830** |
| 2^30 | 176,375 | 174,573 |

(The ordering is not monotonic in SCALE: it depends on where the near-threshold `(n, d)` pairs fall relative to each rounding.) 2^20 — the obvious "millionths of a bit" choice — disagrees at an overlap of 2,575 bases, which is not comfortably out of reach. 2^24 needs an overlap of 32,830, two orders of magnitude past any Illumina read. So no read-length cap is needed anywhere in the tool: fixed point is the specification and is exact at every length, and `int64` overflows only past 2.8×10¹¹ bases.

### 2.2 The argmax is a total order

> **Maximize `score`; among ties, minimize `s`** — the leftmost, most read-through alignment.

This is a specification, not an artifact of iteration order (`merge_core.hpp:20-27`, `_pymerge.py:92-96`), and it is what lets a rewritten kernel be tested for exact equality instead of the weaker "returned *an* argmax". Two facts make the two keys sufficient. Ties across *different* overlap lengths need `Δn * match_q == Δd * step_q`, whose minimal solution is `step_q / gcd(match_q, step_q)`; that gcd is 1 for the shipped weights and for every `e` in {0.001, 0.005, 0.01, 0.02, 0.05}, so `Δn` would have to exceed 1.3×10⁸. Ties therefore only ever occur at equal `n`, and `n` never discriminates.

The visiting order realizes the rule directly. The scan walks the plateau first — every shift achieving the maximal overlap `nmax`, ascending in `s` — and then the two flanks in lockstep at decreasing `n`, read-through side (`s = n - len2`) before the normal side (`s = len1 - n`). At any fixed `n` that is ascending in `s`, and improvement is strict `>`, so the first-visited member of a tie keeps the win. **Do not reorder these loops.** Swapping the two flanks changes the winner on every tied pair — and with it which bases a merged read is built from — while leaving every random-sequence test green: random sequence essentially never ties (a sweep of 7,000 random and adversarial pairs produced exactly zero), so the tie fixtures have to be constructed deliberately.

### 2.3 Pruning is exact in integers

The score is monotone decreasing in `d` alone, which gives two bounds.

- **Outer.** Shifts are visited in decreasing `n`, so once `n * match_q <= best` no remaining shift can win and the whole scan stops.
- **Inner.** `score_q > best` iff `d * step_q < ceiling - best` iff `d <= (ceiling - best - 1) / step_q`, where `ceiling = n * match_q`. That budget is computed once per shift and the comparison loop bails the moment it is exceeded.

Both sides compute the same integer because `ceiling > best` is checked first, so the numerator is non-negative and `step_q` is positive — the range where C++ truncating division and Python floor division agree. Computing the budget as `int((ceiling - best) / step)` in floating point would put a platform-dependent boundary back into control flow, which is the one place §2.1's guarantee cannot survive it.

Pruning is consistent with the tie rule for free. Rejecting on `ceiling <= best` discards only shifts that could at best *tie* the incumbent, and a tie always loses under strict `>`. The incumbent's initial value is `T_trim - 1`, which is what makes a sub-threshold best shift come back as "no overlap" (§1.4).

### 2.4 The inner loop compares raw bytes

```c
const uint8_t* a = s1   + (s > 0 ?  s : 0);      // a shift is just a pointer offset
const uint8_t* b = s2rc + (s < 0 ? -s : 0);
for (; k + 32 <= n; k += 32) {
    d += neq16(a + k, b + k) + neq16(a + k + 16, b + k + 16);
    if (d > dmax) return REJECT;                 // bail every 32 bases
}
for (; k < n; ++k) d += (a[k] != b[k]);          // scalar tail
```

`neq16` counts differing bytes in a 16-byte window: `vceqq_u8` + `vaddvq_u8` on NEON, `_mm_cmpeq_epi8` + `movemask` + popcount on SSE2. Sixteen bytes is baseline on every target — NEON on aarch64, SSE2 as part of the x86-64 base ISA — so there is no `-march` flag and no runtime dispatch, and where neither is available the vector block compiles out and the scalar tail does the whole overlap for the same result. Because a shift is a pointer offset there is no packing step, no cross-word realignment, no guard word and no lookup table, and the guard `k + 32 <= n` is the only thing keeping the loads inside the record.

**Byte comparison is not an approximation of the reference semantics — it *is* the reference semantics.** `N` vs `N` earns a full match, `N` vs `A` is a mismatch, IUPAC codes compare as themselves (`reverse_complement` passes anything outside ACGTN through uncomplemented, deliberately), and any byte at all is defined. There is no purity test, no dispatch, no second code path to keep in sync, and no class of input on which the fast path and the oracle can disagree. A 2-bit packed kernel cannot have that property — `N` would collapse onto a real base and start manufacturing evidence — and would need a packer, bit realignment and a purity dispatch to fake it. (The container's 2-bit alphabet, §4.1, is a separate matter: the merger never packs, and `N` policy is applied at encode time.)

Measured over 50,000 real pairs, full pruned scans, every variant checked against the shipped kernel's `(s, n, d)` on every pair, zero mismatches throughout (Apple M-series, aarch64/NEON, Apple clang, `-O3 -std=c++17`, min-of-N):

| variant | µs/pair |
|---|---:|
| numba JIT (pure-Python-source baseline) | 2.633 |
| scalar C++, direct transliteration | 1.075 |
| 2-bit packed, SWAR packer, fused flanks | 0.535 |
| byte SIMD, 16 B vectors, bail every 16 | 0.534 |
| **byte SIMD, 16 B vectors, bail every 32** | **0.470** |
| byte SIMD, 16 B vectors, bail every 64 | 0.515 |
| byte SIMD, bail 32, fused flanks | 0.555 |

Bail granularity, not vector width, is the variable that matters: the scan is rejection-dominated (roughly 0.32·n comparisons per rejected shift), so 32 bases beats both 16 (too much branching) and 64 (too much work past the point of rejection). Fusing the two flank shifts helps a packed kernel and hurts this one — it doubles register pressure and defeats the early-exit ordering.

The reference backend runs the same search with mismatches accumulated branchlessly in blocks of 8 and the bail tested once per block; the result is bit-for-bit identical either way, since overshooting `dmax` inside a block still rejects. A stride bug there (`k += 7`, a wrong `lim`) changes 6.34% of scores and 0.26% of merge/trim/keep decisions on real data while passing every test that uses clean overlaps, so both loops have a dedicated test that sweeps a single mismatch across every position of a block, across the vector boundary and into the tail.

### 2.5 Two backends, one specification

`src/zna/merge/_pymerge.py` is the **reference oracle, not a fallback**. It is plain Python — no JIT, nothing between the reader and the algorithm — about 50× slower than the compiled backend, and that is the correct trade: speed there would only buy the ability to be wrong in the same way as the thing it checks. It is never deleted and never optimized at the cost of clarity, because it is what the accelerated backend is *defined* to agree with, so it has to stay readable enough to be believed.

`src/zna/merge/backend.py` selects between them by name, validating that a backend exposes the full required set (`scan`, `process_pair`, `merge_chunk`, `split_records`) before returning it, and resolving on first use so importing the package costs nothing. `src/zna/merge/overlap.py` is a pure dispatcher over it. The CLI's `--backend {auto,accel,python}` makes either side first-class, and `auto` **refuses to run on the reference kernel**: a silently correct 50× slowdown at cluster scale is indistinguishable from a slow node and burns the allocation before anyone looks. The library entry point does not refuse, and the JSON stats record which kernel actually ran.

What the arrangement buys is enforced by differential tests at all three levels: identical `(shift, score_q, olen, diff)` from `scan`; identical record tuples and counters from `process_pair`; and from `merge_chunk`, an identical blob, identical consumed byte counts, identical counters and identical histograms. Those are equality assertions rather than "an argmax" assertions, which is only possible because of the integer score and the specified order. Two details keep the suite from being circular or blind: tests that go through the public API are parametrized over *every* available backend, since otherwise the oracle goes unchecked in exactly the environment that ships; and the tie fixtures build both plateau ties (periodic content on unequal-length mates) and flank ties (equal-length periodic mates held mutually out of phase), the latter being the only construction that separates the two flanks' visiting order. The differential catches defects of a kind review reliably misses — a compiled backend that threw on reads over 1024 bp where the reference merged them, and a merged-name buffer sized from the read arena, which overflowed the heap on any FASTQ whose headers outran its reads.

### 2.6 Layering and the consumption protocol

Three layers, and the boundaries are load-bearing.

- **Core** (`merge_core.hpp`) — `scan` and `process_pair` over `const uint8_t*` + length, writing into a caller-owned scratch arena. No I/O, no Python objects, no FASTQ. It is header-only and dependency-free so it can also be compiled into a sanitizer driver that cannot link against Python, and so a second consumer can reuse it without a second copy of the algorithm.
- **Adapter** (`fastq_chunk.hpp`) — `merge_chunk`: raw FASTQ text in, formatted FASTQ text out, plus counters and three histograms. It knows nothing about the scoring rule, and the core knows nothing about FASTQ.
- **Orchestration** (`src/zna/merge/cli.py`) — reading, chunking, threading, writing, statistics.

All three levels are exposed across the Python boundary; the fine-grained ones exist so failures localize, and `merge_chunk` is what production calls, because it captures the parse, the scan, the consensus, the construction, the formatting and the histograms in one GIL-released call. Only `int64`s and the byte table cross that boundary, so no float ever does, and `merge_chunk` is a pure function of its arguments — no globals, no clock, no RNG — which is what makes it fuzzable, safely callable with the GIL released, and byte-deterministic.

`merge_chunk` takes an explicit `(buffer, start, end)` range for each stream plus a `base_index` (`_pymerge.py:401`, `_accel.cpp:90`). The explicit ranges, rather than sliced buffers, are what let worker threads share one immutable buffer per stream with no copy (§2.7).

**The consumption protocol: whole pairs only, with a consumed count per stream.** `merge_chunk` merges pairs until either range stops yielding a complete record, and reports how many bytes it took from R1 and from R2 *separately*, so the two streams may carry different leftovers and Python never has to scan for record boundaries. A partial record at the end of a buffer is not an error — it is simply not consumed. At EOF the caller asserts that **both** buffers came out empty.

That protocol exists to make one specific silent data loss structurally impossible. A merger that consumed `min(records1, records2)` pairs and returned success would drop the trailing R2 records whenever R1 ran out first — and the R1/R2 name-sync check cannot catch that, because comparing base names cannot see records that were never read. Here the dropped records are still sitting in the caller's buffer, and the caller looks. The two ways that happens are diagnosed apart, since a wrong message sends the next person to the wrong file (`cli.py:185-198`, pinned at `tests/test_merge.py:1890` and `:1897`):

```
truncated FASTQ record at the end of <stream>        # leftover bytes, no complete record
<other> exhausted before <stream> (unequal read counts)   # leftover whole records
```

The threaded driver cuts chunks at `min(k1, k2)` whole records per stream, which is the same hazard shape — and is safe for the same reason, because whatever it does not hand to a worker stays in the buffer for the EOF check.

### 2.7 Deterministic output under threading

Because `merge_chunk` releases the GIL for its whole duration, `--threads` are real worker threads rather than processes: no fork context, no pickling of chunks and blobs, no per-worker globals, no broken-pool handling. Each stream also gets one prefetch thread, so a blocking read overlaps merge work instead of stalling it. Chunks are cut to whole records with `split_records` and handed to workers as `(buffer, start, end)` into the stream's immutable `bytes` buffer (`cli.py:82-129`, `:260-263`), so a chunk costs no copy and workers share the buffer safely. Memory is therefore one buffer per stream plus the submit window of 2×threads chunks in flight — not a per-thread copy of the input.

**Blobs are written in submission order** (`cli.py:238-243`), unconditionally and with no flag: futures are drained from the head of a bounded window, so chunk *i*'s output always precedes chunk *i+1*'s. Chunks are independent, so ordering costs only head-of-line blocking bounded by the variance in per-chunk compute time, a few percent for fixed-size chunks of Illumina reads. What it buys is that the output file is a pure function of the input and the parameters — identical at any thread count and any chunk size — which is the actual requirement for training data, and makes any future corpus defect bisectable.

The test suite asserts something weaker than that write guarantee, and the difference is worth knowing: it compares the *set* of emitted records at 1 and 3 threads (`tests/test_merge.py:1708`) and the complete statistics dict — every counter and every histogram, timing fields removed — across `--threads`/`--chunk-size` combinations of (1, 7), (3, 7), (4, 13) and (2, 1000) (`tests/test_merge.py:1809`). No test byte-compares whole output files across thread counts.

---

## 3. Fragment geometry: what the flags mean

A ZNA record is a stored sequence, and the question every downstream consumer asks of it is the same one: **which edges of this sequence are true fragment boundaries?** A model that learns transcript termini learns them from where fragment-end markers pile up, so a marker placed on a read-length cutoff is not a missing label — it is noise injected into the signal. The format therefore carries geometry per record rather than leaving it to be inferred.

### 3.1 The invariant

Every Illumina read begins at a fragment boundary and reads inward, so base 0 of a read as sequenced is a true fragment end. ZNA stores a record either in that frame or reverse-complemented, and records which in the `IS_RC` bit:

| `IS_RC` | true fragment boundary is at |
|---|---|
| clear | the **left** edge of the stored sequence |
| set | the **right** edge |

This holds per record with no knowledge of pairing, mate number, or strandedness, and it holds for un-normalized files too — there `IS_RC` is clear on everything and every record's left edge is its boundary, which is correct, because both mates of an untouched pair start at a fragment end.

The premise is that nothing removed *template* bases from a read's 5' end. Removing *synthetic* prefix — UMI, spacer, adapter — is not merely safe, it is required to recover the true boundary. What breaks the invariant is blind 5' clipping of template (`fastp --trim_front1`, 5'-side quality trimming). `zna merge` never does it: all three construction paths take from index 0 (`s1[:take1]`, `s1[:keep1]`, `s2[:keep2]`), the symmetric trim cuts only 3' ends (§1.6), and the overlap consensus changes base *values*, never positions (§1.7). Measured against the benchmark's simulated truth: base 0 of every one of 1,416,630 emitted records was a true fragment end — 0 records with a shortened 5' end, and 0 with R1 shortened.

### 3.2 The flag byte

| bit | name | meaning |
|---|---|---|
| 0 | `IS_READ1` | record is mate 1 of a pair |
| 1 | `IS_READ2` | record is mate 2 of a pair |
| 2 | `IS_PAIRED` | record has a mate, adjacent in the file |
| 3 | `IS_RC` | strand normalization reverse-complemented this record on the way in |
| 4 | `IS_FULL_FRAGMENT` | the record spans its whole fragment: both edges are boundaries |

Bits 5–7 are reserved (`src/zna/core.py:106`). Bit 4 occupies what was an unused bit, so it needs no format version bump: old readers ignore it, new readers see it clear on old files, which is the safe reading.

### 3.3 What `IS_RC` records

**`IS_RC` is a storage fact, not a mate-number fact.** It says the bases on disk were flipped relative to the bases submitted — nothing more. It cannot be recovered from the sequence, because reverse-complementing a stored record reproduces the submitted read exactly; there is no residue in the bases to test. It cannot be derived from `IS_READ1`/`IS_READ2` either: under unstranded normalization the encoder flips **one mate per pair at random**, so which mate carries the bit is a coin, not a protocol (§4.3 gives the rules). Under strand-specific normalization the bit does track the mate (and an unpaired or merged record follows the R1 rule), but a consumer that relies on that is reading a header setting through a proxy. `ENDS_BY_FLAG` (`src/zna/core.py:130`) is public precisely so that `blocks()` consumers, who get the raw flag column, never re-derive this.

### 3.4 `IS_FULL_FRAGMENT`

`IS_RC` names one edge and cannot say there is another. When the insert is at or below the read length — every overlap-merged read, and any pair whose mates cover the same interval — the record spans the whole fragment and both edges are real. That case is bit 4.

Who supplies it:

- **Full-overlap pairs: detected.** Two mates covering the identical interval are exact reverse complements, so `_is_reverse_complement(r1, r2)` on the submitted frames identifies them without knowing the insert size (`src/zna/cli.py:994`). The test probes the two end characters first, rejecting about 15/16 of ordinary pairs before the O(n) comparison runs. Measured 10051/10051 detected with 0/9949 false positives. Both mates get the bit.
- **Unpaired records: declared,** with `zna encode --treat-unpaired-as-merged`. A merged read and a genuine single-end read are byte-identical, so the encoder cannot tell them apart. Default off: an unpaired record is assumed one-ended, which under-labels rather than placing a marker at an interior position.
- **The library writer never derives it.** `ZnaWriter.write_record(is_full_fragment=…)` is always caller-supplied and says so (`src/zna/core.py:354-358`); detection lives in the `zna encode` path, which groups the stream into fragment units (`[R1, R2]` or `[single]`) and decides per unit (`cli.py:1012`, `cli.py:1040`). The codec backends never touch bit 4.
- **Copy paths carry it** rather than re-deriving it, in both `zna encode` on a `.zna` input and `zna shuffle` (§3.5, §4.5).

What the bit does **not** guarantee is that the merger was right. It asserts that whoever produced the record believes it spans its fragment. On the benchmark, 5,591 of 582,866 merged records (0.96% of merges; 5,591 per million input pairs) are not the true fragment. None is a scan failure: 0 of 848 wrong-shift merges picked a lower-scoring shift than the truth, 56% carry no sequencing error at all, and 25% survive a threshold 3.6× the default. They are fragments whose two ends are genuinely homologous — median 88.1% identity over a median 79-base overlap, concentrated in pericentromeric satellite. This is the same population §1.3 prices when it declines to raise `T_merge`: a model consuming `has_start and has_end` is trusting the merger's assertion, and that assertion is wrong about once per hundred merges on a genome that is ~50% repeat.

### 3.5 The projection: read views and copy paths

`records(with_ends=True)` yields `(seq, is_paired, is_read1, is_read2, has_start, has_end)` — six fields, or seven with a trailing `labels` on a labeled file — computed by table lookup on the flag byte:

```
has_start = is_full_fragment or not is_rc
has_end   = is_full_fragment or is_rc
```

This is the form a consumer placing fragment-end supervision wants, and it is genuinely all it wants. It is also **not invertible**:

| `IS_RC` | `IS_FULL_FRAGMENT` | byte | `(has_start, has_end)` |
|---|---|---|---|
| clear | clear | 0 | `(True, False)` |
| set | clear | 8 | `(False, True)` |
| clear | set | 16 | `(True, True)` |
| **set** | **set** | **24** | `(True, True)` |

Four reachable states, three outcomes (`src/zna/core.py:139-156`). Byte 24 is not a corner case — it is what the encoder writes whenever strand normalization reverse-complements a merged read, which under an unstranded coin is half of them. `ENDS_BY_FLAG` maps 16 and 24 alike to `(True, True)`, correctly, because both records do have two real ends; the cost is that no inverse function can exist. An inverse has to guess, and a guess silently clears `IS_RC` on every full-fragment record that passes through it: on a 200-record normalized file, 101 records carry `IS_RC` before a shuffle and 0 after, at which point `--restore-strand` returns half the corpus in the wrong orientation. The same loss runs through any ZNA → ZNA re-encode built on the projection.

Hence the rule the code enforces:

- **Reading** takes a view — `(is_paired, is_read1, is_read2)`, or `is_rc`, or `(has_start, has_end)`. Each is a projection chosen for a consumer, and none is complete. `restore_strand` is mutually exclusive with both `with_rc` and `with_ends`, since it *consumes* `IS_RC` to undo the flip; `with_rc` and `with_ends` are mutually exclusive with each other.
- **Copying** takes the byte. `ZnaReader.copy_records()` yields `ZnaRecord(seq, flags, labels)` — the stored byte verbatim, including combinations this version does not interpret and bits a future version may define — and `ZnaWriter.write_copy()` writes it unexamined. `ZnaRecord` is deliberately three fields wide and named, so code that confuses it with a `records()` tuple fails at the unpack instead of reading `flags` as `is_paired`. `write_records()` refuses the 6-tuples of `with_ends=True` outright (`core.py:499-507`) rather than reconstructing flags from them.
- **Producing** passes `is_full_fragment` to `write_record()`, which builds the byte.

Both copy halves require `ZnaWriter(preserve_normalization=True)`, and this is not bookkeeping: normalization is not idempotent (§4.3). Without that flag the codec would OR its own `IS_RC` into the flag column while reverse-complementing bases that already carry it, returning the file to an un-normalized state while the header still claims otherwise.

### 3.6 Where records land

| record | flags | geometry |
|---|---|---|
| merged read (`zna merge`, encoded with `--treat-unpaired-as-merged`) | unpaired; `IS_FULL_FRAGMENT`, plus `IS_RC` if normalization flipped it | `(True, True)` — the record is the fragment |
| merged read without the flag | unpaired, bit 4 clear | `(True, False)` — under-labeled, never wrong |
| unmerged mate | `IS_PAIRED` + `IS_READ1`/`IS_READ2`; exactly one mate of a normalized pair carries `IS_RC` | `(True, False)` and `(False, True)` — the two flagged edges are the two ends of the fragment |
| full-overlap pair reaching the encoder unmerged | both mates `IS_FULL_FRAGMENT`; one also `IS_RC` | `(True, True)` on both |
| genuine single-end read | unpaired, bit 4 clear | `(True, False)` |

A merged read is built from R1's 5' end outward with its mate suffix stripped (§1.5), so it arrives in R1's frame and the interleaved parser classifies it as a single. An unmerged pair is emitted all-or-nothing (§1.8), and pair integrity is preserved all-or-nothing (§1.8), because a lone `IS_PAIRED` record would be read downstream as half a molecule. `--npolicy` never breaks a pair: `trim3` reshapes a record and `random` substitutes within it, so neither can remove one mate of two.

`zna inspect --counts` tallies the flag column without decoding sequence and prints a (mate × `IS_RC`) cross-tabulation alongside the full-fragment count. It is the cheapest way to confirm a file's geometry before training on it: a stranded file should show `IS_RC` on essentially all of one mate and none of the other, an unstranded one about 50% of each.

---

## 4. The container and the codec

### 4.1 2-bit packing

Bases are stored two bits each — `A=00`, `C=01`, `G=10`, `T=11` — four to a byte, most significant pair first. Encoding is a 256-entry ASCII→code table (case-folding: `acgt` map to the same codes as `ACGT`); decoding is a 256-entry table mapping each packed byte to its four-character string, so decode is one table lookup and one 4-byte store per input byte. Each record's sequence is padded to a whole byte — `ceil(len/4)` bytes, padding bits zero — so every record starts byte-aligned and its offset follows from the lengths column alone. Both backends zero the pad, which is why the pure-Python and compiled codecs emit byte-identical streams.

The alphabet has exactly four symbols, so **`N` cannot be stored**, and neither can `U` or any IUPAC ambiguity code. `--npolicy` decides what happens to an `N` on the way in:

| Policy | Effect |
|---|---|
| `drop` (default) | The whole **fragment** is discarded. Applied before the writer, per fragment rather than per record, so a mate is never orphaned. |
| `random` | Each `N` is replaced by a uniformly drawn base. |
| `A` / `C` / `G` / `T` | Each `N` is replaced by that base. |

The consequence is worth stating plainly: **a decoded ZNA record can never contain an `N`.** Nothing downstream needs to handle one, and nothing downstream can tell a substituted base from a called one. Two further losses follow from the same table: decode always emits uppercase, so soft-masking does not survive a round trip, and quality scores are not part of the format at all — they are consumed by the merger (§1.7) and discarded at encode.

Two implementation notes. The compiled backend applies the policy to *any* unencodable character while the pure-Python backend substitutes only `N`/`n` and raises on the rest, so an IUPAC code such as `R` encodes on one backend and raises on the other; the divergence is pinned by a test (`tests/test_fuzz_roundtrip.py:691`) rather than silently tolerated, and it is listed below as a live inconsistency. And under `random`, the compiled encoder draws in *output* order and applies the policy *after* complementing, so a reverse-complemented record packs identically whether the two steps are fused or done separately.

### 4.2 Blocks and the columnar payload

A file is a 15-byte fixed header, then variable-length read group and description, then one 89-byte definition per label, then a chain of blocks.

```
file header   magic "ZNA\x1A" | version=3 | seq_len_bytes | header flags |
              compression method | compression level | label count |
              read-group length | description length            (15 bytes, LE)
label def     name (16 B, NUL-padded UTF-8) | description (64 B) |
              dtype code (1 B ASCII) | missing value (8 B)      (89 bytes each)
block header  comp_size | uncomp_size | record_count | flags_size |
              lengths_size                                      (20 bytes, 5x uint32 LE)
block payload one zstd frame:  flags ‖ labels ‖ lengths ‖ sequences
```

Compression is zstd at `DEFAULT_ZSTD_LEVEL = 9` unless `--level` says otherwise (`src/zna/core.py:81`, `src/zna/cli.py:1965`). `seq_len_bytes` (1, 2, or 4; CLI default 2) fixes the width of the lengths column and hence the maximum record length — 255, 65535, or 4294967295 bases. A longer record raises rather than truncating.

**The file header stores no record or block count.** Every block header carries its own, so totals come from walking the chain and seeking over each payload — measured at 2.3 µs per block, or 1.4 ms for a 38 MB / 611-block file, against 89 ms to reach the same counts by decoding. That walk is `block_index()`, which returns per-block offsets and record counts without decompressing anything (`core.py:818`); `blocks(indices=…)` then seeks past unselected blocks. Random access is therefore **block-granular** — it is per-*record* random access that the format does not provide. Blocks are flushed on an *estimated* byte size (`len/4 + 1 + seq_len_bytes` per record, default 4 MiB), so their record counts are near-uniform for fixed-length reads and vary with variable-length ones.

**A block holds whole fragments.** A fragment's reads are stored consecutively, R1 immediately followed by R2, and a fragment never straddles a block boundary. So a block is a self-contained set of molecules and any subset of blocks decodes independently of the rest of the file — which is the whole basis for `blocks(stride=…)` sharding and for block-parallel consumers: a worker handed a block never sees one mate of a pair whose other mate went to a different worker.

`ZnaWriter` enforces both halves on every write path (`OPENS_FRAGMENT` / `CLOSES_FRAGMENT`, `src/zna/core.py`). It needs no state beyond "did the last record open a fragment": that bit must equal `CLOSES_FRAGMENT[flag]` on the next record, which rejects an R1 followed by anything but its R2 *and* an R2 with no R1 in front of it, in one comparison. The size test is then asked only where a fragment ends, so the only record a block may not end on is the one whose mate has not arrived. Enforcement is what bounds the hold: a fragment is at most two records, so a block overruns its target by at most one record.

Version 2 held a weaker rule — the flush was deferred only under unstranded normalization, where the codec's pair detection requires it — and split roughly half of all block boundaries mid-fragment on every other configuration. Fixed-length reads hid it, because a constant per-record size estimate lands every boundary on an even record count. Deferring *without* enforcing was also unbounded: a run of consecutive paired R1s buffered the entire stream and wrote nothing. Version 3 files satisfy the rule by construction, and version 2 files are not readable.

The payload is columnar, not row-oriented: all flags, then all label values, then all lengths, then all packed sequence. Column order is what makes the block compress. Each stream is homogeneous, so zstd sees long runs instead of an interleave of unrelated byte types. Measured on a 25.4M-read, 150 bp library (10.76 GB FASTQ, zstd level 9, 16.5× overall): the flags stream compresses 500–1000×, the lengths stream ~1000× on uniform 150 bp reads, the label columns compress as contiguous numeric arrays, and only the sequence stream — already near-incompressible after 2-bit packing — compresses just 3–5×. Row-oriented storage would smear the three cheap columns through the expensive one.

Column order also makes partial reads possible. Flags come first, so a consumer can decompress a *prefix* of the frame and count flags without touching sequence; the shuffle's counting pass does exactly this (§4.5). Note the one hard dependency: **the block header records `flags_size` and `lengths_size` but not the label columns' width.** That is recovered as `record_count × sum(dtype sizes)` from the header schema, so a reader without the schema splits the payload at the wrong offsets — and does so silently, into plausible garbage.

### 4.3 Strand normalization

Four header bits describe the library: `STRAND_SPECIFIC`, `READ1_ANTISENSE`, `READ2_ANTISENSE`, `STRAND_NORMALIZED`. From them the writer derives three encode-time rules:

- `rc_r1 = normalized and specific and read1_antisense`
- `rc_r2 = normalized and specific and read2_antisense`
- `random = normalized and not specific`

**Stranded.** Each record is reverse-complemented iff its mate number matches an antisense rule. A record that is neither R1 nor R2 — a single-end or overlap-merged read — takes the R1 rule (`src/zna/_accel.cpp:610-614`), because a merged read is R1 followed by the reverse complement of R2's tail and is therefore R1-oriented.

**Unstranded.** There is no sense strand to normalize to, so the rule is symmetry instead. An FR pair covers the two ends of one fragment pointing inward, so as sequenced the mates sit in opposite frames. For each paired R1, one fair coin decides which of the two mates to reverse-complement; **exactly one always is**, which puts both into a single common frame. Unpaired records get an independent coin.

The codec finds a fragment's two mates by looking at adjacent records within one block, which it may do because §4.2 guarantees they are there. It raises rather than guessing if they are not, so a caller reaching the backend directly gets an error instead of a half-normalized fragment.

The resulting geometry is the reason `IS_RC` is stored per record rather than inferred: in the common frame the reverse-complemented mate has its **right** edge on the real fragment boundary and its left edge at a read-length cutoff, and the other mate is the mirror image. §3.3 gives the semantics and §3.5 the read/copy rule that follows from them.

**Normalization is not idempotent.** Applying it a second time returns the data to an un-normalized state while the header still reports `strand_normalized`. Every ZNA→ZNA path therefore *copies* orientation rather than re-deriving it: the writer's `preserve_normalization=True` disables all three rules, and `copy_records()` → `write_copy()` carries the flag byte verbatim, including bits this version does not interpret.

### 4.4 Labels

A label is a fixed-width numeric column, one value per record, defined once in the file header. Types are the SAM/BAM auxiliary type codes `A c C s S i I f`, plus three 8-byte ZNA extensions `d` (float64), `q` (int64), `Q` (uint64). Each definition stores a name (≤16 UTF-8 bytes), a description (≤64 bytes), the dtype code, and a missing-value sentinel packed in the dtype's native little-endian form and zero-padded to 8 bytes. The sentinel defaults to zero for every dtype.

Values are extracted from the FASTQ header by a SAM-tag rule:

1. Split the raw header on **any** whitespace. Splitting on space alone is not enough — tools such as fastp append bare space-separated tokens like `merged_150_87`, which would otherwise be glued onto the last tag's value.
2. **Skip token 0**: it is the read name.
3. Match each remaining token as `KEY:TYPE:VALUE`. `KEY` is everything before the first colon and must be non-empty; the byte after it must be followed by another colon. Keys may be any length, so both two-character SAM tags and long custom keys work. `TYPE` is parsed for position and then **ignored** — the schema's dtype governs conversion, not the header.
4. Convert by dtype: `A` takes the single byte's ordinal, `f`/`d` go through `float()`, everything else through `int()`. The compiled fast path handles only a short run of ASCII digits and hands everything else to CPython, so both the value and the exception are the reference implementation's by construction.
5. Any tag absent from the header gets that label's missing value.

Matching uses the label's `tag` when one is defined, otherwise its `name` — which lets a file parse `AS` from the input and store the column as `aln_score`. `tag` is an encode-time field only and is **not** written to the file, so a decoded file reports names.

On disk the columns sit between the flags and the lengths, in header order, each a contiguous run of `record_count` fixed-width values. A reader can therefore skip them (`blocks(labels=False)`), receive them columnar (`labels=True`), or get them per record — but never by accident: `blocks()` on a labeled file raises unless the caller says which. Skipping does not avoid inflating the bytes, since the block is one zstd frame; it avoids unpacking them into Python objects, which is where the cost is (about 1.9× on a three-column file).

### 4.5 The shuffle

Training data wants arbitrary record order, and `blocks()` sharding hands each worker contiguous runs, so any grouping in the input becomes a biased sample. `zna shuffle` is a two-pass bucket shuffle with bounded memory.

The unit of shuffling is a **fragment**, not a record: an unpaired record, or a paired R1 together with the R2 that follows it. Pairs must stay adjacent and in order, so a pair may never be split across buckets — otherwise the mates end up arbitrarily far apart in the output and the pairing is unrecoverable. A unit start is any flag byte with `(f & 0x05) != 0x04`, i.e. anything that is not a paired non-R1.

**Pass 0 — count.** Walk the block chain. For each block, decompress only the first `flags_size` bytes of the frame (a zstd frame reads as a prefix) and translate the flag column through a 256-entry table to count unit starts; record counts come free from the block headers. No base is decoded.

**Pass 1 — distribute.** `K` is chosen so a loaded bucket is about `buffer_bytes`: `bytes_per_unit = max(file_size / n_units, 64)`, `K = ceil(n_units / (buffer_bytes / bytes_per_unit))`, clamped to `min(n_units, 4096)`. Each unit is assigned to a uniformly random bucket; an R1 pins the bucket so its R2 lands in the same one. Buckets are ZNA files written at zstd level 1 — they are read exactly once and deleted, so spending the source file's compression level on them buys a ratio nothing benefits from. The output file keeps the user's level.

**Pass 2 — collect.** Bucket order is itself shuffled. Each bucket is read whole into memory, its units Fisher-Yates shuffled, and appended to the output. Uniform independent assignment, a uniform permutation within each bucket, and a uniform bucket order compose to a uniform permutation over all units.

Every writer in both passes runs `preserve_normalization=True` and copies each record's flag byte verbatim. A shuffle is a pure permutation, so orientation must be carried, not re-derived (§4.3): re-deriving would reverse-complement once on the bucket pass and again on the output pass, leaving sequences that still look plausible but whose `IS_RC` no longer corresponds to the real fragment boundary — and the same applies to `IS_FULL_FRAGMENT`, which no projection can reconstruct (§3.5).

What survives a shuffle: every record's full flag byte, its label values, R1-before-R2 adjacency within each pair, and the file header — read group, description, `seq_len_bytes`, strand flags, compression settings, and label schema. Given the same seed, input, and buffer size, the output is byte-reproducible. What does not survive is block boundaries; records are re-blocked at the output writer's block size.

---
