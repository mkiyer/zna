# Read-merge redesign: one score, two calibrated thresholds

Status: **implemented** in `src/zna/merge/`, shipped as `zna merge`. Written
2026-08-11 against `lib/hulkrna/merge/`, where the tool was developed; ported into
zna unmodified on 2026-08-12 (see [READ_MERGE_PORT_TO_ZNA.md](READ_MERGE_PORT_TO_ZNA.md)).
Supersedes the original fastp-derived design `READ_MERGE_TOOL_DESIGN.md`, which was
retired 2026-08-12 and lives in hulkrna's git history.

**The design is settled. Do not re-derive it.** The three implementation decisions it
left open are recorded immediately below; ideas that look good and were measured and
rejected are in [MERGE_TOOL_AUDIT.md](MERGE_TOOL_AUDIT.md) §(c), each with the number
that closed it.

Three places where the implementation had to decide something this document leaves
open or approximate. Each is pinned by a test.

1. **The match weight is the exact `log2((1-e)/0.25) = 1.9855`, not 2.0.** §5b's
   "matches needed" table requires it (at T=18, d=0 the answer is 10, which 2.0 would
   make 9). Consequence: `T_trim = 8` is described in §5b as "4 perfect matches", but
   4 matches score 7.94 bits and fall 0.06 short — the smallest overlap that can be
   trimmed is **5** bases. Measured detection and redundancy are unaffected; the
   4-base case was already the coin-flip §5 argues against.
2. **The tie margin is not reported; the score distribution is** (§8b.5). ~~§8's margin
   suggestion is incompatible with §8's pruning.~~ **Corrected 2026-08-12: that claim was
   wrong.** Setting the pruning floor to `max(t_trim, best − margin_req)` instead of
   `best` yields an exact bounded margin — two independent implementations matched a
   fully unpruned scan's margin decision with **0 disagreements** at 4/6/8/10/16 bits,
   for +4.2% of the kernel at margin 10. The margin is still not computed, but for a
   different and better reason: the audit measured what it would be used for and found
   ambiguity is a *microsatellite detector*, not an identifiability problem (§10).
   The per-pair score histogram remains the shipped diagnostic.
3. **The inner loop accumulates mismatches branchlessly in blocks of 8.** On a wrong
   shift each comparison is a coin flip, so a per-position `if` mispredicts
   constantly; hoisting the bail out of the block is worth **3.7x** on the scan and is
   bit-for-bit identical. Without it the redesign cost 5.8x the old kernel; with it,
   1.6x (and **+14%** on the whole tool, end to end — see §8b note).

---

## 1. Why redesign

The current rule is a port of fastp's `OverlapAnalysis::analyze`: slide R1 against
revcomp(R2) and accept the **first** offset whose overlap has few enough
mismatches. Accepting it requires four thresholds that interact:

| parameter | value | what it is trying to decide |
|---|---|---|
| `merge_overlap_require` | 20 | is the overlap long enough to merge? |
| `trim_overlap_require` | 3 | …long enough to trim? |
| `merge_diff_limit` | 3–4 | absolute mismatch cap |
| `merge_diff_percent_limit` | 20% | fractional mismatch cap |
| *(proposed)* `exact_below` | 10 | short overlaps must match exactly |

Plus `correction` with its own two quality cutoffs (Q30 / Q14). That is seven or
eight knobs, and three of them (`diff_limit`, `diff_percent`, `exact_below`) are
three different attempts at one question — **is this overlap real?** — none of
which asks it directly.

The symptoms are what you would expect from proxies rather than the question:

- `floor(0.20 × 4) = 0` makes a 4-base overlap require an exact match, which
  sounds strict but fires by chance 1 time in 256 — and there are ~146 offsets.
  **5.27%** of genuinely unrelated read pairs get an overlap "detected".
- `floor(0.20 × 5) = 1` lets a 5-base overlap through *with* an error, which is 4×
  more likely by chance than the 4-base exact case. The budget is not too tight;
  for short overlaps it is too loose.
- First-accept means a spurious short hit in the forward scan preempts the reverse
  (read-through) scan, so **~1.7% of read-through pairs** are emitted unmerged with
  adapter retained.
- On a realistic insert distribution the current rule picks the **wrong shift for
  2.77%** of pairs.

Adding `exact_below` would patch the second symptom and make the knob count worse.

## 2. What we are actually inferring

Given R1 and revcomp(R2), there is exactly one unknown: the fragment length `L`,
equivalently the offset `s` of revcomp(R2) relative to R1.

Write both reads on one axis. revcomp(R2)'s fragment portion always ends at its
own end, so for reads of length `len2`:

```
s = L - len2
```

uniformly — `s ≥ 0` is a normal overlap, `s < 0` is read-through, and `s = 0` is
exact full overlap. **One axis, one scan.** That alone deletes the two-scan
structure, the forward-before-reverse ordering artifact, and the `require`
parameter (which only existed to bound each scan).

For a candidate `s` we observe `n` aligned base comparisons containing `d`
mismatches. Two hypotheses:

- **H₁, real overlap** — mismatches are sequencing errors, rate `e ≈ 2ε` since
  either read can be wrong; ~1% for typical data.
- **H₀, chance alignment** — unrelated bases agree 1/4 of the time, so the
  mismatch rate is 3/4.

## 3. The score

The log-likelihood ratio, in bits:

```
score(s) = matches × log₂((1−e) / 0.25)  +  mismatches × log₂(e / 0.75)
         ≈ 2 × matches  −  w × mismatches            w = log₂(0.75/e) ≈ 6.2 at e = 1%
```

Every matching base is worth exactly **2 bits** — that is `log₂ 4`, the
information in agreeing on one of four bases. Every mismatch costs ~6 bits. No
tuning: both weights fall out of the error rate.

**The threshold is derived, not chosen.** The score *is* the log-odds against
chance, and for any likelihood ratio `E_{H₀}[LR] = 1`, so `P(score > T) ≤ 2^−T`
per shift. Over `N ≈ 2 × readlen` candidate shifts, a union bound gives

> **Read this as what it is (added 2026-08-12).** `α` is the tolerated rate of
> detection **against chance alignment of uniform random sequence**, and it holds
> there with ~14 bits to spare: 0 spurious merges in 40,000 uniform-random pairs at
> T = 28, at every read length from 50 to 300. It is **not** a bound on detections
> against *real* sequence, where reads share genuine homology and repeat content:
> on shuffled real pairs from one cfRNA library the detection rate is 1.7% at T = 28
> and still 1.3% at T = 100, because 89% of those involve a >90% two-base read.
> Raising T does not buy what the formula suggests it buys — the score is measuring
> real similarity at that point, not chance. See docs/MERGE_TOOL_AUDIT.md.

```
T = log₂(N / α)
```

for a target spurious-merge rate `α`. At `readlen` 150 and `α = 10⁻⁶`,
`T = 28 bits` — about 14 perfect matching bases, or 17 with one mismatch.

Because evidence accumulates **linearly** while the false-positive rate falls
**exponentially**, being very conservative is nearly free: tightening α from 10⁻³
to 10⁻⁹ costs only 10 extra matching bases.

## 4. The decision

```
for every shift s:  compute score(s)
s* = argmax score(s)
if score(s*) ≥ T:  merge at s*        else:  keep both reads unchanged
```

`argmax` instead of first-accept is what removes the preemption bug — the scan can
no longer be captured by a spurious short hit before it reaches the real offset.

**Merge and trim keep separate thresholds — because their errors cost different
amounts.** A wrong *merge* produces a chimera: actively false sequence. A wrong
*trim* removes a few real bases from a read tail, which is cheap. Same evidence,
different consequences, so:

```
score(s*) ≥ T_merge   ->  merge at s*                (one full-fragment record)
T_trim ≤ score < T_merge -> keep both reads, trim the redundant overlap
                            off R2's 3' end so the fragment is tiled once
score(s*) < T_trim    ->  keep both reads unchanged
```

Both thresholds are in **bits** and both mean the same thing — tolerated spurious
rate — so this is one calibrated scale read at two points, not two unrelated
proxies. `T_trim` sits far lower than `T_merge` purely because the downside is
smaller.

This also preserves the boundary invariant ZNA depends on: trimming only ever
removes bases from R2's **3'** end, so base 0 of each mate stays a true fragment
boundary.

**The critical structural point:** because the decision is `argmax` rather than
first-accept, a low `T_trim` no longer costs merges. A spurious 4-base hit scores
8 bits; a real 40-base overlap scores 80 and wins the argmax outright. In the old
design the two were in competition — the scan stopped at whichever it met first —
which is why a permissive trim threshold silently suppressed ~1.7% of
read-through merges. Aggressive trimming and reliable merging stop trading off
against each other.

## 5. Measured against the current rule

4000 pairs, insert ~ N(200, 70) truncated to [50, 400], 2×150, 0.5% per-base
error, real adapters with read-through:

| rule | wrong shift | redundant bases/pair | FP on unrelated pairs |
|---|---|---|---|
| current (3 thresholds, first-accept) | **2.77%** | 0.09 | **5.27%** |
| LR, α = 10⁻³ (T = 18 b) | **0.00%** | 0.12 | 0.00% |
| LR, α = 10⁻⁶ (T = 28 b) | **0.00%** | 0.31 | 0.00% |
| LR, α = 10⁻⁹ (T = 38 b) | **0.00%** | 0.52 | 0.00% |

One threshold beats three on every axis at once. **Recommended default
α = 10⁻⁶**: redundancy costs 0.31 bases per ~300-base pair (0.1%), against
essentially zero chimeric merges.

Note the redundancy column is the honest cost of statistical rigour, and it is
small for a reason worth stating: **short overlaps are rare.** Overlap = 2·readlen
− insert, so 4–11 base overlaps come from inserts of 289–296 — the far tail of any
real library. Aggressively "trimming" 4-base overlaps is not aggressive trimming;
it is coin-flipping on ~1.7% of pairs, and it is what produces the 5.27% false
positive rate.

## 5b. Choosing the two thresholds

**What the score requires, concretely** (e = 1%, so a mismatch costs 6.23 bits):

| mismatches in the overlap | matches needed at T=18 | at T=21 | at T=28 |
|---|---|---|---|
| 0 | 10 | 11 | 15 |
| 1 | 13 | 14 | 18 |
| 2 | 16 | 17 | 21 |
| 3 | 19 | 20 | 24 |

**To reproduce the current tool's *merge* behaviour, use T_merge ≈ 21 bits.** The
current merge rule (≥20 overlapping bases, ≤4 mismatches) corresponds to a
per-shift chance probability of ~4.6 × 10⁻⁷, i.e. α ≈ 1.4 × 10⁻⁴ over ~300 shifts,
which is `log₂(300 / 1.4e-4)` ≈ 21 bits ≈ 11 perfect matching bases. Recommended
default is slightly stricter, **T_merge = 28 bits (α = 10⁻⁶)**, since a chimera is
the expensive error and the extra cost is four matching bases.

**T_trim, swept across literature-realistic insert distributions** (2×150, 0.5%
error, T_merge = 28):

| T_trim | redundant bases/pair | wrongly-cut bases/pair |
|---|---|---|
| none (always trim at argmax) | 0.00 | 0.08 – **0.70** |
| **8 bits (4 perfect matches)** | **0.02 – 0.04** | **0.01** |
| 12 bits | 0.04 – 0.09 | 0.00 |
| 28 bits (no separate trim) | **0.22 – 0.62** | 0.00 |

Ranges span N(120,40) degraded/FFPE, N(200,70) typical mRNA, N(300,80) long
insert, and N(400,100) very long. Merge rate varies enormously across them — 100%,
88%, 44%, 14% respectively — but the threshold behaviour does not, which is the
point: the score is calibrated in evidence, so it does not need retuning per
library.

**Recommended default: T_trim = 8 bits.** This vindicates the existing
`trim_overlap_require: 3` instinct — accepting 4-base overlaps *for trimming* is
correct, because the measured cost is 0.01 wrongly-cut bases per pair. The old
design's mistake was applying that same permissive rule to *merging*, and letting
first-accept turn it into a merge-suppressor.

## 6. Options considered

- **A. Keep the threshold stack and tune it.** Rejected: three proxies for one
  question, and every fix adds a knob. This is the status quo plus `exact_below`.
- **B. Likelihood-ratio score with a derived threshold.** *Recommended.* Replaces
  five thresholds with one interpretable rate.
- **C. Full local alignment (Smith-Waterman).** More general — it handles indels —
  but O(n²) per pair with gap-open and gap-extend as *new* parameters. Illumina
  errors are substitution-dominated, so the generality buys little. Worth keeping
  as the escape hatch if a chemistry with indel errors ever matters.
- **D. k-mer seed-and-extend to propose shifts, then score.** Not a different
  decision rule — a *speed* optimisation that composes with B. Hold until the scan
  is measured to be a bottleneck; §8's pruning is likely enough.
- **E. Do not merge; let ZNA/khorana cope.** Rejected: it maximises the redundancy
  we are trying to eliminate and discards the fragment-end information that only
  the merger can recover.
- **F. Quality-aware LR.** B, refined — see §7. The natural end state.

## 7. Extension: quality-aware scoring (removes the last knob)

`e` is the only remaining tunable, and we can read it off the data instead. Per
aligned position, from the two Phred scores, `eᵢ ≈ p₁ + p₂`, so

```
match at i     -> + log₂((1 − eᵢ) / p_null)
mismatch at i  -> + log₂(eᵢ / (1 − p_null))
```

Precompute a table indexed by `(Q₁, Q₂)` into fixed-point integers: two lookups
per position, numba-friendly, no logs in the inner loop.

This also **eliminates `correction` and its Q30/Q14 cutoffs**: at a mismatch, the
consensus base is simply the one with the higher posterior, which is what the
score is already computing. Three more parameters gone, and the overlap consensus
becomes principled rather than "R1 always wins".

**Low-complexity guard, free.** `p_null` above is the probability two unrelated
bases agree — 1/4 for uniform sequence, but ~1 inside a polyA run. Computing it
from the two reads' actual base composition (`p_null = Σ_b f₁(b)·f₂(b)`) makes a
matching base inside a homopolymer worth ~0 bits, so low-complexity regions stop
generating confident spurious overlaps. Still parameter-free. **The current design
has no protection here at all**, and polyA is common in RNA-seq — this is the one
respect in which the redesign fixes something we have not yet measured.

> **Measured 2026-08-12, three times independently: do not build this.** It does not
> fire on the case it was designed for. The dominant low-complexity motif here is a
> *dinucleotide* repeat, where `p_null = 0.500` exactly and a match is still 0.986
> bits — a 74-base perfect CA overlap still scores 72.9 against T = 28. Net effect
> measured: **−0.55 pp merge rate overall and ~2 pp on exactly the low-complexity
> reads it was meant to protect**, 2.26% of argmaxes moved, wrong-shift rate
> unchanged. Cost +3.4% end to end when fused into `_scan`. One honest caveat: the
> ground-truth subset under-samples low-complexity reads 15x, so for the affected
> population the verdict is *unmeasured*, not *measured null*. Also note the formula
> as written is orientation-ambiguous (p90 = 0.462 on R1 vs revcomp(R2), 0.269 on raw
> R2 — opposite signals), which would need settling first. See docs/MERGE_TOOL_AUDIT.md.

## 8. Implementation notes

- **Cost is unchanged**, O(readlen²) worst case, with better pruning available:
  evaluate shifts in decreasing overlap length and skip any shift whose ceiling
  `2n` cannot beat the best score so far; within a shift, bail when
  `current + 2 × remaining < best`. A real overlap is usually found early and then
  prunes almost everything, which is strictly better than first-accept's bail.
- Reads of unequal length (post quality-trim) need no special case: `s` is defined
  by the offset, and the compared region is the intersection.
- Ties: if the best and second-best shifts are within a bit or two, the pair is
  genuinely ambiguous. Reporting the margin in the JSON stats is cheap and worth
  doing; declining on a small margin would be another knob, so leave it out until
  the data says otherwise.
- `min_read_length` stays. It is an orthogonal QC filter, not part of overlap
  inference, and the all-or-nothing pair rule (the `if paired:` branch at the end of
  `process_pair`) must be kept exactly as is — it is what keeps orphan reads out of
  ZNA.
- **Trim guard, no new parameter.** If a trim would leave R2 below
  `min_read_length`, keep both reads *untrimmed* instead. Measurement says this is
  unreachable in practice — the longest overlap that can land in the trim band is
  `13 + 4.1·d`, and across three insert distributions the worst observed trim was
  22 bases, leaving 128 of 150 — but the guard costs one comparison and turns a
  would-be whole-fragment discard into a no-op. Trimming also only ever removes
  from R2's **3'** end, so base 0 stays a true fragment boundary and ZNA's
  boundary invariant holds.

## 8b. Implementation checklist

1. `overlap.py`: replace `_scan` with a single-axis scorer over
   `s ∈ [-(len2-1), len1-1]` returning `(best_s, best_score, overlap_len)`. Keep
   it `@njit`-able — no Python objects in the loop, integer/float locals only.
   Consider fixed-point integer scoring (hundredths of a bit) if faster under
   numba; measure rather than assume.
2. Branch-and-bound pruning as in §8, so the common case stays at or below the
   current cost. **Benchmark against the current kernel on the same pairs** — this
   loop runs over every read pair in every library.

   *Measured* (min-of-9, in-memory, both kernels interleaved in one process, 20k
   pairs per cell, 2x150):

   | cell | old | new | ratio |
   |---|---|---|---|
   | realistic N(200,70) | 1.34 us | 2.49 us | 1.85x |
   | long insert N(300,80) | 2.08 us | 3.58 us | 1.72x |
   | short insert N(120,40) | 1.41 us | 1.71 us | 1.21x |
   | quality-trimmed mix | 1.38 us | 2.49 us | 1.80x |
   | unrelated (worst case) | 3.03 us | 4.17 us | 1.37x |

   The scan is genuinely more work: the LR rule tolerates ~24% mismatches on a long
   overlap, so *rejecting* a shift costs ~0.32n comparisons against the old rule's ~7.
   That is the price of asking the right question, and pruning cannot remove it —
   the branchless block above is what keeps it to 1.6x instead of 5.8x. End to end
   the merge tool goes **7.4 -> 8.5 us/pair (+14%)**, since the kernel is under a
   third of the per-pair cost. A k-mer seed pass (option D) was implemented and
   measured: **no gain** over the block loop (1.58x vs 1.57x), so it was dropped.
3. `pairs.py`: replace the merge/trim/keep cascade with the three-way threshold
   test. `_build_merged` is unchanged (it already takes a direction and shift).
   The `keep2 <= 0` special case disappears — it is subsumed by `score ≥ T_merge`.
4. `MergeParams`: drop `merge_require`, `trim_require`, `diff_limit`, `diff_pct`;
   add `t_merge=28.0`, `t_trim=8.0`, `err_rate=0.01`. Keep `min_read_length` and
   (for now) `correction`; §7 removes the latter.
5. `cli.py`: `--t-merge` / `--t-trim` replace the four overlap flags. Add the
   score distribution and merge/trim/keep split to the JSON stats.
6. `config.yaml` and `zna.smk` params updated to match.

## 8c. Tests

1. **Threshold arithmetic** — pin "matches needed at d mismatches" for T=18/21/28
   so the score cannot silently drift.
2. **False positives** — unrelated random pairs; assert the spurious detection
   rate is under a bound (< 0.5% at n=20 000). This is the test that would have
   caught the current 5.17%.
3. **Detection** — synthetic pairs with known overlaps 4…40 at 0.5% error; the
   recovered shift must equal the true shift whenever score ≥ T_trim.
4. **Read-through** — inserts below the read length must merge to exactly the
   fragment with adapter removed, for inserts from `min_read_length` upward.
5. **Boundary invariant** — for every emitted record base 0 is a true fragment
   boundary, and every merged record equals its fragment exactly. Build reads at
   full cycle length with filler after the adapter: a shorter-than-cycle read
   makes the true shift `L - len(read2)`, not `L - readlen`, and that discrepancy
   produced two false alarms during the ZNA audit.
6. **Trim guard** — a trim that would leave R2 under `min_read_length` keeps both
   reads untrimmed rather than discarding the fragment.
7. **Parity where it should hold** — on clean pairs with overlap ≥ 30 the new and
   old rules must merge at the same shift; divergence is expected only in the
   short/noisy band.

## 9. Final parameter count

| | before | after |
|---|---|---|
| overlap inference | `merge_require`, `trim_require`, `diff_limit`, `diff_pct`, (`exact_below`) | **`T_merge`, `T_trim`** (bits; one calibrated scale) |
| consensus | `correction`, `GOOD_QUAL`, `BAD_QUAL` | *(none — falls out of the score)* |
| QC | `min_read_length` | `min_read_length` |

**From eight knobs to three**, and all of them have units. `T_merge` and `T_trim`
are read off one evidence scale in bits, each expressible as a tolerated spurious
rate `T = log₂(N/α)`; `min_read_length` is an orthogonal QC filter. Defaults:
`T_merge = 28` (α = 10⁻⁶), `T_trim = 8`.
