# `zna merge` against fastp on simulated ground truth — results

Run 2026-08-12 against zna at `710957d` (0.4.0-unreleased), on **1,000,000 pairs**
simulated from hg38. This is the measurement
the 0.4.0 benchmark plan specified, and it gates the 0.4.0
tag. Everything before it proved that `zna merge` matches *itself*; this is the first
thing that asks whether its answers are **true**.

> **Re-verified 2026-08-13 at the 0.4.0 release commit.** Every number in this document
> was reproduced exactly — `summary.json` and `report.md` diff clean against the run
> above apart from output paths and the throughput line, which is run-to-run noise on a
> shared machine (zna 1.86 → 1.66 µs/pair, fastp 2.83 → 2.93). Per-record provenance
> changes headers only, never a merge decision, so nothing here was expected to move and
> nothing did. The re-run also confirmed the two structural properties the corpus rests
> on: flag byte 24 survives `--shuffle` on 291,810 records, and `--restore-strand`
> reproduces the original base multiset over all 1,416,630.
>
> The simulator emits no `N`, so the N policy is exercised by injecting one no-call into
> ~1.5% of reads, 3'-biased. At that rate, of 30,000 injected: **5,879 rescued from the
> mate** under *both* policies (rescue precedes the policy, so this must agree), then
> `random` substituted exactly the remaining **24,121**, and `trim3` removed 913,823
> bases — 0.43% of emitted bases, below the 1% warning. `trim3` also costs ~1,900 merges
> against `random` (580,971 vs 582,777), which is D4a's coverage retry declining to
> merge a pair that no longer tiles its fragment.

Reproduce with [`scripts/merge_bench/simulate.py`](../scripts/merge_bench/simulate.py)
and [`compare.py`](../scripts/merge_bench/compare.py) — see that directory's README. The
simulator is deterministic in its seed (verified byte-identical across runs), so the
input is regenerable rather than stored.

```
simulate.py --genome hg38.fa --n-pairs 1000000 --read-length 150 \
            --frag-min 60 --frag-max 450 --error-rate 0.002 \
            --quality-model novaseq --seed 42
compare.py --sim-prefix sim --out results/ --threads 4
```

Fragment lengths are **uniform on [60, 450]**, which populates all five geometric
regimes at equal density: 382,657 pairs with no overlap at all, 381,379 partial,
230,815 read-through, 2,538 exactly full, 2,611 exactly at the zero-overlap boundary.
**The overall merge rate from this input is not a library number** — read the per-bin
curve instead.

---

## The headline

| | zna merge | fastp 1.1.0 |
|---|---:|---:|
| **C1 boundary violations** (base 0 of an emitted read is a true fragment end) | **0** of 1,416,630 records | 0 of 1,458,910 |
| **C3 orphans** (one mate emitted without the other) | **0** | 0 |
| Duplicated fragment bases left in unmerged pairs | **28,074** (2.25 per 10k) | 1,181,978 (85.8 per 10k) |
| Merge sensitivity at true overlap ≥ 15 | **99.83%** | 92.98% |
| Merges at true overlap 0 (chimeras) | 4,743 / 385,268 = 1.23% | 2,393 = 0.62% |
| …at `--threshold-merge 60`, which matches fastp's rate | 0.60% at 92.63% sensitivity | — |
| Merged records equal to the true fragment | **86.59%** | 85.90% |
| Overlap errors recovered, vs the best possible | **90.35%** | 74.07% |
| Wrong merges where the scan picked a *lower*-scoring shift than the truth | **0** of 848 | 0 of 1,109 |
| Throughput, same input and thread count | **1.53 µs/pair** | 2.83 µs/pair |

**C1 holds exactly, over every emitted record, against independent truth.** That is the
property the ZNA fragment-boundary contract rests on, and it had only ever been checked
against the tool's own inferences. It is now checked against the genome.

**C2 — a merged record is its fragment exactly — does not hold universally**, and the
rest of this document is about why: 5,591 merges (0.96% of merges) are at the wrong
shift, and every one of them is the genome repeating rather than the merger
malfunctioning.

---

## 1. Sensitivity: the predicted step at 15 bases is real

| true overlap | pairs | zna merged | fastp merged |
|---|---:|---:|---:|
| 0 | 385,268 | 4,245 (1.10%) | 2,082 (0.54%) |
| 1–4 | 10,156 | 60 (0.59%) | 36 (0.35%) |
| 5–9 | 12,973 | 106 (0.82%) | 67 (0.52%) |
| 10–14 | 12,762 | 122 (0.96%) | 69 (0.54%) |
| 15–19 | 12,747 | **11,797 (92.55%)** | 53 (0.42%) |
| 20–29 | 25,601 | 25,540 (99.76%) | 146 (0.57%) |
| 30–49 | 51,467 | 51,466 (100.00%) | 48,973 (95.15%) |
| 50–99 | 127,343 | 127,343 (100.00%) | 127,343 (100.00%) |
| 100–149 | 128,330 | 128,330 (100.00%) | 128,330 (100.00%) |
| 150 (mates are exact reverse complements) | 2,538 | 2,538 (100.00%) | 2,538 (100.00%) |

The cliff sits exactly where the derivation puts it: `ceil(28 / 1.9855) = 15` clean bases
is the shortest overlap that can reach 28 bits, and the tool merges essentially nothing
below it and essentially everything above. The 7.45% shortfall in the 15–19 bin is not
slack in the rule — it is arithmetic: a 15-base overlap scores 29.8 bits clean, and a
single mismatch costs 8.2, so **any** error in a minimal overlap puts it under
threshold. At a 0.4% per-position disagreement rate, 5.8% of 15-base overlaps carry one;
observed 7.45% across a bin that also contains 16–19.

The nonzero rates in the 1–14 bins are **not** detections of the true overlap — they are
wrong merges at some other shift, and they are counted again below. This table counts
only merges that were **emitted**; 504 more were merged and then filtered below
`--min-read-length`, which is why the chimera total above (4,743) exceeds this table's
first row (4,245). See §4.

fastp's curve reads its own defaults directly: `--overlap_len_require 30` puts its cliff
at 30, and `--overlap_diff_limit 5` / `--overlap_diff_percent_limit 20` cost it 4.85% of
the 30–49 bin.

**Read-through fragments** (`L < readlen`, both mates running into adapter — the geometry
that used to truncate 0.271% of merges) reconstruct at **100.00%** in all four length
bins, 230,815 pairs, both tools.

---

## 2. The chimeras are the genome, and they are not a threshold problem

4,743 of 385,268 no-overlap pairs merged (1.23%), plus 848 merges at the wrong shift
where a real overlap existed: 5,591 wrong merges in total. Every one was examined by
re-running the shipped kernel on the pair:

- **median identity 88.1%**, minimum 78.2%, over a **median 79-base** overlap; 5,413 of
  5,591 are ≥80% identical.
- **56% of them carry no sequencing error at all**, so they are structural, not noise.
- their evidence has a **median of 63 bits** against a 28-bit threshold; a quarter are
  over 101 bits.
- the top 5-Mb hotspots are `chr5:45–50`, `chr18:15–20`, `chr1:120–125`, `chr16:35–40`,
  `chr15:15–20`, `chr19:25–30` — every one pericentromeric or an acrocentric short arm,
  i.e. satellite.

So the fragment's two ends genuinely align. **Raising the threshold does not fix this**,
which reproduces an earlier finding
against a completely different measurement:

| `--threshold-merge` | wrong merges surviving |
|---:|---:|
| 28 (default) | 5,591 (100%) |
| 40 | 4,321 (77%) |
| 60 | 2,937 (53%) |
| 100 | 1,403 (25%) |
| 150 | 620 (11%) |

Reaching a tenfold reduction means 150 bits — five times the default — and the
sensitivity that buys it back is not free. The evidence is real; the sequence is
ambiguous.

**And the search is not at fault.** For all 848 wrong-shift merges where a true overlap
existed, the kernel was also asked to score the *true* shift. It picked a lower-scoring
shift **zero** times: the alternative alignment beat the true one by a median of **44
bits**. The argmax is maximising correctly; the objective is simply maximised somewhere
else. That distinction is what separates "fix the code" from "the input does not
determine the answer", and it is the reason no change is proposed here.

fastp merges the same class of pair at roughly half the rate — for the same reason its
sensitivity is lower, not because it is more discerning: its 30-base floor and 20%
mismatch budget reject short and imperfect alignments of both kinds.

---

## 3. The quality-aware consensus, measured against its own ceiling

Every merged record of the correct length is scored three ways against the true
fragment: what the tool emitted, what **R1-wins** (no quality model) would have emitted,
and the **oracle floor** — the best any consensus could do, which is not zero, because a
position where *both* mates are wrong is unrecoverable.

| | zna merge | fastp |
|---|---:|---:|
| records scored | 577,275 | 537,269 |
| errors, R1-wins baseline | 182,184 | 161,158 |
| errors, actual | 84,469 | 82,610 |
| errors, oracle floor | 74,037 | 55,117 |
| **recovered, of what was recoverable** | **90.35%** | 74.07% |
| records left worse than doing nothing | 895 (0.16%) | 166 |

`zna merge` recovers 97,715 of 108,147 recoverable errors and costs 895 bases elsewhere,
a ratio of 109:1. This is the first *controlled* test of the 35.2% claim in
the redesign — that number was against fastp's
Q14/Q30 cutoffs on real reads with no ground truth, and it survives contact with truth.

**The control that makes this a real measurement.** A separate 200,000-pair run with
`--quality-model flat` (constant `I`, the obvious way to write this simulator) gives a
recovery of **0.0% for both tools**: with equal qualities every mismatch is a tie, R1
always wins, and the consensus is not exercised at all. Anyone re-running this must use
`novaseq`, or metric 4 measures nothing. The simulator therefore draws the quality first
from a four-bin position-degrading profile and the error from it at `10^(-Q/10)`, so Q
is monotone in the real error probability — which is the only property the posterior
uses.

---

## 4. Contracts, in full

| | zna merge | fastp |
|---|---:|---:|
| records checked for C1 | 1,416,630 | 1,458,910 |
| **C1: base 0 not a true fragment boundary** | **0** | 0 |
| C1: correct-length merge not framed at offset 0 | 0 | 0 |
| C2: merged record not the fragment (wrong length) | 5,591 | 3,502 |
| **C3: orphaned mates** | **0** | 0 |
| same molecule emitted twice (merged *and* paired) | 0 | 0 |
| **C4: R1 shortened** (a 3' end moved on the read assumed intact) | **0** | 0 |
| C4: duplicated fragment bases surviving | 28,074 | 1,181,978 |

C1 is checked on **every** emitted record — merged, kept and trimmed — by comparing its
first 24 bases against the fragment end it claims to start at. A real 5' end mismatches
in 0 or 1 of those; a shifted one mismatches in ~18. Checking the *wrongly merged*
records too is the point: they are where a 5' shift would hide, and they are clean.

With full-length mates, `L = shift + len2` makes an emitted length equal to the true
fragment length **equivalent** to inferring the true shift, so the C2 count is exactly
the count of wrong inferences — and a pair with no true overlap can never reach the true
length by accident, since it cannot emit more than `2 × readlen − 1` bases.

C4 — *an unmerged pair tiles its fragment exactly once* — is the contract the trim band
exists to hold, and §5 is entirely about it. Its hard half (nothing may move a 5' end,
and R1 is never shortened) is exact; its soft half is a bases-count trade, priced there.

One smaller finding, not a boundary violation: **504 pairs vanish entirely** (fastp:
319). These are wrong merges that produced a record under `--min-read-length 40` and were
filtered away, so the fragment leaves the corpus rather than being emitted as a pair.
They are invisible in the output file; the scorer recovers them from the pairs that
produced no records at all, and they are included in the chimera counts above. This is
`--min-read-length` behaving as documented, on input where the merge was wrong.

---

## 5. The trim band (contract C4) — what the *unmerged* pairs cost

A pair that does not merge is still encoded, so this half of the tool reaches the corpus
just as directly as merging does. `zna merge` splits the redundant overlap between the
two mates' **3'** ends, leaving R1 and R2 tiling the fragment exactly once at equal
length. Two ways to get it wrong, and they are opposite, so they are scored over
**every** kept pair rather than only the overlapping ones — an analysis restricted to
pairs that share sequence cannot see a trim applied to a pair that shares none.

**The trim was asymmetric until 2026-08-12** — the whole overlap came off R2, leaving one
full-length read beside a short one. The overlap sits at the 3' end of *both* mates, so
splitting it tiles the fragment exactly as well, and every number in this section is
**identical under both rules**: same merges, same trims, same duplicated and deleted
bases, same zero boundary violations. What changed is only where the cut falls:

| trimmed pairs, 2×150 | asymmetric | symmetric |
|---|---:|---:|
| mean `\|len(R1) − len(R2)\|` | 10.8 (the overlap) | **0.51** |
| max `\|len(R1) − len(R2)\|` | 110 | **1** |

The residual 1 is an odd overlap, whose extra base stays on R1.

| | zna merge | fastp |
|---|---:|---:|
| kept pairs — no true overlap / overlapping / read-through | 380,525 / 36,609 / 0 | 382,875 / 76,354 / 0 |
| shortest overlap the threshold can act on | 5 bases | – |
| trimmed **exactly** the true overlap | 25,971 | 0 |
| over-trimmed / under-trimmed | 222 / 3 | 0 / 0 |
| not trimmed, overlap below the threshold *(by design)* | 9,963 | 10,120 |
| not trimmed, overlap above the threshold | 450 | 66,234 |
| **false trims** (a pair sharing nothing, trimmed anyway) | 2,733 | 0 |
| an emitted read longer than its input | **0** | 0 |

### Exact-trim rate by true overlap, over each tool's own kept pairs

| true overlap | zna kept | zna trimmed exactly |
|---|---:|---:|
| 1–4 | 10,095 | 0 (0.0%) — below the threshold, by design |
| 5–9 | 12,865 | 12,352 (96.0%) |
| 10–14 | 12,640 | 12,613 (99.8%) |
| 15–19 | 949 | 946 (99.7%) |
| 20–29 | 60 | 60 (100.0%) |

The band is bounded on both sides by the thresholds and the curve shows both edges: below
5 bases nothing can clear `--threshold-trim 8`, and above ~15 the pair merges instead, so
only a handful reach the trim path at all. The 4.0% miss at 5–9 is the same arithmetic as
the merge cliff — a 5-base overlap scores 9.93 bits and one mismatch costs 8.21, so a
single error puts it under.

### In bases, which is what a model actually sees

| | zna merge | fastp |
|---|---:|---:|
| bases emitted in unmerged pairs | 124,826,899 | 137,768,700 |
| duplicated fragment bases, **if nothing were trimmed** | 283,838 | 1,181,978 |
| duplicated bases removed correctly | **255,764** | 0 |
| duplicated bases surviving into the corpus | 28,074 | 1,181,978 |
| — per 10,000 unmerged bases | **2.25** | 85.79 |
| real fragment bases deleted | 57,537 | 0 |
| **benefit ratio** (duplicates removed per base deleted) | **4.45** | – |

fastp is not a competitor here: it has no redundant-overlap trim at all — it trims
*adapter*, not shared sequence — so its column is the size of the problem, not a rival
result. Leaving 85.8 duplicated bases per 10k against 2.25 is the difference the trim
band makes to a corpus.

### Is `--threshold-trim 8` the right value? Now measurable

Sweeping the shipped tool over the same input. The merge count is **identical** in all
four runs (582,866), so this isolates the trim band cleanly:

| `--threshold-trim` | exact trims | false trims | duplicates removed | real bases deleted | duplicates left | ratio |
|---:|---:|---:|---:|---:|---:|---:|
| **8** (default) | 25,971 | 2,733 | 255,764 | 57,537 | 28,074 | 4.45 |
| 12 | 20,847 | 1,177 | 226,272 | 39,487 | 57,566 | 5.73 |
| 16 | 15,454 | 690 | 184,232 | 27,792 | 99,606 | 6.63 |
| 20 | 10,231 | 372 | 132,517 | 16,124 | 151,321 | 8.22 |

Raising the threshold improves the *ratio* and shrinks the *benefit*, which is the shape
you would expect and not by itself an argument either way. The decision needs the two
harms added: duplicated-left plus deleted is **85,611 / 97,053 / 127,398 / 167,445** at
8 / 12 / 16 / 20. **The default minimises it**, and it keeps doing so unless a deleted
base is judged more than **1.63×** as harmful as a duplicated one — which is the
threshold khorana would have to cross to justify a change. Given that a duplicated base
is *false* evidence (the same molecule counted twice with nothing marking it) while a
deleted base is merely absent, the asymmetry runs the other way if anything.

### What the false trims actually are

The 55,565 deleted bases split cleanly into the designed cost and the same repeat
phenomenon as the chimeras:

| | trims | bases | mean | error-free |
|---|---:|---:|---:|---:|
| chance-compatible (found overlap ≤ 14 bases) | 1,724 | 12,167 | 7.1 | 55% |
| repeat-driven (found overlap > 14 bases) | 1,009 | 43,398 | 43.0 | 56% |

The first row is what an 8-bit threshold is *for*: it buys short trims at a chance rate
of ~0.45% of no-overlap pairs, which is the price the redesign knowingly paid. The second
row is **78% of the lost sequence** and is not a threshold artifact — hotspots are
`chr11:50–55`, `chrX:60–65`, `chr17:20–25`, `chr15:15–20`, `chr20:25–30` Mb, every one
pericentromeric, exactly as in §2.

One correction this forced: `args.py` and the redesign notes both said a wrong
trim "costs a few bases". It costs a median of 9, a mean of 20, and up to 110 — a third
of a read. The asymmetry that justifies a low `T_trim` still holds, but it is now stated
with a number instead of an adjective.

**No contract is violated by any of this.** A trim only ever removes from a 3' end, so
base 0 of every mate stays a true fragment boundary — verified directly, 0 of 1,416,630
records — and no emitted read is longer than its input. A false trim loses signal; it
never invents any.

### Two defects the symmetric change introduced, and the benchmark caught

Worth recording, because neither was visible from the test suite and both were found by
re-running this benchmark against the modified tool.

1. **237 corrupted 5' ends.** Splitting the overlap means R2 keeps part of it, so the
   consensus has to be written into R2 as well. Applied unconditionally, that rewrites
   R2's *5'* bases on a **read-through** geometry, where the overlap lands there — on an
   inference too weak to merge on. All 237 were pairs with no true overlap, 228 of them
   at a negative inferred shift — and, re-measured later, **9 at a non-negative one**.
   That matters because the fix was originally justified as "the trim branch requires
   shift >= 0, which excludes that geometry". It does not: the bound is
   `min written R2 index = max(0, len2 - len1 + shift)`, which reaches index 0 at
   non-negative shifts once the mates differ in length. What actually keeps trims clear
   of R2's 5' end is `trim_is_allowed`, which forces `shift >= min_read_length`. The
   same 9 records are the ones item 2 below attributes to dropping that guard. Fixed by writing R2 only on the trim branch, which
   requires `shift ≥ 0` and therefore cannot reach that geometry. Merged and kept records
   are byte-identical to before.
2. **9 more, plus 133 extra false trims and 17,214 extra deleted bases**, from a guard
   that looked redundant and was not. The old rule `len2 − olen ≥ min_read_length`
   reads as "don't trim R2 below the length filter", but it also **capped the overlap the
   trim band may act on**: 145 clean bases score 288 bits and would merge outright, so an
   overlap that long arriving in the 8–28 bit band is carrying a pile of mismatches and
   is almost certainly spurious. Restated symmetrically as *each mate must reach at least
   `min_read_length` past the other's 3' end*, it recovers both properties — the emitted
   lengths are bounded below by the split's own clamp, and near-total overlaps are
   refused. `test_a_near_total_overlap_in_the_trim_band_is_refused` pins it.

---

## 6. Head to head with fastp: matching its false-positive rate, and beating it

The two tools sit at different points on one trade-off, so "which is better" is only
answerable at a matched operating point. Sweeping `--threshold-merge` with everything
else held fixed puts them on the same curve. **Merge is the only thing that moves here**
— `--threshold-trim` stays at 8, so the trim band absorbs whatever stops merging.

| setting | chimera rate | wrong merges | sensitivity ≥15 | missed | FP+FN | precision | exact |
|---|---:|---:|---:|---:|---:|---:|---:|
| **`--threshold-merge 28`** (default) | 1.231% | 5,591 | **99.83%** | 1,012 | **6,603** | 99.041% | 86.59% |
| `--threshold-merge 34` | 1.060% | 4,883 | 98.50% | 8,701 | 13,584 | 99.150% | 87.02% |
| `--threshold-merge 40` | 0.924% | 4,321 | 97.13% | 16,637 | 20,958 | 99.237% | 87.46% |
| `--threshold-merge 50` | 0.738% | 3,533 | 94.89% | 29,557 | 33,090 | 99.360% | 88.17% |
| **`--threshold-merge 60`** | **0.597%** | 2,933 | 92.63% | 42,683 | 45,616 | 99.455% | **88.87%** |
| `--threshold-merge 80` | 0.385% | 2,023 | 88.08% | 68,978 | 71,001 | 99.604% | 90.19% |
| `--threshold-merge 100` | 0.245% | 1,403 | 83.57% | 95,109 | 96,512 | 99.711% | 91.39% |
| **fastp 1.1.0 defaults** | **0.621%** | 3,502 | 92.98% | 40,643 | 44,145 | 99.352% | 85.90% |

*chimera rate* — merged, of the 385,268 pairs with no true overlap: the false-positive
rate. *wrong merges* — chimeras plus merges at the wrong shift, i.e. C2 violations.
*missed* — pairs with a true overlap ≥ 15 that did not merge, `578,841 × (1 − sens)`.
*precision* — merges that are geometrically right. *exact* — merged records equal to the
true fragment base for base.

### To match fastp's false-positive rate: `--threshold-merge 60`

0.597% against fastp's 0.621%, and sensitivity 92.63% against 92.98% — the same
operating point to within a rounding error. **The agreement is mechanical, not luck:**
60 bits is `ceil(60 / 1.9855) = 31` clean bases, and fastp's `--overlap_len_require`
defaults to **30**. Read the other way, this prices fastp's headline knob on a
calibrated scale — *fastp's 30-base floor is worth about 60 bits of evidence.*

At that matched point zna is still the more accurate tool, on the two axes the threshold
does not control:

- **reconstruction** — 88.87% of merged records equal the true fragment against 85.90%,
  because the posterior consensus recovers 90.35% of recoverable overlap errors against
  fastp's 74.07% (§3);
- **precision** — 99.455% against 99.352%, because the decision is an `argmax` over all
  shifts rather than a first-accept, so a spurious short hit cannot preempt the real
  offset.

### For best overall accuracy: keep the default 28

It minimises false positives plus false negatives, and not marginally: **6,603 errors
per million pairs against 44,145** at the fastp-equivalent setting. The curve is steeply
asymmetric because sensitivity falls much faster than the chimera rate does.

| raising the threshold | wrong merges prevented | missed merges caused | exchange rate |
|---|---:|---:|---:|
| 28 → 34 | 708 | 7,689 | 10.9x |
| 28 → 40 | 1,270 | 15,625 | 12.3x |
| 28 → 60 | 2,658 | 41,671 | 15.7x |
| 28 → 100 | 4,188 | 94,097 | 22.5x |

So raising it pays only if a wrong merge costs more than about **11x** what a missed
merge costs. Weigh that knowing what each actually is: a wrong merge injects a fragment
that does not exist, with a false 3' boundary; a missed merge emits the same molecule as
two correctly-bounded reads with the redundant overlap trimmed off — no false
information, just less consolidation. For a corpus whose whole point is fragment-end
supervision the ratio is not obviously below 11, which is why this is stated as a
trade-off with a number rather than a recommendation to raise it.

### Tuning cannot reach zero, and that is a property of the genome

At 100 bits — 3.6x the default — **1,403 wrong merges per million remain**. §2 is the
explanation: every one is a fragment whose two ends are genuinely homologous, the scan
never picks a lower-scoring alignment than the true one (0 of 848 checked), and 56% carry
no sequencing error at all. Anyone quoting a chimera rate from an overlap merger without
saying what genome the fragments came from is quoting the input, not the tool.

---

## 7. Throughput

| tool | wall s | µs/pair |
|---|---:|---:|
| `zna merge --threads 4` | 1.53 | **1.528** |
| `fastp -w 4` | 2.83 | 2.826 |

Same input, same thread count, one Apple M3 Max (16 cores, aarch64). fastp also writes
three output files and computes its own statistics, so this compares the two *commands*,
not their merge kernels.

---

## 8. What was fixed, and what was not

**Fixed before the run**, because the benchmark's own histograms had to be trustworthy at
any read length: the three histograms in `fastq_chunk.hpp` were fixed `uint32_t[1025]`
arrays with every index clamped to the last bin, so past 1024 bp the length and insert
distributions silently aggregated. They now grow with the `Scratch` arena
(`2 × cap + 1` bins, resized O(log L) times per chunk, not per record), so indexing is
still one instruction and there is no cap. `insert_size_censoring.floor` — a second copy
of `params.min_read_length` — was dropped.

**Corrected by the run:** the justification for `--threshold-trim 8` in `args.py` and
the redesign notes said a wrong trim "costs a few
bases". Measured, it costs a median of 9, a mean of 20 and up to 110. The asymmetry it
was arguing for survives — the band still removes 4.45 duplicated bases per base it
deletes, and 8 bits minimises the two harms added — but the claim now carries a number
instead of an adjective.

**Not fixed, deliberately:** the 5,591 wrong merges. They are not a defect in the scan
(0/848 argmax failures), not a threshold artifact (25% survive a 100-bit threshold), and
not noise (56% are error-free). They are fragments whose two ends are genuinely
homologous, drawn uniformly from a genome that is ~50% repeat. The residue is a
*corpus* question, not a merge-tool one, and the production audit already reached the
same conclusion from the other direction: gating on shift
ambiguity costs more than it saves. The same phenomenon accounts for 78% of the bases
lost to false trims, which is why raising `--threshold-trim` does not buy what it looks
like it should.

Three things a reader should carry away rather than a number:

1. **Uniform fragment lengths over the whole genome is a stress test, not a library.**
   38.5% of these pairs have no overlap at all, and pericentromeric satellite appears at
   its genomic frequency. A real RNA library has neither property, and the production
   figures measured on real libraries remain the ones to quote for rates.
2. **The chimera rate is a property of the input, not of the tool.** Both tools find the
   same class of pair, with different rules and different thresholds. If chimeras matter
   for a downstream model, they have to be addressed with the genome in hand — which the
   merger, by design, does not have.
3. **Trimming is half the tool, and it reaches the corpus by the same road.** Merging is
   where the interesting failures are, but most pairs in a real library do not merge, and
   what happens to those is 2.25 duplicated bases per 10k here against 85.79 with no trim
   at all. Any change to `--threshold-trim` should be argued on the sweep in §5, not on
   the ratio alone — the ratio improves monotonically while the benefit falls off a
   cliff.
