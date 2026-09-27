# zna merge 0.6: accuracy by derivation, not tuning

**Status:** **implemented** as zna 0.6.0 (branch `merge-policy-0.6`; release pending).
Design agreed 2026-09-26 (error rate revised the same day: a documented, user-set
parameter, not an estimate — §3); Stages A–C done the same week, and the §8 acceptance
criteria, fixed before Stage C, held at full scale except where the adapter declaration
was not honest (§11). Results and residual limits are §11; the algorithm as built is
`docs/METHODS.md` §1–§3, and the measurements `docs/MERGE_BENCHMARK_RESULTS.md` §9.
**Supersedes:** the 2026-09-25 version of this file (a tuned policy with read-through
concordance constants, a trim mismatch cap and displaced-winner clauses), which was
rejected as overfit-shaped. Its measurements remain valid evidence and are cited below.
**Motivation:** khorana's chr22 merge review (2026-09-25; kept with khorana's evidence, not in this repository), and a reviewer's critique of this plan's first draft.
**Evidence:** kept outside the repository with khorana's simulation data (§10).
**Compatibility:** none. ZNA is beta; 0.6.0 breaks behaviour and format freely, every
existing `.zna` is regenerated, and old releases are reproduced with a pinned environment.

---

## 1. What changes, in one screen

| | 0.5.3 | 0.6 |
|---|---|---|
| merge threshold | `--threshold-merge 28` bits, fixed | **derived per pair**: `T = log2((len1+len2−1)/α)`, `--alpha` (default 1e-6) |
| sequencing error model | `err_rate = 0.01`, hidden constant | **`--error-rate`** (default 0.01): a documented, user-set parameter; zna reports the rate its detected overlaps actually show and warns when, at that rate, the gate is expected to refuse more than 0.1% of true overlaps |
| wrong-but-strong alignments (repeats) | merged | **refused** when the mismatch count is implausible under the error model at the same α |
| read-through | always inferred | **`--adapter-trimmed`**: a declaration zna trusts (read-through becomes impossible), checked on the first 100k pairs |
| trim band (`--threshold-trim 8`) | trims 5–28 bp overlaps | **removed**: a pair is merged or kept whole |
| consensus | merge and trim paths | merge path only (unchanged) |
| file provenance | `writer_version`, `merged_in_process` | + a **merge record** written at merge time and carried through shuffle and re-encode |

User-facing parameters: `α` (chance-merge tolerance) and `--error-rate` (the error model).
Declaration: `--adapter-trimmed`. Everything else is derived or a documented definition.

## 2. The policy

For each pair with read lengths `len1`, `len2`:

```
N      = len1 + len2 − 1                          # candidate shifts on the signed axis
T      = log2(N / α)                              # bits; exact fixed point, table over N
e      = --error-rate (default 0.01; §3)           # one number per run
weights: match = log2((1−e)/0.25), mismatch = log2(0.75/e)     # unchanged formulas
dfit[n]= max d with P(Binom(n, e) ≥ d) ≥ α        # exact integers, table over n

eligible shifts:
    all of s ∈ [−(len2−1), len1−1]                             (default)
    only s ≥ max(0, len1 − len2), i.e. L̂ = s+len2 ≥ max(len1,len2)   (--adapter-trimmed)

W = argmax over eligible shifts of the score, ties → smaller s   (0.5.3's total order)
    with floor T: nothing ≥ T  ⇒  no overlap
d_inf  = mismatches at W minus positions where exactly one base is N
if W exists and d_inf > dfit[n_W]:  no overlap          # "implausible": abstain, never re-place

W exists → MERGE (0.5.3's construction and posterior consensus, R1 wins ties);
           a merged record shorter than --min-read-length is dropped (as today)
otherwise → KEEP both mates unchanged (N policy applied; pair emitted all-or-nothing)
```

**Why each piece is there, and where it came from**

- **T from α.** Under unrelated sequence a shift reaches T bits with probability ≤ 2^−T
  (E[LR] = 1); a union over the N shifts bounds a chance merge per pair by N·2^−T = α.
  The old 28 bits was this at 2×150 (N = 299). Per pair, 2×50 → 26.6 bits, 2×300 → 29.2;
  the shortest mergeable perfect overlap is 14–15 bp. α bounds merges of *unrelated*
  sequence only; measured chance merges sit far below it (0 in 40,000 random pairs).
- **The plausibility gate uses the same α.** Only one shift can be the true one, so under
  the error model a true overlap is refused with probability ≤ α. One tolerance, two
  errors: at most α chance merges and at most α true overlaps refused, per pair. It
  catches the report's dominant failure — a long divergent repeat outscoring a short true
  overlap (C01: 19 mismatches in 122 bases where `dfit[122] = 9` at e = 0.01).
- **Abstain, don't re-place.** An implausible best alignment marks a repetitive context;
  searching for a runner-up re-placed 37–38% of caught wrong merges onto *another* wrong
  shift on the gene-disjoint holdout (9% on dev). Abstaining removes that class by
  construction. Cost: pairs whose true overlap was outscored are kept whole instead of
  re-merged (measured on dev: the gate removes ~10% more wrong merges than re-placing
  and adds no lost fragments; it forgoes 145–271 re-merges per 500k–1M pairs).
- **N carries no information** about agreement, so it does not count against plausibility.
  Without this, reads with N runs had ~1,500 correct merges per 200k refused (the only
  correctness regression the robustness review found).
- **Read-through as a declaration.** After adapter trimming (or for reads clipped to the
  molecule, as khorana's simulator produces), no read extends past its fragment, so an
  alignment implying one does is impossible. No threshold. On input where the
  declaration is true it measured better than every tuned read-through rule (tx-fastp:
  −715 wrong merges and +332 correct merges vs −632 / +249). Without the declaration zna
  keeps 0.5.3's behaviour, which *is* overlap-based adapter removal: a read contains
  adapter only when its insert is shorter than the read, so the mates overlap fully and
  the merged record `[0, L)` excludes the adapter. Declaring it wrongly on raw reads is
  catastrophic (−231k correct merges on hg38-noisy), hence §4's verification.
- **No trim band.** khorana now trains on one randomly chosen mate of an unmerged pair and
  wants it whole; trimming removed ~0.55 duplicated bases per pair (the review's §3: 550,656 redundant copies per 1M chr22 pairs) while causing every
  wrong trim (chr22: 1,454 per 1M pairs, 23,823 unique positions deleted, 1,755
  substitutions). Other uses of overlap trimming (variant calling, coverage) are handled
  after alignment by the tools that need them. Removing it deletes `--threshold-trim`,
  the balanced split, the trim guard, trim-path consensus and `PROV_TRIMMED`.

## 3. The error rate: a parameter, reported against the data

`--error-rate e` (default 0.01) is the expected fraction of disagreeing bases between the
two mates where they truly overlap — about twice the per-base sequencing error, since
either read can be wrong. It is the same quantity 0.5.3 hid as `err_rate`, now documented
and settable. It enters two places, with different roles:

- **The score** (`match = log2((1−e)/0.25)`, `mismatch = log2(0.75/e)`) — how strongly a
  mismatch counts against an overlap. The α bound on chance merges does **not** depend on
  `e`: under unrelated sequence the likelihood ratio averages exactly 1 whatever `e` the
  overlap hypothesis assumes. So `e` here sets strictness, not validity. A larger `e` is
  more tolerant of mismatches — sequencing errors and repeat divergence alike.
- **The plausibility gate** (`dfit`) — how many mismatches a true overlap can carry. Its
  guarantee (a true overlap is refused with probability ≤ α) holds only if `e` is at least
  the library's real disagreement rate. If it is set too low, true overlaps are refused and
  kept whole (never merged wrongly); the closed form is in §8.

**Why a parameter and not an estimate.** An estimate was implemented and measured
(`patches/stageA-with-estimator.patch` in the evidence directory, §10). Feeding a
library's own rate into the score made noisy libraries admit *more* divergent repeats and
short false read-throughs (3′-ramp sets: wrong merges 560 → 741 with `e` = 0.035, vs
560 → 240 at a fixed 1%; lost fragments rose on 1%-error hg38), and clean ones reject a few
short true overlaps with two or more errors. The estimator also needed a buffer, a prior
and a sorted-input caveat. A fixed, documented `e` is simpler and was safer on the panel;
the data are used to *check* it instead.

**How the data check it.** Over every overlap the scan detects (score ≥ T, before the
gate, informative positions only), zna reports `detected_overlap_mismatch_rate`. It is a
check, not an estimate, and it errs both ways: repeats inflate it on clean libraries
(hg38 at full scale: 0.0051 detected vs 0.0041 true), while on degraded ones the most
divergent true overlaps never reach T, so it reads low (3′-ramp: 0.0355 vs 0.0385). From it and the
histogram of detected overlap lengths zna computes what the rate *costs* — the closed
form of §8, `Σₙ count[n]·P(Binom(n, rate) > dfit[n]) / Σₙ count[n]`, reported as
`expected_refused_true_overlap_fraction` — and warns when that exceeds 0.1% (a Stage C
decision: the first trigger, "detected rate above `--error-rate`", fired on 12 of the 46
dev benches, ten of them expecting under 0.02% refused; the new one fires only on the
5% 3′-ramp set, at 0.40%). The warning gives the refused count and a suggested value,
and states the trade: raising `e` also softens the score, so more short or divergent
overlaps merge, false ones included (on the 3′-ramp sets following the suggestion took
wrong merges 240 → 745). Production libraries
measured ~0.009 on 0.5.3's comparable statistic (hulkrna); poor or 3′-degraded libraries
(2–5%) whose runs refuse many pairs are the ones to weigh raising it.

Exact derivations: `e` is parsed from its decimal string into an exact rational; weights,
`T(N)` and `dfit` are computed once per run in Python with `decimal` (software, correctly
rounded) and `Fraction`, and cross into both backends as integers — no libm anywhere in a
decision.

## 4. Diagnostics and warnings

All JSON values finite and of a fixed type (hulkrna's cohort gather rejects `Infinity`
and type changes). Nothing here changes a decision.

| key | meaning | warning |
|---|---|---|
| `error_rate` | the `e` used (`--error-rate`) | — |
| `detected_overlap_mismatch_rate`, `detected_overlap_bases` | informative mismatches / positions over every detected overlap, before the gate | — |
| `expected_refused_true_overlap_fraction`, `detected_overlap_length_histogram` | the share of detected overlaps the gate is expected to refuse were they true and disagreeing at the detected rate (§8's closed form over the histogram) | `> 0.001`: true overlaps are being refused; states the trade and suggests a value |
| `adapter_trimmed`, `readthrough_check_pairs`, `readthrough_check_strong_fraction` | over the first 100,000 pairs (a diagnostic scan run alongside the real one, deterministic in input order): fraction whose unrestricted best alignment is a read-through scoring ≥ T | under `--adapter-trimmed`, `> 0.01`: the reads appear to contain adapter (measured: ~0.1% honest, 2–23% false) |
| `implausible_refused` | pairs refused by the gate | — |
| `overlap_mismatch_rate` | over merged overlaps (post-gate) | documented as post-admission |

## 5. Provenance

`encode --merge-pairs` writes a `merge` object into the prologue: `policy` (`"zna-merge-0.6"`),
`zna_version` at merge time, `alpha`, `error_rate`, `adapter_trimmed`, `min_read_length`,
`npolicy`. `shuffle` and re-encode copy it
unchanged (as `merged_in_process` already is); a file whose origin is unknown carries no
`merge` object, and consumers treat "absent" as unknown, never as a policy. `zna inspect`
prints it. Round-trip tests: encode → shuffle → re-encode preserves it byte for byte.

## 6. API

`find_overlap` returns the **authoritative decision**: the alignment actually used and a
verdict (`merge`, `none`, `implausible`). The unrestricted scan stays available under an
explicit diagnostic name. `pairs.process_pair` and the chunk functions take the derived
parameters (`T` table, `dfit` table, weights, contract) as integers/bytes, like the score
weights today.

## 7. Implementation (ordered; each stage leaves the suite green in its configuration)

*All three stages are done. The stage notes below are the plan as written; where the
build departed from it, §3 and §11 say so — chiefly the `--error-rate` warning, whose
trigger changed in Stage C from "detected rate above `--error-rate`" to the expected
refused share.*

**Stage A — reference implementation (Python), CLI, provenance, tests.**
`params.py` (α, e, exact `decimal`/`Fraction` derivations, `T` and `dfit` tables grown by
doubling); `_pymerge.py` (drop the trim path; per-pair floor `T`; contract range; gate
with informative mismatches; counters incl. detected-overlap totals; the read-through
diagnostic scan for the first 100k pairs); `cli.py` and `encode_stream.py` (new stats keys
and warnings, shared between `zna merge` and `encode --merge-pairs`);
`args.py` (`--alpha`, `--error-rate`, `--adapter-trimmed`; remove both thresholds);
`overlap.py`, `pairs.py`, `backend.py` (API above); `core.py`, `_shuffle.py`, `cli.py`
(prologue merge record, carried through). Tests: delete trim tests; add exact goldens for
`T(N)` and `dfit`, the diagnostics (warning fires when detected disagreement exceeds
`--error-rate`; read-through check; independence of threads and chunk size), contract ranges (`len1 < len2`, `len1 > len2`), N-informative
gating, the report's 19 cases as fixtures, provenance round trips. Run with
`--backend python` until Stage B.

**Stage B — compiled kernel.** `merge_core.hpp` (floor from the `T` table, scan range
under the contract, post-scan gate — a single table lookup, no new loop state),
`fastq_chunk.hpp`, `_accel.cpp`; byte-identical to Stage A on the full suite, the panel and
fuzzed pairs; all three test configurations (see the merge test-environment notes).
Timing on aarch64 against 0.5.3 (the floor rises from 8 bits to ~28, so the scan prunes
earlier) and later on Linux/x86.

**Stage C — qualification and documentation.** §8's panel evaluation, full-scale
confirmation, the sealed holdout-2 panel once; METHODS, README, CHANGELOG (0.6.0,
breaking), ROADMAP, MERGE_BENCHMARK_RESULTS; handoff docs for hulkrna (flags, config, env
pin, the N-policy question) and khorana (require the prologue merge record).

## 8. How it is tested: three layers, seconds each

1. **Exact tests (no data).** `T(N, α)`, `dfit(n, e, α)` goldens and monotonicity; the
   gate's sensitivity as a closed form: P(true overlap of length n refused | true rate
   e′) = P(Binom(n, e′) > dfit[n]). Example (rule built for e = 0.01): at e′ = 1% it is
   ~1e-6; at 3%, 0.07–0.6% (n = 50–150); at 5%, 1–13% — which is why a run is warned
   when this, summed over its detected overlaps at their detected rate, passes 0.1%.
2. **The truth panel.** Mined in 22 s of wall time (6 processes) by zna's compiled backend from 14.1M simulated pairs
   over 46 substrates: every pair 0.5.3 gets wrong (57,595 wrong merges, 24,264 wrong
   trims, 2,077 lost fragments, 20,244 wrong alignments kept, N-trim drops) plus a
   weighted stratified sample of correct pairs — 245k pairs, 62 MB; plus a sealed 34k-pair
   panel from holdout-2 (genes never sampled before, and an unseen hg38 draw). A change is
   judged by a **pair-level ledger** against truth: which pairs changed outcome, from what
   to what. Scripts: `scripts/merge_bench/mine_panel.py` (mines, under a pinned 0.5.3) and
   `panel_eval.py` (evaluates the current tree); data under §10.
3. **Full-scale confirmation, at release only.** Regenerate substrates from their recipes
   (deterministic seeds), run the compiled merger (seconds per million pairs), delete.

**Acceptance (fixed before Stage C runs; every bench, vs 0.5.3 on the same pairs):**
wrong merges retained and lost fragments do not increase on any bench, and fall overall;
baseline-correct merges forgone ≤ 0.01% of baseline correct merges on every bench where
the error model holds; wrong trims and kept-mate substitutions are zero by construction;
adapter bases emitted under an honest declaration do not increase; overall affected
pairs fall; compiled per-pair time within +5% of 0.5.3. Report class *and* overall
counts, per-bench locus concentration of the residual, and the holdout-2 result once.

## 9. Known limits, stated up front

- **Near-identical repeats** (a perfect 15-bp repeat, C07; 7/64 mismatches, C14) are
  plausible under any error model and remain wrong merges. On the gene-disjoint holdout
  they were 52–56% of 0.5.3's wrong merges. Distinguishing them needs information the
  merger does not have.
- **Raw reads without the declaration** keep 0.5.3's false read-through merges and the
  resulting lost fragments (hg38-noisy: 495 per 1M, 0.5.3: 505). The remedy is to trim and declare.
- **`--error-rate` set below a library's real disagreement** refuses true overlaps (kept
  whole, never merged wrongly); the warning in §4 says so and suggests a value.
- **Quality scores** play no part in the decision (only in the merge consensus).
- **The declaration is only as honest as the trimmer** (found in Stage C, §11). hulkrna's
  fastp pass leaves 1–3 bp of adapter on 42 per 500k transcriptome pairs (35 on R1
  alone, its `--cut_tail` class; 7 on both mates); the same pass without `--cut_tail`
  leaves the 7. Declared, those pairs are kept with their adapter or (4 per 500k) merged
  at a wrong length. zna does not trim adapters — that stays fastp's job — and its
  read-through check (§4) is the per-library signal, though at this rate (≤ 0.01%) it
  cannot fire.
- **The residual concentrates.** Wrong merges roughly halve, and what remains sits on a
  few repeat-rich transcripts: on the clean transcriptome set one transcript carries 200
  of 589 (0.5.3: 180 loci, top five 40%; 0.6: 99 loci, top five 55%). A consumer that
  weights by locus should know it.
- **Very long reads pay table building once per run**: a 10 kb read costs ~1.6 s before
  the first pair (the exact floor table, in `decimal`); nothing at Illumina lengths.

## 10. Evidence

The evidence directory (outside the repository): `panel/` (the
truth panel, `mine.py`, per-substrate summaries; `sealed/` is holdout-2 — do not inspect
before Stage C), `bench/` (the exploration harness: `candidate.py` is byte-identical to
0.5.3, `evaluate.py`, the rule-family modules and their result summaries), `wf3/` and
`wf4/` (combination, adversarial reviews, first holdout), `data-recipes/` and
`hold2-recipe/` (regenerate any substrate), `drafts-unverified/` (an interrupted
prototype of this policy — reference only), `patches/stageA-with-estimator.patch` (the
Stage A implementation with the error-rate estimator, before it was dropped),
`fastp_probe/` (the synthetic fastp residual-adapter probe: `gen.py`, `run.sh`,
`run2.sh`, `score.py`), `fastp_single_pass/` (fastp variants at full scale on
tx-dev-noisy, fastp 1.1.0 and 1.3.6: `run.sh`, `run136.sh`, `score.py`, logs and each
variant's fastp JSON), `stageC-results/` (the full-scale qualification and the one
sealed-holdout run: per-bench JSON, pair-level ledgers, timings, the scripts that made
them, and the regeneration check; `tx-dev-fastp3.*` is the recommended single fastp
pass, measured afterwards — its README addendum).

## 11. Results (Stage C, 2026-09-26)

**Full scale.** Seven benches regenerated byte-identically from their recipes (chr22
pilot/val/train, the three transcriptome dev sets, hg38) plus a new two-pass-fastp set:
eight benches, 4.11M pairs, every pair classified against truth under 0.5.3 and 0.6,
each bench under its primary contract (`docs/MERGE_BENCHMARK_RESULTS.md` §9 has every
row):

| | 0.5.3 | 0.6 |
|---|---:|---:|
| wrong merges retained | 11,875 | 5,151 |
| lost fragments | 719 | 544 |
| wrong trims | 7,225 | 0 |
| kept-mate substitutions | 14,366 | 0 |
| affected pairs | 19,819 | 5,695 |
| correct merges | 2,937,253 | 2,937,748 (7 forgone, 502 gained) |
| gate refusals at the true length | — | 0 of 5,632 |
| compiled CPU, `merge --threads 1` | — | 9–29% lower |

**Acceptance.** Every criterion holds on every bench and contract except on *declared*
fastp output, where the declaration is not fully honest: hulkrna's fastp pass (with
`--cut_tail`) forgoes 42 correct merges (0.0111%, over the 0.01% limit) and emits 47
adapter bases; two-pass fastp — whose primary contract is declared — forgoes 7 (0.0019%,
within the limit) but emits 8 adapter bases in 3 pairs where 0.5.3 emitted none.
Undeclared, both pass everything. The error model held on every bench (true
disagreement ≤ 0.0049), and no run warned. Measured afterwards with the same scripts,
the recommended single pass without `--cut_tail` (tx-dev-fastp3) matches two-pass:
declared, the same 7 pairs forgone (0.0019%) and 8 adapter bases in 3 pairs, every
other criterion held, time included (−13.8% CPU); undeclared, everything passes.

**The sealed holdout-2, once.** Wrong merges 10,673 → 5,359, lost 575 → 518, wrong
trims → 0, no bench warns; 0 forgone except tx-hold2-fastp (hulkrna's pass with
`--cut_tail`, declared) at 63 weighted (0.0083%). The panel cannot resolve the 0.01%
criterion — one sampled correct-merge row weighs 0.05–0.065% of baseline — so that
criterion rests on full scale.

**Decisions made in Stage C.**

- *The `--error-rate` warning* fires on `expected_refused_true_overlap_fraction > 0.001`
  (§3, §4), replacing "detected rate above `--error-rate`", which fired on 12 of the 46
  panel benches. Its message states the trade (240 wrong merges at 0.01, 745 at 0.036 on
  the 3′-ramp sets).
- *fastp.* Leftover adapter is fastp's to remove, not zna's (division of labour). The
  synthetic probe traced hulkrna's fastp residual to `--cut_tail`, which runs before
  fastp's overlap-based adapter detection and, by trimming R2's last (adapter) base,
  hides R1's 1–2 bp overhang. **Full scale** (tx-dev-noisy, fastp 1.1.0 and 1.3.6
  identical): 42 pairs keep adapter with `--cut_tail`, 7 without it in one pass, 7 with
  two passes (adapters, then `-A --cut_tail`); `--overlap_diff_limit` made it worse
  (10: 53 with `--cut_tail`, 18 without; 1: 46; 0: 368 — against 7 at the default 5).
  The last 7 keep 1–3 bp on both mates and are missed by fastp's overlap detection
  itself (hypothesis, not verified: low-complexity inserts on which fastp accepts a
  full-length overlap). **Recommendation to hulkrna: one fastp pass
  without `--cut_tail`, declared** — 664 wrong merges against 750 undeclared, 0 lost
  fragments against 15, 7 correct merges forgone (0.0019%) and 8 adapter bases emitted
  (tx-dev-fastp3, `docs/MERGE_BENCHMARK_RESULTS.md` §9). A second pass only if hulkrna
  wants 3' quality trimming: same adapter result, and here it removed 0.9% of sequencing
  errors. Output of a pass with `--cut_tail` is documented as not honest for
  `--adapter-trimmed`.
- *File format.* `_VERSION` stays 3 and `PROLOGUE_SCHEMA` 1: the merge record is an
  optional prologue key nothing gates on, absence already means unknown, and the retired
  provenance bit reaches a file only as a user-declared label value.

**Downstream.** hulkrna drops `--threshold-merge`/`--threshold-trim`, drops fastp's
`--cut_tail` and passes `--adapter-trimmed`, and its cohort gather meets the new JSON
keys (the CHANGELOG lists them; all type-stable, 0 / 0.0 / `{}` on an empty input). khorana requires the prologue's `merge.policy` rather than a writer
version, treats a missing record as unknown, and refuses a record without
`adapter_trimmed: true`.
