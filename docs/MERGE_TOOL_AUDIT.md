# Merge tool: what to do next

> **What this document is, and how to read it here (added 2026-08-12, on the port).**
> A six-lens audit of the read-merge tool, with every claim adversarially verified and
> the disagreements resolved in favour of the verifier. It is a **record of a completed
> audit**, kept verbatim, and it is the answer to "why didn't you do X" — **§(c) is a
> list of ideas that look good and were measured and rejected, each with the number
> that closed it. Check it before proposing an improvement.** Everything in §(a), and
> R1–R7 of the roadmap, landed before the port; §(b) B1#4 (raw-blob IPC) and the §A4
> decode gate did not.
>
> It was written against **hulkrna's** tree, before the tool moved into zna, so its
> paths are hulkrna's. The ones that moved:
>
> | in this document | now |
> |---|---|
> | `lib/hulkrna/merge/*.py` | `src/zna/merge/*.py` |
> | `tests/test_read_merge.py` | `tests/test_merge.py` |
> | `tests/test_merge_zna_e2e.py` | `tests/test_merge_encode.py` |
> | `hulkrna-merge` | `zna merge` |
>
> Everything under `workflow/`, `config/`, `resources/` and `gather/` is pipeline
> configuration and **stays in hulkrna** — those references are still live, just in the
> other repo. Line numbers throughout are as of the audit and have since moved.

**Verdict tally across 51 candidates from 6 lenses:** 24 CONFIRMED, 18 OVERSTATED (kept at the verifier's reduced size), 7 REFUTED (dropped, with one labelling exception noted below), 2 ALREADY_HANDLED. Every number below is the verifier's where the two disagreed; I flag corrections inline.

**One-sentence answer:** the geometry defect you already know about is the only corpus-quality item that clears the bar, four lenses independently landed on the same 374 records, and everything else that looked bigger — the 12–16% microsatellite "ambiguity" — was measured to be a length error inside a repeat that causes no sequence damage and whose proposed fix costs 4.68% corpus duplication. The real leverage after the geometry fix is (i) throughput, which is 1.65x available for four small changes, and (ii) two *class* instruments that would have caught three of your four historical bugs.

---

## Where independent lenses converged (this is the signal)

| Finding | Lenses that found it independently | Agreement |
|---|---|---|
| `_build_merged` truncation | geometry, robustness, observability, testing | Four separate BAM extractions all produced **374 truncated merges / 137,796** |
| Microsatellite shift ambiguity | statistics, geometry, robustness (+ statistics' false-positive analysis) | All four verifiers independently concluded: **report it, do not gate on it** |
| §7 composition null | statistics, geometry, robustness | All three measured it and all three said don't build |
| `zna inspect --counts` cannot verify labels | observability (verified from zna source) | Single lens, but established from code, not demonstration |

Convergence on the truncation count from four different extraction pipelines is worth more than any single lens's argument.

---

## (a) Fix before the re-encode

A re-encode is expensive; these are the things that get baked into the ZNA and can't be patched afterward.

### A1. Rewrite `_build_merged` from the fragment span — `pairs.py:119-140`

**What.** Replace the per-direction construction with the signed-axis one:

```python
s = shift if direction == FORWARD else -shift
L = s + len2
take1 = min(len1, L)
seq = s1[:take1] + s2rc[take1 - s : L - s]
qual = q1[:take1] + q2rev[take1 - s : L - s]
```

`take1 > s` always holds because the scan never returns `s > len1-1`, so the tail slice can never go negative. Delete the strict xfail at `tests/test_read_merge.py:724-743`.

**Why / size (verifier-measured).** 374 of 137,796 merged records (0.2714%; 0.1897% of the 197,131 emitted records) are emitted shorter than the fragment length the tool itself inferred, and `--treat-unpaired-as-merged` then stamps `IS_FULL_FRAGMENT` on them — an interior EOS. 307 via read-through plateau, 67 via FORWARD shift==0. Median displacement is 1 base (247/374 lose exactly 1), mean 2.92, max 33.

**Structural fact none of the auditors stated up front, and the one that should drive your fixture design:** truncation requires `len1 < len2` **strictly** in every case. The read-through *flank* always has `L = n < len1` and can never truncate; only the plateau can, and the plateau contains negative shifts only when `len1 < len2`. Exposure is therefore bounded by the R1-shorter-than-R2 rate, which is 2.05% here — not the 3.78% unequal-length rate, and not the overall trim rate. That is why the production rate is 0.27% and not the several percent the redesign's simulation predicted: fastp adapter-trims PE symmetrically, so 96.2% of pairs have `len1 == len2` and mean |len1−len2| = 0.43.

**Sensitivity (corrected).** Under a *symmetric* independent 3'-trim model (which is what `--cut_tail --cut_window_size 1` at `config/config.yaml:73` actually does): **2.46% / 4.52% / 9.32% / 12.85%** of merges truncated at 5/10/25/50% trim rates. The auditor's 2.52/4.80/11.77/23.07% came from an R1-only model; the verifier reproduced both and the symmetric column is the realistic one. This still matters: the platform spans older 4-channel chemistry where the trim rate is much higher than this library's.

**Cost.** +0.024 µs per merged record (0.2478 → 0.2716 µs, min-of-9, interleaved, in-memory) ≈ +0.2% end to end. Corrected upward from the auditor's +0.013 µs.

**Risk.** Very low. Output is byte-identical on all 137,422 correctly-merged records (seq *and* qual) and a strict prefix-superset on all 374 broken ones. One behavioural consequence: merged records get slightly longer, so a handful of read-through fragments previously dropped at `min_read_length=40` now pass, and `mean_emitted_length` moves. Both correct.

**Land one integer compare beside it** (`span == len(seq)`, one JSON counter). Be honest about what it is: after the fix it is an algebraic tautology, not a metric — `len(s1[:take1]) + len(s2rc[take1-s:L-s]) = take1 + (L - take1) = L` identically. Its value is regression insurance on a 10-line function, and the observability auditor's claim that it is *free* was corrected: it costs ~+1% of the merge loop (~+0.4% end to end) across four interleaved runs. Take it anyway.

### A2. The 4-line e2e fixture change — `tests/test_merge_zna_e2e.py:81-82`

`test_every_edge_reported_is_a_true_fragment_boundary` already contains the right assertion and has never fired because `write_fastqs` builds both mates at exactly `READLEN`. Draw the two lengths independently:

```python
l1 = READLEN - (rng.randrange(1, 60) if rng.random() < 0.4 else 0)
l2 = READLEN - (rng.randrange(1, 60) if rng.random() < 0.4 else 0)
```

Verified both directions: **2 failed / 4 passed** on current code (`AssertionError: ('IU', 11, 100, 125)` — a 100 bp merged record for a 125 bp fragment, caught through merge → `zna encode` → `ZnaReader(with_ends=True)`), and **6 passed** against the fixed `_build_merged`. 0.46 s, 90 fragments. This is the cheapest instrument in the entire audit and the only test in either repo that crosses the merge/zna seam. Ship it with A1.

*(The testing lens also proposed a ~70-line property test for the same defect. Same detector, 17x the code. Take this one.)*

### A3. Decide the SAM-tag payload now — `workflow/rules/zna.smk:32`

`-T "*"` on `rigel/annotated.bam` carries **261.7 bytes of aux tags per record** against 129.4 bytes of SEQ. The tag set is NH HI AS NM MD RG MQ MC ms XS nM uT SA **plus thirteen rigel Z tags**: ZB ZC ZF ZG ZH ZI ZJ ZL ZN ZR ZS ZT ZW. Trimming to `-T "ZI,ZJ,ZF"` leaves 26.3 B: record 567.1 → 331.7 B (1.71x), gzipped 21.3 → 12.7 MB per mate, **production rule 3.370 → 2.612 µs/pair = 1.29x from one string.**

This was CONFIRMED and the auditor's magnitude was corrected *upward* (they claimed 1.09x): the performance lens benchmarked `bam/star.srt.rmdup.collate.bam` rather than `rigel/annotated.bam`, which `rigel_bam_to_fastq` actually consumes. **All of that lens's absolute numbers understate production by 2.4x for the same reason.**

Two things to settle before you trim, because they are re-encode decisions:
- The dropped Z tags carry transcript ID (ZT), gene ID (ZG), symbol (ZR), assignment class (ZC), splice status (ZS), weight (ZW). Grep the pipeline *and* khorana for all ten before dropping them. `zna inspect --counts` after the change is necessary but not sufficient (see A4).
- **Type mismatch:** rigel emits `ZF:i:3` while `resources/zna/label_defs.yaml` declares ZF as type C. Resolve this before you re-encode, or you will discover it in the corpus.

### A4. A G3 gate that decodes the produced ZNA — new rule after `zna.smk:96-128`

This is the highest-value **class fix** in the set and the answer to "what would have caught the bug that worried me."

Established from zna source, not from a demo: `_scan_flag_counts` (`src/zna/cli.py:1554-1592`) tallies only the per-record flags byte, and `_inspect_json` (`cli.py:1500-1550`) emits the label *schema* from the header but never a label *value*. **`zna inspect --counts --json` output is byte-identical between a correctly-labelled file and one where every label sits at its missing sentinel.** The G3 gate as written in `SOS_EOS_ENCODING_PLAN.md` §11 could not have caught the `key:`/`tag:` bug — and the note just added to `resources/zna/label_defs.yaml` ("Verify with `zna inspect --counts` after encoding") is wrong and will send the next person down the same dead end. Fix that line as part of this.

The gate: iterate `ZnaReader.records(with_ends=True)`, assert the five legal `(paired, r1, r2, has_start, has_end)` classes and nothing else, assert `count(R1)==count(R2)`, `count(R1 fwd)==count(R2 rc)`, `count(R1 rc)==count(R2 fwd)`, the RC coin is ~50/50 when unstranded, labels are not all at sentinel, and totals reconcile against `read_merge.json` within the N-drop bound. Emit `zna/geometry_check.json` so gather carries the verdict as a cohort field, not just an exit code.

**Unverified:** the verifier read the zna source and confirmed the central claim, but did **not** run the prototype — the 5-class whitelist, the four identities and the 0.13 s/library figure are the auditor's numbers only. Two cautions carried forward: gate the label check on labels actually being declared (a pre-rigel BAM legitimately has no ZI/ZJ/ZF), and fail loudly on an unrecognised encode flag combination rather than passing by default.

### A5. Fail fast instead of silently succeeding — `cli.py:310-345`, `fastqio.py:36-46`

Three cheap guards, all CONFIRMED, all in the "a library silently vanishes from the corpus" class:

- **Empty input.** `gzip.compress(b"")` as both inputs → rc=0, 20-byte valid gzip, JSON with `input_pairs=0`, then `zna encode` writes a 22-byte `.zna` with 0 records, rc=0. Nothing anywhere flags it. Three lines in `run()`: `if acc[0] == 0: raise SystemExit("no input pairs")`, behind `--allow-empty`. Note the two I/O backends disagree on a *zero-byte* `.gz`: pigz exits 1, stdlib gzip returns rc=0 — pin which backend any regression test targets.
- **Reader validation.** `read_fastq` checks only `if not q`. A FASTQ truncated inside its final quality line yields `len(seq) != len(qual)` and is emitted as a structurally invalid record with rc=0. Offset sweep: cuts of 2..100 bytes are accepted (99 offsets), 101..209 raise ValueError (109). Add `if len(s) != len(q): raise` and `if h[0] != 64: raise`. Cost re-measured at **+1.33% on the reader = ~0.13% end to end** (the auditor's +3.87% was corrected down). Exposure is narrow — `samtools fastq` exits nonzero when killed and Snakemake deletes failed outputs — so justify this on the fix being free, not on frequency.
- **Argument floors.** `--threshold-merge 0.5 --threshold-trim 0.5` is accepted: 93.93% merged vs 82.12%, no warning. `--min-read-length -5`, `--processes 0`, `--chunk-size 0` all accepted. ~6 lines. Derive the `t_merge` floor from `threshold_bits(read_len, alpha)` — that function exists at `overlap.py:86` precisely so the bound is executable. The realistic path is a config typo propagating into a full-cohort re-encode with nothing but `merged_pct` to signal it.

### A6. Optional, if you have appetite: max-posterior overlap consensus

The only item here that improves *base-level* corpus quality, and it can only be captured by re-encoding.

Derived analytically from the real joint (Q1,Q2) mismatch distribution over 33k real merges (4,915,715 merged bases, 28,001 overlap mismatches), under "exactly one mate is wrong", P(R1 wrong) = p1/(p1+p2): residual base errors per 100k merged bases are **R1-wins 196.2, the deployed fastp Q14/Q30 correction 77.5, max-posterior 51.1 — a 34% reduction over what is running today.** The auditor claimed 12.2%; the verifier's number is larger and better grounded. The driver is the Q25 band fastp's hard cutoffs never touch: Q11/Q25 (2.7% of mismatches) and Q25/Q37 (2.4%), where R1-wins is ~95% wrong and the posterior is 15–25:1 the other way.

**Do it for the consensus, not the merge rate.** Two of the statistics auditor's load-bearing claims are false: the quality-aware argmax differs from the constant-weight argmax on 3.78% of real pairs (not 0), and a winner-rescore captures ~40% of a full quality-aware scan's effect, not 100%. The merge-rate gain is +0.095pp — negligible. Note also Q11's true per-base error back-fits to 0.141, **1.8x its nominal Phred value**, so a nominal table under-penalises exactly the bin that matters; calibrate from the data.

**What it deletes:** `_correct_r1_overlap` (pairs.py:45-75), `--correction`, `_GOOD_QUAL`, `_BAD_QUAL`, and the `config.yaml` key — the last three knobs §9 says haven't been removed. It also closes a live hazard: the correction gate at `pairs.py:66` is **effectively untested and currently on** (`config/config.yaml:145` sets `correction: true`; 27,003 bases corrected on 167,784 pairs). A Phred sweep against the existing suite shows `_BAD_QUAL` stays green for Q10..Q39 and `_GOOD_QUAL` for Q0..Q40 — the low bar can rise 25 points and the high bar fall 30 with the suite fully green, and a widened gate changes emitted sequence on 3.06% of pairs. This is the one candidate whose real-world exposure the verifier revised *upward*. If you don't do A6, write the 6-case boundary test instead (~25 lines).

**Cost, and an important cross-lens correction.** The statistics verifier priced a fused njit scan+rescore at +0.366 µs/pair = +7.9% end to end against a 4.65 µs/pair single-process steady state — 3x the auditor's claim. **But in the production `--processes 8` config that cost is ~zero:** the performance lens measured p4 1.923 ≈ p8 1.938 ≈ p12 1.965 µs/pair, i.e. the workers are idle and the entire merge computation contributes 0.005 µs/pair to wall time once IPC is fixed. Neither lens could make this join. The throughput objection to A6 is much weaker than it reads in isolation.

**Risk.** Merged quality strings must be recomputed on the consensus or downstream Q values become dishonest. Non-binned chemistry (MiSeq/HiSeq) touches ~1600 table cells instead of 16 — still 13 KB, but the timing was measured on 4-level NovaSeq data.

---

## (b) Worth doing next

### B1. The throughput stack — 1.65x measured, in this order

All four measured by the verifier on the **correct** input (`rigel/annotated.bam`, 1,006,704 real pairs, interleaved, min-of-5/7). Cumulative: 2.272 → 1.377 µs/pair at 4M pairs.

| # | Change | Verified win | Cost | Note |
|---|---|---|---|---|
| 1 | `-T "ZI,ZJ,ZF"` (A3) | **1.29x** | one string | Auditor said 1.09x — corrected *up* |
| 2 | Reader pigz `-p 1` | **1.13x** (1.20x at t4/c2000) | ~2 lines | See labelling note below |
| 3 | `--chunk-size` 50000 → 2000 | **1.33x** at 1M, **1.24x** at 4M, ~1.21x asymptotic | one line | Auditor's 1.31x headline is a 1M-pair artifact |
| 4 | Raw-blob IPC | **1.20x** at 1M, **1.25x** at 4M (matched chunk) | ~60 lines | Auditor's 1.65x headline was the whole bundle credited to one change |

**Labelling note on #2:** the pigz-threads candidate carries verdict REFUTED, but what was refuted is the *auditor's negative conclusion* ("no wall-time win on this machine"). The verifier patched only `fastqio._open_binary_read`'s threads argument and measured 3.353 → 2.962 µs/pair on the production rule config. The mechanism is oversubscription, not inflate throughput: the rule spawns 8 workers + 2 reader pigz × 8 threads + 1 writer pigz × 8 threads ≈ 32 threads against an 8-core allocation. This is a real 1.13x for two lines and should not be dropped along with the genuinely-refuted items. Corrected on the writer side too: `--threads 1` costs 1.21x (not 1.8x), and t4 vs t8 are within 1%.

**Corrections you should carry.** The chunk-size mechanism on the ticket ("allocator pressure") is wrong — running the unmodified `cli.run()` under `gc.disable()` moves c50000 by 1.5% and c2000 by 0.3%. The "hard 2.05x `--processes` ceiling caused by pickling" is a chunk-size artifact: at c2000 the ceiling is 2.73x. Both matter because the right recommendation — re-run the sweep on the actual cluster node type — needs a correct model to interpret.

**Four defects in the IPC prototype**, three of them not in the stated risks, one silent and disqualifying until fixed:
- truncated R1 mid-record → `IndexError` instead of `"truncated FASTQ record in <path>"`
- R2 shorter than R1 → `IndexError` instead of `"R2 exhausted before R1"`
- **R1 shorter than R2 → returns OK and silently drops the trailing R2 records.** No error, no stat. In an unattended platform tool this is the one failure mode you cannot ship. The claim that "the desync check moves entirely into the worker (it already is there today)" is wrong: `base_name(h1) != base_name(h2)` cannot see records that were never read.
- CRLF input yields `seq=b'ACGT...\r'`, which the kernel would score as a base and write into the corpus.

All fixable inside the 60-line budget. Sequence this **last** of the four.

### B2. `ProcessPoolExecutor` swap — `cli.py:246-263`

`mp.Pool` does not detect abrupt worker death. Verified independently: SIGKILL one worker → parent **alive at 150 s at 0.0% CPU** with four respawned children, no output, no JSON, no error, unreadable partial `.fq.gz`. `_maintain_pool` respawns, but `IMapUnorderedIterator.next()` blocks on `_index == _length` forever. PPE with a `2*P` bounded submit window and `cf.wait(FIRST_COMPLETED)` raises `BrokenProcessPool` in **0.97 s**, at **+2.7%** throughput (3.694 → 3.794 s on 2M pairs, min-of-7 interleaved — inside noise). ~25 lines in `_run_parallel`, no change to `_process_chunk` or `_init_worker`.

Two corrections to the auditor's framing: the tool cannot plausibly OOM itself (590 MB peak parent RSS on a 4M-pair library at p8, against `24000 + 3000*threads` = 48 GB requested; memory is bounded by `--chunk-size`, not library size), so the trigger is node-level/sibling OOM or a native crash in JIT/OpenMP code. And `cli.py:262` does *not* relabel the reader's ValueError — `SystemExit(str(e))` preserves the message, only the traceback is lost. It is still the only unbounded-time failure in the tool and it is silent.

### B3. Make the JSON reach the cohort, and make it readable

- **`gather/tools/read_merge.py:39-46`** keeps only top-level `int`/`float`, so `score_histogram_bits` — the headline observability addition of this rewrite — **reaches the cohort store nowhere.** 21 keys survive today; `length_histogram` survives only because it has its own TABLE OutputSpec. Widen the filter to `(int, float, str)` (matching the `params` branch two lines below), declare TABLE outputs for any new histogram, update `docs/cohort_schema.md:47,142-143`. Note the str branch is currently dead code — there are zero string values in the JSON — so it only becomes load-bearing once provenance lands, which is why the two are coupled.
- **Provenance.** No `__version__` anywhere under `lib/hulkrna/`, no `setup.py`, no `pyproject.toml`. Add `tool`, `tool_version`/git SHA, `numba`, `elapsed_s`, `pairs_per_second`. Smaller than sold: `gather/tools/pipeline.py` already flattens `pipeline_config.yaml` into `pipeline/params`, so `read_merge.threshold_merge`, `err_rate` etc. are already cohort-queryable — **only code identity is missing**. `pairs_per_second` as a cohort field gives free detection of nodes where numba failed to load or pigz was missing.
- **Rebin the score histogram or replace it.** `score = olen*1.9855 − diff*8.2143`, so a perfect n-base overlap floors to bin `2n−1` for all n<69: bins 8–60 hold **8,052 odd counts vs 637 even = 12.6:1**, and bins 1–7 are permanently dead (`find_overlap` returns 0.0 below `t_trim`). Rebinning at 2 bits collapses the ratio to 1.08 (one character). Better: emit `overlap_length_histogram` instead — a uniform-random control (60k pairs, real length profile) predicts ~330 chance detections at olen==5 against **326 observed**, i.e. essentially 100% of olen==5 detections are chance, and the cliff position reads the merge threshold directly. Do it now: nothing consumes `score_histogram_bits`, so there is nothing to break.
- **Insert-size distribution, split by outcome**, with **both** censoring bounds published. p5/p25/p50/p75/p95 = 61/101/141/195/262, TLEN-validated exact for **99.64%** of merges. Hard-censored at exactly 285 = 150+150−ceil(28/1.9855) (269 counts at 285, zero at 286+) — *and* hard-floored at 40 by `min_read_length`, with 270 counts sitting exactly at 40 and a flat ~270/bp density through 40–50, i.e. the real distribution continues below the cut and 641 merged singles are dropped. The auditor warned about the upper bound and then reproduced their own error at the lower one. Also: `insert_size_max_detectable` is not one number — this library has 111 distinct `len1` values and the bound is per-pair `len1+len2−15`.
- **`overlap_mismatch_rate`** (accumulate the `diff`/`olen` `find_overlap` already returns): 0.00911 all-detected, 0.00866 merged-only, against `err_rate = 0.01`. Worth having — but **not** for the reason given. Under `n*1.9855 − d*8.2143 ≥ T` the observable rate is hard-capped at ~0.22, so "> 0.10 means chance alignments are being accepted" is wrong; it can never approach 0.75. Injecting per-base errors into R2: 0.0091 → 0.0234 (2%) → 0.0805 (10%) → 0.2006 (35%), saturating at the predicted ceiling while `merged_pct` collapses 82.2 → 8.8. Its actual value is that at 5% degradation it moves 5x while `merged_pct` moves 1% — it is the **sensitive early-warning channel**, `merged_pct` is the catastrophic one.

### B4. Two small test items with real teeth

- **A 12-line block-loop test.** Sweep a single mismatch across every position of a 40 bp overlap and assert `score == score_of(39, 1)`. This kills `block8-stride-7` (`k += 7` in the hot unrolled loop), which passes the entire existing suite and changes 6.34% of scores, 2.79% of emitted records and 0.26% of merge/trim/keep decisions on 167,784 real pairs. The testing lens proposed a 90-line brute-force differential for this; the verifier showed 5 of the 6 mutants it claimed to kill are **already killed** by the existing suite, so its unique contribution is this one mutant. If you do want a kernel differential, write it **order-independently** — assert the returned shift is *an* argmax with the true n/d/score — because ties are not rare: 20.1% of accepted scans at the auditor's own exhaustive-arm threshold and 0.845% on real pairs at the production floor. Their order-replicating reference would break on exactly the SIMD/vectorisation rewrite it exists to guard.
- **`bases_trimmed`, 2 lines.** `cli.py:123-124` and `:208-209` use `next((r for r in records if r[0] != h1), None)`; when headers match, `r2` is None and the whole of R2 is charged. Replace with `records[1]` — guaranteed to exist in that branch (`pairs.py:206`), and simpler than what it replaces. It is **correct today** (8,678 = ground truth, both paths) because `zna.smk:32` passes `-N`; the 12.1x error is latent, for any future non-BAM input path. Free, so take it.

### B5. Two-line hygiene

`os.path.realpath(in1) == realpath(in2)` at startup — `--in1 X --in2 X` currently gives rc=0 with every read emitted twice as a mate pair (merged collapses 82.12% → 0.05%, nothing reads it). Zero false positives. The `--min-merged-pct` floor the auditor also proposed is the only thing that catches swapped mates (statistics are byte-identical, and `--strand-normalize` *is* wired from `salmon_infer_libtype` at `zna.smk:124-125`, so a swap inverts the whole library's strand) — but a genuinely long-insert library merges at a low rate, so warn-only or per-cohort.

---

## (c) Deliberately not worth doing

**Do not gate merges on shift ambiguity — in any form.** Four lenses found it; all four verifiers said no. The statistic is real and reproduces (12.09–16.4% of merges have a rival above threshold; margin < 10 bits on 12.55%), but:

- **It is a microsatellite detector, not an identifiability test.** Every shift's score is an LLR against the same null, so the posterior over shifts is exactly `2^score` normalised — computable, and nobody computed it. Posterior mass >50% off the argmax: **1.60% of merges**. Off argmax±4: **0.15%**. Off ±9: **0.04%**. 83% of merges have <1e-12 mass off the argmax. A rival at 28 bits against a winner at 250 bits is 2^222:1, not ambiguity.
- **Against STAR-derived ground truth the argmax is exactly right for 90.4% of merges** (73,604 pairs, soft-clip-corrected). Of the 6,161 disagreements, 73% emit a read ≥80% period≤6 periodic, 56% are 100% periodic, 48% differ by a pure indel. On the 1,699 non-repeat disagreements the *data prefers the kernel*: in 58.3% the STAR shift scores ≥40 bits worse, in 14.8% it isn't a valid alignment at all. Merges with a genuinely competitive alternative: **458/66,745 = 0.69%**, and those inspect as impure CA and higher-order satellites.
- **A wrong shift produces a read of the wrong *length*, not the wrong *sequence*.** Error-free GT simulation, CA120 cell with 25.35% wrong shifts: **0.0157%** of 24-mers in wrongly-shifted merged reads are absent from the true fragment, and the emitted read **never** overruns the true fragment end (0.00 bases/read).
- **The fix is net-negative.** 89.9% of flagged merges are REVERSE geometry, where `pairs.py:200` requires FORWARD — so **99.9% of declined pairs fall through to keep-both-whole**, not to trim (13 of 16,659 would actually trim). That re-emits the overlap twice: **+2,034,354 duplicated bases = +4.68% of the library's input bases**, which is precisely the redundancy bias objective #1 exists to remove, plus a 10pp merge-rate loss. On its own ground truth the gate destroys 1.3 correct full-fragment records per wrong one caught (1.0:1 on non-repeat reads) and still misses 87.
- **Cost:** +0.356 µs/pair, +18.3% of the kernel, ~+7.7% end to end.

Two useful residues. (i) The redesign's preamble item 2 — "every runner-up is pinned just below the winner, so reporting the margin honestly would mean disabling the pruning" — **is wrong**: set the pruning floor to `max(t_trim, best − margin_req)` and the margin decision is exact (0 disagreements vs a fully unpruned scan at margins 4/6/8/10/16). Correct the doc. But 99.43% of ambiguous merges are period≤6 repeats covering >80% of R1, only 98 records in 137,796 are ambiguous without being a microsatellite, and an unpruned offline scan over all 167,784 pairs runs in 2.4 s — so sample it offline, don't pay ~+1.7% forever. (ii) The genuinely actionable observation is a **corpus** one, not a merge-tool one: **13.7% of merged reads in this cfRNA library are ≥80% tandem repeat**, and their emitted length is arbitrary to within a repeat period. That is a khorana filtering question. Note it is one library: cross-library replication gives 12.09% / 0.177% / 0.058% ambiguous merges across LBX0190 / LBX0588 / VCaP.

Also on the "don't build" list, each with the measurement that closes it:

| Item | Why not |
|---|---|
| **§7 composition-aware null** | Three lenses, all negative. For CA-vs-CA `p_null = 0.500` exactly, a match is still 0.986 bits, a 74-base perfect overlap still scores 72.9 against T=28 — it does not fire on the case it was designed for. Costs **−0.55pp merge rate overall and ~2pp on exactly the low-complexity reads it was meant to protect**, moves 2.26% of argmaxes, leaves wrong-shift rates unchanged. Fused njit cost +3.4% e2e (the auditor's +16% was a strawman implementation). Honest caveat: on the TLEN-truthed subset it is +2 net correct out of 44,670, and that subset under-samples low-complexity reads 15.2x — so for the affected population the verdict is *unmeasured*, not *measured null*. If ever revisited: the redesign's formula is orientation-ambiguous (p_null p90 = 0.462 on R1 vs revcomp(R2), 0.269 on raw R2 — opposite signals). |
| **Per-pair threshold from read length** | Nothing in the code proposes it (`threshold_bits` at `overlap.py:86-94` is documented as derivation-only and unused). Uniform-null p99.9 is 9.9 bits flat across L=50..300, 0/20,000 above 28 at every length. Real content is a **doc correction**: on shuffled real pairs the spurious rate is 1.710% at T=28 and still 1.323% at T=100 — raising T by 72 bits removes a quarter of them, because 89.3% involve a >90% two-base read. The T→α mapping in `READ_MERGE_REDESIGN.md:106-116`, `cli.py:41-44` and `pairs.py:26-31` is off by 26 orders of magnitude against real repeat content and must say so. |
| **Estimating `err_rate`** | `e_hat = 0.00866` from 13,577,477 aligned overlap positions vs 0.01 assumed — already right to 15%, adopting it costs −0.04pp merge rate. But do not repeat the auditor's "insensitive everywhere": on **real** pairs a 0.002–0.05 sweep moves the merge rate 0.98pp (~8x their simulated 0.13pp) and changes up to 3% of argmaxes. The right statement is narrower: 0.01 sits inside the flat region. |
| **N special-casing, indels/Smith-Waterman** | N: 284/167,784 pairs (0.169%), 298 N bases, **every run length exactly 1**, 139 pairs (0.083%) with N in the winning overlap. Free as a side effect if A6 lands (Q2 → e clamps to 0.75 → both weights exactly 0.00 bits); never worth a branch in the hot loop. Chemistry-dependent, so prefer the structural route. Indels: 3.93% of merges carry anomalous mismatch load, 55.7% of that population is >90% two-base — consistent with wrong shifts in repeats. **The one-gap model fit ("median reduction of 1 mismatch, 4 pairs materially explained") is single-sourced and was not re-run — treat as unverified.** |
| **SWAR 2-bit packing** (`overlap.py:123-139`) | The 4.8x is a microbenchmark of a no-bail sweep that does ~5x the real pruned scan's work. Corrected: ~1.9x on `_scan` compute, 1.80x on the call, **1.19x on the whole `--processes 1` tool, 0 in production**. The prototype's own equivalence assertion prints `agree=False` in all three regimes. Plus an N-semantics change (N collapses to A) in a corpus-quality tool. |
| **Inlining `_shift_score`, hoisting `score_weights`** | 1.107x (0.181 µs/pair) and ~0.08–0.13 µs/pair respectively. 3.4% and 2.5% at p1, **zero at p8** (p4 1.923 ≈ p8 1.938 ≈ p12 1.965 — workers idle). Useful control the verifier added: a null njit kernel with the same signature costs **0.221 µs/pair of pure numba dispatch**, which caps what *any* kernel change can return through this interface. Both correctly self-assessed as "don't take it." |
| **Capped mismatch budget; record-aligned block reading** | Correctly-recorded negatives. The cap cannot be made exact (re-check fires on 18–55% of pairs; exact two-phase is 0.717x) and as an approximation D=16 loses 0.68% of merges. Block reading is 0.76x — `read1()` or ≤64 KB alternation is mandatory for any future multi-stream reader, worth a code comment. |
| **Per-record `ZE:i:` geometry tag** | ALREADY_HANDLED. Every load-bearing fact is already the documented justification at `zna.smk:81-89`. Bytes are cheap (+0.44% gzipped) but it is a two-repo coordinated change for a hole A1 closes inside one tool. |
| **Gating `--correction` on `t_merge`** | REFUTED. 1,117 of 27,003 corrections (4.14%) are sub-threshold, but 96.6% of those bases sit on the shift STAR independently agrees with — ~38 bases per library across 25M R1 bases are on a genuinely wrong alignment. And 965 of the 1,117 are in the FORWARD trim band, where R2's overlap copy is discarded from the output, so correcting R1 from R2 there is the *only* way R2's evidence survives. Keep instead: split `bases_corrected` into merged/unmerged in the JSON (4 lines). |
| **`--count-n` for exact ZNA reconciliation; `--self-check` byte invariants; widening the parallel-vs-single test** | All REFUTED. The N counter doesn't deliver exactness — `--npolicy drop` is pair-atomic, so a per-record counter reports 284 against zna's 316. The self-check's invariants are tautologies against the current construction (`_build_merged` returns literal prefixes of s1; the min-length assert restates `pairs.py:218-221`; the trim-tiling identity holds by construction) — the one non-tautological piece is seq/qual length parity, which belongs in the reader (A5). The parallel test found nothing on 167,784 real pairs at processes 1/3/4/8 × chunk 1/137/7777/50000/1000000 × correction on/off. |
| **The mutation harness as a standing instrument** | Of 6 survivors: 1 provably equivalent, 1 *is* the proposed fix, 1 (merged-qual-from-r1, 27.2% of records) has exactly zero impact because zna discards quality (`zna/cli.py:292,316`), 1 touches 0.047% of pairs. Two real survivors: `block8-stride-7` and `correction-quality-gate` — both covered above by targeted tests an order of magnitude smaller. |

---

## The class fixes, named explicitly

You asked because one bug worried you. These are the ones that close families rather than instances:

1. **Build geometry from the fragment span, once** (A1). The class is *per-direction construction carrying an implicit `len1 ≥ len2` assumption*. The one-line repair is half of it; the other half is **A2/fixtures that draw `len1` and `len2` independently** — `make_pair` and `cycle_pair` structurally cannot express the failure, which is why 135 tests were green over a live defect. Target `len1 < len2` specifically; that is the entire exposure surface.
2. **Decode what you actually wrote** (A4). The class is *nothing in the chain reads the produced artifact*. It covers the `key:`/`tag:` bug, the missing-`--interleaved` bug, orphaned mates, broken co-orientation, biased RC coin, and count divergence — three of your four historical defects in one gate. And it retires a verification instruction that provably cannot work.
3. **Fail instead of silently succeeding** (A5). The class is *a library vanishes or degrades with every stage exiting 0*: empty input, structurally invalid FASTQ, a threshold typo, a hung pool. Each individually is low-probability; collectively they are the only way a platform-scale corpus loses data invisibly.
4. **Provenance** (B3). The class is *"which build produced this file?"* — the question every future corpus defect will open with. Cheap now, unanswerable later.
5. **Make the tool's numbers reach the cohort** (B3). A statistic nobody can query is not observability. The rewrite's headline addition currently reaches nowhere.

A smaller standing hazard worth noting: `_account` (`cli.py:115-137`) and `_process_chunk` (`cli.py:187-222`) duplicate ~35 lines of bookkeeping, and `_run_parallel` hand-merges every histogram. They are line-for-line equivalent today (verified on 167,784 real pairs across every process/chunk combination), but every statistic you add doubles that surface. Collapsing them — `_run_single` calling `_process_chunk` on one chunk — is ~30 lines and removes the class rather than testing around it; it needs its own throughput measurement because the single path writes incrementally rather than accumulating a blob.

---

## Doc and comment corrections to land alongside

- `workflow/rules/zna.smk:83` — "That is exact, not a guess" is false for 1 in 370 merged reads until A1 lands.
- `resources/zna/label_defs.yaml` — the new "Verify with `zna inspect --counts`" note cannot detect the bug it was added for.
- `docs/READ_MERGE_REDESIGN.md` preamble item 2 — runner-ups are *not* pinned; the margin is exactly computable under a relaxed pruning floor.
- `docs/READ_MERGE_REDESIGN.md` §7 — measured three times, negative; record it rather than leaving it as a roadmap item.
- `docs/READ_MERGE_REDESIGN.md:106-116`, `cli.py:41-44`, `pairs.py:26-31` — the T→α mapping holds against uniform sequence with ~14 bits of margin and fails by 26 orders of magnitude against real repeat content.
- `overlap.py:53` — "anything else → N" is false; `bytes.maketrans` passes unlisted bytes through, so `rc(b'RYKMSWBDHVN') = b'NVHDBWSMKYR'`. Fix the comment, **not** the behaviour: mapping non-ACGTN to N would change the kernel's N-vs-N semantics (currently +1.9855 bits/position) and needs measuring. Exposure is zero — of 167,784 real pairs, every non-ACGT byte is N.
- `overlap.py` module docstring — the scan is O(L²) on non-overlapping pairs (L=150: 3.88 µs, L=1000: 116 µs, L=20000: 40 ms). Not reachable from a STAR BAM, and 2×150 → 2×300 is ~1.8x for 2x the bases at the whole-tool level (the "cliff" framing was corrected), but a ~5-line `--max-read-length` fail-fast makes a pathological input distinguishable from the pool hang.

---

## What the auditors got wrong, and what is still unchecked

**Errors worth knowing about, because they change how you read the raw reports:**
- The performance lens benchmarked `bam/star.srt.rmdup.collate.bam` instead of `rigel/annotated.bam`. Every absolute µs/pair figure in that lens understates production by **2.4x**, and its enumerated tag set missed all thirteen rigel Z tags.
- The robustness lens counted the truncation defect at 67 records (0.049%) — it counted only the FORWARD path while describing both. True count 374; the read-through path is 82% of it. Its sequencing advice ("the span rewrite is a small win, the rival gate is the large one") inverts once corrected.
- The observability lens built an entire candidate on misreading `tests/test_read_merge.py:931`, which does assert `score_histogram_bits` equality across process counts.
- The statistics lens's simulator absolute magnitudes are ~4x inflated (signed shift error −0.03/−0.26/−4.48 vs their −1.10/−3.42/−17.22); its "quality weights never change the argmax" is false (3.78% of real pairs); its "full scan == winner rescore" is false (rescore captures ~40%).
- The geometry lens's TLEN ground truth was contaminated (|TLEN| up to 99 Mb among pure-M pairs) until filtered to ≤500.
- The testing lens's simulated truncation rates (1.44%/5.48%) disagree ~2x with the repo's own recorded figures (0.83%/2.68%) for the nominally identical experiment, and neither was reproducible. **Treat every simulated truncation column as unsubstantiated; the real-data column (374 records) is the one that holds.**

**Genuinely unchecked:**
- The G3 gate's runtime behaviour — 5-class whitelist, four identities, 0.13 s/library. The *reason* to build it (that `zna inspect` cannot see labels) was established from zna source; the gate itself was not run.
- The one-gap indel model.
- "Two of the four historical defects were which-build questions" — asserted, unverifiable.

**Unresolved and not resolvable from inside the tool:** against aligner ground truth, ~0.18–0.21% of merged records have a fragment length wrong by >20 nt; ~36% of those are pairs the aligner says do not overlap at all (inverted-repeat/hairpin self-complementarity), the rest are argmax losses inside repeats. A margin guard recovers 11–20% of them. The only real lever would be carrying the aligner's own fragment length through as a SAM tag, which needs an upstream stage to stamp it (`samtools fastq -T` copies aux tags but not TLEN). At 0.18%, not worth a pipeline stage.
