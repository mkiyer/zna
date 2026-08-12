# Read-merge tool: status board

The standing status board for `src/zna/merge/` (`zna merge`). Written 2026-08-12 in
hulkrna, where the tool was developed, as the roadmap to a portable state; it became
zna's board when the port landed the same day. Everything above "Open" is history —
kept because each entry records the measurement that justified the change.

Design: [READ_MERGE_REDESIGN.md](READ_MERGE_REDESIGN.md) (implemented).
Audit: [MERGE_TOOL_AUDIT.md](MERGE_TOOL_AUDIT.md) — §(c) is the "measured and rejected"
list; check it before proposing an improvement.
The port: [READ_MERGE_PORT_TO_ZNA.md](READ_MERGE_PORT_TO_ZNA.md).
Cross-repo contract: [SOS_EOS_ENCODING_PLAN.md](SOS_EOS_ENCODING_PLAN.md).

**Objectives, in priority order.** (1) Corpus quality and correct fragment-end
supervision — a special token at an interior position is noise injected into the signal
the downstream model is learning. (2) Throughput at platform scale. Everything below is
ranked against those two, and (1) beats (2) whenever they conflict.

---

## Where we are

| Milestone | State |
|---|---|
| M0 — LR redesign implemented (one score, two thresholds, argmax) | **done** |
| M1 — geometry correct on every read-length combination | **done** (2026-08-12) |
| M2 — fails loudly instead of silently; no unbounded-time failure | **done** (2026-08-12) |
| M3 — statistics reach the cohort; provenance recorded | **done** (2026-08-12) |
| M4 — intermediate size and pool scaling | **done**; raw-blob IPC deferred, see R7 |
| M5 — freeze, then MOVE to the zna package | **done** (2026-08-12) |
| M6 — numba → C++ | **done** (0.4.0); [MERGE_CPP_DESIGN.md](MERGE_CPP_DESIGN.md) |
| M7 — `zna encode --merge-pairs` | **next** |

**M6 and M7 swapped order (2026-08-12).** The C++ design was measured first and the
measurement changed the sequencing: at `--processes 8` the merge computation already
contributes almost nothing to wall time, so the win from native code is *collapsing
eight cores into one*, and the remaining time is gzip I/O — which is exactly what
`--merge-pairs` deletes. Doing C++ first also means `--merge-pairs` is written once,
against the final backend, instead of twice.

Historical paths: everything below `lib/hulkrna/merge/` is now `src/zna/merge/`, and
`tests/test_read_merge.py` / `tests/test_merge_zna_e2e.py` are `tests/test_merge.py` /
`tests/test_merge_encode.py`. References to `workflow/`, `config/`, `resources/` and
`gather/` are pipeline configuration and are still hulkrna's.

---

## Done — landed 2026-08-12

### M1.1 `_build_merged` builds from the fragment span — the geometry defect

`pairs.py` case-analysed the direction and only ever emitted R2's tail on the forward
path with `shift > 0`, which carried an implicit `len1 >= len2` assumption. Whenever R1
was quality-trimmed shorter than the fragment — forward with `shift == 0` and
`len2 > len1`, or **any** read-through with `len1 < L` — the merged read was truncated
to R1 while `--treat-unpaired-as-merged` still stamped `IS_FULL_FRAGMENT` on it. An
interior EOS token, which is exactly what the ZNA contract exists to prevent.

- **Measured, real data:** 374 of 137,796 merged records (**0.271%**; 0.19% of emitted
  records) on 167,784 pairs from `mctp_LBX0190`. 307 read-through, 67 forward-shift-0.
  66% lose exactly one base; 1,091 bases lost library-wide; max 33.
- **Structural bound:** truncation requires `len1 < len2` *strictly*. Only 2.05% of
  pairs in that library qualify (fastp adapter-trims PE symmetrically), and 10.9% of
  those truncate. Under a symmetric 3'-trim model the rate is 2.5 / 4.5 / 9.3 / 12.9%
  at 5 / 10 / 25 / 50% trim rates — **so this scales with chemistry**, and the platform
  spans older 4-channel data where trimming is far heavier.
- **Fix:** build from the one quantity the scan actually infers. `s = ±shift`,
  `L = s + len2`, `take1 = min(len1, L)`, `seq = s1[:take1] + s2rc[take1-s : L-s]`.
  `len(seq) == L` identically, for every geometry.
- **Regression evidence:** byte-identical `seq` *and* `qual` on all 137,422 correctly
  merged records; strict prefix-superset on all 374 broken ones; zero non-supersets.
  Independent ground truth: inferred `L` matches BAM `|TLEN|` on 99.635% of pure-M
  proper pairs. Cost +0.024 µs/merged record (~+0.2% end to end).

### M1.2 Fixtures that can express the failure

`make_pair` / `cycle_pair` and the e2e `write_fastqs` all built both mates at one
length, so 135 tests plus the whole seam suite were structurally blind. Now:
`TestUnequalReadLengths` (both mechanisms + a span invariant swept over five length
combinations) and independently-drawn mate lengths in `test_merge_zna_e2e.py`.
**Verified to fail on the pre-fix builder** (8 failures) and pass after.

### M1.3 A test for the unrolled block loop

`_shift_score` accumulates mismatches branchlessly in blocks of 8 and tests the bail
once per block. A stride bug there (`k += 7`) mis-scores every overlap containing an
error — 6.34% of scores and 0.26% of merge/trim/keep decisions on real data — and
**passed all 135 tests**. `test_block_loop_sees_every_position` sweeps one mismatch
across all 40 positions of a full overlap. Verified: it is the *only* test that kills
that mutant. This matters most for the C++ rewrite, which will touch exactly this loop.

### M1.4 Documentation that was wrong

- `zna.smk` — the `--treat-unpaired-as-merged` justification now names what it rests on
  and which tests pin it, instead of asserting "exact, not a guess".
- `READ_MERGE_REDESIGN.md` §8 note 2 — the claim that the tie margin is uncomputable
  under pruning was **wrong**; `max(t_trim, best − margin_req)` gives an exact margin.
- `READ_MERGE_REDESIGN.md` §3 / `cli.py` epilog — `T = log₂(N/α)` bounds detection
  against *chance alignment of uniform sequence* (holds, 0/40,000 at every read length),
  not against real repeat-rich sequence (1.7% at T=28, still 1.3% at T=100).
- `READ_MERGE_REDESIGN.md` §7 — the composition-aware `p_null` guard is now recorded as
  measured-negative three times, not left as a roadmap item.
- `overlap.py` — "anything else → N" was false; IUPAC codes pass through uncomplemented.
  Comment fixed; behaviour deliberately not changed.
- `resources/zna/label_defs.yaml` — the verification note pointed at `zna inspect
  --counts`, which provably cannot detect the bug it was added for.

---

## Done — M2/M3/M4, landed 2026-08-12

All of the below are landed and pinned by tests. What follows each item is what was
actually verified, not what was proposed.

### R1. `multiprocessing.Pool` hangs forever on abrupt worker death — **FIXED**

`mp.Pool` does not detect a SIGKILLed worker. Verified: kill one worker →
parent **alive at 150 s at 0.0% CPU**, children respawned, no output, no JSON, no error,
unreadable partial `.fq.gz`. `_maintain_pool` respawns the process but
`IMapUnorderedIterator.next()` blocks on `_index == _length` forever.

At platform scale this does not look like a bug, it looks like a slow node, and it burns
the job's entire wall-clock allocation. It is the **only unbounded-time failure** in the
tool. Trigger is node-level or sibling OOM, or a native crash in JIT/OpenMP code — not
self-inflicted (peak parent RSS 590 MB on a 4M-pair library against 48 GB requested).

**Fix:** `ProcessPoolExecutor(mp_context=fork)` with a `2*P` bounded submit window and
`cf.wait(FIRST_COMPLETED)`. Raises `BrokenProcessPool` in **0.97 s**, costs **+2.7%**
(inside noise). ~25 lines in `_run_parallel`; `_process_chunk` and `_init_worker`
unchanged.

### R2. Silent successes — **FIXED**

- **Empty input.** Two empty `.gz` inputs → rc=0, `input_pairs=0`, then `zna encode`
  writes a 22-byte 0-record `.zna`, rc=0. Nothing anywhere flags it. 3 lines behind
  `--allow-empty`. (Note the backends disagree on a *zero-byte* `.gz`: pigz exits 1,
  stdlib gzip returns 0 — pin which one any test targets.)
- **Reader validation.** `read_fastq` checks only that line 4 is non-empty. A FASTQ
  truncated inside its final quality line yields `len(seq) != len(qual)` and is emitted
  as a structurally invalid record, rc=0. `if len(s) != len(q): raise` plus a `@` check
  on the header: **+1.33% of the reader ≈ +0.13% end to end**. Narrow exposure (gzip
  truncation already raises, Snakemake deletes failed outputs) — justify it on being
  free, not on frequency.
- **Argument floors.** `--threshold-merge 0.5 --threshold-trim 0.5` is accepted and
  merges 93.9% instead of 82.1%, silently. `--min-read-length -5`, `--processes 0`,
  `--chunk-size 0` likewise. ~6 lines. Derive the `t_merge` floor from
  `threshold_bits(read_len, alpha)` — that function exists precisely so the bound is
  executable rather than folklore. The realistic path here is a config typo propagating
  into a full-cohort re-encode with nothing but `merged_pct` to signal it.

### R3. `bases_trimmed` wrong on identical mate headers — **FIXED**

`cli.py` (both accounting paths) uses `next((r for r in records if r[0] != h1), None)`;
identical headers make that `None` and charge the whole of R2. Correct **today** only
because `zna.smk` passes `samtools fastq -N`; latent for any other input path, where it
overstates by 12x. Use `records[1]` — guaranteed to exist in that branch, and simpler.

### R4. Duplicated accounting — **FIXED**

Line-for-line equivalent today (verified on 167,784 real pairs across every
process/chunk combination), but every statistic added doubles the surface, and the
parallel-vs-single test enumerates a hardcoded key list. Collapse them —
`_run_single` calling `_process_chunk` on a single chunk — and make the test
`assert s_stats == p_stats`. ~30 lines; needs its own throughput check because the
single path writes incrementally rather than accumulating a blob.

---

## Done — M3/M4

### R5. Make the tool's numbers reach the cohort — **FIXED**

`gather/tools/read_merge.py` keeps only top-level `int`/`float`, so
`score_histogram_bits` — the headline observability addition of the redesign — reaches
the cohort store **nowhere**. Widen to `(int, float, str)`, matching the `params` branch
two lines below, and declare TABLE outputs for histograms. A statistic nobody can query
is not observability.

### R6. Provenance — **FIXED**

There is no `__version__` anywhere under `lib/hulkrna/`. Add `tool`, `tool_version`/git
SHA, `numba` (a numba-less run is ~50x slower and silently correct), `elapsed_s`,
`pairs_per_second`. Config parameters are already cohort-queryable via
`gather/tools/pipeline.py` — **only code identity is missing**, and it is the question
every future corpus defect will open with. `pairs_per_second` as a cohort field gives
free detection of nodes where numba failed to load.

Also in this bucket: **rebin or replace the score histogram**. `score = olen·1.9855 −
diff·8.2143`, so a clean n-base overlap floors to bin `2n−1` for all n < 69 — bins 8–60
hold 8,052 odd counts against 637 even, a **12.6:1 comb**, and bins 1–7 are permanently
dead. Rebinning at 2 bits collapses it to 1.08. Better: emit
`overlap_length_histogram` instead, whose bins are the natural quantum and whose short
tail reads the spurious-detection rate directly. Nothing consumes the score histogram
today, so there is nothing to break — do it before there is.

Free companions in the same edit: `overlap_mismatch_rate` (`sum(diff)/sum(olen)` =
0.00911 against `err_rate` 0.01 — a **sensitive** early-warning channel that moves 5x
under 5% degradation while `merged_pct` moves 1%), and the insert-size distribution
split by outcome, published **with both censoring bounds** (hard-censored at
`len1+len2−15`, hard-floored at `min_read_length`, with 270 counts sitting exactly at 40).

### R7. Throughput — **3 of 4 DONE**, one deliberately deferred

**Done.**

1. **`-T "*"` → the tags `label_defs.yaml` declares** (`zna.smk`). The list is now
   *derived* from the label defs by `zna_label_tags()`, not hardcoded, so the tag set
   and the label set cannot drift, and a label with no `tag:` is a hard error instead of
   a silent column of sentinels. Measured on the production BAM: mean header **153.6 →
   41.5 bytes**, R1+R2 intermediate **30.5 → 23.1 MB (−24%)**, merged output **21.9 →
   17.1 MB (−22%)**. Nothing is lost permanently — the full tag set is archived in
   `rigel/annotated.cram`, so a future khorana label re-encodes from there.
   *Honest correction to the audit:* its "1.29x" was for the `rigel_bam_to_fastq` rule
   (samtools + pigz), not for the merge tool. Merge wall time did not move measurably;
   the win here is bytes through two large intermediates, and it is real.
2. **Reader pigz threads = 1 when workers are resident** (`fastqio`, chosen in
   `cli._run_parallel`). pigz cannot parallelise inflate, so reader threads only contend
   with the workers. Single-process runs still get the caller's thread count.
3. **`--chunk-size` default 50000 → 2000.**

Together with the sparse-histogram fix below, the pool now scales flat: **4.57 / 4.58 /
4.63 µs/pair at 2 / 4 / 8 processes** (was 9.50 / 5.70 / 10.14 — p8 *worse* than p4).

**A regression this work introduced and then fixed, worth recording.** Adding two
histograms while cutting `--chunk-size` 25x made the per-chunk result tuple carry three
dense 1025-element lists through the pickle, and p8 went to 10.14 µs/pair. Histograms
are now shipped **sparse** (`{bin: count}` for nonzero bins only) and folded by
`.items()`. The general lesson for the port: anything added to the `_process_chunk`
return tuple is paid per chunk, and `--chunk-size` multiplies it.

**Deferred: raw-blob IPC** (~1.2x, ~60 lines). The audit's prototype had four defects,
one disqualifying — **R1 shorter than R2 returns OK and silently drops the trailing R2
records**, which the desync check cannot catch because `base_name(h1) != base_name(h2)`
cannot see records that were never read. It is the only remaining item that touches the
reader/worker boundary, where a mistake loses data silently rather than loudly, and it
is exactly the code the C++ rewrite will replace. Doing it now spends the risk twice.
Revisit in zna, with the differential test in place.

---

## Open — pipeline side, not tool side

### R8. A G3 gate that decodes the produced ZNA

> **Single source of truth: [`ZNA_PIPELINE.md`](ZNA_PIPELINE.md) T4.** This item and
> that one are the same instrument, specified twice. Build it from T4, which carries
> the corrected requirement (the gate must DECODE — `zna inspect --counts` provably
> cannot distinguish populated labels from sentinels). Kept here only as the
> merge-tool-side pointer.

The highest-value class fix in the audit, and it belongs with the re-encode rather than
with the tool. `zna inspect --counts --json` output is **byte-identical** between a
correctly-labelled file and one where every label sits at its missing sentinel — it
reports the label *schema* and the flags byte, never a label *value*. So the gate as
specified in `SOS_EOS_ENCODING_PLAN.md` §11 could not have caught the `key:`/`tag:` bug.

The gate must decode: iterate `records(with_ends=True)`, assert the five legal
`(paired, r1, r2, has_start, has_end)` classes and nothing else, assert
`count(R1) == count(R2)`, `count(R1 fwd) == count(R2 rc)` and its mirror, a ~50/50 RC
coin when unstranded, labels not all at sentinel, and totals reconciling against
`read_merge.json` within the N-drop bound. Emit `zna/geometry_check.json` so gather
carries the verdict as a cohort field. ~0.5 µs/record. **Unverified:** the reason to
build it is established from zna source; the gate itself was never run.

Covers in one instrument: the `key:`/`tag:` bug, the missing-`--interleaved` bug,
orphaned mates, broken co-orientation, a biased RC coin, and count divergence — three
of the four historical defects in this chain.

---

## Decided against, with the measurement that closed it

Recorded so nobody re-opens them without new evidence. Full detail in
[MERGE_TOOL_AUDIT.md](MERGE_TOOL_AUDIT.md) §(c).

| Item | Why not |
|---|---|
| **Gate merges on shift ambiguity** | The 12–16% headline is a microsatellite detector, not an identifiability test. Posterior mass off the argmax is **1.60%** of merges (0.15% off ±4). Against STAR ground truth the argmax is right for 90.4%, and on non-repeat disagreements the data prefers the kernel. A wrong shift gives the wrong *length*, not the wrong *sequence* (0.016% of 24-mers absent from the true fragment; the read never overruns the fragment end). The gate would re-emit **+4.68% duplicated bases** — the exact redundancy objective (1) exists to remove. |
| **Per-record `ZE:` geometry tag** | ALREADY_HANDLED. Seam measured sound on 70,351 records: every merged single 2 ends, every mate exactly 1, correct side, 0 broken pairs, 0 lone mates. A two-repo coordinated change for a hole M1.1 closes inside one tool. |
| **§7 composition-aware `p_null`** | Three lenses, all negative. −0.55 pp merge rate, and ~2 pp on exactly the low-complexity reads it was meant to protect. |
| **Per-pair threshold from read length** | Uniform-null p99.9 is flat at ~10 bits from L=50 to L=300. Nothing to gain; the real content was a doc correction. |
| **Estimating `err_rate` from data** | `ê = 0.00866` vs 0.01 assumed. Adopting it costs −0.04 pp merge rate. 0.01 sits inside the flat region. |
| **N special-casing; Smith-Waterman / indels** | N: 0.169% of pairs, every run length exactly 1, 0 merge decisions changed. Indels: the anomalous-mismatch population is 56% two-base repeat, i.e. wrong shifts in repeats, not indels. |
| **SWAR 2-bit packing in numba** | 1.19x on the whole single-process tool, **0 in production**. A null njit kernel with the same signature costs 0.221 µs/pair of pure dispatch, which caps what any numba-side kernel change can return. **Revisit in C++, where that ceiling does not exist.** |
| **Inlining `_shift_score`; hoisting `score_weights`** | 3.4% and 2.5% at `--processes 1`, zero at 8 (workers are idle). |
| **Gating `--correction` on `t_merge`** | REFUTED. 96.6% of sub-threshold corrections sit on the shift STAR independently agrees with, and 86% are in the trim band where R2's copy is discarded — correcting R1 there is the only way R2's evidence survives. Split `bases_corrected` into merged/unmerged instead. |
| **`--self-check` byte invariants; `--count-n`** | The invariants are tautologies against the current construction; the N counter cannot deliver exactness because `--npolicy drop` is pair-atomic. |

---

## Decisions taken (2026-08-12)

1. **The `-T` trim is safe and is now derived from `label_defs.yaml`** (R7.1). Nothing
   consumes the other ten rigel Z tags *today*; they exist for future khorana labels, and
   they survive in `rigel/annotated.cram`, so re-encoding recovers them.
2. **Max-posterior overlap consensus: adopted**, replacing `--correction`. See
   "Overlap consensus" below.
3. **R5/R6/R7 land here, before the port.**
4. Subsumed by (2): `--correction` and its Q14/Q30 cutoffs no longer exist.

## Overlap consensus — landed

Where the mates overlap and disagree, exactly one is wrong, and the two Phred scores
already say which: `P(R1 wrong) = p1(1-p2) / (p1(1-p2) + p2(1-p1))`. The consensus takes
the better-supported call, and writes the *posterior* quality of that call — always
worse than the winner's own Q, because a contested base is less certain than an
uncontested one. Parameter-free.

Measured on 167,784 real pairs (124,248 overlap disagreements), scoring each rule by its
expected residual errors — no ground truth needed, since the posterior gives the
probability that the chosen base is the wrong one:

| rule | expected residual errors | bases rewritten |
|---|---|---|
| R1 always wins | 43,827 | 0 |
| fastp Q14/Q30 (what was deployed) | 16,948 | 27,003 |
| **max-posterior (landed)** | **10,975** | 33,583 |

**35.2% fewer residual errors than the rule that was running.** The gain is entirely in
the band the cutoffs never reached: on NovaSeq binned quality the Q11-vs-Q25 and
Q25-vs-Q37 cells are 5.3% of disagreements with `P(R1 wrong)` of 0.96 and 0.94, and
fastp's gate does not fire on either — so R1's error was kept ~95% of the time.

Quality is recomputed only where the base changed. The agreement case would justify a
*higher* Q, but ZNA stores no quality at all and the merged FASTQ has no other consumer,
so a full per-position recompute would cost throughput for a number nothing reads. If
that ever changes, the agreement case is `Q_out = Q1 + Q2` and belongs in
`_consensus_r1_overlap`.

Deletes `correction`, `_GOOD_QUAL`, `_BAD_QUAL` and the config key — the last of the
knobs §9 of the redesign wanted gone. The tool now has three: `t_merge`, `t_trim`,
`min_read_length`.

---

## After the freeze

### M5 — the port — **done 2026-08-12**

`lib/hulkrna/merge/` is now `src/zna/merge/`, shipped as `zna merge`, with the six
modules unmodified apart from the `hulkrna` → `zna` renames. numba is the optional
`merge` extra; `import zna.merge` works without it and the CLI refuses to run that way
unless given `--allow-slow`. 177 + 6 tests, green with and without numba. Two deviations
from the port plan, both forced by measurement:

- **The argparse half lives in `zna/merge/args.py` and `zna/merge/__init__.py` resolves
  its exports lazily (PEP 562).** The port plan's trap #1 priced eager import at the
  ~4 ms `_DISAGREE_Q` build; the real cost is **~170 ms**, because `overlap.py` imports
  numba. Registering the subcommand eagerly took `import zna.cli` from 40 ms to 210 ms,
  which is absurd for `zna inspect --json`. `zna/cli.py` imports only `args` (argparse
  and nothing else) and reaches for the runtime half after dispatch.
- **The numba refusal sits in the CLI entry points, not in `run()`.** `run()` is the
  in-process API — what the test suite calls, and what `--merge-pairs` will call — so
  making it refuse would have meant either a numba-less CLI test suite that cannot run
  or 30 edited call sites. The failure being guarded (a silently-correct 50x slowdown
  that looks like a slow node) is a job-submission failure, which is a CLI concern.

**Still to do in hulkrna**, and only after a zna release ships the tool: repin
`workflow/envs/read_merge.yaml` to `zna >= <release>` + `pigz`, switch the
`rigel_read_merge` rule to `zna merge`, delete `scripts/hulkrna-merge`,
`lib/hulkrna/merge/` and the two test files. Keep a thin test that the *rule* invokes
the tool correctly — the seam moves, it does not vanish. A brief period of duplication
is much cheaper than a pipeline that cannot build its environment.

### M6 — the C++ backend — **done, shipped in 0.4.0**

8.34 -> 1.40 us/pair (2 threads), 6.0x single-threaded and 1.4x faster than the old
8-process configuration on a quarter of the cores. numba is gone. Output is
byte-identical to the numba implementation on 200,000 real pairs, from either backend,
at any thread count. See [MERGE_CPP_DESIGN.md](MERGE_CPP_DESIGN.md) §16 for what the
measurements corrected along the way.

### Next

1. **`zna encode --merge-pairs`** — merging becomes an input strategy, feeding the
   `(seq, is_paired, is_read1, is_read2, has_start, has_end)` tuple path that already
   exists for ZNA→ZNA re-encode. Deletes the FASTQ intermediate, the `/1`,`/2`-suffix
   re-pairing, and `--treat-unpaired-as-merged` itself. Keep `zna merge` (FASTQ→FASTQ)
   afterwards: other pipelines want the intermediate, and it is the easiest way to A/B
   the merger against fastp/bbmerge.
2. **numba → C++.** Keep the pure-Python scorer as the *reference backend* (zna already
   has this pattern: `_pycodec` vs `_accel`, pinned by cross-backend equivalence tests),
   not as a fallback. Move to fixed-point integer scoring so the `argmax` is exactly
   reproducible across compilers and platforms — this is training data. Any differential
   test must be **order-independent** (assert the returned shift is *an* argmax with
   matching n/d/score): ties occur on 0.845% of real pairs, and an order-replicating
   reference would break on exactly the SIMD rewrite it exists to guard.
