# ZNA roadmap

What has already been tried and closed by measurement — the section most worth reading
before proposing something. Open work, whether scheduled, considered or an open question,
is in [GitHub Issues](https://github.com/mkiyer/zna/issues).

Shipped work is not described here; see [CHANGELOG.md](../CHANGELOG.md). How the current
algorithms work is [METHODS.md](METHODS.md).

---

## Moved to GitHub Issues

On 2026-11-13 the open items this file held became issues, each with its text and a line
naming the section it came from:

| Section here | Issue |
|---|---|
| Fix — an exact reverse-complement pair is not proof of a full overlap | [#1](https://github.com/mkiyer/zna/issues/1), fixed in 0.7.0 |
| Next — `zna encode`'s per-record driver into C++, which also lets `pigz` go | [#3](https://github.com/mkiyer/zna/issues/3) |
| Faster table building for very long reads | [#4](https://github.com/mkiyer/zna/issues/4) |
| `zna sample --fraction` | [#5](https://github.com/mkiyer/zna/issues/5) |
| String labels | [#6](https://github.com/mkiyer/zna/issues/6) |
| Categorical / dictionary-encoded labels | [#7](https://github.com/mkiyer/zna/issues/7) |
| Per-column compression | [#8](https://github.com/mkiyer/zna/issues/8) |
| Label indexes | [#9](https://github.com/mkiyer/zna/issues/9) |
| Label arithmetic | [#10](https://github.com/mkiyer/zna/issues/10) |
| Arrow / Parquet export | [#11](https://github.com/mkiyer/zna/issues/11) |
| Open questions in the label model | [#12](https://github.com/mkiyer/zna/issues/12), [#13](https://github.com/mkiyer/zna/issues/13) |

---

## Superseded

### Stored block index (sidecar) — SUPERSEDED by the 0.5.0 trailer

The block index now lives **inside the file**: the 0.5.0 trailer stores it, the footer
makes it discoverable O(1) from EOF, and a suffix-range GET of a remote file's tail
indexes it with no full download — the remote-corpus case this item existed for.
`docs/TRAILER_PLAN.md` §14 records how the footer landed *without* a format-version
bump (count-0 pseudo-blocks; 0.4.1 readers decode them as valid empty blocks), which
dissolves the "a footer would be misread" objection that justified a sidecar. No
sidecar will be built.

---

## Known divergence between the codec backends — RESOLVED in 0.4.0

This section used to record a live inconsistency: the compiled backend substituted the
policy base for **any** unencodable character while the reference substituted only
`N`/`n` and raised on the rest. 0.4.0's alphabet strictness (`docs/archive/NPOLICY_PLAN.md`
§8.2) resolved it in the strict direction — **both backends now raise, naming the
character and its offset**, and only `N`/`n` is ever substitutable, under every policy.
Verified directly against both backends (`encode_block` with an `R` under
`--npolicy random`: both raise) while auditing for the 0.5.0 trailer work, which found
this section still claiming the divergence was open. No divergence between the codec
backends is currently known.

---

## Decided against, with the measurement that closed it

Recorded so nobody reopens these without new evidence.

| Proposal | Why not |
|---|---|
| **Gate merges on shift ambiguity** | Posterior mass off the argmax is 1.60% of merges (0.15% off ±4). A wrong shift gives the wrong *length*, not the wrong *sequence* — 0.016% of 24-mers absent from the true fragment, and the read never overruns the fragment end. The gate would re-emit **+4.68% duplicated bases**. (0.6's plausibility gate is a different test: it asks whether the winner's *mismatches* are believable, not whether its score is well separated.) |
| **Per-record geometry tag in the FASTQ** | The seam was measured sound on 70,351 records: every merged single reports two ends, every mate exactly one, correct side, zero broken pairs, zero lone mates. A two-repo coordinated change for a hole the tool already closes. |
| **Composition-aware null model** | Three independent lenses, all negative: −0.55 pp merge rate overall, and ~2 pp on exactly the low-complexity reads it was meant to protect. |
| **Special-casing N in the score; Smith-Waterman / indels** | N appears in 0.169% of pairs, every run length exactly 1, and changed **0** merge decisions. The anomalous-mismatch population is 56% two-base repeat — wrong shifts in repeats, not indels. (0.6 discounts a one-sided N in the plausibility gate only, where it did matter — see below; the score is unchanged.) |

**Where these numbers come from.** They were measured during a 2026-08 audit against
**production libraries that are not distributable**, so the figures are recorded here
rather than reproducible from this repository. Treat them as evidence that the question
was asked and answered, not as something you can re-run — and if you have a reason to
reopen one, measure it on your own data rather than arguing from these.

### The 0.6 merge policy: what was measured and not built

The alternatives weighed while designing 0.6's merge decision (`docs/archive/MERGE_ACCURACY_PLAN.md`
has the design; its §10 evidence directory has every run). **These numbers are from
simulation** — the 46-substrate truth panel (14.1M simulated pairs: chr22 and
transcriptome sets from khorana's simulator and ZNA's own, hg38, read-through and
3′-degradation stress sets) and its gene-disjoint holdouts — except where a row cites the
2026-08 audit, because only simulation knows which merges are wrong. They say which rule is better on those substrates; they
are not production rates.

| Proposal | Why not |
|---|---|
| **Estimate `--error-rate` from each library** | Built and measured, then dropped. Fed into the score, a noisy library's own higher rate cheapened mismatches and admitted *more* divergent repeats and short false read-throughs: 3′-degraded sets 560 → **741** wrong merges at the estimated 0.035, against 560 → **240** at a fixed 0.01; lost fragments rose on a 1%-error hg38 set; clean libraries estimated low and refused short true overlaps carrying two errors. It also needed a buffered sample, a prior and a sorted-input caveat. The 2026-08 audit reached the same place from production data (`ê = 0.00866` against 0.01, −0.04 pp merge rate: 0.01 sits in the flat region). The rate is now a documented setting that each run *checks* (`expected_refused_true_overlap_fraction`). |
| **Re-place an implausible winner** (the best plausible shift, or an argmax constrained to plausible shifts) | An implausible winner marks a repetitive context, and the runner-up there is usually another repeat: re-placing moved **37–38%** of the wrong merges it caught onto *another* wrong shift on the gene-disjoint holdout (9% on dev). Abstaining — keep the pair whole — removed ~10% more wrong merges than re-placing on dev and added no lost fragments. Its cost is 145–271 re-merges forgone per 500k–1M pairs where the truth was outscored. |
| **Read-through concordance rule** (refuse an `s < 0` shift under 80 bits unless the two 8-base adapter overhangs agree within 2 mismatches) | Four tuned constants standing in for a fact the user already knows. Where the input is honestly trimmed, declaring it (`--adapter-trimmed`) beat the best rule: transcriptome fastp set **−715** wrong merges / **+332** correct against −632 / +249. Where adapters are not shared between mates (small-RNA kits) the rule is wrong by construction: it turned 18,201 of 100k raw small-RNA pairs into kept pairs carrying 4.45M adapter bases, and its strict form refused 73,150 correct merges. |
| **A capped trim band** (keep trimming, but only overlaps with at most `dtrim[n]` mismatches, `dtrim` at α = 10⁻³) | Cut wrong trims by only 17–25% (hg38 2,948 → 2,267, still deleting 19,013 real positions), refused 245 exact trims on a 3′-degraded set, and interacted with the gate: re-placed wrong merges landed in the band (74 new wrong trims on hg38). Removing the band gives 0 wrong trims and 0 kept-mate substitutions by construction, at ~0.55 duplicated bases per pair — and khorana wants unmerged mates whole. |
| **Displaced-winner clauses** (refuse a trim or read-through that won only because the true winner was masked) | Every variant that refused a displaced read-through winner gave *identical* merge-side results (−3,140 wrong merges on hg38); the variants differed only in trading wrong trims against duplicated bases (e.g. hg38 −681 wrong trims for +428 duplicate positions). With the trim band gone there is nothing left for them to decide. |
| **The "lazy" kernel form** for those clauses (a capped main scan plus verification re-scans) | Measured at **2.4×** 0.5.3's scan time on hg38 and **3.6×** on read-through-rich input, where 92% of pairs re-scanned; the mismatch caps it relied on saved only 8% of shifts once bail is per 64-base group. An exact one-pass form reached 0.91–0.98×, but the 0.6 gate needs neither: it is one table lookup after an ordinary scan (`docs/METHODS.md` §2.5). |
| **Warn when the detected rate exceeds `--error-rate`** | The first trigger fired on 12 of the 46 panel benches, ten of them expecting under 0.02% of true overlaps refused — repeats nudge the detected rate past 0.01 on clean libraries while costing nothing. Replaced by the expected refused share itself, warned above 0.1%: only the extreme 5% 3′-ramp set crosses it (0.40%). |
| **A faster `--adapter-trimmed` check** | The check's second, unrestricted scan costs +5.4% on the first 100,000 pairs of a declared run only — about 20 ms. A cheaper existence test would need its own proof of equivalence and a matching reference implementation. |
| **Raising the merge floor** (`--alpha` below 10⁻⁶; 0.5.x's `--threshold-merge` above 28) | On the 0.5.x benchmark 28 bits minimized false plus missed merges — 6,603 per million against 45,616 at the fastp-equivalent 60 bits (fastp itself: 44,145) — and raising it traded ~11 missed merges for each wrong merge prevented, reaching 22.5:1 at 100 bits. The gate has since removed about half of those wrong merges. See [MERGE_BENCHMARK_RESULTS.md](MERGE_BENCHMARK_RESULTS.md) §6: a trade-off with a number, not a recommendation. |
| **A per-pair floor *tuned* to read length** | The 2026-08 audit found the uniform-null p99.9 flat at ~10 bits from 50 to 300 bp, so nothing to gain by tuning. 0.6 does vary the floor with read length, but by derivation, not tuning: `log₂(N/α)` is the same union bound at every geometry, 26.6 bits at 2×50 to 29.2 at 2×300. |

---

## Closed by measurement in 0.6

### The merge-accuracy work — DONE in 0.6.0

khorana's chr22 merge review (2026-09-25) is answered by the 0.6 policy: a floor derived
per pair from one tolerance, a plausibility gate at the same tolerance, an
`--adapter-trimmed` declaration, and no trim band (`docs/METHODS.md` §1). Qualified at
full scale against 0.5.3, pair by pair against simulated truth over 4.11M pairs: wrong
merges 11,875 → 5,151, lost fragments 719 → 544, wrong trims and kept-mate
substitutions to 0, 7 correct merges forgone and 502 gained, compiled CPU 9–29% lower
(`docs/MERGE_BENCHMARK_RESULTS.md` §9). Every acceptance criterion fixed before the run
held, except where the adapter declaration was not fully honest: declaring the output
of hulkrna's fastp pass forgoes 0.011% of correct merges and emits 47 adapter bases per
500k pairs, because its `--cut_tail` hides a 1–2 bp overhang from fastp's adapter
detection and leaves 1–3 bp of adapter on 42 pairs. The recommended trimmer is the same
single pass without `--cut_tail` (fastp 1.1.0 and 1.3.6 identical): it leaves 1–3 bp on
7 pairs per 500k and, declared, emits 8 adapter bases in 3 of them (forgone 0.002%,
within the limit). Hence the declaration is documented as for honest input only, with
that residual stated.

**Left open, by design.** Near-identical repeats still merge wrongly (a perfect 15-base
repeat is plausible under any error model), and the residual is more concentrated than
0.5.3's: on transcriptome input one transcript carries a third of it. Distinguishing
those needs information the merger does not have.

---

## Closed by measurement in 0.5.3

### ~~0.5.3 — re-measure the merge kernel on aarch64~~ — DONE

The debt 0.5.2 filed, paid on an Apple M3 Max against 50k real pairs
(`scripts/merge_bench/bench_simd.cpp`, which now carries the shipped kernel shape beside
the 0.5.1 "before"). Both questions answered:

**1. Did the grouped reduction regress NEON?** No — it won. Per-vector kernel at bail 32
(0.5.1): 0.497 µs/pair; grouped at bail 32: 0.473. `vaddvq_u8` really is the expensive
instruction, and folding once per group pays even on a core where it is cheap.

**2. Does `BAIL_VECTORS` move?** Yes, further than x86's: **2 → 4** (64 bases). The
shipped 2 was the worst of {2,3,4,5}:

| bail (bases) | µs/pair | vs scalar |
|---|---|---|
| 32 (0.5.2 shipped) | 0.506 | 2.33× |
| 48 (x86's value) | 0.439 | 2.69× |
| **64** | **0.416** | **2.84×** |
| 80 | 0.415 | 2.75–2.84× |
| 96 | 0.528 | 2.39× |
| 128 | 0.675 | 1.87× |

The cliff past 80 is the step loop's limit: a 96+ base group fits a 150 bp shift at most
once, so the group loop stops contributing. 64 and 80 are within noise across repeats;
64 keeps clear of the cliff down to ~80 bp reads. Byte-identical output at every value
(a pruning schedule cannot change an answer); 1.22× over shipped on the kernel, within
noise end to end because the kernel is now ~25% of a merge's wall time.

The Windows leg 0.5.2 could not build locally is green in CI on its commit, so the
`psadbw` portability claim is verified, not assumed.

## Closed by measurement in 0.5.2

The x86 assessment this ROADMAP had been carrying since 0.4.1. Kept in full because
the headline finding is a trap someone will otherwise fall into again: **the SSE2 path
was running the whole time, and was slow for a reason that had nothing to do with
vector width.**

### ~~0.4.2 — x86 tuning for the merge kernel~~ — DONE in 0.5.2

Taken on a Xeon E5-2680 v3 (Haswell: AVX2, no AVX-512), RHEL 8, g++ 8.5 and 13.2,
against the 1M-pair library from `scripts/merge_bench/gen_library.py`. The three
questions this item asked, and their answers:

**1. Does the SSE2 path execute?** Yes — `VECTOR_WIDTH` reported 16 and `pcmpeqb` was in
the disassembly. It was also **crippled**, which is not what anyone was looking for.
`neq16` reduced each vector compare with `pmovmskb` + `__builtin_popcount`, and POPCNT is
SSE4.2-era while this extension is compiled for baseline x86-64 — so GCC emitted
`callq __popcountdi2@plt`, a PLT call with the caller's registers spilled around it,
**twice per 32 bases, in the innermost loop of the whole kernel**. GCC 8.5 and GCC 13.2
both do this; the comment claiming GCC would "emit a table or a SWAR sequence" instead
was simply wrong. aarch64 never saw it because NEON reduces with `vaddvq_u8`.

Measured cost: the byte-SIMD kernel was only **1.33x** the scalar one on x86, against the
**2.3x** recorded on aarch64. Nearly half the SIMD win was going to a function call.

**2. Does AVX2 win?** Barely, and not enough to build. The fix for (1) is `psadbw` —
baseline SSE2, no `-march`, no dispatch, and the exact analogue of `vaddvq_u8`. With the
reduction fixed, 32-byte AVX2 vectors are **within noise of 16-byte SSE2**:

| kernel | flags | µs/pair |
|---|---|---|
| shipped 0.5.1: `pmovmskb` + popcount, bail 32 | none | 3.140 |
| `psadbw`, bail 48, + 16-byte step | **none** | **1.539** |
| `pmovmskb` + popcount, bail 32 | `-mpopcnt` | 1.739 |
| AVX2 32-byte, `pmovmskb` + popcount, bail 32 | `-march=native` | 1.556 |
| AVX2 32-byte, `psadbw`, bail 32 | `-march=native` | 1.776 |

The baseline-SSE2 kernel at 1.539 **beats the AVX2 one at 1.556**. So the ROADMAP's
instinct was right and the reason is now measured: the win was never in the vector width,
it was in the reduction. No AVX2 path, no CPUID, no runtime dispatch, no second kernel —
and note the last row, where AVX2 *loses* to SSE2 because a 256-bit `psadbw` needs a
cross-lane fold that POPCNT does not.

**3. Re-tune the bail granularity.** It is per-ISA now (`BAIL_VECTORS` in
`merge_core.hpp`): **3 vectors / 48 bases on x86**, against 1.754 at 2 and 1.542 at 4;
**2 / 32 bases on aarch64**, left at its measured Apple-silicon optimum because no
aarch64 machine was available this round.

A wide interval only pays off with the second loop that came out of this: group step,
then a **16-byte step**, then the byte tail. Without it a 48-base group leaves up to 47
bases to the scalar loop — a third of a 150 bp read — and the cheaper reduction is spent
on a longer tail. With it the tail is under 16 bases at any interval, worth a further
1.08x (1.664 → 1.539).

End to end on 1M pairs, with the merged FASTQ byte-identical to 0.5.1's and every
compressed block payload of the `.zna` bit-identical (the files themselves differ only in
the prologue's `writer_version`). Same decompressor on both sides, output to `/dev/null`,
so this row is the kernel and nothing else: `merge --threads 1` **7.61 s → 5.84 s** wall
and 7.77 s → 6.06 s CPU; `merge_chunk` itself 4.698 → 2.978 µs/pair (1.58x).

### ~~Faster inflate~~ — DONE in 0.5.2, as an optional extra

Now that the kernel is 2x, this became the largest single cost of a merge: a profile of
`zna merge --threads 1` put **a quarter of all cycles inside pigz's libz**, `crc32_z`
alone at 8.6%.

`pip install 'zna[fast]'` pulls in `isal` (Intel ISA-L), and the readers prefer it.
Raw inflate on a 106 MB member: stdlib 208 MB/s, `pigz -dc -p1` 193 MB/s, **ISA-L 448
MB/s**. It is strictly optional — absent, 0.5.1's behaviour is unchanged.

**Which decompressor wins depends on the caller, which was not obvious and had to be
measured** (`zna/_gzip.py` holds both numbers):

| 1M pairs, wall / CPU | ISA-L | pigz subprocess |
|---|---|---|
| `zna merge --threads 1` | **4.04 / 4.28 s** | 5.84 / 6.06 s |
| `zna encode` from `.gz` | 6.54 / 6.31 s | **5.50 / 7.26 s** |

`zna merge` prefetches on its own thread and ISA-L releases the GIL, so in-process
inflate overlaps exactly as a subprocess would and being 2.3x faster it wins on **both**
axes. `zna encode` drives `readline` from the main thread, so its inflate competes with a
GIL-bound parse loop and another *core* beats a cheaper *implementation* — pigz stays
first there, with ISA-L ahead of the stdlib as the fallback. All three paths are pinned
byte-identical, and all three still reject a truncated or CRC-damaged member
(`tests/test_gzip_backends.py`).

Cumulative for 0.5.2, 1M pairs, `merge --threads 1`, output to `/dev/null`: **7.61 →
4.04 s** wall (1.88x) and **7.77 → 4.28 s** CPU (1.82x). Writing the 427 MB intermediate
FASTQ to network scratch instead: 9.75 → 6.88 s wall, 7.71 → 4.32 s CPU. `encode
--merge-pairs`, which writes no intermediate, 9.63 → 6.23 s.
