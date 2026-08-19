# ZNA roadmap

What is scheduled, what is being considered, and — the section most worth reading before
proposing something — what has already been tried and closed by measurement.

Shipped work is not described here; see [CHANGELOG.md](../CHANGELOG.md). How the current
algorithms work is [METHODS.md](METHODS.md).

---

## Scheduled

### 0.5.3 — re-measure the merge kernel on aarch64

**0.5.2 changed the NEON code path without being able to measure it, and that is a debt,
not a claim.** The x86 work replaced `neq16` with `neq16x<V>`, which accumulates the
per-lane equality flags in a vector register and folds **once per group** instead of once
per vector. On NEON that turns `V` × `vaddvq_u8` into `V-1` × `vaddq_u8` plus one
`vaddvq_u8`, and `vaddvq_u8` is the more expensive instruction on every ARM core — so it
should be a win or a wash, but "should be" is exactly the reasoning this ROADMAP exists to
distrust. Correctness is not in question: the outputs are byte-identical on 1M pairs, the
suite passes under UBSAN, and `asan_scan` is clean.

Two things to take on an Apple-silicon or Graviton box:

1. **Confirm the grouped reduction did not regress NEON.** `scripts/merge_bench/bench_simd.cpp`
   against the shipped kernel.
2. **Re-tune `BAIL_VECTORS` for aarch64.** It is still **2** (32 bases), the optimum
   measured for the *old* per-vector reduction. x86 moved from 2 to 3 once the reduction
   got cheaper and the 16-byte step loop shortened the tail; aarch64 may well move too.

### 0.6.0 — `zna encode`'s per-record driver into C++, which also lets `pigz` go

Two items in one because the second falls out of the first for free, and the second is
the reason not to treat 0.5.2's gzip work as finished.

**ZNA now has four ways to inflate a `.gz`** — ISA-L, a `pigz` subprocess, a `gzip`
subprocess, and the stdlib — which is one more than 0.5.1 had and three more than anyone
wants. The end state is **two**: ISA-L, and the stdlib behind it. What stands in the way
is not the decompressor, it is *how `zna encode` reads*.

`zna merge` hands raw buffers to a C++ chunk parser and reads them in big blocks on a
prefetch thread, so an in-process inflate that releases the GIL overlaps perfectly and
ISA-L wins outright. `zna encode` drives `readline` from the main thread, so its inflate
competes with a GIL-bound parse loop and a separate *process* beats a faster
*implementation* — 5.50 s against 6.54 s on 1M records. That is why `pigz` is still first
there and cannot simply be deleted.

**A prefetch thread does not rescue it, and this was measured rather than assumed.**
Wrapping ISA-L in a Python prefetching reader recovers only 1.08x (3.05 s → 2.85 s) and
is still 0.79x of the `pigz` pipe: the `readinto` hop and the buffer copy give back what
the overlap wins, and the GIL is still in the loop. There is no arrangement of threads
that fixes a main-thread `readline` parse.

So the ordering is: **move the encode driver to C++ first** — which is worth doing on its
own evidence, since the compiled codec is ~3% of `zna encode` and the Python per-record
work is most of the rest (see the numbers below) — and *then* encode reads chunks like
merge does, ISA-L wins there too, and the two subprocess paths, `_DecompressPipe`,
`_GUNZIP_COMMANDS`, the child reaping and the SIGPIPE special-casing all delete. 0.5.1
had to fix a leaked `pigz` child and prefetch thread; that whole class of bug goes with
them.

Until then 0.5.2's four paths are a deliberate way station, not the design.

#### Why the codec itself is not the target



Recorded here because the obvious guess is wrong and cost this round a detour. Profiling
`zna encode` found the **compiled codec at ~3% of the time**: zstd is another 3%
(`--uncompressed` 4.95 s against `--level 9` 5.10 s on 1M records), the Python FASTQ
parse is ~27%, and the rest is the per-record driver — `_read_key`, `write_record`,
`_trim3`, `_pair_interleaved`, `_fragment_units`, one Python call each per record.

**And the codec kernels have no SIMD headroom left**, which is worth knowing before
someone repeats the merge fix on them. Measured on 30 Mbase, baseline flags:

| kernel | shipped scalar | hand-written SSE2 |
|---|---|---|
| decode (packed → ASCII, `DECODE_LUT4`) | **5465 Mbase/s** | 2679 (0.49x) |
| encode (ASCII → packed, table) | **1723 Mbase/s** | 1527 (0.89x) |

Decode at 5465 Mbase/s is ~5.5 GB/s of stores — already memory-bandwidth-bound on this
part, so no kernel can help. Encode is limited by four dependent table lookups per output
byte, and beating that needs `pshufb` (SSSE3, not baseline) for a ~1.5% CLI win. Both
left alone.

The real fix is moving the per-record loop into C++, which is what `--merge-pairs`
already does for the paired path — a 0.6.0-shaped change, not a patch release.


## Considered, not scheduled

### `zna sample --fraction`

A command writing a new ZNA holding a random fraction of an input's blocks or records.

```bash
zna sample --fraction 0.05 big.zna -o small.zna
zna sample --records 1000000 big.zna -o small.zna
```

**Why it is not the first answer for training.** Training iterates the same corpus for
many epochs and wants *different* reads each pass; a static subset freezes one draw. The
runtime path already covers that — `block_index()` to size the file, `blocks(indices=…)`
to draw a fresh block subset per epoch — without materialising anything.

It is still worth having for publishing a small public excerpt, cheap CI fixtures, and
pinning an evaluation set that must not vary between runs.

**Sketch.** Block-granular sampling is nearly free: copy whole compressed payloads
through without decoding, rewriting only block headers. Record granularity needs a
decode/re-encode pass. Start with blocks, and treat `--records N` as rounding to whole
blocks. Both require a shuffled input to be statistically meaningful.

Since format version 3 a block holds whole fragments, so the block-granular path cannot
orphan a mate and needs no fragment logic of its own — it copies payloads and is done.
The record-granular path would have to group mates itself, which is the second reason to
start with blocks.

### Stored block index (sidecar) — SUPERSEDED by the 0.5.0 trailer

The block index now lives **inside the file**: the 0.5.0 trailer stores it, the footer
makes it discoverable O(1) from EOF, and a suffix-range GET of a remote file's tail
indexes it with no full download — the remote-corpus case this item existed for.
`docs/TRAILER_PLAN.md` §14 records how the footer landed *without* a format-version
bump (count-0 pseudo-blocks; 0.4.1 readers decode them as valid empty blocks), which
dissolves the "a footer would be misread" objection that justified a sidecar. No
sidecar will be built.

### String labels

The fixstr (`Z`) dtype was dropped to simplify the initial format. Three ways back:

| Option | Shape | Cost |
|---|---|---|
| **A — quoted strings** | `MD:Z:"50 A2"`; unquoted values still split on whitespace | minimal format change, backward-compatible |
| B — binary-safe fixed-width | programmatic API only, null-padded; headers stay numeric-only | no parser change, no header support |
| C — variable-length (`V`) | 2-byte length + UTF-8 | needs a different block layout (offset table); most flexible, most complex |

Start with **A** when string labels are actually needed.

### Categorical / dictionary-encoded labels

Tags like `tp:A` have a tiny alphabet. A per-block dictionary plus per-record indices
compresses ~8× for a 4-value categorical against raw `uint8`. New dtype code `E`,
dictionary written before the column, falling back to raw storage above a cardinality
threshold.

### Per-column compression

Today a block is one ZSTD frame over flags ‖ labels ‖ lengths ‖ sequences. Compressing
each column independently would allow per-column strategies (RLE for flags, delta for
sorted labels, plain ZSTD for near-random sequences) and — the real win — **selective
column reads**, skipping the sequence payload entirely when only labels are wanted. The
block header would carry per-column compressed sizes.

### Label indexes

For "all reads where `AS > 100`": store `(min, max)` per label per block in the block
header and skip blocks whose range cannot match. Zero cost at write time.

### Label arithmetic

Some labels are derivable rather than stored — `insert_size`, an alignment-score gap, an
is-mapped predicate. A small expression engine could compute them during encoding.

### Arrow / Parquet export

ZNA already stores labels columnar, so exporting them is mostly a metadata mapping —
enabling zero-copy pandas/polars and SQL over label columns.

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

## Open questions in the label model

Two cases are currently *defined but not validated*, and both concern paired input:

- **R1 and R2 with different tag subsets.** Each mate is read independently, so a tag
  present on only one resolves to that label's `missing` value on the other. Nothing
  reports it.
- **Conflicting tag values on a merged read.** `zna merge` builds the merged record's
  name from R1's header, so R1's tags survive and R2's are discarded. That is right for a
  per-*fragment* tag (the production case — transcript and gene assignment agree across
  mates by construction) and lossy for a genuinely per-*mate* one.

  The subtle part: **whether a fragment contributes one label row or two is decided by the
  overlap score**, so for a per-mate tag the column's meaning varies with sequence content
  and with `--threshold-merge`. Nothing records which branch a fragment took.

  `zna merge` cannot validate this — it is a FASTQ tool with no concept of label
  definitions, and comparing all tags would false-positive on legitimately per-mate ones.
  The check belongs in `--merge-pairs`, where the encoder holds both headers *and* knows
  which tags were declared as labels. See [MERGE_PAIRS_PLAN.md](archive/MERGE_PAIRS_PLAN.md) §3.5.

---

## Decided against, with the measurement that closed it

Recorded so nobody reopens these without new evidence.

| Proposal | Why not |
|---|---|
| **Gate merges on shift ambiguity** | Posterior mass off the argmax is 1.60% of merges (0.15% off ±4). A wrong shift gives the wrong *length*, not the wrong *sequence* — 0.016% of 24-mers absent from the true fragment, and the read never overruns the fragment end. The gate would re-emit **+4.68% duplicated bases**, which is the exact redundancy the trim band exists to remove. |
| **Per-record geometry tag in the FASTQ** | The seam was measured sound on 70,351 records: every merged single reports two ends, every mate exactly one, correct side, zero broken pairs, zero lone mates. A two-repo coordinated change for a hole the tool already closes. |
| **Composition-aware null model** | Three independent lenses, all negative: −0.55 pp merge rate overall, and ~2 pp on exactly the low-complexity reads it was meant to protect. |
| **Per-pair threshold from read length** | The uniform-null p99.9 is flat at ~10 bits from L=50 to L=300. Nothing to gain. |
| **Estimating `err_rate` from the data** | Measured `ê = 0.00866` against the assumed 0.01; adopting it costs −0.04 pp merge rate. 0.01 sits inside the flat region. |
| **Special-casing N; Smith-Waterman / indels** | N appears in 0.169% of pairs, every run length exactly 1, and changed **0** merge decisions. The anomalous-mismatch population is 56% two-base repeat — wrong shifts in repeats, not indels. |
| **Raising `--threshold-trim` above 8** | The two harms added (duplicated bases left + real bases deleted) are minimised at 8: 85,611 against 97,053 / 127,398 / 167,445 at 12 / 16 / 20. It keeps minimising unless a deleted base is judged >1.63× as harmful as a duplicated one — and a duplicated base is *false* evidence while a deleted one is merely absent. |
| **Raising `--threshold-merge` above 28** | 28 minimises false positives plus false negatives, and not marginally: 6,603 errors per million against 44,145 at the fastp-equivalent setting. Raising it trades ~11 missed merges per wrong merge prevented. See [MERGE_BENCHMARK_RESULTS.md](MERGE_BENCHMARK_RESULTS.md) §6 for the full curve — it is a trade-off with a number, not a recommendation. |

**Where these numbers come from.** The last two rows — the two thresholds — are from the
simulated ground-truth benchmark and are fully reproducible; the derivations and the
input are in [MERGE_BENCHMARK_RESULTS.md](MERGE_BENCHMARK_RESULTS.md), and
[`scripts/merge_bench/`](../scripts/merge_bench/) regenerates them. The rest were
measured during a 2026-08 audit against **production libraries that are not
distributable**, so the figures are recorded here rather than reproducible from this
repository. Treat them as evidence that the question was asked and answered, not as
something you can re-run — and if you have a reason to reopen one, measure it on your own
data rather than arguing from these.

---

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
