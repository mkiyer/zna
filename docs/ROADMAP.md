# ZNA roadmap

What is scheduled, what is being considered, and — the section most worth reading before
proposing something — what has already been tried and closed by measurement.

Shipped work is not described here; see [CHANGELOG.md](../CHANGELOG.md). How the current
algorithms work is [METHODS.md](METHODS.md).

---

## Scheduled

### 0.4.2 — x86 tuning for the merge kernel

*(Was scheduled for 0.4.1, which went to the format-version-3 break instead.)*

The scan uses a byte-wise 16-byte vector comparison, which is baseline on every target:
NEON on aarch64, SSE2 on x86-64 — no `-march` flag, no runtime dispatch. That is
deliberate: one kernel, correct and fast everywhere, with no ISA machinery.

But **every measurement behind it was taken on aarch64/NEON.** Three things are therefore
unverified, and all three want a Linux/x86 box:

1. **Confirm the SSE2 path executes.** It is written (`neq16` under `#ifdef`) and has
   never run. `scripts/merge_bench/bench_simd.cpp` builds it.
2. **Evaluate AVX2 (32-byte vectors) before writing any dispatch.** Do not assume it
   wins. The scan is *rejection-dominated* — roughly 0.32·n comparisons per rejected
   shift — so the variable that matters is the bail interval, not the vector width. On
   NEON, bailing every 32 bases beat both 16 and 64, and a packed kernel's 128-base bail
   was measurably worse than its 64. A 32-byte vector may need to bail every single
   vector just to compete.
3. **Re-tune the bail granularity.** It is a hardware property, not an algorithm
   property. 32 bases won on Apple silicon; Zen and Sapphire Rapids may not agree.

Reproduce with [`scripts/merge_bench/README.md`](../scripts/merge_bench/README.md).

### 0.5.0 — `zna encode --merge-pairs`

Merging in-process, deleting the FASTQ intermediate between `zna merge` and
`zna encode`. Fully specified in [MERGE_PAIRS_PLAN.md](MERGE_PAIRS_PLAN.md), including
the three blockers found while auditing the design against the code.

Worth knowing before prioritising it: **the corpus it produces is flag-identical to the
two-step path** on well-formed input — verified over 1,416,630 records, same histogram
either way. So this is a speed and robustness change, not a correctness one, and nothing
downstream needs rebuilding on account of it.

### Faster inflate

Gzip decompression is the floor once compute is native. `libdeflate` or ISA-L inside the
extension would move it — but `--merge-pairs` may moot the question first by deleting the
intermediate that is being decompressed.

---

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

### Stored block index (sidecar)

`block_index()` walks the block-header chain at ~2.3 µs per block — ~1.4 ms for a
611-block file — so a stored index buys nothing locally.

It would matter for **remote** corpora. A sidecar `sample.zna.idx` of
`(offset, n_records)` per block would let a consumer read record counts for corpus
balancing without downloading the data, and issue ranged reads for only the blocks it
sampled — turning a 5% sample into 5% of the transfer.

A sidecar beats a footer: no format-version bump, old readers ignore it, and it can be
generated lazily. A footer would be misread, because today's readers parse block headers
until EOF.

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

## Known divergence between the codec backends

The project's central discipline is that the two codec backends agree exactly. One case
does not, and it is pinned by a test that documents it rather than blessing it
(`tests/test_fuzz_roundtrip.py::test_npolicy_rejects_non_n_ambiguity_codes_in_python_backend`):

With any `--npolicy` set, the compiled backend substitutes the policy base for **any**
unencodable character, while the pure-Python backend substitutes only `N`/`n` and raises
on the rest. So an IUPAC ambiguity code — `R`, `Y`, `S`, … — encodes on the compiled
backend and raises on the reference one.

Neither behaviour is obviously right, which is why it is still open. `--npolicy` is
named for N, which argues for rejecting the rest; but silently rejecting a whole library
because one read carries an `R` is a harsh failure. Resolving it is a behaviour change
either way, so it wants a release boundary and a decision, not a quiet fix.

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
  which tags were declared as labels. See [MERGE_PAIRS_PLAN.md](MERGE_PAIRS_PLAN.md) §3.5.

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
