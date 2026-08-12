# ZNA Future Extensions

*Ideas for future development beyond the current performance optimization roadmap.*

---

## Static Subsampled Files (`zna sample --fraction`)

A CLI that writes a new ZNA holding a random fraction of an input's blocks or
records:

```bash
zna sample --fraction 0.05 big.zna -o small.zna
zna sample --records 1000000 big.zna -o small.zna
```

The motivating case is corpus balancing for LLM training, where a manifest mixes
files spanning several orders of magnitude in size and the large ones dominate
the sample. Building balanced static subsets offline would make the training
run's data pipeline trivial.

**Why it is not the first answer.** Training iterates over the same corpus for
many epochs and wants *different* reads each pass. A static subset freezes one
draw, so every epoch sees the same reads — the opposite of what is wanted. The
runtime path (`block_index()` to size the file, `blocks(indices=...)` to draw a
fresh block subset per epoch) covers that case without materialising anything.

`zna sample` is still worth having for: publishing a small public excerpt of a
large dataset, cheap CI fixtures, and pinning an exact evaluation set that must
not vary between runs.

**Design sketch.** Block-granular sampling is nearly free — copy whole
compressed payloads through without decoding, rewriting only the block headers.
Record-granular sampling requires a decode/re-encode pass. Start with block
granularity and `--records N` as an approximation that rounds to whole blocks.
Both require the input to be shuffled to be statistically meaningful, so the
command should read the file's own history if the format ever records it, and
warn otherwise.

---

## Stored Block Index (sidecar or footer)

`ZnaReader.block_index()` walks the block-header chain, seeking over each
payload: about 2.3 microseconds per block, so ~1.4 ms for a 611-block file.
That is cheap enough that a *stored* index buys nothing for local files.

It would matter for remote corpora. With a sidecar `sample.zna.idx` holding
`(offset, n_records)` per block, a consumer could:

- read record counts for corpus balancing **without downloading the data file**,
  which is the real cost when files arrive over Globus or S3;
- issue ranged reads for only the blocks it sampled, turning a 5% sample into 5%
  of the transfer rather than 100%.

A sidecar is preferable to a footer: it needs no format-version bump, old
readers ignore it, and it can be generated lazily and cached. A footer would
have to be skipped by readers that predate it, and today's readers parse block
headers until EOF — so trailing bytes would be misread as a block header.

## String Label Support (Re-introduce with proper design)

The fixstr (`Z`) dtype was removed in V2 to simplify the initial implementation. Future versions could re-add string support with one of these approaches:

### Option A: Quoted strings in headers
- Values containing whitespace must be double-quoted: `MD:Z:"50 A2"`
- Unquoted values are split on whitespace as usual
- Simple to parse, familiar convention

### Option B: Binary-safe fixed-width strings
- No header parsing changes — strings are only accepted via programmatic API
- Fixed-width, null-padded (same as original fixstr)
- Header parsing only supports numeric types; string labels are set programmatically

### Option C: Variable-length strings (`V` dtype)
- Length-prefixed encoding: 2-byte length + UTF-8 payload
- Per-record variable size means label columns become variable-width
- Requires a different block layout (offset table or length-prefixed runs)
- Most flexible but most complex

**Recommendation:** Start with Option A (quoted strings) when string labels become needed. It requires minimal format changes and is backward-compatible.

---

## Categorical / Dictionary-Encoded Labels

Tags like `tp:A` (alignment type) have a small alphabet (e.g., `P`, `S`, `I`, `*`). A dictionary-encoded column would:

- Store a dictionary of unique values per block
- Store per-record indices (1-2 bits for small alphabets)
- Compress 8x for 4-value categoricals vs raw uint8

This is especially valuable for labels with low cardinality repeated across millions of records.

### Design sketch
- New dtype code `E` (enum)
- Block-level dictionary written before column data
- Per-record values are indices into the dictionary
- Falls back to raw storage if cardinality exceeds threshold (e.g., 256)

---

## Per-Column Compression

Instead of compressing the entire block (flags + labels + lengths + seqs) as a single ZSTD frame, compress each column independently. Benefits:

- Different columns have different entropy profiles (flags are highly redundant, sequences are near-random)
- Per-column compression allows different strategies (RLE for flags, delta for sorted labels, raw ZSTD for sequences)
- Enables selective column reads (skip decompressing sequences when only reading labels)

### Approach
- Block header lists per-column compressed sizes
- Each column is an independent ZSTD frame (or uncompressed if below threshold)
- Reader can seek past columns it doesn't need

---

## Label Indexes

For query patterns like "give me all reads where `AS > 100`", a block-level min/max index per label column would enable block skipping:

- Store `(min_value, max_value)` per label per block in the block header
- Reader can skip entire blocks whose range doesn't intersect the query
- Zero-cost during write (just track running min/max per batch)

---

## Multi-File Paired Labels

The current `stream_inputs_labeled` only supports single-end FASTQ. Future work:

- Support paired files (`r1.fq.gz` + `r2.fq.gz`) with labels from both
- Handle the case where R1/R2 have different tag subsets (e.g., `ts:A` only in R1)
- Interleaved labeled FASTQ support

---

## Label Arithmetic / Derived Labels

Some useful labels are derivable from existing ones:

- `insert_size = merged_end - merged_start` (from fastp merged suffix)
- `mapping_quality = AS - s2` (alignment score gap)
- `is_mapped = tp != '*'`

A lightweight expression engine could compute derived labels during encoding without storing them on disk.

---

## Arrow / Parquet Export

For downstream analysis, exporting ZNA label columns to Apache Arrow or Parquet would enable:

- Zero-copy integration with pandas/polars
- SQL queries over label columns
- Integration with bioinformatics analysis pipelines

Since ZNA already stores labels in columnar format, the conversion is mostly a metadata mapping exercise.
