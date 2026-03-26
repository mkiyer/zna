# ZNA Future Extensions

*Ideas for future development beyond the current performance optimization roadmap.*

---

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
