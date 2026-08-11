# Changelog

All notable changes to the ZNA project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.3.3] - 2026-08-11

### Fixed
- **Labeled encoding no longer discards pairing.** `stream_inputs_labeled` had
  its own single-file FASTQ reader and flagged every record unpaired, so
  `zna encode --interleaved --label-defs …` silently ignored `--interleaved` and
  wrote `is_paired=False` on every record — losing the R1/R2 flags and, because
  the encoder then drew an independent reverse-complement coin per record
  instead of one per pair, mate co-orientation as well. Labeled encoding now
  dispatches over the same input strategies as unlabeled (two paired files,
  interleaved, single-end) and shares one implementation of the interleaved
  pairing rule. A consumer cannot otherwise distinguish a mate (one true
  fragment boundary) from a merged read (two), which is the difference between
  correct fragment-end supervision and endpoint tokens on interior positions.
  **`is_paired` is not recoverable from affected files; re-encode them from
  FASTQ.** Labeled FASTA input is now rejected explicitly rather than parsed as
  FASTQ — labels come from SAM tags in the read header.
- **Copying records between ZNA files no longer re-applies strand
  normalization.** `strand_normalized` in a header describes what was already
  done to the data, but the writer read it as an instruction to do it, so any
  read-then-write path normalized a second time. Orientation is not idempotent:
  applying it twice returns the data to an un-normalized state while the header
  still claims otherwise. Affected `zna encode` on a `.zna` input, `zna shuffle`,
  `zna encode --shuffle`, and `zna._shuffle.shuffle_zna`. Encoding from FASTQ was
  never affected, and no on-disk format or existing file's bytes change.
- **`zna shuffle` now preserves each record's `IS_RC` flag.** It writes through
  two `ZnaWriter`s (buckets, then output), so it applied two further random
  reverse-complements. Mate co-orientation survived by parity, leaving the
  sequences valid and the damage invisible, but `IS_RC` ended up uncorrelated
  with which mate had actually been flipped.
- **Stranded re-encodes are stable.** Under `--strand-specific` the second pass
  reverse-complemented read1 (and merged/single reads) back to their original
  orientation, silently un-normalizing a strand-specific file.
- `zna encode --label ... input.zna` no longer exits 0 having written an empty
  file. A ZNA input routed to the labeled encode path was parsed as FASTQ,
  matched no records, and produced a valid, correctly-labeled, empty output.
- **Affected files, and what to do with them.**
  - A file re-encoded an **odd** number of times has lost mate co-orientation.
    Its stored `IS_RC` records exactly the last pass's reverse-complements, so
    `zna decode --restore-strand` undoes that pass and hands back the correctly
    co-oriented normalized sequences — byte-identical to the file before the
    damage, verified on both backends. What is *not* recoverable from the file
    is the `IS_RC` provenance and the original as-sequenced orientation, since
    the first pass's flags were overwritten.
  - A shuffled file is still correctly co-oriented and remains usable exactly as
    before, but its `IS_RC` flags no longer identify which mate was
    reverse-complemented, so it should not be read with `with_rc` below. Two
    passes leave each pair with either both mates flipped or neither, so
    re-normalizing it cannot recover the geometry.
  - In both cases, re-derive from the original FASTQ if the fragment-boundary
    geometry is needed.

- **`--npolicy drop` is now pair-atomic.** It filtered records *after* pairing,
  so an N in one mate dropped that mate and left its partner behind as a lone
  `IS_PAIRED` record — indistinguishable downstream from a genuine single-end
  read, and the only route to the encoder's mate-detection window mis-pairing
  two unrelated R1s. The whole fragment is now dropped together.

### Added
- **`IS_FULL_FRAGMENT` record flag (bit `0x10`)** — the record spans its entire
  fragment, so **both** of its edges are true fragment boundaries. `IS_RC` can
  only name one edge, which under-marks merged reads and over-marks genuine
  single-end reads, and nothing in the format previously distinguished them.
  Full-overlap **pairs are detected automatically** (mates spanning the same
  interval are exact reverse complements); unpaired records are declared with
  `zna encode --treat-unpaired-as-merged`, which defaults off — the conservative
  reading. Bits 4–7 of the record flag byte were unused and existing readers
  extract specific bits, so this needs **no format version bump**: old readers
  read new files unchanged, and new readers see the bit clear on old files.
- `ZnaReader.records(with_ends=True)` yields `has_start, has_end` — whether the
  left and right edge of the stored sequence is a true fragment boundary. This
  is the form a consumer placing fragment-end supervision wants: it combines
  `IS_RC` with `IS_FULL_FRAGMENT` so the geometry never has to be re-derived.
  The pair round-trips losslessly back to both flags, so it is also the shape
  for a faithful copy — `write_records` accepts it directly under
  `preserve_normalization=True`.
- `zna inspect --counts` reports a full-fragment tally and a **(mate × IS_RC)
  cross-tabulation**, with the expected pattern for the file's strandedness and
  a warning when R1 and R2 counts differ (orphaned mates). A bare RC total is
  equally consistent with a healthy file and a broken one; this is the tally
  that actually verifies a file's geometry before anything trains on it.
- `zna decode` warns when the input is strand-normalized and `--restore-strand`
  was not given: the output is then in normalized orientation, and re-encoding
  it with `--strand-normalize` applies orientation a second time — silently,
  since both files' headers are identical.
- `ZnaReader.records(with_rc=True)` yields the per-record `IS_RC` flag:
  `(seq, is_paired, is_read1, is_read2, is_rc)`, or
  `(seq, is_paired, is_read1, is_read2, is_rc, labels)` for labeled files. This
  is the only way to recover which edge of a mate is a real fragment boundary —
  it cannot be derived from the sequence, and `restore_strand=True` is not a
  substitute because it *consumes* the flag to undo the reverse-complement. The
  two are mutually exclusive and raise `ValueError` together. Opt-in: the
  default tuple widths (4, and 5 labeled) are unchanged and now pinned by test.
- `ZnaWriter(..., preserve_normalization=True)` with an optional trailing
  `is_rc` on `write_record`, plus `write_records` acceptance of the 5-tuple
  `records(with_rc=True)` yields, making a lossless ZNA → ZNA copy expressible
  for the first time. `write_zna` takes the same flag. Supplying `is_rc` without
  it raises `ValueError`.

### Changed
- `zna encode` refuses, rather than silently mangling, two re-encode requests it
  cannot honour: changing `--strand-specific` on an already-normalized input
  (orientation cannot be re-derived after the fact), and re-encoding a labeled
  ZNA file (label columns are not carried through that path).

### Documentation
- README gains "Unstranded Normalization and Fragment Geometry", stating the
  invariant the `IS_RC` flag exists to record: whichever mate was
  reverse-complemented ends up at the right of the common frame, so its right
  edge is the real fragment boundary and its left edge is a read-length cutoff;
  for the other mate it is the mirror image. The record format spec now lists
  `IS_RC` (bit 3), which it had omitted.

### Internal
- 20 tests covering the fragment-boundary contract, re-encode and shuffle
  record-identity, the default tuple widths, the stranded mirror case, and
  cross-backend `IS_RC` equality (both decoders, labeled and unlabeled). The
  pass-through encode path is deterministic, so it is additionally asserted
  byte-identical between the Python and C++ backends — unlike the deriving
  unstranded path, whose coin comes from a different PRNG in each.

## [0.3.1] - 2026-07-13

### Fixed
- Strand normalization now treats single-end / merged reads as **read1** under
  `--strand-specific`. Previously such reads were never reverse-complemented, leaving
  merged reads (e.g. from fastp's mixed interleaved output) on the opposite strand from
  the normalized paired read1 records. Decoding and already-encoded files are unaffected;
  only new strand-specific encodes of single reads change.

### Added
- `zna inspect --counts` reports per-flag record counts (paired R1, paired R2,
  single/merged, and reverse-complemented). Useful for validating mixed
  paired-end + single-end interleaved streams. The default header-only scan is
  unchanged; `--counts` reads block payloads (decoding only the flags column).
- Encode-time warning when `--strand-normalize` is combined with
  `--strand-specific` but neither read1 nor read2 is antisense — a configuration
  in which no read would ever be reverse-complemented (suppressed by `-q`).

### Internal
- Added a cross-backend lockstep test suite asserting the Python (`_pycodec`)
  and C++ (`_accel`) encoders produce byte-identical output across a battery of
  flag / strand / N-policy / length combinations, guarding against drift in the
  duplicated encode logic.

## [0.3.0] - 2026-03-25

### Added
- **Per-sequence label columns** for storing SAM-tag metadata (alignment scores,
  edit distances, etc.) alongside compressed sequences
  - Labels defined via `--label NAME:TYPE` on the CLI (repeatable)
  - 11 numeric dtypes: `A`, `c`/`C`, `s`/`S`, `i`/`I`, `f`, `d`, `q`/`Q`
  - Per-label `missing` value for reads where a tag is absent
  - Columnar storage: each label column packed contiguously per block
- **Decoupled label name and tag** — label `name` (stored in ZNA, up to 16
  chars) can differ from the SAM `tag` parsed at encode time
  - CLI 3-part format: `--label edit_dist:C:NM` (name:type:tag)
  - YAML: optional `tag` field per label definition
  - Tag is encode-time only; not stored in the ZNA file
  - If only name or tag is given, the other defaults to the same value
- **YAML label specification files** (`--label-defs labels.yaml`)
  - Define labels, descriptions, dtypes, and missing values in a single file
  - CLI `--label` / `--label-desc` flags override YAML definitions
  - Sample template: `examples/labels.yaml`
- **C++ accelerated label encode/decode** — labeled files now use the C++
  backend for both encoding and decoding, eliminating the Python fallback
  - `encode_block_labeled`: encodes sequences + pre-packed label columns in C++
  - `decode_block_labeled`: decodes sequences and unpacks label columns by dtype
  - `extract_labels_fast`: zero-allocation C++ SAM-tag extraction from headers
- **Optimized Python header parsing** (fallback path)
  - Eliminated `bytes.decode()` calls — `int()`/`float()` accept bytes directly
  - Replaced lambda dispatch with integer conv-codes
  - Early-exit when all labels found

### Performance
- **Labeled encode**: 336K rec/s (was 80K rec/s) — **4.2× faster**
- **Labeled decode**: 2.23M rec/s (was 610K rec/s) — **3.7× faster**
- Label overhead reduced from ~39% to within typical I/O variance

### Changed
- Binary format: label definition on-disk size is now 89 bytes
  (16 name + 64 desc + 1 dtype_code + 8 missing)
- `LabelDef.param` replaced by `LabelDef.missing` (explicit missing sentinel)
- Header field splitting uses any whitespace (tabs and spaces) so fastp
  `merged_XX_YY` suffixes no longer corrupt trailing tag values
- Label tag parsing now supports variable-length keys in `KEY:TYPE:VALUE`
  header fields (not limited to 2-character SAM tags)
  - Python fallback parser and C++ `extract_labels_fast` use dynamic key parsing
  - Supports mixed key lengths in the same header (e.g. `NM` + `edit_distance`)
- Dropped fixstr (`Z`) label type — all labels are numeric

### Documentation
- Expanded README strand guidance to explicitly document
  `--strand-normalize` behavior and its interaction with `--strand-specific`
- Updated label docs and `examples/labels.yaml` to show `name` + `tag`
  decoupling and custom non-SAM, variable-length header keys

## [0.2.0] - 2026-03-24

### Added
- **Strand normalization** (`--strand-normalize`) for consistent strand orientation in ZNA files
  - Unstranded libraries: randomly reverse-complements one read per pair (or SE reads) with per-record `IS_RC` flag
  - Strand-specific libraries: deterministically RCs antisense reads to sense orientation
  - Orthogonal to library metadata (`--strand-specific`, `--read1-antisense`, etc.)
  - `--restore-strand` on decode reverses RC operations using per-record `IS_RC` flags
- New `IS_RC` record flag (bit 3) tracks which records were reverse-complemented
- New `STRAND_NORMALIZED` header flag (bit 3) indicates strand normalization was applied

### Fixed
- C++ accelerator bug: paired R2 reads could receive an unintended second random RC during unstranded normalization

### Changed
- Bioconda recipe: added `osx-arm64` and `linux-aarch64` to `additional-platforms`

## [0.1.8] - 2026-02-11

### Fixed
- macOS conda build: removed `llvm-tools` from host dependencies

## [0.1.7] - 2026-02-06

- **Documentation changes**

## [0.1.6] - 2026-02-06

### Added
- **Built-in shuffle command** (`zna shuffle`) for random shuffling of ZNA files with bounded memory
  - Bucket-shuffle algorithm with configurable memory budget (default 1GB)
  - Preserves paired-end read associations (R1+R2 stay together)
  - Deterministic shuffling with `--seed` option for reproducibility
  - Suitable for ML training data preparation
- **Encode with shuffle** (`zna encode --shuffle`) to shuffle during encoding
  - Eliminates need for separate shuffle step
  - Uses same memory-bounded algorithm as standalone shuffle command
- New `_shuffle.py` module with clean `shuffle_zna()` API for programmatic use

### Changed
- Refactored CLI: shuffle logic extracted from `cli.py` to `_shuffle.py` for maintainability
- Consolidated `--block-size` argument (now shared across encode and shuffle commands)

### Performance
- **Decode optimization**: Added `yield from` fast path when `restore_strand=False` (10% faster)
- **Encode optimization**: New `write_records()` batch method caches attribute lookups (20% faster)
- Overall throughput remains: 165 MB/s encode, 241 MB/s decode with 9.10x compression

### Documentation
- Updated README.md with shuffle command documentation and examples
- Updated PERFORMANCE.md with recent optimization notes
- All 86 tests passing

## [0.1.5] - 2026-02-04

### Changed
- Columnar block storage format for 2.5x better compression (9.10x vs 3.68x)
- Optimized flag operations (removed enum overhead for 60-135% speedup)
- Default block size: 128 KB (512 KB recommended for archival)

### Performance
- Compression: 9.10x (9.44x with 512 KB blocks)
- Speed: 165 MB/s encode, 241 MB/s decode (pure Python)
- File size competitive with BAM format

## [0.2.0] - 2026-02-03

### Added
- C++ acceleration with nanobind (9.5x speedup over pure Python)
- Block-based architecture for streaming I/O
- Strand-specific library support (dUTP, TruSeq protocols)
- N-policy handling (drop, random, replace)

### Performance
- Short reads: 189.5 MB/s encode, 668.8 MB/s decode
- Long reads: 2,824 MB/s encode, 3,393 MB/s decode

## [0.1.0] - 2026-01-15

### Added
- Initial ZNA format implementation
- 2-bit nucleotide encoding
- Zstd compression support
- CLI tools: encode, decode, inspect
- Basic FASTQ/FASTA support
- Paired-end and interleaved read handling

[0.1.6]: https://github.com/mkiyer/zna/compare/v0.1.5...v0.1.6
[0.1.5]: https://github.com/mkiyer/zna/releases/tag/v0.1.5
