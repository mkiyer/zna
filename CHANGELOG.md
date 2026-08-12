# Changelog

All notable changes to the ZNA project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.4.0] - 2026-08-12

Adds `zna merge`, an overlap merger for paired-end reads, with a compiled
C++ kernel. No on-disk format change, and nothing existing changes
behaviour.

### Added
- **`zna merge` — overlap-merge paired-end reads into one mixed interleaved
  FASTQ**, ready for `zna encode --interleaved`. Replaces fastp's PE-merge step.

  Every pair is scored **once**: R1 is slid against `revcomp(R2)` over the single
  axis of candidate fragment lengths, and each shift is scored as a
  log-likelihood ratio in **bits** against chance alignment — `+1.9855` per
  matching base (`log2 4`, the information in agreeing on one of four bases),
  `-6.2288` per mismatch at a 1% error rate. Both weights fall out of the error
  rate; neither is tuned. The `argmax` shift (not fastp's first-accept) is then
  read at two thresholds: at or above `--threshold-merge` (28 bits) the pair
  becomes one full-fragment record; between that and `--threshold-trim` (8 bits)
  both reads are kept and the redundant overlap is cut off R2's **3'** end; below
  it both are kept untouched. Three parameters in total, all with units, and read
  length is not one of them.

  Measured against the fastp-derived rule it replaces: spurious detection on
  genuinely unrelated pairs **5.17% → 0.26%**, zero spurious merges in 20,000
  uniform-random pairs at 28 bits, merge rate **85.7% → 88.8%**.

  Where the mates overlap and disagree, the consensus takes the better-supported
  call by posterior from the two Phred scores and derates its quality, because a
  contested base is less certain than an uncontested one. On 167,784 real pairs
  that is **35.2% fewer expected residual errors** than fastp's Q14/Q30 cutoffs,
  entirely in the quality band those cutoffs never reach.

  Why it belongs here: base 0 of every emitted read is a true fragment boundary
  (nothing is ever cut from a 5' end), a merged record *is* its fragment exactly,
  and unmerged pairs are emitted all-or-nothing so a lone mate never becomes a
  spurious "full molecule with both endpoints". Those are precisely the
  properties `IS_RC`, `IS_FULL_FRAGMENT` and `--treat-unpaired-as-merged` rest
  on, and they were previously asserted across a repo boundary. See
  [docs/READ_MERGE_REDESIGN.md](docs/READ_MERGE_REDESIGN.md) for the derivation
  and [docs/MERGE_TOOL_AUDIT.md](docs/MERGE_TOOL_AUDIT.md) for what was measured
  and rejected.

  Also available in-process as `zna.merge.process_pair` / `zna.merge.find_overlap`,
  and as `python -m zna.merge`.

### Performance
- **A second compiled extension, `zna.merge._accel`**, carries the whole per-pair
  path: FASTQ parsing, the overlap scan, the posterior consensus, record
  construction, formatting and the histograms. It releases the GIL, so
  `--threads` are real worker threads.

  | | µs/pair |
  |---|---:|
  | `--threads 1` | 2.78 |
  | **`--threads 2`** | **1.40** |
  | `--threads 4` | 1.43 |

  The scan alone went from 2.633 µs/pair (the numba kernel this replaces) to
  0.470 in C++ — 5.6x — measured on 50,000 real pairs with full pruned scans.

  The tool is now I/O bound, and measurably so: with gzip removed from both ends
  it runs at **0.42 µs/pair**, while decompression costs 1.42 and compression
  0.52. Further speed needs a faster inflate or removing the FASTQ intermediate
  altogether, not a faster kernel.

- **The kernel compares raw bytes, 16 at a time**, through one vector primitive:
  NEON on aarch64, SSE2 on x86-64. Both are baseline, so there is no `-march`
  flag, no runtime ISA dispatch, and a scalar fallback everywhere else. Byte
  comparison is not a concession to portability — it *is* the reference
  semantics, for ACGT and equally for N, IUPAC codes and lowercase, so no fast
  path can disagree with the reference on any input. A 2-bit packed kernel
  measured slower (0.535 vs 0.470) *and* would have needed a purity dispatch to
  keep those semantics.

### Reproducibility
- **The score is computed in integers**, at 2²⁴ units per bit. No float takes
  part in a comparison, a bail bound or the argmax, so a given FASTQ produces
  byte-identical output on every platform and compiler. The scale is chosen by
  exhaustive enumeration rather than taste: 2²⁰ first disagrees with float
  arithmetic at an overlap of 2,575 bases, 2²⁴ at 32,830.
- **The argmax is a specified total order** — maximise score, then minimise the
  shift — rather than an artifact of iteration order, verified against an
  unpruned exhaustive scan including deliberately constructed ties.
- **Output is written in submission order**, so it is byte-identical at any
  `--threads`.

### Packaging
- **No new runtime dependencies.** The optional `numba` dependency that the
  merge tool used before it was ported here is gone; the reference kernel is
  plain Python and the fast one is compiled with the rest of the package.
- `zna merge` **refuses to run on the reference kernel** unless asked by name
  with `--backend python`. It is correct and ~50x slower, and a silently-correct
  50x slowdown at cluster scale is indistinguishable from a slow node.

## [0.3.5] - 2026-08-12

Two fixes found by reviewing the pipelines that write and read ZNA files
(hulkrna and khorana). No on-disk format change.

### Fixed
- **`zna encode --shuffle-buffer-size` did nothing.** The flag was declared and
  parsed, but `encode_command` passed a hardcoded `1 << 30` to `shuffle_zna`, so
  the bucket budget was always 1 GiB no matter what was asked for.
  (`zna shuffle`, the standalone command, honoured its `--buffer-size`
  correctly; only the encode path was affected.) The default still parses to
  1 GiB, so no existing invocation changes behaviour — the flag simply works now.

### Added
- **`ZnaReader.blocks(labels=...)`** — `blocks()` previously refused labeled
  files outright, on the reasoning that handing back sequences while silently
  dropping the label columns was worse than not offering the API. That reasoning
  still holds for the *default*, but it made the batch API unusable on any
  labeled corpus, which is what the pipeline it was built for actually produces.

  The default still raises, now naming the way out. `labels=False` skips the
  columns — not silent, the caller asked — and `labels=True` yields a third
  element holding one value-tuple per label column, in header order.
  `len(label_columns)` always equals `len(header.labels)`, so an unlabeled file
  yields `()` rather than erroring.

  Measured on a three-column file (int32, int32, uint8) of 200k x 150 bp:
  `records()` 32.7 ms, `blocks(labels=True)` 15.5 ms (2.12x),
  `blocks(labels=False)` 9.5 ms (3.44x).

  Note this does not avoid *decompressing* the label bytes — a block payload is
  a single zstd frame holding all four columns — it avoids unpacking them into
  Python objects, which is where the cost is.

### Packaging
- **CPython 3.14 wheels are now built.** `CIBW_BUILD` had listed `cp314-*` since
  0.3.3, but the workflow pinned `pypa/cibuildwheel@v2.23`, which predates 3.14
  entirely — the entry matched nothing and was silently dropped, so 0.3.3 and
  0.3.4 shipped no cp314 wheels despite asking for them. Verified: cibuildwheel
  2.23.3 selects cp310-cp313 for this project, 4.2.0 selects cp310-cp314.
- **The manylinux baseline is pinned to `manylinux2014`.** cibuildwheel 3.0
  changed the default image to `manylinux_2_28`, which would have raised the
  glibc floor from 2.17 to 2.28 and dropped RHEL/CentOS 7-era HPC clusters as a
  side effect of adding 3.14. `manylinux2014` still builds 3.14, and `zstandard`
  — zna's own runtime dependency — publishes its cp314 wheel as
  `manylinux_2_17`, so taking the new default would have made zna harder to
  install than the package it depends on.
- **Wheels are now smoke-tested before publication.** They were previously
  uploaded without ever being imported. `scripts/wheel_smoketest.py` runs as
  `CIBW_TEST_COMMAND` on every wheel: it asserts the compiled extension actually
  loaded (a wheel that silently falls back to the pure-Python backend is a
  broken wheel) and exercises `records()`, `blocks()` with and without labels,
  `block_index()`, `restore_strand`, and the N-policy error path.

### Fixed (latent)
- `blocks()` split the block payload as `flags | lengths | sequences`, which is
  only correct without labels. It never ran on a labeled file because of the
  guard above, but lifting that guard without fixing the offset would have
  decoded sequences from the wrong position — plausible garbage, no error. The
  split now accounts for the label columns whether or not their values are
  wanted, and the fuzz harness checks it against `records()` across dtypes,
  compression, sharding, `indices=`, and `restore_strand`.

## [0.3.4] - 2026-08-11

Performance release. Nothing about the on-disk format changes, and no existing
file needs re-encoding.

### Fixed
- **The pure-Python backend silently corrupted sequences containing an invalid
  base at an index congruent to 3 mod 4.** `_pycodec.encode_sequence` validated
  each 4-base group with `val > 255`, applied *after* OR-ing the four 2-bit
  values together. A 255 (invalid) in the group's last slot ORs to exactly 255
  and never exceeds it, so instead of raising, the group was packed as `0xFF`
  and decoded back as `TTTT` — losing the three valid bases beside it. This was
  reachable on real data: the IUPAC ambiguity codes (`R`, `Y`, `S`, `W`, ...)
  occur in FASTQ, and with an N-policy set they hit this path rather than being
  substituted. The C++ backend was never affected. Validation now happens once
  up front in a single scan, which is also faster than the per-group test.

  Found by the new round-trip fuzz harness on its first run.

### Added
- **`tests/test_fuzz_roundtrip.py`** — a round-trip fuzz harness asserting
  bit-exact recovery across sequence content x length (including 0, 1, and every
  residue mod 4) x `seq_len_bytes` x N-policy x strand configuration x
  compression x record layout x backend, against a reference model written from
  the specification rather than from `_pycodec`. It also checks that the two
  backends produce byte-identical columnar streams, that re-encoding is a
  faithful copy, that every label dtype round-trips, and that the C++ decoder
  neither leaks nor over-frees.

  This exists because the suite could not previously catch codec corruption:
  it round-trips fixed inputs and checks API contracts, and two changes proposed
  during the 0.3.4 sweep silently corrupted data with all 282 tests passing.

- **`ZnaReader.blocks(stride=1, offset=0, restore_strand=False)`** — a columnar
  batch API yielding `(list_of_sequences, flags_bytes)` per block, for consumers
  that process a batch at a time. `stride`/`offset` shard by block and seek past
  blocks the shard does not want, so a sharded reader decodes its share instead
  of the whole file: 1.8x for one worker of 2, 9.4x for one worker of 16.
  Requires arbitrary record order (use `zna shuffle`) and many more blocks than
  shards; warns when a shard matches no blocks. Raises on labeled files.
- **`ZnaReader.block_index()`** — one `BlockInfo(index, offset, n_records,
  comp_size, uncomp_size)` per block, built by walking the block-header chain
  and seeking over each payload. The *file* header stores no record or block
  count; every *block* header carries its own count, so the totals cost no
  decompression: 2.3 us per block, about 1.4 ms for a 38 MB / 611-block file,
  against 89 ms to reach the same counts by decoding. Use it to size a
  subsample, or to catalogue a corpus.
- **`blocks(indices=...)`** — decode an explicit set of block numbers instead of
  a stride. Needed when the sampling fraction is not a unit fraction, and when
  repeated passes over one file should see different blocks: `stride` admits
  only `stride` distinct phases, so several epochs at `stride=4` would revisit
  the same four subsets. Mutually exclusive with `stride`/`offset`.
- **`zna inspect --json`** (with `--blocks` and `--counts`) — machine-readable
  output including `n_blocks` and `n_records`, for building manifests that carry
  record counts and can therefore weight a balanced corpus sample without
  opening the data files.
- **`FLAG_FIELDS`** — `flag byte -> (is_paired, is_read1, is_read2)`, so
  `blocks()` consumers need not re-derive the bit layout.
- **`ENDS_BY_FLAG`** — `flag byte -> (has_start, has_end)`, now public. This is
  the only correct way to recover fragment geometry from the raw flags column:
  under unstranded normalization ZNA reverse-complements *one mate per pair, at
  random*, so which edge is a true fragment boundary is a per-record fact
  carried by `IS_RC`, not a property of R1 versus R2. A consumer that assigns
  endpoints by mate number instead is right about half the time.
- `decode_block_sequences` on both codec backends.
- `ZNA_NO_EXTERNAL_GZIP=1` forces gzip input through the Python `gzip` module.

### Changed
- **The C++ codec no longer copies every base twice around the work that
  matters.** nanobind was materialising a `std::vector<std::string>` on the way
  in and a `std::vector<std::tuple<...>>` on the way out. Decode now writes bases
  straight into each `str`'s own buffer, four at a time from a lookup table, and
  returns a hand-built list of GC-untracked tuples; encode reads each sequence's
  buffer in place and folds reverse-complement into the packing loop, so a
  normalized record needs no intermediate string. The decoder also emits the
  caller's tuple width and can undo strand normalization itself, which removes
  the Python-side rebuild of every record. `encode_block` and
  `encode_block_labeled` were hand-copied duplicates and now share one core.

  Verified byte-identical to 0.3.3 across 90 encode/decode configurations,
  including `npolicy="random"` and every reverse-complement path.

- Gzip FASTQ input is decompressed by `pigz`/`gzip` when available, so it runs
  in its own process rather than contending for the GIL. Falls back to the
  `gzip` module. A truncated or corrupt `.gz` is now an error rather than a
  silently short input — previously the partial stream was indistinguishable
  from a small file, and `zna encode` would write a truncated ZNA and report
  success.
- `shuffle`'s counting pass reads flag columns instead of decoding every
  sequence, and its bucket files are written at zstd level 1 rather than the
  source file's level.
- `CMakeLists.txt` no longer passes `STABLE_ABI`. It was silently a no-op —
  nanobind honours it only when CMake finds `Development.SABIModule`, and the
  `find_package` asks for `Development.Module` — and the module could never have
  been limited-API anyway.

### Performance

Measured interleaved, min-of-N, in-process (see `scripts/bench_ab.py`); the C++
figures are cross-process, since the two arms are different builds.

| path | change |
|---|---|
| `shuffle` counting pass | +2899% |
| raw `encode_block`, reverse-complementing | +163% |
| gzip FASTQ input | +123% |
| raw `encode_block` | +113% |
| raw `decode_block` | +109% |
| `records(restore_strand=True)` | +147% |
| `records()` | +97% |
| one sharded worker of 16, via `blocks()` | +837% |
| single-end encode loop | +52% |
| gzip module read-ahead | +51% |
| FASTQ sequence-line trimming | +31% |
| `reverse_complement` | +31% |
| interleaved pairing, bytes-native | +29% |
| N-drop filter | +26% |
| decode output f-strings | +19% |
| `_flags_from_ends` in a write loop | +11% |

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
