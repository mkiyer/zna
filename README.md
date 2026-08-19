# ZNA: Compressed Nucleic Acid Format

**ZNA** (Compressed **Z**-Nucleic **N**-Acid **A**) is a high-performance binary format for storing DNA/RNA sequences with exceptional compression and I/O speed.

## Performance

- **1.5 GB/s decode** via `records()`, **3.6 GB/s** via `blocks()`, on 150 bp
  reads (0.5.1, in-memory, single core, Apple Silicon)
- **455 MB/s encode** on 150 bp reads — the delta from 0.3.x's 726 MB/s is the
  price of the self-describing container: per-block stats for the trailer and a
  content checksum in every frame
- **~4x compression** from 2-bit packing alone, more with Zstd on duplicated data
- **C++ acceleration** with a pure Python fallback

## Features

- **High Compression**: 2-bit encoding (4 bases per byte) + optional Zstd compression
- **Ultra-Fast I/O**: C++ accelerated encode/decode with block-based architecture
- **Minimal Dependencies**: `zstandard` and `pyyaml` only (C++ extensions ship prebuilt)
- **Flexible**: Single-end, paired-end, and interleaved reads
- **Block-Parallel**: fragments never split across blocks, so any subset of blocks decodes independently — shard a file by block without splitting a pair
- **Self-Describing**: every file carries its provenance (writer version, shuffled?) up front and its exact contents (counts, length histograms, block index, checksums) in a trailer — `zna inspect --verify` certifies a file with no sidecar
- **In-Process Merging**: `zna encode --merge-pairs R1.fq R2.fq` merges and encodes in one step, no FASTQ intermediate
- **Overlap Merging**: `zna merge` collapses overlapping pairs into full-fragment reads on one calibrated likelihood-ratio score, with a compiled kernel and byte-identical output on any platform
- **Strand-Specific Support**: dUTP, TruSeq, and custom strand protocols
- **Built-in Shuffle**: Memory-bounded random shuffling for training data preparation
- **Metadata Rich**: Read groups, descriptions, and custom flags
- **Unix-Friendly**: Pipe-compatible CLI for seamless workflow integration
- **Streaming**: Memory-efficient block-based processing

## Installation

```bash
pip install zna
# or, if your inputs are gzipped -- see below
pip install 'zna[fast]'
# or
conda install -c bioconda zna
```

Both ship prebuilt binaries with the C++ extensions already compiled — Linux, macOS
(Intel and Apple Silicon) and Windows, on CPython 3.10–3.14. Nothing needs a compiler.

**`zna[fast]` adds one optional dependency, `isal` (Intel ISA-L), for gzip input.**
Inflate is the largest single cost of `zna merge` — a profile puts a quarter of all
cycles inside `pigz`'s `libz` — and ISA-L reads a gzip member at **448 MB/s against
`pigz`'s 193 and the stdlib's 208**, which on a 1M-pair library takes `zna merge
--threads 1` from 5.84 s to 4.04 s while *also* using less CPU. Output is byte-identical
either way, and a truncated or CRC-damaged member is still rejected. Without the extra,
ZNA behaves exactly as it did before: `pigz` when it is on `PATH`, then stdlib `gzip`.

**Verify both extensions loaded.** ZNA has two: the codec, and `zna merge`'s overlap
scan. Either can be absent without an error — the pure-Python fallbacks are correct and
about 50x slower, so a broken install looks like a slow machine rather than a mistake.
Check it rather than assume it:

```bash
python -c "
import zna
from zna.merge.backend import available_merge_backends
from zna._gzip import inflate_backend_name
print('zna', zna.__version__)
print('codec accelerated:', zna.is_accelerated())
print('merge backends:   ', available_merge_backends())   # want 'accel' in here
print('gzip input via:   ', inflate_backend_name('x.fq.gz'))   # 'isal' with zna[fast]
"
```

From source, for development (needs a C++17 compiler and CMake ≥ 3.15):

```bash
git clone https://github.com/mkiyer/zna.git
cd zna
pip install -e .
```

**Requirements:**
- Python ≥3.10
- C++ compiler (for optimal performance)
- CMake ≥3.15 (auto-installed via pip)

## Quick Start

```bash
# Encode FASTQ to compressed ZNA (default: Zstd level 9)
zna encode sample.fastq.gz -o sample.zna

# Encode with shuffle (for ML training data)
zna encode sample.fastq.gz --shuffle -o shuffled.zna

# Encode with shuffle and explicit memory cap per bucket
zna encode sample.fastq.gz --shuffle --shuffle-buffer-size 512M -o shuffled.zna

# Shuffle an existing ZNA file
zna shuffle input.zna -o shuffled.zna

# Decode back to FASTA
zna decode sample.zna -o sample.fasta

# Inspect file statistics
zna inspect sample.zna

# Merge overlapping pairs and encode, one step, no intermediate
zna encode --merge-pairs R1.fq.gz R2.fq.gz -o sample.zna

# ...or the two-step path (still supported)
zna merge --in1 R1.fq.gz --in2 R2.fq.gz --out merged.fq.gz
zna encode --interleaved --treat-unpaired-as-merged merged.fq.gz -o sample.zna

# Certify a file: integrity, structure, stats vs its trailer; exit 0 iff good
zna inspect --verify sample.zna

# Pipe-friendly workflows
cat reads.fastq | zna encode -o reads.zna
zna decode reads.zna | head -n 1000
```

## Performance Benchmarks

### Throughput by Read Length

Measured on ZNA 0.3.5, Apple Silicon, single core, in-memory (`BytesIO`),
min-of-7. Throughput is sequence bases in/out per second; compression is bases
per stored byte at Zstd level 9 on random (worst-case, incompressible) sequence.

> **Re-measured at 0.5.1** (same method, 150 bp): encode 455 MB/s, decode
> 1,486 MB/s via `records()` and 3,621 MB/s via `blocks()`. The encode delta
> is the self-describing container (trailer stats + frame checksums); decode
> is within noise of the table. The *scaling with read length* below still
> holds; the absolute numbers are 0.3.5-era.
>
> **Every figure on this page is Apple silicon.** The same method at 0.5.2 on a Xeon
> E5-2680 v3 (Haswell, 2.5 GHz, 2014) gives encode 179 MB/s, decode 452 MB/s via
> `records()` and 811 MB/s via `blocks()`, at the identical 3.95x compression — roughly
> 2.5x lower across the board, which is the per-core difference between the two parts
> rather than anything ZNA does differently. **The codec was profiled on x86 in 0.5.2 and
> left alone:** its two kernels measure 5,465 Mbase/s (decode, already memory-bandwidth
> bound) and 1,723 Mbase/s (encode), and hand-written SSE2 was *slower* than both. What
> `zna encode` is actually limited by is its per-record Python driver — the compiled codec
> is ~3% of the command — see `docs/ROADMAP.md`.

| Read Type | Encode (MB/s) | Decode (MB/s) | Encode (rec/s) | Decode (rec/s) | Compression |
|---|---|---|---|---|---|
| Short (Illumina, 150 bp) | 726 | 1,718 | 4.8 M | 11.5 M | 3.95x |
| Medium (300 bp) | 1,073 | 2,443 | 3.6 M | 8.1 M | 4.00x |
| Long (PacBio, 1 kb) | 1,775 | 2,924 | 1.8 M | 2.9 M | 4.00x |
| Very Long (5 kb) | 2,610 | 5,770 | 0.5 M | 1.2 M | 4.00x |
| Ultra Long (15 kb) | **2,701** | **6,621** | 0.2 M | 0.4 M | 4.00x |

**Key Insights:**
- Throughput scales with read length: per-record overhead dominates at 150 bp
  and vanishes by 5 kb.
- 4x is the 2-bit packing floor. Real libraries with duplicate reads compress
  further; unique reads do not, because packed DNA is near-incompressible.
- `blocks()` is faster still for batch consumers — see
  [Batch Reading](#batch-reading-with-blocks).

See [docs/PERFORMANCE.md](docs/PERFORMANCE.md) for detailed benchmarking.

---

## Documentation

This README is the **user manual** — installation, usage, the file format, the command
reference and the Python API are all below.

| | |
|---|---|
| [CHANGELOG.md](CHANGELOG.md) | what changed in each release, and why |
| [docs/METHODS.md](docs/METHODS.md) | the algorithms: the overlap score and its two thresholds, the quality-aware consensus, fragment geometry and what the flags mean, the codec |
| [docs/MERGE_BENCHMARK_RESULTS.md](docs/MERGE_BENCHMARK_RESULTS.md) | `zna merge` scored against known ground truth and head to head with fastp. Read this before changing a threshold |
| [docs/PERFORMANCE.md](docs/PERFORMANCE.md) | compression ratios, throughput, and tuning |
| [docs/ROADMAP.md](docs/ROADMAP.md) | what is scheduled, what is being considered, and what was tried and closed by measurement |
| [docs/RELEASING.md](docs/RELEASING.md) | publishing to PyPI and Bioconda *(maintainers)* |
| [docs/TRAILER_PLAN.md](docs/TRAILER_PLAN.md) | the self-describing container: prologue, trailer, footer, `inspect --verify` — executed in 0.5.0; §14 records the amendments |
| [docs/archive/MERGE_PAIRS_PLAN.md](docs/archive/MERGE_PAIRS_PLAN.md) | `zna encode --merge-pairs` — executed in 0.5.0 |
| [docs/archive/NPOLICY_PLAN.md](docs/archive/NPOLICY_PLAN.md) | the `--npolicy` design and what remains of it |
| [docs/HANDOFF_0.4.0.md](docs/HANDOFF_0.4.0.md) | a record of the 0.4.0 release; its build traps and ground-truth notes are still current *(maintainers)* |

---

## File Format Specification

### Overview

ZNA files use a binary format optimized for nucleic acid sequences:

- **File Extension**: `.zna` (for both compressed and uncompressed files)
- **Default Compression**: Zstd level 9 (use `--uncompressed` to disable)
- **Magic Number**: `ZNA\x1A` (4 bytes)
- **Version**: 3 (1 byte)
- **2-bit Encoding**: A=00, C=01, G=10, T=11
- **Block Structure**: Columnar blocks, compressed as one Zstd frame each
- **Fragment-complete blocks**: a fragment's reads are consecutive and never split
- **Self-describing (0.5.0)**: a provenance prologue after the header, a stats trailer + footer at EOF; zstd frames carry content checksums
- **Metadata**: Read groups, descriptions, and custom information

### File Structure

```
┌─────────────────────────────────────────┐
│  File Header (15 bytes fixed)           │
│   - Magic "ZNA\x1A" (4 bytes)           │
│   - Version = 3 (1 byte)                │
│   - Sequence length width (1 byte)      │
│   - Flags (1 byte)                      │
│   - Compression method (1 byte)         │
│   - Compression level (1 byte)          │
│   - Label count (2 bytes)               │
│   - Read-group / description lens (4 B) │
│  + read group, description (variable)   │
│  + one 89-byte definition per label     │
├─────────────────────────────────────────┤
│  Block 0                                │
│   Block Header (20 bytes)               │
│    * Compressed size    (4 bytes)       │
│    * Uncompressed size  (4 bytes)       │
│    * Record count       (4 bytes)       │
│    * Flags column size  (4 bytes)       │
│    * Lengths column size(4 bytes)       │
│   Payload — ONE Zstd frame, COLUMNAR:   │
│    ┌───────────────────────────────┐    │
│    │ flags    (1 byte  per record) │    │
│    │ labels   (per schema, if any) │    │
│    │ lengths  (1/2/4 B per record) │    │
│    │ sequences (2-bit packed)      │    │
│    └───────────────────────────────┘    │
├─────────────────────────────────────────┤
│  Block 1 ...                            │
├─────────────────────────────────────────┤
│  Sentinel (20-byte header, count = 0)   │
│  Trailer — canonical JSON, derived      │
│   facts only: record/base/pair counts,  │
│   length histograms, per-flag counts,   │
│   the stored block index, prologue CRC  │
├─────────────────────────────────────────┤
│  Footer (32 bytes at EOF)               │
│   magic | ver | crc32(trailer) |        │
│   sentinel offset — O(1) discovery      │
└─────────────────────────────────────────┘
```

A **provenance prologue** — another count-0 pseudo-block — sits between the
file header and Block 0, carrying the facts known at encode *start*:
`writer_version`, `shuffled`, `merged_in_process`. A streaming consumer, pipes
included, learns what wrote the file before decoding a record
(`ZnaReader.provenance`); the trailer (`ZnaReader.trailer`) carries only what
requires the *complete* encode. The trailer is the writer's signature that the
encode finished — an aborted encode writes none, so the file still reads but
`zna inspect --verify` refuses to certify it. Files written by zna < 0.5 have
neither and still read. The format version byte is unchanged: no data block
ever carries `count == 0`, so the pseudo-blocks are invisible to older
readers.

The payload is **columnar, not row-oriented**: all flags come first, then all
label values, then all lengths, then all packed sequence. That is what lets
`zna inspect` tally flags without touching sequence, and `blocks(labels=False)`
skip label columns.

**The file header stores no record or block count.** Each *block* header carries
its own count, so totals come from walking the block chain — see
[`block_index()`](#sizing-a-file-before-reading-it-block_index).

**A block holds whole fragments.** A fragment's reads are stored consecutively,
R1 immediately followed by R2, and never span a block boundary — so a block is a
self-contained set of molecules, and any subset of blocks decodes independently
of the rest of the file. That is what makes block sharding
([`blocks(stride=…)`](#batch-reading-with-blocks)) and every block-parallel
consumer sound: a worker never receives one mate of a pair whose other mate went
to a different worker. `ZnaWriter` enforces it on write, so it is a property of
every ZNA file rather than one that happens to hold. Blocks are flushed on an
estimated byte size, so a block overruns `--block-size` by at most one record.

### Record Format

A record's fields live in separate columns of the block, not adjacent to each
other. Per record:
- **Flags** (1 byte): IS_READ1 (bit 0), IS_READ2 (bit 1), IS_PAIRED (bit 2),
  IS_RC (bit 3 — set when strand normalization reverse-complemented this record),
  IS_FULL_FRAGMENT (bit 4 — the record spans its whole fragment, so *both* edges
  are true fragment boundaries). Bits 5-7 are reserved.
- **Length** (1-4 bytes): Sequence length (configurable)
- **Sequence** (variable): 2-bit encoded bases

### Compression

- **Method 0**: Uncompressed
- **Method 1**: Zstd, levels 1-22 (default 9)
- **Block Size**: Default 4 MiB (`--block-size`, accepts K/M/G suffixes)

Both use the `.zna` extension; compression is recorded in the header, not the
filename. Smaller blocks cost compression ratio only on duplicate-rich data —
on unique reads the packed sequence is incompressible, so block size is free.
A batch consumer holding a decoded block at a time (see `blocks()`) may want
`--block-size 1M` to bound its memory.

---

## Usage Guide

### Encoding

#### Single-End Reads

```bash
# From FASTQ file
zna encode sample.fastq -o sample.zna

# From FASTA file  
zna encode sample.fasta -o sample.zna

# From gzipped input
zna encode sample.fastq.gz -o sample.zna

# Lower level = faster encode, larger file (default is 9)
zna encode sample.fastq --level 5 -o sample.zna

# Uncompressed (rarely needed)
zna encode sample.fastq --uncompressed -o sample.zna

# From stdin (format inferred from content; extensions are used for files)
cat sample.fastq | zna encode -o sample.zna
```

#### Paired-End Reads

```bash
# Separate R1/R2 files
zna encode R1.fastq.gz R2.fastq.gz -o paired.zna

# Interleaved file (strict alternating R1/R2 pairs)
zna encode interleaved.fastq --interleaved -o paired.zna

# Interleaved from stdin
cat interleaved.fastq | zna encode --interleaved -o paired.zna
```

#### Mixed Paired-End and Single-End Reads (Interleaved)

The `--interleaved` mode intelligently detects both paired-end and single-end reads in the same file by analyzing read names. This is useful for output from tools like **fastp** that produce mixed merged (single) and unmerged (paired) reads.

**How it works:**
- Reads with matching base names (e.g., `read1/1` and `read1/2`) are paired
- Reads without matching pairs are treated as single-end
- Read names are used to determine pairing (not just alternating order)

```bash
# Mixed interleaved input (fastp output with merged + unmerged reads)
zna encode fastp_output.fastq --interleaved -o mixed.zna

# Example input structure:
#   @read1/1         →  paired with next read
#   @read1/2
#   @merged1         →  single-end (no pair)
#   @read2/1         →  paired with next read
#   @read2/2
#   @merged2         →  single-end (no pair)
```

**Read name formats supported:**
- `/1` and `/2` suffixes: `read1/1`, `read1/2`
- No suffix: treated as single-end unless next read has matching base name
- Comments ignored: `read1/1 merged_length:150` extracts `read1/1`

**Strand normalization of merged/single reads:** single-end reads (including merged
reads with no mate) are treated as **read1** for strand normalization. Under
`--strand-specific`, a single read is reverse-complemented exactly when read1 is
antisense, so merged reads end up on the same strand as normalized paired R1 reads.

#### Advanced Options

```bash
# Custom metadata
zna encode sample.fastq \
  --read-group "Sample_01" \
  --description "Experiment XYZ" \
  -o sample.zna

# Strand-specific library (default: R1 antisense, R2 sense)
zna encode R1.fastq.gz R2.fastq.gz \
  --strand-specific \
  -o stranded.zna

# Custom strand orientation (e.g., fr-secondstrand protocol)
zna encode R1.fastq.gz R2.fastq.gz \
  --strand-specific --read1-sense --read2-antisense \
  -o stranded.zna

# Handle sequences with N nucleotides
zna encode sample.fastq --npolicy trim3 -o clean.zna      # Cut each read at its first N (default)
zna encode sample.fastq --npolicy random -o clean.zna     # Substitute N from a seeded stream

# Shuffle during encoding (for ML training data preparation)
zna encode sample.fastq --shuffle -o shuffled.zna
zna encode R1.fastq.gz R2.fastq.gz --shuffle --seed 12345 -o shuffled.zna

# Control compression
zna encode sample.fastq \
  --level 9 \
  --block-size 262144 \
  -o sample.zna

# Uncompressed (rarely needed, for maximum I/O speed)
zna encode sample.fastq --uncompressed -o sample.zna

# Sequence length encoding (max sequence length)
zna encode sample.fastq \
  --seq-len-bytes 1 -o short_reads.zna    # max 255 bp

zna encode sample.fastq \
  --seq-len-bytes 2 -o sample.zna         # max 65,535 bp (the default)

zna encode sample.fastq \
  --seq-len-bytes 4 -o long_reads.zna     # max 4.2 billion bp
```

### Decoding

#### Basic Decoding

```bash
# To FASTA file
zna decode sample.zna -o output.fasta

# To gzipped FASTA
zna decode sample.zna -o output.fasta.gz

# To stdout (pipe-friendly)
zna decode sample.zna | head -n 1000

# From stdin
cat sample.zna | zna decode -o output.fasta
```

#### Paired-End Decoding

```bash
# Interleaved output (default)
zna decode paired.zna -o interleaved.fasta

# Split to R1/R2 files (use # placeholder)
zna decode paired.zna -o reads#.fasta
# Creates: reads_1.fasta and reads_2.fasta

# Split with gzip
zna decode paired.zna -o reads#.fasta.gz
# Creates: reads_1.fasta.gz and reads_2.fasta.gz

# Restore original strand for strand-specific libraries
zna decode stranded.zna --restore-strand -o reads.fasta
```

#### Piping Examples

```bash
# Extract first 1M reads
zna decode large.zna | head -n 2000000 > subset.fasta

# Count sequences
zna decode sample.zna | grep -c "^>"

# Convert to gzipped output via pipe
zna decode sample.zna --gzip > output.fasta.gz

# Chain operations
zna decode sample.zna | seqtk seq -r - | gzip > reversed.fasta.gz
```

#### Batch Reading with `blocks()`

`records()` yields one tuple per record. A consumer that works a whole batch at
a time — a training data loader, say — can instead take a block at a time and
skip the per-record tuple entirely:

```python
from zna import ZnaReader, FLAG_FIELDS

with open("sample.zna", "rb") as fh:
    for sequences, flags in ZnaReader(fh).blocks():
        # sequences: list[str];  flags: bytes, one per record, same order
        for seq, fl in zip(sequences, flags):
            is_paired, is_read1, is_read2 = FLAG_FIELDS[fl]
            ...
```

`ENDS_BY_FLAG[fl]` gives `(has_start, has_end)` from the same byte — whether each
edge of the stored sequence is a true fragment boundary. Use it rather than
inferring from the mate number: under unstranded normalization ZNA
reverse-complements one mate per pair *at random*, so the boundary edge is a
per-record fact, not a property of R1 versus R2.

**A block holds whole fragments.** Paired reads sit consecutively, R1 then R2,
and a fragment never straddles a block boundary — so a worker handed a block
never sees one mate of a pair whose other mate went to a different worker. This
is a guarantee of the format, enforced by `ZnaWriter` on every write path, not a
property that happens to hold for a given file. It is what makes block-parallel
consumers safe to write at all.

`stride`/`offset` shard **by block**, and — the point — seek past the blocks
this shard does not want instead of decoding and discarding them:

```python
# Worker 3 of 8: decodes ~1/8 of the file, not all of it.
for sequences, flags in ZnaReader(fh).blocks(stride=8, offset=3):
    ...
```

That is worth 1.8x at 2 workers and 9.4x at 16, compared with striding over
`records()`. Two conditions come with it:

- **Record order must already be arbitrary.** Shards get contiguous runs, not an
  interleave, so a file grouped by anything meaningful hands each worker a
  biased sample. Use `zna shuffle` first.
- **The file needs many more blocks than shards.** Shares are whole blocks, so
  a small file split many ways is lopsided, and past the block count some shards
  get nothing — which `blocks()` warns about rather than passing off as an empty
  file. The default 4 MiB block gives a few hundred blocks per GB; write with a
  smaller `block_size` if you need finer shards.

`blocks()` also takes `restore_strand=True`.

On a **labeled** file it needs an explicit `labels=`, because quietly discarding
label columns is not a decision it should make for you:

```python
reader.blocks()                 # labeled file -> raises
reader.blocks(labels=False)     # skip the columns  -> (sequences, flags)
reader.blocks(labels=True)      # -> (sequences, flags, label_columns)
```

`label_columns` holds one value-tuple per column in header order, each as long as
`sequences`; `len(label_columns)` always equals `len(header.labels)`, so an
unlabeled file yields `()`. On a three-column file, `labels=False` is 3.4x faster
than `records()` and `labels=True` 2.1x. (Both still inflate the label bytes — a
block is one zstd frame — what is saved is unpacking them into Python objects.)

#### Sizing a file before reading it: `block_index()`

**The ZNA file header stores no record or block count** — only the format
version, sequence-length width, strand flags, compression settings and label
schema. Each *block* header does carry its own record count, so the totals are
recovered by walking the block chain, seeking over each payload:

```python
reader = ZnaReader(fh)
index = reader.block_index()          # list[BlockInfo]
total = sum(b.n_records for b in index)
```

This decompresses nothing. Measured at 2.3 µs per block — 1.4 ms for a 38 MB,
611-block, 1M-record file, against 89 ms to reach the same counts by decoding.
Cheap enough to run at open time, or across a whole corpus to build a manifest.

That makes proportional subsampling straightforward: use the counts to decide
how much of each file you want, then decode only those blocks.

```python
import random

index = reader.block_index()
want = round(len(index) * target_fraction)
keep = random.sample([b.index for b in index], want)

for sequences, flags in reader.blocks(indices=keep):
    ...
```

`indices` is mutually exclusive with `stride`/`offset`. Prefer it when the
fraction is not a unit fraction, or when repeated passes over one file should see
*different* blocks — `stride` admits only `stride` distinct phases, so training
several epochs at `stride=4` would revisit the same four subsets.

Blocks are flushed on an estimated byte size, so record counts per block are
near-uniform for fixed-length reads and vary for variable-length ones. That is
why `block_index()` returns per-block counts rather than an average, and why
sampling *k* of *n* blocks gives approximately, not exactly, *k/n* of the records.

#### Cataloguing a corpus: `zna inspect --json`

```bash
zna inspect sample.zna --json
zna inspect sample.zna --json --blocks     # include the per-block array
zna inspect sample.zna --json --counts     # add per-flag record tallies
```

Emits header fields plus `n_blocks` and `n_records`, read from block headers
without decompressing. Fast enough to sweep thousands of files, so a manifest can
record record counts once and weight a balanced sample later without opening any
of them.

Batching alone (without sharding) is worth about 24% for a loader doing real
per-record work, and it fades with read length: ~24% at 150 bp, ~8% at 1 kb, and
nothing measurable at 10 kb, where the sequence dominates the record overhead.

### Inspecting Files

```bash
# Show file statistics
zna inspect sample.zna
```

**Example Output:**
```
File: sample.zna
Total Size: 45.32 MB

--- Header Metadata ---
Read Group:       Sample_01
Description:      Experiment XYZ
Seq Length:       2 bytes (Max: 65535 bp)
Strand Specific:  True
R1 Antisense:     True
R2 Antisense:     False
Strand Normalized:True
Compression:      ZSTD (Level 9)

--- Provenance ---
Writer:           zna 0.5.0
Shuffled:         True

--- Content Statistics ---
Total Blocks:       356
Total Records:      1000000
Total Bases:        142053188
Pairs / Unpaired:   417134 / 165732
Compressed Payload: 42.15 MB
Uncompressed Data:  125.50 MB
Compression Ratio:  2.98x
```

---

## Command Reference

### `zna encode`

Convert FASTQ/FASTA to ZNA format.

**Usage:**

```
zna encode [FILE1] [FILE2] [OPTIONS]

Positional Arguments:
  FILE1 [FILE2]          Input files (0=stdin, 1=single/interleaved, 2=paired R1 R2)

Options:
  --interleaved          Treat input as interleaved (auto-detects mixed paired/single reads)
  --shuffle              Shuffle records after encoding (for ML training data)
  --seed N               One seed for every random decision: --shuffle order,
                         --npolicy random substitutions, unstranded
                         normalization's per-pair coin (default: 42)
  --shuffle-buffer-size N
                         Max memory per bucket for encode --shuffle (default: 1G).
                         Accepts K/M/G suffixes.
Overlap merging (--merge-pairs):
  --merge-pairs          Merge overlapping pairs in process while encoding
                         (takes exactly two FASTQ files, R1 R2, adjacent on the
                         command line; no FASTQ intermediate). Mate names are
                         checked per pair — input that silently mispaired under
                         the plain two-file reader now fails loudly.
  --threshold-merge N    Overlap score (bits) to merge (default: 28)
  --threshold-trim N     Overlap score (bits) to trim the redundant overlap
                         (default: 8)
  --min-read-length N    Drop emitted reads shorter than this (default: 40)
  --no-sync-check        Skip the per-pair mate-name check
  --merge-json PATH      Write the merge run's statistics as JSON (same block
                         as `zna merge --json`)
  --merge-backend {auto,accel,python}
                         Merge kernel (default: auto)
  --allow-empty          Permit an input with no read pairs

Metadata:
  --read-group TEXT      Read group ID (default: "Unknown")
  --description TEXT     Description string
  --strand-specific      Flag library as strand-specific (default: R1 antisense, R2 sense)
  --strand-normalize     Enable strand normalization (RC reads to consistent strand).
                         With --strand-specific: deterministic (antisense reads RC'd).
                         Without: random RC (for unstranded data).
  --read1-sense          Read 1 represents sense strand
  --read1-antisense      Read 1 represents antisense strand (default when --strand-specific)
  --read2-sense          Read 2 represents sense strand (default when --strand-specific)
  --read2-antisense      Read 2 represents antisense strand
  --npolicy {trim3,random}
                         Policy for no-call (N) bases (default: trim3):
                         - trim3: cut the read at its first N, keeping [0, N)
                         - random: substitute a base from a seeded stream
                           (--seed), reproducibly

Format Options:
  -o, --output FILE      Output file (default: stdout)
  --seq-len-bytes N      Bytes for sequence length: 1, 2, or 4 (default: 2)
  --block-size SIZE      Block size; accepts K/M/G suffixes (default: 4M)
  --zstd                 Force Zstd compression
  --uncompressed         Force uncompressed
  --level N              Zstd compression level 1-22 (default: 9)
```

### `zna decode`

Convert ZNA to FASTA format.

**Usage:**

```
zna decode [FILE] [OPTIONS]

Positional Arguments:
  FILE                   Input ZNA file (default: stdin)

Options:
  -o, --output FILE      Output FASTA file. Use '#' for split R1/R2
  -q, --quiet            Suppress progress messages
  --gzip                 Force gzip compression for stdout
  --restore-strand       Restore original strand orientation for antisense reads
  --labels               Append each record's label values as SAM-style tags
                         (labeled files only)
```

### `zna inspect`

Display ZNA file statistics.

**Usage:**

```
zna inspect FILE [--verify] [--json] [--counts] [--blocks]

  input FILE             Input ZNA file to inspect
  --verify               Fully certify the file: decode every block (zstd
                         content checksums verify in the same motion), recount
                         every stat, compare against the trailer, check
                         fragment adjacency. Exit 0 iff the file passes.
  --json                 Machine-readable output; splices the provenance and
                         trailer payloads verbatim, and the verify verdict
                         under --verify.
  --counts               Also report per-flag record counts (paired R1, paired R2,
                         single/merged, reverse-complemented). Reads block payloads,
                         so slower than the default scan.
  --blocks               (--json only) include the per-block index.
```

A bare `inspect` already runs the cheap structural pass on any 0.5 file —
footer CRC, stored block index against a walk of the actual headers — and
exits non-zero on disagreement. `--verify` is the paid tier; it is what a
corpus stager should run after every fetch. A file written by zna < 0.5 (or
an encode that crashed mid-write) has no trailer: it reads normally, but
`--verify` refuses it — the trailer is the writer's signature that the encode
finished.

### `zna shuffle`

Randomly shuffle records in a ZNA file with bounded memory usage. Preserves paired-end read associations.

**Usage:**

```
zna shuffle INPUT -o OUTPUT [OPTIONS]

Positional Arguments:
  INPUT                  Input ZNA file to shuffle

Options:
  -o, --output FILE      Output ZNA file (required)
  -s, --seed N           Random seed for reproducibility (default: 42)
  -b, --buffer-size SIZE Maximum memory per bucket (default: 1G)
                         Accepts K/M/G suffixes (e.g., 512M, 2G)
  --block-size SIZE      Block size for output ZNA (default: 4M)
  --tmp-dir DIR          Directory for temporary files (default: system temp)
  -q, --quiet            Suppress progress messages
```

**Algorithm**: Uses bucket shuffle with bounded memory:
1. Randomly distributes records into K temporary bucket files on disk
2. Shuffles each bucket in memory using Fisher-Yates algorithm
3. Concatenates shuffled buckets to produce uniform random permutation

**Examples:**

```bash
# Shuffle with default settings (1GB memory, seed 42)
zna shuffle input.zna -o shuffled.zna

# Shuffle with custom seed for reproducibility
zna shuffle input.zna -o shuffled.zna --seed 12345

# Shuffle with limited memory (512MB buffer)
zna shuffle input.zna -o shuffled.zna --buffer-size 512M

# Shuffle paired-end data (pairs stay together)
zna shuffle paired.zna -o shuffled_paired.zna
```

**Note**: Paired-end reads (R1+R2) are kept together as a single shuffle unit.

### `zna merge`

Overlap-merge paired-end reads into one mixed interleaved FASTQ, ready for
`zna encode --interleaved`. Replaces fastp's PE-merge step.

Each pair is scored **once**: R1 is slid against `revcomp(R2)` over the single axis of
candidate fragment lengths, and every shift gets a log-likelihood ratio in **bits** —
`+1.99` per matching base (that is `log2 4`, the information in agreeing on one of four
bases), `-6.23` per mismatch at a 1% error rate. Both weights fall out of the error
rate; neither is tuned. The best-scoring shift (`argmax`, not fastp's first-accept) is
then read at two thresholds:

| | condition | action |
|---|---|---|
| **merge** | `score ≥ --threshold-merge` | emit one full-fragment record |
| **trim** | `--threshold-trim ≤ score < merge` | keep both; split the redundant overlap between their **3'** ends |
| **keep** | `score < --threshold-trim` | keep both, untouched |

Three parameters, all with units. Both thresholds read one calibrated scale, so `T` bits
tolerates a spurious rate of about `N · 2^-T` over the `N ≈ 2 · readlen` candidate
shifts — the default 28 is one spurious merge in 10⁶ pairs *against chance alignment*
(measured: 0 in 40,000 uniform-random pairs, at every read length from 50 to 300). It is
not a bound against real sequence, where reads share genuine homology and repeat
content. Trim sits far lower only because a wrong trim deletes bases from a read tail
while a wrong merge invents sequence.

The overlap sits at the 3' end of **both** mates — each read starts at a fragment end
and reads inward — so a trim splits it between them. The emitted pair tiles the fragment
exactly once and comes out at equal length, and where the mates disagree both carry the
consensus call.

#### Choosing `--threshold-merge`, measured against ground truth

The defaults are not a guess. On 1,000,000 simulated pairs from hg38 with the true
fragment length known exactly ([docs/MERGE_BENCHMARK_RESULTS.md](docs/MERGE_BENCHMARK_RESULTS.md)),
against fastp 1.1.0 at its own defaults:

| setting | chimera rate¹ | sensitivity² | merges that are wrong | reconstructed exactly³ |
|---|---:|---:|---:|---:|
| **`--threshold-merge 28`** (default) | 1.231% | **99.83%** | 0.96% | 86.59% |
| `--threshold-merge 60` | 0.597% | 92.63% | 0.55% | 88.87% |
| `--threshold-merge 100` | 0.245% | 83.57% | 0.29% | 91.39% |
| fastp defaults | 0.621% | 92.98% | 0.65% | 85.90% |

¹ fraction of pairs with **no true overlap** that were merged anyway — the false-positive
rate. ² fraction of pairs with a true overlap ≥ 15 bases that merged. ³ merged records
equal to the true fragment, base for base.

**If you want fastp's false-positive rate, use `--threshold-merge 60.`** That is not a
coincidence: 60 bits is 31 clean bases, which is essentially fastp's
`--overlap_len_require 30`. At that matched operating point zna's sensitivity is the
same (92.63% vs 92.98%) and its reconstruction is better (88.87% vs 85.90% exact),
because the overlap consensus recovers 90.4% of recoverable overlap errors against
fastp's 74.1%.

**For best overall accuracy, keep the default 28.** It minimises false positives plus
false negatives by a wide margin — 6,603 total errors per million pairs against 44,145
for the fastp-equivalent setting. Raising the threshold trades ~10.9 extra missed merges
for every wrong merge it prevents at 28→34, worsening to 15.7 at 28→60, so it only pays
if a chimera costs you more than ~11x what a missed merge does. A missed merge is not
lost data: the pair is still emitted, correctly bounded and with its redundant overlap
trimmed.

**Tuning cannot reach zero.** At 100 bits — 3.6x the default — 1,403 wrong merges per
million remain. Every one is a fragment whose two ends are genuinely homologous (median
88% identity over 79 bases, hotspots entirely pericentromeric), and the scan never
picks a lower-scoring alignment than the true one. That residue is a property of the
genome, not of the threshold.

**Usage:**

```
zna merge --in1 R1.fq --in2 R2.fq --out OUT.fq [OPTIONS]

Required:
  --in1 FILE             R1 FASTQ (optionally .gz)
  --in2 FILE             R2 FASTQ (optionally .gz), positionally synced with --in1
  --out FILE             Output mixed interleaved FASTQ (.gz to gzip)

Options:
  --json FILE            Write run statistics as JSON (counts, histograms, provenance)
  --threshold-merge BITS Score at or above this merges the pair (default: 28.0)
  --threshold-trim BITS  Score at or above this (but below merge) keeps both mates
                         and splits the redundant overlap between their 3' ends
                         (default: 8.0)
  --npolicy {trim3,random}
                         No-call policy for an N the overlap could not rescue
                         from the mate (default: trim3; same flag and values
                         as zna encode)
  --seed N               Seed for --npolicy random (default: 42)
  --min-read-length N    Drop emitted reads shorter than this (default: 40)
  --threads N            Merge worker threads (default: min(4, cpu count))
  --io-threads N         pigz threads for the gzipped output (default: 4)
  --chunk-size N         Read pairs per work unit (default: 2000)
  --compress-level N     pigz level for --out (default: 1 — it is an intermediate)
  --backend NAME         auto (default), accel, or python
  --no-sync-check        Skip the per-pair R1/R2 read-name consistency check
  --allow-empty          Exit 0 on an input with no read pairs
  -q, --quiet            Suppress progress logging
```

**Speed.** The merge kernel is compiled C++ and releases the GIL, so `--threads` are
real worker threads. It is not usually the bottleneck — inflate is — so **2 threads
saturate and more does nothing**:

| | µs/pair |
|---|---:|
| `--threads 1` | 2.78 |
| `--threads 2` | **1.40** |
| `--threads 4` | 1.43 |

With gzip removed from both ends the tool runs at 0.42 µs/pair, so on compressed input
it is I/O bound. (Apple silicon, 8 cores.)

> **Install `zna[fast]` if your input is gzipped.** Inflate is the largest single cost of
> a merge — a profile at `--threads 1` puts a quarter of all cycles inside `pigz`'s
> `libz`. The `fast` extra pulls in `isal` (Intel ISA-L), which the reader prefers when
> present: **448 MB/s of inflate against `pigz`'s 193**, and on a 1M-pair library
> `--threads 1` goes from 5.84 s to 4.04 s wall while also using less CPU. It is
> optional — without it, `pigz` is used when on `PATH`, then stdlib `gzip`, exactly as
> before. `zna merge` logs which one ran and records it in `--json`.
>
> **0.5.2 also made the scan itself 2.04x faster on x86-64.** Every earlier number here
> was measured on aarch64; on Linux/x86 the kernel had been reducing each vector compare
> through a PLT call to libgcc's `__popcountdi2`. See CHANGELOG and `docs/METHODS.md` §2.4.

**Determinism.** The score is computed in fixed-point integers and the argmax has a
specified tie-break, so a given FASTQ produces **byte-identical output on any platform,
compiler and thread count**. `--backend python` selects the pure-Python reference
implementation, which exists as an oracle for the compiled one; it is ~50x slower and
is never chosen for you.

**Output** is one stream mixing both shapes: merged reads as single records with the
`/1`,`/2` suffix stripped, unmerged pairs as adjacent `/1`,`/2` records. Feed it to
`zna encode --interleaved --treat-unpaired-as-merged`, which is exact here because
merged records span their fragment and unmerged pairs are emitted all-or-nothing —
never a lone mate.

#### Per-record provenance

The run summary says what happened to a *library*. These say what happened to a *read*.
Existing header fields are always passed through untouched — provenance is **appended**,
never substituted — so `--label` reads the same tags off an emitted record that it would
have read off the input.

A record that nothing happened to is emitted unchanged, so on a clean library this costs
nothing.

```
@SRR1.7  ZI:i:42 ZN:i:6 trim3_12 rescued_1 merged_90_0
          ^tag    ^bits  ^bases cut  ^no-calls   ^fastp-style,
          yours          by trim3    recovered    stays last
                                     from the mate
```

| token | meaning |
|---|---|
| `trim3_<n>` / `subn_<n>` | bases removed by `--npolicy trim3`, or substituted by `--npolicy random` |
| `rescued_<n>` | no-calls this record recovered from its mate, which cost nothing |
| `merged_<n1>_<n2>` | bases contributed by R1 and R2 to a merged record |

The word tokens are for reading; **`ZN:i:<bits>` is the one that survives encoding.** ZNA
does not store headers, so it is the only per-record provenance that reaches a corpus:

```bash
zna encode --interleaved --treat-unpaired-as-merged \
           --label provenance:C:ZN -o reads.zna merged.fq.gz
```

That is an ordinary label column — declare it and you get one byte per record, omit it
and nothing changes. There is no provenance-specific code in the encoder.

| bit | | set when |
|---:|---|---|
| 1 | trimmed | the pair's redundant overlap was split between its mates |
| 2 | rescued | ≥1 no-call was recovered from the mate |
| 4 | N-trimmed | ≥1 base was removed by `--npolicy trim3` |
| 8 | N-substituted | ≥1 base was invented by `--npolicy random` |

There is deliberately **no "merged" bit**: that fact already has two homes, the
`merged_` token here and `IS_FULL_FRAGMENT` in the corpus. Every bit above is one with
nowhere else to live — a trimmed pair in particular is emitted as an ordinary pair, and
nothing in the ZNA flag byte distinguishes it from one kept whole.

**Examples:**

```bash
# Defaults suit 2x150 bp data; you normally set nothing
zna merge --in1 R1.fq.gz --in2 R2.fq.gz --out merged.fq.gz

# With run statistics for a pipeline to collect
zna merge --in1 R1.fq.gz --in2 R2.fq.gz --out merged.fq.gz --json merge.json

# Straight into a training corpus
zna merge --in1 R1.fq.gz --in2 R2.fq.gz --out merged.fq.gz
zna encode --interleaved --treat-unpaired-as-merged --strand-normalize \
           --shuffle merged.fq.gz -o reads.zna
```

**Boundary guarantee.** Base 0 of every emitted read is a true fragment boundary —
nothing is ever removed from a read's 5' end, and trimming only ever cuts 3' ends. A
merged record is the tool's assertion of its fragment. This is what makes `IS_RC` and
`IS_FULL_FRAGMENT` honest for merged input; verified at 0 violations over 1,416,630
records against genome truth. See [docs/METHODS.md](docs/METHODS.md) for the derivation
and [docs/MERGE_BENCHMARK_RESULTS.md](docs/MERGE_BENCHMARK_RESULTS.md) for the
verification.

---

## Performance Characteristics

### Compression Ratios

Typical compression ratios compared to raw FASTQ:

| Format | Size | Ratio | Notes |
|--------|------|-------|-------|
| FASTQ (uncompressed) | 100% | 1.0x | Baseline |
| FASTQ.gz (gzip -6) | 25-30% | 3-4x | Standard |
| ZNA (uncompressed) | 12-15% | 6-8x | 2-bit encoding only |
| ZNA (Zstd L3) | 8-10% | 10-12x | Fast compression (`--level 3`) |
| ZNA (Zstd L9) | 6-8% | 12-16x | **Default** (`DEFAULT_ZSTD_LEVEL = 9`) |

*Results vary based on sequence complexity and redundancy*

### Speed

- **Encoding**: ~4.8M reads/second at 150 bp (single thread)
- **Decoding**: ~11.5M reads/second at 150 bp (single thread)
- **Block-based**: shards and subsamples without decoding what it skips

### Memory Usage

- **Streaming I/O**: constant memory for `records()`; the writer buffers one
  block at a time
- **Default block size**: 4 MiB (`--block-size`)
- **`blocks()` holds one decoded block per open reader**, so block size sets the
  memory of a batch consumer: at 150 bp, a 4 MiB block is ~100k records (~20 MB
  of Python strings) against ~25k (~5 MB) for 1 MiB
- **No stored index required**: `block_index()` walks block headers in ~2.3 µs
  per block, so counts and offsets are recovered without one

---

## Technical Details

### 2-Bit Encoding

DNA bases are encoded in 2 bits:

```
A = 00 = 0
C = 01 = 1
G = 10 = 2
T = 11 = 3
```

Four bases pack into one byte:
```
Byte: [B1][B2][B3][B4]
      76 54 32 10  (bit positions)
```

### Lookup Tables

Pre-computed lookup tables provide O(1) encoding/decoding:

- **Encoding**: 256-element array mapping ASCII → 2-bit
- **Decoding**: 256-element tuple mapping byte → 4-character string

### Block-Based Architecture

Data is organized in independently compressed blocks:

- **Advantages**: block-granular sharding and sampling; per-block counts without
  decompression
- **Overhead**: 20 bytes of header per block
- **Choosing a size**: 4 MiB (the default) maximises compression on
  duplicate-rich data. On unique reads, packed sequence is incompressible and
  block size costs nothing, so prefer smaller blocks (1 MiB) when a consumer
  reads with `blocks()` or shards by block

### Compression Strategy

- **Zstd**: Modern compression algorithm (Facebook)
- **Reusable compressor**: Amortizes initialization cost
- **Memoryview parsing**: Zero-copy decompression
- **Pre-sized buffers**: Eliminates reallocations

---

## Strand-Specific Libraries

ZNA supports strand-specific RNA-seq libraries by normalizing all reads to sense strand orientation during encoding. This enables consistent downstream analysis while preserving the ability to restore original strand information.

### How It Works

1. **Encoding**: reads are reverse-complemented into one common frame. *Which* reads,
   and by what rule, depends on the mode — see **Strand Normalization** below: with
   `--strand-specific` the protocol decides, without it one mate per pair is chosen at
   random.
2. **Storage**: each record's `IS_RC` flag records whether it was flipped. That flag is
   the only record of it — it cannot be recovered from the sequence.
3. **Decoding**: `--restore-strand` consumes `IS_RC` to recover the original orientation.

### Strand Normalization

The `--strand-normalize` flag controls whether reads are reverse-complemented
to a consistent strand during encoding:

- **With `--strand-specific`**: Deterministic normalization — antisense reads
  are reverse-complemented to sense orientation based on the library protocol.
  Each read's IS_RC flag records whether it was flipped.
- **Without `--strand-specific`**: Random reverse-complementing for unstranded
  data.  Useful for data augmentation in ML training.
- **Without `--strand-normalize`**: Reads are stored in their original
  orientation (no reverse-complementing).

```bash
# Strand-normalized encoding (most common for stranded RNA-seq)
zna encode R1.fq.gz R2.fq.gz --strand-specific --strand-normalize -o lib.zna

# Decode with original strand orientation restored
zna decode lib.zna --restore-strand -o original.fasta

# Decode with sense-normalized sequences (for alignment)
zna decode lib.zna -o normalized.fasta
```

### Unstranded Normalization and Fragment Geometry

Unstranded normalization does more than augment the data: it carries information
about the molecule that cannot be reconstructed afterwards.

A fastp-style FR pair covers the two ends of one fragment, pointing inward:

```
    fragment, length L
    |------------------------------------------------|
    |>>>>>>>>>>>|                        |<<<<<<<<<<<<|
     R1 as sequenced                      R2 as sequenced
     = F[0:l1]                            = revcomp(F[L-l2:L])
```

As sequenced the mates are in *opposite* frames.  Normalization
reverse-complements **exactly one** of them so both land in one common frame,
and records which one in that record's `IS_RC` flag:

```
    common frame after normalization
    |------------------------------------------------|
    |<<<<<<<<<<<|                        |<<<<<<<<<<<<|
     not RC'd                             RC'd
     LEFT edge  = real fragment boundary  RIGHT edge = real fragment boundary
     right edge = read-length cutoff      left edge  = read-length cutoff
```

**The invariant:** whichever mate was reverse-complemented ends up at the right
of the common frame, so its **right** edge is the real fragment boundary and its
left edge is a read-length cutoff.  For the other mate it is the mirror image.

`IS_RC` is the only thing that distinguishes the two cases, and **it cannot be
recovered from the sequence**.  Reverse-complementing the right-hand mate
reproduces the fragment-frame sequence exactly, because that mate was stored
reverse-complemented to begin with — there is no residue in the bases to test.
The coin is also independent of the mate number, so `is_read1` is not a
substitute for it.

**Reading the geometry.** Use `records(with_ends=True)`, which answers the
question directly instead of making you re-derive it:

```python
with open("lib.zna", "rb") as f:
    reader = ZnaReader(f)
    for seq, is_paired, is_read1, is_read2, has_start, has_end in \
            reader.records(with_ends=True):
        # has_start: the LEFT edge of seq is a true fragment boundary
        # has_end:   the RIGHT edge is
        ...
```

`records(with_rc=True)` exposes the raw `IS_RC` flag instead, if you want the
orientation itself rather than the boundary geometry.

**A record can have two real ends.** When the insert is at or below the read
length — every overlap-merged read, and any pair after adapter trimming — the
record spans the *whole* fragment and both edges are true boundaries. `IS_RC`
names only one edge, so that case is carried by a separate flag,
`IS_FULL_FRAGMENT`, which `with_ends` folds in for you. Full-overlap **pairs**
are detected automatically at encode time (mates covering the same interval are
exact reverse complements); for **unpaired** records the encoder cannot tell a
merged read from a genuine single-end read, so declare it:

```bash
# reads from an overlap merger: unpaired records span their whole fragment
zna encode --interleaved --treat-unpaired-as-merged -o out.zna merged.fq.gz
```

Without the flag an unpaired record is assumed to have one real edge, which is
the safe reading — a tool marking fragment ends will under-label rather than
place a marker at an interior position.

`restore_strand=True` is not a substitute: it *consumes* the flag to undo the
reverse-complement and hand back original-orientation reads.  A caller that
wants the normalized frame *and* the boundary geometry needs `with_rc`, and the
two options are mutually exclusive.

**Normalization happens once, at encode time, and is not idempotent.** Applying
it a second time returns the data to an un-normalized state while the header
still reports `strand_normalized`.  So anything that copies records between ZNA
files — `zna encode` on a `.zna` input, `zna shuffle` — copies the existing
orientation rather than re-deriving it.

**A view is for reading; the flag byte is for copying.**  `records()` returns
views — each of them chosen for a consumer, and none of them able to carry the
whole flag byte back to a writer.  Copying uses `copy_records()`:

```python
# A lossless ZNA -> ZNA copy.
with open("in.zna", "rb") as fin, open("out.zna", "wb") as fout:
    reader = ZnaReader(fin)
    with ZnaWriter(fout, reader.header, preserve_normalization=True) as writer:
        for rec in reader.copy_records():
            writer.write_copy(rec)
```

`copy_records()` yields `ZnaRecord(seq, flags, labels)` — the stored
`ZnaRecordFlags` byte, verbatim — so a copy carries every bit, including ones
this version does not interpret.  This example used to read
`records(with_ends=True)` into `write_records()`, and that was wrong:
`(has_start, has_end)` has three states where `(IS_RC, IS_FULL_FRAGMENT)` has
four, so every full-fragment record came out of the copy with `IS_RC` cleared.
`write_records()` now refuses that shape rather than accepting it.

### Strand Flags

| Flag | Description |
|------|-------------|
| `--strand-specific` | Enable strand-specific mode (default: R1 antisense, R2 sense) |
| `--read1-sense` | Read 1 represents sense strand |
| `--read1-antisense` | Read 1 represents antisense strand |
| `--read2-sense` | Read 2 represents sense strand |
| `--read2-antisense` | Read 2 represents antisense strand |

### Common Library Protocols

| Protocol | R1 | R2 | ZNA Flags |
|----------|----|----|-----------|
| **dUTP / TruSeq Stranded** | antisense | sense | `--strand-specific` (default) |
| **Illumina Stranded mRNA** | antisense | sense | `--strand-specific` |
| **fr-firststrand** | antisense | sense | `--strand-specific` |
| **fr-secondstrand** | sense | antisense | `--strand-specific --read1-sense --read2-antisense` |
| **Ligation (ScriptSeq)** | sense | antisense | `--strand-specific --read1-sense --read2-antisense` |

### Examples

```bash
# dUTP/TruSeq protocol (most common - this is the default)
zna encode R1.fastq.gz R2.fastq.gz --strand-specific -o library.zna

# fr-secondstrand protocol
zna encode R1.fastq.gz R2.fastq.gz \
  --strand-specific --read1-sense --read2-antisense \
  -o library.zna

# Decode with sense-normalized sequences (for alignment)
zna decode library.zna -o normalized.fasta

# Decode with original strand orientation restored
zna decode library.zna --restore-strand -o original.fasta
```

---

## Per-Sequence Labels

ZNA can store numeric metadata as compact columnar label columns alongside
each sequence.  Labels are parsed from key-value tags in FASTQ headers
(e.g. output from `samtools fastq -T`), but the tag format is not limited
to SAM — any `KEY:TYPE:VALUE` field in the header will work, and keys
can be any length.

### Defining Labels on the CLI

```bash
# Two labels with descriptions
zna encode reads.fq.gz -o reads.zna \
  --label NH:C --label AS:i \
  --label-desc NH:"Number of hits" --label-desc AS:"Alignment score"
```

The `--label` format is `NAME:TYPE` where TYPE is one of:
`A`, `c`, `C`, `s`, `S`, `i`, `I`, `f`, `d`, `q`, `Q`.
Use the smallest type that fits your data to minimize file size.

#### Decoupled Name and Tag

By default, the label **name** (stored in the ZNA header) is also used as
the **tag** to parse from input.  You can decouple these with the
3-part format `NAME:TYPE:TAG`:

```bash
# Store as "edit_dist" in ZNA, but parse "NM" tag from input headers
zna encode reads.fq.gz -o reads.zna \
  --label edit_dist:C:NM --label aln_score:i:AS

# Custom long-form tags (not SAM) work too
zna encode reads.fq.gz -o reads.zna \
  --label score:i:alignment_score --label edits:C:edit_distance
```

The `tag` is only used at encode time and is **not stored** in the ZNA file.
When decoding, the label `name` is used in the output.

### Defining Labels with a YAML File

For many labels, define them in a YAML file instead of many CLI flags:

```bash
zna encode reads.fq.gz -o reads.zna --label-defs labels.yaml
```

```yaml
# labels.yaml
labels:
  - name: NM
    type: C
    description: Edit distance
    missing: 255
  - name: aln_score
    type: i
    tag: AS                    # parse "AS" from input, store as "aln_score"
    description: Alignment score
    missing: -1
```

CLI flags `--label` and `--label-desc` **override** values from the YAML
file, so you can keep a YAML base and tweak individual labels per run.

See [`examples/labels.yaml`](examples/labels.yaml) for a fully-commented
template.

### Decoding Labeled Files

```bash
# Include labels as SAM-style tags in the output
zna decode reads.zna --labels > output.fq

# Inspect to see label definitions
zna inspect reads.zna
```

### Python API with Labels

```python
from zna.core import ZnaHeader, ZnaWriter, ZnaReader
from zna.dtypes import LabelDef, parse_dtype

defs = (
    LabelDef(0, "NM", "Edit distance", parse_dtype("C"), missing=255),
    LabelDef(1, "AS", "Alignment score", parse_dtype("i"), missing=-1),
)
header = ZnaHeader(read_group="sample", labels=defs)

with open("out.zna", "wb") as f:
    with ZnaWriter(f, header) as w:
        w.write_record("ACGT", is_paired=False,
                        is_read1=False, is_read2=False,
                        labels=(3, 280))

with open("out.zna", "rb") as f:
    reader = ZnaReader(f)
    for seq, is_paired, is_r1, is_r2, labels in reader.records():
        print(seq, labels)  # ACGT (3, 280)
```

Labeled files yield a 5-tuple ending in `labels`.  With `with_rc=True` the
`is_rc` flag is inserted *before* it — `(seq, is_paired, is_read1, is_read2,
is_rc, labels)` — so that the unlabeled and labeled tuples agree on where
`is_rc` lives.

---

## Use Cases

### Recommended For

- ✅ **Long-term archival**: High compression with fast retrieval
- ✅ **Data transfer**: Reduced bandwidth requirements
- ✅ **Cloud storage**: Lower storage costs
- ✅ **Pipeline integration**: Unix-friendly streaming
- ✅ **Reference storage**: Efficient genome/transcriptome storage

### Not Recommended For

- ⚠️ **Record-level random access**: there is no record index. Access is *block*
  granular — `block_index()` gives per-block offsets and counts without decompressing,
  and `blocks(indices=...)` reads only the blocks you ask for. That is what sharded
  training uses; seeking to record *n* is what is not supported.
- ❌ **Quality scores**: Sequences only (use CRAM/BAM for qualities)
- ❌ **Small files**: Overhead outweighs benefits (<10K reads)
- ❌ **Real-time streaming**: Use case requires quality scores

---

## Comparison with Other Formats

| Feature | ZNA | FASTA | FASTQ | CRAM | FASTA.gz |
|---------|-----|-------|-------|------|----------|
| Compression | Excellent | None | None | Excellent | Good |
| Speed | Fast | Fastest | Fast | Slow | Medium |
| Quality Scores | ❌ | ❌ | ✅ | ✅ | ❌ |
| Paired-End | ✅ | ❌ | ❌ | ✅ | ❌ |
| Random Access | ❌ | ✅ | ✅ | ✅ | ❌ |
| Streaming | ✅ | ✅ | ✅ | Limited | ✅ |
| Dependencies | 2 | 0 | 0 | Many | 0 |

---

## Python API

In addition to the CLI, ZNA provides a Python API:

```python
from zna import ZnaHeader, ZnaWriter, ZnaReader, COMPRESSION_ZSTD

# Writing
header = ZnaHeader(
    read_group="Sample_01",
    compression_method=COMPRESSION_ZSTD,
    compression_level=5
)

with open("output.zna", "wb") as f:
    with ZnaWriter(f, header) as writer:
        writer.write_record("ACGTACGT", is_paired=False, 
                          is_read1=False, is_read2=False)
        writer.write_record("TGCATGCA", is_paired=False,
                          is_read1=False, is_read2=False)

# Paired reads: write the fragment whole, R1 immediately followed by R2.
with open("paired.zna", "wb") as f:
    with ZnaWriter(f, header) as writer:
        writer.write_record("ACGTACGT", is_paired=True,
                          is_read1=True, is_read2=False)     # R1
        writer.write_record("TGCATGCA", is_paired=True,
                          is_read1=False, is_read2=True)     # its R2

# Reading
with open("output.zna", "rb") as f:
    reader = ZnaReader(f)
    print(f"Read Group: {reader.header.read_group}")
    
    for seq, is_paired, is_read1, is_read2 in reader.records():
        print(seq)

# Copying to another ZNA file: carry the flag byte, not a view of it
with open("output.zna", "rb") as fin, open("copy.zna", "wb") as fout:
    reader = ZnaReader(fin)
    with ZnaWriter(fout, reader.header, preserve_normalization=True) as writer:
        for rec in reader.copy_records():
            writer.write_copy(rec)
```

**The writer requires whole fragments.** A paired R1 must be followed immediately
by its R2 — on `write_record`, `write_records`, and `write_copy` alike — and the
stream may not end on an R1. Anything else raises `ValueError` naming the record.
This is what lets the writer keep every fragment inside one block, which the
`blocks()` sharding above depends on. Unpaired and merged reads are one-record
fragments and need nothing special.

`records()` yields a 4-tuple, or a 5-tuple ending in `labels` for labeled files.
Two options change what it yields:

| Option | Yields | Purpose |
|--------|--------|---------|
| *(default)* | `(seq, is_paired, is_read1, is_read2)` | stored orientation |
| `restore_strand=True` | same 4-tuple | undoes strand normalization, returning original-orientation reads |
| `with_rc=True` | `(seq, is_paired, is_read1, is_read2, is_rc)` | stored orientation plus the per-record `IS_RC` flag |
| `with_ends=True` | `(seq, is_paired, is_read1, is_read2, has_start, has_end)` | which edges are true fragment boundaries |

The options are mutually exclusive: `restore_strand` consumes the orientation,
`with_rc` returns it raw, and `with_ends` returns what it means.

All three are **views**, for consumers.  None of them round-trips back into a
writer — `with_ends` in particular folds `IS_RC` and `IS_FULL_FRAGMENT` into two
booleans that cannot distinguish all four reachable states.  To copy records
between ZNA files use [`copy_records()` / `write_copy()`](#unstranded-normalization-and-fragment-geometry),
which carry the flag byte itself.  See
[Unstranded Normalization and Fragment Geometry](#unstranded-normalization-and-fragment-geometry)
for what `is_rc` means and why it cannot be derived from the sequence.

---

## Development

### Running Tests

```bash
# All tests
PYTHONPATH=src pytest -v

# Specific test suite
PYTHONPATH=src pytest tests/test_cli.py -v
PYTHONPATH=src pytest tests/test_core.py -v

# With coverage
PYTHONPATH=src pytest --cov=zna tests/
```

### Code Quality

```bash
# Format code
black src/ tests/

# Type checking
mypy src/zna/
```

---

## Limitations

1. **No quality scores or read names**: sequences, flags, and declared numeric
   labels only. (Per-record numeric metadata IS supported — see
   [Per-Sequence Labels](#per-sequence-labels).)
2. **Block-granular access**: random access is by block, not by record —
   `blocks(indices=…)` seeks straight to any block; there is no per-record
   index.
3. **DNA/RNA only**: A, C, G, T stored; `N` handled by `--npolicy` at encode
   time; other IUPAC codes are rejected.
4. **Case insensitive**: lowercase is uppercased on encode.
5. ~~No index~~ Since 0.5.0 the trailer stores record/base counts, length
   histograms, and the block index — `zna inspect` reads them in O(1).

---

## Future Enhancements

See [docs/ROADMAP.md](docs/ROADMAP.md) — what is scheduled, what is under consideration,
and what has already been tried and closed by measurement.

---

## License

GNU General Public License v3.0 **or later** (`GPL-3.0-or-later`). See
[LICENSE](LICENSE).

---

## Citation

If you use ZNA in your research, please cite:

```
Iyer, M. (2026). ZNA: A compressed binary format for nucleic acid sequences.
GitHub: https://github.com/mkiyer/zna
```

---

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Ensure all tests pass
5. Submit a pull request

---

## Contact

- **Author**: Matthew Iyer
- **Email**: mkiyer@umich.edu
- **Issues**: https://github.com/mkiyer/zna/issues
