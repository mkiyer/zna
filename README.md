# ZNA: Compressed Nucleic Acid Format

**ZNA** (Compressed **Z**-Nucleic **N**-Acid **A**) is a high-performance binary format for storing DNA/RNA sequences with exceptional compression and I/O speed.

## Performance

- **135 MB/s roundtrip** throughput (9.5x faster than Python baseline)
- **2.8+ GB/s** encoding/decoding for long reads
- **3.7-4.0x compression** ratio with Zstd
- **C++ acceleration** with pure Python fallback

## Features

- **High Compression**: 2-bit encoding (4 bases per byte) + optional Zstd compression
- **Ultra-Fast I/O**: C++ accelerated encode/decode with block-based architecture
- **Minimal Dependencies**: `zstandard` only (C++ extension auto-compiled)
- **Flexible**: Single-end, paired-end, and interleaved reads
- **Strand-Specific Support**: dUTP, TruSeq, and custom strand protocols
- **Built-in Shuffle**: Memory-bounded random shuffling for training data preparation
- **Metadata Rich**: Read groups, descriptions, and custom flags
- **Unix-Friendly**: Pipe-compatible CLI for seamless workflow integration
- **Streaming**: Memory-efficient block-based processing

## Installation

```bash
# From source (recommended - includes C++ acceleration)
git clone https://github.com/mkiyer/zna.git
cd zna
pip install -e .

# Check if C++ acceleration is available
python -c "from zna.core import is_accelerated; print(f'Accelerated: {is_accelerated()}')"
```

**Requirements:**
- Python ≥3.10
- C++ compiler (for optimal performance)
- CMake ≥3.15 (auto-installed via pip)

## Quick Start

```bash
# Encode FASTQ to compressed ZNA (default: Zstd level 3)
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

# Pipe-friendly workflows
cat reads.fastq | zna encode -o reads.zna
zna decode reads.zna | head -n 1000
```

## Performance Benchmarks

### Throughput by Read Length

| Read Type | Encode (MB/s) | Decode (MB/s) | Compression |
|-----------|---------------|---------------|-------------|
| Short (Illumina, 100-150bp) | 189.5 | 668.8 | 3.68x |
| Medium (300-500bp) | 540.5 | 1,280.9 | 3.87x |
| Long (PacBio, 1-5kb) | 1,921.5 | 2,864.6 | 3.98x |
| Very Long (Nanopore, 5-15kb) | **2,824.7** | **3,392.7** | 3.99x |

**Key Insights:**
- Performance scales dramatically with read length
- Compression ratio remains consistent across workloads
- C++ acceleration provides 9.5x speedup over pure Python

See [docs/PERFORMANCE.md](docs/PERFORMANCE.md) for detailed benchmarking.

---

## Documentation

- **[docs/RELEASING.md](docs/RELEASING.md)** - Publishing to PyPI and Bioconda
- **[docs/PERFORMANCE.md](docs/PERFORMANCE.md)** - Benchmarks and tuning

---

## File Format Specification

### Overview

ZNA files use a binary format optimized for nucleic acid sequences:

- **File Extension**: `.zna` (for both compressed and uncompressed files)
- **Default Compression**: Zstd level 3 (use `--uncompressed` flag to disable)
- **Magic Number**: `ZNA\x1A` (4 bytes)
- **Version**: 1 (1 byte)
- **2-bit Encoding**: A=00, C=01, G=10, T=11
- **Block Structure**: Data organized in compressed/uncompressed blocks
- **Metadata**: Read groups, descriptions, and custom information

### File Structure

```
┌─────────────────────────────────────┐
│         File Header                 │
│  - Magic (4 bytes)                  │
│  - Version (1 byte)                 │
│  - Sequence length encoding (1 byte)│
│  - Flags (1 byte)                   │
│  - Compression method (1 byte)      │
│  - Compression level (1 byte)       │
│  - Metadata lengths (6 bytes)       │
│  - Variable metadata strings        │
├─────────────────────────────────────┤
│         Block 1                     │
│  - Block Header (12 bytes)          │
│    * Compressed size (4 bytes)      │
│    * Uncompressed size (4 bytes)    │
│    * Record count (4 bytes)         │
│  - Compressed/Raw Payload           │
│    * Record 1: flags, length, seq   │
│    * Record 2: flags, length, seq   │
│    * ...                            │
├─────────────────────────────────────┤
│         Block 2                     │
│  ...                                │
└─────────────────────────────────────┘
```

### Record Format

Each record in a block contains:
- **Flags** (1 byte): IS_READ1 (bit 0), IS_READ2 (bit 1), IS_PAIRED (bit 2),
  IS_RC (bit 3 — set when strand normalization reverse-complemented this record),
  IS_FULL_FRAGMENT (bit 4 — the record spans its whole fragment, so *both* edges
  are true fragment boundaries). Bits 5-7 are reserved.
- **Length** (1-4 bytes): Sequence length (configurable)
- **Sequence** (variable): 2-bit encoded bases

### Compression

- **Method 0**: Uncompressed (`.zna`)
- **Method 1**: Zstd compression (`.zzna`, levels 1-22)
- **Block Size**: Default 128KB (configurable)

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
zna encode sample.fastq.gz -o sample.zzna

# With high compression (default is level 3)
zna encode sample.fastq --level 5 -o sample.zna

# Uncompressed (rarely needed)
zna encode sample.fastq --uncompressed -o sample.zna

# From stdin
cat sample.fastq | zna encode -o sample.zna

# Force format (when extension detection fails)
cat data.txt | zna encode --fastq -o sample.zna
```

#### Paired-End Reads

```bash
# Separate R1/R2 files
zna encode R1.fastq.gz R2.fastq.gz -o paired.zna

# Interleaved file (strict alternating R1/R2 pairs)
zna encode interleaved.fastq --interleaved -o paired.zna

# Interleaved from stdin
cat interleaved.fastq | zna encode --interleaved -o paired.zzna
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
zna encode sample.fastq --npolicy drop -o clean.zna       # Skip sequences with N
zna encode sample.fastq --npolicy random -o clean.zna     # Replace N with random base
zna encode sample.fastq --npolicy A -o clean.zna          # Replace N with A

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
  --seq-len-bytes 1 \  # Max 255 bp
  -o short_reads.zna

zna encode sample.fastq \
  --seq-len-bytes 2 \  # Max 65,535 bp (default)
  -o sample.zna

zna encode sample.fastq \
  --seq-len-bytes 4 \  # Max 4.2 billion bp
  -o long_reads.zna
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

`blocks()` also takes `restore_strand=True`. It raises on labeled files — the
label columns would have to come back too, and dropping them silently is worse
than not offering the API — so use `records()` there.

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
Compression:      ZSTD (Level 3)

--- Content Statistics ---
Total Blocks:       356
Total Records:      1000000
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
  --seed N               Random seed for --shuffle (default: 42)
  --shuffle-buffer-size N
                         Max memory per bucket for encode --shuffle (default: 1G).
                         Accepts K/M/G suffixes.
  --fasta                Force FASTA format (overrides extension detection)
  --fastq                Force FASTQ format (overrides extension detection)

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
  --npolicy {drop,random,A,C,G,T}
                         Policy for handling 'N' nucleotides:
                         - drop: skip sequences containing N
                         - random: replace N with random base (A/C/G/T)
                         - A/C/G/T: replace N with specific base

Format Options:
  -o, --output FILE      Output file (default: stdout)
  --seq-len-bytes N      Bytes for sequence length: 1, 2, or 4 (default: 2)
  --block-size N         Block size in bytes (default: 131072)
  --zstd                 Force Zstd compression
  --uncompressed         Force uncompressed
  --level N              Zstd compression level 1-22 (default: 3)
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
```

### `zna inspect`

Display ZNA file statistics.

**Usage:**

```
zna inspect FILE [--counts]

  input FILE             Input ZNA file to inspect
  --counts               Also report per-flag record counts (paired R1, paired R2,
                         single/merged, reverse-complemented). Reads block payloads,
                         so slower than the default header-only scan.
```

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

---

## Performance Characteristics

### Compression Ratios

Typical compression ratios compared to raw FASTQ:

| Format | Size | Ratio | Notes |
|--------|------|-------|-------|
| FASTQ (uncompressed) | 100% | 1.0x | Baseline |
| FASTQ.gz (gzip -6) | 25-30% | 3-4x | Standard |
| ZNA (uncompressed) | 12-15% | 6-8x | 2-bit encoding only |
| ZNA (Zstd L3) | 8-10% | 10-12x | Fast compression (default) |
| ZNA (Zstd L9) | 6-8% | 12-16x | High compression |

*Results vary based on sequence complexity and redundancy*

### Speed

- **Encoding**: ~5-10M reads/second (single thread)
- **Decoding**: ~8-15M reads/second (single thread)
- **Block-based**: Enables parallel processing (future)

### Memory Usage

- **Streaming I/O**: Constant memory usage
- **Default block size**: 128KB buffer
- **No index required**: Sequential scan

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

- **Advantages**: Random access, parallel processing potential
- **Overhead**: ~12 bytes per block
- **Optimal size**: 128KB balances compression ratio and I/O

### Compression Strategy

- **Zstd**: Modern compression algorithm (Facebook)
- **Reusable compressor**: Amortizes initialization cost
- **Memoryview parsing**: Zero-copy decompression
- **Pre-sized buffers**: Eliminates reallocations

---

## Strand-Specific Libraries

ZNA supports strand-specific RNA-seq libraries by normalizing all reads to sense strand orientation during encoding. This enables consistent downstream analysis while preserving the ability to restore original strand information.

### How It Works

1. **Encoding**: Antisense reads are reverse-complemented to sense strand
2. **Storage**: All reads stored in sense orientation
3. **Decoding**: Use `--restore-strand` to recover original orientation

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
orientation rather than re-deriving it.  In the Python API that is
`ZnaWriter(..., preserve_normalization=True)` fed from `records(with_rc=True)`:

```python
# A lossless ZNA -> ZNA copy.
with open("in.zna", "rb") as fin, open("out.zna", "wb") as fout:
    reader = ZnaReader(fin)
    with ZnaWriter(fout, reader.header, preserve_normalization=True) as writer:
        writer.write_records(reader.records(with_ends=True))
```

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
zna encode R1.fastq.gz R2.fastq.gz --strand-specific -o library.zzna

# fr-secondstrand protocol
zna encode R1.fastq.gz R2.fastq.gz \
  --strand-specific --read1-sense --read2-antisense \
  -o library.zzna

# Decode with sense-normalized sequences (for alignment)
zna decode library.zzna -o normalized.fasta

# Decode with original strand orientation restored
zna decode library.zzna --restore-strand -o original.fasta
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

- ❌ **Random access**: Sequential format (no index)
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
| Dependencies | 1 | 0 | 0 | Many | 0 |

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

with open("output.zzna", "wb") as f:
    with ZnaWriter(f, header) as writer:
        writer.write_record("ACGTACGT", is_paired=False, 
                          is_read1=False, is_read2=False)
        writer.write_record("TGCATGCA", is_paired=False,
                          is_read1=False, is_read2=False)

# Reading
with open("output.zzna", "rb") as f:
    reader = ZnaReader(f)
    print(f"Read Group: {reader.header.read_group}")
    
    for seq, is_paired, is_read1, is_read2 in reader.records():
        print(seq)
```

`records()` yields a 4-tuple, or a 5-tuple ending in `labels` for labeled files.
Two options change what it yields:

| Option | Yields | Purpose |
|--------|--------|---------|
| *(default)* | `(seq, is_paired, is_read1, is_read2)` | stored orientation |
| `restore_strand=True` | same 4-tuple | undoes strand normalization, returning original-orientation reads |
| `with_rc=True` | `(seq, is_paired, is_read1, is_read2, is_rc)` | stored orientation plus the per-record `IS_RC` flag |
| `with_ends=True` | `(seq, is_paired, is_read1, is_read2, has_start, has_end)` | which edges are true fragment boundaries; also the lossless form for copying |

The options are mutually exclusive: `restore_strand` consumes the orientation,
`with_rc` returns it raw, and `with_ends` returns what it means.  See
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

1. **Sequences only**: No quality scores, headers, or annotations
2. **Sequential access**: No random access without full scan
3. **DNA/RNA only**: A, C, G, T bases (N or IUPAC codes not supported)
4. **Case insensitive**: Lowercase converted to uppercase
5. **No index**: Full file scan required for record counting

---

## Future Enhancements

- [ ] Parallel compression/decompression
- [ ] Optional index for random access
- [ ] Support for IUPAC ambiguity codes
- [ ] Memory-mapped I/O for large files
- [ ] Streaming statistics (GC content, length distribution)

---

## License

GNU General Public License v3.0 (GPLv3)

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
