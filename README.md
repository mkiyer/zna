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

Comprehensive documentation is available in the [docs/](docs/) directory:

- **[docs/RELEASING.md](docs/RELEASING.md)** - Quick release guide
- **[docs/PERFORMANCE.md](docs/PERFORMANCE.md)** - Performance benchmarks
- **[docs/PUBLISHING.md](docs/PUBLISHING.md)** - PyPI publishing guide

See [docs/README.md](docs/README.md) for a complete list.

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
- **Flags** (1 byte): IS_READ1, IS_READ2, IS_PAIRED
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
  --fasta                Force FASTA format (overrides extension detection)
  --fastq                Force FASTQ format (overrides extension detection)

Metadata:
  --read-group TEXT      Read group ID (default: "Unknown")
  --description TEXT     Description string
  --strand-specific      Flag library as strand-specific (default: R1 antisense, R2 sense)
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
zna inspect FILE

  input FILE             Input ZNA file to inspect
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
