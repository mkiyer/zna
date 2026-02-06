# Changelog

All notable changes to the ZNA project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
