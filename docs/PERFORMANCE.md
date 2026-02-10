# ZNA Performance

Benchmarks on Apple Silicon, Python 3.14, February 2026.
Test dataset: 1M paired-end Illumina reads (150 bp) from chr22 transcriptome simulation, 645.8 MB uncompressed FASTQ.

## Summary

| Metric | Value |
|--------|-------|
| **Encode** | 165 MB/s |
| **Decode** | 241 MB/s |
| **Compression** | 9.10x (default) / 9.44x (512 KB blocks) |

## Compression Comparison

| Format | File Size | Ratio |
|--------|-----------|-------|
| FASTQ (uncompressed) | 645.8 MB | 1.0x |
| FASTQ.gz | 175 MB | 3.7x |
| **ZNA (128 KB blocks)** | **71.0 MB** | **9.10x** |
| **ZNA (512 KB blocks)** | **68.4 MB** | **9.44x** |
| BAM (typical) | 60–80 MB | 8–10x |

## Compression Level Tuning

| Level | Encode (MB/s) | Decode (MB/s) | Ratio | Notes |
|-------|---------------|---------------|-------|-------|
| 1 | 166 | 242 | 9.09x | Fastest encode |
| **3** | **165** | **242** | **9.10x** | **Default** |
| 5 | 159 | 241 | 9.13x | |
| 9 | 159 | 243 | 9.16x | |
| 15 | 124 | 237 | 9.29x | Best compression |

Recommendation: level 3 is the sweet spot; level 15 gains only 2% compression at 25% slower encode.

## Block Size Tuning

| Block Size | Encode (MB/s) | Decode (MB/s) | Ratio |
|------------|---------------|---------------|-------|
| 32 KB | 168 | 243 | 9.07x |
| **128 KB** | **165** | **240** | **9.10x** |
| **512 KB** | **162** | **243** | **9.44x** |

Recommendation: 128 KB for streaming, 512 KB for archival (+3.7% compression, negligible speed cost).

## Why It's Fast

The columnar format groups homogeneous data streams together, giving ZSTD much better pattern matching:
- **Flags stream**: paired/single-end patterns compress 500–1000x
- **Lengths stream**: uniform 150 bp reads compress ~1000x
- **Sequences stream**: 2-bit DNA data compresses 3–5x

Remaining hot-path bottleneck is input I/O (~75% of encode time is reading gzipped FASTQ).

## Reproducing Benchmarks

```bash
python scripts/benchmark_columnar.py --all          # Full suite (~5 min)
python scripts/benchmark_columnar.py --compression-levels
python scripts/benchmark_columnar.py --block-sizes
```
