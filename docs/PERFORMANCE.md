# ZNA Performance

Benchmarks on Apple Silicon (M-series), Python 3.12, March 2026.

Test dataset: 25.4M paired-end Illumina reads (150 bp), interleaved FASTQ
from a simulated HeLa transcriptome (minimap2 aligned), 10.76 GB uncompressed.

## Summary

| Metric | Value |
|--------|-------|
| **Encode** | 135 MB/s (312 K rec/s) |
| **Decode** | 626 MB/s (1.44 M rec/s) |
| **Compression** | 16.5× (default) |

## Compression Comparison

| Format | Size | Ratio |
|--------|------|-------|
| FASTQ (uncompressed) | 10.76 GB | 1.0× |
| FASTQ.gz | 1.15 GB | 9.4× |
| ZNA (uncompressed) | 1.05 GB | 10.2× |
| **ZNA (ZSTD L9, 4 MB blocks)** | **666.8 MB** | **16.5×** |
| ZNA (ZSTD L9 + 10 labels) | 718.7 MB | 15.3× |

## Compression Level Tuning

| Level | Encode (MB/s) | Decode (MB/s) | Size | Ratio | Notes |
|-------|---------------|---------------|------|-------|-------|
| 1 | 150 | 651 | 716.0 MB | 15.4× | Fastest encode |
| 3 | 146 | 683 | 696.4 MB | 15.8× | |
| 5 | 142 | 662 | 674.0 MB | 16.4× | |
| **9** | **137** | **635** | **666.8 MB** | **16.5×** | **Default** |
| 15 | 97 | 635 | 664.1 MB | 16.6× | Best compression |

Level 9 is the default.  Level 15 gains < 1% compression at 30% slower
encode.  Level 1 is 10% faster to encode with 7% larger files.

## Block Size Tuning

| Block Size | Encode (MB/s) | Decode (MB/s) | Size | Ratio |
|------------|---------------|---------------|------|-------|
| 512 KB | 136 | 644 | 680.7 MB | 16.2× |
| 1 MB | 135 | 667 | 675.7 MB | 16.3× |
| **4 MB** | **137** | **614** | **666.8 MB** | **16.5×** |
| 8 MB | 133 | 576 | 664.6 MB | 16.6× |

4 MB blocks (default) balance compression and decode throughput.
Smaller blocks decode slightly faster but compress less.

## Labeled Encode/Decode (10 SAM tags)

When storing per-sequence labels (NM, ms, AS, nn, tp, cm, s1, s2, de, rl):

| Metric | Plain | Labeled (10 tags) | Overhead |
|--------|-------|-------------------|----------|
| Encode | 135 MB/s (312 K rec/s) | 134 MB/s (309 K rec/s) | ~1% |
| Decode | 626 MB/s (1.44 M rec/s) | 159 MB/s (367 K rec/s) | Labels emitted as text |
| Size | 666.8 MB | 718.7 MB | +7.8% |

Label encode overhead is negligible thanks to C++ acceleration.  Decode is
slower when `--labels` is used because SAM-style tag strings are formatted
and written to the output; decode without `--labels` runs at full speed.

## Why It's Fast

The columnar format groups homogeneous data streams together, giving ZSTD
much better pattern matching:
- **Flags stream**: paired/single-end patterns compress 500–1000×
- **Lengths stream**: uniform 150 bp reads compress ~1000×
- **Sequences stream**: 2-bit DNA data compresses 3–5×
- **Label columns**: contiguous numeric arrays compress efficiently

Remaining hot-path bottleneck is input I/O (~75% of encode time is
reading gzipped FASTQ).

## Reproducing Benchmarks

```bash
python scripts/benchmark_perf.py
```
