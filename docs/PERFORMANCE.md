# ZNA Performance Report

## Executive Summary

ZNA achieves **135.7 MB/s roundtrip throughput** with C++ acceleration, representing a **9.5x speedup** over the pure Python baseline. For long reads, performance reaches **2.8+ GB/s** for both encoding and decoding, making it suitable for high-throughput genomic applications.

| Metric | Baseline (Python) | Optimized (C++) | Speedup |
|--------|-------------------|-----------------|---------|
| **Roundtrip** | 14.3 MB/s | 135.7 MB/s | **9.5x** |
| **Encode (Short)** | 24 MB/s | 189.5 MB/s | **7.9x** |
| **Decode (Short)** | 42 MB/s | 668.8 MB/s | **15.9x** |
| **Encode (Long)** | 25 MB/s | 2,824.7 MB/s | **113x** |
| **Decode (Long)** | 262 MB/s | 3,392.7 MB/s | **13x** |

---

## Current Performance Metrics

### Workload Performance

| Workload | Size | Encode (MB/s) | Decode (MB/s) | Compression |
|----------|------|---------------|---------------|-------------|
| **Short Reads (Illumina)** | 1.19 MB | 189.5 | 668.8 | 3.68x |
| **Medium Reads** | 1.91 MB | 540.5 | 1,280.9 | 3.87x |
| **Long Reads (PacBio)** | 2.79 MB | 1,921.5 | 2,864.6 | 3.98x |
| **Very Long Reads (Nanopore)** | 4.77 MB | 2,824.7 | 3,392.7 | 3.99x |

### Block Size Analysis

| Block Size | Encode (MB/s) | Decode (MB/s) | File Size | Ratio |
|------------|---------------|---------------|-----------|-------|
| 32 KB | 173.1 | 631.0 | 169,841 | 3.68x |
| 64 KB | 170.3 | 677.9 | 169,742 | 3.68x |
| **128 KB** | **181.6** | **685.2** | **169,740** | **3.68x** ✅ |
| 256 KB | 155.0 | 815.6 | 173,093 | 3.61x |
| 512 KB | 144.0 | 671.8 | 173,093 | 3.61x |
| 1 MB | 158.7 | 776.4 | 173,093 | 3.61x |

**Optimal**: 128 KB (current default)
- Best compression ratio (3.68x)
- Balanced encode/decode performance
- Optimal for streaming workflows

### Compression Level Analysis

| Level | Encode (MB/s) | Decode (MB/s) | File Size | Ratio |
|-------|---------------|---------------|-----------|-------|
| 1 | 138.0 | 539.8 | 169,740 | 3.68x |
| **3** | **166.8** | **670.0** | **169,740** | **3.68x** ✅ |
| 5 | 146.5 | 602.8 | 169,763 | 3.68x |
| 9 | 164.8 | 545.7 | 169,763 | 3.68x |
| 15 | 81.3 | 647.7 | 173,112 | 3.61x |
| 19 | 52.8 | 572.1 | 169,612 | 3.68x |

**Optimal**: Level 3 (current default)
- Best throughput vs compression trade-off
- Minimal file size difference across levels 1-9
- Decode performance remains excellent

---

## Key Performance Insights

1. **Read Length Scaling**: Performance improves dramatically with longer reads
   - Short (100-150bp): 189.5 MB/s encode, 668.8 MB/s decode
   - Medium (300-500bp): 540.5 MB/s encode, 1,280.9 MB/s decode  
   - Long (1-5kb): 1,921.5 MB/s encode, 2,864.6 MB/s decode
   - Very long (5-15kb): **2.8 GB/s encode, 3.4 GB/s decode**

2. **Compression Characteristics**:
   - Consistent 3.7-4.0x compression across all workloads
   - Minimal performance penalty vs uncompressed (5-10%)
   - Level 3 provides optimal speed/compression balance

3. **Memory Efficiency**:
   - Block-based processing (128 KB default)
   - Streaming-friendly architecture
   - Minimal memory overhead

4. **CPU Utilization**:
   - C++ hot loops minimize Python interpreter overhead
   - Block-level decoding reduces per-record costs
   - Near-optimal CPU cache utilization

---

## Optimization History

### Phase 1: Pure Python Optimizations
1. **List comprehension for decode** (7% improvement)
2. **bytes.translate() for encode** (5% improvement)
3. **Combined**: 12% improvement (14.3 → 16.0 MB/s)

### Phase 2: nanobind C++ Extension ✅
- Moved encode/decode hot loops to C++
- Block-level decoding (processes entire blocks in C++)
- Optimized string allocation with reserve()
- **Result**: 9.5x overall speedup, up to 113x for long read encoding

### Recent Optimizations
- Added `constexpr` to lookup tables
- Improved string memory allocation (reserve + resize)
- Optimized decode loop (removed redundant length checks)
- Added `noexcept` specifiers for better optimization

---

## Technical Implementation

### C++ Extension
- **Framework**: nanobind (lightweight Python bindings)
- **Build System**: scikit-build-core + CMake
- **Optimization**: -O3 with modern C++17
- **Fallback**: Pure Python implementation if C++ unavailable

### Algorithmic Optimizations
- 2-bit base encoding (4 bases per byte)
- Lookup table-based encoding/decoding
- Block-based compression for I/O efficiency
- Reused compressor instances

### Memory Management
- Pre-allocated buffers
- Minimal copying
- String reserve/resize pattern
- Aligned lookup tables (64-byte alignment)

---

## Benchmark Reproduction

```bash
# Quick benchmark (~10 seconds)
python scripts/benchmark.py --quick

# Workload analysis (~30 seconds)  
python scripts/benchmark.py --workloads

# Block size testing (~1 minute)
python scripts/benchmark.py --blocks

# Compression levels (~1 minute)
python scripts/benchmark.py --compression

# Run all benchmarks (~3 minutes)
python scripts/benchmark.py --all
```

### System Information
Benchmarks run on: Apple M-series (ARM64) with Python 3.12
- C++ compiled with -O3 optimization
- nanobind 2.x with stable ABI
- zstandard 0.23.0

---

## Future Optimization Opportunities

### High Priority
1. **Parallel block compression** - Expected 2-4x on multi-core systems
2. **Buffered I/O** - Expected 10-20% for many small records
3. **Batch record writing** - Reduce Python→C++ call overhead

### Medium Priority
1. **SIMD vectorization** - AVX2/NEON for encoding (potential 2-10x)
2. **Memory-mapped I/O** - For large file processing
3. **Custom zstd dictionaries** - Trained on genomic data

### Future Exploration
1. **GPU acceleration** - For massive parallel processing
2. **Indexed format** - Random access support
3. **Adaptive block sizing** - Based on record length distribution
python benchmark.py --compression

# Everything
python benchmark.py --all
```

## System Info
- Date: February 3, 2026
- Python: 3.12.9
- Platform: macOS (Apple Silicon)
