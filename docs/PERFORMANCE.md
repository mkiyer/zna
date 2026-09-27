# ZNA Performance

> **These figures predate ZNA 0.3.4 and understate current performance.**
> They were measured in March 2026 against the 0.3.x codec, before the 0.3.4
> rewrite roughly doubled decode and cut encode substantially. The dataset is a
> real 10.76 GB library that is not reproducible here, so the numbers are left
> as recorded rather than guessed at; the compression figures still hold. Format
> version 3 changed only *which* record lands in which block — a block now ends
> where a fragment does, so it overruns its size target by at most one record —
> and left the block payload's layout untouched.
>
> For current throughput see the table in [`../README.md`](../README.md) — its
> "Re-measured at 0.5.1" note carries fresh 150 bp figures — and
> [`../CHANGELOG.md`](../CHANGELOG.md) for the per-path deltas. The merger's numbers
> are current as of 0.6.0: see [Overlap merging](#overlap-merging-060) at the end.

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

## Overlap merging (0.6.0)

Measured for 0.6.0 on an Apple M3 Max (clang -O3), against a 0.5.3 wheel built from its
tag. For the kernel rows both versions were compiled into one driver binary, each with
its own weights (identical at `e = 0.01`); runs are on a quiet machine and a repeat agreed
within ~2%.

### Kernel, µs/pair

| input | contract | scan 0.5.3 → 0.6 | `process_pair` 0.5.3 → 0.6 | `merge_chunk` 0.5.3 → 0.6 |
|---|---|---|---|---|
| chr22, 1M pairs | `--adapter-trimmed` | 0.455 → 0.301 (−34%) | 0.667 → 0.330 (−50%) | 0.795 → 0.458 (−42%) |
| chr22, 1M pairs | undeclared | 0.443 → 0.376 (−15%) | 0.667 → 0.509 (−24%) | 0.792 → 0.633 (−20%) |
| truth panel, trimmed benches | `--adapter-trimmed` | 0.559 → 0.340 (−39%) | 0.819 → 0.359 (−56%) | 0.929 → 0.471 (−49%) |
| truth panel, raw benches | undeclared | 0.519 → 0.435 (−16%) | 0.790 → 0.611 (−23%) | 0.873 → 0.689 (−21%) |

Where it comes from, with instrumented counts on 1M chr22 pairs:

- **The floor.** 0.5.3 started the scan's incumbent at 8 bits, 0.6 at the pair's floor
  (~28 bits). Work falls only slightly — 191.7 → 187.5 shifts and 10,277 → 10,104
  compared bases per pair — yet the scan is 15% faster (0.448 → 0.379): the shifts
  pruned are mostly short, scalar, branchy tail shifts (an inference; the counts are
  measured). At the same floor, 0.6's undeclared scan does exactly 0.5.3's work, with 0
  result mismatches and the same time, so the undeclared gain is the floor alone.
- **The declaration** halves the work: 94.2 shifts and 5,096 bases per pair.
- **The contract is a template parameter.** A runtime `bool` in the flank loop measured
  28% slower undeclared and 12% slower declared (non-inlined wrappers); inlined, it made
  the undeclared scan 6% slower than 0.5.3's despite the higher floor.
- **The `--adapter-trimmed` check** runs one extra unrestricted scan on the first
  100,000 pairs of a declared run: 0.703 µs/pair against 0.5.3's 0.667 (+5.4%) for those
  pairs only, about 20 ms per run. With it on, a declared chr22 chunk is 0.463 µs/pair.

### End to end

`zna merge --threads 1`, plain FASTQ in, three interleaved repetitions, medians, CPU
seconds (qualification run, `docs/MERGE_BENCHMARK_RESULTS.md` §9):

| input | 0.5.3 | 0.6 undeclared | 0.6 `--adapter-trimmed` |
|---|---:|---:|---:|
| chr22, 1M pairs | 1.08 | 0.94 (−13.0%) | 0.77 (−28.7%) |
| tx-dev-clean, 500k | 0.56 | 0.50 (−10.7%) | 0.46 (−17.9%) |
| tx-dev-noisy, 500k | 0.60 | 0.54 (−10.0%) | — |
| tx-dev-fastp, 500k | 0.56 | 0.51 (−8.9%) | 0.45 (−19.6%) |
| hg38, 1M | 1.06 | 0.90 (−15.1%) | — |

On 1M chr22 pairs with 5 repetitions: `--threads 1` 1.10 s (0.5.3), 0.96 s (0.6
undeclared), 0.78 s (declared), at 19.7G, 17.9G and 13.5G instructions; `--threads 4`
is ~0.50 s wall for all three, bound by I/O and Python (CPU 1.45, 1.28, 1.07 s);
`encode --merge-pairs`, which is single-threaded, 2.37, 2.23 and 2.06 s.

### Memory

Peak footprint (macOS, dirty memory) on 1M chr22 pairs is 26–34 MB at 1, 4 and 16
threads, against 0.5.3's 27–33 MB, with byte-identical output at every thread count.
The first 0.6 port measured 30 / 69 / 194 MB at 1 / 4 / 16 threads: its submit window
held each chunk's input buffers until the chunk was written, and the reader makes a new
~2 MB buffer about once per chunk. Only a chunk that stopped for table capacity keeps its
input now (`docs/METHODS.md` §2.8). `encode --merge-pairs` peaks at ~60 MB in both versions.

### Long reads

The policy tables are built once per run for the longest read seen, and grown by
doubling. A 2×150 library never regrows past the initial 256-base capacity. One pair of
long reads among 2,000 chr22 pairs:

| read length | 0.5.3 | 0.6 as first ported | 0.6.0 |
|---|---:|---:|---:|
| 5 kb | 0.09 s | 7.12 s | 0.76 s |
| 10 kb | 0.09 s | 43.9 s | 1.66 s |

The port rebuilt each `dfit` tail from `bⁿ` down, 6.7× per doubling; 0.6.0 advances
the exact tail one step per overlap length (`binom_cap`, `params.py`), capacity 16,384 in
0.41 s. What remains is the floor table, one 50-digit `Decimal` logarithm per entry,
~1.2 s to capacity 16,384 — on the main thread, once per run. `docs/ROADMAP.md` records
the shortcut that would remove it. The scan itself is O(L²) per pair and takes
milliseconds for such a read.

### Sanitizers

The core (`merge_core.hpp`, `fastq_chunk.hpp`) was run under AddressSanitizer and
UBSan with every policy table allocated at exactly its capacity, so an index one past it
lands in a redzone: 445,015 real pairs (truth panel and chr22) through `process_pair` with
tables at capacity 150 (every 150 bp read at the edge), 372,101 fuzzed pairs over capacities 1–1,024, and
34,599 pairs through the chunk path with mid-batch capacity stops. Clean.
`scripts/merge_bench/asan_scan.cpp` is the standing form of that check.

