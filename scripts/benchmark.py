#!/usr/bin/env python3
"""
ZNA Performance Benchmarking Suite
Run with: python benchmark.py
"""
import sys
import time
import random
import tempfile
import statistics
from pathlib import Path
from io import BytesIO
from typing import List, Tuple

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from zna.core import (
    ZnaWriter, ZnaReader, ZnaHeader,
    COMPRESSION_NONE, COMPRESSION_ZSTD,
    DEFAULT_ZSTD_LEVEL
)


# --- Test Data Generation ---

def generate_random_sequence(length: int) -> str:
    """Generate random DNA sequence."""
    return ''.join(random.choice('ACGT') for _ in range(length))


def generate_test_sequences(count: int, min_len: int, max_len: int) -> List[str]:
    """Generate test sequences with varying lengths."""
    return [generate_random_sequence(random.randint(min_len, max_len)) 
            for _ in range(count)]


def generate_paired_sequences(count: int, read_len: int) -> List[Tuple[str, bool, bool, bool]]:
    """Generate paired-end reads."""
    records = []
    for _ in range(count):
        r1 = generate_random_sequence(read_len)
        r2 = generate_random_sequence(read_len)
        records.append((r1, True, True, False))
        records.append((r2, True, False, True))
    return records


# --- Benchmark Functions ---

def benchmark_encode(sequences: List[str], 
                    compression_method: int = COMPRESSION_NONE,
                    block_size: int = 131072,
                    iterations: int = 3) -> dict:
    """Benchmark encoding performance."""
    header = ZnaHeader(
        read_group="benchmark",
        seq_len_bytes=2,
        compression_method=compression_method,
        compression_level=DEFAULT_ZSTD_LEVEL
    )
    
    times = []
    sizes = []
    
    for _ in range(iterations):
        output = BytesIO()
        
        start = time.perf_counter()
        with ZnaWriter(output, header, block_size=block_size) as writer:
            for seq in sequences:
                writer.write_record(seq, False, False, False)
        elapsed = time.perf_counter() - start
        
        times.append(elapsed)
        sizes.append(output.tell())
    
    total_bases = sum(len(s) for s in sequences)
    
    return {
        'time_mean': statistics.mean(times),
        'time_stdev': statistics.stdev(times) if len(times) > 1 else 0,
        'throughput_mb_s': (total_bases / statistics.mean(times)) / (1024 * 1024),
        'throughput_records_s': len(sequences) / statistics.mean(times),
        'file_size': statistics.mean(sizes),
        'compression_ratio': total_bases / statistics.mean(sizes) if statistics.mean(sizes) > 0 else 0
    }


def benchmark_decode(sequences: List[str],
                    compression_method: int = COMPRESSION_NONE,
                    block_size: int = 131072,
                    iterations: int = 3) -> dict:
    """Benchmark decoding performance."""
    # First encode
    header = ZnaHeader(
        read_group="benchmark",
        seq_len_bytes=2,
        compression_method=compression_method,
        compression_level=DEFAULT_ZSTD_LEVEL
    )
    
    encoded = BytesIO()
    with ZnaWriter(encoded, header, block_size=block_size) as writer:
        for seq in sequences:
            writer.write_record(seq, False, False, False)
    
    encoded_data = encoded.getvalue()
    total_bases = sum(len(s) for s in sequences)
    
    # Now benchmark decoding
    times = []
    
    for _ in range(iterations):
        input_buf = BytesIO(encoded_data)
        
        start = time.perf_counter()
        reader = ZnaReader(input_buf)
        decoded = list(reader.records())
        elapsed = time.perf_counter() - start
        
        times.append(elapsed)
    
    return {
        'time_mean': statistics.mean(times),
        'time_stdev': statistics.stdev(times) if len(times) > 1 else 0,
        'throughput_mb_s': (total_bases / statistics.mean(times)) / (1024 * 1024),
        'throughput_records_s': len(sequences) / statistics.mean(times),
        'file_size': len(encoded_data)
    }


def benchmark_roundtrip(sequences: List[str],
                       compression_method: int = COMPRESSION_NONE,
                       block_size: int = 131072,
                       iterations: int = 3) -> dict:
    """Benchmark full encode + decode cycle."""
    times = []
    
    for _ in range(iterations):
        # Encode
        header = ZnaHeader(
            read_group="benchmark",
            seq_len_bytes=2,
            compression_method=compression_method,
            compression_level=DEFAULT_ZSTD_LEVEL
        )
        
        output = BytesIO()
        start = time.perf_counter()
        
        with ZnaWriter(output, header, block_size=block_size) as writer:
            for seq in sequences:
                writer.write_record(seq, False, False, False)
        
        # Decode
        output.seek(0)
        reader = ZnaReader(output)
        decoded = list(reader.records())
        
        elapsed = time.perf_counter() - start
        times.append(elapsed)
    
    total_bases = sum(len(s) for s in sequences)
    
    return {
        'time_mean': statistics.mean(times),
        'time_stdev': statistics.stdev(times) if len(times) > 1 else 0,
        'throughput_mb_s': (total_bases / statistics.mean(times)) / (1024 * 1024),
        'throughput_records_s': len(sequences) / statistics.mean(times)
    }


# --- Test Suites ---

def run_workload_benchmarks():
    """Run benchmarks for different workload types."""
    print("=" * 80)
    print("WORKLOAD BENCHMARKS")
    print("=" * 80)
    
    workloads = [
        ("Short Reads (Illumina)", 10000, 100, 150),
        ("Medium Reads", 5000, 300, 500),
        ("Long Reads (PacBio)", 1000, 1000, 5000),
        ("Very Long Reads (Nanopore)", 500, 5000, 15000),
    ]
    
    for name, count, min_len, max_len in workloads:
        print(f"\n{name}: {count} sequences, {min_len}-{max_len} bp")
        sequences = generate_test_sequences(count, min_len, max_len)
        total_mb = sum(len(s) for s in sequences) / (1024 * 1024)
        print(f"  Total data: {total_mb:.2f} MB")
        
        # Uncompressed
        print("  Uncompressed:")
        enc_result = benchmark_encode(sequences, COMPRESSION_NONE)
        print(f"    Encode: {enc_result['throughput_mb_s']:.1f} MB/s, {enc_result['throughput_records_s']:.0f} rec/s")
        
        dec_result = benchmark_decode(sequences, COMPRESSION_NONE)
        print(f"    Decode: {dec_result['throughput_mb_s']:.1f} MB/s, {dec_result['throughput_records_s']:.0f} rec/s")
        
        # Compressed
        print("  Compressed (ZSTD):")
        enc_result = benchmark_encode(sequences, COMPRESSION_ZSTD)
        print(f"    Encode: {enc_result['throughput_mb_s']:.1f} MB/s, {enc_result['throughput_records_s']:.0f} rec/s")
        print(f"    Compression: {enc_result['compression_ratio']:.2f}x")
        
        dec_result = benchmark_decode(sequences, COMPRESSION_ZSTD)
        print(f"    Decode: {dec_result['throughput_mb_s']:.1f} MB/s, {dec_result['throughput_records_s']:.0f} rec/s")


def run_block_size_benchmarks():
    """Benchmark different block sizes."""
    print("\n" + "=" * 80)
    print("BLOCK SIZE BENCHMARKS")
    print("=" * 80)
    
    sequences = generate_test_sequences(5000, 100, 150)
    block_sizes = [32768, 65536, 131072, 262144, 524288, 1048576]  # 32KB to 1MB
    
    print("\nBlock Size | Encode (MB/s) | Decode (MB/s) | File Size | Comp Ratio")
    print("-" * 80)
    
    for block_size in block_sizes:
        enc_result = benchmark_encode(sequences, COMPRESSION_ZSTD, block_size=block_size, iterations=5)
        dec_result = benchmark_decode(sequences, COMPRESSION_ZSTD, block_size=block_size, iterations=5)
        
        print(f"{block_size:>9d} | {enc_result['throughput_mb_s']:>13.1f} | "
              f"{dec_result['throughput_mb_s']:>13.1f} | {enc_result['file_size']:>9d} | "
              f"{enc_result['compression_ratio']:>10.2f}x")


def run_compression_level_benchmarks():
    """Benchmark different compression levels."""
    print("\n" + "=" * 80)
    print("COMPRESSION LEVEL BENCHMARKS")
    print("=" * 80)
    
    sequences = generate_test_sequences(5000, 100, 150)
    
    print("\nLevel | Encode (MB/s) | Decode (MB/s) | File Size | Comp Ratio")
    print("-" * 80)
    
    for level in [1, 3, 5, 9, 15, 19]:
        header = ZnaHeader(
            read_group="benchmark",
            seq_len_bytes=2,
            compression_method=COMPRESSION_ZSTD,
            compression_level=level
        )
        
        output = BytesIO()
        start = time.perf_counter()
        with ZnaWriter(output, header) as writer:
            for seq in sequences:
                writer.write_record(seq, False, False, False)
        enc_time = time.perf_counter() - start
        
        file_size = output.tell()
        total_bases = sum(len(s) for s in sequences)
        
        # Decode
        output.seek(0)
        start = time.perf_counter()
        reader = ZnaReader(output)
        list(reader.records())
        dec_time = time.perf_counter() - start
        
        enc_throughput = (total_bases / enc_time) / (1024 * 1024)
        dec_throughput = (total_bases / dec_time) / (1024 * 1024)
        comp_ratio = total_bases / file_size
        
        print(f"{level:>5d} | {enc_throughput:>13.1f} | {dec_throughput:>13.1f} | "
              f"{file_size:>9d} | {comp_ratio:>10.2f}x")


def run_quick_benchmark():
    """Quick benchmark for development."""
    print("=" * 80)
    print("QUICK BENCHMARK")
    print("=" * 80)
    
    sequences = generate_test_sequences(10000, 100, 150)
    total_mb = sum(len(s) for s in sequences) / (1024 * 1024)
    print(f"\nTest data: {len(sequences)} sequences, {total_mb:.2f} MB")
    
    print("\nUncompressed:")
    result = benchmark_roundtrip(sequences, COMPRESSION_NONE, iterations=5)
    print(f"  Roundtrip: {result['throughput_mb_s']:.1f} MB/s, {result['throughput_records_s']:.0f} rec/s")
    
    print("\nCompressed (ZSTD level 3):")
    result = benchmark_roundtrip(sequences, COMPRESSION_ZSTD, iterations=5)
    print(f"  Roundtrip: {result['throughput_mb_s']:.1f} MB/s, {result['throughput_records_s']:.0f} rec/s")


# --- Main ---

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description="ZNA Performance Benchmarks")
    parser.add_argument("--quick", action="store_true", help="Run quick benchmark")
    parser.add_argument("--workloads", action="store_true", help="Run workload benchmarks")
    parser.add_argument("--blocks", action="store_true", help="Run block size benchmarks")
    parser.add_argument("--compression", action="store_true", help="Run compression level benchmarks")
    parser.add_argument("--all", action="store_true", help="Run all benchmarks")
    
    args = parser.parse_args()
    
    # Set random seed for reproducibility
    random.seed(42)
    
    if args.quick or not any([args.workloads, args.blocks, args.compression, args.all]):
        run_quick_benchmark()
    
    if args.workloads or args.all:
        run_workload_benchmarks()
    
    if args.blocks or args.all:
        run_block_size_benchmarks()
    
    if args.compression or args.all:
        run_compression_level_benchmarks()
    
    print("\n" + "=" * 80)
    print("Benchmark complete!")
    print("=" * 80)


if __name__ == "__main__":
    main()
