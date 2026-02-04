#!/usr/bin/env python3
"""
Comprehensive benchmarks for columnar ZNA format.

Tests compression ratios, encoding/decoding speed, and profiling
with focus on short read (Illumina) data optimization.
"""
import sys
import time
import tempfile
import os
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple

# Test data paths
R1_PATH = "/Users/mkiyer/University of Michigan Dropbox/MED-mctp-iyerlab/projects/transcriptor/models/chr22/chr22_test_R1.fq.gz"
R2_PATH = "/Users/mkiyer/University of Michigan Dropbox/MED-mctp-iyerlab/projects/transcriptor/models/chr22/chr22_test_R2.fq.gz"


def get_uncompressed_size(fastq_gz_path: str) -> int:
    """Get uncompressed FASTQ size in bytes."""
    result = subprocess.run(
        ["gzcat", fastq_gz_path],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    return len(result.stdout)


def run_encode(r1: str, r2: str, output: str, **kwargs) -> Tuple[float, int]:
    """
    Run ZNA encode and return (duration, output_size).
    
    kwargs can include:
    - sorted: bool
    - level: int
    - uncompressed: bool
    """
    cmd = ["zna", "encode", r1, r2, "-o", output, "--quiet"]
    
    if kwargs.get("sorted"):
        cmd.append("--sorted")
    
    if kwargs.get("uncompressed"):
        cmd.append("--uncompressed")
    elif "level" in kwargs:
        cmd.extend(["--level", str(kwargs["level"])])
    
    start = time.time()
    result = subprocess.run(cmd, capture_output=True)
    duration = time.time() - start
    
    if result.returncode != 0:
        print(f"ERROR: Encode failed: {result.stderr.decode()}", file=sys.stderr)
        sys.exit(1)
    
    output_size = os.path.getsize(output)
    return duration, output_size


def run_decode(input_path: str, output: str, **kwargs) -> float:
    """
    Run ZNA decode and return duration.
    
    kwargs can include:
    - no_reshuffle: bool
    """
    cmd = ["zna", "decode", input_path, "-o", output, "--quiet"]
    
    if kwargs.get("no_reshuffle"):
        cmd.append("--no-reshuffle")
    
    start = time.time()
    result = subprocess.run(cmd, capture_output=True)
    duration = time.time() - start
    
    if result.returncode != 0:
        print(f"ERROR: Decode failed: {result.stderr.decode()}", file=sys.stderr)
        sys.exit(1)
    
    return duration


def format_size(bytes_val: int) -> str:
    """Format bytes as human-readable."""
    for unit in ['B', 'KB', 'MB', 'GB']:
        if bytes_val < 1024.0:
            return f"{bytes_val:.2f} {unit}"
        bytes_val /= 1024.0
    return f"{bytes_val:.2f} TB"


def format_speed(bytes_val: int, seconds: float) -> str:
    """Format throughput as MB/s."""
    mb_per_sec = bytes_val / seconds / (1024 * 1024)
    return f"{mb_per_sec:.1f} MB/s"


def benchmark_compression_levels(r1: str, r2: str, input_size: int):
    """Benchmark different compression levels."""
    print("\n" + "="*80)
    print("COMPRESSION LEVEL ANALYSIS")
    print("="*80)
    print(f"\nInput: {format_size(input_size)} (uncompressed FASTQ)")
    print(f"\n{'Level':<8} {'Encode':<12} {'Decode':<12} {'Size':<12} {'Ratio':<8} {'Notes'}")
    print("-" * 80)
    
    with tempfile.TemporaryDirectory() as tmpdir:
        results = []
        
        for level in [1, 3, 5, 9, 15, 19]:
            zna_path = f"{tmpdir}/test_L{level}.zna"
            out_path = f"{tmpdir}/decoded_L{level}.fa"
            
            # Encode
            encode_time, zna_size = run_encode(r1, r2, zna_path, level=level)
            
            # Decode
            decode_time = run_decode(zna_path, out_path)
            
            ratio = input_size / zna_size
            encode_speed = format_speed(input_size, encode_time)
            decode_speed = format_speed(input_size, decode_time)
            
            notes = ""
            if level == 3:
                notes = "← Default"
            
            print(f"{level:<8} {encode_speed:<12} {decode_speed:<12} "
                  f"{format_size(zna_size):<12} {ratio:.2f}x{' ':<4} {notes}")
            
            results.append({
                'level': level,
                'encode_time': encode_time,
                'decode_time': decode_time,
                'size': zna_size,
                'ratio': ratio
            })
        
        # Find best options
        best_encode = min(results, key=lambda x: x['encode_time'])
        best_decode = min(results, key=lambda x: x['decode_time'])
        best_ratio = max(results, key=lambda x: x['ratio'])
        
        print(f"\nBest encode speed: Level {best_encode['level']}")
        print(f"Best decode speed: Level {best_decode['level']}")
        print(f"Best compression: Level {best_ratio['level']} ({best_ratio['ratio']:.2f}x)")


def benchmark_sorted_mode(r1: str, r2: str, input_size: int):
    """Benchmark sorted vs unsorted encoding."""
    print("\n" + "="*80)
    print("SORTED MODE ANALYSIS")
    print("="*80)
    print(f"\nInput: {format_size(input_size)} (uncompressed FASTQ)")
    print(f"\n{'Mode':<15} {'Encode':<12} {'Decode':<12} {'Size':<12} {'Ratio':<8} {'Improvement'}")
    print("-" * 80)
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Unsorted (baseline)
        unsorted_path = f"{tmpdir}/unsorted.zna"
        unsorted_out = f"{tmpdir}/unsorted.fa"
        
        unsorted_enc_time, unsorted_size = run_encode(r1, r2, unsorted_path)
        unsorted_dec_time = run_decode(unsorted_path, unsorted_out)
        unsorted_ratio = input_size / unsorted_size
        
        print(f"{'Unsorted':<15} {format_speed(input_size, unsorted_enc_time):<12} "
              f"{format_speed(input_size, unsorted_dec_time):<12} "
              f"{format_size(unsorted_size):<12} {unsorted_ratio:.2f}x{' ':<4} (baseline)")
        
        # Sorted
        sorted_path = f"{tmpdir}/sorted.zna"
        sorted_out = f"{tmpdir}/sorted.fa"
        
        sorted_enc_time, sorted_size = run_encode(r1, r2, sorted_path, sorted=True)
        sorted_dec_time = run_decode(sorted_path, sorted_out)
        sorted_ratio = input_size / sorted_size
        
        improvement = (unsorted_size - sorted_size) / unsorted_size * 100
        
        print(f"{'Sorted':<15} {format_speed(input_size, sorted_enc_time):<12} "
              f"{format_speed(input_size, sorted_dec_time):<12} "
              f"{format_size(sorted_size):<12} {sorted_ratio:.2f}x{' ':<4} "
              f"{improvement:+.1f}%")
        
        # Sorted + no reshuffle (fastest decode)
        sorted_dec_no_reshuffle = run_decode(sorted_path, sorted_out, no_reshuffle=True)
        
        print(f"{'Sorted (no-RS)':<15} {format_speed(input_size, sorted_enc_time):<12} "
              f"{format_speed(input_size, sorted_dec_no_reshuffle):<12} "
              f"{format_size(sorted_size):<12} {sorted_ratio:.2f}x{' ':<4} "
              f"(faster decode)")
        
        print(f"\nCompression improvement from sorting: {improvement:+.1f}%")
        print(f"Decode speedup without reshuffle: "
              f"{sorted_dec_time / sorted_dec_no_reshuffle:.2f}x")


def benchmark_block_sizes(r1: str, r2: str, input_size: int):
    """Benchmark different block sizes."""
    print("\n" + "="*80)
    print("BLOCK SIZE ANALYSIS")
    print("="*80)
    print(f"\nInput: {format_size(input_size)} (uncompressed FASTQ)")
    print(f"\n{'Block Size':<12} {'Encode':<12} {'Decode':<12} {'Size':<12} {'Ratio':<8} {'Notes'}")
    print("-" * 80)
    
    with tempfile.TemporaryDirectory() as tmpdir:
        results = []
        
        for block_size in [32768, 65536, 131072, 262144, 524288]:
            zna_path = f"{tmpdir}/test_B{block_size}.zna"
            out_path = f"{tmpdir}/decoded_B{block_size}.fa"
            
            # Encode with custom block size
            cmd = ["zna", "encode", r1, r2, "-o", zna_path, 
                   "--block-size", str(block_size), "--quiet"]
            
            start = time.time()
            subprocess.run(cmd, capture_output=True, check=True)
            encode_time = time.time() - start
            zna_size = os.path.getsize(zna_path)
            
            # Decode
            decode_time = run_decode(zna_path, out_path)
            
            ratio = input_size / zna_size
            encode_speed = format_speed(input_size, encode_time)
            decode_speed = format_speed(input_size, decode_time)
            
            block_kb = block_size // 1024
            notes = ""
            if block_size == 131072:
                notes = "← Default"
            
            print(f"{block_kb} KB{' ':<7} {encode_speed:<12} {decode_speed:<12} "
                  f"{format_size(zna_size):<12} {ratio:.2f}x{' ':<4} {notes}")
            
            results.append({
                'block_size': block_size,
                'size': zna_size,
                'ratio': ratio
            })
        
        best = max(results, key=lambda x: x['ratio'])
        print(f"\nOptimal block size: {best['block_size'] // 1024} KB "
              f"({best['ratio']:.2f}x compression)")


def benchmark_uncompressed(r1: str, r2: str, input_size: int):
    """Benchmark uncompressed ZNA format."""
    print("\n" + "="*80)
    print("UNCOMPRESSED FORMAT ANALYSIS")
    print("="*80)
    print(f"\nInput: {format_size(input_size)} (uncompressed FASTQ)")
    print(f"\n{'Format':<15} {'Encode':<12} {'Decode':<12} {'Size':<12} {'Ratio':<8}")
    print("-" * 80)
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Compressed (default)
        comp_path = f"{tmpdir}/compressed.zna"
        comp_out = f"{tmpdir}/comp_decoded.fa"
        
        comp_enc_time, comp_size = run_encode(r1, r2, comp_path)
        comp_dec_time = run_decode(comp_path, comp_out)
        comp_ratio = input_size / comp_size
        
        print(f"{'Compressed':<15} {format_speed(input_size, comp_enc_time):<12} "
              f"{format_speed(input_size, comp_dec_time):<12} "
              f"{format_size(comp_size):<12} {comp_ratio:.2f}x")
        
        # Uncompressed
        uncomp_path = f"{tmpdir}/uncompressed.zna"
        uncomp_out = f"{tmpdir}/uncomp_decoded.fa"
        
        uncomp_enc_time, uncomp_size = run_encode(r1, r2, uncomp_path, uncompressed=True)
        uncomp_dec_time = run_decode(uncomp_path, uncomp_out)
        uncomp_ratio = input_size / uncomp_size
        
        print(f"{'Uncompressed':<15} {format_speed(input_size, uncomp_enc_time):<12} "
              f"{format_speed(input_size, uncomp_dec_time):<12} "
              f"{format_size(uncomp_size):<12} {uncomp_ratio:.2f}x")
        
        # Calculate overhead
        overhead = (uncomp_size - (input_size / 2)) / (input_size / 2) * 100
        print(f"\nZNA overhead (uncompressed): {overhead:.1f}% of raw 2-bit encoding")
        print(f"Compression speedup: {uncomp_dec_time / comp_dec_time:.2f}x decode, "
              f"{uncomp_enc_time / comp_enc_time:.2f}x encode")


def profile_encode(r1: str, r2: str, n_reads: int = 100000):
    """Profile encode operation with cProfile."""
    print("\n" + "="*80)
    print(f"PROFILING ENCODE ({n_reads} reads)")
    print("="*80)
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Extract subset of reads
        subset_r1 = f"{tmpdir}/subset_R1.fq"
        subset_r2 = f"{tmpdir}/subset_R2.fq"
        
        print(f"\nExtracting {n_reads} reads for profiling...")
        subprocess.run(
            f'gzcat "{r1}" | head -n {n_reads * 4} > "{subset_r1}"',
            shell=True, check=True
        )
        subprocess.run(
            f'gzcat "{r2}" | head -n {n_reads * 4} > "{subset_r2}"',
            shell=True, check=True
        )
        
        zna_out = f"{tmpdir}/profile.zna"
        prof_out = f"{tmpdir}/encode.prof"
        
        # Run with profiling
        cmd = [
            "python", "-m", "cProfile", "-o", prof_out,
            "-m", "zna.cli", "encode",
            subset_r1, subset_r2, "-o", zna_out, "--quiet"
        ]
        
        print("Running profiler...")
        subprocess.run(cmd, check=True)
        
        # Analyze profile
        print("\nTop 20 time-consuming functions:")
        print("-" * 80)
        
        subprocess.run([
            "python", "-c",
            f"""import pstats
p = pstats.Stats('{prof_out}')
p.strip_dirs()
p.sort_stats('cumulative')
p.print_stats(20)
"""
        ])


def profile_decode(r1: str, r2: str, n_reads: int = 100000):
    """Profile decode operation with cProfile."""
    print("\n" + "="*80)
    print(f"PROFILING DECODE ({n_reads} reads)")
    print("="*80)
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Extract subset and encode
        subset_r1 = f"{tmpdir}/subset_R1.fq"
        subset_r2 = f"{tmpdir}/subset_R2.fq"
        zna_path = f"{tmpdir}/profile.zna"
        
        print(f"\nPreparing {n_reads} reads for profiling...")
        subprocess.run(
            f'gzcat "{r1}" | head -n {n_reads * 4} > "{subset_r1}"',
            shell=True, check=True
        )
        subprocess.run(
            f'gzcat "{r2}" | head -n {n_reads * 4} > "{subset_r2}"',
            shell=True, check=True
        )
        
        # Encode
        run_encode(subset_r1, subset_r2, zna_path)
        
        # Profile decode
        prof_out = f"{tmpdir}/decode.prof"
        decoded_out = f"{tmpdir}/decoded.fa"
        
        cmd = [
            "python", "-m", "cProfile", "-o", prof_out,
            "-m", "zna.cli", "decode",
            zna_path, "-o", decoded_out, "--quiet"
        ]
        
        print("Running profiler...")
        subprocess.run(cmd, check=True)
        
        # Analyze profile
        print("\nTop 20 time-consuming functions:")
        print("-" * 80)
        
        subprocess.run([
            "python", "-c",
            f"""import pstats
p = pstats.Stats('{prof_out}')
p.strip_dirs()
p.sort_stats('cumulative')
p.print_stats(20)
"""
        ])


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description="Benchmark columnar ZNA format")
    parser.add_argument("--compression-levels", action="store_true",
                       help="Test different compression levels")
    parser.add_argument("--sorted", action="store_true",
                       help="Test sorted vs unsorted mode")
    parser.add_argument("--block-sizes", action="store_true",
                       help="Test different block sizes")
    parser.add_argument("--uncompressed", action="store_true",
                       help="Compare compressed vs uncompressed")
    parser.add_argument("--profile-encode", action="store_true",
                       help="Profile encode operation")
    parser.add_argument("--profile-decode", action="store_true",
                       help="Profile decode operation")
    parser.add_argument("--all", action="store_true",
                       help="Run all benchmarks (except profiling)")
    parser.add_argument("--profile-reads", type=int, default=100000,
                       help="Number of reads for profiling (default: 100000)")
    
    args = parser.parse_args()
    
    # Check if test data exists
    if not os.path.exists(R1_PATH):
        print(f"ERROR: Test data not found: {R1_PATH}", file=sys.stderr)
        sys.exit(1)
    
    print("ZNA Columnar Format Benchmark")
    print("=" * 80)
    print(f"Test Data: chr22_test (1M paired-end reads, 150bp)")
    print(f"Date: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Calculate input size
    print("\nCalculating input size...")
    input_size = get_uncompressed_size(R1_PATH) + get_uncompressed_size(R2_PATH)
    print(f"Total uncompressed FASTQ size: {format_size(input_size)}")
    
    # Run requested benchmarks
    run_any = False
    
    if args.all or args.compression_levels:
        benchmark_compression_levels(R1_PATH, R2_PATH, input_size)
        run_any = True
    
    if args.all or args.sorted:
        benchmark_sorted_mode(R1_PATH, R2_PATH, input_size)
        run_any = True
    
    if args.all or args.block_sizes:
        benchmark_block_sizes(R1_PATH, R2_PATH, input_size)
        run_any = True
    
    if args.all or args.uncompressed:
        benchmark_uncompressed(R1_PATH, R2_PATH, input_size)
        run_any = True
    
    if args.profile_encode:
        profile_encode(R1_PATH, R2_PATH, args.profile_reads)
        run_any = True
    
    if args.profile_decode:
        profile_decode(R1_PATH, R2_PATH, args.profile_reads)
        run_any = True
    
    if not run_any:
        print("\nNo benchmarks selected. Use --help to see options.")
        print("Quick start: python benchmark_columnar.py --all")
    
    print("\n" + "="*80)
    print("Benchmark complete!")
    print("="*80)


if __name__ == "__main__":
    main()
