#!/usr/bin/env python3
"""
ZNA Compression Comparison: Row vs Columnar vs BAM-style

This script compares compression approaches using more realistic DNA patterns.
"""
from __future__ import annotations

import struct
import random
import zstandard
from dataclasses import dataclass, field
from typing import List, Tuple
import os

# --- ENCODING ---
_BASES = ('A', 'C', 'G', 'T')

def make_dna_encode_table() -> bytes:
    table = bytearray([255] * 256)
    for i, base in enumerate(_BASES):
        table[ord(base)] = i
        table[ord(base.lower())] = i
    return bytes(table)

_TRANS_TABLE = make_dna_encode_table()

def encode_sequence(seq: str) -> bytes:
    """Encode DNA sequence to 2-bit packed bytes."""
    seq_bytes = seq.upper().encode('ascii')
    translated = seq_bytes.translate(_TRANS_TABLE)
    
    seq_len = len(seq)
    full_bytes = seq_len // 4
    remainder = seq_len % 4
    total_bytes = full_bytes + (1 if remainder else 0)
    
    result = bytearray(total_bytes)
    
    for i in range(full_bytes):
        offset = i * 4
        result[i] = (
            (translated[offset] << 6) |
            (translated[offset + 1] << 4) |
            (translated[offset + 2] << 2) |
            translated[offset + 3]
        )
    
    if remainder:
        offset = full_bytes * 4
        byte_val = 0
        for j in range(remainder):
            byte_val |= translated[offset + j] << (6 - j * 2)
        result[full_bytes] = byte_val
    
    return bytes(result)


# --- REALISTIC DNA GENERATION ---

def generate_genome_chunk(length: int) -> str:
    """Generate a genome-like chunk with realistic GC content and patterns."""
    # Realistic GC content ~40-60%
    gc_content = random.uniform(0.4, 0.6)
    
    result = []
    for _ in range(length):
        if random.random() < gc_content:
            result.append(random.choice('GC'))
        else:
            result.append(random.choice('AT'))
    
    return ''.join(result)


def generate_reads_from_genome(genome: str, num_reads: int, read_length: int) -> List[str]:
    """Generate reads by sampling from a genome (simulates real sequencing)."""
    reads = []
    genome_len = len(genome)
    
    for _ in range(num_reads):
        start = random.randint(0, genome_len - read_length)
        read = genome[start:start + read_length]
        reads.append(read)
    
    return reads


def add_sequencing_errors(read: str, error_rate: float = 0.001) -> str:
    """Add random sequencing errors."""
    result = list(read)
    for i in range(len(result)):
        if random.random() < error_rate:
            result[i] = random.choice('ACGT')
    return ''.join(result)


# --- ENCODING FORMATS ---

def encode_row_format(sequences: List[str], flags_list: List[int], seq_len_bytes: int = 2) -> bytes:
    """Row-oriented format (current ZNA v1)."""
    buffer = bytearray()
    if seq_len_bytes == 1:
        len_struct = struct.Struct("<B")
    else:
        len_struct = struct.Struct("<H")
    
    for seq, flags in zip(sequences, flags_list):
        buffer.append(flags)
        buffer.extend(len_struct.pack(len(seq)))
        buffer.extend(encode_sequence(seq))
    
    return bytes(buffer)


def encode_columnar_format(sequences: List[str], flags_list: List[int], seq_len_bytes: int = 2) -> bytes:
    """Columnar format (proposed ZNA v2)."""
    flags_buf = bytearray()
    lengths_buf = bytearray()
    seqs_buf = bytearray()
    
    if seq_len_bytes == 1:
        len_struct = struct.Struct("<B")
    else:
        len_struct = struct.Struct("<H")
    
    for seq, flags in zip(sequences, flags_list):
        flags_buf.append(flags)
        lengths_buf.extend(len_struct.pack(len(seq)))
        seqs_buf.extend(encode_sequence(seq))
    
    return bytes(flags_buf) + bytes(lengths_buf) + bytes(seqs_buf)


# --- ADVANCED COMPRESSION TECHNIQUES ---

def encode_with_rle_flags(sequences: List[str], flags_list: List[int]) -> bytes:
    """Columnar with RLE-encoded flags."""
    # Run-length encode flags
    rle_flags = bytearray()
    if flags_list:
        current_flag = flags_list[0]
        count = 1
        for f in flags_list[1:]:
            if f == current_flag and count < 255:
                count += 1
            else:
                rle_flags.append(current_flag)
                rle_flags.append(count)
                current_flag = f
                count = 1
        rle_flags.append(current_flag)
        rle_flags.append(count)
    
    # Regular lengths and sequences
    lengths_buf = bytearray()
    seqs_buf = bytearray()
    len_struct = struct.Struct("<H")
    
    for seq in sequences:
        lengths_buf.extend(len_struct.pack(len(seq)))
        seqs_buf.extend(encode_sequence(seq))
    
    return bytes(rle_flags) + bytes(lengths_buf) + bytes(seqs_buf)


def encode_with_delta_lengths(sequences: List[str], flags_list: List[int]) -> bytes:
    """Columnar with delta-encoded lengths."""
    flags_buf = bytearray(flags_list)
    
    # Delta encode lengths
    lengths = [len(seq) for seq in sequences]
    delta_lengths = bytearray()
    if lengths:
        delta_lengths.extend(struct.pack("<H", lengths[0]))  # First value as-is
        for i in range(1, len(lengths)):
            delta = lengths[i] - lengths[i-1]
            # Store as signed byte if possible, else 2 bytes
            if -128 <= delta <= 127:
                delta_lengths.append(0)  # Marker for 1-byte delta
                delta_lengths.append(delta & 0xFF)
            else:
                delta_lengths.append(1)  # Marker for 2-byte value
                delta_lengths.extend(struct.pack("<H", lengths[i]))
    
    seqs_buf = bytearray()
    for seq in sequences:
        seqs_buf.extend(encode_sequence(seq))
    
    return bytes(flags_buf) + bytes(delta_lengths) + bytes(seqs_buf)


# --- BENCHMARKING ---

def run_realistic_benchmark():
    """Benchmark with realistic genome-derived reads."""
    print("=" * 70)
    print("ZNA Compression: Realistic DNA Patterns")
    print("=" * 70)
    
    # Generate a "genome" chunk to sample from
    genome_size = 1_000_000  # 1MB genome
    print(f"\nGenerating {genome_size:,} bp genome chunk...")
    genome = generate_genome_chunk(genome_size)
    
    compressor_fast = zstandard.ZstdCompressor(level=3)
    compressor_best = zstandard.ZstdCompressor(level=19)
    
    # Test configurations
    configs = [
        (50000, 150, "Illumina 150bp, 50K reads"),
        (50000, 300, "Illumina 300bp, 50K reads"),
        (10000, 1000, "Long reads 1kb, 10K reads"),
    ]
    
    for num_reads, read_length, description in configs:
        print(f"\n{'─' * 70}")
        print(f"Test: {description}")
        print(f"{'─' * 70}")
        
        # Generate reads from genome (simulates real sequencing)
        reads = generate_reads_from_genome(genome, num_reads, read_length)
        reads = [add_sequencing_errors(r, 0.001) for r in reads]  # 0.1% error rate
        
        # All paired-end
        flags_list = [0x05 if i % 2 == 0 else 0x06 for i in range(num_reads)]
        
        # Encode in different formats
        row_data = encode_row_format(reads, flags_list)
        col_data = encode_columnar_format(reads, flags_list)
        
        # Also test just raw sequences (no metadata)
        raw_seqs = bytearray()
        for seq in reads:
            raw_seqs.extend(encode_sequence(seq))
        raw_seqs = bytes(raw_seqs)
        
        # Compress with different levels
        print(f"\nRaw 2-bit encoded size: {len(raw_seqs):,} bytes")
        print(f"With metadata (row):    {len(row_data):,} bytes")
        print(f"With metadata (col):    {len(col_data):,} bytes")
        
        row_fast = len(compressor_fast.compress(row_data))
        col_fast = len(compressor_fast.compress(col_data))
        raw_fast = len(compressor_fast.compress(raw_seqs))
        
        row_best = len(compressor_best.compress(row_data))
        col_best = len(compressor_best.compress(col_data))
        raw_best = len(compressor_best.compress(raw_seqs))
        
        print(f"\nZstd Level 3 (fast):")
        print(f"  Raw sequences:  {raw_fast:>10,} bytes ({len(raw_seqs)/raw_fast:.2f}x)")
        print(f"  Row format:     {row_fast:>10,} bytes ({len(row_data)/row_fast:.2f}x)")
        print(f"  Columnar:       {col_fast:>10,} bytes ({len(col_data)/col_fast:.2f}x)")
        print(f"  ► Columnar vs Row: {(row_fast-col_fast)/row_fast*100:.1f}% smaller")
        
        print(f"\nZstd Level 19 (best):")
        print(f"  Raw sequences:  {raw_best:>10,} bytes ({len(raw_seqs)/raw_best:.2f}x)")
        print(f"  Row format:     {row_best:>10,} bytes ({len(row_data)/row_best:.2f}x)")
        print(f"  Columnar:       {col_best:>10,} bytes ({len(col_data)/col_best:.2f}x)")
        print(f"  ► Columnar vs Row: {(row_best-col_best)/row_best*100:.1f}% smaller")
        
        # Show what would be needed for text FASTA equivalent
        fasta_size = sum(len(r) + 2 for r in reads)  # +2 for header/newline estimate
        print(f"\nComparison (Zstd L3):")
        print(f"  Uncompressed FASTA: ~{fasta_size:>10,} bytes")
        print(f"  ZNA Columnar:        {col_fast:>10,} bytes ({fasta_size/col_fast:.1f}x smaller)")


def run_overlap_analysis():
    """Analyze how much sequences overlap (important for compression)."""
    print("\n" + "=" * 70)
    print("Sequence Overlap Analysis")
    print("=" * 70)
    
    # High overlap scenario (sorted reads from same region)
    genome = generate_genome_chunk(10000)
    
    # Sample overlapping reads
    reads = []
    for i in range(1000):
        start = i * 2  # High overlap
        if start + 150 <= len(genome):
            reads.append(genome[start:start + 150])
    
    # Encode
    seqs_buf = bytearray()
    for seq in reads:
        seqs_buf.extend(encode_sequence(seq))
    seqs_data = bytes(seqs_buf)
    
    compressor = zstandard.ZstdCompressor(level=3)
    compressed = compressor.compress(seqs_data)
    
    print(f"\nHighly overlapping reads (sliding window):")
    print(f"  Reads: {len(reads)}")
    print(f"  Uncompressed: {len(seqs_data):,} bytes")
    print(f"  Compressed:   {len(compressed):,} bytes")
    print(f"  Ratio: {len(seqs_data)/len(compressed):.1f}x")
    
    # Random sampling (low overlap)
    random_reads = generate_reads_from_genome(genome, 1000, 150)
    seqs_buf = bytearray()
    for seq in random_reads:
        seqs_buf.extend(encode_sequence(seq))
    seqs_data = bytes(seqs_buf)
    compressed = compressor.compress(seqs_data)
    
    print(f"\nRandomly sampled reads:")
    print(f"  Uncompressed: {len(seqs_data):,} bytes")
    print(f"  Compressed:   {len(compressed):,} bytes")
    print(f"  Ratio: {len(seqs_data)/len(compressed):.1f}x")


def compare_with_gzip():
    """Compare Zstd with gzip (used by BAM)."""
    import gzip
    import io
    
    print("\n" + "=" * 70)
    print("Compression Algorithm Comparison")
    print("=" * 70)
    
    genome = generate_genome_chunk(500000)
    reads = generate_reads_from_genome(genome, 50000, 150)
    
    # Just sequences (fair comparison)
    seqs_buf = bytearray()
    for seq in reads:
        seqs_buf.extend(encode_sequence(seq))
    data = bytes(seqs_buf)
    
    print(f"\nInput: 50K reads × 150bp = {len(data):,} bytes")
    
    # Gzip (what BAM uses internally)
    gzip_buf = io.BytesIO()
    with gzip.GzipFile(fileobj=gzip_buf, mode='wb', compresslevel=6) as f:
        f.write(data)
    gzip_size = len(gzip_buf.getvalue())
    
    # Zstd levels
    for level in [1, 3, 9, 19]:
        comp = zstandard.ZstdCompressor(level=level)
        zstd_size = len(comp.compress(data))
        print(f"  Zstd L{level:2d}: {zstd_size:>10,} bytes ({len(data)/zstd_size:.2f}x)")
    
    print(f"  Gzip L6:  {gzip_size:>10,} bytes ({len(data)/gzip_size:.2f}x)")


if __name__ == "__main__":
    random.seed(42)
    run_realistic_benchmark()
    run_overlap_analysis()
    compare_with_gzip()
    print("\n" + "=" * 70)
    print("Analysis complete!")
    print("=" * 70)
