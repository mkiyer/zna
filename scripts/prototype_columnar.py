#!/usr/bin/env python3
"""
ZNA Columnar Format Prototype

This script tests the columnar storage approach to validate compression improvements.
Run with: python scripts/prototype_columnar.py
"""
from __future__ import annotations

import struct
import random
import zstandard
from dataclasses import dataclass, field
from typing import List, Tuple
import time

# --- ENCODING (from core.py) ---
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
    idx = 0
    
    for i in range(full_bytes):
        offset = i * 4
        result[i] = (
            (translated[offset] << 6) |
            (translated[offset + 1] << 4) |
            (translated[offset + 2] << 2) |
            translated[offset + 3]
        )
        idx = i + 1
    
    if remainder:
        offset = full_bytes * 4
        byte_val = 0
        for j in range(remainder):
            byte_val |= translated[offset + j] << (6 - j * 2)
        result[idx] = byte_val
    
    return bytes(result)


def generate_random_sequence(length: int) -> str:
    """Generate random DNA sequence."""
    return ''.join(random.choice('ACGT') for _ in range(length))


def generate_test_data(num_records: int, seq_length: int, paired: bool = True) -> List[Tuple[str, int, int]]:
    """Generate test records: (sequence, flags, length)"""
    records = []
    for i in range(num_records):
        seq = generate_random_sequence(seq_length)
        if paired:
            flags = 0x05 if i % 2 == 0 else 0x06  # IS_READ1|IS_PAIRED or IS_READ2|IS_PAIRED
        else:
            flags = 0x00
        records.append((seq, flags, seq_length))
    return records


# --- ROW-ORIENTED FORMAT (Current V1) ---

def encode_row_format(records: List[Tuple[str, int, int]], seq_len_bytes: int = 2) -> bytes:
    """Encode records in row-oriented format (current ZNA v1)."""
    buffer = bytearray()
    if seq_len_bytes == 1:
        len_struct = struct.Struct("<B")
    elif seq_len_bytes == 2:
        len_struct = struct.Struct("<H")
    else:
        len_struct = struct.Struct("<I")
    
    for seq, flags, seq_len in records:
        # [Flag][Length][EncodedSeq]
        buffer.append(flags)
        buffer.extend(len_struct.pack(seq_len))
        buffer.extend(encode_sequence(seq))
    
    return bytes(buffer)


# --- COLUMNAR FORMAT (Proposed V2) ---

@dataclass
class ColumnarBlock:
    """Columnar block buffer."""
    flags: bytearray = field(default_factory=bytearray)
    lengths: bytearray = field(default_factory=bytearray)
    sequences: bytearray = field(default_factory=bytearray)
    count: int = 0
    seq_len_bytes: int = 2
    
    def __post_init__(self):
        if self.seq_len_bytes == 1:
            self._len_struct = struct.Struct("<B")
        elif self.seq_len_bytes == 2:
            self._len_struct = struct.Struct("<H")
        else:
            self._len_struct = struct.Struct("<I")
    
    def add(self, seq: str, flags: int, seq_len: int) -> None:
        self.flags.append(flags)
        self.lengths.extend(self._len_struct.pack(seq_len))
        self.sequences.extend(encode_sequence(seq))
        self.count += 1
    
    def to_bytes(self) -> bytes:
        """Return columnar payload."""
        return bytes(self.flags) + bytes(self.lengths) + bytes(self.sequences)
    
    @property
    def flags_size(self) -> int:
        return len(self.flags)
    
    @property
    def lengths_size(self) -> int:
        return len(self.lengths)
    
    @property
    def sequences_size(self) -> int:
        return len(self.sequences)


def encode_columnar_format(records: List[Tuple[str, int, int]], seq_len_bytes: int = 2) -> bytes:
    """Encode records in columnar format (proposed ZNA v2)."""
    block = ColumnarBlock(seq_len_bytes=seq_len_bytes)
    
    for seq, flags, seq_len in records:
        block.add(seq, flags, seq_len)
    
    return block.to_bytes()


# --- BENCHMARKING ---

def benchmark_compression(name: str, data: bytes, compressor: zstandard.ZstdCompressor) -> Tuple[int, int, float]:
    """Compress data and return (original_size, compressed_size, ratio)."""
    compressed = compressor.compress(data)
    original = len(data)
    comp_size = len(compressed)
    ratio = original / comp_size if comp_size > 0 else float('inf')
    return original, comp_size, ratio


def run_benchmark():
    """Run compression benchmark comparing row vs columnar formats."""
    print("=" * 70)
    print("ZNA Columnar Format Compression Benchmark")
    print("=" * 70)
    
    # Test configurations
    configs = [
        (10000, 150, True, "Illumina PE 150bp"),
        (10000, 150, False, "Illumina SE 150bp"),
        (10000, 300, True, "Illumina PE 300bp"),
        (5000, 1000, False, "Long reads 1kb"),
        (1000, 10000, False, "Very long reads 10kb"),
    ]
    
    compressor = zstandard.ZstdCompressor(level=3)
    
    for num_records, seq_length, paired, description in configs:
        print(f"\n{'─' * 70}")
        print(f"Test: {description}")
        print(f"Records: {num_records:,}, Seq Length: {seq_length}bp, Paired: {paired}")
        print(f"{'─' * 70}")
        
        # Generate test data
        records = generate_test_data(num_records, seq_length, paired)
        
        # Encode in both formats
        row_data = encode_row_format(records)
        col_data = encode_columnar_format(records)
        
        # Compress
        row_orig, row_comp, row_ratio = benchmark_compression("Row", row_data, compressor)
        col_orig, col_comp, col_ratio = benchmark_compression("Columnar", col_data, compressor)
        
        # Breakdown for columnar
        block = ColumnarBlock(seq_len_bytes=2)
        for seq, flags, seq_len in records:
            block.add(seq, flags, seq_len)
        
        flags_comp = len(compressor.compress(bytes(block.flags)))
        lengths_comp = len(compressor.compress(bytes(block.lengths)))
        seqs_comp = len(compressor.compress(bytes(block.sequences)))
        
        # Results
        print(f"\nUncompressed sizes:")
        print(f"  Row format:      {row_orig:>12,} bytes")
        print(f"  Columnar format: {col_orig:>12,} bytes")
        
        print(f"\nCompressed sizes (Zstd level 3):")
        print(f"  Row format:      {row_comp:>12,} bytes (ratio: {row_ratio:.2f}x)")
        print(f"  Columnar format: {col_comp:>12,} bytes (ratio: {col_ratio:.2f}x)")
        
        improvement = (row_comp - col_comp) / row_comp * 100 if row_comp > 0 else 0
        print(f"\n  ► Columnar is {improvement:.1f}% smaller")
        
        print(f"\nColumnar stream breakdown (compressed):")
        print(f"  Flags:     {block.flags_size:>8,} bytes → {flags_comp:>8,} bytes ({block.flags_size/flags_comp:.0f}x)")
        print(f"  Lengths:   {block.lengths_size:>8,} bytes → {lengths_comp:>8,} bytes ({block.lengths_size/lengths_comp:.0f}x)")
        print(f"  Sequences: {block.sequences_size:>8,} bytes → {seqs_comp:>8,} bytes ({block.sequences_size/seqs_comp:.2f}x)")
        
        # Metadata overhead analysis
        metadata_orig = block.flags_size + block.lengths_size
        metadata_comp = flags_comp + lengths_comp
        print(f"\n  Metadata overhead: {metadata_orig:,} bytes → {metadata_comp:,} bytes ({metadata_orig/metadata_comp:.0f}x compression)")


def run_realistic_test():
    """Test with more realistic data patterns."""
    print("\n" + "=" * 70)
    print("Realistic Data Patterns Test")
    print("=" * 70)
    
    compressor = zstandard.ZstdCompressor(level=3)
    
    # Simulate 1M reads worth of blocks
    num_blocks = 100
    records_per_block = 10000
    seq_length = 150
    
    total_row_comp = 0
    total_col_comp = 0
    total_orig = 0
    
    print(f"\nSimulating {num_blocks} blocks × {records_per_block} records = {num_blocks * records_per_block:,} total reads")
    print(f"Sequence length: {seq_length}bp\n")
    
    for block_idx in range(num_blocks):
        records = generate_test_data(records_per_block, seq_length, paired=True)
        
        row_data = encode_row_format(records)
        col_data = encode_columnar_format(records)
        
        row_comp = len(compressor.compress(row_data))
        col_comp = len(compressor.compress(col_data))
        
        total_row_comp += row_comp
        total_col_comp += col_comp
        total_orig += len(row_data)
        
        if (block_idx + 1) % 20 == 0:
            print(f"  Processed {block_idx + 1} blocks...")
    
    print(f"\nTotal Results:")
    print(f"  Original size:     {total_orig / 1024 / 1024:.2f} MB")
    print(f"  Row compressed:    {total_row_comp / 1024 / 1024:.2f} MB")
    print(f"  Columnar compressed: {total_col_comp / 1024 / 1024:.2f} MB")
    
    savings = (total_row_comp - total_col_comp) / 1024 / 1024
    pct = (total_row_comp - total_col_comp) / total_row_comp * 100
    print(f"\n  ► Columnar saves {savings:.2f} MB ({pct:.1f}% smaller)")


if __name__ == "__main__":
    random.seed(42)  # Reproducibility
    run_benchmark()
    run_realistic_test()
    print("\n" + "=" * 70)
    print("Benchmark complete!")
    print("=" * 70)
