################################################################################
# ZNA: compressed nucleic acid data format
################################################################################
from __future__ import annotations

import struct
import zstandard
from dataclasses import dataclass
from enum import IntFlag
from typing import BinaryIO, Iterable, Iterator, Tuple

# Try to import C++ accelerated functions, fall back to pure Python
try:
    from zna._accel import encode_sequence as _accel_encode_sequence
    from zna._accel import decode_block_records as _accel_decode_block_records
    _USE_ACCEL = True
except ImportError:
    _USE_ACCEL = False

def is_accelerated() -> bool:
    """Check if C++ acceleration is available."""
    return _USE_ACCEL

# --- CONSTANTS ---
# File Magic and Version
# 4-byte magic number to identify ZNA files
# fourth byte is 0x1A for 26 (format created in 2026)
_MAGIC = b"ZNA\x1A"
# Use entire byte for versioning to allow future expansion
_VERSION = 1

# Struct Format: Magic(4s), Ver(B), Len(B), Flags(B), CompMethod(B), CompLevel(B), 2xStrLen(H)
# < = Little Endian
_FILE_HEADER_FMT = "<4sBBBBBHH"
_FILE_HEADER_SIZE = struct.calcsize(_FILE_HEADER_FMT)

# Block Header: CompressedSize (I), UncompressedSize (I), RecordCount (I)
_BLOCK_HEADER_FMT = "<III"
_BLOCK_HEADER_SIZE = struct.calcsize(_BLOCK_HEADER_FMT)

# --- FORMAT LIMITS ---
MIN_SEQ_LEN_BYTES = 1
MAX_SEQ_LEN_BYTES = 4
# Derived from unsigned short (H) in struct header
MAX_METADATA_LENGTH = 65535

# --- COMPRESSION TYPES ---
COMPRESSION_NONE = 0
COMPRESSION_ZSTD = 1

# --- FLAGS ---
class ZnaHeaderFlags(IntFlag):
    STRAND_SPECIFIC = 1
    READ1_ANTISENSE = 2  # Read1 is antisense (needs reverse complement)
    READ2_ANTISENSE = 4  # Read2 is antisense (needs reverse complement)

class ZnaRecordFlags(IntFlag):
    IS_READ1 = 1
    IS_READ2 = 2
    IS_PAIRED = 4


# --- REVERSE COMPLEMENT ---
# Complement table: A<->T, C<->G
_COMPLEMENT_TABLE = str.maketrans('ACGTacgt', 'TGCAtgca')

def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(_COMPLEMENT_TABLE)[::-1]


# --- LOOKUP TABLES ---
_BASES = ('A', 'C', 'G', 'T')

def make_dna_encode_table() -> bytes:
    """
    Create translation table for bytes.translate() encoding.
    
    Maps ASCII bytes to 2-bit integers: A→0, C→1, G→2, T→3
    Invalid characters map to 255 (triggers error during bit packing).
    """
    table = bytearray([255] * 256)
    for i, base in enumerate(_BASES):
        table[ord(base)] = i
        table[ord(base.lower())] = i
    return bytes(table)

# Encode Table: Translation table for vectorized encoding
_TRANS_TABLE = make_dna_encode_table()

# Decode Table: Maps byte value to 4-character string
# Each byte represents 4 bases packed as 2-bit values
_DECODE_TABLE = tuple(
    _BASES[(val >> 6) & 0b11] + _BASES[(val >> 4) & 0b11] + 
    _BASES[(val >> 2) & 0b11] + _BASES[(val >> 0) & 0b11]
    for val in range(256)
)

# Pre-compiled struct formats for record length (only powers of 2 supported)
_SEQ_LEN_STRUCTS = {
    1: struct.Struct("<B"),
    2: struct.Struct("<H"),
    4: struct.Struct("<I"),
}
_VALID_SEQ_LEN_BYTES = frozenset(_SEQ_LEN_STRUCTS.keys())

# Default compression level for zstd (1-22, lower = faster)
DEFAULT_ZSTD_LEVEL = 3
# Default block size for compression (128KB for optimal Zstd performance)
DEFAULT_BLOCK_SIZE = 131072


@dataclass(slots=True)
class ZnaHeader:
    """Header metadata for ZNA files."""
    read_group: str
    description: str = ""
    seq_len_bytes: int = MIN_SEQ_LEN_BYTES
    strand_specific: bool = False
    read1_antisense: bool = False  # True if read1 is antisense strand
    read2_antisense: bool = False  # True if read2 is antisense strand
    compression_method: int = COMPRESSION_NONE
    compression_level: int = DEFAULT_ZSTD_LEVEL

    def __post_init__(self):
        if self.seq_len_bytes not in _VALID_SEQ_LEN_BYTES:
            raise ValueError(f"seq_len_bytes must be one of {_VALID_SEQ_LEN_BYTES}, got {self.seq_len_bytes}")
        if len(self.read_group.encode('utf-8')) > MAX_METADATA_LENGTH:
            raise ValueError(f"read_group too long: {len(self.read_group)} > {MAX_METADATA_LENGTH}")
        if len(self.description.encode('utf-8')) > MAX_METADATA_LENGTH:
            raise ValueError(f"description too long: {len(self.description)} > {MAX_METADATA_LENGTH}")
        if self.compression_method not in (COMPRESSION_NONE, COMPRESSION_ZSTD):
            raise ValueError(f"compression_method must be 0 or 1, got {self.compression_method}")
        if not (1 <= self.compression_level <= 22):
            raise ValueError(f"compression_level must be 1-22, got {self.compression_level}")


class ZnaWriter:
    """Writer for ZNA files, handling block-based compression."""
    def __init__(self, fh: BinaryIO, header: ZnaHeader, block_size: int = DEFAULT_BLOCK_SIZE):
        self._fh = fh
        self._header = header
        self._seq_len_bytes = header.seq_len_bytes
        self._max_len = (1 << (8 * header.seq_len_bytes)) - 1
        self._block_size = block_size
        
        # Pre-size buffer to avoid reallocations
        self._buffer = bytearray(block_size)
        self._buffer_pos = 0
        self._records_in_block = 0
        
        # Pre-compiled struct for sequence length
        self._seq_len_struct = _SEQ_LEN_STRUCTS.get(header.seq_len_bytes)
        
        # Reuse compressor instance
        if header.compression_method == COMPRESSION_ZSTD:
            self._compressor = zstandard.ZstdCompressor(level=header.compression_level)
        else:
            self._compressor = None
            
        self._write_header()

    def __enter__(self) -> ZnaWriter:
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self.close()

    def close(self) -> None:
        self._flush_block()    

    def _write_header(self) -> None:
        """Write file header with metadata."""
        rg_bytes = self._header.read_group.encode("utf-8")
        desc_bytes = self._header.description.encode("utf-8")
        
        flags = 0
        if self._header.strand_specific:
            flags |= ZnaHeaderFlags.STRAND_SPECIFIC
        if self._header.read1_antisense:
            flags |= ZnaHeaderFlags.READ1_ANTISENSE
        if self._header.read2_antisense:
            flags |= ZnaHeaderFlags.READ2_ANTISENSE
        
        header_fixed = struct.pack(
            _FILE_HEADER_FMT,
            _MAGIC,
            _VERSION,
            self._header.seq_len_bytes,
            flags,
            self._header.compression_method,
            self._header.compression_level,
            len(rg_bytes),
            len(desc_bytes)
        )
        
        self._fh.write(header_fixed)
        self._fh.write(rg_bytes)
        self._fh.write(desc_bytes)

    def write_record(self, seq: str, is_paired: bool, is_read1: bool, is_read2: bool) -> None:
        """Write a single sequence record to the buffer.
        
        For strand-specific libraries, antisense reads are automatically
        reverse complemented to normalize all reads to the sense strand.
        """
        # Strand normalization: reverse complement antisense reads
        if self._header.strand_specific:
            if is_read1 and self._header.read1_antisense:
                seq = reverse_complement(seq)
            elif is_read2 and self._header.read2_antisense:
                seq = reverse_complement(seq)
        
        seq_len = len(seq)
        if seq_len > self._max_len:
            raise ValueError(
                f"Sequence length {seq_len} exceeds maximum {self._max_len} "
                f"allowed by header (seq_len_bytes={self._seq_len_bytes})"
            )
        
        flags = 0
        if is_read1:
            flags |= ZnaRecordFlags.IS_READ1
        if is_read2:
            flags |= ZnaRecordFlags.IS_READ2
        if is_paired:
            flags |= ZnaRecordFlags.IS_PAIRED
        
        try:
            encoded_seq = _encode_sequence(seq)
        except ValueError as e:
            raise ValueError(f"Invalid characters in sequence (only A,C,G,T allowed). Violation in: {seq[:20]}...") from e
        
        # Calculate record size
        encoded_len = len(encoded_seq)
        record_size = 1 + self._seq_len_bytes + encoded_len
        
        # Flush before adding if record won't fit (avoids buffer growth in normal cases)
        new_pos = self._buffer_pos + record_size
        if new_pos > len(self._buffer):
            # Flush current block first
            self._flush_block()
            new_pos = record_size
            
            # If single record exceeds block_size, grow buffer (rare edge case)
            if new_pos > len(self._buffer):
                self._buffer = bytearray(new_pos)
        
        # Write flag byte
        pos = self._buffer_pos
        self._buffer[pos] = flags
        pos += 1
        
        # Write sequence length (always power-of-2, no None check needed)
        self._buffer[pos:pos + self._seq_len_bytes] = self._seq_len_struct.pack(seq_len)
        pos += self._seq_len_bytes
        
        # Write encoded sequence
        self._buffer[pos:pos + encoded_len] = encoded_seq
        self._buffer_pos = new_pos
        
        self._records_in_block += 1
        
        # Check if block limit reached
        if self._buffer_pos >= self._block_size:
            self._flush_block()
            
    def _flush_block(self) -> None:
        """Compress and write current block to file."""
        if self._buffer_pos == 0:
            return

        # Get only the used portion of the buffer
        raw_data = memoryview(self._buffer)[:self._buffer_pos]
        uncompressed_size = self._buffer_pos
        count = self._records_in_block
        
        # Compress if needed (reuse compressor instance)
        if self._compressor is not None:
            final_data = self._compressor.compress(bytes(raw_data))
        else:
            final_data = bytes(raw_data)
            
        compressed_size = len(final_data)
        
        # Write block header and data
        block_header = struct.pack(_BLOCK_HEADER_FMT, compressed_size, uncompressed_size, count)
        
        self._fh.write(block_header)
        self._fh.write(final_data)
        
        # Reset buffer position (reuse pre-allocated buffer)
        self._buffer_pos = 0
        self._records_in_block = 0


class ZnaReader:
    """Reader for ZNA files, handling block-based decompression."""
    def __init__(self, fh: BinaryIO):
        self._fh = fh
        self._header = self._read_header()

    @property
    def header(self) -> ZnaHeader:
        return self._header
    
    def __iter__(self) -> Iterator[Tuple[str, bool, bool, bool]]:
        return self.records()

    def _read_exact(self, n: int) -> bytes:
        data = self._fh.read(n)
        if len(data) != n:
            raise EOFError("Unexpected EOF while reading zna")
        return data

    def _read_header(self) -> ZnaHeader:
        # struct.unpack to parse header (clarity)
        fixed_data = self._read_exact(_FILE_HEADER_SIZE)
        fixed_fields = struct.unpack(_FILE_HEADER_FMT, fixed_data)
        (magic, ver, len_bytes, flags, compression_method, compression_level, rg_len, desc_len) = fixed_fields
        
        if magic != _MAGIC:
            raise ValueError("Not a ZNA file")
        if ver != _VERSION:
            raise ValueError(f"Unsupported ZNA version: {ver}")
            
        read_group = self._read_exact(rg_len).decode("utf-8")
        description = self._read_exact(desc_len).decode("utf-8")
        
        return ZnaHeader(
            read_group=read_group,
            description=description,
            seq_len_bytes=len_bytes,
            strand_specific=bool(flags & ZnaHeaderFlags.STRAND_SPECIFIC),
            read1_antisense=bool(flags & ZnaHeaderFlags.READ1_ANTISENSE),
            read2_antisense=bool(flags & ZnaHeaderFlags.READ2_ANTISENSE),
            compression_method=compression_method,
            compression_level=compression_level
        )

    def records(self, restore_strand: bool = False) -> Iterator[Tuple[str, bool, bool, bool]]:
        """Iterate over all records in the file.
        
        Args:
            restore_strand: If True and library is strand-specific, restore
                           original strand by reverse complementing antisense reads.
        """
        read_exact = self._read_exact
        fh_read = self._fh.read
        len_bytes = self._header.seq_len_bytes
        compression_method = self._header.compression_method
        use_accel = _USE_ACCEL
        
        # Strand restoration: reverse complement to recover original antisense reads
        do_restore_r1 = restore_strand and self._header.strand_specific and self._header.read1_antisense
        do_restore_r2 = restore_strand and self._header.strand_specific and self._header.read2_antisense

        if compression_method == COMPRESSION_ZSTD:
            dctx = zstandard.ZstdDecompressor()
        
        # Pure Python decode table (only used if C++ extension not available)
        if not use_accel:
            decode_table = _DECODE_TABLE
        
        while True:
            # 1. Read Block Header
            block_header_data = fh_read(_BLOCK_HEADER_SIZE)
            if not block_header_data:
                break
            if len(block_header_data) < _BLOCK_HEADER_SIZE:
                 raise EOFError(f"Incomplete block header. Expected {_BLOCK_HEADER_SIZE} bytes, got {len(block_header_data)}")
                 
            comp_size, uncomp_size, count = struct.unpack(_BLOCK_HEADER_FMT, block_header_data)
            
            # 2. Read Block Data
            block_payload = read_exact(comp_size)
            
            # 3. Decompress if needed
            if compression_method == COMPRESSION_ZSTD:
                block_data = dctx.decompress(block_payload, max_output_size=uncomp_size)
            else:
                block_data = block_payload
            
            # 4. Decode records - use C++ extension if available
            if use_accel:
                # C++ decodes entire block at once
                for seq, is_paired, is_read1, is_read2 in _accel_decode_block_records(block_data, len_bytes, count):
                    # Restore original strand if requested
                    if do_restore_r1 and is_read1:
                        seq = reverse_complement(seq)
                    elif do_restore_r2 and is_read2:
                        seq = reverse_complement(seq)
                    yield seq, is_paired, is_read1, is_read2
            else:
                # Pure Python fallback
                mv = memoryview(block_data)
                offset = 0
                for _ in range(count):
                    flags = mv[offset]
                    offset += 1
                    
                    is_read1 = bool(flags & ZnaRecordFlags.IS_READ1)
                    is_read2 = bool(flags & ZnaRecordFlags.IS_READ2)
                    is_paired = bool(flags & ZnaRecordFlags.IS_PAIRED)
                    
                    seq_len = int.from_bytes(mv[offset:offset + len_bytes], "little")
                    offset += len_bytes
                    
                    enc_len = (seq_len + 3) >> 2
                    seq_bytes = mv[offset:offset + enc_len]
                    offset += enc_len
                    
                    # Decode: Build list of 4-mers, join once, then trim
                    chunks = [decode_table[b] for b in seq_bytes]
                    full_seq = ''.join(chunks)
                    seq = full_seq[:seq_len]
                    
                    # Restore original strand if requested
                    if do_restore_r1 and is_read1:
                        seq = reverse_complement(seq)
                    elif do_restore_r2 and is_read2:
                        seq = reverse_complement(seq)
                    
                    yield seq, is_paired, is_read1, is_read2


def _encode_sequence(seq: str) -> bytes:
    """
    Encode DNA sequence to 2-bit packed bytes.
    
    Each base (A, C, G, T) is encoded as 2 bits (00, 01, 10, 11).
    Four bases are packed into each byte.
    
    Args:
        seq: DNA sequence string (A, C, G, T)
        
    Returns:
        Packed bytes with 4 bases per byte
        
    Raises:
        ValueError: If sequence contains invalid characters
    """
    # Use C++ implementation if available
    if _USE_ACCEL:
        return _accel_encode_sequence(seq)
    
    # Pure Python fallback with vectorized character lookup
    bseq = seq.encode('ascii')
    mapped = bseq.translate(_TRANS_TABLE)
    length = len(mapped)
    
    # Pre-allocate output buffer
    out = bytearray((length + 3) // 4)
    
    i = 0
    idx = 0
    full_chunks = length // 4
    
    # Pack full 4-base chunks into bytes
    while idx < full_chunks:
        val = (
            (mapped[i] << 6) | 
            (mapped[i+1] << 4) | 
            (mapped[i+2] << 2) | 
            mapped[i+3]
        )
        out[idx] = val
        idx += 1
        i += 4

    # Tail: Handle remaining 1-3 bases
    if i < length:
        val = 0
        shift = 6
        while i < length:
            val |= (mapped[i] << shift)
            i += 1
            shift -= 2
        out[idx] = val

    return bytes(out)


def write_zna(
    fh: BinaryIO,
    header: ZnaHeader,
    records: Iterable[Tuple[str, bool, bool, bool]],
) -> None:
    """Write ZNA file with header and records."""
    with ZnaWriter(fh, header) as writer:
        for seq, is_paired, is_read1, is_read2 in records:
            writer.write_record(seq, is_paired, is_read1, is_read2)


def read_zna(fh: BinaryIO) -> Tuple[ZnaHeader, Iterator[Tuple[str, bool, bool, bool]]]:
    """Read ZNA file header and return records iterator."""
    reader = ZnaReader(fh)
    return reader.header, reader.records()