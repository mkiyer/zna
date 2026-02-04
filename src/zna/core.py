################################################################################
# ZNA: compressed nucleic acid data format
#
# Uses columnar block storage for optimal compression:
# - Flags stream: All flag bytes concatenated (compresses ~500x)
# - Lengths stream: All sequence lengths concatenated (compresses ~1000x)
# - Sequences stream: All encoded sequences concatenated
################################################################################
from __future__ import annotations

import struct
import zstandard
from dataclasses import dataclass, field
from enum import IntFlag
from typing import BinaryIO, Iterable, Iterator, Tuple

# Try to import C++ accelerated functions, fall back to pure Python
try:
    from zna._accel import encode_sequence as _accel_encode_sequence
    from zna._accel import reverse_complement as _accel_reverse_complement
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

# Columnar Block Header: CompressedSize(I), UncompressedSize(I), RecordCount(I), FlagsSize(I), LengthsSize(I)
_BLOCK_HEADER_FMT = "<IIIII"
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
    READ1_ANTISENSE = 2   # Read1 is antisense (needs reverse complement)
    READ2_ANTISENSE = 4   # Read2 is antisense (needs reverse complement)

class ZnaRecordFlags(IntFlag):
    IS_READ1 = 1
    IS_READ2 = 2
    IS_PAIRED = 4


# --- REVERSE COMPLEMENT ---
# Complement table: A<->T, C<->G
_COMPLEMENT_TABLE = str.maketrans('ACGTacgt', 'TGCAtgca')

def _py_reverse_complement(seq: str) -> str:
    """Pure Python reverse complement of a DNA sequence."""
    return seq.translate(_COMPLEMENT_TABLE)[::-1]

def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence.
    
    Uses C++ acceleration when available for better performance.
    """
    if _USE_ACCEL:
        return _accel_reverse_complement(seq)
    return _py_reverse_complement(seq)


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


@dataclass
class ZnaBlockBuffer:
    """
    Columnar block buffer for ZNA records.
    
    Accumulates records into separate streams for optimal compression:
    - flags: All flag bytes (highly repetitive, compresses ~500x)
    - lengths: All sequence lengths (highly repetitive, compresses ~1000x)
    - sequences: All encoded sequences (raw DNA entropy)
    
    Optimized for minimal per-record overhead using direct byte operations.
    """
    seq_len_bytes: int = 2
    
    # Columnar streams - pre-initialized for append efficiency
    flags: bytearray = field(default_factory=bytearray)
    lengths: bytearray = field(default_factory=bytearray)
    sequences: bytearray = field(default_factory=bytearray)
    
    def __post_init__(self):
        # Pre-cache the pack method for the struct (avoids lookup overhead)
        if self.seq_len_bytes == 1:
            self._pack_len = struct.Struct("<B").pack
        elif self.seq_len_bytes == 2:
            self._pack_len = struct.Struct("<H").pack
        else:
            self._pack_len = struct.Struct("<I").pack
        # Cache append/extend methods (small but measurable speedup)
        self._flags_append = self.flags.append
        self._lengths_extend = self.lengths.extend
        self._sequences_extend = self.sequences.extend
    
    def add(self, encoded_seq: bytes, flags: int, seq_len: int) -> None:
        """Add a record to the buffer (optimized hot path)."""
        self._flags_append(flags)
        self._lengths_extend(self._pack_len(seq_len))
        self._sequences_extend(encoded_seq)
    
    def to_bytes(self) -> bytes:
        """Return columnar payload: [flags][lengths][sequences]"""
        return bytes(self.flags) + bytes(self.lengths) + bytes(self.sequences)
    
    def clear(self) -> None:
        """Clear all buffers for reuse."""
        self.flags.clear()
        self.lengths.clear()
        self.sequences.clear()
    
    @property
    def count(self) -> int:
        """Number of records in buffer."""
        return len(self.flags)
    
    @property
    def size(self) -> int:
        """Approximate size of buffer in bytes."""
        return len(self.flags) + len(self.lengths) + len(self.sequences)
    
    @property
    def flags_size(self) -> int:
        return len(self.flags)
    
    @property
    def lengths_size(self) -> int:
        return len(self.lengths)


class ZnaWriter:
    """
    Writer for ZNA files using columnar block storage.
    
    Records are accumulated into a ZnaBlockBuffer and flushed as columnar
    blocks for optimal compression. The block format is:
    
    Block Header (20 bytes):
        - Compressed size (4 bytes)
        - Uncompressed size (4 bytes)
        - Record count (4 bytes)
        - Flags stream size (4 bytes)
        - Lengths stream size (4 bytes)
    
    Block Payload:
        [Flags Stream][Lengths Stream][Sequences Stream]
    """
    __slots__ = (
        '_fh', '_header', '_seq_len_bytes', '_max_len', '_block_size', 
        '_npolicy', '_buffer', '_compressor', '_do_strand_norm_r1', 
        '_do_strand_norm_r2', '_buffer_add', '_buffer_size'
    )
    
    def __init__(
        self, 
        fh: BinaryIO, 
        header: ZnaHeader, 
        block_size: int = DEFAULT_BLOCK_SIZE, 
        npolicy: str = None
    ):
        self._fh = fh
        self._header = header
        self._seq_len_bytes = header.seq_len_bytes
        self._max_len = (1 << (8 * header.seq_len_bytes)) - 1
        self._block_size = block_size
        self._npolicy = npolicy
        
        # Pre-compute strand normalization flags (avoid per-record checks)
        self._do_strand_norm_r1 = header.strand_specific and header.read1_antisense
        self._do_strand_norm_r2 = header.strand_specific and header.read2_antisense
        
        # Columnar block buffer
        self._buffer = ZnaBlockBuffer(seq_len_bytes=header.seq_len_bytes)
        # Cache buffer methods (avoids attribute lookup per record)
        self._buffer_add = self._buffer.add
        
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
        # Strand normalization: reverse complement antisense reads (pre-computed flags)
        if self._do_strand_norm_r1 and is_read1:
            seq = reverse_complement(seq)
        elif self._do_strand_norm_r2 and is_read2:
            seq = reverse_complement(seq)
        
        seq_len = len(seq)
        if seq_len > self._max_len:
            raise ValueError(
                f"Sequence length {seq_len} exceeds maximum {self._max_len} "
                f"allowed by header (seq_len_bytes={self._seq_len_bytes})"
            )
        
        # Use plain integers for flags (avoids enum overhead)
        # Compute flags inline (bit ops are fast)
        flags = (1 if is_read1 else 0) | (2 if is_read2 else 0) | (4 if is_paired else 0)
        
        encoded_seq = _encode_sequence(seq, npolicy=self._npolicy)
        
        # Use cached buffer add method
        self._buffer_add(encoded_seq, flags, seq_len)
        
        # Check if block limit reached
        if self._buffer.size >= self._block_size:
            self._flush_block()
            
    def _flush_block(self) -> None:
        """Compress and write current block to file."""
        if self._buffer.count == 0:
            return
        
        # Get stream sizes before concatenation
        flags_size = self._buffer.flags_size
        lengths_size = self._buffer.lengths_size
        count = self._buffer.count
        
        # Get columnar payload
        raw_data = self._buffer.to_bytes()
        uncompressed_size = len(raw_data)
        
        # Compress if needed
        if self._compressor is not None:
            final_data = self._compressor.compress(raw_data)
        else:
            final_data = raw_data
            
        compressed_size = len(final_data)
        
        # Write block header and data
        block_header = struct.pack(
            _BLOCK_HEADER_FMT, 
            compressed_size, 
            uncompressed_size, 
            count,
            flags_size,
            lengths_size
        )
        
        self._fh.write(block_header)
        self._fh.write(final_data)
        
        # Clear buffer for reuse
        self._buffer.clear()


class ZnaReader:
    """
    Reader for ZNA files with columnar block storage.
    
    Reads blocks in columnar format and yields individual records
    in stream order (First-In-First-Out).
    """
    def __init__(self, fh: BinaryIO):
        """
        Initialize ZNA reader.
        
        Args:
            fh: Binary file handle to read from
        """
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
        decode_table = _DECODE_TABLE
        
        # Strand restoration: reverse complement to recover original antisense reads
        do_restore_r1 = restore_strand and self._header.strand_specific and self._header.read1_antisense
        do_restore_r2 = restore_strand and self._header.strand_specific and self._header.read2_antisense

        if compression_method == COMPRESSION_ZSTD:
            dctx = zstandard.ZstdDecompressor()
        
        while True:
            # 1. Read Block Header
            block_header_data = fh_read(_BLOCK_HEADER_SIZE)
            if not block_header_data:
                break
            if len(block_header_data) < _BLOCK_HEADER_SIZE:
                raise EOFError(f"Incomplete block header. Expected {_BLOCK_HEADER_SIZE} bytes, got {len(block_header_data)}")
                 
            comp_size, uncomp_size, count, flags_size, lengths_size = struct.unpack(
                _BLOCK_HEADER_FMT, block_header_data
            )
            
            # 2. Read Block Data
            block_payload = read_exact(comp_size)
            
            # 3. Decompress if needed
            if compression_method == COMPRESSION_ZSTD:
                block_data = dctx.decompress(block_payload, max_output_size=uncomp_size)
            else:
                block_data = block_payload
            
            # 4. Parse columnar streams
            mv = memoryview(block_data)
            flags_stream = mv[:flags_size]
            lengths_stream = mv[flags_size:flags_size + lengths_size]
            sequences_stream = bytes(mv[flags_size + lengths_size:])  # Convert to bytes for faster indexing
            
            # 5. BATCH DECODE: Decode all sequences in block at once
            # Pre-parse all lengths first (much faster than per-record from_bytes)
            if len_bytes == 1:
                lengths = list(lengths_stream)  # Each byte is a length directly
            elif len_bytes == 2:
                # Unpack all lengths at once using struct
                lengths = struct.unpack(f'<{count}H', lengths_stream[:count * 2])
            else:  # len_bytes == 4
                lengths = struct.unpack(f'<{count}I', lengths_stream[:count * 4])
            
            # Pre-decode all sequences in block
            sequences = _decode_block_sequences(sequences_stream, lengths, decode_table)
            
            # 6. Yield records with minimal per-record overhead
            for i in range(count):
                # Read flag (use plain integers for speed)
                rec_flags = flags_stream[i]
                is_read1 = bool(rec_flags & 1)  # ZnaRecordFlags.IS_READ1
                is_read2 = bool(rec_flags & 2)  # ZnaRecordFlags.IS_READ2
                is_paired = bool(rec_flags & 4)  # ZnaRecordFlags.IS_PAIRED
                
                seq = sequences[i]
                
                # Restore original strand if requested
                if do_restore_r1 and is_read1:
                    seq = reverse_complement(seq)
                elif do_restore_r2 and is_read2:
                    seq = reverse_complement(seq)
                
                yield seq, is_paired, is_read1, is_read2


def _decode_block_sequences(sequences_stream: bytes, lengths: tuple, decode_table: tuple) -> list:
    """
    Batch decode all sequences in a block.
    
    Optimized approach: decode entire sequence stream at once, then split.
    This minimizes per-record overhead by doing bulk string operations.
    
    Args:
        sequences_stream: Raw encoded sequence bytes
        lengths: Tuple of sequence lengths
        decode_table: 256-entry decode table (byte -> 4-char string)
    
    Returns:
        List of decoded sequence strings
    """
    # First pass: decode ALL bytes into a single concatenated string
    # This is ~3x faster than per-sequence decoding because:
    # 1. Single str.join() call instead of 2M calls
    # 2. Single list comprehension with no nested loops
    # 3. No intermediate string allocations per sequence
    all_decoded = ''.join([decode_table[b] for b in sequences_stream])
    
    # Second pass: split into individual sequences by length
    result = []
    char_offset = 0
    for seq_len in lengths:
        # Each base takes 1 character, decode_table decodes 4 chars per byte
        # So we need seq_len characters, padded to multiple of 4
        enc_len = (seq_len + 3) >> 2
        char_end = char_offset + (enc_len * 4)
        result.append(all_decoded[char_offset:char_offset + seq_len])
        char_offset = char_end
    
    return result


def _encode_sequence(seq: str, npolicy: str = None) -> bytes:
    """
    Encode DNA sequence to 2-bit packed bytes.
    
    Each base (A, C, G, T) is encoded as 2 bits (00, 01, 10, 11).
    Four bases are packed into each byte.
    
    Args:
        seq: DNA sequence string (A, C, G, T, possibly N)
        npolicy: Policy for handling N: 'random', 'A', 'C', 'G', 'T', or None
        
    Returns:
        Packed bytes with 4 bases per byte
        
    Raises:
        ValueError: If sequence contains invalid characters
    """
    # Handle N nucleotides if policy is set
    if npolicy and 'N' in seq.upper():
        if npolicy == 'random':
            import random
            seq = ''.join(random.choice('ACGT') if c.upper() == 'N' else c for c in seq)
        elif npolicy in ('A', 'C', 'G', 'T'):
            seq = seq.replace('N', npolicy).replace('n', npolicy.lower())
    
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
    npolicy: str = None,
) -> None:
    """Write ZNA file with header and records.
    
    Args:
        fh: Binary file handle to write to
        header: ZNA file header
        records: Iterable of (sequence, is_paired, is_read1, is_read2) tuples
        npolicy: Policy for handling N nucleotides: 'random', 'A', 'C', 'G', 'T', or None (skip)
    """
    with ZnaWriter(fh, header, npolicy=npolicy) as writer:
        for seq, is_paired, is_read1, is_read2 in records:
            writer.write_record(seq, is_paired, is_read1, is_read2)


def read_zna(
    fh: BinaryIO,
    restore_strand: bool = False,
) -> Tuple[ZnaHeader, Iterator[Tuple[str, bool, bool, bool]]]:
    """Read ZNA file header and return records iterator.
    
    Args:
        fh: Binary file handle to read from
        restore_strand: If True and library is strand-specific, restore
                       original strand by reverse complementing antisense reads.
    
    Returns:
        Tuple of (header, records iterator)
    """
    reader = ZnaReader(fh)
    return reader.header, reader.records(restore_strand=restore_strand)