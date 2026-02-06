"""
ZNA: compressed nucleic acid data format.

This module defines the binary file format — magic bytes, headers, block
framing, and compression.  All sequence encoding/decoding is delegated to
a pluggable *codec backend* (see :mod:`zna.codec`).

Binary layout::

    [File Header]
    [Block 0: Header + Payload]
    [Block 1: Header + Payload]
    ...

Each block payload contains three concatenated columnar streams
(flags, lengths, sequences) compressed as a single ZSTD frame.
"""
from __future__ import annotations

import struct
import zstandard
from dataclasses import dataclass
from enum import IntFlag
from typing import BinaryIO, Iterable, Iterator, Tuple

from .codec import get_backend as _get_backend

# Resolve the best available backend once at import time.
_codec = _get_backend()

# Re-export ``reverse_complement`` so existing imports keep working.
reverse_complement = _codec.reverse_complement

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_MAGIC = b"ZNA\x1A"

_VERSION = 1
_SUPPORTED_VERSIONS = frozenset({1})

# File header: Magic(4s) Ver(B) LenBytes(B) Flags(B) CompMethod(B) CompLevel(B) RGLen(H) DescLen(H)
_FILE_HEADER_FMT = "<4sBBBBBHH"
_FILE_HEADER_SIZE = struct.calcsize(_FILE_HEADER_FMT)

# Block header: CompSize(I) UncompSize(I) RecordCount(I) FlagsSize(I) LengthsSize(I)
_BLOCK_HEADER_FMT = "<IIIII"
_BLOCK_HEADER_SIZE = struct.calcsize(_BLOCK_HEADER_FMT)

# Sequence-length field widths
MIN_SEQ_LEN_BYTES = 1
MAX_SEQ_LEN_BYTES = 4
MAX_METADATA_LENGTH = 65535

# Compression types
COMPRESSION_NONE = 0
COMPRESSION_ZSTD = 1

# Defaults
DEFAULT_ZSTD_LEVEL = 9
DEFAULT_BLOCK_SIZE = 4 * 1024 * 1024  # 4 MB

_VALID_SEQ_LEN_BYTES = frozenset({1, 2, 4})

# ---------------------------------------------------------------------------
# Flags
# ---------------------------------------------------------------------------

class ZnaHeaderFlags(IntFlag):
    STRAND_SPECIFIC = 1
    READ1_ANTISENSE = 2
    READ2_ANTISENSE = 4


class ZnaRecordFlags(IntFlag):
    IS_READ1 = 1
    IS_READ2 = 2
    IS_PAIRED = 4


# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------

@dataclass(slots=True)
class ZnaHeader:
    """Metadata stored in the ZNA file header."""

    read_group: str
    description: str = ""
    seq_len_bytes: int = MIN_SEQ_LEN_BYTES
    strand_specific: bool = False
    read1_antisense: bool = False
    read2_antisense: bool = False
    compression_method: int = COMPRESSION_NONE
    compression_level: int = DEFAULT_ZSTD_LEVEL

    def __post_init__(self) -> None:
        if self.seq_len_bytes not in _VALID_SEQ_LEN_BYTES:
            raise ValueError(
                f"seq_len_bytes must be one of {_VALID_SEQ_LEN_BYTES}, "
                f"got {self.seq_len_bytes}"
            )
        if len(self.read_group.encode('utf-8')) > MAX_METADATA_LENGTH:
            raise ValueError("read_group too long")
        if len(self.description.encode('utf-8')) > MAX_METADATA_LENGTH:
            raise ValueError("description too long")
        if self.compression_method not in (COMPRESSION_NONE, COMPRESSION_ZSTD):
            raise ValueError(
                f"compression_method must be 0 or 1, got {self.compression_method}"
            )
        if not 1 <= self.compression_level <= 22:
            raise ValueError(
                f"compression_level must be 1-22, got {self.compression_level}"
            )


# ---------------------------------------------------------------------------
# Writer
# ---------------------------------------------------------------------------

class ZnaWriter:
    """Write records to a ZNA file.

    Records are buffered and flushed in blocks.  All sequence
    encoding is handled by the active codec backend.

    Block layout (V1)::

        [20-byte header]
        [ZSTD-compressed payload: flags ‖ lengths ‖ sequences]
    """

    __slots__ = (
        "_fh",
        "_header",
        "_seq_len_bytes",
        "_max_len",
        "_block_size",
        "_npolicy",
        "_compressor",
        "_do_strand_norm_r1",
        "_do_strand_norm_r2",
        "_batch_seqs",
        "_batch_flags",
        "_size_estimate",
    )

    def __init__(
        self,
        fh: BinaryIO,
        header: ZnaHeader,
        block_size: int = DEFAULT_BLOCK_SIZE,
        npolicy: str | None = None,
    ) -> None:
        self._fh = fh
        self._header = header
        self._seq_len_bytes = header.seq_len_bytes
        self._max_len = (1 << (8 * header.seq_len_bytes)) - 1
        self._block_size = block_size
        self._npolicy = npolicy or ""

        self._do_strand_norm_r1 = header.strand_specific and header.read1_antisense
        self._do_strand_norm_r2 = header.strand_specific and header.read2_antisense

        self._batch_seqs: list[str] = []
        self._batch_flags = bytearray()
        self._size_estimate = 0

        if header.compression_method == COMPRESSION_ZSTD:
            self._compressor = zstandard.ZstdCompressor(level=header.compression_level)
        else:
            self._compressor = None

        self._write_file_header()

    # -- context manager -----------------------------------------------------

    def __enter__(self) -> ZnaWriter:
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:  # noqa: ANN001
        self.close()

    def close(self) -> None:
        self._flush_block()

    # -- public --------------------------------------------------------------

    def write_record(
        self, seq: str, is_paired: bool, is_read1: bool, is_read2: bool
    ) -> None:
        """Buffer a single record.  Flushes automatically when the block is full."""
        seq_len = len(seq)
        if seq_len > self._max_len:
            raise ValueError(
                f"Sequence length {seq_len} exceeds maximum {self._max_len} "
                f"allowed by header (seq_len_bytes={self._seq_len_bytes})"
            )

        self._batch_seqs.append(seq)
        flag = (1 if is_read1 else 0) | (2 if is_read2 else 0) | (4 if is_paired else 0)
        self._batch_flags.append(flag)

        self._size_estimate += (seq_len // 4) + 1 + self._seq_len_bytes
        if self._size_estimate >= self._block_size:
            self._flush_block()

    def write_records(
        self, records: Iterable[Tuple[str, bool, bool, bool]]
    ) -> None:
        """Buffer many records at once.

        This is faster than calling :meth:`write_record` in a Python loop
        because it caches attribute lookups and avoids per-call method
        dispatch overhead.
        """
        append_seq = self._batch_seqs.append
        append_flag = self._batch_flags.append
        max_len = self._max_len
        seq_len_bytes = self._seq_len_bytes
        block_size = self._block_size
        flush = self._flush_block
        size_est = self._size_estimate

        for seq, is_paired, is_read1, is_read2 in records:
            seq_len = len(seq)
            if seq_len > max_len:
                raise ValueError(
                    f"Sequence length {seq_len} exceeds maximum {max_len} "
                    f"allowed by header (seq_len_bytes={seq_len_bytes})"
                )
            append_seq(seq)
            append_flag(
                (1 if is_read1 else 0)
                | (2 if is_read2 else 0)
                | (4 if is_paired else 0)
            )
            size_est += (seq_len >> 2) + 1 + seq_len_bytes
            if size_est >= block_size:
                self._size_estimate = size_est
                flush()
                size_est = self._size_estimate  # reset after flush

        self._size_estimate = size_est

    # -- private -------------------------------------------------------------

    def _write_file_header(self) -> None:
        rg_bytes = self._header.read_group.encode("utf-8")
        desc_bytes = self._header.description.encode("utf-8")

        flags = 0
        if self._header.strand_specific:
            flags |= ZnaHeaderFlags.STRAND_SPECIFIC
        if self._header.read1_antisense:
            flags |= ZnaHeaderFlags.READ1_ANTISENSE
        if self._header.read2_antisense:
            flags |= ZnaHeaderFlags.READ2_ANTISENSE

        self._fh.write(
            struct.pack(
                _FILE_HEADER_FMT,
                _MAGIC,
                _VERSION,
                self._header.seq_len_bytes,
                flags,
                self._header.compression_method,
                self._header.compression_level,
                len(rg_bytes),
                len(desc_bytes),
            )
        )
        self._fh.write(rg_bytes)
        self._fh.write(desc_bytes)

    def _flush_block(self) -> None:
        """Encode the current batch and write one block."""
        if not self._batch_seqs:
            return

        count = len(self._batch_seqs)

        # 1. Encode via backend
        flags_bytes, lengths_bytes, seqs_bytes = _codec.encode_block(
            self._batch_seqs,
            list(self._batch_flags),
            self._seq_len_bytes,
            self._npolicy,
            self._do_strand_norm_r1,
            self._do_strand_norm_r2,
        )

        # 2. Compress
        uncompressed = flags_bytes + lengths_bytes + seqs_bytes
        compressed = (
            self._compressor.compress(uncompressed)
            if self._compressor is not None
            else uncompressed
        )

        # 3. Write block header + payload
        self._fh.write(
            struct.pack(
                _BLOCK_HEADER_FMT,
                len(compressed),
                len(uncompressed),
                count,
                len(flags_bytes),
                len(lengths_bytes),
            )
        )
        self._fh.write(compressed)

        # 4. Reset
        self._batch_seqs.clear()
        self._batch_flags.clear()
        self._size_estimate = 0


# ---------------------------------------------------------------------------
# Reader
# ---------------------------------------------------------------------------

class ZnaReader:
    """Read records from a ZNA file.

    Yields ``(sequence, is_paired, is_read1, is_read2)`` tuples.
    """

    def __init__(self, fh: BinaryIO) -> None:
        self._fh = fh
        self._header = self._read_file_header()

    @property
    def header(self) -> ZnaHeader:
        return self._header

    def __iter__(self) -> Iterator[Tuple[str, bool, bool, bool]]:
        return self.records()

    # -- public --------------------------------------------------------------

    def records(
        self, restore_strand: bool = False
    ) -> Iterator[Tuple[str, bool, bool, bool]]:
        """Yield every record in file order.

        Parameters
        ----------
        restore_strand
            If ``True`` and the library is strand-specific, reverse-complement
            antisense reads back to their original orientation.
        """
        fh_read = self._fh.read
        read_exact = self._read_exact
        len_bytes = self._header.seq_len_bytes
        compression_method = self._header.compression_method

        do_restore_r1 = (
            restore_strand
            and self._header.strand_specific
            and self._header.read1_antisense
        )
        do_restore_r2 = (
            restore_strand
            and self._header.strand_specific
            and self._header.read2_antisense
        )

        needs_restore = do_restore_r1 or do_restore_r2
        rc = _codec.reverse_complement
        decode = _codec.decode_block

        if compression_method == COMPRESSION_ZSTD:
            dctx = zstandard.ZstdDecompressor()

        while True:
            block_header_data = fh_read(_BLOCK_HEADER_SIZE)
            if not block_header_data:
                break
            if len(block_header_data) < _BLOCK_HEADER_SIZE:
                raise EOFError(
                    f"Incomplete block header. Expected {_BLOCK_HEADER_SIZE} "
                    f"bytes, got {len(block_header_data)}"
                )

            comp_size, uncomp_size, count, flags_size, lengths_size = struct.unpack(
                _BLOCK_HEADER_FMT, block_header_data
            )

            block_payload = read_exact(comp_size)

            if compression_method == COMPRESSION_ZSTD:
                block_data = dctx.decompress(block_payload, max_output_size=uncomp_size)
            else:
                block_data = block_payload

            # Split columnar streams
            flags_end = flags_size
            lengths_end = flags_end + lengths_size
            flags_stream = block_data[:flags_end]
            lengths_stream = block_data[flags_end:lengths_end]
            seqs_stream = block_data[lengths_end:]

            # Decode via backend
            block_records = decode(
                flags_stream, lengths_stream, seqs_stream, len_bytes, count
            )

            if needs_restore:
                for seq, is_paired, is_read1, is_read2 in block_records:
                    if do_restore_r1 and is_read1:
                        seq = rc(seq)
                    elif do_restore_r2 and is_read2:
                        seq = rc(seq)
                    yield seq, is_paired, is_read1, is_read2
            else:
                yield from block_records

    # -- private -------------------------------------------------------------

    def _read_exact(self, n: int) -> bytes:
        data = self._fh.read(n)
        if len(data) != n:
            raise EOFError("Unexpected EOF while reading ZNA")
        return data

    def _read_file_header(self) -> ZnaHeader:
        fixed = self._read_exact(_FILE_HEADER_SIZE)
        (
            magic,
            ver,
            len_bytes,
            flags,
            comp_method,
            comp_level,
            rg_len,
            desc_len,
        ) = struct.unpack(_FILE_HEADER_FMT, fixed)

        if magic != _MAGIC:
            raise ValueError("Not a ZNA file")
        if ver not in _SUPPORTED_VERSIONS:
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
            compression_method=comp_method,
            compression_level=comp_level,
        )


# ---------------------------------------------------------------------------
# Convenience helpers
# ---------------------------------------------------------------------------


def write_zna(
    fh: BinaryIO,
    header: ZnaHeader,
    records: Iterable[Tuple[str, bool, bool, bool]],
    npolicy: str | None = None,
) -> None:
    """Write a complete ZNA file from an iterable of records."""
    with ZnaWriter(fh, header, npolicy=npolicy) as writer:
        writer.write_records(records)


def read_zna(
    fh: BinaryIO,
    restore_strand: bool = False,
) -> Tuple[ZnaHeader, Iterator[Tuple[str, bool, bool, bool]]]:
    """Read a ZNA file header and return a record iterator."""
    reader = ZnaReader(fh)
    return reader.header, reader.records(restore_strand=restore_strand)
