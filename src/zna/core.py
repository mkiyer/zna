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
import warnings

import zstandard
from dataclasses import dataclass
from enum import IntFlag
from typing import BinaryIO, Iterable, Iterator, Tuple

from .codec import get_backend as _get_backend, get_backend_name as _get_backend_name
from .dtypes import LabelDef, LabelDtype, DTYPE_BY_CODE, label_bytes_per_record, resolve_missing, pack_missing, unpack_missing

# Resolve the best available backend once at import time.
_codec = _get_backend()

# Keep a reference to the Python backend for label operations.
from . import _pycodec as _pycodec_mod

# Check if C++ backend has label support
_accel_mod = None
if _get_backend_name() == "accel":
    try:
        from . import _accel as _accel_mod
        # Verify the new functions exist
        if not hasattr(_accel_mod, 'encode_block_labeled'):
            _accel_mod = None
    except ImportError:
        _accel_mod = None

# Re-export ``reverse_complement`` so existing imports keep working.
reverse_complement = _codec.reverse_complement

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_MAGIC = b"ZNA\x1A"

_VERSION = 2
_SUPPORTED_VERSIONS = frozenset({2})

# File header: Magic(4s) Ver(B) LenBytes(B) Flags(B) CompMethod(B) CompLevel(B) NumLabels(H) RGLen(H) DescLen(H)
_FILE_HEADER_FMT = "<4sBBBBBHHH"
_FILE_HEADER_SIZE = struct.calcsize(_FILE_HEADER_FMT)

# On-disk size of one label definition
_LABEL_DEF_SIZE = 16 + 64 + 1 + 8  # name + desc + dtype_code + missing = 89 bytes

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
    STRAND_NORMALIZED = 8


class ZnaRecordFlags(IntFlag):
    IS_READ1 = 1
    IS_READ2 = 2
    IS_PAIRED = 4
    IS_RC = 8
    #: Both edges of this record are true fragment boundaries — the record spans
    #: its entire fragment (a merged read, or a pair whose insert was shorter
    #: than the read).  Without it, ``IS_RC`` names the single edge that is a
    #: boundary; see :meth:`ZnaReader.records` with ``with_ends=True``.
    IS_FULL_FRAGMENT = 16


def _ends_from_flags(is_rc: bool, is_full_fragment: bool) -> tuple[bool, bool]:
    """``(is_rc, is_full_fragment)`` → ``(has_start, has_end)``.

    A read begins at a fragment boundary and runs inward, so base 0 is a true
    boundary; storing it reverse-complemented moves that boundary to the right
    edge.  A full-fragment record has both.
    """
    if is_full_fragment:
        return True, True
    return (not is_rc), bool(is_rc)


#: ``flag byte -> (has_start, has_end)``, precomputed so the decode hot loop is a
#: single index instead of a branch and a function call.
_ENDS_BY_FLAG = tuple(
    (True, True) if (f & 16) else (not (f & 8), bool(f & 8))
    for f in range(256)
)


def _flags_from_ends(has_start: bool, has_end: bool) -> tuple[bool, bool]:
    """Inverse of :func:`_ends_from_flags` — the mapping is a bijection over the
    three reachable states, so ``(has_start, has_end)`` round-trips losslessly.

    This is the readable reference form.  The write loops that run it on every
    record index :data:`_RC_FULL_BY_ENDS` or :data:`_FLAG_BITS_BY_ENDS` instead,
    to avoid a Python call frame per record; the tables are generated from this
    function below, so the three cannot drift apart.
    """
    return (bool(has_end) and not has_start), (bool(has_start) and bool(has_end))


#: ``flag byte -> (is_paired, is_read1, is_read2)`` — the three booleans
#: :meth:`ZnaReader.records` yields, for consumers of :meth:`ZnaReader.blocks`
#: that get the raw flag column instead and should not re-derive the bit layout.
FLAG_FIELDS = tuple(
    (bool(f & 4), bool(f & 1), bool(f & 2)) for f in range(256)
)

#: ``(has_start << 1 | has_end) -> (is_rc, is_full_fragment)``.
_RC_FULL_BY_ENDS = tuple(
    _flags_from_ends(bool(i & 2), bool(i & 1)) for i in range(4)
)

#: ``(has_start << 1 | has_end) -> IS_RC|IS_FULL_FRAGMENT flag bits``, for the
#: write paths that want the packed byte rather than the pair.
_FLAG_BITS_BY_ENDS = tuple(
    (8 if is_rc else 0) | (16 if is_full else 0) for is_rc, is_full in _RC_FULL_BY_ENDS
)


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
    strand_normalized: bool = False
    compression_method: int = COMPRESSION_NONE
    compression_level: int = DEFAULT_ZSTD_LEVEL
    labels: tuple[LabelDef, ...] = ()

    @property
    def num_labels(self) -> int:
        return len(self.labels)

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
        # Validate labels
        for i, ldef in enumerate(self.labels):
            if ldef.label_id != i:
                raise ValueError(
                    f"Label at index {i} has label_id={ldef.label_id}, expected {i}"
                )
            if len(ldef.name.encode('utf-8')) > 16:
                raise ValueError(f"Label name too long (max 16 UTF-8 bytes): {ldef.name!r}")
            if len(ldef.description.encode('utf-8')) > 64:
                raise ValueError(f"Label description too long (max 64 UTF-8 bytes): {ldef.description!r}")


# ---------------------------------------------------------------------------
# Writer
# ---------------------------------------------------------------------------

class ZnaWriter:
    """Write records to a ZNA file.

    Records are buffered and flushed in blocks.  All sequence
    encoding is handled by the active codec backend.

    ``strand_normalized`` in the header being written describes the **output**.
    By default the writer *derives* orientation: it reverse-complements records
    as the header's strand settings dictate and records what it did in each
    record's ``IS_RC`` flag.  If the input is already normalized — re-encoding
    or copying an existing ZNA file — pass ``preserve_normalization=True`` and
    supply each record's ``is_rc`` verbatim instead.  Orientation is **not
    idempotent**: applying it twice returns the data to an un-normalized state
    while leaving the header still claiming otherwise.

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
        "_do_random_norm",
        "_preserve_normalization",
        "_batch_seqs",
        "_batch_flags",
        "_size_estimate",
        "_pending_r1",
        "_label_defs",
        "_batch_labels",
    )

    def __init__(
        self,
        fh: BinaryIO,
        header: ZnaHeader,
        block_size: int = DEFAULT_BLOCK_SIZE,
        npolicy: str | None = None,
        preserve_normalization: bool = False,
    ) -> None:
        self._fh = fh
        self._header = header
        self._seq_len_bytes = header.seq_len_bytes
        self._max_len = (1 << (8 * header.seq_len_bytes)) - 1
        self._block_size = block_size
        self._npolicy = npolicy or ""
        self._preserve_normalization = preserve_normalization

        if preserve_normalization:
            # Pass-through: the caller supplies both the frame and the IS_RC
            # bit.  Re-deriving orientation here would apply a *second*
            # reverse-complement and silently un-normalize the file.
            self._do_strand_norm_r1 = False
            self._do_strand_norm_r2 = False
            self._do_random_norm = False
        else:
            self._do_strand_norm_r1 = (
                header.strand_normalized and header.strand_specific and header.read1_antisense
            )
            self._do_strand_norm_r2 = (
                header.strand_normalized and header.strand_specific and header.read2_antisense
            )
            self._do_random_norm = header.strand_normalized and not header.strand_specific

        self._batch_seqs: list[str] = []
        self._batch_flags = bytearray()
        self._size_estimate = 0
        self._pending_r1 = False  # True when last buffered record is a paired R1

        self._label_defs = header.labels
        self._batch_labels: list[list] = [[] for _ in self._label_defs]

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
        self, seq: str, is_paired: bool, is_read1: bool, is_read2: bool,
        labels: tuple | None = None, is_rc: bool = False,
        is_full_fragment: bool = False,
    ) -> None:
        """Buffer a single record.  Flushes automatically when the block is full.

        ``is_rc`` records that *seq* is already in its normalized frame and was
        reverse-complemented to get there.  It requires
        ``preserve_normalization=True`` on the writer; otherwise the writer
        derives ``IS_RC`` itself and an externally supplied value would conflict
        with it.

        ``is_full_fragment`` records that *seq* spans its entire fragment, so
        **both** of its edges are true fragment boundaries — a merged read, or a
        mate whose insert was shorter than the read.  The encoder never derives
        this (it cannot know the insert size), so it is always caller-supplied
        and needs no pass-through mode.
        """
        # Validate before mutating any batch state: a caller that catches this
        # must not be left with seqs and flags out of sync.
        if is_rc and not self._preserve_normalization:
            raise ValueError(
                "is_rc=True requires ZnaWriter(preserve_normalization=True). "
                "Without it the writer derives orientation itself, and supplying "
                "IS_RC would double-normalize the record."
            )
        seq_len = len(seq)
        if seq_len > self._max_len:
            raise ValueError(
                f"Sequence length {seq_len} exceeds maximum {self._max_len} "
                f"allowed by header (seq_len_bytes={self._seq_len_bytes})"
            )

        if self._label_defs:
            if labels is None or len(labels) != len(self._label_defs):
                raise ValueError(
                    f"Expected {len(self._label_defs)} label values, "
                    f"got {0 if labels is None else len(labels)}"
                )
            for i, val in enumerate(labels):
                self._batch_labels[i].append(val)

        self._batch_seqs.append(seq)
        # Single expression, no follow-up statements: this runs on every record.
        self._batch_flags.append(
            (1 if is_read1 else 0) | (2 if is_read2 else 0) | (4 if is_paired else 0)
            | (8 if is_rc else 0) | (16 if is_full_fragment else 0)
        )

        self._size_estimate += (seq_len // 4) + 1 + self._seq_len_bytes
        # Defer flush when the last record is a paired R1 to keep pairs together
        if self._size_estimate >= self._block_size:
            if is_paired and is_read1 and self._do_random_norm:
                self._pending_r1 = True
            else:
                self._pending_r1 = False
                self._flush_block()
        elif self._pending_r1 and not (is_paired and is_read1):
            # R2 (or SE) arrived after a deferred R1 — safe to flush now
            self._pending_r1 = False
            if self._size_estimate >= self._block_size:
                self._flush_block()

    def write_records(
        self, records: Iterable[Tuple[str, bool, bool, bool]]
    ) -> None:
        """Buffer many records at once.

        This is faster than calling :meth:`write_record` in a Python loop
        because it caches attribute lookups and avoids per-call method
        dispatch overhead.

        Accepts ``(seq, is_paired, is_read1, is_read2)`` tuples.  Under
        ``preserve_normalization=True`` each tuple carries a trailing ``is_rc``
        — i.e. exactly what ``ZnaReader.records(with_rc=True)`` yields — making a
        lossless ZNA → ZNA copy a one-liner.  Plain 4-tuples are then accepted
        only if the output header is not ``strand_normalized``; otherwise they
        would silently write the orientation away.

        Note: this method does not support labels. Use :meth:`write_record`
        for labeled files.
        """
        if self._label_defs:
            raise TypeError(
                "write_records() does not support labels. "
                "Use write_record() with the labels parameter."
            )
        append_seq = self._batch_seqs.append
        append_flag = self._batch_flags.append
        max_len = self._max_len
        seq_len_bytes = self._seq_len_bytes
        block_size = self._block_size
        flush = self._flush_block
        size_est = self._size_estimate

        if self._preserve_normalization:
            # Pass-through. Accepts either shape the reader produces:
            #   5-tuple  (..., is_rc)                  <- records(with_rc=True)
            #   6-tuple  (..., has_start, has_end)     <- records(with_ends=True)
            # A plain 4-tuple is only acceptable when the output is not
            # normalized, since there is then no orientation to carry; otherwise
            # accepting it would clear every IS_RC bit and lose the orientation
            # silently.
            require_rc = self._header.strand_normalized
            flag_bits_by_ends = _FLAG_BITS_BY_ENDS
            for rec in records:
                seq = rec[0]
                is_paired = rec[1]
                is_read1 = rec[2]
                n_fields = len(rec)
                if n_fields > 5:
                    rc_bit = flag_bits_by_ends[
                        (2 if rec[4] else 0) | (1 if rec[5] else 0)
                    ]
                elif n_fields > 4:
                    rc_bit = 8 if rec[4] else 0
                elif require_rc:
                    raise ValueError(
                        "preserve_normalization=True on a strand-normalized "
                        "header requires each record to carry its orientation, "
                        f"but got a {n_fields}-tuple. Read the source with "
                        "ZnaReader.records(with_ends=True)."
                    )
                else:
                    rc_bit = 0
                seq_len = len(seq)
                if seq_len > max_len:
                    raise ValueError(
                        f"Sequence length {seq_len} exceeds maximum {max_len} "
                        f"allowed by header (seq_len_bytes={seq_len_bytes})"
                    )
                append_seq(seq)
                append_flag(
                    (1 if is_read1 else 0)
                    | (2 if rec[3] else 0)
                    | (4 if is_paired else 0)
                    | rc_bit
                )
                size_est += (seq_len >> 2) + 1 + seq_len_bytes
                if size_est >= block_size:
                    self._size_estimate = size_est
                    flush()
                    size_est = self._size_estimate  # reset after flush
            self._size_estimate = size_est
            return

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
                if is_paired and is_read1 and self._do_random_norm:
                    pass  # defer flush — keep pair together
                else:
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
        if self._header.strand_normalized:
            flags |= ZnaHeaderFlags.STRAND_NORMALIZED

        self._fh.write(
            struct.pack(
                _FILE_HEADER_FMT,
                _MAGIC,
                _VERSION,
                self._header.seq_len_bytes,
                flags,
                self._header.compression_method,
                self._header.compression_level,
                len(self._header.labels),
                len(rg_bytes),
                len(desc_bytes),
            )
        )
        self._fh.write(rg_bytes)
        self._fh.write(desc_bytes)

        # Write label definitions (89 bytes each)
        for ldef in self._header.labels:
            self._fh.write(ldef.name.encode('utf-8').ljust(16, b'\x00')[:16])
            self._fh.write(ldef.description.encode('utf-8').ljust(64, b'\x00')[:64])
            self._fh.write(struct.pack('B', ord(ldef.dtype.code)))
            self._fh.write(pack_missing(ldef))

    def _flush_block(self) -> None:
        """Encode the current batch and write one block."""
        if not self._batch_seqs:
            return

        count = len(self._batch_seqs)

        # 1. Encode via backend
        if self._label_defs:
            # Pre-pack label columns for C++ or use Python backend
            if _accel_mod is not None:
                # Pack each label column into bytes with struct
                label_col_data = []
                label_col_sizes = []
                for col_values, ldef in zip(self._batch_labels, self._label_defs):
                    fmt = f"<{len(col_values)}{ldef.dtype.struct_ch}"
                    label_col_data.append(struct.pack(fmt, *col_values))
                    label_col_sizes.append(ldef.dtype.size)

                flags_bytes, labels_bytes, lengths_bytes, seqs_bytes = _accel_mod.encode_block_labeled(
                    self._batch_seqs,
                    list(self._batch_flags),
                    self._seq_len_bytes,
                    self._npolicy,
                    self._do_strand_norm_r1,
                    self._do_strand_norm_r2,
                    self._do_random_norm,
                    label_col_data,
                    label_col_sizes,
                )
            else:
                flags_bytes, labels_bytes, lengths_bytes, seqs_bytes = _pycodec_mod.encode_block(
                    self._batch_seqs,
                    list(self._batch_flags),
                    self._seq_len_bytes,
                    self._npolicy,
                    self._do_strand_norm_r1,
                    self._do_strand_norm_r2,
                    self._do_random_norm,
                    label_values=self._batch_labels,
                    label_defs=list(self._label_defs),
                )
        else:
            flags_bytes, lengths_bytes, seqs_bytes = _codec.encode_block(
                self._batch_seqs,
                list(self._batch_flags),
                self._seq_len_bytes,
                self._npolicy,
                self._do_strand_norm_r1,
                self._do_strand_norm_r2,
                self._do_random_norm,
            )
            labels_bytes = b""

        # 2. Compress — column order: flags, labels, lengths, sequences
        uncompressed = flags_bytes + labels_bytes + lengths_bytes + seqs_bytes
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
        for col in self._batch_labels:
            col.clear()


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
        self, restore_strand: bool = False, with_rc: bool = False,
        with_ends: bool = False,
    ) -> Iterator[Tuple[str, bool, bool, bool] | Tuple[str, bool, bool, bool, tuple]]:
        """Yield every record in file order.

        Parameters
        ----------
        restore_strand
            If ``True`` and the file was strand-normalized, reverse-complement
            records that were RC'd during encoding (using the per-record
            ``IS_RC`` flag) to restore original orientation.
        with_ends
            If ``True``, yield ``has_start, has_end`` — whether the **left** and
            **right** edge of the stored sequence is a true fragment boundary.
            This is the form a consumer placing fragment-end supervision wants:
            it already combines ``IS_RC`` with ``IS_FULL_FRAGMENT``, so a caller
            never has to re-derive the geometry (and a full-fragment record
            correctly reports *both* edges, which ``with_rc`` alone cannot
            express).  The pair round-trips losslessly back to the two flags.
        with_rc
            If ``True``, include the per-record ``IS_RC`` flag in each tuple.
            This is the only way to recover which edge of a mate is a real
            fragment boundary — it cannot be derived from the sequence (see the
            "Unstranded strand normalization" section of the README).  It is
            mutually exclusive with ``restore_strand``, which *consumes* the
            flag to undo the reverse-complement.

        Yields ``(seq, is_paired, is_read1, is_read2)`` for unlabeled files,
        or ``(seq, is_paired, is_read1, is_read2, labels)`` for labeled files.
        With ``with_rc=True``, ``is_rc`` is inserted **before** ``labels``:
        ``(seq, is_paired, is_read1, is_read2, is_rc)`` and
        ``(seq, is_paired, is_read1, is_read2, is_rc, labels)``.  The default
        widths (4, and 5 labeled) are a compatibility promise and never change.
        """
        if restore_strand and (with_rc or with_ends):
            raise ValueError(
                "restore_strand=True is mutually exclusive with with_rc/with_ends: "
                "restore_strand consumes IS_RC to undo the reverse-complement, so "
                "the orientation they describe has already been undone."
            )
        if with_rc and with_ends:
            raise ValueError(
                "with_rc=True and with_ends=True are mutually exclusive; "
                "with_ends already encodes IS_RC (has_end and not has_start)."
            )
        return self._iter_records(restore_strand, with_rc, with_ends)

    def blocks(
        self, stride: int = 1, offset: int = 0, restore_strand: bool = False,
    ) -> Iterator[Tuple[list, bytes]]:
        """Yield ``(sequences, flags)`` once per block, columnar.

        *sequences* is a ``list[str]``; *flags* is the block's raw flag column,
        one byte per record, positionally aligned with *sequences*.  Decode a
        byte with :class:`ZnaRecordFlags`, or index :data:`FLAG_FIELDS` for the
        ``(is_paired, is_read1, is_read2)`` triple the record API yields.

        This is the batch form for consumers that process a whole block at a
        time — a training loader, say.  It skips the per-record tuple entirely,
        which is most of what is left in :meth:`records` once the sequence itself
        is cheap to produce.

        ``stride``/``offset`` shard **by block**: with ``stride=8, offset=3``
        this yields every eighth block starting at the fourth, and — the point —
        seeks past the others without decompressing or decoding them.  A sharded
        consumer that strides over :meth:`records` instead pays full decode cost
        for every record it then discards, which is ``stride``x more work than it
        needs.  Measured on 200k reads, one worker of 16 costs 9.4x less this way.

        Two things this asks of the caller:

        * **Record order must already be arbitrary.**  Shards get contiguous
          runs, not an interleave, so a file grouped by anything meaningful will
          hand each worker a biased sample.  ``zna shuffle`` exists for this.
        * **The file needs many more blocks than shards.**  Shares are whole
          blocks, so with ``stride`` close to the block count the split is
          lopsided, and past it some shards get nothing at all.  A shard that
          yields no blocks warns rather than looking like an empty file.
          ``DEFAULT_BLOCK_SIZE`` gives a few hundred blocks per GB; a small file
          split many ways wants a smaller ``block_size`` at write time.

        Raises on labeled files: the label columns would have to come back too,
        and silently dropping them is worse than not offering the API.
        """
        if stride < 1:
            raise ValueError(f"stride must be >= 1, got {stride}")
        if not 0 <= offset < stride:
            raise ValueError(f"offset must be in [0, {stride}), got {offset}")
        if self._header.labels:
            raise NotImplementedError(
                "blocks() does not support labeled files "
                f"({len(self._header.labels)} label column(s) defined). "
                "Use records(), which returns labels with each record."
            )
        return self._iter_blocks(stride, offset, restore_strand)

    def _iter_blocks(
        self, stride: int, offset: int, restore_strand: bool
    ) -> Iterator[Tuple[list, bytes]]:
        fh = self._fh
        fh_read = fh.read
        len_bytes = self._header.seq_len_bytes
        compression_method = self._header.compression_method
        needs_restore = restore_strand and self._header.strand_normalized
        decode_seqs = _codec.decode_block_sequences

        if compression_method == COMPRESSION_ZSTD:
            dctx = zstandard.ZstdDecompressor()

        index = 0
        yielded = 0
        while True:
            block_header_data = fh_read(_BLOCK_HEADER_SIZE)
            if not block_header_data:
                if stride > 1 and yielded == 0 and index > 0:
                    # Silence here is indistinguishable from an empty file, and
                    # in a training loader it is an idle worker nobody notices.
                    warnings.warn(
                        f"blocks(stride={stride}, offset={offset}) matched none "
                        f"of this file's {index} block(s), so this shard has no "
                        f"data. Shards are whole blocks: write the file with a "
                        f"smaller block_size, or use fewer shards.",
                        RuntimeWarning, stacklevel=3,
                    )
                break
            if len(block_header_data) < _BLOCK_HEADER_SIZE:
                raise EOFError(
                    f"Incomplete block header. Expected {_BLOCK_HEADER_SIZE} "
                    f"bytes, got {len(block_header_data)}"
                )
            comp_size, uncomp_size, count, flags_size, lengths_size = struct.unpack(
                _BLOCK_HEADER_FMT, block_header_data
            )

            index += 1
            if (index - 1) % stride != offset:
                # Not this shard's block: step over the payload without paying
                # to decompress it.  Seek when the stream allows it; a pipe has
                # to be read through.  ``_read_exact`` raises on a short read.
                try:
                    fh.seek(comp_size, 1)
                except (AttributeError, OSError):
                    self._read_exact(comp_size)
                continue

            block_payload = self._read_exact(comp_size)
            if compression_method == COMPRESSION_ZSTD:
                block_data = dctx.decompress(block_payload, max_output_size=uncomp_size)
            else:
                block_data = block_payload

            lengths_end = flags_size + lengths_size
            flags_stream = block_data[:flags_size]
            lengths_stream = block_data[flags_size:lengths_end]
            seqs_stream = block_data[lengths_end:]

            yielded += 1
            yield (
                decode_seqs(flags_stream, lengths_stream, seqs_stream,
                            len_bytes, count, needs_restore),
                flags_stream,
            )

    def _iter_records(
        self, restore_strand: bool, with_rc: bool, with_ends: bool = False
    ) -> Iterator[tuple]:
        fh_read = self._fh.read
        read_exact = self._read_exact
        len_bytes = self._header.seq_len_bytes
        compression_method = self._header.compression_method
        label_defs = self._header.labels
        has_labels = len(label_defs) > 0

        # Pre-compute per-record label byte sizes
        if has_labels:
            label_col_sizes = [label_bytes_per_record(ld) for ld in label_defs]
            total_label_bytes_per_rec = sum(label_col_sizes)
        else:
            label_col_sizes = []
            total_label_bytes_per_rec = 0

        needs_restore = restore_strand and self._header.strand_normalized

        # What the backend should emit.  The three consumers below want three
        # different widths, and asking the backend for the right one is much
        # cheaper than unpacking and rebuilding every record tuple in Python:
        #   with_rc      -> yield the backend's tuples untouched
        #   needs_restore-> backend reverse-complements and drops IS_RC
        #   otherwise    -> backend drops IS_RC
        # ``with_ends`` reads IS_FULL_FRAGMENT straight off the flags column, so
        # it needs the narrow width too.
        want_rc = with_rc
        decode_kwargs = {"with_rc": want_rc, "restore_strand": needs_restore}

        # Choose decode backend
        use_accel_labels = has_labels and _accel_mod is not None
        if use_accel_labels:
            label_dtype_codes = "".join(ld.dtype.code for ld in label_defs)
        decode = _codec.decode_block if not has_labels else None  # set per block below

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

            # Split columnar streams: flags → labels → lengths → sequences
            flags_end = flags_size

            if has_labels:
                total_label_bytes = total_label_bytes_per_rec * count
                labels_end = flags_end + total_label_bytes
                # Split individual label columns
                offset = flags_end
                label_columns: list[bytes] = []
                for col_size in label_col_sizes:
                    col_bytes = col_size * count
                    label_columns.append(block_data[offset:offset + col_bytes])
                    offset += col_bytes
            else:
                labels_end = flags_end
                label_columns = []

            lengths_end = labels_end + lengths_size
            flags_stream = block_data[:flags_end]
            lengths_stream = block_data[labels_end:lengths_end]
            seqs_stream = block_data[lengths_end:]

            # Decode via backend
            if has_labels:
                if use_accel_labels:
                    block_records = _accel_mod.decode_block_labeled(
                        flags_stream, lengths_stream, seqs_stream, len_bytes, count,
                        label_columns, label_col_sizes, label_dtype_codes,
                        **decode_kwargs,
                    )
                else:
                    block_records = _pycodec_mod.decode_block(
                        flags_stream, lengths_stream, seqs_stream, len_bytes, count,
                        label_columns=label_columns, label_defs=list(label_defs),
                        **decode_kwargs,
                    )
            else:
                block_records = decode(
                    flags_stream, lengths_stream, seqs_stream, len_bytes, count,
                    **decode_kwargs,
                )

            if with_ends:
                # IS_FULL_FRAGMENT is read straight off the flags column rather
                # than widening the backends' record tuples, so both decoders
                # (and the compiled extension) are untouched by this feature.
                # zip over the flags bytes and a table lookup keep this close to
                # the cost of the plain path.
                ends = _ENDS_BY_FLAG
                if has_labels:
                    for rec, fl in zip(block_records, flags_stream):
                        has_start, has_end = ends[fl]
                        yield rec[0], rec[1], rec[2], rec[3], has_start, has_end, rec[4]
                else:
                    for rec, fl in zip(block_records, flags_stream):
                        has_start, has_end = ends[fl]
                        yield rec[0], rec[1], rec[2], rec[3], has_start, has_end
            else:
                # The backend was asked for exactly the shape this call yields:
                # IS_RC present only under with_rc, and the reverse-complement
                # already undone under restore_strand.  Nothing left to rebuild.
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
            num_labels,
            rg_len,
            desc_len,
        ) = struct.unpack(_FILE_HEADER_FMT, fixed)

        if magic != _MAGIC:
            raise ValueError("Not a ZNA file")
        if ver not in _SUPPORTED_VERSIONS:
            raise ValueError(f"Unsupported ZNA version: {ver}")

        read_group = self._read_exact(rg_len).decode("utf-8")
        description = self._read_exact(desc_len).decode("utf-8")

        # Read label definitions
        label_defs: list[LabelDef] = []
        for i in range(num_labels):
            name = self._read_exact(16).rstrip(b'\x00').decode('utf-8')
            desc = self._read_exact(64).rstrip(b'\x00').decode('utf-8')
            dtype_code = chr(struct.unpack('B', self._read_exact(1))[0])
            missing_bytes = self._read_exact(8)
            if dtype_code not in DTYPE_BY_CODE:
                raise ValueError(f"Unknown label dtype code: {dtype_code!r}")
            dtype = DTYPE_BY_CODE[dtype_code]
            missing_val = unpack_missing(dtype, missing_bytes)
            label_defs.append(LabelDef(
                label_id=i, name=name, description=desc,
                dtype=dtype, missing=missing_val,
            ))

        return ZnaHeader(
            read_group=read_group,
            description=description,
            seq_len_bytes=len_bytes,
            strand_specific=bool(flags & ZnaHeaderFlags.STRAND_SPECIFIC),
            read1_antisense=bool(flags & ZnaHeaderFlags.READ1_ANTISENSE),
            read2_antisense=bool(flags & ZnaHeaderFlags.READ2_ANTISENSE),
            strand_normalized=bool(flags & ZnaHeaderFlags.STRAND_NORMALIZED),
            compression_method=comp_method,
            compression_level=comp_level,
            labels=tuple(label_defs),
        )


# ---------------------------------------------------------------------------
# Convenience helpers
# ---------------------------------------------------------------------------


def write_zna(
    fh: BinaryIO,
    header: ZnaHeader,
    records: Iterable[Tuple[str, bool, bool, bool]],
    npolicy: str | None = None,
    preserve_normalization: bool = False,
) -> None:
    """Write a complete ZNA file from an iterable of records.

    Pass ``preserve_normalization=True`` when *records* already carry their
    normalized frame and a trailing ``is_rc`` (as from
    :meth:`ZnaReader.records` with ``with_rc=True``); see :class:`ZnaWriter`.
    """
    with ZnaWriter(
        fh, header, npolicy=npolicy, preserve_normalization=preserve_normalization
    ) as writer:
        writer.write_records(records)


def read_zna(
    fh: BinaryIO,
    restore_strand: bool = False,
    with_rc: bool = False,
) -> Tuple[ZnaHeader, Iterator[Tuple[str, bool, bool, bool]]]:
    """Read a ZNA file header and return a record iterator."""
    reader = ZnaReader(fh)
    return reader.header, reader.records(restore_strand=restore_strand, with_rc=with_rc)
