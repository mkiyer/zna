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

**Blocks are fragment-complete.**  A fragment's reads are stored consecutively,
R1 immediately followed by R2, and never span a block boundary.  So a block is a
self-contained unit of molecules, and any set of blocks can be decoded
independently of the rest of the file — which is what makes :meth:`ZnaReader.blocks`
usable as a parallelism and sharding primitive.  :class:`ZnaWriter` enforces both
halves of this on every write path; see :data:`OPENS_FRAGMENT`.
"""
from __future__ import annotations

import json
import struct
import warnings
import zlib
from itertools import compress

import zstandard
from dataclasses import dataclass
from enum import IntFlag
from typing import BinaryIO, Iterable, Iterator, NamedTuple, Tuple

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

_VERSION = 3
_SUPPORTED_VERSIONS = frozenset({3})

# File header: Magic(4s) Ver(B) LenBytes(B) Flags(B) CompMethod(B) CompLevel(B) NumLabels(H) RGLen(H) DescLen(H)
_FILE_HEADER_FMT = "<4sBBBBBHHH"
_FILE_HEADER_SIZE = struct.calcsize(_FILE_HEADER_FMT)

# On-disk size of one label definition
_LABEL_DEF_SIZE = 16 + 64 + 1 + 8  # name + desc + dtype_code + missing = 89 bytes

# Block header: CompSize(I) UncompSize(I) RecordCount(I) FlagsSize(I) LengthsSize(I)
_BLOCK_HEADER_FMT = "<IIIII"
_BLOCK_HEADER_SIZE = struct.calcsize(_BLOCK_HEADER_FMT)

# A block header whose record count is 0 is never a data block -- ``_flush_block``
# returns without writing on an empty batch, so no real block has ever carried
# count == 0 and the value is reserved (from 0.5.0) for the two METADATA
# pseudo-blocks: the provenance PROLOGUE right after the file header, and the
# trailer SENTINEL after the last data block.  Position disambiguates them, and
# only at open: ``ZnaReader.__init__`` consumes the prologue, so to every walk
# loop count == 0 has exactly one meaning -- end of data.
#
# Trailer footer, fixed 32 bytes at EOF:
#   magic | footer_version | reserved | crc32 of the trailer payload AS WRITTEN
#   (compressed form, checkable without inflating) | sentinel_offset | 8 pad bytes
# ``sentinel_offset`` is the byte offset of the sentinel block header, i.e. the
# end of the data blocks -- the value a block-packing consumer needs as its
# final fence now that "the last block ends exactly at EOF" is no longer true.
_FOOTER_FMT = "<8sHHIQ8x"
_FOOTER_SIZE = struct.calcsize(_FOOTER_FMT)  # 32
_FOOTER_MAGIC = b"ZNATRLR\x1A"
_FOOTER_VERSION = 1

#: Schema versions of the two JSON payloads, independent of the format byte.
PROLOGUE_SCHEMA = 1
TRAILER_SCHEMA = 1

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


#: ``flag byte -> (has_start, has_end)`` — whether the left and right edge of the
#: stored sequence is a true fragment boundary.  Precomputed so the decode hot
#: loop is a single index instead of a branch and a function call.
#:
#: Public because it is the only correct way for a :meth:`ZnaReader.blocks`
#: consumer to recover fragment geometry from the raw flags column.  It cannot be
#: inferred from the mate number: under unstranded normalization ZNA picks *one
#: mate per pair at random* to reverse-complement, so which edge is the real
#: boundary is a per-record fact carried by ``IS_RC``, not a property of R1 vs R2.
ENDS_BY_FLAG = tuple(
    (True, True) if (f & 16) else (not (f & 8), bool(f & 8))
    for f in range(256)
)

#: Back-compat alias for the private name used before 0.3.4.
_ENDS_BY_FLAG = ENDS_BY_FLAG


# There is deliberately no inverse of :func:`_ends_from_flags`, and there used to be.
#
# ``_flags_from_ends`` claimed to be "a bijection over the three reachable states".
# There are **four**: byte 24 (``IS_RC | IS_FULL_FRAGMENT``) is what the encoder writes
# whenever strand normalization reverse-complements a merged read.  ``ENDS_BY_FLAG``
# maps 16 and 24 alike to ``(True, True)`` — correctly, because both records do have two
# real fragment ends — so no inverse can exist, and the one that existed silently
# cleared ``IS_RC`` on every full-fragment record that passed through it.  Measured on a
# 200-record normalized file: 101 records carried ``IS_RC`` before ``zna shuffle`` and 0
# after, at which point ``--restore-strand`` returned half the corpus in the wrong
# orientation.  The same loss ran through ZNA -> ZNA re-encode.
#
# **The rule this file now keeps: a view is for reading; the flag byte is for copying.**
# ``(has_start, has_end)`` is a decode-side projection *for consumers* — it answers
# "which edges of this record are true fragment boundaries", which is what a downstream
# model wants and is genuinely all it wants.  It must never travel back into a writer.
# Copying carries the byte (:meth:`ZnaReader.copy_records` -> :meth:`ZnaWriter.write_copy`);
# producing carries ``is_full_fragment`` (:meth:`ZnaWriter.write_record`).


#: ``flag byte -> (is_paired, is_read1, is_read2)`` — the three booleans
#: :meth:`ZnaReader.records` yields, for consumers of :meth:`ZnaReader.blocks`
#: that get the raw flag column instead and should not re-derive the bit layout.
FLAG_FIELDS = tuple(
    (bool(f & 4), bool(f & 1), bool(f & 2)) for f in range(256)
)


#: ``flag byte -> True`` for a paired R1: the record that *opens* a two-read fragment.
#: Its mate must be the very next record, and must land in the same block — so this is
#: exactly the record a block may not end on.
OPENS_FRAGMENT = tuple((f & 0x05) == 0x05 for f in range(256))

#: ``flag byte -> True`` for a paired R2: the record that *closes* a two-read fragment.
#: It is the only record that may follow a paired R1, and may not appear anywhere else.
CLOSES_FRAGMENT = tuple((f & 0x06) == 0x06 for f in range(256))

# Together these are the whole fragment contract, and the writer needs no other state to
# enforce it than "did the last record open a fragment".  ``CLOSES_FRAGMENT[flag] must
# equal that bit`` rejects both halves of a broken stream in one comparison -- an R1
# followed by anything but its R2, and an R2 with no R1 in front of it -- and deciding
# the flush on ``OPENS_FRAGMENT[flag]`` keeps every fragment whole, because the only
# record a block is forbidden to end on is the one whose mate has not arrived yet.
#
# That bound is why no deferral bookkeeping is needed: a fragment is at most two
# records, so a block overruns ``block_size`` by at most one record.  The older code
# deferred a flush *without* enforcing the contract, which had no such bound -- a run of
# consecutive paired R1s buffered the entire stream in memory and wrote nothing.


class ZnaRecord(NamedTuple):
    """One record as it is **stored**, for copying between ZNA files.

    ``flags`` is the raw :class:`ZnaRecordFlags` byte, verbatim — not a projection of
    it — so a copy carries every bit, including combinations this version of the library
    does not interpret and bits a future format version may add.

    Three fields, and named: no ``records()`` shape is three wide, so code that confuses
    the two fails at the unpack instead of silently reading ``flags`` as ``is_paired``.
    """

    seq: str
    flags: int
    labels: tuple = ()


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


#: ``flag byte -> 1`` when the record is unpaired (single or merged read), for
#: selecting unpaired lengths out of a block's columns at C speed via
#: ``bytes.translate`` + ``itertools.compress``.
_UNPAIRED_SELECTOR = bytes(0 if (f & 0x04) else 1 for f in range(256))


def _count_into(hist: dict, values) -> None:
    """Histogram update at C speed (the same mechanism Counter uses)."""
    from collections import _count_elements  # type: ignore[attr-defined]
    _count_elements(hist, values)


def _canonical_json(obj) -> bytes:
    """The one serialization both metadata payloads use.

    Sorted keys, no whitespace, ASCII-only: identical input must produce
    byte-identical files, so nothing here may depend on dict order, locale,
    or anything derived from the environment.
    """
    return json.dumps(
        obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")


@dataclass(slots=True)
class ZnaProvenance:
    """The prologue: facts known at encode START, readable before any record.

    Written by :class:`ZnaWriter` immediately after the file header, so a
    consumer -- including one on a pipe, which can never seek to the trailer --
    can answer "what wrote this, and may it enter training?" before decoding a
    single block.  Facts that require the whole encode live in
    :class:`ZnaTrailer` instead; the split is principled, not incidental:
    position in the file follows when the information is knowable.
    """

    writer_version: str
    shuffled: bool
    merged_in_process: bool
    #: crc32 of the prologue payload as written; :class:`ZnaTrailer` re-states
    #: it, so the end of the file attests the start.
    crc32: int
    #: The parsed payload verbatim, including fields a future schema may add.
    raw: dict


@dataclass(slots=True)
class ZnaTrailer:
    """The trailer: facts derivable only from the COMPLETE encode.

    Tallied in ``_flush_block`` from the codec's returned columns -- never from
    the writer's input arguments, which differ from the file whenever strand
    normalization sets ``IS_RC`` or an N-policy changes stored lengths.
    """

    n_records: int
    n_bases: int
    n_pairs: int
    n_unpaired: int
    #: flag byte -> record count, zero entries omitted.
    flag_counts: dict
    #: stored record length -> count, over all records.
    length_histogram: dict
    #: same, over unpaired records only (merged reads and singles).
    length_histogram_unpaired: dict
    #: Byte offset of block 0's header (after header and prologue).
    data_start: int
    #: Per-block stored index; offset[i] = data_start + sum(20 + comp_sizes[:i]).
    block_comp_sizes: list
    block_uncomp_sizes: list
    block_records: list
    #: crc32 of the prologue payload, closing the integrity gap for
    #: uncompressed files whose prologue carries no zstd content checksum.
    prologue_crc32: int
    #: The parsed payload verbatim, including fields a future schema may add.
    raw: dict

    def block_offsets(self) -> list[int]:
        """Absolute byte offset of every block header, by running sum."""
        out = []
        off = self.data_start
        for c in self.block_comp_sizes:
            out.append(off)
            off += _BLOCK_HEADER_SIZE + c
        return out


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

    Block layout::

        [20-byte header]
        [compressed payload: flags ‖ labels ‖ lengths ‖ sequences]

    The labels column is present only for a labeled file, and the column order is
    the compression order -- flags first because they are the most redundant.

    **Fragments are written whole.**  Every write path here requires a paired R1 to be
    followed immediately by its R2, and refuses to end a block between them, so a block
    always holds entire molecules.  Callers that group their input already satisfy this;
    one that does not gets a :exc:`ValueError` naming the record, rather than a file
    whose blocks cannot be processed independently.
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
        "_expect_mate",
        "_rng_seed",
        "_records_written",
        "_label_defs",
        "_batch_labels",
        "_shuffled",
        "_merged_in_process",
        "_len_char",
        "_flag_counts",
        "_len_hist",
        "_len_hist_unpaired",
        "_n_bases",
        "_block_meta",
        "_data_start_offset",
        "_data_end",
        "_prologue_crc32",
        "_closed",
    )

    def __init__(
        self,
        fh: BinaryIO,
        header: ZnaHeader,
        block_size: int = DEFAULT_BLOCK_SIZE,
        npolicy: str | None = None,
        preserve_normalization: bool = False,
        rng_seed: int = 0,
        shuffled: bool = False,
        merged_in_process: bool = False,
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
        # True when the last buffered record was a paired R1, i.e. a fragment is open
        # and the next record must be the R2 that closes it.
        self._expect_mate = False
        # Both random decisions -- which mate unstranded normalization flips, and the
        # base substituted for a no-call -- are derived from `rng_seed` and the record's
        # GLOBAL index, never from a running stream. That is what makes the output
        # independent of `block_size`: batching cannot shift a position-derived draw.
        self._rng_seed = rng_seed
        self._records_written = 0

        self._label_defs = header.labels
        self._batch_labels: list[list] = [[] for _ in self._label_defs]

        if header.compression_method == COMPRESSION_ZSTD:
            # write_checksum embeds an xxhash64 content checksum in every frame,
            # verified automatically on decompress -- the format's first
            # protection against bit rot, at negligible cost, and invisible to
            # any reader (checksummed frames are ordinary frames).
            self._compressor = zstandard.ZstdCompressor(
                level=header.compression_level, write_checksum=True,
            )
        else:
            self._compressor = None

        # Trailer accumulators.  Every stat is tallied in ``_flush_block`` from
        # the columns the CODEC returned, never from the caller's arguments --
        # strand normalization ORs IS_RC into the flags column and an N-policy
        # can change stored lengths, so the input is simply not what the file
        # holds.  The fuzz suite recounts every generated file against these.
        self._shuffled = shuffled
        self._merged_in_process = merged_in_process
        self._len_char = {1: "B", 2: "H", 4: "I"}[header.seq_len_bytes]
        self._flag_counts = [0] * 256
        self._len_hist: dict[int, int] = {}
        self._len_hist_unpaired: dict[int, int] = {}
        self._n_bases = 0
        self._block_meta: list[tuple[int, int, int]] = []
        self._closed = False

        header_len = self._write_file_header()
        self._write_prologue(header_len)

    # -- context manager -----------------------------------------------------

    def __enter__(self) -> ZnaWriter:
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:  # noqa: ANN001
        # An exception on the way out already explains why the stream stopped mid
        # fragment; raising the pairing error on top of it would replace the cause
        # with its symptom.
        if exc_type is not None:
            self.abort()
            return
        self.close()

    def abort(self) -> None:
        """Close the writer WITHOUT signing the file: no sentinel, no trailer.

        For error paths that cannot use the context manager (``zna shuffle``
        holds many writers at once).  What was buffered is flushed so every
        block on disk stays whole, but no trailer is written: the trailer is
        the writer's signature that the encode finished, and a crashed encode
        must not produce a file that certifies.  ``records()`` still reads the
        blocks -- EOF termination is retained for exactly this shape of file --
        while ``zna inspect --verify`` refuses it.

        An open fragment cannot be flushed whole: if the buffer ends on a
        paired R1 whose mate never arrived, that one record is dropped rather
        than written -- a block ending mid-fragment would violate the very
        invariant "blocks stay whole" names.  Aborting also never raises the
        pairing error, so the exception that caused the abort keeps the stage.
        """
        if self._closed:
            return
        if self._expect_mate and self._batch_seqs:
            self._batch_seqs.pop()
            self._batch_flags.pop()
            for col in self._batch_labels:
                if col:
                    col.pop()
        self._expect_mate = False
        self._flush_block()
        self._closed = True

    def close(self) -> None:
        """Flush the final block, then sign the file: sentinel, trailer, footer.

        Raises if the stream ended on a paired R1 — its mate never arrived, so the
        file would close with half a fragment in its last block.  The buffered tail is
        *not* written in that case: every block already on disk is whole, and the
        caller has an exception telling it the file is short.

        Idempotent: a second call is a no-op, so ``close()`` inside a ``with``
        block does not write the trailer twice.
        """
        if self._closed:
            return
        if self._expect_mate:
            # The writer is finished either way: a second close() must be a
            # no-op, not a flush of the orphaned R1 -- that would write a
            # block ending mid-fragment AND sign it with a trailer whose
            # pair counts cannot balance.  The buffered tail dies with the
            # writer, exactly as this docstring promises, and the file on
            # disk -- whole blocks, no trailer -- reads but does not certify.
            self._expect_mate = False
            self._closed = True
            raise ValueError(
                "the record stream ended on a paired R1 whose R2 never arrived. "
                "Every fragment must be written whole, R1 immediately followed by R2; "
                "write the mate, or write the read unpaired (is_paired=False)."
            )
        self._flush_block()
        self._closed = True
        self._write_trailer()

    def _write_trailer(self) -> None:
        """Append sentinel + trailer payload + footer.  Append-only: works on pipes.

        The sentinel's ``comp_size`` covers the payload AND the 32-byte footer,
        deliberately: a 0.4.1 reader that walks past the last data block then
        consumes everything to EOF as one valid empty block and terminates
        cleanly, instead of parsing the footer as a garbage block header.
        (python-zstandard's one-shot ``decompress`` tolerates trailing bytes
        after a complete frame, so the footer riding inside ``comp_size`` does
        not break the old reader's inflate; ``uncomp_size`` stays honest, which
        is what bounds its ``max_output_size``.)  New readers never use the
        sentinel's sizes -- discovery is footer-first, or count==0 on pipes.
        """
        fc = self._flag_counts
        n_records = sum(fc)
        n_paired = sum(fc[f] for f in range(256) if f & 0x04)
        n_pairs, odd = divmod(n_paired, 2)
        n_unpaired = n_records - n_paired
        # The writer enforces R1-R2 adjacency, so an odd paired count cannot
        # happen; assert the derivation rather than trust it.
        assert odd == 0 and n_unpaired + 2 * n_pairs == n_records

        trailer: dict = {
            "trailer_schema": TRAILER_SCHEMA,
            "n_records": n_records,
            "n_bases": self._n_bases,
            "n_pairs": n_pairs,
            "n_unpaired": n_unpaired,
            "flag_counts": {str(f): c for f, c in enumerate(fc) if c},
            "length_histogram": {str(k): v for k, v in self._len_hist.items()},
            "length_histogram_unpaired": {
                str(k): v for k, v in self._len_hist_unpaired.items()},
            "blocks": {
                "data_start": self._data_start_offset,
                "comp_sizes": [m[0] for m in self._block_meta],
                "uncomp_sizes": [m[1] for m in self._block_meta],
                "records": [m[2] for m in self._block_meta],
            },
            "prologue_crc32": self._prologue_crc32,
        }

        raw = _canonical_json(trailer)
        payload = (self._compressor.compress(raw)
                   if self._compressor is not None else raw)
        self._fh.write(struct.pack(
            _BLOCK_HEADER_FMT,
            len(payload) + _FOOTER_SIZE, len(raw), 0, 0, 0,
        ))
        self._fh.write(payload)
        self._fh.write(struct.pack(
            _FOOTER_FMT,
            _FOOTER_MAGIC, _FOOTER_VERSION, 0,
            zlib.crc32(payload) & 0xFFFFFFFF,
            self._data_end,
        ))

    # -- public --------------------------------------------------------------

    def write_record(
        self, seq: str, is_paired: bool, is_read1: bool, is_read2: bool,
        labels: tuple | None = None, is_rc: bool = False,
        is_full_fragment: bool = False,
    ) -> None:
        """Buffer a single record.  Flushes at the next fragment boundary once full.

        A paired R1 must be followed immediately by its R2 — see :class:`ZnaWriter` —
        and anything else raises :exc:`ValueError`.

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

        if self._closed:
            raise ValueError(
                "this ZnaWriter is closed; a record written now would land "
                "after the trailer and be invisible to every reader."
            )
        flag = (
            (1 if is_read1 else 0) | (2 if is_read2 else 0) | (4 if is_paired else 0)
            | (8 if is_rc else 0) | (16 if is_full_fragment else 0)
        )
        if CLOSES_FRAGMENT[flag] != self._expect_mate:
            raise self._bad_fragment(flag)

        if self._label_defs:
            if labels is None or len(labels) != len(self._label_defs):
                raise ValueError(
                    f"Expected {len(self._label_defs)} label values, "
                    f"got {0 if labels is None else len(labels)}"
                )
            for i, val in enumerate(labels):
                self._batch_labels[i].append(val)

        self._batch_seqs.append(seq)
        self._batch_flags.append(flag)

        self._size_estimate += (seq_len // 4) + 1 + self._seq_len_bytes
        # A block may only end where a fragment does, so the size test is asked at
        # fragment boundaries and nowhere else.
        self._expect_mate = OPENS_FRAGMENT[flag]
        if not self._expect_mate and self._size_estimate >= self._block_size:
            self._flush_block()

    def write_copy(self, rec: ZnaRecord) -> None:
        """Buffer a record read from another ZNA file, flag byte verbatim.

        The write half of the copy path; see :meth:`ZnaReader.copy_records`.  Nothing
        here interprets the byte, so every bit survives — including
        ``IS_RC | IS_FULL_FRAGMENT`` together, which no boolean-shaped API can express,
        and bits a later format version may define.

        Requires ``preserve_normalization=True``.  Without it the codec ORs ``IS_RC``
        into the flag column as it reverse-complements (see ``_pycodec.encode_block``),
        which would both corrupt the copied byte and normalize already-normalized bases
        a second time.
        """
        if not self._preserve_normalization:
            raise ValueError(
                "write_copy() requires ZnaWriter(preserve_normalization=True). "
                "A copy carries each record's stored orientation in its flag byte; "
                "without the flag the writer would apply normalization again, on top "
                "of bases that already have it."
            )
        seq = rec[0]
        seq_len = len(seq)
        if seq_len > self._max_len:
            raise ValueError(
                f"Sequence length {seq_len} exceeds maximum {self._max_len} "
                f"allowed by header (seq_len_bytes={self._seq_len_bytes})"
            )
        flag = rec[1]
        if self._closed:
            raise ValueError(
                "this ZnaWriter is closed; a record written now would land "
                "after the trailer and be invisible to every reader."
            )
        if CLOSES_FRAGMENT[flag] != self._expect_mate:
            raise self._bad_fragment(flag)
        if self._label_defs:
            labels = rec[2]
            if len(labels) != len(self._label_defs):
                raise ValueError(
                    f"Expected {len(self._label_defs)} label values, got {len(labels)}"
                )
            for i, val in enumerate(labels):
                self._batch_labels[i].append(val)

        self._batch_seqs.append(seq)
        self._batch_flags.append(flag)
        self._size_estimate += (seq_len // 4) + 1 + self._seq_len_bytes
        self._expect_mate = OPENS_FRAGMENT[flag]
        if not self._expect_mate and self._size_estimate >= self._block_size:
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
        if self._closed:
            raise ValueError(
                "this ZnaWriter is closed; records written now would land "
                "after the trailer and be invisible to every reader."
            )
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
        # Bound in the loop, not looked up per record: the fragment contract is checked
        # on every record here, and this is the bulk path.
        opens_fragment = OPENS_FRAGMENT
        closes_fragment = CLOSES_FRAGMENT
        expect_mate = self._expect_mate

        if self._preserve_normalization:
            # Pass-through. Accepts:
            #   5-tuple  (..., is_rc)   <- records(with_rc=True)
            # A plain 4-tuple is only acceptable when the output is not
            # normalized, since there is then no orientation to carry; otherwise
            # accepting it would clear every IS_RC bit and lose the orientation
            # silently.
            #
            # A 6-tuple from records(with_ends=True) is REFUSED. It used to be
            # accepted, and this docstring used to advertise it as "a lossless ZNA ->
            # ZNA copy in one line", which was false: (has_start, has_end) has three
            # states and the flag pair it was reconstructing has four, so every
            # full-fragment record came out with IS_RC cleared. See the note above
            # ENDS_BY_FLAG. Copying carries the byte; use copy_records/write_copy.
            require_rc = self._header.strand_normalized
            for rec in records:
                seq = rec[0]
                is_paired = rec[1]
                is_read1 = rec[2]
                n_fields = len(rec)
                if n_fields > 5:
                    raise ValueError(
                        "write_records() no longer accepts the 6-tuples of "
                        "records(with_ends=True): (has_start, has_end) cannot express "
                        "IS_RC on a full-fragment record, so this silently dropped the "
                        "orientation of every merged read. Copy with "
                        "ZnaReader.copy_records() -> ZnaWriter.write_copy() instead, "
                        "which carries the flag byte verbatim."
                    )
                if n_fields > 4:
                    rc_bit = 8 if rec[4] else 0
                elif require_rc:
                    raise ValueError(
                        "preserve_normalization=True on a strand-normalized "
                        "header requires each record to carry its orientation, "
                        f"but got a {n_fields}-tuple. Copy the source with "
                        "ZnaReader.copy_records()."
                    )
                else:
                    rc_bit = 0
                seq_len = len(seq)
                if seq_len > max_len:
                    raise ValueError(
                        f"Sequence length {seq_len} exceeds maximum {max_len} "
                        f"allowed by header (seq_len_bytes={seq_len_bytes})"
                    )
                flag = (
                    (1 if is_read1 else 0)
                    | (2 if rec[3] else 0)
                    | (4 if is_paired else 0)
                    | rc_bit
                )
                if closes_fragment[flag] != expect_mate:
                    self._expect_mate = expect_mate
                    self._size_estimate = size_est
                    raise self._bad_fragment(flag)
                append_seq(seq)
                append_flag(flag)
                size_est += (seq_len >> 2) + 1 + seq_len_bytes
                expect_mate = opens_fragment[flag]
                if not expect_mate and size_est >= block_size:
                    self._size_estimate = size_est
                    flush()
                    size_est = self._size_estimate  # reset after flush
            self._size_estimate = size_est
            self._expect_mate = expect_mate
            return

        for seq, is_paired, is_read1, is_read2 in records:
            seq_len = len(seq)
            if seq_len > max_len:
                raise ValueError(
                    f"Sequence length {seq_len} exceeds maximum {max_len} "
                    f"allowed by header (seq_len_bytes={seq_len_bytes})"
                )
            flag = (
                (1 if is_read1 else 0)
                | (2 if is_read2 else 0)
                | (4 if is_paired else 0)
            )
            if closes_fragment[flag] != expect_mate:
                self._expect_mate = expect_mate
                self._size_estimate = size_est
                raise self._bad_fragment(flag)
            append_seq(seq)
            append_flag(flag)
            size_est += (seq_len >> 2) + 1 + seq_len_bytes
            expect_mate = opens_fragment[flag]
            if not expect_mate and size_est >= block_size:
                self._size_estimate = size_est
                flush()
                size_est = self._size_estimate  # reset after flush

        self._size_estimate = size_est
        self._expect_mate = expect_mate

    # -- private -------------------------------------------------------------

    def _bad_fragment(self, flag: int) -> ValueError:
        """The error for a record that cannot legally follow the previous one."""
        if self._expect_mate:
            return ValueError(
                f"expected the paired R2 completing the preceding R1, got flags="
                f"0x{flag:02x}. A fragment's reads must be written consecutively, R1 "
                f"immediately followed by R2, so that no block ever splits a fragment."
            )
        return ValueError(
            f"paired R2 (flags=0x{flag:02x}) with no R1 in front of it. A fragment's "
            f"reads must be written consecutively, R1 immediately followed by R2."
        )

    def _write_file_header(self) -> int:
        """Write the file header; return its byte length.

        The length is arithmetic, never ``tell()`` -- the writer must work on
        pipes, where every offset in the prologue/trailer bookkeeping is
        derived by running sums from this value.
        """
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

        return (_FILE_HEADER_SIZE + len(rg_bytes) + len(desc_bytes)
                + _LABEL_DEF_SIZE * len(self._header.labels))

    def _write_prologue(self, header_len: int) -> None:
        """Write the provenance prologue: facts known at encode START.

        A count-0 pseudo-block right after the header, so a streaming consumer
        learns what wrote the file -- and whether its order is shuffled --
        before decoding a single record.  Facts that need the whole encode go
        in the trailer instead; parameter echoes (thresholds, npolicy) go
        nowhere, because they are QC material, not file facts.

        A 0.4.1 reader decodes this as a valid empty block and moves on.
        """
        from . import __version__ as writer_version

        prov: dict = {
            "prologue_schema": PROLOGUE_SCHEMA,
            "writer_version": writer_version,
            "shuffled": bool(self._shuffled),
        }
        if self._merged_in_process:
            prov["merged_in_process"] = True

        raw = _canonical_json(prov)
        payload = (self._compressor.compress(raw)
                   if self._compressor is not None else raw)
        self._prologue_crc32 = zlib.crc32(payload) & 0xFFFFFFFF
        self._fh.write(struct.pack(
            _BLOCK_HEADER_FMT, len(payload), len(raw), 0, 0, 0,
        ))
        self._fh.write(payload)

        # Where block 0's header will sit; advanced by every flush, so at close
        # it is the sentinel offset the footer promises.
        self._data_start_offset = header_len + _BLOCK_HEADER_SIZE + len(payload)
        self._data_end = self._data_start_offset

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
                    self._rng_seed,
                    self._records_written,
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
                    rng_seed=self._rng_seed,
                    base_index=self._records_written,
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
                self._rng_seed,
                self._records_written,
            )
            labels_bytes = b""

        # 2. Compress — column order: flags, labels, lengths, sequences
        uncompressed = flags_bytes + labels_bytes + lengths_bytes + seqs_bytes
        compressed = (
            self._compressor.compress(uncompressed)
            if self._compressor is not None
            else uncompressed
        )

        # 2b. Tally trailer statistics from the columns the codec RETURNED.
        # ``flags_bytes`` carries IS_RC bits the caller never supplied, and
        # ``lengths_bytes`` is post-npolicy -- the file's truth, not the input's.
        # Everything below is C-speed: set/count over a column with few distinct
        # bytes, one struct.unpack, dict-of-int updates via Counter mechanics.
        assert count > 0, "count == 0 is reserved for metadata pseudo-blocks"
        fc = self._flag_counts
        for b in set(flags_bytes):
            fc[b] += flags_bytes.count(b)

        lengths = struct.unpack(f"<{count}{self._len_char}", lengths_bytes)
        self._n_bases += sum(lengths)
        lo, hi = min(lengths), max(lengths)
        if lo == hi:
            # Uniform read length -- the overwhelmingly common shape of real
            # libraries -- collapses the histogram update to one dict add.
            self._len_hist[lo] = self._len_hist.get(lo, 0) + count
        else:
            _count_into(self._len_hist, lengths)
        selector = flags_bytes.translate(_UNPAIRED_SELECTOR)
        n_unpaired = selector.count(1)
        if n_unpaired:
            if lo == hi:
                self._len_hist_unpaired[lo] = (
                    self._len_hist_unpaired.get(lo, 0) + n_unpaired)
            elif n_unpaired == count:
                _count_into(self._len_hist_unpaired, lengths)
            else:
                _count_into(self._len_hist_unpaired,
                            compress(lengths, selector))

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
        self._block_meta.append((len(compressed), len(uncompressed), count))
        self._data_end += _BLOCK_HEADER_SIZE + len(compressed)

        # 4. Reset
        # Advance the global record index BEFORE clearing, so the next block's
        # position-derived draws continue where this one left off.
        self._records_written += count
        self._batch_seqs.clear()
        self._batch_flags.clear()
        self._size_estimate = 0
        for col in self._batch_labels:
            col.clear()


# ---------------------------------------------------------------------------
# Reader
# ---------------------------------------------------------------------------

class BlockInfo(NamedTuple):
    """One entry of :meth:`ZnaReader.block_index`.

    *offset* is the byte position of the block header, so a future reader can
    seek straight to it.  *n_records* comes from that header, which is why an
    index costs no decompression to build.
    """

    index: int
    offset: int
    n_records: int
    comp_size: int
    uncomp_size: int


class ZnaReader:
    """Read records from a ZNA file.

    Yields ``(sequence, is_paired, is_read1, is_read2)`` tuples.
    """

    def __init__(self, fh: BinaryIO) -> None:
        self._fh = fh
        self._header = self._read_file_header()

        # Consume the provenance prologue, if there is one, so that every walk
        # loop after this point is prologue-blind and ``count == 0`` has
        # exactly one meaning to them: the trailer sentinel, i.e. end of data.
        # On a 0.4.1-era file the probed 20 bytes are the first data block's
        # header; a seekable stream rewinds, a pipe stashes them for the first
        # walk read.
        self._provenance: ZnaProvenance | None = None
        self._pending_block_header = b""
        self._trailer: ZnaTrailer | None = None
        self._trailer_loaded = False
        probe = fh.read(_BLOCK_HEADER_SIZE)
        if len(probe) == _BLOCK_HEADER_SIZE:
            comp_size, uncomp_size, count, _f, _l = struct.unpack(
                _BLOCK_HEADER_FMT, probe)
            if count == 0:
                self._provenance = self._parse_prologue(comp_size, uncomp_size)
            else:
                self._pending_block_header = probe
        elif probe:
            # 1-19 bytes where a block header should be: stash them so the
            # walks raise their usual incomplete-header error instead of this
            # constructor swallowing the damage.
            self._pending_block_header = probe

        #: Position of the first data block, so :meth:`block_index` can rewind
        #: to it no matter how far a caller has already read.
        try:
            if self._pending_block_header:
                fh.seek(-len(self._pending_block_header), 1)
                self._pending_block_header = b""
            self._data_start = fh.tell()
        except (AttributeError, OSError):
            self._data_start = None

    def _parse_prologue(self, comp_size: int, uncomp_size: int) -> ZnaProvenance:
        payload = self._read_exact(comp_size)
        try:
            if self._header.compression_method == COMPRESSION_ZSTD:
                raw = zstandard.ZstdDecompressor().decompress(
                    payload, max_output_size=uncomp_size)
            else:
                raw = payload
            d = json.loads(raw)
            return ZnaProvenance(
                writer_version=d["writer_version"],
                shuffled=bool(d["shuffled"]),
                merged_in_process=bool(d.get("merged_in_process", False)),
                crc32=zlib.crc32(payload) & 0xFFFFFFFF,
                raw=d,
            )
        except (ValueError, KeyError, TypeError, zstandard.ZstdError) as e:
            raise ValueError(
                f"Corrupt provenance prologue: {e}. The file's data blocks may "
                f"be intact, but its provenance cannot be trusted."
            ) from e

    @property
    def provenance(self) -> ZnaProvenance | None:
        """Facts known at encode start, or ``None`` for files written by zna < 0.5.

        Available on pipes too — the prologue sits before the first block, which
        is the point of putting it there.
        """
        return self._provenance

    @property
    def header(self) -> ZnaHeader:
        return self._header

    @property
    def trailer(self) -> ZnaTrailer | None:
        """Facts derived from the complete encode, or ``None`` when absent.

        ``None`` means the stream is a pipe (no seek to the footer), or the
        file carries no trailer — written by zna < 0.5, an aborted encode, or
        truncated; ``zna inspect --verify`` treats all of those as uncertified.
        A trailer that is *present but corrupt* (bad footer version, CRC
        mismatch, unparseable payload) raises ``ValueError`` instead: damage
        must not be indistinguishable from age.

        Leaves the read position unchanged.  Cached after the first call.
        """
        if self._trailer_loaded:
            return self._trailer
        if self._data_start is None:
            self._trailer_loaded = True
            return None

        fh = self._fh
        resume = fh.tell()
        try:
            fh.seek(0, 2)
            file_size = fh.tell()
            if file_size < self._data_start + _FOOTER_SIZE:
                self._trailer_loaded = True
                return None
            fh.seek(file_size - _FOOTER_SIZE)
            magic, footer_version, _reserved, crc, sentinel_offset = (
                struct.unpack(_FOOTER_FMT, fh.read(_FOOTER_SIZE))
            )
            if magic != _FOOTER_MAGIC:
                self._trailer_loaded = True
                return None
            if footer_version != _FOOTER_VERSION:
                raise ValueError(
                    f"Unknown trailer footer version {footer_version}; this "
                    f"build reads version {_FOOTER_VERSION}. Written by a "
                    f"newer zna?"
                )
            payload_start = sentinel_offset + _BLOCK_HEADER_SIZE
            payload_end = file_size - _FOOTER_SIZE
            if not (self._data_start <= sentinel_offset
                    and sentinel_offset + _BLOCK_HEADER_SIZE <= payload_end):
                raise ValueError(
                    f"Trailer footer points at sentinel offset "
                    f"{sentinel_offset}, outside the file's data region — "
                    f"corrupt trailer."
                )
            fh.seek(sentinel_offset)
            _c, uncomp_size, count, _f, _l = struct.unpack(
                _BLOCK_HEADER_FMT, fh.read(_BLOCK_HEADER_SIZE))
            if count != 0:
                raise ValueError(
                    f"No sentinel at offset {sentinel_offset} where the "
                    f"footer promised one — corrupt trailer."
                )
            payload = self._read_exact(payload_end - payload_start)
            if zlib.crc32(payload) & 0xFFFFFFFF != crc:
                raise ValueError(
                    "Trailer payload CRC mismatch — the file is corrupt."
                )
            if self._header.compression_method == COMPRESSION_ZSTD:
                raw = zstandard.ZstdDecompressor().decompress(
                    payload, max_output_size=uncomp_size)
            else:
                raw = payload
            try:
                d = json.loads(raw)
                blocks = d["blocks"]
                if blocks["data_start"] != self._data_start:
                    raise ValueError(
                        f"trailer says data starts at {blocks['data_start']}, "
                        f"but the reader found it at {self._data_start}"
                    )
                self._trailer = ZnaTrailer(
                    n_records=d["n_records"],
                    n_bases=d["n_bases"],
                    n_pairs=d["n_pairs"],
                    n_unpaired=d["n_unpaired"],
                    flag_counts={int(k): v for k, v in d["flag_counts"].items()},
                    length_histogram={
                        int(k): v for k, v in d["length_histogram"].items()},
                    length_histogram_unpaired={
                        int(k): v
                        for k, v in d["length_histogram_unpaired"].items()},
                    data_start=blocks["data_start"],
                    block_comp_sizes=blocks["comp_sizes"],
                    block_uncomp_sizes=blocks["uncomp_sizes"],
                    block_records=blocks["records"],
                    prologue_crc32=d["prologue_crc32"],
                    raw=d,
                )
            except (KeyError, TypeError, ValueError) as e:
                raise ValueError(f"Corrupt trailer payload: {e}") from e
            self._trailer_loaded = True
            return self._trailer
        finally:
            fh.seek(resume)

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
            express).

            It is a **projection, not an encoding**: ``IS_FULL_FRAGMENT`` and
            ``IS_RC|IS_FULL_FRAGMENT`` both report ``(True, True)``, correctly,
            because both records do have two real fragment ends.  So the pair
            does *not* round-trip back to the flags and must never be fed to a
            writer — use :meth:`copy_records` to copy records between files.
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
                "for a record that is not a full fragment, with_ends already encodes "
                "IS_RC as (has_end and not has_start). For a full-fragment record it "
                "does not -- the ends are (True, True) either way -- so if you need "
                "the orientation of those, read the flag byte with copy_records()."
            )
        return self._iter_records(restore_strand, with_rc, with_ends)

    def copy_records(self) -> Iterator[ZnaRecord]:
        """Yield every record as stored, for writing to another ZNA file.

        The lossless counterpart to :meth:`records`.  ``records`` returns *views* —
        ``(is_paired, is_read1, is_read2)``, ``is_rc``, ``(has_start, has_end)`` — each
        of which is a projection chosen for a consumer, and none of which can carry the
        whole flag byte back to a writer.  This yields the byte.

        Pair it with :meth:`ZnaWriter.write_copy`::

            with open(src, "rb") as fin, open(dst, "wb") as fout:
                reader = ZnaReader(fin)
                with ZnaWriter(fout, reader.header, preserve_normalization=True) as w:
                    for rec in reader.copy_records():
                        w.write_copy(rec)

        Sequences come back in their **stored** orientation, which is the only
        orientation a copy may write: re-deriving it would apply normalization a second
        time.  Hence ``preserve_normalization=True`` on the writer, which
        :meth:`~ZnaWriter.write_copy` requires.
        """
        if self._header.labels:
            for seqs, flags, cols in self.blocks(labels=True):
                for rec in zip(seqs, flags, zip(*cols)):
                    yield ZnaRecord(*rec)
        else:
            for seqs, flags in self.blocks(labels=False):
                for seq, flag in zip(seqs, flags):
                    yield ZnaRecord(seq, flag, ())

    def block_index(self) -> list[BlockInfo]:
        """Return one :class:`BlockInfo` per block, without decompressing any.

        **The file header stores no record or block count** — it holds only the
        format version, sequence-length width, strand flags, compression
        settings and label schema.  Every *block* header does carry its own
        record count, so the totals are recovered by stepping through the block
        chain, seeking over each payload rather than reading it.

        That walk is cheap enough to do at open time: measured 2.3 microseconds
        per block, so 1.4 ms for a 38 MB / 611-block file, against 89 ms to
        reach the same counts by decoding.  A stored index would only be worth
        it to avoid touching the file at all — which is a remote-storage
        problem, not a local one.

        Use it to size a subsample before reading:

            idx = reader.block_index()
            total = sum(b.n_records for b in idx)
            keep = [b.index for b in random.sample(idx, k)]
            for seqs, flags in reader.blocks(indices=keep):
                ...

        Blocks are flushed on an estimated byte size, so their record counts are
        near-uniform for fixed-length reads and vary with variable-length ones —
        which is exactly why this returns per-block counts rather than an
        average.  The final block is a partial and is usually much smaller.

        Requires a seekable stream, and leaves the read position unchanged.
        """
        if self._data_start is None:
            raise TypeError(
                "block_index() requires a seekable stream; this reader was "
                "opened on a pipe or socket."
            )

        # A trailer stores the index the scan would recompute, so a
        # trailer-bearing file answers in one seek instead of one per block.
        # The scan is retained forever: it is what ``zna inspect --verify``
        # cross-checks the stored index against.
        trailer = self.trailer
        if trailer is not None:
            return [
                BlockInfo(index=i, offset=off, n_records=n,
                          comp_size=c, uncomp_size=u)
                for i, (off, c, u, n) in enumerate(zip(
                    trailer.block_offsets(),
                    trailer.block_comp_sizes,
                    trailer.block_uncomp_sizes,
                    trailer.block_records,
                ))
            ]
        return self.scan_block_index()

    def scan_block_index(self) -> list[BlockInfo]:
        """:meth:`block_index` by walking the block headers, trailer ignored.

        One 20-byte read and one seek per block.  This is the ground truth the
        stored index is checked against — a trailer could describe a different
        file than the one it rides in, and only a walk can say so.
        """
        if self._data_start is None:
            raise TypeError(
                "scan_block_index() requires a seekable stream; this reader "
                "was opened on a pipe or socket."
            )
        fh = self._fh
        resume = fh.tell()
        try:
            fh.seek(0, 2)
            file_size = fh.tell()
            fh.seek(self._data_start)
            index: list[BlockInfo] = []
            while True:
                offset = fh.tell()
                block_header_data = fh.read(_BLOCK_HEADER_SIZE)
                if not block_header_data:
                    break
                if len(block_header_data) < _BLOCK_HEADER_SIZE:
                    raise EOFError(
                        f"Incomplete block header at offset {offset}. Expected "
                        f"{_BLOCK_HEADER_SIZE} bytes, got {len(block_header_data)}"
                    )
                comp_size, uncomp_size, count, _flags_size, _lengths_size = (
                    struct.unpack(_BLOCK_HEADER_FMT, block_header_data)
                )
                if count == 0:
                    break  # the trailer sentinel: end of data
                if offset + _BLOCK_HEADER_SIZE + comp_size > file_size:
                    # This used to be a silent seek past EOF followed by a clean
                    # empty read: a truncated or garbage-tailed file returned a
                    # phantom BlockInfo and success.
                    raise ValueError(
                        f"block header at offset {offset} claims a "
                        f"{comp_size}-byte payload, but the file ends at "
                        f"{file_size} — truncated, or trailing garbage."
                    )
                index.append(BlockInfo(
                    index=len(index), offset=offset, n_records=count,
                    comp_size=comp_size, uncomp_size=uncomp_size,
                ))
                fh.seek(comp_size, 1)
            return index
        finally:
            fh.seek(resume)

    def blocks(
        self, stride: int = 1, offset: int = 0, restore_strand: bool = False,
        indices: Iterable[int] | None = None, labels: bool | None = None,
    ) -> Iterator[Tuple[list, bytes] | Tuple[list, bytes, tuple]]:
        """Yield ``(sequences, flags)`` once per block, columnar.

        *sequences* is a ``list[str]``; *flags* is the block's raw flag column,
        one byte per record, positionally aligned with *sequences*.  Decode a
        byte with :class:`ZnaRecordFlags`, or index :data:`FLAG_FIELDS` for the
        ``(is_paired, is_read1, is_read2)`` triple the record API yields.

        This is the batch form for consumers that process a whole block at a
        time — a training loader, say.  It skips the per-record tuple entirely,
        which is most of what is left in :meth:`records` once the sequence itself
        is cheap to produce.

        **A block holds whole fragments.**  Paired reads appear consecutively, R1
        then R2, and a fragment never straddles a block boundary — so a worker
        handed a block never sees half a molecule whose other half went to some
        other worker.  That is what makes sharding here safe to do at all, and it
        is a guarantee of the format, enforced by :class:`ZnaWriter` on write, not
        a property that happens to hold for a particular file.

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

        ``indices`` selects an explicit set of block numbers instead, as
        produced by :meth:`block_index`, and is mutually exclusive with
        ``stride``/``offset``.  Use it when the shard fraction is not a unit
        fraction, or when successive passes over one file should see *different*
        blocks: ``stride`` admits only ``stride`` distinct phases, so training
        several epochs over a file at ``stride=4`` would revisit the same four
        subsets.  Non-selected blocks are seeked over either way.

        ``labels`` decides what happens to a labeled file's label columns:

        * unset (the default) — a labeled file **raises**.  Handing back
          sequences while quietly discarding the label columns is not a decision
          this method should make on the caller's behalf.
        * ``labels=False`` — the columns are skipped.  Not silent: the caller
          asked for it.  This is what a consumer that only wants sequence should
          pass; it is worth about 1.9x on a three-column file.
        * ``labels=True`` — yields a third element, one value-tuple per label
          column in header order, each ``count`` long.  ``len(label_columns)``
          always equals ``len(header.labels)``, so an unlabeled file yields
          ``()`` rather than erroring.

        Note what ``labels=False`` does and does not save.  A block payload is a
        single zstd frame holding all four columns, so the label bytes are
        inflated either way; what is skipped is unpacking them into Python
        objects, which is where the cost actually is.
        """
        if labels is None and self._header.labels:
            raise NotImplementedError(
                "blocks() needs an explicit labels= on a labeled file "
                f"({len(self._header.labels)} label column(s) defined). "
                "Pass labels=False to skip the label columns, labels=True to "
                "receive them columnar, or use records() to get them per record."
            )
        if indices is not None:
            if stride != 1 or offset != 0:
                raise ValueError(
                    "indices is mutually exclusive with stride/offset: pass the "
                    "block numbers you want, or a stride, not both."
                )
            selected = frozenset(indices)
            for i in selected:
                if not isinstance(i, (int,)) or isinstance(i, bool) or i < 0:
                    raise ValueError(
                        f"indices must be non-negative block numbers, got {i!r}"
                    )
        else:
            selected = None
            if stride < 1:
                raise ValueError(f"stride must be >= 1, got {stride}")
            if not 0 <= offset < stride:
                raise ValueError(f"offset must be in [0, {stride}), got {offset}")
        return self._iter_blocks(stride, offset, restore_strand, selected,
                                 bool(labels))

    def _iter_blocks(
        self, stride: int, offset: int, restore_strand: bool,
        selected: frozenset | None, want_labels: bool,
    ) -> Iterator[Tuple[list, bytes] | Tuple[list, bytes, tuple]]:
        fh = self._fh
        fh_read = fh.read
        len_bytes = self._header.seq_len_bytes
        compression_method = self._header.compression_method
        needs_restore = restore_strand and self._header.strand_normalized
        decode_seqs = _codec.decode_block_sequences

        # The payload is flags | labels | lengths | sequences.  The label column
        # widths are needed even when the caller does not want the values,
        # because without them the lengths and sequence streams start at the
        # wrong offset -- which decodes silently, into garbage.
        label_defs = self._header.labels
        label_col_sizes = [label_bytes_per_record(ld) for ld in label_defs]
        label_bytes_per_rec = sum(label_col_sizes)
        label_formats = [ld.dtype.struct_ch for ld in label_defs]

        if compression_method == COMPRESSION_ZSTD:
            dctx = zstandard.ZstdDecompressor()

        index = 0
        yielded = 0
        # On a pipe over a 0.4.1-era file, __init__'s prologue probe consumed
        # the first block header; it is replayed here.
        pending = self._pending_block_header
        self._pending_block_header = b""
        while True:
            block_header_data = pending or fh_read(_BLOCK_HEADER_SIZE)
            pending = b""
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
            if count == 0:
                # The trailer sentinel: end of data.  Deliberately not counted
                # as a block, so the empty-shard warning at loop exit reports
                # the number of DATA blocks a stride could have matched.  The
                # payload is consumed so a second walk finds a clean EOF.
                fh_read(comp_size)
                break

            index += 1
            this_block = index - 1
            wanted = (this_block in selected if selected is not None
                      else this_block % stride == offset)
            if not wanted:
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

            labels_end = flags_size + label_bytes_per_rec * count
            lengths_end = labels_end + lengths_size
            flags_stream = block_data[:flags_size]
            lengths_stream = block_data[labels_end:lengths_end]
            seqs_stream = block_data[lengths_end:]

            yielded += 1
            sequences = decode_seqs(flags_stream, lengths_stream, seqs_stream,
                                    len_bytes, count, needs_restore)
            if not want_labels:
                yield sequences, flags_stream
                continue

            columns = []
            offset_ = flags_size
            for col_size, fmt in zip(label_col_sizes, label_formats):
                col_bytes = col_size * count
                columns.append(struct.unpack(
                    f"<{count}{fmt}",
                    block_data[offset_:offset_ + col_bytes],
                ))
                offset_ += col_bytes
            yield sequences, flags_stream, tuple(columns)

        if yielded == 0 and index > 0 and (stride > 1 or selected is not None):
            # Silence here is indistinguishable from an empty file, and in a
            # training loader it is an idle worker nobody notices.  Checked at
            # loop exit so both terminations -- EOF on a 0.4.1-era file, the
            # trailer sentinel on a 0.5 one -- report it.
            what = (f"indices={sorted(selected)[:8]}..."
                    if selected is not None
                    else f"stride={stride}, offset={offset}")
            warnings.warn(
                f"blocks({what}) matched none of this file's {index} "
                f"block(s), so this shard has no data. Shards are whole "
                f"blocks: write the file with a smaller block_size, or "
                f"select fewer.",
                RuntimeWarning, stacklevel=3,
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

        # On a pipe over a 0.4.1-era file, __init__'s prologue probe consumed
        # the first block header; it is replayed here.
        pending = self._pending_block_header
        self._pending_block_header = b""

        while True:
            block_header_data = pending or fh_read(_BLOCK_HEADER_SIZE)
            pending = b""
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
            if count == 0:
                # The trailer sentinel: end of data.  (EOF termination above is
                # retained for 0.4.1-era files and aborted encodes.)  The
                # payload is consumed -- comp_size runs to EOF -- so a second
                # walk on this reader finds a clean end, not trailer bytes
                # misread as a block header.  Short reads are fine here: the
                # records were all delivered, and trailer damage is
                # ``inspect --verify``'s business, not iteration's.
                fh_read(comp_size)
                break

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
            raise ValueError(
                f"Unsupported ZNA version: {ver} (this build reads version "
                f"{_VERSION}). Version 3 made blocks fragment-complete — a "
                f"fragment's reads are consecutive and never span a block "
                f"boundary — which version 2 files do not satisfy, so they "
                f"cannot be read as if they did. Re-encode from FASTQ."
            )

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
