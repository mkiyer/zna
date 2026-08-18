"""Bucket-shuffle implementation for ZNA files.

Algorithm — *bucket shuffle* with bounded memory:

1.  **Distribute:** Read the input ZNA sequentially.  Each *shuffle
    unit* (an unpaired record, or a consecutive R1+R2 pair) is
    randomly assigned to one of *K* temporary bucket files on disk.
2.  **Collect:** Read each bucket into memory, Fisher-Yates shuffle
    its units in-place, and append to the output ZNA.

Because the assignment in step 1 is random, the concatenation of
independently shuffled buckets yields a statistically uniform
permutation over all units.

*K* is derived from ``buffer_bytes`` so that each bucket fits
comfortably in memory during step 2.

A shuffle is a pure permutation, so orientation must be copied rather than
re-derived: every writer here runs in pass-through mode
(``preserve_normalization=True``) and carries each record's ``IS_RC`` flag
through verbatim.  Re-deriving it would reverse-complement a second time on the
bucket pass and a third time on the output pass, leaving the sequences plausible
but ``IS_RC`` uncorrelated with the real fragment boundary.
"""
from __future__ import annotations

import os
import random
import struct
import sys
import tempfile
from dataclasses import replace
from typing import List, Tuple

import zstandard

from .core import (
    COMPRESSION_ZSTD, ZnaHeader, ZnaWriter, ZnaReader, FLAG_FIELDS,
    _BLOCK_HEADER_FMT, _BLOCK_HEADER_SIZE,
)

#: ZSTD level for the throwaway bucket files.  They exist for one read.
BUCKET_ZSTD_LEVEL = 1

#: ``flag byte -> 1 if the record starts a shuffle unit``.  A unit is an unpaired
#: record or a paired R1; a paired R2 continues the unit its R1 opened.
_IS_UNIT_START = bytes(
    0 if (f & 0x05) == 0x04 else 1 for f in range(256)
)


def _scan_counts(fh, header) -> Tuple[int, int]:
    """Return ``(n_units, n_records)`` without decoding a single base.

    *fh* must be positioned at the first block.

    The counting pass only needs each record's flag byte.  Flags are the first
    column in the block payload, and a zstd frame can be read as a prefix, so
    this decompresses ``flags_size`` bytes per block and stops — rather than
    inflating the whole payload and building a Python string and tuple for every
    record just to look at three bits of each.  Record counts come free from the
    block headers.

    Translating the flag column through a 256-entry table and counting the
    result keeps the per-record work in C as well.
    """
    dctx = (zstandard.ZstdDecompressor()
            if header.compression_method == COMPRESSION_ZSTD else None)
    n_units = 0
    n_records = 0

    while True:
        block_header = fh.read(_BLOCK_HEADER_SIZE)
        if not block_header or len(block_header) < _BLOCK_HEADER_SIZE:
            break
        comp_size, _uncomp_size, count, flags_size, _lengths_size = struct.unpack(
            _BLOCK_HEADER_FMT, block_header
        )
        if count == 0:
            break  # the trailer sentinel: end of data
        n_records += count

        if dctx is not None:
            payload = fh.read(comp_size)
            if len(payload) != comp_size:
                raise EOFError("Unexpected EOF while scanning ZNA blocks")
            flags = dctx.stream_reader(payload).read(flags_size)
        else:
            flags = fh.read(flags_size)
            fh.seek(comp_size - flags_size, 1)

        if len(flags) != flags_size:
            raise EOFError("Unexpected EOF while reading a block's flags column")
        n_units += flags.translate(_IS_UNIT_START).count(1)

    return n_units, n_records


# Type alias for a single record tuple (without or with labels).
Record = tuple


def shuffle_zna(
    input_path: str,
    output_path: str,
    *,
    seed: int = 42,
    buffer_bytes: int = 1 << 30,  # 1 GiB
    block_size: int = 4 * 1024 * 1024,  # 4 MiB
    tmp_dir: str | None = None,
    quiet: bool = False,
) -> None:
    """Shuffle the records in a ZNA file with bounded memory.

    Parameters
    ----------
    input_path:
        Path to the source ``.zna`` file.
    output_path:
        Path for the shuffled ``.zna`` output.
    seed:
        Random seed for reproducibility.
    buffer_bytes:
        Approximate maximum memory consumed per bucket when loaded.
    block_size:
        Block size (in bytes) for the output ZNA writer.
    tmp_dir:
        Directory for temporary bucket files.  Defaults to the
        system temp directory.
    quiet:
        Suppress progress messages on *stderr*.
    """
    if not os.path.isfile(input_path):
        raise FileNotFoundError(f"Input file not found: {input_path}")
    if os.path.abspath(input_path) == os.path.abspath(output_path):
        raise ValueError("Input and output files must be different.")

    rng = random.Random(seed)

    # ── Pass 0: scan for record count and header ──────────────────────
    if not quiet:
        print(f"[ZNA] Scanning {input_path} ...", file=sys.stderr)

    with open(input_path, "rb") as f:
        reader = ZnaReader(f)
        in_header = reader.header
        # ZnaReader has consumed the file header and the provenance
        # prologue (if any), so f is at the first data block.
        n_units, n_records = _scan_counts(f, in_header)

    if n_units == 0:
        raise ValueError("Input file contains no records.")

    # Estimate bytes per unit (rough: each record ≈ avg_bases/4 + overhead)
    input_size = os.path.getsize(input_path)
    bytes_per_unit = max(input_size / n_units, 64)

    # Choose K so each bucket ≈ buffer_bytes when loaded
    units_per_bucket = max(1, int(buffer_bytes / bytes_per_unit))
    n_buckets = max(1, (n_units + units_per_bucket - 1) // units_per_bucket)
    n_buckets = min(n_buckets, n_units, 4096)

    if not quiet:
        print(f"[ZNA] {n_records:,} records, {n_units:,} shuffle units", file=sys.stderr)
        print(
            f"[ZNA] Using {n_buckets} buckets "
            f"(buffer {buffer_bytes // (1024 * 1024)} MB)",
            file=sys.stderr,
        )

    # ── Pass 1: distribute units into bucket files ────────────────────
    if not quiet:
        print("[ZNA] Distributing to buckets ...", file=sys.stderr)

    out_header = ZnaHeader(
        read_group=in_header.read_group,
        description=in_header.description,
        seq_len_bytes=in_header.seq_len_bytes,
        strand_specific=in_header.strand_specific,
        read1_antisense=in_header.read1_antisense,
        read2_antisense=in_header.read2_antisense,
        strand_normalized=in_header.strand_normalized,
        compression_method=in_header.compression_method,
        compression_level=in_header.compression_level,
        labels=in_header.labels,
    )

    # Bucket files are written once, read back immediately, and deleted with the
    # temp directory, so compressing them at the *source file's* level spends
    # real CPU buying a ratio nothing ever benefits from.  Level 1 keeps the
    # temp footprint in the same ballpark for a fraction of the cost.  The
    # output header is untouched: the file the user keeps still gets their level.
    if out_header.compression_method == COMPRESSION_ZSTD:
        bucket_header = replace(out_header, compression_level=BUCKET_ZSTD_LEVEL)
    else:
        bucket_header = out_header

    has_labels = len(in_header.labels) > 0
    # Bound once: both passes below index this per record to recover the pairing
    # triple from the flag byte, and a global lookup per record is what it replaces.
    flag_fields = FLAG_FIELDS

    with tempfile.TemporaryDirectory(dir=tmp_dir, prefix="zna_shuffle_") as tmp_path:
        bucket_paths = [
            os.path.join(tmp_path, f"bucket_{i:04d}.zna") for i in range(n_buckets)
        ]
        bucket_fhs = [open(p, "wb") for p in bucket_paths]
        bucket_writers = [
            ZnaWriter(fh, bucket_header, block_size=block_size,
                      preserve_normalization=True)
            for fh in bucket_fhs
        ]

        aborted = False
        try:
            with open(input_path, "rb") as f:
                reader = ZnaReader(f)
                distributed = 0
                _pending_bucket = 0

                for rec in reader.copy_records():
                    is_paired, is_read1, is_read2 = flag_fields[rec.flags]

                    if is_paired and is_read1:
                        # Hold the bucket so R2 lands in it too: a pair must not be
                        # split across buckets or it cannot be kept adjacent.
                        _pending_bucket = rng.randrange(n_buckets)
                        bucket_writers[_pending_bucket].write_copy(rec)
                    elif is_paired and is_read2:
                        bucket_writers[_pending_bucket].write_copy(rec)
                        distributed += 1
                    else:
                        bucket_writers[rng.randrange(n_buckets)].write_copy(rec)
                        distributed += 1

                    if (
                        distributed % 1_000_000 == 0
                        and distributed > 0
                        and not quiet
                    ):
                        print(
                            f"      Distributed {distributed // 1_000_000}M units ...",
                            end="\r",
                            file=sys.stderr,
                        )
        except BaseException:
            aborted = True
            raise
        finally:
            # On the error path, abort() rather than close(): a close() here
            # would (a) raise the dangling-R1 error over the REAL exception if
            # the interrupt landed between an R1 and its R2, and (b) sign a
            # half-distributed bucket with a trailer.  Buckets die with the
            # temp directory either way, but exceptions must keep the stage.
            for w in bucket_writers:
                (w.abort if aborted else w.close)()
            for fh in bucket_fhs:
                fh.close()

        if not quiet:
            print("\n[ZNA] Shuffling and writing output ...", file=sys.stderr)

        # ── Pass 2: read each bucket, shuffle, write to output ────────
        bucket_order = list(range(n_buckets))
        rng.shuffle(bucket_order)

        written = 0
        with open(output_path, "wb") as out_fh:
            # ``shuffled=True`` goes into the output's provenance prologue:
            # this function IS the attestation khorana's C1 used to take on
            # operator word.  Bucket writers stay unstamped -- they are
            # deleted with the temp directory.
            #
            # The writer is driven as a context manager, NOT close()d in a
            # finally: a KeyboardInterrupt mid-pass would otherwise leave a
            # partial file at the USER'S OUTPUT PATH carrying a complete,
            # self-consistent trailer and shuffled=True -- a certified file
            # silently missing most of its records, which is precisely the
            # failure class the trailer exists to eliminate.  __exit__ writes
            # no trailer on an exception, so an aborted shuffle leaves a file
            # that reads but refuses certification.
            out_writer = ZnaWriter(out_fh, out_header, block_size=block_size,
                                   preserve_normalization=True, shuffled=True)
            try:
                for bi in bucket_order:
                    bp = bucket_paths[bi]
                    if not os.path.isfile(bp) or os.path.getsize(bp) == 0:
                        continue

                    # Read entire bucket into memory
                    units: List[List[Record]] = []
                    with open(bp, "rb") as bfh:
                        breader = ZnaReader(bfh)
                        current_unit: List[Record] = []
                        for rec in breader.copy_records():
                            is_paired, is_read1, is_read2 = flag_fields[rec.flags]
                            if is_paired and is_read1:
                                if current_unit:
                                    units.append(current_unit)
                                current_unit = [rec]
                            elif is_paired and is_read2:
                                current_unit.append(rec)
                                units.append(current_unit)
                                current_unit = []
                            else:
                                if current_unit:
                                    units.append(current_unit)
                                units.append([rec])
                                current_unit = []
                        if current_unit:
                            units.append(current_unit)

                    rng.shuffle(units)

                    for unit in units:
                        for rec in unit:
                            out_writer.write_copy(rec)
                    written += len(units)

                    if not quiet:
                        print(
                            f"      Written {written:,} / {n_units:,} units ...",
                            end="\r",
                            file=sys.stderr,
                        )
            except BaseException:
                out_writer.abort()
                raise
            else:
                out_writer.close()

    return written, n_records
