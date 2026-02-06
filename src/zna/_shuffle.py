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
"""
from __future__ import annotations

import os
import random
import sys
import tempfile
from typing import List, Tuple

from .core import ZnaHeader, ZnaWriter, ZnaReader


# Type alias for a single record tuple.
Record = Tuple[str, bool, bool, bool]


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

        n_units = 0
        n_records = 0
        for _seq, is_paired, is_read1, _is_read2 in reader.records():
            n_records += 1
            if not (is_paired and not is_read1):
                n_units += 1

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
        compression_method=in_header.compression_method,
        compression_level=in_header.compression_level,
    )

    with tempfile.TemporaryDirectory(dir=tmp_dir, prefix="zna_shuffle_") as tmp_path:
        bucket_paths = [
            os.path.join(tmp_path, f"bucket_{i:04d}.zna") for i in range(n_buckets)
        ]
        bucket_fhs = [open(p, "wb") for p in bucket_paths]
        bucket_writers = [
            ZnaWriter(fh, out_header, block_size=block_size) for fh in bucket_fhs
        ]

        try:
            with open(input_path, "rb") as f:
                reader = ZnaReader(f)
                distributed = 0
                _pending_bucket = 0

                for seq, is_paired, is_read1, is_read2 in reader.records():
                    if is_paired and is_read1:
                        _pending_bucket = rng.randrange(n_buckets)
                        bucket_writers[_pending_bucket].write_record(
                            seq, is_paired, is_read1, is_read2
                        )
                    elif is_paired and is_read2:
                        bucket_writers[_pending_bucket].write_record(
                            seq, is_paired, is_read1, is_read2
                        )
                        distributed += 1
                    else:
                        bucket_idx = rng.randrange(n_buckets)
                        bucket_writers[bucket_idx].write_record(
                            seq, is_paired, is_read1, is_read2
                        )
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
        finally:
            for w in bucket_writers:
                w.close()
            for fh in bucket_fhs:
                fh.close()

        if not quiet:
            print("\n[ZNA] Shuffling and writing output ...", file=sys.stderr)

        # ── Pass 2: read each bucket, shuffle, write to output ────────
        bucket_order = list(range(n_buckets))
        rng.shuffle(bucket_order)

        written = 0
        with open(output_path, "wb") as out_fh:
            out_writer = ZnaWriter(out_fh, out_header, block_size=block_size)
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
                        for rec in breader.records():
                            _seq, is_paired, is_read1, _is_read2 = rec
                            if is_paired and is_read1:
                                if current_unit:
                                    units.append(current_unit)
                                current_unit = [rec]
                            elif is_paired and _is_read2:
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
                        for seq, is_paired, is_read1, is_read2 in unit:
                            out_writer.write_record(
                                seq, is_paired, is_read1, is_read2
                            )
                    written += len(units)

                    if not quiet:
                        print(
                            f"      Written {written:,} / {n_units:,} units ...",
                            end="\r",
                            file=sys.stderr,
                        )
            finally:
                out_writer.close()

    return written, n_records
