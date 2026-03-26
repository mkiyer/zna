#!/usr/bin/env python3
"""
Profile and validate ZNA labeled encoding against real FASTQ files.

Tests:
  1. Performance profiling (encode time, decode time, file sizes)
  2. Label parsing accuracy (roundtrip)
  3. Sequence roundtrip accuracy
  4. ZNA label file creation/parsing
"""
import gzip
import os
import sys
import tempfile
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from zna.core import (
    ZnaHeader, ZnaWriter, ZnaReader,
    COMPRESSION_ZSTD, COMPRESSION_NONE,
)
from zna.dtypes import LabelDef, parse_dtype
from zna.cli import (
    parse_fastq_with_headers,
    build_tag_extractor,
    extract_labels_from_header,
)

BASE = "/Users/mkiyer/Downloads/rigel_runs/sim_ccle_hela_salmon/gdna_high_ss_0.90_nrna_default/align_minimap2"
R1_PATH = f"{BASE}/r1.fq.gz"
R2_PATH = f"{BASE}/r2.fq.gz"
MERGED_PATH = f"{BASE}/merged.fq.gz"


# ── Label definitions ─────────────────────────────────────────────
# Tags present in both R1 and R2 (safe for paired encoding):
# Using only tags guaranteed present in all mapped reads.
# Unmapped reads (only rl tag) get defaults.
COMMON_LABEL_DEFS = (
    LabelDef(0, "NM", "Edit distance", parse_dtype("C")),                       # uint8
    LabelDef(1, "ms", "Minimap2 score", parse_dtype("S")),                       # uint16
    LabelDef(2, "AS", "Alignment score", parse_dtype("s")),                      # int16
    LabelDef(3, "nn", "Ambiguous bases", parse_dtype("C")),                      # uint8
    LabelDef(4, "tp", "Alignment type", parse_dtype("A"), missing=ord('*')),     # char -> uint8, unmapped marker
    LabelDef(5, "cm", "Chaining minimizers", parse_dtype("S")),                  # uint16
    LabelDef(6, "s1", "Chaining score", parse_dtype("S")),                       # uint16
    LabelDef(7, "s2", "Second best score", parse_dtype("S")),                    # uint16
    LabelDef(8, "de", "Divergence", parse_dtype("f")),                           # float32
    LabelDef(9, "rl", "Read length", parse_dtype("S")),                          # uint16
)


def make_tag_map(label_defs):
    return build_tag_extractor(label_defs)


def timed(func, *args, **kwargs):
    t0 = time.perf_counter()
    result = func(*args, **kwargs)
    elapsed = time.perf_counter() - t0
    return elapsed, result


# ── Test 1: Header parsing from all three files ──────────────────

def test_header_parsing():
    """Verify that SAM tags can be extracted from all three FASTQ files."""
    print("=" * 70)
    print("TEST 1: Header parsing correctness")
    print("=" * 70)

    tag_map = make_tag_map(COMMON_LABEL_DEFS)
    num_labels = len(COMMON_LABEL_DEFS)

    for name, path in [("R1", R1_PATH), ("R2", R2_PATH), ("Merged", MERGED_PATH)]:
        print(f"\n--- {name}: {path} ---")
        errors = 0
        parsed = 0
        merged_count = 0
        with gzip.open(path, "rb") as fh:
            for raw_header, seq in parse_fastq_with_headers(fh):
                try:
                    labels = extract_labels_from_header(raw_header, tag_map, num_labels, label_defs=COMMON_LABEL_DEFS)
                    parsed += 1
                    if b'merged_' in raw_header:
                        merged_count += 1
                except Exception as e:
                    errors += 1
                    if errors <= 3:
                        print(f"  ERROR: {e}")
                        print(f"  Header: {raw_header[:120]}")
                if parsed >= 100000:  # sample first 100K for speed
                    break

        print(f"  Parsed: {parsed}, Errors: {errors}, Merged: {merged_count}")
        if errors > 0:
            print(f"  *** FAILED — {errors} parsing errors ***")
        else:
            print(f"  PASSED")


# ── Test 2: Single-file encode + decode roundtrip (R1 only) ──────

def test_roundtrip_labeled_single(max_records=500_000):
    """Encode R1 with labels, decode, verify sequences and labels match."""
    print("\n" + "=" * 70)
    print(f"TEST 2: Labeled roundtrip (R1, {max_records:,} records)")
    print("=" * 70)

    label_defs = COMMON_LABEL_DEFS
    tag_map = make_tag_map(label_defs)
    num_labels = len(label_defs)

    header = ZnaHeader(
        read_group="roundtrip_test",
        description="R1 labeled roundtrip",
        seq_len_bytes=2,
        compression_method=COMPRESSION_ZSTD,
        labels=label_defs,
    )

    # Read & encode
    original_seqs = []
    original_labels = []
    with tempfile.NamedTemporaryFile(suffix=".zna", delete=False) as tmp:
        tmp_path = tmp.name

    try:
        t_encode_start = time.perf_counter()
        with open(tmp_path, "wb") as fh:
            with ZnaWriter(fh, header, npolicy='A') as w:
                with gzip.open(R1_PATH, "rb") as fq:
                    for raw_header, seq in parse_fastq_with_headers(fq):
                        labels = extract_labels_from_header(raw_header, tag_map, num_labels, label_defs=COMMON_LABEL_DEFS)
                        w.write_record(seq, False, False, False, labels=labels)
                        original_seqs.append(seq)
                        original_labels.append(labels)
                        if len(original_seqs) >= max_records:
                            break
        t_encode = time.perf_counter() - t_encode_start

        file_size = os.path.getsize(tmp_path)
        print(f"  Encoded {len(original_seqs):,} records in {t_encode:.2f}s")
        print(f"  ZNA file size: {file_size / (1024*1024):.2f} MB")
        print(f"  Encode speed: {len(original_seqs)/t_encode:,.0f} rec/s")

        # Decode & verify
        t_decode_start = time.perf_counter()
        decoded_seqs = []
        decoded_labels = []
        with open(tmp_path, "rb") as fh:
            reader = ZnaReader(fh)
            assert reader.header.num_labels == num_labels
            for rec in reader.records():
                seq, is_paired, is_r1, is_r2, labels = rec
                decoded_seqs.append(seq)
                decoded_labels.append(labels)
        t_decode = time.perf_counter() - t_decode_start

        print(f"  Decoded {len(decoded_seqs):,} records in {t_decode:.2f}s")
        print(f"  Decode speed: {len(decoded_seqs)/t_decode:,.0f} rec/s")

        # Verify
        assert len(decoded_seqs) == len(original_seqs), \
            f"Count mismatch: {len(decoded_seqs)} vs {len(original_seqs)}"

        seq_mismatches = 0
        label_mismatches = 0
        float_tolerance_mismatches = 0
        for i in range(len(original_seqs)):
            if decoded_seqs[i] != original_seqs[i]:
                seq_mismatches += 1
                if seq_mismatches <= 3:
                    print(f"  SEQ MISMATCH at {i}: {original_seqs[i][:50]}... vs {decoded_seqs[i][:50]}...")
            if decoded_labels[i] != original_labels[i]:
                # Check if it's just float32 precision loss
                all_close = True
                for j, (a, b) in enumerate(zip(original_labels[i], decoded_labels[i])):
                    if isinstance(a, float) or isinstance(b, float):
                        if abs(float(a) - float(b)) > 1e-5:
                            all_close = False
                    elif a != b:
                        all_close = False
                if all_close:
                    float_tolerance_mismatches += 1
                else:
                    label_mismatches += 1
                    if label_mismatches <= 3:
                        print(f"  LABEL MISMATCH at {i}: {original_labels[i]} vs {decoded_labels[i]}")

        if seq_mismatches == 0 and label_mismatches == 0:
            print(f"  PASSED — all sequences and labels match")
            if float_tolerance_mismatches > 0:
                print(f"  (Note: {float_tolerance_mismatches} float32 precision differences within tolerance)")
        else:
            print(f"  *** FAILED — {seq_mismatches} seq mismatches, {label_mismatches} label mismatches ***")

    finally:
        os.unlink(tmp_path)


# ── Test 3: Merged file encode + decode roundtrip ─────────────────

def test_roundtrip_merged(max_records=500_000):
    """Encode merged FASTQ with labels, verify roundtrip."""
    print("\n" + "=" * 70)
    print(f"TEST 3: Merged FASTQ labeled roundtrip ({max_records:,} records)")
    print("=" * 70)

    label_defs = COMMON_LABEL_DEFS
    tag_map = make_tag_map(label_defs)
    num_labels = len(label_defs)

    header = ZnaHeader(
        read_group="merged_test",
        seq_len_bytes=2,
        compression_method=COMPRESSION_ZSTD,
        labels=label_defs,
    )

    original_seqs = []
    original_labels = []
    merged_count = 0
    with tempfile.NamedTemporaryFile(suffix=".zna", delete=False) as tmp:
        tmp_path = tmp.name

    try:
        t_encode_start = time.perf_counter()
        with open(tmp_path, "wb") as fh:
            with ZnaWriter(fh, header, npolicy='A') as w:
                with gzip.open(MERGED_PATH, "rb") as fq:
                    for raw_header, seq in parse_fastq_with_headers(fq):
                        labels = extract_labels_from_header(raw_header, tag_map, num_labels, label_defs=COMMON_LABEL_DEFS)
                        w.write_record(seq, False, False, False, labels=labels)
                        original_seqs.append(seq)
                        original_labels.append(labels)
                        if b'merged_' in raw_header:
                            merged_count += 1
                        if len(original_seqs) >= max_records:
                            break
        t_encode = time.perf_counter() - t_encode_start

        file_size = os.path.getsize(tmp_path)
        print(f"  Encoded {len(original_seqs):,} records ({merged_count} merged) in {t_encode:.2f}s")
        print(f"  ZNA file size: {file_size / (1024*1024):.2f} MB")
        print(f"  Encode speed: {len(original_seqs)/t_encode:,.0f} rec/s")

        # Decode & verify
        t_decode_start = time.perf_counter()
        decoded_seqs = []
        decoded_labels = []
        with open(tmp_path, "rb") as fh:
            reader = ZnaReader(fh)
            for rec in reader.records():
                seq, ip, r1, r2, labels = rec
                decoded_seqs.append(seq)
                decoded_labels.append(labels)
        t_decode = time.perf_counter() - t_decode_start

        print(f"  Decoded {len(decoded_seqs):,} records in {t_decode:.2f}s")
        print(f"  Decode speed: {len(decoded_seqs)/t_decode:,.0f} rec/s")

        seq_mismatches = 0
        label_mismatches = 0
        float_tolerance_mismatches = 0
        for i in range(len(original_seqs)):
            if decoded_seqs[i] != original_seqs[i]:
                seq_mismatches += 1
            if decoded_labels[i] != original_labels[i]:
                all_close = True
                for j, (a, b) in enumerate(zip(original_labels[i], decoded_labels[i])):
                    if isinstance(a, float) or isinstance(b, float):
                        if abs(float(a) - float(b)) > 1e-5:
                            all_close = False
                    elif a != b:
                        all_close = False
                if all_close:
                    float_tolerance_mismatches += 1
                else:
                    label_mismatches += 1

        if seq_mismatches == 0 and label_mismatches == 0:
            print(f"  PASSED — all sequences and labels match")
            if float_tolerance_mismatches > 0:
                print(f"  (Note: {float_tolerance_mismatches} float32 precision differences within tolerance)")
        else:
            print(f"  *** FAILED — {seq_mismatches} seq, {label_mismatches} label mismatches ***")

    finally:
        os.unlink(tmp_path)


# ── Test 4: Full encode performance profile (R1 all records, no labels) ──

def test_encode_performance_baseline():
    """Encode R1 without labels to establish baseline performance."""
    print("\n" + "=" * 70)
    print("TEST 4: Encode performance — BASELINE (no labels)")
    print("=" * 70)

    from zna.cli import parse_fastq

    header = ZnaHeader(
        read_group="perf_baseline",
        compression_method=COMPRESSION_ZSTD,
    )

    with tempfile.NamedTemporaryFile(suffix=".zna", delete=False) as tmp:
        tmp_path = tmp.name

    try:
        count = 0
        t0 = time.perf_counter()
        with open(tmp_path, "wb") as fh:
            with ZnaWriter(fh, header, npolicy='A') as w:
                with gzip.open(R1_PATH, "rb") as fq:
                    for seq in parse_fastq(fq):
                        w.write_record(seq, False, False, False)
                        count += 1
        t_total = time.perf_counter() - t0
        file_size = os.path.getsize(tmp_path)

        print(f"  Records: {count:,}")
        print(f"  Time: {t_total:.2f}s")
        print(f"  Speed: {count/t_total:,.0f} rec/s")
        print(f"  ZNA size: {file_size / (1024*1024):.2f} MB")
    finally:
        os.unlink(tmp_path)


def test_encode_performance_labeled():
    """Encode R1 with labels to measure label overhead."""
    print("\n" + "=" * 70)
    print("TEST 5: Encode performance — WITH LABELS (10 labels)")
    print("=" * 70)

    label_defs = COMMON_LABEL_DEFS
    tag_map = make_tag_map(label_defs)
    num_labels = len(label_defs)

    header = ZnaHeader(
        read_group="perf_labeled",
        compression_method=COMPRESSION_ZSTD,
        labels=label_defs,
    )

    with tempfile.NamedTemporaryFile(suffix=".zna", delete=False) as tmp:
        tmp_path = tmp.name

    try:
        count = 0
        t0 = time.perf_counter()
        with open(tmp_path, "wb") as fh:
            with ZnaWriter(fh, header, npolicy='A') as w:
                with gzip.open(R1_PATH, "rb") as fq:
                    for raw_header, seq in parse_fastq_with_headers(fq):
                        labels = extract_labels_from_header(raw_header, tag_map, num_labels, label_defs=COMMON_LABEL_DEFS)
                        w.write_record(seq, False, False, False, labels=labels)
                        count += 1
        t_total = time.perf_counter() - t0
        file_size = os.path.getsize(tmp_path)

        print(f"  Records: {count:,}")
        print(f"  Time: {t_total:.2f}s")
        print(f"  Speed: {count/t_total:,.0f} rec/s")
        print(f"  ZNA size: {file_size / (1024*1024):.2f} MB")
    finally:
        os.unlink(tmp_path)


# ── Test 6: Decode performance ────────────────────────────────────

def test_decode_performance():
    """Measure decode speed for labeled and unlabeled files."""
    print("\n" + "=" * 70)
    print("TEST 6: Decode performance comparison")
    print("=" * 70)

    # First encode both files
    label_defs = COMMON_LABEL_DEFS
    tag_map = make_tag_map(label_defs)
    num_labels = len(label_defs)

    with tempfile.TemporaryDirectory() as tmpdir:
        # Unlabeled
        unlabeled_path = f"{tmpdir}/unlabeled.zna"
        header_no = ZnaHeader(read_group="perf", compression_method=COMPRESSION_ZSTD)
        count = 0
        from zna.cli import parse_fastq
        with open(unlabeled_path, "wb") as fh:
            with ZnaWriter(fh, header_no, npolicy='A') as w:
                with gzip.open(R1_PATH, "rb") as fq:
                    for seq in parse_fastq(fq):
                        w.write_record(seq, False, False, False)
                        count += 1
        print(f"  Encoded {count:,} records (unlabeled + labeled)")

        # Labeled
        labeled_path = f"{tmpdir}/labeled.zna"
        header_lbl = ZnaHeader(read_group="perf", compression_method=COMPRESSION_ZSTD, labels=label_defs)
        with open(labeled_path, "wb") as fh:
            with ZnaWriter(fh, header_lbl, npolicy='A') as w:
                with gzip.open(R1_PATH, "rb") as fq:
                    for raw_header, seq in parse_fastq_with_headers(fq):
                        labels = extract_labels_from_header(raw_header, tag_map, num_labels, label_defs=COMMON_LABEL_DEFS)
                        w.write_record(seq, False, False, False, labels=labels)

        # Decode unlabeled
        t0 = time.perf_counter()
        n = 0
        with open(unlabeled_path, "rb") as fh:
            for _ in ZnaReader(fh).records():
                n += 1
        t_unlabeled = time.perf_counter() - t0
        print(f"  Decode unlabeled: {n:,} records in {t_unlabeled:.2f}s ({n/t_unlabeled:,.0f} rec/s)")

        # Decode labeled
        t0 = time.perf_counter()
        n = 0
        with open(labeled_path, "rb") as fh:
            for rec in ZnaReader(fh).records():
                n += 1
        t_labeled = time.perf_counter() - t0
        print(f"  Decode labeled:   {n:,} records in {t_labeled:.2f}s ({n/t_labeled:,.0f} rec/s)")

        overhead = (t_labeled - t_unlabeled) / t_unlabeled * 100
        print(f"  Label decode overhead: {overhead:+.1f}%")


# ── Test 7: cProfile breakdown ────────────────────────────────────

def test_profile_encode():
    """Detailed cProfile of labeled encode (first 500K records)."""
    print("\n" + "=" * 70)
    print("TEST 7: cProfile of labeled encode (500K records)")
    print("=" * 70)

    import cProfile
    import pstats
    from io import StringIO

    label_defs = COMMON_LABEL_DEFS
    tag_map = make_tag_map(label_defs)
    num_labels = len(label_defs)

    header = ZnaHeader(
        read_group="profile",
        compression_method=COMPRESSION_ZSTD,
        labels=label_defs,
    )

    def do_encode():
        with tempfile.NamedTemporaryFile(suffix=".zna", delete=False) as tmp:
            tmp_path = tmp.name
        try:
            count = 0
            with open(tmp_path, "wb") as fh:
                with ZnaWriter(fh, header, npolicy='A') as w:
                    with gzip.open(R1_PATH, "rb") as fq:
                        for raw_header, seq in parse_fastq_with_headers(fq):
                            labels = extract_labels_from_header(raw_header, tag_map, num_labels, label_defs=COMMON_LABEL_DEFS)
                            w.write_record(seq, False, False, False, labels=labels)
                            count += 1
                            if count >= 500_000:
                                break
            return count
        finally:
            os.unlink(tmp_path)

    pr = cProfile.Profile()
    pr.enable()
    n = do_encode()
    pr.disable()

    s = StringIO()
    ps = pstats.Stats(pr, stream=s).sort_stats('cumulative')
    ps.print_stats(30)
    print(s.getvalue())


# ── Test 8: ZNA inspect labeled file ──────────────────────────────

def test_zna_inspect():
    """Create a labeled ZNA file and verify the header / inspect output."""
    print("\n" + "=" * 70)
    print("TEST 8: ZNA labeled file creation + header inspection")
    print("=" * 70)

    label_defs = COMMON_LABEL_DEFS
    tag_map = make_tag_map(label_defs)
    num_labels = len(label_defs)

    header = ZnaHeader(
        read_group="inspect_test",
        description="Test file for inspection",
        compression_method=COMPRESSION_ZSTD,
        labels=label_defs,
    )

    with tempfile.NamedTemporaryFile(suffix=".zna", delete=False) as tmp:
        tmp_path = tmp.name

    try:
        with open(tmp_path, "wb") as fh:
            with ZnaWriter(fh, header, npolicy='A') as w:
                with gzip.open(R1_PATH, "rb") as fq:
                    count = 0
                    for raw_header, seq in parse_fastq_with_headers(fq):
                        labels = extract_labels_from_header(raw_header, tag_map, num_labels, label_defs=COMMON_LABEL_DEFS)
                        w.write_record(seq, False, False, False, labels=labels)
                        count += 1
                        if count >= 10000:
                            break

        # Read back and verify header
        with open(tmp_path, "rb") as fh:
            reader = ZnaReader(fh)
            h = reader.header
            print(f"  Read group:   {h.read_group}")
            print(f"  Description:  {h.description}")
            print(f"  Compression:  {'ZSTD' if h.compression_method else 'None'}")
            print(f"  Num labels:   {h.num_labels}")
            for ldef in h.labels:
                print(f"    [{ldef.label_id}] {ldef.name:<4} {ldef.dtype.name:<8} ({ldef.dtype.code}) = {ldef.description!r}")

            # Verify first few records
            recs = []
            for rec in reader.records():
                recs.append(rec)
                if len(recs) >= 5:
                    break

            print(f"\n  First 5 records:")
            for i, rec in enumerate(recs):
                seq, ip, r1, r2, labels = rec
                print(f"    [{i}] seq={seq[:30]}... labels={labels}")

        print(f"\n  PASSED — header and records read successfully")
    finally:
        os.unlink(tmp_path)


# ── Main ──────────────────────────────────────────────────────────

if __name__ == "__main__":
    test_header_parsing()
    test_roundtrip_labeled_single()
    test_roundtrip_merged()
    test_encode_performance_baseline()
    test_encode_performance_labeled()
    test_decode_performance()
    test_profile_encode()
    test_zna_inspect()
