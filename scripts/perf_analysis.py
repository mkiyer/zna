#!/usr/bin/env python3
"""
Performance analysis of ZNA with C++ label optimizations.

Targets: merged.fq.gz (all records)
Tests:
  1. Baseline encode (no labels)
  2. Labeled encode (10 labels, C++ path)
  3. Baseline decode (no labels)
  4. Labeled decode (10 labels, C++ path)
  5. Roundtrip correctness (500K records)
  6. cProfile breakdown (500K records)
"""
import gzip
import os
import sys
import tempfile
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from zna.core import (
    ZnaHeader, ZnaWriter, ZnaReader,
    COMPRESSION_ZSTD, _accel_mod,
)
from zna.dtypes import LabelDef, parse_dtype, resolve_missing
from zna.cli import (
    parse_fastq, parse_fastq_with_headers,
    build_tag_extractor, extract_labels_from_header,
)

# Try C++ fast extractor
try:
    from zna._accel import extract_labels_fast as cpp_extract_labels_fast
    HAS_CPP_EXTRACT = True
except ImportError:
    HAS_CPP_EXTRACT = False

from zna import is_accelerated

MERGED_PATH = "/Users/mkiyer/Downloads/rigel_runs/sim_ccle_hela_salmon/gdna_high_ss_0.90_nrna_default/align_minimap2/merged.fq.gz"

COMMON_LABEL_DEFS = (
    LabelDef(0, "NM", "Edit distance", parse_dtype("C")),
    LabelDef(1, "ms", "Minimap2 score", parse_dtype("S")),
    LabelDef(2, "AS", "Alignment score", parse_dtype("s")),
    LabelDef(3, "nn", "Ambiguous bases", parse_dtype("C")),
    LabelDef(4, "tp", "Alignment type", parse_dtype("A"), missing=ord('*')),
    LabelDef(5, "cm", "Chaining minimizers", parse_dtype("S")),
    LabelDef(6, "s1", "Chaining score", parse_dtype("S")),
    LabelDef(7, "s2", "Second best score", parse_dtype("S")),
    LabelDef(8, "de", "Divergence", parse_dtype("f")),
    LabelDef(9, "rl", "Read length", parse_dtype("S")),
)


def make_tag_map(label_defs):
    return build_tag_extractor(label_defs)


def build_cpp_tag_specs(label_defs):
    """Build the C++ tag specs list and missing values tuple."""
    CONV_INT, CONV_FLOAT, CONV_ORD = 0, 1, 2
    tag_specs = []
    for ldef in label_defs:
        tag = ldef.name.encode('ascii')
        if ldef.dtype.code == 'A':
            conv = CONV_ORD
        elif ldef.dtype.code in ('f', 'd'):
            conv = CONV_FLOAT
        else:
            conv = CONV_INT
        tag_specs.append((tag, conv))
    missing_values = tuple(resolve_missing(ld) for ld in label_defs)
    return tag_specs, missing_values


def hr(title):
    print(f"\n{'=' * 70}")
    print(f"  {title}")
    print(f"{'=' * 70}")


# ---------------------------------------------------------------------------

def test_baseline_encode():
    """Encode merged.fq.gz without labels — baseline."""
    hr("TEST 1: Baseline encode (no labels)")
    header = ZnaHeader(
        read_group="perf_baseline",
        seq_len_bytes=2,
        compression_method=COMPRESSION_ZSTD,
    )

    with tempfile.NamedTemporaryFile(suffix=".zna", delete=False) as tmp:
        tmp_path = tmp.name

    try:
        count = 0
        t0 = time.perf_counter()
        with open(tmp_path, "wb") as fh:
            with ZnaWriter(fh, header, npolicy='A') as w:
                with gzip.open(MERGED_PATH, "rb") as fq:
                    for seq in parse_fastq(fq):
                        w.write_record(seq, False, False, False)
                        count += 1
        elapsed = time.perf_counter() - t0
        file_size = os.path.getsize(tmp_path)

        print(f"  Records:  {count:,}")
        print(f"  Time:     {elapsed:.2f}s")
        print(f"  Speed:    {count/elapsed:,.0f} rec/s")
        print(f"  ZNA size: {file_size / (1024*1024):.2f} MB")
        return tmp_path, count, elapsed, file_size
    except Exception as e:
        os.unlink(tmp_path)
        raise


def test_labeled_encode():
    """Encode merged.fq.gz with 10 labels — C++ path."""
    hr("TEST 2: Labeled encode (10 labels, C++ path)")
    label_defs = COMMON_LABEL_DEFS
    num_labels = len(label_defs)

    header = ZnaHeader(
        read_group="perf_labeled",
        seq_len_bytes=2,
        compression_method=COMPRESSION_ZSTD,
        labels=label_defs,
    )

    # Use C++ fast extractor if available
    if HAS_CPP_EXTRACT:
        tag_specs, missing_values = build_cpp_tag_specs(label_defs)
        print(f"  Using: C++ extract_labels_fast")
    else:
        tag_map = make_tag_map(label_defs)
        print(f"  Using: Python extract_labels_from_header")

    with tempfile.NamedTemporaryFile(suffix=".zna", delete=False) as tmp:
        tmp_path = tmp.name

    try:
        count = 0
        t0 = time.perf_counter()
        with open(tmp_path, "wb") as fh:
            with ZnaWriter(fh, header, npolicy='A') as w:
                with gzip.open(MERGED_PATH, "rb") as fq:
                    if HAS_CPP_EXTRACT:
                        for raw_header, seq in parse_fastq_with_headers(fq):
                            labels = cpp_extract_labels_fast(raw_header, tag_specs, num_labels, missing_values)
                            w.write_record(seq, False, False, False, labels=labels)
                            count += 1
                    else:
                        for raw_header, seq in parse_fastq_with_headers(fq):
                            labels = extract_labels_from_header(raw_header, tag_map, num_labels, label_defs=label_defs)
                            w.write_record(seq, False, False, False, labels=labels)
                            count += 1
        elapsed = time.perf_counter() - t0
        file_size = os.path.getsize(tmp_path)

        print(f"  Records:  {count:,}")
        print(f"  Time:     {elapsed:.2f}s")
        print(f"  Speed:    {count/elapsed:,.0f} rec/s")
        print(f"  ZNA size: {file_size / (1024*1024):.2f} MB")
        return tmp_path, count, elapsed, file_size
    except Exception as e:
        os.unlink(tmp_path)
        raise


def test_baseline_decode(zna_path, expected_count):
    """Decode unlabeled ZNA file."""
    hr("TEST 3: Baseline decode (no labels)")

    t0 = time.perf_counter()
    n = 0
    with open(zna_path, "rb") as fh:
        for _ in ZnaReader(fh).records():
            n += 1
    elapsed = time.perf_counter() - t0

    print(f"  Records:  {n:,}")
    print(f"  Time:     {elapsed:.2f}s")
    print(f"  Speed:    {n/elapsed:,.0f} rec/s")
    assert n == expected_count
    return elapsed


def test_labeled_decode(zna_path, expected_count):
    """Decode labeled ZNA file."""
    hr("TEST 4: Labeled decode (10 labels, C++ path)")

    t0 = time.perf_counter()
    n = 0
    with open(zna_path, "rb") as fh:
        for rec in ZnaReader(fh).records():
            n += 1
    elapsed = time.perf_counter() - t0

    print(f"  Records:  {n:,}")
    print(f"  Time:     {elapsed:.2f}s")
    print(f"  Speed:    {n/elapsed:,.0f} rec/s")
    assert n == expected_count
    return elapsed


def test_roundtrip_correctness(max_records=500_000):
    """Verify roundtrip correctness of labeled encode/decode."""
    hr(f"TEST 5: Roundtrip correctness ({max_records:,} records)")
    label_defs = COMMON_LABEL_DEFS
    num_labels = len(label_defs)

    header = ZnaHeader(
        read_group="roundtrip_test",
        seq_len_bytes=2,
        compression_method=COMPRESSION_ZSTD,
        labels=label_defs,
    )

    if HAS_CPP_EXTRACT:
        tag_specs, missing_values = build_cpp_tag_specs(label_defs)
    else:
        tag_map = make_tag_map(label_defs)

    original_seqs = []
    original_labels = []
    with tempfile.NamedTemporaryFile(suffix=".zna", delete=False) as tmp:
        tmp_path = tmp.name

    try:
        with open(tmp_path, "wb") as fh:
            with ZnaWriter(fh, header, npolicy='A') as w:
                with gzip.open(MERGED_PATH, "rb") as fq:
                    for raw_header, seq in parse_fastq_with_headers(fq):
                        if HAS_CPP_EXTRACT:
                            labels = cpp_extract_labels_fast(raw_header, tag_specs, num_labels, missing_values)
                        else:
                            labels = extract_labels_from_header(raw_header, tag_map, num_labels, label_defs=label_defs)
                        w.write_record(seq, False, False, False, labels=labels)
                        original_seqs.append(seq)
                        original_labels.append(labels)
                        if len(original_seqs) >= max_records:
                            break

        # Decode & verify
        decoded_seqs = []
        decoded_labels = []
        with open(tmp_path, "rb") as fh:
            reader = ZnaReader(fh)
            assert reader.header.num_labels == num_labels
            for rec in reader.records():
                seq, _, _, _, labels = rec
                decoded_seqs.append(seq)
                decoded_labels.append(labels)

        assert len(decoded_seqs) == len(original_seqs), \
            f"Count mismatch: {len(decoded_seqs)} vs {len(original_seqs)}"

        seq_mismatches = 0
        label_mismatches = 0
        float_tol = 0
        for i in range(len(original_seqs)):
            if decoded_seqs[i] != original_seqs[i]:
                seq_mismatches += 1
            if decoded_labels[i] != original_labels[i]:
                all_close = True
                for a, b in zip(original_labels[i], decoded_labels[i]):
                    if isinstance(a, float) or isinstance(b, float):
                        if abs(float(a) - float(b)) > 1e-5:
                            all_close = False
                    elif a != b:
                        all_close = False
                if all_close:
                    float_tol += 1
                else:
                    label_mismatches += 1

        if seq_mismatches == 0 and label_mismatches == 0:
            print(f"  PASSED — all {len(original_seqs):,} sequences and labels match")
            if float_tol:
                print(f"  (Note: {float_tol} float32 precision differences within tolerance)")
        else:
            print(f"  *** FAILED — {seq_mismatches} seq, {label_mismatches} label mismatches ***")

    finally:
        os.unlink(tmp_path)


def test_profile_encode(max_records=500_000):
    """cProfile breakdown of labeled encode."""
    hr(f"TEST 6: cProfile of labeled encode ({max_records:,} records)")

    import cProfile
    import pstats
    from io import StringIO

    label_defs = COMMON_LABEL_DEFS
    num_labels = len(label_defs)

    header = ZnaHeader(
        read_group="profile",
        seq_len_bytes=2,
        compression_method=COMPRESSION_ZSTD,
        labels=label_defs,
    )

    if HAS_CPP_EXTRACT:
        tag_specs, missing_values = build_cpp_tag_specs(label_defs)
    else:
        tag_map = make_tag_map(label_defs)

    def do_encode():
        with tempfile.NamedTemporaryFile(suffix=".zna", delete=False) as tmp:
            tmp_path = tmp.name
        try:
            count = 0
            with open(tmp_path, "wb") as fh:
                with ZnaWriter(fh, header, npolicy='A') as w:
                    with gzip.open(MERGED_PATH, "rb") as fq:
                        if HAS_CPP_EXTRACT:
                            for raw_header, seq in parse_fastq_with_headers(fq):
                                labels = cpp_extract_labels_fast(raw_header, tag_specs, num_labels, missing_values)
                                w.write_record(seq, False, False, False, labels=labels)
                                count += 1
                                if count >= max_records:
                                    break
                        else:
                            for raw_header, seq in parse_fastq_with_headers(fq):
                                labels = extract_labels_from_header(raw_header, tag_map, num_labels, label_defs=label_defs)
                                w.write_record(seq, False, False, False, labels=labels)
                                count += 1
                                if count >= max_records:
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
    ps.print_stats(25)
    print(s.getvalue())


# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print(f"Backend: {'C++ accel' if is_accelerated() else 'Python'}")
    print(f"C++ label encode: {'yes' if _accel_mod else 'no'}")
    print(f"C++ label decode: {'yes' if _accel_mod else 'no'}")
    print(f"C++ label extract: {'yes' if HAS_CPP_EXTRACT else 'no'}")
    print(f"Input: {MERGED_PATH}")
    if not os.path.exists(MERGED_PATH):
        sys.exit(f"File not found: {MERGED_PATH}")
    print(f"File size: {os.path.getsize(MERGED_PATH) / (1024*1024):.1f} MB (gzipped)")

    # Test 1: Baseline encode
    baseline_path, baseline_count, t_baseline_enc, baseline_size = test_baseline_encode()

    # Test 2: Labeled encode
    labeled_path, labeled_count, t_labeled_enc, labeled_size = test_labeled_encode()

    # Test 3: Baseline decode
    t_baseline_dec = test_baseline_decode(baseline_path, baseline_count)

    # Test 4: Labeled decode
    t_labeled_dec = test_labeled_decode(labeled_path, labeled_count)

    # Cleanup temp files
    os.unlink(baseline_path)
    os.unlink(labeled_path)

    # Test 5: Roundtrip correctness
    test_roundtrip_correctness()

    # Test 6: cProfile
    test_profile_encode()

    # Summary
    hr("SUMMARY")
    print(f"  {'Metric':<30} {'Baseline':>15} {'Labeled':>15} {'Overhead':>12}")
    print(f"  {'-'*30} {'-'*15} {'-'*15} {'-'*12}")
    print(f"  {'Records':<30} {baseline_count:>15,} {labeled_count:>15,}")
    print(f"  {'Encode time':<30} {t_baseline_enc:>14.2f}s {t_labeled_enc:>14.2f}s {(t_labeled_enc/t_baseline_enc - 1)*100:>+10.1f}%")
    print(f"  {'Encode speed (rec/s)':<30} {baseline_count/t_baseline_enc:>15,.0f} {labeled_count/t_labeled_enc:>15,.0f}")
    print(f"  {'Decode time':<30} {t_baseline_dec:>14.2f}s {t_labeled_dec:>14.2f}s {(t_labeled_dec/t_baseline_dec - 1)*100:>+10.1f}%")
    print(f"  {'Decode speed (rec/s)':<30} {baseline_count/t_baseline_dec:>15,.0f} {labeled_count/t_labeled_dec:>15,.0f}")
    print(f"  {'ZNA file size (MB)':<30} {baseline_size/(1024*1024):>15.2f} {labeled_size/(1024*1024):>15.2f} {(labeled_size/baseline_size - 1)*100:>+10.1f}%")
    print()
    print(f"  Previous baseline (from profiling): encode 80,220 rec/s (labeled), decode 610,030 rec/s (labeled)")
    print()
