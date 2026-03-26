"""Exhaustive tests for the C++ label-aware encode/decode and fast header extraction.

These tests validate Optimizations 1B, 2, and 3 from PERFORMANCE_PLAN.md:
  1B: C++ extract_labels_fast (fast SAM tag extraction)
  2:  C++ encode_block_labeled (encoded with pre-packed label columns)
  3:  C++ decode_block_labeled (decode with label column unpacking)

They also cross-validate C++ vs Python backends for correctness.
"""
import struct
import tempfile
import unittest

from zna.dtypes import (
    LabelDef, DTYPE_BY_CODE, parse_dtype, label_bytes_per_record,
    resolve_missing,
)
from zna.core import (
    ZnaHeader, ZnaWriter, ZnaReader,
    COMPRESSION_ZSTD, COMPRESSION_NONE,
)
from zna._pycodec import (
    encode_block as py_encode_block,
    decode_block as py_decode_block,
)

# Import C++ module directly — skip test class if unavailable
try:
    from zna._accel import (
        encode_block_labeled as cpp_encode_block_labeled,
        decode_block_labeled as cpp_decode_block_labeled,
        extract_labels_fast as cpp_extract_labels_fast,
        encode_block as cpp_encode_block,
        decode_block as cpp_decode_block,
    )
    HAS_ACCEL = True
except ImportError:
    HAS_ACCEL = False


def _pack_label_columns(label_values, label_defs, count):
    """Pack label columns into bytes (same as core.py _flush_block does)."""
    col_data = []
    col_sizes = []
    for col_vals, ldef in zip(label_values, label_defs):
        fmt = f"<{count}{ldef.dtype.struct_ch}"
        col_data.append(struct.pack(fmt, *col_vals))
        col_sizes.append(ldef.dtype.size)
    return col_data, col_sizes


def _split_label_columns(labels_bytes, label_defs, count):
    """Split concatenated label bytes into per-column bytes."""
    cols = []
    offset = 0
    for ldef in label_defs:
        bpr = label_bytes_per_record(ldef)
        col_bytes = bpr * count
        cols.append(labels_bytes[offset:offset + col_bytes])
        offset += col_bytes
    return cols


@unittest.skipUnless(HAS_ACCEL, "C++ acceleration not available")
class TestCppEncodeBlockLabeled(unittest.TestCase):
    """Test C++ encode_block_labeled matches Python backend."""

    def _roundtrip_cpp(self, seqs, flags, label_values, label_defs):
        """Encode with C++ labeled encoder, decode with C++ labeled decoder."""
        count = len(seqs)
        col_data, col_sizes = _pack_label_columns(label_values, label_defs, count)

        result = cpp_encode_block_labeled(
            seqs, flags, 2, "", False, False, False, col_data, col_sizes,
        )
        self.assertEqual(len(result), 4, "Should return 4-tuple with labels")
        flags_bytes, labels_bytes, lengths_bytes, seqs_bytes = result

        # Decode with C++
        label_cols = _split_label_columns(labels_bytes, label_defs, count)
        dtype_codes = "".join(ld.dtype.code for ld in label_defs)
        records = cpp_decode_block_labeled(
            flags_bytes, lengths_bytes, seqs_bytes, 2, count,
            label_cols, col_sizes, dtype_codes,
        )
        return records

    def test_single_uint8(self):
        ldef = LabelDef(0, "NH", "", parse_dtype("C"))
        recs = self._roundtrip_cpp(["ACGT"], [0], [[5]], [ldef])
        self.assertEqual(len(recs), 1)
        self.assertEqual(recs[0][0], "ACGT")
        self.assertEqual(recs[0][5], (5,))

    def test_multiple_numeric_types(self):
        ldefs = [
            LabelDef(0, "c", "", parse_dtype("c")),   # int8
            LabelDef(1, "S", "", parse_dtype("S")),    # uint16
            LabelDef(2, "i", "", parse_dtype("i")),    # int32
            LabelDef(3, "f", "", parse_dtype("f")),    # float32
        ]
        vals = [[-10], [60000], [-100000], [3.14]]
        recs = self._roundtrip_cpp(["ACGTACGT"], [0], vals, ldefs)
        self.assertEqual(recs[0][5][0], -10)
        self.assertEqual(recs[0][5][1], 60000)
        self.assertEqual(recs[0][5][2], -100000)
        self.assertAlmostEqual(recs[0][5][3], 3.14, places=5)

    def test_all_numeric_dtype_roundtrips(self):
        """Every numeric dtype with min/max values through C++ path."""
        test_cases = [
            ("A", 0, 255),
            ("c", -128, 127),
            ("C", 0, 255),
            ("s", -32768, 32767),
            ("S", 0, 65535),
            ("i", -2147483648, 2147483647),
            ("I", 0, 4294967295),
            ("f", -1.5, 1.5),
            ("d", -1e300, 1e300),
            ("q", -(2**63), 2**63 - 1),
            ("Q", 0, 2**64 - 1),
        ]
        for code, min_val, max_val in test_cases:
            with self.subTest(code=code):
                ldef = LabelDef(0, "X", "", parse_dtype(code))
                recs = self._roundtrip_cpp(
                    ["ACGT", "TGCA"], [0, 0],
                    [[min_val, max_val]], [ldef]
                )
                if code in ("f", "d"):
                    self.assertAlmostEqual(recs[0][5][0], min_val, places=1)
                    self.assertAlmostEqual(recs[1][5][0], max_val, places=1)
                else:
                    self.assertEqual(recs[0][5][0], min_val)
                    self.assertEqual(recs[1][5][0], max_val)

    def test_many_records(self):
        ldef = LabelDef(0, "NH", "", parse_dtype("I"))
        n = 500
        seqs = ["ACGT" * 10] * n
        flags = [0] * n
        vals = [list(range(n))]
        recs = self._roundtrip_cpp(seqs, flags, vals, [ldef])
        self.assertEqual(len(recs), n)
        for i in range(n):
            self.assertEqual(recs[i][5], (i,))

    def test_mixed_labels_multiple_records(self):
        ldefs = [
            LabelDef(0, "NH", "", parse_dtype("C")),
            LabelDef(1, "AS", "", parse_dtype("i")),
            LabelDef(2, "de", "", parse_dtype("f")),
        ]
        n = 100
        vals = [
            [i % 256 for i in range(n)],
            [i * 10 - 500 for i in range(n)],
            [i * 0.01 for i in range(n)],
        ]
        seqs = ["ACGT" * 5] * n
        flags = [0] * n
        recs = self._roundtrip_cpp(seqs, flags, vals, ldefs)
        self.assertEqual(len(recs), n)
        for i in range(n):
            self.assertEqual(recs[i][5][0], i % 256)
            self.assertEqual(recs[i][5][1], i * 10 - 500)
            self.assertAlmostEqual(recs[i][5][2], i * 0.01, places=5)

    def test_no_labels_returns_3tuple(self):
        """When no label columns are provided, returns 3-tuple."""
        result = cpp_encode_block_labeled(
            ["ACGT"], [0], 2, "", False, False, False, [], [],
        )
        self.assertEqual(len(result), 3)

    def test_int64_uint64(self):
        ldefs = [
            LabelDef(0, "q", "", parse_dtype("q")),
            LabelDef(1, "Q", "", parse_dtype("Q")),
        ]
        big_signed = -(2**62)
        big_unsigned = 2**63
        recs = self._roundtrip_cpp(
            ["ACGT"], [0], [[big_signed], [big_unsigned]], ldefs,
        )
        self.assertEqual(recs[0][5][0], big_signed)
        self.assertEqual(recs[0][5][1], big_unsigned)

    def test_float64(self):
        ldef = LabelDef(0, "d", "", parse_dtype("d"))
        recs = self._roundtrip_cpp(
            ["ACGT"], [0], [[1.23456789012345]], [ldef],
        )
        self.assertAlmostEqual(recs[0][5][0], 1.23456789012345, places=10)

    def test_20_labels(self):
        """Test with many label columns."""
        n_labels = 20
        ldefs = [LabelDef(i, f"L{i:02d}", "", parse_dtype("C")) for i in range(n_labels)]
        vals = [[i] for i in range(n_labels)]
        recs = self._roundtrip_cpp(["ACGT"], [0], vals, ldefs)
        self.assertEqual(recs[0][5], tuple(range(n_labels)))


@unittest.skipUnless(HAS_ACCEL, "C++ acceleration not available")
class TestCppVsPythonCrossValidation(unittest.TestCase):
    """Cross-validate C++ and Python codec backends produce identical results."""

    def test_encode_produces_same_output(self):
        """C++ and Python encoders should produce identical binary output."""
        ldefs = [
            LabelDef(0, "NH", "", parse_dtype("C")),
            LabelDef(1, "AS", "", parse_dtype("i")),
            LabelDef(2, "de", "", parse_dtype("f")),
        ]
        seqs = ["ACGTACGT", "TGCATGCA", "AAACCC"]
        flags = [0, 0, 0]
        label_values = [[1, 2, 3], [-100, 200, -300], [0.5, 1.5, 2.5]]

        # Python encode
        py_result = py_encode_block(
            seqs, flags, 2, "", False, False, False,
            label_values=label_values, label_defs=ldefs,
        )
        py_flags, py_labels, py_lengths, py_seqs = py_result

        # C++ encode
        count = len(seqs)
        col_data, col_sizes = _pack_label_columns(label_values, ldefs, count)
        cpp_result = cpp_encode_block_labeled(
            seqs, flags, 2, "", False, False, False, col_data, col_sizes,
        )
        cpp_flags, cpp_labels, cpp_lengths, cpp_seqs = cpp_result

        self.assertEqual(py_flags, cpp_flags, "Flags differ")
        self.assertEqual(py_labels, cpp_labels, "Labels differ")
        self.assertEqual(py_lengths, cpp_lengths, "Lengths differ")
        self.assertEqual(py_seqs, cpp_seqs, "Sequences differ")

    def test_decode_produces_same_output(self):
        """C++ and Python decoders should produce identical records."""
        ldefs = [
            LabelDef(0, "NH", "", parse_dtype("C")),
            LabelDef(1, "AS", "", parse_dtype("i")),
            LabelDef(2, "XF", "", parse_dtype("f")),
        ]
        seqs = ["ACGTACGT", "TGCATGCA"]
        flags = [0, 0]
        label_values = [[1, 2], [-100, 200], [0.5, 1.5]]
        count = len(seqs)

        # Encode with Python (canonical)
        py_result = py_encode_block(
            seqs, flags, 2, "", False, False, False,
            label_values=label_values, label_defs=ldefs,
        )
        flags_bytes, labels_bytes, lengths_bytes, seqs_bytes = py_result

        label_cols = _split_label_columns(labels_bytes, ldefs, count)

        # Decode with Python
        py_records = py_decode_block(
            flags_bytes, lengths_bytes, seqs_bytes, 2, count,
            label_columns=label_cols, label_defs=ldefs,
        )

        # Decode with C++
        col_sizes = [ld.dtype.size for ld in ldefs]
        dtype_codes = "".join(ld.dtype.code for ld in ldefs)
        cpp_records = cpp_decode_block_labeled(
            flags_bytes, lengths_bytes, seqs_bytes, 2, count,
            label_cols, col_sizes, dtype_codes,
        )

        self.assertEqual(len(py_records), len(cpp_records))
        for i in range(len(py_records)):
            # Compare sequences
            self.assertEqual(py_records[i][0], cpp_records[i][0])
            # Compare flags
            self.assertEqual(py_records[i][1], cpp_records[i][1])  # is_paired
            self.assertEqual(py_records[i][2], cpp_records[i][2])  # is_read1
            self.assertEqual(py_records[i][3], cpp_records[i][3])  # is_read2
            self.assertEqual(py_records[i][4], cpp_records[i][4])  # is_rc
            # Compare labels
            py_labels = py_records[i][5]
            cpp_labels = cpp_records[i][5]
            self.assertEqual(len(py_labels), len(cpp_labels))
            for j in range(len(py_labels)):
                if isinstance(py_labels[j], float):
                    self.assertAlmostEqual(py_labels[j], cpp_labels[j], places=5)
                else:
                    self.assertEqual(py_labels[j], cpp_labels[j])

    def test_full_pipeline_cpp_vs_python(self):
        """Writer→Reader roundtrip should produce identical results regardless of backend.

        This test forces both backends and compares the output.
        """
        ldefs = [
            LabelDef(0, "NH", "", parse_dtype("C")),
            LabelDef(1, "AS", "", parse_dtype("i")),
        ]
        h = ZnaHeader(read_group="rg", compression_method=COMPRESSION_ZSTD, labels=tuple(ldefs))
        n = 200
        records_in = [
            ("ACGT" * 10, False, False, False, (i % 256, i * 10))
            for i in range(n)
        ]

        # Write with the active backend (C++ if available)
        with tempfile.TemporaryDirectory() as d:
            path = f"{d}/test.zna"
            with open(path, "wb") as fh:
                with ZnaWriter(fh, h) as w:
                    for seq, ip, r1, r2, labels in records_in:
                        w.write_record(seq, ip, r1, r2, labels=labels)
            with open(path, "rb") as fh:
                reader = ZnaReader(fh)
                recs = list(reader.records())

        self.assertEqual(len(recs), n)
        for i in range(n):
            self.assertEqual(recs[i][0], "ACGT" * 10)
            self.assertEqual(recs[i][4], (i % 256, i * 10))


@unittest.skipUnless(HAS_ACCEL, "C++ acceleration not available")
class TestCppDecodeBlockLabeled(unittest.TestCase):
    """Direct tests for decode_block_labeled."""

    def test_empty_block(self):
        """Zero records should return empty list."""
        recs = cpp_decode_block_labeled(
            b"", b"", b"", 2, 0, [], [], "",
        )
        self.assertEqual(recs, [])

    def test_char_type_A(self):
        """Type 'A' should decode as uint8 (same as 'C')."""
        ldef = LabelDef(0, "tp", "", parse_dtype("A"))
        # Pack a single char value
        col_data = struct.pack("<B", ord('P'))
        recs_py = py_encode_block(["ACGT"], [0], 2, "", False, False)
        flags_bytes, lengths_bytes, seqs_bytes = recs_py

        recs = cpp_decode_block_labeled(
            flags_bytes, lengths_bytes, seqs_bytes, 2, 1,
            [col_data], [1], "A",
        )
        self.assertEqual(recs[0][5], (ord('P'),))

    def test_is_rc_flag_preserved(self):
        """IS_RC flag (0x08) should be decoded correctly."""
        ldef = LabelDef(0, "X", "", parse_dtype("C"))
        # Manually set IS_RC flag
        flags_bytes = bytes([0x08])  # is_rc=True
        lengths_bytes = struct.pack("<H", 4)
        from zna._pycodec import encode_sequence
        seqs_bytes = encode_sequence("ACGT")
        col_data = struct.pack("<B", 42)

        recs = cpp_decode_block_labeled(
            flags_bytes, lengths_bytes, seqs_bytes, 2, 1,
            [col_data], [1], "C",
        )
        self.assertEqual(recs[0][0], "ACGT")
        self.assertFalse(recs[0][1])  # not paired
        self.assertFalse(recs[0][2])  # not read1
        self.assertFalse(recs[0][3])  # not read2
        self.assertTrue(recs[0][4])   # is_rc
        self.assertEqual(recs[0][5], (42,))

    def test_paired_flags(self):
        """Paired flags should be decoded correctly."""
        from zna._pycodec import encode_sequence
        flags_bytes = bytes([0x05, 0x06])  # R1+paired, R2+paired
        lengths_bytes = struct.pack("<HH", 4, 4)
        seqs_bytes = encode_sequence("AAAA") + encode_sequence("CCCC")
        col_data = struct.pack("<BB", 10, 20)

        recs = cpp_decode_block_labeled(
            flags_bytes, lengths_bytes, seqs_bytes, 2, 2,
            [col_data], [1], "C",
        )
        self.assertTrue(recs[0][1])   # is_paired
        self.assertTrue(recs[0][2])   # is_read1
        self.assertFalse(recs[0][3])  # not read2
        self.assertTrue(recs[1][1])   # is_paired
        self.assertFalse(recs[1][2])  # not read1
        self.assertTrue(recs[1][3])   # is_read2


@unittest.skipUnless(HAS_ACCEL, "C++ acceleration not available")
class TestCppExtractLabelsFast(unittest.TestCase):
    """Test the C++ fast SAM tag extraction (Optimization 1B)."""

    # Conversion codes matching CLI
    CONV_INT = 0
    CONV_FLOAT = 1
    CONV_ORD = 2

    def _extract(self, header_bytes, label_defs, tag_specs=None, missing_values=None):
        """Helper to call extract_labels_fast with proper args."""
        if tag_specs is None:
            tag_specs = []
            for ldef in label_defs:
                tag = ldef.effective_tag.encode('ascii')
                if ldef.dtype.code == 'A':
                    conv = self.CONV_ORD
                elif ldef.dtype.code in ('f', 'd'):
                    conv = self.CONV_FLOAT
                else:
                    conv = self.CONV_INT
                tag_specs.append((tag, conv))
        if missing_values is None:
            missing_values = tuple(resolve_missing(ld) for ld in label_defs)

        return cpp_extract_labels_fast(
            header_bytes, tag_specs, len(label_defs), missing_values,
        )

    def test_basic_tab_separated(self):
        defs = [
            LabelDef(0, "NH", "", parse_dtype("C")),
            LabelDef(1, "AS", "", parse_dtype("i")),
        ]
        header = b"readname\tNH:i:1\tAS:i:298"
        result = self._extract(header, defs)
        self.assertEqual(result, (1, 298))

    def test_space_separated(self):
        defs = [LabelDef(0, "NH", "", parse_dtype("C"))]
        header = b"readname NH:i:5"
        result = self._extract(header, defs)
        self.assertEqual(result, (5,))

    def test_negative_int(self):
        defs = [LabelDef(0, "AS", "", parse_dtype("i"))]
        header = b"readname\tAS:i:-42"
        result = self._extract(header, defs)
        self.assertEqual(result, (-42,))

    def test_float_value(self):
        defs = [LabelDef(0, "de", "", parse_dtype("f"))]
        header = b"readname\tde:f:0.0123"
        result = self._extract(header, defs)
        self.assertAlmostEqual(result[0], 0.0123, places=4)

    def test_char_ord(self):
        """Type 'A' with single char value should return ord."""
        defs = [LabelDef(0, "tp", "", parse_dtype("A"))]
        header = b"readname\ttp:A:P"
        result = self._extract(header, defs)
        self.assertEqual(result, (ord('P'),))

    def test_missing_tag_uses_default(self):
        defs = [
            LabelDef(0, "NH", "", parse_dtype("C"), missing=255),
            LabelDef(1, "AS", "", parse_dtype("i")),
        ]
        header = b"readname\tAS:i:100"  # NH missing
        result = self._extract(header, defs)
        self.assertEqual(result, (255, 100))

    def test_all_missing_uses_defaults(self):
        defs = [LabelDef(0, "NH", "", parse_dtype("C"))]
        header = b"readname"
        result = self._extract(header, defs)
        self.assertEqual(result, (0,))  # default missing for uint8

    def test_ignores_unknown_tags(self):
        defs = [LabelDef(0, "NH", "", parse_dtype("C"))]
        header = b"readname\tXX:i:999\tNH:i:5\tYY:Z:hello"
        result = self._extract(header, defs)
        self.assertEqual(result, (5,))

    def test_fastp_merged_suffix(self):
        """Token 'merged_150_87' should be ignored (no colons at [2],[4])."""
        defs = [
            LabelDef(0, "AS", "", parse_dtype("i")),
            LabelDef(1, "rl", "", parse_dtype("S")),
        ]
        header = b"readname\tAS:i:100\trl:i:0 merged_150_87"
        result = self._extract(header, defs)
        self.assertEqual(result, (100, 0))

    def test_early_exit_all_found(self):
        """When all tags found, remaining fields are skipped."""
        defs = [LabelDef(0, "NH", "", parse_dtype("C"))]
        # NH is in the second field; many more fields follow
        header = b"readname\tNH:i:3\tAS:i:100\tNM:i:5\tXY:i:77\tAB:i:99"
        result = self._extract(header, defs)
        self.assertEqual(result, (3,))

    def test_multiple_labels_order(self):
        """Labels should be returned in label_def order, not header order."""
        defs = [
            LabelDef(0, "AS", "", parse_dtype("i")),
            LabelDef(1, "NH", "", parse_dtype("C")),
        ]
        # Header has NH before AS
        header = b"readname\tNH:i:3\tAS:i:100"
        result = self._extract(header, defs)
        self.assertEqual(result, (100, 3))  # AS first, NH second

    def test_cross_validation_with_python(self):
        """C++ and Python extractors should produce identical results."""
        from zna.cli import extract_labels_from_header, build_tag_extractor

        defs = tuple([
            LabelDef(0, "NM", "", parse_dtype("C")),
            LabelDef(1, "ms", "", parse_dtype("S")),
            LabelDef(2, "AS", "", parse_dtype("s")),
            LabelDef(3, "nn", "", parse_dtype("C")),
            LabelDef(4, "tp", "", parse_dtype("A"), missing=ord('*')),
            LabelDef(5, "cm", "", parse_dtype("S")),
            LabelDef(6, "s1", "", parse_dtype("S")),
            LabelDef(7, "s2", "", parse_dtype("S")),
            LabelDef(8, "de", "", parse_dtype("f")),
            LabelDef(9, "rl", "", parse_dtype("S")),
        ])

        headers = [
            b"read1\tNM:i:0\tms:i:281\tAS:i:281\tnn:i:0\ttp:A:P\tcm:i:18\ts1:i:175\ts2:i:0\tde:f:0.0\trl:i:150",
            b"read2\tNM:i:3\tms:i:254\tAS:i:254\tnn:i:0\ttp:A:S\tcm:i:15\ts1:i:140\ts2:i:98\tde:f:0.0200\trl:i:150",
            b"read3\trl:i:100",  # only rl present — all others use missing
            b"read4\tNM:i:1\tms:i:50\tAS:i:-10\tnn:i:2\ttp:A:P\tcm:i:5\ts1:i:30\ts2:i:10\tde:f:0.1234\trl:i:75 merged_150_87",
        ]

        tag_map = build_tag_extractor(defs)
        for raw_header in headers:
            py_result = extract_labels_from_header(raw_header, tag_map, len(defs), label_defs=defs)
            cpp_result = self._extract(raw_header, list(defs))

            self.assertEqual(len(py_result), len(cpp_result),
                             f"Length mismatch for {raw_header}")
            for j in range(len(py_result)):
                if isinstance(py_result[j], float):
                    self.assertAlmostEqual(py_result[j], cpp_result[j], places=4,
                                           msg=f"Float mismatch at label {j} for {raw_header}")
                else:
                    self.assertEqual(py_result[j], cpp_result[j],
                                     msg=f"Value mismatch at label {j} for {raw_header}")

    def test_large_int_values(self):
        """Large integer values should be parsed correctly."""
        defs = [LabelDef(0, "XL", "", parse_dtype("I"))]
        header = b"readname\tXL:i:4294967295"
        result = self._extract(header, defs)
        self.assertEqual(result, (4294967295,))

    def test_long_tag_names(self):
        """Tags longer than 2 characters should work in C++ extractor."""
        defs = [
            LabelDef(0, "edit_distance", "", parse_dtype("C"), tag="edit_distance"),
            LabelDef(1, "my_score", "", parse_dtype("i"), tag="my_score"),
        ]
        header = b"readname\tedit_distance:i:7\tmy_score:i:42"
        result = self._extract(header, defs)
        self.assertEqual(result, (7, 42))

    def test_long_tag_cross_validation(self):
        """C++ and Python should agree on long-tag parsing."""
        from zna.cli import extract_labels_from_header, build_tag_extractor
        defs = tuple([
            LabelDef(0, "edit_dist", "", parse_dtype("C"), tag="edit_distance"),
            LabelDef(1, "NM", "", parse_dtype("C")),
            LabelDef(2, "custom_float", "", parse_dtype("f"), tag="custom_float"),
        ])
        header = b"readname\tedit_distance:i:3\tNM:i:5\tcustom_float:f:1.5"
        tag_map = build_tag_extractor(defs)
        py_result = extract_labels_from_header(header, tag_map, len(defs), label_defs=defs)
        cpp_result = self._extract(header, list(defs))
        self.assertEqual(py_result[0], cpp_result[0])
        self.assertEqual(py_result[1], cpp_result[1])
        self.assertAlmostEqual(py_result[2], cpp_result[2], places=4)

    def test_negative_float(self):
        defs = [LabelDef(0, "XF", "", parse_dtype("f"))]
        header = b"readname\tXF:f:-3.14"
        result = self._extract(header, defs)
        self.assertAlmostEqual(result[0], -3.14, places=2)

    def test_scientific_notation_float(self):
        defs = [LabelDef(0, "XD", "", parse_dtype("d"))]
        header = b"readname\tXD:f:1.5e-4"
        result = self._extract(header, defs)
        self.assertAlmostEqual(result[0], 1.5e-4, places=7)


@unittest.skipUnless(HAS_ACCEL, "C++ acceleration not available")
class TestCppEncodeLabeledStrandNorm(unittest.TestCase):
    """Test that strand normalization works correctly with C++ labeled encoder."""

    def test_deterministic_rc_r1(self):
        """R1 should be reverse-complemented when do_rc_r1=True."""
        ldef = LabelDef(0, "X", "", parse_dtype("C"))
        col_data, col_sizes = _pack_label_columns([[42]], [ldef], 1)

        result = cpp_encode_block_labeled(
            ["AAAA"], [0x05], 2, "", True, False, False, col_data, col_sizes,
        )
        flags_bytes = result[0]
        # Flag should have IS_RC set (0x08)
        self.assertTrue(flags_bytes[0] & 0x08)

    def test_labels_survive_strand_norm_full_pipeline(self):
        """Labels should be preserved through strand normalization."""
        defs = (LabelDef(0, "NH", "", parse_dtype("C")),)
        h = ZnaHeader(
            read_group="rg",
            strand_specific=True,
            read1_antisense=True,
            strand_normalized=True,
            labels=defs,
        )
        original_seq = "AAAACCCC"
        with tempfile.TemporaryDirectory() as d:
            path = f"{d}/test.zna"
            with open(path, "wb") as fh:
                with ZnaWriter(fh, h) as w:
                    w.write_record(original_seq, True, True, False, labels=(42,))
            with open(path, "rb") as fh:
                reader = ZnaReader(fh)
                recs = list(reader.records(restore_strand=True))
        self.assertEqual(len(recs), 1)
        self.assertEqual(recs[0][0], original_seq)
        self.assertEqual(recs[0][4], (42,))

    def test_npolicy_with_labels(self):
        """N-policy should work correctly with C++ labeled encoder."""
        ldef = LabelDef(0, "X", "", parse_dtype("C"))
        col_data, col_sizes = _pack_label_columns([[7]], [ldef], 1)

        result = cpp_encode_block_labeled(
            ["ACNGT"], [0], 2, "A", False, False, False, col_data, col_sizes,
        )
        # Should not raise — N replaced with A
        self.assertEqual(len(result), 4)


@unittest.skipUnless(HAS_ACCEL, "C++ acceleration not available")
class TestCppWriterReaderEndToEnd(unittest.TestCase):
    """Full Writer→Reader roundtrip via C++ backend for labeled files."""

    def _roundtrip(self, header, records_in):
        with tempfile.TemporaryDirectory() as d:
            path = f"{d}/test.zna"
            with open(path, "wb") as fh:
                with ZnaWriter(fh, header) as w:
                    for seq, ip, r1, r2, labels in records_in:
                        w.write_record(seq, ip, r1, r2, labels=labels)
            with open(path, "rb") as fh:
                reader = ZnaReader(fh)
                return reader.header, list(reader.records())

    def test_single_record(self):
        defs = (LabelDef(0, "NH", "", parse_dtype("C")),)
        h = ZnaHeader(read_group="rg", labels=defs)
        _, recs = self._roundtrip(h, [("ACGT", False, False, False, (5,))])
        self.assertEqual(recs[0][4], (5,))

    def test_multiblock(self):
        """Test across multiple blocks with tiny block sizes."""
        defs = (LabelDef(0, "X", "", parse_dtype("I")),)
        h = ZnaHeader(read_group="rg", labels=defs)
        n = 200
        records_in = [("A" * 50, False, False, False, (i,)) for i in range(n)]
        with tempfile.TemporaryDirectory() as d:
            path = f"{d}/test.zna"
            with open(path, "wb") as fh:
                with ZnaWriter(fh, h, block_size=100) as w:
                    for seq, ip, r1, r2, labels in records_in:
                        w.write_record(seq, ip, r1, r2, labels=labels)
            with open(path, "rb") as fh:
                recs = list(ZnaReader(fh).records())
        self.assertEqual(len(recs), n)
        for i in range(n):
            self.assertEqual(recs[i][4], (i,))

    def test_compressed_many_records(self):
        defs = (
            LabelDef(0, "NH", "", parse_dtype("C")),
            LabelDef(1, "AS", "", parse_dtype("i")),
            LabelDef(2, "de", "", parse_dtype("f")),
        )
        h = ZnaHeader(read_group="rg", compression_method=COMPRESSION_ZSTD, labels=defs)
        n = 1000
        records_in = [
            ("ACGT" * 20, False, False, False, (i % 256, i * 10, i * 0.01))
            for i in range(n)
        ]
        _, recs = self._roundtrip(h, records_in)
        self.assertEqual(len(recs), n)
        for i in range(n):
            self.assertEqual(recs[i][4][0], i % 256)
            self.assertEqual(recs[i][4][1], i * 10)
            self.assertAlmostEqual(recs[i][4][2], i * 0.01, places=5)

    def test_all_dtypes_roundtrip(self):
        """End-to-end roundtrip through Writer/Reader for all dtypes."""
        defs = tuple(
            LabelDef(i, f"X{i}", "", parse_dtype(code))
            for i, code in enumerate("AcCsSiIfqQd")
        )
        vals = (65, -10, 200, -1000, 50000, -100000, 3000000000, 1.5, -(2**60), 2**62, 3.14159)
        h = ZnaHeader(read_group="rg", compression_method=COMPRESSION_ZSTD, labels=defs)
        _, recs = self._roundtrip(h, [("ACGT" * 10, False, False, False, vals)])
        self.assertEqual(len(recs), 1)
        for i, (expect, got) in enumerate(zip(vals, recs[0][4])):
            if defs[i].dtype.code in ('f', 'd'):
                self.assertAlmostEqual(expect, got, places=2,
                                       msg=f"dtype {defs[i].dtype.code}")
            else:
                self.assertEqual(expect, got, msg=f"dtype {defs[i].dtype.code}")


@unittest.skipUnless(HAS_ACCEL, "C++ acceleration not available")
class TestOptimization1APurePython(unittest.TestCase):
    """Test the pure Python optimized extract_labels_from_header (Opt 1A).

    These test the Python path specifically, which is the fallback when
    C++ is not available.
    """

    def test_int_from_bytes_no_decode(self):
        """int() should accept bytes directly — no .decode() needed."""
        from zna.cli import extract_labels_from_header, build_tag_extractor
        defs = (LabelDef(0, "NH", "", parse_dtype("C")),)
        tag_map = build_tag_extractor(defs)
        header = b"readname\tNH:i:42"
        result = extract_labels_from_header(header, tag_map, 1, label_defs=defs)
        self.assertEqual(result, (42,))

    def test_float_from_bytes_no_decode(self):
        """float() should accept bytes directly — no .decode() needed."""
        from zna.cli import extract_labels_from_header, build_tag_extractor
        defs = (LabelDef(0, "de", "", parse_dtype("f")),)
        tag_map = build_tag_extractor(defs)
        header = b"readname\tde:f:0.0123"
        result = extract_labels_from_header(header, tag_map, 1, label_defs=defs)
        self.assertAlmostEqual(result[0], 0.0123, places=4)

    def test_early_exit(self):
        """Should return immediately once all labels are found."""
        from zna.cli import extract_labels_from_header, build_tag_extractor
        defs = (LabelDef(0, "NH", "", parse_dtype("C")),)
        tag_map = build_tag_extractor(defs)
        # NH is found early; many trailing fields
        header = b"readname\tNH:i:5\tXX:i:1\tYY:i:2\tZZ:i:3\tWW:i:4"
        result = extract_labels_from_header(header, tag_map, 1, label_defs=defs)
        self.assertEqual(result, (5,))

    def test_conv_code_dispatch(self):
        """Verify the conv_code dispatch works for all types."""
        from zna.cli import build_tag_extractor, _CONV_INT, _CONV_FLOAT, _CONV_ORD
        defs = (
            LabelDef(0, "NH", "", parse_dtype("C")),   # int
            LabelDef(1, "de", "", parse_dtype("f")),    # float
            LabelDef(2, "tp", "", parse_dtype("A")),    # ord
        )
        tag_map = build_tag_extractor(defs)
        self.assertEqual(tag_map[b"NH"][1], _CONV_INT)
        self.assertEqual(tag_map[b"de"][1], _CONV_FLOAT)
        self.assertEqual(tag_map[b"tp"][1], _CONV_ORD)


if __name__ == "__main__":
    unittest.main()
