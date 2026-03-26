"""Exhaustive test suite for the labeled sequences feature (ZNA V2)."""
import struct
import tempfile
import unittest

from zna.dtypes import (
    LabelDtype, LabelDef, DTYPE_BY_CODE, DTYPE_BY_NAME,
    parse_dtype, label_bytes_per_record,
    resolve_missing, pack_missing, unpack_missing,
)
from zna.core import (
    ZnaHeader, ZnaWriter, ZnaReader, write_zna, read_zna,
    COMPRESSION_ZSTD, COMPRESSION_NONE,
    _FILE_HEADER_FMT, _FILE_HEADER_SIZE, _MAGIC, _VERSION,
)
from zna._pycodec import encode_block, decode_block


# ---------------------------------------------------------------------------
# §1  Dtype registry
# ---------------------------------------------------------------------------

class TestDtypeRegistry(unittest.TestCase):
    """Test the dtype system in dtypes.py."""

    def test_all_numeric_codes_present(self):
        for code in "AcCsSiIfqdQ":
            self.assertIn(code, DTYPE_BY_CODE, f"Missing code {code!r}")

    def test_no_fixstr(self):
        """Z (fixstr) should no longer be in the registry."""
        self.assertNotIn("Z", DTYPE_BY_CODE)
        self.assertNotIn("fixstr", DTYPE_BY_NAME)

    def test_name_lookup(self):
        self.assertIs(DTYPE_BY_NAME["uint8"], DTYPE_BY_CODE["C"])
        self.assertIs(DTYPE_BY_NAME["int32"], DTYPE_BY_CODE["i"])
        self.assertIs(DTYPE_BY_NAME["float32"], DTYPE_BY_CODE["f"])

    def test_parse_dtype_by_code(self):
        dt = parse_dtype("C")
        self.assertEqual(dt.code, "C")
        self.assertEqual(dt.name, "uint8")

    def test_parse_dtype_by_name(self):
        dt = parse_dtype("uint16")
        self.assertEqual(dt.code, "S")

    def test_parse_dtype_unknown_raises(self):
        with self.assertRaises(ValueError):
            parse_dtype("X")

    def test_parse_dtype_fixstr_raises(self):
        with self.assertRaises(ValueError):
            parse_dtype("Z")

    def test_label_bytes_per_record_numeric(self):
        ldef = LabelDef(0, "x", "", parse_dtype("i"))
        self.assertEqual(label_bytes_per_record(ldef), 4)

    def test_struct_sizes(self):
        """Verify struct_ch sizes match declared sizes for all numeric types."""
        for dt in DTYPE_BY_CODE.values():
            self.assertEqual(
                struct.calcsize(dt.struct_ch), dt.size,
                f"Size mismatch for {dt.code}"
            )

    def test_frozen_dtype(self):
        dt = parse_dtype("C")
        with self.assertRaises(AttributeError):
            dt.code = "X"


# ---------------------------------------------------------------------------
# §1b  Missing value system
# ---------------------------------------------------------------------------

class TestMissingValues(unittest.TestCase):
    """Test the missing value resolution and serialization."""

    def test_resolve_missing_default_int(self):
        ldef = LabelDef(0, "NH", "", parse_dtype("C"))
        self.assertEqual(resolve_missing(ldef), 0)

    def test_resolve_missing_default_float(self):
        ldef = LabelDef(0, "X", "", parse_dtype("f"))
        self.assertEqual(resolve_missing(ldef), 0.0)

    def test_resolve_missing_default_float64(self):
        ldef = LabelDef(0, "X", "", parse_dtype("d"))
        self.assertEqual(resolve_missing(ldef), 0.0)

    def test_resolve_missing_explicit(self):
        ldef = LabelDef(0, "NH", "", parse_dtype("C"), missing=255)
        self.assertEqual(resolve_missing(ldef), 255)

    def test_resolve_missing_explicit_float(self):
        ldef = LabelDef(0, "X", "", parse_dtype("f"), missing=-1.0)
        self.assertEqual(resolve_missing(ldef), -1.0)

    def test_resolve_missing_explicit_zero(self):
        """Explicit 0 should be returned (not confused with None)."""
        ldef = LabelDef(0, "NH", "", parse_dtype("C"), missing=0)
        self.assertEqual(resolve_missing(ldef), 0)

    def test_pack_unpack_roundtrip_int(self):
        for code in "AcCsSiIqQ":
            with self.subTest(code=code):
                dt = parse_dtype(code)
                ldef = LabelDef(0, "X", "", dt, missing=42)
                packed = pack_missing(ldef)
                self.assertEqual(len(packed), 8)
                unpacked = unpack_missing(dt, packed)
                self.assertEqual(unpacked, 42)

    def test_pack_unpack_roundtrip_float(self):
        for code in "fd":
            with self.subTest(code=code):
                dt = parse_dtype(code)
                ldef = LabelDef(0, "X", "", dt, missing=-3.14)
                packed = pack_missing(ldef)
                self.assertEqual(len(packed), 8)
                unpacked = unpack_missing(dt, packed)
                self.assertAlmostEqual(unpacked, -3.14, places=2)

    def test_pack_default_is_zero(self):
        ldef = LabelDef(0, "X", "", parse_dtype("i"))
        packed = pack_missing(ldef)
        self.assertEqual(len(packed), 8)
        unpacked = unpack_missing(parse_dtype("i"), packed)
        self.assertEqual(unpacked, 0)

    def test_all_dtypes_default_zero(self):
        """Every dtype should have 0 or 0.0 as its default missing."""
        for code, dt in DTYPE_BY_CODE.items():
            ldef = LabelDef(0, "X", "", dt)
            val = resolve_missing(ldef)
            if code in ('f', 'd'):
                self.assertEqual(val, 0.0)
            else:
                self.assertEqual(val, 0)


# ---------------------------------------------------------------------------
# §2  Header serialisation
# ---------------------------------------------------------------------------

class TestLabeledHeader(unittest.TestCase):
    """Test ZnaHeader with label definitions."""

    def test_header_no_labels(self):
        h = ZnaHeader(read_group="test")
        self.assertEqual(h.num_labels, 0)
        self.assertEqual(h.labels, ())

    def test_header_with_labels(self):
        defs = (
            LabelDef(0, "NH", "Number of hits", parse_dtype("C")),
            LabelDef(1, "AS", "Alignment score", parse_dtype("i")),
        )
        h = ZnaHeader(read_group="rg1", labels=defs)
        self.assertEqual(h.num_labels, 2)
        self.assertEqual(h.labels[0].name, "NH")
        self.assertEqual(h.labels[1].dtype.code, "i")

    def test_header_label_id_mismatch_raises(self):
        defs = (LabelDef(1, "NH", "", parse_dtype("C")),)
        with self.assertRaises(ValueError):
            ZnaHeader(read_group="rg1", labels=defs)

    def test_header_label_name_too_long(self):
        defs = (LabelDef(0, "A" * 17, "", parse_dtype("C")),)
        with self.assertRaises(ValueError):
            ZnaHeader(read_group="rg1", labels=defs)

    def test_header_label_desc_too_long(self):
        defs = (LabelDef(0, "NH", "X" * 65, parse_dtype("C")),)
        with self.assertRaises(ValueError):
            ZnaHeader(read_group="rg1", labels=defs)

    def test_header_roundtrip_no_labels(self):
        h = ZnaHeader(read_group="rg1", description="desc")
        with tempfile.TemporaryDirectory() as d:
            path = f"{d}/test.zna"
            with open(path, "wb") as fh:
                with ZnaWriter(fh, h) as w:
                    w.write_record("ACGT", False, False, False)
            with open(path, "rb") as fh:
                r = ZnaReader(fh)
                self.assertEqual(r.header.num_labels, 0)
                self.assertEqual(r.header.read_group, "rg1")

    def test_header_roundtrip_one_label(self):
        defs = (LabelDef(0, "NH", "Hits count", parse_dtype("C")),)
        h = ZnaHeader(read_group="rg1", labels=defs)
        with tempfile.TemporaryDirectory() as d:
            path = f"{d}/test.zna"
            with open(path, "wb") as fh:
                with ZnaWriter(fh, h) as w:
                    w.write_record("ACGT", False, False, False, labels=(1,))
            with open(path, "rb") as fh:
                r = ZnaReader(fh)
                self.assertEqual(r.header.num_labels, 1)
                self.assertEqual(r.header.labels[0].name, "NH")
                self.assertEqual(r.header.labels[0].description, "Hits count")
                self.assertEqual(r.header.labels[0].dtype.code, "C")

    def test_header_roundtrip_many_labels(self):
        defs = (
            LabelDef(0, "NH", "Hits", parse_dtype("C")),
            LabelDef(1, "AS", "Score", parse_dtype("i")),
            LabelDef(2, "NM", "", parse_dtype("S")),
            LabelDef(3, "XF", "", parse_dtype("f")),
        )
        h = ZnaHeader(read_group="exp1", description="4 labels", labels=defs)
        with tempfile.TemporaryDirectory() as d:
            path = f"{d}/test.zna"
            with open(path, "wb") as fh:
                with ZnaWriter(fh, h) as w:
                    w.write_record(
                        "ACGT", False, False, False,
                        labels=(1, -100, 500, 3.14),
                    )
            with open(path, "rb") as fh:
                r = ZnaReader(fh)
                self.assertEqual(r.header.num_labels, 4)
                self.assertEqual(r.header.labels[3].dtype.code, "f")

    def test_header_roundtrip_with_missing(self):
        """Missing values should survive the header roundtrip."""
        defs = (
            LabelDef(0, "NH", "", parse_dtype("C"), missing=255),
            LabelDef(1, "AS", "", parse_dtype("i"), missing=-1),
            LabelDef(2, "de", "", parse_dtype("f"), missing=-1.0),
        )
        h = ZnaHeader(read_group="rg1", labels=defs)
        with tempfile.TemporaryDirectory() as d:
            path = f"{d}/test.zna"
            with open(path, "wb") as fh:
                with ZnaWriter(fh, h) as w:
                    w.write_record("ACGT", False, False, False, labels=(1, 100, 0.5))
            with open(path, "rb") as fh:
                r = ZnaReader(fh)
                self.assertEqual(r.header.labels[0].missing, 255)
                self.assertEqual(r.header.labels[1].missing, -1)
                self.assertAlmostEqual(r.header.labels[2].missing, -1.0, places=5)

    def test_header_roundtrip_default_missing(self):
        """Labels with no explicit missing should read back as 0."""
        defs = (LabelDef(0, "NH", "", parse_dtype("C")),)
        h = ZnaHeader(read_group="rg1", labels=defs)
        with tempfile.TemporaryDirectory() as d:
            path = f"{d}/test.zna"
            with open(path, "wb") as fh:
                with ZnaWriter(fh, h) as w:
                    w.write_record("ACGT", False, False, False, labels=(1,))
            with open(path, "rb") as fh:
                r = ZnaReader(fh)
                self.assertEqual(r.header.labels[0].missing, 0)

    def test_version_is_2(self):
        """The on-disk version byte should be 2."""
        h = ZnaHeader(read_group="v2")
        with tempfile.TemporaryDirectory() as d:
            path = f"{d}/test.zna"
            with open(path, "wb") as fh:
                with ZnaWriter(fh, h) as w:
                    w.write_record("ACGT", False, False, False)
            with open(path, "rb") as fh:
                fixed = fh.read(_FILE_HEADER_SIZE)
                vals = struct.unpack(_FILE_HEADER_FMT, fixed)
                self.assertEqual(vals[0], _MAGIC)
                self.assertEqual(vals[1], 2)  # version byte


# ---------------------------------------------------------------------------
# §3  Codec roundtrips
# ---------------------------------------------------------------------------

class TestCodecLabels(unittest.TestCase):
    """Test encode_block / decode_block with label columns."""

    def _roundtrip(self, seqs, flags, label_values, label_defs):
        """Encode then decode records with labels."""
        result = encode_block(
            seqs, flags, 2, "", False, False, False,
            label_values=label_values,
            label_defs=label_defs,
        )
        self.assertEqual(len(result), 4)
        flags_bytes, labels_bytes, lengths_bytes, seqs_bytes = result

        # Split label columns
        count = len(seqs)
        lbl_cols = []
        offset = 0
        for ldef in label_defs:
            bpr = label_bytes_per_record(ldef)
            col_size = bpr * count
            lbl_cols.append(labels_bytes[offset:offset+col_size])
            offset += col_size

        records = decode_block(
            flags_bytes, lengths_bytes, seqs_bytes, 2, count,
            label_columns=lbl_cols, label_defs=label_defs,
        )
        return records

    def test_single_uint8(self):
        ldef = LabelDef(0, "NH", "", parse_dtype("C"))
        recs = self._roundtrip(
            ["ACGT"], [0], [[5]], [ldef]
        )
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
        recs = self._roundtrip(["ACGTACGT"], [0], vals, ldefs)
        self.assertEqual(recs[0][5][0], -10)
        self.assertEqual(recs[0][5][1], 60000)
        self.assertEqual(recs[0][5][2], -100000)
        self.assertAlmostEqual(recs[0][5][3], 3.14, places=5)

    def test_multiple_records(self):
        ldef = LabelDef(0, "NH", "", parse_dtype("C"))
        n = 100
        seqs = ["ACGT"] * n
        flags = [0] * n
        vals = [list(range(n))]
        recs = self._roundtrip(seqs, flags, vals, [ldef])
        self.assertEqual(len(recs), n)
        for i in range(n):
            self.assertEqual(recs[i][5], (i,))

    def test_no_labels_returns_3tuple(self):
        """Without label args, encode_block returns 3-tuple."""
        result = encode_block(["ACGT"], [0], 2, "", False, False)
        self.assertEqual(len(result), 3)

    def test_all_numeric_dtype_roundtrips(self):
        """Test every numeric dtype with min/max values."""
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
                recs = self._roundtrip(
                    ["ACGT", "TGCA"], [0, 0],
                    [[min_val, max_val]], [ldef]
                )
                if code in ("f", "d"):
                    self.assertAlmostEqual(recs[0][5][0], min_val, places=1)
                    self.assertAlmostEqual(recs[1][5][0], max_val, places=1)
                else:
                    self.assertEqual(recs[0][5][0], min_val)
                    self.assertEqual(recs[1][5][0], max_val)

    def test_mixed_numeric_types(self):
        ldefs = [
            LabelDef(0, "NH", "", parse_dtype("C")),
            LabelDef(1, "AS", "", parse_dtype("i")),
            LabelDef(2, "de", "", parse_dtype("f")),
        ]
        vals = [[10, 20], [-42, 999], [0.5, 1.5]]
        recs = self._roundtrip(
            ["ACGT", "TGCA"], [0, 0], vals, ldefs
        )
        self.assertEqual(recs[0][5][0], 10)
        self.assertEqual(recs[0][5][1], -42)
        self.assertAlmostEqual(recs[0][5][2], 0.5, places=5)
        self.assertEqual(recs[1][5][0], 20)
        self.assertEqual(recs[1][5][1], 999)
        self.assertAlmostEqual(recs[1][5][2], 1.5, places=5)


# ---------------------------------------------------------------------------
# §4  Writer / Reader end-to-end
# ---------------------------------------------------------------------------

class TestWriterReaderLabels(unittest.TestCase):
    """Full roundtrip tests through ZnaWriter -> ZnaReader with labels."""

    def _roundtrip(self, header, records_in):
        """Write records then read them back.

        records_in: list of (seq, is_paired, is_r1, is_r2, labels_tuple)
        """
        with tempfile.TemporaryDirectory() as d:
            path = f"{d}/test.zna"
            with open(path, "wb") as fh:
                with ZnaWriter(fh, header) as w:
                    for seq, ip, r1, r2, labels in records_in:
                        w.write_record(seq, ip, r1, r2, labels=labels)
            with open(path, "rb") as fh:
                reader = ZnaReader(fh)
                out_header = reader.header
                out_records = list(reader.records())
        return out_header, out_records

    def test_single_record_single_label(self):
        defs = (LabelDef(0, "NH", "", parse_dtype("C")),)
        h = ZnaHeader(read_group="rg", labels=defs)
        _, recs = self._roundtrip(h, [("ACGT", False, False, False, (5,))])
        self.assertEqual(len(recs), 1)
        seq, ip, r1, r2, labels = recs[0]
        self.assertEqual(seq, "ACGT")
        self.assertEqual(labels, (5,))

    def test_many_records(self):
        defs = (
            LabelDef(0, "NH", "", parse_dtype("C")),
            LabelDef(1, "AS", "", parse_dtype("i")),
        )
        h = ZnaHeader(read_group="rg", labels=defs)
        n = 500
        records_in = [
            ("ACGT" * 10, False, False, False, (i % 256, i * 10))
            for i in range(n)
        ]
        _, recs = self._roundtrip(h, records_in)
        self.assertEqual(len(recs), n)
        for i in range(n):
            self.assertEqual(recs[i][4], (i % 256, i * 10))

    def test_missing_labels_raises(self):
        defs = (LabelDef(0, "NH", "", parse_dtype("C")),)
        h = ZnaHeader(read_group="rg", labels=defs)
        with tempfile.TemporaryDirectory() as d:
            path = f"{d}/test.zna"
            with open(path, "wb") as fh:
                with self.assertRaises(ValueError):
                    with ZnaWriter(fh, h) as w:
                        w.write_record("ACGT", False, False, False)

    def test_wrong_label_count_raises(self):
        defs = (
            LabelDef(0, "NH", "", parse_dtype("C")),
            LabelDef(1, "AS", "", parse_dtype("i")),
        )
        h = ZnaHeader(read_group="rg", labels=defs)
        with tempfile.TemporaryDirectory() as d:
            path = f"{d}/test.zna"
            with open(path, "wb") as fh:
                with self.assertRaises(ValueError):
                    with ZnaWriter(fh, h) as w:
                        w.write_record("ACGT", False, False, False, labels=(1,))

    def test_write_records_rejects_labels(self):
        defs = (LabelDef(0, "NH", "", parse_dtype("C")),)
        h = ZnaHeader(read_group="rg", labels=defs)
        with tempfile.TemporaryDirectory() as d:
            path = f"{d}/test.zna"
            with open(path, "wb") as fh:
                with self.assertRaises(TypeError):
                    with ZnaWriter(fh, h) as w:
                        w.write_records([("ACGT", False, False, False)])

    def test_multiblock_labels(self):
        """Test label correctness across multiple blocks."""
        defs = (LabelDef(0, "X", "", parse_dtype("I")),)
        h = ZnaHeader(read_group="rg", labels=defs)
        n = 200
        records_in = [
            ("A" * 50, False, False, False, (i,)) for i in range(n)
        ]
        with tempfile.TemporaryDirectory() as d:
            path = f"{d}/test.zna"
            with open(path, "wb") as fh:
                with ZnaWriter(fh, h, block_size=100) as w:
                    for seq, ip, r1, r2, labels in records_in:
                        w.write_record(seq, ip, r1, r2, labels=labels)
            with open(path, "rb") as fh:
                reader = ZnaReader(fh)
                recs = list(reader.records())
        self.assertEqual(len(recs), n)
        for i in range(n):
            self.assertEqual(recs[i][4], (i,))

    def test_compressed_labels(self):
        defs = (
            LabelDef(0, "NH", "", parse_dtype("C")),
            LabelDef(1, "AS", "", parse_dtype("i")),
        )
        h = ZnaHeader(read_group="rg", compression_method=COMPRESSION_ZSTD, labels=defs)
        records_in = [
            ("ACGT" * 5, False, False, False, (i, i * 100))
            for i in range(50)
        ]
        with tempfile.TemporaryDirectory() as d:
            path = f"{d}/test.zna"
            with open(path, "wb") as fh:
                with ZnaWriter(fh, h) as w:
                    for seq, ip, r1, r2, labels in records_in:
                        w.write_record(seq, ip, r1, r2, labels=labels)
            with open(path, "rb") as fh:
                reader = ZnaReader(fh)
                recs = list(reader.records())
        self.assertEqual(len(recs), 50)
        self.assertEqual(recs[0][4], (0, 0))
        self.assertEqual(recs[49][4], (49, 4900))

    def test_paired_with_labels(self):
        defs = (LabelDef(0, "AS", "", parse_dtype("i")),)
        h = ZnaHeader(read_group="rg", labels=defs)
        records_in = [
            ("AAAA", True, True, False, (100,)),   # R1
            ("CCCC", True, False, True, (200,)),   # R2
        ]
        _, recs = self._roundtrip(h, records_in)
        self.assertEqual(len(recs), 2)
        self.assertEqual(recs[0][0], "AAAA")
        self.assertTrue(recs[0][2])   # is_read1
        self.assertEqual(recs[0][4], (100,))
        self.assertEqual(recs[1][0], "CCCC")
        self.assertTrue(recs[1][3])   # is_read2
        self.assertEqual(recs[1][4], (200,))

    def test_strand_normalized_with_labels(self):
        """Labels survive strand normalization (deterministic)."""
        defs = (LabelDef(0, "NH", "", parse_dtype("C")),)
        h = ZnaHeader(
            read_group="rg",
            strand_specific=True,
            read1_antisense=True,
            strand_normalized=True,
            labels=defs,
        )
        original_seq = "AAAACCCC"
        records_in = [
            (original_seq, True, True, False, (42,)),
        ]
        with tempfile.TemporaryDirectory() as d:
            path = f"{d}/test.zna"
            with open(path, "wb") as fh:
                with ZnaWriter(fh, h) as w:
                    for seq, ip, r1, r2, labels in records_in:
                        w.write_record(seq, ip, r1, r2, labels=labels)
            with open(path, "rb") as fh:
                reader = ZnaReader(fh)
                recs = list(reader.records(restore_strand=True))
        self.assertEqual(len(recs), 1)
        self.assertEqual(recs[0][0], original_seq)
        self.assertEqual(recs[0][4], (42,))

    def test_int64_uint64_labels(self):
        """Test 8-byte integer labels."""
        defs = (
            LabelDef(0, "q", "", parse_dtype("q")),
            LabelDef(1, "Q", "", parse_dtype("Q")),
        )
        h = ZnaHeader(read_group="rg", labels=defs)
        big_signed = -(2**62)
        big_unsigned = 2**63
        _, recs = self._roundtrip(h, [
            ("ACGT", False, False, False, (big_signed, big_unsigned)),
        ])
        self.assertEqual(recs[0][4], (big_signed, big_unsigned))

    def test_float64_label(self):
        defs = (LabelDef(0, "d", "", parse_dtype("d")),)
        h = ZnaHeader(read_group="rg", labels=defs)
        _, recs = self._roundtrip(h, [
            ("ACGT", False, False, False, (1.23456789012345,)),
        ])
        self.assertAlmostEqual(recs[0][4][0], 1.23456789012345, places=10)

    def test_char_label(self):
        """Type 'A' stores a single byte -- test with ASCII code."""
        defs = (LabelDef(0, "A", "", parse_dtype("A")),)
        h = ZnaHeader(read_group="rg", labels=defs)
        _, recs = self._roundtrip(h, [
            ("ACGT", False, False, False, (ord("+"),)),
        ])
        self.assertEqual(recs[0][4], (ord("+"),))


# ---------------------------------------------------------------------------
# §5  Edge cases
# ---------------------------------------------------------------------------

class TestLabelEdgeCases(unittest.TestCase):
    """Edge case tests for the label system."""

    def test_empty_file_with_labels(self):
        """Writing zero records to a labeled file should not error."""
        defs = (LabelDef(0, "NH", "", parse_dtype("C")),)
        h = ZnaHeader(read_group="rg", labels=defs)
        with tempfile.TemporaryDirectory() as d:
            path = f"{d}/test.zna"
            with open(path, "wb") as fh:
                with ZnaWriter(fh, h) as w:
                    pass
            with open(path, "rb") as fh:
                reader = ZnaReader(fh)
                self.assertEqual(reader.header.num_labels, 1)
                recs = list(reader.records())
                self.assertEqual(len(recs), 0)

    def test_many_labels(self):
        """Test with 20 label columns."""
        n_labels = 20
        defs = tuple(
            LabelDef(i, f"L{i:02d}", "", parse_dtype("C")) for i in range(n_labels)
        )
        h = ZnaHeader(read_group="rg", labels=defs)
        values = tuple(range(n_labels))
        with tempfile.TemporaryDirectory() as d:
            path = f"{d}/test.zna"
            with open(path, "wb") as fh:
                with ZnaWriter(fh, h) as w:
                    w.write_record("ACGT", False, False, False, labels=values)
            with open(path, "rb") as fh:
                reader = ZnaReader(fh)
                recs = list(reader.records())
                self.assertEqual(recs[0][4], values)

    def test_zero_value_numeric(self):
        defs = (
            LabelDef(0, "X", "", parse_dtype("i")),
            LabelDef(1, "Y", "", parse_dtype("f")),
        )
        h = ZnaHeader(read_group="rg", labels=defs)
        with tempfile.TemporaryDirectory() as d:
            path = f"{d}/test.zna"
            with open(path, "wb") as fh:
                with ZnaWriter(fh, h) as w:
                    w.write_record("ACGT", False, False, False, labels=(0, 0.0))
            with open(path, "rb") as fh:
                reader = ZnaReader(fh)
                recs = list(reader.records())
                self.assertEqual(recs[0][4][0], 0)
                self.assertEqual(recs[0][4][1], 0.0)


# ---------------------------------------------------------------------------
# §6  CLI label spec parsing
# ---------------------------------------------------------------------------

class TestCliLabelParsing(unittest.TestCase):
    """Test CLI label argument parsing."""

    def test_parse_label_spec_numeric(self):
        from zna.cli import parse_label_spec
        name, dtype_str, tag = parse_label_spec("NH:C")
        self.assertEqual(name, "NH")
        self.assertEqual(dtype_str, "C")
        self.assertIsNone(tag)

    def test_parse_label_spec_name(self):
        from zna.cli import parse_label_spec
        name, dtype_str, tag = parse_label_spec("AS:int32")
        self.assertEqual(name, "AS")
        self.assertEqual(dtype_str, "int32")
        self.assertIsNone(tag)

    def test_parse_label_spec_three_part(self):
        """NAME:TYPE:TAG format should return tag."""
        from zna.cli import parse_label_spec
        name, dtype_str, tag = parse_label_spec("edits:C:NM")
        self.assertEqual(name, "edits")
        self.assertEqual(dtype_str, "C")
        self.assertEqual(tag, "NM")

    def test_parse_label_spec_three_part_empty_tag(self):
        """NAME:TYPE: with empty tag should be rejected."""
        import argparse
        from zna.cli import parse_label_spec
        with self.assertRaises(argparse.ArgumentTypeError):
            parse_label_spec("edits:C:")

    def test_parse_label_spec_four_colons_rejected(self):
        import argparse
        from zna.cli import parse_label_spec
        with self.assertRaises(argparse.ArgumentTypeError):
            parse_label_spec("a:b:c:d")

    def test_parse_label_spec_invalid_no_colon(self):
        import argparse
        from zna.cli import parse_label_spec
        with self.assertRaises(argparse.ArgumentTypeError):
            parse_label_spec("NH")

    def test_parse_label_spec_fixstr_rejected(self):
        """Z (fixstr) should be rejected as unknown dtype."""
        import argparse
        from zna.cli import parse_label_spec
        with self.assertRaises(argparse.ArgumentTypeError):
            parse_label_spec("MD:Z32")

    def test_parse_label_spec_unknown_dtype(self):
        import argparse
        from zna.cli import parse_label_spec
        with self.assertRaises(argparse.ArgumentTypeError):
            parse_label_spec("X:blah")

    def test_build_label_defs(self):
        from zna.cli import build_label_defs
        defs = build_label_defs(
            ["NH:C", "AS:i"],
            ["NH:Number of hits"],
        )
        self.assertEqual(len(defs), 2)
        self.assertEqual(defs[0].name, "NH")
        self.assertEqual(defs[0].dtype.code, "C")
        self.assertEqual(defs[0].description, "Number of hits")
        self.assertIsNone(defs[0].tag)
        self.assertEqual(defs[0].effective_tag, "NH")
        self.assertEqual(defs[1].name, "AS")
        self.assertEqual(defs[1].dtype.code, "i")

    def test_build_label_defs_with_tag(self):
        from zna.cli import build_label_defs
        defs = build_label_defs(
            ["edit_dist:C:NM", "score:i:AS"],
            ["edit_dist:Number of edits"],
        )
        self.assertEqual(len(defs), 2)
        self.assertEqual(defs[0].name, "edit_dist")
        self.assertEqual(defs[0].tag, "NM")
        self.assertEqual(defs[0].effective_tag, "NM")
        self.assertEqual(defs[0].description, "Number of edits")
        self.assertEqual(defs[1].name, "score")
        self.assertEqual(defs[1].tag, "AS")
        self.assertEqual(defs[1].effective_tag, "AS")

    def test_extract_with_tag_different_from_name(self):
        """When tag differs from name, the tag is used for input parsing."""
        from zna.cli import extract_labels_from_header, build_tag_extractor
        defs = (
            LabelDef(0, "edit_dist", "", parse_dtype("C"), tag="NM"),
            LabelDef(1, "aln_score", "", parse_dtype("i"), tag="AS"),
        )
        tag_map = build_tag_extractor(defs)
        # Input has SAM tags NM and AS (not edit_dist/aln_score)
        raw = b"readname\tNM:i:3\tAS:i:298"
        labels = extract_labels_from_header(raw, tag_map, 2, label_defs=defs)
        self.assertEqual(labels, (3, 298))

    def test_extract_labels_from_header_tab_separated(self):
        from zna.cli import extract_labels_from_header, build_tag_extractor
        defs = (
            LabelDef(0, "NH", "", parse_dtype("C")),
            LabelDef(1, "AS", "", parse_dtype("i")),
        )
        tag_map = build_tag_extractor(defs)
        raw = b"readname\tNH:i:1\tAS:i:298"
        labels = extract_labels_from_header(raw, tag_map, 2, label_defs=defs)
        self.assertEqual(labels, (1, 298))

    def test_extract_labels_from_header_space_separated(self):
        """Fields separated by spaces should also work."""
        from zna.cli import extract_labels_from_header, build_tag_extractor
        defs = (LabelDef(0, "NH", "", parse_dtype("C")),)
        tag_map = build_tag_extractor(defs)
        raw = b"readname NH:i:5"
        labels = extract_labels_from_header(raw, tag_map, 1, label_defs=defs)
        self.assertEqual(labels, (5,))

    def test_extract_labels_fastp_merged_suffix(self):
        """fastp merged_XX_YY suffix should not corrupt tag values."""
        from zna.cli import extract_labels_from_header, build_tag_extractor
        defs = (
            LabelDef(0, "AS", "", parse_dtype("i")),
            LabelDef(1, "rl", "", parse_dtype("S")),
        )
        tag_map = build_tag_extractor(defs)
        raw = b"readname\tAS:i:100\trl:i:0 merged_150_87"
        labels = extract_labels_from_header(raw, tag_map, 2, label_defs=defs)
        self.assertEqual(labels, (100, 0))

    def test_extract_labels_missing_tag_with_label_defs(self):
        """Missing tags should be filled with the label's missing value."""
        from zna.cli import extract_labels_from_header, build_tag_extractor
        defs = (
            LabelDef(0, "NH", "", parse_dtype("C"), missing=255),
            LabelDef(1, "AS", "", parse_dtype("i")),
        )
        tag_map = build_tag_extractor(defs)
        raw = b"readname\tAS:i:100"  # NH is missing
        labels = extract_labels_from_header(raw, tag_map, 2, label_defs=defs)
        self.assertEqual(labels, (255, 100))

    def test_extract_labels_missing_tag_default_zero(self):
        """Missing tags with no explicit missing should default to 0."""
        from zna.cli import extract_labels_from_header, build_tag_extractor
        defs = (LabelDef(0, "NH", "", parse_dtype("C")),)
        tag_map = build_tag_extractor(defs)
        raw = b"readname\tAS:i:100"  # NH is missing
        labels = extract_labels_from_header(raw, tag_map, 1, label_defs=defs)
        self.assertEqual(labels, (0,))

    def test_extract_labels_missing_tag_raises_without_defs(self):
        """Without label_defs, missing tags should raise."""
        from zna.cli import extract_labels_from_header, build_tag_extractor
        defs = (LabelDef(0, "NH", "", parse_dtype("C")),)
        tag_map = build_tag_extractor(defs)
        raw = b"readname\tAS:i:100"
        with self.assertRaises(ValueError):
            extract_labels_from_header(raw, tag_map, 1)

    def test_extract_labels_float_tag(self):
        from zna.cli import extract_labels_from_header, build_tag_extractor
        defs = (LabelDef(0, "XF", "", parse_dtype("f")),)
        tag_map = build_tag_extractor(defs)
        raw = b"readname\tXF:f:3.14"
        labels = extract_labels_from_header(raw, tag_map, 1, label_defs=defs)
        self.assertAlmostEqual(labels[0], 3.14, places=2)

    def test_extract_labels_ignores_unknown_tags(self):
        from zna.cli import extract_labels_from_header, build_tag_extractor
        defs = (LabelDef(0, "NH", "", parse_dtype("C")),)
        tag_map = build_tag_extractor(defs)
        raw = b"readname\tXX:i:9999\tNH:i:5\tYY:i:42"
        labels = extract_labels_from_header(raw, tag_map, 1, label_defs=defs)
        self.assertEqual(labels, (5,))

    def test_extract_labels_ignores_non_tag_fields(self):
        """Tokens like 'merged_150_87' have no colons and should be skipped."""
        from zna.cli import extract_labels_from_header, build_tag_extractor
        defs = (LabelDef(0, "NH", "", parse_dtype("C")),)
        tag_map = build_tag_extractor(defs)
        raw = b"readname\tNH:i:3\tmerged_150_87"
        labels = extract_labels_from_header(raw, tag_map, 1, label_defs=defs)
        self.assertEqual(labels, (3,))

    def test_extract_labels_long_tag_names(self):
        """Tags longer than 2 characters should work."""
        from zna.cli import extract_labels_from_header, build_tag_extractor
        defs = (
            LabelDef(0, "edit_distance", "", parse_dtype("C"), tag="edit_distance"),
            LabelDef(1, "score", "", parse_dtype("i"), tag="score"),
        )
        tag_map = build_tag_extractor(defs)
        raw = b"readname\tedit_distance:i:3\tscore:i:280"
        labels = extract_labels_from_header(raw, tag_map, 2, label_defs=defs)
        self.assertEqual(labels, (3, 280))

    def test_extract_labels_mixed_tag_lengths(self):
        """Mix of 2-char SAM tags and long custom tags."""
        from zna.cli import extract_labels_from_header, build_tag_extractor
        defs = (
            LabelDef(0, "NM", "", parse_dtype("C")),
            LabelDef(1, "my_score", "", parse_dtype("i"), tag="my_score"),
        )
        tag_map = build_tag_extractor(defs)
        raw = b"readname\tNM:i:5\tmy_score:i:42"
        labels = extract_labels_from_header(raw, tag_map, 2, label_defs=defs)
        self.assertEqual(labels, (5, 42))

    def test_extract_labels_long_tag_with_decoupled_name(self):
        """Long tag with different name should parse correctly."""
        from zna.cli import extract_labels_from_header, build_tag_extractor
        defs = (
            LabelDef(0, "edits", "", parse_dtype("C"), tag="edit_distance"),
        )
        tag_map = build_tag_extractor(defs)
        raw = b"readname\tedit_distance:i:7"
        labels = extract_labels_from_header(raw, tag_map, 1, label_defs=defs)
        self.assertEqual(labels, (7,))


# ---------------------------------------------------------------------------
# §7  Shuffle with labels
# ---------------------------------------------------------------------------

class TestShuffleLabels(unittest.TestCase):
    """Test that shuffle preserves labels."""

    def test_shuffle_preserves_labels(self):
        from zna._shuffle import shuffle_zna

        defs = (
            LabelDef(0, "NH", "", parse_dtype("C")),
            LabelDef(1, "AS", "", parse_dtype("i")),
        )
        h = ZnaHeader(read_group="rg", labels=defs)
        n = 50
        records_in = [
            ("ACGT" * 5, False, False, False, (i % 256, i * 10))
            for i in range(n)
        ]

        with tempfile.TemporaryDirectory() as d:
            in_path = f"{d}/in.zna"
            out_path = f"{d}/out.zna"

            # Write input file
            with open(in_path, "wb") as fh:
                with ZnaWriter(fh, h) as w:
                    for seq, ip, r1, r2, labels in records_in:
                        w.write_record(seq, ip, r1, r2, labels=labels)

            # Shuffle
            shuffle_zna(in_path, out_path, seed=123, quiet=True)

            # Read shuffled output
            with open(out_path, "rb") as fh:
                reader = ZnaReader(fh)
                recs = list(reader.records())

            # All records should be present (possibly in different order)
            self.assertEqual(len(recs), n)
            # Check that label values are consistent with sequences
            original_set = set()
            for seq, ip, r1, r2, labels in records_in:
                original_set.add((seq, labels))
            shuffled_set = set()
            for rec in recs:
                shuffled_set.add((rec[0], rec[4]))
            self.assertEqual(original_set, shuffled_set)


if __name__ == "__main__":
    unittest.main()
