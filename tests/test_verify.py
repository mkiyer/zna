"""``zna inspect --verify``: full certification, and the damage matrix.

One command owns "look at a file" (TRAILER_PLAN §14 A3): bare ``inspect`` runs
the cheap structural pass and reports provenance; ``--verify`` decodes every
block, recounts every stat, and compares against the trailer — exit 0 iff the
file certifies.  Each row of the damage matrix below is a distinct failure
mode with a distinct detecting mechanism, and each must fail loudly.
"""
import io
import json
import struct
import subprocess
import sys
import tempfile
import unittest
import zlib
from pathlib import Path

import zstandard

from zna.core import (
    COMPRESSION_ZSTD,
    ZnaHeader,
    ZnaReader,
    ZnaWriter,
    _BLOCK_HEADER_SIZE,
    _FOOTER_FMT,
    _FOOTER_SIZE,
)


def build(n_pairs=60, n_single=8, block_size=1024) -> bytes:
    import random
    rng = random.Random(7)
    buf = io.BytesIO()
    h = ZnaHeader(read_group="verify-test", seq_len_bytes=2,
                  compression_method=COMPRESSION_ZSTD)
    with ZnaWriter(buf, h, block_size=block_size) as w:
        for _ in range(n_pairs):
            w.write_record("".join(rng.choices("ACGT", k=rng.randint(40, 120))),
                           True, True, False)
            w.write_record("".join(rng.choices("ACGT", k=rng.randint(40, 120))),
                           True, False, True)
        for _ in range(n_single):
            w.write_record("".join(rng.choices("ACGT", k=rng.randint(40, 200))),
                           False, False, False, is_full_fragment=True)
    return buf.getvalue()


def run_inspect(data: bytes, *flags: str) -> subprocess.CompletedProcess:
    with tempfile.TemporaryDirectory() as d:
        path = Path(d) / "t.zna"
        path.write_bytes(data)
        return subprocess.run(
            [sys.executable, "-m", "zna.cli", "inspect", *flags, str(path)],
            capture_output=True, text=True)


class TestVerifyPasses(unittest.TestCase):
    def test_a_healthy_file_certifies(self):
        r = run_inspect(build(), "--verify")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("PASSED", r.stdout)

    def test_a_zero_record_file_certifies(self):
        buf = io.BytesIO()
        with ZnaWriter(buf, ZnaHeader(read_group="e", seq_len_bytes=2,
                                      compression_method=COMPRESSION_ZSTD)):
            pass
        r = run_inspect(buf.getvalue(), "--verify")
        self.assertEqual(r.returncode, 0, r.stderr)

    def test_json_carries_the_payloads_and_the_verdict(self):
        r = run_inspect(build(), "--json", "--verify")
        self.assertEqual(r.returncode, 0, r.stderr)
        out = json.loads(r.stdout)
        self.assertEqual(out["provenance"]["shuffled"], False)
        self.assertIn("writer_version", out["provenance"])
        self.assertEqual(out["trailer"]["n_records"], 128)
        self.assertTrue(out["verify"]["passed"])

    def test_bare_inspect_reports_provenance(self):
        r = run_inspect(build())
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("Provenance", r.stdout)
        self.assertIn("Shuffled:         False", r.stdout)


class TestDamageMatrix(unittest.TestCase):
    """Every row: a distinct failure, a distinct mechanism, a loud refusal."""

    @classmethod
    def setUpClass(cls):
        cls.data = build()
        (cls.magic, cls.fver, _r, cls.crc, cls.sentinel) = struct.unpack(
            _FOOTER_FMT, cls.data[-_FOOTER_SIZE:])

    def assert_fails(self, blob, *needles, flags=("--verify",)):
        r = run_inspect(blob, *flags)
        self.assertNotEqual(r.returncode, 0, "damage was certified:\n" + r.stdout)
        text = r.stdout + r.stderr
        for n in needles:
            self.assertIn(n, text, text)

    def test_truncated_at_a_block_boundary(self):
        """The previously SILENT case: the file reads back clean, just shorter.
        Only the trailer's absence betrays it."""
        self.assert_fails(self.data[:self.sentinel], "no trailer")

    def test_truncated_mid_block(self):
        self.assert_fails(self.data[:self.sentinel - 100], "file ends at")

    def test_flipped_trailer_payload_byte(self):
        b = bytearray(self.data)
        b[self.sentinel + _BLOCK_HEADER_SIZE + 6] ^= 0xFF
        self.assert_fails(bytes(b), "CRC mismatch")

    def test_flipped_data_block_byte(self):
        """D6's content checksums: a single flipped byte inside a compressed
        frame fails the decompress itself."""
        b = bytearray(self.data)
        b[len(self.data) // 3] ^= 0xFF
        self.assert_fails(bytes(b), "block")

    def test_tampered_stats(self):
        payload = self.data[self.sentinel + _BLOCK_HEADER_SIZE:-_FOOTER_SIZE]
        raw = zstandard.ZstdDecompressor().decompress(
            payload, max_output_size=1 << 20)
        d = json.loads(raw)
        d["n_records"] -= 1
        raw2 = json.dumps(d, sort_keys=True, separators=(",", ":")).encode()
        pay2 = zstandard.ZstdCompressor(
            level=9, write_checksum=True).compress(raw2)
        blob = (self.data[:self.sentinel]
                + struct.pack("<IIIII", len(pay2) + _FOOTER_SIZE, len(raw2),
                              0, 0, 0)
                + pay2
                + struct.pack(_FOOTER_FMT, self.magic, 1, 0,
                              zlib.crc32(pay2) & 0xFFFFFFFF, self.sentinel))
        self.assert_fails(blob, "does not match a recount")

    def test_aborted_encode_refuses(self):
        buf = io.BytesIO()
        with self.assertRaises(RuntimeError):
            with ZnaWriter(buf, ZnaHeader(read_group="a", seq_len_bytes=2,
                                          compression_method=COMPRESSION_ZSTD)) as w:
                w.write_record("ACGT", False, False, False)
                raise RuntimeError("crash")
        self.assert_fails(buf.getvalue(), "no trailer")

    def test_bare_inspect_catches_index_mismatch(self):
        """The free structural tier: no --verify needed for gross damage."""
        # Rewrite the stored index to claim a wrong comp_size.
        payload = self.data[self.sentinel + _BLOCK_HEADER_SIZE:-_FOOTER_SIZE]
        raw = zstandard.ZstdDecompressor().decompress(
            payload, max_output_size=1 << 20)
        d = json.loads(raw)
        d["blocks"]["comp_sizes"][0] += 4
        raw2 = json.dumps(d, sort_keys=True, separators=(",", ":")).encode()
        pay2 = zstandard.ZstdCompressor(
            level=9, write_checksum=True).compress(raw2)
        blob = (self.data[:self.sentinel]
                + struct.pack("<IIIII", len(pay2) + _FOOTER_SIZE, len(raw2),
                              0, 0, 0)
                + pay2
                + struct.pack(_FOOTER_FMT, self.magic, 1, 0,
                              zlib.crc32(pay2) & 0xFFFFFFFF, self.sentinel))
        self.assert_fails(blob, "disagrees with the file", flags=())


class TestVerifyEngineDirect(unittest.TestCase):
    """The engine, without a subprocess: check list and warning semantics."""

    def test_checks_enumerate_in_order(self):
        from zna.cli import _verify_zna
        data = build()
        fh = io.BytesIO(data)
        v = _verify_zna(fh, ZnaReader(fh))
        self.assertTrue(v["passed"])
        self.assertEqual(len(v["checks"]), 4)
        self.assertIsNone(v["failure"])

    def test_both_antisense_is_a_warning_not_a_failure(self):
        """The fuzz matrix writes this config; verify must not fail legal
        files (TRAILER_PLAN §14 A4)."""
        buf = io.BytesIO()
        h = ZnaHeader(read_group="w", seq_len_bytes=2,
                      compression_method=COMPRESSION_ZSTD,
                      strand_specific=True, read1_antisense=True,
                      read2_antisense=True, strand_normalized=True)
        with ZnaWriter(buf, h) as w:
            w.write_record("ACGTACGT", True, True, False)
            w.write_record("TTGGCCAA", True, False, True)
        from zna.cli import _verify_zna
        fh = io.BytesIO(buf.getvalue())
        v = _verify_zna(fh, ZnaReader(fh))
        self.assertTrue(v["passed"])
        self.assertEqual(len(v["warnings"]), 1)


if __name__ == "__main__":
    unittest.main()
