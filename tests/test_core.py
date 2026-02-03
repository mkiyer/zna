import sys
import tempfile
import unittest
from pathlib import Path
from zna.core import ZnaHeader, write_zna, read_zna, ZnaWriter, COMPRESSION_ZSTD, COMPRESSION_NONE



class TestZnaIO(unittest.TestCase):
    def roundtrip(self, sequences, header):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = f"{tmpdir}/sample.zna"
            with open(path, "wb") as out_fh:
                # Mock records logic for test:
                # i=0: paired=True, r1=True, r2=False (Paired R1)
                # i=1: paired=True, r1=False, r2=True (Paired R2)
                # i=2: paired=False, r1=False, r2=False (Singleton)
                records = []
                for i, seq in enumerate(sequences):
                    mode = i % 3
                    if mode == 0:
                        records.append((seq, True, True, False))
                    elif mode == 1:
                        records.append((seq, True, False, True))
                    else:
                        records.append((seq, False, False, False))
                        
                write_zna(out_fh, header, records)
            with open(path, "rb") as in_fh:
                read_header, rec_iter = read_zna(in_fh)
                read_records = list(rec_iter)
        return read_header, read_records, records

    def test_roundtrip_header(self):
        header = ZnaHeader(
            read_group="rg1",
            description="test",
            seq_len_bytes=1,
            strand_specific=True,
        )
        sequences = ["A", "AC", "ACG", "ACGT"]
        read_header, _, _ = self.roundtrip(sequences, header)
        self.assertEqual(read_header, header)

    def test_roundtrip_sequences_varlen(self):
        header = ZnaHeader(read_group="", description="", seq_len_bytes=1)
        sequences = [
            "A",
            "AC",
            "ACG",
            "ACGT",
            "ACGTA",
            "ACGTAC",
            "ACGTACG",
            "ACGTACGT",
            "ACGTACGTA",
        ]
        _, records, original_records = self.roundtrip(sequences, header)
        
        # Check sequences match
        decoded_seqs = [r[0] for r in records]
        self.assertEqual(decoded_seqs, sequences)
        
        # Check flags match original inputs
        self.assertEqual(records, original_records)

    def test_flags(self):
        header = ZnaHeader(read_group="", description="")
        sequences = ["ACGT", "TGCA", "AAAA"]
        read_header, records, _ = self.roundtrip(sequences, header)
        self.assertEqual(read_header.read_group, "")
        
        # i=0: paired=True, r1=True, r2=False 
        self.assertTrue(records[0][1]) # is_paired
        self.assertTrue(records[0][2]) # is_read1
        self.assertFalse(records[0][3]) # is_read2
        
        # i=1: paired=True, r1=False, r2=True
        self.assertTrue(records[1][1])
        self.assertFalse(records[1][2])
        self.assertTrue(records[1][3])
        
        # i=2: paired=False, r1=False, r2=False
        self.assertFalse(records[2][1])
        self.assertFalse(records[2][2])
        self.assertFalse(records[2][3])

    def test_invalid_base_raises(self):
        header = ZnaHeader(read_group="", description="")
        with tempfile.TemporaryDirectory() as tmpdir:
            path = f"{tmpdir}/bad.zna"
            with open(path, "wb") as out_fh:
                with self.assertRaises(ValueError):
                    write_zna(out_fh, header, [("ACNG", True, True, False)])

    def test_max_len_exceeded(self):
        # 1 byte length = max 255 bases
        header = ZnaHeader(read_group="", description="", seq_len_bytes=1)
        long_seq = "A" * 256
        with tempfile.TemporaryDirectory() as tmpdir:
            path = f"{tmpdir}/bad_len.zna"
            with open(path, "wb") as out_fh:
                with self.assertRaises(ValueError):
                    write_zna(out_fh, header, [(long_seq, True, True, False)])

    def test_empty_sequence(self):
        header = ZnaHeader(read_group="", description="", seq_len_bytes=1)
        sequences = [""]
        _, records, _ = self.roundtrip(sequences, header)
        self.assertEqual(records[0][0], "")

    def test_header_validation(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = f"{tmpdir}/bad_header.zna"
            with open(path, "wb") as out_fh:
                with self.assertRaises(ValueError):
                    # seq_len_bytes must be 1-4
                    bad_header = ZnaHeader(read_group="", seq_len_bytes=5)
                    write_zna(out_fh, bad_header, [])

    def test_compression_zstd_single_block(self):
        """Test basic Zstandard compression with a single block."""
        header = ZnaHeader(
            read_group="rg_zstd", 
            description="compressed", 
            compression_method=COMPRESSION_ZSTD
        )
        sequences = ["ACGT", "TGCA"] * 50
        read_header, records, original_records = self.roundtrip(sequences, header)
        
        self.assertEqual(read_header.compression_method, COMPRESSION_ZSTD)
        self.assertEqual(len(records), 100)
        self.assertEqual([r[0] for r in records], sequences)

    def test_compression_zstd_multi_block(self):
        """Test Zstandard compression with multiple blocks by forcing small block size."""
        header = ZnaHeader(
            read_group="rg_multi", 
            compression_method=COMPRESSION_ZSTD
        )
        # Create sequences
        # 50 bases -> ~13 bytes packed + 1 flag + 1 len = 15 bytes per record
        # A block size of 40 bytes will fit approx 2 records.
        sequences = ["A" * 50] * 20
        
        with tempfile.TemporaryDirectory() as tmpdir:
            path = f"{tmpdir}/multiblock_zstd.zna"
            
            # Write with small block size
            with open(path, "wb") as fh:
                with ZnaWriter(fh, header, block_size=40) as writer:
                    for i, seq in enumerate(sequences):
                        writer.write_record(seq, False, False, False)
            
            # Read back
            with open(path, "rb") as fh:
                read_header, r_iter = read_zna(fh)
                read_records = list(r_iter)
            
            self.assertEqual(read_header.compression_method, COMPRESSION_ZSTD)
            self.assertEqual(len(read_records), 20)
            self.assertEqual(read_records[0][0], sequences[0])
            self.assertEqual(read_records[-1][0], sequences[-1])

    def test_mixed_compression_switch(self):
        """Ensure we can handle explicitly setting NONE as well."""
        header = ZnaHeader(read_group="none", compression_method=COMPRESSION_NONE)
        sequences = ["A"] * 10
        read_header, _, _ = self.roundtrip(sequences, header)
        self.assertEqual(read_header.compression_method, COMPRESSION_NONE)



if __name__ == "__main__":
    unittest.main()
