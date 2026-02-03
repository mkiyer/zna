import sys
import tempfile
import unittest
from pathlib import Path
from zna.core import ZnaHeader, write_zna, read_zna, ZnaWriter, ZnaReader, COMPRESSION_ZSTD, COMPRESSION_NONE
from zna import reverse_complement



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


class TestReverseComplement(unittest.TestCase):
    """Test the reverse complement function."""
    
    def test_simple_sequence(self):
        """Test basic reverse complement."""
        self.assertEqual(reverse_complement("ACGT"), "ACGT")  # palindrome
        self.assertEqual(reverse_complement("AAAA"), "TTTT")
        self.assertEqual(reverse_complement("CCCC"), "GGGG")
        self.assertEqual(reverse_complement("GGGG"), "CCCC")
        self.assertEqual(reverse_complement("TTTT"), "AAAA")
    
    def test_longer_sequence(self):
        """Test longer sequences."""
        self.assertEqual(reverse_complement("ACGTACGT"), "ACGTACGT")  # palindrome
        self.assertEqual(reverse_complement("AACCGGTT"), "AACCGGTT")  # palindrome
        self.assertEqual(reverse_complement("AAAACCCC"), "GGGGTTTT")
    
    def test_mixed_case(self):
        """Test that mixed case is handled correctly."""
        self.assertEqual(reverse_complement("AcGt"), "aCgT")
        self.assertEqual(reverse_complement("acgt"), "acgt")
    
    def test_empty_sequence(self):
        """Test empty sequence."""
        self.assertEqual(reverse_complement(""), "")
    
    def test_single_base(self):
        """Test single base."""
        self.assertEqual(reverse_complement("A"), "T")
        self.assertEqual(reverse_complement("C"), "G")
        self.assertEqual(reverse_complement("G"), "C")
        self.assertEqual(reverse_complement("T"), "A")


class TestStrandSpecific(unittest.TestCase):
    """Test strand-specific functionality in ZnaWriter and ZnaReader."""
    
    def test_header_antisense_flags(self):
        """Test that antisense flags are properly stored in header."""
        header = ZnaHeader(
            read_group="strand_test",
            description="dUTP protocol",
            strand_specific=True,
            read1_antisense=True,
            read2_antisense=False,
        )
        
        with tempfile.TemporaryDirectory() as tmpdir:
            path = f"{tmpdir}/strand.zna"
            
            # Write
            with open(path, "wb") as fh:
                with ZnaWriter(fh, header) as writer:
                    writer.write_record("ACGT", True, True, False)  # R1
            
            # Read and check header
            with open(path, "rb") as fh:
                reader = ZnaReader(fh)
                self.assertTrue(reader.header.strand_specific)
                self.assertTrue(reader.header.read1_antisense)
                self.assertFalse(reader.header.read2_antisense)
    
    def test_r1_antisense_normalization(self):
        """Test that R1 antisense reads are normalized (reverse complemented) during encoding."""
        header = ZnaHeader(
            read_group="dutp",
            strand_specific=True,
            read1_antisense=True,  # R1 is antisense
            read2_antisense=False,
        )
        
        r1_seq = "AAAACCCC"  # Original antisense R1
        r1_rc = "GGGGTTTT"   # Should be stored as this (sense)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            path = f"{tmpdir}/strand.zna"
            
            # Write R1 (should be reverse complemented)
            with open(path, "wb") as fh:
                with ZnaWriter(fh, header) as writer:
                    writer.write_record(r1_seq, True, True, False)  # R1
            
            # Read without strand restoration
            with open(path, "rb") as fh:
                reader = ZnaReader(fh)
                records = list(reader.records(restore_strand=False))
                self.assertEqual(len(records), 1)
                # Stored as reverse complement (sense strand)
                self.assertEqual(records[0][0], r1_rc)
    
    def test_r2_antisense_normalization(self):
        """Test that R2 antisense reads are normalized during encoding."""
        header = ZnaHeader(
            read_group="fr-secondstrand",
            strand_specific=True,
            read1_antisense=False,
            read2_antisense=True,  # R2 is antisense
        )
        
        r2_seq = "AAAACCCC"
        r2_rc = "GGGGTTTT"
        
        with tempfile.TemporaryDirectory() as tmpdir:
            path = f"{tmpdir}/strand.zna"
            
            # Write R2 (should be reverse complemented)
            with open(path, "wb") as fh:
                with ZnaWriter(fh, header) as writer:
                    writer.write_record(r2_seq, True, False, True)  # R2
            
            # Read without strand restoration
            with open(path, "rb") as fh:
                reader = ZnaReader(fh)
                records = list(reader.records(restore_strand=False))
                self.assertEqual(records[0][0], r2_rc)
    
    def test_sense_reads_not_modified(self):
        """Test that sense reads are not reverse complemented."""
        header = ZnaHeader(
            read_group="dutp",
            strand_specific=True,
            read1_antisense=True,  # R1 is antisense
            read2_antisense=False, # R2 is sense
        )
        
        r2_seq = "AAAACCCC"  # R2 is sense, should NOT be modified
        
        with tempfile.TemporaryDirectory() as tmpdir:
            path = f"{tmpdir}/strand.zna"
            
            # Write R2 (should NOT be reverse complemented)
            with open(path, "wb") as fh:
                with ZnaWriter(fh, header) as writer:
                    writer.write_record(r2_seq, True, False, True)  # R2
            
            # Read - should be unchanged
            with open(path, "rb") as fh:
                reader = ZnaReader(fh)
                records = list(reader.records(restore_strand=False))
                self.assertEqual(records[0][0], r2_seq)
    
    def test_restore_strand_roundtrip(self):
        """Test that restore_strand=True recovers original sequences."""
        header = ZnaHeader(
            read_group="dutp",
            strand_specific=True,
            read1_antisense=True,
            read2_antisense=False,
        )
        
        r1_original = "AAAACCCC"
        r2_original = "TTTTGGGG"
        
        with tempfile.TemporaryDirectory() as tmpdir:
            path = f"{tmpdir}/strand.zna"
            
            # Write paired reads
            with open(path, "wb") as fh:
                with ZnaWriter(fh, header) as writer:
                    writer.write_record(r1_original, True, True, False)  # R1
                    writer.write_record(r2_original, True, False, True)  # R2
            
            # Read with strand restoration
            with open(path, "rb") as fh:
                reader = ZnaReader(fh)
                records = list(reader.records(restore_strand=True))
                
                self.assertEqual(len(records), 2)
                # R1 should be restored to original (reverse complemented back)
                self.assertEqual(records[0][0], r1_original)
                # R2 should be unchanged (it was sense, not modified)
                self.assertEqual(records[1][0], r2_original)
    
    def test_dutp_protocol_paired_roundtrip(self):
        """Test complete dUTP protocol roundtrip with paired-end reads."""
        header = ZnaHeader(
            read_group="dutp_full",
            strand_specific=True,
            read1_antisense=True,  # dUTP: R1 is antisense
            read2_antisense=False, # dUTP: R2 is sense
        )
        
        # Create test sequences
        test_pairs = [
            ("AAACCCGGG", "TTTAAACCC"),  # Pair 1
            ("GGGCCCAAA", "CCCGGGTTT"),  # Pair 2
        ]
        
        with tempfile.TemporaryDirectory() as tmpdir:
            path = f"{tmpdir}/dutp.zna"
            
            # Write all records
            with open(path, "wb") as fh:
                with ZnaWriter(fh, header) as writer:
                    for r1, r2 in test_pairs:
                        writer.write_record(r1, True, True, False)
                        writer.write_record(r2, True, False, True)
            
            # Read with restore_strand=True - should get original sequences
            with open(path, "rb") as fh:
                reader = ZnaReader(fh)
                records = list(reader.records(restore_strand=True))
                
                self.assertEqual(len(records), 4)
                
                # Check all sequences match original
                for i, (r1, r2) in enumerate(test_pairs):
                    rec_r1 = records[i * 2]
                    rec_r2 = records[i * 2 + 1]
                    self.assertEqual(rec_r1[0], r1)
                    self.assertEqual(rec_r2[0], r2)
                    self.assertTrue(rec_r1[2])  # is_r1
                    self.assertTrue(rec_r2[3])  # is_r2
    
    def test_non_strand_specific_unchanged(self):
        """Test that non-strand-specific libraries don't modify sequences."""
        header = ZnaHeader(
            read_group="unstranded",
            strand_specific=False,  # NOT strand specific
        )
        
        test_seq = "AAACCCGGG"
        
        with tempfile.TemporaryDirectory() as tmpdir:
            path = f"{tmpdir}/unstranded.zna"
            
            with open(path, "wb") as fh:
                with ZnaWriter(fh, header) as writer:
                    writer.write_record(test_seq, True, True, False)
            
            with open(path, "rb") as fh:
                reader = ZnaReader(fh)
                self.assertFalse(reader.header.strand_specific)
                self.assertFalse(reader.header.read1_antisense)
                self.assertFalse(reader.header.read2_antisense)
                
                records = list(reader.records(restore_strand=True))
                self.assertEqual(records[0][0], test_seq)


if __name__ == "__main__":
    unittest.main()
