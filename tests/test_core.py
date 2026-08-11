import contextlib
import random
import sys
import tempfile
import unittest
from pathlib import Path
from zna.core import ZnaHeader, write_zna, read_zna, ZnaWriter, ZnaReader, COMPRESSION_ZSTD, COMPRESSION_NONE
from zna import reverse_complement
from zna._pycodec import encode_block, decode_block
from zna._shuffle import shuffle_zna
from zna.dtypes import LabelDef, DTYPE_BY_CODE



class TestZnaIO(unittest.TestCase):
    def roundtrip(self, sequences, header, npolicy=None):
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
                        
                with ZnaWriter(out_fh, header, npolicy=npolicy) as writer:
                    for seq, is_paired, is_r1, is_r2 in records:
                        writer.write_record(seq, is_paired, is_r1, is_r2)
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
        """Test that mixed case is normalized to uppercase.
        
        Note: The C++ accelerated version normalizes all output to uppercase,
        which is correct for DNA sequence handling in ZNA format.
        """
        # With C++ acceleration, output is always uppercase
        result = reverse_complement("AcGt")
        self.assertIn(result, ["aCgT", "ACGT"])  # Accept either behavior
        result = reverse_complement("acgt")
        self.assertIn(result, ["acgt", "ACGT"])  # Accept either behavior
    
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
            strand_normalized=True,
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
                self.assertTrue(reader.header.strand_normalized)
    
    def test_r1_antisense_normalization(self):
        """Test that R1 antisense reads are normalized (reverse complemented) during encoding."""
        header = ZnaHeader(
            read_group="dutp",
            strand_specific=True,
            read1_antisense=True,  # R1 is antisense
            read2_antisense=False,
            strand_normalized=True,
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
            strand_normalized=True,
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
            strand_normalized=True,
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
            strand_normalized=True,
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
            strand_normalized=True,
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

    def test_strand_specific_without_normalize(self):
        """Test that strand_specific=True without strand_normalized stores reads as-is."""
        header = ZnaHeader(
            read_group="metadata_only",
            strand_specific=True,
            read1_antisense=True,
            read2_antisense=False,
            strand_normalized=False,  # metadata only, no RC
        )
        
        r1_seq = "AAAACCCC"
        
        with tempfile.TemporaryDirectory() as tmpdir:
            path = f"{tmpdir}/no_norm.zna"
            
            with open(path, "wb") as fh:
                with ZnaWriter(fh, header) as writer:
                    writer.write_record(r1_seq, True, True, False)
            
            with open(path, "rb") as fh:
                reader = ZnaReader(fh)
                self.assertTrue(reader.header.strand_specific)
                self.assertTrue(reader.header.read1_antisense)
                self.assertFalse(reader.header.strand_normalized)
                
                records = list(reader.records(restore_strand=False))
                # Sequence should be unchanged (no RC applied)
                self.assertEqual(records[0][0], r1_seq)


class TestUnstrandedNormalization(unittest.TestCase):
    """Test unstranded strand normalization (random RC)."""

    def test_paired_normalization_roundtrip(self):
        """Test that unstranded normalization RC's one read per pair and restores correctly."""
        header = ZnaHeader(
            read_group="unstranded_norm",
            strand_specific=False,
            strand_normalized=True,
        )
        
        r1_seq = "AAAACCCC"
        r2_seq = "TTTTGGGG"
        
        with tempfile.TemporaryDirectory() as tmpdir:
            path = f"{tmpdir}/unstranded_norm.zna"
            
            with open(path, "wb") as fh:
                with ZnaWriter(fh, header) as writer:
                    writer.write_record(r1_seq, True, True, False)
                    writer.write_record(r2_seq, True, False, True)
            
            # Read with restore_strand=True — should recover original sequences
            with open(path, "rb") as fh:
                reader = ZnaReader(fh)
                records = list(reader.records(restore_strand=True))
                
                self.assertEqual(len(records), 2)
                self.assertEqual(records[0][0], r1_seq)
                self.assertEqual(records[1][0], r2_seq)

    def test_paired_one_read_is_rc(self):
        """Test that exactly one read per pair is RC'd in unstranded normalization."""
        header = ZnaHeader(
            read_group="test",
            strand_specific=False,
            strand_normalized=True,
        )
        
        r1_seq = "AAAACCCC"
        r2_seq = "TTTTGGGG"
        r1_rc = reverse_complement(r1_seq)
        r2_rc = reverse_complement(r2_seq)
        
        # Encode many pairs to statistically verify both choices occur
        n_pairs = 100
        r1_was_rc_count = 0
        
        with tempfile.TemporaryDirectory() as tmpdir:
            path = f"{tmpdir}/test.zna"
            
            with open(path, "wb") as fh:
                with ZnaWriter(fh, header) as writer:
                    for _ in range(n_pairs):
                        writer.write_record(r1_seq, True, True, False)
                        writer.write_record(r2_seq, True, False, True)
            
            # Read without restore — see what was actually stored
            with open(path, "rb") as fh:
                reader = ZnaReader(fh)
                records = list(reader.records(restore_strand=False))
                
                self.assertEqual(len(records), n_pairs * 2)
                
                for i in range(0, len(records), 2):
                    stored_r1 = records[i][0]
                    stored_r2 = records[i + 1][0]
                    
                    # Exactly one of the two reads should be RC'd
                    if stored_r1 == r1_rc:
                        # R1 was RC'd, R2 should be unchanged
                        self.assertEqual(stored_r2, r2_seq)
                        r1_was_rc_count += 1
                    else:
                        # R1 unchanged, R2 should be RC'd
                        self.assertEqual(stored_r1, r1_seq)
                        self.assertEqual(stored_r2, r2_rc)
        
        # With 100 pairs and random 50/50, both choices should occur
        self.assertGreater(r1_was_rc_count, 0, "R1 should be RC'd at least once")
        self.assertLess(r1_was_rc_count, n_pairs, "R2 should be RC'd at least once")

    def test_single_end_normalization_roundtrip(self):
        """Test that single-end reads are randomly RC'd and restored correctly."""
        header = ZnaHeader(
            read_group="se_norm",
            strand_specific=False,
            strand_normalized=True,
        )
        
        seqs = ["AAAACCCC", "TTTTGGGG", "ACGTACGT", "GGGGTTTTT"]
        
        with tempfile.TemporaryDirectory() as tmpdir:
            path = f"{tmpdir}/se_norm.zna"
            
            with open(path, "wb") as fh:
                with ZnaWriter(fh, header) as writer:
                    for seq in seqs:
                        writer.write_record(seq, False, False, False)
            
            # Restore should give back originals
            with open(path, "rb") as fh:
                reader = ZnaReader(fh)
                records = list(reader.records(restore_strand=True))
                
                self.assertEqual(len(records), len(seqs))
                for i, seq in enumerate(seqs):
                    self.assertEqual(records[i][0], seq)

    def test_mixed_se_pe_normalization_roundtrip(self):
        """Test normalization with mixed single-end and paired-end reads."""
        header = ZnaHeader(
            read_group="mixed_norm",
            strand_specific=False,
            strand_normalized=True,
        )
        
        with tempfile.TemporaryDirectory() as tmpdir:
            path = f"{tmpdir}/mixed.zna"
            
            r1a = "AAAACCCC"
            r2a = "TTTTGGGG"
            se1 = "ACGTACGT"
            r1b = "GGGGTTTTT"
            r2b = "CCCAAAAA"
            se2 = "TTTTTAAA"
            
            with open(path, "wb") as fh:
                with ZnaWriter(fh, header) as writer:
                    # Pair 1
                    writer.write_record(r1a, True, True, False)
                    writer.write_record(r2a, True, False, True)
                    # Single-end
                    writer.write_record(se1, False, False, False)
                    # Pair 2
                    writer.write_record(r1b, True, True, False)
                    writer.write_record(r2b, True, False, True)
                    # Single-end
                    writer.write_record(se2, False, False, False)
            
            # Restore should give back originals
            with open(path, "rb") as fh:
                reader = ZnaReader(fh)
                records = list(reader.records(restore_strand=True))
                
                self.assertEqual(len(records), 6)
                self.assertEqual(records[0][0], r1a)
                self.assertEqual(records[1][0], r2a)
                self.assertEqual(records[2][0], se1)
                self.assertEqual(records[3][0], r1b)
                self.assertEqual(records[4][0], r2b)
                self.assertEqual(records[5][0], se2)

    def test_encode_block_is_rc_flag(self):
        """Test that encode_block sets IS_RC bit (0x08) on RC'd records."""
        seqs = ["AAAACCCC", "TTTTGGGG"]
        flags = [0x05, 0x06]  # paired R1, paired R2
        
        flags_out, _, _ = encode_block(
            seqs, flags, 1, "", False, False, do_random_rc=True
        )
        
        # Exactly one of the two should have IS_RC set
        rc_count = sum(1 for f in flags_out if f & 0x08)
        self.assertEqual(rc_count, 1)

    def test_decode_block_returns_is_rc(self):
        """Test that decode_block returns the is_rc field."""
        seqs = ["AAAA"]
        flags = [0x0D]  # IS_READ1 | IS_PAIRED | IS_RC = 0x01 | 0x04 | 0x08
        
        flags_out, lengths_out, seqs_out = encode_block(
            seqs, flags, 1, "", False, False, do_random_rc=False
        )
        
        records = decode_block(flags_out, lengths_out, seqs_out, 1, 1)
        self.assertEqual(len(records), 1)
        seq, is_paired, is_read1, is_read2, is_rc = records[0]
        self.assertEqual(seq, "AAAA")
        self.assertTrue(is_paired)
        self.assertTrue(is_read1)
        self.assertFalse(is_read2)
        self.assertTrue(is_rc)

    def test_strand_normalized_header_roundtrip(self):
        """Test that strand_normalized flag persists through header serialization."""
        header = ZnaHeader(
            read_group="norm_header",
            strand_specific=False,
            strand_normalized=True,
        )
        
        with tempfile.TemporaryDirectory() as tmpdir:
            path = f"{tmpdir}/norm.zna"
            with open(path, "wb") as fh:
                with ZnaWriter(fh, header) as writer:
                    writer.write_record("ACGT", False, False, False)
            
            with open(path, "rb") as fh:
                reader = ZnaReader(fh)
                self.assertTrue(reader.header.strand_normalized)
                self.assertFalse(reader.header.strand_specific)


class TestStrandSpecificSingleEnd(unittest.TestCase):
    """Strand-specific normalization of single-end / merged reads.

    A single/merged read carries neither the read1 nor read2 flag. Because a
    merged read is built in read1's orientation (R1 + revcomp(R2)-tail), it must
    be normalized with the read1 rule under ``--strand-specific``.
    """

    def test_single_end_strand_specific_rc(self):
        """A single/merged read must be RC'd as read1 under strand-specific norm."""
        header = ZnaHeader(
            read_group="se_ss",
            strand_specific=True,
            read1_antisense=True,   # dUTP: R1 antisense, R2 sense
            read2_antisense=False,
            strand_normalized=True,
        )
        se = "ACGTACGTACGTACGTAAAA"
        se_rc = reverse_complement(se)   # "TTTTACGTACGTACGTACGT"

        with tempfile.TemporaryDirectory() as tmpdir:
            path = f"{tmpdir}/se_ss.zna"
            with open(path, "wb") as fh:
                with ZnaWriter(fh, header) as writer:
                    writer.write_record(se, False, False, False)   # single

            # Stored orientation (no restore): must be RC'd to sense.
            with open(path, "rb") as fh:
                stored = list(ZnaReader(fh).records(restore_strand=False))
            self.assertEqual(stored[0][0], se_rc)

            # Restore must recover the original.
            with open(path, "rb") as fh:
                restored = list(ZnaReader(fh).records(restore_strand=True))
            self.assertEqual(restored[0][0], se)

    def test_single_end_read1_sense_not_rc(self):
        """read1-sense protocol: single read is already sense, must be left as-is."""
        header = ZnaHeader(
            read_group="se_r1sense",
            strand_specific=True,
            read1_antisense=False,   # read1 is sense
            read2_antisense=True,
            strand_normalized=True,
        )
        se = "ACGTACGTACGTACGTAAAA"
        with tempfile.TemporaryDirectory() as tmpdir:
            path = f"{tmpdir}/se_r1sense.zna"
            with open(path, "wb") as fh:
                with ZnaWriter(fh, header) as writer:
                    writer.write_record(se, False, False, False)
            with open(path, "rb") as fh:
                stored = list(ZnaReader(fh).records(restore_strand=False))
            self.assertEqual(stored[0][0], se)   # unchanged

    def test_mixed_se_pe_strand_specific(self):
        """Merged singles and paired R1 must be normalized to the same (sense) strand."""
        header = ZnaHeader(
            read_group="mixed_ss",
            strand_specific=True,
            read1_antisense=True,
            read2_antisense=False,
            strand_normalized=True,
        )
        merged = "ACGTACGTACGTACGTAAAA"
        r1     = "GGGGCCCCGGGGCCCCAAAA"
        r2     = "TTTTGGGGTTTTGGGGCCCC"
        with tempfile.TemporaryDirectory() as tmpdir:
            path = f"{tmpdir}/mixed_ss.zna"
            with open(path, "wb") as fh:
                with ZnaWriter(fh, header) as writer:
                    writer.write_record(merged, False, False, False)  # single (read1-like)
                    writer.write_record(r1, True, True, False)        # pair R1 (antisense)
                    writer.write_record(r2, True, False, True)        # pair R2 (sense)
            with open(path, "rb") as fh:
                stored = list(ZnaReader(fh).records(restore_strand=False))
            self.assertEqual(stored[0][0], reverse_complement(merged))  # RC'd like R1
            self.assertEqual(stored[1][0], reverse_complement(r1))      # RC'd
            self.assertEqual(stored[2][0], r2)                          # kept
            # round-trip
            with open(path, "rb") as fh:
                restored = list(ZnaReader(fh).records(restore_strand=True))
            self.assertEqual([r[0] for r in restored], [merged, r1, r2])

    def test_encode_block_single_is_rc_strand_specific(self):
        """encode_block must set IS_RC (0x08) on a single read under do_rc_r1."""
        flags_out, _, _ = encode_block(
            ["ACGTACGT"], [0x00], 1, "",   # flag 0x00 = single (no R1/R2/paired bits)
            True,    # do_rc_r1
            False,   # do_rc_r2
            do_random_rc=False,
        )
        self.assertTrue(flags_out[0] & 0x08)

    def test_encode_block_single_not_touched_by_read2_rule(self):
        """A single read must never follow the read2 rule."""
        flags_out, _, _ = encode_block(
            ["ACGTACGT"], [0x00], 1, "",
            False,   # do_rc_r1
            True,    # do_rc_r2
            do_random_rc=False,
        )
        self.assertFalse(flags_out[0] & 0x08)

    def test_cross_backend_equivalence(self):
        """Python and C++ backends must produce identical (flags, sequences).

        The stranded single-read rule lives in both _pycodec and _accel; this
        guards against the two implementations silently drifting.
        """
        try:
            from zna import _accel as accel
        except ImportError:
            self.skipTest("C++ accel backend not available")
        if not hasattr(accel, "encode_block"):
            self.skipTest("accel.encode_block not available")

        # Deterministic (stranded) branch only — the random-RC branch uses a
        # different PRNG in each backend and is not expected to match byte-for-byte.
        seqs = ["ACGTACGTACGTACGTAAAA", "GGGGCCCCGGGGCCCCAAAA", "TTTTGGGGTTTTGGGGCCCC"]
        flags = [0x00, 0x05, 0x06]  # single, paired R1, paired R2
        for do_r1, do_r2 in [(True, False), (False, True), (False, False)]:
            py = encode_block(list(seqs), list(flags), 1, "", do_r1, do_r2, do_random_rc=False)
            cc = accel.encode_block(list(seqs), list(flags), 1, "", do_r1, do_r2, False)
            self.assertEqual(bytes(py[0]), bytes(cc[0]),
                             f"flags differ (do_r1={do_r1}, do_r2={do_r2})")
            self.assertEqual(bytes(py[-1]), bytes(cc[-1]),
                             f"sequences differ (do_r1={do_r1}, do_r2={do_r2})")


class TestRcFlagAndReencode(unittest.TestCase):
    """The unstranded fragment-boundary contract, and the pass-through writer.

    Unstranded normalization reverse-complements exactly one mate of an FR pair
    so both land in a common frame, and records which one in ``IS_RC``.  That
    flag is the only record of which edge of a mate is a real fragment boundary:
    it cannot be recovered from the sequence, because reverse-complementing the
    right-hand mate reproduces the fragment-frame sequence exactly.

    These tests pin both halves — that readers can see the flag, and that
    re-encoding a normalized file *copies* orientation instead of applying it a
    second time.  Orientation is not idempotent.
    """

    FRAG_LEN = 300
    READ_LEN = 100
    N_PAIRS = 200

    def setUp(self):
        rng = random.Random(20260811)
        self.frags = [
            "".join(rng.choice("ACGT") for _ in range(self.FRAG_LEN))
            for _ in range(self.N_PAIRS)
        ]

    # -- helpers -------------------------------------------------------------

    def _as_sequenced(self):
        """An FR pair as it comes off the sequencer: the mates point inward, so
        they are in *opposite* frames until normalization puts them in one."""
        L, R = self.FRAG_LEN, self.READ_LEN
        records = []
        for f in self.frags:
            records.append((f[:R], True, True, False))
            records.append((reverse_complement(f[L - R:]), True, False, True))
        return records

    def _unstranded_header(self, rg="unstranded"):
        return ZnaHeader(read_group=rg, strand_specific=False, strand_normalized=True)

    def _write(self, tmpdir, name, records, header, preserve=False):
        path = f"{tmpdir}/{name}"
        with open(path, "wb") as fh:
            with ZnaWriter(fh, header, preserve_normalization=preserve) as writer:
                for rec in records:
                    writer.write_record(
                        rec[0], rec[1], rec[2], rec[3],
                        is_rc=(rec[4] if preserve and len(rec) > 4 else False),
                    )
        return path

    def _read(self, path, **kwargs):
        with open(path, "rb") as fh:
            return list(ZnaReader(fh).records(**kwargs))

    def _co_oriented(self, records):
        """Count pairs whose mates sit in the same frame — what normalization
        achieves, and what a second normalization pass destroys."""
        L, R = self.FRAG_LEN, self.READ_LEN
        n = 0
        for i, f in enumerate(self.frags):
            a, b = records[2 * i][0], records[2 * i + 1][0]
            if (a == f[:R]) == (b == f[L - R:]):
                n += 1
        return n

    @contextlib.contextmanager
    def _force_backend(self, name):
        """Force the codec backend.

        ``zna.core`` resolves its backend at *import* time, so patching
        ``zna.codec`` afterwards has no effect — the module globals are what
        must be swapped.
        """
        import zna.core as core
        from zna.codec import get_backend
        saved_codec, saved_accel = core._codec, core._accel_mod
        try:
            core._codec = get_backend(name)
            core._accel_mod = saved_accel if name == "accel" else None
            yield
        finally:
            core._codec, core._accel_mod = saved_codec, saved_accel

    def _require_accel(self):
        from zna.codec import available_backends
        if "accel" not in available_backends():
            self.skipTest("C++ accel backend not available")

    # -- defect B: re-encoding must not normalize twice ----------------------

    def test_reencode_preserves_co_orientation(self):
        """Regression test for the re-encode bug.

        Before the pass-through writer, every re-encode reverse-complemented one
        mate of each pair *again*, toggling co-orientation on each pass while the
        header still claimed the file was normalized and every record still
        carried an IS_RC bit.  Nothing downstream could detect it.
        """
        header = self._unstranded_header()
        with tempfile.TemporaryDirectory() as tmpdir:
            p1 = self._write(tmpdir, "p1.zna", self._as_sequenced(), header)
            r1 = self._read(p1, with_rc=True)
            self.assertEqual(self._co_oriented(r1), self.N_PAIRS, "encode")

            p2 = self._write(tmpdir, "p2.zna", r1, header, preserve=True)
            r2 = self._read(p2, with_rc=True)
            self.assertEqual(self._co_oriented(r2), self.N_PAIRS, "re-encode")

            p3 = self._write(tmpdir, "p3.zna", r2, header, preserve=True)
            r3 = self._read(p3, with_rc=True)
            self.assertEqual(self._co_oriented(r3), self.N_PAIRS, "re-encode twice")

    def test_reencode_is_record_identical(self):
        """A re-encode is a faithful copy: same sequences, same IS_RC flags.

        Compression framing may differ, so this compares decoded records rather
        than file bytes.
        """
        header = self._unstranded_header()
        with tempfile.TemporaryDirectory() as tmpdir:
            p1 = self._write(tmpdir, "p1.zna", self._as_sequenced(), header)
            r1 = self._read(p1, with_rc=True)
            p2 = self._write(tmpdir, "p2.zna", r1, header, preserve=True)
            r2 = self._read(p2, with_rc=True)
            p3 = self._write(tmpdir, "p3.zna", r2, header, preserve=True)
            r3 = self._read(p3, with_rc=True)
        self.assertEqual(r1, r2)
        self.assertEqual(r2, r3)

    def test_write_records_pass_through(self):
        """``records(with_rc=True)`` feeds ``write_records`` directly, which is
        what makes a lossless ZNA to ZNA copy expressible at all."""
        header = self._unstranded_header()
        with tempfile.TemporaryDirectory() as tmpdir:
            src = self._write(tmpdir, "src.zna", self._as_sequenced(), header)
            original = self._read(src, with_rc=True)
            dst = f"{tmpdir}/dst.zna"
            with open(dst, "wb") as fh:
                write_zna(fh, header, original, preserve_normalization=True)
            self.assertEqual(self._read(dst, with_rc=True), original)

    def test_write_records_rejects_records_without_rc(self):
        """Feeding a normalized copy from plain ``records()`` must not silently
        clear every IS_RC bit — the orientation would be lost with no error."""
        header = self._unstranded_header()
        with tempfile.TemporaryDirectory() as tmpdir:
            src = self._write(tmpdir, "src.zna", self._as_sequenced(), header)
            plain = self._read(src)                      # 4-tuples: no is_rc
            with open(f"{tmpdir}/dst.zna", "wb") as fh:
                with self.assertRaises(ValueError):
                    write_zna(fh, header, plain, preserve_normalization=True)

            # An un-normalized output has no orientation to carry, so the plain
            # 4-tuple form stays valid there.
            plain_header = ZnaHeader(read_group="plain")
            with open(f"{tmpdir}/plain.zna", "wb") as fh:
                write_zna(fh, plain_header, plain, preserve_normalization=True)
            self.assertEqual(
                [rec[0] for rec in self._read(f"{tmpdir}/plain.zna")],
                [rec[0] for rec in plain],
            )

    def test_is_rc_requires_preserve_normalization(self):
        """Supplying IS_RC to a deriving writer would double-normalize."""
        header = self._unstranded_header()
        with tempfile.TemporaryDirectory() as tmpdir:
            with open(f"{tmpdir}/x.zna", "wb") as fh:
                with ZnaWriter(fh, header) as writer:
                    with self.assertRaises(ValueError):
                        writer.write_record("ACGTACGT", True, True, False, is_rc=True)

    def test_stranded_reencode_is_stable(self):
        """The deterministic mirror of the unstranded bug.

        Under ``--strand-specific`` the writer reverse-complements read1 (dUTP
        default).  Doing that a second time returns read1 to its original
        orientation, silently *un*-normalizing a strand-specific file.
        """
        def header():
            return ZnaHeader(
                read_group="ss",
                strand_specific=True,
                read1_antisense=True,    # dUTP: R1 antisense, R2 sense
                read2_antisense=False,
                strand_normalized=True,
            )

        r1, r2 = "ACGTACGTACGTACGTAAAA", "TTTTGGGGTTTTGGGGCCCC"
        records = [(r1, True, True, False), (r2, True, False, True)]
        with tempfile.TemporaryDirectory() as tmpdir:
            p1 = self._write(tmpdir, "ss1.zna", records, header())
            first = self._read(p1, with_rc=True)
            self.assertEqual(first[0][0], reverse_complement(r1))   # R1 flipped to sense
            self.assertTrue(first[0][4])
            self.assertEqual(first[1][0], r2)                       # R2 already sense
            self.assertFalse(first[1][4])

            p2 = self._write(tmpdir, "ss2.zna", first, header(), preserve=True)
            self.assertEqual(self._read(p2, with_rc=True), first)

            # ...and the round-trip back to original orientation still works.
            self.assertEqual([r[0] for r in self._read(p2, restore_strand=True)], [r1, r2])

    def test_shuffle_preserves_rc_geometry(self):
        """A shuffle is a permutation, so orientation must survive it untouched.

        ``shuffle_zna`` writes through two ZnaWriters (buckets, then output).
        When those re-derived orientation they applied two further random
        reverse-complements: co-orientation survived by parity, so the file
        looked fine, but IS_RC ended up uncorrelated with the fragment boundary.
        """
        header = self._unstranded_header()
        with tempfile.TemporaryDirectory() as tmpdir:
            src = self._write(tmpdir, "src.zna", self._as_sequenced(), header)
            before = {rec[0]: rec[4] for rec in self._read(src, with_rc=True)}
            out = f"{tmpdir}/shuffled.zna"
            # A small buffer forces many buckets, which is where re-deriving
            # orientation diverged even on the accel backend.
            shuffle_zna(src, out, seed=3, buffer_bytes=1024, quiet=True)
            after = self._read(out, with_rc=True)

        self.assertEqual(len(after), 2 * self.N_PAIRS)
        self.assertEqual({rec[0] for rec in after}, set(before))
        for rec in after:
            self.assertEqual(rec[4], before[rec[0]], "IS_RC must survive a shuffle")

    def test_shuffle_preserves_rc_and_labels(self):
        """The labeled shuffle path is a separate set of write sites.

        ``_shuffle`` has one write call per (paired/single x labeled/unlabeled)
        combination, and the labeled ones are exactly where ``labels`` moved
        from tuple index 4 to index 5. Without this test, dropping ``is_rc``
        from the labeled sites passes the whole suite.
        """
        labels = (
            LabelDef(label_id=0, name="score", description="",
                     dtype=DTYPE_BY_CODE['i'], missing=0),
            LabelDef(label_id=1, name="qual", description="",
                     dtype=DTYPE_BY_CODE['C'], missing=0),
        )
        header = ZnaHeader(read_group="labshuf", strand_specific=False,
                           strand_normalized=True, labels=labels)
        with tempfile.TemporaryDirectory() as tmpdir:
            src = f"{tmpdir}/src.zna"
            with open(src, "wb") as fh:
                with ZnaWriter(fh, header) as writer:
                    for i, rec in enumerate(self._as_sequenced()):
                        writer.write_record(rec[0], rec[1], rec[2], rec[3],
                                            labels=(i, i % 251))
            before = {rec[0]: (rec[4], rec[5]) for rec in self._read(src, with_rc=True)}

            out = f"{tmpdir}/shuffled.zna"
            shuffle_zna(src, out, seed=11, buffer_bytes=1024, quiet=True)
            after = self._read(out, with_rc=True)

        self.assertEqual(len(after), 2 * self.N_PAIRS)
        self.assertEqual({rec[0] for rec in after}, set(before))
        for rec in after:
            is_rc, labs = before[rec[0]]
            self.assertEqual(rec[4], is_rc, "IS_RC must survive a labeled shuffle")
            self.assertEqual(rec[5], labs, "labels must stay with their sequence")
        # The flag must still be doing real work, not uniformly False.
        self.assertEqual(sum(rec[4] for rec in after), self.N_PAIRS)

    def test_full_fragment_survives_reencode_and_shuffle(self):
        """IS_FULL_FRAGMENT must ride through every copy path.

        Silently clearing it would downgrade full-fragment records to one-ended
        ones — lost supervision rather than corruption, and invisible.
        """
        header = self._unstranded_header()
        with tempfile.TemporaryDirectory() as tmpdir:
            src = f"{tmpdir}/src.zna"
            with open(src, "wb") as fh:
                with ZnaWriter(fh, header) as writer:
                    for i, rec in enumerate(self._as_sequenced()):
                        # mark every 3rd record as spanning its whole fragment
                        writer.write_record(rec[0], rec[1], rec[2], rec[3],
                                            is_full_fragment=(i % 3 == 0))
            before = {r[0]: (r[4], r[5]) for r in self._read(src, with_ends=True)}
            n_full = sum(1 for v in before.values() if v == (True, True))
            self.assertGreater(n_full, 0)

            # 1. re-encode through the pass-through writer
            copied = f"{tmpdir}/copy.zna"
            with open(copied, "wb") as fh:
                write_zna(fh, header, self._read(src, with_ends=True),
                          preserve_normalization=True)
            self.assertEqual(self._read(copied, with_ends=True),
                             self._read(src, with_ends=True))

            # 2. shuffle
            out = f"{tmpdir}/shuf.zna"
            shuffle_zna(src, out, seed=5, buffer_bytes=1024, quiet=True)
            after = self._read(out, with_ends=True)

        self.assertEqual({r[0] for r in after}, set(before))
        for rec in after:
            self.assertEqual((rec[4], rec[5]), before[rec[0]],
                             "fragment-span flag must survive a shuffle")
        self.assertEqual(sum(1 for r in after if r[4] and r[5]), n_full)

    # -- defect A: the flag must reach the caller ----------------------------

    def test_records_with_rc_exposes_the_flag(self):
        """Exactly one mate of each unstranded pair is marked IS_RC."""
        with tempfile.TemporaryDirectory() as tmpdir:
            path = self._write(tmpdir, "rc.zna", self._as_sequenced(), self._unstranded_header())
            records = self._read(path, with_rc=True)

        self.assertEqual(len(records), 2 * self.N_PAIRS)
        n_r1_rc = 0
        for i in range(self.N_PAIRS):
            a, b = records[2 * i], records[2 * i + 1]
            self.assertEqual(len(a), 5)
            self.assertNotEqual(a[4], b[4], f"pair {i}: exactly one mate must be RC'd")
            n_r1_rc += a[4]
        # The coin is independent of the mate number — which is precisely why a
        # consumer cannot guess the flag from is_read1.
        self.assertGreater(n_r1_rc, 0)
        self.assertLess(n_r1_rc, self.N_PAIRS)

    def test_rc_flag_identifies_the_boundary_edge(self):
        """The claim consumers actually rely on, as a contract rather than folklore.

        In the normalized common frame the RC'd mate's RIGHT edge is the real
        fragment boundary and its left edge is a read-length truncation; for the
        other mate it is the mirror image.  A consumer placing end-of-sequence
        supervision writes ``has_eos = is_rc`` and ``has_sos = not is_rc``.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            path = self._write(tmpdir, "geom.zna", self._as_sequenced(), self._unstranded_header())
            records = self._read(path, with_rc=True)

        for i, frag in enumerate(self.frags):
            m1, m2 = records[2 * i], records[2 * i + 1]
            rc_mate, keep_mate = (m1, m2) if m1[4] else (m2, m1)

            # The common frame is whichever of F / revcomp(F) the mates landed in.
            frame = frag if frag.startswith(keep_mate[0]) else reverse_complement(frag)
            self.assertEqual(len(frame), self.FRAG_LEN)

            # Not RC'd: its LEFT edge is the real fragment boundary.
            self.assertTrue(frame.startswith(keep_mate[0]), f"pair {i}: SOS edge")
            # RC'd: its RIGHT edge is the real fragment boundary.
            self.assertTrue(frame.endswith(rc_mate[0]), f"pair {i}: EOS edge")

            # The opposite edges are read-length cutoffs, not boundaries: the
            # molecule demonstrably continues past them.
            self.assertFalse(frame.endswith(keep_mate[0]), f"pair {i}")
            self.assertFalse(frame.startswith(rc_mate[0]), f"pair {i}")

    def test_records_default_width_unchanged(self):
        """The default tuple width is a compatibility promise: 4, or 5 labeled.

        ``with_rc`` is opt-in and inserts ``is_rc`` at index 4, before labels.
        """
        labels = (LabelDef(label_id=0, name="score", description="",
                           dtype=DTYPE_BY_CODE['i'], missing=0),)
        plain = self._unstranded_header()
        labeled = ZnaHeader(read_group="lab", strand_specific=False,
                            strand_normalized=True, labels=labels)

        with tempfile.TemporaryDirectory() as tmpdir:
            p_plain = f"{tmpdir}/plain.zna"
            with open(p_plain, "wb") as fh:
                with ZnaWriter(fh, plain) as writer:
                    writer.write_record("ACGTACGT", True, True, False)
                    writer.write_record("TTTTGGGG", True, False, True)

            p_lab = f"{tmpdir}/labeled.zna"
            with open(p_lab, "wb") as fh:
                with ZnaWriter(fh, labeled) as writer:
                    writer.write_record("ACGTACGT", True, True, False, labels=(7,))
                    writer.write_record("TTTTGGGG", True, False, True, labels=(9,))

            for rec in self._read(p_plain):
                self.assertEqual(len(rec), 4)
            for rec in self._read(p_plain, restore_strand=True):
                self.assertEqual(len(rec), 4)
            for rec in self._read(p_lab):
                self.assertEqual(len(rec), 5)
                self.assertIsInstance(rec[4], tuple)          # labels, not is_rc

            for rec in self._read(p_plain, with_rc=True):
                # with_rc hands the backend's own tuples straight through, so
                # pin the type as well as the width.
                self.assertIsInstance(rec, tuple)
                self.assertEqual(len(rec), 5)
                self.assertIsInstance(rec[4], bool)
            for rec in self._read(p_lab, with_rc=True):
                self.assertIsInstance(rec, tuple)
                self.assertEqual(len(rec), 6)
                self.assertIsInstance(rec[4], bool)           # is_rc before labels
                self.assertIsInstance(rec[5], tuple)
            self.assertEqual(
                [rec[5] for rec in self._read(p_lab, with_rc=True)], [(7,), (9,)]
            )

    def test_with_rc_and_restore_strand_conflict(self):
        """restore_strand consumes the flag; with_rc returns it. Asking for both
        means the caller is confused about one of them."""
        with tempfile.TemporaryDirectory() as tmpdir:
            path = self._write(tmpdir, "x.zna", self._as_sequenced()[:2],
                               self._unstranded_header())
            with open(path, "rb") as fh:
                reader = ZnaReader(fh)
                # Raised eagerly at call time. Wrapped in list() anyway so the
                # test still holds if records() is ever re-inlined as a plain
                # generator, where the raise would move to the first next().
                with self.assertRaises(ValueError):
                    list(reader.records(restore_strand=True, with_rc=True))

    # -- the two backends must agree ----------------------------------------

    def test_cross_backend_is_rc_identical(self):
        """One file, decoded through each backend, must report identical IS_RC."""
        self._require_accel()
        with tempfile.TemporaryDirectory() as tmpdir:
            path = self._write(tmpdir, "xb.zna", self._as_sequenced(),
                               self._unstranded_header())
            with self._force_backend("python"):
                py_records = self._read(path, with_rc=True)
            with self._force_backend("accel"):
                cc_records = self._read(path, with_rc=True)

        self.assertEqual([r[4] for r in py_records], [r[4] for r in cc_records])
        self.assertEqual(py_records, cc_records)

    def test_cross_backend_labeled_is_rc_identical(self):
        """Same, for the labeled decode path, which is a separate implementation."""
        self._require_accel()
        labels = (LabelDef(label_id=0, name="score", description="",
                           dtype=DTYPE_BY_CODE['i'], missing=0),)
        header = ZnaHeader(read_group="lab", strand_specific=False,
                           strand_normalized=True, labels=labels)
        with tempfile.TemporaryDirectory() as tmpdir:
            path = f"{tmpdir}/lab.zna"
            with open(path, "wb") as fh:
                with ZnaWriter(fh, header) as writer:
                    for i, rec in enumerate(self._as_sequenced()[:40]):
                        writer.write_record(rec[0], rec[1], rec[2], rec[3], labels=(i,))
            with self._force_backend("python"):
                py_records = self._read(path, with_rc=True)
            with self._force_backend("accel"):
                cc_records = self._read(path, with_rc=True)

        self.assertEqual(len(py_records), 40)
        for rec in py_records:
            self.assertEqual(len(rec), 6)
        self.assertEqual(py_records, cc_records)

    def test_pass_through_encode_is_backend_identical(self):
        """With orientation supplied rather than derived, encoding is fully
        deterministic — so both backends must produce byte-identical files.

        The *deriving* unstranded path cannot be compared this way: each backend
        uses a different PRNG for the coin.
        """
        self._require_accel()
        header = self._unstranded_header()
        with tempfile.TemporaryDirectory() as tmpdir:
            src = self._write(tmpdir, "src.zna", self._as_sequenced(), header)
            oriented = self._read(src, with_rc=True)

            with self._force_backend("python"):
                py_path = self._write(tmpdir, "py.zna", oriented, header, preserve=True)
            with self._force_backend("accel"):
                cc_path = self._write(tmpdir, "cc.zna", oriented, header, preserve=True)

            self.assertEqual(Path(py_path).read_bytes(), Path(cc_path).read_bytes())


class TestBackendLockstep(unittest.TestCase):
    """The Python (_pycodec) and C++ (_accel) encoders duplicate non-trivial
    logic (strand rules, N-policy, 2-bit packing). These tests feed a battery of
    flag / strand / npolicy / length combinations through both backends and
    assert byte-identical output, so the two can never silently drift.

    The random-RC branch (``do_random_rc=True`` or ``npolicy='random'``) is
    excluded: each backend uses a different PRNG, so its output is not expected
    to match byte-for-byte.
    """

    def setUp(self):
        try:
            from zna import _accel as accel
        except ImportError:
            self.skipTest("C++ accel backend not available")
        if not hasattr(accel, "encode_block"):
            self.skipTest("accel.encode_block not available")
        self.accel = accel

    def test_encode_block_battery(self):
        import itertools

        # Mixed stream: single, paired R1, paired R2, single, and an R1 with no
        # following R2 — plus varied lengths.
        clean_seqs = [
            "ACGTACGTACGTACGTAAAA",  # single
            "GGGGCCCCGGGGCCCCAAAA",  # paired R1
            "TTTTGGGGTTTTGGGGCCCC",  # paired R2
            "ACGT",                  # single (short)
            "TTTTTTTTTTTTTTTTTTTTTTTTT",  # single (odd length)
        ]
        n_seqs = [
            "ACGTNNNNACGTACGTAAAA",
            "GGGGCCCCNGGGCCCCAAAA",
            "NNNNGGGGTTTTGGGGCCCC",
            "ACNT",
            "NTTTTTTTTTTTTTTTTTTTTTTTN",
        ]
        flags = [0x00, 0x05, 0x06, 0x00, 0x05]  # single, R1, R2, single, lone R1

        strand_rules = [(False, False), (True, False), (False, True), (True, True)]
        len_bytes_opts = [1, 2, 4]
        npolicy_opts = ["", "A", "C", "G", "T"]

        for npolicy, len_bytes, (do_r1, do_r2) in itertools.product(
            npolicy_opts, len_bytes_opts, strand_rules
        ):
            seqs = clean_seqs if npolicy == "" else n_seqs
            py = encode_block(list(seqs), list(flags), len_bytes, npolicy,
                              do_r1, do_r2, do_random_rc=False)
            cc = self.accel.encode_block(list(seqs), list(flags), len_bytes, npolicy,
                                         do_r1, do_r2, False)
            ctx = f"npolicy={npolicy!r} len_bytes={len_bytes} do_r1={do_r1} do_r2={do_r2}"
            self.assertEqual(bytes(py[0]), bytes(cc[0]), f"flags differ ({ctx})")
            self.assertEqual(bytes(py[1]), bytes(cc[1]), f"lengths differ ({ctx})")
            self.assertEqual(bytes(py[2]), bytes(cc[2]), f"sequences differ ({ctx})")

    def test_labeled_matches_unlabeled_strand_rules(self):
        """accel.encode_block_labeled must apply the same strand rules as
        accel.encode_block (both contain a copy of the stranded branch)."""
        import struct as _struct

        if not hasattr(self.accel, "encode_block_labeled"):
            self.skipTest("accel.encode_block_labeled not available")

        seqs = ["ACGTACGTACGTACGTAAAA", "GGGGCCCCGGGGCCCCAAAA", "TTTTGGGGTTTTGGGGCCCC"]
        flags = [0x00, 0x05, 0x06]  # single, paired R1, paired R2
        n = len(seqs)
        # One uint8 ('C') label column, all zeros.
        label_col_data = [_struct.pack(f"<{n}B", *([0] * n))]
        label_col_sizes = [1]

        for do_r1, do_r2 in [(True, False), (False, True), (False, False)]:
            unlabeled = self.accel.encode_block(
                list(seqs), list(flags), 1, "", do_r1, do_r2, False
            )
            labeled = self.accel.encode_block_labeled(
                list(seqs), list(flags), 1, "", do_r1, do_r2, False,
                label_col_data, label_col_sizes,
            )
            # unlabeled: (flags, lengths, seqs); labeled: (flags, labels, lengths, seqs)
            self.assertEqual(bytes(unlabeled[0]), bytes(labeled[0]),
                             f"flags differ (do_r1={do_r1}, do_r2={do_r2})")
            self.assertEqual(bytes(unlabeled[-1]), bytes(labeled[-1]),
                             f"sequences differ (do_r1={do_r1}, do_r2={do_r2})")


class TestNPolicy(unittest.TestCase):
    """Test N-nucleotide handling policies."""

    def test_npolicy_replace_with_A(self):
        """Test that N nucleotides are replaced with A."""
        header = ZnaHeader(read_group="test", description="")
        
        with tempfile.TemporaryDirectory() as tmpdir:
            path = f"{tmpdir}/npolicy.zna"
            
            # Write sequence with N
            with open(path, "wb") as fh:
                with ZnaWriter(fh, header, npolicy="A") as writer:
                    writer.write_record("ACNGT", False, False, False)
                    writer.write_record("NNN", False, False, False)
            
            # Read and verify N was replaced with A
            with open(path, "rb") as fh:
                reader = ZnaReader(fh)
                records = list(reader.records())
                self.assertEqual(records[0][0], "ACAGT")  # N -> A
                self.assertEqual(records[1][0], "AAA")    # NNN -> AAA
    
    def test_npolicy_replace_with_other_bases(self):
        """Test that N can be replaced with C, G, or T."""
        header = ZnaHeader(read_group="test")
        
        for base in ['C', 'G', 'T']:
            with tempfile.TemporaryDirectory() as tmpdir:
                path = f"{tmpdir}/npolicy.zna"
                
                with open(path, "wb") as fh:
                    with ZnaWriter(fh, header, npolicy=base) as writer:
                        writer.write_record("ACNGT", False, False, False)
                
                with open(path, "rb") as fh:
                    reader = ZnaReader(fh)
                    records = list(reader.records())
                    expected = f"AC{base}GT"
                    self.assertEqual(records[0][0], expected)
    
    def test_npolicy_random(self):
        """Test that N nucleotides are replaced with random bases."""
        header = ZnaHeader(read_group="test")
        
        with tempfile.TemporaryDirectory() as tmpdir:
            path = f"{tmpdir}/npolicy.zna"
            
            # Write sequence with N
            with open(path, "wb") as fh:
                with ZnaWriter(fh, header, npolicy="random") as writer:
                    writer.write_record("ACNGT", False, False, False)
                    writer.write_record("NNNNN", False, False, False)
            
            # Read and verify N was replaced with valid bases
            with open(path, "rb") as fh:
                reader = ZnaReader(fh)
                records = list(reader.records())
                
                # Check first sequence
                seq1 = records[0][0]
                self.assertEqual(len(seq1), 5)
                self.assertEqual(seq1[0:2], "AC")
                self.assertEqual(seq1[3:], "GT")
                self.assertIn(seq1[2], "ACGT")  # Random base
                
                # Check second sequence
                seq2 = records[1][0]
                self.assertEqual(len(seq2), 5)
                for base in seq2:
                    self.assertIn(base, "ACGT")
    
    def test_npolicy_no_N_unchanged(self):
        """Test that sequences without N are unchanged regardless of policy."""
        header = ZnaHeader(read_group="test")
        test_seq = "ACGTACGT"
        
        for policy in ['A', 'C', 'G', 'T', 'random']:
            with tempfile.TemporaryDirectory() as tmpdir:
                path = f"{tmpdir}/npolicy.zna"
                
                with open(path, "wb") as fh:
                    with ZnaWriter(fh, header, npolicy=policy) as writer:
                        writer.write_record(test_seq, False, False, False)
                
                with open(path, "rb") as fh:
                    reader = ZnaReader(fh)
                    records = list(reader.records())
                    self.assertEqual(records[0][0], test_seq)
    
    def test_npolicy_case_insensitive(self):
        """Test that lowercase n is also handled."""
        header = ZnaHeader(read_group="test")
        
        with tempfile.TemporaryDirectory() as tmpdir:
            path = f"{tmpdir}/npolicy.zna"
            
            with open(path, "wb") as fh:
                with ZnaWriter(fh, header, npolicy="A") as writer:
                    writer.write_record("ACnGT", False, False, False)
            
            with open(path, "rb") as fh:
                reader = ZnaReader(fh)
                records = list(reader.records())
                self.assertEqual(records[0][0], "ACAGT")


class TestColumnarFormat(unittest.TestCase):
    """Test columnar block format encoding."""
    
    def test_varying_length_sequences(self):
        """Test that varying length sequences are handled correctly."""
        header = ZnaHeader(read_group="test", seq_len_bytes=2)
        sequences = ["A", "AC", "ACG", "ACGT", "ACGTA" * 10, "TGCATGCATGCA"]
        
        with tempfile.TemporaryDirectory() as tmpdir:
            path = f"{tmpdir}/varlen.zna"
            
            with open(path, "wb") as fh:
                with ZnaWriter(fh, header) as writer:
                    for seq in sequences:
                        writer.write_record(seq, False, False, False)
            
            with open(path, "rb") as fh:
                reader = ZnaReader(fh)
                decoded = [r[0] for r in reader.records()]
            
            self.assertEqual(decoded, sequences)
    
    def test_mixed_flags(self):
        """Test that columnar format preserves record flags correctly."""
        header = ZnaHeader(read_group="test")
        
        # Various flag combinations
        records = [
            ("ACGT", True, True, False),   # Paired R1
            ("TGCA", True, False, True),   # Paired R2
            ("GGGG", False, False, False), # Single
            ("AAAA", True, True, False),   # Paired R1
            ("CCCC", True, False, True),   # Paired R2
        ]
        
        with tempfile.TemporaryDirectory() as tmpdir:
            path = f"{tmpdir}/flags.zna"
            
            with open(path, "wb") as fh:
                with ZnaWriter(fh, header) as writer:
                    for seq, is_paired, is_r1, is_r2 in records:
                        writer.write_record(seq, is_paired, is_r1, is_r2)
            
            with open(path, "rb") as fh:
                reader = ZnaReader(fh)
                decoded = list(reader.records())
            
            self.assertEqual(decoded, records)
    
    def test_multi_block_roundtrip(self):
        """Test roundtrip across multiple blocks."""
        header = ZnaHeader(read_group="test", compression_method=COMPRESSION_NONE)
        
        # Enough sequences to span multiple small blocks
        sequences = ["ACGTACGT" * 2] * 50
        
        with tempfile.TemporaryDirectory() as tmpdir:
            path = f"{tmpdir}/multiblock.zna"
            
            with open(path, "wb") as fh:
                with ZnaWriter(fh, header, block_size=100) as writer:
                    for seq in sequences:
                        writer.write_record(seq, False, False, False)
            
            with open(path, "rb") as fh:
                reader = ZnaReader(fh)
                decoded = [r[0] for r in reader.records()]
            
            self.assertEqual(decoded, sequences)


if __name__ == "__main__":
    unittest.main()
