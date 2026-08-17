"""
Unit tests for CLI functionality (encode, decode, inspect).
"""
import gzip
import io
import random
import tempfile
import struct
from pathlib import Path
from io import BytesIO

import pytest

from zna import reverse_complement
from zna.cli import (
    parse_fasta, parse_fastq, choose_parser,
    stream_inputs, encode_command, decode_command, inspect_command,
    get_base_name, get_read_suffix_number, parse_fastq_with_names,
    parse_block_size, shuffle_command,
)
from zna._shuffle import shuffle_zna
from zna.core import (
    ZnaHeader, ZnaWriter, ZnaReader,
    COMPRESSION_ZSTD, COMPRESSION_NONE,
    _BLOCK_HEADER_FMT, _FILE_HEADER_SIZE, _BLOCK_HEADER_SIZE,
)


# --- Test Data ---

FASTA_DATA = b""">seq1
ACGTACGT
>seq2
TGCATGCA
>seq3 with description
AAAA
CCCC
GGGG
"""

FASTQ_DATA = b"""@read1
ACGTACGT
+
IIIIIIII
@read2
TGCATGCA
+
IIIIIIII
@read3
AAAACCCCGGGG
+
IIIIIIIIIIII
"""

INTERLEAVED_FASTQ = b"""@read1/1
ACGTACGT
+
IIIIIIII
@read1/2
TGCATGCA
+
IIIIIIII
@read2/1
AAAA
+
IIII
@read2/2
CCCC
+
IIII
"""


# --- Parser Tests ---

class TestParsers:
    def test_parse_fasta_simple(self):
        """Test FASTA parser with simple sequences."""
        fh = BytesIO(FASTA_DATA)
        sequences = list(parse_fasta(fh))
        assert len(sequences) == 3
        assert sequences[0] == "ACGTACGT"
        assert sequences[1] == "TGCATGCA"
        assert sequences[2] == "AAAACCCCGGGG"

    def test_parse_fasta_empty(self):
        """Test FASTA parser with empty input."""
        fh = BytesIO(b"")
        sequences = list(parse_fasta(fh))
        assert len(sequences) == 0

    def test_parse_fasta_single_line(self):
        """Test FASTA with single sequence."""
        fh = BytesIO(b">seq1\nACGT\n")
        sequences = list(parse_fasta(fh))
        assert sequences == ["ACGT"]

    def test_parse_fastq_simple(self):
        """Test FASTQ parser with simple sequences."""
        fh = BytesIO(FASTQ_DATA)
        sequences = list(parse_fastq(fh))
        assert len(sequences) == 3
        assert sequences[0] == "ACGTACGT"
        assert sequences[1] == "TGCATGCA"
        assert sequences[2] == "AAAACCCCGGGG"

    def test_parse_fastq_empty(self):
        """Test FASTQ parser with empty input."""
        fh = BytesIO(b"")
        sequences = list(parse_fastq(fh))
        assert len(sequences) == 0
    
    def test_parse_fastq_with_names(self):
        """Test FASTQ parser that returns names and sequences."""
        fh = BytesIO(FASTQ_DATA)
        entries = list(parse_fastq_with_names(fh))
        assert len(entries) == 3
        assert entries[0] == ("read1", "ACGTACGT")
        assert entries[1] == ("read2", "TGCATGCA")
        assert entries[2] == ("read3", "AAAACCCCGGGG")

    def test_get_parser_fasta(self):
        """Test parser selection for FASTA files."""
        assert choose_parser("test.fasta") == parse_fasta
        assert choose_parser("test.fa") == parse_fasta
        assert choose_parser("test.fna") == parse_fasta
        assert choose_parser("test.fa.gz") == parse_fasta
        assert choose_parser("test.fasta.gz") == parse_fasta
        assert choose_parser("test.fna.gz") == parse_fasta

    def test_get_parser_fastq_default(self):
        """Test parser defaults to FASTQ."""
        assert choose_parser("test.fastq") == parse_fastq
        assert choose_parser("test.fq") == parse_fastq
        assert choose_parser("test.fastq.gz") == parse_fastq
        assert choose_parser("test.fq.gz") == parse_fastq
        # Unknown extensions default to FASTQ with warning
        assert choose_parser("test.txt") == parse_fastq
        assert choose_parser(None) == parse_fastq
    
    def test_get_parser_format_override(self):
        """Test format override with --fasta/--fastq flags."""
        # Override to FASTA even for .fastq file
        assert choose_parser("test.fastq", format_override='fasta') == parse_fasta
        # Override to FASTQ even for .fasta file
        assert choose_parser("test.fasta", format_override='fastq') == parse_fastq
        # Override works with stdin
        assert choose_parser(None, format_override='fasta') == parse_fasta
        assert choose_parser(None, format_override='fastq') == parse_fastq


# --- Read Name Helper Tests ---

class TestReadNameHelpers:
    def test_get_base_name_with_suffix(self):
        """Test extracting base name from reads with /1 or /2 suffix."""
        assert get_base_name("read1/1") == "read1"
        assert get_base_name("read1/2") == "read1"
        assert get_base_name("INSTRUMENT:123:FLOWCELL:1:1:1234:5678/1") == "INSTRUMENT:123:FLOWCELL:1:1:1234:5678"
    
    def test_get_base_name_without_suffix(self):
        """Test extracting base name from reads without suffix."""
        assert get_base_name("read1") == "read1"
        assert get_base_name("single_read") == "single_read"
    
    def test_get_base_name_with_comment(self):
        """Test extracting base name ignores comments."""
        assert get_base_name("read1/1 comment text") == "read1"
        assert get_base_name("read1/2 merged") == "read1"
        assert get_base_name("read1 merged comment") == "read1"
    
    def test_get_read_suffix_number(self):
        """Test extracting read number from suffix."""
        assert get_read_suffix_number("read1/1") == 1
        assert get_read_suffix_number("read1/2") == 2
        assert get_read_suffix_number("read1") == 0
        assert get_read_suffix_number("single") == 0
    
    def test_get_read_suffix_number_with_comment(self):
        """Test suffix number extraction ignores comments."""
        assert get_read_suffix_number("read1/1 comment") == 1
        assert get_read_suffix_number("read1/2 merged") == 2
        assert get_read_suffix_number("read1 comment") == 0


class TestParseBlockSize:
    def test_plain_integer(self):
        assert parse_block_size("524288") == 524288

    def test_kilobytes(self):
        assert parse_block_size("512K") == 512 * 1024
        assert parse_block_size("512KB") == 512 * 1024
        assert parse_block_size("512k") == 512 * 1024

    def test_megabytes(self):
        assert parse_block_size("4M") == 4 * 1024 * 1024
        assert parse_block_size("4MB") == 4 * 1024 * 1024
        assert parse_block_size("4m") == 4 * 1024 * 1024

    def test_whitespace(self):
        assert parse_block_size("  4M  ") == 4 * 1024 * 1024

    def test_invalid_raises(self):
        import argparse
        with pytest.raises(argparse.ArgumentTypeError):
            parse_block_size("abc")
        with pytest.raises(argparse.ArgumentTypeError):
            parse_block_size("4X")


# --- Encode Tests ---

class TestEncode:
    def test_encode_fasta_uncompressed(self):
        """Test encoding FASTA to uncompressed ZNA."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create input FASTA
            fasta_path = f"{tmpdir}/input.fasta"
            with open(fasta_path, "wb") as f:
                f.write(FASTA_DATA)

            # Encode
            zna_path = f"{tmpdir}/output.zna"
            
            class Args:
                files = [fasta_path]
                interleaved = False
                fasta = False
                fastq = False
                read_group = "test_rg"
                description = "test_desc"
                strand_specific = False
                output = zna_path
                seq_len_bytes = 1
                block_size = 131072
                compress_flag = False  # Uncompressed
                level = 3

            encode_command(Args())

            # Verify output
            assert Path(zna_path).exists()
            with open(zna_path, "rb") as f:
                reader = ZnaReader(f)
                assert reader.header.read_group == "test_rg"
                assert reader.header.description == "test_desc"
                assert reader.header.compression_method == COMPRESSION_NONE
                
                records = list(reader.records())
                assert len(records) == 3
                assert records[0][0] == "ACGTACGT"
                assert records[1][0] == "TGCATGCA"
                assert records[2][0] == "AAAACCCCGGGG"

    def test_encode_fastq_compressed(self):
        """Test encoding FASTQ to compressed ZNA."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create input FASTQ
            fastq_path = f"{tmpdir}/input.fastq"
            with open(fastq_path, "wb") as f:
                f.write(FASTQ_DATA)

            # Encode with compression
            zna_path = f"{tmpdir}/output.zzna"
            
            class Args:
                files = [fastq_path]
                interleaved = False
                fasta = False
                fastq = False
                read_group = "compressed_test"
                description = ""
                strand_specific = False
                output = zna_path
                seq_len_bytes = 1
                block_size = 131072
                compress_flag = True  # Compressed
                level = 3

            encode_command(Args())

            # Verify output
            with open(zna_path, "rb") as f:
                reader = ZnaReader(f)
                assert reader.header.compression_method == COMPRESSION_ZSTD
                assert reader.header.compression_level == 3
                
                records = list(reader.records())
                assert len(records) == 3

    def test_encode_paired_end_files(self):
        """Test encoding paired-end reads from separate files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create R1 and R2 files
            r1_path = f"{tmpdir}/R1.fastq"
            r2_path = f"{tmpdir}/R2.fastq"
            
            r1_data = b"@read1\nACGT\n+\nIIII\n@read2\nTGCA\n+\nIIII\n"
            r2_data = b"@read1\nAAAA\n+\nIIII\n@read2\nCCCC\n+\nIIII\n"
            
            with open(r1_path, "wb") as f:
                f.write(r1_data)
            with open(r2_path, "wb") as f:
                f.write(r2_data)

            # Encode
            zna_path = f"{tmpdir}/paired.zna"
            
            class Args:
                files = [r1_path, r2_path]
                interleaved = False
                fasta = False
                fastq = False
                read_group = "paired"
                description = ""
                strand_specific = False
                output = zna_path
                seq_len_bytes = 1
                block_size = 131072
                compress_flag = False
                level = 3

            encode_command(Args())

            # Verify paired-end structure
            with open(zna_path, "rb") as f:
                reader = ZnaReader(f)
                records = list(reader.records())
                
                assert len(records) == 4
                # First pair
                assert records[0] == ("ACGT", True, True, False)  # R1
                assert records[1] == ("AAAA", True, False, True)  # R2
                # Second pair
                assert records[2] == ("TGCA", True, True, False)  # R1
                assert records[3] == ("CCCC", True, False, True)  # R2

    def test_encode_interleaved_input(self):
        """Test encoding interleaved paired-end input."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create interleaved file
            interleaved_path = f"{tmpdir}/interleaved.fastq"
            with open(interleaved_path, "wb") as f:
                f.write(INTERLEAVED_FASTQ)

            # Encode
            zna_path = f"{tmpdir}/interleaved.zna"
            
            class Args:
                files = [interleaved_path]
                interleaved = True
                fasta = False
                fastq = False
                read_group = "interleaved"
                description = ""
                strand_specific = False
                output = zna_path
                seq_len_bytes = 1
                block_size = 131072
                compress_flag = False
                level = 3

            encode_command(Args())

            # Verify
            with open(zna_path, "rb") as f:
                reader = ZnaReader(f)
                records = list(reader.records())
                
                assert len(records) == 4
                assert records[0] == ("ACGTACGT", True, True, False)
                assert records[1] == ("TGCATGCA", True, False, True)
                assert records[2] == ("AAAA", True, True, False)
                assert records[3] == ("CCCC", True, False, True)

    def test_encode_gzipped_input(self):
        """Test encoding from gzipped input."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create gzipped FASTQ
            fastq_gz_path = f"{tmpdir}/input.fastq.gz"
            with gzip.open(fastq_gz_path, "wb") as f:
                f.write(FASTQ_DATA)

            # Encode
            zna_path = f"{tmpdir}/output.zna"
            
            class Args:
                files = [fastq_gz_path]
                interleaved = False
                fasta = False
                fastq = False
                read_group = "gzipped"
                description = ""
                strand_specific = False
                output = zna_path
                seq_len_bytes = 1
                block_size = 131072
                compress_flag = False
                level = 3

            encode_command(Args())

            # Verify
            with open(zna_path, "rb") as f:
                reader = ZnaReader(f)
                records = list(reader.records())
                assert len(records) == 3


# --- Decode Tests ---

class TestDecode:
    def test_decode_roundtrip_uncompressed(self):
        """Test roundtrip: encode then decode, verify sequences match."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create test sequences
            test_seqs = ["ACGTACGT", "TGCATGCA", "AAAACCCCGGGG"]
            
            # Encode
            zna_path = f"{tmpdir}/test.zna"
            header = ZnaHeader(
                read_group="roundtrip",
                compression_method=COMPRESSION_NONE
            )
            with open(zna_path, "wb") as f:
                with ZnaWriter(f, header) as writer:
                    for seq in test_seqs:
                        writer.write_record(seq, False, False, False)

            # Decode
            fasta_path = f"{tmpdir}/output.fasta"
            
            class Args:
                input = zna_path
                output = fasta_path
                quiet = True
                gzip = False

            decode_command(Args())

            # Verify sequences
            with open(fasta_path, "r") as f:
                lines = f.readlines()
            
            # Extract sequences (skip headers)
            decoded_seqs = [lines[i].strip() for i in range(1, len(lines), 2)]
            assert decoded_seqs == test_seqs

    def test_decode_roundtrip_compressed(self):
        """Test roundtrip with compressed ZNA."""
        with tempfile.TemporaryDirectory() as tmpdir:
            test_seqs = ["ACGT" * 100, "TGCA" * 100, "AAAA" * 100]
            
            # Encode with compression (use seq_len_bytes=2 for longer sequences)
            zna_path = f"{tmpdir}/test.zzna"
            header = ZnaHeader(
                read_group="compressed_roundtrip",
                seq_len_bytes=2,  # Support up to 65535 bp
                compression_method=COMPRESSION_ZSTD,
                compression_level=5
            )
            with open(zna_path, "wb") as f:
                with ZnaWriter(f, header) as writer:
                    for seq in test_seqs:
                        writer.write_record(seq, False, False, False)

            # Decode
            fasta_path = f"{tmpdir}/output.fasta"
            
            class Args:
                input = zna_path
                output = fasta_path
                quiet = True
                gzip = False

            decode_command(Args())

            # Verify
            with open(fasta_path, "r") as f:
                lines = f.readlines()
            
            decoded_seqs = [lines[i].strip() for i in range(1, len(lines), 2)]
            assert decoded_seqs == test_seqs

    def test_decode_paired_end_split(self):
        """Test decoding paired-end to split files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create paired-end ZNA
            zna_path = f"{tmpdir}/paired.zna"
            header = ZnaHeader(read_group="split_test", compression_method=COMPRESSION_NONE)
            
            with open(zna_path, "wb") as f:
                with ZnaWriter(f, header) as writer:
                    writer.write_record("ACGT", True, True, False)  # R1
                    writer.write_record("AAAA", True, False, True)  # R2
                    writer.write_record("TGCA", True, True, False)  # R1
                    writer.write_record("CCCC", True, False, True)  # R2

            # Decode to split files
            output_pattern = f"{tmpdir}/out#.fasta"
            
            class Args:
                input = zna_path
                output = output_pattern
                quiet = True
                gzip = False

            decode_command(Args())

            # Verify R1 file
            r1_path = f"{tmpdir}/out_1.fasta"
            with open(r1_path, "r") as f:
                lines = f.readlines()
            r1_seqs = [lines[i].strip() for i in range(1, len(lines), 2)]
            assert r1_seqs == ["ACGT", "TGCA"]

            # Verify R2 file
            r2_path = f"{tmpdir}/out_2.fasta"
            with open(r2_path, "r") as f:
                lines = f.readlines()
            r2_seqs = [lines[i].strip() for i in range(1, len(lines), 2)]
            assert r2_seqs == ["AAAA", "CCCC"]

    def test_decode_gzipped_output(self):
        """Test decoding to gzipped output."""
        with tempfile.TemporaryDirectory() as tmpdir:
            test_seqs = ["ACGT", "TGCA"]
            
            # Create ZNA
            zna_path = f"{tmpdir}/test.zna"
            header = ZnaHeader(read_group="gzip_out", compression_method=COMPRESSION_NONE)
            with open(zna_path, "wb") as f:
                with ZnaWriter(f, header) as writer:
                    for seq in test_seqs:
                        writer.write_record(seq, False, False, False)

            # Decode to gzipped output
            fasta_gz_path = f"{tmpdir}/output.fasta.gz"
            
            class Args:
                input = zna_path
                output = fasta_gz_path
                quiet = True
                gzip = False  # Inferred from filename

            decode_command(Args())

            # Verify gzipped output
            assert Path(fasta_gz_path).exists()
            with gzip.open(fasta_gz_path, "rt") as f:
                lines = f.readlines()
            
            decoded_seqs = [lines[i].strip() for i in range(1, len(lines), 2)]
            assert decoded_seqs == test_seqs

    def test_decode_record_count_consistency(self):
        """Test that record counts are preserved through encode/decode."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create large dataset
            n_records = 1000
            test_seqs = [f"ACGT{'A' * (i % 100)}" for i in range(n_records)]
            
            # Encode
            zna_path = f"{tmpdir}/large.zna"
            header = ZnaHeader(read_group="count_test", compression_method=COMPRESSION_NONE)
            with open(zna_path, "wb") as f:
                with ZnaWriter(f, header) as writer:
                    for seq in test_seqs:
                        writer.write_record(seq, False, False, False)

            # Decode
            fasta_path = f"{tmpdir}/output.fasta"
            
            class Args:
                input = zna_path
                output = fasta_path
                quiet = True
                gzip = False

            decode_command(Args())

            # Count records
            with open(fasta_path, "r") as f:
                lines = [l for l in f if l.startswith(">")]
            
            assert len(lines) == n_records


# --- Inspect Tests ---

class TestInspect:
    def test_inspect_uncompressed(self):
        """Test inspect on uncompressed ZNA file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create test file
            zna_path = f"{tmpdir}/test.zna"
            header = ZnaHeader(
                read_group="inspect_test",
                description="test description",
                seq_len_bytes=2,
                strand_specific=True,
                compression_method=COMPRESSION_NONE
            )
            
            test_seqs = ["ACGT" * 10] * 50  # 50 sequences
            
            with open(zna_path, "wb") as f:
                with ZnaWriter(f, header, block_size=500) as writer:
                    for seq in test_seqs:
                        writer.write_record(seq, False, False, False)

            # Inspect
            class Args:
                input = zna_path

            # Capture would require redirecting stdout, so just verify it runs
            inspect_command(Args())

    def test_inspect_compressed(self):
        """Test inspect on compressed ZNA file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            zna_path = f"{tmpdir}/compressed.zzna"
            header = ZnaHeader(
                read_group="compressed_inspect",
                compression_method=COMPRESSION_ZSTD,
                compression_level=5
            )
            
            test_seqs = ["ACGTACGT"] * 100
            
            with open(zna_path, "wb") as f:
                with ZnaWriter(f, header) as writer:
                    for seq in test_seqs:
                        writer.write_record(seq, False, False, False)

            # Inspect
            class Args:
                input = zna_path

            inspect_command(Args())

    def test_inspect_block_statistics(self):
        """Test that inspect correctly counts blocks and records."""
        with tempfile.TemporaryDirectory() as tmpdir:
            zna_path = f"{tmpdir}/blocks.zna"
            header = ZnaHeader(read_group="block_test", compression_method=COMPRESSION_NONE)
            
            # Force multiple blocks with small block size
            n_records = 100
            test_seqs = ["ACGT"] * n_records
            
            with open(zna_path, "wb") as f:
                with ZnaWriter(f, header, block_size=50) as writer:  # Small block
                    for seq in test_seqs:
                        writer.write_record(seq, False, False, False)

            # Manually verify block count
            with open(zna_path, "rb") as f:
                reader = ZnaReader(f)
                h = reader.header
                
                # Skip to blocks
                f.seek(_FILE_HEADER_SIZE + len(h.read_group) + len(h.description))
                
                block_count = 0
                total_records = 0
                
                while True:
                    b_header = f.read(_BLOCK_HEADER_SIZE)  # Block header (20 bytes)
                    if not b_header:
                        break
                    
                    c_size, u_size, n_recs, flags_size, lengths_size = struct.unpack(
                        _BLOCK_HEADER_FMT, b_header
                    )
                    block_count += 1
                    total_records += n_recs
                    f.seek(c_size, 1)  # Skip payload
                
                assert total_records == n_records
                assert block_count > 1  # Should have multiple blocks

            # Run inspect command
            class Args:
                input = zna_path

            inspect_command(Args())


# --- Integration Tests ---

class TestIntegration:
    def test_full_workflow_single_end(self):
        """Test complete workflow: FASTQ -> ZNA -> FASTA."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create input FASTQ
            fastq_path = f"{tmpdir}/input.fastq"
            with open(fastq_path, "wb") as f:
                f.write(FASTQ_DATA)

            # Encode
            zna_path = f"{tmpdir}/encoded.zzna"
            
            class EncArgs:
                files = [fastq_path]
                interleaved = False
                fasta = False
                fastq = False
                read_group = "workflow_test"
                description = "full workflow"
                strand_specific = False
                output = zna_path
                seq_len_bytes = 1
                block_size = 131072
                compress_flag = True
                level = 3

            encode_command(EncArgs())

            # Inspect
            class InspArgs:
                input = zna_path

            inspect_command(InspArgs())

            # Decode
            fasta_path = f"{tmpdir}/decoded.fasta"
            
            class DecArgs:
                input = zna_path
                output = fasta_path
                quiet = True
                gzip = False

            decode_command(DecArgs())

            # Verify roundtrip
            original_seqs = []
            fh = BytesIO(FASTQ_DATA)
            for seq in parse_fastq(fh):
                original_seqs.append(seq)

            with open(fasta_path, "r") as f:
                lines = f.readlines()
            decoded_seqs = [lines[i].strip() for i in range(1, len(lines), 2)]

            assert decoded_seqs == original_seqs

    def test_full_workflow_paired_end(self):
        """Test complete workflow with paired-end data."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create paired files
            r1_path = f"{tmpdir}/R1.fastq.gz"
            r2_path = f"{tmpdir}/R2.fastq.gz"
            
            r1_data = b"@r1\nACGT\n+\nIIII\n@r2\nTGCA\n+\nIIII\n"
            r2_data = b"@r1\nAAAA\n+\nIIII\n@r2\nCCCC\n+\nIIII\n"
            
            with gzip.open(r1_path, "wb") as f:
                f.write(r1_data)
            with gzip.open(r2_path, "wb") as f:
                f.write(r2_data)

            # Encode
            zna_path = f"{tmpdir}/paired.zzna"
            
            class EncArgs:
                files = [r1_path, r2_path]
                interleaved = False
                fasta = False
                fastq = False
                read_group = "paired_workflow"
                description = ""
                strand_specific = False
                output = zna_path
                seq_len_bytes = 1
                block_size = 131072
                compress_flag = True
                level = 5

            encode_command(EncArgs())

            # Decode to split gzipped files
            output_pattern = f"{tmpdir}/decoded#.fasta.gz"
            
            class DecArgs:
                input = zna_path
                output = output_pattern
                quiet = True
                gzip = False

            decode_command(DecArgs())

            # Verify R1
            with gzip.open(f"{tmpdir}/decoded_1.fasta.gz", "rt") as f:
                lines = f.readlines()
            r1_decoded = [lines[i].strip() for i in range(1, len(lines), 2)]
            assert r1_decoded == ["ACGT", "TGCA"]

            # Verify R2
            with gzip.open(f"{tmpdir}/decoded_2.fasta.gz", "rt") as f:
                lines = f.readlines()
            r2_decoded = [lines[i].strip() for i in range(1, len(lines), 2)]
            assert r2_decoded == ["AAAA", "CCCC"]


# --- Strand-Specific CLI Tests ---

class TestStrandProtocol:
    """Test CLI strand-specific functionality."""
    
    def test_encode_with_strand_specific_defaults(self):
        """Test encoding with --strand-specific uses dUTP defaults (R1 antisense, R2 sense)."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create paired FASTQ files
            r1_path = f"{tmpdir}/R1.fastq"
            r2_path = f"{tmpdir}/R2.fastq"
            
            # R1: AAAACCCC, R2: TTTTGGGG
            r1_data = b"@read1\nAAAACCCC\n+\nIIIIIIII\n"
            r2_data = b"@read1\nTTTTGGGG\n+\nIIIIIIII\n"
            
            with open(r1_path, "wb") as f:
                f.write(r1_data)
            with open(r2_path, "wb") as f:
                f.write(r2_data)
            
            # Encode with strand-specific (uses dUTP defaults)
            zna_path = f"{tmpdir}/strand.zna"
            
            class Args:
                files = [r1_path, r2_path]
                interleaved = False
                fasta = False
                fastq = False
                read_group = "dutp_test"
                description = ""
                strand_specific = True  # Enable strand-specific
                strand_normalize = True  # Enable normalization
                read1_sense = False     # Default: R1 is antisense
                read2_antisense = False # Default: R2 is sense
                output = zna_path
                seq_len_bytes = 1
                block_size = 131072
                compress_flag = False
                level = 3
            
            encode_command(Args())
            
            # Read and verify header
            with open(zna_path, "rb") as f:
                reader = ZnaReader(f)
                assert reader.header.strand_specific == True
                assert reader.header.read1_antisense == True
                assert reader.header.read2_antisense == False
                assert reader.header.strand_normalized == True
    
    def test_decode_restore_strand(self):
        """Test decoding with --restore-strand flag."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create paired FASTQ files
            r1_path = f"{tmpdir}/R1.fastq"
            r2_path = f"{tmpdir}/R2.fastq"
            
            r1_original = "AAAACCCC"
            r2_original = "TTTTGGGG"
            
            r1_data = f"@read1\n{r1_original}\n+\nIIIIIIII\n".encode()
            r2_data = f"@read1\n{r2_original}\n+\nIIIIIIII\n".encode()
            
            with open(r1_path, "wb") as f:
                f.write(r1_data)
            with open(r2_path, "wb") as f:
                f.write(r2_data)
            
            # Encode with strand-specific (dUTP defaults)
            zna_path = f"{tmpdir}/strand.zna"
            
            class EncArgs:
                files = [r1_path, r2_path]
                interleaved = False
                fasta = False
                fastq = False
                read_group = "dutp_test"
                description = ""
                strand_specific = True
                strand_normalize = True
                read1_sense = False      # R1 is antisense
                read2_antisense = False  # R2 is sense
                output = zna_path
                seq_len_bytes = 1
                block_size = 131072
                compress_flag = False
                level = 3
            
            encode_command(EncArgs())
            
            # Decode WITH restore_strand=True (should get original sequences)
            output_restored = f"{tmpdir}/restored.fasta"
            
            class DecArgsRestored:
                input = zna_path
                output = output_restored
                quiet = True
                gzip = False
                restore_strand = True
            
            decode_command(DecArgsRestored())
            
            # Verify sequences are restored to original
            with open(output_restored, "r") as f:
                lines = f.readlines()
            
            seqs = [lines[i].strip() for i in range(1, len(lines), 2)]
            assert seqs[0] == r1_original  # R1 restored
            assert seqs[1] == r2_original  # R2 unchanged
    
    def test_decode_without_restore_strand(self):
        """Test decoding without --restore-strand gives sense-normalized sequences."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create paired FASTQ files
            r1_path = f"{tmpdir}/R1.fastq"
            r2_path = f"{tmpdir}/R2.fastq"
            
            r1_original = "AAAACCCC"
            r1_revcomp = "GGGGTTTT"  # This is what should be stored
            r2_original = "TTTTGGGG"
            
            r1_data = f"@read1\n{r1_original}\n+\nIIIIIIII\n".encode()
            r2_data = f"@read1\n{r2_original}\n+\nIIIIIIII\n".encode()
            
            with open(r1_path, "wb") as f:
                f.write(r1_data)
            with open(r2_path, "wb") as f:
                f.write(r2_data)
            
            # Encode with strand-specific (dUTP defaults)
            zna_path = f"{tmpdir}/strand.zna"
            
            class EncArgs:
                files = [r1_path, r2_path]
                interleaved = False
                fasta = False
                fastq = False
                read_group = "dutp_test"
                description = ""
                strand_specific = True
                strand_normalize = True
                read1_sense = False      # R1 is antisense
                read2_antisense = False  # R2 is sense
                output = zna_path
                seq_len_bytes = 1
                block_size = 131072
                compress_flag = False
                level = 3
            
            encode_command(EncArgs())
            
            # Decode WITHOUT restore_strand (default) - should get normalized sequences
            output_normalized = f"{tmpdir}/normalized.fasta"
            
            class DecArgsNormal:
                input = zna_path
                output = output_normalized
                quiet = True
                gzip = False
                restore_strand = False
            
            decode_command(DecArgsNormal())
            
            # Verify sequences are sense-normalized
            with open(output_normalized, "r") as f:
                lines = f.readlines()
            
            seqs = [lines[i].strip() for i in range(1, len(lines), 2)]
            assert seqs[0] == r1_revcomp  # R1 is stored as reverse complement
            assert seqs[1] == r2_original  # R2 is unchanged
    
    def test_inspect_shows_strand_info(self):
        """Test that inspect command shows strand-specific information."""
        import io
        import sys
        
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create a simple strand-specific file
            r1_path = f"{tmpdir}/R1.fastq"
            with open(r1_path, "wb") as f:
                f.write(b"@read1\nACGT\n+\nIIII\n")
            
            zna_path = f"{tmpdir}/strand.zna"
            
            class EncArgs:
                files = [r1_path]
                interleaved = False
                fasta = False
                fastq = False
                read_group = "inspect_test"
                description = ""
                strand_specific = True
                strand_normalize = True
                read1_sense = False      # R1 is antisense (default)
                read2_antisense = False  # R2 is sense (default)
                output = zna_path
                seq_len_bytes = 1
                block_size = 131072
                compress_flag = False
                level = 3
            
            encode_command(EncArgs())
            
            # Capture inspect output
            captured_output = io.StringIO()
            sys.stdout = captured_output
            
            class InspArgs:
                input = zna_path
            
            inspect_command(InspArgs())
            
            sys.stdout = sys.__stdout__
            output = captured_output.getvalue()
            
            # Verify strand info is displayed
            assert "Strand Specific:  True" in output
            assert "R1 Antisense:     True" in output
            assert "R2 Antisense:     False" in output


# --- ZNA -> ZNA Re-encode Tests ---

class TestReencodeKeepsFullFragment:
    """`zna encode` of a ZNA file is a COPY, and must carry IS_FULL_FRAGMENT.

    The normalized case was handled; the *un*-normalized one was not, and fell through
    a gap nobody had a test for. `preserve_normalization` was
    ``is_reencoding and input_header.strand_normalized``, so re-encoding a file that was
    never normalized took the plain 4-tuple path, where IS_FULL_FRAGMENT has no source
    at all and was simply re-derived as False.

    Measured before the fix: a merged read encoded with `--treat-unpaired-as-merged`
    read back as ends (True, True); re-encoding that .zna gave (True, False). Silent,
    exit 0 — a whole-molecule record downgraded to a one-ended one, which is lost
    supervision rather than visible corruption.
    """

    def _args(self, files, out, **over):
        class Args:
            pass
        a = Args()
        a.files, a.output = files, out
        a.interleaved = a.fasta = False
        a.fastq = True
        a.read_group, a.description = "rg", ""
        a.strand_specific = a.read1_sense = a.read2_antisense = False
        a.strand_normalize = False
        a.npolicy = None
        a.seq_len_bytes, a.block_size = 2, 131072
        a.compress_flag, a.level, a.quiet = False, 3, True
        a.treat_unpaired_as_merged = False
        for k, v in over.items():
            setattr(a, k, v)
        return a

    def _ends(self, path):
        with open(path, "rb") as fh:
            return [(r[4], r[5]) for r in ZnaReader(fh).records(with_ends=True)]

    def test_a_merged_read_stays_a_whole_fragment_across_a_reencode(self, tmp_path):
        fq = tmp_path / "m.fq"
        fq.write_text("".join(f"@merged{i} merged_20_0\nACGTACGTACGTACGTACGT\n+\n"
                              f"{'I' * 20}\n" for i in range(6)))
        first, second = tmp_path / "a.zna", tmp_path / "b.zna"
        encode_command(self._args([str(fq)], str(first),
                                  treat_unpaired_as_merged=True))
        assert self._ends(first) == [(True, True)] * 6

        encode_command(self._args([str(first)], str(second)))
        assert self._ends(second) == [(True, True)] * 6, \
            "re-encoding a non-normalized file dropped IS_FULL_FRAGMENT"

    def test_reencoding_is_idempotent_over_several_generations(self, tmp_path):
        fq = tmp_path / "m.fq"
        fq.write_text("".join(f"@merged{i} merged_20_0\nACGTACGTACGTACGTACGT\n+\n"
                              f"{'I' * 20}\n" for i in range(6)))
        path = tmp_path / "gen0.zna"
        encode_command(self._args([str(fq)], str(path),
                                  treat_unpaired_as_merged=True))
        expected = self._ends(path)
        for gen in range(1, 4):
            nxt = tmp_path / f"gen{gen}.zna"
            encode_command(self._args([str(path)], str(nxt)))
            assert self._ends(nxt) == expected, f"generation {gen} diverged"
            path = nxt


class TestReencode:
    """Re-encoding an already strand-normalized ZNA must copy its orientation.

    Orientation is not idempotent: the writer used to re-derive it from the
    header, reverse-complementing one mate of every pair a second time and
    silently un-normalizing the file while the header still claimed otherwise.
    """

    FRAG_LEN = 300
    READ_LEN = 100
    N_PAIRS = 120

    def _fragments(self):
        rng = random.Random(4242)
        return ["".join(rng.choice("ACGT") for _ in range(self.FRAG_LEN))
                for _ in range(self.N_PAIRS)]

    def _write_normalized(self, path, frags, **header_kwargs):
        """Write FR pairs as sequenced, letting the writer normalize them."""
        L, R = self.FRAG_LEN, self.READ_LEN
        kwargs = dict(read_group="reenc", strand_specific=False,
                      strand_normalized=True, seq_len_bytes=1)
        kwargs.update(header_kwargs)
        header = ZnaHeader(**kwargs)
        with open(path, "wb") as fh:
            with ZnaWriter(fh, header) as writer:
                for f in frags:
                    writer.write_record(f[:R], True, True, False)
                    writer.write_record(reverse_complement(f[L - R:]), True, False, True)

    def _read(self, path, **kwargs):
        with open(path, "rb") as fh:
            return list(ZnaReader(fh).records(**kwargs))

    def _args(self, in_path, out_path, **overrides):
        class Args:
            files = [in_path]
            interleaved = False
            fasta = False
            fastq = False
            read_group = "Unknown"      # sentinel: inherit from input header
            description = ""            # sentinel: inherit
            strand_specific = False
            strand_normalize = False
            read1_sense = False
            read2_antisense = False
            output = out_path
            seq_len_bytes = 2           # sentinel: inherit
            block_size = 512            # deliberately unlike the input's blocks
            compress_flag = False
            level = 3
            quiet = True
        for k, v in overrides.items():
            setattr(Args, k, v)
        return Args()

    def test_reencode_preserves_orientation(self, tmp_path):
        """A re-encode is a faithful copy: same sequences, same IS_RC flags.

        The block size differs from the input's on purpose. The accel backend
        seeds its coin to a constant per block, so with identical blocks a
        second pass re-flips the same mate and the damage shows up only in the
        sequences — the IS_RC flags still round-trip. Changing the block
        boundaries makes the second pass pick different mates, so the flags go
        wrong too and the assertion covers strictly more.
        """
        frags = self._fragments()
        src = str(tmp_path / "norm.zna")
        out = str(tmp_path / "reencoded.zna")
        self._write_normalized(src, frags)

        before = self._read(src, with_rc=True)
        encode_command(self._args(src, out))
        after = self._read(out, with_rc=True)

        assert after == before

        # And the pairs are still co-oriented: both mates in one frame.
        L, R = self.FRAG_LEN, self.READ_LEN
        for i, f in enumerate(frags):
            a, b = after[2 * i][0], after[2 * i + 1][0]
            assert (a == f[:R]) == (b == f[L - R:]), f"pair {i} lost co-orientation"

    def test_reencode_twice_is_stable(self, tmp_path):
        """Repeated re-encodes must converge, not toggle."""
        frags = self._fragments()
        p0 = str(tmp_path / "a.zna")
        p1 = str(tmp_path / "b.zna")
        p2 = str(tmp_path / "c.zna")
        self._write_normalized(p0, frags)

        encode_command(self._args(p0, p1))
        encode_command(self._args(p1, p2, block_size=4096))

        assert self._read(p2, with_rc=True) == self._read(p0, with_rc=True)

    def test_reencode_restores_original_orientation(self, tmp_path):
        """--restore-strand on a re-encoded file must still recover the input reads."""
        frags = self._fragments()
        L, R = self.FRAG_LEN, self.READ_LEN
        src = str(tmp_path / "norm.zna")
        out = str(tmp_path / "reencoded.zna")
        self._write_normalized(src, frags)
        encode_command(self._args(src, out))

        restored = [rec[0] for rec in self._read(out, restore_strand=True)]
        expected = []
        for f in frags:
            expected.append(f[:R])
            expected.append(reverse_complement(f[L - R:]))
        assert restored == expected

    def test_reencode_labeled_input_refused(self, tmp_path):
        """A labeled input would lose its label columns — refuse, don't crash."""
        from zna.dtypes import LabelDef, DTYPE_BY_CODE

        src = str(tmp_path / "labeled.zna")
        header = ZnaHeader(
            read_group="lab", seq_len_bytes=1,
            labels=(LabelDef(label_id=0, name="score", description="",
                             dtype=DTYPE_BY_CODE['i'], missing=0),),
        )
        with open(src, "wb") as fh:
            with ZnaWriter(fh, header) as writer:
                writer.write_record("ACGTACGT", False, False, False, labels=(1,))

        with pytest.raises(SystemExit):
            encode_command(self._args(src, str(tmp_path / "out.zna")))

    def test_reencode_with_label_flag_refused(self, tmp_path):
        """--label on a ZNA input used to parse the binary as FASTQ and write
        an empty file with a zero exit status."""
        src = str(tmp_path / "norm.zna")
        self._write_normalized(src, self._fragments()[:4])

        with pytest.raises(SystemExit):
            encode_command(self._args(src, str(tmp_path / "out.zna"),
                                      label=["NH:C"], label_desc=[]))

    def test_reencode_cannot_change_strand_specificity(self, tmp_path):
        """Orientation already applied cannot be re-derived under new rules."""
        src = str(tmp_path / "norm.zna")
        self._write_normalized(src, self._fragments()[:4])

        with pytest.raises(SystemExit):
            encode_command(self._args(src, str(tmp_path / "out.zna"),
                                      strand_specific=True, strand_normalize=True))

    def test_reencode_unnormalized_input_still_normalizes(self, tmp_path):
        """Pass-through must be conditional: an un-normalized input asked to be
        normalized still gets orientation derived for it."""
        L, R = self.FRAG_LEN, self.READ_LEN
        frags = self._fragments()[:20]
        src = str(tmp_path / "plain.zna")
        header = ZnaHeader(read_group="plain", seq_len_bytes=1)
        with open(src, "wb") as fh:
            with ZnaWriter(fh, header) as writer:
                for f in frags:
                    writer.write_record(f[:R], True, True, False)
                    writer.write_record(reverse_complement(f[L - R:]), True, False, True)
        assert not any(rec[4] for rec in self._read(src, with_rc=True))

        out = str(tmp_path / "out.zna")
        encode_command(self._args(src, out, strand_normalize=True))

        after = self._read(out, with_rc=True)
        with open(out, "rb") as fh:
            assert ZnaReader(fh).header.strand_normalized
        for i in range(len(frags)):
            assert after[2 * i][4] != after[2 * i + 1][4], "exactly one mate must be RC'd"


# --- Labeled encoding must preserve pairing ---

class TestLabeledPairing:
    """Labels must not cost pairing information.

    ``stream_inputs_labeled`` used to have its own single-file FASTQ reader and
    flag every record unpaired, so ``--interleaved`` was silently ignored and the
    R1/R2 flags were discarded. A consumer then cannot tell a mate (one true
    fragment end) from a merged read (two), which is the difference between
    correct endpoint supervision and endpoint tokens on interior positions.
    """

    # A mixed interleaved stream as a read merger emits it: adjacent /1,/2 pairs
    # plus merged singles whose pair suffix has been stripped.
    MIXED = (
        b"@readA/1\tNH:i:3\nAAAACCCCGGGGTTTT\n+\nIIIIIIIIIIIIIIII\n"
        b"@readA/2\tNH:i:3\nTTTTGGGGCCCCAAAA\n+\nIIIIIIIIIIIIIIII\n"
        b"@readB merged_20_8\tNH:i:7\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n"
        b"@readC/1\tNH:i:1\nGGGGAAAACCCCTTTT\n+\nIIIIIIIIIIIIIIII\n"
        b"@readC/2\tNH:i:1\nCCCCTTTTGGGGAAAA\n+\nIIIIIIIIIIIIIIII\n"
    )
    EXPECTED_FLAGS = [
        (True, True, False),    # readA R1
        (True, False, True),    # readA R2
        (False, False, False),  # readB merged single
        (True, True, False),    # readC R1
        (True, False, True),    # readC R2
    ]

    def _args(self, files, out, **overrides):
        class Args:
            interleaved = False
            fasta = False
            fastq = True
            read_group = "lab"
            description = ""
            strand_specific = False
            strand_normalize = True
            read1_sense = False
            read2_antisense = False
            output = out
            seq_len_bytes = 1
            block_size = 131072
            compress_flag = False
            level = 3
            quiet = True
            label = None
            label_desc = None
            label_defs = None
        Args.files = files
        for k, v in overrides.items():
            setattr(Args, k, v)
        return Args()

    def _flags(self, path):
        with open(path, "rb") as fh:
            return [(r[1], r[2], r[3]) for r in ZnaReader(fh).records()]

    def test_labeled_interleaved_preserves_pairing(self, tmp_path):
        """The regression test: labels must not change the pairing flags."""
        src = tmp_path / "mixed.fq"
        src.write_bytes(self.MIXED)
        plain = str(tmp_path / "plain.zna")
        labeled = str(tmp_path / "labeled.zna")

        encode_command(self._args([str(src)], plain, interleaved=True))
        encode_command(self._args([str(src)], labeled, interleaved=True, label=["NH:i"]))

        assert self._flags(plain) == self.EXPECTED_FLAGS
        assert self._flags(labeled) == self.EXPECTED_FLAGS

        # ...and the labels still line up with their sequences.
        with open(labeled, "rb") as fh:
            recs = list(ZnaReader(fh).records())
        assert [r[4] for r in recs] == [(3,), (3,), (7,), (1,), (1,)]

    def test_labeled_interleaved_normalizes_per_pair(self, tmp_path):
        """With pairing restored, unstranded normalization RC's exactly one mate
        per pair again instead of drawing an independent coin per record.

        Uses many pairs on purpose: with only a couple, an independent-coin
        encoder passes this by chance. At 60 pairs that is 2**-60.
        """
        n_pairs = 60
        rec = []
        for i in range(n_pairs):
            rec.append(b"@read%d/1\tNH:i:1\nAAAACCCCGGGGTTTT\n+\nIIIIIIIIIIIIIIII\n" % i)
            rec.append(b"@read%d/2\tNH:i:1\nTTTTGGGGCCCCAAAA\n+\nIIIIIIIIIIIIIIII\n" % i)
        src = tmp_path / "pairs.fq"
        src.write_bytes(b"".join(rec))
        out = str(tmp_path / "labeled.zna")
        encode_command(self._args([str(src)], out, interleaved=True, label=["NH:i"]))

        with open(out, "rb") as fh:
            recs = list(ZnaReader(fh).records(with_rc=True))
        assert len(recs) == 2 * n_pairs
        for i in range(n_pairs):
            assert recs[2 * i][4] != recs[2 * i + 1][4], f"pair {i}: not exactly one RC"

    def test_labeled_paired_files_preserves_pairing(self, tmp_path):
        """The two-file labeled path must emit R1/R2, not two singles."""
        r1 = tmp_path / "r1.fq"
        r2 = tmp_path / "r2.fq"
        r1.write_bytes(b"@readA/1\tNH:i:2\nAAAACCCCGGGGTTTT\n+\nIIIIIIIIIIIIIIII\n")
        r2.write_bytes(b"@readA/2\tNH:i:2\nTTTTGGGGCCCCAAAA\n+\nIIIIIIIIIIIIIIII\n")
        out = str(tmp_path / "pe.zna")

        encode_command(self._args([str(r1), str(r2)], out, label=["NH:i"]))

        assert self._flags(out) == [(True, True, False), (True, False, True)]
        with open(out, "rb") as fh:
            assert [r[4] for r in ZnaReader(fh).records()] == [(2,), (2,)]

    def test_labeled_single_end_still_unpaired(self, tmp_path):
        """Genuine single-end labeled input is unchanged by the refactor."""
        src = tmp_path / "se.fq"
        src.write_bytes(b"@readA\tNH:i:5\nAAAACCCCGGGGTTTT\n+\nIIIIIIIIIIIIIIII\n")
        out = str(tmp_path / "se.zna")
        encode_command(self._args([str(src)], out, label=["NH:i"]))
        assert self._flags(out) == [(False, False, False)]

    def test_labeled_fasta_refused(self, tmp_path):
        """Labels come from SAM tags in FASTQ headers; FASTA cannot carry them."""
        src = tmp_path / "in.fa"
        src.write_bytes(b">readA\nACGTACGT\n")
        with pytest.raises(SystemExit):
            encode_command(self._args([str(src)], str(tmp_path / "o.zna"),
                                      fasta=True, fastq=False, label=["NH:i"]))


# --- Fragment-span flag and pair-atomic filtering ---

class TestFullFragment:
    """A record that spans its whole fragment has BOTH edges as true boundaries.

    ``IS_RC`` can only name one edge, so without this flag a consumer placing
    fragment-end supervision either under-marks merged reads or over-marks
    genuine single-end reads — and it cannot tell the two apart.
    """

    def _args(self, files, out, **overrides):
        class Args:
            interleaved = True
            fasta = False
            fastq = True
            read_group = "ff"
            description = ""
            strand_specific = False
            strand_normalize = True
            read1_sense = False
            read2_antisense = False
            output = out
            seq_len_bytes = 1
            block_size = 131072
            compress_flag = False
            level = 3
            quiet = True
            label = None
            label_desc = None
            label_defs = None
            npolicy = "trim3"
            treat_unpaired_as_merged = False
        Args.files = files
        for k, v in overrides.items():
            setattr(Args, k, v)
        return Args()

    def _write_fq(self, path, records):
        with open(path, "w") as fh:
            for name, seq in records:
                fh.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")

    def _ends(self, path):
        with open(path, "rb") as fh:
            return [(r[4], r[5]) for r in ZnaReader(fh).records(with_ends=True)]

    def test_full_overlap_pair_detected(self, tmp_path):
        """Mates whose insert was <= the read length span the same interval, so
        they are exact reverse complements — detectable with no insert size."""
        frag = "ACGTACGTAAGGCCTTACGTACGT"
        src = tmp_path / "in.fq"
        self._write_fq(src, [
            ("full/1", frag), ("full/2", reverse_complement(frag)),   # full overlap
            ("ord/1", "AAAACCCCGGGGTTTT"), ("ord/2", "TTTTGGGGCCCCAAAA"),
        ])
        out = str(tmp_path / "o.zna")
        encode_command(self._args([str(src)], out))
        ends = self._ends(out)
        assert ends[0] == (True, True) and ends[1] == (True, True), "full-overlap pair"
        # the ordinary pair gets exactly one true edge per mate
        assert ends[2] != (True, True) and ends[3] != (True, True)
        assert ends[2][0] != ends[3][0]

    def test_treat_unpaired_as_merged_flag(self, tmp_path):
        """An unpaired record's span cannot be inferred, so it is declared."""
        src = tmp_path / "in.fq"
        self._write_fq(src, [("merged1 merged_16_0", "ACGTACGTAAGGCCTT")])

        off = str(tmp_path / "off.zna")
        encode_command(self._args([str(src)], off))
        assert self._ends(off)[0] != (True, True), "default: one real edge"

        on = str(tmp_path / "on.zna")
        encode_command(self._args([str(src)], on, treat_unpaired_as_merged=True))
        assert self._ends(on)[0] == (True, True), "declared: both edges real"

    def test_trim3_never_orphans_a_mate(self, tmp_path):
        """trim3 reshapes a record; it never removes one, so a pair stays a pair.

        This replaces a pair-atomicity test for the old `drop` policy. Dropping one mate
        would have left a lone paired record, which downstream cannot distinguish from a
        genuine single — the reason `drop` was fragment-atomic. trim3 has no such hazard
        by construction, and it salvages the rest of the read instead of discarding the
        whole fragment.
        """
        src = tmp_path / "in.fq"
        self._write_fq(src, [
            ("keep/1", "AAAACCCCGGGGTTTT"), ("keep/2", "TTTTGGGGCCCCAAAA"),
            ("bad/1", "AAAACCCCGGGGTTTN"), ("bad/2", "TTTTGGGGCCCCAAAA"),  # N in R1 only
        ])
        out = str(tmp_path / "o.zna")
        encode_command(self._args([str(src)], out))
        with open(out, "rb") as fh:
            recs = list(ZnaReader(fh).records())
        assert len(recs) == 4, "trim3 must keep both mates of both fragments"
        assert [(r[1], r[2], r[3]) for r in recs] == [
            (True, True, False), (True, False, True),
            (True, True, False), (True, False, True)]
        assert recs[2][0] == "AAAACCCCGGGGTTT", "R1 must be cut at its N, 3' only"
        assert recs[3][0] == "TTTTGGGGCCCCAAAA", "the clean mate must be untouched"
        assert all("N" not in r[0] for r in recs)

    def test_trim3_treats_each_record_independently(self, tmp_path):
        """A no-call in one single must not disturb its neighbours."""
        src = tmp_path / "in.fq"
        self._write_fq(src, [
            ("a", "AAAACCCCGGGGTTTT"), ("b", "AAAACCCCGGGGTTTN"),
            ("c", "TTTTGGGGCCCCAAAA"),
        ])
        out = str(tmp_path / "o.zna")
        encode_command(self._args([str(src)], out))
        with open(out, "rb") as fh:
            seqs = [r[0] for r in ZnaReader(fh).records()]
        # This fixture normalizes, so a record may be stored reverse-complemented.
        # What trim3 owes us is the LENGTH and the absence of no-calls.
        assert [len(x) for x in seqs] == [16, 15, 16], seqs
        assert all("N" not in x for x in seqs)
        assert seqs[0] != seqs[2], "the two clean records should not collide"

    def _norm_zna(self, path):
        header = ZnaHeader(read_group="g", strand_specific=False,
                           strand_normalized=True, seq_len_bytes=1)
        with open(path, "wb") as fh:
            with ZnaWriter(fh, header) as w:
                for i in range(6):
                    w.write_record("ACGTACGTAAGG", True, i % 2 == 0, i % 2 == 1)
                w.write_record("TTTTGGGGCCCC", False, False, False,
                               is_full_fragment=True)

    def test_decode_warns_when_reencode_would_double_normalize(self, tmp_path, capsys):
        """`zna decode` emits the normalized frame by default; re-encoding that
        with --strand-normalize applies orientation twice, silently."""
        src = str(tmp_path / "n.zna")
        self._norm_zna(src)

        class Args:
            input = src
            output = str(tmp_path / "out.fa")
            quiet = False
            gzip = False
            labels = False
            restore_strand = False

        decode_command(Args())
        err = capsys.readouterr().err
        assert "not idempotent" in err and "--restore-strand" in err

        Args.restore_strand = True
        Args.output = str(tmp_path / "out2.fa")
        decode_command(Args())
        assert "not idempotent" not in capsys.readouterr().err

    def test_inspect_counts_reports_fragment_geometry(self, tmp_path, capsys):
        """--counts must expose the (mate x IS_RC) cross-tab and full-fragment
        tally; a bare RC total is equally consistent with a healthy and a broken
        file."""
        src = str(tmp_path / "n.zna")
        self._norm_zna(src)

        class InspArgs:
            input = src
            counts = True

        inspect_command(InspArgs())
        out = capsys.readouterr().out
        assert "Full-fragment:" in out
        assert "Fragment Geometry" in out
        assert "unstranded" in out


# --- N-Policy CLI Tests ---

class TestNPolicyCLI:
    """Test CLI N-policy functionality."""
    
    def test_npolicy_trim3(self):
        """`--npolicy trim3` (the default) cuts each record at its first no-call."""
        with tempfile.TemporaryDirectory() as tmpdir:
            fastq_path = f"{tmpdir}/test.fastq"
            with open(fastq_path, "w") as f:
                f.write("@r1\nACGTACGT\n+\nIIIIIIII\n")
                f.write("@r2\nACGNACGT\n+\nIIIIIIII\n")
                f.write("@r3\nNACGTACG\n+\nIIIIIIII\n")
            zna_path = f"{tmpdir}/out.zna"

            class Args:
                files = [fastq_path]
                interleaved = False
                fasta = False
                fastq = True
                read_group = "test"
                description = ""
                strand_specific = False
                read1_sense = False
                read2_antisense = False
                npolicy = "trim3"
                output = zna_path
                seq_len_bytes = 2
                block_size = 131072
                compress_flag = False
                level = 3
                quiet = True

            encode_command(Args())
            with open(zna_path, "rb") as f:
                seqs = [r[0] for r in ZnaReader(f).records()]
            # Cut at the N, keeping the 5' side — including down to nothing.
            assert seqs == ["ACGTACGT", "ACG", ""]

    def test_npolicy_random(self):
        """Test that N nucleotides are replaced with random bases."""
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = f"{tmpdir}/input.fasta"
            with open(fasta_path, "w") as f:
                f.write(">seq_with_N\nACNGT\n")
            
            zna_path = f"{tmpdir}/output.zna"
            
            class Args:
                files = [fasta_path]
                interleaved = False
                fasta = True
                fastq = False
                read_group = "test"
                description = ""
                strand_specific = False
                read1_sense = False
                read2_antisense = False
                npolicy = "random"
                output = zna_path
                seq_len_bytes = 1
                block_size = 131072
                compress_flag = False
                level = 3
            
            encode_command(Args())
            
            # Decode and verify N was replaced with a valid base
            with open(zna_path, "rb") as f:
                reader = ZnaReader(f)
                records = list(reader.records())
                
                assert len(records) == 1
                seq = records[0][0]
                assert len(seq) == 5
                assert seq[0:2] == "AC"
                assert seq[3:] == "GT"
                assert seq[2] in "ACGT"  # Random replacement


# --- Mixed Interleaved Tests ---

class TestZnaInputIsOnlyValidAlone:
    """A .zna is input only in single-file re-encode mode.

    ``is_reencoding`` is ``len(files) == 1 and is_zna_file(...)``, so with two .zna files
    the check never fired, the binary went to the FASTQ parser, and encode wrote a valid
    0-record file at exit 0 — warning only about *format inference*, never about the real
    problem.
    """

    def _args(self, files, out):
        class Args:
            pass
        a = Args()
        a.files, a.output = files, out
        a.interleaved = a.fasta = a.fastq = False
        a.read_group, a.description = "rg", ""
        a.strand_specific = a.read1_sense = a.read2_antisense = False
        a.npolicy = None
        a.seq_len_bytes, a.block_size = 2, 131072
        a.compress_flag, a.level, a.quiet = False, 3, True
        return a

    def _make_zna(self, tmp_path, name):
        src = tmp_path / "in.fq"
        src.write_text("@r1\nACGTACGT\n+\nIIIIIIII\n")
        out = tmp_path / name
        a = self._args([str(src)], str(out))
        a.fastq = True
        encode_command(a)
        return out

    @pytest.mark.parametrize("second", ["zna", "fastq"])
    def test_a_zna_alongside_another_file_is_refused(self, tmp_path, second):
        z = self._make_zna(tmp_path, "a.zna")
        if second == "zna":
            other = self._make_zna(tmp_path, "b.zna")
        else:
            other = tmp_path / "r2.fq"
            other.write_text("@r1\nACGTACGT\n+\nIIIIIIII\n")
        with pytest.raises(SystemExit) as exc:
            encode_command(self._args([str(z), str(other)], str(tmp_path / "o.zna")))
        assert "already a ZNA file" in str(exc.value)

    def test_a_lone_zna_still_reencodes(self, tmp_path):
        z = self._make_zna(tmp_path, "a.zna")
        out = tmp_path / "b.zna"
        encode_command(self._args([str(z)], str(out)))
        with open(out, "rb") as fh:
            assert len(list(ZnaReader(fh).records())) == 1


class TestPairedFilesMustHaveEqualCounts:
    """Two-file input with unequal record counts is an error, not a silent truncation.

    The loop was ``for s1, s2 in zip(p1, p2)``, which stops at the shorter file: an R2
    with fewer records than R1 dropped every remaining R1 read and exited 0.  A whole
    half-library could vanish with every stage green.

    Note ``zip`` also makes the bug unfixable from outside — by the time it stops it has
    already pulled and discarded the extra value from the longer iterator, so the
    off-by-one case looks like a clean finish.  Hence the explicit two-``next`` loop.
    """

    def _write(self, path, n, mate):
        with open(path, "w") as f:
            for i in range(n):
                f.write(f"@r{i}/{mate}\nACGTACGT\n+\nIIIIIIII\n")

    def _args(self, files, out):
        class Args:
            pass
        a = Args()
        a.files, a.output = files, out
        a.interleaved = a.fasta = False
        a.fastq = True
        a.read_group, a.description = "rg", ""
        a.strand_specific = a.read1_sense = a.read2_antisense = False
        a.npolicy = None
        a.seq_len_bytes, a.block_size = 2, 131072
        a.compress_flag, a.level, a.quiet = False, 3, True
        return a

    @pytest.mark.parametrize("n1,n2", [(3, 2), (2, 3)])
    def test_unequal_counts_are_refused(self, tmp_path, n1, n2):
        r1, r2 = tmp_path / "r1.fq", tmp_path / "r2.fq"
        self._write(r1, n1, 1)
        self._write(r2, n2, 2)
        with pytest.raises(SystemExit) as exc:
            encode_command(self._args([str(r1), str(r2)], str(tmp_path / "o.zna")))
        msg = str(exc.value)
        assert "ran out of records before" in msg
        # The message names the file that ended early, and the one that outlasted it.
        short, long_ = (r2, r1) if n1 > n2 else (r1, r2)
        assert f"{short} ran out of records before {long_}" in msg

    def test_equal_counts_still_encode(self, tmp_path):
        r1, r2 = tmp_path / "r1.fq", tmp_path / "r2.fq"
        self._write(r1, 3, 1)
        self._write(r2, 3, 2)
        out = tmp_path / "o.zna"
        encode_command(self._args([str(r1), str(r2)], str(out)))
        with open(out, "rb") as fh:
            recs = list(ZnaReader(fh).records())
        assert len(recs) == 6
        assert [(r[1], r[2], r[3]) for r in recs[:2]] == [
            (True, True, False), (True, False, True)]


class TestMergedReadsAreNeverMates:
    """A merged read is a whole molecule, so it must never be paired with a neighbour.

    Interleaved pairing matches on the read NAME, and names collide.  Two independently
    merged molecules that happened to share a name were paired into one fragment: both
    lost ``IS_FULL_FRAGMENT``, both gained ``IS_PAIRED``/``IS_READ1``/``IS_READ2``, and
    two different molecules were encoded as each other's mate — with a zero exit status
    and no warning.  Duplicate read names are ordinary in real data (concatenated lanes,
    some SRA dumps, UMI-collapsed files) and nothing upstream forbids them.

    `zna merge`, fastp and khorana's `parse_merged_fastq` all already agree on the
    ``merged_<n1>_<n2>`` token, so :func:`_read_key` reads it and reports
    :data:`MERGED_SINGLE`.  Nothing errors: the records simply come out right.
    """

    def test_read_key_reports_a_merged_read(self):
        from zna.cli import _read_key, MERGED_SINGLE
        assert _read_key(b"f12 merged_90_0") == (b"f12", MERGED_SINGLE)
        assert _read_key(b"f12\tZI:i:7 merged_90_0") == (b"f12", MERGED_SINGLE)
        # A mate keeps its mate number and never pays for the check.
        assert _read_key(b"f12/1\tZI:i:7") == (b"f12", 1)
        assert _read_key(b"f12/2\tZI:i:7") == (b"f12", 2)
        # An ordinary single is still an ordinary single.
        assert _read_key(b"f12 some comment") == (b"f12", 0)
        assert _read_key(b"f12") == (b"f12", 0)

    def test_two_merged_reads_sharing_a_name_stay_separate(self):
        from zna.cli import _pair_interleaved, _read_key
        recs = [(*_read_key(h), h) for h in
                (b"dup merged_90_0", b"dup merged_95_0", b"uniq merged_100_0")]
        out = list(_pair_interleaved(recs))
        assert len(out) == 3
        for _payload, is_paired, is_read1, is_read2 in out:
            assert not is_paired and not is_read1 and not is_read2

    def test_genuine_mates_sharing_a_name_with_a_merged_read_still_pair(self):
        """The fix must not cost a real pair its pairing."""
        from zna.cli import _pair_interleaved, _read_key
        recs = [(*_read_key(h), h) for h in
                (b"dup merged_90_0", b"dup/1", b"dup/2")]
        out = list(_pair_interleaved(recs))
        assert [(r[1], r[2], r[3]) for r in out] == [
            (False, False, False),        # the merged read, alone
            (True, True, False),          # ...and the pair that shares its name
            (True, False, True),
        ]


class TestMixedInterleaved:
    """Tests for interleaved mode with mixed paired-end and single-end reads."""

    def test_mixed_paired_and_single(self):
        """Test interleaved mode with both paired and single-end reads."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create test file with mixed reads
            fastq_path = f"{tmpdir}/mixed.fastq"
            with open(fastq_path, "w") as f:
                f.write("@read1/1\nACGTACGT\n+\nIIIIIIII\n")
                f.write("@read1/2\nTGCATGCA\n+\nIIIIIIII\n")
                f.write("@single1\nGGGGAAAA\n+\nIIIIIIII\n")
                f.write("@read2/1\nCCCCTTTT\n+\nIIIIIIII\n")
                f.write("@read2/2\nAAAACCCC\n+\nIIIIIIII\n")
                f.write("@single2\nTTTTGGGG\n+\nIIIIIIII\n")
            
            zna_path = f"{tmpdir}/output.zna"
            
            class Args:
                files = [fastq_path]
                interleaved = True
                fasta = False
                fastq = True
                read_group = "mixed"
                description = "test"
                strand_specific = False
                read1_sense = False
                read2_antisense = False
                npolicy = None
                output = zna_path
                seq_len_bytes = 2
                block_size = 131072
                compress_flag = False
                level = 3
                quiet = True
            
            encode_command(Args())
            
            # Decode and verify structure
            with open(zna_path, "rb") as f:
                reader = ZnaReader(f)
                records = list(reader.records())
                
                # Should have 6 records: 2 pairs + 2 singles
                assert len(records) == 6
                
                # Check first pair
                assert records[0] == ("ACGTACGT", True, True, False)
                assert records[1] == ("TGCATGCA", True, False, True)
                
                # Check first single
                assert records[2] == ("GGGGAAAA", False, False, False)
                
                # Check second pair
                assert records[3] == ("CCCCTTTT", True, True, False)
                assert records[4] == ("AAAACCCC", True, False, True)
                
                # Check second single
                assert records[5] == ("TTTTGGGG", False, False, False)
    
    def test_all_paired_reads(self):
        """Test interleaved mode with only paired reads."""
        with tempfile.TemporaryDirectory() as tmpdir:
            fastq_path = f"{tmpdir}/paired.fastq"
            with open(fastq_path, "w") as f:
                f.write("@read1/1\nACGT\n+\nIIII\n")
                f.write("@read1/2\nTGCA\n+\nIIII\n")
                f.write("@read2/1\nGGGG\n+\nIIII\n")
                f.write("@read2/2\nCCCC\n+\nIIII\n")
            
            zna_path = f"{tmpdir}/output.zna"
            
            class Args:
                files = [fastq_path]
                interleaved = True
                fasta = False
                fastq = True
                read_group = "paired"
                description = ""
                strand_specific = False
                read1_sense = False
                read2_antisense = False
                npolicy = None
                output = zna_path
                seq_len_bytes = 1
                block_size = 131072
                compress_flag = False
                level = 3
                quiet = True
            
            encode_command(Args())
            
            with open(zna_path, "rb") as f:
                reader = ZnaReader(f)
                records = list(reader.records())
                
                assert len(records) == 4
                assert all(rec[1] for rec in records)  # All paired
    
    def test_all_single_reads(self):
        """Test interleaved mode with only single-end reads."""
        with tempfile.TemporaryDirectory() as tmpdir:
            fastq_path = f"{tmpdir}/singles.fastq"
            with open(fastq_path, "w") as f:
                f.write("@single1\nACGT\n+\nIIII\n")
                f.write("@single2\nTGCA\n+\nIIII\n")
                f.write("@single3\nGGGG\n+\nIIII\n")
            
            zna_path = f"{tmpdir}/output.zna"
            
            class Args:
                files = [fastq_path]
                interleaved = True
                fasta = False
                fastq = True
                read_group = "singles"
                description = ""
                strand_specific = False
                read1_sense = False
                read2_antisense = False
                npolicy = None
                output = zna_path
                seq_len_bytes = 1
                block_size = 131072
                compress_flag = False
                level = 3
                quiet = True
            
            encode_command(Args())
            
            with open(zna_path, "rb") as f:
                reader = ZnaReader(f)
                records = list(reader.records())
                
                assert len(records) == 3
                assert all(not rec[1] for rec in records)  # All single
    
    def test_read_name_without_suffix(self):
        """Test reads without /1 or /2 suffix are treated as single."""
        with tempfile.TemporaryDirectory() as tmpdir:
            fastq_path = f"{tmpdir}/no_suffix.fastq"
            with open(fastq_path, "w") as f:
                f.write("@read1\nACGT\n+\nIIII\n")
                f.write("@read2\nTGCA\n+\nIIII\n")
            
            zna_path = f"{tmpdir}/output.zna"
            
            class Args:
                files = [fastq_path]
                interleaved = True
                fasta = False
                fastq = True
                read_group = "nosuffix"
                description = ""
                strand_specific = False
                read1_sense = False
                read2_antisense = False
                npolicy = None
                output = zna_path
                seq_len_bytes = 1
                block_size = 131072
                compress_flag = False
                level = 3
                quiet = True
            
            encode_command(Args())
            
            with open(zna_path, "rb") as f:
                reader = ZnaReader(f)
                records = list(reader.records())
                
                assert len(records) == 2
                assert all(not rec[1] for rec in records)  # All treated as single
    
    def test_fastp_merged_format(self):
        """Test handling of fastp output with merged notation."""
        with tempfile.TemporaryDirectory() as tmpdir:
            fastq_path = f"{tmpdir}/fastp.fastq"
            with open(fastq_path, "w") as f:
                # Merged read (no suffix)
                f.write("@read1 merged\nACGTACGT\n+\nIIIIIIII\n")
                # Unmerged pair
                f.write("@read2/1\nCCCCTTTT\n+\nIIIIIIII\n")
                f.write("@read2/2\nAAAACCCC\n+\nIIIIIIII\n")
                # Another merged
                f.write("@read3 merged\nGGGGAAAA\n+\nIIIIIIII\n")
            
            zna_path = f"{tmpdir}/output.zna"
            
            class Args:
                files = [fastq_path]
                interleaved = True
                fasta = False
                fastq = True
                read_group = "fastp"
                description = ""
                strand_specific = False
                read1_sense = False
                read2_antisense = False
                npolicy = None
                output = zna_path
                seq_len_bytes = 2
                block_size = 131072
                compress_flag = False
                level = 3
                quiet = True
            
            encode_command(Args())
            
            with open(zna_path, "rb") as f:
                reader = ZnaReader(f)
                records = list(reader.records())
                
                # Should have 4 records: 2 singles + 1 pair
                assert len(records) == 4
                
                # First merged/single
                assert records[0][1] == False  # not paired
                
                # Unmerged pair
                assert records[1] == ("CCCCTTTT", True, True, False)
                assert records[2] == ("AAAACCCC", True, False, True)
                
                # Second merged/single
                assert records[3][1] == False  # not paired


# --- SHUFFLE TESTS ---

class TestShuffle:
    """Tests for the ``zna shuffle`` command."""

    @staticmethod
    def _make_zna(path, records, compressed=True, strand_specific=False,
                  read1_antisense=False):
        """Helper: write a list of (seq, is_paired, is_read1, is_read2) to a ZNA file."""
        header = ZnaHeader(
            read_group="test",
            seq_len_bytes=2,
            compression_method=COMPRESSION_ZSTD if compressed else COMPRESSION_NONE,
            compression_level=3,
            strand_specific=strand_specific,
            read1_antisense=read1_antisense,
        )
        with open(path, "wb") as fh:
            with ZnaWriter(fh, header) as writer:
                for seq, is_paired, is_read1, is_read2 in records:
                    writer.write_record(seq, is_paired, is_read1, is_read2)

    @staticmethod
    def _read_zna(path):
        """Helper: read all records from a ZNA file."""
        with open(path, "rb") as fh:
            reader = ZnaReader(fh)
            return list(reader.records())

    @staticmethod
    def _make_args(**kwargs):
        """Build a namespace mimicking argparse output for shuffle_command."""
        defaults = dict(
            input=None,
            output=None,
            seed=42,
            buffer_size="1M",
            block_size="4M",
            tmp_dir=None,
            quiet=True,
        )
        defaults.update(kwargs)

        class Args:
            pass

        a = Args()
        for k, v in defaults.items():
            setattr(a, k, v)
        return a

    def test_single_end_shuffle(self, tmp_path):
        """Single-end records are shuffled; all records are preserved."""
        in_path = str(tmp_path / "in.zna")
        out_path = str(tmp_path / "out.zna")

        # Generate 100 distinct DNA sequences
        rng = __import__("random").Random(0)
        bases = "ACGT"
        seqs = ["".join(rng.choices(bases, k=40)) for _ in range(100)]
        records = [(s, False, False, False) for s in seqs]
        self._make_zna(in_path, records)

        args = self._make_args(input=in_path, output=out_path)
        shuffle_command(args)

        out_records = self._read_zna(out_path)
        assert len(out_records) == len(records)

        # All sequences should be present (set comparison)
        in_seqs = sorted(r[0] for r in records)
        out_seqs = sorted(r[0] for r in out_records)
        assert in_seqs == out_seqs

        # Order should differ (with 100 records and seed=42, vanishingly unlikely to match)
        in_order = [r[0] for r in records]
        out_order = [r[0] for r in out_records]
        assert in_order != out_order

    def test_paired_end_shuffle_preserves_pairs(self, tmp_path):
        """Paired-end R1+R2 stay adjacent after shuffling."""
        in_path = str(tmp_path / "in.zna")
        out_path = str(tmp_path / "out.zna")

        # 50 pairs: R1 and R2 share a common tag to verify pairing
        rng = __import__("random").Random(1)
        bases = "ACGT"
        records = []
        tags = ["".join(rng.choices(bases, k=20)) for _ in range(50)]
        for tag in tags:
            records.append(("AAAA" + tag + "AAAA", True, True, False))   # R1
            records.append(("TTTT" + tag + "TTTT", True, False, True))   # R2
        self._make_zna(in_path, records)

        args = self._make_args(input=in_path, output=out_path)
        shuffle_command(args)

        out_records = self._read_zna(out_path)
        assert len(out_records) == 100  # 50 pairs × 2 records

        # Verify every R1 is immediately followed by its matching R2
        for i in range(0, len(out_records), 2):
            r1_seq, r1_paired, r1_is_r1, r1_is_r2 = out_records[i]
            r2_seq, r2_paired, r2_is_r1, r2_is_r2 = out_records[i + 1]

            assert r1_paired and r1_is_r1 and not r1_is_r2, f"Record {i} should be R1"
            assert r2_paired and not r2_is_r1 and r2_is_r2, f"Record {i+1} should be R2"

            # The tag (positions 4:24) must match between R1 and R2
            r1_tag = r1_seq[4:24]
            r2_tag = r2_seq[4:24]
            assert r1_tag == r2_tag, f"Pair broken at {i}: R1 tag={r1_tag}, R2 tag={r2_tag}"

    def test_mixed_paired_and_single(self, tmp_path):
        """Mix of paired and single-end records: all preserved, pairs intact."""
        in_path = str(tmp_path / "in.zna")
        out_path = str(tmp_path / "out.zna")

        rng = __import__("random").Random(2)
        bases = "ACGT"
        records = []
        pair_tags = []
        # 20 pairs
        for _ in range(20):
            tag = "".join(rng.choices(bases, k=20))
            pair_tags.append(tag)
            records.append(("AA" + tag + "AA", True, True, False))
            records.append(("TT" + tag + "TT", True, False, True))
        # 30 singles
        for _ in range(30):
            records.append(("CC" + "".join(rng.choices(bases, k=20)) + "CC", False, False, False))

        self._make_zna(in_path, records)

        args = self._make_args(input=in_path, output=out_path)
        shuffle_command(args)

        out_records = self._read_zna(out_path)
        assert len(out_records) == 70  # 20*2 + 30

        # Check that all sequences are present
        in_seqs = sorted(r[0] for r in records)
        out_seqs = sorted(r[0] for r in out_records)
        assert in_seqs == out_seqs

        # Check every paired R1 is followed by its R2 with matching tag
        i = 0
        while i < len(out_records):
            seq, is_paired, is_r1, is_r2 = out_records[i]
            if is_paired and is_r1:
                r2 = out_records[i + 1]
                assert r2[1] and r2[3], f"Expected R2 after R1 at position {i}"
                # Tags at positions 2:22 must match
                assert seq[2:22] == r2[0][2:22], "Pair tag mismatch"
                i += 2
            else:
                assert not is_paired, f"Unexpected isolated paired record at {i}"
                i += 1

    def test_deterministic_with_seed(self, tmp_path):
        """Same seed produces identical output."""
        in_path = str(tmp_path / "in.zna")
        out1 = str(tmp_path / "out1.zna")
        out2 = str(tmp_path / "out2.zna")

        rng = __import__("random").Random(3)
        bases = "ACGT"
        records = [("".join(rng.choices(bases, k=30)), False, False, False) for _ in range(200)]
        self._make_zna(in_path, records)

        args1 = self._make_args(input=in_path, output=out1, seed=12345)
        shuffle_command(args1)

        args2 = self._make_args(input=in_path, output=out2, seed=12345)
        shuffle_command(args2)

        recs1 = self._read_zna(out1)
        recs2 = self._read_zna(out2)
        assert [r[0] for r in recs1] == [r[0] for r in recs2]

    def test_different_seed_different_order(self, tmp_path):
        """Different seeds produce different orderings."""
        in_path = str(tmp_path / "in.zna")
        out1 = str(tmp_path / "out1.zna")
        out2 = str(tmp_path / "out2.zna")

        rng = __import__("random").Random(4)
        bases = "ACGT"
        records = [("".join(rng.choices(bases, k=30)), False, False, False) for _ in range(200)]
        self._make_zna(in_path, records)

        args1 = self._make_args(input=in_path, output=out1, seed=1)
        shuffle_command(args1)

        args2 = self._make_args(input=in_path, output=out2, seed=2)
        shuffle_command(args2)

        order1 = [r[0] for r in self._read_zna(out1)]
        order2 = [r[0] for r in self._read_zna(out2)]
        assert order1 != order2

    def test_header_preserved(self, tmp_path):
        """Shuffle preserves the ZNA header metadata."""
        in_path = str(tmp_path / "in.zna")
        out_path = str(tmp_path / "out.zna")

        header = ZnaHeader(
            read_group="my_sample",
            description="test description",
            seq_len_bytes=2,
            strand_specific=True,
            read1_antisense=True,
            compression_method=COMPRESSION_ZSTD,
            compression_level=3,
        )
        records = [("ACGTACGT", False, False, False)] * 10
        with open(in_path, "wb") as fh:
            with ZnaWriter(fh, header) as w:
                for r in records:
                    w.write_record(*r)

        args = self._make_args(input=in_path, output=out_path)
        shuffle_command(args)

        with open(out_path, "rb") as fh:
            reader = ZnaReader(fh)
            h = reader.header
            assert h.read_group == "my_sample"
            assert h.description == "test description"
            assert h.strand_specific is True
            assert h.read1_antisense is True
            assert h.compression_method == COMPRESSION_ZSTD

    def test_multi_bucket_with_small_buffer(self, tmp_path):
        """A small buffer forces multiple buckets; result is still correct."""
        in_path = str(tmp_path / "in.zna")
        out_path = str(tmp_path / "out.zna")

        rng = __import__("random").Random(5)
        bases = "ACGT"
        # 500 distinct records
        records = [("".join(rng.choices(bases, k=40)), False, False, False) for _ in range(500)]
        self._make_zna(in_path, records)

        args = self._make_args(input=in_path, output=out_path, buffer_size="1K")
        shuffle_command(args)

        out_records = self._read_zna(out_path)
        assert len(out_records) == 500
        assert sorted(r[0] for r in records) == sorted(r[0] for r in out_records)

    def test_parse_block_size_gigabytes(self):
        """parse_block_size handles G/GB suffix."""
        assert parse_block_size("1G") == 1024 * 1024 * 1024
        assert parse_block_size("2GB") == 2 * 1024 * 1024 * 1024

    def test_shuffle_zna_directly(self, tmp_path):
        """shuffle_zna can be called as a library function."""
        in_path = str(tmp_path / "in.zna")
        out_path = str(tmp_path / "out.zna")

        rng = __import__("random").Random(6)
        records = [("" .join(rng.choices("ACGT", k=30)), False, False, False) for _ in range(80)]
        self._make_zna(in_path, records)

        written, n_records = shuffle_zna(in_path, out_path, seed=99, quiet=True)
        assert written == 80
        assert n_records == 80

        out_records = self._read_zna(out_path)
        assert sorted(r[0] for r in records) == sorted(r[0] for r in out_records)

    def test_encode_with_shuffle(self, tmp_path):
        """encode --shuffle produces a shuffled ZNA file."""
        fq_path = str(tmp_path / "input.fastq")
        out_path = str(tmp_path / "out.zna")

        # Write 100 distinct FASTQ reads
        rng = __import__("random").Random(7)
        seqs = ["" .join(rng.choices("ACGT", k=40)) for _ in range(100)]
        with open(fq_path, "w") as f:
            for i, s in enumerate(seqs):
                f.write(f"@read{i}\n{s}\n+\n{'I' * len(s)}\n")

        # Encode with --shuffle
        class Args:
            pass
        a = Args()
        a.files = [fq_path]
        a.output = out_path
        a.interleaved = False
        a.fasta = False
        a.fastq = True
        a.quiet = True
        a.read_group = "test"
        a.description = ""
        a.seq_len_bytes = 2
        a.strand_specific = False
        a.read1_sense = False
        a.read2_antisense = False
        a.compress_flag = None
        a.level = 3
        a.block_size = "4M"
        a.npolicy = "trim3"
        a.shuffle = True
        a.seed = 42
        a.shuffle_buffer_size = "1G"

        encode_command(a)

        out_records = self._read_zna(out_path)
        assert len(out_records) == 100

        # All sequences present
        assert sorted(seqs) == sorted(r[0] for r in out_records)

        # Order should differ from input
        assert seqs != [r[0] for r in out_records]


class TestGzipInput:
    """Gzip input is served by an external decompressor when one is available,
    which is worth ~2.2x but puts a second process in the data path.  These pin
    the properties that matters most: the two paths agree, and a corrupt input
    is rejected rather than silently truncated."""

    RECORDS = 500

    def _make_gz(self, path):
        rng = random.Random(4242)
        seqs = ["".join(rng.choices("ACGT", k=150)) for _ in range(self.RECORDS)]
        with gzip.open(path, "wb") as f:
            for i, s in enumerate(seqs):
                f.write(f"@read{i} comment\n{s}\n+\n{'I' * 150}\n".encode())
        return seqs

    def test_both_paths_agree(self, monkeypatch):
        from zna.cli import get_input_handle, parse_fastq
        with tempfile.TemporaryDirectory() as tmp:
            p = str(Path(tmp) / "in.fastq.gz")
            seqs = self._make_gz(p)

            with get_input_handle(p) as fh:
                via_default = list(parse_fastq(fh))

            monkeypatch.setenv("ZNA_NO_EXTERNAL_GZIP", "1")
            with get_input_handle(p) as fh:
                assert isinstance(fh, io.BufferedReader)
                via_module = list(parse_fastq(fh))

            assert via_default == seqs
            assert via_module == seqs

    def test_truncated_gzip_is_rejected(self, monkeypatch):
        """A partial stream plus a non-zero exit must not look like a short file.

        Without the exit-status check this is the worst kind of failure: encode
        would write a truncated ZNA and report success.
        """
        from zna.cli import get_input_handle, parse_fastq
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "in.fastq.gz"
            self._make_gz(str(p))
            raw = p.read_bytes()
            bad = Path(tmp) / "trunc.fastq.gz"
            bad.write_bytes(raw[: len(raw) // 2])

            with pytest.raises((OSError, EOFError, gzip.BadGzipFile)):
                fh = get_input_handle(str(bad))
                try:
                    list(parse_fastq(fh))
                finally:
                    fh.close()

            # ...and the pure-Python path must reject it too.
            monkeypatch.setenv("ZNA_NO_EXTERNAL_GZIP", "1")
            with pytest.raises((OSError, EOFError, gzip.BadGzipFile)):
                fh = get_input_handle(str(bad))
                try:
                    list(parse_fastq(fh))
                finally:
                    fh.close()

    def test_encode_matches_across_paths(self, monkeypatch):
        """End to end: the same .gz must encode to the same records either way."""
        from zna.cli import encode_command

        def run(path, out):
            class Args:
                files = [path]
                output = out
                interleaved = False
                read_group = "gz"
                description = ""
                seq_len_bytes = 2
                strand_specific = False
                strand_normalize = False
                npolicy = "trim3"
                compress_flag = False
                level = 3
                block_size = "64K"
                quiet = True
                fasta = False
                fastq = True
                shuffle = False
            encode_command(Args())
            with open(out, "rb") as f:
                return [r[0] for r in ZnaReader(f).records()]

        with tempfile.TemporaryDirectory() as tmp:
            p = str(Path(tmp) / "in.fastq.gz")
            seqs = self._make_gz(p)
            a = run(p, str(Path(tmp) / "a.zna"))
            monkeypatch.setenv("ZNA_NO_EXTERNAL_GZIP", "1")
            b = run(p, str(Path(tmp) / "b.zna"))
            assert a == seqs
            assert a == b


class TestEncodeShuffleBufferSize:
    """`--shuffle-buffer-size` was parsed and then ignored: `encode --shuffle`
    hardcoded a 1 GiB bucket budget, so the documented flag did nothing."""

    def _args(self, fq_path, out_path, buffer_size):
        # NB: the parameters must not be named `fastq`/`output`; those are also
        # class attributes below, which makes the class body's LOAD_NAME skip the
        # enclosing function scope.
        class Args:
            files = [fq_path]
            output = out_path
            interleaved = False
            read_group = "rg"
            description = ""
            seq_len_bytes = 2
            strand_specific = False
            strand_normalize = False
            npolicy = "trim3"
            compress_flag = False
            level = 3
            block_size = "64K"
            quiet = True
            fasta = False
            fastq = True
            shuffle = True
            seed = 7
            shuffle_buffer_size = buffer_size
        return Args()

    def _write_fastq(self, path, n=200):
        rng = random.Random(1234)
        with open(path, "w") as fh:
            for i in range(n):
                s = "".join(rng.choices("ACGT", k=150))
                fh.write(f"@r{i}\n{s}\n+\n{'I' * 150}\n")

    def test_flag_reaches_shuffle_zna(self, monkeypatch):
        from zna import cli

        seen = {}
        real = cli.shuffle_zna

        def spy(*a, **kw):
            seen.update(kw)
            return real(*a, **kw)

        monkeypatch.setattr(cli, "shuffle_zna", spy)
        with tempfile.TemporaryDirectory() as tmp:
            fq = str(Path(tmp) / "in.fastq")
            self._write_fastq(fq)
            cli.encode_command(
                self._args(fq, str(Path(tmp) / "a.zna"), "8M"))
        assert seen.get("buffer_bytes") == 8 * 1024 * 1024, (
            f"--shuffle-buffer-size did not reach shuffle_zna: {seen!r}")

    def test_default_is_still_one_gib(self, monkeypatch):
        """The default must parse to the value that was previously hardcoded,
        so no existing invocation changes behaviour."""
        from zna import cli

        seen = {}
        real = cli.shuffle_zna

        def spy(*a, **kw):
            seen.update(kw)
            return real(*a, **kw)

        monkeypatch.setattr(cli, "shuffle_zna", spy)
        with tempfile.TemporaryDirectory() as tmp:
            fq = str(Path(tmp) / "in.fastq")
            self._write_fastq(fq)
            args = self._args(fq, str(Path(tmp) / "b.zna"), None)
            cli.encode_command(args)
        assert seen.get("buffer_bytes") == 1 << 30

    def test_small_buffer_still_produces_a_valid_shuffled_file(self):
        """A tiny budget means many buckets; the output must still round-trip."""
        from zna import cli

        with tempfile.TemporaryDirectory() as tmp:
            fq = str(Path(tmp) / "in.fastq")
            self._write_fastq(fq, n=400)
            out = str(Path(tmp) / "c.zna")
            cli.encode_command(self._args(fq, out, "16K"))
            with open(out, "rb") as fh:
                got = [r[0] for r in ZnaReader(fh).records()]
            expected = []
            with open(fq) as fh:
                for i, line in enumerate(fh):
                    if i % 4 == 1:
                        expected.append(line.strip())
            assert sorted(got) == sorted(expected)


class TestFormatErrorsReachTheUser:
    """A file this build cannot read must produce a message, not a traceback.

    The version-3 break makes this the *first* thing anyone with older files sees, and
    it used to be a stack trace from `main()` on three of the four commands — `inspect`
    caught its own, `decode`, `encode` and `shuffle` did not.
    """

    #: argv tail per command; only `inspect` takes no output path.
    _ARGV = {
        "decode": lambda src, out: [src, "-o", out],
        "encode": lambda src, out: [src, "-o", out],
        "shuffle": lambda src, out: [src, "-o", out],
        "inspect": lambda src, out: [src],
    }

    def _v2_file(self, tmp) -> str:
        """A well-formed file with its version byte set back to 2."""
        path = str(Path(tmp) / "old.zna")
        buf = BytesIO()
        header = ZnaHeader(read_group="old", seq_len_bytes=2,
                           compression_method=COMPRESSION_ZSTD)
        with ZnaWriter(buf, header) as w:
            w.write_record("ACGTACGT", True, True, False)
            w.write_record("TTTTGGGG", True, False, True)
        data = bytearray(buf.getvalue())
        data[4] = 2                      # version byte
        Path(path).write_bytes(bytes(data))
        return path

    @pytest.mark.parametrize("command", ["decode", "encode", "inspect", "shuffle"])
    def test_an_unreadable_file_exits_cleanly(self, command, monkeypatch, capsys):
        from zna import cli

        with tempfile.TemporaryDirectory() as tmp:
            src = self._v2_file(tmp)
            out = str(Path(tmp) / "out.zna")
            argv = ["zna", command] + self._ARGV[command](src, out)
            monkeypatch.setattr("sys.argv", argv)
            with pytest.raises(SystemExit) as exc:
                cli.main()
            captured = capsys.readouterr()

        code = exc.value.code
        # argparse also exits non-zero; make sure we failed on the FILE, not the argv.
        assert code != 2, f"{command}: argparse rejected {argv!r}"
        assert code not in (0, None), f"{command} exited successfully"
        text = f"{code}{captured.out}{captured.err}"
        assert "version" in text.lower(), (
            f"{command} did not name the version as the problem: {text!r}"
        )
        assert "Traceback" not in text

    def test_a_readable_file_is_unaffected(self, monkeypatch):
        """The handler must not swallow a successful run."""
        from zna import cli

        with tempfile.TemporaryDirectory() as tmp:
            src = str(Path(tmp) / "good.zna")
            buf = BytesIO()
            with ZnaWriter(buf, ZnaHeader(read_group="x", seq_len_bytes=2)) as w:
                w.write_record("ACGTACGT", False, False, False)
            Path(src).write_bytes(buf.getvalue())
            monkeypatch.setattr("sys.argv", ["zna", "inspect", src])
            cli.main()          # must not raise
