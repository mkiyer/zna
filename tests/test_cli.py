"""
Unit tests for CLI functionality (encode, decode, inspect).
"""
import gzip
import tempfile
import struct
from pathlib import Path
from io import BytesIO

import pytest

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


# --- N-Policy CLI Tests ---

class TestNPolicyCLI:
    """Test CLI N-policy functionality."""
    
    def test_npolicy_drop(self):
        """Test that sequences with N are dropped when using --npolicy drop."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create FASTA with some sequences containing N
            fasta_path = f"{tmpdir}/input.fasta"
            with open(fasta_path, "w") as f:
                f.write(">seq1\nACGT\n")
                f.write(">seq2_with_N\nACNGT\n")
                f.write(">seq3\nTGCA\n")
                f.write(">seq4_with_N\nNNNN\n")
            
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
                npolicy = "drop"
                output = zna_path
                seq_len_bytes = 1
                block_size = 131072
                compress_flag = False
                level = 3
            
            encode_command(Args())
            
            # Decode and verify only sequences without N remain
            with open(zna_path, "rb") as f:
                reader = ZnaReader(f)
                records = list(reader.records())
                
                assert len(records) == 2  # Only seq1 and seq3
                assert records[0][0] == "ACGT"
                assert records[1][0] == "TGCA"
    
    def test_npolicy_replace_A(self):
        """Test that N nucleotides are replaced with A."""
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
                npolicy = "A"
                output = zna_path
                seq_len_bytes = 1
                block_size = 131072
                compress_flag = False
                level = 3
            
            encode_command(Args())
            
            # Decode and verify N was replaced with A
            with open(zna_path, "rb") as f:
                reader = ZnaReader(f)
                records = list(reader.records())
                
                assert len(records) == 1
                assert records[0][0] == "ACAGT"  # N -> A
    
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
        a.npolicy = "drop"
        a.shuffle = True
        a.seed = 42

        encode_command(a)

        out_records = self._read_zna(out_path)
        assert len(out_records) == 100

        # All sequences present
        assert sorted(seqs) == sorted(r[0] for r in out_records)

        # Order should differ from input
        assert seqs != [r[0] for r in out_records]
