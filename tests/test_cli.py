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
    stream_inputs, encode_command, decode_command, inspect_command
)
from zna.core import (
    ZnaHeader, ZnaWriter, ZnaReader,
    COMPRESSION_ZSTD, COMPRESSION_NONE,
    _BLOCK_HEADER_FMT, _FILE_HEADER_SIZE
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
                    b_header = f.read(12)  # Block header size
                    if not b_header:
                        break
                    
                    c_size, u_size, n_recs = struct.unpack(_BLOCK_HEADER_FMT, b_header)
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
