import sys
import argparse
import gzip
import statistics
import time
from pathlib import Path
from typing import BinaryIO, Iterator, Tuple, Optional, List
from .core import (
    ZnaHeader, ZnaWriter, ZnaReader, 
    COMPRESSION_ZSTD, COMPRESSION_NONE, 
    DEFAULT_ZSTD_LEVEL, _BLOCK_HEADER_SIZE, _BLOCK_HEADER_FMT, _MAGIC, _VERSION, _FILE_HEADER_SIZE
)
import struct

# --- Helper: Lightweight FASTQ/FASTA Parser ---

def open_smart(filepath: str) -> BinaryIO:
    """Opens file as binary, handling .gz extension automatically."""
    if filepath.endswith(".gz"):
        return gzip.open(filepath, "rb")
    return open(filepath, "rb")

def parse_fastq(fh: BinaryIO) -> Iterator[Tuple[str, str]]:
    """Yields (header, sequence) from FASTQ stream."""
    # 4-line format: @Header, Seq, +, Qual
    while True:
        line = fh.readline()
        if not line: break
        header = line.strip().decode('ascii', errors='ignore')
        if not header.startswith('@'): continue # Skip garbage
        
        seq = fh.readline().strip().decode('ascii', errors='ignore')
        fh.readline() # +
        fh.readline() # Qual
        yield header[1:], seq

def parse_fasta(fh: BinaryIO) -> Iterator[Tuple[str, str]]:
    """Yields (header, sequence) from FASTA stream."""
    header = None
    seq_parts = []
    
    for line in fh:
        line = line.strip().decode('ascii', errors='ignore')
        if not line: continue
        
        if line.startswith(">"):
            if header:
                yield header, "".join(seq_parts)
            header = line[1:]
            seq_parts = []
        else:
            seq_parts.append(line)
            
    if header:
        yield header, "".join(seq_parts)

def get_parser(filepath: str, format_type: str = 'auto'):
    if format_type == 'fastq' or (format_type == 'auto' and ('fastq' in filepath or 'fq' in filepath)):
        return parse_fastq
    return parse_fasta

# --- Command: ENCODE ---

def encode_command(args):
    start_time = time.time()
    
    # 1. Input Validation
    if args.read2 and not args.read1:
        sys.exit("Error: --read2 requires --read1")
        
    inputs = []
    if args.read1:
        inputs.append(args.read1)
        if args.read2:
            inputs.append(args.read2)
    elif args.interleaved:
        inputs.append(args.interleaved)
    else:
        # Default to stdin? Or require input? Let's require input for clarity
        sys.exit("Error: Must provide input via --read1/--read2 or --interleaved")

    # 2. Setup Header
    comp_method = COMPRESSION_ZSTD if args.compress else COMPRESSION_NONE
    header = ZnaHeader(
        read_group=args.read_group,
        description=args.description,
        seq_len_bytes=args.seq_len_bytes,
        strand_specific=args.strand_specific,
        compression_method=comp_method,
        compression_level=args.level
    )

    print(f"[ZNA] Encoding...", file=sys.stderr)
    print(f"       Inputs: {inputs}", file=sys.stderr)
    print(f"       Compression: {'ZSTD (L' + str(args.level) + ')' if args.compress else 'None'}", file=sys.stderr)

    # 3. Processing Loop
    count = 0
    
    with open(args.output, "wb") as f_out:
        with ZnaWriter(f_out, header, block_size=args.block_size) as writer:
            
            # Case A: Paired Files (Read 1 + Read 2)
            if args.read1 and args.read2:
                with open_smart(args.read1) as f1, open_smart(args.read2) as f2:
                    p1 = get_parser(args.read1)(f1)
                    p2 = get_parser(args.read2)(f2)
                    
                    for (h1, s1), (h2, s2) in zip(p1, p2):
                        writer.write_record(s1, is_paired=True, is_read1=True, is_read2=False)
                        writer.write_record(s2, is_paired=True, is_read1=False, is_read2=True)
                        count += 1
            
            # Case B: Interleaved (Read 1, Read 2, Read 1, Read 2...)
            elif args.interleaved:
                with open_smart(args.interleaved) as f:
                    parser = get_parser(args.interleaved)(f)
                    try:
                        while True:
                            h1, s1 = next(parser)
                            h2, s2 = next(parser)
                            writer.write_record(s1, is_paired=True, is_read1=True, is_read2=False)
                            writer.write_record(s2, is_paired=True, is_read1=False, is_read2=True)
                            count += 1
                    except StopIteration:
                        pass

            # Case C: Single End (Read 1 only)
            else:
                with open_smart(args.read1) as f:
                    for h, s in get_parser(args.read1)(f):
                        writer.write_record(s, is_paired=False, is_read1=False, is_read2=False)
                        count += 1

    duration = time.time() - start_time
    print(f"[ZNA] Done. Wrote {count} records (or pairs) in {duration:.2f}s.", file=sys.stderr)


# --- Command: DECODE ---

def decode_command(args):
    # Output is always FASTA for simplicity, or raw sequences if requested
    out_fh = sys.stdout
    
    with open(args.input, "rb") as f_in:
        reader = ZnaReader(f_in)
        
        if not args.quiet:
            print(f"[ZNA] Decoding {args.input} (RG: {reader.header.read_group})...", file=sys.stderr)

        counter = 1
        for seq, is_paired, is_read1, is_read2 in reader.records():
            # Construct a FASTA-like header
            flags = []
            if is_paired: flags.append("paired")
            if is_read1: flags.append("read1")
            if is_read2: flags.append("read2")
            flag_str = ",".join(flags) if flags else "single"
            
            # Write FASTA
            out_fh.write(f">record_{counter} {flag_str}\n{seq}\n")
            counter += 1


# --- Command: INSPECT ---

def inspect_command(args):
    f_path = Path(args.input)
    file_size = f_path.stat().st_size
    
    print(f"File: {args.input}")
    print(f"Total Size: {file_size / (1024*1024):.2f} MB")
    
    with open(args.input, "rb") as f:
        # 1. Read File Header
        try:
            reader = ZnaReader(f)
            h = reader.header
        except Exception as e:
            sys.exit(f"Error reading header: {e}")

        print("\n--- Header Metadata ---")
        print(f"Read Group:      {h.read_group}")
        print(f"Description:     {h.description}")
        print(f"Seq Length Bytes:{h.seq_len_bytes}")
        print(f"Strand Specific: {h.strand_specific}")
        print(f"Compression:     {'ZSTD' if h.compression_method == COMPRESSION_ZSTD else 'None'} (Level {h.compression_level})")

        # 2. Scan Blocks for Statistics
        block_count = 0
        total_records = 0
        compressed_payload = 0
        uncompressed_payload = 0
        seq_lengths = []
        
        # We need to manually iterate blocks to get stats without decoding everything
        # We start AFTER the file header
        f.seek(_FILE_HEADER_SIZE + len(h.read_group) + len(h.description) + len(h.extra_info))
        
        try:
            while True:
                # Read Block Header
                b_header = f.read(_BLOCK_HEADER_SIZE)
                if not b_header: break
                
                c_size, u_size, n_recs = struct.unpack(_BLOCK_HEADER_FMT, b_header)
                
                block_count += 1
                total_records += n_recs
                compressed_payload += c_size
                uncompressed_payload += u_size
                
                # Skip the payload data
                f.seek(c_size, 1)
                
                # Sampling logic for lengths (optional, might require decoding blocks to be accurate)
                # To be fast, inspect usually just reads block headers. 
                # If we want exact seq length stats, we have to decode.
                # Let's skip deep decoding for speed unless requested.
                
        except struct.error:
            print("\n[Warning] File appears truncated.")

        print("\n--- Content Statistics ---")
        print(f"Total Blocks:       {block_count}")
        print(f"Total Records:      {total_records}")
        print(f"Compressed Data:    {compressed_payload / (1024*1024):.2f} MB")
        print(f"Uncompressed Data:  {uncompressed_payload / (1024*1024):.2f} MB")
        
        ratio = uncompressed_payload / compressed_payload if compressed_payload > 0 else 0
        print(f"Compression Ratio:  {ratio:.2f}x")
        
        overhead = file_size - compressed_payload
        print(f"Format Overhead:    {overhead / (1024*1024):.2f} MB ({(overhead/file_size)*100:.1f}%)")


# --- MAIN ---

def main():
    parser = argparse.ArgumentParser(description="zna: ZNA Processing Toolkit")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # --- ENCODE ---
    enc = subparsers.add_parser("encode", help="Convert FASTQ/FASTA to ZNA")
    input_group = enc.add_argument_group("Input Options")
    input_group.add_argument("--read1", help="Read 1 file (FASTA/FASTQ, can be .gz)")
    input_group.add_argument("--read2", help="Read 2 file (FASTA/FASTQ, can be .gz)")
    input_group.add_argument("--interleaved", help="Interleaved Paired-End file")
    
    meta_group = enc.add_argument_group("Metadata")
    meta_group.add_argument("--read-group", default="Unknown", help="Read Group ID")
    meta_group.add_argument("--description", default="", help="Description string")
    meta_group.add_argument("--strand-specific", action="store_true", help="Flag library as strand-specific")
    
    fmt_group = enc.add_argument_group("Format Options")
    fmt_group.add_argument("-o", "--output", required=True, help="Output .zna file")
    fmt_group.add_argument("--seq-len-bytes", type=int, choices=[1, 2, 4], default=2, help="Bytes for seq len (1=255bp, 2=65kb)")
    fmt_group.add_argument("--no-compress", dest="compress", action="store_false", help="Disable ZSTD compression")
    fmt_group.add_argument("--level", type=int, default=DEFAULT_ZSTD_LEVEL, help="ZSTD compression level (1-22)")
    fmt_group.add_argument("--block-size", type=int, default=131072, help="Internal block size (bytes)")
    enc.set_defaults(compress=True)

    # --- DECODE ---
    dec = subparsers.add_parser("decode", help="Convert ZNA to FASTA")
    dec.add_argument("-i", "--input", required=True, help="Input .zna file")
    dec.add_argument("-q", "--quiet", action="store_true", help="Suppress stderr logs")

    # --- INSPECT ---
    insp = subparsers.add_parser("inspect", help="Show ZNA file statistics")
    insp.add_argument("input", help="Input .zna file")

    args = parser.parse_args()
    
    if args.command == "encode":
        encode_command(args)
    elif args.command == "decode":
        decode_command(args)
    elif args.command == "inspect":
        inspect_command(args)

if __name__ == "__main__":
    main()