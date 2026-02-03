import sys
import argparse
import gzip
import time
from pathlib import Path
from contextlib import ExitStack
from typing import BinaryIO, Iterator, Tuple, Optional, IO
import struct

# Import from core
from .core import (
    ZnaHeader, ZnaWriter, ZnaReader, 
    COMPRESSION_ZSTD, COMPRESSION_NONE, 
    DEFAULT_ZSTD_LEVEL, _BLOCK_HEADER_SIZE, _BLOCK_HEADER_FMT, _FILE_HEADER_SIZE,
    ZnaHeaderFlags, reverse_complement
)

# --- I/O HELPERS ---

def get_input_handle(filepath: Optional[str]) -> BinaryIO:
    """Returns a file handle for reading (supports gzip and stdin)."""
    if filepath is None or filepath == "-":
        return sys.stdin.buffer
    if filepath.endswith(".gz"):
        return gzip.open(filepath, "rb")
    return open(filepath, "rb")

def get_output_handle(filepath: Optional[str]) -> BinaryIO:
    """Returns a binary file handle for writing ZNA (supports stdout)."""
    if filepath is None or filepath == "-":
        return sys.stdout.buffer
    return open(filepath, "wb")

def open_text_output(filepath: str, compress: bool = False) -> IO[str]:
    """
    Opens a file for text writing (FASTA), optionally with gzip.
    """
    if compress or filepath.endswith(".gz"):
        return gzip.open(filepath, "wt")
    return open(filepath, "w")


# --- PARSERS ---

def get_base_name(full_name: str) -> str:
    """
    Extracts base read ID for pairing verification.
    Handles headers like: @ID/1 merged_... or @ID comment
    Returns 'ID' without /1 or /2 suffix and without comments.
    """
    # Split on whitespace to remove comments
    read_id = full_name.split()[0]
    # Split on slash to remove /1 or /2 pair indicators
    if "/" in read_id:
        return read_id.rsplit("/", 1)[0]
    return read_id


def get_read_suffix_number(full_name: str) -> int:
    """
    Returns 1 if name ends in /1, 2 if /2, else 0.
    Considers only the ID part before whitespace.
    """
    read_id = full_name.split()[0]
    if read_id.endswith("/1"):
        return 1
    if read_id.endswith("/2"):
        return 2
    return 0


def parse_fastq(fh: BinaryIO) -> Iterator[str]:
    """Yields sequence only from FASTQ stream."""
    while True:
        header = fh.readline()
        if not header: break
        if not header.startswith(b'@'): continue
        seq = fh.readline().strip().decode('ascii', errors='ignore')
        plus = fh.readline()  # +
        qual = fh.readline()  # Quality
        if seq:  # Only yield if we got a sequence
            yield seq


def parse_fastq_with_names(fh: BinaryIO) -> Iterator[Tuple[str, str]]:
    """Yields (read_name, sequence) tuples from FASTQ stream."""
    while True:
        header = fh.readline()
        if not header: break
        if not header.startswith(b'@'): continue
        read_name = header[1:].strip().decode('ascii', errors='ignore')
        seq = fh.readline().strip().decode('ascii', errors='ignore')
        plus = fh.readline()  # +
        qual = fh.readline()  # Quality
        if seq:  # Only yield if we got a sequence
            yield read_name, seq

def parse_fasta(fh: BinaryIO) -> Iterator[str]:
    """Yields sequence only from FASTA stream."""
    seq_parts = []
    for line in fh:
        line = line.strip()
        if not line: continue
        if line.startswith(b">"):
            if seq_parts:
                yield b"".join(seq_parts).decode('ascii', errors='ignore')
            seq_parts = []
        else:
            seq_parts.append(line)
    if seq_parts:
        yield b"".join(seq_parts).decode('ascii', errors='ignore')


def choose_parser(filepath: Optional[str], format_override: Optional[str] = None):
    """Returns appropriate parser based on format flag or file extension.
    
    Args:
        filepath: Path to input file or None for stdin
        format_override: 'fasta', 'fastq', or None to infer from extension
    
    Returns:
        Parser function (parse_fasta or parse_fastq)
    
    Raises:
        ValueError: If format cannot be determined
    """
    # 1. Format explicitly specified via command line flag
    if format_override:
        if format_override == 'fasta':
            return parse_fasta
        elif format_override == 'fastq':
            return parse_fastq
        else:
            raise ValueError(f"Unknown format: {format_override}")
    
    # 2. Infer from file extension
    if filepath and filepath != "-":
        # Remove .gz extension if present to check underlying format
        lower = filepath.lower()
        if lower.endswith('.gz'):
            lower = lower[:-3]  # Remove .gz suffix
        
        # Check for FASTA extensions
        if lower.endswith(('.fasta', '.fa', '.fna')):
            return parse_fasta
        
        # Check for FASTQ extensions
        if lower.endswith(('.fastq', '.fq')):
            return parse_fastq
        
        # Extension doesn't match known formats
        print(f"[Warning] Cannot determine format from filename '{filepath}'. "
              f"Use --fasta or --fastq to specify format explicitly. Defaulting to FASTQ.",
              file=sys.stderr)
        return parse_fastq
    
    # 3. Reading from stdin without format specified
    print("[Warning] Reading from stdin without format specified (--fasta or --fastq). "
          "Defaulting to FASTQ.", file=sys.stderr)
    return parse_fastq


# --- INPUT STRATEGY ---

def is_zna_file(filepath: Optional[str]) -> bool:
    """Check if a file is in ZNA format by reading magic bytes."""
    if not filepath or filepath == "-":
        # Can't easily check stdin without consuming bytes
        return False
    try:
        with open(filepath, "rb") as f:
            magic = f.read(4)
            from .core import _MAGIC
            return magic == _MAGIC
    except (IOError, OSError):
        return False


def stream_inputs(args) -> Iterator[Tuple[str, bool, bool, bool]]:
    """
    Uniform generator yielding (sequence, is_paired, is_read1, is_read2).
    Handles Files and Stdin (Single or Interleaved).
    
    Input modes:
    - 0 files: read from stdin (single or interleaved)
    - 1 file: read from file (single or interleaved, or ZNA for reencoding)
    - 2 files: paired-end (read1, read2)
    """
    # Determine format override from command line flags
    format_override = None
    if hasattr(args, 'fasta') and args.fasta:
        format_override = 'fasta'
    elif hasattr(args, 'fastq') and args.fastq:
        format_override = 'fastq'
    
    files = args.files if args.files else []
    
    # Special case: single ZNA file = reencoding mode
    if len(files) == 1 and is_zna_file(files[0]):
        with open(files[0], "rb") as f:
            reader = ZnaReader(f)
            for record in reader.records():
                yield record
        return
    
    with ExitStack() as stack:
        # 1. Two Files = Paired-End
        if len(files) == 2:
            f1 = stack.enter_context(get_input_handle(files[0]))
            f2 = stack.enter_context(get_input_handle(files[1]))
            p1 = choose_parser(files[0], format_override)(f1)
            p2 = choose_parser(files[1], format_override)(f2)
            
            for s1, s2 in zip(p1, p2):
                yield s1, True, True, False
                yield s2, True, False, True

        # 2. One File or Stdin + Interleaved Flag = Interleaved (Mixed Paired/Single-End)
        elif args.interleaved:
            src = files[0] if len(files) == 1 else None
            f = stack.enter_context(get_input_handle(src))
            
            # For interleaved mode, we need read names to detect pairing
            # Currently only support FASTQ with names for smart interleaved mode
            format_type = format_override
            if not format_type:
                # Infer format
                if src and src != "-":
                    lower = src.lower()
                    if lower.endswith('.gz'):
                        lower = lower[:-3]
                    if lower.endswith(('.fasta', '.fa', '.fna')):
                        format_type = 'fasta'
                    elif lower.endswith(('.fastq', '.fq')):
                        format_type = 'fastq'
                    else:
                        format_type = 'fastq'  # default
                else:
                    format_type = 'fastq'  # default for stdin
            
            if format_type == 'fastq':
                # Use smart interleaved mode with read name tracking
                parser = parse_fastq_with_names(f)
                prev_entry = None
                
                for curr_name, curr_seq in parser:
                    if prev_entry is None:
                        prev_entry = (curr_name, curr_seq)
                        continue
                    
                    prev_name, prev_seq = prev_entry
                    
                    # Check if consecutive reads belong to the same pair
                    if get_base_name(prev_name) == get_base_name(curr_name):
                        # Found a pair (R1 and R2)
                        n1 = get_read_suffix_number(prev_name)
                        n2 = get_read_suffix_number(curr_name)
                        
                        # Validation
                        if n1 == 2:
                            print(f"[Warning] Found Read 2 before Read 1: {prev_name} -> {curr_name}", file=sys.stderr)
                        if n2 == 1:
                            print(f"[Warning] Found Read 1 after Read 1: {prev_name} -> {curr_name}", file=sys.stderr)
                        
                        # Yield as paired reads
                        yield prev_seq, True, True, False   # R1
                        yield curr_seq, True, False, True   # R2
                        prev_entry = None
                    else:
                        # prev_entry was a singleton (single-end read)
                        yield prev_seq, False, False, False
                        prev_entry = (curr_name, curr_seq)
                
                # Handle the final remaining entry
                if prev_entry is not None:
                    yield prev_entry[1], False, False, False
            
            else:
                # FASTA format: fallback to simple interleaved pairing
                parser = parse_fasta(f)
                while True:
                    try:
                        s1 = next(parser)
                        try:
                            s2 = next(parser)
                        except StopIteration:
                            print("[Warning] Interleaved input ended with orphan read.", file=sys.stderr)
                            break
                        yield s1, True, True, False
                        yield s2, True, False, True
                    except StopIteration:
                        break

        # 3. One File or Stdin = Single-End
        else:
            src = files[0] if len(files) == 1 else None
            f = stack.enter_context(get_input_handle(src))
            for s in choose_parser(src, format_override)(f):
                yield s, False, False, False


# --- COMMAND: ENCODE ---

def encode_command(args):
    start_time = time.time()
    
    # Validation
    files = args.files if args.files else []
    if len(files) > 2:
        sys.exit("Error: Maximum 2 input files allowed")
    if len(files) == 2 and args.interleaved:
        sys.exit("Error: Cannot use --interleaved with 2 input files (already paired-end)")
    
    # Check if input is ZNA (reencoding mode)
    input_header = None
    is_reencoding = len(files) == 1 and is_zna_file(files[0])
    
    if is_reencoding:
        # Read existing header to use as defaults
        with open(files[0], "rb") as f:
            reader = ZnaReader(f)
            input_header = reader.header
        
        if not args.quiet:
            print(f"[ZNA] Reencoding {files[0]}...", file=sys.stderr)
    
    # Determine Output Mode
    # If -o is missing or "-", writing to stdout
    is_stdout = (args.output is None or args.output == "-")

    # Compression Logic
    should_compress = True 
    if args.compress_flag is not None:
        should_compress = args.compress_flag
    elif not is_stdout:
        # Auto-detect from filename if available
        ext = Path(args.output).suffix.lower()
        if ext == '.zna': should_compress = False
    
    comp_method = COMPRESSION_ZSTD if should_compress else COMPRESSION_NONE
    comp_level = args.level

    # Handle strand-specific flags
    strand_specific_flag = args.strand_specific
    
    # Determine R1/R2 antisense settings
    # Default for strand-specific: R1 antisense, R2 sense (dUTP protocol)
    if strand_specific_flag:
        # Check explicit flags, default to dUTP (R1 antisense, R2 sense)
        read1_antisense = not getattr(args, 'read1_sense', False)  # Default: antisense
        read2_antisense = getattr(args, 'read2_antisense', False)  # Default: sense
    else:
        read1_antisense = False
        read2_antisense = False

    # Build header: use input header as defaults if reencoding, otherwise use CLI args
    if is_reencoding and input_header:
        # Start with existing header values
        read_group = args.read_group if args.read_group != "Unknown" else input_header.read_group
        description = args.description if args.description != "" else input_header.description
        seq_len_bytes = args.seq_len_bytes if args.seq_len_bytes != 2 else input_header.seq_len_bytes
        strand_specific = strand_specific_flag if strand_specific_flag else input_header.strand_specific
        # Inherit antisense settings from input header if not explicitly overridden
        if not strand_specific_flag:
            read1_antisense = input_header.read1_antisense
            read2_antisense = input_header.read2_antisense
        
        # Always use new compression settings (that's usually why you reencode)
        # unless they match defaults, then keep original
        if args.compress_flag is None and comp_level == DEFAULT_ZSTD_LEVEL:
            comp_method = input_header.compression_method
            comp_level = input_header.compression_level
    else:
        # New encoding: use CLI args
        read_group = args.read_group
        description = args.description
        seq_len_bytes = args.seq_len_bytes
        strand_specific = strand_specific_flag

    # Header Setup
    header = ZnaHeader(
        read_group=read_group,
        description=description,
        seq_len_bytes=seq_len_bytes,
        strand_specific=strand_specific,
        read1_antisense=read1_antisense,
        read2_antisense=read2_antisense,
        compression_method=comp_method,
        compression_level=comp_level
    )

    if not is_stdout:
        c_str = f"ZSTD (L{comp_level})" if should_compress else "None"
        if not is_reencoding:
            print(f"[ZNA] Encoding to {args.output} [{c_str}]...", file=sys.stderr)
        else:
            print(f"[ZNA] Reencoding to {args.output} [{c_str}]...", file=sys.stderr)

    # Encode Loop
    count = 0
    quiet = hasattr(args, 'quiet') and args.quiet
    # Use ExitStack to safely close output file (or leave stdout open)
    npolicy = getattr(args, 'npolicy', None)
    with ExitStack() as stack:
        f_out = stack.enter_context(get_output_handle(args.output))
        writer = stack.enter_context(ZnaWriter(f_out, header, block_size=args.block_size, npolicy=npolicy))
        
        for seq, is_paired, is_r1, is_r2 in stream_inputs(args):
            # Skip sequence if it contains N and policy is 'drop'
            if npolicy == 'drop' and 'N' in seq.upper():
                continue
            writer.write_record(seq, is_paired, is_r1, is_r2)
            count += 1
            
            if count % 1_000_000 == 0 and not is_stdout and not quiet:
                print(f"      Processed {count//1_000_000}M records...", end='\r', file=sys.stderr)

    duration = time.time() - start_time
    if not is_stdout and not quiet:
        print(f"\n[ZNA] Done. Wrote {count} records in {duration:.2f}s.", file=sys.stderr)


# --- COMMAND: DECODE ---

def decode_command(args):
    # Determine Input
    input_file = args.input if args.input else None
    f_in_handle = get_input_handle(input_file)
    
    # Determine Output Config
    mode = "interleaved" # default
    outs = [] # List of file handles
    
    try:
        reader = ZnaReader(f_in_handle)
        rg = reader.header.read_group
        
        if not args.quiet:
            src_name = input_file if input_file else "stdin"
            print(f"[ZNA] Decoding {src_name} (RG: {rg})...", file=sys.stderr)

        with ExitStack() as stack:
            # --- 1. File Output Mode ---
            if args.output and args.output != "-":
                if "#" in args.output:
                    # Split Mode (e.g. "out#.fa.gz" -> "out_1.fa.gz", "out_2.fa.gz")
                    mode = "split"
                    fn1 = args.output.replace("#", "_1")
                    fn2 = args.output.replace("#", "_2")
                    
                    # Compression inferred from extension
                    f1 = stack.enter_context(open_text_output(fn1))
                    f2 = stack.enter_context(open_text_output(fn2))
                    outs = [f1, f2]
                else:
                    # Interleaved File Mode
                    mode = "interleaved"
                    f = stack.enter_context(open_text_output(args.output))
                    outs = [f]

            # --- 2. Stdout Mode ---
            else:
                mode = "interleaved"
                if args.gzip:
                    # Compressed Stdout
                    f = stack.enter_context(gzip.open(sys.stdout.buffer, "wt"))
                    outs = [f]
                else:
                    # Standard Stdout
                    outs = [sys.stdout]

            # --- Decode Loop ---
            counter = 0
            restore_strand = getattr(args, 'restore_strand', False)

            for seq, is_paired, is_r1, is_r2 in reader.records(restore_strand=restore_strand):
                
                # --- Strict ID Logic ---
                # Requirement: Read 1 (or Unpaired) always precedes Read 2.
                # Increment ID only on new template start.
                if is_r1 or (not is_paired):
                    counter += 1
                
                # --- Formatting ---
                suffix = ""
                if is_r1: suffix = "/1"
                elif is_r2: suffix = "/2"
                
                header_str = f">{rg}:{counter}{suffix}"
                
                if mode == "split":
                    if is_r2:
                        outs[1].write(f"{header_str}\n{seq}\n")
                    else:
                        # R1 or Unpaired goes to File 1
                        outs[0].write(f"{header_str}\n{seq}\n")
                else:
                    # Interleaved (File or Stdout)
                    outs[0].write(f"{header_str}\n{seq}\n")

    except BrokenPipeError:
        sys.stderr.close()
    finally:
        if input_file and input_file != "-":
            f_in_handle.close()


# --- COMMAND: INSPECT ---

def inspect_command(args):
    f_path = Path(args.input)
    if not f_path.exists():
        sys.exit(f"Error: File {args.input} not found.")

    file_size = f_path.stat().st_size
    print(f"File: {args.input}")
    print(f"Total Size: {file_size / (1024*1024):.2f} MB")
    
    with open(args.input, "rb") as f:
        try:
            reader = ZnaReader(f)
            h = reader.header
        except Exception as e:
            sys.exit(f"Error reading header: {e}")

        print("\n--- Header Metadata ---")
        print(f"Read Group:       {h.read_group}")
        print(f"Description:      {h.description}")
        print(f"Seq Length:       {h.seq_len_bytes} bytes (Max: {(1<<(8*h.seq_len_bytes))-1} bp)")
        print(f"Strand Specific:  {h.strand_specific}")
        if h.strand_specific:
            print(f"R1 Antisense:     {h.read1_antisense}")
            print(f"R2 Antisense:     {h.read2_antisense}")
        
        method = "None"
        if h.compression_method == COMPRESSION_ZSTD:
            method = f"ZSTD (Level {h.compression_level})"
        print(f"Compression:      {method}")

        # Scan Blocks
        block_count = 0
        total_records = 0
        compressed_payload = 0
        uncompressed_payload = 0
        
        # Seek past header
        f.seek(_FILE_HEADER_SIZE + len(h.read_group) + len(h.description))
        
        while True:
            b_header = f.read(_BLOCK_HEADER_SIZE)
            if not b_header: break
            
            c_size, u_size, n_recs = struct.unpack(_BLOCK_HEADER_FMT, b_header)
            
            block_count += 1
            total_records += n_recs
            compressed_payload += c_size
            uncompressed_payload += u_size
            
            f.seek(c_size, 1) # Skip payload

        print("\n--- Content Statistics ---")
        print(f"Total Blocks:       {block_count}")
        print(f"Total Records:      {total_records}")
        
        if compressed_payload > 0:
             print(f"Compressed Payload: {compressed_payload / (1024*1024):.2f} MB")
             print(f"Uncompressed Data:  {uncompressed_payload / (1024*1024):.2f} MB")
             ratio = uncompressed_payload / compressed_payload
             print(f"Compression Ratio:  {ratio:.2f}x")
        else:
             print(f"Data Size:          {uncompressed_payload / (1024*1024):.2f} MB")
             print(f"Compression Ratio:  1.00x")


# --- MAIN ---

def main():
    parser = argparse.ArgumentParser(description="zna: ZNA Processing Toolkit")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # --- ENCODE ---
    enc = subparsers.add_parser("encode", help="Convert FASTQ/FASTA to ZNA")
    enc.add_argument("files", nargs="*", help="Input files (0=stdin, 1=single/interleaved, 2=paired R1 R2)")
    
    input_group = enc.add_argument_group("Input Options")
    input_group.add_argument("--interleaved", action="store_true", help="Treat input as interleaved paired-end")
    input_group.add_argument("-q", "--quiet", action="store_true", help="Suppress progress messages")
    
    format_group = input_group.add_mutually_exclusive_group()
    format_group.add_argument("--fasta", action="store_true", help="Force FASTA format (overrides extension detection)")
    format_group.add_argument("--fastq", action="store_true", help="Force FASTQ format (overrides extension detection)")
    
    meta_group = enc.add_argument_group("Metadata")
    meta_group.add_argument("--read-group", default="Unknown", help="Read Group ID")
    meta_group.add_argument("--description", default="", help="Description string")
    meta_group.add_argument("--strand-specific", action="store_true", 
                           help="Flag library as strand-specific (default: R1 antisense, R2 sense)")
    meta_group.add_argument("--npolicy", choices=["drop", "random", "A", "C", "G", "T"],
                           help="Policy for handling 'N' nucleotides: drop (skip sequences), random (replace with random base), or A/C/G/T (replace with specific base)")
    
    # R1 strand orientation (mutually exclusive)
    r1_strand = meta_group.add_mutually_exclusive_group()
    r1_strand.add_argument("--read1-sense", dest="read1_sense", action="store_true",
                          help="Read 1 represents sense strand")
    r1_strand.add_argument("--read1-antisense", dest="read1_sense", action="store_false",
                          help="Read 1 represents antisense strand (default for --strand-specific)")
    
    # R2 strand orientation (mutually exclusive)
    r2_strand = meta_group.add_mutually_exclusive_group()
    r2_strand.add_argument("--read2-sense", dest="read2_antisense", action="store_false",
                          help="Read 2 represents sense strand (default for --strand-specific)")
    r2_strand.add_argument("--read2-antisense", dest="read2_antisense", action="store_true",
                          help="Read 2 represents antisense strand")
    
    fmt_group = enc.add_argument_group("Format Options")
    fmt_group.add_argument("-o", "--output", help="Output file. Defaults to stdout (-).")
    fmt_group.add_argument("--seq-len-bytes", type=int, choices=[1, 2, 4], default=2, help="Bytes for seq len")
    fmt_group.add_argument("--block-size", type=int, default=131072, help="Block size")
    
    comp_group = fmt_group.add_mutually_exclusive_group()
    comp_group.add_argument("--zstd", dest="compress_flag", action="store_true", default=None, help="Force ZSTD")
    comp_group.add_argument("--uncompressed", dest="compress_flag", action="store_false", help="Force uncompressed")
    fmt_group.add_argument("--level", type=int, default=DEFAULT_ZSTD_LEVEL, help="ZSTD level")

    # --- DECODE ---
    dec = subparsers.add_parser("decode", help="Convert ZNA to FASTA")
    dec.add_argument("input", nargs="?", help="Input .zna file (default: stdin)")
    dec.add_argument("-q", "--quiet", action="store_true", help="Suppress logs")
    
    # Combined Output Logic
    dec.add_argument("-o", "--output", help="Output filename. Use '#' for split files (e.g. 'out#.fa'). Ends in .gz for gzip.")
    dec.add_argument("--gzip", action="store_true", help="Force gzip compression (useful for stdout)")
    dec.add_argument("--restore-strand", action="store_true", 
                    help="Restore original strand orientation for antisense reads (reverse complement)")

    # --- INSPECT ---
    insp = subparsers.add_parser("inspect", help="Show ZNA file statistics")
    insp.add_argument("input", help="Input .zna/.zzna file")

    args = parser.parse_args()
    
    if args.command == "encode":
        encode_command(args)
    elif args.command == "decode":
        decode_command(args)
    elif args.command == "inspect":
        inspect_command(args)

if __name__ == "__main__":
    main()