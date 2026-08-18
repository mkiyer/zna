import sys
import argparse
import gzip
import io
import os
import subprocess
import tempfile
import time
from pathlib import Path
from contextlib import ExitStack
from typing import BinaryIO, Iterator, Tuple, Optional, IO
import struct

# Import from core
from .core import (
    ZnaHeader, ZnaWriter, ZnaReader, 
    COMPRESSION_ZSTD, COMPRESSION_NONE, 
    DEFAULT_ZSTD_LEVEL, DEFAULT_BLOCK_SIZE,
    _FILE_HEADER_FMT, _FILE_HEADER_SIZE,
    _BLOCK_HEADER_FMT, _BLOCK_HEADER_SIZE, _VERSION,
    ZnaHeaderFlags, ZnaRecordFlags, reverse_complement,
    FLAG_FIELDS,
)
from .dtypes import LabelDef, parse_dtype, label_bytes_per_record
from ._shuffle import shuffle_zna

# Try to import the C++ fast label extractor
_accel_extract = None
try:
    from ._accel import extract_labels_fast as _accel_extract
except (ImportError, AttributeError):
    pass


#: Read-ahead in front of the gzip *module*, whose own readline buffer is small.
_GZIP_READ_BUFFER = 1 << 20

#: Read-ahead in front of a decompressor pipe.  Measured flat from 32 KiB to
#: 1 MiB, so this just matches what subprocess would have chosen.
_PIPE_READ_BUFFER = io.DEFAULT_BUFFER_SIZE

#: External decompressors to try, fastest first.  ``pigz`` is multithreaded.
_GUNZIP_COMMANDS = (("pigz", "-dc"), ("gzip", "-dc"))


class _DecompressPipe(io.BufferedReader):
    """A read stream fed by an external decompressor.

    :meth:`close` reaps the child and checks its exit status.  That check is the
    whole reason this exists rather than handing back ``proc.stdout`` directly:
    a truncated or corrupt ``.gz`` makes the decompressor emit a partial stream
    and *then* exit non-zero.  Without the check that is indistinguishable from
    a short input file, and ``zna encode`` would write a silently truncated ZNA
    and report success.

    The check is not free.  CPython's buffered ``readline`` has a fast path that
    applies only to an exact ``BufferedReader``, so subclassing costs about 18%
    of the read (93 ms against 79 ms on 300k records).  Interposing a
    ``RawIOBase`` below an exact ``BufferedReader`` was measured too, on the
    theory that a per-buffer-fill ``readinto`` would be cheaper than a
    per-``readline`` slow path; it is worse, at 114 ms.  So: subclass, and pay
    18% of a 2.3x win for not silently truncating people's data.
    """

    def __init__(self, proc, path, buffer_size=None):
        super().__init__(proc.stdout, buffer_size or _PIPE_READ_BUFFER)
        self._proc = proc
        self._path = path

    def close(self) -> None:
        if self.closed:
            return
        try:
            super().close()
        finally:
            # Reap the child however the close above went, so a failure here can
            # never leave the process behind.
            returncode = self._proc.wait()
        # A negative code is a signal.  SIGPIPE means *we* stopped reading early
        # (a downstream `head`, a BrokenPipeError), not a decompression failure.
        if returncode not in (0, -13):
            raise OSError(
                f"Decompressing {self._path} failed: "
                f"{' '.join(self._proc.args)} exited with {returncode}. "
                f"The file may be truncated or corrupt."
            )


def _open_gzip(filepath: str):
    """Open a gzip file for reading, preferring an external decompressor.

    Handing decompression to a separate process lets it run in parallel with
    Python's parsing instead of competing with it for the GIL; measured ~2.75x
    over the ``gzip`` module on 300k records.  ``pigz`` is tried first because it
    threads the decompression too.

    Falls back to the ``gzip`` module when no external tool is present, and
    ``ZNA_NO_EXTERNAL_GZIP=1`` forces that path for debugging.
    """
    if not os.environ.get("ZNA_NO_EXTERNAL_GZIP"):
        for cmd in _GUNZIP_COMMANDS:
            try:
                proc = subprocess.Popen(
                    [*cmd, filepath], stdout=subprocess.PIPE, bufsize=0,
                )
            except OSError:
                # No such binary, or it is not executable — try the next one.
                continue
            return _DecompressPipe(proc, filepath)
    # A FASTQ record is four readline calls, and GzipFile serves each from its
    # own modest internal buffer; a wide BufferedReader in front amortises that
    # across many records instead of paying it per line.
    return io.BufferedReader(gzip.open(filepath, "rb"), buffer_size=_GZIP_READ_BUFFER)


def parse_block_size(value) -> int:
    """Parse a human-readable block size string (e.g. '512K', '4M', '8M').
    
    Accepts plain integers (bytes) or suffixed values:
        K/KB = kilobytes, M/MB = megabytes
    """
    if isinstance(value, int):
        return value
    value = str(value).strip().upper()
    multipliers = {'K': 1024, 'KB': 1024, 'M': 1024 * 1024, 'MB': 1024 * 1024,
                   'G': 1024 * 1024 * 1024, 'GB': 1024 * 1024 * 1024}
    for suffix, mult in sorted(multipliers.items(), key=lambda x: -len(x[0])):
        if value.endswith(suffix):
            try:
                return int(value[:-len(suffix)]) * mult
            except ValueError:
                raise argparse.ArgumentTypeError(
                    f"Invalid block size: '{value}'. Use integers with optional K/M suffix."
                )
    try:
        return int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"Invalid block size: '{value}'. Use integers with optional K/M suffix (e.g. 512K, 4M)."
        )

# --- I/O HELPERS ---

def get_input_handle(filepath: Optional[str]) -> BinaryIO:
    """Returns a file handle for reading (supports gzip and stdin)."""
    if filepath is None or filepath == "-":
        return sys.stdin.buffer
    if filepath.endswith(".gz"):
        return _open_gzip(filepath)
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


#: ``suffix_number`` for a record that is a merged read rather than a mate.  Negative so
#: it can never collide with a real mate number, and so ``if suffix`` still reads as
#: "this record names a mate position".
MERGED_SINGLE = -1

#: The token `zna merge` and fastp append to a merged read's name.  See :func:`_read_key`.
_MERGED_MARKER = b" merged_"


def _read_key(raw_header: bytes) -> Tuple[bytes, int]:
    """``(base_name, suffix_number)`` from a raw FASTQ header line, as bytes.

    The bytes-native equivalent of ``get_base_name`` + ``get_read_suffix_number``,
    computed in one pass and without decoding.  Those two helpers each ran
    ``full_name.split()[0]``, which tokenizes the *entire* header — comments and
    all — and they ran twice per record between them, on a string that had
    already been decoded from ASCII for no other reason.

    *base_name* is the ID up to the first whitespace with any ``/suffix``
    removed; *suffix_number* is 1 for ``/1``, 2 for ``/2``,
    :data:`MERGED_SINGLE` for a merged read, else 0.
    """
    sp = raw_header.find(b" ")
    tab = raw_header.find(b"\t")
    if sp < 0:
        end = tab
    elif tab < 0:
        end = sp
    else:
        end = sp if sp < tab else tab
    read_id = raw_header if end < 0 else raw_header[:end]

    slash = read_id.rfind(b"/")
    if slash >= 0:
        suffix = read_id[slash + 1:]
        if suffix == b"1":
            return read_id[:slash], 1
        if suffix == b"2":
            return read_id[:slash], 2
        read_id = read_id[:slash]

    # No mate number.  A merged read is a whole molecule and can never be a mate, so say
    # so explicitly rather than leaving pairing to guess from the name -- names collide.
    # Two independently merged molecules sharing a read name were silently paired into
    # one fragment: both lost IS_FULL_FRAGMENT, both were flagged IS_PAIRED/IS_READ1/
    # IS_READ2, and two different molecules were encoded as each other's mate, with a
    # zero exit status.  Duplicate names are ordinary in real data (concatenated lanes,
    # some SRA dumps, UMI-collapsed files) and nothing upstream forbids them.
    #
    # The trailing " merged_<n1>_<n2>" token is what marks one.  It is not a heuristic:
    # `zna merge` writes it (merge/merge_core.hpp), fastp writes it, and khorana's
    # `parse_merged_fastq` requires it -- it is already load-bearing in this contract,
    # so reading it here makes an existing promise checked rather than adding a new one.
    #
    # Tested only when there is no /1,/2 suffix, so it costs one `rfind` per merged read
    # and nothing per mate.
    if end >= 0 and raw_header.rfind(_MERGED_MARKER, end) >= 0:
        return read_id, MERGED_SINGLE
    return read_id, 0


def parse_fastq(fh: BinaryIO) -> Iterator[str]:
    """Yields sequence only from FASTQ stream.
    
    Optimized for minimal overhead: reads 4 lines at a time and
    only decodes the sequence line.
    """
    readline = fh.readline  # Cache method lookup
    while True:
        header = readline()
        if not header: 
            break
        # Skip non-record lines (shouldn't happen in valid FASTQ)
        if header[0] != 64:  # ord('@') = 64, faster than startswith
            continue
        seq_line = readline()
        readline()  # + line (discard)
        readline()  # Quality line (discard)
        if seq_line:
            yield seq_line.rstrip(b"\r\n").decode('ascii')


def parse_fastq_with_names(fh: BinaryIO) -> Iterator[Tuple[str, str]]:
    """Yields (read_name, sequence) tuples from FASTQ stream.
    
    Optimized for minimal overhead.
    """
    readline = fh.readline  # Cache method lookup
    while True:
        header = readline()
        if not header: 
            break
        if header[0] != 64:  # ord('@') = 64
            continue
        # Extract read name (skip @ and strip)
        read_name = header[1:].rstrip(b"\r\n").decode('ascii')
        
        seq_line = readline()
        readline()  # + line
        readline()  # Quality line
        
        if seq_line:
            yield read_name, seq_line.rstrip(b"\r\n").decode('ascii')


def parse_fastq_keyed(fh: BinaryIO) -> Iterator[Tuple[bytes, int, str]]:
    """Yields ``(base_name_bytes, suffix_number, sequence)`` from FASTQ.

    The pairing rule only ever compares base names for equality and reads the
    ``/1``,``/2`` suffix, so neither needs to be a ``str``.  Producing the key
    here — once, from the raw header bytes — keeps the interleaved path from
    decoding a name it never prints and re-deriving it per comparison.
    """
    readline = fh.readline
    while True:
        header = readline()
        if not header:
            break
        if header[0] != 64:  # ord('@')
            continue
        base, suffix = _read_key(header[1:].rstrip(b"\r\n"))

        seq_line = readline()
        readline()  # + line
        readline()  # quality line

        if seq_line:
            yield base, suffix, seq_line.rstrip(b"\r\n").decode('ascii')


def parse_fastq_with_headers(fh: BinaryIO) -> Iterator[Tuple[bytes, str]]:
    """Yields (raw_header_bytes, sequence) tuples from FASTQ stream.

    The header is the full line after ``@`` (as bytes, without the ``@`` or
    the trailing newline).  The caller is responsible for parsing SAM tags
    from the header.
    """
    readline = fh.readline
    while True:
        header = readline()
        if not header:
            break
        if header[0] != 64:  # ord('@')
            continue
        raw_header = header[1:].rstrip(b"\r\n")

        seq_line = readline()
        readline()  # +
        readline()  # quality

        if seq_line:
            yield raw_header, seq_line.rstrip(b"\r\n").decode('ascii')

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


# --- LABEL HELPERS ---

def parse_label_spec(spec: str) -> Tuple[str, str, Optional[str]]:
    """Parse a ``--label NAME:TYPE`` or ``--label NAME:TYPE:TAG`` spec string.

    Returns ``(name, dtype_code_or_name, tag_or_none)``.

    Examples::

        "NH:C"      -> ("NH", "C", None)
        "AS:i"      -> ("AS", "i", None)
        "edits:C:NM" -> ("edits", "C", "NM")
    """
    if ':' not in spec:
        raise argparse.ArgumentTypeError(
            f"Invalid --label spec {spec!r}: expected NAME:TYPE or NAME:TYPE:TAG (e.g. NH:C)"
        )
    parts = spec.split(':')
    if len(parts) == 2:
        name, type_str = parts
        tag: Optional[str] = None
    elif len(parts) == 3:
        name, type_str, tag = parts
        if not tag:
            raise argparse.ArgumentTypeError(f"Empty tag in --label {spec!r}")
    else:
        raise argparse.ArgumentTypeError(
            f"Invalid --label spec {spec!r}: expected NAME:TYPE or NAME:TYPE:TAG"
        )
    if not name:
        raise argparse.ArgumentTypeError(f"Empty label name in --label {spec!r}")

    # Validate dtype
    try:
        parse_dtype(type_str)
    except ValueError:
        raise argparse.ArgumentTypeError(f"Unknown dtype {type_str!r} in --label {spec!r}")

    return name, type_str, tag


def build_label_defs(
    label_specs: list[str],
    label_descs: list[str],
) -> tuple[LabelDef, ...]:
    """Build a tuple of :class:`LabelDef` from CLI ``--label`` and ``--label-desc`` args.

    *label_specs* are ``"NAME:TYPE"`` or ``"NAME:TYPE:TAG"`` strings.
    *label_descs* are ``"NAME:Description text"`` strings (optional,
    may be shorter than *label_specs*).
    """
    # Index descriptions by label name
    desc_map: dict[str, str] = {}
    for desc_spec in label_descs:
        if ':' in desc_spec:
            name, desc = desc_spec.split(':', 1)
            desc_map[name] = desc

    defs: list[LabelDef] = []
    for i, spec in enumerate(label_specs):
        name, type_str, tag = parse_label_spec(spec)
        dtype = parse_dtype(type_str)
        description = desc_map.get(name, "")
        defs.append(LabelDef(
            label_id=i, name=name, description=description,
            dtype=dtype, tag=tag,
        ))
    return tuple(defs)


def load_label_file(path: str) -> tuple[LabelDef, ...]:
    """Load label definitions from a YAML file.

    Each entry in the top-level ``labels`` list must have at least ``name``
    and ``type``.  Optional keys: ``description``, ``missing``.

    Example YAML::

        labels:
          - name: NM
            type: C
            description: Edit distance
          - name: AS
            type: i
            description: Alignment score
            missing: -1

    Returns a tuple of :class:`LabelDef` objects.  Raises
    ``SystemExit`` with a user-friendly message on errors.
    """
    # Imported here rather than at module scope on purpose: `import yaml` costs ~15 ms,
    # and `--label-defs` is the only thing in zna that wants it.  Every `zna decode`,
    # `zna inspect` and unlabeled `zna encode` would otherwise pay it.  Keep it local.
    try:
        import yaml
    except ImportError:
        sys.exit(
            "Error: PyYAML is missing.  It is a required dependency of zna, so this "
            "install is incomplete rather than merely missing an extra -- reinstall "
            "with `pip install zna` or `conda install -c bioconda zna`."
        )

    try:
        with open(path) as fh:
            data = yaml.safe_load(fh)
    except FileNotFoundError:
        sys.exit(f"Error: label file not found: {path}")
    except yaml.YAMLError as e:
        sys.exit(f"Error: invalid YAML in {path}: {e}")

    if not isinstance(data, dict) or "labels" not in data:
        sys.exit(
            f"Error: {path} must contain a top-level 'labels' key "
            f"with a list of label definitions."
        )

    entries = data["labels"]
    if not isinstance(entries, list):
        sys.exit(f"Error: 'labels' in {path} must be a list.")

    defs: list[LabelDef] = []
    for i, entry in enumerate(entries):
        if not isinstance(entry, dict):
            sys.exit(f"Error: label entry {i} in {path} must be a mapping.")

        name = entry.get("name")
        if not name:
            sys.exit(f"Error: label entry {i} in {path} is missing 'name'.")

        type_str = entry.get("type")
        if not type_str:
            sys.exit(f"Error: label '{name}' in {path} is missing 'type'.")

        type_str = str(type_str)
        try:
            dtype = parse_dtype(type_str)
        except ValueError:
            sys.exit(
                f"Error: unknown dtype {type_str!r} for label '{name}' in {path}.\n"
                f"Valid types: A, c, C, s, S, i, I, f, d, q, Q"
            )

        description = str(entry.get("description", ""))
        missing_raw = entry.get("missing")
        missing: int | float | None = None
        if missing_raw is not None:
            if isinstance(missing_raw, str) and len(missing_raw) == 1:
                # Single character → ord (for type A)
                missing = ord(missing_raw)
            else:
                try:
                    missing = int(missing_raw) if dtype.code not in ('f', 'd') else float(missing_raw)
                except (ValueError, TypeError):
                    sys.exit(
                        f"Error: invalid missing value {missing_raw!r} "
                        f"for label '{name}' in {path}."
                    )

        tag_str = entry.get("tag")
        if tag_str is not None:
            tag_str = str(tag_str)

        defs.append(LabelDef(
            label_id=i, name=name, description=description,
            dtype=dtype, missing=missing, tag=tag_str,
        ))

    return tuple(defs)


def merge_label_defs(
    yaml_defs: tuple[LabelDef, ...],
    cli_specs: list[str],
    cli_descs: list[str],
) -> tuple[LabelDef, ...]:
    """Merge CLI ``--label`` / ``--label-desc`` overrides onto YAML base defs.

    CLI flags override YAML values by tag name.  Any CLI ``--label`` that
    does not appear in the YAML file is appended at the end.
    """
    # Build override maps from CLI (keyed by label name)
    cli_type_map: dict[str, str] = {}
    cli_tag_map: dict[str, Optional[str]] = {}
    for spec in cli_specs:
        name, type_str, tag = parse_label_spec(spec)
        cli_type_map[name] = type_str
        cli_tag_map[name] = tag

    cli_desc_map: dict[str, str] = {}
    for desc_spec in cli_descs:
        if ':' in desc_spec:
            name, desc = desc_spec.split(':', 1)
            cli_desc_map[name] = desc

    # Pass 1: update existing YAML defs with CLI overrides
    seen_names: set[str] = set()
    merged: list[LabelDef] = []
    for ldef in yaml_defs:
        name = ldef.name
        seen_names.add(name)

        dtype = ldef.dtype
        if name in cli_type_map:
            dtype = parse_dtype(cli_type_map[name])

        desc = ldef.description
        if name in cli_desc_map:
            desc = cli_desc_map[name]

        tag = ldef.tag
        if name in cli_tag_map and cli_tag_map[name] is not None:
            tag = cli_tag_map[name]

        merged.append(LabelDef(
            label_id=len(merged), name=name, description=desc,
            dtype=dtype, missing=ldef.missing, tag=tag,
        ))

    # Pass 2: append CLI-only labels not in YAML
    for spec in cli_specs:
        name, type_str, tag = parse_label_spec(spec)
        if name not in seen_names:
            dtype = parse_dtype(type_str)
            desc = cli_desc_map.get(name, "")
            merged.append(LabelDef(
                label_id=len(merged), name=name, description=desc,
                dtype=dtype, tag=tag,
            ))
            seen_names.add(name)

    return tuple(merged)


# Conversion type codes for tag_map (avoids closure/lambda dispatch)
_CONV_INT = 0
_CONV_FLOAT = 1
_CONV_ORD = 2


def build_tag_extractor(label_defs: tuple[LabelDef, ...]) -> dict[bytes, Tuple[int, int, LabelDef]]:
    """Build a lookup dict for fast SAM tag extraction from FASTQ headers.

    Returns ``{b'NH': (label_index, conv_code, ldef), ...}`` where
    *conv_code* is one of ``_CONV_INT``, ``_CONV_FLOAT``, ``_CONV_ORD``.
    """
    tag_map: dict[bytes, Tuple[int, int, LabelDef]] = {}
    for i, ldef in enumerate(label_defs):
        tag_bytes = ldef.effective_tag.encode('ascii')
        if ldef.dtype.code == 'A':
            conv_code = _CONV_ORD
        elif ldef.dtype.code in ('f', 'd'):
            conv_code = _CONV_FLOAT
        else:
            conv_code = _CONV_INT
        tag_map[tag_bytes] = (i, conv_code, ldef)

    return tag_map


def extract_labels_from_header(
    raw_header: bytes,
    tag_map: dict[bytes, Tuple[int, int, object]],
    num_labels: int,
    label_defs: tuple[LabelDef, ...] | None = None,
) -> tuple:
    """Extract label values from a FASTQ header with key-value tags.

    The header format is::

        READNAME<ws>KEY:TYPE:VALUE<ws>KEY:TYPE:VALUE...

    Keys can be any length (e.g. standard 2-char SAM tags like ``NM``
    or longer custom keys like ``edit_distance``).  TYPE is a single
    character (ignored by the parser — the dtype from the label
    definition is used).

    Fields are split on **any whitespace** (tabs and spaces) so that
    tools like fastp that append space-separated tokens (e.g.
    ``merged_150_87``) do not corrupt the last tag value.

    Missing tags are filled with the label's ``missing`` value (from
    :func:`resolve_missing`).  If *label_defs* is not provided, a
    ``ValueError`` is raised for any absent tag.
    """
    from .dtypes import resolve_missing as _resolve_missing

    fields = raw_header.split()  # split on any whitespace
    values: list = [None] * num_labels
    remaining = num_labels

    for field in fields[1:]:  # skip read name
        # KEY:TYPE:VALUE — find first colon to extract key
        colon1 = field.find(b':')
        if colon1 < 1:
            continue
        # TYPE char after first colon, second colon separates value
        if len(field) < colon1 + 3 or field[colon1 + 2:colon1 + 3] != b':':
            continue
        tag = field[:colon1]
        entry = tag_map.get(tag)
        if entry is not None:
            idx, conv_code, _ldef = entry
            raw_val = field[colon1 + 3:]
            if conv_code == _CONV_INT:
                values[idx] = int(raw_val)
            elif conv_code == _CONV_FLOAT:
                values[idx] = float(raw_val)
            else:  # _CONV_ORD
                values[idx] = raw_val[0] if len(raw_val) == 1 else int(raw_val)
            remaining -= 1
            if remaining == 0:
                return tuple(values)

    # Fill in missing values from label definitions
    for i in range(num_labels):
        if values[i] is None:
            if label_defs is not None:
                values[i] = _resolve_missing(label_defs[i])
            else:
                raise ValueError(
                    f"Missing SAM tag for label index {i} in header: "
                    f"{raw_header[:80].decode('ascii', errors='replace')!r}"
                )

    return tuple(values)


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


# --- INPUT STRATEGIES ---
# Each strategy is a focused generator for a specific input mode.

def _stream_zna_reencode(
    filepath: str, with_ends: bool = False
) -> Iterator[Tuple[str, bool, bool, bool]]:
    """Stream records from an existing ZNA file for reencoding.

    With *with_ends*, each record carries ``has_start, has_end`` — the lossless
    form of ``IS_RC`` plus ``IS_FULL_FRAGMENT`` — so the writer can copy the
    existing orientation and fragment-span verbatim instead of re-deriving them.
    """
    with open(filepath, "rb") as f:
        reader = ZnaReader(f)
        for record in reader.records(with_ends=with_ends):
            yield record


def _stream_paired_files(f1: BinaryIO, f2: BinaryIO,
                         path1: Optional[str], path2: Optional[str],
                         format_override: Optional[str]) -> Iterator[Tuple[str, bool, bool, bool]]:
    """Stream paired-end reads from two separate files.

    Pairing is **positional** — record *i* of R1 with record *i* of R2 — which is what
    the two-file mode has always meant.  Read names are not compared; use ``zna merge``
    or ``--interleaved`` if you need that.

    The two files must therefore hold the same number of records, and this checks it.
    It used to be ``for s1, s2 in zip(p1, p2)``, which stops at the shorter input: an R2
    file with fewer records than R1 truncated the library silently and exited 0, every
    remaining R1 read gone with nothing in the output to say so.

    ``zip`` cannot be repaired by inspecting the iterators afterwards, either — it has
    already pulled and discarded the extra value from the longer one, so the very
    off-by-one case that matters most looks like a clean finish.  Pull both explicitly.
    """
    p1 = choose_parser(path1, format_override)(f1)
    p2 = choose_parser(path2, format_override)(f2)
    end = _END

    while True:
        s1 = next(p1, end)
        s2 = next(p2, end)
        if s1 is end:
            if s2 is end:
                return
            longer, shorter = path2, path1
            break
        if s2 is end:
            longer, shorter = path1, path2
            break
        yield s1, True, True, False
        yield s2, True, False, True

    sys.exit(
        f"Error: {_name_of(shorter)} ran out of records before {_name_of(longer)}. "
        f"Paired-end input must hold the same number of records in both files; "
        f"encoding what was read would silently drop the rest of the library."
    )


#: Sentinel for "this stream is exhausted", so a legitimately empty sequence cannot be
#: mistaken for the end of the file.
_END = object()


def _name_of(path: Optional[str]) -> str:
    return path if path else "<stdin>"


def _describe_read(base: bytes, suffix: int) -> str:
    """Render a keyed read name for a warning message."""
    name = base.decode('ascii', errors='replace')
    return f"{name}/{suffix}" if suffix else name


def _pair_interleaved(records) -> Iterator[Tuple[object, bool, bool, bool]]:
    """Assign pairing flags to an interleaved stream of
    ``(base_name, suffix_number, payload)``.

    Consecutive records whose base names match are emitted as an R1/R2 pair;
    a record whose neighbour has a different base name is a single.

    A record :func:`_read_key` identified as a merged read (:data:`MERGED_SINGLE`) is
    **never** paired, whatever its name.  Matching on the name alone was not enough: two
    independently merged molecules that share a read name were paired into one fragment
    and both silently lost ``IS_FULL_FRAGMENT``.

    Yields ``(payload, is_paired, is_read1, is_read2)``.  *payload* is opaque —
    the sequence for unlabeled input, or ``(sequence, labels)`` for labeled —
    so that labeled and unlabeled encoding share one implementation of the
    pairing rule instead of drifting apart.

    The base name and suffix arrive pre-computed (see :func:`_read_key`) because
    this comparison is the whole per-record cost of interleaved input: deriving
    them here meant tokenizing each header twice and decoding a name that is only
    ever compared, never printed.
    """
    prev_base = None
    prev_suffix = 0
    prev_payload = None
    have_prev = False

    for base, suffix, payload in records:
        if not have_prev:
            prev_base, prev_suffix, prev_payload = base, suffix, payload
            have_prev = True
            continue

        if (prev_base == base
                and prev_suffix != MERGED_SINGLE and suffix != MERGED_SINGLE):
            if prev_suffix == 2:
                print(f"[Warning] Found Read 2 before Read 1: "
                      f"{_describe_read(prev_base, prev_suffix)} -> "
                      f"{_describe_read(base, suffix)}", file=sys.stderr)
            if suffix == 1:
                print(f"[Warning] Found Read 1 after Read 1: "
                      f"{_describe_read(prev_base, prev_suffix)} -> "
                      f"{_describe_read(base, suffix)}", file=sys.stderr)

            yield prev_payload, True, True, False   # R1
            yield payload, True, False, True        # R2
            have_prev = False
        else:
            # prev was a singleton (single-end or merged read)
            yield prev_payload, False, False, False
            prev_base, prev_suffix, prev_payload = base, suffix, payload

    if have_prev:
        yield prev_payload, False, False, False


def _stream_interleaved_fastq(f: BinaryIO) -> Iterator[Tuple[str, bool, bool, bool]]:
    """Stream interleaved FASTQ with smart paired/single detection based on read names.

    Returns the pairing generator directly rather than re-yielding through it: for
    unlabeled input the payload *is* the sequence, so the shapes already match and
    a wrapper layer would cost a generator resume per record.
    """
    return _pair_interleaved(parse_fastq_keyed(f))


def _stream_interleaved_fasta(f: BinaryIO) -> Iterator[Tuple[str, bool, bool, bool]]:
    """Stream strict interleaved FASTA (alternating R1/R2 pairs)."""
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


def _stream_single_end(f: BinaryIO, filepath: Optional[str], 
                       format_override: Optional[str]) -> Iterator[Tuple[str, bool, bool, bool]]:
    """Stream single-end reads from a file or stdin."""
    for s in choose_parser(filepath, format_override)(f):
        yield s, False, False, False


def _infer_format(filepath: Optional[str], format_override: Optional[str]) -> str:
    """Infer format from filepath or return override."""
    if format_override:
        return format_override
    if filepath and filepath != "-":
        lower = filepath.lower()
        if lower.endswith('.gz'):
            lower = lower[:-3]
        if lower.endswith(('.fasta', '.fa', '.fna')):
            return 'fasta'
        elif lower.endswith(('.fastq', '.fq')):
            return 'fastq'
    return 'fastq'  # default


def stream_inputs(args, with_ends: bool = False) -> Iterator[Tuple[str, bool, bool, bool]]:
    """
    Uniform generator yielding (sequence, is_paired, is_read1, is_read2).
    Dispatches to appropriate strategy based on input configuration.

    Input modes:
    - 0 files: read from stdin (single or interleaved)
    - 1 file: read from file (single or interleaved, or ZNA for reencoding)
    - 2 files: paired-end (read1, read2)

    *with_ends* applies to the ZNA re-encode mode only, where it appends each
    record's boundary geometry: ``(seq, is_paired, is_read1, is_read2,
    has_start, has_end)``.  Every other input mode is producing fresh records
    that have no orientation history, so they are unaffected.
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
        yield from _stream_zna_reencode(files[0], with_ends=with_ends)
        return
    
    with ExitStack() as stack:
        # Strategy 1: Two Files = Paired-End
        if len(files) == 2:
            f1 = stack.enter_context(get_input_handle(files[0]))
            f2 = stack.enter_context(get_input_handle(files[1]))
            yield from _stream_paired_files(f1, f2, files[0], files[1], format_override)

        # Strategy 2: Interleaved Mode
        elif args.interleaved:
            src = files[0] if len(files) == 1 else None
            f = stack.enter_context(get_input_handle(src))
            
            format_type = _infer_format(src, format_override)
            
            if format_type == 'fastq':
                yield from _stream_interleaved_fastq(f)
            else:
                yield from _stream_interleaved_fasta(f)

        # Strategy 3: Single-End
        else:
            src = files[0] if len(files) == 1 else None
            f = stack.enter_context(get_input_handle(src))
            yield from _stream_single_end(f, src, format_override)


_COMPLEMENT_TABLE = str.maketrans("ACGTacgt", "TGCAtgca")
_COMPLEMENT_CHAR = {"A": "T", "C": "G", "G": "C", "T": "A",
                    "a": "t", "c": "g", "g": "c", "t": "a"}


def _is_reverse_complement(a: str, b: str) -> bool:
    """True when *b* is the reverse complement of *a*.

    Two mates of a pair whose insert was at or below the read length span the
    identical fragment interval, so they are exact reverse complements — that is
    how a full-fragment pair is recognised without knowing the insert size.

    The end probes reject ~15/16 of ordinary pairs in two character
    comparisons, so the O(n) translate only runs on genuine candidates.
    """
    n = len(a)
    if n == 0 or n != len(b):
        return False
    if _COMPLEMENT_CHAR.get(a[-1]) != b[0] or _COMPLEMENT_CHAR.get(a[0]) != b[-1]:
        return False
    return b == a.translate(_COMPLEMENT_TABLE)[::-1]


def _fragment_units(records) -> Iterator[list]:
    """Group a record stream into fragments: ``[R1, R2]`` or ``[single]``.

    A paired R1 is grouped with the immediately following paired R2; anything
    else stands alone.  Grouping lets record-level policies (the N-drop filter,
    full-overlap detection) act on a whole fragment, so they can never leave a
    lone mate behind — a lone ``IS_PAIRED`` record would otherwise be encoded as
    half a fragment and mis-read downstream as a full molecule.
    """
    pending = None
    for rec in records:
        if pending is not None:
            if rec[1] and rec[3]:          # paired R2 completes the held R1
                yield [pending, rec]
                pending = None
            else:
                yield [pending]
                pending = rec if (rec[1] and rec[2]) else None
                if pending is None:
                    yield [rec]
        elif rec[1] and rec[2]:            # paired R1 — hold for its mate
            pending = rec
        else:
            yield [rec]
    if pending is not None:
        yield [pending]


def _full_fragment_flags(unit: list, treat_unpaired_as_merged: bool) -> list:
    """Which records of *unit* span their entire fragment (both edges real)."""
    if len(unit) == 1:
        return [treat_unpaired_as_merged]
    full = _is_reverse_complement(unit[0][0], unit[1][0])
    return [full, full]


def _make_label_extractor(label_defs: tuple[LabelDef, ...], tag_map: dict):
    """Return ``extract(raw_header) -> labels``, using the C++ path when built."""
    from .dtypes import resolve_missing as _resolve_missing

    num_labels = len(label_defs)
    if _accel_extract is not None:
        tag_specs = []
        for ldef in label_defs:
            tag_bytes = ldef.effective_tag.encode('ascii')
            if ldef.dtype.code == 'A':
                conv_code = _CONV_ORD
            elif ldef.dtype.code in ('f', 'd'):
                conv_code = _CONV_FLOAT
            else:
                conv_code = _CONV_INT
            tag_specs.append((tag_bytes, conv_code))
        missing_values = tuple(_resolve_missing(ld) for ld in label_defs)

        def extract(raw_header):
            return _accel_extract(raw_header, tag_specs, num_labels, missing_values)
    else:
        def extract(raw_header):
            return extract_labels_from_header(
                raw_header, tag_map, num_labels, label_defs=label_defs
            )
    return extract


def _labeled_seqs(f: BinaryIO, extract) -> Iterator[Tuple[str, tuple]]:
    """Yield ``(sequence, labels)`` from a labeled FASTQ stream.

    Deliberately does not parse the read name: only interleaved input needs it,
    and this runs on every read.
    """
    for raw_header, seq in parse_fastq_with_headers(f):
        yield seq, extract(raw_header)


def _labeled_named(f: BinaryIO, extract) -> Iterator[Tuple[bytes, int, Tuple[str, tuple]]]:
    """Yield ``(base_name, suffix_number, (sequence, labels))``.

    For interleaved pairing only, and keyed the same way as
    :func:`parse_fastq_keyed` so both feed one :func:`_pair_interleaved`.
    """
    for raw_header, seq in parse_fastq_with_headers(f):
        base, suffix = _read_key(raw_header)
        yield base, suffix, (seq, extract(raw_header))


def stream_inputs_labeled(
    args, label_defs: tuple[LabelDef, ...], tag_map: dict
) -> Iterator[Tuple[str, bool, bool, bool, tuple]]:
    """Yield ``(seq, is_paired, is_read1, is_read2, labels)`` for labeled encoding.

    Dispatches over the same input modes as :func:`stream_inputs` — two paired
    files, interleaved, or single-end — so that labeled encoding preserves
    pairing.  It previously assumed single-end and flagged every record
    unpaired, which silently discarded the R1/R2 flags of interleaved input.

    Labels are read from SAM tags in the FASTQ header, so FASTA input is not
    supported.  Uses the C++ fast label extractor when available.
    """
    files = args.files if args.files else []
    extract = _make_label_extractor(label_defs, tag_map)

    format_override = None
    if getattr(args, 'fasta', False):
        format_override = 'fasta'
    elif getattr(args, 'fastq', False):
        format_override = 'fastq'

    with ExitStack() as stack:
        # Strategy 1: two files = paired-end
        if len(files) == 2:
            if _infer_format(files[0], format_override) != 'fastq':
                sys.exit("Error: labeled encoding requires FASTQ input "
                         "(labels are parsed from SAM tags in the read header).")
            f1 = stack.enter_context(get_input_handle(files[0]))
            f2 = stack.enter_context(get_input_handle(files[1]))
            for (s1, l1), (s2, l2) in zip(_labeled_seqs(f1, extract),
                                          _labeled_seqs(f2, extract)):
                yield s1, True, True, False, l1
                yield s2, True, False, True, l2
            return

        src = files[0] if len(files) == 1 else None
        if _infer_format(src, format_override) != 'fastq':
            sys.exit("Error: labeled encoding requires FASTQ input "
                     "(labels are parsed from SAM tags in the read header).")
        f = stack.enter_context(get_input_handle(src))

        # Strategy 2: interleaved — share the pairing rule with stream_inputs
        if getattr(args, 'interleaved', False):
            for (seq, labels), is_paired, is_read1, is_read2 in _pair_interleaved(
                    _labeled_named(f, extract)):
                yield seq, is_paired, is_read1, is_read2, labels
            return

        # Strategy 3: single-end
        for seq, labels in _labeled_seqs(f, extract):
            yield seq, False, False, False, labels


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

    # A .zna is only ever valid as input on its own, in re-encode mode.  Given two of
    # them the check above never fires, so the binary went to the FASTQ parser, which
    # found no records and wrote a valid 0-record file with a zero exit status --
    # warning only about *format inference*, never about the real problem.  A library
    # could disappear from a corpus with every stage green.
    if not is_reencoding:
        for path in files:
            if is_zna_file(path):
                sys.exit(
                    f"Error: {path} is already a ZNA file. Re-encoding takes exactly "
                    f"one .zna input and no other files; encoding it as FASTQ would "
                    f"write an empty output and report success."
                )

    merge_pairs = getattr(args, 'merge_pairs', False)
    if merge_pairs:
        # Each rejection says WHY, not just that (MERGE_PAIRS_PLAN.md §5).
        if len(files) != 2:
            sys.exit("Error: --merge-pairs is a two-file mode: give R1 and R2 "
                     "as exactly two positional files. (argparse requires the "
                     "two paths to be adjacent: put option flags before or "
                     "after both, not between them.)")
        if "-" in files:
            sys.exit("Error: --merge-pairs cannot read stdin; the merge "
                     "reader takes paths and drives its own decompression.")
        if getattr(args, 'interleaved', False):
            sys.exit("Error: --merge-pairs takes two separate mate files, "
                     "not an interleaved stream.")
        if getattr(args, 'fasta', False):
            sys.exit("Error: --merge-pairs needs FASTQ: the merge decision "
                     "is quality-aware, and FASTA has no qualities.")
        if getattr(args, 'treat_unpaired_as_merged', False):
            sys.exit("Error: --merge-pairs makes --treat-unpaired-as-merged "
                     "meaningless: the merger hands over each record's exact "
                     "geometry, which is what that blanket assertion "
                     "approximates.")
        if args.seq_len_bytes == 1:
            sys.exit("Error: --merge-pairs needs --seq-len-bytes >= 2: a "
                     "merged record can be up to twice a read long, and with "
                     "1-byte lengths (max 255) a 2x150 merge would fail "
                     "mid-stream, leaving a partial file.")

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
    # Default to compression (Zstd level 3) unless explicitly disabled
    should_compress = True 
    if args.compress_flag is not None:
        should_compress = args.compress_flag
    
    comp_method = COMPRESSION_ZSTD if should_compress else COMPRESSION_NONE
    comp_level = args.level

    # Handle strand-specific flags
    strand_specific_flag = args.strand_specific
    strand_normalize_flag = getattr(args, 'strand_normalize', False)
    
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
        strand_normalized = strand_normalize_flag if strand_normalize_flag else input_header.strand_normalized
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
        strand_normalized = strand_normalize_flag

    # An already-normalized input carries its orientation in the per-record
    # IS_RC flags. Copy it verbatim: orientation is not idempotent, and
    # re-deriving it here would reverse-complement a second time, silently
    # un-normalizing the file while the header still claimed otherwise.
    preserve_normalization = bool(
        is_reencoding and input_header is not None and input_header.strand_normalized
    )
    if preserve_normalization and strand_specific != input_header.strand_specific:
        sys.exit(
            "Error: cannot change strand-specificity when re-encoding an "
            "already strand-normalized file. Its orientation was applied at "
            "encode time and cannot be re-derived. Decode with --restore-strand "
            "and re-encode from the original reads instead."
        )
    # Build label definitions (if any)
    label_defs: tuple[LabelDef, ...] = ()
    tag_map: dict = {}
    label_defs_path = getattr(args, 'label_defs', None)
    label_specs = getattr(args, 'label', None) or []
    label_descs = getattr(args, 'label_desc', None) or []
    if label_defs_path:
        yaml_defs = load_label_file(label_defs_path)
        if label_specs or label_descs:
            label_defs = merge_label_defs(yaml_defs, label_specs, label_descs)
        else:
            label_defs = yaml_defs
        tag_map = build_tag_extractor(label_defs)
    elif label_specs:
        label_defs = build_label_defs(label_specs, label_descs)
        tag_map = build_tag_extractor(label_defs)

    # Labels are not carried through the ZNA -> ZNA re-encode path. Refuse
    # loudly: a labeled input widens the record tuple and dies mid-stream, and
    # --label routes the input through the FASTQ label parser, which finds no
    # records in a binary ZNA and writes an empty file with a zero exit status.
    if is_reencoding:
        if input_header is not None and input_header.labels:
            sys.exit(
                "Error: re-encoding a labeled ZNA file is not supported "
                f"({files[0]} defines {len(input_header.labels)} label(s): "
                f"{', '.join(ld.name for ld in input_header.labels)}). "
                "Its labels would be dropped."
            )
        if label_defs:
            sys.exit(
                "Error: --label/--label-defs cannot be applied when re-encoding "
                "a ZNA file. Labels are read from FASTQ headers at encode time; "
                f"{files[0]} has none to read."
            )

    # Header Setup
    header = ZnaHeader(
        read_group=read_group,
        description=description,
        seq_len_bytes=seq_len_bytes,
        strand_specific=strand_specific,
        read1_antisense=read1_antisense,
        read2_antisense=read2_antisense,
        strand_normalized=strand_normalized,
        compression_method=comp_method,
        compression_level=comp_level,
        labels=label_defs,
    )

    # Warn on a strand-normalization config that would normalize nothing:
    # strand-specific + normalize but neither read1 nor read2 flagged antisense
    # means do_rc_r1, do_rc_r2 (and the single/read1 rule) are all off, so no
    # read is ever reverse-complemented. Usually a misconfigured protocol.
    quiet = hasattr(args, 'quiet') and args.quiet
    if (header.strand_normalized and header.strand_specific
            and not header.read1_antisense and not header.read2_antisense
            and not quiet):
        print(
            "[Warning] --strand-normalize with --strand-specific but neither read1 "
            "nor read2 is antisense: no reads will be reverse-complemented. "
            "Did you mean to keep the default --read1-antisense (dUTP), or to drop "
            "--strand-specific for unstranded random normalization?",
            file=sys.stderr,
        )

    if not is_stdout:
        c_str = f"ZSTD (L{comp_level})" if should_compress else "None"
        if not is_reencoding:
            print(f"[ZNA] Encoding to {args.output} [{c_str}]...", file=sys.stderr)
        else:
            print(f"[ZNA] Reencoding to {args.output} [{c_str}]...", file=sys.stderr)

    # Encode Loop
    count = 0
    # Use ExitStack to safely close output file (or leave stdout open)
    npolicy = getattr(args, 'npolicy', None)
    # trim3 reshapes records; it is not an encoding policy, so the codec never sees it.
    # Keeping the two apart is what stopped `drop` reaching the codec and being silently
    # read as "substitute A" -- see docs/NPOLICY_PLAN.md 8.1.
    trim3 = (npolicy == 'trim3')
    codec_npolicy = 'random' if npolicy == 'random' else ''
    # What the N policy did, for the summary below. `zna merge` has always reported
    # this and `zna encode` did not, so the same policy on the same library was
    # accountable on one side of the seam and silent on the other.
    #
    # trim3's count is free -- `_trim3` is already finding the first no-call. `random`
    # substitutes inside the codec, which returns no count, so encode counts the
    # no-calls itself rather than widening the codec ABI across two backends for a
    # statistic. That costs one extra C-level scan per record, under `random` only.
    npolicy_bases = 0
    npolicy_records = 0
    npolicy_total_bases = 0
    block_size = parse_block_size(args.block_size)
    if preserve_normalization and not quiet:
        print(
            "[ZNA] Input is strand-normalized; copying orientation verbatim.",
            file=sys.stderr,
        )

    # ZNA -> ZNA re-encode copies records in order, so the input's shuffled
    # attestation survives verbatim into the output's provenance prologue.
    # Every other input produces records in source order: not shuffled.  (The
    # post-encode --shuffle pass rewrites the file through shuffle_zna, which
    # stamps its own True.)
    shuffled_in = False
    merged_in_process_in = False
    if is_reencoding:
        with open(files[0], "rb") as _f:
            _prov = ZnaReader(_f).provenance
        shuffled_in = bool(_prov is not None and _prov.shuffled)
        # A copy preserves record content, so the merged-in-process
        # attestation survives re-encode exactly as `shuffled` does.
        merged_in_process_in = bool(_prov is not None and _prov.merged_in_process)

    with ExitStack() as stack:
        f_out = stack.enter_context(get_output_handle(args.output))
        writer = stack.enter_context(ZnaWriter(
            f_out, header, block_size=block_size, npolicy=codec_npolicy,
            preserve_normalization=preserve_normalization,
            rng_seed=getattr(args, 'seed', 0) or 0,
            shuffled=shuffled_in,
            merged_in_process=merge_pairs or merged_in_process_in,
        ))

        def _trim3(seq):
            """Cut a record at its first no-call, keeping ``[0, first N)``.

            3' only. Base 0 is a true fragment boundary -- a read starts at a fragment
            end and runs inward -- so cutting from the far end never disturbs it, and the
            record stays honestly anchored however short it gets. There is deliberately
            no minimum length here: a core read-processing tool reports what it did and
            leaves the length decision to the consumer.

            Returns ``(seq, bases_removed)``; the count is what the summary reports, and
            it is "bases removed", not "no-calls found" -- the same quantity `zna merge`
            reports, so the two halves of the pipeline are directly comparable.
            """
            i = seq.find('N')
            if i < 0:
                i = seq.find('n')
            return (seq, 0) if i < 0 else (seq[:i], len(seq) - i)

        def _count_subs(seq):
            """No-calls the codec is about to substitute under ``--npolicy random``."""
            return seq.count('N') + seq.count('n')
        treat_unpaired_as_merged = getattr(args, 'treat_unpaired_as_merged', False)

        merge_acc = None
        if merge_pairs:
            # Dispatched BEFORE label_defs, deliberately: stream_inputs_labeled
            # also knows the two-file mode and would silently emit an UNMERGED
            # stream, writing an unmerged corpus at exit 0
            # (MERGE_PAIRS_PLAN.md §0.3 -- the audit's third blocker).
            from .merge.encode_stream import stream_merge_pairs
            from .merge.cli import _new_acc
            from .merge.params import MergeParams
            merge_acc = _new_acc()
            merge_params = MergeParams(
                t_merge=getattr(args, 't_merge', 28.0),
                t_trim=getattr(args, 't_trim', 8.0),
                min_read_length=getattr(args, 'min_read_length', 40),
                npolicy=npolicy or "trim3",
                rng_seed=getattr(args, 'seed', 0) or 42,
            )
            stream = stream_merge_pairs(
                files, merge_params,
                check_sync=not getattr(args, 'no_sync_check', False),
                merge_backend=getattr(args, 'merge_backend', 'auto'),
                label_defs=label_defs, tag_map=tag_map,
                allow_empty=getattr(args, 'allow_empty', False),
                acc=merge_acc, quiet=quiet)
            carries_ends = True
        elif label_defs:
            stream = stream_inputs_labeled(args, label_defs, tag_map)
            carries_ends = False
        else:
            stream = stream_inputs(args)
            carries_ends = False

        # Single-end input can never yield a paired record: only two-file,
        # interleaved, and ZNA re-encode modes set IS_PAIRED.  Fragment grouping
        # would then wrap every record in its own one-element list — a generator
        # resume and a list allocation per read — to decide something already
        # known, and every full-fragment decision would return the same constant.
        single_end_input = not (len(files) == 2
                                or getattr(args, 'interleaved', False)
                                or is_reencoding)
        show_progress = not is_stdout and not quiet

        if is_reencoding:
            # ZNA -> ZNA is a COPY, so carry each record's flag byte verbatim rather
            # than any view of it.  It used to go through records(with_ends=True) and
            # rebuild the flags from (has_start, has_end), which loses information in
            # both directions and did so silently:
            #
            #   * a strand-normalized source lost IS_RC on every full-fragment record,
            #     because (True, True) cannot say whether the record was flipped -- so
            #     --restore-strand then returned those records in the wrong orientation;
            #   * a NON-normalized source took the 4-tuple path instead and re-derived
            #     IS_FULL_FRAGMENT from scratch, i.e. dropped it altogether.
            #
            # No fragment grouping and no N-drop here, and neither is an oversight: ZNA
            # stores two bits per base, so a decoded record cannot contain an N for the
            # filter to act on.
            with open(files[0], "rb") as f_re:
                src = ZnaReader(f_re)
                if preserve_normalization:
                    write_copy = writer.write_copy
                    for rec in src.copy_records():
                        write_copy(rec)
                        count += 1
                        if show_progress and count % 1_000_000 == 0:
                            print(f"      Processed {count//1_000_000}M records...",
                                  end='\r', file=sys.stderr)
                else:
                    # The source is not normalized, so the writer owns IS_RC -- it may
                    # be about to apply normalization for the first time.  Hand over
                    # only the bit it cannot derive.
                    write_record = writer.write_record
                    fields = FLAG_FIELDS
                    for seq, flags, _labels in src.copy_records():
                        is_paired, is_read1, is_read2 = fields[flags]
                        write_record(seq, is_paired, is_read1, is_read2,
                                     is_full_fragment=bool(flags & 16))
                        count += 1
                        if show_progress and count % 1_000_000 == 0:
                            print(f"      Processed {count//1_000_000}M records...",
                                  end='\r', file=sys.stderr)
        elif single_end_input:
            is_full = treat_unpaired_as_merged
            write_record = writer.write_record
            for rec in stream:
                seq = rec[0]
                # Membership on the sequence itself: the old test uppercased a
                # copy of every read (150 bytes each) purely to look for one
                # character.
                if trim3:
                    seq, cut = _trim3(seq)
                    if cut:
                        npolicy_bases += cut
                        npolicy_records += 1
                elif codec_npolicy:
                    subs = _count_subs(seq)
                    if subs:
                        npolicy_bases += subs
                        npolicy_records += 1
                npolicy_total_bases += len(seq)
                if label_defs:
                    write_record(seq, rec[1], rec[2], rec[3],
                                 labels=rec[4], is_full_fragment=is_full)
                else:
                    write_record(seq, rec[1], rec[2], rec[3],
                                 is_full_fragment=is_full)
                count += 1
                if show_progress and count % 1_000_000 == 0:
                    print(f"      Processed {count//1_000_000}M records...",
                          end='\r', file=sys.stderr)
        elif carries_ends:
            # The stream hands over each record's geometry: 6-tuples
            # (seq, is_paired, is_read1, is_read2, has_start, has_end), labels
            # seventh when declared.  No N-policy work here -- under
            # --merge-pairs the policy already ran inside the kernel, so this
            # path sees no no-calls at all -- and no _full_fragment_flags:
            # re-deriving what was handed over is the failure mode this seam
            # exists to delete.  IS_RC is never supplied; orientation belongs
            # to the writer's strand-normalization settings, which compose
            # with this stream unchanged (plan §0.2).
            #
            # _fragment_units still groups mates -- the kernel emits MATE1
            # immediately followed by MATE2, and grouping is what the writer's
            # fragment contract expects to see arrive together (plan §6.8).
            for unit in _fragment_units(stream):
                npolicy_total_bases += sum(len(rec[0]) for rec in unit)
                for rec in unit:
                    writer.write_record(
                        rec[0], rec[1], rec[2], rec[3],
                        labels=rec[6] if label_defs else None,
                        is_full_fragment=rec[4] and rec[5])
                count += len(unit)
                if (count % 1_000_000 < len(unit) and count >= 1_000_000
                        and show_progress):
                    print(f"      Processed {count//1_000_000}M records...",
                          end='\r', file=sys.stderr)
        else:
            for unit in _fragment_units(stream):
                if trim3:
                    # Per record, not per fragment: trimming never orphans a mate, so
                    # there is nothing here for the old pair-atomic drop to protect.
                    cut_unit = []
                    for rec in unit:
                        seq, cut = _trim3(rec[0])
                        if cut:
                            npolicy_bases += cut
                            npolicy_records += 1
                        cut_unit.append((seq,) + tuple(rec[1:]))
                    unit = cut_unit
                elif codec_npolicy:
                    for rec in unit:
                        subs = _count_subs(rec[0])
                        if subs:
                            npolicy_bases += subs
                            npolicy_records += 1
                npolicy_total_bases += sum(len(rec[0]) for rec in unit)
                full = _full_fragment_flags(unit, treat_unpaired_as_merged)
                if label_defs:
                    for rec, is_full in zip(unit, full):
                        writer.write_record(rec[0], rec[1], rec[2], rec[3],
                                            labels=rec[4], is_full_fragment=is_full)
                else:
                    for rec, is_full in zip(unit, full):
                        writer.write_record(rec[0], rec[1], rec[2], rec[3],
                                            is_full_fragment=is_full)
                count += len(unit)
                if (count % 1_000_000 < len(unit) and count >= 1_000_000
                        and show_progress):
                    print(f"      Processed {count//1_000_000}M records...",
                          end='\r', file=sys.stderr)

    if merge_pairs and merge_acc is not None:
        # The records the merger emitted ARE the records in the file -- there
        # is no dropper between the kernel and the writer.  Not an input
        # error: if these ever disagree the adapter or the write loop lost
        # records, which must be a crash, not a statistic
        # (MERGE_PAIRS_PLAN.md §10.2; subsumed at read time by
        # `zna inspect --verify` check 4).
        emitted = merge_acc[0][4]
        if emitted != count:
            raise AssertionError(
                f"merge emitted {emitted} records but {count} were written")
        if not quiet and not is_stdout:
            c = merge_acc[0]
            print(f"[ZNA] merge: {c[0]} pairs -> {c[1]} merged, "
                  f"{c[2]} trimmed, {c[3]} kept; {c[5]} dropped short",
                  file=sys.stderr)
        merge_json_path = getattr(args, 'merge_json', None)
        if merge_json_path:
            import json as _json
            from .merge.cli import _assemble_stats
            from .merge.backend import get_merge_backend_name
            stats = _assemble_stats(merge_acc, merge_params)
            # _assemble_stats reads the process-global backend selection, which
            # --merge-backend deliberately never mutates; report the kernel
            # that actually ran.
            stats["backend"] = get_merge_backend_name(
                None if getattr(args, 'merge_backend', 'auto') == 'auto'
                else args.merge_backend)
            # Same block `zna merge --json` writes; "tool" stays "zna-merge"
            # deliberately -- hulkrna's gather consumes these blocks and may
            # key on it, and khorana no longer reads merge.json at all
            # (TRAILER_PLAN §8).
            with open(merge_json_path, "w") as jf:
                _json.dump(stats, jf, indent=2)
                jf.write("\n")

    # ── Optional shuffle pass ─────────────────────────────────────────
    if getattr(args, 'shuffle', False) and not is_stdout:
        if not quiet:
            print(f"\n[ZNA] Shuffling ...", file=sys.stderr)
        # Shuffle in-place via a temp file
        tmp_fd, tmp_shuffle = tempfile.mkstemp(
            suffix=".zna", dir=os.path.dirname(args.output) or "."
        )
        os.close(tmp_fd)
        try:
            shuffle_zna(
                args.output,
                tmp_shuffle,
                seed=getattr(args, 'seed', 42),
                # --shuffle-buffer-size was parsed and then ignored here, so the
                # documented flag did nothing and the budget was always 1 GiB.
                # Its default still parses to 1 GiB, so this changes no existing
                # invocation's behaviour -- it only makes the flag work.
                buffer_bytes=parse_block_size(
                    getattr(args, 'shuffle_buffer_size', None) or "1G"
                ),
                block_size=block_size,
                quiet=True,
            )
            os.replace(tmp_shuffle, args.output)
        except Exception:
            if os.path.exists(tmp_shuffle):
                os.unlink(tmp_shuffle)
            raise

    duration = time.time() - start_time
    if not is_stdout and not quiet:
        print(f"\n[ZNA] Done. Wrote {count} records in {duration:.2f}s.", file=sys.stderr)
        # Always say what the N policy did, exactly as `zna merge` does. The failure
        # this guards against is silent: one dark cycle can make a policy consume most
        # of a library while the run still ends with "Done".
        #
        # Skipped when re-encoding a .zna, and that is not an oversight -- ZNA stores
        # two bits per base, so a decoded record cannot contain a no-call and no policy
        # ran. Reporting "removed 0 bases" there would imply one had been applied.
        if merge_pairs and merge_acc is not None:
            # Under --merge-pairs the policy ran inside the kernel, after
            # overlap rescue; encode's own counters are legitimately zero and
            # reporting them would claim the policy did nothing.  Report the
            # kernel's numbers -- the same quantities `zna merge` prints.
            npolicy_bases = merge_acc[0][13]
            npolicy_records = 0        # the kernel counts bases, not records
        if npolicy and not is_reencoding and merge_pairs and merge_acc is not None:
            verb = "removed" if npolicy == "trim3" else "substituted"
            rescued = merge_acc[0][14]
            print(f"[ZNA] no-calls: --npolicy {npolicy} {verb} {npolicy_bases} "
                  f"base{'' if npolicy_bases == 1 else 's'} in the merge kernel; "
                  f"{rescued} rescued from the mate at no cost.",
                  file=sys.stderr)
        elif npolicy and not is_reencoding:
            verb = "removed" if npolicy == "trim3" else "substituted"
            print(f"[ZNA] no-calls: --npolicy {npolicy} {verb} {npolicy_bases} "
                  f"base{'' if npolicy_bases == 1 else 's'} "
                  f"in {npolicy_records} record{'' if npolicy_records == 1 else 's'}.",
                  file=sys.stderr)
            if npolicy_total_bases and npolicy_bases / npolicy_total_bases > 0.01:
                pct = 100.0 * npolicy_bases / npolicy_total_bases
                print(f"[ZNA] WARNING: --npolicy {npolicy} affected {pct:.1f}% of "
                      f"encoded bases. That is high enough to be a run problem (a dark "
                      f"cycle, a failed tile) rather than ordinary no-calls -- check "
                      f"the input before using this library.", file=sys.stderr)


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

        # Decoding a normalized file without --restore-strand emits the STORED
        # (normalized) frame.  Re-encoding that output with --strand-normalize
        # applies orientation a second time and desynchronizes IS_RC from the
        # bases — silently, because the two files' headers are identical.
        if (reader.header.strand_normalized
                and not getattr(args, 'restore_strand', False)
                and not args.quiet):
            print(
                "[Warning] This file is strand-normalized and --restore-strand was "
                "not given, so the output is in NORMALIZED orientation. Do not "
                "re-encode it with --strand-normalize: orientation is not idempotent, "
                "and applying it twice un-normalizes the data while the header still "
                "claims otherwise. Use --restore-strand to recover the original read "
                "orientation.",
                file=sys.stderr,
            )

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
            show_labels = getattr(args, 'labels', False)
            has_labels = len(reader.header.labels) > 0

            if has_labels and show_labels:
                for rec in reader.records(restore_strand=restore_strand):
                    seq, is_paired, is_r1, is_r2, labels = rec
                    if is_r1 or (not is_paired):
                        counter += 1
                    suffix = ""
                    if is_r1: suffix = "/1"
                    elif is_r2: suffix = "/2"

                    # Build SAM-style tag string
                    tag_parts = []
                    for ldef, val in zip(reader.header.labels, labels):
                        if ldef.dtype.code in ('f', 'd'):
                            tag_parts.append(f"{ldef.name}:f:{val}")
                        elif ldef.dtype.code == 'A':
                            tag_parts.append(f"{ldef.name}:A:{chr(val)}")
                        else:
                            tag_parts.append(f"{ldef.name}:i:{val}")
                    tag_str = "\t".join(tag_parts)
                    # One f-string, not a header f-string interpolated into a
                    # second one: this runs on every record.
                    record = f">{rg}:{counter}{suffix}\t{tag_str}\n{seq}\n"

                    if mode == "split" and is_r2:
                        outs[1].write(record)
                    else:
                        outs[0].write(record)
            else:
                for rec in reader.records(restore_strand=restore_strand):
                    if has_labels:
                        seq, is_paired, is_r1, is_r2, _labels = rec
                    else:
                        seq, is_paired, is_r1, is_r2 = rec
                    if is_r1 or (not is_paired):
                        counter += 1
                    suffix = ""
                    if is_r1: suffix = "/1"
                    elif is_r2: suffix = "/2"
                    record = f">{rg}:{counter}{suffix}\n{seq}\n"

                    if mode == "split" and is_r2:
                        outs[1].write(record)
                    else:
                        outs[0].write(record)

    except BrokenPipeError:
        sys.stderr.close()
    finally:
        if input_file and input_file != "-":
            f_in_handle.close()


# --- COMMAND: INSPECT ---

def _verify_zna(fh, reader) -> dict:
    """Full certification of one open ZNA file (``zna inspect --verify``).

    Decodes every block and checks, in order, stopping at the first failure:

      1. footer + trailer + prologue present and parse (the reader raises on a
         corrupt trailer; ABSENCE is the two-cause failure: pre-0.5 writer, or
         truncation);
      2. the stored block index matches a walk of the actual headers, and the
         sentinel sits where the footer says;
      3. every block decompresses (zstd content checksums verify implicitly)
         and its columns split cleanly -- including the lengths-vs-sequence
         cross-check ``sum(ceil(len/4)) == len(seq_stream)``, which validates
         the lengths column against the sequence bytes without materializing
         a single string;
      4. recounted stats equal the trailer: flag counts, both histograms,
         n_records, n_bases, n_pairs, n_unpaired, and the prologue's crc;
      5. structural invariants: every paired R1 immediately followed by its
         R2, no block ending mid-fragment (the 0.4.1 guarantee).

    "Stranded implies exactly one antisense mate" is reported as a WARNING,
    not a failure: the fuzz matrix deliberately writes both-antisense configs
    and they are legal files (TRAILER_PLAN §14 A4).

    Returns ``{"passed", "checks", "warnings", "failure"}``; certifies
    integrity, structure, and stats fidelity.  It cannot prove ``shuffled``
    (a permutation attests itself no better than a flag) nor that strand
    normalization was applied correctly; those stay stamped facts.
    """
    import zstandard
    from collections import Counter
    from zna.core import (CLOSES_FRAGMENT, OPENS_FRAGMENT,
                          _FOOTER_SIZE as FOOTER_SIZE)

    h = reader.header
    checks: list = []
    warnings_: list = []
    result = {"passed": False, "checks": checks, "warnings": warnings_,
              "failure": None}

    def fail(check: str, detail: str) -> dict:
        result["failure"] = {"check": check, "detail": detail}
        return result

    # -- 1. trailer + prologue presence -----------------------------------
    try:
        trailer = reader.trailer
        prov = reader.provenance
    except ValueError as e:
        return fail("trailer", str(e))
    if trailer is None or prov is None:
        return fail(
            "trailer",
            "no trailer: written by zna < 0.5, or truncated — re-encode "
            "from source either way.",
        )
    checks.append("trailer and prologue present, CRC valid, payloads parse")

    # -- 2. stored index vs the walked file --------------------------------
    try:
        scanned = reader.scan_block_index()
    except (ValueError, EOFError) as e:
        return fail("index", f"walking the block chain failed: {e}")
    stored = [(b.offset, b.comp_size, b.uncomp_size, b.n_records)
              for b in reader.block_index()]
    walked = [(b.offset, b.comp_size, b.uncomp_size, b.n_records)
              for b in scanned]
    if stored != walked:
        return fail(
            "index",
            f"stored block index disagrees with the file: "
            f"{len(stored)} stored vs {len(walked)} walked block(s), first "
            f"difference at block "
            f"{next(i for i, (a, b) in enumerate(zip(stored + [None], walked + [None])) if a != b)}.",
        )
    fh.seek(0, 2)
    file_size = fh.tell()
    sentinel_offset = (walked[-1][0] + _BLOCK_HEADER_SIZE + walked[-1][1]
                       if walked else trailer.data_start)
    fh.seek(sentinel_offset)
    sent = fh.read(_BLOCK_HEADER_SIZE)
    sc, su, sn, _sf, _sl = struct.unpack(_BLOCK_HEADER_FMT, sent)
    if sn != 0:
        return fail("index", f"no sentinel at offset {sentinel_offset} after "
                             f"the last data block.")
    if sentinel_offset + _BLOCK_HEADER_SIZE + sc != file_size:
        return fail(
            "index",
            f"the sentinel's payload does not reach EOF: sentinel at "
            f"{sentinel_offset}, comp_size {sc}, file size {file_size} — "
            f"trailing bytes after the footer, or a mis-sized trailer.",
        )
    checks.append(f"stored index matches the walk ({len(walked)} blocks); "
                  f"sentinel and footer agree")

    # -- 3 + 4 + 5. decode every block, recount, check structure -----------
    dctx = (zstandard.ZstdDecompressor()
            if h.compression_method == COMPRESSION_ZSTD else None)
    slb = h.seq_len_bytes
    len_ch = {1: "B", 2: "H", 4: "I"}[slb]
    label_bytes_per_rec = sum(ld.dtype.size for ld in h.labels)
    flag_counts: Counter = Counter()
    len_hist: Counter = Counter()
    len_hist_unpaired: Counter = Counter()
    n_records = n_bases = 0
    for b in scanned:
        fh.seek(b.offset)
        _c, u_size, _n, flags_size, lengths_size = struct.unpack(
            _BLOCK_HEADER_FMT, fh.read(_BLOCK_HEADER_SIZE))
        payload = fh.read(b.comp_size)
        if len(payload) != b.comp_size:
            return fail("blocks", f"block {b.index}: payload truncated.")
        if dctx is not None:
            try:
                data = dctx.decompress(payload, max_output_size=b.uncomp_size)
            except zstandard.ZstdError as e:
                return fail("blocks", f"block {b.index}: {e}")
        else:
            data = payload
        count = b.n_records
        if len(data) != u_size:
            return fail("blocks", f"block {b.index}: decompressed to "
                                  f"{len(data)} bytes, header says {u_size}.")
        if flags_size != count or lengths_size != count * slb:
            return fail("blocks", f"block {b.index}: column sizes do not "
                                  f"match its record count.")
        flags = data[:flags_size]
        lengths_off = flags_size + count * label_bytes_per_rec
        lengths = struct.unpack(
            f"<{count}{len_ch}", data[lengths_off:lengths_off + lengths_size])
        seq_stream_len = len(data) - lengths_off - lengths_size
        if sum((L + 3) >> 2 for L in lengths) != seq_stream_len:
            return fail("blocks", f"block {b.index}: lengths column does not "
                                  f"tile its sequence bytes.")
        expect_mate = False
        for fl in flags:
            if CLOSES_FRAGMENT[fl] != expect_mate:
                return fail("structure", f"block {b.index}: paired R1/R2 "
                                         f"adjacency violated.")
            expect_mate = OPENS_FRAGMENT[fl]
        if expect_mate:
            return fail("structure", f"block {b.index} ends on an unmatched "
                                     f"paired R1 (fragment split across blocks).")
        n_records += count
        for fl in set(flags):
            flag_counts[fl] += flags.count(fl)
        n_bases += sum(lengths)
        for L, fl in zip(lengths, flags):
            len_hist[L] += 1
            if not (fl & 0x04):
                len_hist_unpaired[L] += 1
    checks.append(f"{len(walked)} block(s) decode cleanly; columns tile; "
                  f"fragments whole")

    n_paired = sum(c for f, c in flag_counts.items() if f & 0x04)
    recount = {
        "n_records": n_records, "n_bases": n_bases,
        "n_pairs": n_paired // 2, "n_unpaired": n_records - n_paired,
        "flag_counts": dict(flag_counts),
        "length_histogram": dict(len_hist),
        "length_histogram_unpaired": dict(len_hist_unpaired),
    }
    stated = {
        "n_records": trailer.n_records, "n_bases": trailer.n_bases,
        "n_pairs": trailer.n_pairs, "n_unpaired": trailer.n_unpaired,
        "flag_counts": trailer.flag_counts,
        "length_histogram": trailer.length_histogram,
        "length_histogram_unpaired": trailer.length_histogram_unpaired,
    }
    for key in stated:
        if stated[key] != recount[key]:
            return fail("stats", f"trailer's {key} does not match a recount "
                                 f"of the file.")
    if trailer.prologue_crc32 != prov.crc32:
        return fail("stats", "trailer's prologue_crc32 does not match the "
                             "prologue as stored.")
    checks.append("recounted stats equal the trailer")

    if (h.strand_specific and h.strand_normalized
            and h.read1_antisense == h.read2_antisense):
        warnings_.append(
            "stranded and normalized, but both mates share one antisense "
            "setting — legal, but unusual for a real protocol."
        )

    result["passed"] = True
    return result


def _inspect_json(args, fh, reader, h, file_size) -> None:
    """Emit one machine-readable object describing the file.

    Built for cataloguing a corpus: a manifest that records ``n_records`` per
    file can weight a balanced sample without opening any of them at read time.
    The block walk this uses decompresses nothing — it reads each 20-byte block
    header and seeks over the payload — so it is fast enough to run across
    thousands of files.

    ``--counts`` additionally tallies the per-record flags, which does have to
    decompress each block's flags column.
    """
    import json

    index = reader.block_index()
    out = {
        "path": str(args.input),
        "file_size": file_size,
        "format_version": _VERSION,
        "header": {
            "read_group": h.read_group,
            "description": h.description,
            "seq_len_bytes": h.seq_len_bytes,
            "max_seq_len": (1 << (8 * h.seq_len_bytes)) - 1,
            "strand_specific": h.strand_specific,
            "read1_antisense": h.read1_antisense,
            "read2_antisense": h.read2_antisense,
            "strand_normalized": h.strand_normalized,
            "compression": ("zstd" if h.compression_method == COMPRESSION_ZSTD
                            else "none"),
            "compression_level": h.compression_level,
            "labels": [
                {"id": ld.label_id, "name": ld.name, "dtype": ld.dtype.code,
                 "dtype_name": ld.dtype.name, "description": ld.description,
                 "missing": ld.missing}
                for ld in h.labels
            ],
        },
        "n_blocks": len(index),
        "n_records": sum(b.n_records for b in index),
        "compressed_payload": sum(b.comp_size for b in index),
        "uncompressed_payload": sum(b.uncomp_size for b in index),
    }
    comp = out["compressed_payload"]
    out["compression_ratio"] = (
        out["uncompressed_payload"] / comp if comp else 1.0
    )

    # The two metadata payloads, verbatim -- khorana's indexer reads these.
    prov = reader.provenance
    trailer = reader.trailer          # raises on a corrupt trailer, on purpose
    out["provenance"] = prov.raw if prov is not None else None
    out["trailer"] = trailer.raw if trailer is not None else None

    if getattr(args, 'verify', False):
        out["verify"] = _verify_zna(fh, reader)

    if getattr(args, 'blocks', False):
        out["blocks"] = [
            {"index": b.index, "offset": b.offset, "n_records": b.n_records,
             "comp_size": b.comp_size, "uncomp_size": b.uncomp_size}
            for b in index
        ]

    if getattr(args, 'counts', False):
        out["counts"] = _scan_flag_counts(fh, reader, h)

    json.dump(out, sys.stdout, indent=2)
    sys.stdout.write("\n")
    if getattr(args, 'verify', False) and not out["verify"]["passed"]:
        sys.exit(1)


def _scan_flag_counts(fh, reader, h) -> dict:
    """Per-flag record tallies, decompressing only each block's flags column."""
    import zstandard

    dctx = (zstandard.ZstdDecompressor()
            if h.compression_method == COMPRESSION_ZSTD else None)
    tally = {"paired_r1": 0, "paired_r2": 0, "single_or_merged": 0,
             "reverse_complemented": 0, "full_fragment": 0,
             "rc_by_mate": {"R1": 0, "R2": 0, "single_or_merged": 0}}

    fh.seek(reader._data_start)
    while True:
        b_header = fh.read(_BLOCK_HEADER_SIZE)
        if not b_header or len(b_header) < _BLOCK_HEADER_SIZE:
            break
        c_size, _u, _n, flags_size, _l = struct.unpack(_BLOCK_HEADER_FMT, b_header)
        if _n == 0:
            break  # the trailer sentinel: end of data
        if dctx is not None:
            flags_bytes = dctx.stream_reader(fh.read(c_size)).read(flags_size)
        else:
            flags_bytes = fh.read(flags_size)
            fh.seek(c_size - flags_size, 1)
        for fl in flags_bytes:
            rc = fl & 8
            if rc:
                tally["reverse_complemented"] += 1
            if fl & 16:
                tally["full_fragment"] += 1
            if fl & 1:
                tally["paired_r1"] += 1
                if rc:
                    tally["rc_by_mate"]["R1"] += 1
            elif fl & 2:
                tally["paired_r2"] += 1
                if rc:
                    tally["rc_by_mate"]["R2"] += 1
            else:
                tally["single_or_merged"] += 1
                if rc:
                    tally["rc_by_mate"]["single_or_merged"] += 1
    return tally


def inspect_command(args):
    f_path = Path(args.input)
    if not f_path.exists():
        sys.exit(f"Error: File {args.input} not found.")

    as_json = getattr(args, 'json', False)
    file_size = f_path.stat().st_size
    if not as_json:
        print(f"File: {args.input}")
        print(f"Total Size: {file_size / (1024*1024):.2f} MB")

    with open(args.input, "rb") as f:
        try:
            reader = ZnaReader(f)
            h = reader.header
        except Exception as e:
            sys.exit(f"Error reading header: {e}")

        if as_json:
            _inspect_json(args, f, reader, h, file_size)
            return

        print("\n--- Header Metadata ---")
        print(f"Read Group:       {h.read_group}")
        print(f"Description:      {h.description}")
        print(f"Seq Length:       {h.seq_len_bytes} bytes (Max: {(1<<(8*h.seq_len_bytes))-1} bp)")
        print(f"Strand Specific:  {h.strand_specific}")
        if h.strand_specific:
            print(f"R1 Antisense:     {h.read1_antisense}")
            print(f"R2 Antisense:     {h.read2_antisense}")
        print(f"Strand Normalized:{h.strand_normalized}")
        
        method = "None"
        if h.compression_method == COMPRESSION_ZSTD:
            method = f"ZSTD (Level {h.compression_level})"
        print(f"Compression:      {method}")

        # Label schema
        if h.labels:
            print(f"\nLabels: {len(h.labels)}")
            for ldef in h.labels:
                type_info = f"{ldef.dtype.name:<8} ({ldef.dtype.code})"
                desc_str = f'  "{ldef.description}"' if ldef.description else ""
                miss_str = f"  [missing={ldef.missing}]" if ldef.missing is not None else ""
                print(f"  [{ldef.label_id}] {ldef.name:<4}{type_info}{desc_str}{miss_str}")

        # Corrupt metadata is an error, not a shrug: a trailer that is
        # PRESENT but damaged must never read as "absent".
        try:
            prov = reader.provenance
            trailer = reader.trailer
        except ValueError as e:
            sys.exit(f"Error: {e}")

        print("\n--- Provenance ---")
        if prov is not None:
            print(f"Writer:           zna {prov.writer_version}")
            print(f"Shuffled:         {prov.shuffled}")
            if prov.merged_in_process:
                print(f"Merged in-process:True")
        else:
            print("(none: written by zna < 0.5, or truncated)")

        # Block stats come from block_index(): the stored trailer index on a
        # 0.5 file, one 20-byte read per block otherwise.  The hand-rolled walk
        # this replaces seeked to a hand-computed header offset that measured
        # the read group in CHARACTERS, not UTF-8 bytes -- wrong on any
        # non-ASCII read group -- and broke silently on a short read, so damage
        # under-reported instead of erroring (TRAILER_PLAN T6).
        try:
            index = reader.block_index()
        except (ValueError, EOFError) as e:
            sys.exit(f"Error walking blocks: {e}")
        block_count = len(index)
        total_records = sum(b.n_records for b in index)
        compressed_payload = sum(b.comp_size for b in index)
        uncompressed_payload = sum(b.uncomp_size for b in index)

        if trailer is not None:
            # The cheap structural pass every inspect gets for free: the
            # stored index against a walk of the actual headers.
            walked = reader.scan_block_index()
            if [tuple(b) for b in index] != [tuple(b) for b in walked]:
                sys.exit(
                    f"Error: stored block index disagrees with the file "
                    f"({len(index)} stored vs {len(walked)} walked) — "
                    f"corrupt or truncated. Run --verify for detail."
                )

        if getattr(args, 'verify', False):
            v = _verify_zna(f, reader)
            print("\n--- Verification ---")
            for c in v["checks"]:
                print(f"  ok  {c}")
            for wmsg in v["warnings"]:
                print(f"  [Warning] {wmsg}")
            if not v["passed"]:
                fl = v["failure"]
                sys.exit(f"  FAILED [{fl['check']}] {fl['detail']}")
            print("  PASSED: integrity, structure, and stats certified")

        count_flags = getattr(args, 'counts', False)
        if count_flags:
            tally = _scan_flag_counts(f, reader, h)
            n_paired_r1 = tally["paired_r1"]
            n_paired_r2 = tally["paired_r2"]
            n_single = tally["single_or_merged"]
            n_rc = tally["reverse_complemented"]
            n_full = tally["full_fragment"]
            rc_by_mate = {"R1": tally["rc_by_mate"]["R1"],
                          "R2": tally["rc_by_mate"]["R2"],
                          "single": tally["rc_by_mate"]["single_or_merged"]}

        print("\n--- Content Statistics ---")
        print(f"Total Blocks:       {block_count}")
        print(f"Total Records:      {total_records}")
        if trailer is not None:
            print(f"Total Bases:        {trailer.n_bases}")
            print(f"Pairs / Unpaired:   {trailer.n_pairs} / {trailer.n_unpaired}")

        if count_flags:
            print(f"  Paired R1:        {n_paired_r1}")
            print(f"  Paired R2:        {n_paired_r2}")
            print(f"  Single/merged:    {n_single}")
            print(f"  Reverse-comp'd:   {n_rc}")
            print(f"  Full-fragment:    {n_full}   (both edges are true fragment boundaries)")
            print("\n--- Fragment Geometry (mate x IS_RC) ---")
            for mate, total in (("R1", n_paired_r1), ("R2", n_paired_r2),
                                ("single/merged", n_single)):
                key = {"R1": "R1", "R2": "R2", "single/merged": "single"}[mate]
                rc_n = rc_by_mate[key]
                pct = (100.0 * rc_n / total) if total else 0.0
                print(f"  {mate:<14s} {total:>10d} records, {rc_n:>10d} RC'd ({pct:5.1f}%)")
            if h.strand_normalized and h.strand_specific:
                print("  expect: stranded — RC'd should be ~all of one mate and ~none of the other")
            elif h.strand_normalized:
                print("  expect: unstranded — RC'd should be ~50% of each mate, one per pair")
            if n_paired_r1 != n_paired_r2:
                print(f"  [Warning] R1 and R2 counts differ by {abs(n_paired_r1 - n_paired_r2)}: "
                      "the file contains orphaned mates.")

        if compressed_payload > 0:
             print(f"Compressed Payload: {compressed_payload / (1024*1024):.2f} MB")
             print(f"Uncompressed Data:  {uncompressed_payload / (1024*1024):.2f} MB")
             ratio = uncompressed_payload / compressed_payload
             print(f"Compression Ratio:  {ratio:.2f}x")
        else:
             print(f"Data Size:          {uncompressed_payload / (1024*1024):.2f} MB")
             print(f"Compression Ratio:  1.00x")


def get_zna_version():
    try:
        from . import __version__
        return f"zna {__version__}"
    except ImportError:
        return "zna (unknown version)"


# --- COMMAND: SHUFFLE ---

def shuffle_command(args):
    """CLI wrapper for the bucket-shuffle algorithm."""
    start_time = time.time()
    buffer_bytes = parse_block_size(args.buffer_size)
    block_size = parse_block_size(args.block_size)

    try:
        written, n_records = shuffle_zna(
            args.input,
            args.output,
            seed=args.seed,
            buffer_bytes=buffer_bytes,
            block_size=block_size,
            tmp_dir=args.tmp_dir,
            quiet=args.quiet,
        )
    except (FileNotFoundError, ValueError) as exc:
        sys.exit(f"Error: {exc}")

    duration = time.time() - start_time
    if not args.quiet:
        print(
            f"\n[ZNA] Done. Shuffled {written:,} units ({n_records:,} records) "
            f"in {duration:.2f}s.",
            file=sys.stderr,
        )

# --- MAIN ---

def main():
    parser = argparse.ArgumentParser(description="zna: ZNA Processing Toolkit")
    parser.add_argument('--version', action='version', version=get_zna_version(), help="Show zna version and exit")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # --- ENCODE ---
    enc = subparsers.add_parser("encode", help="Convert FASTQ/FASTA to ZNA")
    enc.add_argument("files", nargs="*", help="Input files (0=stdin, 1=single/interleaved, 2=paired R1 R2)")
    
    input_group = enc.add_argument_group("Input Options")
    input_group.add_argument("--interleaved", action="store_true", help="Treat input as interleaved paired-end")
    input_group.add_argument("--shuffle", action="store_true", help="Shuffle records after encoding")
    input_group.add_argument("--seed", type=int, default=42, help="Seed for every random decision this command makes: "
                                    "--shuffle, unstranded strand normalization, and "
                                    "--npolicy random. All three are derived from it and "
                                    "the record index, never from a running stream, so "
                                    "the output does not depend on --block-size "
                                    "(default: 42)")
    input_group.add_argument(
        "--shuffle-buffer-size",
        type=str,
        default="1G",
        help="Max memory per bucket for encode --shuffle (default: 1G). Accepts K/M/G suffixes.",
    )
    input_group.add_argument("-q", "--quiet", action="store_true", help="Suppress progress messages")
    
    meta_group = enc.add_argument_group("Metadata")
    meta_group.add_argument("--read-group", default="Unknown", help="Read Group ID")
    meta_group.add_argument("--description", default="", help="Description string")
    meta_group.add_argument("--strand-specific", action="store_true",
                           help="Flag library as strand-specific (default: R1 antisense, R2 sense)")
    meta_group.add_argument("--strand-normalize", action="store_true",
                           help="Enable strand normalization (RC reads to consistent strand). "
                                "With --strand-specific: deterministic (antisense reads RC'd). "
                                "Without: random RC (for unstranded data).")
    meta_group.add_argument("--npolicy", choices=["trim3", "random"], default="trim3",
                           help="What to do with a base call that is not A/C/G/T. ZNA "
                                "stores two bits per base, so an Illumina 'N' no-call "
                                "cannot be stored. trim3 (default): cut the read at it, "
                                "keeping [0, first N) -- 3' only, so base 0 stays a true "
                                "fragment boundary however short the read gets. random: "
                                "substitute a base drawn from a seeded stream (--seed), "
                                "reproducible and identical on both backends. Any OTHER "
                                "non-ACGT character (an IUPAC code, punctuation) is an "
                                "error under both policies.")
    meta_group.add_argument("--treat-unpaired-as-merged", dest="treat_unpaired_as_merged",
                           action="store_true",
                           help="Unpaired records span their whole fragment, so BOTH edges "
                                "are true fragment boundaries (IS_FULL_FRAGMENT). Use for "
                                "overlap-merged reads. Default off: an unpaired record is "
                                "assumed to have one real edge, like a single-end read. "
                                "Full-overlap PAIRS are detected automatically and need no flag.")
    
    # R1 strand orientation (mutually exclusive)
    # Default: R1 is antisense (dUTP protocol)
    r1_strand = meta_group.add_mutually_exclusive_group()
    r1_strand.add_argument("--read1-sense", dest="read1_sense", action="store_true",
                          help="Read 1 represents sense strand")
    r1_strand.add_argument("--read1-antisense", dest="read1_sense", action="store_false",
                          help="Read 1 represents antisense strand (default for --strand-specific)")
    enc.set_defaults(read1_sense=False)  # Default: antisense (dUTP)
    
    # R2 strand orientation (mutually exclusive)
    # Default: R2 is sense (dUTP protocol)
    r2_strand = meta_group.add_mutually_exclusive_group()
    r2_strand.add_argument("--read2-sense", dest="read2_antisense", action="store_false",
                          help="Read 2 represents sense strand (default for --strand-specific)")
    r2_strand.add_argument("--read2-antisense", dest="read2_antisense", action="store_true",
                          help="Read 2 represents antisense strand")
    enc.set_defaults(read2_antisense=False)  # Default: sense (dUTP)
    
    mp_group = enc.add_argument_group(
        "Overlap merging (--merge-pairs)",
        "Merge overlapping pairs in process while encoding, deleting the "
        "FASTQ intermediate between `zna merge` and `zna encode`. Takes "
        "exactly two positional FASTQ files (R1 R2, adjacent on the command "
        "line). Pairing is handed over by the merger rather than re-derived "
        "from read names -- and mate names are CHECKED per pair, so two "
        "files that do not correspond now fail loudly where the two-file "
        "reader silently mispaired them (see --no-sync-check).")
    mp_group.add_argument("--merge-pairs", action="store_true",
                          help="merge overlapping pairs in process while encoding")
    from .merge.args import add_merge_algorithm_arguments
    add_merge_algorithm_arguments(mp_group)
    mp_group.add_argument("--merge-json", metavar="PATH",
                          help="write the merge run's statistics as JSON -- the "
                               "same block `zna merge --json` writes")
    mp_group.add_argument("--merge-backend", choices=("auto", "accel", "python"),
                          default="auto",
                          help="merge kernel. `auto` uses the compiled backend "
                               "and fails if it is missing; `python` is the "
                               "~50x-slower reference oracle, never chosen "
                               "for you.")
    mp_group.add_argument("--allow-empty", action="store_true",
                          help="permit an input with no read pairs. Off by "
                               "default: a 0-record .zna otherwise vanishes "
                               "from a corpus with every stage green.")

    fmt_group = enc.add_argument_group("Format Options")
    fmt_group.add_argument("-o", "--output", help="Output file. Defaults to stdout (-).")
    fmt_group.add_argument("--seq-len-bytes", type=int, choices=[1, 2, 4], default=2, help="Bytes for seq len")
    fmt_group.add_argument("--block-size", type=str, default="4M",
                           help="Block size (e.g. 512K, 4M, 8M). Larger = better compression. Default: 4M")
    
    comp_group = fmt_group.add_mutually_exclusive_group()
    comp_group.add_argument("--zstd", dest="compress_flag", action="store_true", default=None, help="Force ZSTD")
    comp_group.add_argument("--uncompressed", dest="compress_flag", action="store_false", help="Force uncompressed")
    fmt_group.add_argument("--level", type=int, default=DEFAULT_ZSTD_LEVEL, help="ZSTD level")

    label_group = enc.add_argument_group("Labels")
    label_group.add_argument("--label-defs", metavar="FILE",
                             help="Load label definitions from a YAML file")
    label_group.add_argument("--label", action="append", metavar="NAME:TYPE[:TAG]",
                             help="Define a label column (repeatable). NAME is stored in ZNA, "
                                  "TYPE is the dtype (A, c, C, s, S, i, I, f, d, q, Q). "
                                  "Optional TAG is the SAM tag to parse from input (defaults to NAME).")
    label_group.add_argument("--label-desc", action="append", metavar="TAG:TEXT",
                             help="Optional description for label TAG (max 64 chars)")

    # --- DECODE ---
    dec = subparsers.add_parser("decode", help="Convert ZNA to FASTA")
    dec.add_argument("input", nargs="?", help="Input .zna file (default: stdin)")
    dec.add_argument("-q", "--quiet", action="store_true", help="Suppress logs")
    
    # Combined Output Logic
    dec.add_argument("-o", "--output", help="Output filename. Use '#' for split files (e.g. 'out#.fa'). Ends in .gz for gzip.")
    dec.add_argument("--gzip", action="store_true", help="Force gzip compression (useful for stdout)")
    dec.add_argument("--restore-strand", action="store_true", 
                    help="Restore original strand orientation for antisense reads (reverse complement)")
    dec.add_argument("--labels", action="store_true",
                    help="Include label values in output as SAM-style tags")

    # --- INSPECT ---
    insp = subparsers.add_parser("inspect", help="Show ZNA file statistics")
    insp.add_argument("input", help="Input .zna file")
    insp.add_argument("--verify", action="store_true",
                      help="Fully certify the file: decode every block, recount "
                           "every stat, and compare against its trailer. Exit 0 "
                           "iff the file passes. A bare inspect already runs the "
                           "cheap structural checks; this is the paid tier.")
    insp.add_argument("--json", action="store_true",
                      help="Emit machine-readable JSON instead of text. Includes "
                           "n_records and n_blocks, read from block headers "
                           "without decompressing anything -- fast enough to "
                           "catalogue a whole corpus for manifest building.")
    insp.add_argument("--blocks", action="store_true",
                      help="With --json, include a per-block array (offset, "
                           "n_records, sizes).")
    insp.add_argument("--counts", action="store_true",
                      help="Report per-flag record counts (paired R1/R2, single/merged, "
                           "reverse-complemented). Reads block payloads, so slower than the "
                           "default header-only scan.")

    # --- SHUFFLE ---
    shuf = subparsers.add_parser("shuffle", help="Shuffle records in a ZNA file")
    shuf.add_argument("input", help="Input .zna file")
    shuf.add_argument("-o", "--output", required=True, help="Output .zna file")
    shuf.add_argument("-s", "--seed", type=int, default=42,
                      help="Random seed for reproducibility (default: 42)")
    shuf.add_argument("-b", "--buffer-size", type=str, default="1G",
                      help="Max memory per bucket (default: 1G). Accepts K/M/G suffixes.")
    shuf.add_argument("--block-size", type=str, default="4M",
                      help="Block size for output ZNA (default: 4M).")
    shuf.add_argument("--tmp-dir", type=str, default=None,
                      help="Directory for temporary bucket files (default: system temp)")
    shuf.add_argument("-q", "--quiet", action="store_true", help="Suppress progress messages")

    # --- MERGE ---
    # Registered from zna.merge.args, which imports nothing but argparse. The runtime
    # half (zna.merge.cli) pulls in the kernel, the extension module and a 64 KiB
    # lookup table -- paying that on `zna inspect` would be absurd, so it is imported
    # below, only once the user has actually asked for `merge`.
    from .merge.args import add_merge_parser
    add_merge_parser(subparsers)

    args = parser.parse_args()
    try:
        if args.command == "encode":
            encode_command(args)
        elif args.command == "decode":
            decode_command(args)
        elif args.command == "inspect":
            inspect_command(args)
        elif args.command == "shuffle":
            shuffle_command(args)
        elif args.command == "merge":
            from .merge.cli import run_command
            sys.exit(run_command(args))
    except ValueError as e:
        # Everything the format rejects -- an unsupported version, a base outside ACGTN,
        # a record longer than seq_len_bytes allows, a half-written fragment -- arrives
        # as a ValueError carrying a message written for a person to read.  A traceback
        # buries that message under a stack the user cannot act on.  `inspect` and
        # `shuffle` already caught their own; this is what `decode` and `encode` were
        # missing, and the version error most of all, since that is the one every pre-3
        # file produces on first contact.
        sys.exit(f"Error: {e}")


if __name__ == "__main__":
    main()