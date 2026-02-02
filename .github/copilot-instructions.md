# ZNA (Compressed Nucleic Acid) Project Instructions

## Project Overview
ZNA is a specialized binary format project for compressing and storing nucleic acid sequences. It focuses on high-performance I/O and storage efficiency using 2-bit encoding for DNA bases (A, C, G, T). The project is a pure Python implementation requiring no external runtime dependencies.

## Architecture & Core Components

### 1. Core Logic (`src/zna/core.py`)
This is the heart of the library, handling all binary serialization and deserialization.
-   **`ZnaHeader`**: definition of the file header metadata using `@dataclass(slots=True)`.
-   **`ZnaWriter`**: Handles encoding text sequences into the binary format.
-   **`ZnaReader`**: Handles parsing the binary format back into text sequences.
-   **Lookup Tables**: Uses pre-computed tables (`_CHAR_TO_2BIT`, `_DECODE_TABLE`) for fast translation between characters and 2-bit integers.

### 2. CLI (`src/zna/cli.py`)
-   Entry point for the application (`zna` command).
-   **`encode`**: Stream-based processing reading raw sequences from stdin.
-   **`decode`**: Reads .zna files and outputs raw sequences to stdout.

## Binary Format Specification
When working on `core.py`, strictly adhere to the defined binary structure:
-   **Magic Bytes**: `b"ZNA\x00"` (4 bytes)
-   **Endpoints**: Little-endian (`<`) is used for all struct packing.
-   **Header Structure**:
    -   Format: `<4sBBBHHH` (Magic, Ver, Len, Flags, 3xStrLen)
    -   Followed by variable-length strings for `read_group`, `description`, `extra_info`.
-   **Encoding**: Bases map to 2-bit integers (A=0, C=1, G=2, T=3). These are bit-packed into bytes.

## Development Workflows

### Build & Dependency Management
-   **Build System**: `flit`. The project uses `pyproject.toml` exclusively (no `setup.py`).
-   **Dependencies**: Zero runtime dependencies. Keep it that way.
-   **Dev Dependencies**: `pytest`, `black`, `mypy`.

### Testing
-   Run tests using `pytest`.
-   Tests are located in `tests/`.

## Coding Conventions
-   **Type Hinting**: Extensive use of `typing` (`BinaryIO`, `Iterable`, `Iterator`) and `__future__.annotations` is required.
-   **Performance**: Prefer pre-computed lookup tables over runtime logic for hot paths (like base encoding).
-   **Dataclasses**: Use `slots=True` for data structures to minimize memory overhead.
-   **I/O**: Always use binary mode (`"rb"`, `"wb"`) when dealing with ZNA files.

## Common Patterns
### bit-packing
When implementing new encoding logic, ensure 2-bit values are packed into integers before writing to disk to maximize compression.

### CLI piping
The CLI is designed for unix-style pipes.
-   Input: `cat sequences.txt |zna encode -o out.zna`
-   Output: `zna decode -o out.zna > decoded.txt`
