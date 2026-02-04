# ZNA CLI Test Summary

## Overview
Comprehensive test suite created for the ZNA CLI tool with 22 tests covering all major functionality.

## Issues Fixed in CLI Code

### 1. Parser Selection Logic ([cli.py](src/zna/cli.py))
**Issue**: The `get_parser()` function was checking for 'fa' in filepath, which incorrectly matched 'fastq' files.

**Fix**: Changed to use precise extension matching:
```python
def get_parser(filepath: Optional[str]):
    if filepath:
        lower = filepath.lower()
        if lower.endswith('.fasta') or lower.endswith('.fa') or \
           lower.endswith('.fasta.gz') or lower.endswith('.fa.gz'):
            return parse_fasta
    return parse_fastq
```

### 2. FASTQ Parser ([cli.py](src/zna/cli.py))
**Issue**: Parser wasn't properly extracting sequence lines and would include FASTQ headers/quality lines.

**Fix**: Added proper header validation and sequence extraction:
```python
def parse_fastq(fh: BinaryIO) -> Iterator[str]:
    while True:
        header = fh.readline()
        if not header: break
        if not header.startswith(b'@'): continue
        seq = fh.readline().strip().decode('ascii', errors='ignore')
        plus = fh.readline()  # +
        qual = fh.readline()  # Quality
        if seq:
            yield seq
```

## Test Coverage

### Parser Tests (7 tests)
✅ FASTA parsing (simple, empty, single-line)
✅ FASTQ parsing (simple, empty)
✅ Parser selection (FASTA vs FASTQ)

### Encode Tests (5 tests)
✅ FASTA to uncompressed ZNA
✅ FASTQ to compressed ZNA (Zstd)
✅ Paired-end from separate files
✅ Interleaved paired-end input
✅ Gzipped input files

### Decode Tests (5 tests)
✅ Roundtrip uncompressed (sequences match)
✅ Roundtrip compressed (sequences match)
✅ Paired-end split output (R1/R2 to separate files)
✅ Gzipped output
✅ Record count consistency (1000 records)

### Inspect Tests (3 tests)
✅ Uncompressed file inspection
✅ Compressed file inspection
✅ Block statistics accuracy

### Integration Tests (2 tests)
✅ Full single-end workflow (FASTQ → ZNA → FASTA)
✅ Full paired-end workflow (gzipped inputs → compressed ZNA → split gzipped outputs)

## Test Results
```
32 tests total: 22 CLI tests + 10 core tests
All tests PASSED ✓
```

## Key Test Features

### Roundtrip Verification
Tests verify that sequences encoded and then decoded match the original input exactly.

### Compression Testing
- Tests both compressed (Zstd) and uncompressed modes
- Verifies compression levels are properly stored and retrieved
- Tests with various compression levels (3, 5)

### Format Support
- FASTA input parsing
- FASTQ input parsing
- Gzipped input (.gz)
- Gzipped output (.gz)
- Split paired-end files (# pattern)
- Interleaved paired-end

### Data Integrity
- Sequence accuracy verification
- Record count consistency
- Paired-end read pairing correctness (R1/R2 flags)
- Block statistics validation

## Code Quality Improvements

### CLI Robustness
- Fixed parser selection logic
- Improved FASTQ parsing reliability
- Better handling of edge cases (empty files, orphaned reads)

### Test Structure
- Organized into logical test classes
- Clear test names describing functionality
- Comprehensive edge case coverage
- Uses tempfile for clean test isolation

## Running Tests

```bash
# Run all tests
PYTHONPATH=src pytest -v

# Run only CLI tests
PYTHONPATH=src pytest tests/test_cli.py -v

# Run specific test class
PYTHONPATH=src pytest tests/test_cli.py::TestEncode -v

# Run with coverage
PYTHONPATH=src pytest --cov=zna tests/
```

## Notes

- No modifications needed to core.py (binary format logic)
- CLI code changes were minimal and focused on bug fixes
- All tests use temporary directories for isolation
- Tests cover realistic bioinformatics workflows
- Integration tests simulate complete user workflows
