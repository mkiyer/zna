# ZNA project instructions

ZNA is a binary container for DNA/RNA sequence: 2-bit packed bases in a columnar,
block-framed file, with optional Zstd compression and per-record numeric labels. It ships
two C++ extensions with pure-Python fallbacks, and a CLI (`zna encode|decode|inspect|shuffle|merge`).

Deeper references: [`docs/METHODS.md`](../docs/METHODS.md) for how the format and the
merge kernel actually work, [`docs/ROADMAP.md`](../docs/ROADMAP.md) for what has already
been tried and closed by measurement, [`CHANGELOG.md`](../CHANGELOG.md) for what shipped.

## Layout

| Path | What it is |
|---|---|
| `src/zna/core.py` | The format: headers, block framing, `ZnaWriter`, `ZnaReader` |
| `src/zna/_pycodec.py` | Reference codec — the specification, in Python |
| `src/zna/_accel.cpp` | Compiled codec (`zna._accel`), must match `_pycodec` byte for byte |
| `src/zna/codec.py` | Backend selection (`python` / `accel`) |
| `src/zna/cli.py` | CLI for encode/decode/inspect/shuffle |
| `src/zna/_shuffle.py` | Bucket shuffle with bounded memory |
| `src/zna/merge/` | `zna merge`: overlap merger, its own C++ kernel + `_pymerge.py` reference |
| `src/zna/dtypes.py` | Label dtypes |

## Build and test

```bash
pip install -e . --no-build-isolation     # scikit-build-core + nanobind + CMake
python -m pytest tests -q
```

`--no-build-isolation` is load-bearing: without it scikit-build-core rewrites
`CMakeCache` to point at a directory pip then deletes, and alternating between the two
forms costs a full rebuild each way. **A CMake failure ends with an unremarkable pip
banner and leaves the old `.so` installed**, so the suite can run green against stale
code — grep the output for `Successfully built`.

Runtime dependencies are `zstandard` and `pyyaml`; both are declared in
`pyproject.toml` and in `conda/meta.yaml`'s `run:`. Adding a third is a real decision,
not a convenience — this is a library other pipelines install.

The compiled extensions are optional at runtime and about 50x faster. Either can be
absent without an error, so a broken build looks like a slow machine:

```python
import zna; from zna.merge.backend import available_merge_backends
zna.is_accelerated()          # codec extension
available_merge_backends()    # want 'accel' in here
```

## Format

Little-endian throughout. Magic `b"ZNA\x1A"`, **format version 3**.

```
file header    <4sBBBBBHHH>  magic | version | seq_len_bytes | header flags |
                             compression method | level | label count |
                             read-group length | description length
               then the read group, the description, and one 89-byte label
               definition per label
block header   <IIIII>       comp_size | uncomp_size | record_count |
                             flags_size | lengths_size
block payload  one zstd frame:  flags ‖ labels ‖ lengths ‖ sequences
```

Per-record flags (`ZnaRecordFlags`): `IS_READ1 1`, `IS_READ2 2`, `IS_PAIRED 4`,
`IS_RC 8`, `IS_FULL_FRAGMENT 16`.

The payload is columnar because that is what makes it compress: each stream is
homogeneous, so zstd sees runs instead of an interleave. Flags come first so a consumer
can decompress a *prefix* of the frame and read the flag column without touching
sequence.

**No record or block count is stored in the file header.** Every block header carries
its own; totals come from walking the chain (`block_index()`).

## Invariants that must not be broken

1. **Blocks are fragment-complete.** A fragment's reads are consecutive, R1 immediately
   followed by R2, and never span a block boundary. `ZnaWriter` enforces this on every
   write path and raises otherwise. It is what makes `blocks(stride=…)` sharding sound.
2. **The two backends produce byte-identical output.** `_pycodec.py` is the
   specification; `_accel.cpp` must match it. `tests/test_fuzz_roundtrip.py` runs the
   whole matrix against both. One divergence is currently open and pinned by a test —
   see ROADMAP.
3. **Output must not depend on `--block-size`.** Both random streams — the unstranded RC
   coin and `--npolicy random`'s substitution — are derived from the record's *global*
   index via `_mix64`, never from a running per-block state.
4. **A view is for reading; the flag byte is for copying.** `records(with_ends=…)` and
   `with_rc=` are lossy projections for consumers and must never travel back into a
   writer. Copy with `copy_records()` → `write_copy()`, which carries the byte verbatim.
5. **Strand normalization is not idempotent.** Applying it twice un-normalizes the data
   while the header still claims otherwise. Every ZNA→ZNA path copies orientation rather
   than re-deriving it (`preserve_normalization=True`).

## Conventions

- `from __future__ import annotations`, and type hints throughout.
- Lookup tables over branches on the hot path (`ENDS_BY_FLAG`, `FLAG_FIELDS`,
  `OPENS_FRAGMENT`, `_DECODE_TABLE`) — these are per-record code paths.
- `@dataclass(slots=True)` for data structures; `__slots__` on the writer/reader.
- Binary mode (`"rb"`/`"wb"`) for ZNA files, always.
- Comments explain *why*, especially where a defect was fixed — several in `core.py`
  record the measurement that closed a question. Do not delete them to tidy up.
