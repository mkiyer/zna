"""
Label dtype system for ZNA V2.

Provides the type registry for per-sequence label columns, using
SAM/BAM auxiliary tag type codes as the specification basis.

Supported numeric types:
- BAM-compatible: A, c, C, s, S, i, I, f
- ZNA extensions: d, q, Q  (8-byte)
"""
from __future__ import annotations

import struct as _struct
from dataclasses import dataclass


@dataclass(slots=True, frozen=True)
class LabelDtype:
    """Immutable descriptor for a label storage type."""

    code: str        # single ASCII char — type code on disk
    name: str        # human-readable name, e.g. "uint8", "float32"
    size: int        # bytes per value
    struct_ch: str   # Python struct pack/unpack char
    is_numeric: bool


@dataclass(slots=True)
class LabelDef:
    """Definition of one label column in the ZNA file header."""

    label_id: int       # 0-based index (== position in the labels array)
    name: str           # label name stored in ZNA, max 16 UTF-8 bytes
    description: str    # human-readable description, max 64 UTF-8 bytes
    dtype: LabelDtype   # from the dtype registry
    missing: int | float | None = None  # explicit missing sentinel (None → dtype default)
    tag: str | None = None  # input file tag to match at encode time (not stored in ZNA)

    @property
    def effective_tag(self) -> str:
        """Return the tag used for parsing input headers.

        Falls back to ``name`` when ``tag`` is not set, which preserves
        backward compatibility with label definitions that do not
        distinguish between name and tag.
        """
        return self.tag if self.tag is not None else self.name


# ---------------------------------------------------------------------------
# Dtype registry
# ---------------------------------------------------------------------------

_DTYPES: list[LabelDtype] = [
    # BAM-compatible numeric types
    LabelDtype("A", "char",    1, "B", True),
    LabelDtype("c", "int8",    1, "b", True),
    LabelDtype("C", "uint8",   1, "B", True),
    LabelDtype("s", "int16",   2, "h", True),
    LabelDtype("S", "uint16",  2, "H", True),
    LabelDtype("i", "int32",   4, "i", True),
    LabelDtype("I", "uint32",  4, "I", True),
    LabelDtype("f", "float32", 4, "f", True),
    # ZNA extensions (8-byte)
    LabelDtype("d", "float64", 8, "d", True),
    LabelDtype("q", "int64",   8, "q", True),
    LabelDtype("Q", "uint64",  8, "Q", True),
]

DTYPE_BY_CODE: dict[str, LabelDtype] = {dt.code: dt for dt in _DTYPES}
DTYPE_BY_NAME: dict[str, LabelDtype] = {dt.name: dt for dt in _DTYPES}

# Default missing value per dtype code (all zeros)
_DTYPE_DEFAULT_MISSING: dict[str, int | float] = {
    "A": 0, "c": 0, "C": 0,
    "s": 0, "S": 0, "i": 0, "I": 0,
    "f": 0.0, "d": 0.0,
    "q": 0, "Q": 0,
}

# Struct for packing the 8-byte missing slot in the file header
_MISSING_STRUCT: dict[str, _struct.Struct] = {
    code: _struct.Struct(f"<{dt.struct_ch}") for code, dt in DTYPE_BY_CODE.items()
}


def resolve_missing(ldef: LabelDef) -> int | float:
    """Return the effective missing value for *ldef*.

    Uses the explicit ``missing`` field if set, otherwise the dtype default
    (zero for all types).
    """
    if ldef.missing is not None:
        return ldef.missing
    return _DTYPE_DEFAULT_MISSING[ldef.dtype.code]


def pack_missing(ldef: LabelDef) -> bytes:
    """Pack the resolved missing value into an 8-byte slot (little-endian).

    The value is packed in the dtype's native format and zero-padded to 8 bytes.
    """
    val = resolve_missing(ldef)
    native = _MISSING_STRUCT[ldef.dtype.code].pack(val)
    return native.ljust(8, b'\x00')


def unpack_missing(dtype: LabelDtype, data: bytes) -> int | float:
    """Unpack a missing value from an 8-byte header slot."""
    st = _MISSING_STRUCT[dtype.code]
    return st.unpack(data[:st.size])[0]


def parse_dtype(s: str) -> LabelDtype:
    """Resolve a dtype string.

    Accepts type codes (``'C'``, ``'i'``) or human-readable
    names (``'uint8'``, ``'int32'``).

    Raises ``ValueError`` if *s* is not recognised.
    """
    if s in DTYPE_BY_CODE:
        return DTYPE_BY_CODE[s]
    if s in DTYPE_BY_NAME:
        return DTYPE_BY_NAME[s]
    raise ValueError(f"Unknown label dtype: {s!r}")


def label_bytes_per_record(ldef: LabelDef) -> int:
    """Return the per-record byte width for *ldef*."""
    return ldef.dtype.size
