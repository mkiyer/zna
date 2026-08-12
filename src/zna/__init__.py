__version__ = "0.4.0"


from .core import (
    ZnaHeader, ZnaWriter, ZnaReader, BlockInfo,
    read_zna, write_zna,
    ZnaHeaderFlags, ZnaRecordFlags, FLAG_FIELDS, ENDS_BY_FLAG,
    reverse_complement,
)
from .dtypes import LabelDef, LabelDtype, parse_dtype, resolve_missing
from .codec import get_backend, get_backend_name, available_backends


def is_accelerated() -> bool:
    """Check if the C++ acceleration backend is active."""
    return get_backend_name() == "accel"


__all__ = [
    "ZnaHeader", "ZnaWriter", "ZnaReader", "BlockInfo",
    "read_zna", "write_zna",
    "ZnaHeaderFlags", "ZnaRecordFlags", "FLAG_FIELDS", "ENDS_BY_FLAG",
    "reverse_complement",
    "LabelDef", "LabelDtype", "parse_dtype", "resolve_missing",
    "is_accelerated",
    "get_backend", "get_backend_name", "available_backends",
    "__version__",
]
