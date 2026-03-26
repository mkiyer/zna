__version__ = "0.3.0"


from .core import (
    ZnaHeader, ZnaWriter, ZnaReader,
    read_zna, write_zna,
    ZnaHeaderFlags, ZnaRecordFlags,
    reverse_complement,
)
from .dtypes import LabelDef, LabelDtype, parse_dtype, resolve_missing
from .codec import get_backend, get_backend_name, available_backends


def is_accelerated() -> bool:
    """Check if the C++ acceleration backend is active."""
    return get_backend_name() == "accel"


__all__ = [
    "ZnaHeader", "ZnaWriter", "ZnaReader",
    "read_zna", "write_zna",
    "ZnaHeaderFlags", "ZnaRecordFlags",
    "reverse_complement",
    "LabelDef", "LabelDtype", "parse_dtype", "resolve_missing",
    "is_accelerated",
    "get_backend", "get_backend_name", "available_backends",
    "__version__",
]
