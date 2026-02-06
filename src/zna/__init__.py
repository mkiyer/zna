__version__ = "0.1.6"


from .core import (
    ZnaHeader, ZnaWriter, ZnaReader,
    read_zna, write_zna,
    ZnaHeaderFlags, ZnaRecordFlags,
    reverse_complement,
)
from .codec import get_backend, get_backend_name, available_backends


def is_accelerated() -> bool:
    """Check if the C++ acceleration backend is active."""
    return get_backend_name() == "accel"


__all__ = [
    "ZnaHeader", "ZnaWriter", "ZnaReader",
    "read_zna", "write_zna",
    "ZnaHeaderFlags", "ZnaRecordFlags",
    "reverse_complement",
    "is_accelerated",
    "get_backend", "get_backend_name", "available_backends",
    "__version__",
]
