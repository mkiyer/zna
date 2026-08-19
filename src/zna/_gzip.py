"""One place that decides how a ``.gz`` gets inflated.

Three ways to inflate:

1. **ISA-L**, via the optional :mod:`isal` package. In process, releasing the GIL.
2. **An external ``pigz``/``gzip``**, in its own process.
3. **The stdlib :mod:`gzip` module**, always available.

Raw inflate rate, this repository's 1M-pair library (106 MB gzip member, 327 MB out),
Xeon E5-2680 v3, one core:

===================================  ========  ============
path                                 wall      inflate rate
===================================  ========  ============
stdlib ``gzip.open``                 1.52 s    208 MB/s
``pigz -dc -p1`` in a subprocess      1.63 s    193 MB/s
``isal.igzip.open``                  0.70 s    448 MB/s
===================================  ========  ============

``pigz`` looks bad there and is still worth having: it cannot thread *inflate* (only the
CRC and read-ahead), so on one gzip member it is plain single-stream zlib — but it runs
in **another process**, so its cost overlaps with the caller's instead of adding to it.

**Which of those two matters depends on the caller, and the two ZNA read paths differ.**
That is not a preference; it was measured after assuming otherwise:

============================  ==============  ==============  ===================
1M pairs, wall / user         ISA-L           pigz            faster
============================  ==============  ==============  ===================
``zna merge --threads 1``     **4.04 / 4.28** 5.84 / 6.06     ISA-L, on both axes
``zna encode`` from ``.gz``   6.54 / 6.31     **5.50 / 7.26** pigz on wall
============================  ==============  ==============  ===================

``zna merge`` reads through :class:`zna.merge.cli._RawStream`, which prefetches on its
own thread; ISA-L releases the GIL, so an in-process inflate overlaps there exactly as a
subprocess would, and being 2.3x faster it wins outright — 1.45x wall *and* 1.42x CPU.
``zna encode`` has no reader thread: it drives ``readline`` from the main thread, so its
inflate competes with a GIL-bound parse loop, and moving that work to another core beats
making it cheaper. ISA-L still uses 1.15x less CPU there, which is why it is the choice
on a machine with no ``pigz`` — where its real competition is the 208 MB/s stdlib.

So: :func:`prefer_isal` takes the one fact that decides it, and each caller passes what
is true of itself.

ISA-L is **optional and never required**. It is not imported until a ``.gz`` is actually
opened, so ``import zna`` costs nothing for it, and a machine without it behaves exactly
as 0.5.1 did. ``pip install 'zna[fast]'`` or ``pip install isal``.

Two environment escapes, for debugging a suspected difference between decompressors:
``ZNA_NO_ISAL=1`` skips ISA-L and ``ZNA_NO_EXTERNAL_GZIP=1`` skips the subprocess;
both set forces the stdlib.
"""
from __future__ import annotations

import gzip
import io
import os
from typing import Any, Optional

#: Cached result of probing for :mod:`isal`. ``None`` until the first probe, then the
#: module or ``False``. Probed once per process, on first use, never at import.
_ISAL: Any = None


def isal_module() -> Optional[Any]:
    """``isal.igzip`` if it is importable and not disabled, else ``None``."""
    global _ISAL
    if _ISAL is None:
        if os.environ.get("ZNA_NO_ISAL"):
            _ISAL = False
        else:
            try:
                from isal import igzip
            except Exception:
                # Any failure at all -- not installed, wrong ABI, a broken build --
                # means "use the next path", never a traceback out of an open().
                _ISAL = False
            else:
                _ISAL = igzip
    return _ISAL or None


def external_gzip_allowed() -> bool:
    """Whether an external decompressor may be spawned."""
    return not os.environ.get("ZNA_NO_EXTERNAL_GZIP")


def prefer_isal(*, own_read_thread: bool) -> bool:
    """Should ISA-L be tried before spawning an external decompressor?

    *own_read_thread* says whether the caller reads on a dedicated thread. If it does,
    an in-process inflate that releases the GIL overlaps with the caller's work just as
    a subprocess would, so the faster implementation wins on both wall and CPU. If it
    does not — a ``readline`` loop on the main thread — the subprocess's separate core
    beats ISA-L's lower cost on wall time, and ISA-L is better kept as the fallback for
    when no external tool exists. The module docstring has both measurements.

    ISA-L is always used ahead of the stdlib either way; this only orders it against
    the subprocess.
    """
    if not own_read_thread and external_gzip_allowed():
        return False
    return True


def inflate_backend_name(path: str, *, own_read_thread: bool = True) -> str:
    """Which path a caller will take for *path*: isal, external, gzip, or plain.

    For logs and provenance: a throughput number is not interpretable without knowing
    which decompressor produced it.
    """
    if not str(path).endswith(".gz"):
        return "plain"
    if prefer_isal(own_read_thread=own_read_thread) and isal_module() is not None:
        return "isal"
    if external_gzip_allowed():
        return "external"
    if isal_module() is not None:
        return "isal"
    return "gzip"


def open_isal(path: str, buffer_size: int):
    """Open *path* through ISA-L, wrapped so reads come in large blocks.

    Returns ``None`` if ISA-L is unavailable, which is the caller's signal to try the
    next path.

    ``igzip.open`` takes no ``buffering`` argument -- it mirrors ``gzip.open``, which
    does not either -- so it is wrapped the same way the stdlib path is: the callers
    here read either whole 256 KiB blocks or four lines per record, and both want one
    wide buffer rather than ``IGzipFile``'s own general-purpose one.
    """
    igzip = isal_module()
    if igzip is None:
        return None
    return io.BufferedReader(igzip.open(path, "rb"), buffer_size=buffer_size)


def open_stdlib_gzip(path: str, buffer_size: int):
    """Open *path* through the stdlib, with a wide buffer in front.

    A FASTQ record is four ``readline`` calls and ``GzipFile`` serves each from its own
    modest internal buffer; the wrapper amortises that across many records instead of
    paying it per line.
    """
    return io.BufferedReader(gzip.open(path, "rb"), buffer_size=buffer_size)
