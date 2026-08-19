"""FASTQ I/O: paired reader and a single-stream writer.

Gzip is handled by ``pigz`` (fast, multithreaded) when available, falling back to the
stdlib ``gzip`` module, or plain files for non-``.gz`` paths. Records are ``bytes``
triples ``(header, seq, qual)`` — header without the leading ``@`` or trailing newline;
sequences upper-cased for correct overlap comparison.
"""
from __future__ import annotations

import gzip
import shutil
import subprocess

from .. import _gzip

_BUF = 1 << 20


class InputError(ValueError):
    """Malformed or inconsistent input: truncation, desync, unequal read counts.

    Distinct from a bare ``ValueError`` so the CLI can turn *these* into a clean
    ``SystemExit`` with a message while letting a genuine bug keep its traceback.
    """


def _open_binary_read(path: str, threads: int, *, own_read_thread: bool):
    """Return (stream, proc). ``proc`` is the pigz process to reap, or None.

    *own_read_thread* must say whether the caller reads on a dedicated thread, because
    that is what decides between ISA-L and a subprocess and **the two callers of this
    function differ**. :class:`zna.merge.cli._RawStream` prefetches on its own thread and
    wants ISA-L; :func:`read_fastq` below drives ``readline`` from the calling thread and
    wants the subprocess. Getting that backwards is not subtle — it cost `read_fastq`
    2.836 -> 4.007 us/pair when this function briefly hardcoded ``True``. The policy and
    both measurements are in :mod:`zna._gzip`.

    ISA-L is worth wanting where it fits: 2.3x the pigz pipe on this repository's
    1M-pair library, and a profile of ``zna merge --threads 1`` puts a quarter of all
    cycles inside pigz's ``libz``, which makes inflate the largest single cost of a merge.

    pigz cannot parallelise inflate (it threads only the CRC and read-ahead), so a
    reader gains almost nothing from extra threads — and measured, ``-p4`` per stream is
    *worse* than ``-p1`` because the threads contend with the merge workers. The reader
    therefore always passes 1. ``--io-threads`` sizes the *writer* instead, where
    deflate genuinely parallelises.
    """
    if path.endswith(".gz"):
        if _gzip.prefer_isal(own_read_thread=own_read_thread):
            stream = _gzip.open_isal(path, _BUF)
            if stream is not None:
                return stream, None
        pigz = shutil.which("pigz") if _gzip.external_gzip_allowed() else None
        if pigz:
            proc = subprocess.Popen(
                [pigz, "-dc", "-p", str(max(1, threads)), path],
                stdout=subprocess.PIPE, bufsize=_BUF,
            )
            return proc.stdout, proc
        stream = _gzip.open_isal(path, _BUF)
        if stream is not None:
            return stream, None
        return gzip.open(path, "rb"), None
    return open(path, "rb", buffering=_BUF), None


def read_fastq(path, threads: int = 1):
    """Yield ``(header, seq, qual)`` bytes triples from a FASTQ (optionally gzipped).

    Not on the merge path — the backend parses raw buffers itself (see
    ``merge_chunk``). This is here for the development tools under
    ``scripts/merge_bench/``, which want records rather than bytes.
    """
    path = str(path)
    # This loop is the reader: `readline` runs on the calling thread, so an in-process
    # inflate would compete with it rather than overlap. See `_open_binary_read`.
    stream, proc = _open_binary_read(path, threads, own_read_thread=False)
    readline = stream.readline
    try:
        while True:
            h = readline()
            if not h:
                break
            s = readline()
            readline()          # '+' separator
            q = readline()
            if not q:
                raise InputError(f"truncated FASTQ record in {path}")
            if h[0] != 0x40:    # b'@'
                raise InputError(f"malformed FASTQ header in {path}: {h[:60]!r}")
            s = s.rstrip(b"\r\n").upper()
            q = q.rstrip(b"\r\n")
            # A file truncated inside its LAST quality line still satisfies `if not q`,
            # so without this the record is emitted with a short quality string and the
            # tool exits 0. Every earlier truncation is already caught above.
            if len(s) != len(q):
                raise InputError(
                    f"FASTQ record in {path} has {len(s)} bases but {len(q)} quality "
                    f"scores (truncated or malformed): {h[1:60]!r}")
            yield h[1:].rstrip(b"\r\n"), s, q
    finally:
        try:
            stream.close()
        except Exception:
            pass
        if proc is not None:
            proc.wait()
            # Positive exit = real pigz error (missing/corrupt). Negative = killed by a
            # signal (e.g. -13 SIGPIPE when the consumer stops early) — not an error.
            if proc.returncode and proc.returncode > 0:
                raise IOError(f"pigz failed reading {path} (exit {proc.returncode})")


def read_pairs(path1, path2, threads: int = 1):
    """Yield ``(r1, r2)`` record pairs from two positionally-synced FASTQ files."""
    it2 = read_fastq(path2, threads)
    for r1 in read_fastq(path1, threads):
        r2 = next(it2, None)
        if r2 is None:
            raise InputError("R2 exhausted before R1 (unequal read counts)")
        yield r1, r2
    if next(it2, None) is not None:
        raise InputError("R1 exhausted before R2 (unequal read counts)")


class FastqWriter:
    """Write FASTQ records to one stream (gzipped via pigz when the path ends .gz)."""

    def __init__(self, path, threads: int = 1, level: int = 1):
        self.path = str(path)
        self._proc = None
        self._sink = None
        if self.path.endswith(".gz"):
            pigz = shutil.which("pigz")
            if pigz:
                self._sink = open(self.path, "wb")
                self._proc = subprocess.Popen(
                    [pigz, "-c", "-p", str(max(1, threads)), f"-{level}"],
                    stdin=subprocess.PIPE, stdout=self._sink, bufsize=_BUF,
                )
                self._fh = self._proc.stdin
            else:
                self._fh = gzip.open(self.path, "wb", compresslevel=level)
        else:
            self._fh = open(self.path, "wb", buffering=_BUF)

    def write(self, header: bytes, seq: bytes, qual: bytes) -> None:
        self._fh.write(b"@%b\n%b\n+\n%b\n" % (header, seq, qual))

    def write_raw(self, blob: bytes) -> None:
        """Write a pre-formatted FASTQ byte blob (used by the parallel path)."""
        self._fh.write(blob)

    def close(self) -> None:
        self._fh.close()
        if self._proc is not None:
            self._proc.wait()
            if self._sink is not None:
                self._sink.close()
            if self._proc.returncode:
                raise IOError(f"pigz failed writing {self.path} (exit {self._proc.returncode})")

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
