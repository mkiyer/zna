"""Read-name helpers.

Their own module because both the reference backend and :mod:`zna.merge.pairs` need
them, and neither should have to import the other to get them. Nothing heavy is
imported here, so this stays free to pull in from anywhere.

Headers arrive without the leading ``@`` and without the trailing newline.
"""
from __future__ import annotations


def id_end(header: bytes) -> int:
    """Index of the first whitespace (space or tab), or ``len(header)`` if none.

    Uses C-level ``bytes.find`` rather than a Python per-char loop; this runs on every
    read.
    """
    sp = header.find(b" ")
    tab = header.find(b"\t")
    if sp == -1:
        return tab if tab != -1 else len(header)
    if tab == -1:
        return sp
    return sp if sp < tab else tab


def base_name(header: bytes) -> bytes:
    """Read ID without the ``/1``/``/2`` pair suffix and without comments/tags.

    Mirrors ZNA's ``get_base_name`` so the R1/R2 sync check matches ZNA's pairing.
    """
    idtok = header[:id_end(header)]
    slash = idtok.rfind(b"/")
    return idtok[:slash] if slash != -1 else idtok


def strip_pair_suffix(header: bytes) -> bytes:
    """Drop a trailing ``/1`` or ``/2`` from the ID token, preserving tags/separator."""
    cut = id_end(header)
    idtok, rest = header[:cut], header[cut:]
    if idtok.endswith(b"/1") or idtok.endswith(b"/2"):
        idtok = idtok[:-2]
    return idtok + rest
