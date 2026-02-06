"""
Pure Python codec backend for ZNA.

Implements DNA sequence encoding/decoding using 2-bit packing.
All functions are stateless — this module serves as a backend
that can be swapped with the C++ accelerator or future backends.
"""
from __future__ import annotations


# ---------------------------------------------------------------------------
# Lookup tables
# ---------------------------------------------------------------------------

_BASES = ("A", "C", "G", "T")

# Complement table: A<->T, C<->G (for str.translate)
_COMPLEMENT_TABLE = str.maketrans("ACGTacgt", "TGCAtgca")


def _make_encode_table() -> bytes:
    """Create translation table: ASCII byte → 2-bit integer (A=0 C=1 G=2 T=3).

    Invalid characters map to 255.
    """
    table = bytearray([255] * 256)
    for i, base in enumerate(_BASES):
        table[ord(base)] = i
        table[ord(base.lower())] = i
    return bytes(table)


_TRANS_TABLE = _make_encode_table()

# Decode table: maps a packed byte to its 4-character string.
_DECODE_TABLE = tuple(
    _BASES[(v >> 6) & 3]
    + _BASES[(v >> 4) & 3]
    + _BASES[(v >> 2) & 3]
    + _BASES[v & 3]
    for v in range(256)
)


# ---------------------------------------------------------------------------
# Public API — every backend must expose these four functions.
# ---------------------------------------------------------------------------


def encode_sequence(seq: str) -> bytes:
    """Encode a single DNA sequence to 2-bit packed bytes.

    Four bases are packed into each byte (A=00, C=01, G=10, T=11).
    Raises ``ValueError`` on characters outside ACGT/acgt.
    """
    bseq = seq.encode("ascii")
    mapped = bseq.translate(_TRANS_TABLE)
    length = len(mapped)

    out = bytearray((length + 3) // 4)

    i = 0
    idx = 0
    full_chunks = length // 4

    while idx < full_chunks:
        val = (mapped[i] << 6) | (mapped[i + 1] << 4) | (mapped[i + 2] << 2) | mapped[i + 3]
        if val > 255:
            # One of the mapped values was 255 (invalid character)
            raise ValueError("Invalid character in sequence")
        out[idx] = val
        idx += 1
        i += 4

    if i < length:
        val = 0
        shift = 6
        while i < length:
            if mapped[i] == 255:
                raise ValueError("Invalid character in sequence")
            val |= mapped[i] << shift
            i += 1
            shift -= 2
        out[idx] = val

    return bytes(out)


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of *seq* (A↔T, C↔G)."""
    return seq.translate(_COMPLEMENT_TABLE)[::-1]


def encode_block(
    seqs: list[str],
    flags: list[int],
    len_bytes_fmt: int,
    npolicy: str,
    do_rc_r1: bool,
    do_rc_r2: bool,
) -> tuple[bytes, bytes, bytes]:
    """Batch-encode sequences into columnar byte streams.

    Returns ``(flags_bytes, lengths_bytes, sequences_bytes)``.

    Handles strand normalisation (reverse-complement) and N-policy
    replacement before 2-bit packing.
    """
    import struct as _struct

    _len_struct = {1: _struct.Struct("<B"), 2: _struct.Struct("<H"), 4: _struct.Struct("<I")}[len_bytes_fmt]
    pack_len = _len_struct.pack

    flags_out = bytes(flags)
    lengths_buf = bytearray()
    seqs_buf = bytearray()

    for i, seq in enumerate(seqs):
        flag = flags[i]
        is_read1 = bool(flag & 1)
        is_read2 = bool(flag & 2)

        if do_rc_r1 and is_read1:
            seq = reverse_complement(seq)
        elif do_rc_r2 and is_read2:
            seq = reverse_complement(seq)

        # N-policy
        if npolicy and "N" in seq.upper():
            if npolicy == "random":
                import random
                seq = "".join(random.choice("ACGT") if c.upper() == "N" else c for c in seq)
            elif npolicy in ("A", "C", "G", "T"):
                seq = seq.replace("N", npolicy).replace("n", npolicy.lower())

        lengths_buf.extend(pack_len(len(seq)))
        seqs_buf.extend(encode_sequence(seq))

    return flags_out, bytes(lengths_buf), bytes(seqs_buf)


def decode_block(
    flags_data: bytes,
    lengths_data: bytes,
    seqs_data: bytes,
    len_bytes: int,
    count: int,
) -> list[tuple[str, bool, bool, bool]]:
    """Batch-decode columnar byte streams into record tuples.

    Returns a list of ``(sequence, is_paired, is_read1, is_read2)``.
    """
    import struct as _struct

    # Parse lengths
    if len_bytes == 1:
        lengths: list[int] | tuple[int, ...] = list(lengths_data)
    elif len_bytes == 2:
        lengths = _struct.unpack(f"<{count}H", lengths_data[: count * 2])
    else:
        lengths = _struct.unpack(f"<{count}I", lengths_data[: count * 4])

    # Bulk-decode all sequence bytes at once, then slice per record
    decode_table = _DECODE_TABLE
    all_decoded = "".join(decode_table[b] for b in seqs_data)

    results: list[tuple[str, bool, bool, bool]] = []
    char_offset = 0
    for i in range(count):
        seq_len = lengths[i]
        enc_len = (seq_len + 3) >> 2
        char_end = char_offset + enc_len * 4
        seq = all_decoded[char_offset : char_offset + seq_len]
        char_offset = char_end

        flag = flags_data[i]
        results.append((
            seq,
            bool(flag & 4),  # is_paired
            bool(flag & 1),  # is_read1
            bool(flag & 2),  # is_read2
        ))

    return results
