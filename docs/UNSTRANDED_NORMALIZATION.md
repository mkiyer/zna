# Proposal: Unified Strand Normalization

**Status:** Draft / Under Review  
**Date:** 2026-03-24 (rev 2)

---

## 1. Motivation

ZNA currently normalizes strand-specific (stranded) paired-end data by
reverse-complementing antisense reads at encode time, placing all stored
sequences on the "sense" strand. This is controlled by the `--strand-specific`
flag and the header flags `READ1_ANTISENSE` / `READ2_ANTISENSE`. For
unstranded data, no normalization is performed.

### Why normalize unstranded data?

For unstranded paired-end data, read1 and read2 are on opposite strands
(Illumina convention), but we do not know which is sense or antisense. If we
want all stored reads on a consistent strand (useful for kmer analysis,
deduplication, compression), we must reverse-complement one read per pair.
Since we cannot know which read is "correct," we choose at random. This
guarantees both mates are on the **same** strand, without guaranteeing it is
the sense strand.

For unstranded single-end data, the strand is similarly unknown. Randomly
reverse-complementing each read (with a recorded per-record flag) keeps the
code path uniform and the operation fully reversible.

### Design insight: two orthogonal axes

The current design conflates two independent concepts:

1. **Library metadata** — Is the library strand-specific? If so, which
   reads are sense/antisense? These are properties of the sequencing
   protocol and should always be stored in the header.

2. **Strand normalization** — Should sequences be reverse-complemented
   at encode time to place them on a consistent strand? This is an
   encoding action that can apply to both stranded and unstranded data.

**Proposal:** Separate these into two orthogonal settings:

| | Normalize OFF | Normalize ON |
|---|---|---|
| **Stranded** | Reads stored as-is. Header records sense/antisense metadata. | Antisense reads RC'd to sense strand (current behavior). |
| **Unstranded** | Reads stored as-is (current behavior). | One paired read RC'd at random; SE reads RC'd at random. |

A new per-record flag (`IS_RC`) records which reads were reverse-complemented,
making the operation perfectly reversible regardless of mode.

---

## 2. Current Architecture Summary

### 2.1 Header Flags (1 byte, `core.py`)

```
Bit 0 (0x01): STRAND_SPECIFIC
Bit 1 (0x02): READ1_ANTISENSE
Bit 2 (0x04): READ2_ANTISENSE
Bits 3-7:     Unused (free for new flags)
```

### 2.2 Per-Record Flags (1 byte per record, stored in block payload)

```
Bit 0 (0x01): IS_READ1
Bit 1 (0x02): IS_READ2
Bit 2 (0x04): IS_PAIRED
Bits 3-7:     Unused (free for new flags)
```

### 2.3 Current Strand Normalization Path (stranded data only)

**Encode:** `ZnaWriter._flush_block()` → `encode_block(do_rc_r1, do_rc_r2)` →
if a read matches the antisense flag it is reverse-complemented before 2-bit
encoding. Both `_pycodec.py` and `_accel.cpp` implement this.

**Decode:** `ZnaReader.records(restore_strand=True)` → after decoding,
the reader checks `header.strand_specific` and `header.read1_antisense` /
`read2_antisense` to determine which reads to reverse-complement back to
their original orientation. There is no per-record flag; the decision is
based entirely on header metadata + read identity (R1 vs R2).

### 2.4 Mixed Single-End / Paired-End Support

Yes — the format supports it. Record flags distinguish unpaired (`flag=0x00`),
paired-R1 (`flag=0x05`), and paired-R2 (`flag=0x06`) on a per-record basis.
The CLI's interleaved-FASTQ parser (`_stream_interleaved_fastq`) can emit
both unpaired and paired records in the same stream when consecutive read
names do not match.

---

## 3. Design

### 3.1 Core Principle: Orthogonal Axes

**Library metadata** (stored in header, always available):
- `strand_specific: bool` — Is this a stranded library?
- `read1_antisense: bool` — R1 represents the antisense strand (meaningful
  only when `strand_specific=True`).
- `read2_antisense: bool` — R2 represents the antisense strand (meaningful
  only when `strand_specific=True`).

**Encoding action** (stored in header + per-record flags):
- `strand_normalized: bool` — Were reads strand-normalized during encoding?

**Per-record state**:
- `IS_RC` (bit 3 of per-record flag byte) — Was this specific record
  reverse-complemented during encoding?

### 3.2 New Header Flag

Add a new header-level flag:

```
Bit 3 (0x08): STRAND_NORMALIZED
```

This flag is independent of `STRAND_SPECIFIC`. They represent different
things: one is a library property, the other is an encoding action.

Add a corresponding field to `ZnaHeader`:

```python
strand_normalized: bool = False
```

### 3.3 New Per-Record Flag: `IS_RC`

Add a per-record flag to track whether a record was reverse-complemented:

```
Bit 3 (0x08): IS_RC
```

This bit is set on every record that was reverse-complemented during
encoding — regardless of whether the normalization was stranded (deterministic)
or unstranded (random). This makes decode trivially simple: check `IS_RC`,
reverse-complement if set.

### 3.4 Encode-Time Behavior

When `strand_normalized=True`, `encode_block` applies reverse-complement and
sets the `IS_RC` flag. The logic differs by library type:

#### 3.4.1 Stranded + Normalized

Deterministic. The header's `read1_antisense` / `read2_antisense` tells us
exactly which reads to RC:

```
For each record:
    if (read1_antisense AND record is R1) → RC, set IS_RC
    elif (read2_antisense AND record is R2) → RC, set IS_RC
    else → store as-is, IS_RC = 0
```

This is equivalent to the current behavior but now additionally records the
action in `IS_RC`. The header flags (`read1_antisense`, `read2_antisense`)
become **purely metadata** — they describe the library, while `IS_RC`
records what was actually done.

#### 3.4.2 Unstranded + Normalized

Random. We do not have antisense metadata to guide us, so:

**Paired reads:**
```
For each consecutive R1 + R2 pair:
    coin = random_bit()
    if coin == 1: RC read1, set IS_RC on R1
    else:         RC read2, set IS_RC on R2
```

**Single-end reads:**
```
For each unpaired record:
    coin = random_bit()
    if coin == 1: RC the read, set IS_RC
    else:         store as-is, IS_RC = 0
```

Rationale for normalizing single-end reads: for unstranded data, a single-end
read's strand is unknown. Randomly RC'ing (with a recorded flag) keeps the
code path uniform and the operation fully reversible. It does not achieve
strand "consistency" for SE reads (there is no mate to be consistent with),
but it is harmless, and it avoids special-casing single-end reads in an
otherwise uniform normalization pipeline.

#### 3.4.3 Not Normalized (either library type)

No reverse-complement is applied. All `IS_RC` bits are `0`. Sequences are
stored verbatim. Library metadata (strand_specific, antisense flags) is
still stored in the header for informational purposes.

### 3.5 Decode-Time Behavior

The decode path becomes simple and **uniform** across all modes:

```python
if restore_strand and header.strand_normalized:
    for record in block_records:
        if record.is_rc:
            seq = reverse_complement(seq)
        yield seq, ...
else:
    yield from block_records
```

No need to check `strand_specific`, `read1_antisense`, `read2_antisense`,
or read identity at decode time — `IS_RC` carries all the information needed
to reverse the operation. This is a significant simplification over the
current reader, which must reconstruct the RC decision from header flags.

### 3.6 Handling Mixed Single-End + Paired-End

When `strand_normalized=True`, the encode loop handles three record types:

| Record Type | Stranded + Normalized | Unstranded + Normalized |
|---|---|---|
| Paired R1 + R2 | RC the antisense read (deterministic) | RC one at random |
| Unpaired (SE) | RC if antisense (deterministic) | RC at random |

For the stranded case with single-end reads: if `read1_antisense` is set
and the read is conceptually "R1" (or the library default), it gets RC'd.
In practice, single-end reads in a stranded library flow through the same
per-record check: if the read's identity matches the antisense flag, RC it.
Single-end reads have `IS_READ1=0, IS_READ2=0`, so they won't match the
R1/R2 antisense checks and will be stored as-is. If a user wants SE reads
normalized for stranded data, they would need to tag them as R1 or R2 at
input time (which the CLI already does not do for SE input). This is correct
behavior — we only normalize reads whose strand relationship we know.

### 3.7 Block Boundary Safety

Paired records (R1 then R2) are always emitted consecutively by all CLI
input streaming functions. However, the `write_record` method checks the
flush threshold after each record, which could split a pair across blocks.

**Mitigation:** When `strand_normalized=True`, defer the flush check when
the last record written is a paired R1. Flush only after its R2 arrives.
Add a safety assertion in `_flush_block` that no dangling paired R1 exists.

The batch `write_records` method has its own inline flush logic and needs
the same deferral.

### 3.8 CLI Interface

**Encode:**
```
--strand-normalize       Enable strand normalization
--strand-specific        Library is strand-specific (sets library metadata)
--read1-sense            R1 is on the sense strand (override dUTP default)
--read2-antisense        R2 is on the antisense strand
```

`--strand-normalize` can be used alone (unstranded normalization) or together
with `--strand-specific` (stranded normalization). Without `--strand-normalize`,
reads are stored as-is regardless of `--strand-specific`.

**Decode:**
```
--restore-strand         Reverse-complement IS_RC reads to original orientation
```

The existing `--restore-strand` flag works unchanged — it simply checks
`strand_normalized` in the header and reverses any record with `IS_RC` set.

### 3.9 Backward Compatibility

#### Reading Old Files (no `STRAND_NORMALIZED` flag)

Old files have `STRAND_NORMALIZED = 0` and no `IS_RC` bits set. The new
reader sees `strand_normalized=False` and skips all restoration → identical
to current behavior.

However, old files that used stranded normalization (the current encode
path) also have `STRAND_NORMALIZED = 0` because the flag didn't exist.
For these files, the reader must fall back to the **legacy restore path**:
check `strand_specific`, `read1_antisense`, `read2_antisense`, and restore
based on read identity (R1/R2). This can be handled by:

```python
if header.strand_normalized:
    # New path: use IS_RC per-record flag
    ...
elif restore_strand and header.strand_specific:
    # Legacy path: use header antisense flags + read identity
    ...
else:
    yield from block_records
```

#### Old Readers on New Files

Old readers ignore bits 3-7 in both header and per-record flags. A new
file with `STRAND_NORMALIZED` set will be read by old code as if the flag
doesn't exist. The sequences are valid 2-bit encoded data, so reads will
decode correctly — they just won't be restored to original orientation.
The old `restore_strand` logic would also not trigger for unstranded files
(since `strand_specific=False`). For stranded + normalized new files, the
old reader's restore logic would still work correctly because the
deterministic R1/R2 antisense flags are preserved. **Safe.**

---

## 4. File Format Changes

### 4.1 Header Flags Byte

| Bit | Name | Existing? |
|-----|------|-----------|
| 0 | STRAND_SPECIFIC | Yes |
| 1 | READ1_ANTISENSE | Yes |
| 2 | READ2_ANTISENSE | Yes |
| **3** | **STRAND_NORMALIZED** | **New** |
| 4-7 | Reserved | — |

### 4.2 Per-Record Flags Byte

| Bit | Name | Existing? |
|-----|------|-----------|
| 0 | IS_READ1 | Yes |
| 1 | IS_READ2 | Yes |
| 2 | IS_PAIRED | Yes |
| **3** | **IS_RC** | **New** |
| 4-7 | Reserved | — |

### 4.3 Version Compatibility

The file format version remains `1`. Older readers only test bits 0-2 in
per-record flags via bitmask (`flag & 1`, `flag & 2`, `flag & 4`), so the
`IS_RC` bit (bit 3) is invisible to them. **Verify this holds in both
`_pycodec.py` and `_accel.cpp` decode paths before relying on it.** If any
path reads the full byte in a way that would break, bump to version `2`.

---

## 5. Changes by File

### `src/zna/core.py`

- Add `STRAND_NORMALIZED = 8` to `ZnaHeaderFlags`.
- Add `IS_RC = 8` to `ZnaRecordFlags`.
- Add `strand_normalized: bool = False` to `ZnaHeader`.
- Serialize/deserialize the new header flag in `_write_file_header` /
  `_read_file_header`.
- In `ZnaWriter.__init__`, compute:
  - `self._do_stranded_norm_r1 = strand_normalized and strand_specific and read1_antisense`
  - `self._do_stranded_norm_r2 = strand_normalized and strand_specific and read2_antisense`
  - `self._do_random_norm = strand_normalized and not strand_specific`
- Pass these to `encode_block`.
- In `write_record` / `write_records`: defer flush when a paired-R1 was
  just buffered and `strand_normalized` is on (prevent pair splitting).
- In `ZnaReader.records`:
  - If `restore_strand` and `strand_normalized`: reverse-complement any
    record with `IS_RC` set (new path).
  - Elif `restore_strand` and `strand_specific` and not `strand_normalized`:
    use legacy header-based restore (backward compat for old files).

### `src/zna/_pycodec.py`

- Update `encode_block` signature: add `do_random_rc: bool` parameter.
- When `do_random_rc=True`:
  - For each paired R1+R2 pair: random bit → RC one, set `IS_RC` on its flag.
  - For each unpaired record: random bit → RC or not, set `IS_RC` accordingly.
- When `do_rc_r1` or `do_rc_r2` is `True` (stranded normalization):
  existing logic, but now additionally set `IS_RC` on the record's flag
  byte (mutate `flags[i] |= 0x08`).
- `decode_block`: expose the `IS_RC` bit. Either:
  - Return it as part of the tuple: `(seq, is_paired, is_read1, is_read2, is_rc)`.
  - Or return the raw flag byte and let `records()` extract it.
  Preferred: return the raw flag byte alongside the tuple. The reader
  extracts `is_rc = bool(flag & 8)`.

### `src/zna/_accel.cpp`

- Mirror the `_pycodec.py` changes in C++.
- Add `do_random_rc` parameter to `encode_block`.
- Use a fast PRNG (`std::mt19937` seeded once) for random bit generation.
- Set `IS_RC` (bit 3) on the flag byte for RC'd records.

### `src/zna/codec.py`

- Update the codec dispatch / wrapper if it wraps `encode_block`.

### `src/zna/cli.py`

- Replace `--strand-specific` encode behavior:
  - `--strand-specific` remains: sets library metadata only.
  - Add `--strand-normalize`: enables normalization (can combine with
    `--strand-specific` or stand alone).
- When `--strand-normalize` without `--strand-specific`:
  sets `strand_normalized=True`, `strand_specific=False`.
- When `--strand-normalize` with `--strand-specific`:
  sets `strand_normalized=True`, `strand_specific=True`,
  plus `read1_antisense` / `read2_antisense` per protocol flags.
- When `--strand-specific` without `--strand-normalize`:
  sets `strand_specific=True`, `strand_normalized=False` (metadata only,
  no RC applied).
- `--restore-strand` at decode time: works for both modes via `IS_RC`.

### `tests/test_core.py`

- Test: stranded + normalized → deterministic RC + `IS_RC` set on correct reads.
- Test: unstranded + normalized → one read per pair RC'd + `IS_RC` set.
- Test: unstranded + normalized → SE reads randomly RC'd + `IS_RC` set.
- Test: `restore_strand=True` recovers original sequences (both modes).
- Test: mixed SE + PE with unstranded normalization: pairs normalized,
  singletons randomly normalized, all reversible.
- Test: `strand_specific=True` without normalization → reads stored as-is,
  library metadata preserved.
- Test: backward compat — old-style stranded file (no `STRAND_NORMALIZED`)
  still restores correctly via legacy path.

### `tests/test_cli.py`

- Test: `--strand-normalize` alone (unstranded) round-trips correctly.
- Test: `--strand-normalize --strand-specific` round-trips correctly.
- Test: `--strand-specific` without `--strand-normalize` stores metadata
  but does not RC.
- Test: `--restore-strand` at decode reverses both modes.

---

## 6. Edge Cases & Risks

### 6.1 Pair Splitting Across Blocks

**Risk:** A paired-R1 record triggers a block flush, and its R2 mate lands
in the next block. The `encode_block` lookahead would not find the mate.

**Mitigation:** Defer the flush threshold check when the last buffered
record is a paired R1 and `strand_normalized=True`. Flush only after the
R2 arrives. Add a safety assertion in `_flush_block`.

### 6.2 Single-End Reads with Unstranded Normalization

**Behavior:** Randomly RC each SE read. Set `IS_RC` accordingly. This keeps
the pipeline uniform and is fully reversible, though it does not achieve
inter-read strand consistency (no mate to be consistent with).

### 6.3 Single-End Reads with Stranded Normalization

**Behavior:** SE reads have `IS_READ1=0, IS_READ2=0`, so neither the R1
nor R2 antisense check matches. They pass through un-RC'd with `IS_RC=0`.
This is correct — we only normalize reads whose strand we know.

### 6.4 Reproducibility / Determinism

Random normalization means encoding the same input twice produces different
binary output. The decoded content (with `restore_strand=True`) is always
identical. If determinism is desired, an optional `--seed` argument could
be added in a follow-up.

### 6.5 Backward Compatibility of Per-Record `IS_RC` Bit

**Risk:** Older readers may interpret the `IS_RC` bit unexpectedly.

**Analysis:** The current decode path in `_pycodec.py` reads flags as:
```python
flag = flags_data[i]
bool(flag & 4)  # is_paired — masks bit 2 only
bool(flag & 1)  # is_read1  — masks bit 0 only
bool(flag & 2)  # is_read2  — masks bit 1 only
```
Bit 3 is never tested. Older code ignores `IS_RC`. **Safe.**

The C++ decode path should be verified similarly.

### 6.6 `write_records` Batch API

The batch `write_records` method has its own inline flush logic and must
also defer flushes mid-pair, mirroring the `write_record` change.

### 6.7 Legacy Stranded Files

Old files that were stranded-normalized **before this change** have
`STRAND_NORMALIZED = 0` (the bit didn't exist). The reader must detect this
and use the legacy restore path (header flags + read identity). See
Section 3.9 for the fallback logic.

---

## 7. Implementation Order

1. **`core.py` flags & header** — Add `STRAND_NORMALIZED`, `IS_RC`, header
   field, serialization. Low risk, testable in isolation.
2. **`_pycodec.py` encode** — Add `do_random_rc` parameter. Implement
   paired-lookahead random RC + SE random RC. Set `IS_RC` on affected
   records. Also set `IS_RC` in the existing stranded RC path.
3. **`core.py` writer** — Compute new normalization booleans, pass them
   to `encode_block`. Add flush-deferral for pairs.
4. **`core.py` reader** — Add new `IS_RC`-based restore path + legacy
   fallback. Test roundtrips.
5. **`_accel.cpp`** — Port Python codec changes to C++.
6. **`cli.py`** — Add `--strand-normalize`, update `--strand-specific`
   semantics, wire through to header.
7. **Tests** — Full suite for both modes, mixed SE/PE, backward compat.
8. **Docs** — Update README and binary format spec.

---

## 8. Open Questions

1. **Should we add `--seed` for reproducible random normalization?** Useful
   for testing and debugging. Could be a follow-up.
2. **Should `zna info` display the `STRAND_NORMALIZED` and `IS_RC` flags?**
   Probably yes.
3. **Should `--strand-normalize` without `--strand-specific` warn if input
   is single-end only?** The operation is valid (random RC of each read),
   but the user may not realize there is no pair-based consistency. A
   one-line stderr note seems appropriate.
4. **Migration path for existing stranded workflows:** Users currently
   passing `--strand-specific` get normalization automatically. After this
   change, `--strand-specific` alone stores metadata only. Users must add
   `--strand-normalize` to get the RC behavior. Should we:
   (a) Require explicit `--strand-normalize` (breaking change, cleaner),
   (b) Have `--strand-specific` imply `--strand-normalize` for backward
       compat, with a deprecation warning, or
   (c) Keep current behavior: `--strand-specific` implies normalization,
       add `--no-strand-normalize` to opt out?
   **Recommend (a)** for clarity, with a clear migration note in the
   changelog. The old behavior was implicit; making it explicit is better.
