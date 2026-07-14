# Fix: strand-normalize single/merged reads as read1

**Status:** proposed fix (2026-07-13)
**Component:** strand normalization in the encode path (`_pycodec.py`, `_accel.cpp`)
**Severity:** correctness — merged/single reads are left in the wrong strand under
`--strand-specific --strand-normalize`.
**Backward compatibility:** changes what *new* encodes store for single reads in
strand-specific mode; already-written ZNA files and all decoding are unaffected.

---

## 1. Summary

Under `--strand-specific --strand-normalize`, ZNA reverse-complements a record **only**
if it is flagged read1 (and read1 is antisense) or read2 (and read2 is antisense). A
**single-end / merged read** carries neither flag (`is_read1 = is_read2 = False`), so it
is **never reverse-complemented** — it stays in its native (read1) orientation while the
paired read1 records in the same file *are* flipped to the sense strand. The result is a
strand-**inconsistent** file.

This matters for the "mixed paired-end + single-end interleaved" feature (README §"Mixed
Paired-End and Single-End Reads"), which is explicitly meant to ingest **fastp**-style
output where merged reads are singles. Those merged reads are read1-oriented (a merged
read is built as `R1 + revcomp(R2)-tail`, entirely in R1's frame), so they should follow
the **read1** normalization rule.

**Fix:** in the deterministic (stranded) branch, treat a single read (neither read1 nor
read2) as **read1**.

---

## 2. Reproduction

```bash
printf '@merged1 merged_20_0\nACGTACGTACGTACGTAAAA\n+\nIIIIIIIIIIIIIIIIIIII\n'\
'@readP/1\nGGGGCCCCGGGGCCCCAAAA\n+\nIIIIIIIIIIIIIIIIIIII\n'\
'@readP/2\nTTTTGGGGTTTTGGGGCCCC\n+\nIIIIIIIIIIIIIIIIIIII\n' > mixed.fastq

zna encode --interleaved --strand-specific --strand-normalize -o ss.zna mixed.fastq
zna decode ss.zna -o -            # stored orientation (no --restore-strand)
```

Observed (v0.3.0):

```
>...:1                              # merged1  (single)
ACGTACGTACGTACGTAAAA                # UNCHANGED  ← bug: should be RC'd to sense
>...:2/1                            # readP/1  (pair R1, antisense)
TTTTGGGGCCCCGGGGCCCC                # RC'd to sense  ✓
>...:2/2                            # readP/2  (pair R2, sense)
TTTTGGGGTTTTGGGGCCCC                # kept sense     ✓
```

The merged single is left antisense; the paired R1 is normalized to sense. After the fix
the merged single becomes `TTTTACGTACGTACGTACGT` (RC), consistent with the pairs.

Unstranded mode (`--strand-normalize` without `--strand-specific`) is **not** affected:
singles already get random RC there.

---

## 3. Root cause

The `do_rc_r1 / do_rc_r2 / do_random_rc` flags are derived in `core.py:204-210`:

```python
self._do_strand_norm_r1 = header.strand_normalized and header.strand_specific and header.read1_antisense
self._do_strand_norm_r2 = header.strand_normalized and header.strand_specific and header.read2_antisense
self._do_random_norm    = header.strand_normalized and not header.strand_specific
```

They are consumed in three places, all with the **same** stranded-branch rule that
ignores single reads:

1. `src/zna/_pycodec.py` — `encode_block`, lines **171-178**
2. `src/zna/_accel.cpp` — `encode_block`, line **377**
3. `src/zna/_accel.cpp` — `encode_block_labeled`, line **641** (identical block)

Decoding is fine: `core.py:576-590` restores strand purely from the stored `IS_RC` bit,
so once encode sets `IS_RC` on RC'd singles, `--restore-strand` round-trips them.

---

## 4. The fix

Define a single read as one carrying neither the read1 nor read2 flag, and normalize it
with the **read1** rule.

### 4.1 `src/zna/_pycodec.py` (`encode_block`)

Replace the stranded branch (lines 171-178):

```python
        else:
            # Stranded (deterministic) normalization
            if do_rc_r1 and is_read1:
                seq = reverse_complement(seq)
                flags[i] |= IS_RC_BIT
            elif do_rc_r2 and is_read2:
                seq = reverse_complement(seq)
                flags[i] |= IS_RC_BIT
```

with:

```python
        else:
            # Stranded (deterministic) normalization.
            # A single/merged read carries neither read1 nor read2. It is
            # read1-oriented (a merged read is R1 + revcomp(R2)-tail), so
            # normalize it with the read1 rule.
            is_single = not is_read1 and not is_read2
            if do_rc_r1 and (is_read1 or is_single):
                seq = reverse_complement(seq)
                flags[i] |= IS_RC_BIT
            elif do_rc_r2 and is_read2:
                seq = reverse_complement(seq)
                flags[i] |= IS_RC_BIT
```

### 4.2 `src/zna/_accel.cpp` (`encode_block`, ~line 375-381)

Replace:

```cpp
        } else {
            // Deterministic stranded normalization
            needs_rc = (do_rc_r1 && is_read1) || (do_rc_r2 && is_read2);
            if (needs_rc) {
                flags_out[i] = static_cast<char>(flag | IS_RC_BIT);
            }
        }
```

with:

```cpp
        } else {
            // Deterministic stranded normalization.
            // A single/merged read (neither R1 nor R2) is R1-oriented:
            // normalize it with the read1 rule.
            bool is_single = !is_read1 && !is_read2;
            needs_rc = (do_rc_r1 && (is_read1 || is_single)) || (do_rc_r2 && is_read2);
            if (needs_rc) {
                flags_out[i] = static_cast<char>(flag | IS_RC_BIT);
            }
        }
```

### 4.3 `src/zna/_accel.cpp` (`encode_block_labeled`, ~line 641)

The `encode_block_labeled` stranded branch has the identical line; apply the same change:

```cpp
            bool is_single = !is_read1 && !is_read2;
            needs_rc = (do_rc_r1 && (is_read1 || is_single)) || (do_rc_r2 && is_read2);
```

> Keep all three in lockstep — see the cross-backend equivalence test in §6.4.

---

## 5. Why "single → read1" is the correct rule

- **Merged reads** (the motivating case) are constructed in read1's orientation, so they
  are antisense exactly when read1 is antisense. RC-ing them iff `read1_antisense` puts
  them on the same (sense) strand as the normalized paired R1 reads. ✓
- **Genuine single-end strand-specific libraries**: the lone read *is* read1 by
  convention, so the read1 rule is already the correct one. ✓
- **read1-sense protocols** (`--read1-sense`, so `read1_antisense = False`): `do_rc_r1`
  is `False`, so singles are left as-is — correct, they are already sense. ✓
- **Unstranded** (`do_random_rc`): unchanged — the random branch already handles singles.

There is no case where a single read should follow the read2 rule (a single is never
read2-oriented), so anchoring on read1 is unambiguous.

---

## 6. Tests to add (`tests/test_core.py`, unittest style)

These lock the fixed behavior. They mirror the existing
`test_single_end_normalization_roundtrip` / `test_mixed_se_pe_normalization_roundtrip`.

### 6.1 Single read is RC'd under strand-specific (the regression test)

```python
def test_single_end_strand_specific_rc(self):
    """A single/merged read must be RC'd as read1 under strand-specific norm."""
    header = ZnaHeader(
        read_group="se_ss",
        strand_specific=True,          # default: read1 antisense, read2 sense
        strand_normalized=True,
    )
    se = "ACGTACGTACGTACGTAAAA"
    se_rc = reverse_complement(se)     # "TTTTACGTACGTACGTACGT"

    with tempfile.TemporaryDirectory() as tmpdir:
        path = f"{tmpdir}/se_ss.zna"
        with open(path, "wb") as fh:
            with ZnaWriter(fh, header) as writer:
                writer.write_record(se, False, False, False)   # single

        # Stored orientation (no restore): must be RC'd to sense.
        with open(path, "rb") as fh:
            stored = list(ZnaReader(fh).records(restore_strand=False))
        self.assertEqual(stored[0][0], se_rc)

        # Restore must recover the original.
        with open(path, "rb") as fh:
            restored = list(ZnaReader(fh).records(restore_strand=True))
        self.assertEqual(restored[0][0], se)
```

### 6.2 Single read NOT RC'd under a read1-sense protocol

```python
def test_single_end_read1_sense_not_rc(self):
    """read1-sense protocol: single read is already sense, must be left as-is."""
    header = ZnaHeader(
        read_group="se_r1sense",
        strand_specific=True,
        read1_antisense=False,         # read1 is sense
        read2_antisense=True,
        strand_normalized=True,
    )
    se = "ACGTACGTACGTACGTAAAA"
    with tempfile.TemporaryDirectory() as tmpdir:
        path = f"{tmpdir}/se_r1sense.zna"
        with open(path, "wb") as fh:
            with ZnaWriter(fh, header) as writer:
                writer.write_record(se, False, False, False)
        with open(path, "rb") as fh:
            stored = list(ZnaReader(fh).records(restore_strand=False))
        self.assertEqual(stored[0][0], se)     # unchanged
```

### 6.3 Mixed stream under strand-specific: singles and R1 end up on the same strand

```python
def test_mixed_se_pe_strand_specific(self):
    """Merged singles and paired R1 must be normalized to the same (sense) strand."""
    header = ZnaHeader(
        read_group="mixed_ss", strand_specific=True, strand_normalized=True,
    )
    merged = "ACGTACGTACGTACGTAAAA"
    r1     = "GGGGCCCCGGGGCCCCAAAA"
    r2     = "TTTTGGGGTTTTGGGGCCCC"
    with tempfile.TemporaryDirectory() as tmpdir:
        path = f"{tmpdir}/mixed_ss.zna"
        with open(path, "wb") as fh:
            with ZnaWriter(fh, header) as writer:
                writer.write_record(merged, False, False, False)   # single (read1-like)
                writer.write_record(r1, True, True, False)         # pair R1 (antisense)
                writer.write_record(r2, True, False, True)         # pair R2 (sense)
        with open(path, "rb") as fh:
            stored = list(ZnaReader(fh).records(restore_strand=False))
        self.assertEqual(stored[0][0], reverse_complement(merged))  # RC'd like R1
        self.assertEqual(stored[1][0], reverse_complement(r1))      # RC'd
        self.assertEqual(stored[2][0], r2)                          # kept
        # round-trip
        with open(path, "rb") as fh:
            restored = list(ZnaReader(fh).records(restore_strand=True))
        self.assertEqual([r[0] for r in restored], [merged, r1, r2])
```

### 6.4 encode_block: IS_RC set on a strand-specific single

```python
def test_encode_block_single_is_rc_strand_specific(self):
    """encode_block must set IS_RC (0x08) on a single read under do_rc_r1."""
    flags_out, _, _ = encode_block(
        ["ACGTACGT"], [0x00], 1, "",       # flag 0x00 = single (no R1/R2/paired bits)
        True,   # do_rc_r1
        False,  # do_rc_r2
        do_random_rc=False,
    )
    self.assertTrue(flags_out[0] & 0x08)
```

### 6.5 (Recommended) cross-backend equivalence

Because the rule lives in both `_pycodec.py` and `_accel.cpp`, add a test that runs the
same records through **both** backends and asserts identical `(flags, sequences)` output,
so the two implementations can never silently drift. Force the Python backend (e.g. by
monkeypatching `_accel_mod`/`_pycodec_mod` selection or calling each `encode_block`
directly) and compare against the C++ backend for the §6.1/§6.3 inputs.

---

## 7. Documentation updates

### 7.1 README — "Mixed Paired-End and Single-End Reads (Interleaved)"

Add a strand-normalization note (there is none today):

> **Strand normalization of merged/single reads:** single-end reads (including merged
> reads with no mate) are treated as **read1** for strand normalization. Under
> `--strand-specific`, a single read is reverse-complemented exactly when read1 is
> antisense, so merged reads end up on the same strand as normalized paired R1 reads.

### 7.2 CHANGELOG

```markdown
## [0.3.1] - 2026-07-XX

### Fixed
- Strand normalization now treats single-end / merged reads as **read1** under
  `--strand-specific`. Previously such reads were never reverse-complemented, leaving
  merged reads (e.g. from fastp's mixed interleaved output) on the opposite strand from
  the normalized paired read1 records. Decoding and already-encoded files are unaffected;
  only new strand-specific encodes of single reads change.
```

---

## 8. Additional suggestions (optional, lower priority)

1. **`zna inspect`: report single vs paired counts.** When ingesting mixed interleaved
   streams it is hard to confirm the paired/single split. Surfacing per-flag record
   counts (paired R1, paired R2, single) in `inspect` would make mixed-stream debugging
   and pipeline validation straightforward.
2. **Keep `_pycodec` and `_accel` provably in lockstep.** Beyond §6.5, consider a small
   parametrized suite that feeds a battery of flag/strand combinations to both backends
   and asserts byte-identical output — the two encoders duplicate non-trivial logic
   (strand rules, N-policy, 2-bit packing) and are a natural drift risk.
3. **Header validation.** Optionally warn if `strand_normalized` is set without
   `strand_specific` and without any random-norm intent, or if both `read1_antisense`
   and `read2_antisense` are false under `--strand-specific` (no read would ever be
   normalized), which usually indicates a misconfigured protocol string.
