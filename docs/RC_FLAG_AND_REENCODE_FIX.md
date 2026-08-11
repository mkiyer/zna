# Fix: expose `IS_RC` to readers, and stop re-normalizing normalized files

Status: **proposed**. Written 2026-08-10 against `baf902a` (v0.3.1-1).
Both defects below were reproduced against the installed 0.3.1 codec.

---

## 1. Summary

Two independent defects, one root cause.

`ZnaHeaderFlags.STRAND_NORMALIZED` is a *record of something that already
happened* — "these sequences were oriented at encode time". But:

- the **writer** reads it as an *instruction*: "orient these sequences", and
- the **reader** discards the per-record `IS_RC` evidence of what was done.

So the fact of normalization survives, and the content of it does not.

| | defect | severity |
|---|---|---|
| **A** | `ZnaReader.records()` decodes `IS_RC` and throws it away. No consumer can recover which end of a mate is a real fragment boundary, and it is not derivable from the sequence. | information loss at the API boundary |
| **B** | Re-encoding an unstranded normalized file applies a **second** random reverse-complement, breaking mate co-orientation and making `IS_RC` meaningless. | **silent data corruption** |

**A has a measured downstream cost.** The `khorana` trainer must decide which
edge of a mate is a real fragment boundary in order to place its `SOS`/`EOS`
supervision. Unable to read `IS_RC`, it guesses by mate number. Measured on
4000 synthetic FR pairs through the real codec: **49.0% of endpoint tokens land
on a read-length truncation point** — ordinary interior sequence where the
molecule demonstrably continues. With the flag exposed and used, **0.0%**. All
of that project's training data is unstranded, so this is the supervision its
shipped 25B-parameter model was trained on.

**B must be fixed with or before A.** In a file that has been re-encoded even
once, `IS_RC` no longer corresponds to the fragment boundary, so exposing it
would hand consumers a flag that is right for fresh files and wrong for
re-encoded ones — with nothing to distinguish them.

---

## 2. Background: what normalization is for

A fastp-style FR pair covers the two ends of one fragment, pointing inward:

```
    fragment F, length L
    |------------------------------------------------|
    |>>>>>>>>>>>|                        |<<<<<<<<<<<<|
     R1 as sequenced                      R2 as sequenced
     = F[0:l1]                            = revcomp(F[L-l2:L])
```

As sequenced the mates are in *opposite* frames. Unstranded normalization
reverse-complements **exactly one** of them so both land in one common frame,
and records which one in `IS_RC`. That is correct and lossless — the encoder
does its job.

The consequence that matters downstream is geometric:

> whichever mate was RC'd ends up at the **right** of the common frame, so its
> **right** edge is the real fragment boundary and its left edge is a
> read-length cutoff. For the other mate it is the mirror image.

```
    common frame after normalization
    |------------------------------------------------|
    |<<<<<<<<<<<|                        |<<<<<<<<<<<<|
     not RC'd                             RC'd
     LEFT edge = real boundary            RIGHT edge = real boundary
     right edge = truncation              left edge  = truncation
```

`IS_RC` is the only thing that distinguishes the two cases, and **it cannot be
recovered from the sequence.** RC'ing the right-hand mate reproduces the
fragment-frame sequence exactly, because that mate was stored
reverse-complemented to begin with. There is no residue in the bases to test.

---

## 3. Reproduction

### 3.1 Defect B — re-encode breaks co-orientation

Self-contained; needs only `zna`.

```python
import io, random
from zna import ZnaWriter, ZnaReader, reverse_complement as rc
from zna.core import ZnaHeader

rng = random.Random(0)
L, READ, N = 300, 100, 500
frags = ["".join(rng.choice("ACGT") for _ in range(L)) for _ in range(N)]

def roundtrip(recs):
    h = ZnaHeader(read_group="p", strand_specific=False, strand_normalized=True)
    buf = io.BytesIO(); w = ZnaWriter(buf, h)
    for r in recs:
        w.write_record(*r)
    w.close()
    return list(ZnaReader(io.BytesIO(buf.getvalue())).records())

# an FR pair AS SEQUENCED: mates point inward, so NOT co-oriented
src = []
for f in frags:
    src.append((f[:READ], True, True, False))           # left mate
    src.append((rc(f[L - READ:]), True, False, True))   # right mate

def co_oriented(out):
    """Both mates in the SAME fragment frame -- what normalization achieves."""
    n = 0
    for i, f in enumerate(frags):
        a, b = out[2 * i][0], out[2 * i + 1][0]
        if (a == f[:READ]) == (b == f[L - READ:]):
            n += 1
    return n

p1 = roundtrip(src)
p2 = roundtrip([tuple(x) for x in p1])      # what `zna` re-encode does
p3 = roundtrip([tuple(x) for x in p2])
print(co_oriented(p1), co_oriented(p2), co_oriented(p3))
```

Observed:

```
    input (as sequenced), co-oriented:     0/500     correct
    after 1 encode,       co-oriented:   500/500     correct
    after re-encode,      co-oriented:     0/500     <-- BUG
    after 2nd re-encode,  co-oriented:   500/500     toggles
```

Every re-encode flips co-orientation. 50% of record *sequences* change on each
pass. The header still says `strand_normalized=True` and every record still
carries an `IS_RC` bit, so nothing downstream can detect that this happened.

### 3.2 Defect A — the flag never reaches the caller

```python
recs = list(ZnaReader(fh).records())
len(recs[0])        # 4  (or 5 with labels) -- IS_RC is not among them
```

The end-to-end cost is reproduced by, in the `khorana` repo:

```
python scripts/debug/check_unstranded_endpoint_geometry.py
```

which reports `49.0%` of special tokens on a non-boundary, and `0.0%` under the
proposed fix. It is self-checking and carries a stranded control.

---

## 4. Root cause

### 4.1 Defect A — `core.py`, `ZnaReader.records()` (~line 474)

Both decode backends already produce the flag: `_pycodec.decode_block` yields
5-tuples ending in `is_rc`, and `_accel.cpp` sets `IS_RC_BIT` (line ~293). The
public generator is the only thing that drops it:

```python
            if has_labels:
                if needs_restore:
                    for seq, is_paired, is_read1, is_read2, is_rc, labels in block_records:
                        if is_rc:
                            seq = rc(seq)
                        yield seq, is_paired, is_read1, is_read2, labels
                else:
                    for seq, is_paired, is_read1, is_read2, _is_rc, labels in block_records:
                        yield seq, is_paired, is_read1, is_read2, labels      # <-- dropped
```

`restore_strand=True` is not a substitute: it *consumes* the flag to undo the
RC, returning original-orientation sequences. A consumer that wants the
normalized frame **and** the boundary geometry — which is what a trainer wants —
has no way to get both.

### 4.2 Defect B — three places conspire

| where | what it does |
|---|---|
| `cli.py:730-732` | a single ZNA input is a supported *re-encode* mode |
| `cli.py:866` | `strand_normalized = flag if flag else input_header.strand_normalized` — carries the bit through |
| `core.py:210` | `self._do_random_norm = header.strand_normalized and not header.strand_specific` |

So the writer normalizes again. The stranded path has the same shape
(`_do_strand_norm_r1/r2`, lines 204-208) and the same latent issue, though it is
deterministic rather than random: RC'ing R1 a second time returns it to its
original orientation, silently un-normalizing the file.

---

## 5. The fix

### 5.1 Writer — never normalize input that is already normalized

This is the load-bearing change and should land first.

The writer cannot currently be told "this record is already oriented, and here
is its flag" — `write_record(seq, is_paired, is_read1, is_read2, labels)` has
nowhere to put it. Two options:

**Option 1 (recommended): a pass-through mode.** Add a writer-level flag, set by
the CLI when the input is a ZNA file whose header is already
`strand_normalized`, meaning *do not orient; the caller supplies the frame and
the `IS_RC` bit verbatim.* Requires `write_record` / `write_records` to accept an
optional trailing `is_rc`.

```python
    ZnaWriter(fh, header, preserve_normalization=True)
    writer.write_record(seq, is_paired, is_read1, is_read2, labels=None, is_rc=False)
```

Pairs naturally with §5.2: `records(with_rc=True)` produces exactly the tuple
`write_record` then consumes, so re-encode becomes a faithful copy.

**Option 2 (smaller, less good): make re-encode refuse.** Have the CLI error on
re-encoding an already-normalized file unless given an explicit override. Cheap
and safe, but it turns a corruption bug into an ergonomics one, and leaves the
library API unable to express a lossless copy.

Whichever is chosen, the invariant to state in the docstring is:

> `strand_normalized` in a header being written describes the **output**. If the
> input is already normalized, the writer must copy orientation rather than
> re-derive it. Orientation is not idempotent: applying it twice returns the
> data to an un-normalized state while leaving the header claiming otherwise.

### 5.2 Reader — expose `IS_RC`, opt-in

**It must be opt-in.** The default tuple width is load-bearing inside this repo
(§6). Suggested:

```python
    def records(self, restore_strand: bool = False, with_rc: bool = False):
        """...
        with_rc
            Append the per-record ``IS_RC`` flag to each tuple. This is the only
            way to recover which edge of a mate is a real fragment boundary; it
            cannot be derived from the sequence. Mutually exclusive in spirit
            with ``restore_strand``, which consumes the flag to undo the RC.
        """
```

Yielding `(seq, is_paired, is_read1, is_read2, is_rc)` and, for labeled files,
`(seq, is_paired, is_read1, is_read2, labels, is_rc)`.

Note the labeled ordering question: appending after `labels` keeps the existing
5-tuple prefix stable for `rec[:5]`-style slicing, which is what
`khorana/data.py` does. Inserting `is_rc` before `labels` would be tidier but
breaks that. **Recommend appending**, and saying so explicitly in the docstring
because it is the kind of thing that gets "cleaned up" later.

Consider raising if both `restore_strand=True` and `with_rc=True` are passed —
after restoration the flag describes work already undone, and every caller that
wants both is confused about one of them.

### 5.3 Longer term — a record type

Four positional booleans plus two optional trailing fields is at the limit of
what a tuple should carry, and this document exists partly because adding a
field is scary. A `NamedTuple` (`ZnaRecord`) would be tuple-compatible for
existing positional unpacking of the first four fields, make `record.is_rc`
self-documenting, and make the next field free. Worth costing separately; not
required for this fix.

---

## 6. Compatibility analysis

Every in-repo consumer of `records()`, checked:

| site | unpacking | safe under an opt-in change? |
|---|---|---|
| `_shuffle.py:80,141,217` | `rec[1]`, `rec[2]` — indexed | yes |
| `cli.py:617` (`_stream_zna_reencode`) | yields `record` whole into the writer path | yes, but see below |
| `cli.py:1052` | `seq, is_paired, is_r1, is_r2, labels = rec` — **fixed width 5** | breaks if the default width changes |
| `cli.py:1078` | same shape | same |
| `core.py:677` (`read_zna`) | passes the generator through | yes |
| `tests/test_cli.py:1526` | `seq, is_paired, is_r1, is_r2 = out_records[i]` — **fixed width 4** | breaks if the default width changes |

External: `khorana/data.py:531` does `record[:4]`, which is width-tolerant.

Conclusion: **do not change the default tuple width.** `cli.py:1052/1078` and
`test_cli.py:1526` are the tripwires.

Also note `cli.py:617` yields records straight into `stream_inputs`, whose
consumers at `cli.py:949/957` unpack fixed-width 4- and 5-tuples. If §5.1
Option 1 is taken, that path needs to carry one more field end to end.

### Format compatibility

None of this changes the on-disk format. `IS_RC` is already written, already
read, already tested. Existing files stay readable and their bytes stay valid.
§5.1 changes the *output* of a re-encode, which is a bug fix, not a format
change — but it is a behavioural change to a shipped code path and belongs in
the changelog as such.

**Files already re-encoded an odd number of times are corrupt and cannot be
repaired**, because the information needed to undo it was never stored. Worth a
sentence in the release notes so anyone with such a file can re-derive it from
FASTQ rather than trusting it.

---

## 7. Implementation plan

Ordered so each step is independently testable.

1. **Regression test for B first**, from §3.1, asserting co-orientation survives
   a re-encode. It should fail on `baf902a`.
2. **Writer pass-through** (§5.1). Decide Option 1 vs 2 before writing code —
   Option 1 is more work and is the one that makes a lossless `zna` → `zna` copy
   expressible at all.
3. **`records(with_rc=...)`** (§5.2), python backend path first.
4. **Verify the accel backend needs nothing.** `_accel.cpp` already emits
   `is_rc` in its record tuples; the change should be confined to `core.py`. If
   that turns out to be wrong, the `.so` needs a rebuild, not just an edit.
5. **Cross-backend equivalence test** — same file, `get_backend` forced each
   way, identical `is_rc` per record. The precedent doc's §6.5 does this.
6. **CLI**: thread the flag through `_stream_zna_reencode` → `stream_inputs` →
   writer; make `zna inspect --counts` report RC'd counts if it does not already
   (the 0.3.1 changelog suggests it does).
7. **Docs**: README section on unstranded normalization stating the geometric
   invariant from §2 — that is the part no reader currently knows; CHANGELOG;
   version bump.

Semver: §5.2 is additive (minor). §5.1 changes existing behaviour, and since
the old behaviour was corrupting data, "fixed" is the honest label.

---

## 8. Tests to add

Matching the house style of `docs/STRAND_NORM_MERGED_READ_FIX.md` §6
(`unittest`, in `tests/test_core.py`).

1. **`test_reencode_preserves_co_orientation`** — §3.1, the regression test for B.
2. **`test_reencode_is_byte_identical`** — stronger and worth trying: encode,
   re-encode, assert the record stream is unchanged. Compression framing may
   differ, so compare decoded records rather than file bytes.
3. **`test_records_with_rc_exposes_the_flag`** — exactly one mate of each
   unstranded pair has `is_rc=True`.
4. **`test_rc_flag_identifies_the_boundary_edge`** — the claim consumers will
   actually rely on: for a known fragment, the RC'd mate's right edge and the
   other mate's left edge coincide with the fragment ends. This is the test that
   makes the geometry a contract rather than folklore.
5. **`test_records_default_width_unchanged`** — pins the tuple width at 4 (and 5
   labeled) so the compatibility promise in §6 cannot regress silently.
6. **`test_with_rc_and_restore_strand_conflict`** — if §5.2's guard is adopted.
7. **`test_stranded_reencode_is_stable`** — the deterministic mirror of B; RC'ing
   R1 twice must not silently un-normalize a strand-specific file.
8. **Cross-backend equivalence** for 3 and 4.

---

## 9. Downstream: what `khorana` does with this

Recorded so the API shape is chosen with its consumer in mind, and so this
document is enough context on its own.

`khorana`'s `BaseRNASeqDataset._assign_endpoints` decides whether a read carries
the biological 5′ origin (`SOS`) or 3′ terminus (`EOS`). For **stranded** data it
derives this correctly from the header. For **unstranded** data it cannot see
`IS_RC`, so it hard-codes `R1 → EOS, R2 → SOS` — and since the RC coin is
independent of the mate number, half of those labels land on the wrong edge.

With `with_rc` available the fix there is two lines:

```python
    has_eos = is_rc            # the RC'd mate: its RIGHT edge is the boundary
    has_sos = not is_rc        # the other mate: its LEFT edge is
```

verified to take non-boundary special tokens from 49.0% to 0.0%.

That project tracks this as **M0** in `docs/OPEN_ISSUES.md`, currently marked
*blocked on a `zna` reader-API change*. It only benefits their model after a
retrain, so there is no urgency beyond getting the API right — but the API is on
their critical path, and a released `zna` with `with_rc` unblocks it.

---

## 10. Open questions for the implementer

1. **Is the re-encode path used in practice on unstranded files?** If yes,
   Defect B has already corrupted data and the release notes need to say so
   plainly. If it is effectively dead, Option 2 in §5.1 becomes reasonable.
2. **Should `strand_normalized` be split** into "orientation was applied" (a
   provenance bit, immutable) and "apply orientation" (an encode-time request)?
   That is the clean fix to the root cause, and it is a header-semantics change
   — bigger than this document's scope, but everything here is a symptom of the
   two meanings sharing one bit.
3. **Does anything outside these two repos consume `records()`?** The
   compatibility argument in §6 assumes not.
