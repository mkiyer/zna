# Cross-repo plan: correct fragment-boundary supervision, end to end

Status: **proposed**. Written 2026-08-11, revised after auditing the full chain.
Target: **zna 0.3.3** (all zna fixes bundled), with coordinated changes in
`hulkrna` and `khorana`. All numbers were measured, not estimated.

**Why one document.** Every defect found here was a *cross-repo assumption that
nobody checked*: `khorana` assumed pairing survived encoding, `zna` assumed labeled
input was single-end, `hulkrna` assumed `zna encode` honoured `--interleaved`. Each
repo's code was self-consistent and the suite was green in all three. Splitting
this into three documents would rebuild exactly the seam the bugs lived in. So the
contract, the ordering, and the verification gates live here, in `zna`, because
`zna` defines the format contract that the other two depend on. §8–§10 are
per-repo checklists an implementer can work from directly.

The chain under audit:

```
raw FASTQ -> fastp (trim, UMI, adapter) -> dedup -> STAR/rigel -> BAM
          -> samtools fastq -> hulkrna-merge -> zna encode -> ZNA -> khorana
```

> **Note (2026-08-12).** `hulkrna-merge` has since moved into this repo and is now the
> `zna merge` subcommand (`src/zna/merge/`); read every mention of it below as that.
> Its own documents are [READ_MERGE_REDESIGN.md](READ_MERGE_REDESIGN.md) (the design),
> [MERGE_TOOL_AUDIT.md](MERGE_TOOL_AUDIT.md) (what was measured and rejected) and
> [READ_MERGE_ROADMAP.md](READ_MERGE_ROADMAP.md) (status). The §4 properties this plan
> verified by inspection are now stated as a contract, C1–C7, in
> [READ_MERGE_PORT_TO_ZNA.md](READ_MERGE_PORT_TO_ZNA.md) §4, and pinned by
> `tests/test_merge.py` and `tests/test_merge_encode.py`. Line references into
> `pairs.py` throughout this document predate the geometry rewrite and have moved.

---

## 1. The contract

`khorana` marks the two ends of the sequenced fragment with `SOS` and `EOS`:

> A special token must never land at an interior position. The model infers true
> transcript termini from where fragment-end tokens *pile up*, so an interior token
> is not a missing label — it is noise injected into the signal being learned.

Polarity is secondary: strand-specific data must place `SOS` at the transcript-5'
end and `EOS` at the 3' end; for unstranded data the two are equivalent and a
random assignment is correct.

So `zna` must answer one question per record, exactly: **which edge, or edges, of
this stored sequence are true fragment boundaries?**

## 2. The invariant, and the one thing missing from it

Every Illumina read begins at a fragment boundary and reads inward, so base 0 of a
read is a true fragment end. `zna` stores the read either as submitted (`IS_RC`
clear, base 0 leftmost) or reverse-complemented (`IS_RC` set, base 0 rightmost):

| `IS_RC` | true fragment boundary is at |
|---|---|
| clear | the **left** edge |
| set | the **right** edge |

This holds per record with no knowledge of pairing, mate number, or strandedness,
and survived every degenerate input tested (zero-length and 1-base sequences,
palindromes, orphan mates, pairs split across blocks, both codec backends).

The premise is that nothing removed *template* bases from the read's 5' end.
Removing *synthetic* prefix — UMI, spacer, adapter — is not merely safe, it is
required to recover the true boundary. `hulkrna` handles UMI and adapter removal
upstream of the merge, so base 0 is correct by the time reads reach `zna`. What
would break it is blind 5' clipping of template bases (`fastp --trim_front1`,
5'-side quality trimming); the pipeline does none.

**What the format cannot express is the two-ended record.** A read spanning its
entire fragment has **both** edges as true boundaries. `IS_RC` names one and there
is no way to say there is another. In this pipeline that is the common case:
**~80% of emitted records are full-fragment merged reads.**

## 3. What is broken

Ranked by severity. Only §3.1 and §3.2 can place an interior token.

### 3.1 CATASTROPHIC, in production: `--label-defs` destroys all pairing

`zna encode --interleaved --label-defs …` **silently ignores `--interleaved`** and
flags every record unpaired. `stream_inputs_labeled` (`cli.py:760`) has its own
single-file FASTQ reader and yields `(seq, False, False, False, labels)`
unconditionally; the interleaved, two-file, and ZNA-input strategies in
`stream_inputs` are never reached. Measured on a mixed interleaved stream of two
pairs plus one merged read:

| | `paired` | `read1` | `read2` |
|---|---|---|---|
| without `--label` | T, T, F, T, T | correct | correct |
| **with `--label`** | **False ×5** | **False ×5** | **False ×5** |

`workflow/rules/zna.smk:108` passes `--label-defs`, so **every ZNA file this
pipeline has produced has `is_paired=False` on every record.** Consequences:

1. `khorana` returns `(True, True)` for every unpaired record (`data.py:249-251`),
   so **every unmerged mate gets one interior token** — its far edge is a
   read-length truncation. Unmerged mates are ~20% of records, so ~10% of endpoint
   tokens are interior, on top of §3.2.
2. Pair co-orientation is destroyed: with no `IS_PAIRED` bit the encoder draws an
   independent coin per record instead of one per pair.
3. `is_paired` is **not recoverable**. Affected files must be re-encoded.

### 3.2 Unstranded endpoint assignment guesses by mate number

`khorana`'s unstranded branch keys on mate number, which is independent of `zna`'s
coin, so ~half of pairs get *both* tokens interior: **792/1600 (49.5%)** measured.
Strand-specific paired-end is **correct today** (0/1600) — its header-driven rule
reproduces `IS_RC` exactly. Fixed by `with_rc` (already implemented) plus §10.

### 3.3 `--npolicy drop` is not pair-atomic

`--npolicy drop` (the CLI default, and in effect for this pipeline) filters records
*after* pairing, so it can drop one mate and leave a lone `IS_PAIRED` record —
exactly the breakage `hulkrna-merge` works to avoid (`pairs.py:147-152`). A lone
`IS_PAIRED` record is also the only route to the encoder's mis-pairing window
(§6.2). Must be fixed in `zna`, and the flag passed explicitly in `hulkrna`.

### 3.4 Spurious overlap detection in `hulkrna-merge`

`trim_overlap_require: 3` (`config.yaml:141`) lets the scan accept a **4-base**
overlap, and the per-offset budget `floor(0.20 × 4)` is **zero** — so a chance
4-mer match is accepted. Each candidate overlap length is tried exactly once
across ~146 offsets, and 4 exact bases is 1 in 256.

Measured false-positive rate on genuinely unrelated 2×150 read pairs, 20 000 pairs
per cell, using the real kernel:

| `--trim-overlap` | false-positive rate |
|---|---|
| **3 (current)** | **5.17%** |
| 6 | 0.51% |
| 9 | 0.18% |
| **12** | **0.00%** |
| 20 | 0.00% |

The false-positive cost that matters is **not** the wrongly-trimmed bases — those
are a few bases off a long read. It is that a spurious *forward* hit preempts the
reverse (read-through) scan, so **~1.7% of read-through pairs are emitted unmerged
with adapter retained**. Adapter sequence in a pretraining corpus is worse than a
few duplicated bases. Boundary geometry survives either way (base 0 is never
touched), so this is corpus quality, not token placement.

**The objective is to minimise redundancy**, since duplicated sequence biases the
model — so a *missed* real overlap costs more than a spurious trim. That rules out
raising the threshold much. Measured against both axes, with a realistic 0.5%
per-base sequencing error rate:

| rule | FP rate | redundant bases/pair |
|---|---|---|
| current (`require=3`, floor budget) | 5.17% | 0.02 |
| **`require=3` + exact match for olen < 10** | **1.27%** | **0.10** |
| `require=3` + exact match for olen < 20 | 1.20% | 0.80 |
| `require=5` + exact match for olen < 10 | 0.18% | 0.49 |
| `require=12`, floor budget | 0.02% | ~3.3 |

**Fix: keep `trim_require = 3`; require a zero-mismatch match for overlaps shorter
than 10 bases.** Three lines in `_scan`: after computing `dl`, add
`if olen < exact_below: dl = 0`, with `exact_below = 10` (expose as
`--overlap-exact-below` for tuning). That cuts false positives 4× for +0.08
redundant bases per pair, and 4- and 5-base overlaps are still detected at 99%, so
trimming stays exactly as aggressive as it is today.

Two rejected alternatives, for the record:

- **"Round instead of floor" makes it 3.5× worse** — 5.17% → 17.92% — because it
  grants a 4-base overlap one free mismatch. The budget is not too tight; it is
  too *loose* for short overlaps. `floor(0.20 × 5) = 1` lets a 5-base overlap
  match with an error, which fires 4× more often by chance than an exact 4-mer.
- **Raising `trim_require` to 12** costs ~3.3 redundant bases per pair, which is
  the wrong trade for this corpus.

Note the exact-match rule slightly *improves* true-overlap detection (4-base
overlaps 97% → 99%), because the scan is less likely to stop early on a spurious
hit before reaching the real offset. The fix is a strict improvement on both axes
at `exact_below = 10`.

### 3.5 Existing ZNA files have unreliable `IS_RC`

`zna.smk:102` passes `--shuffle`, and before this release shuffle re-derived
orientation, leaving `IS_RC` uncorrelated with which mate was flipped. Combined
with §3.1, **every existing ZNA file must be re-encoded.** Sequences are
unaffected; only the flags are wrong.

### 3.6 Minor

- **Dead parameter.** `zna.smk:93` defines `npolicy = "drop"` but never passes
  `--npolicy`. Behaviour is unchanged (same default) but it is misleading.
- **BAM round-trip, unverified.** Soft clips are retained in `SEQ` and
  secondary/supplementary records are skipped by default, but a **hard-clipped
  primary** alignment would truncate the read and could remove base 0. Worth
  confirming the aligner emits no hard-clipped primaries — this is the one
  remaining unverified link in the chain.

## 4. What the chain already gets right

Verified against the real `hulkrna.merge` code, 20 000 simulated fragments,
inserts 45–300 nt, read lengths 100/150, with adapter read-through:

| Property | Result |
|---|---|
| **Merged read is exactly the fragment** (both edges true boundaries) | **15892 / 15892** |
| **Paired mate's base 0 is still a true boundary** | **4002 / 4002** |

And by inspection:

- **No 5' trimming in the merge tool.** All three construction paths take from the
  start — `s1[:olen]`, `s1[:len1_take]`, `s2[:keep2]` — so base 0 always survives.
  Overlap base correction (`pairs.py:34`) changes base *values*, never positions.
- **Merged reads are in R1's frame**, exactly what `zna`'s single/merged read1
  normalization rule assumes (the v0.3.1 fix).
- **No orphans are produced.** `process_pair` is all-or-nothing for unmerged pairs
  (`pairs.py:196-198`), and upstream `samtools fastq -s /dev/null -0 /dev/null`
  discards singletons. `pairs.py:147-152` documents precisely the failure mode this
  avoids.
- **Pair adjacency survives `--processes N`**: `_iter_chunks` splits on whole pairs
  and each worker's blob is written atomically.
- **Merged reads have the `/1`,`/2` suffix stripped**, so `zna`'s interleaved
  parser classifies them as singles.
- **The salmon → zna strand mapping is correct**: `salmon.smk:132-145` maps
  ISF → `--read1-sense --read2-antisense`, ISR → `--read1-antisense --read2-sense`.
- `--seq-len-bytes 2` is correct for merged reads over 255 bp; `min_read_length:
  40` is applied identically by fastp and the merge tool.

**Consequence for the plan.** Because the merge tool emits *only* full-fragment
merged singles and intact pairs, **every unpaired record in this pipeline is
genuinely a full-fragment read.** The "genuine single-end read" hazard does not
arise, so the §5.2 declaration is exact rather than a guess, and `khorana`'s
`(True, True)` for unpaired records is correct — once §3.1 lands.

## 5. The zna fixes (0.3.3)

Already implemented and green (270 tests): `records(with_rc=True)`,
`ZnaWriter(preserve_normalization=True)`, the shuffle and ZNA→ZNA re-encode
orientation fixes, and the re-encode guards. The rest of this section is new work.

### 5.1 Labeled encoding must honour the input strategies — blocking

`stream_inputs_labeled` must not have its own reader. Refactor so labeled encoding
reuses `stream_inputs`' strategies (interleaved, two-file, single-end, ZNA) and
attaches labels to the records they yield. Until that lands, `--label-defs` with
`--interleaved` or two files must be a **hard error** rather than silent
corruption: a loud wrong answer is recoverable, this is not.

### 5.2 `IS_FULL_FRAGMENT` record flag (bit `0x10`)

```python
class ZnaRecordFlags(IntFlag):
    IS_READ1 = 1
    IS_READ2 = 2
    IS_PAIRED = 4
    IS_RC = 8
    IS_FULL_FRAGMENT = 16      # both edges are true fragment boundaries
```

Bits 4–7 are unused and existing readers extract specific bits, so this is
forward-compatible: **no format version bump**, old readers read new files
unchanged, new readers see the bit clear on old files — the safe reading.

- **Paired records: detected automatically.** A full-overlap pair's mates span the
  identical interval, so they are exact reverse complements. Testing
  `len(s1) == len(s2) and s2 == revcomp(s1)` detected **10051/10051** with
  **0/9949** false positives — one comparison per pair on the equal-length path.
  (With `hulkrna-merge` upstream these are already merged, so this protects other
  pipelines, including fastp-merge users.)
- **Unpaired records: declared** via `--treat-unpaired-as-merged`, default off. For this
  pipeline it is always correct to pass it (§4).

### 5.3 `records(with_ends=True)` — hand over the answer, not the ingredients

```python
    for seq, is_paired, is_read1, is_read2, has_start, has_end in \
            reader.records(with_ends=True):
```

with `has_start = is_full_fragment or not is_rc` and
`has_end = is_full_fragment or is_rc`. The derivation is where consumers go wrong —
the two-liner in `RC_FLAG_AND_REENCODE_FIX.md` §9 is itself incomplete for
full-fragment records — so `zna` should own it. `with_rc` stays for raw access;
default widths (4, 5 labeled) stay pinned.

### 5.4 Carry the new bit through every copy path

`preserve_normalization` carries `IS_RC` through re-encode and shuffle;
`IS_FULL_FRAGMENT` must ride along identically — `write_record`, `_shuffle.py`'s
read and write sites, and the CLI re-encode path. Silently clearing it would
downgrade full-fragment records to one-ended ones.

### 5.5 Pair-atomic `--npolicy drop`

Consume the record stream in units (a pair, or a single) and apply the drop
decision to the whole unit, so a dropped mate takes its partner with it. This
matches `hulkrna-merge` and fastp, and removes the only route to §6.2.

### 5.6 Warn on the `decode` → `encode` round trip

`zna decode` emits the **normalized** frame unless `--restore-strand` is given;
re-encoding that with `--strand-normalize` normalizes twice and desynchronizes
`IS_RC` from the bases — measured **200/400** records with their indicated edge on
a truncation point, byte-identical headers, exit 0, silent. Warn on decode.

## 6. Lower priority (zna)

Verified real; none can place an interior token.

1. **Orphan mates land in the wrong frame.** `_stream_single_end` (`cli.py:694`)
   emits every unpaired read as neither R1 nor R2, so the read1 rule applies —
   wrong for a `--unpaired2` file. Does not arise in this pipeline (§4); inverts
   polarity only, never the boundary edge.
2. **The unstranded pair window tests `IS_PAIRED` but not `IS_READ2`**
   (`_pycodec.py:143`, `_accel.cpp:346`), so two consecutive R1s can be
   co-oriented as mates, and the backends then diverge on identical input
   (`[T,F,T,F,T]` python vs `[T,F,F,T,T]` accel). §5.5 removes the only route to it.
3. **`read1_antisense` and `read2_antisense` both true** is silently accepted and
   reverse-complements both mates into opposite frames; the existing warning covers
   only the both-false mirror. `khorana` already rejects this header.
4. **The pure-Python backend silently corrupts an invalid base** at string position
   3, 7, 11 … into `TTTT` instead of raising (`_pycodec.py:67-70`); the C++ backend
   raises correctly.

## 7. Compatibility

- On-disk format unchanged in shape; bit `0x10` occupies an unused bit. No version
  bump, old readers unaffected.
- Default tuple widths unchanged (4, or 5 labeled), still pinned by test.
  `with_ends` is opt-in and independent of `with_rc`.
- Existing files stay readable; they simply never report `has_start and has_end`,
  which is the safe, incomplete reading. But per §3.1 and §3.5 they carry wrong
  `is_paired` and unreliable `IS_RC`, so **they must be re-encoded regardless.**
- §5.1 changes what labeled interleaved encoding *produces*. That is a bug fix, but
  files already written cannot be repaired — re-encode from the merged FASTQ.

## 8. zna checklist (0.3.3)

Ordered; each step independently testable.

1. **§5.1 labeled + interleaved.** Blocking — everything downstream is moot until
   pairing survives encoding. Land the hard error first if the refactor slips.
2. `IS_FULL_FRAGMENT` in `ZnaRecordFlags`; `write_record(..., is_full_fragment=)`.
3. Full-overlap auto-detection in `_pycodec.encode_block`, `_accel.cpp
   encode_block`, and `_accel.cpp encode_block_labeled` — all three in lockstep.
   **Needs a rebuild** (`pip install -e .`); no other step does.
4. `records(with_ends=True)`.
5. Thread the bit through `preserve_normalization`, `_shuffle.py`, CLI re-encode.
6. `--treat-unpaired-as-merged`; pair-atomic `--npolicy drop` (§5.5).
7. `zna decode` round-trip warning; `zna inspect --counts` reports a full-fragment
   tally and a (mate × is_rc) cross-tabulation — the only tally that verifies a
   file's geometry before training on it.
8. README, CHANGELOG, version already at 0.3.3.

Tests, alongside `TestRcFlagAndReencode`:

1. `test_labeled_interleaved_preserves_pairing` — the §3.1 regression test: a mixed
   interleaved stream with labels must produce the same flags as without labels.
   **Must fail before the fix.**
2. `test_labeled_paired_files_preserves_pairing`.
3. `test_full_overlap_pair_detected` — insert ≤ read sets the bit on both mates;
   insert > read does not.
4. `test_full_fragment_flag_survives_reencode_and_shuffle`.
5. `test_with_ends_geometry_contract` — the load-bearing test: across dUTP,
   fr-secondstrand and unstranded, every edge reported by `with_ends=True` is a
   true fragment boundary **and** every true boundary is reported. Both directions.
6. `test_with_ends_default_width_unchanged` — pins 4 / 5 / 6.
7. `test_npolicy_drop_is_pair_atomic`.
8. `test_decode_reencode_warns`.
9. Cross-backend equivalence for 3 and 5.

## 9. hulkrna checklist

1. **§3.4 overlap fix**: in `lib/hulkrna/merge/overlap.py` `_scan`, add an
   `exact_below` parameter (default 10) forcing `dl = 0` for short overlaps.
   Keep `trim_overlap_require: 3`. Expose as `--overlap-exact-below`.
2. Add a merge-tool regression test asserting the spurious-overlap rate on
   unrelated random pairs stays below a threshold (e.g. < 0.1% at 20 000 pairs).
   This is the guard that would have caught §3.4, and it is cheap.
3. `zna.smk`: pass `--treat-unpaired-as-merged`; pass `--npolicy {params.npolicy}` or
   delete the dead param (§3.6).
4. Verify no hard-clipped primary alignments reach `samtools fastq` (§3.6).
5. Add an end-to-end test: merge → encode → read back, asserting every merged
   record reports both ends and every paired record exactly one. **This is the test
   that would have caught §3.1**, and it belongs here because it is the only place
   that sees both sides of the seam.
6. Re-run the pipeline and re-encode all libraries once zna 0.3.3 is released
   (§3.5).

## 10. khorana checklist

1. `records()` → `records(with_ends=True)` at `data.py:526` and `data.py:640`;
   `record[:4]` unpacking stays valid.
2. `_assign_endpoints` takes `has_start, has_end` and returns them directly. Its
   body collapses; `header_info` is no longer needed for endpoint assignment,
   though its validation (`data.py:205-217`) is worth keeping — it already rejects
   `strand_normalized=False` and the both/neither-antisense headers.
3. Do not re-derive endpoints from `is_paired`. Keeping `(True, True)` for unpaired
   records is correct *for this pipeline* (§4) but wrong in general — for any
   genuine single-end or orphan record — and `with_ends` makes the distinction
   explicit.
4. Add an assertion that the file was written by zna ≥ 0.3.3, or that
   `with_ends=True` is available, so a stale-file run fails loudly rather than
   training on wrong labels.
5. Interim, if training must proceed before 0.3.3: `has_sos = not is_rc;
   has_eos = is_rc` applied to *every* record removes the catastrophic case using
   `with_rc` alone, at the cost of the `bookend` crop (it never yields both ends,
   so `data.py:352` stops firing).

## 11. Ordering and gates

The three repos have a hard dependency order. Do not overlap these.

| Gate | Condition to pass |
|---|---|
| **G1** | zna 0.3.3 released with §8.1–§8.7; §8 test 1 fails on 0.3.2 and passes after. |
| **G2** | hulkrna §9.1–§9.5 merged; the end-to-end test (§9.5) green against 0.3.3. |
| **G3** | Pipeline re-run complete; `zna inspect --counts` on a sample of new files shows sane (mate × is_rc) and full-fragment tallies. |
| **G4** | khorana §10 merged, reading only files that passed G3. Then retrain. |

G3 is the gate that matters most: it is the first point at which anyone can look
at a produced file and see whether the geometry is right. Every defect in §3 was
invisible precisely because no such check existed.

Note for anyone validating geometry by locating sequences in a reference: carry a
fragment id through as a ZNA **label** instead. An earlier audit produced a
spurious 97.5% failure rate that turned out to be birthday collisions on fragment
start coordinates, not a codec defect.
