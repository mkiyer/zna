# `--npolicy`: `trim3` and `random` — implementation plan

Status: **built.** Written 2026-08-13 after the consensus/N-rescue work landed and an
audit of `--npolicy` against both codec backends found four defects, one of them on the
default path; completed the same day with §5's encode summary and D5's per-record
provenance.

| section | state |
|---|---|
| D1 policy set `{trim3, random}` | **built** — same flag, same values, on `zna merge` and `zna encode` |
| D2 no `--max-n` | **decided**, not built by design |
| D3 stumps allowed | **built** |
| D4 / D4a option C + the coverage retry | **built**, both backends byte-identical on N-bearing input |
| §8.2 alphabet strictness | **built** — only `N`/`n` substitutable, everything else raises with character and offset |
| §8.3 seeded position-derived streams | **built** — covers `--npolicy random` *and* the unstranded RC coin; verified reproducible, backend-identical, block- and chunk-size independent |
| §5 reporting | **built** — run-level counters and a >1% warning on both sides of the seam; `zna encode` now reports its policy the way `zna merge` always has |
| D5 per-record provenance | **built** — header tokens for reading, `ZN:i:<bits>` for the corpus |

A defect the audit found that is not in §1: `_accel` measured record length in UTF-8
*bytes*, so a non-ASCII character corrupted the lengths column (`ACGéT` stored as length
6, decoding to `ACGAAT`). Fixed by §8.2's strictness; both backends now raise.

The goal is a flag a user can predict without reading the source: **two policies, no
tuning parameters, and a report of what happened.**

---

## 1. What is wrong today

Measured directly against both backends (`encode_block` with a single `ACGN` record):

| policy | `_pycodec` | `_accel` | |
|---|---|---|---|
| *(unset)* | `ValueError` | `ValueError` | agree |
| **`drop`** *(the default)* | `ValueError` | **encodes → `ACGA`** | **diverge** |
| `random` | `ACGT`, **non-deterministic** | `ACGT`, deterministic | **diverge, and unreproducible** |
| `A`/`C`/`G`/`T` | substitutes `N` only | substitutes **any** unencodable char | diverge on IUPAC |

Four defects:

1. **`drop` silently substitutes `A` on the compiled backend.** Root cause is one
   reasoning slip in a dispatch chain ([_accel.cpp:553-558](../../src/zna/_accel.cpp#L553-L558)):
   ```cpp
   const bool has_npolicy = !npolicy.empty();       // "drop" is non-empty -> true
   if (npolicy == "C") … else if ("G") … else if ("T") … else if ("random") …
   // else: A (default) or empty (will throw on N)  <-- "drop" lands here too
   ```
   It is masked in the CLI because `drop_n` removes N-containing fragments first — but
   `drop_n` tests only `N`/`n`, so an IUPAC code slips past it, and any library caller
   using `ZnaWriter(npolicy="drop")` gets silent substitution. Every unrecognised policy
   string behaves the same way.
2. **`random` is not reproducible.** `_pycodec` uses the unseeded global `random`, so the
   same input produces a different file on every run — in a tool whose stated promise is
   byte-identical output.
3. **The two backends disagree on `random`** even when both are deterministic.
4. **`--npolicy` is not an N policy.** `_accel` implements an *anything-unencodable*
   policy; `_pycodec` an N-only one.

---

## 2. Decisions

### D1. The policy set becomes `{trim3, random}`

`drop` and `A`/`C`/`G`/`T` are removed.

- **`drop` is subsumed by `trim3`.** Trimming at the first `N` and letting the existing
  length filter discard what is left is the same outcome for a read that is mostly N, and
  strictly better for a read with one N near its 3' end — which is where instrument no-calls
  actually sit.
- **`A`/`C`/`G`/`T` inject a composition bias** correlated with low quality, and buy
  nothing `random` does not: every N becomes the same base, so the corpus gains a spurious
  A-enrichment exactly at low-quality positions. 15 test call-sites use `npolicy='A'` only
  because it is deterministic; once `random` is properly seeded they can use that.

### D2. No `--max-n`

An audit proposed bounding the number of no-calls in a retained read, as fastp and
cutadapt do. **Rejected**, and the reason is that it is a *substitution* knob:

- Under `random` it means something — "at most K of 150 bases are guesses."
- Under `trim3` it means nothing. The output contains **zero** N by construction, and what
  determines the loss is the *position* of the first N, not the count. A read with one N at
  position 3 and a read with fifty N starting at position 3 produce the identical 3-base
  output.

There is no intuition for how a user would set it, so it does not exist.

### D3. Stumps are allowed. The tool does not second-guess the consumer

`trim3` may cut a read to one base, or to none. Verified that ZNA stores both:

```
lengths written -> read back:  [0, 1, 3, 144]
```

`zna merge --min-read-length` (default 40) already discards stumps **pair-atomically**, so
the sensible default exists without a new knob. `zna encode` standalone has **no length
floor at all**, so `zna encode --npolicy trim3` will write 0-base records rather than
dropping them. That is documented, not fixed: a core read-processing tool should report
what it did and let the consumer decide.

### D4. A pair with a surviving `N` is not merged; it is trimmed and reconsidered

Settled by a worked example rather than by argument. Take a 220 bp fragment, 2x150 reads,
and an `N` at R1[40] — outside the 80-base overlap, so rescue cannot reach it:

| | emitted | `IS_FULL_FRAGMENT` | fragment bases kept |
|---|---|---|---|
| A. trim the merged record, keep the flag | 40 bp | 1 — **false**: the 3' terminus was at 219 | 40 of 220 |
| B. trim the merged record, clear the flag | 40 bp | 0 — honest | 40 of 220 |
| **C. don't merge; trim and keep the pair** | 40 bp + 150 bp | 0, 0 — honest | **190 of 220** |

**C wins on both axes**, and the second one was the surprise: merging and *then* trimming
throws away the entire tail of the fragment, including all 150 bases R2 independently
observed, because the merged record is one indivisible sequence and the `N` sits early in
R1. C keeps 4.75x the sequence *and* leaves both fragment termini intact.

The geometric reason A is wrong is worth stating once, because it is not obvious: **a
merged record's 3' end is R2's 5' end.** Verified in both geometries —
`M[L-1] == complement(R2[0])` for a normal overlap and for read-through. So "trim the
merged record's 3' end" means "discard R2's 5' anchor", which is exactly the anchor the
trim-3'-only rule exists to protect. On an unmerged *pair* the two anchors live in two
separate records and trimming inward touches neither; that is why trim3 is safe there and
only there.

### D4a. After trimming, retry the merge — using the original overlap, not a re-scan

Trimming removes *interior* bases and leaves both 5' anchors, so R1' covers `[0, k1)` and
R2' covers `[L-k2, L)`. The pair therefore still tiles the fragment iff

```
k1 + k2 >= L
```

and when it does, the reconstruction **is** the fragment, exactly and N-free — verified:

```
 N at R1[]   k1    k2   k1+k2   >=L?   outcome
        40   40   150     190  False   keep pair (30-base gap in the middle)
        80   80   150     230   True   MERGE -> 220 bp, == fragment: True, N-free: True
       120  120   150     270   True   MERGE -> 220 bp, == fragment: True, N-free: True
       145  145   150     295   True   MERGE -> 220 bp, == fragment: True, N-free: True
```

**The retry reuses the original shift. It does not re-scan.** trim3 cuts 3' ends, and for
a normal overlap the overlap *lives* at the 3' ends, so trimming preferentially destroys
the evidence a re-scan would need:

| N at R1[] | pre-trim overlap | post-trim overlap |
|---:|---:|---:|
| 40 | 80 | 0 |
| 80 | 80 | 10 |
| 120 | 80 | 50 |
| 145 | 80 | 75 |

At N@80 the overlap collapses from 80 bases to 10, which scores ~19.9 bits — under the
28-bit threshold. A re-scan would refuse a merge there was 80 bases of evidence for a
moment earlier. Removing a no-call must not make the tool less confident about a fragment
length it has already established.

The shift is a claim about **fragment length**, and trimming interior bases does not change
a fragment's length. So the retry is a *coverage test*, not an inference:

```
1. scan                -> shift s, score, L = s + len2
2. rescue              -> an N opposite a real call takes that call
3. trim3               -> cut each read at its first surviving N, giving k1, k2
4. decide, reusing s:
     k1 + k2 >= L  and score >= t_merge      -> MERGE, full fragment, IS_FULL_FRAGMENT=1
     k1 + k2 >  L  and score in [t_trim, t_merge) -> symmetric trim on the new lengths
     otherwise                               -> keep the pair
5. min_read_length, pair-atomic
```

**A knob deliberately not added:** a minimum *residual* overlap at step 4. It sounds
prudent, but the evidence for `s` was already gathered; demanding it again on less data is
the same mistake as re-scanning. The pre-trim score is the evidence, `k1 + k2 >= L` is the
geometry — two questions, each asked once.

### D5. Transparency goes in two places, because they serve different readers

`zna merge` writes a FASTQ; **ZNA does not store headers.** So a header token is
transparency for the *intermediate* only, and under `zna encode --merge-pairs` (0.5.0)
there is no intermediate at all.

- **Header tokens** on the merge FASTQ, for human inspection: `trim3_<n>` for bases cut
  (`subn_<n>` under `random`), `rescued_<n>` for no-calls recovered from the mate,
  beside the existing `merged_<n1>_<n2>`. Safe against label extraction — a token with no
  colon is skipped by the tag parser, exactly as `merged_150_87` is:
  ```
  f12  ZI:i:7 ZN:i:6 trim3_12 rescued_1 merged_90_0   ->   labels (7,)   unaffected
  ```
  They are emitted on **every** outcome, not only the merged one. A kept pair whose R1
  lost 12 bases to `trim3` is emitted shorter, and without a token nothing on the record
  says why — which is the gap this exists to close. A record nothing happened to is
  emitted byte-unchanged and stays zero-copy in the compiled backend: 1,332,353 of
  1,418,525 records on the 1M library at a 1.5% no-call rate.
- **A provenance label** for the corpus, since that is the only thing that survives
  encoding: one `C` (uint8) column, one byte per record.

  It is carried as an ordinary SAM tag, `ZN:i:<bits>`, and read by the `--label`
  machinery that already ships — `--label provenance:C:ZN`. That settles the question
  §7 left open ("it adds a column to every labeled file"): it does not. Declaring the
  column is opt-in per file, an absent tag resolves to 0 through the label path's own
  missing-value handling, and there is no provenance-specific code in the encoder.

  **Bits: trimmed 1, rescued 2, N-trimmed 4, N-substituted 8 — and deliberately no
  "merged".** The plan above listed one. It was dropped on measurement: `merged` already
  has two homes (the `merged_` token, and `IS_FULL_FRAGMENT` in the flag byte), and
  spending a bit on it put ` ZN:i:1` on ~82% of emitted records to repeat what two other
  places already said. Every remaining bit is a fact with nowhere else to live —
  `trimmed` above all, since a trimmed pair is emitted as an ordinary pair and no ZNA
  flag distinguishes it from one kept whole. The consequence worth knowing:
  `IS_FULL_FRAGMENT` is set only under `zna encode --treat-unpaired-as-merged`, so an
  encode that omits that flag records neither — which is what omitting it means.

  The vocabulary is shared with `zna encode --merge-pairs` (0.5.0), which computes the
  same `PairResult` and writes the same bits with no FASTQ in between.
- **Run-level counters** in `--merge-json` and the encode summary (§5).

---

## 3. The specification

`--npolicy` governs exactly one thing: **what happens to `N` and `n`.**

**The encodable alphabet** is `A C G T a c g t`, folded to uppercase. `N`/`n` are the
ambiguity characters. **Every other character is invalid under every policy** — all IUPAC
codes (`R Y S W K M B D H V`), `.`, `-`, digits, whitespace, non-ASCII — and raises

```
ValueError: invalid character 'R' at offset 37 in sequence
```

in **both** backends. No policy can cause a non-`N` character to be silently encoded. This
is the fix for defects 1 and 4.

| policy | behaviour |
|---|---|
| **`trim3`** *(default)* | Cut each read at its first `N`, keeping `[0, first_N)`. 3' only, never 5', so base 0 remains a true fragment boundary. Applied **after** overlap rescue, as a last resort. |
| **`random`** | Every `N` becomes a base drawn from a seeded PRNG (§4), identical in both backends. |
| *(unset — library API only)* | `ValueError` on any `N`. The CLI always passes a policy. |

The policy set `{None, "", "trim3", "random"}` is closed; anything else raises at
`ZnaWriter` construction rather than defaulting.

### Ordering — this is the part that matters

```
1. scan for the overlap
2. rescue: an N opposite a real call takes that call          (shipped in 0.4.0)
3. decide merge / trim / keep, and build the records
4. trim3: cut each emitted record at its first SURVIVING N
5. length filter (min_read_length), pair-atomic
```

Rescue precedes trimming because an N recovered from the mate costs nothing, while
trimming discards every base after it. Step 4 is on the **emitted records**, not the input
reads: trimming before the scan would forfeit the rescue, and trimming between 3 and 4
would change the geometry the merge decision was made on.

This is also why `trim3` belongs in the **merge** path and cannot live in the codec: the
codec sees one record and knows nothing about a mate.

---

## 4. Making `random` reproducible

Harder than it looks, and the reason is worth stating: `_accel` substitutes during
**packing**, advancing a `xorshift32` once per substituted base, in *stored* (post-RC)
order, with the state **carried across records within a block**
([_accel.cpp:439-443](../../src/zna/_accel.cpp#L439-L443), seeded `0xDEADBEEF` per
`encode_core` call). `_pycodec` substitutes on the string before packing, per record,
using the global `random`.

The orders happen to agree — both are post-RC, left to right — so mirroring is possible.
Two options:

**(a) Mirror the existing stream.** Give `_pycodec` the same `xorshift32`, threaded
through its per-record loop. Smallest change; achieves determinism and backend agreement.
Leaves one wart: the stream is seeded per *block*, so the substituted bases depend on
`--block-size`. Document it.

**(b) Make it stateless and position-derived.** `base = mix(seed, record_index, offset) & 3`,
evaluated in stored order. Independent of block size, of batching, and of the
`do_random_rc` coin stream; trivially identical across backends. Needs a base record index
threaded into `encode_block`, as `merge_chunk` already does with `base_index`, plus wiring
`--seed` (which `zna encode` already has, currently used only by `--shuffle`).

**Recommendation: (b).** If `random` is one of only two policies it should not change
because someone tuned `--block-size`. (a) is the fallback if the signature change proves
disruptive.

---

## 5. Reporting

Add to the encode summary and to `--merge-json`, for every run:

```
[ZNA] N policy 'trim3': 5,823 reads trimmed (2.91%), 1.75 Mbase removed;
                        11,204 no-calls rescued from the mate;
                        412 reads dropped below --min-read-length
```

Warn above a threshold (say 1% of reads), in the style of the existing
strand-normalization warning. The motivating case: a dark cycle — one failed tile, an
ordinary Illumina occurrence — took a 400,000-read library to 918 records under the old
`drop` policy, with stderr reporting `Done.` A default whose worst case is silent deletion
of a library is a bad default however good its median case. `trim3` largely removes that
failure mode; the report removes the *silence*, which is the part that actually bit.

---

## 6. Build order

1. **Alphabet strictness, both backends.** Substitute only for `N`/`n`; raise on every
   other invalid character; close the policy set and reject unknown strings at
   `ZnaWriter` construction. Fixes defects 1 and 4. *Test:* the behaviour matrix of §1,
   asserted cell for cell across both backends.
2. **`random`, per §4.** *Test:* same input encodes byte-identically across runs, across
   backends, and — under (b) — across `--block-size`.
3. **`trim3` in the merge path**, at step 4 of the ordering. *Test:* a read with an N at
   position *k* emits exactly *k* bases; no emitted record contains `N`; base 0 unchanged;
   pair atomicity preserved through the length filter.
4. **Remove `drop` and `A`/`C`/`G`/`T`**; migrate the 15 test call-sites to `random`.
   Update the help text and README.
5. **Reporting** (§5).
6. **Transparency** (D5): header tokens, then the optional provenance label.

Steps 1–2 are independent of 3–6 and fix shipping bugs; they can land first.

---

## 7. Open

- **D4** — keep `IS_FULL_FRAGMENT` on a trimmed merged record (decided), or adopt the
  reconciliation in D4 and decline to merge pairs with a surviving N.
- **Is `trim3` or `random` the default?** `trim3` never invents a base, which suits a
  corpus whose value is that it reports what the instrument reported. `random` never loses
  a base. The plan assumes `trim3`.
- ~~**Does the provenance label ship in 0.5.0 or later?**~~ **Closed: it ships in
  0.4.0.** The premise was wrong — carried as a `ZN:i:` tag and read through the
  existing `--label` path, it adds a column only to files that ask for one. See D5.

---

## 8. The codec side

`trim3` in the merge kernel is built (§D4a). This section is the encoder.

### 8.1 What is wrong, measured over all 256 byte values

```
policy    input class     _pycodec        _accel
''        N / IUPAC / other   RAISE       RAISE        agree
'A'       N n                 A           A            agree
'A'       IUPAC / other       RAISE       A            DIVERGE
'random'  N n                 non-determ. T            DIVERGE
'random'  IUPAC / other       RAISE       T            DIVERGE
```

And the reproducibility hole is **wider than `--npolicy`**: `_pycodec` uses the unseeded
global `random` for the unstranded strand-normalization coin as well.

```
unstranded --strand-normalize, IS_RC column:
  pycodec   4 runs identical: False        <- non-reproducible
  accel     4 runs identical: True
  pycodec vs accel identical: False        <- backends disagree
```

Worse, on the shipping backend the coin is **block-dependent** — the same 16 records get a
different orientation assignment depending on how they are batched:

```
one block of 16 : [13, 6, 5, 14, 5, 14, 13, 6, 5, 14, 13, 6, 5, 14, 13, 6]
two blocks of 8 : [13, 6, 5, 14, 5, 14, 13, 6, 13, 6, 5, 14, 5, 14, 13, 6]
```

So `--block-size` silently changes a normalized corpus. Both random streams have the same
defect and take the same fix, which is the argument for the position-derived design over
mirroring the existing per-block xorshift.

### 8.2 Target semantics

The encodable alphabet is `A C G T a c g t`, folded to uppercase on storage. `N`/`n` is
the ambiguity character. **Every other byte is invalid under every policy** and raises

```
ValueError: invalid character 'R' at offset 37 in sequence
```

identically in both backends.

| codec policy | `N` / `n` | any other non-ACGT |
|---|---|---|
| `""` / `None` *(strict; library default)* | `ValueError` | `ValueError` |
| `"trim3"` | cut the record at it, keep `[0, first N)` | `ValueError` |
| `"random"` | a base from the seeded stream (§8.3) | `ValueError` |

The set is closed: anything else raises at `ZnaWriter` construction rather than
defaulting to a substitution.

### 8.3 One seeded, position-derived stream, for both random decisions

Stateless, so it cannot depend on block size or batching. Integer-only, so the two
languages cannot drift.

```c
static inline uint64_t zna_mix64(uint64_t x) {          // splitmix64 finalizer
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}
```

```python
_M64 = 0xFFFFFFFFFFFFFFFF

def _mix64(x: int) -> int:                              # the same function
    x = (x + 0x9E3779B97F4A7C15) & _M64
    x = ((x ^ (x >> 30)) * 0xBF58476D1CE4E5B9) & _M64
    x = ((x ^ (x >> 27)) * 0x94D049BB133111EB) & _M64
    return x ^ (x >> 31)
```

Two independent streams, keyed by the record's **global** index `r` and, for
substitution, the offset within the **stored** (post-RC) sequence:

```
rc_coin(r)        = mix64(seed + 0x9E3779B97F4A7C15 * (r + 1))        & 1
sub_base(r, off)  = mix64(seed + 0xBF58476D1CE4E5B9 * (r + 1)
                               + 0x94D049BB133111EB * (off + 1))      & 3
```

Different multipliers, so the two never shift each other. Substitution is drawn **after**
reverse-complement, in stored coordinates, because that is the only coordinate both
backends already share — and it means a substituted base is never itself complemented.

**Plumbing.** `encode_block` gains `rng_seed: int` and `base_index: int`, mirroring the
`base_index` pattern `merge_chunk` already uses. `ZnaWriter` tracks a running
`_records_written` and passes it as `base_index`; `rng_seed` comes from the existing
`zna encode --seed` (today used only by `--shuffle`). Substitution therefore happens
exactly once, at first encode: `copy_records()` / `write_copy()` carry the flag byte and
the packed bases verbatim and never re-substitute.

### 8.4 Removing `drop` and `A`/`C`/`G`/`T`

`drop` was never a codec policy — it is a record-admission filter, and letting the string
reach the codec is the direct cause of the worst defect above. The CLI keeps its
pair-atomic fragment filter; the codec is simply never handed `"drop"` again.

`A`/`C`/`G`/`T` go. Test call sites using `npolicy='A'` do so only because it is
deterministic and migrate to `'random'` with a fixed seed.

### 8.5 Landing order

1. **Alphabet strictness**, both backends. Only `N`/`n` substitutable; everything else
   raises with character and offset. *Suite must stay green.*
2. **Seeded streams**, both backends, both decisions, plus the `ZnaWriter` plumbing.
   *Corpus-affecting: every unstranded normalized file changes.*
3. **`trim3` in the codec** — truncate at the first `N` before packing.
4. **CLI**: `--npolicy {trim3,random}`, default `trim3`; remove `drop`/`A`/`C`/`G`/`T`;
   migrate call sites.
5. **Reporting** and docs.

Steps 1–2 fix shipping bugs and are independent of 3–5.
