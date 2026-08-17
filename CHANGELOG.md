# Changelog

All notable changes to the ZNA project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and
version numbers follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html) —
**except that ZNA is alpha and has broken the on-disk format in a patch release** (0.4.1)
rather than hold the change. Read the notes, not the number: a release that breaks your
files says so in its first paragraph.

## [Unreleased]

### Fixed — `update-conda-sha.sh` silently rewrote nothing on macOS

The rewrite used `sed -E "s|^(\s*)sha256:.*|…|"`. **BSD sed — the macOS default — has no
`\s` in ERE and treats it as a literal `s`**, so the pattern required the line to *begin*
with the letter s and matched nothing on a normally-indented recipe. GNU sed on Linux
accepts `\s`, so it worked everywhere except the machine releases are actually cut from.

BSD *grep* does accept `\s`, which is what made this confusing: every read in the script
(the current hash, the one-sha256-line guard, the post-write verify) found the line
correctly, and only the write failed. Now `[[:space:]]` throughout.

Caught by the script's own post-write verify — the check added in 0.4.0 for exactly this,
"confirm the edit actually landed rather than trusting sed". Without it the recipe would
have gone to bioconda carrying 0.4.0's hash under 0.4.1's version, and bioconda CI would
have been the first thing to notice.

## [0.4.1] - 2026-08-17

> **This is a patch release that breaks the on-disk format.** That is deliberate and it
> is not what the version number would normally tell you: ZNA is alpha, has no legacy
> support to keep, and the release number was chosen over the convention. Read the next
> section before upgrading — **every existing `.zna` file must be re-encoded from
> FASTQ**, and any downstream consumer reading them (hulkrna, khorana) will hit the
> version error on first contact rather than reading anything.

### Changed — blocks are fragment-complete (breaking; format version 3)

**A fragment's reads are now stored consecutively, R1 immediately followed by R2, and
never span a block boundary.** A block is therefore a self-contained set of molecules,
and any subset of blocks decodes independently of the rest of the file. That is the
premise `blocks(stride=…)` sharding and every block-parallel consumer was already
written against; until now the format did not actually provide it.

- **The on-disk version is 3, and version 2 files are refused** with an error that says
  why. They can hold split fragments, so they cannot be read as though they cannot.
  Re-encode from FASTQ.
- **`ZnaWriter` enforces the contract on every write path** — `write_record`,
  `write_records` (both modes), and `write_copy`. A paired R1 not followed by its R2, an
  R2 with no R1 in front of it, or a stream ending on an R1 raises `ValueError` naming
  the record. `zna encode` and `zna shuffle` already grouped their output this way, so
  no CLI behaviour changes; direct API callers writing half fragments now find out.
- Unpaired and merged reads are one-record fragments and need nothing special.

**What was wrong.** The writer held a block open for a pending R1 only when the header
asked for *unstranded* strand normalization, where the codec's pair detection requires
it. Every other configuration — no normalization, stranded normalization, the bulk
`write_records` path, and the `write_copy` path that `zna shuffle` and ZNA → ZNA
re-encode run through — split fragments freely. Measured at the default 4 MiB block
size over 600k pairs of varied length: **44% of blocks ended mid-fragment.** Fixed-length
reads hid it completely, because a constant per-record size estimate lands every
boundary on an even record count.

**Two defects fixed along the way.**

- *Unbounded memory.* The old deferral was not bounded: it held the block open without
  requiring the mate to arrive, so a run of consecutive paired R1s buffered the entire
  stream and wrote nothing. Fifty thousand such records produced a 16-byte file — the
  header alone — with every record still in memory. Enforcing the contract bounds the
  hold at one record, because a fragment is at most two records long.
- *Backend divergence.* For a pair split across blocks under unstranded normalization,
  the Python codec gave the orphaned R2 its own coin flip and the C++ codec gave it
  none, so the two backends would have written different files. Unreachable in practice
  only because that was exactly the case the deferral prevented. Both backends now share
  one rule and raise if a paired R1's mate is not in the block.

Dead code removed: the second arm of the old deferral in `write_record` was unreachable
(the size estimate grows monotonically until a flush zeroes it, so it could not be below
`block_size` while an R1 was pending). Verified across 160k records over the config
matrix before deleting.

## [0.4.0] - 2026-08-13

Adds `zna merge`, an overlap merger for paired-end reads, with a compiled
C++ kernel, and repairs six defects found while auditing the merge → encode
seam — one of them a heap overflow, four of them silent. No on-disk format
change.

**Read the "Fixed (silent)" section before upgrading.** Two of those fixes change
what a previously-"successful" run does: an encode that swallowed a malformed SAM
tag now fails, and a ZNA → ZNA copy through `records(with_ends=True)` now raises
instead of quietly dropping flag bits.

### Changed — `--npolicy` (breaking; affects the corpus)

The policy set is now **`{trim3, random}`**, default `trim3`. `drop` and `A`/`C`/`G`/`T`
are removed.

- **`trim3`** cuts a read at its first no-call, keeping `[0, first N)`. 3' only, so base 0
  stays a true fragment boundary however short the read gets. In `zna merge` it runs
  *after* overlap rescue and is followed by a coverage retry: a pair that still tiles the
  fragment (`k1 + k2 >= L`) merges anyway, and the reconstruction is the fragment exactly.
  Measured on 200k pairs with no-calls in 1.5% of reads: **zero N in the output** for 0.22%
  of bases and 0.31% of merges demoted to pairs.
- **`drop` is subsumed.** Trimming and letting the length filter discard the remainder is
  the same outcome for a read that is mostly N, and strictly better for one with an N near
  its 3' end — which is where instrument no-calls sit. `trim3` also never orphans a mate,
  so the fragment-atomic machinery `drop` needed is gone.
- **`A`/`C`/`G`/`T` are removed**: a fixed base injects a composition bias correlated with
  low quality, and buys nothing `random` does not.

**Reported, not silent.** `--merge-json` and the merge summary now say what the policy
did, and warn above 1% of emitted bases — because the failure being guarded against is a
run that eats most of a library and still ends with `done`. On a library with no-calls in
1.5% of reads the two policies read very differently:

```
no-calls: 1110 rescued from the mate; --npolicy trim3  then removed     182,477 bases
no-calls: 1110 rescued from the mate; --npolicy random then substituted   4,844 bases
```

trim3 discards everything after each surviving no-call, so it costs 37x more sequence than
random invents. That trade-off was previously invisible.

**`zna encode` now reports its policy too**, in the same shape and with the same >1%
warning. The same policy on the same library used to be accountable on one side of the
merge → encode seam and silent on the other:

```
[ZNA] Done. Wrote 300 records in 0.00s.
[ZNA] no-calls: --npolicy trim3 removed 11558 bases in 203 records.
```

Skipped when re-encoding a `.zna`, and that is not an oversight — ZNA stores two bits per
base, so a decoded record cannot contain a no-call and no policy ran. Reporting
"removed 0 bases" there would imply one had been applied.

**Four defects fixed with it.** Measured over all 256 byte values x every policy:

- **`--npolicy drop` — the default — silently wrote `A` on the compiled backend.** One
  reasoning slip: `has_npolicy = !npolicy.empty()` made `drop` truthy, and the
  `C`/`G`/`T`/`random` chain did not match it, so it fell to the `A` default. Every
  unrecognised policy string did the same. The set is now closed and rejects.
- **`--npolicy` was not an N policy.** `_accel` substituted for *any* unencodable byte, so
  an IUPAC code became a real base; `_pycodec` raised. Only `N`/`n` is substitutable now,
  in both, and every other non-ACGT byte raises with the character and offset.
- **`random` was not reproducible.** `_pycodec` used the unseeded global `random`.
- **Unstranded strand normalization was not reproducible either** — the same unseeded
  global chose which mate to reverse-complement, and `_accel`'s per-block xorshift made
  the choice depend on `--block-size`:

  ```
  one block of 16 : [13, 6, 5, 14, 5, 14, 13, 6, 5, 14, 13, 6, 5, 14, 13, 6]
  two blocks of 8 : [13, 6, 5, 14, 5, 14, 13, 6, 13, 6, 5, 14, 5, 14, 13, 6]
  ```

  Both random decisions are now derived from `(--seed, global record index, offset)`
  through a shared splitmix64 finalizer written identically in C++ and Python — stateless,
  so batching cannot shift a draw. Verified reproducible across runs, identical across
  backends, and unchanged across `--block-size`.

  **This changes every existing unstranded-normalized file**: the orientation assignment
  is different under the new scheme. Deliberate.

- Also fixed in passing: `_accel` measured a record's length in UTF-8 *bytes*, so a
  non-ASCII character corrupted the lengths column (`ACGéT` stored as length 6, decoding
  to `ACGAAT`). Both backends now raise.

### Changed — merge consensus (affects the corpus)

- **The overlap consensus is now written only into records whose construction depends on
  the overlap being real**: R1 on a merge, both mates on a trim, **neither on a kept
  pair**. It previously went into R1 on every detected overlap, including kept pairs —
  where both mates are emitted in full, so the corpus carried one corrected copy of the
  region beside one uncorrected copy of it.

  The reason given for that asymmetry was circular ("a kept pair emits R2 untouched",
  offered as the reason for emitting R2 untouched), and measurement says the write was
  net harmful. A detection that lands in *kept* is spurious almost by construction: at
  `shift >= 0` the pair is there only because the trim guard refused it, which requires
  an inferred overlap of more than about 110 bases — and a genuine overlap that long
  scores ~218 bits and would have merged at 28. Over 1M ground-truth pairs, of the 3,068
  kept pairs carrying a detected overlap **none found the true shift** and 97.3% had no
  true overlap at all. Wrong emitted bases in the overlap window: **208** correcting
  neither, 1,509 correcting R1 only, 17,870 correcting both. The old R1-only write turned
  1,379 correct bases wrong in order to fix 78.

  Corpus impact: 1,934 of 1,416,630 emitted records differ (0.137%) — 712 sequence lines
  and 1,934 quality lines. Merge decisions are untouched: the merged / trimmed / kept
  counts are identical.

  A claim that had propagated into three places is retracted with it: `shift >= 0` does
  **not** keep the write clear of R2's 5' end. The real bound is
  `min written R2 index = max(0, len2 - len1 + shift)`, which reaches index 0 at
  non-negative shifts once the mates differ in length — 9 of the recorded 237 corrupted
  5' ends are at `shift >= 0`, not read-through at all. What keeps the trim branch clear
  of that boundary is `trim_is_allowed`, which forces `shift >= min_read_length`.

- **An `N` is now rescued from the mate by design rather than by luck.** An `N` is a
  no-call: it carries no base information, so a real call on the other mate beats it
  whatever the two quality scores say, and the rescued base keeps the surviving mate's
  own quality rather than a contested-base derating — there was no contest. Previously
  the ordinary posterior decided, which rescued an N only because an instrument usually
  assigns one a low quality; a **high-quality N beat a real base and survived into the
  corpus**, where nothing downstream can distinguish it from a genuine ambiguity. Only
  `N` is rescued — the IUPAC codes carry partial information and are left alone — and the
  scan is untouched, so merge decisions do not move.

  Worth recording: the first implementation of this was wrong in the compiled backend
  only, and the cross-backend equality tests passed anyway, because no fixture in the
  suite had an `N` inside an overlap. There is one now.

### Fixed (silent)

Each of these produced a plausible-looking output and a zero exit status. Each is
now covered by a test that fails if the defect is reintroduced.

- **Heap buffer overflow in the compiled merge kernel.** `Scratch::name` was sized
  from the *read-length* arena (minimum 1024 bytes) but written from the *header*:
  a merged record's name is R1's header plus a `merged_<n1>_<n2>` tag. Any FASTQ
  whose headers outran its reads wrote past the end of the block. Reproduced with
  51 bp reads and a 16 KB header: one run aborted under malloc's heap check, an
  identical rerun returned well-formed output — meaning the quiet runs were
  corrupting adjacent heap. The reference backend was never affected, so no
  cross-backend test could catch it without a long-header fixture; there is one
  now. `name` has its own `ensure_name()`, sized from the header.

- **`--shuffle` and ZNA → ZNA re-encode silently dropped flag bits.** Every copy
  path routed the flags byte through `(has_start, has_end)`, which has three
  states where `(IS_RC, IS_FULL_FRAGMENT)` has four: `ENDS_BY_FLAG` maps byte 16
  and byte 24 alike to `(True, True)` — correctly, since both records do have two
  real fragment ends — so no inverse exists, and the one that existed cleared
  `IS_RC` on every full-fragment record. `--shuffle` is in the documented
  merged-read recipe and passed through the conversion twice. Measured on a
  200-record normalized file: `IS_RC` on 101 records before a shuffle, **0**
  after, at which point `--restore-strand` returned half the corpus in the wrong
  orientation. Fragment *geometry* was never affected, so this cost reversibility
  rather than supervision. See **Added** for the copy API that replaces it.

- **Re-encoding a non-normalized ZNA dropped `IS_FULL_FRAGMENT` entirely.** The
  same seam, the other branch: `preserve_normalization` was `is_reencoding and
  input_header.strand_normalized`, so a file that had never been normalized took
  the plain 4-tuple path, where the flag has no source at all and was re-derived
  as `False`. A merged read encoded with `--treat-unpaired-as-merged` read back as
  ends `(True, True)`; re-encoding that file gave `(True, False)`. ZNA → ZNA is
  now a copy in both branches.

- **Two independently merged reads sharing a name were encoded as one fragment.**
  `zna encode --interleaved` infers pairing from read names, so two merged
  molecules that happen to share one were paired: both lost `IS_FULL_FRAGMENT`,
  both gained `IS_PAIRED`/`IS_READ1`/`IS_READ2`, and two different molecules were
  written as each other's mate. Duplicate read names are ordinary in real data
  (concatenated lanes, some SRA dumps, UMI-collapsed files). A merged read is now
  never a mate: `_read_key` reports it from the `merged_<n1>_<n2>` token that
  `zna merge`, fastp and khorana's `parse_merged_fastq` already agree on, and
  pairing skips it. Nothing errors — the records simply come out right.

- **Two-file input with unequal record counts truncated the library.** The loop
  was `zip(p1, p2)`, which stops at the shorter file: every remaining R1 read was
  dropped, at exit 0. (`zip` also made it unfixable from outside — it has already
  pulled and discarded the extra value by the time it stops, so the off-by-one
  case looks like a clean finish.) Both files must now hold the same number of
  records, and the error names the one that ended early. Pairing remains
  positional; use `zna merge` or `--interleaved` if you need names compared.

- **A `.zna` passed alongside any other input wrote an empty file at exit 0.**
  The re-encode check is `len(files) == 1`, so with two inputs the binary went to
  the FASTQ parser, which found no records. It warned only about *format
  inference*. Now refused.

### Fixed

- **The two label extractors disagreed on malformed tag values**, and the
  compiled one produced numbers rather than errors: `NM:i:abc` → `5451`
  (`'a'-'0'` = 49, `'b'-'0'` = 50, …), `NM:i:3.7` → `287`, `NM:i:` → `0`, and a
  value past `int64` was signed overflow — undefined behaviour, not merely a
  wrong answer. The float path used bare `strtod`, so `3.7x` became `3.7` and an
  empty value became `0.0`. A label column is not self-describing, so nothing
  downstream could tell a garbage value from a real one.

  The fast path now handles plain decimal digits and **everything else is
  delegated to CPython's own `int()`/`float()`**, so both the value and the
  exception match the reference by construction rather than by a reimplementation
  kept in step by hand. Arbitrary-precision integers now survive exactly.

  **This is a behaviour change:** an encode that previously produced a corpus with
  garbage label values now fails with `ValueError`. That is the intent, but a
  pipeline relying on the old leniency will stop.

- `scripts/wheel_smoketest.py` never imported `zna.merge._accel`. It asserted the
  codec extension only, so it caught the merge target overwriting the codec but
  not the reverse, nor the merge target failing to build on a platform — in which
  case the wheel shipped with `zna merge` refusing to run. It now checks both
  extensions are present and are the right ones.

- `pytest.importorskip` for both extensions now passes `exc_type=ImportError`, so
  a half-built environment (a stale `.so`, a missing runtime dependency) skips
  rather than erroring. Required from pytest 9.1.

### Added
- **A lossless copy path: `ZnaReader.copy_records()` → `ZnaWriter.write_copy()`.**

  The rule it establishes, and which the code now states: **a view is for reading;
  the flag byte is for copying.** `records()` returns *views* —
  `(is_paired, is_read1, is_read2)`, `is_rc`, `(has_start, has_end)` — each chosen
  for a consumer, and none able to carry the whole flag byte back to a writer.
  `(has_start, has_end)` in particular remains exactly right for its purpose
  (*which edges of this record are true fragment boundaries*, which is what a
  downstream model wants) and is unchanged; it simply must never travel back into
  a writer.

  ```python
  with open(src, "rb") as fin, open(dst, "wb") as fout:
      reader = ZnaReader(fin)
      with ZnaWriter(fout, reader.header, preserve_normalization=True) as w:
          for rec in reader.copy_records():
              w.write_copy(rec)
  ```

  `copy_records()` yields `ZnaRecord(seq, flags, labels)` — the raw
  `ZnaRecordFlags` byte, verbatim. Nothing on either side interprets it, so a copy
  carries every bit, including combinations this version does not interpret and
  bits a later format version may define. `ZnaRecord` is three fields and named,
  and no `records()` shape is three wide, so code that confuses the two fails at
  the unpack instead of reading `flags` as `is_paired`.

  `zna shuffle` and `zna encode` of a `.zna` both use it. Neither codec backend
  changed: both were already the identity on the flags column with normalization
  off, and OR in `IS_RC` with it on — which is why `write_copy` requires
  `preserve_normalization=True` and says so rather than corrupting silently.

- **`zna merge` — overlap-merge paired-end reads into one mixed interleaved
  FASTQ**, ready for `zna encode --interleaved`. Replaces fastp's PE-merge step.

  Every pair is scored **once**: R1 is slid against `revcomp(R2)` over the single
  axis of candidate fragment lengths, and each shift is scored as a
  log-likelihood ratio in **bits** against chance alignment — `+1.9855` per
  matching base (`log2 4`, the information in agreeing on one of four bases),
  `-6.2288` per mismatch at a 1% error rate. Both weights fall out of the error
  rate; neither is tuned. The `argmax` shift (not fastp's first-accept) is then
  read at two thresholds: at or above `--threshold-merge` (28 bits) the pair
  becomes one full-fragment record; between that and `--threshold-trim` (8 bits)
  both reads are kept and the redundant overlap is cut off R2's **3'** end; below
  it both are kept untouched. Three parameters in total, all with units, and read
  length is not one of them.

  Measured against the fastp-derived rule it replaces: spurious detection on
  genuinely unrelated pairs **5.17% → 0.26%**, zero spurious merges in 20,000
  uniform-random pairs at 28 bits, merge rate **85.7% → 88.8%**.

  Where the mates overlap and disagree, the consensus takes the better-supported
  call by posterior from the two Phred scores and derates its quality, because a
  contested base is less certain than an uncontested one. On 167,784 real pairs
  that is **35.2% fewer expected residual errors** than fastp's Q14/Q30 cutoffs,
  entirely in the quality band those cutoffs never reach.

  Why it belongs here: base 0 of every emitted read is a true fragment boundary
  (nothing is ever cut from a 5' end), a merged record *is* its fragment exactly,
  and unmerged pairs are emitted all-or-nothing so a lone mate never becomes a
  spurious "full molecule with both endpoints". Those are precisely the
  properties `IS_RC`, `IS_FULL_FRAGMENT` and `--treat-unpaired-as-merged` rest
  on, and they were previously asserted across a repo boundary. See
  [docs/METHODS.md](docs/METHODS.md) for the derivation, and
  [docs/ROADMAP.md](docs/ROADMAP.md) for what was measured and rejected.

  Also available in-process as `zna.merge.process_pair` / `zna.merge.find_overlap`,
  and as `python -m zna.merge`.

- **Per-record provenance — what happened to a *read*, not just to a library.**

  `zna merge` appends tokens to each emitted record's header, on all three outcomes.
  Existing header fields are always passed through untouched — provenance is *appended*,
  never substituted — so `--label` reads the same tags off an emitted record that it
  would have read off the input:

  ```
  @SRR1.7  ZI:i:42 ZN:i:6 trim3_12 rescued_1 merged_90_0
  ```

  `trim3_<n>` / `subn_<n>` are bases the N policy removed or substituted, `rescued_<n>`
  is no-calls this record recovered from its mate. A record nothing happened to is
  emitted byte-unchanged, so a clean library pays nothing: 1,332,353 of 1,418,525
  records on the 1M benchmark at a 1.5% no-call rate.

  **`ZN:i:<bits>` is the half that survives encoding.** ZNA stores no headers, so it is
  the only per-record provenance that reaches a corpus, and it gets there as an ordinary
  label column — declare it and you get one byte per record, omit it and nothing changes:

  ```bash
  zna encode --interleaved --treat-unpaired-as-merged \
             --label provenance:C:ZN -o reads.zna merged.fq.gz
  ```

  Bits: **trimmed 1, rescued 2, N-trimmed 4, N-substituted 8.** There is deliberately no
  "merged" bit — that fact already has two homes (the `merged_` token, and
  `IS_FULL_FRAGMENT` in the flag byte), and spending a bit on it would have put
  ` ZN:i:1` on ~82% of emitted records to repeat them. Every bit above is one with
  nowhere else to live; `trimmed` above all, since a trimmed pair is emitted as an
  ordinary pair and no ZNA flag distinguishes it from one kept whole. Note
  `IS_FULL_FRAGMENT` is set only under `--treat-unpaired-as-merged`, so an encode that
  omits that flag records neither.

### Performance
- **A second compiled extension, `zna.merge._accel`**, carries the whole per-pair
  path: FASTQ parsing, the overlap scan, the posterior consensus, record
  construction, formatting and the histograms. It releases the GIL, so
  `--threads` are real worker threads.

  | | µs/pair |
  |---|---:|
  | `--threads 1` | 2.78 |
  | **`--threads 2`** | **1.40** |
  | `--threads 4` | 1.43 |

  The scan alone went from 2.633 µs/pair (the numba kernel this replaces) to
  0.470 in C++ — 5.6x — measured on 50,000 real pairs with full pruned scans.

  The tool is now I/O bound, and measurably so: with gzip removed from both ends
  it runs at **0.42 µs/pair**, while decompression costs 1.42 and compression
  0.52. Further speed needs a faster inflate or removing the FASTQ intermediate
  altogether, not a faster kernel.

- **The kernel compares raw bytes, 16 at a time**, through one vector primitive:
  NEON on aarch64, SSE2 on x86-64. Both are baseline, so there is no `-march`
  flag, no runtime ISA dispatch, and a scalar fallback everywhere else. Byte
  comparison is not a concession to portability — it *is* the reference
  semantics, for ACGT and equally for N, IUPAC codes and lowercase, so no fast
  path can disagree with the reference on any input. A 2-bit packed kernel
  measured slower (0.535 vs 0.470) *and* would have needed a purity dispatch to
  keep those semantics.

### Verified against ground truth
- **The full production recipe was re-run end to end on the 1M-pair library**, at
  `zna merge` → `zna encode --interleaved --treat-unpaired-as-merged
  --strand-normalize --shuffle`, and the flag column checked record by record:

  ```
  1,416,630 records
    merged, not flipped      byte 16   288,509  ┐
    merged, reverse-compl.   byte 24   293,853  ┘ 582,362 = every merged read
    mates                  bytes 5,6,13,14      834,268
    fragment geometry     (True,True)  582,362   (True,False)/(False,True) 417,134 each
  ```

  **The 293,853 records at byte 24 are the ones this release stops losing** — every
  one of them previously came back as byte 16 with `IS_RC` cleared. Decoding the
  shuffled, normalized file with `--restore-strand` now reproduces the original
  base multiset **exactly**, all 1,416,630 records.

- **`zna merge` was benchmarked against fastp on 1,000,000 simulated pairs from
  hg38**, with the true fragment length, true overlap and every injected error
  known exactly. Everything measured before this proved the tool matches
  *itself*; this is the first thing that checks whether its answers are true.

  **Base 0 of every emitted read is a true fragment boundary — 0 violations in
  1,416,630 records**, and no orphaned mates. That is the contract `IS_RC`,
  `IS_FULL_FRAGMENT` and `--treat-unpaired-as-merged` rest on, and it had only
  ever been checked against the tool's own inferences.

  Merge sensitivity steps exactly where the derivation puts it (92.6% at a true
  overlap of 15–19, 100% at 30+), read-through fragments reconstruct at 100.00%,
  merged records equal the true fragment 86.6% of the time against fastp's
  85.9%, and the quality-aware consensus recovers **90.4%** of the recoverable
  overlap errors against fastp's 74.1% — with a constant quality string both
  recover 0.0%, which is why that comparison needs a realistic quality model.
  Throughput 1.53 µs/pair against fastp's 2.83 on the same input.

  0.96% of merges are at the wrong shift, all of them fragments whose two ends
  genuinely repeat: the scan never picked a lower-scoring alignment than the true
  one (0 of 848 checked), 56% carry no sequencing error, and the hotspots are
  pericentromeric satellite. See
  [docs/MERGE_BENCHMARK_RESULTS.md](docs/MERGE_BENCHMARK_RESULTS.md) for the
  analysis and `scripts/merge_bench/` for the simulator and scorer.

  **The trim band was measured too**, because a pair that does not merge is still
  encoded and its redundant overlap reaches the corpus just as directly. Scored
  against truth over every kept pair: the overlap comes off exactly on 96.0% of
  pairs with a true overlap of 5–9 bases and 99.8% at 10–14, R1's 3' end is never
  moved, and **2.25 duplicated bases survive per 10,000 unmerged bases against
  85.79 with no trim at all**. Sweeping `--threshold-trim` over 8/12/16/20 with
  the merge decision held fixed shows the default minimises duplicated-plus-
  deleted bases — the first time that parameter has been chosen by measurement
  rather than by argument.

### Changed
- **`ZnaWriter.write_records()` no longer accepts the 6-tuples of
  `records(with_ends=True)`; it raises, naming `copy_records()`.** Its docstring
  advertised that shape as "a lossless ZNA → ZNA copy in one line", which was
  false — see the flag-bit entry above. 4-tuples and the 5-tuples of
  `records(with_rc=True)` are unchanged.

  Removed with it: the private `_flags_from_ends`, `_RC_FULL_BY_ENDS` and
  `_FLAG_BITS_BY_ENDS`. `_flags_from_ends`' docstring claimed "a bijection over
  the three reachable states"; there are four, and every caller was a corruption
  site, so there was nothing to fix. The public read-side tables — `ENDS_BY_FLAG`,
  `FLAG_FIELDS` — are correct and unchanged, and `ENDS_BY_FLAG` is now documented
  as the deliberate many-to-one projection it is.

- **The trim is now symmetric: the redundant overlap is split between the two
  mates' 3' ends** instead of being taken entirely off R2. The overlap sits at the
  3' end of *both* reads — each starts at a fragment end and reads inward — so
  splitting it tiles the fragment exactly once just as before, while leaving the
  two emitted reads the same length, which is what downstream aligners and models
  expect, and discarding the last cycles of both reads rather than one read's whole
  copy. Where the mates disagree, both now carry the consensus call.

  Measured on 1,000,000 simulated pairs, **every decision and every accuracy metric
  is identical** to the old rule — same merges, same trims, same duplicated and
  deleted bases, zero boundary violations — while the mean `|len(R1) − len(R2)|` on
  trimmed pairs falls from the whole overlap (10.8 bases, up to 110) to **0.51**,
  at most 1. Throughput is unchanged.

  **This changes the output bytes for trimmed pairs**, so a corpus re-encoded with
  0.4.0 will not be byte-identical to one built with the pre-release tool. Merged
  and untouched pairs are unaffected.

  Two subtleties, both found by re-running the ground-truth benchmark and both
  pinned by tests. The consensus is written into R2 **only** on the trim branch: on
  a read-through geometry the overlap lands on R2's *5'* end, and rewriting bases
  there on an inference too weak to merge on corrupted 237 fragment boundaries per
  million pairs. And the trim guard, which read as "don't leave R2 under
  `--min-read-length`", was also capping the overlap the trim band may act on; it
  is now stated symmetrically as *each mate must reach at least `--min-read-length`
  past the other's 3' end*, which recovers both properties.
- **`--threshold-trim`'s justification now carries a number.** Both `args.py` and
  the redesign document said a wrong trim "costs a few bases"; measured against
  ground truth it costs a median of 9, a mean of 20 and up to 110. The asymmetry
  the claim was making still holds — the band removes 4.45 duplicated bases for
  every real base it deletes — but it is stated with the measurement now. No
  behaviour change; the default stays 8.
- **The three `zna merge` histograms are no longer capped at 1024 bins.**
  `length_histogram`, `overlap_length_histogram` and `insert_size_histogram`
  clamped every index to bin 1024, so reads past that length silently aggregated
  — while read length itself is uncapped. They now size themselves from the
  scratch arena (`2 × cap + 1` bins), so indexing stays a single dense lookup in
  the per-record path and there is no cap at any read length. Histogram bins with
  no counts were already omitted from the JSON, so nothing downstream changes.
- **`insert_size_censoring.floor` was removed** from the `zna merge` JSON. It was
  a second copy of `params.min_read_length`; `min_mergeable_overlap`, the bound
  the JSON did not otherwise carry, stays. No consumer read it.

### Reproducibility
- **The score is computed in integers**, at 2²⁴ units per bit. No float takes
  part in a comparison, a bail bound or the argmax, so a given FASTQ produces
  byte-identical output on every platform and compiler. The scale is chosen by
  exhaustive enumeration rather than taste: 2²⁰ first disagrees with float
  arithmetic at an overlap of 2,575 bases, 2²⁴ at 32,830.
- **The argmax is a specified total order** — maximise score, then minimise the
  shift — rather than an artifact of iteration order, verified against an
  unpruned exhaustive scan including deliberately constructed ties.
- **Output is written in submission order**, so it is byte-identical at any
  `--threads`.

### Packaging
- **No new runtime dependencies.** The optional `numba` dependency that the
  merge tool used before it was ported here is gone; the reference kernel is
  plain Python and the fast one is compiled with the rest of the package.
- `zna merge` **refuses to run on the reference kernel** unless asked by name
  with `--backend python`. It is correct and ~50x slower, and a silently-correct
  50x slowdown at cluster scale is indistinguishable from a slow node.
- **The Windows wheel did not build.** `neq16`'s SSE2 path counted equal lanes with
  `__builtin_popcount`, a GCC/Clang extension MSVC does not have, so the first Windows
  build of the merge extension failed with `C3861`. It survived to a tag for two
  compounding reasons: `zna merge` is new in 0.4.0, so no MSVC had ever compiled
  `merge_core.hpp`; and an arm64 developer machine takes the NEON path, where the
  popcount does not appear at all.

  Replaced with a portable SWAR fold. MSVC's `__popcnt16` is deliberately *not* used —
  it emits the POPCNT instruction, which is not baseline x86-64, and would have turned a
  build error into an illegal-instruction fault on an older CPU. The fold is compiled on
  every platform even where nothing calls it, and `tests/test_merge.py::TestPopcount`
  checks it over all 65,536 inputs, so the branch is no longer testable only on the one
  build that cannot run the suite.

  `merge_core.hpp` also now includes `<vector>` itself, which it uses and had been
  getting from whichever translation unit included it first.
- **CI now runs on every push, not only on `v*` tags.** The only workflow was
  `publish.yml`, so nothing compiled this project on Windows or Linux until a release was
  already tagged — which is the whole reason the above cost a tag. `.github/workflows/ci.yml`
  builds and tests on Linux, macOS and Windows across the ends of the supported Python
  range, asserts *both* extensions actually loaded, and does not fail-fast, so one
  platform's compiler error cannot hide another's.
- **`scripts/release.sh` had two ways to produce a wrong release**, both found by
  reading it rather than by running it:
  - It pushed a hardcoded `main` and then tagged `HEAD`. Run from a feature
    branch — which is where the work happens — that pushes whatever local `main`
    points at and then tags a commit `main` does not contain. It now refuses
    unless you are on `main`, and says how to get there.
  - Its version bump is a no-op when `src/zna/__init__.py` is already at the
    target, which is the normal case; `git commit` with nothing staged exits
    non-zero, and under `set -e` that aborted the release *after* the
    confirmation prompt and before the tag. Nothing to commit is now a normal
    state.

  `docs/RELEASING.md` gains the merge-to-`main` step, the two-environment
  verification, and a post-publish wheel check for the two-extension collision.

## [0.3.5] - 2026-08-12

Two fixes found by reviewing the pipelines that write and read ZNA files
(hulkrna and khorana). No on-disk format change.

### Fixed
- **`zna encode --shuffle-buffer-size` did nothing.** The flag was declared and
  parsed, but `encode_command` passed a hardcoded `1 << 30` to `shuffle_zna`, so
  the bucket budget was always 1 GiB no matter what was asked for.
  (`zna shuffle`, the standalone command, honoured its `--buffer-size`
  correctly; only the encode path was affected.) The default still parses to
  1 GiB, so no existing invocation changes behaviour — the flag simply works now.

### Added
- **`ZnaReader.blocks(labels=...)`** — `blocks()` previously refused labeled
  files outright, on the reasoning that handing back sequences while silently
  dropping the label columns was worse than not offering the API. That reasoning
  still holds for the *default*, but it made the batch API unusable on any
  labeled corpus, which is what the pipeline it was built for actually produces.

  The default still raises, now naming the way out. `labels=False` skips the
  columns — not silent, the caller asked — and `labels=True` yields a third
  element holding one value-tuple per label column, in header order.
  `len(label_columns)` always equals `len(header.labels)`, so an unlabeled file
  yields `()` rather than erroring.

  Measured on a three-column file (int32, int32, uint8) of 200k x 150 bp:
  `records()` 32.7 ms, `blocks(labels=True)` 15.5 ms (2.12x),
  `blocks(labels=False)` 9.5 ms (3.44x).

  Note this does not avoid *decompressing* the label bytes — a block payload is
  a single zstd frame holding all four columns — it avoids unpacking them into
  Python objects, which is where the cost is.

### Packaging
- **CPython 3.14 wheels are now built.** `CIBW_BUILD` had listed `cp314-*` since
  0.3.3, but the workflow pinned `pypa/cibuildwheel@v2.23`, which predates 3.14
  entirely — the entry matched nothing and was silently dropped, so 0.3.3 and
  0.3.4 shipped no cp314 wheels despite asking for them. Verified: cibuildwheel
  2.23.3 selects cp310-cp313 for this project, 4.2.0 selects cp310-cp314.
- **The manylinux baseline is pinned to `manylinux2014`.** cibuildwheel 3.0
  changed the default image to `manylinux_2_28`, which would have raised the
  glibc floor from 2.17 to 2.28 and dropped RHEL/CentOS 7-era HPC clusters as a
  side effect of adding 3.14. `manylinux2014` still builds 3.14, and `zstandard`
  — zna's own runtime dependency — publishes its cp314 wheel as
  `manylinux_2_17`, so taking the new default would have made zna harder to
  install than the package it depends on.
- **Wheels are now smoke-tested before publication.** They were previously
  uploaded without ever being imported. `scripts/wheel_smoketest.py` runs as
  `CIBW_TEST_COMMAND` on every wheel: it asserts the compiled extension actually
  loaded (a wheel that silently falls back to the pure-Python backend is a
  broken wheel) and exercises `records()`, `blocks()` with and without labels,
  `block_index()`, `restore_strand`, and the N-policy error path.

### Fixed (latent)
- `blocks()` split the block payload as `flags | lengths | sequences`, which is
  only correct without labels. It never ran on a labeled file because of the
  guard above, but lifting that guard without fixing the offset would have
  decoded sequences from the wrong position — plausible garbage, no error. The
  split now accounts for the label columns whether or not their values are
  wanted, and the fuzz harness checks it against `records()` across dtypes,
  compression, sharding, `indices=`, and `restore_strand`.

## [0.3.4] - 2026-08-11

Performance release. Nothing about the on-disk format changes, and no existing
file needs re-encoding.

### Fixed
- **The pure-Python backend silently corrupted sequences containing an invalid
  base at an index congruent to 3 mod 4.** `_pycodec.encode_sequence` validated
  each 4-base group with `val > 255`, applied *after* OR-ing the four 2-bit
  values together. A 255 (invalid) in the group's last slot ORs to exactly 255
  and never exceeds it, so instead of raising, the group was packed as `0xFF`
  and decoded back as `TTTT` — losing the three valid bases beside it. This was
  reachable on real data: the IUPAC ambiguity codes (`R`, `Y`, `S`, `W`, ...)
  occur in FASTQ, and with an N-policy set they hit this path rather than being
  substituted. The C++ backend was never affected. Validation now happens once
  up front in a single scan, which is also faster than the per-group test.

  Found by the new round-trip fuzz harness on its first run.

### Added
- **`tests/test_fuzz_roundtrip.py`** — a round-trip fuzz harness asserting
  bit-exact recovery across sequence content x length (including 0, 1, and every
  residue mod 4) x `seq_len_bytes` x N-policy x strand configuration x
  compression x record layout x backend, against a reference model written from
  the specification rather than from `_pycodec`. It also checks that the two
  backends produce byte-identical columnar streams, that re-encoding is a
  faithful copy, that every label dtype round-trips, and that the C++ decoder
  neither leaks nor over-frees.

  This exists because the suite could not previously catch codec corruption:
  it round-trips fixed inputs and checks API contracts, and two changes proposed
  during the 0.3.4 sweep silently corrupted data with all 282 tests passing.

- **`ZnaReader.blocks(stride=1, offset=0, restore_strand=False)`** — a columnar
  batch API yielding `(list_of_sequences, flags_bytes)` per block, for consumers
  that process a batch at a time. `stride`/`offset` shard by block and seek past
  blocks the shard does not want, so a sharded reader decodes its share instead
  of the whole file: 1.8x for one worker of 2, 9.4x for one worker of 16.
  Requires arbitrary record order (use `zna shuffle`) and many more blocks than
  shards; warns when a shard matches no blocks. Raises on labeled files.
- **`ZnaReader.block_index()`** — one `BlockInfo(index, offset, n_records,
  comp_size, uncomp_size)` per block, built by walking the block-header chain
  and seeking over each payload. The *file* header stores no record or block
  count; every *block* header carries its own count, so the totals cost no
  decompression: 2.3 us per block, about 1.4 ms for a 38 MB / 611-block file,
  against 89 ms to reach the same counts by decoding. Use it to size a
  subsample, or to catalogue a corpus.
- **`blocks(indices=...)`** — decode an explicit set of block numbers instead of
  a stride. Needed when the sampling fraction is not a unit fraction, and when
  repeated passes over one file should see different blocks: `stride` admits
  only `stride` distinct phases, so several epochs at `stride=4` would revisit
  the same four subsets. Mutually exclusive with `stride`/`offset`.
- **`zna inspect --json`** (with `--blocks` and `--counts`) — machine-readable
  output including `n_blocks` and `n_records`, for building manifests that carry
  record counts and can therefore weight a balanced corpus sample without
  opening the data files.
- **`FLAG_FIELDS`** — `flag byte -> (is_paired, is_read1, is_read2)`, so
  `blocks()` consumers need not re-derive the bit layout.
- **`ENDS_BY_FLAG`** — `flag byte -> (has_start, has_end)`, now public. This is
  the only correct way to recover fragment geometry from the raw flags column:
  under unstranded normalization ZNA reverse-complements *one mate per pair, at
  random*, so which edge is a true fragment boundary is a per-record fact
  carried by `IS_RC`, not a property of R1 versus R2. A consumer that assigns
  endpoints by mate number instead is right about half the time.
- `decode_block_sequences` on both codec backends.
- `ZNA_NO_EXTERNAL_GZIP=1` forces gzip input through the Python `gzip` module.

### Changed
- **The C++ codec no longer copies every base twice around the work that
  matters.** nanobind was materialising a `std::vector<std::string>` on the way
  in and a `std::vector<std::tuple<...>>` on the way out. Decode now writes bases
  straight into each `str`'s own buffer, four at a time from a lookup table, and
  returns a hand-built list of GC-untracked tuples; encode reads each sequence's
  buffer in place and folds reverse-complement into the packing loop, so a
  normalized record needs no intermediate string. The decoder also emits the
  caller's tuple width and can undo strand normalization itself, which removes
  the Python-side rebuild of every record. `encode_block` and
  `encode_block_labeled` were hand-copied duplicates and now share one core.

  Verified byte-identical to 0.3.3 across 90 encode/decode configurations,
  including `npolicy="random"` and every reverse-complement path.

- Gzip FASTQ input is decompressed by `pigz`/`gzip` when available, so it runs
  in its own process rather than contending for the GIL. Falls back to the
  `gzip` module. A truncated or corrupt `.gz` is now an error rather than a
  silently short input — previously the partial stream was indistinguishable
  from a small file, and `zna encode` would write a truncated ZNA and report
  success.
- `shuffle`'s counting pass reads flag columns instead of decoding every
  sequence, and its bucket files are written at zstd level 1 rather than the
  source file's level.
- `CMakeLists.txt` no longer passes `STABLE_ABI`. It was silently a no-op —
  nanobind honours it only when CMake finds `Development.SABIModule`, and the
  `find_package` asks for `Development.Module` — and the module could never have
  been limited-API anyway.

### Performance

Measured interleaved, min-of-N, in-process (see `scripts/bench_ab.py`); the C++
figures are cross-process, since the two arms are different builds.

| path | change |
|---|---|
| `shuffle` counting pass | +2899% |
| raw `encode_block`, reverse-complementing | +163% |
| gzip FASTQ input | +123% |
| raw `encode_block` | +113% |
| raw `decode_block` | +109% |
| `records(restore_strand=True)` | +147% |
| `records()` | +97% |
| one sharded worker of 16, via `blocks()` | +837% |
| single-end encode loop | +52% |
| gzip module read-ahead | +51% |
| FASTQ sequence-line trimming | +31% |
| `reverse_complement` | +31% |
| interleaved pairing, bytes-native | +29% |
| N-drop filter | +26% |
| decode output f-strings | +19% |
| `_flags_from_ends` in a write loop | +11% |

## [0.3.3] - 2026-08-11

### Fixed
- **Labeled encoding no longer discards pairing.** `stream_inputs_labeled` had
  its own single-file FASTQ reader and flagged every record unpaired, so
  `zna encode --interleaved --label-defs …` silently ignored `--interleaved` and
  wrote `is_paired=False` on every record — losing the R1/R2 flags and, because
  the encoder then drew an independent reverse-complement coin per record
  instead of one per pair, mate co-orientation as well. Labeled encoding now
  dispatches over the same input strategies as unlabeled (two paired files,
  interleaved, single-end) and shares one implementation of the interleaved
  pairing rule. A consumer cannot otherwise distinguish a mate (one true
  fragment boundary) from a merged read (two), which is the difference between
  correct fragment-end supervision and endpoint tokens on interior positions.
  **`is_paired` is not recoverable from affected files; re-encode them from
  FASTQ.** Labeled FASTA input is now rejected explicitly rather than parsed as
  FASTQ — labels come from SAM tags in the read header.
- **Copying records between ZNA files no longer re-applies strand
  normalization.** `strand_normalized` in a header describes what was already
  done to the data, but the writer read it as an instruction to do it, so any
  read-then-write path normalized a second time. Orientation is not idempotent:
  applying it twice returns the data to an un-normalized state while the header
  still claims otherwise. Affected `zna encode` on a `.zna` input, `zna shuffle`,
  `zna encode --shuffle`, and `zna._shuffle.shuffle_zna`. Encoding from FASTQ was
  never affected, and no on-disk format or existing file's bytes change.
- **`zna shuffle` now preserves each record's `IS_RC` flag.** It writes through
  two `ZnaWriter`s (buckets, then output), so it applied two further random
  reverse-complements. Mate co-orientation survived by parity, leaving the
  sequences valid and the damage invisible, but `IS_RC` ended up uncorrelated
  with which mate had actually been flipped.
- **Stranded re-encodes are stable.** Under `--strand-specific` the second pass
  reverse-complemented read1 (and merged/single reads) back to their original
  orientation, silently un-normalizing a strand-specific file.
- `zna encode --label ... input.zna` no longer exits 0 having written an empty
  file. A ZNA input routed to the labeled encode path was parsed as FASTQ,
  matched no records, and produced a valid, correctly-labeled, empty output.
- **Affected files, and what to do with them.**
  - A file re-encoded an **odd** number of times has lost mate co-orientation.
    Its stored `IS_RC` records exactly the last pass's reverse-complements, so
    `zna decode --restore-strand` undoes that pass and hands back the correctly
    co-oriented normalized sequences — byte-identical to the file before the
    damage, verified on both backends. What is *not* recoverable from the file
    is the `IS_RC` provenance and the original as-sequenced orientation, since
    the first pass's flags were overwritten.
  - A shuffled file is still correctly co-oriented and remains usable exactly as
    before, but its `IS_RC` flags no longer identify which mate was
    reverse-complemented, so it should not be read with `with_rc` below. Two
    passes leave each pair with either both mates flipped or neither, so
    re-normalizing it cannot recover the geometry.
  - In both cases, re-derive from the original FASTQ if the fragment-boundary
    geometry is needed.

- **`--npolicy drop` is now pair-atomic.** It filtered records *after* pairing,
  so an N in one mate dropped that mate and left its partner behind as a lone
  `IS_PAIRED` record — indistinguishable downstream from a genuine single-end
  read, and the only route to the encoder's mate-detection window mis-pairing
  two unrelated R1s. The whole fragment is now dropped together.

### Added
- **`IS_FULL_FRAGMENT` record flag (bit `0x10`)** — the record spans its entire
  fragment, so **both** of its edges are true fragment boundaries. `IS_RC` can
  only name one edge, which under-marks merged reads and over-marks genuine
  single-end reads, and nothing in the format previously distinguished them.
  Full-overlap **pairs are detected automatically** (mates spanning the same
  interval are exact reverse complements); unpaired records are declared with
  `zna encode --treat-unpaired-as-merged`, which defaults off — the conservative
  reading. Bits 4–7 of the record flag byte were unused and existing readers
  extract specific bits, so this needs **no format version bump**: old readers
  read new files unchanged, and new readers see the bit clear on old files.
- `ZnaReader.records(with_ends=True)` yields `has_start, has_end` — whether the
  left and right edge of the stored sequence is a true fragment boundary. This
  is the form a consumer placing fragment-end supervision wants: it combines
  `IS_RC` with `IS_FULL_FRAGMENT` so the geometry never has to be re-derived.
  The pair round-trips losslessly back to both flags, so it is also the shape
  for a faithful copy — `write_records` accepts it directly under
  `preserve_normalization=True`.
- `zna inspect --counts` reports a full-fragment tally and a **(mate × IS_RC)
  cross-tabulation**, with the expected pattern for the file's strandedness and
  a warning when R1 and R2 counts differ (orphaned mates). A bare RC total is
  equally consistent with a healthy file and a broken one; this is the tally
  that actually verifies a file's geometry before anything trains on it.
- `zna decode` warns when the input is strand-normalized and `--restore-strand`
  was not given: the output is then in normalized orientation, and re-encoding
  it with `--strand-normalize` applies orientation a second time — silently,
  since both files' headers are identical.
- `ZnaReader.records(with_rc=True)` yields the per-record `IS_RC` flag:
  `(seq, is_paired, is_read1, is_read2, is_rc)`, or
  `(seq, is_paired, is_read1, is_read2, is_rc, labels)` for labeled files. This
  is the only way to recover which edge of a mate is a real fragment boundary —
  it cannot be derived from the sequence, and `restore_strand=True` is not a
  substitute because it *consumes* the flag to undo the reverse-complement. The
  two are mutually exclusive and raise `ValueError` together. Opt-in: the
  default tuple widths (4, and 5 labeled) are unchanged and now pinned by test.
- `ZnaWriter(..., preserve_normalization=True)` with an optional trailing
  `is_rc` on `write_record`, plus `write_records` acceptance of the 5-tuple
  `records(with_rc=True)` yields, making a lossless ZNA → ZNA copy expressible
  for the first time. `write_zna` takes the same flag. Supplying `is_rc` without
  it raises `ValueError`.

### Changed
- `zna encode` refuses, rather than silently mangling, two re-encode requests it
  cannot honour: changing `--strand-specific` on an already-normalized input
  (orientation cannot be re-derived after the fact), and re-encoding a labeled
  ZNA file (label columns are not carried through that path).

### Documentation
- README gains "Unstranded Normalization and Fragment Geometry", stating the
  invariant the `IS_RC` flag exists to record: whichever mate was
  reverse-complemented ends up at the right of the common frame, so its right
  edge is the real fragment boundary and its left edge is a read-length cutoff;
  for the other mate it is the mirror image. The record format spec now lists
  `IS_RC` (bit 3), which it had omitted.

### Internal
- 20 tests covering the fragment-boundary contract, re-encode and shuffle
  record-identity, the default tuple widths, the stranded mirror case, and
  cross-backend `IS_RC` equality (both decoders, labeled and unlabeled). The
  pass-through encode path is deterministic, so it is additionally asserted
  byte-identical between the Python and C++ backends — unlike the deriving
  unstranded path, whose coin comes from a different PRNG in each.

## [0.3.1] - 2026-07-13

### Fixed
- Strand normalization now treats single-end / merged reads as **read1** under
  `--strand-specific`. Previously such reads were never reverse-complemented, leaving
  merged reads (e.g. from fastp's mixed interleaved output) on the opposite strand from
  the normalized paired read1 records. Decoding and already-encoded files are unaffected;
  only new strand-specific encodes of single reads change.

### Added
- `zna inspect --counts` reports per-flag record counts (paired R1, paired R2,
  single/merged, and reverse-complemented). Useful for validating mixed
  paired-end + single-end interleaved streams. The default header-only scan is
  unchanged; `--counts` reads block payloads (decoding only the flags column).
- Encode-time warning when `--strand-normalize` is combined with
  `--strand-specific` but neither read1 nor read2 is antisense — a configuration
  in which no read would ever be reverse-complemented (suppressed by `-q`).

### Internal
- Added a cross-backend lockstep test suite asserting the Python (`_pycodec`)
  and C++ (`_accel`) encoders produce byte-identical output across a battery of
  flag / strand / N-policy / length combinations, guarding against drift in the
  duplicated encode logic.

## [0.3.0] - 2026-03-25

### Added
- **Per-sequence label columns** for storing SAM-tag metadata (alignment scores,
  edit distances, etc.) alongside compressed sequences
  - Labels defined via `--label NAME:TYPE` on the CLI (repeatable)
  - 11 numeric dtypes: `A`, `c`/`C`, `s`/`S`, `i`/`I`, `f`, `d`, `q`/`Q`
  - Per-label `missing` value for reads where a tag is absent
  - Columnar storage: each label column packed contiguously per block
- **Decoupled label name and tag** — label `name` (stored in ZNA, up to 16
  chars) can differ from the SAM `tag` parsed at encode time
  - CLI 3-part format: `--label edit_dist:C:NM` (name:type:tag)
  - YAML: optional `tag` field per label definition
  - Tag is encode-time only; not stored in the ZNA file
  - If only name or tag is given, the other defaults to the same value
- **YAML label specification files** (`--label-defs labels.yaml`)
  - Define labels, descriptions, dtypes, and missing values in a single file
  - CLI `--label` / `--label-desc` flags override YAML definitions
  - Sample template: `examples/labels.yaml`
- **C++ accelerated label encode/decode** — labeled files now use the C++
  backend for both encoding and decoding, eliminating the Python fallback
  - `encode_block_labeled`: encodes sequences + pre-packed label columns in C++
  - `decode_block_labeled`: decodes sequences and unpacks label columns by dtype
  - `extract_labels_fast`: zero-allocation C++ SAM-tag extraction from headers
- **Optimized Python header parsing** (fallback path)
  - Eliminated `bytes.decode()` calls — `int()`/`float()` accept bytes directly
  - Replaced lambda dispatch with integer conv-codes
  - Early-exit when all labels found

### Performance
- **Labeled encode**: 336K rec/s (was 80K rec/s) — **4.2× faster**
- **Labeled decode**: 2.23M rec/s (was 610K rec/s) — **3.7× faster**
- Label overhead reduced from ~39% to within typical I/O variance

### Changed
- Binary format: label definition on-disk size is now 89 bytes
  (16 name + 64 desc + 1 dtype_code + 8 missing)
- `LabelDef.param` replaced by `LabelDef.missing` (explicit missing sentinel)
- Header field splitting uses any whitespace (tabs and spaces) so fastp
  `merged_XX_YY` suffixes no longer corrupt trailing tag values
- Label tag parsing now supports variable-length keys in `KEY:TYPE:VALUE`
  header fields (not limited to 2-character SAM tags)
  - Python fallback parser and C++ `extract_labels_fast` use dynamic key parsing
  - Supports mixed key lengths in the same header (e.g. `NM` + `edit_distance`)
- Dropped fixstr (`Z`) label type — all labels are numeric

### Documentation
- Expanded README strand guidance to explicitly document
  `--strand-normalize` behavior and its interaction with `--strand-specific`
- Updated label docs and `examples/labels.yaml` to show `name` + `tag`
  decoupling and custom non-SAM, variable-length header keys

## [0.2.0] - 2026-03-24

### Added
- **Strand normalization** (`--strand-normalize`) for consistent strand orientation in ZNA files
  - Unstranded libraries: randomly reverse-complements one read per pair (or SE reads) with per-record `IS_RC` flag
  - Strand-specific libraries: deterministically RCs antisense reads to sense orientation
  - Orthogonal to library metadata (`--strand-specific`, `--read1-antisense`, etc.)
  - `--restore-strand` on decode reverses RC operations using per-record `IS_RC` flags
- New `IS_RC` record flag (bit 3) tracks which records were reverse-complemented
- New `STRAND_NORMALIZED` header flag (bit 3) indicates strand normalization was applied

### Fixed
- C++ accelerator bug: paired R2 reads could receive an unintended second random RC during unstranded normalization

### Changed
- Bioconda recipe: added `osx-arm64` and `linux-aarch64` to `additional-platforms`

## [0.1.8] - 2026-02-11

### Fixed
- macOS conda build: removed `llvm-tools` from host dependencies

## [0.1.7] - 2026-02-06

- **Documentation changes**

## [0.1.6] - 2026-02-06

### Added
- **Built-in shuffle command** (`zna shuffle`) for random shuffling of ZNA files with bounded memory
  - Bucket-shuffle algorithm with configurable memory budget (default 1GB)
  - Preserves paired-end read associations (R1+R2 stay together)
  - Deterministic shuffling with `--seed` option for reproducibility
  - Suitable for ML training data preparation
- **Encode with shuffle** (`zna encode --shuffle`) to shuffle during encoding
  - Eliminates need for separate shuffle step
  - Uses same memory-bounded algorithm as standalone shuffle command
- New `_shuffle.py` module with clean `shuffle_zna()` API for programmatic use

### Changed
- Refactored CLI: shuffle logic extracted from `cli.py` to `_shuffle.py` for maintainability
- Consolidated `--block-size` argument (now shared across encode and shuffle commands)

### Performance
- **Decode optimization**: Added `yield from` fast path when `restore_strand=False` (10% faster)
- **Encode optimization**: New `write_records()` batch method caches attribute lookups (20% faster)
- Overall throughput remains: 165 MB/s encode, 241 MB/s decode with 9.10x compression

### Documentation
- Updated README.md with shuffle command documentation and examples
- Updated PERFORMANCE.md with recent optimization notes
- All 86 tests passing

## [0.1.5] - 2026-02-04

### Changed
- Columnar block storage format for 2.5x better compression (9.10x vs 3.68x)
- Optimized flag operations (removed enum overhead for 60-135% speedup)
- Default block size: 128 KB (512 KB recommended for archival)

### Performance
- Compression: 9.10x (9.44x with 512 KB blocks)
- Speed: 165 MB/s encode, 241 MB/s decode (pure Python)
- File size competitive with BAM format

## [0.2.0] - 2026-02-03

### Added
- C++ acceleration with nanobind (9.5x speedup over pure Python)
- Block-based architecture for streaming I/O
- Strand-specific library support (dUTP, TruSeq protocols)
- N-policy handling (drop, random, replace)

### Performance
- Short reads: 189.5 MB/s encode, 668.8 MB/s decode
- Long reads: 2,824 MB/s encode, 3,393 MB/s decode

## [0.1.0] - 2026-01-15

### Added
- Initial ZNA format implementation
- 2-bit nucleotide encoding
- Zstd compression support
- CLI tools: encode, decode, inspect
- Basic FASTQ/FASTA support
- Paired-end and interleaved read handling

[0.1.6]: https://github.com/mkiyer/zna/compare/v0.1.5...v0.1.6
[0.1.5]: https://github.com/mkiyer/zna/releases/tag/v0.1.5
