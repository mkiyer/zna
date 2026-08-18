# The stats trailer, the stored block index, and `zna verify` — implementation plan

Status: **specified, not built.** Written 2026-08-18 against `main` at `f9697d2` (v0.4.1 + 1
commit, 492 test functions green). Line references were verified at that commit and will drift;
resolve them by symbol, not by line.

Companion documents: [MERGE_PAIRS_PLAN.md](MERGE_PAIRS_PLAN.md) (rides the same release, specified
separately), [RELEASING.md](RELEASING.md) (the release this lands in), and khorana's
`docs/DATA_SERVICE_DESIGN.md` §6.1/§17 (the consumer whose problem this solves).

---

## 0. Why

A `.zna` file cannot answer four questions about itself, and its main consumer (khorana's data
service) needs all four before a file may enter training:

1. **What wrote it?** The header carries a format version byte, not a writer version. The only
   place the writer version exists is `merge.json`'s `tool_version` — a QC sidecar written by
   `zna merge`, one pipeline stage *earlier* than the file it is used to vouch for. khorana's C5
   requirement ("writer ≥ 0.4.1, because `zna shuffle` ≤ 0.3.5 cleared `IS_RC` and pre-0.4.1
   blocks were not fragment-complete") is enforced against that sidecar.
2. **What does it hold?** No record count, base count, or length distribution exists in the file
   (`block_index()` docstring, `core.py:929-934`). khorana's sampling weights need the length
   histogram of all records and of the unpaired records; today both come from `merge.json`.
3. **Was it shuffled?** khorana's C1 ("record order random, or every balance guarantee fails
   silently") is un-verifiable from the file, so `khorana-corpus build --shuffled` is an operator
   attestation. `zna encode --shuffle` *knows*, and the knowledge is thrown away.
4. **Is it intact?** No checksum exists anywhere in the format — not per block (the zstd
   compressor is built without `write_checksum`, `core.py:371`), not per file. A file truncated
   exactly at a block boundary reads back clean, just shorter. Today the only thing that catches
   it is the `merge.json` `emitted_records` cross-check — the same sidecar again.

The sidecar dependency is the single largest operational defect class in khorana's data service:
its first stager fetched `.zna` files without their sidecars, every transfer reported success, and
every file went to `state='error'` at the index step. The fix is not a better sidecar; it is a
file that describes itself.

There is also a fifth question, from khorana's Phase 5 (ranged reads over remote storage): **where
are the blocks?** `block_index()` walks the file, so indexing a remote file requires downloading
it — that is *why* khorana's v1 needs a stager at all. khorana's design doc already concluded
"a footer inside the `.zna` is strictly better [than a `.idx` sidecar] but is a format change"
(its open question 9), and this repo's own ROADMAP proposed a sidecar index only "to avoid a
format-version bump". §1 below gets the footer without the bump, which supersedes both.

---

## 1. What changes on disk

**The format version byte stays 3.** The change is purely additive: every byte of the header and
every data block is written exactly as today, and three new structures are appended after the
final data block:

```
[header, unchanged]
[block 1] ... [block N, unchanged]
[SENTINEL: a 20-byte block header with count = 0]      <- terminates forward walks
[TRAILER PAYLOAD: canonical JSON, zstd iff the file is zstd]
[FOOTER: fixed 32 bytes ending in a magic]             <- O(1) discovery from EOF
```

### D1 — stay on format version 3; the extension is additive

A version bump to 4 was considered and rejected:

- v3 is days old (0.4.1 shipped it as a deliberate format break; CHANGELOG "Read the notes, not
  the number"). Bumping again immediately would invalidate the `zna>=0.4.1` pin ecosystem for a
  change that removes nothing and reinterprets nothing.
- A 0.4.1 reader that meets an unknown version raises a message that says "Re-encode from FASTQ"
  (`core.py:1328-1335`) — wrong advice for a file it could in fact read perfectly.
- Under the sentinel design a 0.4.1 reader reading a 0.5 file **degrades benignly**: the walk
  loops decode the trailer as a valid empty block and yield nothing extra. (Verified reasoning in
  §10 traps — this holds only because the payload is a real zstd frame with an honest
  `uncomp_size`.) The one 0.4.1 behavior that is *not* benign — `block_index()` appending a
  phantom `BlockInfo` with `n_records=0` — inflates `n_blocks` by one and costs a consumer at
  most an empty visit, never a wrong record.

What is given up by not bumping: the format itself cannot declare "a trailer is mandatory", so a
0.5-written file that lost its trailer (truncation) is indistinguishable *at the format level*
from a 0.4.1 file that never had one. The refusal therefore lives in policy: `zna verify` and
khorana's indexer treat a missing trailer as failure, and the message names both causes ("written
by zna < 0.5, or truncated — re-encode from source either way"). Since any truncation also strips
the footer, the detection is equally strong; only the diagnosis is a two-way ambiguity.

### D2 — the sentinel is a count-0 block header

`_flush_block` returns immediately on an empty batch (`core.py:705-706`), so **no real block has
ever had `count == 0`**; the value is free. The sentinel is an ordinary 20-byte block header
(`_BLOCK_HEADER_FMT = "<IIIII"`) with:

| field | value |
|---|---|
| comp_size | byte length of the trailer payload as written |
| uncomp_size | byte length of the payload after decompression (== comp_size when uncompressed) |
| count | **0** — reserved as the sentinel from 0.5.0 on |
| flags_size | 0 |
| lengths_size | 0 |

Forward walks (`_iter_records`, `_iter_blocks`, `block_index` — the three loops that today break
on an empty read) additionally break on `count == 0`. Old files without a sentinel still
terminate on EOF; both conditions stay. This is what keeps `records()` working on **pipes**,
where no seek to a footer is possible.

### D3 — the payload is canonical JSON, compressed like the file

Deterministic serialization: `json.dumps(..., sort_keys=True, separators=(",", ":"),
ensure_ascii=True)`, encoded UTF-8, then zstd-compressed at the file's `compression_level` iff
`header.compression_method == 1` (it must be a genuine frame — see D1). **No timestamps, no
hostnames, no filesystem paths, nothing derived from wall clock**: khorana's fixture cache and
this repo's own determinism stance (`_pycodec.py:239-243`, output independent of block size)
both rely on byte-identical re-encodes of identical input.

Fields, all required unless marked optional:

| field | type | semantics |
|---|---|---|
| `trailer_schema` | int, `1` | schema version of this payload, independent of the format byte |
| `writer_version` | str | `zna.__version__` of the **last writer** of these bytes. `zna shuffle` and ZNA→ZNA re-encode stamp their own. This is khorana's C5 input. |
| `shuffled` | bool | true iff the record order was produced by `zna shuffle` / `encode --shuffle`, or preserved verbatim from an input that was (see §4). khorana's C1 input. |
| `n_records` | int | total records; must equal the sum over the block index |
| `n_bases` | int | total stored sequence bases (post-npolicy, i.e. of the bytes actually in the file) |
| `n_pairs` | int | count of R1+R2 pairs |
| `n_unpaired` | int | count of unpaired records (merged reads and singles). Invariant, asserted at close: `n_unpaired + 2*n_pairs == n_records`. |
| `flag_counts` | {str(int): int} | per-flag-byte record counts, zero bytes omitted. Subsumes the two counts above; kept anyway for legibility and audit. |
| `length_histogram` | {str(int): int} | stored record length → count, over **all** records |
| `length_histogram_unpaired` | {str(int): int} | same, over **unpaired records only**. This is the field khorana today reconstructs from merge.json's `insert_size_histogram` ("a merged read's length IS its insert size") — named for what it is, not for the inference it replaces. |
| `blocks` | object | the stored index: `{"data_start": int, "comp_sizes": [int], "uncomp_sizes": [int], "records": [int]}`. Absolute offsets are reconstructed by running sum: `offset[i] = data_start + Σ_{j<i} (20 + comp_sizes[j])`. |
| `merged_in_process` | bool, optional | present and true when written by `zna encode --merge-pairs` — answers MERGE_PAIRS_PLAN.md open question §10.4 (no other file-level provenance exists) |

Histograms and stats are derived **from the encoded columns** (`flags_bytes`, `lengths_bytes`)
accumulated in `_flush_block`, never from the writer's input arguments — see trap T1.

### D4 — the trailer is the single source of truth; no new header bits

`writer_version` and `shuffled` could also live in the front header, but header flag bits do not
round-trip: a 0.4.1 reader masks the four known bits and drops the rest (`core.py:1360-1363`),
and `_shuffle.py` rebuilds headers field-by-field (`_shuffle.py:172-182`). Two copies that can
disagree are worse than one that cannot. The header wire format is untouched.

### D5 — the footer, fixed 32 bytes at EOF

```
<8s H H I Q 8x   little-endian, 32 bytes total
 magic = b"ZNATRLR\x1A"
 footer_version = 1
 reserved (0)
 crc32 of the payload bytes AS WRITTEN (compressed form) — checkable without inflating
 sentinel_offset = byte offset of the sentinel block header, i.e. THE END OF DATA BLOCKS
 8 reserved bytes (0)
```

`sentinel_offset` is deliberately the value khorana's `pack_block_offsets` needs as its new final
sentinel (today it uses `file_size`, with a docstring asserting "the last block ends exactly at
EOF" — true through 0.4.1, false the moment a trailer exists; see §7). Discovery: seek to
EOF−32, check magic, read `sentinel_offset`/payload. A file smaller than 32 bytes or without the
magic has no trailer.

### D6 — zstd content checksums, on

Every zstd frame (blocks and trailer payload) is now compressed with `write_checksum=True`.
python-zstandard verifies content checksums automatically on decompress, so every block read
becomes an integrity check for free (xxhash64, negligible cost). Frames with checksums remain
ordinary frames to any zstd — 0.4.1 readers are unaffected. This is the format's first
protection against bit rot, and `zna verify` inherits it without writing any comparison code.

---

## 2. Writer changes (`core.py`)

- `ZnaWriter` gains accumulators, updated in `_flush_block` **from the values the codec
  returned** (`flags_bytes`, `lengths_bytes`) between payload assembly and the write: per-flag
  counts, per-length counts overall and for unpaired flag bytes, base total, per-block
  `(comp_size, uncomp_size, count)` list. See trap T1 for why the writer's own inputs are the
  wrong source.
- `ZnaWriter.__init__` gains `shuffled: bool = False` (stored, stamped into the trailer).
- `ZnaWriter.close()` — after the existing dangling-R1 check and final `_flush_block()` — writes
  sentinel, payload, footer. Append-only throughout; `zna encode` writing to stdout keeps
  working. A writer that never wrote a record still writes the trailer (a 0-record file is
  header + sentinel + trailer + footer, and verifies).
- `write_copy`/`copy_records` need no change: they funnel through the same flush path, so the
  accumulators see the copied flag bytes verbatim.

## 3. Reader changes (`core.py`)

- The three walk loops break on `count == 0` in addition to the empty read.
- New cached property `ZnaReader.trailer -> ZnaTrailer | None` (a `slots` dataclass mirroring the
  payload, plus `.raw`): seek EOF−32, magic check, crc check, inflate, parse. `None` on pipes and
  on trailer-less files.
- `block_index()` returns the **stored** index when a trailer is present (reconstructing offsets
  by running sum) and falls back to the scan otherwise. The scan is retained forever: it is
  `zna verify`'s cross-check.
- **Fix the phantom-block bug while in there** (pre-existing, latent): today `block_index()` on a
  file with garbage after the last block appends a garbage `BlockInfo`, seeks past EOF without
  error, and returns success (`core.py:966-986`). After this change, non-sentinel trailing bytes
  ≥ 1 raise a `ValueError` naming the offset. The same silent-`break` exists in
  `inspect_command`'s hand-rolled walk (`cli.py:1862-1875`) — replace it with `block_index()`.

## 4. `zna shuffle` and ZNA→ZNA re-encode

- `shuffle_zna` constructs its output writer with `shuffled=True`; the fresh trailer is automatic
  because the output goes through `ZnaWriter`. Bucket temp files also get trailers (harmless,
  small K, not worth a code path to suppress).
- The ZNA→ZNA re-encode path (`cli.py:1418-1446`) copies records in order; its output propagates
  the *input's* `shuffled` value when the copy preserves record order, else `False`. The
  post-encode `--shuffle` pass replaces the file via `shuffle_zna` and stamps `True` regardless.
- Both paths stamp their own `writer_version` — "last writer" is the C5 contract.

## 5. `zna verify` — the certification command

`zna verify [--json] [--fast] FILE...`, exit 0 iff every file passes.

Full mode decodes every block and checks, in order, stopping at the first failure per file with a
message that names the cause:

1. Footer present, magic and crc valid, payload parses, `trailer_schema` known. Absence: *"no
   trailer: written by zna < 0.5, or truncated — re-encode from source"*.
2. Stored block index matches the walked file exactly (offsets, sizes, counts) and the sentinel
   sits where the footer says.
3. Every block decompresses (zstd content checksums verified implicitly by D6) and its columns
   split cleanly (`flags_size == count`, `lengths_size == count * seq_len_bytes`).
4. Recounted stats equal the trailer: flag counts, both histograms, `n_records`, `n_bases`.
5. Structural invariants: every paired R1 immediately followed by its R2; no block ends on an
   unmatched R1 (fragment-complete blocks, the 0.4.1 guarantee); header strand flags consistent
   (stranded ⇒ exactly one antisense mate).

`--fast` runs check 1 plus an index-vs-block-header walk without decompressing payloads (the
seek-only scan `block_index` already does). What verify certifies — and what it cannot: it proves
integrity, structure, and stats fidelity. It cannot prove `shuffled` (a permutation attests
itself no better than a flag) or that strand normalization was applied correctly; those remain
stamped facts trusted from the writer.

## 6. `zna inspect --json`

Splices the trailer fields in verbatim under a `"trailer"` key (or `"trailer": null`), keeps its
existing header/block reporting, and switches its hand-rolled walk to `block_index()` (§3).

## 7. The consumer contract (khorana)

What khorana reads, where it comes from after 0.5, and which requirement it serves:

| trailer field | khorana consumer | requirement |
|---|---|---|
| `writer_version` | `check_zna_version`, floor becomes `(0, 5, 0)` | C5 |
| `shuffled` | catalogue `files.shuffled`; the `--shuffled` build flag is dropped | C1 |
| `length_histogram` | `FileFacts.histogram` → `usable_fragments`, `usable_bases`, `slide_records`, `bases` | sampling weights |
| `length_histogram_unpaired` | `FileFacts.merged_lengths` (replaces the insert-size inference) | fragment weights |
| `n_unpaired`, `n_pairs`, `n_records` | `FragmentCounts` — now exact by construction, not inferred | fragment weights |
| `blocks` | catalogue `block_offsets`/`block_records`; Phase 5 ranged reads (suffix-range GET of the trailer indexes a remote file with no full download) | serving, Phase 5 |
| footer `sentinel_offset` | `pack_block_offsets`'s final sentinel (**replaces `file_size`** — its "last block ends exactly at EOF" assumption is false once a trailer exists) | serving |
| (header, unchanged) `strand_normalized` | `_open_zna` refusal | C2 |

Acceptance criterion for the khorana side (tested in khorana's repo when it adopts 0.5): a
recipe-built library indexes with **no `merge.json` present**, and the resulting catalogue row
equals, field for field, what today's merge.json path produces. `merge.json` reverts to what it
always was — `zna merge`'s QC report — and khorana deletes `parse_merge_json`, the sidecar-beside-
the-file contract, and the stager's two-file fetch. Staging certification becomes `zna verify`
after every fetch.

## 8. `--merge-pairs` rides the same release

[MERGE_PAIRS_PLAN.md](MERGE_PAIRS_PLAN.md) is the specification; it is not restated here. Three
interactions with this plan, resolving two of its open questions:

- Its output goes through `ZnaWriter`, so the trailer is automatic; it additionally stamps
  `merged_in_process: true` (§10.4's provenance question — answered here).
- Its §10.1 question (should `--merge-json`'s `tool` field change?) is mooted: khorana stops
  reading merge.json entirely, so the field stays `"zna-merge"` for hulkrna compatibility and no
  consumer breaks.
- Its §10.2 request ("assert `emitted_records` equals the `.zna`'s records") is subsumed by
  `zna verify` check 4.

Build the trailer first (§9): it is what khorana blocks on; `--merge-pairs` is a speed and
robustness upgrade whose payoff is largest during the corpus re-encode, which must wait for the
trailer anyway.

## 9. Build order

1. **PR 1 — the format extension.** Writer accumulators + sentinel/payload/footer +
   `ZnaWriter(shuffled=)` + D6 checksums; reader walk termination + `trailer` property +
   stored-index `block_index()` + the phantom-block fix; shuffle/re-encode stamping (§4).
   Tests with it (§10). This PR alone unblocks khorana.
2. **PR 2 — `zna verify` + `inspect --json` trailer splice** (§5, §6). CHANGELOG entry leads
   with the format extension, per RELEASING.md's format-break discipline.
3. **PR 3 — `zna encode --merge-pairs`**, per MERGE_PAIRS_PLAN.md's own build order §4,
   including its three audit blockers (§0.1–§0.3 there).
4. **PR 4 — docs.** Add this file and the new commands to README's documentation index; ROADMAP:
   mark the sidecar block index item **superseded by the trailer**, move shipped items to the
   CHANGELOG; fix the three stale claims found while auditing for this plan (ROADMAP's "Known
   divergence between the codec backends" describes a divergence that 0.4.0 fixed —
   `tests/test_fuzz_roundtrip.py:699` now proves both backends raise; same stale claim at
   `.github/copilot-instructions.md:83` and `docs/METHODS.md:456`; `scripts/merge_bench/README.md:6`
   says x86 tuning is "0.4.1", ROADMAP says 0.4.2).
5. **Release 0.5.0** per RELEASING.md: both-env suites, `release.sh`, conda flow. khorana then
   bumps its pin to `zna>=0.5.0` and begins its consumption pass; the corpus re-encode starts
   only after both.

## 10. Tests

Follow the house conventions: behavioral names, both backends, compiled and extension-less envs
(~56 skips expected there), `_pycodec` as specification.

- **Fuzz integration** (`test_fuzz_roundtrip.py`): every fuzz-generated file now also asserts
  (a) trailer present and parses, (b) recount-from-decode equals trailer exactly — the whole
  matrix (npolicy × strand × compression × layout × backend) becomes a trailer-correctness
  instrument for free. This is the load-bearing test.
- **Determinism**: encode identical input twice (same seed) → byte-identical files including
  trailer and footer.
- **Termination**: `records()`/`blocks()`/`block_index()` on trailer-bearing files, including
  over a pipe (sentinel termination, no seek); trailer-less 0.4.1-written files still read.
- **Damage matrix** for `zna verify`: truncation mid-block, at a block boundary (the previously
  silent case), mid-trailer, footer-only corruption, a flipped payload byte (crc), a flipped
  block byte (zstd checksum), stats tampered — each fails with the intended message; `--fast`
  catches the subset it claims.
- **Shuffle/re-encode stamping**: `shuffled` and `writer_version` propagate per §4; a shuffled
  file's trailer recounts identically (permutation invariance of every stat).
- **Edge**: 0-record file (with `--allow-empty` provenance) verifies; `seq_len_bytes` ∈ {1,2,4};
  labeled files (label columns don't perturb stats); `count == 0` never emitted by any data path
  (assert in `_flush_block`).
- **Back-compat probe**: a 0.5 file read by the *current* code with the new termination logic
  removed (simulated) yields the benign-degradation behavior D1 claims — pin it so the claim
  stays true if the sentinel layout is ever touched.

## 11. Traps

- **T1 — accumulate from the codec's output, not the writer's input.** The codec ORs `IS_RC`
  into the flags column during strand normalization, and `--npolicy trim3` changes stored
  lengths (it can produce empty sequences). Stats tallied from `write_record`'s arguments would
  disagree with a recount of the file. Tally from the returned `flags_bytes`/`lengths_bytes` in
  `_flush_block`, and the fuzz recount test enforces it across every npolicy/strand combination.
- **T2 — the payload must be a real frame with an honest `uncomp_size`.** D1's benign
  degradation for 0.4.1 readers depends on it (`max_output_size` is passed at decompress,
  `core.py:1152`). An uncompressed file's trailer must be raw JSON for the same reason.
- **T3 — `_shuffle.py` rebuilds headers field-by-field.** Any future header field must be added
  there too or it silently vanishes on shuffle (`_shuffle.py:172-182`). This plan adds no header
  field, which is D4's point — but the trap stands for whoever does.
- **T4 — byte-identity is per-environment, not universal.** Compressed frame bytes depend on the
  bundled libzstd; a `zstandard` upgrade may change them for identical input. khorana's fixture
  cache keys on `zna.__version__` — an upgrade of `zstandard` alone can invalidate byte-level
  assumptions without a zna version change. Determinism (D3) means: same input, same seed, same
  environment ⇒ same bytes.
- **T5 — the stale-`.so` build trap** (HANDOFF_0.4.0 traps 1–2): a failed CMake build leaves the
  old extension installed and everything imports fine. `pip install -e . --no-build-isolation`
  and grep for `Successfully built`. The C++ codec is untouched by this plan (all trailer logic
  is Python; C++ never parses block headers), which is what keeps PR 1 small — but the fuzz
  suite must still run against both backends.
- **T6 — `inspect`'s silent `break`.** The hand-rolled walk breaks silently on a short read
  (`cli.py:1865-1866`) and would under-report rather than error on damage; §3/§6 replace it.
  Until then, do not use `inspect` output as verification.

## 12. Considered and not taken

| item | why not now |
|---|---|
| Format version bump to 4 | D1 — additive suffices; bump costs the ecosystem and buys only a cleaner refusal message. Revisit only if a change ever *reinterprets* existing bytes. |
| `writer_version`/`shuffled` as header bits | D4 — bits don't round-trip 0.4.1 and `_shuffle.py` rebuilds headers; one source of truth. |
| Per-column compression, label min/max indexes, string/categorical label dtypes (ROADMAP "considered") | Real format work with no current consumer; none is additive-only; each can ride a future genuine bump. Not this release. |
| x86 merge-kernel tuning (ROADMAP 0.4.2) | Orthogonal to the format; needs a Linux box; unchanged schedule. |
| msgpack/CBOR payload | JSON is stdlib, self-describing, debuggable, and small (a 100k-block file's index is ~2 MB raw, ~hundreds of KB compressed); zero new dependencies. |
| khorana-side facts derivation (decode a bare `.zna` in khorana) | Superseded by this plan: the file describes itself, `zna verify` recounts, and khorana never re-implements format arithmetic. Pre-0.5 files get "re-encode with current zna". |

## 13. Open questions for the owner

1. `zna verify` default depth at khorana's staging call site: full decode (the certification the
   corpus owner asked for; one decode per staged file) with `--fast` as the local-corpus option —
   or the reverse. This plan assumes full-by-default.
2. Should the trailer also carry `read_group`-level free-text provenance (e.g. the merge
   parameter block that today lives only in merge.json's `params`)? Deterministic and small, but
   it duplicates QC material khorana never reads. This plan says no; easy to add under
   `trailer_schema: 2` later.

---

## 14. Amendments — 2026-08-18, at implementation

Recorded before the code was written, so the plan above stays honest about what it got wrong
and what the owner overrode. Where an amendment contradicts a section above, the amendment wins.

### A1 — the sentinel's `comp_size` covers the payload AND the 32-byte footer

D2's table says `comp_size` is "byte length of the trailer payload as written". Implemented that
way, D1's benign-degradation claim is **false**: a 0.4.1 walk, after swallowing the sentinel
block, reads the footer's first 20 bytes as a block header — `comp_size` becomes the magic's
first four bytes (~1.4 GB), `_read_exact` hits EOF, and every `records()`/`blocks()` call on a
0.5 file ends in `EOFError` after the last real record. `block_index()` appends **two** phantoms,
one of them claiming `n_records=1` (the footer's version field). D1 promised one phantom with
`n_records=0` and clean walks; only one layout delivers that: `comp_size = len(payload) + 32`, so
an old reader consumes payload and footer in one gulp and its next read is a clean EOF.

This works because python-zstandard's one-shot `decompress()` tolerates trailing bytes after a
complete frame by default (verified on the installed 0.25.0: `allow_extra_data=True` is the
default), so the footer bytes riding inside the sentinel's `comp_size` do not break the old
reader's decompress of the trailer frame. `uncomp_size` stays honest (the raw JSON length), which
is what bounds `max_output_size`. New readers never use the sentinel's `comp_size` — discovery is
footer-first (seekable) or count==0-terminated (pipes) — and `inspect --verify` checks the
relation `comp_size == (EOF-32) - (sentinel_offset+20) + 32` explicitly.

### A2 — provenance moves to a PROLOGUE; the trailer is derived facts only

The owner's ruling on §13.2, generalized past the question asked: **position in the file follows
when the information is knowable.** `writer_version`, `shuffled`, and `merged_in_process` are
inputs — known before the first record is buffered — so they do not belong in a structure whose
reason to exist is "facts that require the whole encode". D3/D4 put them in the trailer only
because the plan assumed the trailer was the sole additive extension point. It is not: the same
count-0 pseudo-block mechanism works at the *front* of the file.

The layout becomes symmetric:

```
[header, unchanged]
[PROLOGUE: count-0 block header + provenance JSON]     <- written at open; start-known facts
[block 1] ... [block N]
[SENTINEL: count-0 block header + trailer JSON]        <- written at close; derived facts
[FOOTER: 32 bytes]
```

- **Prologue payload** (canonical JSON, compressed like the file): `prologue_schema: 1`,
  `writer_version`, `shuffled`, and `merged_in_process` (present iff true). Nothing else — the
  owner's ruling also excludes parameter echoes (merge thresholds, npolicy): those are QC
  material, not file facts, and stay out entirely.
- **Trailer payload** drops those three fields and keeps only derived facts: `trailer_schema`,
  `n_records`, `n_bases`, `n_pairs`, `n_unpaired`, `flag_counts`, both histograms, `blocks`,
  plus `prologue_crc32` — the end attests the start, closing the integrity gap for uncompressed
  files whose prologue has no zstd checksum.
- **Reader**: `ZnaReader.__init__` probes the first block header. count==0 ⇒ prologue: consume,
  parse, expose as `ZnaReader.provenance`; else ⇒ 0.4.1-era file, the 20 bytes are rewound
  (seekable) or stashed for the first walk read (pipes), `provenance` is None. Every walk
  thereafter is prologue-blind and `count == 0` has exactly one meaning: end of data. A pipe
  consumer therefore learns `writer_version`/`shuffled` **before** decoding a single record —
  something the trailer-only design could never do, and the concrete payoff of the ruling.
- **0.4.1 readers** see the prologue as one more valid empty block at the front (same benign
  degradation as the sentinel; block numbering shifts by one for old readers of new files, which
  is the same class of tolerated skew as D1's phantom).
- khorana's §7 table: `writer_version` and `shuffled` come from the **prologue**; the acceptance
  criterion is unchanged in substance — index a library from the *file* alone, no merge.json.

### A3 — no `zna verify` command; certification is `zna inspect --verify`

The owner's ruling on §13.1: one command owns "look at a file". Bare `zna inspect` gains the
cheap structural pass automatically when a trailer exists (footer magic + CRC, stored index vs
walked headers, prologue parse) and reports provenance and trailer facts; `zna inspect --verify`
is §5's full certification — decode every block (zstd content checksums), recount every stat,
compare to the trailer, check structural invariants — exit 0 iff certified. `--json` composes
with both. §5's `--fast` tier is subsumed by bare `inspect`; khorana's staging call becomes
`zna inspect --verify` after every fetch. Recount reads the flags and lengths columns; the
2-bit sequence column's bytes are covered by the zstd content checksum and the
`sum(ceil(len/4)) == len(seq_stream)` cross-check, which validates the lengths column against
the sequence bytes without materializing strings.

### A4 — §5 check 5's strand-flags test is a lint, not a failure

"Stranded ⇒ exactly one antisense mate" is not an invariant of the format: the fuzz matrix
deliberately writes `(strand_specific, r1_antisense=True, r2_antisense=True, normalized)`
configurations and they are legal files. `inspect --verify` reports it as a warning and does not
fail the file. The fragment-adjacency half of check 5 (every paired R1 immediately followed by
its R2, no block ending mid-fragment) is a real invariant and stays fatal.

### A5 — an aborted encode leaves no trailer, deliberately

`ZnaWriter.__exit__` on an exception flushes the buffered records (blocks on disk stay whole,
as today) but does **not** write sentinel/trailer/footer. The trailer is the writer's signature
that the encode *finished*; a crashed encode now yields a file that `records()` still reads
(EOF termination is retained) but `inspect --verify` refuses — which is exactly the refusal
khorana's staging wants. `close()` is idempotent; only a clean close signs the file.
