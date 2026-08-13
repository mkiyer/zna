# `zna encode --merge-pairs` — implementation plan

Status: **specified, not built.** Revision 2, written 2026-08-13 after auditing the
2026-08-13 draft against the code. Revision 1 was written from the design intent;
this one is written from what `src/zna/cli.py`, `src/zna/core.py` and
`src/zna/merge/` actually do, with the load-bearing claims checked by running them.

Baseline at revision 3: `627 passed, 156 subtests passed` in
`/Users/mkiyer/sw/miniforge3/envs/zna_merge` (compiled backend, `zna.is_accelerated()`
True, `zna.merge` backend `accel`).

---

## 0.0 What 0.4.0 changed under this plan — read this before §0

Revision 2 was audited against the code *earlier on the same day* that the `--npolicy`
and provenance work landed, so several of its load-bearing claims are now false. They
are corrected here rather than quietly edited, because each one is a belief a reader of
rev 2 would carry into the implementation.

- **`--npolicy drop` no longer exists.** The policy set is `{trim3, random}`, default
  `trim3`. Every appearance of `drop` below has been rewritten, but if you are working
  from a rev-2 copy: `trim3` **never orphans a mate**, so the fragment-level atomicity
  `drop` needed is gone, and the N filter no longer needs to sit inside
  `_fragment_units` at all. `zna encode` now trims per record.
- **`zna merge` has a full `--npolicy`.** Rev 2 said "there is no N policy anywhere in
  `zna merge` (a grep for `npolicy` across `src/zna/merge/` returns nothing)". That was
  true when written and is now flatly wrong — the policy runs *inside* the kernel, after
  overlap rescue, with a coverage retry. `--merge-pairs` inherits it through `Params`
  rather than needing an encode-side filter, which removes an ordering hazard rev 2
  did not know it had.
- **Per-record provenance exists**, which answers §8 item 4 below.
- **`scripts/wheel_smoketest.py` now imports `zna.merge._accel`** and checks both
  extensions and their collision case, so §8 item 9 is done.

## 0. Corrections to revision 1 — read this first

Revision 1 had three defects that would each have shipped a silently wrong corpus.
They are recorded here rather than quietly fixed, because each one is a plausible
belief that the code refutes, and the next person will have the same belief.

### 0.1 The geometry table gave R2 the wrong flags (§2.2 of rev 1) — **blocker**

Revision 1 said:

| emitted record | has_start | has_end | ⇒ IS_RC |
|---|:--:|:--:|:--:|
| R2 of a kept/trimmed pair | ✗ | ✔ | **1** |

That is wrong, and rev 1's own note two lines below the table contradicts it ("R2 is
emitted in its own orientation and **its base 0 is the fragment's 3' end**" — i.e. base
0 *is* a true boundary, so `has_start` is True).

`has_start`/`has_end` describe the left and right edge of the **stored** sequence, not
of the fragment ([core.py:110-118](../src/zna/core.py#L110-L118): *"A read begins at a
fragment boundary and runs inward, so base 0 is a true boundary; storing it
reverse-complemented moves that boundary to the right edge"*). `IS_RC` is therefore a
**storage** fact owned by the writer's strand-normalization settings — not a geometry
fact the merger knows. `zna merge` emits R2 in its own sequenced orientation, so R2's
base 0 is a real boundary and `IS_RC` must be **0**.

Measured on the existing two-step path, same input, two flag sets:

```
zna merge -> zna encode --interleaved --treat-unpaired-as-merged
    R1      (has_start=True,  has_end=False)  x20      <- IS_RC 0
    R2      (has_start=True,  has_end=False)  x20      <- IS_RC 0
    merged  (has_start=True,  has_end=True )  x40      <- IS_FULL_FRAGMENT

...the same, plus --strand-normalize --strand-specific --read1-sense --read2-antisense
    R1      (has_start=True,  has_end=False)  x20
    R2      (has_start=False, has_end=True )  x20      <- IS_RC 1, set by the WRITER
    merged  (has_start=True,  has_end=True )  x40
```

Transcribing rev 1's table would have stamped `IS_RC` on every R2 in an
un-normalized corpus while its bases sat forward — the file claiming a
reverse-complement that never happened, for half of every unmerged pair. Invisible to
any test that checks sequences.

**The corrected mapping, and it is smaller than rev 1 implied:**

| emitted record | slot | is_paired | is_read1 | is_read2 | has_start | has_end | IS_FULL_FRAGMENT |
|---|---|:--:|:--:|:--:|:--:|:--:|:--:|
| merged (spans the fragment) | `MERGED` | 0 | 0 | 0 | ✔ | ✔ | **1** |
| R1 of a kept/trimmed pair | `MATE1` | 1 | 1 | 0 | ✔ | ✗ | 0 |
| R2 of a kept/trimmed pair | `MATE2` | 1 | 0 | 1 | ✔ | ✗ | 0 |

`has_start` is **True for every record on this path**. So the geometry transfer is
exactly **one bit per record: "is this record a whole fragment?"** Everything else
the encoder already got right. That is a smaller claim than rev 1 made, and it is the
true one — §1 is rewritten accordingly.

### 0.2 `preserve_normalization=True` would silently disable `--strand-normalize` — **blocker**

Rev 1 §3.4 and trap 2 said to set it, and to "confirm that is harmless for a fresh
encode". It is not harmless. [core.py:284-298](../src/zna/core.py#L284-L298) uses the
same flag to force `_do_strand_norm_r1 = _do_strand_norm_r2 = _do_random_norm = False`,
and those three booleans are the codec's **only** path to reverse-complementing
anything. Meanwhile the `STRAND_NORMALIZED` header bit is still written.

So `zna encode --merge-pairs --strand-normalize` would have produced a file whose header
says normalized, whose bases are not, and whose `IS_RC` column is empty — and
`--strand-normalize` is in the flag set the production pipeline uses
([tests/test_merge_encode.py:148-158](../tests/test_merge_encode.py#L148-L158)).

**Do not set `preserve_normalization` on this path.** It is not needed: with the
corrected mapping the adapter never supplies `is_rc=True`, so `write_record`'s guard
([core.py:349](../src/zna/core.py#L349)) never fires. `IS_FULL_FRAGMENT` and `IS_RC`
are independent bits in the flags byte, and `_ends_from_flags` returns `(True, True)`
for a full-fragment record *regardless* of `IS_RC` — so the caller can own the fragment
span while the writer owns orientation, and the two compose.

Verified by prototype: `--merge-pairs` + `--label-defs` + `--strand-normalize
--strand-specific --read2-antisense`, 81 records over all three outcomes, **100%
geometry-correct against the known fragments**, with `IS_RC` correctly set on R2 by the
writer and `preserve_normalization` never touched. (The prototype ran under the
then-default `--npolicy drop`; the geometry it verified is independent of the policy,
which acts on read *content*, not on flags.)

### 0.3 `--label-defs` wins the stream selection, so the combination would silently not merge — **blocker**

[cli.py:1266](../src/zna/cli.py#L1266) is `if label_defs: stream =
stream_inputs_labeled(...)`, tested **before** the geometry branch.
`stream_inputs_labeled` knows three input modes and, given two files, emits an ordinary
un-merged R1/R2 stream. So `zna encode --merge-pairs --label-defs ...` — the exact
invocation the production pipeline wants — would have parsed the flag, ignored it, and
written an unmerged corpus with a zero exit status.

`--merge-pairs` must be dispatched **before** `label_defs` in that chain.

### 0.4 Smaller corrections

- **`--treat-unpaired-as-merged` is not wrong on `zna merge` output.** Rev 1 §3.3 says
  it "stamps IS_FULL_FRAGMENT on every unpaired record including any that were never
  merged (a mate whose partner was filtered)". `merge_core.hpp` guarantees pair
  atomicity (`n_recs` is 0 or 2 for an unmerged pair,
  [merge_core.hpp:498-502](../src/zna/merge/merge_core.hpp#L498-L502)), and the
  benchmark measured C3 orphans at **0**. So on well-formed merge output the flag is
  *correct*, and §5.1's "assert the difference by name" would have asserted a difference
  that is not there. The real failure mode is the name-based re-pairing in §1, and it is
  worse.
- **Only two tables derive from `_flags_from_ends`**, not three: `_RC_FULL_BY_ENDS` and
  `_FLAG_BITS_BY_ENDS` (which derives from the first). `ENDS_BY_FLAG`
  ([core.py:130](../src/zna/core.py#L130)) is an independent restatement on the decode
  side — a real drift risk rev 1's trap 3 did not name.
- **`_flags_from_ends` is not a bijection**, and its docstring says it is. See §8.1 —
  a pre-existing corruption, not one this change introduces, but adjacent. It has a
  second head: ZNA→ZNA re-encode of a **non**-normalized file silently drops
  `IS_FULL_FRAGMENT` altogether. Verified — `encode m.fq --treat-unpaired-as-merged`
  gives `(has_start, has_end) = (True, True)`; re-encoding that `.zna` gives
  `(True, False)`. Same family, same fix, same rebuild.
- **Trap 5 (the two `_accel` targets colliding) is already fixed** and pinned by
  `TestExtensionsAreDistinct`. Keep it as a caution for a third extension; it is not
  outstanding work.
- **`names.py` is not on the compiled path at all** — the compiled backend
  re-implements `base_name`/`strip_pair_suffix` in C++. Rev 1 §8's deletion list
  overstates what goes away.

---

## 1. What this is, and why now

Today, merging and encoding are two processes joined by a file:

```
R1.fq.gz ─┐
          ├─► zna merge ─► merged.fq.gz ─► zna encode --interleaved ─► out.zna
R2.fq.gz ─┘                    ▲                     ▲
                        format + deflate      inflate + re-parse
```

`--merge-pairs` deletes the middle. Two reasons.

**Speed.** The merge kernel is 0.42 µs/pair of compute against ~1.5 µs/pair end to end.
The intermediate is most of the I/O: one FASTQ format, one deflate, one write, one read,
one inflate, one full re-parse — plus quality strings that ZNA does not store,
formatted and parsed for nothing.

**Correctness — and state this one accurately.** The two-step path is not merely
"inferring what the merger knew"; it is re-deriving the *pairing itself* from read
names, and that inference has a silent failure mode. `zna encode --interleaved` pairs
**consecutive records whose base names match**
([cli.py:761-809](../src/zna/cli.py#L761-L809)). Two independently merged molecules that
share a read name are therefore re-paired into one fragment. Reproduced:

```
input: 4 fragments, all merging; two of them share the read name "dup"

merged FASTQ                       two-step .zna
@dup   merged_90_0      ──►  len= 90 paired=True  R1=True   FULL_FRAGMENT=False
@dup   merged_95_0      ──►  len= 95 paired=True  R2=True   FULL_FRAGMENT=False
@uniqA merged_100_0     ──►  len=100 paired=False           FULL_FRAGMENT=True
@uniqB merged_105_0     ──►  len=105 paired=False           FULL_FRAGMENT=True
```

Two whole molecules lose `IS_FULL_FRAGMENT`, are marked `IS_PAIRED` with `IS_READ1` /
`IS_READ2`, and are encoded as two mates of one fragment. No warning, exit 0.
Duplicate read names occur in real data (concatenated lanes, some SRA dumps,
UMI-collapsed files), and nothing upstream forbids them.

In process, each record's geometry is *handed over* rather than reconstructed, so this
class of error is structurally impossible — and `--treat-unpaired-as-merged`, a blanket
assertion about a whole file, stops being needed at this seam.

**But be honest about the size of the correctness win, because it was measured.** Run on
the real 1M-pair benchmark, the two-step path emits 1,416,630 records (matching
[MERGE_BENCHMARK_RESULTS.md](MERGE_BENCHMARK_RESULTS.md) §4's C1 count) with this
geometry:

```
merged   (has_start=True,  has_end=True )   582,362     <- IS_FULL_FRAGMENT
R1       (has_start=True,  has_end=False)   417,134     <- IS_RC 0
R2       (has_start=True,  has_end=False)   417,134     <- IS_RC 0, NOT 1
```

`IS_FULL_FRAGMENT` sits on exactly the 582,362 merged records and nothing else, and a
correct `--merge-pairs` produces the **identical** histogram. On well-formed,
uniquely-named input **no flag changes at all**. (This is also §0.1 confirmed at scale:
rev 1's table would have written `IS_RC` on 417,134 records per million.) Rev 1 built
its urgency
on "this changes which records carry `IS_FULL_FRAGMENT`, so land it before the rebuild";
that premise is false, and the rebuild argument should not be made on it.

What is actually gained: the speed below, and the removal of two *conditional* failure
modes — the name collision above, and a blanket flag that is only correct because the
merger happens to guarantee pair atomicity. Neither fires on the pipeline's current
input. That is a good reason to build this; it is not a reason to rush it ahead of the
rebuild.

**A caveat on `IS_FULL_FRAGMENT` itself, which rev 1 cited backwards.** Rev 1 §2.2 said
the mapping restates "contracts C1/C2/C4 … verified at 0 violations". C1, C3 and C4 are
0; **C2 is 5,591** per million, and
[MERGE_BENCHMARK_RESULTS.md](MERGE_BENCHMARK_RESULTS.md)'s own headline says C2 "does not
hold universally". So 0.96% of merged records (5,591 of 582,362) are not their true
fragment, and their right
edge is a false boundary. `--merge-pairs` neither creates nor fixes that — the two-step
path stamps the identical records — but the flag's meaning is "the merger asserts this
spans the fragment", not "this spans the fragment", and the plan should not claim
otherwise.

**Speed, measured** on the 1M-pair benchmark (M3 Max): `zna merge --threads 4` **1.62 s**
wall (6.13 s user, 419% CPU); `zna encode --interleaved --treat-unpaired-as-merged`
**2.89 s** wall (4.06 s user). The output deflate costs CPU but ~zero wall — it is
overlapped by the writer's io threads. So the prize is the encode side's
inflate-and-reparse, not the whole 4.5 s; budget the claim accordingly and re-measure
after building.

**Timing.** khorana has committed to a corpus rebuild and to dropping legacy formats.
Landing this before the rebuild is a convenience, not a correctness deadline (see
above). The thing that *does* want fixing before the rebuild is §8.1.

---

## 2. The seam already exists

### 2.1 The encoder already has a geometry-carrying path

`stream_inputs(args, with_ends=True)` yields 6-tuples `(seq, is_paired, is_read1,
is_read2, has_start, has_end)`, and the write loop at
[cli.py:1307-1314](../src/zna/cli.py#L1307-L1314) turns the last two into flags via
`_RC_FULL_BY_ENDS`. This exists for ZNA→ZNA re-encode. `--merge-pairs` is its second
caller — but see §0.2: it is the *shape* that is reused, not `preserve_normalization`.

Note the loop already lives **inside** `for unit in _fragment_units(stream)`. Rev 2 said
to keep it there for `--npolicy drop`'s pair atomicity; `drop` is gone and `trim3` never
orphans a mate, so that reason has expired. Keep it there anyway, for a different and
still-live one: `_fragment_units` is what groups mates so `_full_fragment_flags` can see
a whole fragment, and a write loop that hoisted `carries_ends` above it would break the
grouping (§6 test 8).

### 2.2 The merge tool already computes the one bit that matters

`OUTCOME_MERGED` ⇒ `has_end=True`; everything else ⇒ `has_end=False`. `has_start` is
always True. That is the whole transfer (§0.1).

### 2.3 The C++ core is already I/O-agnostic

`merge_core.hpp` knows nothing about FASTQ; `fastq_chunk.hpp` is one adapter over it.
This work adds the second adapter beside it. One caveat rev 1 missed: the fastp-style
`merged_<n1>_<n2>` name is built **inside the core**
([merge_core.hpp:441-467](../src/zna/merge/merge_core.hpp#L441-L467)), not in the
adapter — so "reuse the core untouched" and "no name construction on this path" cannot
both hold. v1 accepts the wasted snprintf and leaves the core alone; §3.2 says why that
is nearly free.

---

## 3. Design

### 3.1 Scope discipline

`ZnaWriter` buffers Python `str`. Do **not** build a columnar C++→codec fast path here;
that is a separate change with its own measurement burden. Budget: one C++ adapter, one
Python reference, one Python input strategy, a flag group, deletions.

### 3.2 The adapter — `merge_chunk_records`

Beside `merge_chunk` in `fastq_chunk.hpp`, same two raw buffers, same consumption
protocol:

```cpp
enum Slot : uint8_t { SLOT_MERGED = 0, SLOT_MATE1 = 1, SLOT_MATE2 = 2 };

/// One emitted record. `seq_*` index into `seqs`; `hdr_*` index into the ORIGINAL
/// input buffer -- buf1 for MERGED and MATE1, buf2 for MATE2, derivable from `slot`.
struct RecordEnd {
    uint32_t seq_off, seq_len;
    uint32_t hdr_off, hdr_len;   ///< hdr_len == 0 when want_headers is false
    uint8_t  slot;
};

/// Merge every whole pair available, appending RECORDS instead of FASTQ text.
inline void merge_chunk_records(const uint8_t* buf1, size_t n1, size_t& pos1,
                                const uint8_t* buf2, size_t n2, size_t& pos2,
                                const Params& p, bool check_sync, int64_t base_index,
                                bool want_headers,
                                ChunkScratch& sc, std::string& seqs,
                                std::vector<RecordEnd>& ends, ChunkStats& st);
```

**Design notes, each with a reason:**

- **No quality strings.** ZNA does not store them. This is a real part of the win.
- **Headers are offsets into the caller's input buffer, not copies**, and only when
  `want_headers`. This is not just cheaper than emitting a header blob — it is the only
  *correct* option for `MERGED`. `out.recs[0].h` for a merged record points at
  `Scratch::name`, a single per-pair buffer that the next pair overwrites; an adapter
  that recorded that span would hand Python a window onto the wrong record. Pointing at
  `r1.h` in the input buffer sidesteps it, and gives byte-identical label values to the
  two-step path (§3.5).
- **Slot, not flags, crosses the boundary.** Keeps ZNA's flag layout out of the C++ and
  keeps `_flags_from_ends` the single definition.
- **`ChunkStats` shared verbatim** with `merge_chunk`, so `--merge-json` can emit the
  same block. Watch the `ensure_bins` discipline (trap 4).
- **One shared inner loop.** Refactor `merge_chunk` onto a template/`Emit` functor
  rather than copying the loop — it carries the consumption protocol, the desync check
  and the statistics, and two copies will drift.

**The Python binding is a different shape from the C++ helper, and the plan must say
so** — `merge_chunk`'s binding takes **absolute slice ranges**, not a buffer and a
length, and returns seven values:

```python
merge_chunk_records(buf1, start1, end1, buf2, start2, end2,
                    match_q, step_q, t_merge_q, t_trim_q,
                    min_read_length, disagree_q, check_sync, base_index,
                    want_headers)
  -> (seqs, ends, consumed1, consumed2, counters, len_hist, olen_hist, insert_hist)
```

`ends` is a list of `(seq_off, seq_len, hdr_off, hdr_len, slot)` tuples — the same shape
the reference backend returns, so the differential test compares them element for
element. `hdr_off` is **absolute** into the caller's buffer (like `split_records`);
`consumed1`/`consumed2` are **relative** to `start` (like `merge_chunk`). Mirror both
conventions exactly.

**Consumed counts are RELATIVE to `start`; `split_records` returns an ABSOLUTE offset.**
The caller does `pos += c` for one and `pos = o` for the other. Reproduce that exactly;
copying the wrong convention skips records from the second chunk onward.

**Python side decodes once per chunk**, not once per record: `seqs.decode("ascii")` then
slice `str`. Sequences are ASCII so byte offsets are character offsets.

### 3.3 Behaviour changes to state loudly

| today (`zna merge` → `zna encode --interleaved`) | with `--merge-pairs` |
|---|---|
| pairing re-derived from read names; duplicate names silently mis-pair (§1) | pairing handed over; structurally impossible |
| `--treat-unpaired-as-merged` asserts every unpaired record is a full fragment | exact per record; **the flag is rejected** |
| full-overlap *pairs* detected by `_is_reverse_complement` on the sequences | the merge decision (quality-aware, evidence-scored) decides instead |

Rejecting `--merge-pairs --treat-unpaired-as-merged` must say *why* it is gone, not just
that it is. The flag itself stays — the two-step path is still supported, and
`test_merge_encode.py` and `test_cli.py` both exercise it.

**A third change, and it is an improvement that will still surprise someone.**
`--merge-pairs` replaces encode's two-file reader, and that reader does **no name
checking at all**: `_stream_paired_files` is `for s1, s2 in zip(p1, p2)`
([cli.py:743-752](../src/zna/cli.py#L743-L752)), so today two files that do not
correspond are silently mispaired and the longer one is silently truncated. The merge
kernel compares `base_name` per pair and raises `R1/R2 out of sync at pair N`. So
`--merge-pairs` turns two classes of silent corruption into a hard error — good, but it
means an input that "worked" before now fails. Say so in the flag's help, and point at
`--no-sync-check` for input whose names genuinely differ by design.

**One case rev 1 did not list:** `_full_fragment_flags` today stamps `IS_FULL_FRAGMENT`
on **both** mates of a pair whose mates are exact reverse complements
([cli.py:967-972](../src/zna/cli.py#L967-L972)), documented in the README as automatic
full-overlap detection. Under `--merge-pairs` such a pair merges instead (150 clean
bases score ~300 bits against a 28-bit threshold), so it becomes **one** full-fragment
record rather than two. That is a corpus change, it is the more correct one, and it
should be asserted rather than discovered.

### 3.4 The Python seam

Dispatch, at [cli.py:1266](../src/zna/cli.py#L1266) — note `merge_pairs` goes **first**:

```python
if merge_pairs:
    stream = stream_merge_pairs(args, label_defs, tag_map)   # 6- or 7-tuples
    carries_ends = True
elif label_defs:
    stream = stream_inputs_labeled(args, label_defs, tag_map)
    carries_ends = False
elif preserve_normalization:
    stream = stream_inputs(args, with_ends=True)
    carries_ends = True
else:
    stream = stream_inputs(args)
    carries_ends = False
```

`carries_ends` replaces `preserve_normalization` as the write loop's predicate — that is
rev 1's refactor and it is still right, for a better reason: the two things were never
the same, and `--merge-pairs` carries ends while explicitly *not* preserving
normalization. `ZnaWriter(preserve_normalization=...)` keeps its current value
(re-encode only).

**Why a 6-tuple when only one bit is new — and the cheaper alternative.** Since
`has_start` is always True here (§0.1), the minimum viable stream is a 5-tuple
`(seq, is_paired, is_read1, is_read2, is_full)`, which needs no `carries_ends` predicate
and no step 4 at all. The 6-tuple is chosen anyway because it makes the write loop
*identical* to the branch the re-encode path already uses, so one arm serves both
callers instead of two arms drifting apart. That is a judgement call, not a
requirement. But note the 5-tuple does **not** avoid the refactor — it collides with
`stream_inputs_labeled`, which already puts `labels` at index 4, and with
`ZnaWriter.write_records`, which reads `rec[4], rec[5]` as the ends whenever
`len(rec) > 5`. Two of the three existing consumers already fix this shape, so the
6-tuple is the cheaper option as well as the tidier one. The refactor is a
**prerequisite**, not a follow-up — see §4 step 3 for why it must land first.

**Step 4 has a regression surface on users who never pass `--merge-pairs`.** Today
`preserve_normalization = is_reencoding and input_header.strand_normalized`, so
re-encoding a **non**-normalized ZNA takes the 4-tuple path and re-derives
`IS_FULL_FRAGMENT` with `treat_unpaired_as_merged` defaulting False — i.e. it **drops the
flag**. Measured: `encode m.fq --treat-unpaired-as-merged` → flags `{16: 6}`; re-encoding
that `.zna` → `{0: 6}`. Gating on `carries_ends` instead would silently *fix* that. It is
a fix, but it is a behaviour change for existing invocations and it must be a deliberate,
tested one — not a side effect noticed later.

**Tuple shapes.** Labels go **last**, matching what `ZnaReader.records(with_ends=True)`
already yields on a labeled file (labels at index 6):

```
unlabeled + ends:  (seq, is_paired, is_read1, is_read2, has_start, has_end)
labeled   + ends:  (seq, is_paired, is_read1, is_read2, has_start, has_end, labels)
```

**The write loop** — the `carries_ends` branch must gain a labeled arm, or
`write_record` raises `ValueError: Expected N label values, got 0`
([core.py:363-367](../src/zna/core.py#L363-L367)) the first time a labeled stream
carries ends:

```python
for unit in _fragment_units(stream):                     # keep: mate grouping
    # No fragment-level N filter here any more: `drop` is gone, and `trim3` cuts per
    # record because trimming never orphans a mate. Under `--merge-pairs` the policy
    # has already run inside the kernel, so this path sees no no-calls at all.
    if carries_ends:
        for rec in unit:
            is_rc, is_full = _RC_FULL_BY_ENDS[(2 if rec[4] else 0) | (1 if rec[5] else 0)]
            writer.write_record(rec[0], rec[1], rec[2], rec[3],
                                labels=rec[6] if label_defs else None,
                                is_rc=is_rc, is_full_fragment=is_full)
    else:
        ...                                              # unchanged
```

Keep the single-end fast path at [cli.py:1283-1302](../src/zna/cli.py#L1283-L1302); it
exists to avoid a generator resume and a list allocation per read, and `--merge-pairs`
never reaches it (two files ⇒ `single_end_input` is False).

**Do not remove or rename `stream_inputs`** — `tests/test_cli.py` imports it at module
scope, so deleting the name turns a refactor into a collection-time `ImportError`
across the whole file.

**Read every new arg with `getattr(args, name, default)`.** 30 ad-hoc `class Args`
shells drive `encode_command` across 43 call sites and none of them will have the new
attributes.

### 3.5 Labels — **DECIDED: support them**

`--label-defs` is supported with `--merge-pairs`, and each emitted record's labels come
from **its own source header**: `MERGED` → R1's header, `MATE1` → R1's, `MATE2` → R2's.

Four reasons, in order of weight:

1. **The production pipeline already does merging and labels together.**
   [tests/test_merge_encode.py](../tests/test_merge_encode.py) runs `zna merge` →
   `zna encode --interleaved --treat-unpaired-as-merged --strand-normalize
   --label-defs`, and exists because a bug in exactly that combination once flagged
   every record unpaired. Rejecting the combination would delete a working capability
   precisely at the seam this change is supposed to simplify.
2. **The risk rev 1 worried about does not exist.** `--label-defs` is not a pattern
   system — it is a YAML file of column definitions, and extraction is fixed: split the
   header on whitespace, drop token 0, match `KEY` of each `KEY:T:VALUE` token against
   each label's tag ([cli.py:668-693](../src/zna/cli.py#L668-L693)). The appended
   ` merged_<n1>_<n2>` token contains no colon and is rejected by `colon1 < 1` — it is
   *structurally* unmatchable, and the docstring at
   [cli.py:658-660](../src/zna/cli.py#L658-L660) says so by name. R1's tags survive
   byte-identically ahead of it, so the one-step and two-step label values are equal by
   construction.
3. **The adapter cost is one `(offset, length)` pair per record**, into a buffer that
   already exists, gated on a bool. No copy, no allocation, no second loop.
4. **"Labels from R1's header" — rev 1's phrasing — would have been a silent behaviour
   change.** Labels today are per record from that record's own header: the two-file
   path extracts `l1` from R1 and `l2` from R2 independently
   ([cli.py:1054-1057](../src/zna/cli.py#L1054-L1057)). Applying R1's values to `MATE2`
   would overwrite R2's wherever the two mates' tags differ, and nothing would notice.

**State in the docs:** on a merged record the label values are R1's, and R2's are
discarded. For per-fragment tags (`ZI`/`ZJ`/`ZF`) the two are identical by construction;
for a genuinely per-mate tag they are not, and there is no way to keep both in one
record.

### 3.6 Threading

v1: merge inline, single-threaded, and measure. `_drive_threaded`'s ordering guarantee
is what keeps output byte-identical across thread counts; any later reuse must keep it.

Note the reader is threaded regardless: `_RawStream.__init__` starts a
`ThreadPoolExecutor(max_workers=1)` and a pigz child **per stream**. Inside a generator
those must be torn down on early exit — wrap the stream in `contextlib.closing` or drive
it under the `ExitStack` that `encode_command` already holds. A generator abandoned
mid-file otherwise leaks two threads and two pigz processes.

---

## 4. Build order

Do not batch. Each step is independently reviewable; the one ordering constraint that is
*not* negotiable is 3 before 4, for the reason given there.

**Rebuild command** (the plan must state it; a fresh session will not guess the flag):

```bash
/Users/mkiyer/sw/miniforge3/envs/zna_merge/bin/pip install -e . --no-build-isolation
```

`--no-build-isolation` is load-bearing. Without it, scikit-build-core reconfigures the
existing build tree against pip's ephemeral environment and rewrites `CMakeCache` to
point at a directory pip then deletes; alternating between the two forms costs a full
rebuild each way. Measured: ~1.5 s no-op, ~2.3 s after touching a header (CMake's
dependency scanning recompiles only `_merge_accel`), ~9 s on cache invalidation.
`merge_chunk_records` lives in `fastq_chunk.hpp`, already a header of
`src/zna/merge/_accel.cpp`, so **no `CMakeLists.txt` change is needed**.

1. **`merge_chunk_records` in `fastq_chunk.hpp` + binding.** Do **not** add it to
   `backend._REQUIRED_FUNCTIONS` yet. The failure mode is worse than an error: only
   `_load()` validates the set ([backend.py:89-95](../src/zna/merge/backend.py#L89-L95)),
   while `available_merge_backends()` merely *imports*
   ([backend.py:49-57](../src/zna/merge/backend.py#L49-L57)) — so against a stale `.so`,
   `accel` still shows as available while `get_merge_backend()` catches the ImportError
   and **silently falls through to the reference backend**. That is the 50x slowdown
   `args.py`'s `--backend` help says the design exists to prevent: on a cluster it is
   indistinguishable from a slow node. Add the name in step 2, after both backends have
   it, and rebuild. *Test:* `test_record_adapter_agrees_with_merge_chunk` — same input
   through both adapters; sequences equal record for record, and the slot recovered from
   `merge_chunk`'s FASTQ (single ⇒ MERGED, `/1` ⇒ MATE1, `/2` ⇒ MATE2) equals the slot
   the new one reports. Consumed counts and all 13 counters equal.
2. **The reference implementation in `_pymerge.py`.** *Now* add the name to
   `_REQUIRED_FUNCTIONS`. *Test:* extend `_chunk_backends()` to
   `_chunk_record_backends()` and add `test_record_chunks_agree_across_backends` —
   seqs blob, ends list, consumed counts, counters, all three histograms.
3. **The `carries_ends` refactor, landed FIRST and as a pure no-op.** Set
   `carries_ends = preserve_normalization` and add the labeled arm (§3.4). Nothing
   changes behaviourally — the two are still equal — so this is a refactor that can be
   reviewed and reverted on its own. *Test:* the existing suite, unchanged, plus one
   test that a labeled ends-carrying stream reaches `write_record` with its labels
   (constructible today by hand-feeding the write loop; it is unreachable via the CLI
   until step 4).

   **This must precede step 4, not follow it.** Landing `stream_merge_pairs` while the
   write loop still branches on `preserve_normalization` would send merge-pairs records
   down the `else` arm, where `_full_fragment_flags` re-derives the flag and `rec[4]` —
   `has_start` — is read as `labels`. Step 4 alone is not independently landable; this
   is the ordering that makes it so.
4. **`stream_merge_pairs`** in `zna/merge/` (not `zna/cli.py` — it needs `_RawStream`,
   `MergeParams` and `DISAGREE_Q`, none of which belong in encode's startup), plus the
   `--merge-pairs` flag group and the rejections of §5. Flip `carries_ends` to
   `preserve_normalization or merge_pairs`.
5. **Statistics:** `--merge-json`.
6. **Docs:** README `zna encode`, CHANGELOG, mark this plan executed.

Steps 1 and 2 are independent of 3–4 and can land in either order relative to them.

---

## 5. CLI surface

**Added to `zna encode`** (mirroring `zna merge`; rev 1 never enumerated these):

| flag | dest | default | note |
|---|---|---|---|
| `--merge-pairs` | `merge_pairs` | False | requires exactly two input files |
| `--threshold-merge` | `t_merge` | 28.0 | with `--threshold-trim` and `--min-read-length`, the entire `MergeParams` constructor |
| `--threshold-trim` | `t_trim` | 8.0 | |
| `--min-read-length` | `min_read_length` | 40 | |
| `--no-sync-check` | `no_sync_check` | False | |
| `--merge-json` | `merge_json` | None | named to avoid confusion with `zna inspect --json` |
| `--merge-backend` | `merge_backend` | `auto` | **not** `--backend`; and it must not call `backend.use()`, which mutates process-global state for every other consumer. Pass the chosen backend down explicitly. |
| `--allow-empty` | `allow_empty` | False | port it. Without it a zero-pair input writes a 0-record `.zna` and a library vanishes from the corpus with every stage green. |

Do **not** reuse `add_merge_arguments`: it declares `--in1`/`--in2`/`--out` as
`required=True` and adds its own `-q/--quiet`, which collides with encode's.
Factor the scoring knobs into a shared `add_merge_algorithm_arguments(p)` so the two
parsers cannot drift.

`--chunk-size` is **not** mirrored: its help says "read pairs per work unit" but
`_drive_serial` never uses it that way, and v1 is single-threaded.

**Rejected combinations**, each with a message that says why:

| combination | reason |
|---|---|
| `--merge-pairs` with ≠ 2 files | it is a two-file mode |
| `--merge-pairs --interleaved` | reachable today, but the existing message at [cli.py:1088](../src/zna/cli.py#L1088) names the *wrong* flag ("Cannot use --interleaved with 2 input files"). Add a specific one. |
| `--merge-pairs --seq-len-bytes 1` | a merged record is up to 2× a read; see trap 10 |
| `--merge-pairs --treat-unpaired-as-merged` | the geometry is exact; the assertion is what this deletes |
| `--merge-pairs` on ZNA input | re-encode has its own geometry source |
| `--merge-pairs --fasta` | no qualities, so no consensus; and labels need FASTQ |
| `--merge-pairs -` / stdin | the merge reader takes **paths** and drives pigz itself |

`--strand-normalize`, `--strand-specific`, `--shuffle`, `--npolicy`, `--block-size`,
`--label`/`--label-defs`/`--label-desc` all remain **supported** and compose (§0.2).

Argparse note: `files` is `nargs="*"`, so "exactly two" is a hand-written check next to
[cli.py:1085](../src/zna/cli.py#L1085). Also, argparse requires `nargs="*"` positionals
to be contiguous, so `zna encode --merge-pairs R1.fq --quiet R2.fq` fails at the
*top-level* parser with `unrecognized arguments: R2.fq`. Worth one line in the help.

**Error surface.** `zna merge` gets clean messages because `run()` wraps the drive loop
in `except InputError as e: raise SystemExit(str(e))`. `encode_command` must do the
same, and it must run the "both buffers drained" check **before** the output file is
finalized — in `zna merge` today that check sits outside the `with FastqWriter(...)`
block, so a desynced input leaves a complete, plausible-looking output on disk before
exiting non-zero. Do not inherit that ordering.

---

## 6. Tests

In the order they should be written.

1. **`test_merge_pairs_matches_the_two_step_pipeline`** — the headline. Same input
   through (a) `zna merge` → `zna encode --interleaved --treat-unpaired-as-merged` and
   (b) `zna encode --merge-pairs`; compare decoded records **and flags**. Expect
   **exact equality** on well-formed input (§0.4: the blanket flag is not wrong on merge
   output). Compare with `--strand-normalize` off, or with `--strand-specific`, because
   unstranded normalization is a coin flip — the pure-Python codec uses the *unseeded
   global* `random`, so two runs of the same path disagree.
2. **`test_merge_pairs_survives_duplicate_read_names`** — the improvement, pinned. The
   fixture of §1: two merging fragments sharing a read name. The two-step path returns
   two `IS_PAIRED` records with `IS_FULL_FRAGMENT` clear; `--merge-pairs` returns two
   unpaired full-fragment records. Assert both sides, so the difference is a *pinned
   improvement* rather than a tolerated drift.
3. **`test_merge_pairs_geometry_is_exact`** — pairs with known inserts across all five
   regimes (no-overlap / boundary / partial / full / read-through), asserting
   `IS_RC` and `IS_FULL_FRAGMENT` per record against truth. Re-derive the regimes
   locally: `scripts/merge_bench/simulate.py` needs numpy and an indexed FASTA, has no
   `__init__.py`, keeps its classification inline in `main()`, and says in its own
   docstring that it is not part of the suite.
4. **`test_merge_pairs_composes_with_strand_normalize`** — `--merge-pairs
   --label-defs --strand-normalize --strand-specific --read2-antisense`: merged records
   `(True, True)`, mates `has_start != has_end`, every reported edge a real terminus.
   This is the test that would have caught §0.2.
5. **`test_merge_pairs_labels_come_from_each_records_own_header`** — give R1 and R2
   *different* tag values and assert `MATE2` gets R2's. This is the test that would have
   caught §0.4/rev-1-§3.5. Plus: merged-record labels equal the two-step path's.
6. **`test_merge_pairs_rejects_treat_unpaired_as_merged`** — with the message.
7. **`test_merge_pairs_thresholds_reach_the_kernel`** — the cheapest high-value test in
   the feature, and rev 1 had nothing like it. Encode the same input at
   `--threshold-merge 28` and at `60`; assert the merged-record count changes. A
   parameter that is parsed and then never reaches `MergeParams` is invisible to every
   other test here.
8. **`test_merge_pairs_mate_grouping`** — a pair emitted by the kernel is grouped as one
   fragment by `_fragment_units`. **Note what this actually tests.** Rev 1 called it "C3
   at the encode seam" and rev 2 called it `pair_atomicity` and drove it with
   `--npolicy drop`; both are now wrong. `drop` is gone, and `zna merge` *does* have an
   N policy (rev 2 said it did not) — but that policy is `trim3`, which never orphans a
   mate, so there is no atomicity left to test at this seam. C3 is the merger's own
   all-or-nothing `min_read_length` rule at
   [merge_core.hpp:496-502](../src/zna/merge/merge_core.hpp#L496-L502), and it is
   already covered. What this test pins is that the new strategy sets
   `is_paired`/`is_read1`/`is_read2` well enough for `_fragment_units` to group mates —
   which is the real risk, and which a write loop that hoisted `carries_ends` above
   `_fragment_units` would break.
9. **`test_merge_pairs_requires_two_files`**, and the rest of §5's rejection table.
10. **Cross-backend equality** for the new adapter (build step 2).
11. **Mutation test** (rev 1 trap 6): swap `has_start`/`has_end` in the adapter and
    confirm test 3 fails. Note `write_record` validates almost nothing — not
    `is_read1 && is_read2`, not `is_paired` consistency — so the tests are the only
    guard.

Home: [tests/test_merge_encode.py](../tests/test_merge_encode.py) already has the
`encode()` subprocess helper and `read_back()`. `read_back` unpacks 7 values and so
requires a labeled file; add an unlabeled variant rather than forcing labels on every
new test.

Run all three configurations, and know what each actually proves:

| configuration | how | what it proves |
|---|---|---|
| compiled | `envs/zna_merge`, **with its `bin/` prepended to PATH** or pigz is silently skipped | the shipping path |
| reference | *not* a suite-wide switch — there is no `conftest.py` and no pytest option. It exists only via the `any_backend` fixture and the explicit `_chunk_backends()` pairs; every CLI-level test calls `_backend.use("auto")` and so always runs compiled | the oracle, on the tests that route through it |
| extension-less | `envs/zna` — has `zna._accel` but **not** `zna.merge._accel` | the reference merge kernel against the *accelerated* label extractor |

Cross-backend tests **skip** where the merge extension is absent, so a green run in
`envs/zna` proves less than it looks — check the skip count, not just the exit code.
Baseline: 576 passed / 33 subtests in `zna_merge`, 538 passed / 38 skipped in `zna`.

If a `--merge-backend python` suite run is wanted as a real configuration, it needs a
`conftest.py` option; that is new work, not an existing capability.

---

## 7. Verification beyond the suite

The **inputs** exist at `~/proj/zna_merge_bench/` (`refs/hg38.fa`, `run1M/sim_1.fq.gz`,
`run1M/sim_2.fq.gz`, `run1M/sim.truth.tsv`) — 1M pairs, already built. The **comparator
does not**: `scripts/merge_bench/compare.py` never reads `.zna`; it runs the tools
itself and scores merged FASTQ against the truth sidecar. A `.zna`-level comparator
(decode both files with `records(with_ends=True)`, join to truth, score
`IS_FULL_FRAGMENT`) is **new tooling** — budget it, and put it beside `compare.py`.

**Use absolute interpreter paths.** The ambient `PATH` here begins with
`envs/zna/bin`, so a bare `zna merge` runs the extension-less py3.14 install and is
refused by design, and a bare `cmake` is the wrong env's — which is what silently
invalidates the build cache.

```bash
Z=/Users/mkiyer/sw/miniforge3/envs/zna_merge/bin
PATH=$Z:$PATH $Z/zna merge --in1 run1M/sim_1.fq.gz --in2 run1M/sim_2.fq.gz --out m.fq.gz
PATH=$Z:$PATH $Z/zna encode --interleaved --treat-unpaired-as-merged m.fq.gz -o two_step.zna
PATH=$Z:$PATH $Z/zna encode --merge-pairs run1M/sim_1.fq.gz run1M/sim_2.fq.gz -o one_step.zna
```

Then compare both `.zna` files record by record against `sim.truth.tsv`. Every
`IS_FULL_FRAGMENT` claim is checkable against the known fragment length — a stronger
statement than the merge benchmark could make, because it tests the *encoded flag*, not
the FASTQ that preceded it.

Report: records where the two paths disagree and which one truth agrees with; and,
separately, the residual `IS_FULL_FRAGMENT` records that are **not** the true fragment.
Expect ~5,591 per million from the wrong-shift merges of
[MERGE_BENCHMARK_RESULTS.md](MERGE_BENCHMARK_RESULTS.md) §2 — **both** paths carry
those, they are a property of the genome, and `--merge-pairs` neither creates nor fixes
them. Say so, so the number is not read as a regression.

Also report wall time and peak RSS for both paths: the speed claim is the other half of
the case and it has not been measured yet.

---

## 8. Traps

1. **The consumption protocol is not optional, and it has six parts, not one.** Rev 1
   named only the last. From `_drive_serial`
   ([merge/cli.py:201-218](../src/zna/merge/cli.py#L201-L218)):
   1. `fill(target)` **both** streams before each call, `target = max(256 KiB,
      chunk_size * 1024)` — a byte count, not a record count.
   2. Stop if either stream has zero available bytes.
   3. Call the adapter, then **stop if neither stream consumed anything.** Without this,
      a record larger than the refill target spins forever.
   4. Advance by the **relative** consumed counts (`pos += c`), not by assignment.
   5. Close both `_RawStream`s in a `finally` — each owns a prefetch thread and a pigz
      child, and `close()` is what surfaces a pigz error.
   6. Assert both buffers drained at EOF — *before* finalizing the output (§5).
2. **Relative vs absolute offsets.** `merge_chunk` consumed counts are relative;
   `split_records` returns absolute.
3. **`_REQUIRED_FUNCTIONS` is validated at import.** Register the new name only once
   both backends have it (§4).
4. **`ChunkStats` histograms grow with the arena.** Use the same `ensure_bins`
   discipline or long reads silently re-cap.
5. **`Scratch::name` heap overflow — verified, live, in shipped compiled code.**
   [merge_core.hpp:225](../src/zna/merge/merge_core.hpp#L225) sizes `name` as `cap + 64`
   where `cap` is the **read-length** arena (minimum 1024, doubling on read length), but
   [:461-466](../src/zna/merge/merge_core.hpp#L461-L466) writes `idlen + (h.n - cut)`
   bytes plus up to 64 more — all derived from the **header**. Any header longer than
   `arena - 64` writes past the buffer.

   Reproduced against the shipped `zna merge`, 51 bp reads (arena 1024, name buffer
   1088) and a 16 KB FASTQ header: one run aborted with `SIGABRT` from malloc's heap
   check, a second run of the same input survived and returned well-formed output — the
   nondeterministic signature of a heap overflow, which means the quiet runs are
   corrupting adjacent heap rather than being safe. Diagnosis confirmed by holding the
   header at 16 KB and raising the read length to 20,000: the arena grows to 32,768, the
   name buffer to 32,832, and the identical header is then handled cleanly.

   **The reference backend is unaffected**, so no cross-backend test can catch this
   without a long-header fixture. Fix: size `name` from `r1.h.n + 64`, or bound the copy.
   Not caused by this work, but `merge_chunk_records` calls the same `process_pair`, so
   it inherits it — fix it in step 1 and add the fixture.
6. **Do not re-derive `_flags_from_ends`.** Two tables derive from it; `ENDS_BY_FLAG`
   is a third, independent statement of the same rule (§0.4).
7. **The two label extractors disagree on malformed values.** The compiled one produces
   silent garbage where the reference raises `ValueError` (e.g. a Casava comment
   `1:N:0:ATCACG` read as tag `1`). Pre-existing; it means a cross-environment
   comparison of label values is not automatically apples to apples.
8. **A pair can emit zero records** (both mates, or a merged read, below
   `--min-read-length`; counted as `frags_short`). The adapter's `for i < r.n_recs` loop
   handles it naturally. Rev 2 added "distinguish dropped-by-the-merger from
   dropped-by-`--npolicy`" — under `--merge-pairs` there is now only one dropper, since
   the policy runs inside the kernel and its losses are already in `npolicy_bases`.
9. ~~**`scripts/wheel_smoketest.py` never imports `zna.merge._accel`**~~ — **done in
   0.4.0.** It now asserts both extensions load *and* that neither target overwrote the
   other, so a wheel missing the merge kernel fails before publication.
10. **A merged record can be twice a read long, and encode's only length check is a
    *maximum*.** There is no minimum-length filter anywhere in encode — `--min-read-length`
    is the sole floor on this path, which is another reason it must be mirrored. The
    ceiling is the hazard: `write_record` raises `ValueError` when `seq_len >
    (1 << 8*seq_len_bytes) - 1` ([core.py:356-360](../src/zna/core.py#L356-L360)). With
    `--seq-len-bytes 1` (max 255) and 2×150 input, a merged read of 256–300 bases raises
    **mid-stream**, leaving a structurally valid partial `.zna` on disk. Today every
    record on the two-file path is at most one read long, so this cannot happen. Either
    reject `--merge-pairs --seq-len-bytes 1`, or validate `2 * max_read_len` up front.

### 8.1 Prerequisite — **DONE in 0.4.0**

This was the one thing that had to be fixed before `--merge-pairs`, and it has been.
Recorded here because the reasoning still governs how the adapter must hand records over.

`_flags_from_ends` claimed to be "a bijection over the three reachable states". There
are **four**: flag byte 24 (`IS_RC | IS_FULL_FRAGMENT`) is what the encoder writes
whenever strand normalization reverse-complements a merged read. `ENDS_BY_FLAG` maps 16
and 24 alike to `(True, True)` — correctly, since both records do have two real fragment
ends — so no inverse exists, and the one that existed cleared `IS_RC` on every
full-fragment record that passed through it. `--shuffle` was the live vector and it is in
the production recipe, so the loss was reaching the corpus.

**Fixed** by carrying the flag byte instead of a view of it: `ZnaReader.copy_records()`
-> `ZnaWriter.write_copy()`, with `_flags_from_ends`, `_RC_FULL_BY_ENDS` and
`_FLAG_BITS_BY_ENDS` deleted and `write_records()` refusing the 6-tuple shape. Verified
on the 1M-pair library through the full recipe: 293,853 records at byte 24 survive
`--strand-normalize --shuffle`, and `--restore-strand` reproduces the original base
multiset exactly over all 1,416,630 records.

**What this means for `--merge-pairs`:** the law is now written into `core.py` and the
adapter must obey it. *A view is for reading; the flag byte is for copying.* The
merge-pairs stream is **producing**, not copying, so it hands over `is_full_fragment` and
lets the writer own `IS_RC` — which is §0.1 and §0.2, and is why this feature needs no
`preserve_normalization` and composes with `--strand-normalize`.

---

## 9. What this deletes

- the FASTQ intermediate, its pigz write, its gzip read, its full re-parse
- quality-string formatting and parsing on the merge path
- name-based re-pairing at the encode seam, and the duplicate-name corruption of §1
- `--treat-unpaired-as-merged` *at this seam* (the flag itself stays for the two-step
  path, which is still supported)
- the `preserve_normalization` / "stream carries ends" conflation in `encode_command`

If the diff is net-positive by much, the design drifted — most likely into §3.1's
columnar fast path or §3.6's threading.

---

## 10. Open questions for the corpus owner

1. **`--merge-json`'s `"tool"` field.** `_assemble_stats` hardcodes `"tool":
   "zna-merge"`. A block emitted by `zna encode --merge-pairs` naming a tool that did
   not run defeats the provenance the field exists for — but `hulkrna`'s
   `gather/tools/read_merge.py` consumes these blocks and may key on it. Recommend
   keeping every existing key, setting `"tool": "zna-encode --merge-pairs"`, and
   checking `read_merge.py` before landing.
2. **Counting.** `--merge-json` counts records the *merger* emitted. Rev 2's reason for
   the gap (`--npolicy drop` dropping fragments afterwards) is gone with `drop`, and on
   the one-step path the policy runs inside the kernel, so `emitted_records` *should*
   equal the records in the `.zna`. Assert that rather than documenting a difference —
   it is now a real invariant, and a cheap one to check.
3. **Two `.zna` inputs write a 0-record file at exit 0.** Verified: `zna encode a.zna
   b.zna` takes the two-file branch (the ZNA check is `len(files) == 1`), hands binary
   to the FASTQ parser, and succeeds with 0 records. It does warn — but only about
   *format inference*, not about the real problem. Pre-existing and unrelated to this
   change; found while enumerating the reject matrix, and worth a one-line guard beside
   the new ones.
4. **No provenance in the container**, at the *file* level. Still true and still worth
   deciding before the rebuild: the ZNA header has no field distinguishing a
   `--merge-pairs` corpus from a two-step one, so telling them apart after the fact
   needs a header field or a `--description` convention.

   *Per-record* provenance is no longer missing — 0.4.0 added the `ZN:i:<bits>` tag and
   its optional `C` column (`docs/NPOLICY_PLAN.md` D5). `--merge-pairs` must write the
   same `PROV_*` bits directly off `PairResult`, with no tag round-trip, and a test
   should assert the one-step and two-step provenance columns agree — which is the
   sharpest available version of the §7 comparison, since it comes per record rather
   than per run.
