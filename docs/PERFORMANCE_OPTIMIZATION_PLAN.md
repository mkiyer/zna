# Performance optimization plan

Status: **proposed**. Written 2026-08-11 against `ab05813` (0.3.3).
Nothing here is applied. The 0.3.3 tree is the correct, tested baseline these
changes would build on.

Produced by a seven-layer measurement sweep (writer, reader, C++ encode, C++
decode, FASTQ parsing, compression/IO, shuffle/CLI). 31 candidate findings were
raised; the 29 claiming ≥5% went to independent verification, where a second
engineer wrote their own benchmark from scratch, applied the change in a scratch
copy, and ran the full suite. **19 reproduced with the suite green, 7 were
rejected, 4 were lost to infrastructure errors and remain unverified.**

Every number below is the **verifier's** figure, not the finder's. The two
differed constantly and almost always in the same direction:

| finder claimed | verifier measured |
|---|---|
| 348.8% | 274% |
| 260% | 178% |
| 162% | 27.6% |
| 129% | 122.4% |
| 128% | 38.7% *(and incorrect — see §1)* |
| 36% | 5.3% *(rejected)* |

---

## 1. Read this first: the suite does not catch codec corruption

Two verifiers independently reported a change that **silently corrupted data
while all 282 tests passed**. One wrote it plainly:

> the full 282-test suite passes anyway, so "tests green" is worthless evidence
> here.

The NEON 2-bit packing kernel was faster, plausible, and produced wrong bases.
The `encode --shuffle` fusion sketch "silently corrupts data on three real paths"
and also passed.

Our tests check round-trips through *fixed* inputs and API contracts. They do not
assert that arbitrary sequence survives arbitrary codec configuration bit-exactly,
which is precisely the property a codec optimization threatens.

**Prerequisite for every C++ change below:** a round-trip fuzz harness asserting
bit-exact recovery across random sequences × lengths (including non-multiples of
4, 0, and 1) × `seq_len_bytes` ∈ {1,2,4} × N-policies × strand configurations ×
compression on/off, run against both backends. Until that exists, a green suite
is not evidence that a codec patch is correct.

## 2. Verified findings

Ordered by verified speedup **on the path changed**. The "share" column is what
matters for prioritisation — several large wins sit on small slices of real work.

### Tier A — C++ decoder (one root cause, biggest wins)

The decoder allocates a `std::string` per record, nanobind copies it into a
`PyUnicode`, and Python then unpacks and repacks the tuple. Three copies of every
base, plus a tuple with four bool slots, per record.

| verified | change | file |
|---|---|---|
| **+178%** | `decode_block_columnar`: decode straight into `PyUnicode_New` + a 256-entry LUT emitting 4 bases per store, and GC-untrack the record tuples | `_accel.cpp` |
| **+168%** | `decode_block_labeled`: the same, plus a direct `nb::list` | `_accel.cpp` |
| **+78.6%** | decode into `PyUnicode` rather than into a `std::string` that is then copied | `_accel.cpp` |
| **+29%** | GC-untrack emitted record tuples — they can never be in a cycle | `_accel.cpp` |
| **+28.2%** | have the C++ decoder emit the caller's tuple width directly, killing the Python rebuild | `_accel.cpp` |
| **+22%** | return a hand-built `PyList` of `PyTuple` instead of `std::vector<std::tuple<...>>` | `_accel.cpp` |

These overlap heavily — they are facets of one rewrite, not six independent wins.
Note the Python backend has had the 4-bases-at-a-time `_DECODE_TABLE` since the
beginning; the C++ side never did.

**Worth knowing:** one verifier built a decoder that keeps the existing 5-tuple
shape but constructs it with the raw C API, and measured **1.53x on the C++ path
with zero new API surface**. That is roughly 35% of the total available win,
available without any interface change.

### Tier B — C++ encoder and reverse-complement

| verified | change | file |
|---|---|---|
| **+122.4%** | `encode_block` spends about half its time in nanobind's `std::vector<std::string>` conversion; take the sequences without the round-trip | `_accel.cpp` |
| **+87%** | reverse-complement allocates a `std::string` per record inside the hot loop | `_accel.cpp` |
| **+63.5%** | do the `restore_strand` reverse-complement inside the decoder instead of per-record in Python | `_accel.cpp` + `core.py` |

### Tier C — Python and CLI (no C++, low risk)

Four of these are regressions introduced by the 0.3.3 correctness work. They are
the cheapest wins available and should go first.

| verified | change | file |
|---|---|---|
| **+69.4%** | compress blocks on a worker thread pool — zstd releases the GIL | `core.py` |
| **+50.2%** | fuse the shuffle into encode instead of writing the file then reading it three more times. **The sketch as written corrupts data on three paths — see §1. A corrected version was measured equally fast.** | `cli.py`, `_shuffle.py` |
| **+48%** | bytes-native interleaved pairing: stop decoding and re-splitting the read name per record | `cli.py` |
| **+34.1%** | single-end input pays per-record fragment grouping it can never use | `cli.py` ← *0.3.3 regression* |
| **+27.6%** | the default N-drop filter allocates a 150-byte uppercase copy of **every** read; test membership without allocating | `cli.py` ← *pre-existing, preserved in 0.3.3* |
| **+26.8%** | replace hand-rolled newline trimming in the three FASTQ parsers | `cli.py` |
| **+19.7%** | bucket intermediates are compressed at the source file's zstd level, then immediately re-read and deleted | `_shuffle.py` |
| **+18.3%** | gzip FASTQ input goes through `GzipFile.readline` four times per record | `cli.py` |
| **+12.1%** | `_flags_from_ends` is a Python function call on every record in every write loop | `core.py`, `_shuffle.py` ← *0.3.3 regression* |

### Tier D — new API

| verified | change | file |
|---|---|---|
| **+274%** | a columnar `blocks()` batch API yielding `(list_of_sequences, flags_bytes)` per block | `core.py` + `_accel.cpp` |

The largest single number, and the most caveated. The verifier's notes:

- the 89%-of-loader framing assumed a loader that does nothing but append strings;
  with a realistic per-batch step the end-to-end win is **2.78x, not 3.36x**
- `zna encode` is **0% affected**; `zna decode` gains 1.48x and only if `cli.py` is
  rewritten to decode flag bits itself
- consumers needing per-record flags lose a third of the win
- it decays with read length: 3.74x at 150 bp, **1.69x** at 1–20 kb long reads
- **the incremental value over a fixed `decode_block_columnar` is ~2.6x, not 4.5x**
- the sketch raises on labeled files and cannot serve `restore_strand`
- `PyUnicode_New`/`PyList_SET_ITEM` are not in the limited API, and
  `CMakeLists.txt` passes `STABLE_ABI` to `nanobind_add_module`. It does not
  engage today, so this builds — but the change silently depends on that

This is the right shape for khorana's training loader specifically. Design it
**with that loader in hand**, after Tier A lands, since Tier A is most of the win.

## 3. Rejected

| claim | why |
|---|---|
| NEON 2-bit packing (+128%) | kernel win is real (+38.7% on today's code, not 128%) but **the code as written silently corrupts data** |
| `DEFAULT_BLOCK_SIZE` 4→16 MiB (+36%) | decode win is **5.3%**, an artifact of the harness holding its source dataset alive; the stated mechanism is refuted. The *compression-ratio* half does reproduce |
| bulk `write_records_full()` (+16%) | not reproduced |
| `write_records_labeled()` (+34.8%) | not reproduced |
| share one `encode_core` (+400%) | refuted as a performance claim; **confirmed as a correct, behaviour-preserving refactor** — still worth doing for drift, since the two entry points are hand-copied duplicates |
| N-policy in the vector path (+96%) | mechanism real, magnitudes do not reproduce, one stated risk inverted |
| bucket shuffle no-op round-trip (+44.5%) | reproduced but overstated by ~6 points; two supporting claims wrong |

## 4. Unverified

Four verifiers died on API errors. **These are not evidence of anything** and must
be re-measured before use:

- `reverse_complement` str→str with no intermediate `std::string` (claimed 551%)
- shuffle pass-0 counting scan reads block headers only (claimed 95.9%)
- gzip input via an external `gunzip`/`pigz` subprocess (claimed 60%)
- the decode output loop builds two f-strings per record

## 5. Suggested order

1. **The fuzz harness (§1).** Everything else depends on it.
2. **Tier C Python fixes**, starting with the four 0.3.3 regressions. No C++, no
   API change, no rebuild.
3. **Tier A/B C++ decoder + encoder**, as one coherent rewrite rather than six
   patches, since they share a root cause. Ship as 0.3.4.
4. **Re-measure the four unverified items.**
5. **Tier D `blocks()`**, designed against khorana's loader, once Tier A has taken
   most of the win out of it.

## 6. Measurement notes

This machine is noisy — identical code varied 40% between runs on file-based
benchmarks. Everything here used in-memory `BytesIO`, min-of-N with N ≥ 7, and
baseline/candidate interleaved in one process so drift cancels. Reproduce that
discipline or the numbers are meaningless.

Two traps that cost real time during the sweep:

- **The editable install's meta-path finder overrides `PYTHONPATH`.** Pointing
  `sys.path` at a different source tree does *not* load it; you get the working
  tree and a silently invalid A/B comparison. Use `python -S` with an explicit
  path, or swap files in place and restore.
- **Comparing two builds needs distinct `NB_DOMAIN`/module names**, or the second
  import silently returns the first.
