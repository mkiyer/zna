# Porting the read-merge tool from hulkrna into zna

Status: **executed 2026-08-12**, in zna 0.3.6. Written 2026-08-12 against hulkrna at
the commit that carries the finished merge tool and zna v0.3.5. What actually landed,
and the two places it deviated from this plan, are recorded in §14 at the end; the rest
of the document is left as written, because §4 (the contract) and §12 (the traps) are
the parts that stay true and are expensive to rediscover.

Companion reading, now all in this repo's `docs/`:
- [`READ_MERGE_REDESIGN.md`](READ_MERGE_REDESIGN.md) — why the scoring rule is what it
  is, with its three implementation decisions recorded in the header. **The design is
  settled; do not re-derive it.**
- [`MERGE_TOOL_AUDIT.md`](MERGE_TOOL_AUDIT.md) — a six-lens audit, every claim
  adversarially verified. The measurements quoted throughout this document come from
  there or from the session that followed it. §(c) is the list of ideas that look good
  and were measured and rejected — check it before proposing an improvement.
- [`READ_MERGE_ROADMAP.md`](READ_MERGE_ROADMAP.md) — the status board: what was fixed,
  what was deliberately left, what is next.
- [`SOS_EOS_ENCODING_PLAN.md`](SOS_EOS_ENCODING_PLAN.md) — the cross-repo
  fragment-boundary contract the merge tool exists to satisfy.

Still in `/Users/mkiyer/proj/hulkrna/docs/`: `ZNA_PIPELINE.md`, the pipeline work that
surrounds this tool and stays on that side.

---

## 1. What this tool is, in one page

It replaces fastp's paired-end merge step. Input: two positionally-synced FASTQ files
(R1, R2) from `samtools fastq`. Output: **one mixed interleaved FASTQ** — merged reads as
single records with the `/1`,`/2` suffix stripped, unmerged pairs as adjacent `/1`,`/2`
records — which `zna encode --interleaved` consumes.

Per pair it slides R1 against `revcomp(R2)` over the single axis of candidate fragment
lengths and scores every shift as a log-likelihood ratio in **bits** against chance
alignment:

```
score(s) = matches × log2((1−e)/0.25) + mismatches × log2(e/0.75)
         = matches × 1.9855 − mismatches × 6.2288          at e = 0.01
```

It takes the **argmax** over shifts (not fastp's first-accept) and reads that one score
at two thresholds:

| | condition | action |
|---|---|---|
| merge | `score ≥ t_merge` (28) | emit one full-fragment record |
| trim | `t_trim ≤ score < t_merge` (8) | keep both; cut the redundant overlap off R2's **3'** end |
| keep | `score < t_trim` | keep both, untouched |

Three parameters total: `t_merge`, `t_trim`, `min_read_length`. Everything else fell out
of the model. Measured against the fastp-derived rule it replaced: spurious detection on
unrelated pairs **5.17% → 0.26%**, zero spurious merges, merge rate **85.7% → 88.8%**.

**Why it belongs in zna.** Every historical bug in this chain lived at the repo boundary:
the merger computes the fragment geometry, and a *separate process* then re-infers it
from `/1`,`/2` header suffixes plus a global `--treat-unpaired-as-merged` declaration.
Inside zna that inference disappears (§8).

---

## 2. File manifest

Everything that moves. Line counts are as of the port commit.

| Source (hulkrna) | Lines | Purpose |
|---|---|---|
| `lib/hulkrna/merge/__init__.py` | 31 | package docstring + public exports |
| `lib/hulkrna/merge/overlap.py` | 227 | the `@njit` scoring kernel; `find_overlap`, `score_weights`, `threshold_bits`, `reverse_complement` |
| `lib/hulkrna/merge/pairs.py` | 295 | `MergeParams`, `process_pair`, `_build_merged`, the posterior consensus, name helpers |
| `lib/hulkrna/merge/cli.py` | 487 | argparse, JSON stats, single + parallel run paths |
| `lib/hulkrna/merge/fastqio.py` | 139 | pigz/gzip paired reader, single-stream writer, `InputError` |
| `lib/hulkrna/merge/__main__.py` | 3 | `python -m hulkrna.merge` |
| `scripts/hulkrna-merge` | 20 | CLI shim that puts `lib/` on `sys.path` — **does not move**, see §7 |
| `tests/test_read_merge.py` | 1189 | 173 tests |
| `tests/test_merge_zna_e2e.py` | 250 | 6 tests, the merge→zna seam |
| `workflow/envs/read_merge.yaml` | 8 | conda env — **does not move**, see §7 |

**~1180 lines of library code, ~1440 lines of tests.**

### 2.1 Dependencies — the tool is almost free-standing

Across all six modules the imports are: `argparse`, `dataclasses`, `gzip`, `json`,
`logging`, `math`, `platform`, `shutil`, `subprocess`, `sys`, `time`,
`concurrent.futures`, `multiprocessing` — **all stdlib** — plus:

- **`numba`**, imported in `overlap.py` inside a `try/except ImportError` with a
  no-op `njit` fallback. The pure-Python path is correct and ~50x slower.
- **`pigz`**, an external binary, used opportunistically via `shutil.which`; falls back
  to stdlib `gzip`.

**There is no `numpy` dependency.** `workflow/envs/read_merge.yaml` lists numpy, but
nothing imports it — do not carry that into zna's requirements.

### 2.2 The only coupling to hulkrna

`cli.py` does `import hulkrna; hulkrna.__version__` (one line, in `_assemble_stats`, for
the `tool_version` provenance field). Plus cosmetics: the logger name `hulkrna.merge`,
`prog="hulkrna-merge"`, and the `[hulkrna-merge]` log format. That is the entire
coupling. **The package is genuinely liftable.**

`lib/hulkrna/__init__.py` (`__version__ = "0.2.0"`) exists only to satisfy that one
import and becomes dead once the tool leaves. Harmless to keep.

### 2.3 Four things that look like details and are not

1. **`_DISAGREE_Q = _build_disagree_table()` runs at import time** (`pairs.py`), building
   a 127×127 table with two `pow` and a `log10` per cell — **~7.9 ms**. In hulkrna
   nothing imports `pairs` until the CLI runs. **Do not re-export `zna.merge` from
   `zna/__init__.py`**, or every `import zna` anywhere pays that. Import it lazily from
   the subcommand, or make the table build lazily on first use.
2. **The no-op `njit` shim is public surface, not an implementation detail.**
   `test_read_merge.py` imports `njit` from `overlap` to compile its `legacy_scan`
   reference. Keep the symbol exported.
3. **`scripts/hulkrna-merge` exists only because hulkrna has no packaging** — there is no
   `pyproject.toml` or `setup.py` in that repo, so the shim prepends `lib/` to
   `sys.path`. zna is a real package with `[project.scripts] zna = "zna.cli:main"`, so
   the shim has no analogue and must not be carried across.
4. **`mp.get_context("fork")`.** The parallel path forks deliberately, so workers inherit
   the parent's already-JIT-compiled kernel instead of each paying the numba compile.
   fork is deprecated-by-default on macOS and is moving that way in CPython generally.
   Under `spawn` the JIT-warming trick stops working and every worker recompiles — check
   this on zna's supported Python versions and, if fork is unavailable, either accept the
   per-worker compile or set `NUMBA_CACHE_DIR` so the compiled kernel is reused. The
   kernels already pass `cache=True`.

---

## 3. Target layout in zna

zna today (`src/zna/`): `cli.py` 1916, `core.py` 1175, `_accel.cpp` 1291,
`_shuffle.py` 356, `_pycodec.py` 347, `codec.py` 157, `dtypes.py` 131, `__init__.py` 28.

Put the merger in a **subpackage**, not loose modules — it is a coherent unit with its
own kernel, and keeping it separable is what makes the later C++ swap and the
`--merge-pairs` integration tractable:

```
src/zna/merge/__init__.py      exports: MergeParams, PairOutcome, process_pair,
                                        find_overlap, reverse_complement,
                                        score_weights, threshold_bits
src/zna/merge/overlap.py       the kernel (later: a C++ backend beside it)
src/zna/merge/pairs.py         decision + record construction
src/zna/merge/fastqio.py       paired FASTQ reader/writer, InputError
src/zna/merge/cli.py           argparse + stats + run paths  → becomes `zna merge`
tests/test_merge.py            from tests/test_read_merge.py
tests/test_merge_encode.py     from tests/test_merge_zna_e2e.py (see §6.3)
```

Rationale for `zna/merge/` over `zna/_merge.py`: zna's flat modules are all *format*
concerns (codec, shuffle, dtypes). Merging is a data-processing concern with a different
dependency profile (numba) and a different lifecycle (a C++ rewrite is planned). A
subpackage keeps `pip install zna` importable when numba is absent — see §5.

---

## 4. The contract — read this before changing anything

These are the properties the rest of the chain depends on. They are not stylistic. Each
was violated at some point and each violation was silent.

**C1 — base 0 of every emitted read is a true fragment boundary.** Every Illumina read
begins at a fragment end and reads inward. The tool must never remove bases from a
read's **5'** end. All construction paths take from the start (`s1[:take1]`, `s2[:keep2]`);
trimming only ever cuts R2's **3'** end. zna's `IS_RC` flag is a claim about exactly this.

**C2 — a merged record is its fragment, exactly.** `_build_merged` constructs from the
one quantity the scan infers — the offset `s` of `revcomp(R2)` on the shared axis — so
the fragment is `[0, L)` with `L = s + len2` for *every* geometry:

```python
s     = shift if direction == FORWARD else -shift
L     = s + len2
take1 = min(len1, L)
seq   = s1[:take1] + s2rc[take1 - s : L - s]      # len(seq) == L, identically
```

Do **not** reintroduce per-direction case analysis. The previous version did, carried an
implicit `len1 ≥ len2` assumption, and truncated **374 of 137,796 merged records
(0.271%)** on a production library — each one stamped `IS_FULL_FRAGMENT` while missing
bases, i.e. an EOS token at an interior position. Truncation required `len1 < len2`
*strictly*, which is why every equal-length test fixture was blind to it.

**C3 — no orphans, ever.** `process_pair` is all-or-nothing for unmerged pairs: if either
mate falls below `min_read_length`, the whole fragment is discarded. A lone mate would be
encoded as a spurious "single" — a full molecule with both endpoints. Upstream
`samtools fastq -s /dev/null -0 /dev/null` discards singletons for the same reason.

**C4 — merged reads stay in R1's frame.** zna normalizes an unpaired record with the
read1 rule. An R2-framed merger would put an entire library in the antisense frame with
no error raised.

**C5 — the trim guard.** A trim that would leave R2 below `min_read_length` keeps both
reads *untrimmed* rather than discarding the fragment. Counted as
`trim_guard_kept_untrimmed`; expected ~0.

**C6 — pair adjacency survives the parallel path.** `_iter_chunks` splits only on whole
pairs and each worker's output blob is written atomically, so `/1`,`/2` records stay
adjacent for zna's interleaved parser.

**C7 — `--treat-unpaired-as-merged` is a claim about the data, not a formatting flag.**
It asserts every unpaired record spans its whole fragment. True here only because of C2
and C3. §8 removes the need for it.

---

## 5. Packaging: numba must be optional

zna is a format library people should be able to `pip install` casually. numba drags
llvmlite and hard Python-version pins (hulkrna's env already carries
`numba>=0.61 # py3.13 needs >=0.61`).

`pyproject.toml` already has an `[project.optional-dependencies]` table with a `dev`
group, so the mechanism exists:

```toml
[project.optional-dependencies]
merge = ["numba>=0.61"]
```

Three rules:

1. **`overlap.py` keeps its `try/except ImportError` fallback.** It is what makes
   `import zna.merge` work without numba, and it is the reference implementation the C++
   backend will be differentially tested against (§9).
2. **`zna merge` must refuse to run without numba rather than silently running ~50x
   slower.** Today `cli.run()` logs a warning; in zna make it a hard error with an
   `--allow-slow` escape. A silently-correct 50x slowdown on a cluster looks like a slow
   node, which is precisely the failure mode that cost a day in the `mp.Pool` hang.
3. **Do not add numba to zna's core `dependencies`.**

Also verify the conda recipe (`conda/meta.yaml`) — the merge extra needs a decision
there: either a separate `zna-merge` output or numba as a `run_constrained` entry.

---

## 6. Port procedure

Land it as **three commits**, in this order. Every step ends green.

### 6.1 Commit 1 — move the code, unmodified

```bash
mkdir -p src/zna/merge
cp /Users/mkiyer/proj/hulkrna/lib/hulkrna/merge/*.py src/zna/merge/
```

Then the mechanical edits:

- `from .overlap import ...` / `from .pairs import ...` / `from .fastqio import ...`
  already use relative imports — **no change needed**.
- `cli.py`: `import hulkrna` → `import zna`; `hulkrna.__version__` → `zna.__version__`;
  `"tool": "hulkrna-merge"` → `"zna-merge"`; logger `hulkrna.merge` → `zna.merge`;
  `prog="hulkrna-merge"` → `prog="zna merge"`; log format `[hulkrna-merge]`.
- `__main__.py`: keep, so `python -m zna.merge` works.
- Docstring references to `docs/READ_MERGE_REDESIGN.md` need repointing to wherever that
  document lands in zna (§10).

**Do not "clean up" anything in this commit.** Several oddities are load-bearing and
documented in place: the branchless block-8 loop, the `_DISAGREE_Q` table, the sparse
histograms, `InputError`, the `2*P` submit window.

### 6.2 Commit 2 — wire up `zna merge`

zna's CLI is argparse subparsers in `main()` at `src/zna/cli.py:1783`:

```python
subparsers = parser.add_subparsers(dest="command", required=True)
enc  = subparsers.add_parser("encode", ...)     # :1789
dec  = subparsers.add_parser("decode", ...)     # :1863
insp = subparsers.add_parser("inspect", ...)    # :1876
shuf = subparsers.add_parser("shuffle", ...)    # :1892
```

The merge CLI currently owns its own `ArgumentParser` (`build_parser()`). Refactor it to
an `add_merge_parser(subparsers)` that adds the same flags to a subparser, and dispatch
to `zna.merge.cli.run(args)`. Keep `build_parser()` as a thin wrapper so the existing
tests (which call `cli.build_parser().parse_args([...])` ~30 times) keep working — that
is the cheapest way to move 173 tests without touching them.

The full flag set to preserve:

```
--in1 --in2 --out --json --threads --processes --chunk-size --compress-level
--threshold-merge (28.0) --threshold-trim (8.0) --min-read-length/--length-required (40)
--no-sync-check --allow-empty -q/--quiet
```

### 6.3 Commit 3 — tests

`tests/test_read_merge.py` → `tests/test_merge.py`. Edits needed:

- `from hulkrna.merge...` → `from zna.merge...` (3 module-level, 2 function-local at
  `:675` and `:758`, 1 inside the `_KILL_DRIVER` subprocess template at `:1146`).
- The `_KILL_DRIVER` template injects `sys.path.insert(0, <lib dir>)` — in zna the
  package is installed, so that line goes.
- hulkrna's `tests/conftest.py` only provides a `sys.path` shim and an unrelated
  `genome_fa` fixture. **Nothing else is needed**; do not port conftest.

`tests/test_merge_zna_e2e.py` → `tests/test_merge_encode.py`. This one becomes *more*
valuable inside zna, because both sides are now local:

- It runs the real `zna` CLI by `subprocess` and reads
  `hulkrna/resources/zna/label_defs.yaml`. In zna, write a small label-defs fixture in
  the test itself rather than depending on hulkrna.
- Its `pytest.importorskip("zna")` and the `>= 0.3.5` version gate become unnecessary.
- **Keep the independently-drawn mate lengths** (`l1`, `l2` in `write_fastqs`). That two-
  line detail is the only reason the suite can see C2 violations at all; with both mates
  at `READLEN` the geometry assertions are structurally unable to fail.

### 6.4 Verify

```bash
pip install -e ".[merge,dev]"
pytest tests/ -q                      # zna's own suite must stay green
pytest tests/test_merge.py tests/test_merge_encode.py -q     # 179 tests
pip uninstall numba && pytest tests/test_merge.py -q          # pure-Python path
```

The pure-Python run is not optional: production uses numba and CI usually will not, so
that is the configuration where a divergence hides.

---

## 7. What stays in hulkrna, and what changes there

**Stays:**
- `lib/hulkrna/gather/tools/read_merge.py` — parses the merge JSON into the cohort store.
  It is a hulkrna concern and consumes only the JSON contract (§11).
- `config/config.yaml` `read_merge:` block — pipeline configuration.
- `workflow/rules/zna.smk` `rigel_read_merge` rule.
- `workflow/envs/read_merge.yaml` — but it becomes `zna>=X` + `pigz` instead of
  `numba`/`numpy`, and it may merge with `workflow/envs/zna.yaml` entirely once the tool
  ships inside zna.

**Changes:**
- `scripts/hulkrna-merge` is deleted. The rule invokes `zna merge` from the conda env
  instead of `python {params.merge_script}`.
- `tests/test_read_merge.py` and `tests/test_merge_zna_e2e.py` are deleted from hulkrna
  (they live in zna). Consider keeping one thin hulkrna test that the *rule* invokes the
  tool correctly — the seam moves, it does not vanish.
- `lib/hulkrna/merge/` is deleted.

**Sequencing note.** Do not delete from hulkrna until zna has shipped a release
containing the tool and `workflow/envs/read_merge.yaml` can pin it. Until then hulkrna
runs its own copy. A brief period of duplication is much cheaper than a pipeline that
cannot build its environment.

---

## 8. Phase 2 — `zna encode --merge-pairs`

This is the reason the port is worth doing, and zna already has the ingest path.

`encode_command` consumes a record stream. For ZNA→ZNA re-encode it already uses the
**geometry-carrying 6-tuple**:

```python
stream = stream_inputs(args, with_ends=True)     # cli.py:1135
for unit in _fragment_units(stream):
    if preserve_normalization:
        for rec in unit:
            is_rc, is_full = _flags_from_ends(rec[4], rec[5])       # cli.py:1145
            writer.write_record(rec[0], rec[1], rec[2], rec[3],
                                is_rc=is_rc, is_full_fragment=is_full)
```

A merging input strategy yields exactly that shape:
`(seq, is_paired, is_read1, is_read2, has_start, has_end)` — merged record →
`(True, True)`, mate → `(True, False)` or `(False, True)` per its orientation. **No new
writer plumbing.** The merger becomes a second producer for a path that already exists
and is already tested by the re-encode suite.

What this deletes outright:
- the `reads_merged.fastq.gz` intermediate (a pigz write + a gzip read + a full FASTQ
  re-parse per library; measured **17.1 MB** of gzip for a 168k-pair library, and
  production libraries are ~200x that),
- re-pairing records by inspecting `/1`,`/2` suffixes,
- the "strip the suffix so zna classifies merged reads as singles" convention,
- **`--treat-unpaired-as-merged` entirely** — the geometry is passed in-process, so
  C7 stops being a claim anyone can get wrong.

Note `--shuffle` is a **post-pass** (`shuffle_zna` on a temp file, then `os.replace`,
`cli.py:1331`), so it does not constrain where records come from.

Keep `zna merge` (FASTQ→FASTQ) after this lands: other pipelines want the intermediate,
and it is the easiest way to A/B the merger against fastp/bbmerge.

---

## 9. Phase 3 — numba → C++

zna is already a C++-accelerated codebase (nanobind + scikit-build-core + CMake;
`_accel.cpp` is 1291 lines), and it already has **exactly the testing discipline this
needs**: `_pycodec.py` and `_accel.cpp` are two implementations of one codec, pinned
equal by `test_cross_backend_equivalence`, `test_cross_backend_is_rc_identical`,
`test_pass_through_encode_is_backend_identical`, plus `tests/test_fuzz_roundtrip.py`
(1272 lines, which found a real data-corruption bug on its first run).

Adopt that shape:

1. **The pure-Python scorer becomes the reference backend, not a fallback.** Do not
   delete it. It is the oracle.
2. **Move to fixed-point integer scoring.** The redesign floats this (§8b.1, "hundredths
   of a bit"); float64 was kept for clarity. In C++ it is both faster and exactly
   reproducible — no FMA contraction, no `-ffast-math` surprise, no platform-dependent
   tie-breaking in the `argmax`. This is training data; a given FASTQ should produce a
   byte-identical ZNA anywhere. The one place it currently bites is
   `dmax = int((ceiling - best) / step)`, where a float truncation decides control flow.
3. **Any differential test must be order-independent** — assert the returned shift is
   *an* argmax with matching `n`/`d`/`score`. Ties occur on **0.845%** of real pairs, and
   an order-replicating reference would break on exactly the SIMD rewrite it exists to
   guard.
4. **Write the block-loop test first** (it is already in `test_merge.py`:
   `test_block_loop_sees_every_position`). It sweeps one mismatch across all 40 positions
   of a full overlap. A `k += 7` stride bug in the unrolled loop changes 6.34% of scores
   and 0.26% of merge/trim/keep decisions on real data and **passes every other test in
   the suite**. That is the exact class of bug a hand-written SIMD loop introduces.

**What C++ unlocks that numba cannot.** The branchless block-8 loop is a workaround for
what numba's LLVM path will and will not vectorize; it was worth 3.7x. With 2-bit packing
and popcount over 64-bit words it is ~8x fewer operations per position, at the cost of
losing the early bail — which is why it must be measured, not assumed. In numba it
measured **1.19x on the whole single-process tool and ~0 in production**, because a null
njit kernel with the same signature costs **0.221 µs/pair of pure dispatch** — a ceiling
on what *any* numba-side kernel change can return. In C++ that ceiling does not exist.

Second prize: a GIL-releasing kernel makes the parallel path threads instead of
processes, and most of `cli.py`'s chunking, pickling and blob-merging deletes.

---

## 10. Documentation to carry across

Move into `zna/docs/`, adjusting cross-references:

- **`READ_MERGE_REDESIGN.md`** — the design. Its header records three implementation
  decisions the design left open (exact match weight, why the tie margin is not reported,
  the branchless block loop) and its §7 now records the composition-aware `p_null` guard
  as measured-negative three times. Essential.
- **`MERGE_TOOL_AUDIT.md`** — 51 candidates, each adversarially verified, with the
  measurements. This is where "why didn't you do X" is answered; §(c) is a list of things
  deliberately not done, each with the number that closed it.
- **`READ_MERGE_ROADMAP.md`** — fold its "After the freeze" section into this document
  and retire the rest, or keep it as the zna-side status board.

`READ_MERGE_TOOL_DESIGN.md` (the original fastp-derived design) was **retired** on
2026-08-12; anything still true about I/O and the output contract is captured in §1 and
§4 here. Git history has it if the archaeology is ever needed.

---

## 11. The JSON stats contract

`_assemble_stats` emits these. `hulkrna/lib/hulkrna/gather/tools/read_merge.py` and
`hulkrna/docs/cohort_schema.md` consume them, so **renaming a key breaks the cohort
store silently** (the gather tool takes whatever scalars it finds).

Provenance: `tool`, `tool_version`, `numba`, `python`, `elapsed_s`, `pairs_per_second`.
Counts: `input_pairs`, `merged`, `trimmed_pairs`, `kept_pairs`, `merged_pct`,
`trimmed_pct`, `kept_pct`, `emitted_records`, `dropped_below_min_length`,
`fragments_dropped_short_mate`, `bases_trimmed`, `bases_consensus_changed`,
`trim_guard_kept_untrimmed`, `mean_emitted_length`.
Calibration: `overlap_mismatch_rate`, `overlap_bases_compared`.
Params block: `threshold_merge_bits`, `threshold_trim_bits`, `err_rate`, `match_bits`,
`mismatch_bits`, `min_read_length`.
Histograms: `length_histogram`, `overlap_length_histogram`, `insert_size_histogram`,
`insert_size_censoring`.

`overlap_mismatch_rate` is the one to watch: it reads **0.009113** against the assumed
`err_rate` of 0.01 on real data, and it is the *sensitive* degradation channel — 5%
per-base degradation moves it 5x while `merged_pct` moves 1%. It cannot approach the 0.75
of chance alignment; the threshold caps the observable rate near 0.22.

---

## 12. Traps

Things that cost real time to discover. Each is a comment in the code; this is the index.

1. **Equal-length fixtures cannot see the truncation class** (C2). Any new geometry test
   must draw `len1` and `len2` independently.
2. **The block-8 loop needs its own test.** See §9.4.
3. **Anything added to the `_process_chunk` return tuple is paid per chunk, and
   `--chunk-size` multiplies it.** Adding two histograms while cutting chunk size 25x put
   three dense 1025-element lists through the pickle and took p8 from 4.6 to 10.1
   µs/pair. Histograms now ship sparse.
4. **`mp.Pool` cannot detect worker death** — it respawns the process while the parent
   blocks forever on a result that will never come (measured: alive at 150 s at 0% CPU).
   `ProcessPoolExecutor` raises `BrokenProcessPool` in ~1 s. Do not "simplify" back.
5. **The reader must validate `len(seq) == len(qual)`.** A FASTQ truncated inside its
   *final* quality line otherwise satisfies `if not q` and is emitted malformed with
   rc=0.
6. **`T = log2(N/α)` bounds detection against chance alignment of uniform sequence**
   (0/40,000 at every read length from 50 to 300), **not** against real repeat-rich
   sequence (1.7% at T=28, still 1.3% at T=100). Do not quote α as a real-data guarantee.
7. **N scores as an ordinary base** (`N` vs `N` earns a full match). Harmless only
   because `--npolicy drop` discards N-containing fragments; 0.169% of pairs, every run
   length exactly 1, 0 merge decisions changed.
8. **`reverse_complement` passes IUPAC codes through uncomplemented** —
   `rc(b"RYKMSWBDHVN") == b"NVHDBWSMKYR"`. Deliberate: remapping to N would change the
   kernel's N-vs-N semantics. Exposure is nil on real data.
9. **Do not gate merges on shift ambiguity.** It looks compelling — 12–16% of merges have
   a competing shift above threshold — but it is a microsatellite detector, not an
   identifiability test: posterior mass off the argmax is 1.60%, the argmax matches
   aligner ground truth 90.4% of the time, and the gate would re-emit **+4.68% duplicated
   bases**. `MERGE_TOOL_AUDIT.md` §(c) has the full refutation.

---

## 13. Checklist

- [x] `src/zna/merge/` created; six modules moved; `hulkrna` references replaced
- [x] `numba` declared as an optional extra; `import zna.merge` works without it
- [x] `zna merge` subcommand registered and dispatching
- [x] `zna merge` refuses to run without numba (unless `--allow-slow`)
- [x] `tests/test_merge.py` — 173 passing, with and without numba (+4, see §14)
- [x] `tests/test_merge_encode.py` — 6 passing, independent mate lengths preserved
- [x] zna's existing suite still green
- [x] `conda/meta.yaml` decision recorded for the numba extra
- [x] CHANGELOG entry; version bumped
- [x] Design docs moved and cross-references fixed
- [ ] hulkrna: `workflow/envs/read_merge.yaml` repinned, rule switched to `zna merge`,
      old copy deleted — **only after a zna release ships the tool**

---

## 14. What actually landed, 2026-08-12

Three commits, as planned. 177 + 6 tests green with numba (Python 3.12) and without it
(Python 3.14); zna's own 324 unchanged. Both instruments were re-verified against
deliberate mutants after the move: regating the R2 tail on
`direction == FORWARD and shift > 0` (the pre-fix construction, i.e. a C2 violation)
fails 10 tests including **both** parametrizations of the e2e boundary test, and
`k += 7` in the unrolled block loop is killed by `test_block_loop_sees_every_position`
and by nothing else.

**Two deviations, both forced by measurement.**

1. **Trap #1 was much bigger than priced, and the fix is different.** This document
   costs eager import at the ~4 ms `_DISAGREE_Q` build. The real figure is **~170 ms**,
   because `overlap.py` imports numba: `import zna.cli` measured 40 ms, and
   `import zna.cli, zna.merge.cli` measured 210 ms. Registering the subcommand eagerly
   in `main()` would therefore have paid that on *every* `zna` invocation — absurd for
   `zna inspect --json`, which this repo advertises as fast enough to catalogue a
   corpus. So the argparse half was split into **`src/zna/merge/args.py`**, which
   imports nothing but `argparse`, and **`zna/merge/__init__.py` resolves its exports
   lazily (PEP 562)** so that reaching `args` does not execute the package's eager
   imports. `zna/cli.py` imports `args` to register the subparser and reaches for
   `merge.cli.run_command` only after dispatch. Measured after: `import zna.cli` 40 ms,
   `import zna.merge` 5.7 ms with numba untouched. `zna/__init__.py` still does not
   mention the package, as §2.3 requires.

2. **The numba refusal sits in the CLI entry points, not in `run()`.** §5.2 says
   "`zna merge` must refuse to run without numba"; putting the check in `run()` would
   have meant the numba-free test run — the one §6.4 calls non-optional — could not
   execute a single CLI test without `--allow-slow` on 10 call sites. `_require_numba`
   is called from `run_command` (which serves both `zna merge` and
   `python -m zna.merge`); `run()` stays usable in-process, which is also what
   `--merge-pairs` (§8) will call. It still logs the warning.

**Four tests added** beyond the 173, all for behaviour that is new in zna and would
otherwise be unpinned: the CLI refuses without numba, `--allow-slow` gets past it,
`run()` never refuses, and `zna merge` is really reachable through zna's top-level
parser (checked out of process — an in-process check would miss the dispatch).

**Confirmed, not assumed.** `mp.get_context("fork")` (§2.3.4) works on Python 3.12 and
3.14 on macOS, with no `DeprecationWarning` raised under `-W error::DeprecationWarning`;
the parallel path produces a record set byte-identical to the single-process path, and
to the pure-Python path, on the same input.

**Conda (§5).** `conda/meta.yaml` gets numba as a `run_constrained` entry rather than a
separate `zna-merge` output: it keeps one package and one recipe, states the version
floor where a solver can act on it, and does not drag llvmlite into every `zna` install.
`zna merge --help` is added to the recipe's test commands, which exercises the
subcommand registration without needing numba present at test time.
