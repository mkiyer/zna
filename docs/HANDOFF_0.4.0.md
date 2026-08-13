# Handoff: finishing 0.4.0

Written 2026-08-13. Branch `zna-0.4.0-hardening`, two commits, working tree clean,
609 tests green compiled / 562 + 47 skipped extension-less.

This is what is left, why it is left, and the traps between here and the tag.

---

## 1. What remains

### 1.1 Per-record provenance — the only real feature left

Agreed but not built. Run-level reporting landed (`--merge-json`, the merge summary, a
warning above 1% of emitted bases), so a *library* now says what happened to it. A
*read* does not.

Two halves, for two different readers:

- **Header tokens on the merged FASTQ**, for human inspection: `trim3_<n>` for bases cut
  and `rescued_<n>` for no-calls recovered, beside the existing `merged_<n1>_<n2>`.
  Verified safe against label extraction — a token with no colon is skipped by the tag
  parser, exactly as `merged_150_87` is:

  ```
  f12  ZI:i:7 merged_90_0 trim3_12   ->   labels (7,)   unaffected
  ```

  The construction lives in `merge_core.hpp::process_pair` (the `snprintf` into
  `sc.name`), mirrored in `_pymerge.process_pair`. **Both must change together** or the
  cross-backend differential fails — which is the system working.

- **An optional provenance label for the corpus.** ZNA does not store headers, so the
  tokens above vanish at encode time, and under `--merge-pairs` (0.5.0) there is no
  intermediate FASTQ at all. One `C` (uint8) column, one byte per record, with bits for
  merged / trimmed / N-rescued / N-trimmed, is the only per-record provenance that
  survives into the corpus.

  This one needs a decision before it is built: it adds a column to every labeled file.

### 1.2 The encode-side `--npolicy` summary

`zna merge` reports what its policy did; `zna encode --npolicy trim3` does not. The merge
side is `src/zna/merge/cli.py` (search `no-calls:`); encode's summary is the
`[ZNA] Done. Wrote N records` line in `encode_command`. It needs a counter threaded
through the write loop, which is straightforward — `_trim3` is already a local function
there.

### 1.3 Before tagging

- Re-run the 1M benchmark end to end and record the numbers in
  `MERGE_BENCHMARK_RESULTS.md` if any moved.
- `docs/RELEASING.md` for the PyPI/Bioconda steps.
- Decide whether the conda recipe and CI need anything for the changed policy set.

---

## 2. Decisions already taken — do not reopen without new evidence

Each of these was argued and closed with a measurement. `docs/NPOLICY_PLAN.md` and
`docs/METHODS.md` carry the reasoning.

| Decision | Why |
|---|---|
| **`--merge-pairs` deferred to 0.5.0** | The one-step and two-step paths produce the *identical* corpus on the 1M library — same 1,416,630 records, same flag histogram. No correctness reason to rush. `docs/MERGE_PAIRS_PLAN.md` is execution-ready. |
| **No `--max-n`** | It is a *substitution* knob. Under trim3 the output has zero N by construction and the loss is set by the *position* of the first N, not the count. No intuition for how a user would set it. |
| **A pair with a surviving N is not merged** (option C) | Worked example: merging then trimming keeps 40 of 220 fragment bases; not merging keeps 190. C wins on data retention *and* flag honesty. A merged record's 3' end **is R2's 5' end**, so trimming it discards the anchor the 3'-only rule exists to protect. |
| **The post-trim retry reuses the original shift** | trim3 cuts 3' ends, which is where a normal overlap lives. At one fixture the residual overlap collapses 80 → 10 bases; a re-scan would refuse a merge there was 80 bases of evidence for. The shift is a claim about fragment *length*, and trimming interior bases does not change a length. |
| **No consensus on a kept pair** | Of 3,068 kept pairs with a detected overlap, **none** found the true shift and 97.3% had no true overlap. Wrong emitted bases: 208 correcting neither, 1,509 correcting R1 only, 17,870 correcting both. |
| **Stumps are allowed** | ZNA stores 1-base and 0-base records (verified). `--min-read-length` handles them on the merge path. A core read-processing tool reports what it did and leaves the length decision to the consumer. |
| **Breaking backwards compatibility is fine** | Explicitly authorised. Every unstranded-normalized file changes under the new RC coin. |

---

## 3. Traps — these cost real time in the last session

1. **Rebuild command.** `pip install -e . --no-build-isolation`, with that env's pip.
   `--no-build-isolation` is load-bearing: without it scikit-build-core rewrites
   `CMakeCache` to point at a directory pip then deletes, and alternating between the two
   forms costs a full rebuild each way.

2. **Grep the rebuild output for `Successfully built`.** A CMake failure ends with an
   unremarkable pip banner, the old `.so` stays installed, and the whole suite runs green
   against stale code. This happened.

3. **`touch` a `.hpp` after restoring it with `cp`.** CMake's dependency scan can miss a
   `cp`-restored header, so a mutation test runs against the *mutated* `.so` and the fix
   looks broken. This also happened.

4. **The ambient `PATH` starts with `envs/zna/bin`** — the py3.14 extension-less install,
   which refuses to merge by design. Use absolute paths or prepend `envs/zna_merge/bin`.

5. **`envs/zna` is maintained by hand.** `pip install -e .` there builds *both*
   extensions, turning it into a second copy of the compiled environment. Re-disable
   after any rebuild:

   ```bash
   M=/Users/mkiyer/sw/miniforge3/envs/zna/lib/python3.14/site-packages/zna/merge/_accel.cpython-314-darwin.so
   mv "$M" "$M.disabled"
   ```

   Confirm with `available_merge_backends() == ['python']`. Expect ~47 skips.

6. **`_fold` in `merge/cli.py` sums a fixed prefix** of the counter tuple and
   special-cases the maximum. A counter added past that point silently reports **zero** —
   the first version of the N-policy counters did exactly that on a library that was
   visibly being trimmed. Extend the loop.

7. **`backend._REQUIRED_FUNCTIONS` is validated by `_load` but not by
   `available_merge_backends`**, so adding a name there against a stale `.so` leaves
   `accel` listed as available while `get_merge_backend()` silently falls back to the
   50x-slower reference kernel.

---

## 4. Where the ground truth lives

- `~/proj/zna_merge_bench/run1M/` — `sim_1.fq.gz`, `sim_2.fq.gz`, `sim.truth.tsv`,
  1M pairs from hg38. The simulator emits **no N**, so anything about the N policy needs
  N injected (there is a recipe in the last session's scratch; ~1.5% of reads with one
  no-call, biased 3', is realistic).
- `~/proj/zna_merge_bench/refs/hg38.fa`.
- `scripts/merge_bench/` — `simulate.py` and `compare.py` score merged FASTQ against the
  truth sidecar. **Neither reads `.zna`**; a `.zna`-level comparator is new tooling.

---

## 5. The two properties worth re-checking after any change here

```bash
Z=/Users/mkiyer/sw/miniforge3/envs/zna_merge/bin
PATH=$Z:$PATH $Z/zna merge  --in1 sim_1.fq.gz --in2 sim_2.fq.gz --out m.fq.gz --threads 4
PATH=$Z:$PATH $Z/zna encode --interleaved --treat-unpaired-as-merged \
                            --strand-normalize --shuffle --seed 7 -o out.zna m.fq.gz
```

1. **Flag byte 24 survives `--shuffle`** — ~291,810 records. It is the bit that three
   separate copy paths used to destroy.
2. **`--restore-strand` reproduces the original base multiset exactly**, all 1,416,630
   records. That is the property the lossless copy path exists for.
