# Handoff: 0.4.0 is done; what the next session needs

Updated 2026-08-13. Branch `zna-0.4.0-hardening`, working tree clean, **627 tests green
compiled / 574 + 53 skipped extension-less**.

Sections 3–5 are the durable part — the traps, the ground truth, and the two properties
worth re-checking after *any* change to the merge or encode path. They cost real time to
learn and they outlive this release.

---

## 1. What 0.4.0 shipped

Everything the previous handoff listed as remaining is built and measured.

- **Per-record provenance** (`docs/NPOLICY_PLAN.md` D5). Header tokens
  `trim3_<n>` / `subn_<n>` / `rescued_<n>` on **every** outcome, not only merged ones,
  built in `merge_core.hpp::build_name` and mirrored in `_pymerge._prov_name`. A record
  nothing happened to is emitted byte-unchanged and stays zero-copy — 1,332,353 of
  1,418,525 records on the 1M library at a 1.5% no-call rate.
- **The corpus column**, as a `ZN:i:<bits>` SAM tag read through the `--label` machinery
  that already shipped: `--label provenance:C:ZN`. This is what closed the open question
  "does it add a column to every labeled file?" — it does not, because declaring it is
  per file and an absent tag resolves to 0. Bits: trimmed 1, rescued 2, N-trimmed 4,
  N-substituted 8, and **deliberately no "merged"** — see D5 for why that bit was
  dropped after it turned out to put ` ZN:i:1` on 82% of records to repeat
  `IS_FULL_FRAGMENT`.
- **The encode-side `--npolicy` summary**, so both halves of the seam are accountable.
  `trim3`'s count is free; `random` substitutes inside the codec, so encode counts the
  no-calls itself rather than widening the codec ABI across two backends for a statistic.
- **`docs/RELEASING.md`** corrected against two real defects in `scripts/release.sh`
  (below), and the benchmark re-verified at the release commit.

### Two things found by auditing, not by testing

- `scripts/release.sh` hardcoded `git push origin main` while all work was on a feature
  branch — it would have pushed local `main` and then tagged a commit `main` did not
  contain. It now refuses unless you are on `main`.
- The same script's version bump is a no-op when `__init__.py` is already at the target
  (the normal case), and `git commit` with nothing staged exits non-zero, which under
  `set -e` aborted the release *after* the confirmation prompt. Now handled.

### One superstition retired

`merge_core.hpp` claimed khorana's `parse_merged_fastq` required `merged_<n1>_<n2>` to
be the last header token. khorana consumes `.zna` through `ZnaReader`; `seq.py`'s FASTQ
path is stale. The tokens are still ordered to keep `merged_` last, because that is
fastp's convention and costs nothing — but it is a convention now, not a constraint, and
the comment says so.

---

## 2. What is next

`docs/MERGE_PAIRS_PLAN.md` is execution-ready and unchanged: `zna encode --merge-pairs`
for 0.5.0, deferred because the one-step and two-step paths produce the *identical*
corpus on the 1M library. Its provenance story is already settled — it computes the same
`PairResult` and writes the same `PROV_*` bits directly, with no FASTQ in between.

Decisions closed by measurement, **not to be reopened without new evidence** — the table
in `docs/NPOLICY_PLAN.md` §2 and `docs/METHODS.md` carry the reasoning: no `--max-n`; a
pair with a surviving N is not merged; the post-trim retry reuses the original shift; no
consensus on a kept pair; stumps are allowed; breaking backwards compatibility is
authorised.

---

## 3. Traps — these cost real time

1. **Rebuild command.** `pip install -e . --no-build-isolation`, with that env's pip.
   `--no-build-isolation` is load-bearing: without it scikit-build-core rewrites
   `CMakeCache` to point at a directory pip then deletes, and alternating between the two
   forms costs a full rebuild each way.

2. **Grep the rebuild output for `Successfully built`.** A CMake failure ends with an
   unremarkable pip banner, the old `.so` stays installed, and the whole suite runs green
   against stale code. This happened.

3. **`touch` a `.hpp` after restoring it with `cp`.** CMake's dependency scan can miss a
   `cp`-restored header, so a mutation test runs against the *mutated* `.so` and the fix
   looks broken. This also happened — and it matters specifically when you are checking
   that a new fixture actually catches something, which is the one time you are relying
   on a rebuild you did not read.

4. **The ambient `PATH` starts with `envs/zna/bin`** — the py3.14 extension-less install,
   which refuses to merge by design. Use absolute paths or prepend `envs/zna_merge/bin`.

5. **`envs/zna` is maintained by hand.** `pip install -e .` there builds *both*
   extensions, turning it into a second copy of the compiled environment. Re-disable
   after any rebuild:

   ```bash
   M=/Users/mkiyer/sw/miniforge3/envs/zna/lib/python3.14/site-packages/zna/merge/_accel.cpython-314-darwin.so
   mv "$M" "$M.disabled"
   ```

   Confirm with `available_merge_backends() == ['python']`. Expect ~53 skips. Note
   `zna.is_accelerated()` stays **True** there — that is the *codec* extension, which is
   a different build target and is not what this env disables.

6. **`_fold` in `merge/cli.py` sums a fixed prefix** of the counter tuple and
   special-cases the maximum. A counter added past that point silently reports **zero** —
   the first version of the N-policy counters did exactly that on a library that was
   visibly being trimmed. The per-mate provenance counters avoid this entirely by
   summing to the existing totals rather than extending the tuple; keep it that way.

7. **`backend._REQUIRED_FUNCTIONS` is validated by `_load` but not by
   `available_merge_backends`**, so adding a name there against a stale `.so` leaves
   `accel` listed as available while `get_merge_backend()` silently falls back to the
   50x-slower reference kernel.

---

## 4. Where the ground truth lives

- `~/proj/zna_merge_bench/run1M/` — `sim_1.fq.gz`, `sim_2.fq.gz`, `sim.truth.tsv`,
  1M pairs from hg38. The simulator emits **no N**, so anything about the N policy needs
  N injected: one no-call in ~1.5% of reads, placed uniformly over the read's 3' half,
  is realistic and is what the numbers in `MERGE_BENCHMARK_RESULTS.md` were re-verified
  against. That placement is not decoration — an N inside the overlap is rescued from the
  mate, one past it survives to meet the policy, so a 3' bias is what makes both fates
  occur.
- `~/proj/zna_merge_bench/refs/hg38.fa`.
- `scripts/merge_bench/` — `simulate.py` and `compare.py` score merged FASTQ against the
  truth sidecar. **Neither reads `.zna`**; a `.zna`-level comparator is new tooling.
  `compare.py` takes ~4 minutes end to end and diffs cleanly against the recorded run
  apart from paths and the throughput line.

---

## 5. The two properties worth re-checking after any change here

```bash
Z=/Users/mkiyer/sw/miniforge3/envs/zna_merge/bin
PATH=$Z:$PATH $Z/zna merge  --in1 sim_1.fq.gz --in2 sim_2.fq.gz --out m.fq.gz --threads 4
PATH=$Z:$PATH $Z/zna encode --interleaved --treat-unpaired-as-merged \
                            --strand-normalize --shuffle --seed 7 -o out.zna m.fq.gz
```

1. **Flag byte 24 survives `--shuffle`** — 291,810 records. It is the bit that three
   separate copy paths used to destroy.
2. **`--restore-strand` reproduces the original base multiset exactly**, all 1,416,630
   records. That is the property the lossless copy path exists for. Check it as an
   order-independent checksum (XOR of per-record digests), not a set — `--shuffle`
   reorders, and holding 1.4M sequences twice proves nothing extra.

Both verified at the 0.4.0 release commit.

---

## 6. What the cross-backend fixtures are actually for

Two real bugs in the last session, and one in this one, were caught by cross-backend
fixtures that did not previously exist — and in every case the mutation left the *rest*
of the suite green.

The lesson that generalises: **a differential is only as good as the branch its fixture
reaches.** When per-mate rescue counting was added here, charging both directions to
mate 1 left the randomised N-bearing differential green under `--npolicy trim3`, because
random fixtures produce R2-rescues-R1 in bulk and R1-rescues-R2 almost never — and
R1-rescues-R2 only happens on the trim path at all, since a merged pair discards R2's
copy of the overlap. It took a hand-built pair
(`test_a_rescue_is_charged_to_the_mate_it_repaired`) to pin it.

So: when you add a counter or a branch, ask which fixture *reaches* it, and if the
answer is "a random one, sometimes", write the deterministic one too. Then verify the
fixture by mutating the C++ and watching it fail — with trap 3 in mind.
