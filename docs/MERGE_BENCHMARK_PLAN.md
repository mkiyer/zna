# Benchmarking `zna merge` against fastp on simulated ground truth

Status: **specified, not built**. Written 2026-08-12, against zna at `f5db18f`
(0.4.0-unreleased). This is the next piece of work and it **gates the 0.4.0 tag**.

This document is written for a session that has not seen the merge work. It should be
enough to build and run the benchmark without re-deriving anything.

---

## 1. Why this exists

Everything that has been verified about `zna merge` so far proves that it **matches
itself**: byte-identical to the numba implementation it replaced, byte-identical across
the reference and compiled backends, across thread counts, across platforms. That is a
strong *correctness* argument and a **zero** *accuracy* argument.

Nothing has ever asked: when the tool says two reads overlap by 113 bases, is that
true? The redesign's original numbers ([READ_MERGE_REDESIGN.md](READ_MERGE_REDESIGN.md)
§5) came from a simulation that no longer exists in runnable form, and the production
figures in [MERGE_TOOL_AUDIT.md](MERGE_TOOL_AUDIT.md) are self-consistency checks
against an aligner, not controlled ground truth.

So: simulate reads whose true fragment length, true overlap and injected errors are
known exactly, run both mergers, and score each against the truth rather than against
each other.

**The property that matters most** is not the merge rate. It is contract C1/C2 from
[READ_MERGE_PORT_TO_ZNA.md](READ_MERGE_PORT_TO_ZNA.md) §4 — *base 0 of every emitted
read is a true fragment boundary, and a merged record is its fragment exactly*. The
whole ZNA fragment-boundary contract rests on it, and it has only ever been checked
against the tool's own inferences. Independent truth can check it properly.

---

## 2. What already exists

- **`~/proj/khorana/scripts/sim.py`** + `khorana/src/khorana/sim/reads.py` — a
  transcript-level RNA-seq simulator. Useful for its **error-injection model**
  (`_apply_snvs` / `_apply_insertions` / `_apply_deletions`, shuffled per read so errors
  compound), its fragment-length sampling, and its multiprocessing FASTQ writer.
  **It does not simulate adapter read-through at all**, and it needs a GTF. Read it for
  structure; do not try to extend it.
- **`zna/scripts/merge_bench/`** — the benchmark harness from the C++ work:
  `gen_library.py` (a crude 2×150 generator, no ground truth), `bench_breakdown.py`,
  `bench_scan.cpp`, `bench_simd.cpp`, `proto_merge.cpp`, `asan_scan.cpp`, and a README
  that reproduces every number in [MERGE_CPP_DESIGN.md](MERGE_CPP_DESIGN.md). The new
  scripts belong here beside them.
- **fastp** is installed at `/Users/mkiyer/sw/miniforge3/envs/fastp`.
- **The working environment** is `zna_merge`
  (`/Users/mkiyer/sw/miniforge3/envs/zna_merge`), which has zna installed editable plus
  `pigz`. Put its `bin/` on `PATH` or `shutil.which("pigz")` returns None and the tools
  silently fall back to stdlib gzip.

---

## 3. `simulate.py` — the simulator

```
python simulate.py --genome hg38.fa --n-pairs 1000000 --read-length 150 \
                   --frag-min 60 --frag-max 450 --error-rate 0.002 \
                   --out-prefix sim --seed 42
```

**Inputs**

| flag | meaning |
|---|---|
| `--genome` | indexed FASTA (`.fai` beside it; use `pysam.FastaFile` or read the index directly) |
| `--n-pairs` | fragments to draw (1M is the target scale) |
| `--read-length` | cycles per mate, default 150 |
| `--frag-min` / `--frag-max` | **uniform** fragment length, inclusive |
| `--error-rate` | per-base substitution rate, default 0.002 |
| `--adapter1` / `--adapter2` | default Illumina TruSeq (see §3.3) |
| `--seed` | everything derives from it; two runs with one seed are identical |

**Why uniform rather than a realistic insert distribution.** The point is to populate
*every regime* at equal density, not to look like a library:

| fragment length | regime |
|---|---|
| `> 2 × readlen` | no overlap — any merge is a **chimera** |
| `= 2 × readlen` | overlap 0, the boundary case |
| `readlen < L < 2 × readlen` | partial overlap, 1..readlen bases — **the sensitivity curve** |
| `= readlen` | full overlap, mates are exact reverse complements |
| `< readlen` | **read-through**: both mates run into adapter |

`--frag-min 60 --frag-max 450` at 2×150 covers all five. Sample `L` uniformly, then
draw the genomic start uniformly over the concatenated non-N sequence.

### 3.1 Building a pair

```
frag   = genome[chrom][start : start+L]          (revcomp it 50% of the time)
r1_raw = (frag + adapter1)[:readlen]             pad with random bases if still short
r2_raw = (revcomp(frag) + adapter2)[:readlen]
```
Then inject substitution errors **independently per mate**, recording position and
base change. Quality strings: constant `I` is fine for a first pass, but see §7 — the
posterior consensus is quality-driven, so a second mode with realistic binned qualities
(NovaSeq: 2/12/23/37) is worth having, because a constant Q means every mismatch is a
tie and R1 always wins.

**True overlap length** = `max(0, 2 × readlen − L)` clipped to `readlen`, i.e. the
number of fragment bases both mates cover. For `L < readlen` the mates cover the whole
fragment and read through; call that `rt=1` and set `ovl = L`.

### 3.2 Ground truth: sidecar TSV is authoritative

Write **both**:

- **`<prefix>.truth.tsv`** — the authoritative record, one line per pair:
  `read_id, chrom, start, strand, frag_len, true_ovl, read_through, len1, len2, n_err1, n_err2, err_positions`
- a **compact copy in the FASTQ comment**, for eyeballing and for tools that preserve
  it: `@sim0000123/1 fl=187 ov=113 rt=0 e1=2 e2=1`

Use the sidecar for scoring. The comment is *not* reliable: `zna merge` rewrites merged
headers as `<id> merged_<n1>_<n2>` (preserving what follows the first whitespace), and
fastp has its own opinions. The read ID before the first whitespace is the join key, and
both tools preserve it — `zna merge` strips only the `/1`,`/2` suffix, which is exactly
what makes the ID a stable key.

### 3.3 Adapters

Illumina TruSeq, the same constants the existing tests use:

```
adapter1 = AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
adapter2 = AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
```

### 3.4 Output

`<prefix>_1.fq.gz`, `<prefix>_2.fq.gz`, `<prefix>.truth.tsv`, and a `<prefix>.json`
recording every parameter and the seed, so a result can be traced to the input that
produced it.

---

## 4. `compare.py` — the head-to-head

```
python compare.py --sim-prefix sim --out results/ [--threads 4]
```

Runs both tools, joins their output to the truth by read ID, and writes a report.

### 4.1 Invocations, and making them fair

```
zna merge --in1 sim_1.fq.gz --in2 sim_2.fq.gz --out zna.fq.gz \
          --json zna.json --min-read-length 40 --threads 4

fastp --in1 sim_1.fq.gz --in2 sim_2.fq.gz \
      --merge --merged_out fastp.merged.fq.gz \
      --out1 fastp.un1.fq.gz --out2 fastp.un2.fq.gz \
      --length_required 40 --json fastp.json --html /dev/null
```

**Fairness needs care, and the conclusions depend on it. Record exactly what was
passed.**

- fastp does quality filtering, low-complexity filtering, polyG trimming and adapter
  trimming *by default*; `zna merge` does none of those. Disable them
  (`--disable_quality_filtering --disable_length_filtering` and consider `-A`) so the
  comparison is of *merging*, not of preprocessing. But note that fastp's adapter
  trimming interacts with read-through: with `-A`, does fastp still strip adapter from
  merged reads? **Check this empirically before trusting any read-through numbers.**
- Match `--length_required` to `--min-read-length`.
- fastp emits merged reads to a separate file; `zna merge` emits one mixed interleaved
  stream. Normalise both to `{read_id: (merged?, sequence)}` before scoring.
- fastp reorders nothing, but do not rely on order — join by ID.

### 4.2 Metrics, all scored against the sidecar

Per tool, and **binned by true overlap length** (the curve is the deliverable):

1. **Merge sensitivity** — of pairs with true overlap ≥ k, what fraction merged.
2. **Merge specificity / chimera rate** — of pairs with true overlap 0, what fraction
   merged. Any non-zero value here is serious; the redesign's whole claim is that this
   is ~0 at 28 bits.
3. **Reconstruction accuracy** — for merged records: does the emitted sequence equal the
   true fragment? Report exact-match rate, and for mismatches, the length error
   (`len(emitted) − frag_len`) and the edit distance. A wrong *length* inside a repeat
   is a different failure from wrong *bases*.
4. **Base accuracy in the overlap** — the consensus is supposed to *fix* errors. Count
   errors in the emitted record against the true fragment, split by whether the position
   was in the overlap. Compare to the naive "R1 always wins" baseline. This is where the
   35.2% claim in [READ_MERGE_REDESIGN.md](READ_MERGE_REDESIGN.md) gets tested for real.
5. **Boundary violations (C1/C2)** — for every emitted record, is base 0 a true fragment
   end? For merged records, are *both* ends true? **This must be exactly zero.**
6. **Unmerged handling** — of pairs that should not merge, are both mates emitted
   intact? No orphans (C3).
7. **Throughput** — wall time and µs/pair for both, on the same input and thread count.

### 4.3 Output

A `results/report.md` with the tables and curves, `results/summary.json` for machine
consumption, and — most valuable — **`results/zna_errors.tsv`, every pair where
`zna merge` got it wrong**, with the truth, what it emitted, and the score it assigned.
That file is the point of the exercise.

---

## 5. Acceptance criteria for the 0.4.0 tag

Not thresholds to game — a checklist of things to look at before shipping:

- [ ] **Zero boundary violations** (metric 5). This is a contract, not a statistic.
- [ ] **Zero chimeras** at true overlap 0 in 1M pairs, or an explanation of every one.
- [ ] Merge sensitivity at overlap ≥ 15 (the theoretical floor, `ceil(28/1.9855)`)
      matches the predicted step: essentially 0 below 15, essentially 1 above.
- [ ] Read-through fragments (`L < readlen`) reconstruct to the fragment with adapter
      removed — this is the geometry that used to truncate 0.271% of merges.
- [ ] Reconstruction accuracy ≥ fastp's, or the difference is understood.
- [ ] `zna_errors.tsv` has been *read*, not just counted.

---

## 6. Two known open items to fold in

**Histogram bins are clamped at 1024** (`fastq_chunk.hpp: HIST_MAX`). Anything longer
piles into bin 1024, so with long reads the length and insert distributions silently
aggregate. The fix is not a sparse dict — that puts a hash in the per-record path — but
to **grow the three dense arrays with the `Scratch` arena** (a merged record is at most
`len1 + len2`, so `2 × cap + 1` bins suffice), giving O(1) indexing and no cap. Small,
contained, wants a test. Do it before the benchmark so the benchmark's own histograms
are trustworthy at any read length.

**`insert_size_censoring.floor` duplicates `params.min_read_length`.** The block exists
to say the insert histogram is censored at both ends, but only `min_mergeable_overlap`
is information the JSON does not already carry. Drop `floor`, or keep it and say in the
docs that it is a convenience copy. Either is defensible; do not leave it unexplained.

---

## 7. Gotchas worth knowing before you start

1. **Constant quality strings make the consensus untestable.** With every base at `I`,
   a mismatch is always a tie and R1 always wins, so metric 4 measures nothing. Give
   `simulate.py` a `--quality-model {flat,novaseq}` and run the accuracy comparison on
   `novaseq`.
2. **Uniform fragment lengths are not a library.** Do not quote the merge rate from this
   benchmark as if it were a production number — it is a coverage tool. Quote per-bin
   sensitivity instead.
3. **The genome introduces repeats, and that is the point** — but it also means some
   "wrong" merges are the tool finding a *real* alternative alignment. Expect a
   microsatellite tail; [MERGE_TOOL_AUDIT.md](MERGE_TOOL_AUDIT.md) §(c) already measured
   that gating on shift ambiguity costs more than it saves, so treat repeat-driven
   disagreements as expected and separate them out in the report.
4. **Draw fragment starts away from N runs**, or the N-policy interacts with the
   measurement.
5. **1M pairs at 2×150 is ~600 MB of FASTQ.** Generation is the slow part; write it
   once, keep it, and make `simulate.py` deterministic from the seed so it can be
   regenerated rather than stored.
6. **Both tools must see identical input.** Generate once; do not re-simulate per tool.
