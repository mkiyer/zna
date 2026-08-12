# Evidence for `docs/MERGE_CPP_DESIGN.md`

Every number in the design document comes from these scripts. They are kept so the
design can be re-argued against measurements rather than recollection, and so the same
numbers can be taken on a Linux/x86 box — this session's were all aarch64, which matters
for `popcount` (design §10, §14.1).

Nothing here is part of the zna package or the test suite.

## Reproducing

```bash
cd scripts/merge_bench

# 1. a representative library: 2x150, insert ~ N(200,70) in [50,400], 0.5% error.
#    Merges at 88.6% against production's measured 88.8%.
python gen_library.py 200000 r1.fq.gz r2.fq.gz

# 2. where the per-pair time actually goes  (design §1)
python bench_breakdown.py r1.fq.gz r2.fq.gz

# 3. is the argmax tie-break a specifiable total order?  (design §5)
#    Builds ties deliberately -- random sequence essentially never ties.
python verify_tiebreak.py r1.fq.gz r2.fq.gz

# 4. scan kernel variants, full pruned scans, packing inside the timed region,
#    every variant checked against the shipped kernel pair by pair  (design §6.1)
python dump_pairs.py r1.fq.gz r2.fq.gz 50000 pairs.bin
c++ -O3 -std=c++17 -o bench_scan bench_scan.cpp && ./bench_scan pairs.bin

# 5. the whole path in C++, which must be BYTE-IDENTICAL to `zna merge`  (design §3)
c++ -O3 -std=c++17 -o proto_merge proto_merge.cpp
./proto_merge r1.fq.gz r2.fq.gz cpp.fq 40
zna merge --in1 r1.fq.gz --in2 r2.fq.gz --out py.fq --min-read-length 40 -q
cmp cpp.fq py.fq        # must be silent
```

## `proto_merge.cpp` is a prototype, not the implementation

It is the artifact that proves the design lands where it claims: on 200,000 real pairs
it emits a byte-identical file to `zna merge` at 2.79 µs/pair against 8.34. Read it for
the packed kernel, the SWAR packer and the complement-by-XOR-2, all of which carry over.

Do **not** read it as a model for the real backend. It is single-threaded, its
`popen("pigz -dc")` reader is a stand-in for the chunk protocol in design §7.3, it
builds the consensus table in C++ (design §4 says build it once in Python and pass it
across), and it copies several per-pair buffers naively (design §7.4). It also has no
error handling worth the name.

## One number that is easy to get wrong

`bench_scan.cpp` measures **with packing inside the timed loop**. An earlier version
hoisted it, which flattered the packed variant by ~20% and would have hidden the finding
that a table-driven packer costs the entire win (1.027 vs 0.631 µs/pair). The audit made
the same mistake in the other direction — its 4.8x SWAR figure was a no-bail sweep doing
~5x the real pruned scan's work. Measure the whole scan, with everything it needs.
