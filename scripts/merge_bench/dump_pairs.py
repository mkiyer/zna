"""Dump the first N pairs, with the shipped scan's answer, for the C++ scan benches.

    python dump_pairs.py r1.fq.gz r2.fq.gz 50000 pairs.bin

Record layout, little-endian, after a leading uint32 count:
``len1:H len2:H shift:i overlap_len:H diff:H`` then ``len1`` bytes of R1 and ``len2``
bytes of revcomp(R2).  ``shift`` is on the signed single axis of ``overlap.py`` --
negative is read-through -- which is what the C++ kernels take directly.

The expected answer travels with the pair so ``bench_scan``/``bench_simd`` can check
every variant against the shipped kernel pair by pair rather than against each other.

**What is dumped is the SCAN, not the 0.6 decision.** The benches time the scan loop
alone -- every shift, a fixed 8-bit floor, no contract, weights at ``e = 0.01`` -- and
so this calls the shipped backend's ``scan`` with exactly those arguments. The 0.6
decision wraps the same loop in a per-pair floor (~28 bits at 2x150), the
``--adapter-trimmed`` contract and the plausibility gate; none of that changes which
kernel variant is faster or whether it finds the same argmax. The two weights are the
integers the benches derive with ``llround(log2(...) * 2**24)``: 33,311,170 and
137,813,407, pinned by ``tests/test_merge.py``.
"""
import struct
import sys

from zna.merge import backend
from zna.merge.fastqio import read_pairs
from zna.merge.overlap import reverse_complement
from zna.merge.params import SCALE, MergeParams

N = int(sys.argv[3])
FLOOR_Q = 8 * SCALE                 # the benches' FLOOR_Q
params = MergeParams(error_rate="0.01")
scan = backend.active().scan
n = 0
with open(sys.argv[4], "wb") as out:
    out.write(struct.pack("<I", N))
    for r1, r2 in read_pairs(sys.argv[1], sys.argv[2], 1):
        s1, s2rc = r1[1], reverse_complement(r2[1])
        s, _score, olen, diff = scan(s1, s2rc, len(s1), len(s2rc), params.match_q,
                                     params.step_q, FLOOR_Q, 0)
        out.write(struct.pack("<HHiHH", len(s1), len(s2rc), s, olen, diff))
        out.write(s1)
        out.write(s2rc)
        n += 1
        if n >= N:
            break
print("dumped", n, "with the", backend.active_name(), "backend")
