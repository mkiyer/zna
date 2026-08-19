"""Dump the first N pairs, with the shipped kernel's answer, for the C++ scan benches.

    python dump_pairs.py r1.fq.gz r2.fq.gz 50000 pairs.bin

Record layout, little-endian, after a leading uint32 count:
``len1:H len2:H shift:i overlap_len:H diff:H`` then ``len1`` bytes of R1 and ``len2``
bytes of revcomp(R2).  ``shift`` is on the signed single axis of ``overlap.py`` --
negative is read-through -- which is what the C++ kernels take directly.

The expected answer travels with the pair so ``bench_scan``/``bench_simd`` can check
every variant against the shipped kernel pair by pair rather than against each other.
"""
import struct
import sys

from zna.merge.fastqio import read_pairs
from zna.merge.overlap import FORWARD, NO_OVERLAP, find_overlap, reverse_complement
from zna.merge.params import MergeParams

N = int(sys.argv[3])
params = MergeParams()          # the shipped thresholds; the benches derive the same
                                # weights from log2() themselves, in the same scale
n = 0
with open(sys.argv[4], "wb") as out:
    out.write(struct.pack("<I", N))
    for r1, r2 in read_pairs(sys.argv[1], sys.argv[2], 1):
        s1, s2rc = r1[1], reverse_complement(r2[1])
        d, sh, olen, diff, score = find_overlap(s1, s2rc, params)
        s = 0 if d == NO_OVERLAP else (sh if d == FORWARD else -sh)
        out.write(struct.pack("<HHiHH", len(s1), len(s2rc), s, olen, diff))
        out.write(s1)
        out.write(s2rc)
        n += 1
        if n >= N:
            break
print("dumped", n)
