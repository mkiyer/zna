import struct, sys
from zna.merge.fastqio import read_pairs
from zna.merge.overlap import FORWARD, NO_OVERLAP, find_overlap, reverse_complement

N = int(sys.argv[3]); out = open(sys.argv[4], "wb"); n = 0
out.write(struct.pack("<I", N))
for r1, r2 in read_pairs(sys.argv[1], sys.argv[2], 1):
    s1, s2rc = r1[1], reverse_complement(r2[1])
    d, sh, olen, diff, score = find_overlap(s1, s2rc, 8.0, 0.01)
    s = 0 if d == NO_OVERLAP else (sh if d == FORWARD else -sh)
    out.write(struct.pack("<HHiHH", len(s1), len(s2rc), s, olen, diff))
    out.write(s1); out.write(s2rc)
    n += 1
    if n >= N: break
out.close(); print("dumped", n)
