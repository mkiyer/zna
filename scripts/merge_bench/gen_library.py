"""Generate a realistic 2x150 library: insert ~ N(200,70) truncated [50,400],
0.5% per-base error, real adapters with read-through. ~88% merge rate."""
import gzip, random, sys

N = int(sys.argv[1]) if len(sys.argv) > 1 else 200_000
OUT1, OUT2 = sys.argv[2], sys.argv[3]
READLEN = 150
A1 = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
A2 = b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
COMP = bytes.maketrans(b"ACGTN", b"TGCAN")
rng = random.Random(20260812)
BASES = "ACGT"

def draw(n): return bytes("".join(rng.choices(BASES, k=n)), "ascii")
def rc(s): return s.translate(COMP)[::-1]
def mutate(b):
    ba = bytearray(b)
    for i in range(len(ba)):
        if rng.random() < 0.005:
            ba[i] = ord(rng.choice(BASES))
    return bytes(ba)

with gzip.open(OUT1, "wb", compresslevel=1) as f1, gzip.open(OUT2, "wb", compresslevel=1) as f2:
    for i in range(N):
        while True:
            ins = int(rng.gauss(200, 70))
            if 50 <= ins <= 400: break
        frag = draw(ins)
        # fastp trims 3' ends; ~4% of pairs end unequal, as measured in production
        l1 = READLEN - (rng.randrange(1, 40) if rng.random() < 0.02 else 0)
        l2 = READLEN - (rng.randrange(1, 40) if rng.random() < 0.02 else 0)
        r1 = mutate((frag + A1 + draw(READLEN))[:l1])
        r2 = mutate((rc(frag) + A2 + draw(READLEN))[:l2])
        q1 = bytes(rng.choice((70, 70, 70, 58, 44)) for _ in range(len(r1)))
        q2 = bytes(rng.choice((70, 70, 70, 58, 44)) for _ in range(len(r2)))
        f1.write(b"@rd%d/1\n%b\n+\n%b\n" % (i, r1, q1))
        f2.write(b"@rd%d/2\n%b\n+\n%b\n" % (i, r2, q2))
print("wrote", N, "pairs")
