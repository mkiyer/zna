"""Force ties and check WHICH shift the kernel picks.

Ties need two shifts with identical (n, d): score = n*MW - d*STEP and MW/STEP is
not a simple ratio, so equal score implies equal n and equal d. Equal n at
different s happens on the two flanks (s = len1-n and s = n-len2) -- which for a
PERIODIC sequence both align perfectly. Random sequence essentially never ties,
which is why the previous sweep found none.
"""
import random, sys
from zna.merge.overlap import FORWARD, NO_OVERLAP, find_overlap, reverse_complement, score_weights
from zna.merge.fastqio import read_pairs

MW, MMW = score_weights(0.01); STEP = MW + MMW

def exhaustive(s1, s2rc, t=8.0):
    len1, len2 = len(s1), len(s2rc)
    out = []
    for s in range(-(len2 - 1), len1):
        lo, hi = max(s, 0), min(len1, s + len2)
        n = hi - lo
        if n <= 0: continue
        off = lo - s
        d = 0
        for k in range(n):
            d += s1[lo + k] != s2rc[off + k]
        sc = n * MW - d * STEP
        if sc >= t: out.append((s, n, d, sc))
    return out

def kernel_s(s1, s2rc):
    direction, shift, olen, diff, score = find_overlap(s1, s2rc, 8.0, 0.01)
    if direction == NO_OVERLAP: return None
    return (shift if direction == FORWARD else -shift), olen, diff

def probe(s1, s2rc, label, tally):
    allsc = exhaustive(s1, s2rc)
    got = kernel_s(s1, s2rc)
    if not allsc:
        assert got is None, f"{label}: kernel {got}, exhaustive none"; return
    top = max(t[3] for t in allsc)
    ties = sorted([t for t in allsc if abs(t[3] - top) < 1e-9])
    tally["cases"] += 1
    if len(ties) > 1:
        tally["tied"] += 1
        tally["max_tie"] = max(tally["max_tie"], len(ties))
    want = max(ties, key=lambda t: (t[1], -t[0]))          # max n, then min s
    want_maxs = max(ties, key=lambda t: (t[1], t[0]))      # max n, then MAX s (the rival rule)
    if got != (want[0], want[1], want[2]):
        tally["bad_min_s"] += 1
        if tally["bad_min_s"] <= 4:
            print(f"  {label}: kernel {got} | min-s rule {want[:3]} | max-s rule {want_maxs[:3]}"
                  f" | {len(ties)} tied at {top:.3f}")
    if len(ties) > 1 and got == (want_maxs[0], want_maxs[1], want_maxs[2]) and want != want_maxs:
        tally["matches_max_s"] += 1

T = dict(cases=0, tied=0, bad_min_s=0, matches_max_s=0, max_tie=0)

# --- 1. perfect tandem repeats: every flank pair ties exactly ---------------
for period in ("CA", "AT", "CAG", "AAAC", "ACGT"):
    for reps in range(6, 40):
        seq = (period * 100)[: len(period) * reps]
        probe(seq, seq, f"rep{period}x{reps}", T)          # s1 == s2rc, fully periodic
print(f"tandem repeats: {T}")

# --- 2. homopolymers: maximal tie degeneracy -------------------------------
T2 = dict(cases=0, tied=0, bad_min_s=0, matches_max_s=0, max_tie=0)
for L in range(8, 60):
    probe(b"A" * L, b"A" * L, f"polyA{L}", T2)
print(f"homopolymers:   {T2}")

# --- 3. unequal lengths -> a real plateau, plus periodic content -----------
T3 = dict(cases=0, tied=0, bad_min_s=0, matches_max_s=0, max_tie=0)
for l1 in range(20, 46, 3):
    for l2 in range(20, 46, 3):
        probe((b"CA" * 40)[:l1], (b"CA" * 40)[:l2], f"plateau{l1}x{l2}", T3)
        probe((b"ACG" * 30)[:l1], (b"ACG" * 30)[:l2], f"plat3-{l1}x{l2}", T3)
print(f"plateaus:       {T3}")

# --- 4. real full-length 2x150 reads ---------------------------------------
T4 = dict(cases=0, tied=0, bad_min_s=0, matches_max_s=0, max_tie=0)
n = 0
for r1, r2 in read_pairs(sys.argv[1], sys.argv[2], 1):
    probe(r1[1], reverse_complement(r2[1]), f"real{n}", T4)
    n += 1
    if n >= 800: break
print(f"real 2x150 (n={n}): {T4}")

tot_bad = T["bad_min_s"] + T2["bad_min_s"] + T3["bad_min_s"] + T4["bad_min_s"]
tot_tied = T["tied"] + T2["tied"] + T3["tied"] + T4["tied"]
print(f"\nties exercised: {tot_tied}   violations of (max score, max n, min s): {tot_bad}")
print("VERDICT:", "TIE-BREAK RULE HOLDS" if tot_bad == 0 else "RULE VIOLATED")
