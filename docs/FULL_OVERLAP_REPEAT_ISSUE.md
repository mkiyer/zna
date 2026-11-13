# Issue: an exact reverse-complement pair is not proof of a full overlap

Status: **fixed in 0.7.0** (option A, the length gate; [issue #1](https://github.com/mkiyer/zna/issues/1)). Reported 2026-11-13 from khorana's
aberration-library build on khorana3 (zna 0.6.0; `src/zna/cli.py` identical at main
0a19141). The operator's four decisions (§6, §8):

1. **L is declared**, as `zna encode --read-length L`. A read longer than L is an error,
   and a run in which no read reached L warns.
2. **Without L, no pair is marked full.** The summary counts the exact pairs left
   one-ended. That changes the default, hence a minor release.
3. **The merger is unchanged** (policy `zna-merge-0.6`). Zero offset is its most reliable
   merge, and the repeat case there is part of the homologous-ends residual METHODS §3.4
   documents. `--read-length` is refused with `--merge-pairs`.
4. **Files already written are re-encoded from source.** No ZNA → ZNA repair exists. The
   prologue's new `full_fragment` record (METHODS §3.8) is absent on every 0.6 file, so a
   consumer that requires it refuses them.

The same release stops a record `--npolicy trim3` cut from being marked full, which 0.6 did
for declared merged reads. The rule is METHODS §3.4. The rest of this report is kept as it
was written.

## 1. Summary

`zna encode` flags an unmerged pair as a full overlap (`IS_FULL_FRAGMENT` on both mates)
whenever R1 is the exact reverse complement of R2 (`_full_fragment_flags`,
`_is_reverse_complement`, `src/zna/cli.py`). METHODS §3.4 states the premise: "Two mates
covering the identical interval are exact reverse complements." The converse does not hold.
A fragment **longer** than the read whose first and last read-length bases are the same
sequence also gives mates that are exact reverse complements. Tandem repeats do this whenever
the insert minus the read length is a whole number of repeat units.

Both mates then read `(has_start, has_end) = (True, True)`. Each claims a fragment boundary at
its inner edge, where the fragment does not end. A consumer that places fragment-end
supervision from `ENDS_BY_FLAG` learns an end inside the repeat. This is the one direction the
format's own rule forbids: METHODS §3.6 calls under-labelling "never wrong", and this is
over-labelling.

Merging is not involved. The pair is kept as two mates, and the flag alone is wrong. The merge
path has the same ambiguity (§4).

## 2. The condition

For a stranded pair under the dUTP convention, with fragment F and read length L:
R2 = F[0, L) and R1 = rc(F[|F| − L, |F|)). So

    R1 == rc(R2)   ⟺   F[0, L) == F[|F| − L, |F|)

The fragment's two ends carry the same L bases in the same orientation: a **direct** repeat
spanning both ends. It is not an inverted repeat; khorana's FRAGMENT_MODEL.md uses that word
for this case. When |F| ≤ L, the condition is a true full overlap. When |F| > L, it is a
coincidence of the sequence. In a tandem array of period p, it holds for every |F| = L + k·p
that lies inside the array.

## 3. Reproduction

Four synthetic pairs and one real pair, encoded with the flags khorana's simulator uses. The
real pair is from GENCODE's UBC (ENST00000339647.6, fragment 1337–1715 of the transcript,
whose ubiquitin units are 228 bases long).

```python
# python -I repro.py OUTDIR   (zna 0.6.0 on PATH beside this interpreter)
import random, subprocess, sys
from pathlib import Path
from zna import ZnaReader

COMP = str.maketrans('ACGT', 'TGCA')
rc = lambda s: s.translate(COMP)[::-1]

def pair(fragment, L=150):            # dUTP: R1 antisense; reads clipped to the molecule
    n = min(L, len(fragment))
    return rc(fragment[-n:]), fragment[:n]

out = Path(sys.argv[1]); rng = random.Random(7)
seq = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
unit = seq(228)
cases = {'tandem_repeat_378': (unit * 3)[:378], 'short_insert_120': seq(120),
         'ordinary_300': seq(300), 'insert_equals_read_150': seq(150)}
with open(out / 'r1.fq', 'w') as f1, open(out / 'r2.fq', 'w') as f2:
    for name, frag in cases.items():
        r1, r2 = pair(frag)
        f1.write(f'@{name}/1\n{r1}\n+\n{"I" * len(r1)}\n')
        f2.write(f'@{name}/2\n{r2}\n+\n{"I" * len(r2)}\n')
subprocess.run([str(Path(sys.executable).parent / 'zna'), 'encode', '--strand-normalize',
                '--strand-specific', '--seed', '1', '-o', str(out / 'repro.zna'),
                '-q', str(out / 'r1.fq'), str(out / 'r2.fq')], check=True)
with open(out / 'repro.zna', 'rb') as s:
    for i, (sq, paired, is_r1, is_r2, hs, he) in enumerate(ZnaReader(s).records(with_ends=True)):
        print(list(cases)[i // 2], 'R1' if is_r1 else 'R2', (hs, he))
```

zna 0.6.0, with the real UBC pair appended:

| case | insert | read | R1 == rc(R2) | mate | `(has_start, has_end)` | truth |
|---|---:|---:|---|---|---|---|
| tandem_repeat_378 | 378 | 150 | yes | R1 | **(True, True)** | (False, True) |
| | | | | R2 | **(True, True)** | (True, False) |
| short_insert_120 | 120 | 120 | yes | R1, R2 | (True, True) | (True, True) |
| ordinary_300 | 300 | 150 | no | R1 / R2 | (False, True) / (True, False) | same |
| insert_equals_read_150 | 150 | 150 | yes | R1, R2 | (True, True) | (True, True) |
| UBC, real | 378 | 150 | yes | R1 | **(True, True)** | (False, True) |
| | | | | R2 | **(True, True)** | (True, False) |

`zna encode --merge-pairs`, with and without `--adapter-trimmed`, merges both repeat pairs into
one 150-base record with `(True, True)`. That asserts a 150-base molecule where the true one is
378 bases.

## 4. How often, and what a fix trades

Measured on khorana's simulated aberration library: GENCODE v50 basic at NA12878's abundances,
19,245,880 pairs of 2 × 150 stranded reads clipped to the molecule, fragments 250 ± 50 (minimum
50), 0.5% substitutions and 0.03% insertions and deletions. Each pair was tallied by its true
insert and by zna 0.6's test (R1 == rc(R2), equal length):

| true insert | pairs | exact reverse complements (0.6 flags full) |
|---|---:|---:|
| < 150 | 455,223 | 106,929, correct |
| = 150 | 22,018 | 3,994, correct |
| > 150 | 18,768,639 | **1, false** (the UBC pair) |

`zna inspect --counts` on the encoded file shows 221,846 full-fragment records (110,923
pairs). The FASTQ tally gives 110,924; the one-pair gap was not traced (a pair cut by
`--npolicy trim3` at a genome N would explain it).

- **The false case is rare but systematic.** khorana's backlog (`zna-inverted-repeat-overlap`)
  reports 4 of 100M records at genome scale and 6 of 173.6M in another GENCODE v50 basic
  library. Those numbers were not re-measured here. A second uniform-abundance GENCODE library
  has an ACAN pair (ENST00000560601.4, 3691–3955) of the same kind. Both cases seen are in
  tandem-repeat genes. Other repeat genes (UBB, the mucins) would be expected but were not
  checked. UBC is among the most highly expressed genes, so its labels recur across libraries.
- **METHODS §3.4's measurement cannot see it.** "10051/10051 detected with 0/9949 false
  positives" was measured on unrelated sequence, where an exact full-length match does not
  occur by chance.
- **Most true full overlaps are already missed.** With read errors, both mates must be
  error-free to match exactly: 110,923 of 477,241 true inserts of 150 or less (23%) are
  flagged here. Under-labelling is the existing, accepted cost of the exact test.

## 5. Where it comes from

The test uses only the two reads, and a full-length pair cannot tell an insert of exactly L
from a direct repeat. Positive evidence of a full overlap exists only when the reads are
**shorter** than the sequenced read length. Then the insert ended inside the read, and the
read ran into adapter, which was trimmed (or the simulator clipped it). An exact match at full
read length is ambiguous, and on the numbers above it is true 3,994 times for every false one.

## 6. Fix options

**A. A length gate (recommended).** Flag a pair full only when the mates are exact reverse
complements **and shorter than the sequenced read length L**. Leave an exact pair at full
length one-ended, each mate with the edge it was read from, which is always true.

- It removes every false case in the table. A repeat cannot produce mates shorter than L
  unless trimming made them so, and then only rarely.
- Cost: inserts of exactly L lose their second end. Here that is 3,994 of 110,924 flagged
  pairs (3.6%), each falling back to the one-ended labels its 18,024 error-bearing neighbours
  already have. This is METHODS §3.6's "under-labelled, never wrong".
- L must come from somewhere. Options, for the implementer to choose:
  - an explicit `--read-length L`, which is predictable and fits `--npolicy`'s "a flag a user
    can predict";
  - the longest mate in the input, which needs a pass before the first record is flagged;
  - the header, if the format records it.

  Without L, either keep 0.6's rule (no silent change of output) or never flag at full length.
  The default is the decision to make. Report the count either way, as `--npolicy` does:
  "N exact pairs at full read length left one-ended".
- The merge path (`--merge-pairs`, `zna merge`) has the same ambiguity at shift 0 with equal
  full-length mates. Whether the gate applies there too (keep such pairs unmerged) belongs to
  the same decision; METHODS §3.4 already lists homologous fragment ends among the residual
  wrong merges.

**B. Report and document only.** Keep the rule, add the count of full-length exact pairs to
the encode summary, and correct METHODS §3.4's premise. This is cheapest, and the false ends
stay.

**C. Never flag unmerged pairs.** Drop detection, leaving full overlaps to the merger and
`--treat-unpaired-as-merged`. This is simplest, but every short-insert pair loses its second
end (here 110,923 pairs).

## 7. Tests to add

- A tandem-repeat fragment of insert L + p, with full-length mates that are exact reverse
  complements: under A, both mates are one-ended (`(False, True)` and `(True, False)` after
  stranded normalization).
- A true short insert (reads shorter than L): both mates stay `(True, True)`.
- An insert of exactly L: under A, one-ended. The test pins the documented cost.
- The encode summary counts the full-length exact pairs left one-ended.
- The same cases through `--merge-pairs` and with `--adapter-trimmed`, whatever the merge
  decision is.

## 8. Release

The fix changes the flag byte of some records written by `zna encode`. An opt-in flag with
unchanged defaults fits a patch release. A changed default is a behaviour change for every
consumer that places end supervision, and fits a minor release with a CHANGELOG entry. Update
METHODS §3.4 (the premise and the measurement) and §3.6's row for full-overlap pairs.

Nothing in the records says which rule set bit 4, as nothing in them said which merge policy
made them (§3.8). So the fix should write the rule into the prologue the same way: an optional
key naming the full-overlap rule and the L it used, carried verbatim by shuffle and ZNA → ZNA
re-encode, absent meaning unknown (the 0.6 rule). A consumer then checks the key for equality,
not `writer_version`.

Files already written cannot be repaired by copying, because copy paths carry the flag byte
verbatim (§3.5). Either re-encode from FASTQ, or provide a re-derivation for ZNA input. The
re-derivation is possible because both mates' sequences and lengths are stored: given L,
clear bit 4 on paired records at full length. khorana's simulated libraries keep their FASTQ.
Its real libraries, encoded with 0.6, would need the re-derivation or their source FASTQ.

## 9. Not this issue

khorana's audit of the same library also flags 13 reads that cross a base the reference genome
itself spells N (IPMK and MT-RNR2's chrM 3107). `--npolicy trim3` cuts each read at its first N
and keeps `[0, first N)`, exactly as documented. Base 0 stays a true boundary, and no flag is
wrong. That is a khorana audit expectation to fix, not a zna defect.

## 10. References

- `src/zna/cli.py`: `_is_reverse_complement`, `_fragment_units`, `_full_fragment_flags`.
- METHODS §3.4 `IS_FULL_FRAGMENT`, §3.5 the projection, §3.6 where records land, §1.7
  `--adapter-trimmed`.
- khorana: `docs/dev/FRAGMENT_MODEL.md` ("Endpoints without the merger"), backlog item
  `zna-inverted-repeat-overlap`, and the library's audit
  (`/data/khorana/sim/20261113-gencode50-basic-aberrations-na12878/`, on khorana3).
