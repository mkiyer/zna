"""End-to-end: `zna merge` -> `zna encode` -> read back, checking the geometry.

This is the test that sees BOTH sides of the merge/encode seam, which is why the
contract is tested here (docs/METHODS.md, "Fragment geometry"). Two bugs would have been
caught by it and were not caught by either half's own suite:

  * `zna encode --label-defs` silently ignored `--interleaved` and flagged every
    record unpaired, so a downstream model gave every unmerged mate a fragment-end
    token on an edge that is really a read-length truncation;
  * a merged read that is not actually its whole fragment would make
    `--treat-unpaired-as-merged` a lie.

The contract, per record: **which edges of the stored sequence are true fragment
boundaries?** A merged read spans its fragment, so both. A mate of an unmerged pair
reaches a real boundary at base 0 only, so exactly one — which edge depends on whether
strand normalization reverse-complemented it, and that is precisely what
``records(with_ends=True)`` resolves.

Both halves are now local, so nothing here is conditional: the merge runs in process
and `zna encode` runs out of process against this interpreter's own `zna.cli`, with
the flag set a production pipeline uses.
"""
import gzip
import random
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

from zna.core import ZnaReader

from zna.merge import cli
from zna.merge.overlap import reverse_complement as rc

READLEN = 150
MIN_READ_LENGTH = 40
ADAPTER1 = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER2 = b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

# The label block a pipeline hands to `zna encode`, written locally rather than read
# from the pipeline repo. The key that names the SAM tag is `tag:` — an unknown key
# such as `key:` is silently ignored and the label falls back to matching on `name:`,
# which never appears in a header, so every record gets the `missing` sentinel. That
# failure is invisible to `zna inspect --counts` (it reports the label *schema* and the
# flags byte, never a label *value*), which is why
# `test_label_defs_file_really_extracts_the_sam_tags` reads the values back instead.
LABEL_DEFS_YAML = textwrap.dedent("""\
    labels:
      - name: transcript_index
        description: transcript index (-1 if unassigned)
        tag: ZI
        type: i
        missing: -1
      - name: gene_index
        description: gene index (-1 if unassigned)
        tag: ZJ
        type: i
        missing: -1
      - name: assignment_flags
        description: properties of the assigned transcript
        tag: ZF
        type: C
        missing: 0
    """)


def label_defs(tmp_path) -> Path:
    path = tmp_path / "label_defs.yaml"
    if not path.exists():
        path.write_text(LABEL_DEFS_YAML)
    return path


# --------------------------------------------------------------------------- #
# input simulation
# --------------------------------------------------------------------------- #

def _draw(rng, n):
    return bytes("".join(rng.choices("ACGT", k=n)), "ascii")


def simulate(rng, n_each=30):
    """Fragments spanning all three merge outcomes, keyed by fragment id.

    Inserts below the read length read through into adapter (merge), inserts a little
    above it overlap (merge), and long inserts do not overlap at all (kept pair).
    """
    frags = {}
    for i in range(n_each):
        frags[i] = _draw(rng, rng.randrange(60, 140))            # read-through
        frags[n_each + i] = _draw(rng, rng.randrange(170, 280))  # overlapping
        frags[2 * n_each + i] = _draw(rng, rng.randrange(330, 420))  # disjoint
    return frags


def write_fastqs(tmp_path, frags, rng):
    """Write R1/R2 as `samtools fastq -N -T "*"` would: /1,/2 suffix + SAM tags."""
    in1, in2 = tmp_path / "r1.fastq.gz", tmp_path / "r2.fastq.gz"
    with gzip.open(in1, "wb") as f1, gzip.open(in2, "wb") as f2:
        for fid, frag in frags.items():
            # Draw the two mate lengths INDEPENDENTLY, as fastp's 3' quality trim does.
            # With both mates at READLEN the geometry assertions below cannot fail: the
            # merged read is truncated only when len(R1) < len(R2), so an equal-length
            # fixture is structurally blind to the whole class.
            l1 = READLEN - (rng.randrange(1, 60) if rng.random() < 0.4 else 0)
            l2 = READLEN - (rng.randrange(1, 60) if rng.random() < 0.4 else 0)
            r1 = (frag + ADAPTER1 + _draw(rng, READLEN))[:l1]
            r2 = (rc(frag) + ADAPTER2 + _draw(rng, READLEN))[:l2]
            tags = b"\tZI:i:%d\tZJ:i:0\tZF:C:0" % fid
            f1.write(b"@f%d/1%b\n%b\n+\n%b\n" % (fid, tags, r1, b"I" * len(r1)))
            f2.write(b"@f%d/2%b\n%b\n+\n%b\n" % (fid, tags, r2, b"I" * len(r2)))
    return in1, in2


def run_merge(tmp_path, in1, in2):
    out = tmp_path / "merged.fastq.gz"
    args = cli.build_parser().parse_args([
        "--in1", str(in1), "--in2", str(in2), "--out", str(out),
        "--min-read-length", str(MIN_READ_LENGTH), "--threads", "1", "-q",
    ])
    stats = cli.run(args)
    return out, stats


def parse_merged(path):
    """Group the mixed interleaved stream by fragment id -> list of (header, seq)."""
    by_id = {}
    with gzip.open(path, "rb") as fh:
        lines = fh.read().splitlines()
    for i in range(0, len(lines), 4):
        header, seq = lines[i][1:], lines[i + 1]
        idtok = header.split(b"\t", 1)[0]                 # "f12/1", or "f12" when merged
        if idtok.endswith(b"/1") or idtok.endswith(b"/2"):
            idtok = idtok[:-2]
        by_id.setdefault(int(idtok[1:]), []).append((header, seq))
    return by_id


def encode(tmp_path, fastq, strand_args, name="reads.zna"):
    """`zna encode` with exactly the flag set a production pipeline uses.

    Invoked as `python -m zna.cli` rather than via a `zna` on PATH, so the CLI under
    test is unambiguously this interpreter's.
    """
    out = tmp_path / name
    cmd = [sys.executable, "-m", "zna.cli", "encode", "--interleaved",
           "--shuffle", "--seed", "7",
           "--seq-len-bytes", "2",
           "--npolicy", "trim3",
           "--treat-unpaired-as-merged",
           "--read-group", "lib1",
           "--description", "test lib1",
           "--strand-normalize",
           *strand_args,
           "--label-defs", str(label_defs(tmp_path)),
           "-o", str(out), str(fastq)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    return out


def read_back(path):
    """[(fragment_id, seq, is_paired, is_read1, is_read2, has_start, has_end)]"""
    out = []
    with open(path, "rb") as fh:
        for seq, is_paired, is_read1, is_read2, has_start, has_end, labels in \
                ZnaReader(fh).records(with_ends=True):
            out.append((labels[0], seq.encode(), is_paired, is_read1, is_read2,
                        has_start, has_end))
    return out


# --------------------------------------------------------------------------- #

STRANDS = {
    # the two shapes a stranded/unstranded RNA-seq pipeline emits
    "ISF": ["--strand-specific", "--read1-sense", "--read2-antisense"],
    "IU": [],
}


@pytest.fixture(scope="module")
def corpus(tmp_path_factory):
    tmp_path = tmp_path_factory.mktemp("e2e")
    rng = random.Random(20260811)
    frags = simulate(rng)
    in1, in2 = write_fastqs(tmp_path, frags, rng)
    merged_fastq, stats = run_merge(tmp_path, in1, in2)
    by_id = parse_merged(merged_fastq)
    return tmp_path, frags, by_id, merged_fastq, stats


class TestMergeToZna:
    def test_merge_produced_all_three_outcomes(self, corpus):
        """Guard the fixture itself: a corpus with no kept pairs would prove nothing."""
        _tmp, frags, by_id, _fq, stats = corpus
        singles = [f for f, r in by_id.items() if len(r) == 1]
        pairs = [f for f, r in by_id.items() if len(r) == 2]
        assert len(singles) >= 30 and len(pairs) >= 20
        assert stats["merged"] == len(singles)
        assert stats["trimmed_pairs"] + stats["kept_pairs"] >= len(pairs)
        # Every emitted record belongs to a pair or a single — never an orphan.
        assert all(len(r) in (1, 2) for r in by_id.values())

    @pytest.mark.parametrize("libtype", sorted(STRANDS))
    def test_every_edge_reported_is_a_true_fragment_boundary(self, corpus, libtype):
        """The load-bearing assertion: merged records report BOTH ends, mates exactly
        one, and every reported edge really is a fragment terminus."""
        tmp_path, frags, by_id, merged_fastq, _stats = corpus
        znafile = encode(tmp_path, merged_fastq, STRANDS[libtype], f"{libtype}.zna")
        records = read_back(znafile)

        assert len(records) == sum(len(r) for r in by_id.values())

        seen_pairs = {}
        for fid, seq, is_paired, is_read1, is_read2, has_start, has_end in records:
            frag = frags[fid]
            forms = (frag, rc(frag))          # strand normalization may have flipped it
            if not is_paired:
                # A merged read spans its whole fragment: both edges are real.
                assert has_start and has_end, (libtype, fid)
                assert seq in forms, (libtype, fid, len(seq), len(frag))
                assert not is_read1 and not is_read2
            else:
                # A mate reaches a real boundary at base 0 only — exactly one edge.
                assert has_start != has_end, (libtype, fid)
                if has_start:                 # left edge real -> stored 5'->3' from it
                    assert any(f.startswith(seq) for f in forms), (libtype, fid)
                else:                         # right edge real -> stored reverse-complemented
                    assert any(f.endswith(seq) for f in forms), (libtype, fid)
                assert is_read1 != is_read2
                seen_pairs.setdefault(fid, []).append(is_read1)

        # The regression this exists for: with the labeled-encoding bug every record
        # came back unpaired, so this dict would be empty and every mate would have
        # claimed a fragment end on its truncated edge.
        expected_pairs = {f for f, r in by_id.items() if len(r) == 2}
        assert set(seen_pairs) == expected_pairs
        for fid, mates in seen_pairs.items():
            assert sorted(mates) == [False, True], (libtype, fid)   # one R1, one R2

    def test_label_defs_file_really_extracts_the_sam_tags(self, corpus):
        """A label definition must name its tag with `tag:`, not `key:`.

        zna's schema ignores unknown keys and falls back to matching on the label
        `name:`, which never appears in a header — so a typo there does not fail, it
        silently writes the `missing` sentinel into every record.
        """
        tmp_path, frags, _by_id, merged_fastq, _stats = corpus
        znafile = encode(tmp_path, merged_fastq, STRANDS["ISF"], "labels.zna")
        ids = {r[0] for r in read_back(znafile)}
        assert -1 not in ids, "ZI label came back as the missing sentinel"
        assert ids == set(frags)

    def test_full_fragment_flag_is_set_only_on_merged_records(self, corpus):
        """`zna inspect --counts` is the field-level version of this check."""
        tmp_path, _frags, by_id, merged_fastq, _stats = corpus
        znafile = encode(tmp_path, merged_fastq, STRANDS["ISF"], "flags.zna")
        both = sum(1 for r in read_back(znafile) if r[5] and r[6])
        assert both == sum(1 for r in by_id.values() if len(r) == 1)

    def test_without_the_flag_merged_reads_lose_their_second_end(self, corpus):
        """Why a merge-fed encode must pass --treat-unpaired-as-merged: a merged read
        and a genuine single-end read are byte-identical, so zna cannot infer this."""
        tmp_path, _frags, _by_id, merged_fastq, _stats = corpus
        out = tmp_path / "noflag.zna"
        cmd = [sys.executable, "-m", "zna.cli", "encode", "--interleaved",
               "--seq-len-bytes", "2",
               "--npolicy", "trim3", "--strand-normalize", *STRANDS["ISF"],
               "--label-defs", str(label_defs(tmp_path)), "-o", str(out),
               str(merged_fastq)]
        assert subprocess.run(cmd, capture_output=True, text=True).returncode == 0
        singles = [r for r in read_back(out) if not r[2]]
        assert singles and all(r[5] != r[6] for r in singles)   # one end only


# --------------------------------------------------------------------------- #
# per-record provenance: the only thing about a record's history that survives
# encoding, since ZNA does not store headers
# --------------------------------------------------------------------------- #

def write_fastqs_with_no_calls(tmp_path, frags, rng):
    """As :func:`write_fastqs`, but with no-calls injected 3'-biased.

    The simulator emits no ``N`` at all, so nothing downstream of it can exercise the
    N policy or the provenance it produces. Real no-calls cluster at the 3' end (later
    cycles are the worse ones), and that placement matters here rather than being
    decoration: an ``N`` inside the overlap gets rescued from the mate and one past it
    survives to meet the policy, so a 3' bias is what makes both fates occur.
    """
    in1, in2 = tmp_path / "n1.fastq.gz", tmp_path / "n2.fastq.gz"
    with gzip.open(in1, "wb") as f1, gzip.open(in2, "wb") as f2:
        for fid, frag in frags.items():
            r1 = bytearray((frag + ADAPTER1 + _draw(rng, READLEN))[:READLEN])
            r2 = bytearray((rc(frag) + ADAPTER2 + _draw(rng, READLEN))[:READLEN])
            for read in (r1, r2):
                if rng.random() < 0.35:
                    lo = int(len(read) * 0.5)              # 3'-biased, as instruments are
                    read[rng.randrange(lo, len(read))] = ord("N")
            tags = b"\tZI:i:%d" % fid
            f1.write(b"@f%d/1%b\n%b\n+\n%b\n" % (fid, tags, bytes(r1), b"I" * len(r1)))
            f2.write(b"@f%d/2%b\n%b\n+\n%b\n" % (fid, tags, bytes(r2), b"I" * len(r2)))
    return in1, in2


@pytest.fixture(scope="module")
def n_corpus(tmp_path_factory):
    """Merge an N-bearing library, then encode it with the provenance column declared."""
    tmp_path = tmp_path_factory.mktemp("prov")
    rng = random.Random(20260813)
    frags = simulate(rng)
    in1, in2 = write_fastqs_with_no_calls(tmp_path, frags, rng)
    merged, stats = run_merge(tmp_path, in1, in2)
    out = tmp_path / "prov.zna"
    cmd = [sys.executable, "-m", "zna.cli", "encode", "--interleaved",
           "--treat-unpaired-as-merged", "--seq-len-bytes", "2",
           "--npolicy", "trim3",
           "--label", "prov:C:ZN", "--label", "zi:i:ZI",
           "-o", str(out), str(merged)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    recs = []
    with open(out, "rb") as fh:
        for seq, flags, labels in ZnaReader(fh).copy_records():
            recs.append((seq, flags, labels))
    return stats, recs


class TestProvenanceReachesTheCorpus:
    """The `ZN:i:<bits>` tag is provenance that survives `zna encode`.

    Header tokens are for a human reading the intermediate FASTQ; they vanish at encode
    time, and under `--merge-pairs` (0.5.0) there is no intermediate FASTQ at all. A
    declared `C` column is the only per-record provenance that reaches a model, so what
    matters is not that the tag is written but that it arrives — through the ordinary
    `--label` path, with no provenance-specific code in the encoder.
    """

    def test_the_fixture_actually_produced_provenance(self, n_corpus):
        """Guard the fixture: a library with no no-calls would prove nothing here."""
        stats, recs = n_corpus
        assert stats["npolicy_bases"] > 0 and stats["n_rescued_from_mate"] > 0
        assert any(r[2][0] for r in recs), "no record carries any provenance bit"

    def test_an_untouched_record_reports_nothing(self, n_corpus):
        """An absent tag must resolve to 0 through the label machinery's own missing
        path — that is what makes declaring the column optional rather than a schema
        change to every labeled file."""
        _stats, recs = n_corpus
        assert sum(1 for r in recs if r[2][0] == 0) > 0, "every record was touched"

    def test_the_trimmed_bit_counts_the_trimmed_pairs(self, n_corpus):
        """PROV_TRIMMED is the one bit with no other home: a trimmed pair is emitted as
        an ordinary pair, and nothing in the ZNA flag byte distinguishes it from a pair
        that was kept whole. Both mates of every surviving trimmed pair carry it."""
        stats, recs = n_corpus
        from zna.merge._pymerge import PROV_TRIMMED
        n_trimmed_recs = sum(1 for r in recs if r[2][0] & PROV_TRIMMED)
        assert n_trimmed_recs == 2 * stats["trimmed_pairs"], (
            f"{n_trimmed_recs} records carry PROV_TRIMMED but "
            f"{stats['trimmed_pairs']} pairs were trimmed")

    def test_trim3_never_claims_a_substitution(self, n_corpus):
        """The two N policies are mutually exclusive per run, and the bits say which
        one ran. Under trim3 nothing was invented, so PROV_NSUBBED must be absent —
        this is the bit that would lie about a corpus containing made-up bases."""
        _stats, recs = n_corpus
        from zna.merge._pymerge import PROV_NSUBBED, PROV_NTRIMMED
        assert not any(r[2][0] & PROV_NSUBBED for r in recs)
        assert any(r[2][0] & PROV_NTRIMMED for r in recs), "trim3 cut nothing"

    def test_no_bit_outside_the_defined_vocabulary(self, n_corpus):
        """A byte with an undefined bit set means the two sides disagree about what the
        column means, which is worse than an absent column."""
        _stats, recs = n_corpus
        from zna.merge._pymerge import (PROV_NSUBBED, PROV_NTRIMMED, PROV_RESCUED,
                                        PROV_TRIMMED)
        known = PROV_TRIMMED | PROV_RESCUED | PROV_NTRIMMED | PROV_NSUBBED
        assert all(r[2][0] & ~known == 0 for r in recs)

    def test_the_other_labels_are_unaffected(self, n_corpus):
        """Provenance is added beside the existing tags, never instead of them: the ZI
        tag merge passed through must still read back on every record."""
        _stats, recs = n_corpus
        assert all(r[2][1] >= 0 for r in recs)
        assert len({r[2][1] for r in recs}) > 1, "the ZI label collapsed to one value"


class TestEncodeReportsItsNpolicy:
    """`zna merge` has always said what its N policy did; `zna encode` did not.

    The same policy on the same library was accountable on one side of the seam and
    silent on the other, and the failure it guards against is silent by nature: one dark
    cycle can make a policy consume most of a library while the run still ends with
    "Done."
    """

    @staticmethod
    def _encode(tmp_path, fastq, npolicy, name, extra=()):
        out = tmp_path / name
        cmd = [sys.executable, "-m", "zna.cli", "encode", "--npolicy", npolicy,
               *extra, "-o", str(out), *[str(f) for f in fastq]]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        assert proc.returncode == 0, proc.stderr
        return proc.stderr

    def test_random_reports_exactly_the_no_calls_it_substituted(self, tmp_path):
        """The count is checkable against the input rather than against itself: under
        `random` every no-call becomes exactly one substituted base, so the reported
        number must equal the number of Ns in the file."""
        rng = random.Random(11)
        frags = simulate(rng, n_each=10)
        in1, _in2 = write_fastqs_with_no_calls(tmp_path, frags, rng)
        with gzip.open(in1, "rb") as fh:
            n_calls = sum(line.count(b"N")
                          for i, line in enumerate(fh.read().splitlines())
                          if i % 4 == 1)
        assert n_calls > 0, "the fixture injected no no-calls"
        err = self._encode(tmp_path, [in1], "random", "r.zna")
        assert f"substituted {n_calls} bases" in err, err

    def test_trim3_reports_bases_removed_not_no_calls_found(self, tmp_path):
        """trim3 cuts from the first no-call to the 3' end, so it removes at least as
        many bases as there are no-calls, and reporting the no-call count instead would
        understate the loss — which is the whole point of the line."""
        rng = random.Random(12)
        frags = simulate(rng, n_each=10)
        in1, _in2 = write_fastqs_with_no_calls(tmp_path, frags, rng)
        with gzip.open(in1, "rb") as fh:
            lines = fh.read().splitlines()
        seqs = [lines[i] for i in range(1, len(lines), 4)]
        expect = sum(len(s) - s.find(b"N") for s in seqs if b"N" in s)
        n_recs = sum(1 for s in seqs if b"N" in s)
        err = self._encode(tmp_path, [in1], "trim3", "t.zna")
        assert f"removed {expect} bases in {n_recs} records" in err, err
        assert expect > n_recs, "the fixture never cut more than one base"

    def test_the_paired_path_counts_both_mates(self, tmp_path):
        """Two-file input goes through a different write loop than single-end, and an
        untouched loop reports zero on a library that is visibly being trimmed."""
        rng = random.Random(13)
        frags = simulate(rng, n_each=10)
        in1, in2 = write_fastqs_with_no_calls(tmp_path, frags, rng)
        one = self._encode(tmp_path, [in1], "random", "p1.zna")
        both = self._encode(tmp_path, [in1, in2], "random", "p2.zna")
        n_one = int(one.split("substituted ")[1].split()[0])
        n_both = int(both.split("substituted ")[1].split()[0])
        assert n_one > 0 and n_both > n_one, (n_one, n_both)

    def test_a_heavy_policy_warns(self, tmp_path):
        """Above 1% of encoded bases the line stops being informational: that is a run
        problem, not ordinary no-calls."""
        seqs = [b"ACGT" * 20 + b"N" + b"ACGT" * 20 for _ in range(50)]
        fq = tmp_path / "heavy.fastq"
        fq.write_bytes(b"".join(b"@r%d\n%b\n+\n%b\n" % (i, s, b"I" * len(s))
                                for i, s in enumerate(seqs)))
        err = self._encode(tmp_path, [fq], "trim3", "h.zna")
        assert "WARNING" in err and "% of encoded bases" in err, err

    def test_reencoding_a_zna_claims_no_policy_ran(self, tmp_path):
        """ZNA stores two bits per base, so a decoded record cannot contain a no-call
        and no policy ran. Reporting "removed 0 bases" would imply one had been."""
        rng = random.Random(14)
        frags = simulate(rng, n_each=10)
        in1, _in2 = write_fastqs_with_no_calls(tmp_path, frags, rng)
        self._encode(tmp_path, [in1], "trim3", "src.zna")
        err = self._encode(tmp_path, [tmp_path / "src.zna"], "trim3", "re.zna")
        assert "no-calls:" not in err, err
