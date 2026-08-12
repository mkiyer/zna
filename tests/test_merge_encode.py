"""End-to-end: `zna merge` -> `zna encode` -> read back, checking the geometry.

This is the test that sees BOTH sides of the merge/encode seam, which is why the
contract is tested here (docs/SOS_EOS_ENCODING_PLAN.md §9.5). Two bugs would have been
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
           "--npolicy", "drop",
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
               "--npolicy", "drop", "--strand-normalize", *STRANDS["ISF"],
               "--label-defs", str(label_defs(tmp_path)), "-o", str(out),
               str(merged_fastq)]
        assert subprocess.run(cmd, capture_output=True, text=True).returncode == 0
        singles = [r for r in read_back(out) if not r[2]]
        assert singles and all(r[5] != r[6] for r in singles)   # one end only
