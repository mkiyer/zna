"""The inflate-backend policy, and the three paths agreeing on every input.

ZNA can inflate a ``.gz`` three ways -- ISA-L, an external ``pigz``/``gzip``, or the
stdlib -- and which one is preferred depends on whether the caller reads on its own
thread (see :mod:`zna._gzip`). These tests pin two things that matter more than the
speed the choice exists for:

* **All three produce the same bytes.** A decompressor swapped in for throughput must
  not be able to change what gets encoded.
* **All three still refuse damaged input.** The subprocess path exists partly to check
  a decompressor's exit status, because a truncated ``.gz`` otherwise looks exactly
  like a short file and ``zna encode`` would write a truncated ZNA and report success.
  An in-process path has to be held to the same standard.

The ISA-L cases skip when :mod:`isal` is absent -- it is an optional extra -- and the
subprocess cases skip when no ``pigz``/``gzip`` binary is on PATH.
"""
from __future__ import annotations

import gzip
import shutil
import subprocess
import sys

import pytest

from zna import _gzip

HAVE_ISAL = _gzip.isal_module() is not None
HAVE_EXTERNAL = shutil.which("pigz") is not None or shutil.which("gzip") is not None

RECORDS = 4000
READ_LEN = 60


def _fastq_bytes(n=RECORDS, read_len=READ_LEN):
    """An interleaved paired FASTQ, deterministic without needing a seed argument."""
    out = []
    for i in range(n):
        # A rotating but non-degenerate sequence: 2-bit packing and zstd both behave
        # differently on constant input, and a constant file would hide a length bug.
        seq = "".join("ACGT"[(i + k * 7) % 4] for k in range(read_len))
        mate = 1 if i % 2 == 0 else 2
        out.append(f"@rd{i // 2}/{mate}\n{seq}\n+\n{'I' * read_len}\n")
    return "".join(out).encode()


@pytest.fixture(scope="module")
def plain_fastq(tmp_path_factory):
    p = tmp_path_factory.mktemp("gzb") / "in.fq"
    p.write_bytes(_fastq_bytes())
    return p


@pytest.fixture(scope="module")
def gz_fastq(tmp_path_factory):
    p = tmp_path_factory.mktemp("gzb") / "in.fq.gz"
    with gzip.open(p, "wb", compresslevel=1) as fh:
        fh.write(_fastq_bytes())
    return p


@pytest.fixture(scope="module")
def truncated_gz(tmp_path_factory):
    p = tmp_path_factory.mktemp("gzb") / "trunc.fq.gz"
    whole = gzip.compress(_fastq_bytes(), compresslevel=1)
    # Cut inside the deflate stream, so the member has no end-of-stream marker.
    p.write_bytes(whole[: len(whole) // 2])
    return p


@pytest.fixture(scope="module")
def corrupt_gz(tmp_path_factory):
    p = tmp_path_factory.mktemp("gzb") / "corrupt.fq.gz"
    whole = bytearray(gzip.compress(_fastq_bytes(), compresslevel=1))
    # Flip bits in the middle of the deflate stream but leave the trailer, so the
    # member decompresses to the wrong bytes and only the CRC catches it.
    whole[len(whole) // 2] ^= 0xFF
    p.write_bytes(bytes(whole))
    return p


#: The three ways to force one path, as environment overlays.
MODES = {
    "isal": {"ZNA_NO_ISAL": None, "ZNA_NO_EXTERNAL_GZIP": "1"},
    "external": {"ZNA_NO_ISAL": "1", "ZNA_NO_EXTERNAL_GZIP": None},
    "stdlib": {"ZNA_NO_ISAL": "1", "ZNA_NO_EXTERNAL_GZIP": "1"},
}


def _requires(mode):
    if mode == "isal" and not HAVE_ISAL:
        pytest.skip("isal not installed (optional extra)")
    if mode == "external" and not HAVE_EXTERNAL:
        pytest.skip("no pigz/gzip on PATH")


def _run_encode(path, out, env_overlay):
    """Run ``zna encode`` in a subprocess so the env overlay really applies.

    In process would not do: :func:`zna._gzip.isal_module` caches its probe for the
    life of the process, which is what makes it cheap, so ``ZNA_NO_ISAL`` has to be set
    before the interpreter that reads it starts.
    """
    import os
    env = dict(os.environ)
    for k, v in env_overlay.items():
        if v is None:
            env.pop(k, None)
        else:
            env[k] = v
    return subprocess.run(
        [sys.executable, "-m", "zna.cli", "encode", "--interleaved",
         str(path), "-o", str(out), "-q"],
        env=env, capture_output=True, text=True,
    )


class TestPolicy:
    """:func:`zna._gzip.prefer_isal` and the name it reports."""

    def test_own_read_thread_prefers_isal(self):
        assert _gzip.prefer_isal(own_read_thread=True) is True

    def test_main_thread_reader_prefers_subprocess(self, monkeypatch):
        monkeypatch.delenv("ZNA_NO_EXTERNAL_GZIP", raising=False)
        assert _gzip.prefer_isal(own_read_thread=False) is False

    def test_main_thread_reader_takes_isal_when_no_subprocess(self, monkeypatch):
        # With the subprocess ruled out, ISA-L is the best remaining option even for a
        # main-thread reader: its real competition there is the stdlib's 208 MB/s.
        monkeypatch.setenv("ZNA_NO_EXTERNAL_GZIP", "1")
        assert _gzip.prefer_isal(own_read_thread=False) is True

    def test_plain_path_is_never_a_gzip_backend(self):
        assert _gzip.inflate_backend_name("x.fq") == "plain"
        assert _gzip.inflate_backend_name("x.fastq") == "plain"

    def test_backend_name_is_one_of_the_known_paths(self):
        for own in (True, False):
            assert _gzip.inflate_backend_name("x.fq.gz", own_read_thread=own) in {
                "isal", "external", "gzip"}

    def test_isal_probe_is_cached(self):
        # Two calls must return the identical object, not merely an equal one: the
        # probe is meant to happen once per process.
        assert _gzip.isal_module() is _gzip.isal_module()


class TestOpenHelpers:
    """The two direct openers return the same bytes as the stdlib."""

    def test_stdlib_opener_roundtrips(self, gz_fastq):
        with _gzip.open_stdlib_gzip(str(gz_fastq), 1 << 16) as fh:
            assert fh.read() == _fastq_bytes()

    @pytest.mark.skipif(not HAVE_ISAL, reason="isal not installed (optional extra)")
    def test_isal_opener_roundtrips(self, gz_fastq):
        stream = _gzip.open_isal(str(gz_fastq), 1 << 16)
        assert stream is not None
        with stream as fh:
            assert fh.read() == _fastq_bytes()

    @pytest.mark.skipif(not HAVE_ISAL, reason="isal not installed (optional extra)")
    def test_isal_opener_supports_readline(self, gz_fastq):
        # The encode path drives readline, not read; a BufferedReader wrapper that got
        # the layering wrong would still pass a read() test.
        with _gzip.open_isal(str(gz_fastq), 1 << 16) as fh:
            assert fh.readline() == b"@rd0/1\n"
            assert len(fh.readline().rstrip()) == READ_LEN

    def test_isal_opener_returns_none_when_disabled(self, gz_fastq, monkeypatch):
        monkeypatch.setattr(_gzip, "_ISAL", False)
        assert _gzip.open_isal(str(gz_fastq), 1 << 16) is None


class TestAllPathsAgree:
    """Same input, three decompressors, one output."""

    @pytest.mark.parametrize("mode", list(MODES))
    def test_encode_from_gz_matches_encode_from_plain(
            self, mode, gz_fastq, plain_fastq, tmp_path):
        _requires(mode)
        ref = tmp_path / "ref.zna"
        got = tmp_path / f"{mode}.zna"
        r = _run_encode(plain_fastq, ref, {})
        assert r.returncode == 0, r.stderr
        r = _run_encode(gz_fastq, got, MODES[mode])
        assert r.returncode == 0, r.stderr
        assert got.read_bytes() == ref.read_bytes()


class TestDamagedInputIsStillRefused:
    """Throughput must not cost the integrity check."""

    @pytest.mark.parametrize("mode", list(MODES))
    def test_truncated_gz_fails(self, mode, truncated_gz, tmp_path):
        _requires(mode)
        out = tmp_path / f"trunc_{mode}.zna"
        r = _run_encode(truncated_gz, out, MODES[mode])
        assert r.returncode != 0, (
            f"{mode} accepted a truncated .gz and wrote {out} — a silently truncated "
            f"corpus is the failure this check exists for")

    @pytest.mark.parametrize("mode", list(MODES))
    def test_corrupt_gz_fails(self, mode, corrupt_gz, tmp_path):
        _requires(mode)
        out = tmp_path / f"corrupt_{mode}.zna"
        r = _run_encode(corrupt_gz, out, MODES[mode])
        assert r.returncode != 0, f"{mode} accepted a corrupt .gz and wrote {out}"


class TestMergeReadPath:
    """``zna merge``'s reader resolves through the same policy."""

    @pytest.mark.parametrize("mode", list(MODES))
    def test_merge_output_identical_across_backends(self, mode, tmp_path):
        _requires(mode)
        import os
        r1 = tmp_path / "r1.fq.gz"
        r2 = tmp_path / "r2.fq.gz"
        # Two mates that overlap over their whole length, so every pair merges and the
        # scan is actually exercised rather than short-circuited.
        frag = "".join("ACGT"[(i * 5 + 1) % 4] for i in range(READ_LEN))
        comp = str.maketrans("ACGT", "TGCA")
        rc = frag.translate(comp)[::-1]
        recs1 = "".join(f"@rd{i}/1\n{frag}\n+\n{'I' * READ_LEN}\n" for i in range(500))
        recs2 = "".join(f"@rd{i}/2\n{rc}\n+\n{'I' * READ_LEN}\n" for i in range(500))
        with gzip.open(r1, "wb", compresslevel=1) as fh:
            fh.write(recs1.encode())
        with gzip.open(r2, "wb", compresslevel=1) as fh:
            fh.write(recs2.encode())

        env = dict(os.environ)
        for k, v in MODES[mode].items():
            if v is None:
                env.pop(k, None)
            else:
                env[k] = v
        out = tmp_path / f"m_{mode}.fq"
        # This test is about the READER, not the kernel, so it runs on whichever
        # merge backend the environment has.  `zna merge` refuses to fall back to
        # the reference kernel silently (a 50x slowdown looks like a slow node at
        # cluster scale), so the extension-less configuration must ask for it by
        # name -- which is exactly the configuration that first caught this test
        # assuming the compiled backend was always present.
        from zna.merge.backend import available_merge_backends
        backend = "accel" if "accel" in available_merge_backends() else "python"
        r = subprocess.run(
            [sys.executable, "-m", "zna.merge", "--in1", str(r1), "--in2", str(r2),
             "--out", str(out), "--threads", "1", "-q", "--backend", backend],
            env=env, capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
        body = out.read_bytes()
        assert body.count(b"@") >= 500
        # Every pair fully overlaps, so every one merges to a single record.
        assert body.count(b"merged_") == 500
