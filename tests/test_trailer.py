"""The provenance prologue, the stats trailer, and the footer (format 0.5).

The layout is symmetric and the split is principled: position in the file
follows when the information is knowable.  Facts known at encode START --
writer version, shuffled order, in-process merging -- live in a count-0
pseudo-block right after the header, readable before any record, on pipes too.
Facts that need the COMPLETE encode -- counts, histograms, the block index --
live in the trailer at the end, discovered O(1) from the 32-byte footer.

The format version byte stays 3: both structures are additive, and a 0.4.1
reader walks a 0.5 file to a clean stop, decoding each pseudo-block as a valid
empty block (pinned below, in TestOldReaderDegradesBenignly).
"""
import io
import struct
import unittest
import zlib

import zna
from zna.core import (
    COMPRESSION_NONE,
    COMPRESSION_ZSTD,
    ZnaHeader,
    ZnaReader,
    ZnaWriter,
    _BLOCK_HEADER_FMT,
    _BLOCK_HEADER_SIZE,
    _FOOTER_FMT,
    _FOOTER_MAGIC,
    _FOOTER_SIZE,
)
from zna._shuffle import shuffle_zna


def build(comp=COMPRESSION_ZSTD, n_pairs=40, n_single=5, block_size=512,
          seq_len_bytes=2, **writer_kw) -> bytes:
    import random
    rng = random.Random(20260818)
    buf = io.BytesIO()
    h = ZnaHeader(read_group="trailer-test", seq_len_bytes=seq_len_bytes,
                  compression_method=comp)
    with ZnaWriter(buf, h, block_size=block_size, **writer_kw) as w:
        for _ in range(n_pairs):
            w.write_record("".join(rng.choices("ACGT", k=rng.randint(30, 90))),
                           True, True, False)
            w.write_record("".join(rng.choices("ACGT", k=rng.randint(30, 90))),
                           True, False, True)
        for _ in range(n_single):
            w.write_record("".join(rng.choices("ACGT", k=rng.randint(30, 190))),
                           False, False, False, is_full_fragment=True)
    return buf.getvalue()


def strip_to_legacy(data: bytes) -> bytes:
    """Rebuild the 0.4.1-era shape of a 0.5 file: header + data blocks only."""
    r = ZnaReader(io.BytesIO(data))
    t = r.trailer
    _m, _v, _r_, _crc, sentinel_offset = struct.unpack(
        _FOOTER_FMT, data[-_FOOTER_SIZE:])
    # The prologue header sits between the file header and data_start; find it
    # by its own arithmetic: header_len + 20 + payload_len == data_start.
    for plen in range(t.data_start):
        off = t.data_start - _BLOCK_HEADER_SIZE - plen
        if off < 0:
            break
        c, _u, n, f_, l_ = struct.unpack(
            _BLOCK_HEADER_FMT, data[off:off + _BLOCK_HEADER_SIZE])
        if c == plen and n == 0 and f_ == 0 and l_ == 0:
            return data[:off] + data[t.data_start:sentinel_offset]
    raise AssertionError("no prologue found")


class NoSeek:
    """A pipe: read-only, no seek, no tell."""

    def __init__(self, data: bytes) -> None:
        self._b = io.BytesIO(data)

    def read(self, n: int = -1) -> bytes:
        return self._b.read(n)


class TestProvenancePrologue(unittest.TestCase):
    """Start-known facts are readable before the first record."""

    def test_fresh_encode_stamps_version_and_not_shuffled(self):
        r = ZnaReader(io.BytesIO(build()))
        p = r.provenance
        self.assertEqual(p.writer_version, zna.__version__)
        self.assertFalse(p.shuffled)
        self.assertFalse(p.merged_in_process)

    def test_provenance_is_available_on_a_pipe(self):
        """The point of putting it at the front: no seek needed."""
        r = ZnaReader(NoSeek(build()))
        self.assertIsNotNone(r.provenance)
        self.assertEqual(r.provenance.writer_version, zna.__version__)
        self.assertIsNone(r.trailer)  # the trailer, by contrast, needs a seek

    def test_a_legacy_file_has_none(self):
        legacy = strip_to_legacy(build())
        r = ZnaReader(io.BytesIO(legacy))
        self.assertIsNone(r.provenance)
        self.assertIsNone(r.trailer)


class TestTrailer(unittest.TestCase):
    def test_counts_match_a_full_decode(self):
        data = build(n_pairs=40, n_single=5)
        r = ZnaReader(io.BytesIO(data))
        t = r.trailer
        recs = list(ZnaReader(io.BytesIO(data)).records())
        self.assertEqual(t.n_records, 85)
        self.assertEqual(t.n_pairs, 40)
        self.assertEqual(t.n_unpaired, 5)
        self.assertEqual(t.n_bases, sum(len(x[0]) for x in recs))
        self.assertEqual(sum(t.flag_counts.values()), 85)
        self.assertEqual(sum(t.length_histogram.values()), 85)
        self.assertEqual(sum(t.length_histogram_unpaired.values()), 5)

    def test_uncompressed_files_carry_raw_json(self):
        data = build(comp=COMPRESSION_NONE)
        t = ZnaReader(io.BytesIO(data)).trailer
        self.assertEqual(t.n_records, 85)

    def test_stored_index_equals_a_scan(self):
        """The trailer's block index and the legacy walk agree exactly."""
        data = build()
        stored = ZnaReader(io.BytesIO(data)).block_index()
        # Force the scan by stripping the trailer, then map offsets back: the
        # legacy copy has no prologue, so every offset shifts by its length.
        legacy = strip_to_legacy(data)
        scanned = ZnaReader(io.BytesIO(legacy)).block_index()
        shift = stored[0].offset - scanned[0].offset
        self.assertEqual(
            [(b.index, b.offset - shift, b.n_records, b.comp_size, b.uncomp_size)
             for b in stored],
            [tuple(b) for b in scanned],
        )

    def test_deterministic_bytes(self):
        """Identical input, seed, and environment => byte-identical files."""
        self.assertEqual(build(), build())

    def test_corrupt_payload_raises_not_none(self):
        """Damage must not be indistinguishable from age (a missing trailer)."""
        data = bytearray(build())
        _m, _v, _r_, _crc, sentinel = struct.unpack(
            _FOOTER_FMT, bytes(data[-_FOOTER_SIZE:]))
        data[sentinel + _BLOCK_HEADER_SIZE + 4] ^= 0xFF  # flip a payload byte
        with self.assertRaisesRegex(ValueError, "CRC mismatch"):
            ZnaReader(io.BytesIO(bytes(data))).trailer

    def test_missing_footer_reads_as_no_trailer(self):
        data = build()
        r = ZnaReader(io.BytesIO(data[:-_FOOTER_SIZE]))
        self.assertIsNone(r.trailer)
        # ...and the records are still all there (sentinel termination).
        self.assertEqual(len(list(r.records())), 85)

    def test_trailer_leaves_read_position_unchanged(self):
        data = build()
        r = ZnaReader(io.BytesIO(data))
        first_batch = next(r.blocks())
        _ = r.trailer
        rest = sum(len(s) for s, _f in r.blocks())
        self.assertEqual(len(first_batch[0]) + rest, 85)


class TestTermination(unittest.TestCase):
    """records()/blocks()/block_index() stop at the sentinel; pipes included."""

    def test_all_walks_on_a_trailer_bearing_file(self):
        data = build()
        self.assertEqual(len(list(ZnaReader(io.BytesIO(data)).records())), 85)
        self.assertEqual(
            sum(len(s) for s, _f in ZnaReader(io.BytesIO(data)).blocks()), 85)
        idx = ZnaReader(io.BytesIO(data)).block_index()
        self.assertEqual(sum(b.n_records for b in idx), 85)

    def test_records_and_blocks_on_a_pipe(self):
        data = build()
        self.assertEqual(len(list(ZnaReader(NoSeek(data)).records())), 85)
        self.assertEqual(
            sum(len(s) for s, _f in ZnaReader(NoSeek(data)).blocks()), 85)

    def test_legacy_file_on_a_pipe_replays_the_probed_header(self):
        """__init__'s prologue probe reads 20 bytes a pipe cannot give back."""
        legacy = strip_to_legacy(build())
        self.assertEqual(len(list(ZnaReader(NoSeek(legacy)).records())), 85)
        self.assertEqual(
            sum(len(s) for s, _f in ZnaReader(NoSeek(legacy)).blocks()), 85)

    def test_legacy_file_seekable(self):
        legacy = strip_to_legacy(build())
        self.assertEqual(len(list(ZnaReader(io.BytesIO(legacy)).records())), 85)
        idx = ZnaReader(io.BytesIO(legacy)).block_index()
        self.assertEqual(sum(b.n_records for b in idx), 85)

    def test_trailing_garbage_no_longer_yields_a_phantom_block(self):
        """block_index() on a garbage tail raised nothing before 0.5 -- it
        appended a phantom BlockInfo, seeked past EOF, and returned success."""
        legacy = strip_to_legacy(build())
        garbage = legacy + b"\x21" * 24  # a parseable-but-absurd block header
        with self.assertRaisesRegex(ValueError, "trailing garbage|truncated"):
            ZnaReader(io.BytesIO(garbage)).block_index()


class TestStamping(unittest.TestCase):
    def test_shuffle_stamps_shuffled_true(self):
        import tempfile, pathlib
        with tempfile.TemporaryDirectory() as d:
            src = pathlib.Path(d) / "in.zna"
            dst = pathlib.Path(d) / "out.zna"
            src.write_bytes(build(n_pairs=200, n_single=20, block_size=256))
            shuffle_zna(str(src), str(dst), seed=7, buffer_bytes=1 << 16,
                        block_size=256, tmp_dir=d, quiet=True)
            data = dst.read_bytes()
        r = ZnaReader(io.BytesIO(data))
        self.assertTrue(r.provenance.shuffled)
        self.assertEqual(r.provenance.writer_version, zna.__version__)
        # A shuffle is a permutation, so every trailer stat is order-invariant.
        t = r.trailer
        t0 = ZnaReader(io.BytesIO(build(n_pairs=200, n_single=20,
                                        block_size=256))).trailer
        self.assertEqual(t.n_records, t0.n_records)
        self.assertEqual(t.flag_counts, t0.flag_counts)
        self.assertEqual(t.length_histogram, t0.length_histogram)
        self.assertEqual(t.length_histogram_unpaired,
                         t0.length_histogram_unpaired)

    def test_writer_kwarg_round_trips(self):
        data = build(shuffled=True, merged_in_process=True)
        p = ZnaReader(io.BytesIO(data)).provenance
        self.assertTrue(p.shuffled)
        self.assertTrue(p.merged_in_process)


class TestWriterEdges(unittest.TestCase):
    def test_zero_record_file_has_a_trailer(self):
        buf = io.BytesIO()
        with ZnaWriter(buf, ZnaHeader(read_group="empty", seq_len_bytes=2,
                                      compression_method=COMPRESSION_ZSTD)):
            pass
        r = ZnaReader(io.BytesIO(buf.getvalue()))
        self.assertEqual(r.trailer.n_records, 0)
        self.assertEqual(list(r.records()), [])
        self.assertEqual(r.block_index(), [])

    def test_an_aborted_encode_leaves_no_trailer(self):
        """The trailer is the writer's signature that the encode FINISHED."""
        buf = io.BytesIO()
        with self.assertRaises(RuntimeError):
            with ZnaWriter(buf, ZnaHeader(read_group="x", seq_len_bytes=2,
                                          compression_method=COMPRESSION_ZSTD)) as w:
                w.write_record("ACGT", False, False, False)
                raise RuntimeError("simulated crash")
        r = ZnaReader(io.BytesIO(buf.getvalue()))
        self.assertIsNone(r.trailer)                     # not certified...
        self.assertEqual(len(list(r.records())), 1)      # ...but readable

    def test_dangling_r1_poisons_the_writer(self):
        """After the dangling-R1 error a second close() is a no-op -- it must
        not flush the orphan (a block ending mid-fragment) nor sign a trailer
        whose pair counts cannot balance."""
        buf = io.BytesIO()
        w = ZnaWriter(buf, ZnaHeader(read_group="x", seq_len_bytes=2,
                                     compression_method=COMPRESSION_ZSTD))
        w.write_record("ACGTACGT", True, True, False)
        w.write_record("TTTTGGGG", True, False, True)
        w.write_record("CCCCAAAA", True, True, False)      # orphan
        with self.assertRaisesRegex(ValueError, "R2 never arrived"):
            w.close()
        size = len(buf.getvalue())
        w.close()                                          # silent no-op
        self.assertEqual(len(buf.getvalue()), size)
        r = ZnaReader(io.BytesIO(buf.getvalue()))
        self.assertIsNone(r.trailer)                       # not certified
        # The whole buffered tail died with the writer, as close() documents.
        self.assertEqual(list(r.records()), [])

    def test_abort_drops_an_open_fragment_but_flushes_whole_ones(self):
        """The __exit__ abort path may not write a block ending on a lone R1:
        the orphan is dropped, everything complete before it survives."""
        from zna.core import OPENS_FRAGMENT
        buf = io.BytesIO()
        with self.assertRaises(RuntimeError):
            with ZnaWriter(buf, ZnaHeader(read_group="x", seq_len_bytes=2,
                                          compression_method=COMPRESSION_ZSTD)) as w:
                w.write_record("ACGTACGT", True, True, False)
                w.write_record("TTTTGGGG", True, False, True)
                w.write_record("CCCCAAAA", True, True, False)   # open fragment
                raise RuntimeError("crash")
        recs = list(ZnaReader(io.BytesIO(buf.getvalue())).records())
        self.assertEqual(len(recs), 2)
        for seqs, flags in ZnaReader(io.BytesIO(buf.getvalue())).blocks():
            self.assertFalse(OPENS_FRAGMENT[flags[-1]],
                             "a flushed block ends mid-fragment")

    def test_close_is_idempotent(self):
        buf = io.BytesIO()
        w = ZnaWriter(buf, ZnaHeader(read_group="x", seq_len_bytes=2,
                                     compression_method=COMPRESSION_ZSTD))
        w.write_record("ACGT", False, False, False)
        w.close()
        size = len(buf.getvalue())
        w.close()
        self.assertEqual(len(buf.getvalue()), size, "second close wrote bytes")

    def test_a_closed_writer_refuses_writes(self):
        """A record accepted after close() would land behind the sentinel,
        invisible to every reader -- silent data loss, so it raises instead."""
        w = ZnaWriter(io.BytesIO(), ZnaHeader(read_group="x", seq_len_bytes=2,
                                              compression_method=COMPRESSION_ZSTD))
        w.write_record("ACGT", False, False, False)
        w.close()
        with self.assertRaisesRegex(ValueError, "closed"):
            w.write_record("ACGT", False, False, False)
        with self.assertRaisesRegex(ValueError, "closed"):
            w.write_records([("ACGT", False, False, False)])

    def test_a_second_walk_on_an_exhausted_reader_is_empty(self):
        """Sentinel termination consumes the trailer bytes, so a re-walk sees
        a clean EOF -- identical to the 0.4.1 second-walk behavior -- rather
        than trailer bytes misread as a block header."""
        data = build()
        r = ZnaReader(io.BytesIO(data))
        self.assertEqual(len(list(r.records())), 85)
        self.assertEqual(list(r.records()), [])
        r2 = ZnaReader(io.BytesIO(data))
        self.assertEqual(sum(len(s) for s, _f in r2.blocks()), 85)
        self.assertEqual(list(r2.blocks()), [])

    def test_corrupt_prologue_raises_valueerror(self):
        """A flipped prologue byte must surface as the documented ValueError,
        whatever the underlying decoder raised."""
        data = build()
        t = ZnaReader(io.BytesIO(data)).trailer
        b = bytearray(data)
        for plen in range(t.data_start):
            off = t.data_start - _BLOCK_HEADER_SIZE - plen
            if off < 0:
                break
            c, _u, n, f_, l_ = struct.unpack(
                _BLOCK_HEADER_FMT, bytes(b[off:off + _BLOCK_HEADER_SIZE]))
            if c == plen and n == 0 and f_ == 0 and l_ == 0:
                b[off + _BLOCK_HEADER_SIZE + 3] ^= 0xFF
                break
        with self.assertRaisesRegex(ValueError, "prologue"):
            ZnaReader(io.BytesIO(bytes(b)))

    def test_every_seq_len_width(self):
        for slb in (1, 2, 4):
            data = build(seq_len_bytes=slb, n_pairs=10, n_single=2)
            t = ZnaReader(io.BytesIO(data)).trailer
            self.assertEqual(t.n_records, 22, f"seq_len_bytes={slb}")


class TestOldReaderDegradesBenignly(unittest.TestCase):
    """Pin D1/A1: a 0.4.1 reader on a 0.5 file reads every record, decodes each
    pseudo-block as a valid empty block, and stops at a clean EOF.

    These walkers reproduce the 0.4.1 termination conditions exactly (break on
    empty read only; raise on short header or short payload).  If the sentinel
    layout is ever touched -- say ``comp_size`` stops covering the footer --
    this is the test that goes red.
    """

    def _walk_records(self, data, comp):
        import zstandard
        legacy = strip_to_legacy(data)
        # 0.4.1 sat right after the file header; that offset is where the
        # legacy copy and the real file first differ (the prologue header).
        header_len = next(
            i for i in range(min(len(legacy), len(data)))
            if i >= len(legacy) or legacy[i] != data[i]
        )
        fh = io.BytesIO(data)
        fh.seek(header_len)
        dctx = zstandard.ZstdDecompressor() if comp == COMPRESSION_ZSTD else None
        total = empty = 0
        while True:
            hdr = fh.read(_BLOCK_HEADER_SIZE)
            if not hdr:
                return total, empty, "clean-eof"
            if len(hdr) < _BLOCK_HEADER_SIZE:
                return total, empty, "short-header"
            c, u, n, _f, _l = struct.unpack(_BLOCK_HEADER_FMT, hdr)
            payload = fh.read(c)
            if len(payload) != c:
                return total, empty, "short-payload"
            if dctx is not None:
                dctx.decompress(payload, max_output_size=u)
            total += n
            empty += (n == 0)

    def test_both_compression_modes(self):
        for comp in (COMPRESSION_ZSTD, COMPRESSION_NONE):
            data = build(comp=comp)
            total, empty, term = self._walk_records(data, comp)
            self.assertEqual(
                (total, empty, term), (85, 2, "clean-eof"),
                f"comp={comp}: a 0.4.1-style walk must read every record, see "
                f"two empty pseudo-blocks, and stop cleanly; got "
                f"{total} records, {empty} empties, {term}",
            )


class TestFooter(unittest.TestCase):
    def test_footer_shape(self):
        data = build()
        magic, ver, reserved, crc, sentinel = struct.unpack(
            _FOOTER_FMT, data[-_FOOTER_SIZE:])
        self.assertEqual(magic, _FOOTER_MAGIC)
        self.assertEqual(ver, 1)
        self.assertEqual(reserved, 0)
        # sentinel_offset points at a count-0 block header whose comp_size
        # covers payload + footer, i.e. everything to EOF.
        c, u, n, f_, l_ = struct.unpack(
            _BLOCK_HEADER_FMT, data[sentinel:sentinel + _BLOCK_HEADER_SIZE])
        self.assertEqual(n, 0)
        self.assertEqual(sentinel + _BLOCK_HEADER_SIZE + c, len(data))
        # crc covers the payload as written
        payload = data[sentinel + _BLOCK_HEADER_SIZE:-_FOOTER_SIZE]
        self.assertEqual(zlib.crc32(payload) & 0xFFFFFFFF, crc)


if __name__ == "__main__":
    unittest.main()
