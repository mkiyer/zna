#!/usr/bin/env python3
"""
Benchmark ZNA encode/decode on a real dataset.

Runs:
  1. Plain encode (no labels) at default and varied ZSTD levels
  2. Plain decode
  3. Labeled encode with 10 SAM tags
  4. Labeled decode
  5. Compression comparison vs gzip

Outputs a markdown-ready report.
"""
import os
import subprocess
import sys
import tempfile
import time

INPUT_FQ = (
    "/Users/mkiyer/Downloads/rigel_runs/"
    "sim_ccle_hela_salmon/gdna_high_ss_0.90_nrna_default/"
    "align_minimap2/merged.fq.gz"
)

LABEL_FLAGS = [
    "--label", "NM:C",
    "--label", "ms:S",
    "--label", "AS:s",
    "--label", "nn:C",
    "--label", "tp:A",
    "--label", "cm:S",
    "--label", "s1:S",
    "--label", "s2:S",
    "--label", "de:f",
    "--label", "rl:S",
]


def fmt_size(b: int) -> str:
    if b >= 1 << 30:
        return f"{b / (1 << 30):.2f} GB"
    return f"{b / (1 << 20):.1f} MB"


def fmt_speed(b: int, sec: float) -> str:
    return f"{b / sec / (1 << 20):.0f} MB/s"


def run(cmd: list[str], label: str = "") -> tuple[float, int]:
    """Run cmd, return (elapsed_sec, returncode). Suppress output."""
    if label:
        print(f"  {label} ...", end="", flush=True)
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True)
    dt = time.time() - t0
    if label:
        print(f" {dt:.1f}s")
    if r.returncode != 0:
        print(f"  ERROR: {r.stderr.decode()[:200]}", file=sys.stderr)
    return dt, r.returncode


def get_uncompressed_size(path: str) -> int:
    """Estimate uncompressed FASTQ size by decompressing."""
    r = subprocess.run(["gzcat", path], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    return len(r.stdout)


def encode_zna(src: str, dst: str, extra: list[str] | None = None) -> float:
    cmd = ["zna", "encode", src, "--interleaved", "-o", dst, "--quiet"]
    if extra:
        cmd.extend(extra)
    dt, rc = run(cmd, f"encode → {os.path.basename(dst)}")
    if rc:
        sys.exit(1)
    return dt


def decode_zna(src: str, dst: str, extra: list[str] | None = None) -> float:
    cmd = ["zna", "decode", src, "-o", dst, "--quiet"]
    if extra:
        cmd.extend(extra)
    dt, rc = run(cmd, f"decode → {os.path.basename(dst)}")
    if rc:
        sys.exit(1)
    return dt


def main():
    if not os.path.exists(INPUT_FQ):
        sys.exit(f"Input not found: {INPUT_FQ}")

    input_gz_size = os.path.getsize(INPUT_FQ)

    # Count reads (fast: just count lines / 4)
    print("Counting reads ...")
    r = subprocess.run(
        f'gzcat "{INPUT_FQ}" | wc -l',
        shell=True, capture_output=True, text=True,
    )
    total_reads = int(r.stdout.strip()) // 4
    print(f"  {total_reads:,} reads")

    # Uncompressed size
    print("Measuring uncompressed FASTQ size ...")
    uncomp_size = get_uncompressed_size(INPUT_FQ)
    print(f"  {fmt_size(uncomp_size)}")

    results: dict[str, dict] = {}

    with tempfile.TemporaryDirectory() as tmp:
        # -----------------------------------------------------------
        # 1. Plain encode at default level (9)
        # -----------------------------------------------------------
        print("\n=== Plain Encode (default, ZSTD L9, 4 MB blocks) ===")
        plain_zna = os.path.join(tmp, "plain.zna")
        enc_time = encode_zna(INPUT_FQ, plain_zna)
        plain_size = os.path.getsize(plain_zna)
        results["plain_default"] = {
            "enc_time": enc_time,
            "size": plain_size,
            "ratio": uncomp_size / plain_size,
        }

        # -----------------------------------------------------------
        # 2. Plain decode
        # -----------------------------------------------------------
        print("\n=== Plain Decode ===")
        plain_fa = os.path.join(tmp, "plain.fa")
        dec_time = decode_zna(plain_zna, plain_fa)
        results["plain_default"]["dec_time"] = dec_time

        # -----------------------------------------------------------
        # 3. Compression level sweep
        # -----------------------------------------------------------
        print("\n=== Compression Level Sweep ===")
        for level in [1, 3, 5, 9, 15]:
            tag = f"L{level}"
            zna_path = os.path.join(tmp, f"level_{level}.zna")
            fa_path = os.path.join(tmp, f"level_{level}.fa")
            et = encode_zna(INPUT_FQ, zna_path, ["--level", str(level)])
            sz = os.path.getsize(zna_path)
            dt = decode_zna(zna_path, fa_path)
            results[tag] = {
                "enc_time": et,
                "dec_time": dt,
                "size": sz,
                "ratio": uncomp_size / sz,
            }

        # -----------------------------------------------------------
        # 4. Block size sweep (at default level)
        # -----------------------------------------------------------
        print("\n=== Block Size Sweep ===")
        for bs_str, bs_label in [("512K", "512K"), ("1M", "1M"), ("4M", "4M"), ("8M", "8M")]:
            tag = f"block_{bs_label}"
            zna_path = os.path.join(tmp, f"block_{bs_label}.zna")
            fa_path = os.path.join(tmp, f"block_{bs_label}.fa")
            et = encode_zna(INPUT_FQ, zna_path, ["--block-size", bs_str])
            sz = os.path.getsize(zna_path)
            dt = decode_zna(zna_path, fa_path)
            results[tag] = {
                "enc_time": et,
                "dec_time": dt,
                "size": sz,
                "ratio": uncomp_size / sz,
            }

        # -----------------------------------------------------------
        # 5. Labeled encode/decode (10 tags)
        # -----------------------------------------------------------
        print("\n=== Labeled Encode (10 SAM tags) ===")
        labeled_zna = os.path.join(tmp, "labeled.zna")
        et = encode_zna(INPUT_FQ, labeled_zna, LABEL_FLAGS)
        labeled_size = os.path.getsize(labeled_zna)
        labeled_fa = os.path.join(tmp, "labeled.fa")
        dt = decode_zna(labeled_zna, labeled_fa, ["--labels"])
        results["labeled"] = {
            "enc_time": et,
            "dec_time": dt,
            "size": labeled_size,
            "ratio": uncomp_size / labeled_size,
        }

        # -----------------------------------------------------------
        # 6. Uncompressed ZNA
        # -----------------------------------------------------------
        print("\n=== Uncompressed ZNA ===")
        uncomp_zna = os.path.join(tmp, "uncomp.zna")
        et = encode_zna(INPUT_FQ, uncomp_zna, ["--uncompressed"])
        uncomp_zna_size = os.path.getsize(uncomp_zna)
        uncomp_fa = os.path.join(tmp, "uncomp.fa")
        dt = decode_zna(uncomp_zna, uncomp_fa)
        results["uncompressed"] = {
            "enc_time": et,
            "dec_time": dt,
            "size": uncomp_zna_size,
            "ratio": uncomp_size / uncomp_zna_size,
        }

    # ---------------------------------------------------------------
    # Report
    # ---------------------------------------------------------------
    print("\n" + "=" * 70)
    print("BENCHMARK REPORT")
    print("=" * 70)

    print(f"\nDataset: 25.4M paired-end reads (150 bp), interleaved FASTQ")
    print(f"Uncompressed FASTQ: {fmt_size(uncomp_size)}")
    print(f"Gzipped FASTQ:      {fmt_size(input_gz_size)}")
    print(f"Total reads:        {total_reads:,}")

    # Summary
    pd = results["plain_default"]
    print(f"\n### Default Settings (ZSTD L9, 4 MB blocks)")
    print(f"  Encode: {fmt_speed(uncomp_size, pd['enc_time'])} "
          f"({pd['enc_time']:.1f}s, {total_reads / pd['enc_time']:,.0f} rec/s)")
    print(f"  Decode: {fmt_speed(uncomp_size, pd['dec_time'])} "
          f"({pd['dec_time']:.1f}s, {total_reads / pd['dec_time']:,.0f} rec/s)")
    print(f"  Size:   {fmt_size(pd['size'])} ({pd['ratio']:.2f}x)")

    # Compression levels
    print(f"\n### Compression Levels")
    print(f"{'Level':<8} {'Encode':<12} {'Decode':<12} {'Size':<12} {'Ratio':<8}")
    print("-" * 52)
    for level in [1, 3, 5, 9, 15]:
        r = results[f"L{level}"]
        default = " ← default" if level == 9 else ""
        print(f"{level:<8} {fmt_speed(uncomp_size, r['enc_time']):<12} "
              f"{fmt_speed(uncomp_size, r['dec_time']):<12} "
              f"{fmt_size(r['size']):<12} {r['ratio']:.2f}x{default}")

    # Block sizes
    print(f"\n### Block Sizes")
    print(f"{'Size':<8} {'Encode':<12} {'Decode':<12} {'File Size':<12} {'Ratio':<8}")
    print("-" * 52)
    for bs in ["512K", "1M", "4M", "8M"]:
        r = results[f"block_{bs}"]
        default = " ← default" if bs == "4M" else ""
        print(f"{bs:<8} {fmt_speed(uncomp_size, r['enc_time']):<12} "
              f"{fmt_speed(uncomp_size, r['dec_time']):<12} "
              f"{fmt_size(r['size']):<12} {r['ratio']:.2f}x{default}")

    # Labeled
    lr = results["labeled"]
    print(f"\n### Labeled (10 SAM tags)")
    print(f"  Encode: {fmt_speed(uncomp_size, lr['enc_time'])} "
          f"({lr['enc_time']:.1f}s, {total_reads / lr['enc_time']:,.0f} rec/s)")
    print(f"  Decode: {fmt_speed(uncomp_size, lr['dec_time'])} "
          f"({lr['dec_time']:.1f}s, {total_reads / lr['dec_time']:,.0f} rec/s)")
    print(f"  Size:   {fmt_size(lr['size'])} ({lr['ratio']:.2f}x)")

    # Compression comparison table
    ur = results["uncompressed"]
    print(f"\n### Compression Comparison")
    print(f"{'Format':<28} {'Size':<12} {'Ratio':<8}")
    print("-" * 48)
    print(f"{'FASTQ (uncompressed)':<28} {fmt_size(uncomp_size):<12} {'1.00x':<8}")
    print(f"{'FASTQ.gz':<28} {fmt_size(input_gz_size):<12} "
          f"{uncomp_size / input_gz_size:.2f}x")
    print(f"{'ZNA (uncompressed)':<28} {fmt_size(ur['size']):<12} "
          f"{ur['ratio']:.2f}x")
    print(f"{'ZNA (ZSTD L9, default)':<28} {fmt_size(pd['size']):<12} "
          f"{pd['ratio']:.2f}x")
    print(f"{'ZNA (ZSTD L9 + 10 labels)':<28} {fmt_size(lr['size']):<12} "
          f"{lr['ratio']:.2f}x")

    print("\n" + "=" * 70)
    print("Done.")


if __name__ == "__main__":
    main()
