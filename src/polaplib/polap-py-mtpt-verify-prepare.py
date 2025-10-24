#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Populate verify/ for mtpt verification.

Inputs:
  --fasta MT_FASTA
  --reads LONG_READS (ONT/HiFi)
  --mtpt-tsv mtpt.tsv (tab-delimited: mt_contig, mt_start, mt_end, length, ...)
  --flank N (default 800)
  --threads T (default 8)
  --out verify  (output root for per-edge folders)
  --preset map-ont|map-hifi|sr (default map-ont)

For each mtpt row:
  verify/edge_<seqid>_<start>_<end>/
    tmp/
      junctions.edge_<seqid>.<start>.<end>.fa          # >JLEFT, >JRIGHT
      junctions.whole.edge_<seqid>.<start>.<end>.fa    # >JWHOLE
    reads/
      junc.edge_<seqid>.<start>.<end>.bam{,.bai}
      jwhole.edge_<seqid>.<start>.<end>.bam{,.bai}
"""

import argparse, csv, os, sys, subprocess
from typing import Tuple
import pysam


def eprint(*a):
    print(*a, file=sys.stderr)


def must(cmd: list, cwd: str = None):
    p = subprocess.run(cmd, cwd=cwd)
    if p.returncode != 0:
        raise RuntimeError(f"Command failed: {' '.join(cmd)}")


def fetch_circular(fa: pysam.FastaFile, ctg: str, s1: int, e1: int) -> str:
    """
    1-based inclusive coords on a circular contig.
    If s1 <= 0 or e1 > len, wrap around.
    """
    L = fa.get_reference_length(ctg)

    # Normalize range to [1..L], allowing wrap
    def clamp1(x):  # move into [1..L] with wrap
        x = ((x - 1) % L) + 1
        return x

    s = clamp1(s1)
    e = clamp1(e1)
    if s <= e:
        # pysam.fetch is 0-based, end-excl
        return fa.fetch(ctg, s - 1, e)
    else:
        # wrap-around
        part1 = fa.fetch(ctg, s - 1, L)
        part2 = fa.fetch(ctg, 0, e)
        return part1 + part2


def write_fasta(path: str, records: Tuple[Tuple[str, str], ...]):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        for name, seq in records:
            f.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i : i + 80] + "\n")


def build_edge(
    fa: pysam.FastaFile,
    reads: str,
    out_root: str,
    seqid: str,
    start: int,
    end: int,
    flank: int,
    threads: int,
    preset: str,
):
    # folders
    edge_key = f"edge_{seqid}_{start}_{end}"
    edge_dir = os.path.join(out_root, edge_key)
    tmp_dir = os.path.join(edge_dir, "tmp")
    rd_dir = os.path.join(edge_dir, "reads")
    os.makedirs(tmp_dir, exist_ok=True)
    os.makedirs(rd_dir, exist_ok=True)

    # windows (1-based inclusive)
    # JLEFT: center on 'start', width = 2*flank
    jl_s = start - flank
    jl_e = start + flank - 1
    # JRIGHT: center on 'end'
    jr_s = end - flank + 1
    jr_e = end + flank
    # JWHOLE: [start - flank, end + flank]
    jw_s = start - flank
    jw_e = end + flank

    # sequences
    jleft = fetch_circular(fa, seqid, jl_s, jl_e)
    jright = fetch_circular(fa, seqid, jr_s, jr_e)
    jwhole = fetch_circular(fa, seqid, jw_s, jw_e)

    # FASTAs (names expected by your verifier)
    junc_fa = os.path.join(tmp_dir, f"junctions.edge_{seqid}.{start}.{end}.fa")
    jwhole_fa = os.path.join(tmp_dir, f"junctions.whole.edge_{seqid}.{start}.{end}.fa")

    write_fasta(junc_fa, (("JLEFT", jleft), ("JRIGHT", jright)))
    write_fasta(jwhole_fa, (("JWHOLE", jwhole),))

    # Map reads -> BAMs (names expected by your verifier)
    junc_bam = os.path.join(rd_dir, f"junc.edge_{seqid}.{start}.{end}.bam")
    jwhole_bam = os.path.join(rd_dir, f"jwhole.edge_{seqid}.{start}.{end}.bam")

    # Choose minimap2 preset
    mm2_preset = {"map-ont": "map-ont", "map-hifi": "map-hifi", "sr": "sr"}.get(
        preset, "map-ont"
    )

    # junc (JLEFT+JRIGHT)
    must(
        [
            "bash",
            "-lc",
            f"minimap2 -x {mm2_preset} -a -t {threads} {sh(junc_fa)} {sh(reads)}"
            f" | samtools sort -@ {threads} -o {sh(junc_bam)} - && samtools index -@ {threads} {sh(junc_bam)}",
        ]
    )
    # jwhole (JWHOLE)
    must(
        [
            "bash",
            "-lc",
            f"minimap2 -x {mm2_preset} -a -t {threads} {sh(jwhole_fa)} {sh(reads)}"
            f" | samtools sort -@ {threads} -o {sh(jwhole_bam)} - && samtools index -@ {threads} {sh(jwhole_bam)}",
        ]
    )


def sh(path: str) -> str:
    # simple shell quoting
    return "'" + path.replace("'", "'\\''") + "'"


def main():
    ap = argparse.ArgumentParser(
        description="Prepare verify/ (windows + BAMs) for mtpt verification"
    )
    ap.add_argument("--fasta", required=True, help="Mitochondrial FASTA (circular)")
    ap.add_argument(
        "--reads",
        required=True,
        help="Long-read FASTQ/FASTA (gz ok if minimap2 supports)",
    )
    ap.add_argument(
        "--mtpt-tsv",
        required=True,
        help="TSV with columns: mt_contig, mt_start, mt_end, length, ...",
    )
    ap.add_argument(
        "--flank",
        type=int,
        default=800,
        help="Flank size (bp) for windows (default: 800)",
    )
    ap.add_argument(
        "--threads", type=int, default=8, help="Threads for mapping/sorting"
    )
    ap.add_argument(
        "--preset",
        default="map-ont",
        choices=["map-ont", "map-hifi", "sr"],
        help="minimap2 -x preset (default: map-ont)",
    )
    ap.add_argument("--out", default="verify", help="Output verify/ root directory")
    args = ap.parse_args()

    # Open FASTA
    try:
        fa = pysam.FastaFile(args.fasta)
    except Exception as e:
        eprint(f"[error] cannot open FASTA: {args.fasta} :: {e}")
        sys.exit(2)

    # Read mtpt.tsv
    try:
        with open(args.mtpt_tsv, newline="") as f:
            rdr = csv.DictReader(f, delimiter="\t")
            rows = list(rdr)
    except Exception as e:
        eprint(f"[error] cannot read mtpt.tsv: {args.mtpt_tsv} :: {e}")
        sys.exit(2)

    if not rows:
        eprint("[error] mtpt.tsv has no rows")
        sys.exit(2)

    os.makedirs(args.out, exist_ok=True)

    # Process each row
    for i, row in enumerate(rows, 1):
        try:
            seqid = row["mt_contig"]
            start = int(row["mt_start"])
            end = int(row["mt_end"])
        except Exception as e:
            eprint(f"[warn] skip row {i}: bad fields :: {e}")
            continue

        if seqid not in fa.references:
            eprint(f"[warn] skip row {i}: contig {seqid} not in FASTA")
            continue

        eprint(f"[info] building {seqid}:{start}-{end}")
        build_edge(
            fa,
            args.reads,
            args.out,
            seqid,
            start,
            end,
            args.flank,
            args.threads,
            args.preset,
        )

    fa.close()
    eprint("[OK] verify/ populated")


if __name__ == "__main__":
    main()
