#!/usr/bin/env python3
import sys, argparse, gzip, io, os
from multiprocessing import cpu_count

try:
    import pysam  # optional, not required
except Exception:
    pysam = None


def open_any(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--ids",
        required=True,
        help="One OrthoDB gene id per line (e.g., 3702_0:001ABC)",
    )
    ap.add_argument("--fasta", required=True, help="odb12v1_cds_fasta.gz")
    ap.add_argument("--out", required=True, help="Output FASTA .gz")
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()

    # Load ID set
    ids = set()
    with open(args.ids) as f:
        for line in f:
            s = line.strip()
            if s:
                ids.add(s)

    # Stream through FASTA
    # Headers in OrthoDB FASTA start with OrthoDB gene id (per docs), so split first token
    # and test membership.
    with open_any(args.fasta, "rt") as fi, gzip.open(args.out, "wt") as fo:
        keep = False
        for line in fi:
            if line.startswith(">"):
                hdr = line[1:].strip().split()[0]
                keep = hdr in ids
            if keep:
                fo.write(line)


if __name__ == "__main__":
    main()
