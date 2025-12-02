#!/usr/bin/env python3
# File: scripts/polap_py_make_combined_fasta.py
# Version: v0.1.1
import sys, argparse, gzip


def open_any(path, mode="rt"):
    if path == "-" or path is None:
        return sys.stdin if "r" in mode else sys.stdout
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def stream_prefixed(fasta_path, prefix, out_fh):
    if not fasta_path:
        return
    op = open_any(fasta_path, "rt")
    for ln in op:
        if ln.startswith(">"):
            hdr = ln[1:].rstrip("\n")
            out_fh.write(">" + prefix + hdr + "\n")
        else:
            out_fh.write(ln)
    if op is not sys.stdin:
        op.close()


def main():
    ap = argparse.ArgumentParser(description="Prefix FASTA headers and concatenate.")
    ap.add_argument("--target-fasta", required=True)
    ap.add_argument("--target-prefix", required=True)
    ap.add_argument("--other-fasta", default=None)
    ap.add_argument("--other-prefix", default="")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    with open_any(a.out, "wt") as out:
        stream_prefixed(a.target_fasta, a.target_prefix, out)
        if a.other_fasta:
            stream_prefixed(a.other_fasta, a.other_prefix, out)


if __name__ == "__main__":
    main()
