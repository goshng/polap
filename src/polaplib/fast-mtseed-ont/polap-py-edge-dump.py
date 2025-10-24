#!/usr/bin/env python3
# polap-py-edge-dump.py  v0.1.0
# Dump an edge list from a PAF:
#   u  v  alen  ident  weight
# Where:
#   alen  = PAF col11
#   ident = nmatch / alen  (PAF col10/col11)
#   weight= ident * (alen / min(qlen,tlen))
#
# Usage:
#   python polap-py-edge-dump.py --paf in.paf[.gz] --out edges.tsv[.gz] \
#     [--min_olen 0] [--min_ident 0] [--w_floor 0] [--dedup]
#
import sys, os, argparse, gzip
from collections import defaultdict


def openg(path, mode="rt"):
    if path == "-":
        return sys.stdin
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def openw(path):
    return (
        gzip.open(path, "wt")
        if path and path.endswith(".gz")
        else (open(path, "w") if path and path != "-" else sys.stdout)
    )


def parse_args():
    ap = argparse.ArgumentParser(description="PAF -> edge list (u v alen ident weight)")
    ap.add_argument("--paf", required=True, help="PAF file (.gz ok) or '-' for stdin")
    ap.add_argument("--out", default="-", help="output TSV (.gz ok) or '-' [stdout]")
    ap.add_argument("--min_olen", type=int, default=0, help="min aligned length [0]")
    ap.add_argument(
        "--min_ident", type=float, default=0.0, help="min identity in [0,1] [0]"
    )
    ap.add_argument(
        "--w_floor", type=float, default=0.0, help="min weight (ident*alen/shorter) [0]"
    )
    # Dedup is ON by default. User can pass --no-dedup to disable.
    ap.add_argument(
        "--dedup",
        dest="dedup",
        action="store_true",
        default=False,
        help="keep one edge per unordered pair (max weight) [default: on]",
    )
    ap.add_argument(
        "--no-dedup",
        dest="dedup",
        action="store_false",
        help="allow multiple edges per unordered pair (no deduplication)",
    )
    return ap.parse_args()


def dump_edges(args):
    # If dedup: pair -> alen,ident,weight)
    keep = {} if args.dedup else None
    out = openw(args.out)
    wrote_header = False

    def maybe_write(u, v, al, idn, w):
        nonlocal wrote_header
        if al < args.min_olen:
            return
        if idn < args.min_ident:
            return
        if w < args.w_floor:
            return
        if not wrote_header and out is not sys.stdout:
            pass  # header written at the end for sys.stdout parity
        if args.dedup:
            a, b = (u, v) if u <= v else (v, u)
            cur = keep.get((a, b))
            if cur is None or w > cur[2]:
                keep[(a, b)] = (al, idn, w)
        else:
            out.write(f"{u}\t{v}\t{al}\t{idn:.6f}\t{w:.6f}\n")

    with openg(args.paf, "rt") as f:
        for ln in f:
            if not ln or ln[0] == "#":
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 12:
                continue
            u = p[0]
            v = p[5]
            if u == v:
                continue
            try:
                ql = float(p[1])
                tl = float(p[6])
                qs = float(p[2])
                qe = float(p[3])
                ts = float(p[7])
                te = float(p[8])
                nm = float(p[9])
                al = float(p[10])
            except ValueError:
                continue
            if al <= 0 or ql <= 0 or tl <= 0:
                continue
            idn = max(0.0, min(1.0, nm / al))
            shorter = ql if ql < tl else tl
            frac = max(0.0, min(1.0, al / shorter))
            w = idn * frac
            maybe_write(u, v, int(al), idn, w)

    if args.dedup:
        for (a, b), (al, idn, w) in keep.items():
            out.write(f"{a}\t{b}\t{al}\t{idn:.6f}\t{w:.6f}\n")

    if out is not sys.stdout:
        out.close()


def main():
    args = parse_args()
    dump_edges(args)


if __name__ == "__main__":
    main()
