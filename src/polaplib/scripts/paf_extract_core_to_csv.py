#!/usr/bin/env python3
# scripts/paf_extract_core_to_csv.py
# Version: v0.1.0
# Read PAF (optionally .gz), emit CSV with:
# qname,qlen,qstart,qend,strand,tname,tlen,tstart,tend,nmatch,alen,mapq,tp,cm,s1,s2,identity,q_aln_frac
# Filters: --min-alen, --min-ident (identity=nmatch/alen)
import sys, os, gzip, argparse, csv


def openg(path):
    if path == "-":
        return sys.stdin
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def parse_tags(fields):
    tags = {"tp": "", "cm": "", "s1": "", "s2": ""}
    for f in fields[12:]:
        try:
            k, t, v = f.split(":", 2)
        except ValueError:
            continue
        if k in tags:
            tags[k] = v
    return tags


def main():
    ap = argparse.ArgumentParser(description="Extract core PAF fields to CSV")
    ap.add_argument("--paf", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-alen", type=int, default=50)
    ap.add_argument("--min-ident", type=float, default=0.0)
    args = ap.parse_args()

    n_in = n_out = 0
    with openg(args.paf) as f, open(args.out, "w", newline="") as w:
        cw = csv.writer(w)
        cw.writerow(
            [
                "qname",
                "qlen",
                "qstart",
                "qend",
                "strand",
                "tname",
                "tlen",
                "tstart",
                "tend",
                "nmatch",
                "alen",
                "mapq",
                "tp",
                "cm",
                "s1",
                "s2",
                "identity",
                "q_aln_frac",
            ]
        )
        for ln in f:
            if not ln.strip() or ln.startswith("#"):
                continue
            n_in += 1
            F = ln.rstrip("\n").split("\t")
            if len(F) < 12:
                continue
            try:
                qname = F[0]
                qlen = int(F[1])
                qstart = int(F[2])
                qend = int(F[3])
                strand = F[4]
                tname = F[5]
                tlen = int(F[6])
                tstart = int(F[7])
                tend = int(F[8])
                nmatch = int(F[9])
                alen = int(F[10])
                mapq = int(F[11])
            except Exception:
                continue
            if alen <= 0:
                continue
            ident = nmatch / float(alen)
            if alen < args.min_alen or ident < args.min_ident:
                continue
            qfrac = (qend - qstart) / float(qlen) if qlen > 0 else 0.0
            tags = parse_tags(F)
            cw.writerow(
                [
                    qname,
                    qlen,
                    qstart,
                    qend,
                    strand,
                    tname,
                    tlen,
                    tstart,
                    tend,
                    nmatch,
                    alen,
                    mapq,
                    tags["tp"],
                    tags["cm"],
                    tags["s1"],
                    tags["s2"],
                    f"{ident:.6f}",
                    f"{qfrac:.6f}",
                ]
            )
            n_out += 1
    print(f"[paf_extract_core_to_csv] kept {n_out}/{n_in} alignments", file=sys.stderr)


if __name__ == "__main__":
    main()
