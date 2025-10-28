#!/usr/bin/env python3
# Version: v0.1.0
# Name: polap-py-stack-compare-tables.py
# Purpose: Stack multiple compare tables into one with a section label column.
# Usage:
#   python3 polap-py-stack-compare-tables.py LABEL1 file1.tsv LABEL2 file2.tsv ... --out merged.tsv
# Output columns:
#   section  metric  A  B  delta(B-A)  better
import sys, argparse

ap = argparse.ArgumentParser()
ap.add_argument("pairs", nargs="+", help="LABEL and file pairs")
ap.add_argument("--out", required=True)
args = ap.parse_args()

if len(args.pairs) % 2 != 0:
    sys.stderr.write("Error: provide LABEL and file pairs\n")
    sys.exit(2)

rows = []
for i in range(0, len(args.pairs), 2):
    label, path = args.pairs[i], args.pairs[i + 1]
    try:
        with open(path) as f:
            hdr = f.readline()
            for ln in f:
                if not ln.strip():
                    continue
                parts = ln.rstrip("\n").split("\t")
                if len(parts) < 5:
                    continue
                rows.append([label] + parts)
    except FileNotFoundError:
        continue

with open(args.out, "w") as w:
    w.write("section\tmetric\tA\tB\tdelta(B-A)\tbetter\n")
    for r in rows:
        w.write("\t".join(r) + "\n")
