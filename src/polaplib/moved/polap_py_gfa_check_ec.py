#!/usr/bin/env python3
# scripts/polap_py_gfa_check_ec.py
# Version: v0.2.0
# Usage: polap_py_gfa_check_ec.py --gfa file.gfa
import sys, argparse, gzip


def openg(p, mode="rt"):
    return gzip.open(p, mode) if p.endswith(".gz") else open(p, mode)


ap = argparse.ArgumentParser()
ap.add_argument("--gfa", required=True)
a = ap.parse_args()

missing = 0
shown = 0
with openg(a.gfa, "rt") as f:
    for ln in f:
        if not ln or ln[0] in "#\n":
            continue
        cols = ln.rstrip("\n").split("\t")
        if cols[0] != "L":
            continue
        tags = cols[6:] if len(cols) > 6 else []
        if not any(t.startswith("ec:i:") for t in tags):
            missing += 1
            if shown < 10:
                print(ln.rstrip("\n"), file=sys.stderr)
                shown += 1

print(f"L lines missing ec:i = {missing}")
if missing > 0:

    sys.exit(3)
