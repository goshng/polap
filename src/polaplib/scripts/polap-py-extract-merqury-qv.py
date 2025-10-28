#!/usr/bin/env python3
# Version: v0.1.0
# Name: polap-py-extract-merqury-qv.py
# Purpose: Extract Merqury QV from a merqury.sh run folder into a key-value TSV.
# Usage:
#   python3 polap-py-extract-merqury-qv.py --merq-out OUT_PREFIX --out A_qv.tsv
# Output TSV:
#   metric\tvalue
#   merqury_qv\t<value>
import sys, argparse, re, os

ap = argparse.ArgumentParser()
ap.add_argument(
    "--merq-out", required=True, help="Merqury output prefix directory (contains 'qv')"
)
ap.add_argument("--out", required=True)
args = ap.parse_args()

qv_path = os.path.join(args.merq_out, "qv")
qv_val = None


def try_parse(line):
    # Accept lines like: "QV: 42.37" or "Assembly QV: 42.37" or just a number
    m = re.search(r"([0-9]+(?:\.[0-9]+)?)", line)
    return float(m.group(1)) if m else None


if os.path.exists(qv_path):
    with open(qv_path) as f:
        for ln in f:
            v = try_parse(ln)
            if v is not None:
                qv_val = v
else:
    # Some Merqury versions put it in *.qv file(s); scan directory
    for root, _, files in os.walk(args.merq_out):
        for fn in files:
            if fn.endswith(".qv") or fn == "qv":
                with open(os.path.join(root, fn)) as f:
                    for ln in f:
                        v = try_parse(ln)
                        if v is not None:
                            qv_val = v

with open(args.out, "w") as w:
    w.write("metric\tvalue\n")
    if qv_val is None:
        w.write("merqury_qv\tNA\n")
    else:
        w.write(f"merqury_qv\t{qv_val:.2f}\n")
