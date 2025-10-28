#!/usr/bin/env python3
# Version: v0.1.0
# Name: polap-py-merge-kv-tables.py
# Purpose: Merge multiple "metric\tvalue" TSVs into one (later files override earlier on key clash).
#
# Usage:
#   python polap-py-merge-kv-tables.py --tables a.tsv b.tsv c.tsv --out merged.tsv
import sys, argparse
from collections import OrderedDict


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tables", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    return ap.parse_args()


def load_kv(path):
    D = OrderedDict()
    with open(path) as f:
        header = f.readline()
        for ln in f:
            if not ln.strip():
                continue
            k, v = ln.rstrip("\n").split("\t", 1)
            D[k] = v
    return D


def main():
    a = parse_args()
    merged = OrderedDict()
    for t in a.tables:
        kv = load_kv(t)
        merged.update(kv)
    with open(a.out, "w") as w:
        w.write("metric\tvalue\n")
        for k, v in merged.items():
            w.write(f"{k}\t{v}\n")


if __name__ == "__main__":
    main()
