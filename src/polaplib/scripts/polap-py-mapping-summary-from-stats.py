#!/usr/bin/env python3
# Version: v0.2.0
# Name: polap-py-mapping-summary-from-stats.py
# Purpose: Parse `samtools stats` to produce alignment-derived identity & QV summary.
#
# Input:  --stats mapping.stats  (from: samtools stats reads_to_polished.bam)
# Output: --out   TSV with: error_rate, identity, QV, insertion_rate, deletion_rate, mismatch_rate
#
# QV = -10 * log10(error_rate), Identity = 1 - error_rate
import sys, argparse, math


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stats", required=True)
    ap.add_argument("--out", required=True)
    return ap.parse_args()


def as_float(s):
    try:
        return float(s)
    except:
        return float("nan")


def main():
    a = parse_args()
    kv = {}
    with open(a.stats) as f:
        for ln in f:
            if not ln.startswith("SN"):
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            k = parts[1].strip(":").strip()
            v = parts[2].strip()
            kv[k] = v
    er = as_float(kv.get("error rate", "nan"))
    insr = as_float(kv.get("insertion rate", "nan"))
    delr = as_float(kv.get("deletion rate", "nan"))
    mmr = as_float(kv.get("mismatch rate", "nan"))
    ident = (1.0 - er) if (er == er) else float("nan")
    qv = (-10.0 * math.log10(er)) if (er == er and er > 0.0) else float("nan")

    with open(a.out, "w") as w:
        w.write("metric\tvalue\n")
        w.write(f"error_rate\t{er:.8f}\n" if er == er else "error_rate\tNA\n")
        w.write(f"identity\t{ident:.6f}\n" if ident == ident else "identity\tNA\n")
        w.write(f"QV\t{qv:.2f}\n" if qv == qv else "QV\tNA\n")
        w.write(
            f"insertion_rate\t{insr:.8f}\n" if insr == insr else "insertion_rate\tNA\n"
        )
        w.write(
            f"deletion_rate\t{delr:.8f}\n" if delr == delr else "deletion_rate\tNA\n"
        )
        w.write(f"mismatch_rate\t{mmr:.8f}\n" if mmr == mmr else "mismatch_rate\tNA\n")


if __name__ == "__main__":
    main()
