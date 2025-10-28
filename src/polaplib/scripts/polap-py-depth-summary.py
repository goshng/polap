#!/usr/bin/env python3
# Name: polap-py-depth-summary.py v0.2.0
# Purpose: Summarize samtools depth TSV into mean/median/MAD and zero-coverage fraction.
#
# Input:  --depth depth.tsv  (columns: contig  pos  depth)
# Output: --out   TSV with "metric\tvalue"
import sys, argparse, statistics as st


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--depth", required=True)
    ap.add_argument("--out", required=True)
    return ap.parse_args()


def main():
    a = parse_args()
    vals = []
    with open(a.depth) as f:
        for ln in f:
            if not ln.strip():
                continue
            try:
                v = int(ln.split("\t")[2])
                vals.append(v)
            except Exception:
                pass
    with open(a.out, "w") as w:
        w.write("metric\tvalue\n")
        if not vals:
            w.write(
                "mean_depth\tNA\nmedian_depth\tNA\nmad_depth\tNA\nzero_cov_fraction\tNA\n"
            )
            return
        mean = sum(vals) / len(vals)
        med = st.median(vals)
        mad = st.median([abs(v - med) for v in vals])
        zeros = sum(v == 0 for v in vals)
        w.write(f"mean_depth\t{mean:.4f}\n")
        w.write(f"median_depth\t{med:.4f}\n")
        w.write(f"mad_depth\t{mad:.4f}\n")
        w.write(f"zero_cov_fraction\t{zeros/len(vals):.6f}\n")


if __name__ == "__main__":
    main()
