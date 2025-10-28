#!/usr/bin/env python3
# Version: v0.2.0
# Name: polap-py-compare-polish-metrics.py
# Purpose: Compare polishing metrics between two assemblies (A vs B)
#          using key–value summary TSV files produced by POLAP polishing workflows.
#
# Usage:
#   python3 polap-py-compare-polish-metrics.py A.tsv B.tsv out.tsv
#
# Inputs:
#   A.tsv and B.tsv — two summary tables with the format:
#       metric <TAB> value
#
# Output:
#   out.tsv — side-by-side comparison with columns:
#       metric, A, B, delta(B-A), better
#
# Metric interpretation:
#   Higher-is-better  → identity, QV, mean_depth, median_depth
#   Lower-is-better   → error_rate, mad_depth, zero_cov_fraction
#
# Example:
#   python3 polap-py-compare-polish-metrics.py A/summary.tsv B/summary.tsv compare.tsv
#
# Example output:
#   metric    A         B         delta(B-A)    better
#   identity  0.9985    0.9991    0.0006        B
#   QV        32.4      34.8      2.4000        B
#   error_rate 0.0015   0.0009   -0.0006        B
#   mad_depth  4.2      3.7      -0.5           B

import sys

if len(sys.argv) < 4:
    sys.stderr.write("Usage: polap-py-compare-polish-metrics.py A.tsv B.tsv out.tsv\n")
    sys.exit(2)

A, B, OUT = sys.argv[1], sys.argv[2], sys.argv[3]


def read_kv(path):
    """Read metric\tvalue TSV into a dict."""
    d = {}
    with open(path) as f:
        hdr = f.readline()
        for ln in f:
            if not ln.strip():
                continue
            k, v = ln.rstrip("\n").split("\t", 1)
            d[k] = v
    return d


def to_float(s):
    try:
        return float(s)
    except ValueError:
        return None


a = read_kv(A)
b = read_kv(B)
metrics = sorted(set(a.keys()) | set(b.keys()))

with open(OUT, "w") as w:
    w.write("metric\tA\tB\tdelta(B-A)\tbetter\n")
    for m in metrics:
        av = a.get(m, "NA")
        bv = b.get(m, "NA")
        af = to_float(av)
        bf = to_float(bv)
        delta = "NA"
        better = "NA"

        if af is not None and bf is not None:
            d = bf - af
            delta = f"{d:.6f}"
            # Define metric polarity
            if m in ("error_rate", "mad_depth", "zero_cov_fraction"):
                # lower is better
                if af < bf:
                    better = "A"
                elif bf < af:
                    better = "B"
                else:
                    better = "equal"
            else:
                # higher is better
                if af > bf:
                    better = "A"
                elif bf > af:
                    better = "B"
                else:
                    better = "equal"

        w.write(f"{m}\t{av}\t{bv}\t{delta}\t{better}\n")
