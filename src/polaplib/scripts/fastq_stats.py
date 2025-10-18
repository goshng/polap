#!/usr/bin/env python3
# Version: v0.2.0
import sys, gzip, statistics as st


def open_auto(p, m="rt"):
    return gzip.open(p, m) if p.endswith(".gz") else open(p, m)


def lengths(p):
    L = []
    with open_auto(p, "rt") as f:
        i = 0
        for line in f:
            i += 1
            if i % 4 == 2:
                L.append(len(line.rstrip()))
    return L


if len(sys.argv) != 2:
    sys.exit("Usage: fastq_stats.py in.fq[.gz]")
L = lengths(sys.argv[1])
print("reads\tbases\tmean\tmedian\tN50\tp95\tmax")
if not L:
    print("0\t0\t0\t0\t0\t0\t0")
    sys.exit(0)
L.sort()
reads = len(L)
bases = sum(L)


def N50(arr, total):
    half = total / 2
    acc = 0
    for l in sorted(arr, reverse=True):
        acc += l
        if acc >= half:
            return l


print(
    f"{reads}\t{bases}\t{bases/reads:.1f}\t{st.median(L)}\t{N50(L,bases)}\t{int(st.quantiles(L,n=100)[94])}\t{max(L)}"
)
