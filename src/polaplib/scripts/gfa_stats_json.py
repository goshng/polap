#!/usr/bin/env python3
# scripts/gfa_stats_json.py
# Version : v0.1.0
# Author  : Sang Chul Choi (POLAP)
# Date    : 2025-10-05
# License : GPL-3.0+
import sys, argparse, gzip, json


def open_any(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "r")


def n50(v):
    if not v:
        return 0
    s = sorted(v, reverse=True)
    tot = sum(s)
    acc = 0
    for x in s:
        acc += x
        if acc >= tot / 2:
            return x
    return s[-1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gfa", required=True)
    a = ap.parse_args()
    nS = nL = 0
    Ls = []
    circ = False
    with open_any(a.gfa) as f:
        for line in f:
            if not line or line[0] == "#":
                continue
            t = line[0]
            if t == "S":
                nS += 1
                p = line.rstrip().split("\t")
                seq = p[2]
                L = len(seq) if seq != "*" else 0
                if L == 0:
                    for x in p[3:]:
                        if x.startswith("LN:i:"):
                            L = int(x.split(":")[-1])
                            break
                Ls.append(L)
            elif t == "L":
                nL += 1
            elif t == "P":
                circ = True
    js = dict(
        n_segments=nS,
        n_links=nL,
        total_len=sum(Ls),
        N50=n50(Ls),
        max_seg=(max(Ls) if Ls else 0),
        is_circular=int(circ),
    )
    print(json.dumps(js, separators=(",", ":")), end="")


if __name__ == "__main__":
    sys.exit(main())
