#!/usr/bin/env python3
# Version: v0.2.0
import sys, gzip, argparse


def open_auto(p, m="rt"):
    return gzip.open(p, m) if p.endswith(".gz") else open(p, m)


def fq_iter(h):
    while True:
        H = h.readline()
        if not H:
            return
        S = h.readline()
        P = h.readline()
        Q = h.readline()
        if not Q:
            return
        yield H.rstrip(), S.rstrip(), P.rstrip(), Q.rstrip()


ap = argparse.ArgumentParser()
ap.add_argument("--minlen", type=int, default=2000)
ap.add_argument("infq")
ap.add_argument("target_bases", type=int)
ap.add_argument("outfq")
a = ap.parse_args()
Lens = []
with open_auto(a.infq, "rt") as f:
    idx = 0
    for h, s, p, q in fq_iter(f):
        if not h.startswith("@"):
            continue
        L = len(s)
        if L >= a.minlen:
            Lens.append((L, idx, h.split()[0][1:]))
        idx += 1
if not Lens:
    sys.exit("No reads >= --minlen")
Lens.sort(key=lambda x: (-x[0], x[1]))
budget = a.target_bases
chosen = set()
acc = 0
for L, i, n in Lens:
    if acc >= budget:
        break
    chosen.add(n)
    acc += L
with open_auto(a.infq, "rt") as fin, gzip.open(a.outfq, "wt") as out:
    for h, s, p, q in fq_iter(fin):
        if h.split()[0][1:] in chosen:
            out.write(h + "\n")
            out.write(s + "\n")
            out.write(p + "\n")
            out.write(q + "\n")
sys.stderr.write(f"[INFO] selected={len(chosen)} total_bases>={acc}\n")
