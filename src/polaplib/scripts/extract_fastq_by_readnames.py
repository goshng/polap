#!/usr/bin/env python3
# Version: v0.2.0
import sys, gzip


def open_auto(p, m="rt"):
    return gzip.open(p, m) if p.endswith(".gz") else open(p, m)


def fq_iter(h):
    while True:
        h1 = h.readline()
        if not h1:
            return
        s = h.readline()
        p = h.readline()
        q = h.readline()
        if not q:
            return
        yield h1.rstrip(), s.rstrip(), p.rstrip(), q.rstrip()


if len(sys.argv) != 4:
    sys.exit("Usage: extract_fastq_by_readnames.py in.fq[.gz] names.txt out.fq.gz")
infq, names, outfq = sys.argv[1], sys.argv[2], sys.argv[3]
keep = set(x.strip() for x in open(names) if x.strip())
with open_auto(infq, "rt") as fin, gzip.open(outfq, "wt") as fout:
    for h, s, p, q in fq_iter(fin):
        if h.startswith("@") and h.split()[0][1:] in keep:
            fout.write(h + "\n")
            fout.write(s + "\n")
            fout.write(p + "\n")
            fout.write(q + "\n")
