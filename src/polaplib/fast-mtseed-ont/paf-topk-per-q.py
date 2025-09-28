#!/usr/bin/env python3
# paf-topk-per-q.py  v0.1.0
# Keep top-K hits per qname by (MAPQ, then ALEN). Reads PAF on stdin.
import sys, heapq

K = int(sys.argv[1]) if len(sys.argv) > 1 else 0
if K <= 0:
    for ln in sys.stdin:
        sys.stdout.write(ln)
        sys.exit(0)
cur = None
buf = []


def flush():
    if not buf:
        return
    # highest first by (mapq, alen)
    out = sorted(buf, key=lambda x: (x[0], x[1]), reverse=True)[:K]
    for _, __, ln in out:
        sys.stdout.write(ln)


for ln in sys.stdin:
    if not ln.strip() or ln[0] == "#":
        continue
    p = ln.rstrip("\n").split("\t")
    if len(p) < 12:
        continue
    q = p[0]
    al = int(p[10])
    mq = int(p[11])
    if cur is None:
        cur = q
    if q != cur:
        flush()
        buf = []
        cur = q
    buf.append((mq, al, ln))
flush()
