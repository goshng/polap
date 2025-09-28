#!/usr/bin/env python3
# paf-to-edges-stream.py  v0.1.0
# Stream PAF from stdin, apply graph filters, emit:
#   1) edges.tsv (u v weight) one edge per unordered pair (max weight)
#   2) overlapness.tsv to stdout if --overlapness out.tsv passed
# Options:
#   --min-olen N --min-ident F --w-floor F --topk-node K --overlapness out
import sys, argparse, math, collections

ap = argparse.ArgumentParser()
ap.add_argument("--min-olen", type=int, default=1200)
ap.add_argument("--min-ident", type=float, default=0.84)
ap.add_argument("--w-floor", type=float, default=0.12)
ap.add_argument("--topk-node", type=int, default=0)
ap.add_argument("--edges", required=True)
ap.add_argument("--overlapness", required=True)
args = ap.parse_args()
# store max weight per pair; track deg/wdeg per node
pair2w = {}
wdeg = collections.defaultdict(float)
deg = collections.defaultdict(int)


def upd(u, v, w):
    a, b = (u, v) if u <= v else (v, u)
    old = pair2w.get((a, b), 0.0)
    if w > old:
        pair2w[(a, b)] = w


for ln in sys.stdin:
    if not ln.strip() or ln[0] == "#":
        continue
    p = ln.rstrip("\n").split("\t")
    if len(p) < 12:
        continue
    q = p[0]
    ql = int(p[1])
    qs = int(p[2])
    qe = int(p[3])
    t = p[5]
    tl = int(p[6])
    ts = int(p[7])
    te = int(p[8])
    nm = int(p[9])
    al = int(p[10])
    mq = int(p[11])
    if q == t or al < args.min_olen:
        continue
    ident = nm / float(al)
    if ident < args.min_ident:
        continue
    frac = al / float(min(ql, tl)) if min(ql, tl) > 0 else 0.0
    w = ident * frac
    if w < args.w_floor:
        continue
    upd(q, t, w)
# optional Top-K edges per node
if args.topk_node > 0:
    # build adjacency, keep topK by weight per node
    adj = collections.defaultdict(list)
    for (u, v), w in pair2w.items():
        adj[u].append((v, w))
        adj[v].append((u, w))
    pair2w = {}
    for u, lst in adj.items():
        lst.sort(key=lambda x: x[1], reverse=True)
        for v, w in lst[: args.topk_node]:
            a, b = (u, v) if u <= v else (v, u)
            if (a, b) not in pair2w or w > pair2w[(a, b)]:
                pair2w[(a, b)] = w
# overlapness
for (u, v), w in pair2w.items():
    wdeg[u] += w
    wdeg[v] += w
# degree: count unique partners
for u, v in pair2w.keys():
    deg[u] += 1
    deg[v] += 1
# write edges
with open(args.edges, "w") as g:
    g.write("u\tv\tweight\n")
    for (u, v), w in pair2w.items():
        g.write(f"{u}\t{v}\t{w:.6f}\n")
# write overlapness
with open(args.overlapness, "w") as h:
    h.write("#read_id\tdegree\twdegree\n")
    for u in wdeg.keys():
        h.write(f"{u}\t{deg[u]}\t{wdeg[u]:.6f}\n")
