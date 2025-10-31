#!/usr/bin/env python3
# Version: v0.2.0
# Pick ONE greedy high-ec:i path across an (organelle) subgraph (numeric IDs, simple CIGAR).
import sys, argparse, gzip
from collections import defaultdict


def openg(p, mode="rt"):
    return gzip.open(p, mode) if p.endswith(".gz") else open(p, mode)


def parse(gfa):
    adj = defaultdict(list)
    indeg = defaultdict(int)
    outdeg = defaultdict(int)
    nodes = set()
    with openg(gfa) as f:
        for ln in f:
            if not ln or ln[0] in "#\n":
                continue
            c = ln.rstrip("\n").split("\t")
            if c[0] == "S":
                nid = c[1]
                nodes.add((nid, "+"))
                nodes.add((nid, "-"))
            elif c[0] == "L":
                u, uo, v, vo, cg = c[1:6]
                ec = 1
                for t in c[6:]:
                    if t.startswith("ec:i:"):
                        try:
                            ec = int(t.split(":", 2)[2])
                        except:
                            ec = 1
                        break
                # overlap length:
                ol = 0
                if cg.endswith("M"):
                    try:
                        ol = int(cg[:-1])
                    except:
                        ol = 0
                adj[(u, uo)].append((v, vo, ec, ol))
                outdeg[(u, uo)] += 1
                indeg[(v, vo)] += 1
    return nodes, adj, indeg, outdeg


def tips(nodes, indeg, outdeg):
    t = [n for n in nodes if indeg.get(n, 0) == 0 and outdeg.get(n, 0) > 0]
    if t:
        return t
    cand = [n for n in nodes if outdeg.get(n, 0) > 0]
    return sorted(cand, key=lambda x: (outdeg.get(x, 0), indeg.get(x, 0)))[:5]


def walk_greedy(start, adj, indeg):
    p, seen, cur = [], set(), start
    while True:
        if cur in seen:
            break
        seen.add(cur)
        p.append(cur)
        succs = adj.get(cur, [])
        if not succs:
            break
        if len(succs) == 1:
            cur = (succs[0][0], succs[0][1])
            continue
        succs.sort(key=lambda x: (-x[2], -x[3], indeg.get((x[0], x[1]), 0)))
        cur = (succs[0][0], succs[0][1])
    return p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gfa", required=True)
    ap.add_argument("--out-path", required=True)
    ap.add_argument("--min-len", type=int, default=1)
    ap.add_argument("--greedy", action="store_true")
    a = ap.parse_args()

    nodes, adj, indeg, outdeg = parse(a.gfa)
    starts = tips(nodes, indeg, outdeg)
    best = []
    for s in starts:
        p = walk_greedy(s, adj, indeg)
        if len(p) > len(best):
            best = p
    if not best:
        sys.exit("ERR: could not find path")
    toks = [f"{n}{o}" for (n, o) in best]
    with open(a.out_path, "w") as w:
        w.write(",".join(toks) + "\n")


if __name__ == "__main__":
    main()
