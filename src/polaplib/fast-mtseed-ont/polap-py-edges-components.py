#!/usr/bin/env python3
# polap-py-edges-components.py  v0.1.1
# Read an edge list (u v weight) [TSV, optionally .gz], compute:
#   • connected components
#   • per-node degree & weighted degree
#   • summary: n_nodes, n_edges, #comps, top1..top5 sizes, Gini(deg,wdeg)
# Writes:
#   --summary   summary TSV (1 row)
#   --sizes     comp sizes TSV (comp_id, size)
#   --membership node->comp TSV (read_id, comp_id)           [optional]
#   --deg       degree/wdegree TSV (read_id, degree, wdegree) [optional]
#   --giant-nodes file with nodes of the largest component    [optional]
#   --giant-edges edge list restricted to the largest comp    [optional]
# Plus:
#   --print     print a small banner with key stats to stderr

import sys, argparse, gzip, collections, math


def openg(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "rt")


def read_edges(path):
    adj = collections.defaultdict(list)
    wdeg = collections.defaultdict(float)
    deg = collections.Counter()
    m = 0
    with openg(path) as f:
        for ln in f:
            if not ln.strip() or ln[0] == "#":
                continue
            u, v, w = ln.rstrip("\n").split("\t")[:3]
            if u == v:
                continue
            try:
                w = float(w)
            except:
                w = 0.0
            adj[u].append((v, w))
            adj[v].append((u, w))
            wdeg[u] += w
            wdeg[v] += w
            deg[u] += 1
            deg[v] += 1
            m += 1
    return adj, set(adj.keys()), deg, wdeg, m


def comps(adj, nodes):
    comp_id = {}, {}
    comp_id, seen = {}, set()
    sizes = []
    cid = 0
    for s in nodes:
        if s in seen:
            continue
        cid += 1
        q = [s]
        seen.add(s)
        comp_id[s] = cid
        sz = 1
        for u in q:
            for v, _ in adj[u]:
                if v not in seen:
                    seen.add(v)
                    comp_id[v] = cid
                    sz += 1
                    q.append(v)
        sizes.append((cid, sz))
    sizes.sort(key=lambda x: x[1], reverse=True)
    return comp_id, sizes


def gini(vals):
    xs = sorted(vals)
    n = len(xs)
    if n == 0:
        return 0.0
    s = sum(xs)
    if s <= 0:
        return 0.0
    cum = 0.0
    for i, x in enumerate(xs, start=1):
        cum += i * x
    return (2 * cum) / (n * s) - (n + 1) / n


def main():
    ap = argparse.ArgumentParser(description="Components & stats from edges.tsv(.gz)")
    ap.add_argument("edges", help="edges.tsv or edges.tsv.gz (u v weight)")
    ap.add_argument("--summary", required=True)
    ap.add_argument("--sizes", required=True)
    ap.add_argument("--membership")
    ap.add_argument("--deg")
    ap.add_argument("--topk-nodes", type=int, default=0)
    ap.add_argument("--giant-nodes")
    ap.add_argument("--giant-edges")
    ap.add_argument(
        "--print", action="store_true", help="print a brief banner to stderr"
    )
    args = ap.parse_args()

    adj, nodes, deg, wdeg, m = read_edges(args.edges)
    n = len(nodes)
    comp_map, sizes = comps(adj, nodes)
    ncomps = len(sizes)
    top = [s for _, s in sizes[:5]] + [0] * max(0, 5 - len(sizes))
    g_deg = gini([deg[v] for v in nodes])
    g_wdeg = gini([wdeg[v] for v in nodes])

    # summary
    with open(args.summary, "w") as g:
        g.write("#metric\tvalue\n")
        g.write(f"n_nodes\t{n}\n")
        g.write(f"n_edges\t{m}\n")
        g.write(f"n_components\t{ncomps}\n")
        g.write(f"top1\t{top[0]}\n")
        g.write(f"top2\t{top[1]}\n")
        g.write(f"top3\t{top[2]}\n")
        g.write(f"top4\t{top[3]}\n")
        g.write(f"top5\t{top[4]}\n")
        g.write(f"gini_degree\t{g_deg:.5f}\n")
        g.write(f"gini_wdegree\t{g_wdeg:.5f}\n")

    # sizes
    with open(args.sizes, "w") as g:
        g.write("comp_id\tsize\n")
        for cid, sz in sizes:
            g.write(f"{cid}\t{sz}\n")

    # membership
    if args.membership:
        with open(args.membership, "w") as g:
            g.write("read_id\tcomp_id\n")
            for v, cid in comp_map.items():
                g.write(f"{v}\t{cid}\n")

    # degree table
    if args.deg:
        with open(args.deg, "w") as g:
            g.write("read_id\tdegree\twdegree\n")
            for v in nodes:
                g.write(f"{v}\t{deg[v]}\t{wdeg[v]:.6f}\n")

    # top-K nodes by wdegree → stdout (if requested)
    if args.topk_nodes and args.topk_nodes > 0:
        topk = sorted(wdeg.items(), key=lambda kv: kv[1], reverse=True)[
            : args.topk_nodes
        ]
        for v, _ in topk:
            print(v)

    # giant component exports
    if args.giant_nodes or args.giant_edges:
        if sizes:
            giant_id = sizes[0][0]
            giant_set = {v for v, cid in comp_map.items() if cid == giant_id}
            if args.giant_nodes:
                with open(args.giant_nodes, "w") as g:
                    for v in sorted(giant_set):
                        g.write(v + "\n")
            if args.giant_edges:
                with openg(args.edges) as f, open(args.giant_edges, "w") as g:
                    g.write("u\tv\tweight\n")
                    for ln in f:
                        if not ln.strip() or ln[0] == "#":
                            continue
                        u, v, w = ln.rstrip("\n").split("\t")[:3]
                        if u in giant_set and v in giant_set:
                            g.write(f"{u}\t{v}\t{w}\n")

    # banner
    if args.print:
        sys.stderr.write(
            f"[edges-components] nodes={n} edges={m} comps={ncomps} "
            f"top={top[:3]} gini_deg={g_deg:.3f} gini_wdeg={g_wdeg:.3f}\n"
        )


if __name__ == "__main__":
    main()
