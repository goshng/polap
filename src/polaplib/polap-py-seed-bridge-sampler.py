#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
polap-py-seed-bridge-sampler.py

One mapping only:
  1) Build a read–read overlap graph from a single all-vs-all minimap2 PAF.
  2) Take precomputed seed -> reads groupings (no remapping).
  3) For each unordered seed pair, sample random walks between a random read
     from group A and a random read from group B.
  4) Summarize per-read frequencies; high-frequency reads are candidate
     “middle anchors” (often intergenic mtDNA corridors).

Inputs
------
1) allvsall.paf[.gz]    : minimap2 all-vs-all PAF (reads vs reads)
2) seed_groups.tsv      : two-column TSV: seed_id<TAB>read_id (one line per read)

Outputs (in --outdir)
---------------------
- edges.tsv                 : final filtered edges (u v weight)
- graph.gexf                : graph for Gephi (optional)
- pairs/<A>__<B>/walk_*.tsv : sampled successful paths for each seed pair
- pairs/<A>__<B>/summary.tsv: per-pair ranked paths
- read_frequencies.tsv      : read_id, count, fraction
- highfreq_reads.txt        : above thresholds (count or top fraction)
- freq_hist.tsv             : histogram of read counts

References
----------
- Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences.
  Bioinformatics, 34(18), 3094–3100.
- PAF field semantics: minimap2 docs (fields 0..11: qname, qlen, ..., nmatch, alen, mapq).
"""

import sys, os, argparse, gzip, math, random, itertools
from collections import defaultdict, Counter
from typing import Dict, List, Tuple, Optional
import networkx as nx

# ────────────────────────────────────────────────────────────────────────────────
# I/O helpers


def _open(path, mode="rt"):
    if path == "-" and "r" in mode:
        return sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, mode=mode)
    return open(path, mode=mode)


# ────────────────────────────────────────────────────────────────────────────────
# PAF parsing & graph build


def parse_paf_line(line: str):
    p = line.rstrip("\n").split("\t")
    if len(p) < 12:
        return None
    try:
        qname = p[0]
        qlen = int(p[1])
        qstart = int(p[2])
        qend = int(p[3])
        strand = p[4]
        tname = p[5]
        tlen = int(p[6])
        tstart = int(p[7])
        tend = int(p[8])
        nmatch = int(p[9])
        alen = int(p[10])
        mapq = int(p[11])
    except Exception:
        return None
    ident = (nmatch / alen) if alen > 0 else 0.0
    qcov = (qend - qstart) / qlen if qlen > 0 else 0.0
    tcov = (tend - tstart) / tlen if tlen > 0 else 0.0
    return {
        "q": qname,
        "qlen": qlen,
        "qcov": qcov,
        "t": tname,
        "tlen": tlen,
        "tcov": tcov,
        "alen": alen,
        "ident": ident,
        "mapq": mapq,
    }


def edge_weight(rec, formula="ident_x_norm", min_len_norm=True):
    ident = rec["ident"]
    alen = rec["alen"]
    qlen = rec["qlen"]
    tlen = rec["tlen"]
    mapq = rec["mapq"]
    denom = min(qlen, tlen) if min_len_norm else max(1, qlen + tlen)
    frac = alen / denom if denom > 0 else 0.0
    if formula == "ident_x_norm":
        return ident * frac
    elif formula == "ident":
        return ident
    elif formula == "alen_norm":
        return frac
    elif formula == "mapq_x_ident":
        return (mapq / 60.0) * ident
    else:
        return ident


def add_edge_max(G: nx.Graph, a: str, b: str, w: float, data: dict):
    u, v = (a, b) if a <= b else (b, a)
    if u == v:
        return
    if G.has_edge(u, v):
        if w > G[u][v]["weight"]:
            G[u][v]["weight"] = w
            G[u][v]["data"] = data
    else:
        G.add_edge(u, v, weight=w, data=data)


def paf_to_graph(
    paf_path: str, min_alen: int, min_ident: float, min_weight: float, formula: str
) -> Tuple[nx.Graph, int]:
    G = nx.Graph()
    n_in = 0
    with _open(paf_path, "rt") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            rec = parse_paf_line(line)
            if rec is None:
                continue
            n_in += 1
            if rec["alen"] < min_alen:
                continue
            if rec["ident"] < min_ident:
                continue
            w = edge_weight(rec, formula=formula)
            if w < min_weight:
                continue
            add_edge_max(G, rec["q"], rec["t"], w, rec)
    return G, n_in


def write_edges(G: nx.Graph, out_tsv: str):
    with open(out_tsv, "wt") as g:
        print("#u\tv\tweight", file=g)
        for u, v, d in G.edges(data=True):
            print(f"{u}\t{v}\t{d['weight']:.6f}", file=g)


# ────────────────────────────────────────────────────────────────────────────────
# Seed groups


def load_seed_groups(path: str) -> Dict[str, List[str]]:
    """
    seed_groups.tsv: seed_id<TAB>read_id
    Multiple lines per seed, one read per line.
    """
    groups = defaultdict(list)
    with _open(path, "rt") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            seed, read = parts[0], parts[1]
            groups[seed].append(read)
    return groups


# ────────────────────────────────────────────────────────────────────────────────
# Random walk (temperature-biased) from source read to target read


def softmax_probs(weight_list: List[float], T: float) -> List[float]:
    if T <= 0:
        # greedy
        m = max(weight_list)
        return [1.0 if w == m else 0.0 for w in weight_list]
    xs = [w / T for w in weight_list]
    m = max(xs)
    exps = [math.exp(x - m) for x in xs]
    s = sum(exps)
    return [e / s for e in exps] if s > 0 else [1.0 / len(exps)] * len(exps)


def random_walk_once(
    G: nx.Graph,
    src: str,
    dst: str,
    rng: random.Random,
    T: float,
    min_edge_w: float,
    max_hops: int,
    simple: bool,
    p_restart: float,
) -> Tuple[bool, List[str], float, int]:
    if src not in G or dst not in G:
        return (False, [src], float("-inf"), 0)

    path = [src]
    visited = {src} if simple else set()
    cur = src
    logsum = 0.0
    hops = 0

    while True:
        if cur == dst:
            return (True, path, logsum, hops)
        if hops >= max_hops:
            return (False, path, logsum, hops)

        # restart to src?
        if rng.random() < p_restart:
            path = [src]
            visited = {src} if simple else set()
            cur = src
            logsum = 0.0
            hops = 0
            continue

        neighs = []
        ws = []
        for nbr in G.neighbors(cur):
            if simple and nbr in visited:
                continue
            w = G[cur][nbr]["weight"]
            if w >= min_edge_w:
                neighs.append(nbr)
                ws.append(w)
        if not neighs:
            return (False, path, logsum, hops)

        # gentle bias: if destination is adjacent, nudge its weight a bit
        for i, nbr in enumerate(neighs):
            if nbr == dst:
                ws[i] = max(ws[i], min(1.0, ws[i] * 1.05))

        ps = softmax_probs(ws, T)
        # sample next
        r = rng.random()
        acc = 0.0
        idx = len(ps) - 1
        for i, p in enumerate(ps):
            acc += p
            if r <= acc:
                idx = i
                break

        nxt = neighs[idx]
        w = ws[idx]
        path.append(nxt)
        if simple:
            visited.add(nxt)
        logsum += math.log(max(w, 1e-300))
        cur = nxt
        hops += 1


# ────────────────────────────────────────────────────────────────────────────────
# Main driver: per-pair sampling + global frequency summary


def main():
    ap = argparse.ArgumentParser(
        description="Sample random read paths between seed read groups (one all-vs-all PAF only).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument(
        "paf", help="all-vs-all minimap2 PAF of reads vs reads (.paf or .paf.gz)"
    )
    ap.add_argument("seed_groups", help="TSV seed_id<TAB>read_id (one read per line)")
    ap.add_argument("-o", "--outdir", default="seed_bridge", help="output directory")

    # Build graph
    ap.add_argument(
        "--formula",
        choices=["ident_x_norm", "ident", "alen_norm", "mapq_x_ident"],
        default="ident_x_norm",
        help="edge weight formula",
    )
    ap.add_argument(
        "--min-alen", type=int, default=600, help="edge filter: minimum aligned length"
    )
    ap.add_argument(
        "--min-ident", type=float, default=0.86, help="edge filter: minimum identity"
    )
    ap.add_argument(
        "--min-weight", type=float, default=0.12, help="edge filter: minimum weight"
    )

    # Pair sampling
    ap.add_argument(
        "--walks-per-pair",
        type=int,
        default=100,
        help="number of walk attempts per seed pair",
    )
    ap.add_argument("--max-hops", type=int, default=40, help="maximum hops per walk")
    ap.add_argument(
        "--temperature",
        type=float,
        default=0.25,
        help="softmax temperature (lower=greedier)",
    )
    ap.add_argument(
        "--min-edge-w-walk", type=float, default=0.12, help="edge floor during walking"
    )
    ap.add_argument(
        "--simple", action="store_true", help="simple paths only (no node revisits)"
    )
    ap.add_argument(
        "--restart", type=float, default=0.02, help="per-step restart probability"
    )
    ap.add_argument("--seed", type=int, default=13, help="random seed")

    # Frequency thresholding
    ap.add_argument(
        "--min-count", type=int, default=5, help="min count to call high-frequency"
    )
    ap.add_argument(
        "--top-frac", type=float, help="top fraction (0–1) to call high-frequency"
    )

    # Optional: restrict to subset of seed IDs (comma-separated)
    ap.add_argument(
        "--only-seeds", help="Comma-separated seed IDs to include (others ignored)"
    )

    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    rng = random.Random(args.seed)

    # 1) Build overlap graph once
    G, n_in = paf_to_graph(
        args.paf, args.min_alen, args.min_ident, args.min_weight, args.formula
    )
    print(
        f"[info] PAF rows={n_in} edges={G.number_of_edges()} nodes={G.number_of_nodes()}",
        file=sys.stderr,
    )

    edges_tsv = os.path.join(args.outdir, "edges.tsv")
    write_edges(G, edges_tsv)
    try:
        nx.write_gexf(G, os.path.join(args.outdir, "graph.gexf"))
    except Exception as e:
        print(f"[warn] cannot write GEXF: {e}", file=sys.stderr)

    # 2) Load seed -> reads and intersect with graph nodes
    groups_all = load_seed_groups(args.seed_groups)
    if args.only_seeds:
        keep = set([x.strip() for x in args.only_seeds.split(",") if x.strip()])
        groups_all = {k: v for k, v in groups_all.items() if k in keep}

    nodes = set(G.nodes())
    groups = {}
    for sid, reads in groups_all.items():
        in_graph = [r for r in reads if r in nodes]
        if len(in_graph) >= 1:
            groups[sid] = in_graph

    seeds = sorted(groups.keys())
    print(
        f"[info] seeds_with_reads_in_graph={len(seeds)} (of {len(groups_all)})",
        file=sys.stderr,
    )
    if len(seeds) < 2:
        print(
            "[error] need at least two seed groups covering nodes in the graph",
            file=sys.stderr,
        )
        sys.exit(2)

    # 3) For each unordered seed pair, sample walks
    freq_counter = Counter()
    pairs_dir = os.path.join(args.outdir, "pairs")
    os.makedirs(pairs_dir, exist_ok=True)

    total_pairs = 0
    for a, b in itertools.combinations(seeds, 2):
        grpA = groups[a]
        grpB = groups[b]
        if not grpA or not grpB:
            continue
        total_pairs += 1
        pair_tag = f"{a}__{b}"
        outp = os.path.join(pairs_dir, pair_tag)
        os.makedirs(outp, exist_ok=True)

        successes = []
        for i in range(1, args.walks_per_pair + 1):
            # random endpoints from each seed group
            src = rng.choice(grpA)
            dst = rng.choice(grpB)
            ok, nodes_path, logsum, hops = random_walk_once(
                G,
                src,
                dst,
                rng=rng,
                T=args.temperature,
                min_edge_w=args.min_edge_w_walk,
                max_hops=args.max_hops,
                simple=args.simple,
                p_restart=args.restart,
            )
            if ok:
                successes.append((logsum, hops, nodes_path))
                # update global frequencies
                freq_counter.update(nodes_path)
                # write the walk
                with open(os.path.join(outp, f"walk_{i:05d}.tsv"), "wt") as fo:
                    print("#idx\tu\tv\tweight", file=fo)
                    for j in range(len(nodes_path) - 1):
                        u, v = nodes_path[j], nodes_path[j + 1]
                        w = G[u][v]["weight"]
                        print(f"{j+1}\t{u}\t{v}\t{w:.6f}", file=fo)

        # per-pair summary
        successes.sort(key=lambda x: x[0], reverse=True)
        with open(os.path.join(outp, "summary.tsv"), "wt") as so:
            print("rank\tlogsum\tedges\tnodes_csv", file=so)
            for rnk, (logsum, hops, nodes_path) in enumerate(successes, start=1):
                print(f"{rnk}\t{logsum:.6f}\t{hops}\t{','.join(nodes_path)}", file=so)

        print(
            f"[pair] {a} vs {b}: walks={args.walks_per_pair} success={len(successes)}",
            file=sys.stderr,
        )

    if total_pairs == 0:
        print(
            "[error] no valid seed pairs (after filtering) to sample", file=sys.stderr
        )
        sys.exit(3)

    # 4) Global frequency summary
    total_successful_walks = sum(
        1 for _ in freq_counter.elements()
    )  # counts all nodes across walks
    # Build frequency table by read (count == times read appears across all successful paths)
    # We also need number of successful paths to compute fractions by 'presence in paths' if desired.
    # Here, we compute 'fraction by walks': count(read)/sum_over_reads_on_paths? Instead, more interpretable:
    # fraction = count / total_paths_seen_nodes. Keep it simple & monotonic.
    total_node_visits = sum(freq_counter.values())
    rows = []
    for rid, cnt in freq_counter.most_common():
        frac = (cnt / total_node_visits) if total_node_visits > 0 else 0.0
        rows.append((rid, cnt, frac))

    import csv

    with open(os.path.join(args.outdir, "read_frequencies.tsv"), "wt", newline="") as g:
        w = csv.writer(g, delimiter="\t")
        w.writerow(["read_id", "count", "fraction"])
        w.writerows(rows)

    # histogram
    hist = Counter(freq_counter.values())
    with open(os.path.join(args.outdir, "freq_hist.tsv"), "wt") as gh:
        print("count\tn_reads", file=gh)
        for k in sorted(hist.keys()):
            print(f"{k}\t{hist[k]}", file=gh)

    # high-frequency set
    hf_reads = []
    if args.top_frac is not None:
        topN = max(1, int(len(rows) * args.top_frac))
        hf_reads = [rid for rid, _, _ in rows[:topN]]
    else:
        hf_reads = [rid for rid, cnt, _ in rows if cnt >= args.min_count]

    with open(os.path.join(args.outdir, "highfreq_reads.txt"), "wt") as hf:
        for rid in hf_reads:
            print(rid, file=hf)

    print(
        f"[done] pairs={total_pairs} nodes_seen={len(freq_counter)} outdir={args.outdir}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
