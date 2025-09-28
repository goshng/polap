#!/usr/bin/env python3
# polap-py-threshold-from-nuclear.py  v0.0.1
# Input:
#   overlapness.tsv (read_id, degree, wdegree)
#   nuc.ids         (one read id per line)
# Output:
#   shell lines: wdeg_min=..., deg_min=...
#   writes: threshold_from_nuclear.tsv (diagnostic)
import sys, argparse, statistics as st


def load_nucs(path):
    S = set()
    with open(path) as f:
        for ln in f:
            ln = ln.strip()
            if ln:
                S.add(ln)
    return S


def load_ov(path):
    D = {}
    with open(path) as f:
        for ln in f:
            if not ln.strip() or ln[0] == "#":
                continue
            rid, d, w = ln.rstrip("\n").split("\t")
            D[rid] = (int(d), float(w))
    return D


def perc(xs, q):
    if not xs:
        return 0.0
    xs = sorted(xs)
    k = max(0, min(len(xs) - 1, int(round(q * (len(xs) - 1)))))
    return xs[k]


def main():
    ap = argparse.ArgumentParser(
        description="Derive degree/wdegree cutoffs from nuclear labels"
    )
    ap.add_argument("overlapness")
    ap.add_argument("nuc_ids")
    ap.add_argument("--mode", choices=["fpr", "perc", "otsu"], default="fpr")
    ap.add_argument(
        "--fpr", type=float, default=0.01, help="target nuclear FPR (for mode=fpr)"
    )
    ap.add_argument(
        "--q", type=float, default=0.99, help="percentile for mode=perc (e.g., 0.99)"
    )
    ap.add_argument("--diag", default="threshold_from_nuclear.tsv")
    args = ap.parse_args()

    nuc = load_nucs(args.nuc_ids)
    ov = load_ov(args.overlapness)

    n_deg = [ov[r][0] for r in ov if r in nuc]
    n_wdg = [ov[r][1] for r in ov if r in nuc]
    a_deg = [ov[r][0] for r in ov]
    a_wdg = [ov[r][1] for r in ov]

    if args.mode == "fpr":
        # cutoff s.t. fraction of nuclear ≥ cutoff <= fpr
        # i.e., wdeg_min = nuclear q_(1-fpr)
        wdeg_min = perc(n_wdg, 1.0 - args.fpr) if n_wdg else 0.0
        deg_min = int(perc(n_deg, 1.0 - args.fpr)) if n_deg else 0
        # nudge just above to be strict
        wdeg_min = round(wdeg_min + 1e-6, 6)
        deg_min = int(deg_min + 1)
    elif args.mode == "perc":
        wdeg_min = perc(n_wdg, args.q) if n_wdg else 0.0
        deg_min = int(perc(n_deg, args.q)) if n_deg else 0
        wdeg_min = round(wdeg_min + 1e-6, 6)
        deg_min = int(deg_min + 1)
    else:  # otsu over all wdegree, then lift above nuclear tail
        # simple otsu on [min..max] with 200 bins
        mn, mx = (min(a_wdg), max(a_wdg)) if a_wdg else (0.0, 1.0)
        nb = 200
        step = (mx - mn) / nb if mx > mn else 1.0
        hist = [0] * nb
        edges = [mn + i * step for i in range(nb + 1)]
        for x in a_wdg:
            i = min(nb - 1, max(0, int((x - mn) / step))) if step > 0 else 0
            hist[i] += 1
        # otsu
        tot = sum(hist)
        sum_all = sum(((edges[i] + edges[i + 1]) * 0.5) * c for i, c in enumerate(hist))
        wB = 0
        sumB = 0
        varMax = -1
        thr = edges[1] if nb > 1 else mn
        for i, c in enumerate(hist):
            wB += c
            if wB == 0:
                continue
            wF = tot - wB
            if wF == 0:
                break
            mB = (sumB + ((edges[i] + edges[i + 1]) * 0.5) * c) / wB
            mF = (sum_all - sumB - ((edges[i] + edges[i + 1]) * 0.5) * c) / wF
            var = wB * wF * (mB - mF) * (mB - mF)
            if var > varMax:
                varMax = var
                thr = (edges[i] + edges[i + 1]) * 0.5
        wdeg_min = thr
        # ensure below a small fraction of nuclear
        tail = perc(n_wdg, 0.99) if n_wdg else wdeg_min
        wdeg_min = max(wdeg_min, tail + 1e-6)
        # degree: align to nuclear 99th perc
        deg_min = int(perc(n_deg, 0.99)) + 1 if n_deg else 0

    # diagnostics
    with open(args.diag, "w") as g:
        g.write("#metric\tvalue\n")
        g.write(f"n_nuclear\t{len(n_wdg)}\n")
        g.write(f"n_all\t{len(a_wdg)}\n")
        g.write(f"wdeg_min\t{wdeg_min}\n")
        g.write(f"deg_min\t{deg_min}\n")

    print(f"wdeg_min={wdeg_min}")
    print(f"deg_min={deg_min}")


if __name__ == "__main__":
    main()
