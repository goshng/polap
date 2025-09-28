#!/usr/bin/env python3
# polap-py-pt-ident-threshold.py  v0.3.0
# Compute an identity cutoff for plastid read selection from 1+ PAFs,
# using PT-origin, MT-origin and/or nuclear-origin read-id sets.
#
# Decision order (identity is nmatch/alen; MAPQ is IGNORED now):
#   • If PT + MT labels exist (recommended):
#       ident_min = min(  PT 5th percentile,
#                         MT 99th percentile + delta )
#     (protect MT, while still removing ≥~95% PT)
#   • Else if PT + Nuclear exist (legacy):
#       ident_min = max( Nuclear 99th percentile,
#                        PT 5th percentile )
#   • Else if PT only:
#       ident_min = max( PT 5th percentile,
#                        median(PT) - 2*MAD(PT),
#                        0.88 )
#   • Else if Nuclear only:
#       ident_min = Nuclear 99th percentile
#   Clamp to [0.85, 0.99]. Apply optional aligned-length guard (alen_min).
#
# Emits:
#   stdout: ident_min=... \n alen_min=... \n
#   --diag: small TSV with quantiles and mode
#   --emit-pt-ids: union qnames from ALL PAF rows passing (ident>=cut && alen>=alen_min)

import sys, argparse, gzip, statistics as st


def openg(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "rt")


def read_ids(path):
    S = set()
    with open(path) as f:
        for ln in f:
            ln = ln.strip()
            if ln:
                S.add(ln)
    return S


def best_by_mapq_alen(paf_paths):
    # dict[q] = (ident, alen, mapq) but MAPQ is not used for the cutoff
    best = {}
    for paf in paf_paths:
        with openg(paf) as f:
            for ln in f:
                if not ln or ln[0] == "#":
                    continue
                p = ln.rstrip("\n").split("\t")
                if len(p) < 12:
                    continue
                q = p[0]
                nmatch = int(p[9])
                alen = int(p[10])
                mapq = int(p[11])
                if alen <= 0:
                    continue
                ident = nmatch / alen
                if (
                    q not in best
                    or mapq > best[q][2]
                    or (mapq == best[q][2] and alen > best[q][1])
                ):
                    best[q] = (ident, alen, mapq)
    return best


def pct(xs, q):
    xs = sorted(xs)
    if not xs:
        return None
    if q <= 0:
        return xs[0]
    if q >= 1:
        return xs[-1]
    k = int(round(q * (len(xs) - 1)))
    return xs[k]


def mad(vals, center=None):
    if not vals:
        return 0.0
    if center is None:
        center = st.median(vals)
    dev = [abs(x - center) for x in vals]
    return st.median(dev)


def clamp(x, lo=0.15, hi=0.99):
    return max(lo, min(hi, x))


def main():
    ap = argparse.ArgumentParser(
        description="PT identity cutoff from PAFs with PT/MT/Nuclear labels"
    )
    ap.add_argument(
        "--paf", nargs="+", required=True, help="PAF files (e.g., formA.paf formB.paf)"
    )
    ap.add_argument("--pt-ids", required=False, help="PT-origin read IDs")
    ap.add_argument("--mt-ids", required=False, help="MT-origin read IDs")
    ap.add_argument("--nuc-ids", required=False, help="Nuclear-origin read IDs")
    ap.add_argument(
        "--alen-min", type=int, default=0, help="Aligned length guard for selection [0]"
    )
    ap.add_argument(
        "--delta-mt",
        type=float,
        default=0.005,
        help="Margin added to MT 99th (e.g., 0.005)",
    )
    ap.add_argument(
        "--fpr", type=float, default=0.01, help="Legacy: nuclear FPR target [0.01]"
    )
    ap.add_argument(
        "--tpr", type=float, default=0.95, help="Legacy: PT TPR target [0.95]"
    )
    ap.add_argument("--diag", default=None, help="Write diagnostics TSV here")
    ap.add_argument(
        "--emit-pt-ids",
        default=None,
        help="Emit union PT ids (apply gates) to this file",
    )
    args = ap.parse_args()

    best = best_by_mapq_alen(args.paf)

    S_pt = read_ids(args.pt_ids) if args.pt_ids else set()
    S_mt = read_ids(args.mt_ids) if args.mt_ids else set()
    S_nuc = read_ids(args.nuc_ids) if args.nuc_ids else set()

    id_pt = [best[q][0] for q in S_pt if q in best]
    id_mt = [best[q][0] for q in S_mt if q in best]
    id_nuc = [best[q][0] for q in S_nuc if q in best]

    mode = "fallback"
    ident_min = 0.92

    # (1) MT-guided (preferred when both PT and MT exist)
    if id_pt and id_mt:
        pt_p05 = pct(id_pt, 0.05)
        mt_p99 = pct(id_mt, 0.99)
        if pt_p05 is None and mt_p99 is None:
            ident_min = 0.92
            mode = "fallback"
        else:
            lower = (mt_p99 + args.delta_mt) if mt_p99 is not None else 0.0
            upper = pt_p05 if pt_p05 is not None else 0.99
            ident_min = clamp(min(upper, lower))
            mode = "pt+mt"
    # (2) Legacy nuclear-guided (if PT + Nuclear are present, but no MT)
    elif id_pt and id_nuc:
        nuc_p99 = pct(id_nuc, 1.0 - args.fpr)  # ~0.99
        pt_p05 = pct(id_pt, 1.0 - args.tpr)  # 0.05 if tpr=0.95
        lower = nuc_p99 if nuc_p99 is not None else 0.0
        upper = pt_p05 if pt_p05 is not None else 0.99
        ident_min = clamp(max(lower, upper))
        mode = "pt+nuc"
    # (3) PT-only robust
    elif id_pt:
        pt_p05 = pct(id_pt, 0.05)
        pt_med = st.median(sorted(id_pt))
        pt_mad = mad(id_pt, pt_med)
        ident_min = clamp(
            max(pt_p05 if pt_p05 is not None else 0.0, pt_med - 2 * pt_mad, 0.88)
        )
        mode = "pt_only"
    # (4) Nuclear-only
    elif id_nuc:
        nuc_p99 = pct(id_nuc, 1.0 - args.fpr)
        ident_min = clamp(nuc_p99 if nuc_p99 is not None else 0.92)
        mode = "nuc_only"
    else:
        ident_min = 0.92
        mode = "fallback"

    # Emit shell-friendly lines
    print(f"ident_min={ident_min:.3f}")
    print(f"alen_min={args.alen_min}")

    # Diagnostics
    if args.diag:
        with open(args.diag, "w") as g:
            g.write("#metric\tvalue\n")
            g.write(f"mode\t{mode}\n")
            g.write(f"n_best\t{len(best)}\n")
            g.write(f"n_pt_label\t{len(S_pt)}\n")
            g.write(f"n_pt_in_best\t{len(id_pt)}\n")
            g.write(f"n_mt_label\t{len(S_mt)}\n")
            g.write(f"n_mt_in_best\t{len(id_mt)}\n")
            g.write(f"n_nuc_label\t{len(S_nuc)}\n")
            g.write(f"n_nuc_in_best\t{len(id_nuc)}\n")
            if id_pt:
                g.write(f"id_pt_p05\t{(pct(id_pt,0.05) or 0):.5f}\n")
                g.write(f"id_pt_median\t{st.median(sorted(id_pt)):.5f}\n")
                g.write(f"id_pt_mad\t{mad(id_pt):.5f}\n")
            if id_mt:
                g.write(f"id_mt_p99\t{(pct(id_mt,0.99) or 0):.5f}\n")
            if id_nuc:
                g.write(f"id_nuc_p99\t{(pct(id_nuc,0.99) or 0):.5f}\n")
            g.write(f"ident_min\t{ident_min:.5f}\n")
            g.write(f"alen_min\t{args.alen_min}\n")

    # Optionally emit selected PT ids from ALL PAF rows (union)
    if args.emit_pt_ids:
        out_path = args.emit_pt_ids
        keep = set()
        for paf in args.paf:
            with openg(paf) as f:
                for ln in f:
                    if not ln or ln[0] == "#":
                        continue
                    p = ln.rstrip("\n").split("\t")
                    if len(p) < 12:
                        continue
                    q = p[0]
                    nmatch = int(p[9])
                    alen = int(p[10])
                    if alen <= 0:
                        continue
                    ident = nmatch / alen
                    if ident >= ident_min and alen >= args.alen_min:
                        keep.add(q)
        with open(out_path, "w") as h:
            for q in sorted(keep):
                h.write(q + "\n")


if __name__ == "__main__":
    main()
