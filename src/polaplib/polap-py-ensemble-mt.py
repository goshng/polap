#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
polap-py-ensemble-mt.py  v0.3.0

Combine coverage-selected IDs (xcover/xlda) and PPR-selected IDs into a final organelle set.

New:
  --label LABEL         write OUT.<label>.ids, OUT.<label>.nuclear.ids, OUT.<label>.report.tsv
                        (defaults to 'mt'; use 'pt' for plastid pass, etc.)

Inputs:
  --xlda  FILE          coverage-selected IDs (required)
  --ppr   FILE          PPR-selected IDs (optional)
  --mode  union|majority|strict   [majority]

Outputs:
  -o OUTPREFIX
    OUTPREFIX.<label>.ids
    OUTPREFIX.<label>.nuclear.ids
    OUTPREFIX.<label>.report.tsv

Notes:
  * majority with 2 sources is the same as union; kept for parity if you add more sources later.
  * Universe for 'nuclear' is the union of all provided sources (xlda ∪ ppr).
"""

import sys, os, argparse

VERSION = "0.3.0"


def load_ids(path: str):
    S = set()
    if not path:
        return S
    if not os.path.exists(path):
        return S
    with open(path) as fh:
        for line in fh:
            t = line.strip().split()
            if t:
                S.add(t[0])
    return S


def write_ids(path, S):
    with open(path, "w") as w:
        for rid in sorted(S):
            w.write(rid + "\n")


def main():
    ap = argparse.ArgumentParser(
        description="Ensemble combiner for organelle read IDs (coverage + PPR)."
    )
    ap.add_argument("--xlda", required=True, help="coverage-selected IDs (xcover/xlda)")
    ap.add_argument("--ppr", help="PPR-selected IDs (optional)")
    ap.add_argument(
        "--mode",
        choices=["union", "majority", "strict"],
        default="majority",
        help="Ensemble rule: union | majority | strict (intersection)",
    )
    ap.add_argument("-o", "--out", required=True, help="output prefix")
    ap.add_argument(
        "--label", default="mt", help="label for outputs (e.g., mt or pt) [default: mt]"
    )
    ap.add_argument("--version", action="store_true")
    args = ap.parse_args()

    if args.version:
        print(VERSION)
        return

    # Load inputs
    X = load_ids(args.xlda)  # coverage
    P = load_ids(args.ppr) if args.ppr else set()

    # Candidate universe
    U = X | P

    # Decide final set
    if args.mode == "union":
        final = X | P
    elif args.mode == "strict":
        final = X & P if args.ppr else set()  # if no PPR, strict -> empty
    else:
        # majority of 2 sources == union; leave here if you later add a 3rd source
        final = X & P

    nuclear = U - final

    # Outputs (labeled)
    lab = (args.label or "mt").strip()
    out_ids = f"{args.out}.{lab}.ids"
    out_nuc = f"{args.out}.{lab}.nuclear.ids"
    out_rep = f"{args.out}.{lab}.report.tsv"

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    write_ids(out_ids, final)
    write_ids(out_nuc, nuclear)

    # Report
    x_n = len(X)
    p_n = len(P)
    u_n = len(U)
    f_n = len(final)
    n_n = len(nuclear)
    inter = len(X & P)
    uni = len(X | P)
    jacc = (inter / uni) if uni > 0 else 0.0

    with open(out_rep, "w") as w:
        w.write("metric\tvalue\n")
        w.write(f"label\t{lab}\n")
        w.write(f"mode\t{args.mode}\n")
        w.write(f"xlda_count\t{x_n}\n")
        w.write(f"ppr_count\t{p_n}\n")
        w.write(f"intersection\t{inter}\n")
        w.write(f"union\t{uni}\n")
        w.write(f"jaccard\t{jacc:.6f}\n")
        w.write(f"final_count\t{f_n}\n")
        w.write(f"nuclear_count\t{n_n}\n")

    print(f"[done] {out_ids}")
    print(f"[done] {out_nuc}")
    print(f"[done] {out_rep}")


if __name__ == "__main__":
    main()
