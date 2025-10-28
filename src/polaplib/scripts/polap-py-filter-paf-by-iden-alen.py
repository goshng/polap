#!/usr/bin/env python3
# Version: v0.3.0
# polap-py-filter-paf-by-iden-alen.py
import sys, argparse


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--paf", required=True)
    ap.add_argument("--min-ident", type=float, default=0.91)
    ap.add_argument("--min-alen", type=int, default=3000)
    ap.add_argument("--primary-only", action="store_true", default=True)
    ap.add_argument("--out", required=True)
    return ap.parse_args()


def is_primary(cols):
    # minimap2 emits tp:A:P/S; keep only P (primary)
    for t in cols[12:]:
        if t.startswith("tp:A:"):
            return t.endswith(":P")
    return False


def main():
    a = parse_args()
    kept = 0
    with open(a.paf) as f, open(a.out, "w") as w:
        for ln in f:
            if not ln.strip() or ln[0] == "#":
                continue
            c = ln.rstrip("\n").split("\t")
            if len(c) < 12:
                continue
            try:
                nmatch = int(c[9])
                alen = int(c[10])
            except ValueError:
                continue
            if a.primary_only and not is_primary(c):
                continue
            if alen < a.min_alen:
                continue
            ident = (nmatch / alen) if alen > 0 else 0.0
            if ident < a.min_ident:
                continue
            w.write(ln)
            kept += 1
    sys.stderr.write(
        f"[paf-filter] kept {kept} overlaps (primary_only={a.primary_only}, ident>={a.min_ident}, alen>={a.min_alen})\n"
    )


if __name__ == "__main__":
    main()
