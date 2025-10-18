#!/usr/bin/env python3
# scripts/manifest_table_data_min.py
# Version: v0.1.0
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Emit a TSV with columns:
#   code2  species  long_total_bases  short1_total_bases
#
import sys, json, argparse

def fmt_num(x):
    if x is None:
        return ""
    try:
        # normalize 123.0 -> 123
        f = float(x)
        if f.is_integer(): return str(int(f))
        return str(f)
    except Exception:
        return str(x)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--out-tsv",   required=True)
    ap.add_argument("--sort",      default="code2", choices=("code2","species"))
    args = ap.parse_args()

    with open(args.manifest) as fh:
        man = json.load(fh)

    items = man.get("items", [])
    rows = []
    for it in items:
        code2   = it.get("code2", "")
        species = it.get("species", "")
        data    = it.get("data", {}) or {}
        s1      = it.get("short1data", {}) or {}

        long_total   = data.get("total_bases", None)
        short1_total = s1.get("total_bases", None)

        rows.append( (code2, species, long_total, short1_total) )

    if args.sort == "code2":
        rows.sort(key=lambda r: (r[0] if r[0] is not None else ""))
    else:
        rows.sort(key=lambda r: (r[1] if r[1] is not None else ""))

    with open(args.out_tsv, "w", encoding="utf-8") as out:
        out.write("code2\tspecies\tlong_total_bases\tshort1_total_bases\n")
        for code2, sp, lt, s1t in rows:
            out.write(f"{code2}\t{sp}\t{fmt_num(lt)}\t{fmt_num(s1t)}\n")

if __name__ == "__main__":
    sys.exit(main())
