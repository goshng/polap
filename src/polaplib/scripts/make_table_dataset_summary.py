#!/usr/bin/env python3
###############################################################################
# scripts/make_table_dataset_summary.py
#
# Version : v0.1.0
# Author  : Sang Chul Choi (POLAP)
# Date    : 2025-10-05
# License : GPL-3.0+
#
# Purpose :
#   Read the POLAP manifest JSON and export Table S1 (dataset summary)
#   as TSV (and optionally Markdown) for manuscripts.
#
# Input JSON fields (under each item["data"]):
#   total_bases, read_count, mean_length, N50, avg_qual, gc_content
#
# Output TSV columns:
#   species  total_bases  read_count  mean_length  N50  avg_qual  gc_content
#
# Usage :
#   python scripts/make_table_dataset_summary.py \
#       --manifest md/manifest-some.json \
#       --out md/tableS1-dataset-summary.tsv \
#       [--markdown]
###############################################################################
import os, sys, json, csv, argparse
from math import isnan


def safe_float(x):
    try:
        return float(x)
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True, help="Input manifest JSON")
    ap.add_argument("--out", required=True, help="Output TSV (and .md if --markdown)")
    ap.add_argument(
        "--markdown",
        action="store_true",
        default=False,
        help="Also write Markdown table beside TSV",
    )
    a = ap.parse_args()

    if not os.path.exists(a.manifest):
        sys.exit(f"[ERR] manifest not found: {a.manifest}")
    os.makedirs(os.path.dirname(a.out), exist_ok=True)

    with open(a.manifest) as f:
        data = json.load(f)

    items = data.get("items", [])
    if not items:
        sys.exit("[ERR] No 'items' in manifest")

    # collect
    header = [
        "species",
        "total_bases",
        "read_count",
        "mean_length",
        "N50",
        "avg_qual",
        "gc_content",
    ]
    rows = []
    for it in items:
        d = it.get("data", {})
        if not d:  # skip if dataset block missing
            continue
        row = [
            it.get("species", "NA"),
            d.get("total_bases", "NA"),
            d.get("read_count", "NA"),
            d.get("mean_length", "NA"),
            d.get("N50", "NA"),
            d.get("avg_qual", "NA"),
            d.get("gc_content", "NA"),
        ]
        rows.append(row)

    # write TSV
    with open(a.out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)
        w.writerows(rows)
    print(f"[INFO] Wrote TSV: {a.out}")

    # optional Markdown
    if a.markdown:
        mdout = os.path.splitext(a.out)[0] + ".md"
        with open(mdout, "w") as f:
            f.write("| " + " | ".join(header) + " |\n")
            f.write("|" + "|".join(["---"] * len(header)) + "|\n")
            for r in rows:
                f.write("| " + " | ".join(str(x) for x in r) + " |\n")
        print(f"[INFO] Wrote Markdown: {mdout}")


if __name__ == "__main__":
    sys.exit(main())
