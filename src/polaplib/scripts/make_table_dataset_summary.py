#!/usr/bin/env python3
# -*- coding: utf-8 -*-
###############################################################################
# scripts/make_table_dataset_summary.py
#
# Version : v0.2.1
# Author  : Sang Chul Choi (POLAP)
# Date    : 2025-12-04
# Date    : 2025-10-22 v0.2.0
# License : GPL-3.0+
#
# Purpose :
#   Read the POLAP manifest JSON and export Table S1 (dataset summary)
#   as TSV (and optionally Markdown) for manuscripts.
#
# Input JSON schema (per item):
#   {
#     "code2": "Aa01",             # used as Species in the table
#     "species": "Genus_species"   # Species name
#     "data": {
#       "total_bases": <number|string>,   # raw bases
#       "read_count":  <number|string>,   # reads -> delete
#       "mean_length": <number|string>,   # average read length (bases)
#       "N50":         <number|string>,   # N50 (bases)
#       "avg_qual":    <number|string>    # average quality (float)
#       "sra_id":      <number|string>    # SRA
#     }
#   }
#
# Output TSV columns (UNITS APPLIED as requested):
#   SpeciesCode  Species  SRA  Bases  Reads  AvgLen  N50  AvgQ
#   - Species: string
#   - SRA    : string
#   - Bases  : gigabases, 1 decimal (e.g., 7.7)
#   - AvgLen : kilobases, 1 decimal (e.g., 11.4)
#   - N50    : kilobases, 1 decimal (e.g., 20.1)
#   - AvgQ   : integer (rounded, e.g., 9)
#
# Usage :
#   python scripts/make_table_dataset_summary.py \
#       --manifest md/manifest-some.json \
#       --out md/tableS1-dataset-summary.tsv \
#       [--markdown]
###############################################################################
import os
import sys
import json
import csv
import argparse
from typing import Any, Optional


def num(x: Any) -> Optional[float]:
    """Parse x to float, accepting strings with commas; return None if not parseable."""
    if x is None:
        return None
    if isinstance(x, (int, float)):
        return float(x)
    if isinstance(x, str):
        s = x.strip().replace(",", "")
        if s == "":
            return None
        try:
            return float(s)
        except Exception:
            return None
    return None


def fmt_gb(v: Optional[float]) -> str:
    """Bases -> gigabases, 1 decimal. Empty if None."""
    if v is None:
        return ""
    return f"{(v / 1e9):.1f}"


def fmt_kb(v: Optional[float]) -> str:
    """Bases -> kilobases, 1 decimal. Empty if None."""
    if v is None:
        return ""
    return f"{(v / 1e3):.1f}"


def fmt_int(v: Optional[float]) -> str:
    """Rounded integer. Empty if None."""
    if v is None:
        return ""
    # Round half up behavior for non-negative typical quality; Python's round is banker's.
    # Use +0.5 then floor via int() for non-negative; if negative, fall back to round().
    if v >= 0:
        return str(int(v + 0.5))
    return str(int(round(v)))


def write_tsv(path: str, header: list[str], rows: list[list[str]]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def write_markdown(path_md: str, header: list[str], rows: list[list[str]]) -> None:
    with open(path_md, "w") as f:
        f.write("| " + " | ".join(header) + " |\n")
        # left-align all columns (MultiMarkdown/pipe)
        f.write("|" + "|".join([":---"] * len(header)) + "|\n")
        for r in rows:
            f.write(
                "| " + " | ".join(str(x) if x is not None else "" for x in r) + " |\n"
            )


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Produce formatted Table S1 TSV (and optional Markdown)."
    )
    ap.add_argument("--manifest", required=True, help="Input POLAP manifest JSON")
    ap.add_argument(
        "--out",
        required=True,
        help="Output TSV path (Markdown will use same stem if --markdown)",
    )
    ap.add_argument(
        "--markdown",
        action="store_true",
        default=False,
        help="Also write Markdown table beside TSV (same stem, .md extension)",
    )
    args = ap.parse_args()

    if not os.path.exists(args.manifest):
        print(f"[ERR] manifest not found: {args.manifest}", file=sys.stderr)
        return 2

    # Load JSON
    try:
        with open(args.manifest, "r") as fh:
            manifest = json.load(fh)
    except Exception as e:
        print(f"[ERR] failed to read JSON: {e}", file=sys.stderr)
        return 2

    items = manifest.get("items", [])
    if not isinstance(items, list) or not items:
        print("[ERR] No 'items' array in manifest or empty.", file=sys.stderr)
        return 2

    header = [
        "Code",
        "Species",
        "SRA",
        "Bases_Gb",
        # "Reads",
        "AvgLen_kb",
        "N50_kb",
        "AvgQ",
    ]
    out_rows: list[list[str]] = []

    # Build rows with conversion/scaling applied
    for it in items:
        code2 = it.get("code2", "")  # Species code
        species_name = it.get("species", "")  # Species name
        d = it.get("data", {}) or {}

        v_bases = num(d.get("total_bases"))
        v_reads = num(d.get("read_count"))
        v_avlen = num(d.get("mean_length"))
        v_n50 = num(d.get("N50"))
        v_q = num(d.get("avg_qual"))
        v_sra = d.get("sra_id")

        row = [
            str(code2) if code2 is not None else "",
            str(species_name) if species_name is not None else "",
            str(v_sra) if v_sra is not None else "",
            fmt_gb(v_bases),  # Bases (GB, 1 dec)
            # str(int(v_reads)) if v_reads is not None else "",  # Reads (integer display)
            fmt_kb(v_avlen),  # AvgLen (kb, 1 dec)
            fmt_kb(v_n50),  # N50 (kb, 1 dec)
            fmt_int(v_q),  # AvgQ (rounded int)
        ]
        out_rows.append(row)

    # Write TSV
    write_tsv(args.out, header, out_rows)
    print(f"[INFO] Wrote TSV: {args.out}")

    # Optionally write Markdown
    if args.markdown:
        md_out = os.path.splitext(args.out)[0] + ".md"
        write_markdown(md_out, header, out_rows)
        print(f"[INFO] Wrote Markdown: {md_out}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
