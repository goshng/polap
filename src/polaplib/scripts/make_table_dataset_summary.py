#!/usr/bin/env python3
# -*- coding: utf-8 -*-
###############################################################################
# scripts/make_table_dataset_summary.py
#
# Version : v0.2.3
# Author  : Sang Chul Choi (POLAP)
# Date    : 2025-12-15
# Date    : 2025-12-05 v0.2.2
# Date    : 2025-12-04 v0.2.1
# Date    : 2025-10-22 v0.2.0
# License : GPL-3.0+
#
# Purpose :
#   Read the POLAP manifest JSON and export Table S1 (dataset summary)
#   as TSV (and optionally Markdown) for manuscripts.
#
# New in v0.2.3:
#   - Add --columns to select a subset of columns (and order).
#
# Available output columns (use these names in --columns):
#   Code, Species, SRA, Bases (Gb), Reads, AvgLen (kb), N50 (kb), AvgQ
#
# Example:
#   python scripts/make_table_dataset_summary.py \
#       --manifest md/manifest.json \
#       --out md/tableS1.tsv \
#       --columns Code,Species,SRA \
#       --markdown
###############################################################################
import os
import sys
import json
import csv
import argparse
from typing import Any, Optional, Dict, List


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
    if v is None:
        return ""
    return f"{(v / 1e9):.1f}"


def fmt_kb(v: Optional[float]) -> str:
    if v is None:
        return ""
    return f"{(v / 1e3):.1f}"


def fmt_int(v: Optional[float]) -> str:
    if v is None:
        return ""
    if v >= 0:
        return str(int(v + 0.5))
    return str(int(round(v)))


def fmt_species(name: str) -> str:
    """Remove underscores and wrap in Markdown italics."""
    if not name:
        return ""
    clean = name.replace("_", " ").strip()
    return f"*{clean}*"


def write_tsv(path: str, header: List[str], rows: List[List[str]]) -> None:
    outdir = os.path.dirname(path)
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def write_markdown(path_md: str, header: List[str], rows: List[List[str]]) -> None:
    with open(path_md, "w") as f:
        f.write("| " + " | ".join(header) + " |\n")
        f.write("|" + "|".join([":---"] * len(header)) + "|\n")
        for r in rows:
            f.write("| " + " | ".join(r) + " |\n")


def parse_columns_arg(s: Optional[str], available: List[str]) -> List[str]:
    """Parse --columns 'A,B,C' into a validated list."""
    if not s:
        return []
    cols = [c.strip() for c in s.split(",") if c.strip()]
    bad = [c for c in cols if c not in available]
    if bad:
        raise ValueError(
            "Unknown column(s) in --columns: "
            + ", ".join(bad)
            + "\nAvailable: "
            + ", ".join(available)
        )
    return cols


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Produce formatted Table S1 TSV (and optional Markdown)."
    )
    ap.add_argument("--manifest", required=True, help="Input POLAP manifest JSON")
    ap.add_argument(
        "--out",
        required=True,
        help="Output TSV path (Markdown uses same stem if --markdown)",
    )
    ap.add_argument(
        "--markdown",
        action="store_true",
        default=False,
        help="Also write Markdown table beside TSV (same stem, .md extension)",
    )
    ap.add_argument(
        "--columns",
        default="",
        help=(
            "Comma-separated output columns (subset + order). "
            "Example: Code,Species,SRA. "
            "Available: Code,Species,SRA,Bases (Gb),Reads,AvgLen (kb),N50 (kb),AvgQ"
        ),
    )
    args = ap.parse_args()

    if not os.path.exists(args.manifest):
        print(f"[ERR] manifest not found: {args.manifest}", file=sys.stderr)
        return 2

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

    # Define ALL available output columns.
    AVAILABLE = [
        "Code",
        "Species",
        "SRA",
        "Bases (Gb)",
        "Reads",
        "AvgLen (kb)",
        "N50 (kb)",
        "AvgQ",
    ]

    # Your current default output (v0.2.2 style)
    DEFAULT_HEADER = ["Code", "Species", "SRA", "Bases (Gb)", "AvgLen (kb)", "AvgQ"]

    try:
        selected = parse_columns_arg(args.columns, AVAILABLE)
    except ValueError as e:
        print(f"[ERR] {e}", file=sys.stderr)
        return 2

    header = selected if selected else DEFAULT_HEADER

    out_rows: List[List[str]] = []
    for it in items:
        code2 = it.get("code2", "")
        species_name = it.get("species", "")
        d = it.get("data", {}) or {}

        v_bases = num(d.get("total_bases"))
        v_reads = num(d.get("read_count"))
        v_avlen = num(d.get("mean_length"))
        v_n50 = num(d.get("N50"))
        v_q = num(d.get("avg_qual"))
        v_sra = d.get("sra_id")

        # Build a full row dict, then emit only selected columns in order.
        row_map: Dict[str, str] = {
            "Code": str(code2) if code2 is not None else "",
            "Species": fmt_species(species_name) if species_name is not None else "",
            "SRA": str(v_sra) if v_sra is not None else "",
            "Bases (Gb)": fmt_gb(v_bases),
            "Reads": str(int(v_reads)) if v_reads is not None else "",
            "AvgLen (kb)": fmt_kb(v_avlen),
            "N50 (kb)": fmt_kb(v_n50),
            "AvgQ": fmt_int(v_q),
        }

        out_rows.append([row_map.get(col, "") for col in header])

    write_tsv(args.out, header, out_rows)
    print(f"[INFO] Wrote TSV: {args.out}")

    if args.markdown:
        md_out = os.path.splitext(args.out)[0] + ".md"
        write_markdown(md_out, header, out_rows)
        print(f"[INFO] Wrote Markdown: {md_out}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
