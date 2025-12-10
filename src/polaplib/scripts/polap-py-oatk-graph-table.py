#!/usr/bin/env python3
###############################################################################
# scripts/polap-py-oatk-graph-table.py
#
# Version : v0.6.0  (2025-12-09)
# Author  : Sang Chul Choi (POLAP)
# License : GPL-3.0+
#
# Purpose
# -------
# From a POLAP OATK manifest JSON, build a long-format CSV of metrics
# across tools and patterns, plus MT/PT GFA/PNG paths.
#
# Output CSV columns:
#   species
#   code2
#   avg_qual
#   mean_length
#   tool                # "oatk", "himt", "tippo"
#   pattern             # e.g. "h-t","t-h","t-x","x-h","x-x","x" (flat oatk), or "" for others
#   cov                 # coverage number for oatk rows, "" for himt/tippo
#   mem_gb
#   time_hours
#   total_mem_gb
#   total_time_hours
#   completeness_prop   # geneset_completeness_prop
#   num_contigs
#   total_length
#   N50
#   mt_gfa
#   mt_png
#   pt_gfa
#   pt_png
#
# Tools and sources:
#   • OATK (flat):
#       item["oatk"]["c10_mem_gb"], "c10_time_hours", "c10_gfa" (plastid),
#       similarly c20_*, c30_*, ...
#       pattern = "x", cov = 10,20,30,...
#
#   • OATK (patterned):
#       item["oatk-tiara-himt-<pattern>"]["cXX"] fields:
#         total_mem_gb, total_time_hours, mem_gb, time_hours,
#         geneset_completeness_prop, num_contigs, total_length, N50,
#         mito_gfa, mito_png, plastid_gfa, plastid_png
#       where pattern ∈ { "h-t","t-h","t-x","x-h","x-x" } and XX is coverage.
#
#   • HiMT (global tool):
#       item["himt"]["mt_gfa"], "mt_png", "pt_gfa", "pt_png",
#       "mem_gb","time_hours",
#       "geneset_completeness_prop","num_contigs","total_length","N50".
#
#   • TIPPo (global tool):
#       item["tippo"]["mt_gfa"], "mt_png", "pt_gfa", "pt_png",
#       "mem_gb","time_hours",
#       "geneset_completeness_prop","num_contigs","total_length","N50".
###############################################################################

import argparse
import csv
import json
import os
import re
import sys
from typing import Any, Dict, Iterable


def _as_str(v: Any) -> str:
    """Return v as string, or empty string if None."""
    if v is None:
        return ""
    return str(v)


def _get(d: Dict[str, Any], key: str, default: str = "") -> str:
    """Safely get d[key] as string, or default if missing/None."""
    if not isinstance(d, dict):
        return default
    v = d.get(key, default)
    return default if v is None else str(v)


def _guess_png_from_gfa(path: str) -> str:
    """If path ends with .gfa, guess a .png alongside it; else return empty."""
    if not path:
        return ""
    if path.endswith(".gfa"):
        return path[:-4] + ".png"
    return ""


def collect_oatk_flat_rows(item: Dict[str, Any], base_info: Dict[str, str]) -> Iterable[Dict[str, str]]:
    """
    Collect rows for the legacy flat OATK block:
      item["oatk"]["c10_mem_gb"], "c10_time_hours", "c10_gfa", etc.

    pattern is "x"; cov is coverage number after 'c' in keys (10,20,30...).
    """
    oatk = item.get("oatk") or {}
    if not isinstance(oatk, dict):
        return

    cov_set = set()
    for key in oatk.keys():
        m = re.match(r"^c(\d+)_mem_gb$", key)
        if m:
            cov_set.add(m.group(1))

    for cov_str in sorted(cov_set, key=lambda s: int(s)):
        prefix = f"c{cov_str}"

        mem_gb = _get(oatk, f"{prefix}_mem_gb")
        time_hours = _get(oatk, f"{prefix}_time_hours")

        # total_mem_gb / total_time_hours rarely present for flat OATK
        total_mem_gb = _get(oatk, f"{prefix}_total_mem_gb")
        total_time_hours = _get(oatk, f"{prefix}_total_time_hours")

        # geneset metrics: per-cov if they exist, else fall back to global.
        completeness = _get(
            oatk,
            f"{prefix}_geneset_completeness_prop",
            _get(oatk, "geneset_completeness_prop"),
        )
        num_contigs = _get(
            oatk,
            f"{prefix}_num_contigs",
            _get(oatk, "num_contigs"),
        )
        total_length = _get(
            oatk,
            f"{prefix}_total_length",
            _get(oatk, "total_length"),
        )
        N50 = _get(
            oatk,
            f"{prefix}_N50",
            _get(oatk, "N50"),
        )

        # GFA/PNG: flat OATK block generally has plastid-only gfa ("c10_gfa").
        pt_gfa = _get(oatk, f"{prefix}_gfa")
        pt_png = _get(oatk, f"{prefix}_png", _guess_png_from_gfa(pt_gfa))
        mt_gfa = ""
        mt_png = ""

        row = {
            **base_info,
            "tool": "oatk",
            "pattern": "x",
            "cov": cov_str,
            "mem_gb": mem_gb,
            "time_hours": time_hours,
            "total_mem_gb": total_mem_gb,
            "total_time_hours": total_time_hours,
            "completeness_prop": completeness,
            "num_contigs": num_contigs,
            "total_length": total_length,
            "N50": N50,
            "mt_gfa": mt_gfa,
            "mt_png": mt_png,
            "pt_gfa": pt_gfa,
            "pt_png": pt_png,
        }
        yield row


def collect_oatk_pattern_rows(item: Dict[str, Any], base_info: Dict[str, str]) -> Iterable[Dict[str, str]]:
    """
    Collect rows for patterned OATK runs stored in:
      item["oatk-tiara-himt-<pattern>"]["cXX"]

    pattern in { "h-t", "t-h", "t-x", "x-h", "x-x" }.
    """
    patterns = ["h-t", "t-h", "t-x", "x-h", "x-x"]

    for pattern in patterns:
        org_key = f"oatk-tiara-himt-{pattern}"
        block = item.get(org_key)
        if not isinstance(block, dict):
            continue

        for c_key, c_val in block.items():
            if not isinstance(c_key, str):
                continue
            if not c_key.startswith("c"):
                continue
            try:
                cov_int = int(c_key[1:])
                cov_str = str(cov_int)
            except ValueError:
                cov_str = c_key[1:]

            cov_block = c_val if isinstance(c_val, dict) else {}

            mem_gb = _get(cov_block, "mem_gb")
            time_hours = _get(cov_block, "time_hours")
            total_mem_gb = _get(cov_block, "total_mem_gb")
            total_time_hours = _get(cov_block, "total_time_hours")

            completeness = _get(cov_block, "geneset_completeness_prop")
            num_contigs = _get(cov_block, "num_contigs")
            total_length = _get(cov_block, "total_length")
            N50 = _get(cov_block, "N50")

            mt_gfa = _get(cov_block, "mito_gfa")
            mt_png = _get(cov_block, "mito_png")
            pt_gfa = _get(cov_block, "plastid_gfa")
            pt_png = _get(cov_block, "plastid_png", _guess_png_from_gfa(pt_gfa))

            row = {
                **base_info,
                "tool": "oatk",
                "pattern": pattern,
                "cov": cov_str,
                "mem_gb": mem_gb,
                "time_hours": time_hours,
                "total_mem_gb": total_mem_gb,
                "total_time_hours": total_time_hours,
                "completeness_prop": completeness,
                "num_contigs": num_contigs,
                "total_length": total_length,
                "N50": N50,
                "mt_gfa": mt_gfa,
                "mt_png": mt_png,
                "pt_gfa": pt_gfa,
                "pt_png": pt_png,
            }
            yield row


def collect_himt_rows(item: Dict[str, Any], base_info: Dict[str, str]) -> Iterable[Dict[str, str]]:
    """Single-row HiMT tool metrics (no coverage)."""
    himt = item.get("himt") or {}
    if not isinstance(himt, dict):
        return

    has_any = any(k not in ("present",) for k in himt.keys())
    if not has_any:
        return

    mem_gb = _get(himt, "mem_gb")
    time_hours = _get(himt, "time_hours")
    total_mem_gb = _get(himt, "total_mem_gb")
    total_time_hours = _get(himt, "total_time_hours")

    completeness = _get(himt, "geneset_completeness_prop")
    num_contigs = _get(himt, "num_contigs")
    total_length = _get(himt, "total_length")
    N50 = _get(himt, "N50")

    mt_gfa = _get(himt, "mt_gfa")
    mt_png = _get(himt, "mt_png")
    pt_gfa = _get(himt, "pt_gfa")
    pt_png = _get(himt, "pt_png")

    row = {
        **base_info,
        "tool": "himt",
        "pattern": "",
        "cov": "",
        "mem_gb": mem_gb,
        "time_hours": time_hours,
        "total_mem_gb": total_mem_gb,
        "total_time_hours": total_time_hours,
        "completeness_prop": completeness,
        "num_contigs": num_contigs,
        "total_length": total_length,
        "N50": N50,
        "mt_gfa": mt_gfa,
        "mt_png": mt_png,
        "pt_gfa": pt_gfa,
        "pt_png": pt_png,
    }
    yield row


def collect_tippo_rows(item: Dict[str, Any], base_info: Dict[str, str]) -> Iterable[Dict[str, str]]:
    """Single-row TIPPo tool metrics (no coverage)."""
    tippo = item.get("tippo") or {}
    if not isinstance(tippo, dict):
        return

    has_any = any(k not in ("present",) for k in tippo.keys())
    if not has_any:
        return

    mem_gb = _get(tippo, "mem_gb")
    time_hours = _get(tippo, "time_hours")
    total_mem_gb = _get(tippo, "total_mem_gb")
    total_time_hours = _get(tippo, "total_time_hours")

    completeness = _get(tippo, "geneset_completeness_prop")
    num_contigs = _get(tippo, "num_contigs")
    total_length = _get(tippo, "total_length")
    N50 = _get(tippo, "N50")

    mt_gfa = _get(tippo, "mt_gfa")
    mt_png = _get(tippo, "mt_png")
    pt_gfa = _get(tippo, "pt_gfa")
    pt_png = _get(tippo, "pt_png")

    row = {
        **base_info,
        "tool": "tippo",
        "pattern": "",
        "cov": "",
        "mem_gb": mem_gb,
        "time_hours": time_hours,
        "total_mem_gb": total_mem_gb,
        "total_time_hours": total_time_hours,
        "completeness_prop": completeness,
        "num_contigs": num_contigs,
        "total_length": total_length,
        "N50": N50,
        "mt_gfa": mt_gfa,
        "mt_png": mt_png,
        "pt_gfa": pt_gfa,
        "pt_png": pt_png,
    }
    yield row


def iter_rows_for_item(item: Dict[str, Any]) -> Iterable[Dict[str, str]]:
    """Yield all rows (oatk flat, oatk patterns, himt, tippo) for one manifest item."""
    species = _get(item, "species")
    code2 = _get(item, "code2")

    data_block = item.get("data") or {}
    avg_qual = _get(data_block, "avg_qual")
    mean_length = _get(data_block, "mean_length")

    base_info = {
        "species": species,
        "code2": code2,
        "avg_qual": avg_qual,
        "mean_length": mean_length,
    }

    # OATK (flat c10,c20,c30 etc.)
    for row in collect_oatk_flat_rows(item, base_info):
        yield row

    # OATK (patterned h-t, t-h, t-x, x-h, x-x)
    for row in collect_oatk_pattern_rows(item, base_info):
        yield row

    # HiMT
    for row in collect_himt_rows(item, base_info):
        yield row

    # TIPPo
    for row in collect_tippo_rows(item, base_info):
        yield row


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=(
            "Build unified metrics CSV (oatk/himt/tippo) from a POLAP manifest JSON.\n"
            "One row per (tool, pattern, cov) combination, with mt/pt GFA/PNG paths."
        )
    )
    ap.add_argument(
        "--manifest",
        required=True,
        help="Input manifest JSON (e.g., v2-0-auto-manifest.json)",
    )
    ap.add_argument(
        "--out-csv",
        required=True,
        help="Output CSV path (e.g., oatk-metrics-table.csv)",
    )
    args = ap.parse_args(argv)

    with open(args.manifest, "r", encoding="utf-8") as fh:
        manifest = json.load(fh)

    items = manifest.get("items", [])
    if not isinstance(items, list):
        raise SystemExit("manifest['items'] must be a list")

    out_dir = os.path.dirname(os.path.abspath(args.out_csv))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    fieldnames = [
        "species",
        "code2",
        "avg_qual",
        "mean_length",
        "tool",
        "pattern",
        "cov",
        "mem_gb",
        "time_hours",
        "total_mem_gb",
        "total_time_hours",
        "completeness_prop",
        "num_contigs",
        "total_length",
        "N50",
        "mt_gfa",
        "mt_png",
        "pt_gfa",
        "pt_png",
    ]

    with open(args.out_csv, "w", newline="", encoding="utf-8") as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=fieldnames)
        writer.writeheader()
        for item in items:
            for row in iter_rows_for_item(item):
                writer.writerow(row)

    return 0


if __name__ == "__main__":
    sys.exit(main())

