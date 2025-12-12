#!/usr/bin/env python3
###############################################################################
# scripts/polap-py-oatk-graph-table.py
#
# Version : v0.6.1  (2025-12-12)
# Author  : Sang Chul Choi (POLAP)
# License : GPL-3.0+
#
# Purpose
# -------
# From a POLAP OATK manifest JSON, build a long-format CSV of metrics
# across tools and patterns, plus MT/PT GFA/PNG paths.
#
# IMPORTANT CHANGE (v0.6.1)
# ------------------------
# Flat OATK now uses coverage-scoped MT/PT keys:
#   c10_mt_gfa, c10_mt_png, c10_mt_fasta
#   c10_pt_gfa, c10_pt_png, c10_pt_fasta
# (and c20, c30, ...)
###############################################################################

import argparse
import csv
import json
import os
import re
import sys
from typing import Any, Dict, Iterable


def _get(d: Dict[str, Any], key: str, default: str = "") -> str:
    if not isinstance(d, dict):
        return default
    v = d.get(key, default)
    return "" if v is None else str(v)


# ----------------------------------------------------------------------
# OATK (flat: pattern = "x", cov = 10/20/30)
# ----------------------------------------------------------------------
def collect_oatk_flat_rows(
    item: Dict[str, Any], base_info: Dict[str, str]
) -> Iterable[Dict[str, str]]:

    oatk = item.get("oatk") or {}
    if not isinstance(oatk, dict):
        return

    # detect available coverages from cXX_mem_gb
    covs = []
    for k in oatk:
        m = re.match(r"^c(\d+)_mem_gb$", k)
        if m:
            covs.append(m.group(1))

    for cov in sorted(covs, key=int):
        p = f"c{cov}"

        row = {
            **base_info,
            "tool": "oatk",
            "pattern": "x",
            "cov": cov,
            # resources
            "mem_gb": _get(oatk, f"{p}_mem_gb"),
            "time_hours": _get(oatk, f"{p}_time_hours"),
            "disk_gb": _get(oatk, f"{p}_disk_gb"),
            "total_mem_gb": _get(oatk, f"{p}_total_mem_gb"),
            "total_time_hours": _get(oatk, f"{p}_total_time_hours"),
            "total_disk_gb": _get(oatk, f"{p}_total_disk_gb"),
            # assembly quality
            "completeness_prop": _get(oatk, f"{p}_geneset_completeness_prop"),
            "num_contigs": _get(oatk, f"{p}_num_contigs"),
            "total_length": _get(oatk, f"{p}_total_length"),
            "N50": _get(oatk, f"{p}_N50"),
            # MT / PT paths (NEW SCHEMA)
            "mt_gfa": _get(oatk, f"{p}_mt_gfa"),
            "mt_png": _get(oatk, f"{p}_mt_png"),
            "pt_gfa": _get(oatk, f"{p}_pt_gfa"),
            "pt_png": _get(oatk, f"{p}_pt_png"),
        }

        yield row


# ----------------------------------------------------------------------
# OATK (patterned: oatk-tiara-himt-*)
# (unchanged; keys already explicit)
# ----------------------------------------------------------------------
def collect_oatk_pattern_rows(
    item: Dict[str, Any], base_info: Dict[str, str]
) -> Iterable[Dict[str, str]]:

    patterns = ["h-t", "t-h", "t-x", "x-h", "x-x"]

    for pattern in patterns:
        org = f"oatk-tiara-himt-{pattern}"
        block = item.get(org)
        if not isinstance(block, dict):
            continue

        for ckey, cval in block.items():
            if not ckey.startswith("c"):
                continue
            cov = ckey[1:]

            row = {
                **base_info,
                "tool": "oatk",
                "pattern": pattern,
                "cov": cov,
                "mem_gb": _get(cval, "mem_gb"),
                "time_hours": _get(cval, "time_hours"),
                "disk_gb": _get(cval, "disk_gb"),
                "total_mem_gb": _get(cval, "total_mem_gb"),
                "total_time_hours": _get(cval, "total_time_hours"),
                "total_disk_gb": _get(cval, "total_disk_gb"),
                "completeness_prop": _get(cval, "geneset_completeness_prop"),
                "num_contigs": _get(cval, "num_contigs"),
                "total_length": _get(cval, "total_length"),
                "N50": _get(cval, "N50"),
                "mt_gfa": _get(cval, "mito_gfa"),
                "mt_png": _get(cval, "mito_png"),
                "pt_gfa": _get(cval, "plastid_gfa"),
                "pt_png": _get(cval, "plastid_png"),
            }
            yield row


# ----------------------------------------------------------------------
# HiMT / TIPPo (unchanged)
# ----------------------------------------------------------------------
def collect_single_tool_rows(
    item: Dict[str, Any], base_info: Dict[str, str], tool: str
) -> Iterable[Dict[str, str]]:

    blk = item.get(tool) or {}
    if not isinstance(blk, dict):
        return
    if all(k == "present" for k in blk):
        return

    yield {
        **base_info,
        "tool": tool,
        "pattern": "",
        "cov": "",
        "mem_gb": _get(blk, "mem_gb"),
        "time_hours": _get(blk, "time_hours"),
        "disk_gb": _get(blk, "disk_gb"),
        "total_mem_gb": _get(blk, "total_mem_gb"),
        "total_time_hours": _get(blk, "total_time_hours"),
        "total_disk_gb": _get(blk, "total_disk_gb"),
        "completeness_prop": _get(blk, "geneset_completeness_prop"),
        "num_contigs": _get(blk, "num_contigs"),
        "total_length": _get(blk, "total_length"),
        "N50": _get(blk, "N50"),
        "mt_gfa": _get(blk, "mt_gfa"),
        "mt_png": _get(blk, "mt_png"),
        "pt_gfa": _get(blk, "pt_gfa"),
        "pt_png": _get(blk, "pt_png"),
    }


# ----------------------------------------------------------------------
def iter_rows_for_item(item: Dict[str, Any]) -> Iterable[Dict[str, str]]:
    data = item.get("data") or {}
    base = {
        "species": _get(item, "species"),
        "code2": _get(item, "code2"),
        "avg_qual": _get(data, "avg_qual"),
        "mean_length": _get(data, "mean_length"),
    }

    yield from collect_oatk_flat_rows(item, base)
    yield from collect_oatk_pattern_rows(item, base)
    yield from collect_single_tool_rows(item, base, "himt")
    yield from collect_single_tool_rows(item, base, "tippo")


# ----------------------------------------------------------------------
def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--out-csv", required=True)
    args = ap.parse_args(argv)

    with open(args.manifest, encoding="utf-8") as fh:
        manifest = json.load(fh)

    items = manifest.get("items", [])
    if not isinstance(items, list):
        raise SystemExit("manifest['items'] must be a list")

    os.makedirs(os.path.dirname(os.path.abspath(args.out_csv)), exist_ok=True)

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
        "disk_gb",
        "total_mem_gb",
        "total_time_hours",
        "total_disk_gb",
        "completeness_prop",
        "num_contigs",
        "total_length",
        "N50",
        "mt_gfa",
        "mt_png",
        "pt_gfa",
        "pt_png",
    ]

    with open(args.out_csv, "w", newline="", encoding="utf-8") as out:
        w = csv.DictWriter(out, fieldnames=fieldnames)
        w.writeheader()
        for item in items:
            for row in iter_rows_for_item(item):
                w.writerow(row)

    return 0


if __name__ == "__main__":
    sys.exit(main())
