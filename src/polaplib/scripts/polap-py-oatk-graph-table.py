#!/usr/bin/env python3
###############################################################################
# scripts/polap-py-oatk-graph-table.py
#
# Version : v0.4.0  (2025-12-08)
# Author  : Sang Chul Choi (POLAP)
# License : GPL-3.0+
#
# Purpose
# -------
# From a POLAP OATK manifest JSON, build a long-format CSV of metrics
# AND assembly-graph PNG paths for each pattern/coverage combination.
#
# This version v0.4.0 is *compatible* with v0.3.0:
#   - All previous columns are preserved with the same semantics:
#       species, code2, pattern, cov,
#       geneset_completeness_prop, num_contigs, total_length, N50,
#       oatk_mem_gb, oatk_time_hours,
#       himt_mem_gb, himt_time_hours,
#       tiara_mem_gb, tiara_time_hours
#   - NEW columns appended at the end:
#       mito_png, plastid_png
#
# Output CSV columns (in order):
#   species
#   code2
#   pattern
#   cov
#   geneset_completeness_prop
#   num_contigs
#   total_length
#   N50
#   oatk_mem_gb
#   oatk_time_hours
#   himt_mem_gb
#   himt_time_hours
#   tiara_mem_gb
#   tiara_time_hours
#   mito_png
#   plastid_png
#
# JSON layout assumed:
#   • Flat OATK metrics in item["oatk"] as keys:
#       "<tag>_geneset_completeness_prop",
#       "<tag>_num_contigs",
#       "<tag>_total_length",
#       "<tag>_N50"
#     where <tag> looks like "ch-t-02", "ct-h-03", "cx-x-20", etc.
#
#   • Nested pattern blocks:
#       item["oatk-tiara-himt-<pattern>"]["cXX"]["mem_gb"/"time_hours"/"mito_png"/"plastid_png"]
#       item["oatk-tiara-himt-<pattern>"]["himt"]["filter_mem_gb"/"filter_time_hours"]
#       item["oatk-tiara-himt-<pattern>"]["tiara"]["mem_gb"/"time_hours"]
#
#   where:
#       pattern ∈ { "h-t", "t-h", "t-x", "x-h", "x-x" } when present
#       cXX uses zero-padded coverage (e.g. "c02", "c20", "c30")
###############################################################################

import argparse
import csv
import json
import os
import sys
from typing import Any, Dict, List, Tuple


# 4C-ish metrics from the flat "oatk" block
OATK_METRIC_NAMES = (
    "geneset_completeness_prop",
    "num_contigs",
    "total_length",
    "N50",
)


def _as_str(v: Any) -> str:
    """Return v as string, or empty string if None."""
    if v is None:
        return ""
    return str(v)


def collect_oatk_tags(oatk: Dict[str, Any]) -> Dict[str, Dict[str, str]]:
    """
    From the 'oatk' dict, collect metrics by tag.

    Expects keys of the form:
        <tag>_<metric>
    where <metric> is in OATK_METRIC_NAMES.

    Returns:
        { tag: { metric_name: value_str, ... }, ... }
    """
    result: Dict[str, Dict[str, str]] = {}
    if not isinstance(oatk, dict):
        return result

    for key, val in oatk.items():
        if not isinstance(key, str):
            continue
        for metric in OATK_METRIC_NAMES:
            suffix = "_" + metric
            if key.endswith(suffix):
                tag = key[: -len(suffix)]
                if not tag:
                    continue
                tag_block = result.setdefault(tag, {})
                tag_block[metric] = _as_str(val)
                break

    return result


def split_tag(tag: str) -> Tuple[str, str]:
    """
    Split a tag like "ch-t-02" into (pattern, cov):

      "ch-t-02" -> ("h-t", "2")
      "ct-h-02" -> ("t-h", "2")
      "ct-x-02" -> ("t-x", "2")
      "cx-h-02" -> ("x-h", "2")
      "cx-x-20" -> ("x-x", "20")

    If the format is unexpected, returns ("", "").
    """
    if not tag:
        return "", ""

    parts = tag.split("-")
    if len(parts) < 2:
        return "", ""

    cov_str = parts[-1]
    core = "-".join(parts[:-1])  # e.g. "ch-t"

    # remove leading "c" from the core pattern if present
    if core.startswith("c") and len(core) > 1:
        pattern = core[1:]
    else:
        pattern = core

    # normalize cov to integer string (no leading zeros when numeric)
    cov = cov_str
    try:
        cov = str(int(cov_str))
    except ValueError:
        cov = cov_str

    return pattern, cov


def get_pattern_blocks(
    item: Dict[str, Any], pattern: str, cov: str
) -> Tuple[Dict[str, Any], Dict[str, Any], Dict[str, Any], Dict[str, Any]]:
    """
    For a given pattern and coverage, return:
      pattern_block, cov_block, himt_block, tiara_block

    Where:
      pattern_block = item["oatk-tiara-himt-<pattern>"] (or {})
      cov_block     = pattern_block["cXX"] (or {}),
                      where XX is zero-padded cov
      himt_block    = pattern_block["himt"] (or {})
      tiara_block   = pattern_block["tiara"] (or {})
    """
    org = f"oatk-tiara-himt-{pattern}"
    pattern_block = item.get(org) or {}

    try:
        cov_int = int(cov)
        cov_key = f"c{cov_int:02d}"
    except ValueError:
        cov_key = f"c{cov}"

    cov_block = pattern_block.get(cov_key) or {}
    himt_block = pattern_block.get("himt") or {}
    tiara_block = pattern_block.get("tiara") or {}

    return pattern_block, cov_block, himt_block, tiara_block


def iter_rows(item: Dict[str, Any]):
    """
    Yield rows for a single manifest item (species):

      species, code2, pattern, cov,
      geneset_completeness_prop, num_contigs, total_length, N50,
      oatk_mem_gb, oatk_time_hours,
      himt_mem_gb, himt_time_hours,
      tiara_mem_gb, tiara_time_hours,
      mito_png, plastid_png
    """
    species = _as_str(item.get("species", ""))
    code2 = _as_str(item.get("code2", ""))

    oatk_flat = item.get("oatk") or {}
    tags = collect_oatk_tags(oatk_flat)

    for tag, metrics in sorted(tags.items()):
        pattern, cov = split_tag(tag)
        if not pattern or not cov:
            continue

        _, cov_block, himt_block, tiara_block = get_pattern_blocks(item, pattern, cov)

        # OATK mem/time from cov_block
        oatk_mem_gb = _as_str(cov_block.get("mem_gb", ""))
        oatk_time_hours = _as_str(cov_block.get("time_hours", ""))

        # HiMT mem/time (pattern-wide)
        himt_mem_gb = _as_str(himt_block.get("filter_mem_gb", ""))
        himt_time_hours = _as_str(himt_block.get("filter_time_hours", ""))

        # Tiara mem/time (pattern-wide)
        tiara_mem_gb = _as_str(tiara_block.get("mem_gb", ""))
        tiara_time_hours = _as_str(tiara_block.get("time_hours", ""))

        # PNG paths from cov_block
        mito_png = _as_str(cov_block.get("mito_png", ""))
        plastid_png = _as_str(cov_block.get("plastid_png", ""))

        yield {
            "species": species,
            "code2": code2,
            "pattern": pattern,
            "cov": cov,
            "geneset_completeness_prop": metrics.get("geneset_completeness_prop", ""),
            "num_contigs": metrics.get("num_contigs", ""),
            "total_length": metrics.get("total_length", ""),
            "N50": metrics.get("N50", ""),
            "oatk_mem_gb": oatk_mem_gb,
            "oatk_time_hours": oatk_time_hours,
            "himt_mem_gb": himt_mem_gb,
            "himt_time_hours": himt_time_hours,
            "tiara_mem_gb": tiara_mem_gb,
            "tiara_time_hours": tiara_time_hours,
            "mito_png": mito_png,
            "plastid_png": plastid_png,
        }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=(
            "Build long-format metrics CSV for OATK patterns from a POLAP "
            "manifest JSON. v0.4.0 (compatible with v0.3.0) – adds mito_png "
            "and plastid_png columns at the end."
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
        help="Output CSV path (e.g., oatk-graph-metrics.csv)",
    )
    args = ap.parse_args(argv)

    with open(args.manifest, "r", encoding="utf-8") as fh:
        manifest = json.load(fh)

    items = manifest.get("items", [])
    if not isinstance(items, list):
        raise SystemExit("manifest['items'] must be a list")

    os.makedirs(os.path.dirname(os.path.abspath(args.out_csv)), exist_ok=True)

    fieldnames = [
        "species",
        "code2",
        "pattern",
        "cov",
        "geneset_completeness_prop",
        "num_contigs",
        "total_length",
        "N50",
        "oatk_mem_gb",
        "oatk_time_hours",
        "himt_mem_gb",
        "himt_time_hours",
        "tiara_mem_gb",
        "tiara_time_hours",
        "mito_png",
        "plastid_png",
    ]

    with open(args.out_csv, "w", newline="", encoding="utf-8") as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=fieldnames)
        writer.writeheader()
        for item in items:
            for row in iter_rows(item):
                writer.writerow(row)

    return 0


if __name__ == "__main__":
    sys.exit(main())

