#!/usr/bin/env python
# scripts/polap-py-hifi-graph-table.py
# Version : v0.1.0  (2025-12-02)
# Author  : Sang Chul Choi (POLAP)
# License : GPL-3.0+
#
# Read v5-0-auto-manifest.json and emit a long-format CSV:
#   species, pipeline,
#   mem_gb, time_hours, disk_gb,
#   geneset_completeness_prop,
#   num_contigs, N50, fragmentation_index,
#   contig_length_cv, max_contig_prop, total_length
#
# Pipelines:
#   pmat  -> JSON key "pmat2"
#   tippo -> JSON key "tippo"
#   himt  -> JSON key "himt"
#   oatk  -> JSON key "oatk" (c30_* metrics preferred)
#
import argparse
import csv
import json
import sys
from typing import Any, Dict, Optional


def _get(d: Optional[Dict[str, Any]], key: str) -> str:
    """Return d[key] as string, or empty string if missing/None."""
    if not isinstance(d, dict):
        return ""
    v = d.get(key, "")
    if v is None:
        return ""
    return str(v)


def _oatk_metric(oatk: Optional[Dict[str, Any]], suffix: str) -> str:
    """
    Get Oatk metric with c30_* preferred.
    For example suffix='mem_gb' -> c30_mem_gb.
    """
    if not isinstance(oatk, dict):
        return ""
    # c30_* only for now (you can extend to fallback to c20/c10 later)
    key = f"c30_{suffix}"
    v = oatk.get(key, "")
    if v is None:
        return ""
    return str(v)


def build_rows(item: Dict[str, Any]):
    """Yield rows for pmat, tippo, himt, oatk for a single manifest item."""
    species = str(item.get("species", ""))

    pmat2 = item.get("pmat2") or {}
    tippo = item.get("tippo") or {}
    himt = item.get("himt") or {}
    oatk = item.get("oatk") or {}

    # helper to build a row dict
    def row(
        pipeline: str,
        mem_gb: str,
        time_hours: str,
        disk_gb: str,
        geneset_completeness_prop: str,
        num_contigs: str,
        N50: str,
        fragmentation_index: str,
        contig_length_cv: str,
        max_contig_prop: str,
        total_length: str,
    ):
        return {
            "species": species,
            "pipeline": pipeline,
            "mem_gb": mem_gb,
            "time_hours": time_hours,
            "disk_gb": disk_gb,
            "geneset_completeness_prop": geneset_completeness_prop,
            "num_contigs": num_contigs,
            "N50": N50,
            "fragmentation_index": fragmentation_index,
            "contig_length_cv": contig_length_cv,
            "max_contig_prop": max_contig_prop,
            "total_length": total_length,
        }

    # PMAT (pmat2)
    yield row(
        "pmat",
        _get(pmat2, "mem_gb"),
        _get(pmat2, "time_hours"),
        _get(pmat2, "disk_gb"),
        _get(pmat2, "geneset_completeness_prop"),
        _get(pmat2, "num_contigs"),
        _get(pmat2, "N50"),
        _get(pmat2, "fragmentation_index"),
        _get(pmat2, "contig_length_cv"),
        _get(pmat2, "max_contig_prop"),
        _get(pmat2, "total_length"),
    )

    # TIPPo
    yield row(
        "tippo",
        _get(tippo, "mem_gb"),
        _get(tippo, "time_hours"),
        _get(tippo, "disk_gb"),
        _get(tippo, "geneset_completeness_prop"),
        _get(tippo, "num_contigs"),
        _get(tippo, "N50"),
        _get(tippo, "fragmentation_index"),
        _get(tippo, "contig_length_cv"),
        _get(tippo, "max_contig_prop"),
        _get(tippo, "total_length"),
    )

    # HiMT
    yield row(
        "himt",
        _get(himt, "mem_gb"),
        _get(himt, "time_hours"),
        _get(himt, "disk_gb"),
        _get(himt, "geneset_completeness_prop"),
        _get(himt, "num_contigs"),
        _get(himt, "N50"),
        _get(himt, "fragmentation_index"),
        _get(himt, "contig_length_cv"),
        _get(himt, "max_contig_prop"),
        _get(himt, "total_length"),
    )

    # Oatk (c30_* metrics)
    yield row(
        "oatk",
        _oatk_metric(oatk, "mem_gb"),
        _oatk_metric(oatk, "time_hours"),
        _oatk_metric(oatk, "disk_gb"),
        _oatk_metric(oatk, "geneset_completeness_prop"),
        _oatk_metric(oatk, "num_contigs"),
        _oatk_metric(oatk, "N50"),
        _oatk_metric(oatk, "fragmentation_index"),
        _oatk_metric(oatk, "contig_length_cv"),
        _oatk_metric(oatk, "max_contig_prop"),
        _oatk_metric(oatk, "total_length"),
    )


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Build long-format metrics CSV for HiFi mtDNA graphs."
    )
    parser.add_argument(
        "--manifest",
        required=True,
        help="v5-0-auto-manifest.json",
    )
    parser.add_argument(
        "--out-csv",
        required=True,
        help="Output CSV path for metrics table",
    )
    args = parser.parse_args(argv)

    with open(args.manifest, "r", encoding="utf-8") as fh:
        manifest = json.load(fh)

    items = manifest.get("items", [])
    if not isinstance(items, list):
        raise SystemExit("manifest['items'] must be a list")

    fieldnames = [
        "species",
        "pipeline",
        "mem_gb",
        "time_hours",
        "disk_gb",
        "geneset_completeness_prop",
        "num_contigs",
        "N50",
        "fragmentation_index",
        "contig_length_cv",
        "max_contig_prop",
        "total_length",
    ]

    with open(args.out_csv, "w", newline="", encoding="utf-8") as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=fieldnames)
        writer.writeheader()
        for item in items:
            for r in build_rows(item):
                writer.writerow(r)


if __name__ == "__main__":
    main()
