#!/usr/bin/env python3
# scripts/manifest_assemble.py
#
# Version : v0.2.4
# Author  : Sang Chul Choi (POLAP)
# Date    : 2025-10-15
# License : GPL-3.0+
#
# Purpose
# -------
# Read a TSV "facts" table (columns: species, tier, inum, organelle, kind, key, value)
# and build a structured JSON manifest.
#
# What’s new (v0.2.4)
# -------------------
# • Properly marshal short-read metric blocks:
#     short1data  -> {total_bases, read_count, mean_length, N50, avg_qual, gc_content, …}
#     short2data  -> { … }
#   These are handled identically to the long-read “data” block.
# • Keep existing behavior for PT/MT/PTPT files and attributes (ncbi_* to .ref, others to .meta).
# • Optional species two-letter codes via --codes (space-delimited “code species”).
# • Attach GFA stats to pt/ptpt/mt (if scripts/gfa_stats_json.py exists).
#
import os
import sys
import csv
import json
import argparse
import subprocess
import time


# ----------------------------- IO helpers -------------------------------------
def read_facts(path: str):
    """Read TSV into list of dicts (expects header)."""
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def load_codes(path: str):
    """
    Read a space-delimited mapping with header 'code species' or comment lines.
    Returns dict: species -> code
    """
    mapping = {}
    if not path or not os.path.isfile(path):
        return mapping
    with open(path) as fh:
        for ln in fh:
            ln = ln.strip()
            if not ln or ln.startswith("#") or ln.lower().startswith("code species"):
                continue
            parts = ln.split()
            if len(parts) < 2:
                continue
            code, species = parts[0], parts[1]
            mapping[species] = code
    return mapping


# ----------------------------- GFA stats --------------------------------------
def call_gfa_stats_json(gfa_path: str):
    """Call scripts/gfa_stats_json.py --gfa <path> and return parsed JSON (or None)."""
    if not gfa_path or not os.path.isfile(gfa_path):
        return None
    base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    helper = os.path.join(base, "scripts", "gfa_stats_json.py")
    if not os.path.isfile(helper):
        return None
    try:
        out = subprocess.check_output(
            [sys.executable, helper, "--gfa", gfa_path],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
        if out and out != "{}":
            return json.loads(out)
    except Exception:
        return None
    return None


def attach_stats(block: dict):
    """Attach .stats to a block if it has a 'gfa' file and helper is available."""
    gfa = block.get("gfa")
    if not gfa:
        return
    st = call_gfa_stats_json(gfa)
    if st:
        block["stats"] = st


# ----------------------------- value coercion ---------------------------------
def maybe_number(s: str):
    """Convert numeric-looking strings to int/float; keep 'NA'/'', etc. as-is/None."""
    if s is None:
        return None
    s = str(s).strip()
    if s in ("", "NA", "na", "NaN", "nan", "None"):
        return None
    try:
        v = float(s)
        return int(v) if v.is_integer() else v
    except Exception:
        return s


# ----------------------------- main -------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Assemble manifest JSON from harvested TSV facts."
    )
    ap.add_argument("--facts", required=True, help="Input TSV facts")
    ap.add_argument("--set", required=True, help="Set label used for harvesting")
    ap.add_argument("--tier", required=True, help="Tier (e.g., v6)")
    ap.add_argument("--inum", required=True, help="Instance number (stringified)")
    ap.add_argument("--out", required=True, help="Output JSON path")
    ap.add_argument("--pretty", action="store_true", default=False, help="Pretty JSON")
    ap.add_argument(
        "--species-codes",
        default=None,
        help="Optional mapping file (space-delimited: 'code species')",
    )
    args = ap.parse_args()

    rows = read_facts(args.facts)
    codes_map = load_codes(args.species_codes)

    # Root and accumulation
    root = dict(
        version="v0.2.4",
        generated_at=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        set=args.set,
        tier=args.tier,
        inum=str(args.inum),
        items=[],
    )
    # species -> item(dict)
    items_by_species = {}

    # Organelles that are metric/attribute blocks (all keys flattened into dict)
    ATTR_BLOCKS = {"data", "short1data", "short2data"}

    for r in rows:
        sp = r.get("species", "")
        tier = r.get("tier", args.tier)
        inum = str(r.get("inum", args.inum))
        org = r.get("organelle", "")
        kind = r.get("kind", "")
        key = r.get("key", "")
        val = r.get("value", "")

        if sp not in items_by_species:
            item = dict(species=sp, tier=tier, inum=inum)
            if sp in codes_map:
                item["code2"] = codes_map[sp]
            items_by_species[sp] = item

        item = items_by_species[sp]
        # Ensure block container exists
        if org and org not in item:
            # For metric blocks, initialize as plain dict; for others, also dict
            item[org] = {}

        blk = item.get(org, {})

        if kind == "file":
            # file paths
            if key in ("gfa", "png"):
                blk[key] = val
            elif key == "qc_scatter":
                blk.setdefault("qc", {})["scatter_pdf"] = val
            elif key == "gene_count_file":
                blk.setdefault("annotation", {})["gene_count_file"] = val
            elif key == "circular_fa":
                blk.setdefault("circular", {}).setdefault("fasta", []).append(val)

        elif kind == "attr":
            v = maybe_number(val)

            # metric blocks (long reads + short mates)
            if org in ATTR_BLOCKS:
                blk[key] = v

            # organelle-specific attributes
            elif org in {"pt", "mt"}:
                # Normalize ncbi_* under .ref ; others under .meta
                if key.startswith("ncbi_"):
                    blk.setdefault("ref", {})[key] = v
                else:
                    blk.setdefault("meta", {})[key] = v

            elif org == "ptpt":
                # keep as plain attributes
                blk[key] = v

        # write-back (in case blk was created)
        if org:
            item[org] = blk

    # Attach stats where applicable
    for sp, item in items_by_species.items():
        for sect in ("pt", "ptpt", "mt"):
            if sect in item and isinstance(item[sect], dict):
                attach_stats(item[sect])

    # Preserve order of appearance in facts
    seen = set()
    for r in rows:
        sp = r.get("species", "")
        if sp and sp not in seen and sp in items_by_species:
            root["items"].append(items_by_species[sp])
            seen.add(sp)

    # Write JSON
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as f:
        if args.pretty:
            json.dump(root, f, indent=2)
        else:
            json.dump(root, f, separators=(",", ":"))

    return 0


if __name__ == "__main__":
    sys.exit(main())
