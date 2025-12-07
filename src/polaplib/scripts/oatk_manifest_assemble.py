#!/usr/bin/env python3
# scripts/oatk_manifest_assemble.py
#
# Version : v0.3.0
# Author  : Sang Chul Choi (POLAP)
# Date    : 2025-12-07
# License : GPL-3.0+
#
# Purpose
# -------
# Read a TSV "facts" table (columns: species, tier, inum, organelle, kind, key, value)
# and build a structured JSON manifest.
#
# Changes (v0.3.0)
# ----------------
# • Keep v0.2.5 behavior for:
#     - data / short1data / short2data  (metric blocks, numeric coercion)
#     - pt / mt / ptpt                  (file + ncbi_* attrs into .ref/.meta)
#     - species two-letter codes via --species-codes
#     - attach GFA stats for pt/ptpt/mt (if scripts/gfa_stats_json.py exists)
#     - assembly blocks: oatk, tippo, himt, pmat2 (attrs kept as strings)
# • NEW: Organelles of the form "oatk-tiara-himt-<pattern>" are treated as
#   assembly blocks but with internal grouping:
#     - Keys "himt_*"  go under block["himt"][subkey]
#     - Keys "tiara_*" go under block["tiara"][subkey]
#     - Keys "cXX_*"   (e.g. c02_mito_gfa) go under block["c02"][subkey]
#   where "subkey" is the part after the prefix and underscore.
#   Other keys (e.g. "present") remain top-level in that organelle block.
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
    """
    Convert numeric-looking strings to int/float.
    Keep 'NA','',etc. as None; the caller may or may not use this.
    """
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


# ----------------------------- org helpers ------------------------------------
def is_assembly_block(org: str) -> bool:
    """
    Return True if this org name should be treated as an assembly-specific
    block (OATK/TIPPo/HiMT/PMAT2 or derived patterns), i.e. attrs retained
    as strings.
    """
    if not org:
        return False
    org = str(org)
    # Any oatk-* style block (oatk-tiara-himt-*, oatk-foo, etc.)
    if org.startswith("oatk-"):
        return True
    if org in ("oatk", "tippo", "himt", "pmat2"):
        return True
    return False


def is_oatk_tiara_himt_block(org: str) -> bool:
    """Return True if organelle is of the form 'oatk-tiara-himt-<pattern>'."""
    if not org:
        return False
    return str(org).startswith("oatk-tiara-himt-")


def assign_oatk_tiara_himt_key(block: dict, key: str, val: str):
    """
    For oatk-tiara-himt-* blocks, place keys into nested dictionaries:

      himt_*        -> block["himt"][<rest>]
      tiara_*       -> block["tiara"][<rest>]
      cXX_* / cNN_* -> block["cXX"][<rest>] (coverage label is the prefix)
      other keys    -> block[<key>]

    Where <rest> is everything after the first underscore.
    """
    if not key:
        return

    # handle himt_* / tiara_*
    for prefix in ("himt_", "tiara_"):
        if key.startswith(prefix):
            subkey = key[len(prefix) :]
            if not subkey:
                subkey = key  # fallback, though unexpected
            block.setdefault(prefix.rstrip("_"), {})[subkey] = val
            return

    # handle coverage blocks c02_*, c20_*, c3_*, etc.
    if "_" in key and key[0] == "c":
        cov_label, rest = key.split("_", 1)
        # cov_label must be 'c' followed by digits
        if cov_label[1:].isdigit():
            block.setdefault(cov_label, {})[rest] = val
            return

    # fallback: just store at the top level of this organelle block
    block[key] = val


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
        version="v0.3.0",
        generated_at=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        set=args.set,
        tier=args.tier,
        inum=str(args.inum),
        items=[],
    )
    # species -> item(dict)
    items_by_species = {}

    # Organelles that are metric/attribute blocks (all keys flattened into dict, numeric)
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

        if org and org not in item:
            item[org] = {}

        blk = item.get(org, {})

        if kind == "file":
            # Assembly file paths: treat as simple key->path for assembly blocks
            if is_assembly_block(org):
                sval = "" if val is None else str(val)
                if is_oatk_tiara_himt_block(org):
                    assign_oatk_tiara_himt_key(blk, key, sval)
                else:
                    blk[key] = sval

            # regular file paths
            elif key in ("gfa", "png"):
                blk[key] = val
            elif key == "qc_scatter":
                blk.setdefault("qc", {})["scatter_pdf"] = val
            elif key == "gene_count_file":
                blk.setdefault("annotation", {})["gene_count_file"] = val
            elif key == "circular_fa":
                blk.setdefault("circular", {}).setdefault("fasta", []).append(val)

        elif kind == "attr":
            # 1) Long/short metric blocks: coerce to number when possible
            if org in ATTR_BLOCKS:
                blk[key] = maybe_number(val)

            # 2) Assembly blocks: keep values as *strings*
            elif is_assembly_block(org):
                sval = "" if val is None else str(val)
                if is_oatk_tiara_himt_block(org):
                    assign_oatk_tiara_himt_key(blk, key, sval)
                else:
                    blk[key] = sval

            # 3) Organellar pt/mt attributes: use ref/meta
            elif org in {"pt", "mt"}:
                v = maybe_number(val)
                if key.startswith("ncbi_"):
                    blk.setdefault("ref", {})[key] = v
                else:
                    blk.setdefault("meta", {})[key] = v

            # 4) ptpt attributes: plain dict (numeric)
            elif org == "ptpt":
                blk[key] = maybe_number(val)

            # 5) Any other org: just store attrs as plain strings
            else:
                blk[key] = "" if val is None else str(val)

        # write-back (in case blk was created/updated)
        if org:
            item[org] = blk

    # Attach stats where applicable (unchanged)
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

