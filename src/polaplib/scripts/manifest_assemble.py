#!/usr/bin/env python3
# scripts/manifest_assemble.py
#
# Version : v0.2.3
# Author  : Sang Chul Choi (POLAP)
# Date    : 2025-10-13
# License : GPL-3.0+
#
# Purpose:
#   Read a TSV facts table produced by Bash and build a structured JSON manifest.
#   Adds:
#     - GFA stats for each organelle via scripts/gfa_stats_json.py (if present)
#     - Dataset metrics from seqkit stats (data.*)
#     - Organellar attributes (pt.ref.ncbi_accession, ncbi_len_bp, etc.)
#     - (NEW) Optional species 2-letter codes via --codes file ("code species")
#
import os, sys, csv, json, argparse, subprocess, time


def read_facts(path: str):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def call_gfa_stats_json(gfa_path: str):
    if not gfa_path or not os.path.isfile(gfa_path):
        return None
    base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    helper = os.path.join(base, "scripts", "gfa_stats_json.py")
    if not os.path.isfile(helper):
        return None
    try:
        out = subprocess.check_output(
            [sys.executable, helper, "--gfa", gfa_path], text=True
        )
        out = out.strip()
        if out and out != "{}":
            return json.loads(out)
    except Exception:
        return None
    return None


def attach_stats(block: dict):
    gfa = block.get("gfa")
    if not gfa:
        return
    st = call_gfa_stats_json(gfa)
    if st:
        block["stats"] = st


def load_codes(path: str):
    """
    Read a space-delimited mapping file with header 'code species'.
    Returns dict: species -> code
    """
    m = {}
    if not path or not os.path.isfile(path):
        return m
    with open(path) as f:
        for ln in f:
            ln = ln.strip()
            if not ln:
                continue
            # skip header and comments
            if ln.lower().startswith("code species") or ln.startswith("#"):
                continue
            parts = ln.split()
            if len(parts) < 2:
                continue
            code, species = parts[0], parts[1]
            m[species] = code
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--facts", required=True)
    ap.add_argument("--set", required=True)
    ap.add_argument("--tier", required=True)
    ap.add_argument("--inum", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--pretty", action="store_true", default=False)
    # NEW: optional codes file (space-delimited: "code species")
    ap.add_argument(
        "--codes",
        default=None,
        help="Optional mapping file (space-delimited: 'code species')",
    )
    a = ap.parse_args()

    codes_map = load_codes(a.codes)
    rows = read_facts(a.facts)
    data = {}
    root = dict(
        version="v0.2.3",
        generated_at=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        set=a.set,
        tier=a.tier,
        inum=str(a.inum),
        items=[],
    )

    for row in rows:
        sp = row["species"]
        org = row["organelle"]
        kind = row["kind"]
        key = row["key"]
        val = row["value"]

        if sp not in data:
            item = dict(species=sp, tier=a.tier, inum=str(a.inum))
            # inject code2 if available
            if sp in codes_map:
                item["code2"] = codes_map[sp]
            data[sp] = item

        if org not in data[sp]:
            data[sp][org] = {}
        blk = data[sp][org]

        if kind == "file":
            if key in ("gfa", "png"):
                blk[key] = val
            elif key == "qc_scatter":
                blk.setdefault("qc", {})["scatter_pdf"] = val
            elif key == "gene_count_file":
                blk.setdefault("annotation", {})["gene_count_file"] = val
            elif key == "circular_fa":
                blk.setdefault("circular", {}).setdefault("fasta", []).append(val)

        elif kind == "attr":
            v = val
            try:
                if "." in val:
                    v = float(val)
                else:
                    v = int(val)
            except Exception:
                pass

            if key.startswith("circular_"):
                blk.setdefault("circular", {})[key.split("_", 1)[1]] = v
            elif key.startswith("gene_"):
                blk.setdefault("annotation", {})[key.split("_", 1)[1]] = v
            elif org == "data":
                # Dataset metrics & tags (sra_id, data_file_bytes, total_bases, etc.)
                blk[key] = v
            elif org in ("pt", "mt"):
                # Normalize organellar attributes: ncbi_* -> .ref, others -> .meta
                if key.startswith("ncbi_"):
                    blk.setdefault("ref", {})[key] = v
                else:
                    blk.setdefault("meta", {})[key] = v

    # attach GFA stats when gfa present
    for sp, item in data.items():
        for org in ("pt", "ptpt", "mt"):
            if org in item:
                attach_stats(item[org])

    # items[] in order of appearance in facts
    seen = set()
    for row in rows:
        sp = row["species"]
        if sp not in seen:
            root["items"].append(data[sp])
            seen.add(sp)

    with open(a.out, "w") as f:
        if a.pretty:
            json.dump(root, f, indent=2)
        else:
            json.dump(root, f, separators=(",", ":"))


if __name__ == "__main__":
    sys.exit(main())
