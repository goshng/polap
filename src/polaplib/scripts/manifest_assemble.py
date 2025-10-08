#!/usr/bin/env python3
# scripts/manifest_assemble.py
#
# Version : v0.2.0
# Author  : Sang Chul Choi (POLAP)
# Date    : 2025-10-05
# License : GPL-3.0+
#
# Purpose:
#   Read a TSV facts table produced by Bash and build a structured JSON manifest.
#   Adds:
#     - GFA stats for each organelle
#     - Dataset metrics from seqkit stats
#
import os, sys, csv, json, argparse, subprocess, shlex, time
from collections import defaultdict


def read_facts(path):
    rows = []
    with open(path, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            rows.append(row)
    return rows


def call_gfa_stats_json(gfa_path):
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
        if out.strip() and out.strip() != "{}":
            return json.loads(out)
    except Exception:
        return None
    return None


def attach_stats(block):
    gfa = block.get("gfa")
    if not gfa:
        return
    st = call_gfa_stats_json(gfa)
    if st:
        block["stats"] = st


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--facts", required=True)
    ap.add_argument("--set", required=True)
    ap.add_argument("--tier", required=True)
    ap.add_argument("--inum", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--pretty", action="store_true", default=False)
    a = ap.parse_args()

    rows = read_facts(a.facts)
    data = {}
    root = dict(
        version="v0.2.0",
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
            data[sp] = dict(species=sp, tier=a.tier, inum=str(a.inum))
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
            try:
                v = float(val) if "." in val else int(val)
            except:
                v = val
            if key.startswith("circular_"):
                blk.setdefault("circular", {})[key.split("_", 1)[1]] = v
            elif key.startswith("gene_"):
                blk.setdefault("annotation", {})[key.split("_", 1)[1]] = v
            elif org == "data":
                blk[key] = v  # dataset metrics

    for sp, item in data.items():
        for org in ("pt", "ptpt", "mt"):
            if org in item:
                attach_stats(item[org])

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
    print(a.out)


if __name__ == "__main__":
    sys.exit(main())
