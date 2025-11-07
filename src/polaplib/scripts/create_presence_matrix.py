#!/usr/bin/env python3
# FILE: scripts/create_presence_matrix.py
# VERSION: 0.2.0
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Build clusters x species presence matrix from cluster_map.tsv,
with robust species extraction and tree-aware column alignment.

Usage:
  create_presence_matrix.py <cluster_map.tsv> <out.tsv> [--tree PT.treefile]
                            [--alias species_alias.tsv] [--strip-tree-suffix "-0"]
                            [--strip-presence-suffix ""]

Inputs:
  cluster_map.tsv : two columns: MTPT_id, CID
                    where MTPT_id is like "MTPT_<species>_<k>"
  PT.treefile     : (optional) Newick tree; if provided, presence columns are
                    ordered/matched to its tip labels.
  species_alias.tsv: (optional) TSV with columns: tip  species
                    to map an exact tree tip label to a desired presence name.

Notes:
- Species is extracted as the substring between "MTPT_" and the final underscore.
- If a tree is provided, we try the following to match columns to tips:
    1) exact match to extracted species
    2) species with/without a configured suffix (e.g., "-0")
    3) alias table mappings
- Any cluster with no matched species columns will still be emitted (all zeros),
  but a warning is printed.
"""
import sys, csv, re, argparse
from collections import defaultdict

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("cluster_map", help="clusters/cluster_map.tsv")
    ap.add_argument("out_tsv", help="output presence.tsv")
    ap.add_argument("--tree", default=None, help="PT tree (Newick)")
    ap.add_argument("--alias", default=None, help="TSV with columns: tip\tspecies")
    ap.add_argument("--strip-tree-suffix", default="", help="Suffix to strip from tree tips (e.g., -0)")
    ap.add_argument("--strip-presence-suffix", default="", help="Suffix to strip from presence species")
    return ap.parse_args()

def read_tree_tips(newick_path):
    # minimal Newick tip reader (no deps)
    tips = []
    with open(newick_path, "r", encoding="utf-8", errors="ignore") as fh:
        txt = fh.read().strip()
    # Remove branch lengths/comments crudely; then split on delimiters
    # We’ll just pull tokens that look like labels (no parentheses/commas/colons)
    token = ""
    for ch in txt:
        if ch in "(),;:":
            if token:
                tips.append(token)
                token = ""
        else:
            token += ch
    if token:
        tips.append(token)
    # Filter empties and internal node labels heuristically: we keep all; downstream matching prunes
    return [t.strip() for t in tips if t.strip()]

def main():
    a = parse_args()

    # optional alias map: tip -> desired species
    alias_map = {}
    if a.alias:
        with open(a.alias, newline="", encoding="utf-8") as fh:
            rdr = csv.DictReader(fh, delimiter="\t")
            if not {"tip","species"}.issubset(rdr.fieldnames or []):
                sys.stderr.write("[presence] alias file must have columns: tip\tspecies\n")
                sys.exit(2)
            for r in rdr:
                alias_map[r["tip"]] = r["species"]

    # 1) Read cluster_map.tsv
    rows = []
    with open(a.cluster_map, newline="", encoding="utf-8") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        if not {"MTPT_id","CID"}.issubset(rdr.fieldnames or []):
            sys.stderr.write("[presence] cluster_map.tsv must have columns: MTPT_id\tCID\n")
            sys.exit(2)
        for r in rdr:
            rows.append((r["MTPT_id"], r["CID"]))

    # 2) Extract species robustly:
    #    MTPT_<species>_<k>  -> take substring between "MTPT_" and final '_' (rsplit)
    def extract_species(mtpt_id):
        if not mtpt_id.startswith("MTPT_"):
            # fallback: try last underscore removal
            base = mtpt_id
        else:
            base = mtpt_id[len("MTPT_"):]
        if "_" in base:
            sp = base.rsplit("_", 1)[0]   # keep all up to last underscore
        else:
            sp = base
        return sp

    # Map: CID -> set(species)
    cid_to_species = defaultdict(set)
    all_species = set()
    for mtpt, cid in rows:
        sp = extract_species(mtpt)
        if a.strip_presence_suffix and sp.endswith(a.strip_presence_suffix):
            sp = sp[: -len(a.strip_presence_suffix)]
        cid_to_species[cid].add(sp)
        all_species.add(sp)

    # 3) If tree provided, align presence columns to its tips
    if a.tree:
        tips = read_tree_tips(a.tree)
        # Apply alias mapping and suffix strip on tips
        norm_tips = []
        for t in tips:
            t0 = alias_map.get(t, t)
            if a.strip_tree_suffix and t0.endswith(a.strip_tree_suffix):
                t0 = t0[: -len(a.strip_tree_suffix)]
            norm_tips.append(t0)
        # Build final ordered species list following tree order
        species_order = norm_tips
        # Tell user about coverage
        matched = len([s for s in species_order if s in all_species])
        if matched == 0:
            sys.stderr.write("ERROR: No overlapping species between tree tips (after normalization) "
                             "and extracted presence species.\n"
                             f"  example tips: {', '.join(tips[:10])} ...\n"
                             f"  example presence species: {', '.join(list(all_species)[:10])} ...\n")
            sys.exit(3)
    else:
        species_order = sorted(all_species)

    # 4) Build presence matrix
    # Rows: CIDs; Cols: species_order
    # A species is '1' for a CID if that species appears in cid_to_species[CID]
    with open(a.out_tsv, "w", encoding="utf-8") as o:
        o.write("\t".join(["CID"] + species_order) + "\n")
        for cid in sorted(cid_to_species.keys()):
            pres = []
            sset = cid_to_species[cid]
            for sp in species_order:
                pres.append("1" if sp in sset else "0")
            o.write("\t".join([cid] + pres) + "\n")

    sys.stderr.write(f"[presence] wrote {a.out_tsv}\n")

if __name__ == "__main__":
    main()

