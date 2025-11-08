#!/usr/bin/env python3
# FILE: scripts/create_presence_matrix.py
# VERSION: 0.2.1
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Build clusters x species presence matrix from cluster_map.tsv,
with robust species extraction and tree-aware column alignment.

Usage:
  create_presence_matrix.py <cluster_map.tsv> <out.tsv> [--tree PT.treefile]
                            [--alias species_alias.tsv]
                            [--strip-tree-suffix "-0"]
                            [--strip-presence-suffix ""]
"""

import sys, csv, re, argparse
from collections import defaultdict


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("cluster_map", help="clusters/cluster_map.tsv")
    ap.add_argument("out_tsv", help="output presence.tsv")
    ap.add_argument("--tree", default=None, help="PT tree (Newick)")
    ap.add_argument("--alias", default=None, help="TSV with columns: tip\tspecies")
    ap.add_argument(
        "--strip-tree-suffix",
        default="",
        help="Suffix to strip from tree tips (e.g., -0)",
    )
    ap.add_argument(
        "--strip-presence-suffix",
        default="",
        help="Suffix to strip from presence species",
    )
    return ap.parse_args()


def read_tree_tips_biophylo(newick_path):
    try:
        from Bio import Phylo

        tree = Phylo.read(newick_path, "newick")
        # Only leaf (terminal) names
        tips = [t.name for t in tree.get_terminals() if t.name]
        return tips
    except Exception:
        return None


def read_tree_tips_manual(newick_path):
    # Manual, leaf-only: capture labels after '(' or ',' up to delimiters :,()[];
    with open(newick_path, "r", encoding="utf-8", errors="ignore") as fh:
        s = fh.read()

    tips = []
    i = 0
    n = len(s)
    while i < n:
        ch = s[i]
        if ch in "(,":
            i += 1
            # skip whitespace
            while i < n and s[i].isspace():
                i += 1
            # if subtree starts, no label here
            if i < n and s[i] == "(":
                continue
            # quoted label?
            if i < n and s[i] in "'\"":
                q = s[i]
                i += 1
                start = i
                while i < n and s[i] != q:
                    i += 1
                label = s[start:i]
                i += 1
            else:
                start = i
                while i < n and s[i] not in "(),:[];":
                    i += 1
                label = s[start:i].strip()
            if label:
                tips.append(label)
        else:
            i += 1
    return tips


def read_tree_tips(newick_path):
    tips = read_tree_tips_biophylo(newick_path)
    if tips is not None:
        return tips
    return read_tree_tips_manual(newick_path)


def extract_species(mtpt_id):
    # MTPT_<species>_<k>  -> substring between "MTPT_" and final '_'
    base = mtpt_id[len("MTPT_") :] if mtpt_id.startswith("MTPT_") else mtpt_id
    return base.rsplit("_", 1)[0] if "_" in base else base


def main():
    a = parse_args()

    # Auto-harmonize suffix stripping if only tree-side was provided
    strip_tree = a.strip_tree_suffix or ""
    strip_presence = (
        a.strip_presence_suffix if a.strip_presence_suffix != "" else strip_tree
    )

    # optional alias map: tip -> desired species
    alias_map = {}
    if a.alias:
        with open(a.alias, newline="", encoding="utf-8") as fh:
            rdr = csv.DictReader(fh, delimiter="\t")
            if not {"tip", "species"}.issubset(rdr.fieldnames or []):
                sys.stderr.write(
                    "[presence] alias file must have columns: tip\tspecies\n"
                )
                sys.exit(2)
            for r in rdr:
                alias_map[r["tip"]] = r["species"]

    # 1) Read cluster_map.tsv
    rows = []
    with open(a.cluster_map, newline="", encoding="utf-8") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        if not {"MTPT_id", "CID"}.issubset(rdr.fieldnames or []):
            sys.stderr.write(
                "[presence] cluster_map.tsv must have columns: MTPT_id\tCID\n"
            )
            sys.exit(2)
        for r in rdr:
            rows.append((r["MTPT_id"], r["CID"]))

    # 2) Extract normalized species from MTPT ids
    cid_to_species = defaultdict(set)
    all_species = set()
    for mtpt, cid in rows:
        sp = extract_species(mtpt)
        if strip_presence and sp.endswith(strip_presence):
            sp = sp[: -len(strip_presence)]
        cid_to_species[cid].add(sp)
        all_species.add(sp)

    # 3) If tree provided, align presence columns to the tree's (normalized) leaf order
    if a.tree:
        raw_tips = read_tree_tips(a.tree)
        # Apply alias and strip suffix for tree tips
        norm_tips = []
        for t in raw_tips:
            t0 = alias_map.get(t, t)
            if strip_tree and t0.endswith(strip_tree):
                t0 = t0[: -len(strip_tree)]
            norm_tips.append(t0)

        # Order: first those present in matrix and in tree order, then any remaining species (sorted)
        in_both = [t for t in norm_tips if t in all_species]
        if len(in_both) == 0:
            # Show normalized sets to help debugging
            sys.stderr.write(
                "ERROR: No overlapping species between tree tips (after normalization) and extracted presence species.\n"
                f"  example norm tips: {', '.join(norm_tips[:10])} ...\n"
                f"  example presence species: {', '.join(list(all_species)[:10])} ...\n"
            )
            sys.exit(3)
        remaining = sorted(all_species.difference(in_both))
        species_order = in_both + remaining
    else:
        species_order = sorted(all_species)

    # 4) Write presence matrix
    with open(a.out_tsv, "w", encoding="utf-8") as o:
        o.write("\t".join(["CID"] + species_order) + "\n")
        for cid in sorted(cid_to_species.keys()):
            sset = cid_to_species[cid]
            row = ["1" if sp in sset else "0" for sp in species_order]
            o.write("\t".join([cid] + row) + "\n")

    sys.stderr.write(f"[presence] wrote {a.out_tsv}\n")


if __name__ == "__main__":
    main()
