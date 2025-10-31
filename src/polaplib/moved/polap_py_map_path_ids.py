#!/usr/bin/env python3
# scripts/polap_py_map_path_ids.py
# Version: v0.1.0
# Input: 2-col TSV (old_id\tnew_id), and a comma-separated path "old1+,old2-,..."
# Output: "new1+,new2-,..."
import sys, argparse

ap = argparse.ArgumentParser()
ap.add_argument("--map", required=True)
ap.add_argument("--path", required=True)
a = ap.parse_args()

m = {}
with open(a.map) as f:
    for ln in f:
        old, new = ln.rstrip("\n").split("\t")
        m[old] = new

out = []
for tok in [t.strip() for t in a.path.split(",") if t.strip()]:
    if tok[-1] not in "+-":
        sys.exit(f"ERR: bad token (need +/-): {tok}")
    nid, ori = tok[:-1], tok[-1]
    if nid not in m:
        sys.exit(f"ERR: id not in map: {nid}")
    out.append(m[nid] + ori)

print(",".join(out))
