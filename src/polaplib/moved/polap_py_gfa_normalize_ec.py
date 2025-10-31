#!/usr/bin/env python3
# scripts/polap_py_gfa_normalize_ec.py
# Version: v0.5.0
# Input : --gfa in.gfa
# Output: --out out.norm.gfa  (IDs 1..N; L overlap '*'->'0M'; S ensure ll:f; L ensure ec:i)
#         --idmap out.idmap.tsv (old_id \t new_id)
# Options: --ec-mode mean|min|max|const, --ec-const 1, --ec-scale 1.0, --ec-round round|ceil|floor
import sys, argparse, gzip, math


def openg(path, mode="rt"):
    if path == "-":
        return sys.stdin
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gfa", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--idmap", required=True)
    ap.add_argument(
        "--ec-mode", default="mean", choices=["mean", "min", "max", "const"]
    )
    ap.add_argument("--ec-const", default="1")
    ap.add_argument("--ec-scale", default="1.0")
    ap.add_argument("--ec-round", default="round", choices=["round", "ceil", "floor"])
    return ap.parse_args()


def rnd(x, mode):
    if mode == "ceil":
        return math.ceil(x)
    if mode == "floor":
        return math.floor(x)
    # nearest int
    return int(x + 0.5) if x >= 0 else int(x - 0.5)


def build_idmap_and_ll(gfa_path):
    """First pass: assign numeric IDs in S-order; collect ll from existing tags or dp fallback."""
    idmap = {}
    ll = {}  # new_id -> float
    next_id = 1
    with openg(gfa_path, "rt") as f:
        for ln in f:
            if not ln or ln[0] in "#\n":
                continue
            cols = ln.rstrip("\n").split("\t")
            if cols[0] == "S":
                old = cols[1]
                if old not in idmap:
                    idmap[old] = str(next_id)
                    next_id += 1
                # scan tags for ll:f / dp:i / dp:f
                llv = None
                dpv = None
                for tag in cols[3:]:
                    if tag.startswith("ll:f:"):
                        try:
                            llv = float(tag.split(":", 2)[2])
                        except:
                            llv = None
                    elif tag.startswith("dp:i:") or tag.startswith("dp:f:"):
                        try:
                            dpv = float(tag.split(":", 2)[2])
                        except:
                            dpv = None
                new_id = idmap[old]
                if llv is not None:
                    ll[new_id] = llv
                else:
                    ll[new_id] = dpv if dpv is not None else 1.0
    if next_id == 1:
        raise SystemExit("ERR: no S lines found")
    return idmap, ll


def ensure_ec(u, v, existing_tags, ll, mode, cconst, scale, rmode):
    """Return tags list ensuring ec:i:*, preserving existing ec if present."""
    have_ec = any(t.startswith("ec:i:") for t in existing_tags)
    if have_ec:
        return existing_tags
    covu = ll.get(u, 1.0)
    covv = ll.get(v, 1.0)
    if mode == "min":
        ec = min(covu, covv)
    elif mode == "max":
        ec = max(covu, covv)
    elif mode == "mean":
        ec = 0.5 * (covu + covv)
    else:
        ec = float(cconst)
    ec *= float(scale)
    ec_int = max(0, rnd(ec, rmode))
    return existing_tags + [f"ec:i:{ec_int}"]


def main():
    a = parse_args()
    idmap, ll = build_idmap_and_ll(a.gfa)

    with open(a.idmap, "w") as w:
        for old, new in idmap.items():
            w.write(f"{old}\t{new}\n")

    out = open(a.out, "w")
    with openg(a.gfa, "rt") as f:
        for ln in f:
            if not ln or ln[0] == "#":
                continue
            cols = ln.rstrip("\n").split("\t")
            typ = cols[0]
            if typ == "H":
                out.write(ln)  # pass-through
            elif typ == "S":
                old = cols[1]
                if old not in idmap:
                    raise SystemExit(f"ERR: S id not in map: {old}")
                cols[1] = idmap[old]
                # ensure ll:f present, using first-pass ll table
                has_ll = any(t.startswith("ll:f:") for t in cols[3:])
                if not has_ll:
                    cols.append(f"ll:f:{ll[cols[1]]}")
                out.write("\t".join(cols) + "\n")
            elif typ == "L":
                u_old, uo, v_old, vo, cigar = cols[1:6]
                if u_old not in idmap or v_old not in idmap:
                    raise SystemExit(f"ERR: L endpoint not in map: {u_old} {v_old}")
                cols[1] = idmap[u_old]
                cols[3] = idmap[v_old]
                if cigar == "*":
                    cols[5] = "0M"  # simple overlap for gfatk
                # ensure ec:i
                tags = cols[6:] if len(cols) > 6 else []
                tags = ensure_ec(
                    cols[1],
                    cols[3],
                    tags,
                    ll,
                    a.ec_mode,
                    a.ec_const,
                    a.ec_scale,
                    a.ec_round,
                )
                cols = cols[:6] + tags
                out.write("\t".join(cols) + "\n")
            else:
                out.write(ln)
    out.close()


if __name__ == "__main__":
    main()
