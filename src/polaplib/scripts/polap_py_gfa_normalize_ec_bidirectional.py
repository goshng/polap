#!/usr/bin/env python3
# Version: v0.7.0
# Normalize GFA for gfatk (no awk):
#  - Ensure header H	VN:Z:1.0 if missing
#  - Renumber S IDs to 1..N; write idmap
#  - On S: ensure ll:f (from dp:i/f fallback 1.0)
#  - On L: coerce CIGAR to <int>M (sum M; '*'->'0M'), ensure ec:i (mean/min/max/const from endpoints' ll)
#  - Ensure reciprocal L edges exist (with same cleaned CIGAR + ec:i)
import sys, argparse, gzip, math, re


def openg(p, mode="rt"):
    if p == "-":
        return sys.stdin
    return gzip.open(p, mode) if p.endswith(".gz") else open(p, mode)


def rnd(x, mode):
    if mode == "ceil":
        return math.ceil(x)
    if mode == "floor":
        return math.floor(x)
    return int(x + 0.5) if x >= 0 else int(x - 0.5)


cigar_tok = re.compile(r"(\d+)([MIDNSHP=X])")


def clean_cigar(cg: str) -> str:
    if not cg or cg == "*":
        return "0M"
    msum = 0
    for n, op in cigar_tok.findall(cg):
        if op == "M":
            msum += int(n)
    return f"{msum}M"


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


def main():
    a = parse_args()

    # Pass1: assign numeric IDs in S-order; collect ll (node coverage)
    idmap = {}
    ll = {}  # new_id -> float
    order = []
    with openg(a.gfa) as f:
        for ln in f:
            if not ln or ln[0] in "#\n":
                continue
            c = ln.rstrip("\n").split("\t")
            if c[0] == "S":
                old = c[1]
                if old not in idmap:
                    nid = str(len(idmap) + 1)
                    idmap[old] = nid
                    order.append(old)
                # pull ll or dp
                llv, dpv = None, None
                for t in c[3:]:
                    if t.startswith("ll:f:"):
                        try:
                            llv = float(t.split(":", 2)[2])
                        except:
                            pass
                    elif t.startswith("dp:i:") or t.startswith("dp:f:"):
                        try:
                            dpv = float(t.split(":", 2)[2])
                        except:
                            pass
                new = idmap[old]
                ll[new] = llv if llv is not None else (dpv if dpv is not None else 1.0)

    if not idmap:
        sys.exit("ERR: no S lines found")

    # Decide if input has header
    has_header = False
    with openg(a.gfa) as f:
        for ln in f:
            if ln.startswith("H\t"):
                has_header = True
                break

    # Pass2: write normalized and collect seen edges
    out = open(a.out, "w")
    if not has_header:
        out.write("H\tVN:Z:1.0\n")

    edges = set()  # (u,uo,v,vo,cg,ec)
    with openg(a.gfa) as f:
        for ln in f:
            if not ln or ln[0] == "#":
                continue
            c = ln.rstrip("\n").split("\t")
            t = c[0]
            if t == "H":
                out.write(ln)
            elif t == "S":
                old = c[1]
                new = idmap.get(old)
                if not new:
                    sys.exit(f"ERR: S id not in map: {old}")
                c[1] = new
                has_ll = any(x.startswith("ll:f:") for x in c[3:])
                if not has_ll:
                    c.append(f"ll:f:{ll[new]}")
                out.write("\t".join(c) + "\n")
            elif t == "L":
                u_old, uo, v_old, vo, cg = c[1:6]
                u = idmap.get(u_old)
                v = idmap.get(v_old)
                if not u or not v:
                    sys.exit(f"ERR: L endpoint not in map: {u_old} {v_old}")
                c[1], c[3] = u, v
                cg = clean_cigar(cg)
                c[5] = cg
                # ensure ec:i
                tags = c[6:] if len(c) > 6 else []
                ec = None
                for t2 in tags:
                    if t2.startswith("ec:i:"):
                        try:
                            ec = int(t2.split(":", 2)[2])
                        except:
                            pass
                        break
                if ec is None:
                    covu = ll.get(u, 1.0)
                    covv = ll.get(v, 1.0)
                    if a.ec_mode == "min":
                        e = min(covu, covv)
                    elif a.ec_mode == "max":
                        e = max(covu, covv)
                    elif a.ec_mode == "mean":
                        e = 0.5 * (covu + covv)
                    else:
                        e = float(a.ec_const)
                    e *= float(a.ec_scale)
                    ec = max(0, rnd(e, a.ec_round))
                    tags.append(f"ec:i:{ec}")
                out.write("\t".join(c[:6] + tags) + "\n")
                edges.add((u, uo, v, vo, cg, ec))
            else:
                out.write(ln)

    # Add missing reciprocal edges
    def ro(o):
        return "+" if o == "-" else "-"

    existing = set((u, uo, v, vo) for (u, uo, v, vo, _, _) in edges)
    added = 0
    for u, uo, v, vo, cg, ec in list(edges):
        ru, rvo, rv, ruo = v, ro(vo), u, ro(uo)
        if (ru, rvo, rv, ruo) not in existing:
            out.write("\t".join(["L", ru, rvo, rv, ruo, cg, f"ec:i:{ec}"]) + "\n")
            added += 1
    out.close()

    with open(a.idmap, "w") as m:
        for old in order:
            m.write(f"{old}\t{idmap[old]}\n")

    sys.stderr.write(f"[INFO] Added {added} reciprocal L edges\n")


if __name__ == "__main__":
    main()
