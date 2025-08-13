#!/usr/bin/env python3
import sys, argparse, gzip
from collections import defaultdict, namedtuple

Edge = namedtuple("Edge", "to strand ovl iden score")


def open_text(fn):
    if fn in (None, "-"):
        return sys.stdin
    return gzip.open(fn, "rt") if fn.endswith(".gz") else open(fn)


def rc(seq):
    tbl = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(tbl)[::-1]


def read_fx(fn):
    seqs = {}
    with open_text(fn) as f:
        first = f.readline()
        if not first:
            return seqs
        is_fastq = first.startswith("@")
        hdr = first if first[:1] in "@>" else None
        buf = []

        def commit():
            if not hdr:
                return
            name = hdr.split()[0][1:]
            if name not in seqs:
                seqs[name] = "".join(buf)

        for ln in f:
            c = ln[:1]
            if c in "@>":
                commit()
                hdr = ln
                buf = []
            elif is_fastq and c == "+":
                # skip qualities
                need = sum(len(x.strip()) for x in buf)
                got = 0
                for q in f:
                    got += len(q.strip())
                    if got >= need:
                        break
                commit()
                hdr = None
                nxt = f.readline()
                if not nxt:
                    break
                if nxt[:1] in "@>":
                    hdr = nxt
                    buf = []
                else:
                    raise ValueError("FASTQ parsing error")
            else:
                if ln.strip():
                    buf.append(ln.strip())
        if hdr:
            commit()
    return seqs


def parse_paf(paf_path, min_id, min_ovl):
    edges = defaultdict(list)
    with open_text(paf_path) as f:
        for ln in f:
            if not ln or ln[0] == "#":
                continue
            a = ln.rstrip("\n").split("\t")
            if len(a) < 12:
                continue
            qn, ql, qs, qe, st, tn, tl, ts, te, nm, al, mapq = a[:12]
            ql = int(ql)
            qs = int(qs)
            qe = int(qe)
            tl = int(tl)
            ts = int(ts)
            te = int(te)
            nm = int(nm)
            al = int(al)
            if qn == tn:
                continue
            ovl = min(qe - qs, te - ts)
            if ovl < min_ovl:
                continue
            iden = (nm / al) if al > 0 else 0.0
            if iden < min_id:
                continue
            score = iden * ovl
            edges[qn].append(Edge(tn, st, ovl, iden, score))
            # symmetrical listing for reciprocal lookup
            edges[tn].append(Edge(qn, st, ovl, iden, score))
    # sort neighbors by score desc
    for k in edges:
        edges[k].sort(key=lambda e: (e.score, e.iden, e.ovl), reverse=True)
    return edges


def reciprocal_best(a, b, edges, topk=3):
    """Return True if b lists a within its top-k neighbors."""
    nb = edges.get(b, [])
    limit = min(topk, len(nb))
    for i in range(limit):
        if nb[i].to == a:
            return True
    return False


def too_branchy(cur, candidate, edges, delta=0.10):
    """Avoid nodes whose neighbors have many near-tie scores (repeat-ish)."""
    neigh = edges.get(candidate, [])
    if not neigh:
        return False
    best = neigh[0].score
    if best <= 0:
        return False
    # count neighbors within (1-delta)*best
    near = sum(1 for e in neigh if e.score >= (1.0 - delta) * best)
    return near >= 5  # heuristic; tweak if needed


def orient_join(cur_or, st, t_seq):
    # st is '+' if query and target have same orientation.
    t_or = +1 if ((cur_or == +1 and st == "+") or (cur_or == -1 and st == "-")) else -1
    return (t_seq if t_or == +1 else rc(t_seq)), t_or


def choose_next(current, used, edges, topk, branch_delta):
    for e in edges.get(current, []):
        if e.to in used:
            continue
        # reciprocal-best filter
        if not reciprocal_best(current, e.to, edges, topk=topk):
            continue
        # branchiness guard
        if too_branchy(current, e.to, edges, delta=branch_delta):
            continue
        return e
    return None


def build_contig(seed, seqs, edges, target_len, topk, branch_delta):
    used = set([seed])
    path = [(seed, +1)]
    seq = seqs[seed]
    # extend right
    cur = seed
    cur_or = +1
    while len(seq) < target_len:
        e = choose_next(cur, used, edges, topk, branch_delta)
        if e is None:
            break
        t_seq, t_or = orient_join(cur_or, e.strand, seqs[e.to])
        if e.ovl >= len(t_seq) or e.ovl >= len(seq):
            used.add(e.to)
            cur = e.to
            cur_or = t_or
            continue
        seq = seq + t_seq[e.ovl :]
        used.add(e.to)
        path.append((e.to, t_or))
        cur = e.to
        cur_or = t_or

    # extend left from the first element
    left_cur, left_or = path[0]
    while len(seq) < target_len:
        e = choose_next(left_cur, used, edges, topk, branch_delta)
        if e is None:
            break
        t_seq, t_or = orient_join(left_or, e.strand, seqs[e.to])
        if e.ovl >= len(t_seq) or e.ovl >= len(seq):
            used.add(e.to)
            left_cur = e.to
            left_or = t_or
            continue
        seq = t_seq + seq[e.ovl :]
        used.add(e.to)
        path.insert(0, (e.to, t_or))
        left_cur = e.to
        left_or = t_or

    return seq, path, used


def main():
    ap = argparse.ArgumentParser(
        description="Greedy ~100 kb constructor from ONT PAF overlaps."
    )
    ap.add_argument("--reads", required=True)
    ap.add_argument("--paf", required=True)
    ap.add_argument("--min-id", type=float, default=0.90)
    ap.add_argument("--min-ovl", type=int, default=4000)
    ap.add_argument("--target", type=int, default=100000)
    ap.add_argument("--max-contigs", type=int, default=5)
    ap.add_argument("--topk", type=int, default=3, help="reciprocal-best top-K")
    ap.add_argument(
        "--branch_delta",
        type=float,
        default=0.10,
        help="near-tie window for branchiness",
    )
    ap.add_argument("--out", required=True)
    ap.add_argument("--paths", default=None, help="TSV: unitig/read/order/orient")
    args = ap.parse_args()

    print(f"[INFO] Loading reads: {args.reads}", file=sys.stderr)
    seqs = read_fx(args.reads)
    print(f"[INFO] Reads loaded: {len(seqs)}", file=sys.stderr)

    print(f"[INFO] Reading PAF: {args.paf}", file=sys.stderr)
    edges = parse_paf(args.paf, args.min_id, args.min_ovl)
    print(f"[INFO] Nodes with edges: {len(edges)}", file=sys.stderr)

    seeds = sorted(seqs.keys(), key=lambda k: len(seqs[k]), reverse=True)
    used_global = set()
    written = 0

    paths_out = open(args.paths, "w") if args.paths else None
    if paths_out:
        paths_out.write("unitig\torder\tread\torient\n")

    with open(args.out, "w") as out:
        for s in seeds:
            if s in used_global:
                continue
            contig, path, used_local = build_contig(
                s, seqs, edges, args.target, args.topk, args.branch_delta
            )
            if len(contig) >= min(
                args.target, int(0.7 * args.target)
            ):  # accept >=70% of target
                written += 1
                used_global |= {rid for rid, _ in path}
                name = f"unitig{written}_len{len(contig)}_seed={s}_nreads={len(path)}"
                out.write(f">{name}\n")
                for i in range(0, len(contig), 80):
                    out.write(contig[i : i + 80] + "\n")
                if paths_out:
                    for i, (rid, orient) in enumerate(path, start=1):
                        paths_out.write(f"{name}\t{i}\t{rid}\t{orient}\n")
                print(f"[INFO] wrote {name}", file=sys.stderr)
            if written >= args.max_contigs:
                break

    if paths_out:
        paths_out.close()
    if written == 0:
        print(
            "[WARN] No contig reached sufficient length. Relax --min-id/--min-ovl or recruit reads.",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
