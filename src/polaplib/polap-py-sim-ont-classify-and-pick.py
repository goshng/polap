#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
polap-py-classify-and-pick.py
Python replacement for AWK-based classification and greedy selection.
"""
import sys, os, argparse, gzip
from collections import defaultdict, namedtuple

Hit = namedtuple("Hit", "qname qlen tname tlen qstart qend tstart tend nmatch alen mapq")

def open_maybe_gzip(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode=mode)
    return open(path, mode=mode)

def fasta_length_total(path: str) -> int:
    total = 0
    with open(path, "rt") as fh:
        for ln in fh:
            if ln.startswith(">"):
                continue
            total += len(ln.strip())
    return total

def read_fastq_lengths(path: str):
    lens = {}
    ids = []
    with open_maybe_gzip(path, "rt") as fh:
        i = 0
        cur_id = None
        for ln in fh:
            i += 1
            if i % 4 == 1:
                cur_id = ln.strip()
                if cur_id.startswith("@"):
                    cur_id = cur_id[1:]
                ids.append(cur_id)
            elif i % 4 == 2:
                seq = ln.strip()
                if cur_id is not None:
                    lens[cur_id] = len(seq)
    return lens, ids

def parse_paf_hits(path: str):
    hits = []
    with open(path, "rt") as fh:
        for ln in fh:
            if not ln.strip(): 
                continue
            f = ln.rstrip("
").split("	")
            if len(f) < 12:
                continue
            qname, qlen, qstart, qend, strand, tname, tlen, tstart, tend, nmatch, alen, mapq = f[:12]
            try:
                h = Hit(qname=qname, qlen=int(qlen), tname=tname, tlen=int(tlen),
                        qstart=int(qstart), qend=int(qend), tstart=int(tstart), tend=int(tend),
                        nmatch=int(nmatch), alen=int(alen), mapq=int(mapq))
                hits.append(h)
            except ValueError:
                continue
    return hits

def best_class_by_span(hits, mt_name_hint, pt_name_hint):
    by_q = defaultdict(list)
    for h in hits:
        by_q[h.qname].append(h)
    best = {}
    for q, lst in by_q.items():
        best_span = -1
        best_mapq = -1
        best_tuple = None
        for h in lst:
            span = abs(h.tend - h.tstart)
            t = h.tname
            if t == mt_name_hint or ("mt" in t.lower() or "mito" in t.lower() or "mitochond" in t.lower()):
                cls = "mt"
            elif t == pt_name_hint or ("pt" in t.lower() or "plast" in t.lower() or "chloro" in t.lower()):
                cls = "pt"
            else:
                cls = "nuc"
            if span > best_span or (span == best_span and h.mapq > best_mapq):
                best_span = span
                best_mapq = h.mapq
                best_tuple = (cls, h.qlen, span, h.mapq)
        if best_tuple:
            best[q] = best_tuple
    return best

def write_sorted_candidates(best_assign, lens, out_prefix):
    bucket = {"nuc": [], "mt": [], "pt": []}
    for q, (cls, qlen, span, mapq) in best_assign.items():
        L = lens.get(q, qlen)
        bucket[cls].append((q, L))
    for cls in bucket:
        arr = sorted(bucket[cls], key=lambda x: x[1], reverse=True)
        with open(f"{out_prefix}.{cls}.candidates", "wt") as fo:
            for rid, L in arr:
                fo.write(f"{rid}	{L}
")
    return {cls: f"{out_prefix}.{cls}.candidates" for cls in bucket}

def greedy_ids(candidates_path, target_bases, out_ids_path):
    picked = 0
    total = 0
    with open(candidates_path, "rt") as fh, open(out_ids_path, "wt") as out:
        for ln in fh:
            rid, Ls = ln.rstrip("
").split("	")
            L = int(Ls)
            if total >= target_bases:
                break
            out.write(rid + "
")
            total += L
            picked += 1
    return picked, total

def main():
    ap = argparse.ArgumentParser(description="Classify reads (PAF) into nuc/mt/pt and pick IDs to reach target coverage.")
    ap.add_argument("--fastq", required=True, help="Mixed input reads (FASTQ or FASTQ.gz)")
    ap.add_argument("--paf", required=True, help="PAF alignments of FASTQ vs combined refs")
    ap.add_argument("--nuc-fa", required=True, help="nuclear FASTA (for length)")
    ap.add_argument("--mt-fa", required=True, help="mt FASTA (for length/name hint)")
    ap.add_argument("--pt-fa", required=True, help="pt FASTA (for length/name hint)")
    ap.add_argument("--depth-nuc", type=float, default=10.0)
    ap.add_argument("--depth-mt", type=float, default=50.0)
    ap.add_argument("--depth-pt", type=float, default=500.0)
    ap.add_argument("--out-prefix", default="out")
    args = ap.parse_args()

    lens, all_ids = read_fastq_lengths(args.fastq)
    all_id_set = set(all_ids)

    hits = parse_paf_hits(args.paf)

    def first_header_name(path):
        with open(path, "rt") as fh:
            for ln in fh:
                if ln.startswith(">"):
                    return ln[1:].strip().split()[0]
        return None

    mt_name = first_header_name(args.mt_fa) or "mt"
    pt_name = first_header_name(args.pt_fa) or "pt"

    best = best_class_by_span(hits, mt_name, pt_name)  # q -> (cls, qlen, span, mapq)

    mapped_names = set(best.keys())
    unmapped = list(all_id_set - mapped_names)
    for q in unmapped:
        L = lens.get(q, 0)
        best[q] = ("nuc", L, 0, 0)

    with open(f"{args.out_prefix}.best.tsv", "wt") as fo:
        for q, (cls, qlen, span, mapq) in best.items():
            fo.write(f"{q}	{cls}	{qlen}
")

    cand_paths = write_sorted_candidates(best, lens, args.out_prefix)

    Ln = fasta_length_total(args.nuc_fa)
    Lm = fasta_length_total(args.mt_fa)
    Lp = fasta_length_total(args.pt_fa)
    Tn = int(round(args.depth_nuc * Ln))
    Tm = int(round(args.depth_mt * Lm))
    Tp = int(round(args.depth_pt * Lp))

    stats = []
    n_picked, n_bp = greedy_ids(cand_paths["nuc"], Tn, f"{args.out_prefix}.nuc.ids")
    stats.append(("nuc", n_picked, n_bp, Tn))
    m_picked, m_bp = greedy_ids(cand_paths["mt"], Tm, f"{args.out_prefix}.mt.ids")
    stats.append(("mt", m_picked, m_bp, Tm))
    p_picked, p_bp = greedy_ids(cand_paths["pt"], Tp, f"{args.out_prefix}.pt.ids")
    stats.append(("pt", p_picked, p_bp, Tp))

    with open(f"{args.out_prefix}.summary.txt", "wt") as fo:
        for cls, k, bp, tgt in stats:
            fo.write(f"{cls}	reads={k}	bases={bp}	target={tgt}
")

if __name__ == "__main__":
    main()
