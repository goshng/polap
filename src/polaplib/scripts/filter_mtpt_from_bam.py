#!/usr/bin/env python3
# scripts/filter_mtpt_from_bam.py
# Version: v0.1.1
# License: GPL-3.0+
"""
Classify ONT reads as mt.keep / pt.drop / ambig / none by comparing best
minimap2 alignments to mt and pt references using:
  - identity ≈ 1 - NM / query_alignment_length
  - aligned fraction = query_alignment_length / query_length
  - clipping = query_length - query_alignment_length
Best hit per reference is chosen by score = identity * query_alignment_length.
NEW: --primary-only, --min-mapq, --margin (default 0.10) for tie-breaking.
"""

import argparse, sys
from collections import defaultdict
import pysam


def best_hits_from_bam(bam_path, primary_only=False, min_mapq=0):
    best = {}  # qname -> dict(metrics)
    bam = pysam.AlignmentFile(bam_path, "rb")
    for aln in bam.fetch(until_eof=True):
        if aln.is_unmapped:
            continue
        if primary_only and (aln.is_secondary or aln.is_supplementary):
            continue
        if aln.mapping_quality < min_mapq:
            continue

        qn = aln.query_name
        qlen = aln.query_length or 0
        alnlen = aln.query_alignment_length or 0
        if qlen <= 0 or alnlen <= 0:
            continue

        try:
            nm = aln.get_tag("NM")
        except KeyError:
            # no NM -> cannot compute identity robustly; skip this alignment
            continue

        ident = max(0.0, 1.0 - (nm / max(1, alnlen)))
        frac = alnlen / qlen
        clip = qlen - alnlen
        score = ident * alnlen

        prev = best.get(qn)
        if (prev is None) or (score > prev["score"]):
            best[qn] = dict(
                qlen=qlen, alnlen=alnlen, ident=ident, frac=frac, clip=clip, score=score
            )
    bam.close()
    return best


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pt-bam", required=True)
    ap.add_argument("--mt-bam", required=True)
    ap.add_argument("--id-mt", type=float, default=0.75)
    ap.add_argument("--frac-mt", type=float, default=0.40)
    ap.add_argument("--id-pt", type=float, default=0.85)
    ap.add_argument("--frac-pt", type=float, default=0.30)
    ap.add_argument("--max-clip", type=int, default=4000)
    ap.add_argument("--min-read", type=int, default=2000)
    ap.add_argument("--ambig", choices=["drop", "keep", "none"], default="drop")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument(
        "--emit-metrics", default=None, help="write per-read metrics TSV here"
    )
    # NEW options
    ap.add_argument(
        "--primary-only", action="store_true", help="use only primary alignments"
    )
    ap.add_argument(
        "--min-mapq", type=int, default=0, help="ignore alignments with MAPQ < N"
    )
    ap.add_argument(
        "--margin",
        type=float,
        default=0.10,
        help="tie-break margin as fraction (0.10 = 10%)",
    )
    return ap.parse_args()


def main():
    args = parse_args()

    pt = best_hits_from_bam(
        args.pt_bam, primary_only=args.primary_only, min_mapq=args.min_mapq
    )
    mt = best_hits_from_bam(
        args.mt_bam, primary_only=args.primary_only, min_mapq=args.min_mapq
    )
    qnames = set(pt.keys()) | set(mt.keys())

    def pass_mt(m):
        return (
            m
            and m["ident"] >= args.id_mt
            and m["frac"] >= args.frac_mt
            and m["clip"] <= args.max_clip
            and m["qlen"] >= args.min_read
        )

    def pass_pt(m):
        return (
            m
            and m["ident"] >= args.id_pt
            and m["frac"] >= args.frac_pt
            and m["clip"] <= args.max_clip
            and m["qlen"] >= args.min_read
        )

    keep, drop, ambig, none = [], [], [], []
    margin = max(0.0, args.margin)  # ensure non-negative

    metrics_fp = open(args.emit_metrics, "w") if args.emit_metrics else None
    if metrics_fp:
        metrics_fp.write(
            "read\tclass\t"
            "mt_ident\tmt_frac\tmt_clip\tmt_score\t"
            "pt_ident\tpt_frac\tpt_clip\tpt_score\t"
            "qlen\n"
        )

    for q in qnames:
        m_mt = mt.get(q)
        m_pt = pt.get(q)
        mt_ok = pass_mt(m_mt)
        pt_ok = pass_pt(m_pt)
        mt_score = m_mt["score"] if m_mt else 0.0
        pt_score = m_pt["score"] if m_pt else 0.0

        if pt_ok and not mt_ok:
            cls = "pt.drop"
            drop.append(q)
        elif mt_ok and not pt_ok:
            cls = "mt.keep"
            keep.append(q)
        elif mt_ok and pt_ok:
            # tie-break using multiplicative margin
            if pt_score >= mt_score * (1.0 + margin):
                cls = "pt.drop"
                drop.append(q)
            elif mt_score >= pt_score * (1.0 + margin):
                cls = "mt.keep"
                keep.append(q)
            else:
                cls = "ambig"
                ambig.append(q)
        else:
            cls = "none"
            none.append(q)

        if metrics_fp:
            metrics_fp.write(
                f"{q}\t{cls}\t"
                f"{(m_mt['ident'] if m_mt else 'NA')}\t{(m_mt['frac'] if m_mt else 'NA')}\t{(m_mt['clip'] if m_mt else 'NA')}\t{mt_score}\t"
                f"{(m_pt['ident'] if m_pt else 'NA')}\t{(m_pt['frac'] if m_pt else 'NA')}\t{(m_pt['clip'] if m_pt else 'NA')}\t{pt_score}\t"
                f"{(m_mt['qlen'] if m_mt else (m_pt['qlen'] if m_pt else 'NA'))}\n"
            )

    if metrics_fp:
        metrics_fp.close()

    # write lists
    def write_list(path, arr):
        with open(path, "w") as f:
            for x in arr:
                f.write(x + "\n")

    write_list(f"{args.out_prefix}.keep_ids.txt", keep)
    write_list(f"{args.out_prefix}.drop_ids.txt", drop)
    write_list(f"{args.out_prefix}.ambig_ids.txt", ambig)
    write_list(f"{args.out_prefix}.none_ids.txt", none)


if __name__ == "__main__":
    main()
