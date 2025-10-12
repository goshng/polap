#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Two-boundary verifier for mt/pt junctions.

Given an edge verify directory (e.g., verify/edge_1_137249_139540),
expects a JUNC BAM with synthetic contigs JLEFT/JRIGHT.
Finds the best boundary positions by sweeping (W=100, step=10),
chooses W adaptively if needed, counts spanners (any vs primary-only),
decides ACCEPT/BORDERLINE/REJECT, and writes verify/.../twoboundary.tsv.

Outputs (tab-separated):
seqid start end pos_left pos_right w_left w_right any_left pri_left any_right pri_right twob_decision twob_note
"""
import argparse, os, sys, csv, re
from dataclasses import dataclass
from typing import Tuple, Optional
import pysam

# ----------------------- utils -----------------------

def parse_edge_dir_name(basename: str) -> Tuple[str,int,int]:
    # edge_1_137249_139540  -> ("edge_1", 137249, 139540)
    m = re.match(r'^(.*)_(\d+)_(\d+)$', basename)
    if not m:
        raise ValueError(f"Unrecognized edge directory name: {basename}")
    return m.group(1), int(m.group(2)), int(m.group(3))

def pct_len_ge_2W(bam: pysam.AlignmentFile, ctg: str, W: int) -> float:
    ok = tot = 0
    for read in bam.fetch(ctg):
        if read.is_unmapped: 
            continue
        L = (read.reference_end or 0) - (read.reference_start or 0)
        if L >= 2*W:
            ok += 1
        tot += 1
    return (100.0*ok/tot) if tot else 0.0

def count_spanners(bam: pysam.AlignmentFile, ctg: str, pos1: int, W: int, include_supp: bool) -> int:
    """
    Count unique read names that span [pos1-W, pos1+W] on reference (1-based pos).
    """
    names = set()
    for read in bam.fetch(ctg):
        if read.is_unmapped: 
            continue
        if (read.is_secondary or read.is_supplementary) and not include_supp:
            continue
        # BAM is 0-based, half-open
        start1 = (read.reference_start or 0) + 1
        end1   = (read.reference_end   or 0)     # already 1-based exclusive -> end coord as 1-based inclusive is end1
        # We want inclusive end: reference_end is exclusive; subtract 0 because comparison uses >=
        if start1 <= (pos1 - W) and end1 >= (pos1 + W):
            names.add(read.query_name)
    return len(names)

def best_pos(bam: pysam.AlignmentFile, ctg: str, W: int = 100, step: int = 10) -> int:
    """
    Sweep positions (step) to find the POS with max possible spanners.
    Uses only the constraint L >= 2W (ref span).
    """
    tid = bam.get_tid(ctg)
    if tid < 0:
        raise ValueError(f"Contig {ctg} not found in BAM")
    ctg_len = bam.get_reference_length(ctg)
    span = {}
    for read in bam.fetch(ctg):
        if read.is_unmapped:
            continue
        s1 = (read.reference_start or 0) + 1
        e1 = (read.reference_end   or 0)
        if e1 - s1 + 1 < 2*W:
            continue
        pmin = s1 + W
        pmax = e1 - W
        if pmin > pmax:
            continue
        # bin to step
        bmin = (pmin // step) * step
        bmax = (pmax // step) * step
        if bmin < step:
            bmin = step
        if bmax > ctg_len - step:
            bmax = ctg_len - step
        p = bmin
        while p <= bmax:
            span[p] = span.get(p, 0) + 1
            p += step
    if not span:
        # fallback: center
        return ctg_len // 2 or 1
    # argmax
    best_p, best_n = None, -1
    for p, n in span.items():
        if n > best_n:
            best_p, best_n = p, n
    return int(best_p)

@dataclass
class Decision:
    decision: str
    note: str
    any_left: int
    pri_left: int
    any_right: int
    pri_right: int
    pos_left: int
    pos_right: int
    w_left: int
    w_right: int

def decide_for_edge(edge_dir: str,
                    pass_any:int=1000, pass_pri:int=50,
                    border_any:int=200, border_pri:int=10,
                    default_W:int=100, step:int=10) -> Decision:
    reads_dir = os.path.join(edge_dir, "reads")
    junc_bams = [os.path.join(reads_dir, f) for f in os.listdir(reads_dir) if f.startswith("junc.") and f.endswith(".bam")]
    if not junc_bams:
        raise FileNotFoundError(f"No junc.*.bam under {reads_dir}")
    junc_bam = junc_bams[0]

    with pysam.AlignmentFile(junc_bam, "rb") as bam:
        # best POS per side
        posL = best_pos(bam, "JLEFT", W=default_W, step=step)
        posR = best_pos(bam, "JRIGHT", W=default_W, step=step)

        # choose W adaptively (shrink if too few reads could possibly span)
        def choose_W(ctg: str, start_W: int) -> int:
            W = start_W
            pct = pct_len_ge_2W(bam, ctg, W)
            if pct < 95.0:
                W = 50
                pct = pct_len_ge_2W(bam, ctg, W)
            if pct < 95.0:
                W = 25
            return W

        wL = choose_W("JLEFT", default_W)
        wR = choose_W("JRIGHT", default_W)

        # counts
        anyL = count_spanners(bam, "JLEFT",  posL, wL, include_supp=True)
        priL = count_spanners(bam, "JLEFT",  posL, wL, include_supp=False)
        anyR = count_spanners(bam, "JRIGHT", posR, wR, include_supp=True)
        priR = count_spanners(bam, "JRIGHT", posR, wR, include_supp=False)

    # decision logic
    if (anyL >= pass_any and priL >= pass_pri and anyR >= pass_any and priR >= pass_pri):
        dec = "ACCEPT"
    elif (anyL >= border_any and priL >= border_pri and anyR >= border_any and priR >= border_pri):
        dec = "BORDERLINE"
    else:
        dec = "REJECT"

    note = f"JLEFT(pos={posL},W={wL},any={anyL},pri={priL});JRIGHT(pos={posR},W={wR},any={anyR},pri={priR})"
    return Decision(dec, note, anyL, priL, anyR, priR, posL, posR, wL, wR)

def write_twoboundary_tsv(edge_dir: str, dec: Decision) -> None:
    base = os.path.basename(edge_dir)
    seqid, start, end = parse_edge_dir_name(base)
    out_tsv = os.path.join(edge_dir, "twoboundary.tsv")
    header = ["seqid","start","end","pos_left","pos_right","w_left","w_right",
              "twob_left","twob_left_primary","twob_right","twob_right_primary",
              "twob_decision","twob_note"]
    row = [seqid, start, end, dec.pos_left, dec.pos_right, dec.w_left, dec.w_right,
           dec.any_left, dec.pri_left, dec.any_right, dec.pri_right, dec.decision, dec.note]
    with open(out_tsv, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)
        w.writerow(row)

# ----------------------- CLI -----------------------

def main():
    ap = argparse.ArgumentParser(description="Two-boundary verifier (auto-POS, adaptive W, supp-aware).")
    ap.add_argument("--edge-dir", required=True, help="verify/edge_*/ directory")
    ap.add_argument("--pass-any", type=int, default=1000)
    ap.add_argument("--pass-pri", type=int, default=50)
    ap.add_argument("--border-any", type=int, default=200)
    ap.add_argument("--border-pri", type=int, default=10)
    ap.add_argument("--W", type=int, default=100, help="default half-window size")
    ap.add_argument("--step", type=int, default=10, help="sweep step for bestPOS")
    args = ap.parse_args()

    dec = decide_for_edge(args.edge_dir, args.pass_any, args.pass_pri, args.border_any, args.border_pri, args.W, args.step)
    write_twoboundary_tsv(args.edge_dir, dec)
    # Also print a one-line summary for logs
    base = os.path.basename(args.edge_dir)
    print(f"{base}\t{dec.decision}\t{dec.note}")

if __name__ == "__main__":
    main()
