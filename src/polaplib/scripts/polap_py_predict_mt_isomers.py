#!/usr/bin/env python3
# Version: v0.8.0
"""
Build recombination junction templates for each repeat pair.

Input:
  --assembly   FASTA
  --repeats    repeats.tsv (pair_id r_name r_s r_e q_name q_s q_e orient len pid)
  --flank      bp per side to concatenate

Output:
  --out-fasta  junctions.fasta      (>J|TYPE|PAIR|posA|posB|flank=N)
  --out-meta   junctions.tsv        (junction_id pair_id type r_name r_pos q_name q_pos flank)
  --out-gfa    minimal network.gfa  (optional sketch)
"""
from __future__ import annotations
import argparse
from Bio import SeqIO
from Bio.Seq import Seq


def grab(seq: Seq, s: int, e: int) -> Seq:
    s = max(1, s)
    e = min(len(seq), e)
    if e < s:
        return Seq("")
    return seq[s - 1 : e]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--assembly", required=True)
    ap.add_argument("--repeats", required=True)
    ap.add_argument("--flank", type=int, default=500)
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-meta", required=True)
    ap.add_argument("--out-gfa", required=True)
    a = ap.parse_args()

    seqs = {r.id: r.seq.upper() for r in SeqIO.parse(a.assembly, "fasta")}
    rows = []
    with open(a.repeats) as fh:
        fh.readline()
        for ln in fh:
            if not ln.strip():
                continue
            pair, r, rs, re, q, qs, qe, orient, _, _ = ln.rstrip("\n").split("\t")
            rows.append((pair, r, int(rs), int(re), q, int(qs), int(qe), orient))

    ofa = open(a.out_fasta, "w")
    omt = open(a.out_meta, "w")
    gfa = open(a.out_gfa, "w")
    omt.write("junction_id\tpair_id\ttype\tr_name\tr_pos\tq_name\tq_pos\tflank\n")
    gfa.write("H\tVN:Z:1.0\n")
    j = 0
    for pair, r, rs, re, q, qs, qe, orient in rows:
        rseq, qseq = seqs[r], seqs[q]
        rL = grab(rseq, min(rs, re) - a.flank, min(rs, re) - 1)
        rR = grab(rseq, max(rs, re), max(rs, re) + a.flank - 1)
        qL = grab(qseq, min(qs, qe) - a.flank, min(qs, qe) - 1)
        qR = grab(qseq, max(qs, qe), max(qs, qe) + a.flank - 1)
        if orient == "DIR":
            parts = [
                ("DR_del", rL + qR, f"{r}:{min(rs,re)}", f"{q}:{max(qs,qe)}"),
                ("DR_circ", qL + rR, f"{q}:{min(qs,qe)}", f"{r}:{max(rs,re)}"),
            ]
        else:
            parts = [
                (
                    "INV",
                    rL + qL.reverse_complement(),
                    f"{r}:{min(rs,re)}",
                    f"{q}:{min(qs,qe)}",
                ),
                (
                    "INV",
                    rR + qR.reverse_complement(),
                    f"{r}:{max(rs,re)}",
                    f"{q}:{max(qs,qe)}",
                ),
            ]
        for t, seq, apos, bpos in parts:
            if len(seq) == 0:
                continue
            j += 1
            ofa.write(f">J{j}|{t}|{pair}|{apos}|{bpos}|flank={a.flank}\n{seq}\n")
            omt.write(f"J{j}\t{pair}\t{t}\t{r}\t{apos}\t{q}\t{bpos}\t{a.flank}\n")
        # light GFA sketch
        gfa.write(f"S\t{pair}_A\t*\tLN:i:{abs(re-rs)+1}\n")
        gfa.write(f"S\t{pair}_B\t*\tLN:i:{abs(qe-qs)+1}\n")
        gfa.write(f"L\t{pair}_A\t+\t{pair}_B\t{('-' if orient=='INV' else '+')}\t0M\n")

    for fh in (ofa, omt, gfa):
        fh.close()


if __name__ == "__main__":
    main()
