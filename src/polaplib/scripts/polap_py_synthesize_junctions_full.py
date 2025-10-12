#!/usr/bin/env python3
# Version: v0.7.0
"""
Synthesize parental and recombinant junction templates for each repeat pair.

Inputs
------
--assembly    FASTA
--repeats     repeats.tsv  (pair_id r_name r_s r_e q_name q_s q_e orient len pid)
--flank       int  (bp per side to concatenate)
--out-parent-fasta / --out-parent-meta
--out-recomb-fasta / --out-recomb-meta

Outputs
-------
Parent FASTA/meta: two parental adjacencies per pair (PAR_A on r, PAR_B on q)
Recomb FASTA/meta: two recombinant adjacencies per pair:
  DIR: DR_del (r_left + q_right), DR_circ (q_left + r_right)
  INV: INV (r_left + q_left_rc), INV (r_right + q_right_rc)
Meta columns: junction_id pair_id type r_name r_pos q_name q_pos flank
"""
from __future__ import annotations
import argparse, logging
from typing import Dict, Tuple
from Bio import SeqIO
from Bio.Seq import Seq

__version__ = "0.7.0"


def grab(seq: Seq, s: int, e: int) -> Seq:
    s = max(1, s)
    e = min(len(seq), e)
    if e < s:
        return Seq("")
    return seq[s - 1 : e]


def main():
    p = argparse.ArgumentParser(
        description="Synthesize parental/recombinant junction templates."
    )
    p.add_argument("--assembly", required=True)
    p.add_argument("--repeats", required=True)
    p.add_argument("--flank", type=int, default=500)
    p.add_argument("--out-parent-fasta", required=True)
    p.add_argument("--out-parent-meta", required=True)
    p.add_argument("--out-recomb-fasta", required=True)
    p.add_argument("--out-recomb-meta", required=True)
    p.add_argument(
        "-v",
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
    )
    p.add_argument("--version", action="store_true")
    a = p.parse_args()
    if a.version:
        print(__version__)
        return
    logging.basicConfig(
        level=getattr(logging, a.log_level), format="%(levelname)s: %(message)s"
    )

    seqs: Dict[str, Seq] = {
        r.id: r.seq.upper() for r in SeqIO.parse(a.assembly, "fasta")
    }
    rows = []
    with open(a.repeats) as fh:
        fh.readline()
        for ln in fh:
            if not ln.strip():
                continue
            pair, r, rs, re, q, qs, qe, orient, _, _ = ln.rstrip("\n").split("\t")
            rows.append((pair, r, int(rs), int(re), q, int(qs), int(qe), orient))

    pf = open(a.out_parent_fasta, "w")
    pm = open(a.out_parent_meta, "w")
    rf = open(a.out_recomb_fasta, "w")
    rm = open(a.out_recomb_meta, "w")
    pm.write("junction_id\tpair_id\ttype\tr_name\tr_pos\tq_name\tq_pos\tflank\n")
    rm.write("junction_id\tpair_id\ttype\tr_name\tr_pos\tq_name\tq_pos\tflank\n")
    jidP = 0
    jidR = 0

    for pair, r, rs, re, q, qs, qe, orient in rows:
        rseq, qseq = seqs[r], seqs[q]
        rL = grab(rseq, min(rs, re) - a.flank, min(rs, re) - 1)
        rR = grab(rseq, max(rs, re), max(rs, re) + a.flank - 1)
        qL = grab(qseq, min(qs, qe) - a.flank, min(qs, qe) - 1)
        qR = grab(qseq, max(qs, qe), max(qs, qe) + a.flank - 1)

        # parentals
        jidP += 1
        pf.write(
            f">P{jidP}|PAR_A|{pair}|{r}:{min(rs,re)}|{r}:{max(rs,re)}|flank={a.flank}\n{(rL + rR)}\n"
        )
        pm.write(
            f"P{jidP}\t{pair}\tPAR_A\t{r}\t{r}:{min(rs,re)}\t{r}\t{r}:{max(rs,re)}\t{a.flank}\n"
        )
        jidP += 1
        pf.write(
            f">P{jidP}|PAR_B|{pair}|{q}:{min(qs,qe)}|{q}:{max(qs,qe)}|flank={a.flank}\n{(qL + qR)}\n"
        )
        pm.write(
            f"P{jidP}\t{pair}\tPAR_B\t{q}\t{q}:{min(qs,qe)}\t{q}\t{q}:{max(qs,qe)}\t{a.flank}\n"
        )

        # recombinants
        if orient == "DIR":
            jidR += 1
            rf.write(
                f">R{jidR}|DR_del|{pair}|{r}:{min(rs,re)}|{q}:{max(qs,qe)}|flank={a.flank}\n{(rL + qR)}\n"
            )
            rm.write(
                f"R{jidR}\t{pair}\tDR_del\t{r}\t{r}:{min(rs,re)}\t{q}\t{q}:{max(qs,qe)}\t{a.flank}\n"
            )
            jidR += 1
            rf.write(
                f">R{jidR}|DR_circ|{pair}|{q}:{min(qs,qe)}|{r}:{max(rs,re)}|flank={a.flank}\n{(qL + rR)}\n"
            )
            rm.write(
                f"R{jidR}\t{pair}\tDR_circ\t{q}\t{q}:{min(qs,qe)}\t{r}\t{r}:{max(rs,re)}\t{a.flank}\n"
            )
        else:
            jidR += 1
            rf.write(
                f">R{jidR}|INV|{pair}|{r}:{min(rs,re)}|{q}:{min(qs,qe)}|flank={a.flank}\n{(rL + qL.reverse_complement())}\n"
            )
            rm.write(
                f"R{jidR}\t{pair}\tINV\t{r}\t{r}:{min(rs,re)}\t{q}\t{q}:{min(qs,qe)}\t{a.flank}\n"
            )
            jidR += 1
            rf.write(
                f">R{jidR}|INV|{pair}|{r}:{max(rs,re)}|{q}:{max(qs,qe)}|flank={a.flank}\n{(rR + qR.reverse_complement())}\n"
            )
            rm.write(
                f"R{jidR}\t{pair}\tINV\t{r}\t{r}:{max(rs,re)}\t{q}\t{q}:{max(qs,qe)}\t{a.flank}\n"
            )

    for fh in (pf, pm, rf, rm):
        fh.close()


if __name__ == "__main__":
    main()
