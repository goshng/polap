#!/usr/bin/env python3
# Version: v0.1.0
# BED -> TSV with interval length (and contig length if .fai present)
import sys, argparse, os
def read_fai_lengths(fa):
    fai=fa+".fai"
    if not os.path.exists(fai): return {}
    L={}
    with open(fai) as fh:
        for ln in fh:
            f=ln.rstrip("\n").split("\t")
            L[f[0]]=int(f[1])
    return L
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--bed", required=True)
    ap.add_argument("--out", required=True)
    a=ap.parse_args()
    lens=read_fai_lengths(a.fasta)
    n=0; tot=0
    with open(a.bed) as fh, open(a.out,"w") as oh:
        oh.write("contig\tstart\tend\tlength_bp\tcontig_len\n")
        for ln in fh:
            if not ln.strip(): continue
            t,s,e=ln.rstrip("\n").split("\t")[:3]
            s=int(s); e=int(e); L=max(0,e-s); tot+=L; n+=1
            oh.write(f"{t}\t{s}\t{e}\t{L}\t{lens.get(t,'')}\n")
    sys.stderr.write(f"[bed_add_len] {n} intervals, total masked={tot} bp -> {a.out}\n")
if __name__=="__main__": main()

