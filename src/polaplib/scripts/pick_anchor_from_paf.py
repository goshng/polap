#!/usr/bin/env python3
# scripts/pick_anchor_from_paf.py
# Version: v0.1.0
import argparse, sys

def parse_paf(path):
    hits=[]
    with open(path) as f:
        for line in f:
            if not line.strip(): continue
            a=line.rstrip("\n").split("\t")
            qname=a[0]; tname=a[5]
            strand=a[4]  # '+' or '-'
            tstart=int(a[7])+1  # PAF: 0-based, half-open -> 1-based inclusive
            tend=int(a[8])
            mapq=int(a[11]) if len(a)>11 else 0
            alen=tend - (tstart-1)
            # optional tags: look for 'cm:i:' (matches) and 'dv:f:' (divergence)
            score=mapq
            hits.append(dict(src="paf",q=qname,t=tname,ts=tstart,te=tend,strand=strand,mapq=mapq,alen=alen,score=score))
    return hits

def parse_blast_tsv(path):
    hits=[]
    with open(path) as f:
        for line in f:
            if not line.strip(): continue
            qid, tid, sstart, send, sstrand, bits, alen, pident = line.rstrip("\n").split("\t")
            sstart=int(sstart); send=int(send); bits=float(bits); alen=int(alen); pident=float(pident)
            if sstrand in ("plus","Plus"): strand="+"
            elif sstrand in ("minus","Minus"): strand="-"
            else: strand="+"
            # normalize to 1-based increasing
            ts=min(sstart, send); te=max(sstart, send)
            hits.append(dict(src="blast",q=qid,t=tid,ts=ts,te=te,strand=strand,mapq=int(bits),alen=alen,score=float(bits)))
    return hits

ap=argparse.ArgumentParser()
ap.add_argument("--paf")
ap.add_argument("--blast-tsv")
ap.add_argument("--prefer-plus", type=int, default=0)
ap.add_argument("--report-all", type=int, default=0)
ap.add_argument("--out-tsv", required=True)
args=ap.parse_args()

if not args.paf and not args.blast_tsv:
    sys.exit("provide --paf or --blast-tsv")

hits = parse_paf(args.paf) if args.paf else parse_blast_tsv(args.blast_tsv)
if not hits:
    print("", end="")
    sys.exit(0)

# rank
def key(h):
    pref_plus = 1 if (args.prefer_plus and h["strand"]=="+") else 0
    return (h["score"], h["alen"], pref_plus, -h["ts"])

hits.sort(key=key, reverse=True)
best=hits[0]

# write all hits (optional)
with open(args.out_tsv,"w") as o:
    o.write("qname\ttname\ttstart\ttend\tstrand\tscore\talen\n")
    for h in hits:
        o.write(f"{h['q']}\t{h['t']}\t{h['ts']}\t{h['te']}\t{h['strand']}\t{h['score']}\t{h['alen']}\n")

# print the winner to STDOUT for bash caller (single line)
print(f"{best['q']}\t{best['t']}\t{best['ts']}\t{best['te']}\t{best['strand']}\t{best['score']}")
