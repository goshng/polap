#!/usr/bin/env python3
# Version: v0.1.0
# PAF -> TSV of hits with identity and lengths
import sys, argparse
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--paf", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    n=0
    with open(a.paf) as fh, open(a.out,"w") as oh:
        oh.write("qname\tqlen\tqstart\tqend\tstrand\ttname\ttlen\ttstart\ttend\tnmatch\talen\tidentity\n")
        for ln in fh:
            if not ln.strip() or ln[0]=='#': continue
            f=ln.rstrip("\n").split("\t")
            if len(f)<12: continue
            qn,ql,qs,qe,strand,tn,tl,ts,te,nm,al=f[0],int(f[1]),int(f[2]),int(f[3]),f[4],f[5],int(f[6]),int(f[7]),int(f[8]),int(f[9]),int(f[10])
            ident = (nm/ al) if al>0 else 0.0
            oh.write(f"{qn}\t{ql}\t{qs}\t{qe}\t{strand}\t{tn}\t{tl}\t{ts}\t{te}\t{nm}\t{al}\t{ident:.4f}\n")
            n+=1
    sys.stderr.write(f"[paf_to_hits] wrote {n} rows to {a.out}\n")
if __name__=="__main__": main()

