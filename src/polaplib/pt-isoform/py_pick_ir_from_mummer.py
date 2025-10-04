#!/usr/bin/env python3
# py_pick_ir_from_mummer.py  v0.2.1 (wrap-aware, origin-safe pairing)
import sys, argparse
def parse_coords(path):
    out=[]
    with open(path) as f:
        for ln in f:
            ln=ln.strip()
            if not ln: continue
            xs=ln.split()
            if len(xs) < 12: continue
            try:
                s1=int(xs[0]); e1=int(xs[1])
                s2=int(xs[2]); e2=int(xs[3])
                l1=int(xs[4]); l2=int(xs[5]); idy=float(xs[6])
            except ValueError:
                continue
            qn=xs[-2]; tn=xs[-1]
            out.append((s1,e1,s2,e2,l1,l2,idy,qn,tn))
    return out
def is_inv(s1,e1,s2,e2): return (s1<=e1 and s2>=e2) or (s1>=e1 and s2<=e2)
def near_diag(s1,e1,s2,e2,p): return abs(s1-s2)<p and abs(e1-e2)<p
def within(x,a,b): return a <= x <= b
def norm1(x,L): t=(x-1)%L; return t+1
def make_variant(a1,a2,L):
    if a1<=a2: s_d,e_d=a1,a2
    else:      s_d,e_d=a2,a1
    len_d=e_d-s_d+1
    len_w=L-len_d
    s_w=norm1(e_d-(len_w-1),L)
    e_w=norm1(s_w+(len_w-1),L)
    return (s_d,e_d,len_d),(s_w,e_w,len_w)
def circ_overlap(sA,lA,sB,lB,L):
    a0=(sA-1)%L; b0=(sB-1)%L
    def segs(s,l):
        e=(s+l-1)%L
        return [(s,e)] if s<=e else [(s,L-1),(0,e)]
    A=segs(a0,lA); B=segs(b0,lB)
    ov=0
    for x0,x1 in A:
        for y0,y1 in B:
            ov += max(0, min(x1,y1)-max(x0,y0)+1)
    return ov
def pick_pair(rows,L,minL,maxL,reld,diag):
    C=[]
    for s1,e1,s2,e2,l1,l2,idy,qn,tn in rows:
        if qn!=tn: continue
        if not is_inv(s1,e1,s2,e2): continue
        if near_diag(s1,e1,s2,e2,diag): continue
        A_d,A_w=make_variant(s1,e1,L)
        B_d,B_w=make_variant(s2,e2,L)
        C.append((A_d,A_w,B_d,B_w))
    best=None; best_score=(-1,-1)
    for A_d,A_w,B_d,B_w in C:
        for Ad in (A_d,A_w):
            sA,eA,lA=Ad
            if not within(lA,minL,maxL): continue
            for Bd in (B_d,B_w):
                sB,eB,lB=Bd
                if not within(lB,minL,maxL): continue
                mu=(lA+lB)/2.0
                if mu==0: continue
                if abs(lA-lB)/mu>reld: continue
                ov=circ_overlap(sA,lA,sB,lB,L)
                if ov>1000: continue
                score=((lA+lB)/2.0,-ov)
                if score>best_score:
                    best_score=score; best=(sA,eA,lA,sB,eB,lB)
    return best
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("coords"); ap.add_argument("--L",type=int,required=True)
    ap.add_argument("--min_ir_len",type=int,default=5000)
    ap.add_argument("--max_ir_len",type=int,default=80000)
    ap.add_argument("--len_rel_diff",type=float,default=0.20)
    ap.add_argument("--near_diag_pad",type=int,default=200)
    ap.add_argument("--out",required=True)
    a=ap.parse_args()
    rows=parse_coords(a.coords)
    best=pick_pair(rows,a.L,a.min_ir_len,a.max_ir_len,a.len_rel_diff,a.near_diag_pad)
    with open(a.out,'w') as w:
        w.write("s1\te1\ts2\te2\tlen\n")
        if not best: return
        sA,eA,lA,sB,eB,lB=best
        # emit with A first by start
        if sA<=sB:
            w.write(f"{sA}\t{eA}\t{sB}\t{eB}\t{int((lA+lB)/2.0)}\n")
        else:
            w.write(f"{sB}\t{eB}\t{sA}\t{eA}\t{int((lA+lB)/2.0)}\n")
if __name__=="__main__":
    main()
