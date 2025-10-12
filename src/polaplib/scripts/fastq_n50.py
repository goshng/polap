#!/usr/bin/env python3
import sys,gzip
from statistics import median
def openx(p): return gzip.open(p,'rt') if p.endswith('.gz') else open(p)
L=[]
with openx(sys.argv[1]) as f:
    i=0
    for line in f:
        i+=1
        if i%4==2: L.append(len(line.strip()))
L.sort()
tot=sum(L); half=tot/2; acc=0; n50=0
for l in L:
    acc+=l
    if acc>=half: n50=l; break
print(f"n={len(L)} total={tot} n50={n50} median={median(L) if L else 0}")
