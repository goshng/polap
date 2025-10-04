#!/usr/bin/env python3
# py_rotate_fa_to_start.py  v0.1.0
# Usage: py_rotate_fa_to_start.py in.fa start out.fa
import sys
if len(sys.argv)<4: sys.exit(2)
inp, start, outp = sys.argv[1], int(sys.argv[2]), sys.argv[3]
name=None; seq=[]
with open(inp) as f:
    for ln in f:
        if ln.startswith('>'):
            if name is None: name=ln[1:].strip().split()[0]
            else: break
        else: seq.append(ln.strip())
S=''.join(seq)
L=len(S); start=((start-1)%L)+1
rot=S[start-1:]+S[:start-1]
with open(outp,'w') as w:
    w.write('>pt.rot\n'); w.write(rot+'\n')
