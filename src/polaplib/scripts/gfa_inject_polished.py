#!/usr/bin/env python3
import sys, gzip
gfa, polished_fa = sys.argv[1], sys.argv[2]
def open_auto(p): return gzip.open(p,'rt') if p.endswith('.gz') else open(p)
seq={}
with open_auto(polished_fa) as f:
    cur=None; buf=[]
    for line in f:
        if line.startswith('>'):
            if cur: seq[cur]=''.join(buf); buf=[]
            cur=line[1:].strip().split()[0]
        else: buf.append(line.strip())
    if cur: seq[cur]=''.join(buf)
with open_auto(gfa) as f:
    for line in f:
        if not line.strip(): continue
        if line[0]=='S':
            parts=line.rstrip('\n').split('\t')
            sid=parts[1]
            if sid in seq: parts[2]=seq[sid]
            print('\t'.join(parts))
        else:
            print(line, end='')
