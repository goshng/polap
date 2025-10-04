#!/usr/bin/env python3
# py_make_doubled_fa.py  v0.1.0
# Read first record from a FASTA, write a single header (>pt) with seq+seq.
import sys
inp, outp = sys.argv[1], sys.argv[2]
name, seq = None, []
with open(inp) as f:
  for ln in f:
    if ln.startswith('>'):
      if name is None:
        name = ln[1:].strip().split()[0]
      else:
        break
    else:
      seq.append(ln.strip())
S = ''.join(seq)
with open(outp, 'w') as w:
  w.write('>pt\n'); w.write(S+S+'\n')
