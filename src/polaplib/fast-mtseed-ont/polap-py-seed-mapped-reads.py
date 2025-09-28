#!/usr/bin/env python3
# polap-py-seed-mapped-reads.py  v0.0.1
# Read PAF from stdin; emit unique query IDs (after upstream filtering).
import sys
ids = set()
for ln in sys.stdin:
    if not ln.strip() or ln[0] == '#': continue
    p = ln.split('\t')
    if len(p) >= 1: ids.add(p[0])
for x in sorted(ids): print(x)
