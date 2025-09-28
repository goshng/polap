#!/usr/bin/env python3
# polap-py-stop-delta.py  v0.0.1
import sys
delta = float(sys.argv[1]); thr = float(sys.argv[2])
print(1 if delta < thr else 0)
