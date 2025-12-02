#!/usr/bin/env python3
# Version: v0.1.0
import sys, time

if len(sys.argv) < 3:
    print("Usage: pretend_work.py <seconds> <outfile>", file=sys.stderr)
    sys.exit(2)
sec = float(sys.argv[1])
out = sys.argv[2]
time.sleep(sec)
with open(out, "w") as f:
    f.write("ok\n")
print(f"wrote {out}")
