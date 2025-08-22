# Gcov.py
# Usage: python Gcov.py total_bases.txt pivot_nuclear.txt fracNuc
import sys

if len(sys.argv) != 4:
    sys.exit("Usage: python Gcov.py total_bases.txt pivot_nuclear.txt fracNuc")
Tb = float(open(sys.argv[1]).read().strip() or 0)
Dn = float(open(sys.argv[2]).read().strip() or 0)
fracNuc = float(sys.argv[3])
if Dn <= 0 or Tb <= 0:
    print(0)
else:
    G = (fracNuc * Tb) / Dn
    print(f"{G:.3f}")
