#!/usr/bin/env python3
# Version: v0.1.0
# Name: polap-py-paf-filter.py
# Purpose: Filter minimap2 PAF by identity (nmatch/alen) and aligned length.
# Usage:
#   python3 polap-py-paf-filter.py <min_ident> <min_alen> < in.paf > out.paf
# Notes:
#   - Splits on ANY whitespace (tabs or spaces).
#   - Expects minimap2 PAF with columns: 10=nmatch, 11=alen (use minimap2 -c).
import sys


def main():
    if len(sys.argv) < 3:
        sys.stderr.write("usage: polap-py-paf-filter.py <min_ident> <min_alen>\n")
        sys.exit(2)
    try:
        min_ident = float(sys.argv[1])
        min_alen = int(sys.argv[2])
    except Exception as e:
        sys.stderr.write(f"[paf_filter] bad args: {e}\n")
        sys.exit(2)

    seen = kept = 0
    for ln in sys.stdin:
        if not ln.strip() or ln[0] == "#":
            continue
        seen += 1
        c = ln.strip().split()  # ANY whitespace
        if len(c) < 12:
            continue
        try:
            nm = int(c[9])  # nmatch
            al = int(c[10])  # aligned length
        except Exception:
            continue
        if al <= 0 or nm > al:
            continue
        if al >= min_alen and (nm / al) >= min_ident:
            sys.stdout.write(ln)
            kept += 1
    sys.stderr.write(f"[paf_filter] seen={seen} kept={kept}\n")


if __name__ == "__main__":
    main()
