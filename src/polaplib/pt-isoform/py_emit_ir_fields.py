#!/usr/bin/env python3
# py_emit_ir_fields.py  v0.1.0
import sys
p=sys.argv[1] if len(sys.argv)>1 else "-"
with open(p) as f:
    _=next(f,""); row=next(f,"").strip()
    if not row: exit(3)
    s=row.split(); print(s[0],s[1],s[2],s[3],s[4])
