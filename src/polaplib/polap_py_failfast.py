# polaplib/polap_py_failfast.py
# Make Python failures short & useful: print the failing frame, line text,
# a compact traceback, and exit with a code (default 7).

import os, sys, traceback
from datetime import datetime

def polap_py_enable_failsafe(default_code=7):
    """
    Install a global excepthook:
      - prints: [ts (conda) func@file:line] <source line>
      - shows a compact traceback (limit 10)
      - exits with exception.code if present else default_code
    """
    def excepthook(exc_type, exc, tb):
        ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        cenv = os.environ.get("CONDA_DEFAULT_ENV")
        envtag = f" ({cenv})" if cenv else ""
        hdr = f"[{ts}{envtag} PY-CRASH]"

        last = None
        if tb is not None:
            tbe = traceback.extract_tb(tb)
            if tbe:
                last = tbe[-1]

        if last:
            fname = os.path.basename(last.filename)
            sys.stderr.write(f"{hdr} {last.name}@{fname}:{last.lineno}\n")
            if last.line:
                sys.stderr.write(f"{hdr} line: {last.line.strip()}\n")

        traceback.print_exception(exc_type, exc, tb, limit=10, file=sys.stderr)

        code = getattr(exc, "code", default_code)
        try:
            code = int(code)
        except Exception:
            code = default_code
        sys.exit(code)

    sys.excepthook = excepthook

