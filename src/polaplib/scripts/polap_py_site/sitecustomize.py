#!/usr/bin/env python3
# Version: v0.1.0
# Auto-loaded by Python if this directory is on PYTHONPATH.  [oai_citation:2‡Python documentation](https://docs.python.org/3/library/site.html?utm_source=chatgpt.com)
import sys, traceback


def _polap_excepthook(etype, value, tb):
    last = traceback.extract_tb(tb)[-1] if tb else None
    where = f"{last.filename}:{last.lineno} in {last.name}" if last else "?:?"
    sys.stderr.write(f"PY-FAILED {where}: {value}\n")
    traceback.print_exception(etype, value, tb)


sys.excepthook = _polap_excepthook  #  [oai_citation:3‡Python documentation](https://docs.python.org/3/library/sys.html?utm_source=chatgpt.com)
