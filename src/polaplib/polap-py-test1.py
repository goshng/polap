#!/usr/bin/env python3
# Version: v0.1.0
# Python test: import the Python stack logger and raise an error from a function.

import os, sys
from importlib.machinery import SourceFileLoader
from pathlib import Path

BASE = (
    os.environ.get("POLAPLIB_DIR")
    or os.environ.get("_POLAPLIB_DIR")
    or str(Path(__file__).resolve().parent)
)
LOGPY = os.path.join(BASE, "polap-lib-logcallstack.py")
polaplog = SourceFileLoader("polaplog", LOGPY).load_module()  # import file with a dash

log = polaplog.get_logger("INFO")


def inner_py():
    polaplog.error_stack("Python: simulated failure inside inner_py()")
    raise RuntimeError("boom from Python")


def outer_py():
    log.info("Python: entering outer_py()")
    inner_py()


if __name__ == "__main__":
    try:
        outer_py()
    except Exception as e:
        # Log full stack (call frames + exception trace) then exit non-zero
        polaplog.error_stack("Python: unhandled exception bubbled to __main__", exc=e)
        sys.exit(2)
