#!/usr/bin/env python3
# Version: v0.3.0
"""
Uniform log lines + stack trace for Python.
File name has a dash; import via SourceFileLoader:

  import os
  from importlib.machinery import SourceFileLoader
  POLAPLIB_DIR = os.environ.get("POLAPLIB_DIR",".")
  polaplog = SourceFileLoader("polaplog", f"{POLAPLIB_DIR}/polap-lib-logcallstack.py").load_module()
  log = polaplog.get_logger("INFO")
  polaplog.error_stack("message")

Or run directly to self-test:
  python3 polap-lib-logcallstack.py
"""
from __future__ import annotations
import logging, os, sys, time, traceback

_FMT = "[%(asctime)s %(funcName)s@%(filename)s:%(lineno)d] %(message)s"
_DATE = "%Y-%m-%d %H:%M:%S"


def get_logger(level: str = "INFO") -> logging.Logger:
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format=_FMT,
        datefmt=_DATE,
        force=True,
    )
    return logging.getLogger("polap")


def _now() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())


def error_stack(msg: str, *, exc: BaseException | None = None) -> None:
    """
    Emit one ERROR line, then print a call-stack in the same house style.
    If 'exc' is provided, also print its traceback lines (prefixed).
    """
    logging.getLogger("polap").error(msg)
    # Current stack (exclude this frame)
    stk = traceback.extract_stack()[:-1]
    # Avoid printing frames inside this log helper
    me = os.path.basename(__file__)
    frames = [f for f in stk if os.path.basename(f.filename) != me]
    for idx, f in enumerate(reversed(frames)):
        fname = os.path.basename(f.filename)
        func = f.name or "main"
        sys.stderr.write(f"[{_now()} {func}@{fname}:{f.lineno}] #{idx}\n")

    if exc is not None:
        tb = "".join(traceback.format_exception(type(exc), exc, exc.__traceback__))
        for line in tb.rstrip("\n").splitlines():
            sys.stderr.write(f"[{_now()} stack@exception:0] {line}\n")


def log_exceptions(func):
    """
    Decorator: on exception, prints error_stack and re-raises.
    """

    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except BaseException as e:
            error_stack(f"Unhandled exception in {func.__name__}", exc=e)
            raise

    return wrapper


# Self-test
if __name__ == "__main__":
    log = get_logger("INFO")

    def a():
        b()

    def b():
        error_stack("Demo error_stack from Python")
        raise SystemExit(0)

    a()
