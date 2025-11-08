#!/usr/bin/env python3
# polaplib/scripts/test_crash.py
# Version: v0.1.0
#
# PURPOSE
#   Minimal Python crash helper:
#     • prints a one-line context to stderr,
#     • optionally writes to an output file (to prove pre-crash work),
#     • raises a nested exception so the traceback has depth,
#     • exits with your chosen non-zero code (default 7).
#
# USAGE
#   python3 polaplib/scripts/test_crash.py [--message STR] [--exit-code N] [--in FILE] [--out FILE]
#
# EXAMPLES
#   python3 polaplib/scripts/test_crash.py
#   python3 polaplib/scripts/test_crash.py --message "provoking" --exit-code 13 --in in.txt --out out.txt

import argparse
import sys
import traceback
from typing import Optional


# --- robust excepthook that exits with requested code -------------------------
def _excepthook(exc_type, exc, tb):
    # print full traceback to stderr
    traceback.print_exception(exc_type, exc, tb)
    code = getattr(exc, "code", None)
    try:
        code = int(code) if code is not None else 7
    except Exception:
        code = 7
    sys.exit(code)


sys.excepthook = _excepthook


# --- cli ----------------------------------------------------------------------
def _parse_args(argv: Optional[list[str]] = None):
    p = argparse.ArgumentParser(description="polap python crash test")
    p.add_argument(
        "--message",
        default="intentional test failure (python)",
        help="failure message (default: %(default)s)",
    )
    p.add_argument(
        "--exit-code",
        type=int,
        default=7,
        dest="exit_code",
        help="exit code to use on failure (default: %(default)s)",
    )
    p.add_argument("--in", dest="infile", default="", help="optional input file")
    p.add_argument(
        "--out",
        dest="outfile",
        default="",
        help="optional output file (written pre-crash)",
    )
    return p.parse_args(argv)


# --- nested functions to make a readable stack --------------------------------
class ExitCodeError(RuntimeError):
    def __init__(self, msg: str, code: int):
        super().__init__(msg)
        self.code = int(code)


def level3(msg: str, code: int):
    raise ExitCodeError(msg, code)


def level2(msg: str, code: int):
    level3(msg, code)


def level1(msg: str, code: int):
    level2(msg, code)


# --- main ---------------------------------------------------------------------
def main():
    args = _parse_args()

    # print a context line to stderr
    print(
        f"PY test_crash: message='{args.message}' code={args.exit_code} "
        f"infile='{args.infile}' outfile='{args.outfile}'",
        file=sys.stderr,
    )

    # prove we can do work before failing
    if args.outfile:
        with open(args.outfile, "w", encoding="utf-8") as fh:
            fh.write("test_crash.py wrote before failing\n")

    # intentionally fail with a nested stack
    level1(args.message, args.exit_code)


if __name__ == "__main__":
    main()
