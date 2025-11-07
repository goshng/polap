#!/usr/bin/env python3
# ##############################################################################
# This file is part of polap. GPLv3+ (see top-level license)
# ##############################################################################
# test_fail.py
# Intentional Python failure with full traceback and exit code control.
# Usage: test_fail.py [--message STR] [--exit-code N] [INFILE [OUTFILE]]

import argparse, sys, os, traceback

parser = argparse.ArgumentParser(add_help=True)
parser.add_argument("--message", default="intentional test failure (python)")
parser.add_argument("--exit-code", type=int, default=3, dest="exit_code")
parser.add_argument("infile", nargs="?")
parser.add_argument("outfile", nargs="?")
args = parser.parse_args()

# Optional write
if args.outfile:
    with open(args.outfile, "w") as f:
        f.write("test_fail.py writing before failure\n")


class ExitCodeError(RuntimeError):
    def __init__(self, message, code):
        super().__init__(message)
        self.code = int(code)


# Ensure exit code follows the exception
def _excepthook(exc_type, exc, tb):
    traceback.print_exception(exc_type, exc, tb)
    code = getattr(exc, "code", args.exit_code)
    sys.exit(code)


sys.excepthook = _excepthook


def level3():
    raise ExitCodeError(args.message, args.exit_code)


def level2():
    level3()


def level1():
    level2()


# Emit a little context and then fail
print(
    f"PY test_fail: message='{args.message}' code={args.exit_code} infile='{args.infile or ''}' outfile='{args.outfile or ''}'",
    file=sys.stderr,
)
level1()
