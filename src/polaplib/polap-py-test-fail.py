#!/usr/bin/env python3
# polap-py-test-fail.py
# Version: v0.1.0
#
# Minimal crash-only script for testing bash callstack reporters.
# Usage:
#   python3 polap-py-test-fail.py [--message STR] [--exit-code N] [INFILE [OUTFILE]]

import argparse
import sys
import traceback


def parse_args():
    p = argparse.ArgumentParser(add_help=True)
    p.add_argument("--message", default="intentional test failure (python)")
    p.add_argument("--exit-code", type=int, default=7, dest="exit_code")
    p.add_argument("infile", nargs="?")
    p.add_argument("outfile", nargs="?")
    return p.parse_args()


class ExitCodeError(RuntimeError):
    def __init__(self, msg, code):
        super().__init__(msg)
        self.code = int(code)


def main():
    args = parse_args()

    # Optional write to show positional IO got through
    if args.outfile:
        with open(args.outfile, "w") as f:
            f.write("polap-py-test-fail.py wrote this before failing\n")

    # Ensure exit code follows the exception (unhandled)
    def _excepthook(exc_type, exc, tb):
        traceback.print_exception(exc_type, exc, tb)
        code = getattr(exc, "code", args.exit_code)
        sys.exit(code)

    sys.excepthook = _excepthook

    # Emit a small context line (stderr so it’s visible even when stdout is piped)
    print(
        f"PY test_fail: message='{args.message}' code={args.exit_code} "
        f"infile='{args.infile or ''}' outfile='{args.outfile or ''}'",
        file=sys.stderr,
    )

    # Create a few frames so the traceback has depth
    def level3():
        raise ExitCodeError(args.message, args.exit_code)

    def level2():
        level3()

    def level1():
        level2()

    level1()


if __name__ == "__main__":
    main()
