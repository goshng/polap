#!/usr/bin/env python3
# polaplib/polap-py-test-crash.py
# Version: v0.2.1
#
# Minimal Python crash helper:
#  - Prints a one-line context (stderr)
#  - Optionally writes to an output file before failing
#  - Fails via nested calls to produce a meaningful traceback
#  - Exits with the provided non-zero code
#  - --no-trace suppresses the Python traceback (bash failsafe will still report the call site)

import argparse
import sys
import traceback
from typing import Optional


class ExitCodeError(RuntimeError):
    """Carry an exit code through an exception path."""

    def __init__(self, msg: str, code: int):
        super().__init__(msg)
        self.code = int(code)


def _excepthook_with_trace(exc_type, exc, tb):
    """Print full traceback to stderr and exit with code."""
    traceback.print_exception(exc_type, exc, tb)
    code = getattr(exc, "code", 1)
    try:
        code = int(code)
    except Exception:
        code = 1
    sys.exit(code if code != 0 else 1)


def _excepthook_quiet(exc_type, exc, tb):
    """Quiet: print message only, exit with code."""
    msg = str(exc) or exc_type.__name__
    print(msg, file=sys.stderr)
    code = getattr(exc, "code", 1)
    try:
        code = int(code)
    except Exception:
        code = 1
    sys.exit(code if code != 0 else 1)


def parse_args(argv: Optional[list[str]] = None):
    p = argparse.ArgumentParser(
        description="polap crash helper (Python): prints context, optionally writes a file, then fails with a controlled exit code."
    )
    p.add_argument(
        "--message",
        "-m",
        default="intentional test failure (python)",
        help="Failure message (default: %(default)s)",
    )
    p.add_argument(
        "--code",
        "-c",
        type=int,
        default=7,
        help="Non-zero exit code (default: %(default)s)",
    )
    p.add_argument("--exit-code", dest="code_alt", type=int, help="Alias for --code")
    p.add_argument("--in", dest="infile", default="", help="Optional input filename")
    p.add_argument(
        "--out",
        dest="outfile",
        default="",
        help="Optional output filename (written pre-crash)",
    )
    p.add_argument(
        "--no-trace",
        dest="no_trace",
        action="store_true",
        help="Suppress Python traceback (bash failsafe will still report the call site)",
    )
    return p.parse_args(argv)


def level3(msg: str, code: int):
    raise ExitCodeError(msg, code)


def level2(msg: str, code: int):
    level3(msg, code)


def level1(msg: str, code: int):
    level2(msg, code)


def main():
    args = parse_args()

    # normalize code
    code = args.code_alt if args.code_alt is not None else args.code
    if not isinstance(code, int) or code == 0:
        code = 7

    # choose excepthook
    sys.excepthook = _excepthook_quiet if args.no_trace else _excepthook_with_trace

    # Context line (stderr) so bash failsafe can capture/expand the parent call
    print(
        f"PY test_crash: message='{args.message}' code={code} infile='{args.infile}' outfile='{args.outfile}'",
        file=sys.stderr,
        flush=True,
    )

    # Optional side-effect
    if args.outfile:
        try:
            with open(args.outfile, "w", encoding="utf-8") as fh:
                fh.write("polap-py-test-crash.py wrote before failing\n")
        except Exception as e:
            print(
                f"[child] warning: failed to write outfile '{args.outfile}': {e}",
                file=sys.stderr,
            )

    # Trigger nested failure
    level1(args.message, code)


if __name__ == "__main__":
    main()
