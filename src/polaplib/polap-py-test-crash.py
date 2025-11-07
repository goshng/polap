#!/usr/bin/env python3
# polap-py-test-crash.py
# Version: v0.1.0
# Self-contained: no PYTHONPATH/sitecustomize needed; Python tracebacks show file:line.
import argparse


def leaf(msg: str, code: int) -> None:
    raise RuntimeError(f"demo crash (python): {msg} [code={code}]")


def mid(msg: str, code: int) -> None:
    return leaf(msg, code)


def top(msg: str, code: int) -> None:
    return mid(msg, code)


def main():
    ap = argparse.ArgumentParser(description="Trigger a controlled Python crash.")
    ap.add_argument("--msg", default="hello from polap")
    ap.add_argument("--code", type=int, default=1)
    args = ap.parse_args()
    top(args.msg, args.code)


if __name__ == "__main__":
    main()
