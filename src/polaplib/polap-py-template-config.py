#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
polap-py-template.py

Template: load a config and use values in Python.
"""
import sys, os, argparse

# Make sure we can import the loader from ${_POLAPLIB_DIR}
_POLAPLIB_DIR = os.environ.get("_POLAPLIB_DIR")
if _POLAPLIB_DIR and _POLAPLIB_DIR not in sys.path:
    sys.path.insert(0, _POLAPLIB_DIR)

from polap_py_config_load import load_config, resolve_config  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description="Template: load config and do work.")
    ap.add_argument("--path", help="explicit config path")
    ap.add_argument("--config-dir", help="base config dir (default: ~/.polap/profiles)")
    ap.add_argument("--preset", help="profile name (used with --config-dir)")
    args = ap.parse_args()

    # default config dir if needed
    if not args.path and not args.config_dir:
        args.config_dir = os.path.join(os.path.expanduser("~"), ".polap", "profiles")

    cfg_path = resolve_config(args.config_dir, args.preset, args.path)
    cfg = load_config(cfg_path)

    # Example: read fields
    preset = cfg.get("preset", "hifi")
    reads = cfg.get("reads")
    threads = int(cfg.get("threads", 8))
    k = int(cfg.get("k", 121))
    s = int(cfg.get("s", 27))
    hpc = bool(cfg.get("hpc", False))

    print(f"[config] path={cfg_path}")
    print(f"  preset={preset}  threads={threads}  k={k} s={s} hpc={hpc}")
    print(f"  reads={reads}")

    # TODO: put your application logic here using cfg[...] values


if __name__ == "__main__":
    sys.exit(main())
