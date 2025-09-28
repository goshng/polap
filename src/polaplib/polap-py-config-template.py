#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
polap-py-config-template.py  (flat YAML)

Write base YAML templates for HiFi and/or ONT presets into a config dir,
with **all items at the top level only** (no nested sections).

Defaults:
  - Config dir: ~/.polap/profiles
  - Files:      ~/.polap/profiles/hifi.yaml, ~/.polap/profiles/ont.yaml
  - Will NOT overwrite unless --force
  - Creates parent dirs (mkdir -p)

Usage:
  # write both templates (default)
  polap-py-config-template.py

  # only HiFi or only ONT
  polap-py-config-template.py --which hifi
  polap-py-config-template.py --which ont

  # custom config dir
  polap-py-config-template.py --config-dir /project/polap/profiles

  # also print to stdout (preview)
  polap-py-config-template.py --print

  # overwrite existing files
  polap-py-config-template.py --force
"""
import argparse, os, sys, datetime


def _ensure_dir(d: str):
    os.makedirs(d, exist_ok=True)


def _dump_yaml_flat(d: dict, fh):
    """
    Write a flat dict as YAML. Use PyYAML if available; fallback to a simple writer.
    """
    try:
        import yaml  # type: ignore

        yaml.safe_dump(d, fh, default_flow_style=False, sort_keys=False)
    except Exception:
        for k, v in d.items():
            if isinstance(v, bool):
                vv = "true" if v else "false"
            elif v is None:
                vv = "null"
            else:
                s = str(v)
                if any(c in s for c in [":", "#", "'", '"', " "]):
                    s = s.replace("'", "''")
                    vv = f"'{s}'"
                else:
                    vv = s
            fh.write(f"{k}: {vv}\n")


def _template_flat(preset: str) -> dict:
    """
    Return a **flat** YAML dict for a preset ('hifi' or 'ont').
    All keys are top-level.
    """
    ts = datetime.datetime.now(datetime.timezone.utc).isoformat()
    if preset == "ont":
        data = {
            "preset": "ont",
            "wgs_mode": True,
            "reads": "/path/to/reads.fastq[.gz]",
            "anchors": "/path/to/anchors.id.txt",
            "hmm_db": "/path/to/mt_or_pt_genes.hmm",
            "outdir": "/path/to/autotune_out",
            "threads": 16,
            "k": 41,
            "s": 21,
            "hpc": True,
            "final_if": "genes_score < 0.90 || breadth < 0.95",
            "c_radius": "0,10,20",
            "min_shared": 5,
            "jaccard_min": 0.0075,
            "topk_nei": 40,
            "steps": 1,
            "description": "ONT Q20+/duplex presets; HPC + smaller k/s; looser Jaccard",
            "normalized_at": ts,
        }
    else:
        data = {
            "preset": "hifi",
            "wgs_mode": True,
            "reads": "/path/to/reads.fastq[.gz]",
            "anchors": "/path/to/anchors.id.txt",
            "hmm_db": "/path/to/mt_or_pt_genes.hmm",
            "outdir": "/path/to/autotune_out",
            "threads": 16,
            "k": 121,
            "s": 27,
            "hpc": False,
            "final_if": "genes_score < 0.95 || breadth < 0.97",
            "c_radius": "0,10,20",
            "min_shared": 5,
            "jaccard_min": 0.015,
            "topk_nei": 50,
            "steps": 2,
            "description": "PacBio HiFi presets for organelle tuning",
            "normalized_at": ts,
        }
    return data


def write_template(path: str, data: dict, force: bool, also_print: bool) -> int:
    # mkdir -p
    _ensure_dir(os.path.dirname(os.path.abspath(path)))
    if os.path.exists(path):
        if not force:
            print(f"[skip] {path} exists (use --force to overwrite)", file=sys.stderr)
            return 0
        else:
            print(f"[overwrite] {path}", file=sys.stderr)

    try:
        with open(path, "w", encoding="utf-8") as fh:
            _dump_yaml_flat(data, fh)
        print(f"[write] {path}", file=sys.stderr)
        if also_print:
            _dump_yaml_flat(data, sys.stdout)
        return 0
    except Exception as e:
        print(f"[error] failed to write {path}: {e}", file=sys.stderr)
        return 4


def main():
    ap = argparse.ArgumentParser(
        description="Write flat YAML templates for HiFi/ONT presets."
    )
    ap.add_argument(
        "--config-dir",
        default=os.path.join(os.path.expanduser("~"), ".polap", "profiles"),
        help="target config directory (default: ~/.polap/profiles)",
    )
    ap.add_argument(
        "--which",
        choices=["hifi", "ont", "both"],
        default="both",
        help="which template(s) to write (default: both)",
    )
    ap.add_argument("--force", action="store_true", help="overwrite existing files")
    ap.add_argument(
        "--print",
        dest="do_print",
        action="store_true",
        help="also print YAML to stdout",
    )
    args = ap.parse_args()

    cfgdir = os.path.abspath(args.config_dir)
    _ensure_dir(cfgdir)

    todo = []
    if args.which in ("hifi", "both"):
        todo.append(("hifi", os.path.join(cfgdir, "hifi.yaml")))
    if args.which in ("ont", "both"):
        todo.append(("ont", os.path.join(cfgdir, "ont.yaml")))

    rc = 0
    for preset, outp in todo:
        data = _template_flat(preset)
        code = write_template(outp, data, force=args.force, also_print=args.do_print)
        rc = rc or code
    return rc


if __name__ == "__main__":
    sys.exit(main())
