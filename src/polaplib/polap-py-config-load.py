#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
polap-py-config-load.py

Load a flat YAML config into a Python dict or emit env/json for shell consumption.

CLI:
  python3 polap-py-config-load.py --path FILE.yaml [--format env|json] [--prefix PCFG]
  python3 polap-py-config-load.py --config-dir DIR --preset NAME [--format env|json] [--prefix PCFG]

Notes
- Flat YAML (top-level keys) is expected.
- --prefix defaults to 'PCFG' so variables look like PCFG_PRESET, PCFG_READS, ...
"""
import sys, os, argparse, json


def _load_yaml(path):
    try:
        import yaml  # type: ignore

        with open(path, "r", encoding="utf-8") as fh:
            d = yaml.safe_load(fh) or {}
        if not isinstance(d, dict):
            raise ValueError("YAML is not a mapping")
        return d
    except ImportError:
        # fallback minimal parser (flat: key: value)
        d = {}
        with open(path, "r", encoding="utf-8") as fh:
            for ln in fh:
                ln = ln.strip()
                if not ln or ln.startswith("#") or ":" not in ln:
                    continue
                k, v = ln.split(":", 1)
                k = k.strip()
                v = v.strip().strip("'").strip('"')
                # cast simple booleans and numbers
                vl = v.lower()
                if vl in ("true", "false"):
                    d[k] = vl == "true"
                else:
                    try:
                        if "." in v:
                            d[k] = float(v)
                        else:
                            d[k] = int(v)
                    except ValueError:
                        d[k] = v
        return d


def load_config(path):
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    return _load_yaml(path)


def resolve_config(config_dir=None, preset=None, path=None):
    if path:
        return os.path.expanduser(path)
    if not config_dir or not preset:
        raise ValueError("Need --path or both --config-dir and --preset")
    return os.path.join(os.path.expanduser(config_dir), f"{preset}.yaml")


def _emit_env(d, prefix="PCFG"):
    # prefix without underscore; we add the underscore in the formatted var
    lines = []
    pfx = (prefix or "PCFG").upper()
    for k, v in d.items():
        key = pfx + "_" + str(k).upper().replace("-", "_")
        if isinstance(v, bool):
            val = "true" if v else "false"
            lines.append(f"{key}={val}")
        elif v is None:
            lines.append(f"{key}=")
        else:
            s = str(v).replace("'", "'\"'\"'")
            lines.append(f"{key}='{s}'")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description="Load a YAML config and emit env/json.")
    ap.add_argument("--path", help="explicit YAML file")
    ap.add_argument("--config-dir", help="base config dir (default: ~/.polap/profiles)")
    ap.add_argument("--preset", help="profile name (used with --config-dir)")
    ap.add_argument("--format", choices=["env", "json"], default="env")
    ap.add_argument(
        "--prefix",
        default="PCFG",
        help="env var prefix for --format env (default: PCFG)",
    )
    args = ap.parse_args()

    if not args.path and not (args.config_dir and args.preset):
        print("[error] provide --path or (--config-dir and --preset)", file=sys.stderr)
        return 2

    if not args.path and not args.config_dir:
        args.config_dir = os.path.join(os.path.expanduser("~"), ".polap", "profiles")

    try:
        cfg_path = resolve_config(args.config_dir, args.preset, args.path)
        cfg = load_config(cfg_path)
        if args.format == "json":
            print(json.dumps(cfg, ensure_ascii=False))
        else:
            print(_emit_env(cfg, prefix=args.prefix))
        return 0
    except Exception as e:
        print(f"[error] {e}", file=sys.stderr)
        return 4


if __name__ == "__main__":
    sys.exit(main())
