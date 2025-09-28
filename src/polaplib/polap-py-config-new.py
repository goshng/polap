#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
polap-py-config-new.py
Create or update a **flat** YAML config by overlaying user-provided key/values
over an optional template YAML.

File resolution:
  --path PATH.yaml                             # explicit
  OR
  --config-dir DIR --preset NAME               # -> DIR/NAME.yaml

CLI:
  python3 polap-py-config-new.py \
    --config-dir DIR --preset NAME \
    [--path PATH.yaml] \
    [--template BASE.yaml|NAME] \
    [--kv key1=value1 --kv key2=value2 ...] \
    [--bool bkey true|false] \
    [--wgs-mode|--no-wgs-mode|--no-wgs_mode] \
    [--hpc|--no-hpc]

Notes
- Result YAML is flat (top-level keys only).
- Values from flags (--wgs-mode/--hpc) override --kv/--bool/template.
- If PyYAML is not available, falls back to a safe minimal writer/reader.

Exit codes: 0 ok; 2 usage; 3 IO; 4 parse
"""
import argparse, os, sys


# ---------------- YAML I/O ----------------
def _load_yaml(path):
    if not path or not os.path.exists(path):
        return {}
    try:
        import yaml  # type: ignore

        with open(path, "r", encoding="utf-8") as fh:
            doc = yaml.safe_load(fh) or {}
        if not isinstance(doc, dict):
            raise ValueError("YAML is not a mapping")
        # keep only scalars (flat)
        return {str(k): doc[k] for k in doc}
    except ImportError:
        data = {}
        with open(path, "r", encoding="utf-8") as fh:
            for ln in fh:
                ln = ln.strip()
                if not ln or ln.startswith("#") or ":" not in ln:
                    continue
                k, v = ln.split(":", 1)
                data[k.strip()] = v.strip().strip("'").strip('"')
        return data


def _dump_yaml(data, path):
    try:
        import yaml  # type: ignore

        with open(path, "w", encoding="utf-8") as fh:
            yaml.safe_dump(data, fh, default_flow_style=False, sort_keys=False)
    except ImportError:
        with open(path, "w", encoding="utf-8") as fh:
            for k, v in data.items():
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


# ---------------- helpers ----------------
def _coerce_scalar(s):
    """basic type coercion for --kv."""
    vl = s.lower()
    if vl in ("true", "false"):
        return vl == "true"
    try:
        if "." in s:
            return float(s)
        return int(s)
    except ValueError:
        return s


def _resolve_template(template_arg, config_dir):
    """Return template path or None."""
    if not template_arg:
        return None
    cand1 = os.path.join(config_dir, f"{template_arg}.yaml")
    if os.path.exists(cand1):
        return cand1
    # else treat as path
    return template_arg


# ---------------- main ----------------
def main():
    ap = argparse.ArgumentParser(description="Create/update a flat YAML config.")
    # standard addressing
    ap.add_argument("--config-dir", required=True, help="base config directory")
    ap.add_argument(
        "--preset", required=True, help="profile name (used as <dir>/<preset>.yaml)"
    )
    ap.add_argument(
        "--path", help="explicit YAML path (overrides --config-dir/--preset)"
    )

    # template & overlay
    ap.add_argument(
        "--template",
        help="base YAML name or path (resolved as <config-dir>/<name>.yaml if exists)",
    )
    ap.add_argument("--kv", action="append", default=[], help="key=value (repeat)")
    ap.add_argument(
        "--bool",
        nargs=2,
        action="append",
        default=[],
        metavar=("KEY", "BOOL"),
        help="boolean key true|false (repeat)",
    )

    # boolean flags
    ap.add_argument(
        "--wgs-mode",
        dest="flag_wgs_mode",
        action="store_true",
        help="set wgs_mode=true",
    )
    ap.add_argument(
        "--no-wgs-mode",
        dest="flag_wgs_mode",
        action="store_false",
        help="set wgs_mode=false",
    )
    ap.add_argument(
        "--no-wgs_mode",
        dest="flag_wgs_mode",
        action="store_false",
        help="set wgs_mode=false",
    )
    ap.add_argument("--hpc", dest="flag_hpc", action="store_true", help="set hpc=true")
    ap.add_argument(
        "--no-hpc", dest="flag_hpc", action="store_false", help="set hpc=false"
    )

    args = ap.parse_args()

    # Compute target path
    target_path = args.path or os.path.join(args.config_dir, f"{args.preset}.yaml")
    template_path = _resolve_template(args.template, args.config_dir)

    try:
        cfg = {}
        # template overlay first
        if template_path:
            cfg.update(_load_yaml(template_path))

        # --kv overlays
        for item in args.kv:
            if "=" not in item:
                print(f"[error] --kv must be key=value, got {item}", file=sys.stderr)
                return 4
            k, v = item.split("=", 1)
            cfg[k.strip()] = _coerce_scalar(v.strip())

        # --bool overlays
        for k, bv in args.bool or []:
            val = str(bv).lower()
            if val not in ("true", "false"):
                print(
                    f"[error] --bool {k} must be true|false, got {bv}", file=sys.stderr
                )
                return 4
            cfg[k.strip()] = val == "true"

        # boolean flags take precedence
        if args.flag_wgs_mode is not None:
            cfg["wgs_mode"] = bool(args.flag_wgs_mode)
        if args.flag_hpc is not None:
            cfg["hpc"] = bool(args.flag_hpc)

        # ensure directory & save
        os.makedirs(os.path.dirname(os.path.abspath(target_path)), exist_ok=True)
        _dump_yaml(cfg, target_path)
        print(f"[save] {target_path}")
        return 0

    except FileNotFoundError as e:
        print(f"[error] {e}", file=sys.stderr)
        return 3
    except Exception as e:
        print(f"[error] {e}", file=sys.stderr)
        return 4


if __name__ == "__main__":
    sys.exit(main())
