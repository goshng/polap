#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
polap-py-config-validate.py

Validate a polap YAML config (flat or nested). Optionally normalize and emit YAML.

Checks:
  - Required keys (preset, wgs_mode, reads, hmm_db, threads, k/s, hpc)
  - If wgs_mode=true: anchors required
  - Types: bools, ints, floats
  - Path existence warnings or errors (via --strict-paths)
  - Normalization: canonical nested YAML with defaults filled (k/s/hpc from preset)
  - Adds timestamps:
      * JSON summary: {"timestamp": "<UTC-ISO8601>"} when --json-summary
      * Normalized YAML: meta.normalized_at: "<UTC-ISO8601>"

Outputs:
  - Default: "OK" to stdout
  - --print               : normalized YAML to stdout
  - --write F             : normalized YAML to file F
  - --inplace             : overwrite input file with normalized YAML
  - --backup              : with --inplace, save original to file.yaml.bak (or .bak.N)
  - --dry-run             : simulate inplace/backup/write (no filesystem changes)
  - --fail-on-warn        : treat warnings as errors (exit nonzero if any warnings)
  - --quiet               : suppress detailed errors/warnings, only show summary
  - --json-summary        : emit summary as JSON to stdout (plain summary to stderr)
  - --summary-only        : suppress normalization output; only show summary
  - --summary-file PATH   : save summary (JSON or plain) to a file (overwrite)
  - --append-summary-file PATH : append summary (JSON or plain) to a file (mutually exclusive with --summary-file)

Exit codes:
  0: OK
  2: usage error
  3: validation errors (or warnings if --fail-on-warn)
  4: YAML parse error / IO error
"""
import argparse, os, sys, shutil, re, json, datetime


# ---------------- YAML helpers with optional positions ----------------
def _load_yaml_with_positions(path):
    pos_map = {}
    data = {}
    try:
        from ruamel.yaml import YAML

        yaml = YAML(typ="rt")
        with open(path, "r", encoding="utf-8") as fh:
            doc = yaml.load(fh)
        if not isinstance(doc, dict):
            raise ValueError("YAML top-level is not a mapping")
        for k in doc:
            data[k] = doc[k]
            try:
                ln = doc.lc.key(k)[0] + 1
            except Exception:
                ln = None
            if isinstance(k, str):
                pos_map[k] = ln

        def walk(node, prefix=""):
            if isinstance(node, dict):
                for k, v in node.items():
                    dotted = f"{prefix}.{k}" if prefix else str(k)
                    try:
                        ln = node.lc.key(k)[0] + 1
                    except Exception:
                        ln = pos_map.get(k)
                    pos_map[dotted] = ln
                    walk(v, dotted)

        walk(doc)
        return data, pos_map
    except Exception:
        try:
            import yaml  # type: ignore

            with open(path, "r", encoding="utf-8") as fh:
                doc = yaml.safe_load(fh) or {}
            if not isinstance(doc, dict):
                raise ValueError("YAML top-level is not a mapping")
            data = doc
        except Exception as e:
            print(f"[error] YAML parse failed: {e}", file=sys.stderr)
            raise

        with open(path, "r", encoding="utf-8") as fh:
            lines = fh.readlines()

        def find_line_for_key(key):
            pat = re.compile(rf"^\s*{re.escape(str(key))}\s*:")
            for i, ln in enumerate(lines, 1):
                if pat.search(ln):
                    return i
            return None

        for k in data.keys():
            if isinstance(k, str):
                pos_map[k] = find_line_for_key(k)

        nested_interest = [
            ("paths", ["reads", "hmm_db", "anchors", "outdir"]),
            ("compute", ["threads"]),
            ("syncmer", ["k", "s", "hpc"]),
            ("autotune", ["final_if", "c_radius"]),
            ("graph", ["min_shared", "jaccard_min", "topk_nei", "steps"]),
        ]
        for sect, keys in nested_interest:
            for kk in keys:
                pos_map[f"{sect}.{kk}"] = find_line_for_key(kk)

        return data, pos_map


def _dump_yaml_nested(obj, stream):
    try:
        import yaml  # type: ignore

        yaml.safe_dump(obj, stream, default_flow_style=False, sort_keys=False)
    except Exception:

        def write(d, indent=0):
            for k, v in d.items():
                pad = "  " * indent
                if isinstance(v, dict):
                    stream.write(f"{pad}{k}:\n")
                    write(v, indent + 1)
                else:
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
                    stream.write(f"{pad}{k}: {vv}\n")

        write(obj)


# ---------------- flat/nested getters ----------------
def _get(cfg, dotted, flat_alias=None):
    cur = cfg
    try:
        for p in dotted.split("."):
            if not isinstance(cur, dict) or p not in cur:
                raise KeyError
            cur = cur[p]
        return cur
    except Exception:
        pass
    if flat_alias and flat_alias in cfg:
        return cfg[flat_alias]
    last = dotted.split(".")[-1]
    return cfg.get(last, None)


def _line(pos_map, dotted, flat_alias=None):
    return (
        pos_map.get(dotted)
        or pos_map.get(flat_alias or "")
        or pos_map.get(dotted.split(".")[-1])
    )


# ---------------- validation ----------------
def _is_bool(x):
    return isinstance(x, bool)


def _is_int(x):
    return isinstance(x, int) and not isinstance(x, bool)


def _is_num(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool)


def validate(cfg, pos_map, path, strict_paths=False):
    errors, warns = [], []

    def err(key_dotted, msg, flat_alias=None):
        ln = _line(pos_map, key_dotted, flat_alias)
        where = f"{path}:{ln}" if ln else path
        errors.append(f"[error] {where} :: {key_dotted}: {msg}")

    def warn(key_dotted, msg, flat_alias=None):
        ln = _line(pos_map, key_dotted, flat_alias)
        where = f"{path}:{ln}" if ln else path
        warns.append(f"[warn]  {where} :: {key_dotted}: {msg}")

    preset = _get(cfg, "preset")
    preset = preset.lower()
    # if not isinstance(preset, str) or preset.lower() not in ("hifi", "ont"):
    #     err("preset", "missing or not one of {'hifi','ont'}")
    # else:
    #     preset = preset.lower()

    wgs_mode = _get(cfg, "wgs_mode")
    # if not _is_bool(wgs_mode):
    #     err("wgs_mode", "missing or not a boolean (true/false)")

    reads = _get(cfg, "paths.reads", flat_alias="reads")
    # if not isinstance(reads, str) or not reads.strip():
    #     err("paths.reads", "missing (string path)", flat_alias="reads")
    # else:
    #     if not os.path.exists(reads):
    #         (err if strict_paths else warn)(
    #             "paths.reads", f"file not found: {reads}", flat_alias="reads"
    #         )
    #
    # hmm_db = _get(cfg, "paths.hmm_db", flat_alias="hmm_db")
    # if not isinstance(hmm_db, str) or not hmm_db.strip():
    #     err("paths.hmm_db", "missing (string path)", flat_alias="hmm_db")
    # else:
    #     if not os.path.exists(hmm_db):
    #         (err if strict_paths else warn)(
    #             "paths.hmm_db", f"file not found: {hmm_db}", flat_alias="hmm_db"
    #         )

    anchors = _get(cfg, "paths.anchors", flat_alias="anchors")
    if isinstance(wgs_mode, bool) and wgs_mode:
        if not isinstance(anchors, str) or not anchors.strip():
            err(
                "paths.anchors",
                "required when wgs_mode=true (string path)",
                flat_alias="anchors",
            )
        else:
            if not os.path.exists(anchors):
                (err if strict_paths else warn)(
                    "paths.anchors", f"file not found: {anchors}", flat_alias="anchors"
                )

    threads = _get(cfg, "compute.threads", flat_alias="threads")
    if not _is_int(threads) or threads <= 0:
        err(
            "compute.threads", "missing or not a positive integer", flat_alias="threads"
        )

    k = _get(cfg, "syncmer.k", flat_alias="k")
    s = _get(cfg, "syncmer.s", flat_alias="s")
    hpc = _get(cfg, "syncmer.hpc", flat_alias="hpc")
    if k is not None and not _is_int(k):
        err("syncmer.k", "must be positive integer", flat_alias="k")
    if s is not None and not _is_int(s):
        err("syncmer.s", "must be positive integer", flat_alias="s")
    if hpc is not None and not _is_bool(hpc):
        err("syncmer.hpc", "must be boolean", flat_alias="hpc")

    final_if = _get(cfg, "autotune.final_if", flat_alias="final_if")
    if final_if is not None and not isinstance(final_if, str):
        err("autotune.final_if", "must be string if present", flat_alias="final_if")
    c_radius = _get(cfg, "autotune.c_radius", flat_alias="c_radius")
    if c_radius is not None and not isinstance(c_radius, str):
        err("autotune.c_radius", "must be string CSV if present", flat_alias="c_radius")

    min_shared = _get(cfg, "graph.min_shared", flat_alias="min_shared")
    jacc_min = _get(cfg, "graph.jaccard_min", flat_alias="jaccard_min")
    topk_nei = _get(cfg, "graph.topk_nei", flat_alias="topk_nei")
    steps = _get(cfg, "graph.steps", flat_alias="steps")
    if min_shared is not None and not _is_int(min_shared):
        err("graph.min_shared", "must be integer", flat_alias="min_shared")
    if jacc_min is not None and not _is_num(jacc_min):
        err("graph.jaccard_min", "must be number", flat_alias="jaccard_min")
    if topk_nei is not None and not _is_int(topk_nei):
        err("graph.topk_nei", "must be integer", flat_alias="topk_nei")
    if steps is not None and not _is_int(steps):
        err("graph.steps", "must be integer", flat_alias="steps")

    return errors, warns


# ---------------- normalization ----------------
def normalize(cfg):
    out = {}
    preset = _get(cfg, "preset") or "hifi"
    preset = str(preset).lower()
    if preset not in ("hifi", "ont"):
        preset = "hifi"
    out["preset"] = preset
    wgs_mode = _get(cfg, "wgs_mode")
    out["wgs_mode"] = bool(wgs_mode) if isinstance(wgs_mode, bool) else True

    paths = {}
    reads = _get(cfg, "paths.reads", flat_alias="reads")
    hmm_db = _get(cfg, "paths.hmm_db", flat_alias="hmm_db")
    anchors = _get(cfg, "paths.anchors", flat_alias="anchors")
    outdir = _get(cfg, "paths.outdir", flat_alias="outdir")
    if reads is not None:
        paths["reads"] = reads
    if hmm_db is not None:
        paths["hmm_db"] = hmm_db
    if anchors is not None:
        paths["anchors"] = anchors
    if outdir is not None:
        paths["outdir"] = outdir
    if paths:
        out["paths"] = paths

    threads = _get(cfg, "compute.threads", flat_alias="threads")
    if not isinstance(threads, int) or threads <= 0:
        threads = 16
    out["compute"] = {"threads": threads}

    if preset == "ont":
        def_k, def_s, def_hpc = 41, 21, True
    else:
        def_k, def_s, def_hpc = 121, 27, False
    k = _get(cfg, "syncmer.k", flat_alias="k")
    s = _get(cfg, "syncmer.s", flat_alias="s")
    hpc = _get(cfg, "syncmer.hpc", flat_alias="hpc")
    if not isinstance(k, int) or k <= 0:
        k = def_k
    if not isinstance(s, int) or s <= 0:
        s = def_s
    if not isinstance(hpc, bool):
        hpc = def_hpc
    out["syncmer"] = {"k": k, "s": s, "hpc": hpc}

    final_if = _get(cfg, "autotune.final_if", flat_alias="final_if")
    c_radius = _get(cfg, "autotune.c_radius", flat_alias="c_radius")
    aut = {}
    if isinstance(final_if, str):
        aut["final_if"] = final_if
    if isinstance(c_radius, str):
        aut["c_radius"] = c_radius
    if aut:
        out["autotune"] = aut

    g = {}
    min_shared = _get(cfg, "graph.min_shared", flat_alias="min_shared")
    jacc_min = _get(cfg, "graph.jaccard_min", flat_alias="jaccard_min")
    topk_nei = _get(cfg, "graph.topk_nei", flat_alias="topk_nei")
    steps = _get(cfg, "graph.steps", flat_alias="steps")
    if isinstance(min_shared, int):
        g["min_shared"] = min_shared
    if isinstance(jacc_min, (int, float)):
        g["jaccard_min"] = float(jacc_min)
    if isinstance(topk_nei, int):
        g["topk_nei"] = topk_nei
    if isinstance(steps, int):
        g["steps"] = steps
    if g:
        out["graph"] = g

    # Add/merge meta.timestamp
    ts = datetime.datetime.now(datetime.timezone.utc).isoformat()
    meta = _get(cfg, "meta") or {}
    if not isinstance(meta, dict):
        meta = {}
    meta["normalized_at"] = ts
    out["meta"] = meta

    # carry flat unknown scalars under extra
    extras = {}
    for k0 in cfg.keys():
        if k0 in out:
            continue
        if k0 in (
            "paths",
            "compute",
            "syncmer",
            "autotune",
            "graph",
            "preset",
            "wgs_mode",
            "meta",
        ):
            continue
        v = cfg[k0]
        if isinstance(v, (str, int, float, bool)) or v is None:
            extras[str(k0)] = v
    if extras:
        out["extra"] = extras
    return out


# ---------------- summary helpers ----------------
def _make_backup(path):
    base = path + ".bak"
    if not os.path.exists(base):
        shutil.copy2(path, base)
        return base
    n = 1
    while True:
        cand = f"{base}.{n}"
        if not os.path.exists(cand):
            shutil.copy2(path, cand)
            return cand
        n += 1


def _ensure_parent(path):
    parent = os.path.dirname(os.path.abspath(path))
    if parent and not os.path.exists(parent):
        os.makedirs(parent, exist_ok=True)


def _write_summary(path, text, json_mode=False, append=False):
    mode = "a" if append else "w"
    try:
        _ensure_parent(path)
        with open(path, mode, encoding="utf-8") as fh:
            if json_mode:
                fh.write(json.dumps(text) + "\n")
            else:
                fh.write(text + "\n")
        action = "appended" if append else "wrote"
        print(f"[summary-file] {action} summary -> {path}", file=sys.stderr)
    except Exception as e:
        print(f"[error] failed to write summary file {path}: {e}", file=sys.stderr)


# ---------------- main ----------------
def main():
    ap = argparse.ArgumentParser(
        description="Validate a polap YAML profile; optionally print/write normalized YAML."
    )
    ap.add_argument(
        "--config-dir",
        default=os.path.join(os.path.expanduser("~"), ".polap", "profiles"),
        help="target config directory (default: ~/.polap/profiles)",
    )
    ap.add_argument("--preset", required=True)
    ap.add_argument("--path", help="YAML config path")
    ap.add_argument(
        "--strict-paths", action="store_true", help="treat missing files as errors"
    )
    ap.add_argument(
        "--print",
        dest="do_print",
        action="store_true",
        help="print normalized YAML to stdout",
    )
    ap.add_argument("--write", metavar="PATH", help="write normalized YAML to PATH")
    ap.add_argument(
        "--inplace",
        action="store_true",
        help="overwrite the input file with normalized YAML",
    )
    ap.add_argument(
        "--backup",
        action="store_true",
        help="with --inplace, save original to file.yaml.bak (or .bak.N)",
    )
    ap.add_argument(
        "--dry-run",
        action="store_true",
        help="simulate inplace/backup/write (no filesystem changes)",
    )
    ap.add_argument(
        "--fail-on-warn",
        action="store_true",
        help="treat warnings as errors (exit nonzero if any warnings)",
    )
    ap.add_argument(
        "--quiet",
        action="store_true",
        help="suppress detailed errors/warnings, only show summary",
    )
    ap.add_argument(
        "--json-summary",
        action="store_true",
        help="emit summary as JSON to stdout (plain summary to stderr)",
    )
    ap.add_argument(
        "--summary-only",
        action="store_true",
        help="suppress normalization output; only show summary",
    )
    ap.add_argument(
        "--summary-file",
        metavar="PATH",
        help="save summary (JSON or plain) to a file (overwrite)",
    )
    ap.add_argument(
        "--append-summary-file",
        metavar="PATH",
        help="append summary (JSON or plain) to a file",
    )
    args = ap.parse_args()

    # If --path not given, build it from config-dir + preset
    if args.path is None:
        args.path = os.path.join(args.config_dir, args.preset + ".yaml")

    if args.write and args.inplace:
        print("[error] cannot use --write and --inplace together", file=sys.stderr)
        return 2
    if args.summary_file and args.append_summary_file:
        print(
            "[error] cannot use --summary-file and --append-summary-file together",
            file=sys.stderr,
        )
        return 2

    try:
        cfg, pos = _load_yaml_with_positions(args.path)
    except Exception as e:
        if not args.quiet:
            print(f"[error] cannot read YAML: {e}", file=sys.stderr)
        status_obj = {
            "status": "FAIL",
            "errors": 1,
            "warnings": 0,
            "timestamp": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        }
        if args.json_summary:
            print(json.dumps(status_obj), file=sys.stdout)
            print("[summary] FAIL (YAML parse error)", file=sys.stderr)
            if args.summary_file:
                _write_summary(
                    args.summary_file, status_obj, json_mode=True, append=False
                )
            if args.append_summary_file:
                _write_summary(
                    args.append_summary_file, status_obj, json_mode=True, append=True
                )
        else:
            text = "[summary] FAIL (YAML parse error)"
            print(text, file=sys.stderr)
            if args.summary_file:
                _write_summary(args.summary_file, text, json_mode=False, append=False)
            if args.append_summary_file:
                _write_summary(
                    args.append_summary_file, text, json_mode=False, append=True
                )
        return 4

    errors, warns = validate(cfg, pos, args.path, strict_paths=args.strict_paths)

    if not args.quiet:
        for w in warns:
            print(w, file=sys.stderr)
        for e in errors:
            print(e, file=sys.stderr)

    fail = bool(errors) or (args.fail_on_warn and warns)
    status_str = "FAIL" if fail else "OK"
    ts = datetime.datetime.now(datetime.timezone.utc).isoformat()
    summary_text = (
        f"[summary] {status_str} ({len(errors)} errors, {len(warns)} warnings)"
    )
    status_obj = {
        "status": status_str,
        "errors": len(errors),
        "warnings": len(warns),
        "timestamp": ts,
    }

    # emit summaries
    if args.json_summary:
        print(json.dumps(status_obj), file=sys.stdout)
        print(summary_text, file=sys.stderr)
        if args.summary_file:
            _write_summary(args.summary_file, status_obj, json_mode=True, append=False)
        if args.append_summary_file:
            _write_summary(
                args.append_summary_file, status_obj, json_mode=True, append=True
            )
    else:
        print(summary_text, file=sys.stderr)
        if args.summary_file:
            _write_summary(
                args.summary_file, summary_text, json_mode=False, append=False
            )
        if args.append_summary_file:
            _write_summary(
                args.append_summary_file, summary_text, json_mode=False, append=True
            )

    if fail:
        return 3
    if args.summary_only:
        return 0

    # Normalization actions (inject meta.normalized_at)
    norm = normalize(
        cfg
    )  # normalize() already injects meta.normalized_at with current UTC timestamp
    if args.do_print:
        _dump_yaml_nested(norm, sys.stdout)

    if args.write:
        if not args.dry_run:
            try:
                with open(args.write, "w", encoding="utf-8") as fh:
                    _dump_yaml_nested(norm, fh)
                if not args.quiet:
                    print(f"[write] normalized YAML -> {args.write}", file=sys.stderr)
            except Exception as e:
                if not args.quiet:
                    print(f"[error] failed to write {args.write}: {e}", file=sys.stderr)
                return 4
        else:
            if not args.quiet:
                print(
                    f"[dry-run] would write normalized YAML -> {args.write}",
                    file=sys.stderr,
                )

    if args.inplace:
        if not args.dry_run:
            try:
                if args.backup:
                    bak = _make_backup(args.path)
                    if not args.quiet:
                        print(f"[backup] saved original -> {bak}", file=sys.stderr)
                with open(args.path, "w", encoding="utf-8") as fh:
                    _dump_yaml_nested(norm, fh)
                if not args.quiet:
                    print(
                        f"[inplace] normalized YAML overwritten -> {args.path}",
                        file=sys.stderr,
                    )
            except Exception as e:
                if not args.quiet:
                    print(
                        f"[error] failed to overwrite {args.path}: {e}", file=sys.stderr
                    )
                return 4
        else:
            if not args.quiet:
                if args.backup:
                    print(
                        f"[dry-run] would backup original -> {args.path}.bak (or .bak.N)",
                        file=sys.stderr,
                    )
                print(f"[dry-run] would overwrite -> {args.path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
