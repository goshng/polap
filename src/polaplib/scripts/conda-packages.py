#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# conda-packages.py
# Version: v1.1.0 (conda-only, single-process, optional --simple-channel)
#
# Reads a concatenated combined-environment.yml (multiple env blocks with headers),
# extracts every conda dependency under any 'dependencies:' list, collects channels
# (+ any namespaces like 'goshng::pkg'), and prints CSV:
#   package,latest_version,channel,build
#
# Usage:
#   python3 conda-packages.py [--simple-channel] combined-environment.yml > latest.csv

import os
import re
import sys
import json
import shutil
import subprocess
from typing import List, Tuple, Dict, Any

DEFAULT_CHANNELS = ["conda-forge", "bioconda", "defaults"]

# --- Simple YAML-ish scanners (no external deps) ------------------------------
DEPENDENCIES_START = re.compile(r"^[ \t]*dependencies:[ \t]*$")
LIST_ITEM = re.compile(r"^[ \t]*-[ \t]+(.+?)\s*$")
CHANNELS_START = re.compile(r"^[ \t]*channels:[ \t]*$")
CHANNEL_ITEM = re.compile(r"^[ \t]*-[ \t]+([A-Za-z0-9_.+-]+)\s*$")
NS_PREFIX = re.compile(r"^([A-Za-z0-9_.+-]+)::")


def runcmd(cmd: List[str]) -> Tuple[int, str, str]:
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    out, err = p.communicate()
    return p.returncode, out, err


def ensure_conda_in_env():
    """Try to make `conda` usable even if the shell wasn't initialized."""
    if shutil.which("conda"):
        return
    for base in (
        os.environ.get("CONDA_PREFIX", ""),
        os.path.expanduser("~/miniconda3"),
        os.path.expanduser("~/mambaforge"),
        os.path.expanduser("~/anaconda3"),
        "/opt/conda",
    ):
        if not base:
            continue
        sh = os.path.join(base, "etc", "profile.d", "conda.sh")
        if os.path.exists(sh):
            try:
                out = subprocess.check_output(
                    ["bash", "-lc", f'. "{sh}" && conda activate base && env'],
                    text=True,
                )
                for line in out.splitlines():
                    k, _, v = line.partition("=")
                    if k:
                        os.environ[k] = v
                return
            except Exception:
                pass


def read_lines(path: str) -> List[str]:
    return open(path, encoding="utf-8", errors="ignore").read().splitlines(True)


def extract_specs_and_channels(lines: List[str]) -> Tuple[List[str], List[str]]:
    in_deps = False
    in_channels = False
    specs: List[str] = []
    channels: List[str] = []
    for raw in lines:
        line = raw.rstrip("\n")

        if CHANNELS_START.match(line):
            in_channels, in_deps = True, False
            continue
        if in_channels:
            m = CHANNEL_ITEM.match(line)
            if m:
                channels.append(m.group(1))
                continue
            in_channels = False

        if DEPENDENCIES_START.match(line):
            in_deps, in_channels = True, False
            continue
        if in_deps:
            m = LIST_ITEM.match(line)
            if m:
                item = m.group(1).strip()
                if item.lower().startswith("pip:"):
                    continue
                specs.append(item)
                continue
            in_deps = False

    # de-dup preserve order
    def dedup(xs):
        seen = set()
        out = []
        for x in xs:
            if x not in seen:
                seen.add(x)
                out.append(x)
        return out

    return dedup(specs), dedup(channels)


def normalize_pkg(spec: str) -> str:
    """Drop selectors, namespaces, and version pins."""
    s = spec.strip()
    s = s.split("[", 1)[0]
    s = s.split("::", 1)[-1]
    if "==" in s:
        s = s.split("==", 1)[0]
    if "=" in s:
        s = s.split("=", 1)[0]
    return s.strip()


def collect_namespaces(specs: List[str]) -> List[str]:
    ns = []
    seen = set()
    for s in specs:
        m = NS_PREFIX.match(s.strip())
        if m:
            ch = m.group(1)
            if ch not in seen:
                seen.add(ch)
                ns.append(ch)
    return ns


def pick_latest(records: List[Dict[str, Any]]) -> Tuple[str, str, str]:
    """Choose newest primarily by timestamp, then build_number, then version string."""
    if not records:
        return ("NOT_FOUND", "", "")

    def key(rec):
        ts = int(rec.get("timestamp", 0))
        bn = int(rec.get("build_number", 0))
        ver = rec.get("version", "")
        return (ts, bn, ver)

    rec = sorted(records, key=key)[-1]
    return (
        str(rec.get("version", "")),
        str(rec.get("channel", rec.get("subdir", ""))),
        str(rec.get("build", "")),
    )


def simplify_channel(url: str) -> str:
    """Simplify long channel URLs to short forms."""
    if not url:
        return url
    url = url.strip()
    # common patterns
    if url.startswith("https://conda.anaconda.org/"):
        url = url.replace("https://conda.anaconda.org/", "")
    elif url.startswith("http://conda.anaconda.org/"):
        url = url.replace("http://conda.anaconda.org/", "")
    elif url.startswith("https://repo.anaconda.com/pkgs/"):
        url = url.replace("https://repo.anaconda.com/pkgs/", "pkgs/")
    # remove trailing platform dirs like linux-64, noarch, etc.
    url = re.sub(r"/(linux|osx|win|noarch)[-_]?\d*/*$", "", url)
    return url.rstrip("/")


def conda_search_one(
    pkg: str, channels: List[str], simple: bool
) -> Tuple[str, str, str]:
    cmd = ["conda", "search", "--json", "--override-channels", "-q"]
    for ch in channels:
        cmd += ["-c", ch]
    cmd.append(pkg)
    rc, out, err = runcmd(cmd)
    if rc != 0 or not out.strip():
        sys.stderr.write(
            f"# conda QUERY_ERROR for {pkg}: rc={rc} stderr={err.strip()}\n"
        )
        return ("QUERY_ERROR", "", "")
    try:
        j = json.loads(out)
        if "packages" in j:
            recs = list(j["packages"].values())
        elif pkg in j:
            recs = j[pkg]
        else:
            recs = []
        ver, ch, build = pick_latest(recs)
        if simple:
            ch = simplify_channel(ch)
        return (ver, ch, build)
    except Exception as e:
        sys.stderr.write(f"# conda JSON parse error for {pkg}: {e}\n")
        return ("QUERY_ERROR", "", "")


def main():
    if len(sys.argv) < 2:
        print(
            "Usage: conda-packages.py [--simple-channel] combined-environment.yml",
            file=sys.stderr,
        )
        sys.exit(1)

    simple = False
    args = sys.argv[1:]
    if args[0] == "--simple-channel":
        simple = True
        args = args[1:]

    if not args:
        print("ERROR: Missing input YAML path.", file=sys.stderr)
        sys.exit(1)

    ipath = args[0]
    if not os.path.isfile(ipath):
        print(f"ERROR: file not found: {ipath}", file=sys.stderr)
        sys.exit(1)

    ensure_conda_in_env()
    if not shutil.which("conda"):
        print("ERROR: conda not found in PATH.", file=sys.stderr)
        sys.exit(2)

    lines = read_lines(ipath)
    specs, file_channels = extract_specs_and_channels(lines)
    namespaces = collect_namespaces(specs)

    # Build channel list: namespaces first, then file channels (or defaults)
    chs: List[str] = []
    chs.extend(namespaces)
    chs.extend(file_channels or DEFAULT_CHANNELS)
    seen = set()
    uniq: List[str] = []
    for c in chs:
        if c and c not in seen:
            seen.add(c)
            uniq.append(c)
    if not uniq:
        uniq = DEFAULT_CHANNELS[:]

    # Normalize & de-duplicate package names
    pkgs: List[str] = []
    seenp = set()
    for s in specs:
        p = normalize_pkg(s)
        if p and p not in seenp:
            seenp.add(p)
            pkgs.append(p)

    print("package,latest_version,channel,build")
    if not pkgs:
        return

    for pkg in pkgs:
        ver, ch, build = conda_search_one(pkg, uniq, simple)
        print(f"{pkg},{ver},{ch},{build}")


if __name__ == "__main__":
    main()
