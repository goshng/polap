#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
polap-py-syncfilter-cut-auto.py
Config + CLI overrides (CLI wins). Strict-band reclassifier for syncfilter quickview TSV.

Inputs
------
Quickview TSV (columns): read_id, n_syncmer, mean, median, q25, q75, class
X = log10(median), ignoring rows with median <= 0.

Addressing
----------
  --config-path FILE.yaml
  --config-dir DIR --preset NAME     (default DIR=~/.polap/profiles)

YAML keys (flat; CLI with '-' overrides)
----------------------------------------
# IO
input: /path/to/quickview.tsv
cuts_out: /path/to/cuts.tsv
tsv_out: /path/to/reclass.tsv
png: /path/to/overlay.png

# banding
x_method: hybrid           # mad|window|std|valley|hybrid
x_tail: 0.05               # right-tail for nuclear cut
tail_mt: 0.30              # right-tail for mt band (falls back to x_tail if empty)
tail_pt: 0.30              # right-tail for pt band (falls back to x_tail if empty)
x_win: 0.20                # window for std around μ_nuc
band_max: 2.0              # define nuclear band X<x
sigma_floor: 0.03
nbins: 320

# anchors
mt_anchors: /path/to/mt.id.txt
pt_anchors: /path/to/pt.id.txt

Strict policy
-------------
- If X within mt band -> mt
- elif X within pt band -> pt
- elif X <= cut_nuc_log -> nuclear
- else -> nuclear

Outputs
-------
- cuts_out TSV with one line for nuclear and (optional) lines for mt, pt.
- tsv_out TSV with 'class' rewritten to nuclear|mt|pt
- <basename>.nuclear.ids / .mt.ids / .pt.ids
- optional PNG overlay (bands + μ)

Logging
-------
- -v/--verbose prints chatty steps to stderr (repeatable).
"""

from __future__ import annotations
import argparse, os, sys, math, json
from typing import Dict, Any, Tuple, List, Optional
from collections import defaultdict

# ---------------- logging ----------------
VERBOSE = 0


def log(level: int, *args):
    if VERBOSE >= level:
        sys.stderr.write("[cut-auto] " + " ".join(str(a) for a in args) + "\n")


# ---------------- YAML loader ----------------
def _load_yaml(path: str) -> Dict[str, Any]:
    if not path or not os.path.exists(path):
        return {}
    try:
        import yaml  # type: ignore

        with open(path, "r", encoding="utf-8") as fh:
            d = yaml.safe_load(fh) or {}
        if not isinstance(d, dict):
            raise ValueError("YAML root is not a mapping")
        # flat mapping only
        return dict(d)
    except ImportError:
        # minimal fallback: flat "k: v" lines
        d: Dict[str, Any] = {}
        with open(path, "r", encoding="utf-8") as fh:
            for ln in fh:
                ln = ln.strip()
                if not ln or ln.startswith("#") or ":" not in ln:
                    continue
                k, v = ln.split(":", 1)
                k = k.strip()
                v = v.strip().strip("'").strip('"')
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


def _resolve_cfg_path(
    config_dir: Optional[str], preset: Optional[str], config_path: Optional[str]
) -> Optional[str]:
    if config_path:
        return os.path.expanduser(config_path)
    if config_dir and preset:
        return os.path.join(os.path.expanduser(config_dir), f"{preset}.yaml")
    return None


# ---------------- core utils ----------------
def read_ids(path: Optional[str]) -> set[str]:
    s: set[str] = set()
    if not path:
        return s
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            tok = line.strip().split()
            if tok:
                s.add(tok[0])
    return s


def load_quickview(tsv: str) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    with open(tsv, "r", encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}
        # expect: read_id, n_syncmer, mean, median, q25, q75, class
        for ln in fh:
            sp = ln.rstrip("\n").split("\t")
            if len(sp) < 6:
                continue
            try:
                rid = sp[idx.get("read_id", 0)].split()[0]  # normalize
                med = float(sp[idx.get("median", 3)])
                if med <= 0:
                    continue
                X = math.log10(med)
                rows.append({"read_id": rid, "median": med, "X": X})
            except Exception:
                pass
    return rows


def sigma_mad(xs: List[float]) -> Optional[float]:
    if not xs:
        return None
    import statistics as st

    med = st.median(xs)
    mad = st.median([abs(v - med) for v in xs])
    if mad == 0:
        try:
            return float(st.pvariance(xs) ** 0.5) if len(xs) > 1 else None
        except Exception:
            return None
    return float(mad / 0.6744897501960817)


# right-tail z
try:
    from statistics import NormalDist
except Exception:
    NormalDist = None


def z_for_tail(t: float) -> float:
    t = float(t)
    if t <= 0.0:  # zero tail => +inf z
        return float("inf")
    if t >= 1.0:
        return 0.0
    if NormalDist is not None:
        return NormalDist().inv_cdf(1.0 - t)
    # fallback (Abramowitz-Stegun 26.2.23)
    import math

    p = 1.0 - t
    a1, a2, a3, a4, a5, a6 = (
        -39.69683028665376,
        220.9460984245205,
        -275.9285104469687,
        138.3577518672690,
        -30.66479806614716,
        2.506628277459239,
    )
    b1, b2, b3, b4 = (
        -54.47609879822406,
        161.5858368580409,
        -155.6989798598866,
        66.80131188771972,
    )
    c1, c2, c3, c4, c5, c6 = (
        -0.007784894002430293,
        -0.3223964580411365,
        -2.400758277161838,
        -2.549732539343734,
        4.374664141464968,
        2.938163982698783,
    )
    d1, d2, d3 = (0.007784695709041462, 0.3224671290700398, 2.445134137142996)
    plow, phigh = 0.02425, 1.0 - 0.02425
    if p < plow:
        q = (-2.0 * math.log(p)) ** 0.5
        x = (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / (
            (((d1 * q + d2) * q + d3) * q + 1.0)
        )
    elif p > phigh:
        q = (-2.0 * math.log(1.0 - p)) ** 0.5
        x = -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / (
            (((d1 * q + d2) * q + d3) * q + 1.0)
        )
    else:
        q = p - 0.5
        r = q * q
        x = (
            (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6)
            * q
            / (((((b1 * r + b2) * r + b3) * r + b4) * r + 1.0))
        )
    return float(x)


def smooth(y: List[float], rounds: int = 2) -> List[float]:
    x = y[:]
    for _ in range(rounds):
        if len(x) < 3:
            break
        x = (
            [x[0]]
            + [
                0.25 * x[i - 1] + 0.5 * x[i] + 0.25 * x[i + 1]
                for i in range(1, len(x) - 1)
            ]
            + [x[-1]]
        )
    return x


def windowed_mode(logx: List[float], nbins: int = 320, win: float = 0.12) -> float:
    import numpy as np

    if not logx:
        return 0.0
    lo, hi = min(logx), max(logx)
    if lo == hi:
        return lo
    hist, edges = np.histogram(logx, bins=nbins, range=(lo, hi))
    ctrs = 0.5 * (edges[:-1] + edges[1:])
    s = smooth(list(hist), 1)
    bw = (hi - lo) / nbins if nbins > 0 else 1e-6
    half = max(1, int((win / 2) / bw))
    mv = []
    for i in range(len(s)):
        L = max(0, i - half)
        R = min(len(s), i + half + 1)
        mv.append(sum(s[L:R]))
    idx = int(max(range(len(mv)), key=lambda i: mv[i]))
    return float(ctrs[idx])


def sigma_window(
    logx: List[float], mu: float, w: float = 0.20, min_n: int = 50
) -> Optional[float]:
    import numpy as np

    sel = [x for x in logx if (mu - w) <= x <= mu]
    if len(sel) < min_n:
        return None
    return float(np.std(sel)) if len(sel) > 1 else None


# plotting
def plot_bands(
    X_all: List[float], bands: Dict[str, Dict[str, Any]], png_path: str
) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return
    import numpy as np

    plt.figure(figsize=(8.5, 4.8))
    plt.hist(np.array(X_all), bins=140, alpha=0.35, label="all")
    colors = {"nuclear": "tab:blue", "mt": "tab:orange", "pt": "tab:green"}
    for t in ("nuclear", "mt", "pt"):
        if t in bands:
            b = bands[t]
            mu = b["mu"]
            rng = b["range"]
            c = colors.get(t, "k")
            if rng[0] is not None:
                plt.axvline(rng[0], linestyle="--", linewidth=1.8, color=c)
            if rng[1] is not None:
                plt.axvline(rng[1], linestyle="--", linewidth=1.8, color=c)
            plt.axvline(
                mu, linestyle=":", linewidth=1.6, color=c, label=f"{t} μ={mu:.3f}"
            )
    plt.xlabel("log10(median syncmer multiplicity / coverage proxy)")
    plt.ylabel("count (reads)")
    plt.legend(loc="upper right", fontsize=9)
    plt.tight_layout()
    try:
        plt.savefig(png_path, dpi=150)
    except Exception:
        pass
    plt.close()


# ---------------- main ----------------
def main() -> int:
    p = argparse.ArgumentParser(
        description="Strict-band cut & reclass from syncfilter quickview TSV."
    )
    # addressing
    p.add_argument("--config-path", help="explicit YAML config path")
    p.add_argument("--config-dir", help="base config dir (default: ~/.polap/profiles)")
    p.add_argument("--preset", help="profile name (with --config-dir)")

    # verbosity
    p.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="increase verbosity (repeatable)",
    )

    # IO
    p.add_argument("-i", "--input", help="syncfilter quickview TSV")
    p.add_argument("--cuts-out", help="write cuts TSV here")
    p.add_argument("--tsv-out", help="write reclass TSV here")
    p.add_argument("--png", help="overlay PNG")

    # banding knobs
    p.add_argument(
        "--method",
        dest="x_method",
        choices=["mad", "window", "std", "valley", "hybrid"],
    )
    p.add_argument("--tail", type=float, dest="x_tail")
    p.add_argument("--tail-mt", type=float)
    p.add_argument("--tail-pt", type=float)
    p.add_argument("--window-width", type=float, dest="x_win")
    p.add_argument("--band-max", type=float)
    p.add_argument("--sigma-floor", type=float)
    p.add_argument("--nbins", type=int)

    # anchors (strict names)
    p.add_argument("--mt-anchors", help="mt anchor ID list file")
    p.add_argument("--pt-anchors", help="pt anchor ID list file")

    args = p.parse_args()
    global VERBOSE
    VERBOSE = int(args.verbose or 0)

    # load YAML
    cfg_path = _resolve_cfg_path(args.config_dir, args.preset, args.config_path)
    cfg = _load_yaml(cfg_path) if cfg_path else {}
    log(1, f"config={cfg_path or '(none)'}")

    # helper: merge with priority CLI > YAML > default
    def get(key: str, default=None, cli_attr: Optional[str] = None):
        cli_key = cli_attr or key
        val_cli = getattr(args, cli_key.replace("-", "_"), None)
        if val_cli is not None:
            return val_cli
        return cfg.get(key, default)

    # IO resolution
    input_tsv = get("input") or cfg.get("quickview_tsv")
    if not input_tsv:
        log(1, "error: need --input or YAML 'input'/'quickview_tsv'")
        return 2
    cuts_out = (
        get("cuts_out")
        or f"{os.path.splitext(os.path.splitext(input_tsv)[0])[0]}.cuts.tsv"
    )
    tsv_out = (
        get("tsv_out")
        or f"{os.path.splitext(os.path.splitext(input_tsv)[0])[0]}.reclass.tsv"
    )
    png_out = get("png")

    # banding params
    method = get("x_method", "hybrid", "x_method")
    tail_nuc = float(get("x_tail", 0.05, "x_tail") or 0.05)
    tail_mt = get("tail_mt", None)
    tail_pt = get("tail_pt", None)
    xwin = float(get("x_win", 0.20, "x_win") or 0.20)
    band_max = float(get("band_max", 2.0, "band_max") or 2.0)
    sigma_floor = float(get("sigma_floor", 0.03, "sigma_floor") or 0.03)
    nbins = int(get("nbins", 320, "nbins") or 320)

    # anchors
    mt_anchors = get("mt_anchors")
    pt_anchors = get("pt_anchors")
    mt_ids = read_ids(mt_anchors)
    pt_ids = read_ids(pt_anchors)

    # load quickview
    rows = load_quickview(input_tsv)
    if not rows:
        log(1, "error: no usable rows in quickview TSV")
        return 3

    X_all = [r["X"] for r in rows]
    X_nuc_band = [x for x in X_all if x < band_max]
    if not X_nuc_band:
        log(1, f"error: no X < {band_max}")
        return 3

    # μ_nuc (robust mode-ish)
    mu_nuc = windowed_mode(X_nuc_band, nbins=nbins, win=xwin)
    # σ_nuc
    if method == "mad":
        sg_nuc = sigma_mad(X_nuc_band) or sigma_floor
    elif method == "std":
        import statistics as st

        sg_nuc = (
            float(st.pvariance(X_nuc_band) ** 0.5)
            if len(X_nuc_band) > 1
            else sigma_floor
        )
    elif method == "window":
        sg_nuc = (
            sigma_window(X_nuc_band, mu_nuc, w=xwin, min_n=50)
            or sigma_mad(X_nuc_band)
            or sigma_floor
        )
    elif method == "valley":
        # valley doesn't define σ; approximate via window
        sg_nuc = (
            sigma_window(X_nuc_band, mu_nuc, w=xwin, min_n=50)
            or sigma_mad(X_nuc_band)
            or sigma_floor
        )
    else:  # hybrid: use window/MAD as fallback
        sg_nuc = (
            sigma_window(X_nuc_band, mu_nuc, w=xwin, min_n=50)
            or sigma_mad(X_nuc_band)
            or sigma_floor
        )
    sg_nuc = max(sg_nuc, sigma_floor)
    z_nuc = z_for_tail(tail_nuc)
    cut_nuc_log = mu_nuc + z_nuc * sg_nuc
    cut_nuc_lin = 10**cut_nuc_log

    log(
        1,
        f"μ_nuc={mu_nuc:.4f} σ_nuc={sg_nuc:.4f} tail={tail_nuc} -> cut_nuc≈{cut_nuc_lin:.2f}x (log10={cut_nuc_log:.3f})",
    )

    # bands from anchors (if provided)
    import statistics as st

    def learn_band_from_ids(ids: set[str]) -> Optional[Tuple[float, float]]:
        X = [r["X"] for r in rows if r["read_id"] in ids]
        if len(X) < 10:
            return None
        mu = float(st.median(X))
        try:
            sg = float(st.pvariance(X) ** 0.5) if len(X) > 1 else None
        except Exception:
            sg = None
        if sg is None or sg <= 0:
            sg = sigma_mad(X) or sigma_floor
        sg = max(sg, sigma_floor)
        return (mu, sg)

    band_mt = learn_band_from_ids(mt_ids) if mt_ids else None
    band_pt = learn_band_from_ids(pt_ids) if pt_ids else None

    # tail for bands
    if tail_mt is None:
        tail_mt = tail_nuc
    if tail_pt is None:
        tail_pt = tail_nuc
    z_mt = z_for_tail(float(tail_mt))
    z_pt = z_for_tail(float(tail_pt))

    # write cuts.tsv
    with open(cuts_out, "w", encoding="utf-8") as w:
        w.write("f1\tf2\tthree_bins\tmethod\tmu_log10\tsigma_log10\ttail\ttype\n")
        w.write(
            f"{cut_nuc_lin:.8f}\tNA\t0\t{method}\t{mu_nuc:.5f}\t{sg_nuc:.5f}\t{tail_nuc}\tnuclear\n"
        )
        if band_mt:
            mu_mt, sg_mt = band_mt
            f1_mt = 10 ** (mu_mt - z_mt * sg_mt)
            f2_mt = 10 ** (mu_mt + z_mt * sg_mt)
            w.write(
                f"{f1_mt:.8f}\t{f2_mt:.8f}\t1\tstd\t{mu_mt:.5f}\t{sg_mt:.5f}\t{tail_mt}\tmt\n"
            )
            log(
                1,
                f"mt μ={mu_mt:.4f} σ={sg_mt:.4f} tail={tail_mt} band=[{f1_mt:.2f}x..{f2_mt:.2f}x]",
            )
        if band_pt:
            mu_pt, sg_pt = band_pt
            f1_pt = 10 ** (mu_pt - z_pt * sg_pt)
            f2_pt = 10 ** (mu_pt + z_pt * sg_pt)
            w.write(
                f"{f1_pt:.8f}\t{f2_pt:.8f}\t1\tstd\t{mu_pt:.5f}\t{sg_pt:.5f}\t{tail_pt}\tpt\n"
            )
            log(
                1,
                f"pt μ={mu_pt:.4f} σ={sg_pt:.4f} tail={tail_pt} band=[{f1_pt:.2f}x..{f2_pt:.2f}x]",
            )

    log(1, f"[cuts] {cuts_out}")

    # build bands in log space for reclass
    bands: Dict[str, Dict[str, Any]] = {
        "nuclear": {"mu": mu_nuc, "sigma": sg_nuc, "range": (cut_nuc_log, cut_nuc_log)}
    }
    rng_mt = None
    rng_pt = None
    if band_mt:
        mu_mt, sg_mt = band_mt
        rng_mt = (mu_mt - z_mt * sg_mt, mu_mt + z_mt * sg_mt)
        bands["mt"] = {"mu": mu_mt, "sigma": sg_mt, "range": rng_mt}
    if band_pt:
        mu_pt, sg_pt = band_pt
        rng_pt = (mu_pt - z_pt * sg_pt, mu_pt + z_pt * sg_pt)
        bands["pt"] = {"mu": mu_pt, "sigma": sg_pt, "range": rng_pt}

    def in_rng(x: float, rng: Optional[Tuple[float, float]]) -> bool:
        return (
            (rng is not None)
            and (rng[0] is not None)
            and (rng[1] is not None)
            and (rng[0] <= x <= rng[1])
        )

    # reclassify strictly
    for r in rows:
        x = r["X"]
        if rng_mt and in_rng(x, rng_mt):
            r["class"] = "mt"
        elif rng_pt and in_rng(x, rng_pt):
            r["class"] = "pt"
        elif x <= cut_nuc_log:
            r["class"] = "nuclear"
        else:
            r["class"] = "nuclear"

    # write reclass.tsv + ids
    with open(tsv_out, "w", encoding="utf-8") as w:
        # w.write("read_id\tn_syncmer\tmean\tmedian\tq25\tq75\tclass\n")
        w.write("read_id\tn_syncmer\tmean\tmedian\tq25\tq75\tclass\n")
        for r in rows:
            # we only guarantee 'read_id','median','class' here; other columns left blank
            w.write(f"{r['read_id']}\t\t\t{r['median']:.3f}\t\t\t{r['class']}\n")
    log(1, f"[reclass] {tsv_out}")

    base = os.path.splitext(tsv_out)[0]
    buckets: Dict[str, List[str]] = defaultdict(list)
    for r in rows:
        buckets[r["class"]].append(r["read_id"])
    for name in ("nuclear", "mt", "pt"):
        fn = f"{base}.{name}.ids"
        with open(fn, "w", encoding="utf-8") as w:
            for rid in buckets.get(name, []):
                w.write(rid + "\n")
        log(1, f"[ids] {name} -> {fn} (n={len(buckets.get(name, []))})")

    # optional plot
    if png_out:
        plot_bands(X_all, bands, png_out)
        log(1, f"[plot] {png_out}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
