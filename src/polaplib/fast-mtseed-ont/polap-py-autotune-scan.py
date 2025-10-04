#!/usr/bin/env python3
# polap-py-autotune-scan.py  v0.3.0
#
# Autotune Step 2a (scan) for overlapness with *decoupled* mapping/reducing:
#   1) run a fast scan with --secondary=no (no -N) → scan.paf.gz
#   2) infer reducer gates (min_olen/min_ident/w_floor) from scan PAF
#   3) run OLPY on scan.paf.gz → scan.overlapness.tsv
#   4) accept if non-empty & not too sparse; suggest scan_keep_frac/topk_per_q
#
# Emits on success (stdout):
#   scan_k=...
#   scan_w=...
#   scan_secondary=no
#   scan_mask=...
#   scan_minocc=...
#   scan_min_olen=...
#   scan_min_ident=...
#   scan_w_floor=...
#   scan_keep_frac=...
#   topk_per_q=...
#   scan_overlap_tsv=<.../scan.overlapness.tsv>
# On failure: prints NOTHING and exits 1 (caller should fallback to full ava-ont).
#
# Usage:
#   polap-py-autotune-scan.py \
#     --reads R1.fq.gz --threads 32 \
#     --olpy /path/polap-py-overlapness-from-paf.py \
#     --minimap2 minimap2 --outdir out/03-allvsall \
#     [--target-frac 0.15]
import sys, os, argparse, subprocess, statistics as st, gzip, math


def run(cmd, shell=False):
    return subprocess.run(cmd, shell=shell, check=True)


def total_reads(reads):
    try:
        out = subprocess.check_output(["seqkit", "stats", "-T", reads], text=True)
        return int(out.strip().splitlines()[1].split("\t")[3])
    except Exception:
        return 0


def n50(reads):
    try:
        out = subprocess.check_output(["seqkit", "stats", "-T", reads], text=True)
        return int(float(out.strip().splitlines()[1].split("\t")[-1]))
    except Exception:
        return 0


def parse_overlapness(path):
    n = 0
    sdeg = 0
    degs = []
    try:
        with open(path) as f:
            for ln in f:
                if not ln.strip() or ln[0] == "#":
                    continue
                p = ln.rstrip("\n").split("\t")
                if len(p) < 3:
                    continue
                d = int(float(p[1]))
                n += 1
                sdeg += d
                degs.append(d)
    except FileNotFoundError:
        return 0, 0.0, []
    epr = (sdeg / (2 * n)) if n > 0 else 0.0
    return n, epr, degs


def sample_metrics_from_paf(paf_gz, max_lines=2_000_000, step=1):
    """Return lists of alen, identity, weight from a gz PAF (first 12 cols) sampling every 'step' line."""
    alens = []
    idens = []
    wgts = []
    opener = gzip.open if paf_gz.endswith(".gz") else open
    with opener(paf_gz, "rt") as f:
        for i, ln in enumerate(f, 1):
            if i > max_lines:
                break
            if step > 1 and (i % step) != 0:
                continue
            if not ln.strip() or ln[0] == "#":
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 12:
                continue
            try:
                qlen = float(p[1])
                tlen = float(p[6])
                alen = float(p[10])
                nmatch = float(p[9])
                if alen <= 0 or qlen <= 0 or tlen <= 0:
                    continue
                ident = max(0.0, min(1.0, nmatch / alen))
                frac = max(0.0, min(1.0, alen / min(qlen, tlen)))
                w = ident * frac
                alens.append(alen)
                idens.append(ident)
                wgts.append(w)
            except Exception:
                continue
    return alens, idens, wgts


def quantile(xs, p):
    if not xs:
        return None
    xs = sorted(xs)
    k = max(0, min(len(xs) - 1, int(round(p * (len(xs) - 1)))))
    return xs[k]


# ---------------- CLI ----------------
ap = argparse.ArgumentParser(
    description="Autotune scan stage for overlapness (decoupled map/reduce)"
)
ap.add_argument("--reads", required=True, help="PT-depleted reads (R1) FASTQ(.gz)")
ap.add_argument("--threads", type=int, default=16)
ap.add_argument(
    "--olpy", required=True, help="path to polap-py-overlapness-from-paf.py"
)
ap.add_argument("--minimap2", default="minimap2")
ap.add_argument("--outdir", required=True)
ap.add_argument(
    "--target-frac",
    type=float,
    default=0.15,
    help="minimum shortlist fraction (default 0.15)",
)
ap.add_argument(
    "--frac-min", type=float, default=0.15, help="lower clamp (default 0.15)"
)
ap.add_argument(
    "--frac-max", type=float, default=0.60, help="upper clamp (default 0.60)"
)
ap.add_argument(
    "--max-paf-lines", type=int, default=2_000_000, help="safety sampling for metrics"
)
args = ap.parse_args()

# -------- tiers for scan seeding (no -N, --secondary=no) --------
tiers = [
    {"k": 21, "w": 9, "mask": 0.70, "occ": 20},
    {"k": 19, "w": 7, "mask": 0.60, "occ": 10},
    {"k": 17, "w": 5, "mask": 0.55, "occ": 8},
    {"k": 15, "w": 5, "mask": 0.50, "occ": 5},
]

N_total = total_reads(args.reads)
N50 = n50(args.reads)
min_olen_floor = 800 if N50 <= 0 else max(800, int(0.25 * N50))

auto_dir = os.path.join(args.outdir, "01-auto")
os.makedirs(auto_dir, exist_ok=True)
scan_paf = os.path.join(auto_dir, "scan.paf.gz")
scan_ovl = os.path.join(auto_dir, "scan.overlapness.tsv")

accepted = None

for t in tiers:
    k = t["k"]
    w = t["w"]
    mask = t["mask"]
    occ = t["occ"]

    # 1) run fast scan → scan.paf.gz  (no -N; --secondary=no)
    cmd_map = (
        f"{args.minimap2} -t {args.threads} -x ava-ont "
        f"-k {k} -w {w} -N 1 --mask-level {mask} --min-occ-floor {occ} "
        f"{args.reads} {args.reads} | gzip -1 > {scan_paf}"
    )
    try:
        run(cmd_map, shell=True)
    except subprocess.CalledProcessError:
        continue

    # 2) metrics from the small scan PAF
    alens, idens, wgts = sample_metrics_from_paf(
        scan_paf, max_lines=args.max_paf_lines, step=1
    )
    if not alens or not idens or not wgts:
        continue

    # 3) suggest reducer gates from scan PAF (robust, like the R helper)
    mo = max(min_olen_floor, int(quantile(alens, 0.50) or min_olen_floor))
    center_ident = quantile(idens, 0.65) or 0.84
    mi = max(0.78, min(0.95, center_ident - 0.05))
    wf = max(0.05, min(0.20, (quantile(wgts, 0.25) or 0.10)))

    # 4) run OLPY on scan PAF with the inferred gates → scan.overlapness.tsv
    cmd_reduce = f"python {args.olpy} {scan_paf} --min_olen {mo} --min_ident {mi:.3f} --w_floor {wf:.3f} > {scan_ovl}"
    try:
        run(cmd_reduce, shell=True)
    except subprocess.CalledProcessError:
        continue

    # 5) acceptance test (non-empty & not too sparse)
    n, epr, degs = parse_overlapness(scan_ovl)
    ok_n = (N_total == 0) or (n >= 0.30 * N_total)
    ok_e = epr >= 3.0 and epr <= 25.0
    if n == 0 or not ok_n:
        # too sparse → try more sensitive seeding
        continue

    accepted = (t, mo, mi, wf, n, epr, degs)
    break

# Failure: emit nothing, exit 1 → caller falls back to full ava-ont
if accepted is None:
    sys.exit(0)

t, mo, mi, wf, n, epr, degs = accepted


# 6) shortlist fraction: cum-wdegree 80% but never below target-frac
def cum80_frac(path):
    vals = []
    with open(path) as f:
        for ln in f:
            if not ln.strip() or ln[0] == "#":
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            try:
                vals.append(float(p[2]))
            except Exception:
                pass
    if not vals:
        return 0.40
    vals.sort(reverse=True)
    tot = sum(vals)
    if tot <= 0:
        return 0.40
    tgt = 0.80 * tot
    acc = 0.0
    for i, v in enumerate(vals, 1):
        acc += v
        if acc >= tgt:
            return i / len(vals)
    return 0.60


frac80 = cum80_frac(scan_ovl)
scan_keep_frac = max(args.frac_min, min(args.frac_max, max(frac80, args.target_frac)))

# 7) Top-K per query cap suggestion from p95(degree)
degs.sort()
p95 = degs[int(round(0.95 * (len(degs) - 1)))] if degs else 0
topk = int(round(0.75 * p95))
if topk < 15:
    topk = 0

# 8) Emit final parameters (shell-friendly)
print(f"scan_k={t['k']}")
print(f"scan_w={t['w']}")
print(f"scan_secondary=no")
print(f"scan_mask={t['mask']}")
print(f"scan_minocc={t['occ']}")
print(f"scan_min_olen={mo}")
print(f"scan_min_ident={mi:.3f}")
print(f"scan_w_floor={wf:.3f}")
print(f"scan_keep_frac={scan_keep_frac:.4f}")
print(f"topk_per_q={topk}")
print(f"scan_overlap_tsv={scan_ovl}")
