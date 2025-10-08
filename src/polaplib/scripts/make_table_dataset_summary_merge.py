#!/usr/bin/env python3
###############################################################################
# scripts/make_table_dataset_summary_merge.py
#
# Version : v0.2.0
# Author  : Sang Chul Choi (POLAP)
# Date    : 2025-10-05
# License : GPL-3.0+
#
# Purpose :
#   Merge multiple POLAP manifest JSON files into one Table S1 (dataset summary),
#   with optional per-platform / per-tier / platform×tier / overall summary rows.
#
# Output columns:
#   species  platform  tier  inum  total_bases  read_count  mean_length  N50
#   avg_qual  gc_content  source_manifest  stat  N
#
# Usage :
#   python scripts/make_table_dataset_summary_merge.py \
#     --manifest-list lists/manifests.txt \
#     --out md/tableS1-dataset-summary.tsv \
#     --platform-map md/platform.tsv \
#     --dedup species --select max-bases \
#     --summary-scopes platform,tier,overall \
#     --summary-stat both \
#     --markdown
###############################################################################
import os, sys, json, csv, argparse, datetime, statistics

NUM_COLS = ["total_bases", "read_count", "mean_length", "N50", "avg_qual", "gc_content"]
BASE_HEADER = (
    ["species", "platform", "tier", "inum"]
    + NUM_COLS
    + ["source_manifest", "stat", "N"]
)


def parse_iso8601(s):
    try:
        return datetime.datetime.strptime(s, "%Y-%m-%dT%H:%M:%SZ")
    except Exception:
        return datetime.datetime.min


def load_platform_map(path):
    # TSV: species<TAB>platform
    mp = {}
    if not path:
        return mp
    with open(path) as f:
        header = f.readline().rstrip("\n")
        # allow either with/without header; if not matching, assume first line is data
        if header.lower().split("\t")[:2] == ["species", "platform"]:
            it = f
        else:
            # treat first line as data
            it = [header] + list(f)
        for line in it:
            if not line.strip():
                continue
            sp, plat = line.rstrip("\n").split("\t")[:2]
            mp[sp] = plat
    return mp


def load_manifest(path):
    with open(path) as f:
        j = json.load(f)
    gen_at = parse_iso8601(j.get("generated_at", ""))
    tier = j.get("tier", "NA")
    items = j.get("items", [])
    rows = []
    for it in items:
        sp = it.get("species", "NA")
        inum = it.get("inum", "NA")
        d = it.get("data", {})
        if not d:  # skip if dataset block missing
            continue
        rows.append(
            {
                "species": sp,
                "tier": tier,
                "inum": inum,
                "total_bases": d.get("total_bases", "NA"),
                "read_count": d.get("read_count", "NA"),
                "mean_length": d.get("mean_length", "NA"),
                "N50": d.get("N50", "NA"),
                "avg_qual": d.get("avg_qual", "NA"),
                "gc_content": d.get("gc_content", "NA"),
                "source_manifest": os.path.basename(path),
                "_generated_at": gen_at,
                "platform": "NA",  # filled later via --platform-map (if provided)
                "stat": "",
                "N": "",
            }
        )
    return rows


def coerce_num(x):
    try:
        # allow int/float in string form
        if isinstance(x, (int, float)):
            return x
        xs = str(x).strip()
        if xs == "" or xs.upper() == "NA":
            return None
        # Some seqkit numbers may be float (e.g., mean_length, avg_qual, gc)
        if "." in xs:
            return float(xs)
        return int(xs)
    except Exception:
        return None


def choose_row(rows, select_mode):
    if select_mode == "max-bases":

        def tb(r):
            v = coerce_num(r.get("total_bases"))
            return -1 if v is None else v

        return max(rows, key=tb)
    elif select_mode == "latest":
        return max(rows, key=lambda r: r.get("_generated_at"))
    elif select_mode == "first":
        return rows[0]
    return rows[0]


def summarize_group(rows, label_platform, label_tier, which_stat):
    # which_stat ∈ {"mean","median"}
    # produce one summary row; use statistic across numeric columns separately
    out = {
        "species": "__SUMMARY__",
        "platform": label_platform,
        "tier": label_tier,
        "inum": "NA",
        "source_manifest": "SUMMARY",
        "stat": which_stat,
        "N": str(len(rows)),
    }
    for col in NUM_COLS:
        vals = [coerce_num(r.get(col)) for r in rows]
        vals = [v for v in vals if v is not None]
        if not vals:
            out[col] = "NA"
            continue
        try:
            if which_stat == "mean":
                out[col] = f"{statistics.fmean(vals):.3f}"
            else:  # median
                out[col] = f"{statistics.median(vals):.3f}"
        except Exception:
            out[col] = "NA"
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--manifest", action="append", default=[], help="Add one manifest (repeatable)"
    )
    ap.add_argument("--manifest-list", help="File with one manifest path per line")
    ap.add_argument("--platform-map", help="TSV: species<TAB>platform")
    ap.add_argument("--out", required=True, help="Output TSV path")
    ap.add_argument("--markdown", action="store_true", default=False)
    ap.add_argument(
        "--dedup", default="species", choices=["species", "species-tier", "none"]
    )
    ap.add_argument(
        "--select", default="max-bases", choices=["max-bases", "latest", "first"]
    )
    ap.add_argument(
        "--summary-scopes",
        default="",
        help="Comma list of scopes: platform,tier,platform-tier,overall",
    )
    ap.add_argument(
        "--summary-stat", default="median", choices=["median", "mean", "both"]
    )
    a = ap.parse_args()

    # Build manifest list
    manifests = list(a.manifest)
    if a.manifest_list:
        if not os.path.exists(a.manifest_list):
            sys.exit(f"[ERR] manifest-list not found: {a.manifest_list}")
        with open(a.manifest_list) as f:
            for line in f:
                p = line.strip()
                if p:
                    manifests.append(p)
    if not manifests:
        sys.exit("[ERR] no manifests provided")

    # Load platform map (optional)
    platform_map = load_platform_map(a.platform_map) if a.platform_map else {}

    # Load rows from manifests
    all_rows = []
    for m in manifests:
        if not os.path.exists(m):
            sys.exit(f"[ERR] manifest not found: {m}")
        rows = load_manifest(m)
        # fill platform if mapped
        for r in rows:
            r["platform"] = platform_map.get(r["species"], "NA")
        all_rows.extend(rows)

    if not all_rows:
        sys.exit("[ERR] no dataset entries found in given manifests")

    # De-dup
    if a.dedup == "none":
        kept = all_rows
    else:
        groups = {}
        for r in all_rows:
            key = (r["species"],) if a.dedup == "species" else (r["species"], r["tier"])
            groups.setdefault(key, []).append(r)
        kept = []
        for key, rows in groups.items():
            kept.append(choose_row(rows, a.select))

    # Sort per-species section
    kept.sort(key=lambda r: (r["species"], r["tier"], r.get("platform", "NA")))

    # Prepare summaries
    summary_rows = []
    scopes = [x.strip() for x in a.summary_scopes.split(",") if x.strip()]
    stats_to_make = (
        ["median"]
        if a.summary_stat == "median"
        else (["mean"] if a.summary_stat == "mean" else ["median", "mean"])
    )

    # Helper to append summary rows for a slice
    def add_summaries(rows, platform_label, tier_label):
        for st in stats_to_make:
            summary_rows.append(summarize_group(rows, platform_label, tier_label, st))

    if scopes:
        # overall
        if "overall" in scopes:
            add_summaries(kept, "ALL", "ALL")

        # by platform
        if "platform" in scopes:
            byp = {}
            for r in kept:
                byp.setdefault(r.get("platform", "NA"), []).append(r)
            for plat, rows in sorted(byp.items()):
                add_summaries(rows, plat, "ALL")

        # by tier
        if "tier" in scopes:
            byt = {}
            for r in kept:
                byt.setdefault(r.get("tier", "NA"), []).append(r)
            for ti, rows in sorted(byt.items()):
                add_summaries(rows, "ALL", ti)

        # platform × tier
        if "platform-tier" in scopes:
            bypt = {}
            for r in kept:
                key = (r.get("platform", "NA"), r.get("tier", "NA"))
                bypt.setdefault(key, []).append(r)
            for (plat, ti), rows in sorted(bypt.items()):
                add_summaries(rows, plat, ti)

    # Write TSV
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(BASE_HEADER)
        for r in kept:  # species rows
            w.writerow([r.get(h, "") for h in BASE_HEADER])
        for r in summary_rows:  # appended summaries
            w.writerow([r.get(h, "") for h in BASE_HEADER])
    print(f"[INFO] Wrote TSV: {a.out}")

    # Markdown (two blocks: species rows then summary rows)
    if a.markdown:
        mdout = os.path.splitext(a.out)[0] + ".md"
        with open(mdout, "w") as f:
            # species block
            f.write("### Table S1 — Dataset summary (per species)\n\n")
            f.write("| " + " | ".join(BASE_HEADER) + " |\n")
            f.write("|" + "|".join(["---"] * len(BASE_HEADER)) + "|\n")
            for r in kept:
                f.write(
                    "| " + " | ".join(str(r.get(h, "")) for h in BASE_HEADER) + " |\n"
                )
            # summary block (only if any)
            if summary_rows:
                f.write("\n\n### Table S1 — Summary rows\n\n")
                f.write("| " + " | ".join(BASE_HEADER) + " |\n")
                f.write("|" + "|".join(["---"] * len(BASE_HEADER)) + "|\n")
                for r in summary_rows:
                    f.write(
                        "| "
                        + " | ".join(str(r.get(h, "")) for h in BASE_HEADER)
                        + " |\n"
                    )
        print(f"[INFO] Wrote Markdown: {mdout}")


if __name__ == "__main__":
    sys.exit(main())
