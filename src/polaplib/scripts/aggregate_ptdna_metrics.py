#!/usr/bin/env python3
# Version: v0.1.0
# Purpose: Merge recruited/ds FASTQ stats + (optional) mosdepth window stats → single CSV row
# Usage  : aggregate_ptdna_metrics.py recruited.stats.tsv ds.stats.tsv pre.regions.bed.gz post.regions.bed.gz out.csv [key=val ...]
import sys, gzip, statistics as st, csv, os


def read_stats(tsv_path):
    hdr = None
    vals = {}
    with open(tsv_path, "rt") as f:
        for i, line in enumerate(f):
            parts = line.rstrip("\n").split("\t")
            if i == 0:
                hdr = parts
                continue
            for k, v in zip(hdr, parts):
                vals[k] = v
            break
    return vals


def read_mosdepth_regions(bed_gz):
    if bed_gz == "-" or not os.path.exists(bed_gz):
        return None
    depths = []
    with gzip.open(bed_gz, "rt") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            depths.append(float(parts[3]))
    if not depths:
        return None
    depths.sort()

    def q(p):
        idx = int(round((p / 100.0) * (len(depths) - 1)))
        return depths[idx]

    return {
        "min": depths[0],
        "median": st.median(depths),
        "p05": q(5),
        "p25": q(25),
        "p75": q(75),
        "p95": q(95),
        "max": depths[-1],
    }


if len(sys.argv) < 6:
    sys.exit(
        "Usage: aggregate_ptdna_metrics.py rec.tsv ds.tsv pre.bed.gz post.bed.gz out.csv [key=val ...]"
    )

rec_tsv, ds_tsv, pre_bed, post_bed, out_csv = sys.argv[1:6]
kvs = dict(kv.split("=", 1) for kv in sys.argv[6:]) if len(sys.argv) > 6 else {}

rec = read_stats(rec_tsv)
ds = read_stats(ds_tsv)
pre = read_mosdepth_regions(pre_bed)
post = read_mosdepth_regions(post_bed)

row = {
    "outdir": kvs.get("outdir", ""),
    "cov": kvs.get("cov", ""),
    "gsize": kvs.get("gsize", ""),
    "method": kvs.get("method", ""),
    "minlen": kvs.get("minlen", ""),
    "window": kvs.get("window", ""),
    "recruited_reads": rec.get("reads", "0"),
    "recruited_bases": rec.get("bases", "0"),
    "recruited_mean": rec.get("mean", "0"),
    "recruited_median": rec.get("median", "0"),
    "recruited_N50": rec.get("N50", "0"),
    "recruited_p95": rec.get("p95", "0"),
    "recruited_max": rec.get("max", "0"),
    "ds_reads": ds.get("reads", "0"),
    "ds_bases": ds.get("bases", "0"),
    "ds_mean": ds.get("mean", "0"),
    "ds_median": ds.get("median", "0"),
    "ds_N50": ds.get("N50", "0"),
    "ds_p95": ds.get("p95", "0"),
    "ds_max": ds.get("max", "0"),
}
if pre:
    for k, v in pre.items():
        row[f"pre_{k}"] = f"{v:.3f}"
if post:
    for k, v in post.items():
        row[f"post_{k}"] = f"{v:.3f}"

write_header = not os.path.exists(out_csv) or os.stat(out_csv).st_size == 0
with open(out_csv, "a", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(row.keys()))
    if write_header:
        w.writeheader()
    w.writerow(row)
