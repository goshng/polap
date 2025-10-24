#!/usr/bin/env python3
# scripts/report_readassemble_pt.py
# Version: v0.2.0
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Build a JSON-only KPI bundle for PT-1 (plastid read selection & assembly).
# Inputs are discovered under --base-dir (polap-readassemble/annotate-read-pt).
#
# Output:
#   --out <base-dir>/report/pt1-report.json
#
# Fields include:
#  - pt_reads: count, total_bases (from pt-contig-annotation-depth-table.txt where PT>MT)
#  - pt_input_seqkit: num_seqs,sum_len,avg_len,N50,AvgQual,GC(%) from pt.fq.seqkit.stats.ta.txt
#  - coverage pre/post: mean, median from pt1/tmp/pre.cov.csv / post.cov.csv (if present)
#  - genes: pt_gene_count (pt.gene.count), annotated_pt_length (sum Length where PT>MT)
#  - assembly pt0/pt1: n,total_bp,N50,max_bp (prefer contigs_stats.txt; fallback to GFA S lines)
#  - isoforms: circular_path_count from pt1/ptdna/ptdna/circular_path_count.txt
#  - thresholds: mapq, aa_len from pt_thresh.vars (if present)
#  - qc artifact paths (graph, coverage pngs, scatter pdf)
#
import os, sys, json, argparse


def read_first_int(path):
    try:
        with open(path) as f:
            line = f.readline().strip()
            if not line:
                return None
            token = line.split()[0]
            return int(token)
    except Exception:
        return None


def parse_seqkit_ta(path):
    """Return dict with num_seqs,sum_len,avg_len,N50,AvgQual,GC(%) or blanks."""
    res = {
        "num_seqs": "",
        "sum_len": "",
        "avg_len": "",
        "N50": "",
        "AvgQual": "",
        "GC(%)": "",
    }
    if not path or not os.path.isfile(path):
        return res
    try:
        with open(path, encoding="utf-8") as f:
            lines = [ln.rstrip("\n\r") for ln in f]
        if len(lines) < 2:
            return res
        hdr = lines[0].split("\t")
        row = lines[1].split("\t")
        h = {k: i for i, k in enumerate(hdr)}

        def get(key):
            idx = h.get(key)
            return row[idx] if idx is not None and idx < len(row) else ""

        res["num_seqs"] = get("num_seqs")
        res["sum_len"] = get("sum_len")
        res["avg_len"] = get("avg_len")
        res["N50"] = get("N50")
        res["AvgQual"] = get("AvgQual")
        # GC sometimes appears as "GC(%)"
        res["GC(%)"] = (
            get("GC(%)") if "GC(%)" in h else (get("GC") if "GC" in h else "")
        )
        return res
    except Exception:
        return res


def parse_pt_table(path):
    """
    From pt-contig-annotation-depth-table.txt, compute:
      pt_count = #rows where PT>MT
      pt_bases = sum of Length where PT>MT
    Expect header: Contig Length Depth Copy MT PT Edge
    """
    res = {"pt_count": "", "pt_bases": ""}
    if not path or not os.path.isfile(path):
        return res
    try:
        import csv

        with open(path, encoding="utf-8") as f:
            rd = csv.reader(f, delimiter="\t")
            hdr = next(rd, None)
            if not hdr:
                return res
            pos = {name: i for i, name in enumerate(hdr)}
            i_len = pos.get("Length")
            i_mt = pos.get("MT")
            i_pt = pos.get("PT")
            if i_len is None or i_mt is None or i_pt is None:
                # try case-insensitive
                low = {name.lower(): i for i, name in enumerate(hdr)}
                i_len = low.get("length")
                i_mt = low.get("mt")
                i_pt = low.get("pt")
                if i_len is None or i_mt is None or i_pt is None:
                    return res
            cnt = 0
            total = 0
            for row in rd:
                try:
                    ln = int(float(row[i_len]))
                    mt = int(float(row[i_mt]))
                    pt = int(float(row[i_pt]))
                    if pt > mt:
                        cnt += 1
                        total += ln
                except Exception:
                    continue
            res["pt_count"] = str(cnt)
            res["pt_bases"] = str(total)
            return res
    except Exception:
        return res


def parse_cov_csv(path):
    """Return mean/median if available in header, else blanks."""
    res = {"mean": "", "median": ""}
    if not path or not os.path.isfile(path):
        return res
    try:
        import csv

        with open(path, newline="", encoding="utf-8") as f:
            rd = csv.reader(f)
            hdr = next(rd, None)
            row = next(rd, None)
            if not hdr or not row:
                return res
            h = {k.strip().lower(): i for i, k in enumerate(hdr)}
            if "mean" in h:
                res["mean"] = row[h["mean"]]
            if "median" in h:
                res["median"] = row[h["median"]]
            return res
    except Exception:
        return res


def parse_contigs_stats(path):
    """
    Return dict {n,total_bp,N50,max_bp} from Flye contigs_stats.txt.
    Tolerant: if we can’t parse N50 from file, we recompute from lengths.
    """
    res = {"n": "", "total_bp": "", "N50": "", "max_bp": ""}
    if not path or not os.path.isfile(path):
        return res
    try:
        lengths = []
        with open(path, encoding="utf-8") as f:
            for i, ln in enumerate(f):
                if i == 0 and any(x in ln.lower() for x in ("length", "len", "contig")):
                    continue
                parts = ln.strip().split()
                if not parts:
                    continue
                try:
                    L = int(float(parts[0]))
                    if L > 0:
                        lengths.append(L)
                except Exception:
                    continue
        if not lengths:
            return res
        lengths.sort(reverse=True)
        total = sum(lengths)
        n50 = ""
        half = total / 2
        acc = 0
        for L in lengths:
            acc += L
            if acc >= half:
                n50 = L
                break
        res["n"] = str(len(lengths))
        res["total_bp"] = str(total)
        res["N50"] = str(n50)
        res["max_bp"] = str(max(lengths))
        return res
    except Exception:
        return res


def parse_gfa_S(path):
    """Fallback: compute n,total_bp,N50,max_bp from GFA S records with concrete sequence."""
    res = {"n": "", "total_bp": "", "N50": "", "max_bp": ""}
    if not path or not os.path.isfile(path):
        return res
    try:
        lengths = []
        with open(path, encoding="utf-8") as f:
            for ln in f:
                if not ln or ln[0] != "S":
                    continue
                parts = ln.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                seq = parts[2]
                if seq == "*":  # no concrete sequence
                    continue
                L = len(seq)
                if L > 0:
                    lengths.append(L)
        if not lengths:
            return res
        lengths.sort(reverse=True)
        total = sum(lengths)
        half = total / 2
        acc = 0
        n50 = ""
        for L in lengths:
            acc += L
            if acc >= half:
                n50 = L
                break
        res["n"] = str(len(lengths))
        res["total_bp"] = str(total)
        res["N50"] = str(n50)
        res["max_bp"] = str(max(lengths))
        return res
    except Exception:
        return res


def load_thresholds(path):
    """Parse simple key=value lines (case-insensitive) for mapq / aa length."""
    out = {"mapq": "", "aa_len": ""}
    if not path or not os.path.isfile(path):
        return out
    try:
        with open(path, encoding="utf-8") as f:
            for ln in f:
                ln = ln.strip()
                if not ln or ln.startswith("#"):
                    continue
                if "=" not in ln:
                    continue
                k, v = ln.split("=", 1)
                k = k.strip().lower()
                v = v.strip()
                if "mapq" in k and not out["mapq"]:
                    out["mapq"] = v
                if ("aa" in k or "aalen" in k) and not out["aa_len"]:
                    out["aa_len"] = v
        return out
    except Exception:
        return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--base-dir", required=True, help="Path to polap-readassemble/annotate-read-pt"
    )
    ap.add_argument("--out", required=True, help="Output JSON path")
    args = ap.parse_args()

    base = os.path.abspath(args.base_dir)
    pt0 = os.path.join(base, "pt0")
    pt1 = os.path.join(base, "pt1")

    # Files
    pt_id_txt = os.path.join(base, "pt.id.txt")
    pt_input = os.path.join(base, "pt.fq.seqkit.stats.ta.txt")
    pt_table = os.path.join(base, "pt-contig-annotation-depth-table.txt")
    pt_gene = os.path.join(base, "pt.gene.count")

    pre_cov_csv = os.path.join(pt1, "tmp", "pre.cov.csv")
    post_cov_csv = os.path.join(pt1, "tmp", "post.cov.csv")
    pre_cov_png = os.path.join(pt1, "tmp", "pre.cov.png")
    post_cov_png = os.path.join(pt1, "tmp", "post.cov.png")

    pt0_stats = os.path.join(pt0, "30-contigger", "contigs_stats.txt")
    pt1_stats = os.path.join(pt1, "30-contigger", "contigs_stats.txt")
    pt0_gfa = os.path.join(pt0, "30-contigger", "graph_final.gfa")
    pt1_gfa = os.path.join(pt1, "30-contigger", "graph_final.gfa")
    pt1_graph_gfa = os.path.join(pt1, "assembly_graph.gfa")
    pt0_png = os.path.join(pt0, "30-contigger", "graph_final.png")
    pt1_png = os.path.join(pt1, "30-contigger", "graph_final.png")
    pt1_graph_png = os.path.join(pt1, "assembly_graph.png")

    circ_paths = os.path.join(pt1, "ptdna", "ptdna", "circular_path_count.txt")
    scatter_pdf = os.path.join(base, "pt-contig-annotation-depth-table.txt.scatter.pdf")
    thresh_vars = os.path.join(base, "pt_thresh.vars")

    # KPIs
    # A) PT from table (PT>MT)
    pt_tab = parse_pt_table(pt_table)  # pt_count, pt_bases
    pt_count = pt_tab.get("pt_count", "")
    pt_bases = pt_tab.get("pt_bases", "")

    # B) Seqkit input
    seq_in = parse_seqkit_ta(pt_input)

    # C) Coverage
    cov_pre = parse_cov_csv(pre_cov_csv)
    cov_post = parse_cov_csv(post_cov_csv)

    # D) Genes and annotated length
    pt_gene_count = read_first_int(pt_gene)
    annotated_len = pt_bases  # by definition from PT>MT sum(Length)

    # E) Assembly stats
    def best_assembly_stats(stats_txt, gfa_main, gfa_alt=None):
        if os.path.isfile(stats_txt) and os.path.getsize(stats_txt) > 0:
            return parse_contigs_stats(stats_txt)
        if os.path.isfile(gfa_main) and os.path.getsize(gfa_main) > 0:
            return parse_gfa_S(gfa_main)
        if gfa_alt and os.path.isfile(gfa_alt) and os.path.getsize(gfa_alt) > 0:
            return parse_gfa_S(gfa_alt)
        return {"n": "", "total_bp": "", "N50": "", "max_bp": ""}

    asm_pt0 = best_assembly_stats(pt0_stats, pt0_gfa, None)
    asm_pt1 = best_assembly_stats(pt1_stats, pt1_gfa, pt1_graph_gfa)

    # F) Isoforms & thresholds
    circ = read_first_int(circ_paths)
    thr = load_thresholds(thresh_vars)

    # G) QC artifacts (pick best existing)
    graph_png = ""
    for cand in (pt1_png, pt1_graph_png, pt0_png):
        if os.path.isfile(cand):
            graph_png = cand
            break
    scatter = scatter_pdf if os.path.isfile(scatter_pdf) else ""
    pre_png = pre_cov_png if os.path.isfile(pre_cov_png) else ""
    post_png = post_cov_png if os.path.isfile(post_cov_png) else ""

    data = {
        "base_dir": base,
        "kpis": {
            "pt_reads": {"count": pt_count, "total_bases": pt_bases},
            "pt_input_seqkit": {
                "num_seqs": seq_in.get("num_seqs", ""),
                "sum_len": seq_in.get("sum_len", ""),
                "avg_len": seq_in.get("avg_len", ""),
                "N50": seq_in.get("N50", ""),
                "AvgQual": seq_in.get("AvgQual", ""),
                "GC_pct": seq_in.get("GC(%)", ""),
            },
            "coverage": {
                "pre": {
                    "mean": cov_pre.get("mean", ""),
                    "median": cov_pre.get("median", ""),
                },
                "post": {
                    "mean": cov_post.get("mean", ""),
                    "median": cov_post.get("median", ""),
                },
            },
            "genes": {
                "pt_gene_count": "" if pt_gene_count is None else str(pt_gene_count),
                "annotated_pt_length": annotated_len,
            },
            "assembly": {"pt0": asm_pt0, "pt1": asm_pt1},
            "isoforms": {"circular_path_count": "" if circ is None else str(circ)},
            "thresholds": {
                "mapq": thr.get("mapq", ""),
                "aa_len": thr.get("aa_len", ""),
            },
        },
        "qc": {
            "graph_png": graph_png,
            "pre_cov_png": pre_png,
            "post_cov_png": post_png,
            "scatter_pdf": scatter,
        },
    }

    # Ensure output directory exists
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)


if __name__ == "__main__":
    sys.exit(main())
