#!/usr/bin/env python3
# scripts/report_readassemble_mt.py
# Version: v0.1.0
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Build a JSON-only report of the MT pipeline under an mtseed base dir.
#
# JSON structure:
# {
#   "base_dir": ".../mtseed",
#   "steps": {
#     "panel": {
#       "ir": {"pick_tsv": ".../01-panel/ir.pick.tsv", "summary": [...]},
#       "isomers": {"A_len": "...", "B_len": "...", "doubleA_len": "...", "doubleA_ok": true, ...},
#       "log": "01-panel/log/pt_isoform.log"
#     },
#     "map": {
#       "all_ids": {"count": "..."},
#       "pt_ids": {"count": "..."},
#       "removed_pt_frac": "...",
#       "nonpt_ids": {"file": ".../keep.nonpt.ids", "count": "..."},
#       "nonpt_seqkit": {"num_seqs":"", "sum_len":"", "avg_len":"", "N50":"", "AvgQual":""},
#       "thresholds": {"vars": ".../pt_thresh.vars"},
#       "diag": {"pt_thresh_diag": ".../pt_thresh.diag.tsv"}
#     },
#     "allvsall": {
#       "qc": {
#         "degree_hist_pdf": ".../03-allvsall/04-qc/scan_degree_hist.pdf",
#         "wdegree_hist_pdf":".../scan_wdegree_hist.pdf",
#         "cum_wdegree_pdf":".../scan_cum_wdegree.pdf"
#       },
#       "vars": {"overlap_qc_vars": ".../04-qc/overlap_qc.vars"},
#       "overlapness": {"loose_tsv":".../overlapness.tsv", "strict_tsv":".../overlapness_strict.tsv"}
#     },
#     "busco": {
#       "selected_ids": {"file":".../05-round/select_ids.txt", "count":"..."}
#     },
#     "miniasm": {
#       "seeds_fa":".../06-miniasm/m_seeds_raw.fa",
#       "gfa":".../06-miniasm/miniasm.gfa",
#       "stats": {"S": "...", "total_bp":"...", "N50":"...", "max_bp":"..."}
#     },
#     "iterations": [
#        {"name":"mt0", "paths":{...}, "stats":{...}, "genes":{...}},
#        {"name":"mt1", ...}, {"name":"mt2", ...}, {"name":"mt3", ...}
#     ]
#   },
#   "qc": { "panel_log": "...", "allvsall_pdfs": {...} }
# }
#
import os
import sys
import json
import argparse

def read_first_int(path):
    if not path or not os.path.isfile(path):
        return None
    try:
        with open(path, encoding="utf-8") as f:
            ln = f.readline().strip()
            if not ln: return None
            tok = ln.split()[0]
            return int(tok)
    except Exception:
        return None

def count_lines(path):
    if not path or not os.path.isfile(path):
        return None
    try:
        n=0
        with open(path, encoding="utf-8") as f:
            for ln in f:
                if ln.strip(): n+=1
        return n
    except Exception:
        return None

def fasta_total_len(path):
    """Return total sequence length and number of records for a (possibly medium) FASTA."""
    if not path or not os.path.isfile(path):
        return {"n":None,"total_bp":None}
    n=0; total=0; inseq=False
    try:
        with open(path, encoding="utf-8") as f:
            for ln in f:
                if not ln: continue
                if ln.startswith(">"):
                    n+=1; inseq=True
                else:
                    # strip whitespace to guard
                    total += len(ln.strip())
        return {"n":n or None, "total_bp": total or None}
    except Exception:
        return {"n":None,"total_bp":None}

def parse_gfa_S(path):
    """Compute S count, total_bp, N50, max_bp from GFA S records with explicit sequence."""
    res={"S":None,"total_bp":None,"N50":None,"max_bp":None}
    if not path or not os.path.isfile(path):
        return res
    try:
        lens=[]
        with open(path, encoding="utf-8") as f:
            for ln in f:
                if not ln or ln[0]!="S": continue
                parts = ln.rstrip("\n").split("\t")
                if len(parts)<3: continue
                seq = parts[2]
                if seq=="*":  # no explicit sequence
                    continue
                lens.append(len(seq))
        if not lens: return res
        lens.sort(reverse=True)
        total=sum(lens)
        half=total/2; acc=0; n50=None
        for L in lens:
            acc+=L
            if acc>=half:
                n50=L
                break
        res["S"]=len(lens); res["total_bp"]=total; res["N50"]=n50; res["max_bp"]=max(lens)
        return res
    except Exception:
        return res

def parse_seqkit_t(path):
    """Parse a seqkit -T-like table (header + one row)."""
    res={"num_seqs":"","sum_len":"","avg_len":"","N50":"","AvgQual":""}
    if not path or not os.path.isfile(path):
        return res
    try:
        with open(path, encoding="utf-8") as f:
            lines=[ln.rstrip("\n\r") for ln in f]
        if len(lines)<2: return res
        hdr = lines[0].split("\t")
        row = lines[1].split("\t")
        h  = {k:i for i,k in enumerate(hdr)}
        def get(k):
            i=h.get(k); 
            return row[i] if i is not None and i<len(row) else ""
        res["num_seqs"]=get("num_seqs")
        res["sum_len"]= get("sum_len")
        res["avg_len"]= get("avg_len")
        res["N50"]=     get("N50")
        res["AvgQual"]= get("AvgQual")
        return res
    except Exception:
        return res

def read_ir_pick(path):
    """
    Parse 01-panel/ir.pick.tsv (if present). We report lines as a raw list of dicts.
    Expect columns e.g., name, start, end, len, etc. We pass through as-is up to a few rows.
    """
    if not path or not os.path.isfile(path):
        return []
    try:
        import csv
        out=[]
        with open(path, newline="", encoding="utf-8") as f:
            rd=csv.DictReader(f, delimiter="\t")
            for i,row in enumerate(rd):
                out.append(row)
                if i>=9: break  # limit
        return out
    except Exception:
        return []

def fasta_len(path):
    """Return sum length of all sequences in a fasta (small/medium files)."""
    info = fasta_total_len(path)
    return info.get("total_bp")

def build_iteration(mt_root, name):
    """
    Build stats for iteration mt0/mt1/mt2/mt3.
    Prefer 30-contigger/graph_final.fasta; else assembly_graph.gfa S; else 30-contigger/graph_final.gfa S.
    Include gene counts from 50-annotation/mt.gene.count and annotated total length from contig-annotation-depth-table.txt.
    """
    it_dir = os.path.join(mt_root, name)
    if not os.path.isdir(it_dir):
        return {"name":name, "paths":{}, "stats":{}, "genes":{}}

    contig_fa = os.path.join(it_dir, "30-contigger", "graph_final.fasta")
    graph_png = os.path.join(it_dir, "30-contigger", "graph_final.png")
    gfa_gfa   = os.path.join(it_dir, "30-contigger", "graph_final.gfa")
    asmg_gfa  = os.path.join(it_dir, "assembly_graph.gfa")
    gene_file = os.path.join(it_dir, "50-annotation", "mt.gene.count")
    contig_ann= os.path.join(it_dir, "contig-annotation-depth-table.txt")

    # stats
    stats = {}
    if os.path.isfile(contig_fa):
        info=fasta_total_len(contig_fa)
        # compute N50 quickly from fasta lengths again
        lens=[]
        try:
            with open(contig_fa, encoding="utf-8") as f:
                L=0
                for ln in f:
                    if ln.startswith(">"):
                        if L>0: lens.append(L); L=0
                    else:
                        L+=len(ln.strip())
                if L>0: lens.append(L)
            if lens:
                lens.sort(reverse=True)
                total=sum(lens); half=total/2; acc=0; n50=None
                for L in lens:
                    acc+=L
                    if acc>=half: n50=L; break
                stats={"n":len(lens), "total_bp":total, "N50":n50, "max_bp":max(lens)}
        except Exception:
            pass
    elif os.path.isfile(asmg_gfa):
        g=parse_gfa_S(asmg_gfa); 
        stats={"n":g.get("S"), "total_bp":g.get("total_bp"), "N50":g.get("N50"), "max_bp":g.get("max_bp")}
    elif os.path.isfile(gfa_gfa):
        g=parse_gfa_S(gfa_gfa); 
        stats={"n":g.get("S"), "total_bp":g.get("total_bp"), "N50":g.get("N50"), "max_bp":g.get("max_bp")}
    else:
        stats={"n":None,"total_bp":None,"N50":None,"max_bp":None}

    # genes
    genes={}
    mtg = read_first_int(gene_file)
    if mtg is not None: genes["mt_gene_count"]=mtg
    # annotated length
    if os.path.isfile(contig_ann):
        try:
            import csv
            total=0; cnt=0
            with open(contig_ann, newline="", encoding="utf-8") as f:
                rd=csv.reader(f, delimiter="\t")
                hdr=next(rd,None)
                if hdr:
                    pos={k.strip().lower():i for i,k in enumerate(hdr)}
                    i_len=pos.get("length")
                    if i_len is not None:
                        for row in rd:
                            try:
                                L=float(row[i_len]); 
                                if L>0: total+=int(L); cnt+=1
                            except Exception:
                                continue
            genes["annotated_mt_length"]= total if total>0 else None
            genes["annotated_rows"]= cnt if cnt>0 else None
        except Exception:
            pass

    return {
        "name": name,
        "paths": {
            "graph_final_fasta": contig_fa if os.path.isfile(contig_fa) else "",
            "graph_png": graph_png if os.path.isfile(graph_png) else "",
            "assembly_graph_gfa": asmg_gfa if os.path.isfile(asmg_gfa) else (gfa_gfa if os.path.isfile(gfa_gfa) else ""),
        },
        "stats": stats,
        "genes": genes
    }

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base-dir", required=True, help="Path to .../polap-readassemble/mtseed")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    base = os.path.abspath(args.base_dir)
    out_dir = os.path.dirname(os.path.abspath(args.out))
    os.makedirs(out_dir, exist_ok=True)

    # STEP Panel (01-panel)
    panel_dir = os.path.join(base, "01-panel")
    ir_tsv = os.path.join(panel_dir, "ir.pick.tsv")
    iso_log= os.path.join(panel_dir, "log", "pt_isoform.log")

    isomerA = os.path.join(panel_dir, "pt_isomerA.fa")
    isomerB = os.path.join(panel_dir, "pt_isomerB.fa")
    doubleA = os.path.join(panel_dir, "pt_isomerA.double.fa")
    doubleB = os.path.join(panel_dir, "pt_isomerB.double.fa")

    def_len_A = fasta_len(isomerA)
    def_len_B = fasta_len(isomerB)
    dbl_len_A = fasta_len(doubleA)
    dbl_len_B = fasta_len(doubleB)

    panel = {
        "ir": {
            "pick_tsv": ir_tsv if os.path.isfile(ir_tsv) else "",
            "summary": read_ir_pick(ir_tsv)
        },
        "isomers": {
            "A_len": def_len_A,
            "B_len": def_len_B,
            "doubleA_len": dbl_len_A,
            "doubleA_ok": (dbl_len_A == (2 * def_len_A)) if (dbl_len_A and def_len_A) else None,
            "doubleB_len": dbl_len_B,
            "doubleB_ok": (dbl_len_B == (2 * def_len_B)) if (dbl_len_B and def_len_B) else None
        },
        "log": iso_log if os.path.isfile(iso_log) else ""
    }

    # STEP map (02-map)
    map_dir = os.path.join(base, "02-map")
    all_ids = os.path.join(map_dir, "all.ids")
    pt_ids  = os.path.join(map_dir, "pt.ids.sorted")

    all_count = count_lines(all_ids)
    pt_count  = count_lines(pt_ids)

    removed_frac = None
    if all_count and pt_count is not None:
        try: removed_frac = round(int(pt_count)/int(all_count), 6)
        except Exception: removed_frac = None

    keep_nonpt = os.path.join(base, "keep.nonpt.ids")
    keep_nonpt_count = count_lines(keep_nonpt)

    nonpt_seqkit = os.path.join(base, "reads.nonpt.fq.gz.seqkit.stats.T.txt")
    nonpt_stats = parse_seqkit_t(nonpt_seqkit)

    pt_thresh_vars = os.path.join(base, "pt_thresh.vars")
    pt_thresh_diag = os.path.join(base, "pt_thresh.diag.tsv")

    mapstep = {
        "all_ids": {"file": all_ids if os.path.isfile(all_ids) else "", "count": all_count},
        "pt_ids":  {"file": pt_ids if os.path.isfile(pt_ids) else "",   "count": pt_count},
        "removed_pt_frac": removed_frac,
        "nonpt_ids": {"file": keep_nonpt if os.path.isfile(keep_nonpt) else "", "count": keep_nonpt_count},
        "nonpt_seqkit": nonpt_stats,
        "thresholds": {"vars": pt_thresh_vars if os.path.isfile(pt_thresh_vars) else ""},
        "diag": {"pt_thresh_diag": pt_thresh_diag if os.path.isfile(pt_thresh_diag) else ""}
    }

    # STEP allvsall (03-allvsall)
    a_dir = os.path.join(base, "03-allvsall")
    qc_dir = os.path.join(a_dir, "04-qc")
    allvsall = {
        "qc": {
            "degree_hist_pdf": os.path.join(qc_dir, "scan_degree_hist.pdf") if os.path.isfile(os.path.join(qc_dir,"scan_degree_hist.pdf")) else "",
            "wdegree_hist_pdf": os.path.join(qc_dir, "scan_wdegree_hist.pdf") if os.path.isfile(os.path.join(qc_dir,"scan_wdegree_hist.pdf")) else "",
            "cum_wdegree_pdf": os.path.join(qc_dir, "scan_cum_wdegree.pdf") if os.path.isfile(os.path.join(qc_dir,"scan_cum_wdegree.pdf")) else ""
        },
        "vars": {"overlap_qc_vars": os.path.join(qc_dir, "overlap_qc.vars") if os.path.isfile(os.path.join(qc_dir,"overlap_qc.vars")) else ""},
        "overlapness": {
            "loose_tsv": os.path.join(a_dir, "overlapness.tsv") if os.path.isfile(os.path.join(a_dir,"overlapness.tsv")) else "",
            "strict_tsv": os.path.join(a_dir, "overlapness_strict.tsv") if os.path.isfile(os.path.join(a_dir,"overlapness_strict.tsv")) else ""
        }
    }

    # STEP busco (05-round)
    b_dir = os.path.join(base, "05-round")
    selected_ids = os.path.join(b_dir, "select_ids.txt")
    busco = {
        "selected_ids": {
            "file": selected_ids if os.path.isfile(selected_ids) else "",
            "count": count_lines(selected_ids)
        }
    }

    # STEP miniasm (06-miniasm)
    m_dir = os.path.join(base, "06-miniasm")
    seeds_fa = os.path.join(m_dir, "m_seeds_raw.fa")
    miniasm_gfa = os.path.join(m_dir, "miniasm.gfa")
    miniasm_stats = parse_gfa_S(miniasm_gfa) if os.path.isfile(miniasm_gfa) else {"S":None,"total_bp":None,"N50":None,"max_bp":None}
    miniasm = {
        "seeds_fa": seeds_fa if os.path.isfile(seeds_fa) else "",
        "gfa": miniasm_gfa if os.path.isfile(miniasm_gfa) else "",
        "stats": {
            "S": miniasm_stats.get("S"),
            "total_bp": miniasm_stats.get("total_bp"),
            "N50": miniasm_stats.get("N50"),
            "max_bp": miniasm_stats.get("max_bp")
        }
    }

    # STEP iterations mt0..mt3
    its=[]
    for name in ("mt0","mt1","mt2","mt3"):
        its.append(build_iteration(base, name))

    data = {
        "base_dir": base,
        "steps": {
            "panel": panel,
            "map": mapstep,
            "allvsall": allvsall,
            "busco": busco,
            "miniasm": miniasm,
            "iterations": its
        }
    }

    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)

if __name__ == "__main__":
    sys.exit(main())
