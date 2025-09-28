#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
polap-py-eval-numt-nupt-in-assembly.py
Evaluate inclusion of NUMT/NUPT in an organelle assembly using PAF alignments.

Inputs (two modes):
  A) Provide PAFs:
     --paf-nuc asm_vs_nuc.paf --paf-mt asm_vs_mt.paf --paf-pt asm_vs_pt.paf
  B) Or let the script run minimap2:
     --assembly assembly.fa --nuc nuclear.fa --mt mt.fa --pt pt.fa [--minimap2 minimap2] [-x asm5] [--threads 8]

Optional truth:
  --truth-bed spiked_nuc.inserts.bed   # BED of NUMT/NUPT intervals on the nuclear reference

Outputs:
  <out>.contigs.tsv     : per-contig stats and classification
  <out>.summary.tsv     : counts/lengths per label

Classification heuristics (tunable):
  --th-org-frac 0.60    : fraction of contig covered by mt/pt to call organelle
  --th-org-ident 0.90   : minimum best identity to mt/pt
  --th-nuc-minor 0.20   : max allowed nuclear fraction in an organelle call
  --th-split 0.30       : min fractions for split (organelle+nuclear) to flag NUMT/NUPT

Example (PAF given):
  polap-py-eval-numt-nupt-in-assembly.py \
    --assembly flye_spiked/assembly.fasta \
    --paf-nuc spiked_vs_nuc.paf --paf-mt spiked_vs_mt.paf --paf-pt spiked_vs_pt.paf \
    --nuc nuc_spiked.augmented.fa --mt mt.fa --pt pt.fa \
    --truth-bed nuc_spiked.inserts.bed --out eval_spiked

Example (auto-minimap2):
  polap-py-eval-numt-nupt-in-assembly.py \
    --assembly flye_clean/assembly.fasta \
    --nuc nuclear.fa --mt mt.fa --pt pt.fa \
    --out eval_clean --threads 16 -x asm5
"""
import sys, os, argparse, subprocess, shutil, tempfile, math
from collections import defaultdict, namedtuple, Counter
from typing import List, Dict, Tuple

Hit = namedtuple("Hit", "qname qlen tname tlen qstart qend tstart tend nmatch alen mapq strand ident")

def read_fa_sizes(path: str) -> Dict[str,int]:
    sizes: Dict[str,int] = {}
    with open(path) as fh:
        name = None; n = 0
        for ln in fh:
            if ln.startswith(">"):
                if name is not None: sizes[name] = n
                name = ln[1:].strip().split()[0]; n = 0
            else:
                n += len(ln.strip())
        if name is not None: sizes[name] = n
    return sizes

def parse_paf(path: str) -> List[Hit]:
    hits: List[Hit] = []
    with open(path) as fh:
        for ln in fh:
            if not ln.strip(): continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 12: continue
            qname, qlen, qstart, qend, strand, tname, tlen, tstart, tend, nmatch, alen, mapq = f[:12]
            try:
                qlen = int(qlen); qstart = int(qstart); qend = int(qend)
                tlen = int(tlen); tstart = int(tstart); tend = int(tend)
                nmatch = int(nmatch); alen = int(alen); mapq = int(mapq)
            except ValueError:
                continue
            ident = nmatch / alen if alen > 0 else 0.0
            hits.append(Hit(qname, qlen, tname, tlen, qstart, qend, tstart, tend, nmatch, alen, mapq, strand, ident))
    return hits

def merge_intervals(iv: List[Tuple[int,int]]) -> List[Tuple[int,int]]:
    if not iv: return []
    iv = sorted(iv)
    out = [list(iv[0])]
    for s,e in iv[1:]:
        if s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s,e])
    return [(a,b) for a,b in out]

def covered_length(iv: List[Tuple[int,int]]) -> int:
    return sum(b-a for a,b in merge_intervals(iv))

def load_truth_bed(path: str):
    # chrom start end name (score) strand
    truth = defaultdict(lambda: {"NUMT": [], "NUPT": [], "RAW": []})
    if not path: return truth
    with open(path) as fh:
        for ln in fh:
            if not ln.strip() or ln.startswith("#"): continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 6: continue
            chrom, start, end, name, score, strand = f[:6]
            try:
                start = int(start); end = int(end)
            except:
                continue
            tag = "NUMT" if name.upper().startswith("NUMT") else ("NUPT" if name.upper().startswith("NUPT") else "UNK")
            truth[chrom]["RAW"].append((start,end,name))
            if tag in ("NUMT","NUPT"):
                truth[chrom][tag].append((start,end,name))
    return truth

def intersect_bp(iv: List[Tuple[int,int]], truth_list: List[Tuple[int,int,str]]) -> int:
    bp = 0
    for a,b in iv:
        for s,e,_ in truth_list:
            x = max(a,s); y = min(b,e)
            if x < y: bp += (y-x)
    return bp

def summarize_by_target(hits: List[Hit]) -> Dict[str,Dict]:
    per_t = defaultdict(list)
    for h in hits:
        per_t[h.tname].append(h)
    out = {}
    for t, lst in per_t.items():
        qlens = set(h.qlen for h in lst)
        qlen = max(qlens) if qlens else 0
        iv_q = [(h.qstart, h.qend) for h in lst]
        iv_t = [(min(h.tstart,h.tend), max(h.tstart,h.tend)) for h in lst]
        cov_bp = covered_length(iv_q)
        best_ident = max((h.ident for h in lst), default=0.0)
        out[t] = {"qlen": qlen, "cov_bp": cov_bp, "cov_frac": cov_bp/qlen if qlen else 0.0,
                  "best_ident": best_ident, "iv_q": merge_intervals(iv_q), "iv_t": merge_intervals(iv_t)}
    return out

def run_minimap2(minimap2: str, preset: str, threads: int, ref: str, qry: str, out_paf: str):
    cmd = [minimap2, "-x", preset, "-t", str(threads), ref, qry]
    with open(out_paf, "w") as fo:
        subprocess.run(cmd, check=True, stdout=fo)

def classify(qlen: int, nuc_cov: int, mt_cov: int, pt_cov: int,
             best_nuc: float, best_mt: float, best_pt: float,
             th_org_frac: float, th_org_ident: float, th_nuc_minor: float, th_split: float):
    fn = nuc_cov/qlen if qlen else 0.0
    fm = mt_cov/qlen  if qlen else 0.0
    fp = pt_cov/qlen  if qlen else 0.0
    label, reason = "unclassified", []
    if fm >= th_org_frac and best_mt >= th_org_ident and fn < th_nuc_minor and fp < th_nuc_minor:
        label = "organelle_mt"; reason.append("mt-major")
    elif fp >= th_org_frac and best_pt >= th_org_ident and fn < th_nuc_minor and fm < th_nuc_minor:
        label = "organelle_pt"; reason.append("pt-major")
    elif fm >= th_split and fn >= th_split and best_mt >= 0.85 and best_nuc >= 0.85:
        label = "NUMT_candidate"; reason.append("split mt+nuc")
    elif fp >= th_split and fn >= th_split and best_pt >= 0.85 and best_nuc >= 0.85:
        label = "NUPT_candidate"; reason.append("split pt+nuc")
    elif fn >= 0.8:
        label = "nuclear"; reason.append("nuc-major")
    elif fm >= 0.4 and fp >= 0.4:
        label = "chimeric_mt_pt"; reason.append("mt+pt")
    return label, ";".join(reason), dict(frac_nuc=fn, frac_mt=fm, frac_pt=fp)

def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Evaluate NUMT/NUPT inclusion in organelle assemblies.")
    # Assembly & references
    ap.add_argument("--assembly", required=True, help="assembly FASTA (contigs)")
    ap.add_argument("--nuc", help="nuclear reference FASTA (augmented if spiked)")
    ap.add_argument("--mt",  help="mitochondrial FASTA")
    ap.add_argument("--pt",  help="plastid FASTA")
    # PAFs (optional)
    ap.add_argument("--paf-nuc", help="PAF of assembly vs nuclear")
    ap.add_argument("--paf-mt",  help="PAF of assembly vs mt")
    ap.add_argument("--paf-pt",  help="PAF of assembly vs pt")
    # Minimap2 options if PAFs absent
    ap.add_argument("--minimap2", default="minimap2", help="path to minimap2")
    ap.add_argument("-x", "--preset", default="asm5", help="minimap2 preset (asm5 recommended for contigs)")
    ap.add_argument("--threads", type=int, default=8)
    # Truth
    ap.add_argument("--truth-bed", help="BED of NUMT/NUPT regions on nuclear reference (optional)")
    # Thresholds
    ap.add_argument("--th-org-frac", type=float, default=0.60)
    ap.add_argument("--th-org-ident", type=float, default=0.90)
    ap.add_argument("--th-nuc-minor", type=float, default=0.20)
    ap.add_argument("--th-split", type=float, default=0.30)
    # Output
    ap.add_argument("--out", required=True, help="output prefix")
    args = ap.parse_args()

    # map if needed
    tmpdir = None
    try:
        if not (args.paf_nuc and args.paf_mt and args.paf_pt):
            # need refs
            for ref in ("nuc","mt","pt"):
                if getattr(args, ref) is None:
                    raise SystemExit(f"--{ref} is required when PAFs are not supplied")
            tmpdir = tempfile.mkdtemp(prefix="eval_numt_")
            paf_nuc = args.paf_nuc or os.path.join(tmpdir, "asm_vs_nuc.paf")
            paf_mt  = args.paf_mt  or os.path.join(tmpdir, "asm_vs_mt.paf")
            paf_pt  = args.paf_pt  or os.path.join(tmpdir, "asm_vs_pt.paf")
            run_minimap2(args.minimap2, args.preset, args.threads, args.nuc, args.assembly, paf_nuc)
            run_minimap2(args.minimap2, args.preset, args.threads, args.mt,  args.assembly, paf_mt)
            run_minimap2(args.minimap2, args.preset, args.threads, args.pt,  args.assembly, paf_pt)
        else:
            paf_nuc, paf_mt, paf_pt = args.paf_nuc, args.paf_mt, args.paf_pt

        # parse PAFs
        nuc_hits = parse_paf(paf_nuc) if os.path.exists(paf_nuc) else []
        mt_hits  = parse_paf(paf_mt)  if os.path.exists(paf_mt)  else []
        pt_hits  = parse_paf(paf_pt)  if os.path.exists(paf_pt)  else []

        # index by query
        per_q = defaultdict(lambda: {"nuc": [], "mt": [], "pt": []})
        for h in nuc_hits: per_q[h.qname]["nuc"].append(h)
        for h in mt_hits:  per_q[h.qname]["mt"].append(h)
        for h in pt_hits:  per_q[h.qname]["pt"].append(h)

        # truth
        truth = load_truth_bed(args.truth_bed) if args.truth_bed else None

        # sizes
        asm_sizes = read_fa_sizes(args.assembly)

        rows = []
        for q, d in per_q.items():
            qlen = asm_sizes.get(q, max([h.qlen for v in d.values() for h in v], default=0))
            # summarize per target set
            nuc_sum = summarize_by_target(d["nuc"])
            mt_sum  = summarize_by_target(d["mt"])
            pt_sum  = summarize_by_target(d["pt"])

            cov_nuc = sum(v["cov_bp"] for v in nuc_sum.values())
            cov_mt  = sum(v["cov_bp"] for v in mt_sum.values())
            cov_pt  = sum(v["cov_bp"] for v in pt_sum.values())
            best_nuc = max((v["best_ident"] for v in nuc_sum.values()), default=0.0)
            best_mt  = max((v["best_ident"] for v in mt_sum.values()),  default=0.0)
            best_pt  = max((v["best_ident"] for v in pt_sum.values()),  default=0.0)

            label, reason, fr = classify(qlen, cov_nuc, cov_mt, cov_pt,
                                         best_nuc, best_mt, best_pt,
                                         args.th_org_frac, args.th_org_ident, args.th_nuc_minor, args.th_split)

            # truth overlaps on nuclear target coordinates
            numt_bp = nupt_bp = 0
            if truth and d["nuc"]:
                # group by nuclear tname
                per_tn = defaultdict(list)
                for h in d["nuc"]:
                    a = min(h.tstart,h.tend); b = max(h.tstart,h.tend)
                    per_tn[h.tname].append((a,b))
                for tn, iv in per_tn.items():
                    if tn in truth:
                        ivm = merge_intervals(iv)
                        numt_bp += intersect_bp(ivm, truth[tn]["NUMT"])  # type: ignore
                        nupt_bp += intersect_bp(ivm, truth[tn]["NUPT"])  # type: ignore

            rows.append({
                "contig": q, "len": qlen,
                "cov_nuc": cov_nuc, "cov_mt": cov_mt, "cov_pt": cov_pt,
                "frac_nuc": fr["frac_nuc"], "frac_mt": fr["frac_mt"], "frac_pt": fr["frac_pt"],
                "best_nuc": best_nuc, "best_mt": best_mt, "best_pt": best_pt,
                "label": label, "reason": reason,
                "numt_bp_overlap": numt_bp, "nupt_bp_overlap": nupt_bp
            })

        # write per-contig
        out_tsv = f"{args.out}.contigs.tsv"
        with open(out_tsv, "w") as fo:
            fo.write("\t".join(["contig","len","cov_nuc","cov_mt","cov_pt","frac_nuc","frac_mt","frac_pt",
                                "best_nuc","best_mt","best_pt","label","reason","numt_bp_overlap","nupt_bp_overlap"]) + "\n")
            for r in sorted(rows, key=lambda x: (x["label"], -x["len"], x["contig"])):
                fo.write("\t".join([
                    r["contig"], str(r["len"]),
                    str(r["cov_nuc"]), str(r["cov_mt"]), str(r["cov_pt"]),
                    f'{r["frac_nuc"]:.4f}', f'{r["frac_mt"]:.4f}', f'{r["frac_pt"]:.4f}',
                    f'{r["best_nuc"]:.4f}', f'{r["best_mt"]:.4f}', f'{r["best_pt"]:.4f}',
                    r["label"], r["reason"],
                    str(r["numt_bp_overlap"]), str(r["nupt_bp_overlap"])
                ]) + "\n")

        # summary
        by_label = Counter(r["label"] for r in rows)
        len_by_label = Counter()
        for r in rows: len_by_label[r["label"]] += r["len"]
        out_sum = f"{args.out}.summary.tsv"
        with open(out_sum, "w") as fo:
            fo.write("label\tcount\ttotal_len\n")
            for lab in sorted(by_label.keys()):
                fo.write(f"{lab}\t{by_label[lab]}\t{len_by_label[lab]}\n")

        sys.stderr.write(f"[ok] wrote {out_tsv}\n[ok] wrote {out_sum}\n")

    finally:
        if tmpdir and os.path.isdir(tmpdir):
            shutil.rmtree(tmpdir, ignore_errors=True)

if __name__ == "__main__":
    main()
