#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
polap-py-oatk-core-judge-assembly.py v0.0.6

Judge a syncasm k-attempt using GFA(s). Robust under `set -euo pipefail`:
- Always exits 0 (stdout: PASS|CONTINUE; stderr: one-line SUMMARY)
- Accepts GFA-only inputs:
    --gfa-final    cleaned graph (.utg.final.gfa[.gz])   ← used for decision
    --gfa-utg      raw graph (.utg.gfa[.gz])             ← diagnostics
- Auto-extracts unitigs FASTA from whichever GFA has S-sequences (prefers final)
- Optional coverage evenness (minimap2) on ~N sampled reads (non-fatal)
- SUMMARY includes:
    n,total,max,N50,tangle_final,circ,from,final_nodes,final_branch,utg_nodes,utg_branch,reasons
- Optional outputs:
    --diag-file  : write SUMMARY line to file
    --tsv-file   : append a clean TSV row (reasons LAST)
    --quiet-diag : suppress JSON block (keep SUMMARY)

Stdout : PASS | CONTINUE
Stderr : one-line SUMMARY (+ optional JSON)
Exit   : 0 (defensive; never kills the ladder)
"""

import argparse
import sys
import os
import gzip
import shutil
import subprocess
import json
import math
import tempfile
import traceback
import signal
from collections import Counter, defaultdict

VERSION = "0.0.6"


# ---------------- utilities ----------------
def file_ok(path: str) -> bool:
    return bool(path) and os.path.isfile(path) and os.path.getsize(path) > 0


def have_tool(name: str) -> bool:
    return shutil.which(name) is not None


def open_gz(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def safe_print_err(msg: str):
    sys.stderr.write(msg + "\n")
    sys.stderr.flush()


def safe_json_err(obj):
    try:
        safe_print_err(json.dumps(obj, indent=2))
    except Exception as e:
        safe_print_err(f'{{"json_dump_error":"{type(e).__name__}:{e}"}}')


def ignore_sigpipe():
    try:
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    except Exception:
        pass


# ---------------- FASTA & GFA stats ----------------
def fa_stats(path):
    if not file_ok(path):
        return dict(n=0, total=0, max=0, N50=0, names=[])
    n = 0
    total = 0
    lens = []
    names = []
    buf = []
    name = None
    with open_gz(path) as fh:
        for ln in fh:
            if ln.startswith(">"):
                if buf:
                    L = sum(map(len, buf))
                    total += L
                    lens.append(L)
                    buf = []
                    names.append(name)
                name = ln[1:].strip().split()[0]
                n += 1
            else:
                buf.append(ln.strip())
        if buf:
            L = sum(map(len, buf))
            total += L
            lens.append(L)
            names.append(name)

    def n50(vals):
        if not vals:
            return 0
        vals_sorted = sorted(vals, reverse=True)
        half = sum(vals_sorted) / 2
        s = 0
        for L in vals_sorted:
            s += L
            if s >= half:
                return L
        return 0

    return dict(
        n=n, total=total, max=(max(lens) if lens else 0), N50=n50(lens), names=names
    )


def gfa_degrees(path):
    if not file_ok(path):
        return dict(nodes=0, nbranch=0, tips=0)
    deg = Counter()
    nodes = set()
    with open_gz(path) as fh:
        for ln in fh:
            if ln.startswith("S"):
                nodes.add(ln.split("\t", 2)[1])
            elif ln.startswith("L"):
                f = ln.rstrip("\n").split("\t")
                if len(f) >= 4:
                    deg[f[1]] += 1
                    deg[f[3]] += 1
    tips = sum(1 for n in nodes if deg[n] == 1)
    nbranch = sum(1 for n in nodes if deg[n] > 2)
    return dict(nodes=len(nodes), nbranch=nbranch, tips=tips)


def extract_fasta_from_gfa(gfa_path: str, out_fa: str) -> bool:
    """Write FASTA from GFA S-lines with sequences. Return True if >=1 sequence written."""
    if not file_ok(gfa_path):
        return False
    wrote = 0
    try:
        with open_gz(gfa_path) as fi, open(out_fa, "w") as fo:
            for ln in fi:
                if ln.startswith("S\t"):
                    parts = ln.rstrip("\n").split("\t")
                    if len(parts) >= 3 and parts[2] not in ("*", ""):
                        fo.write(">" + parts[1] + "\n" + parts[2] + "\n")
                        wrote += 1
        if wrote == 0:
            try:
                os.remove(out_fa)
            except Exception:
                pass
            return False
        return True
    except Exception:
        try:
            os.remove(out_fa)
        except Exception:
            pass
        return False


# ---------------- circularity heuristic ----------------
def head_tail_overlap(seq, min_ovlp):
    if len(seq) < min_ovlp:
        return 0, 0.0
    K = 31
    head = seq[:min_ovlp]
    tail = seq[-min_ovlp:]
    pos = head.find(tail[:K])
    if pos < 0:
        return 0, 0.0
    window = min_ovlp
    matches = sum(1 for a, b in zip(head[pos : pos + window], tail[:window]) if a == b)
    return window, matches / window


def best_circular_candidate(unitigs_fa, min_ovlp, min_ident):
    if not file_ok(unitigs_fa):
        return None
    seqs = []
    buf = []
    name = None
    with open_gz(unitigs_fa) as fh:
        for ln in fh:
            if ln.startswith(">"):
                if name is not None:
                    seqs.append(("".join(buf), name))
                name = ln[1:].strip().split()[0]
                buf = []
            else:
                buf.append(ln.strip())
        if buf:
            seqs.append(("".join(buf), name))
    seqs.sort(key=lambda x: len(x[0]), reverse=True)
    for s, nm in seqs[:5]:
        ov, ident = head_tail_overlap(s, min_ovlp)
        if ov >= min_ovlp and ident >= min_ident:
            return dict(name=nm, ovlp=ov, ident=round(ident, 4))
    return None


# ---------------- coverage evenness (optional) ----------------
def sample_reads(in_fastx, out_fastq, sample_n):
    if have_tool("seqkit"):
        tmp = out_fastq + ".tmp"
        subprocess.run(
            ["seqkit", "sample", "-n", str(sample_n), "-o", tmp, in_fastx], check=True
        )
        with open(tmp, "rt") as fh:
            first = fh.readline().strip()
        if first.startswith("@"):
            shutil.move(tmp, out_fastq)
        else:
            with open(tmp, "rt") as fi, open(out_fastq, "wt") as fo:
                hdr = None
                for ln in fi:
                    if ln.startswith(">"):
                        hdr = "@" + ln[1:].strip()
                        fo.write(hdr + "\n")
                    else:
                        seq = ln.strip()
                        if not seq:
                            continue
                        fo.write(seq + "\n+\n" + "I" * len(seq) + "\n")
            os.remove(tmp)
        return
    # fallback: first N
    readN = 0
    with open_gz(in_fastx) as fi, open(out_fastq, "wt") as fo:
        if in_fastx.endswith((".fq", ".fastq", ".fq.gz", ".fastq.gz")):
            i = 0
            buf = ["", "", "", ""]
            for ln in fi:
                buf[i % 4] = ln.rstrip("\n")
                i += 1
                if i % 4 == 0:
                    fo.write("\n".join(buf) + "\n")
                    readN += 1
                    if readN >= sample_n:
                        break
        else:
            name = None
            for ln in fi:
                if ln.startswith(">"):
                    name = ln.strip().replace(">", "@")
                else:
                    seq = ln.strip()
                    if not seq:
                        continue
                    fo.write(f"{name}\n{seq}\n+\n{'I'*len(seq)}\n")
                    readN += 1
                    if readN >= sample_n:
                        break


def run_minimap2(unitigs_fa, reads_fastq, threads, paf_out):
    cmd = [
        "minimap2",
        "-x",
        "map-ont",
        "--secondary=no",
        "-N",
        "1",
        "-t",
        str(threads),
        unitigs_fa,
        reads_fastq,
    ]
    with open(paf_out, "wt") as fo:
        subprocess.run(cmd, stdout=fo, check=True)


def paf_cov_evenness(paf_path, unitig_lens, bin_size):
    cov = defaultdict(lambda: defaultdict(int))
    nbins = {u: max(1, (L + bin_size - 1) // bin_size) for u, L in unitig_lens.items()}
    with open(paf_path, "rt") as fh:
        for ln in fh:
            if not ln.strip() or ln.startswith("@"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            t = f[5]
            try:
                ts = int(f[7])
                te = int(f[8])
                alen = int(f[11])
                mapq = int(f[12]) if len(f) > 12 and f[12].isdigit() else 0
            except ValueError:
                continue
            if t not in unitig_lens or alen <= 0 or mapq < 1:
                continue
            span_s = max(0, min(ts, te))
            span_e = min(unitig_lens[t], max(ts, te))
            if span_e <= span_s or (span_e - span_s) < 50:
                continue
            b0 = span_s // bin_size
            b1 = (span_e - 1) // bin_size
            for b in range(b0, b1 + 1):
                cov[t][b] += 1
    vals = []
    for t, bins in cov.items():
        for b in range(nbins[t]):
            vals.append(bins.get(b, 0))
    if not vals:
        return 0.0, 0.0, float("inf")
    mean = sum(vals) / len(vals)
    var = max(0.0, sum((v - mean) ** 2 for v in vals) / len(vals))
    return mean, var**0.5, (var**0.5 / mean if mean > 0 else float("inf"))


# ---------------- emit (SUMMARY, TSV, JSON), reasons LAST ----------------
def emit(
    verdict,
    ust,
    gst_final,
    gst_utg,
    circ,
    tangle_final,
    cov_stat,
    reasons,
    unitigs_used_from,
    args,
    meta=None,
):
    def _get(d, k, default="."):
        return d.get(k, default) if isinstance(d, dict) else default

    circ_s = (
        "none"
        if circ is None
        else f"{circ.get('name','?')}@{circ.get('ovlp','?')}/{circ.get('ident','?')}"
    )
    reasons_flat = ";".join(str(r) for r in reasons) if reasons else "."
    reasons_flat = reasons_flat.replace("\n", " ").replace("\r", " ").replace("\t", " ")

    final_nodes = _get(gst_final, "nodes", "0")
    final_branch = _get(gst_final, "nbranch", "0")
    utg_nodes = _get(gst_utg, "nodes", "0")
    utg_branch = _get(gst_utg, "nbranch", "0")
    tangle_s = "NA" if str(final_nodes) == "0" else f"{float(tangle_final):.4f}"

    summary = (
        f"SUMMARY v{VERSION} verdict={verdict} "
        f"n={ust['n']} total={ust['total']} max={ust['max']} n50={ust['N50']} "
        f"tangle_final={tangle_s} circ={circ_s} from={unitigs_used_from} "
        f"final_nodes={final_nodes} final_branch={final_branch} "
        f"utg_nodes={utg_nodes} utg_branch={utg_branch} "
        f"reasons={reasons_flat}"
    )

    print(verdict)  # stdout: single token for the ladder
    safe_print_err(summary)  # stderr: human-readable one-liner

    if not getattr(args, "quiet_diag", False):
        obj = dict(
            version=VERSION,
            ust=ust,
            final_graph=gst_final,
            utg_graph=gst_utg,
            tangle_final=(
                None if str(final_nodes) == "0" else float(f"{float(tangle_final):.4f}")
            ),
            circ=circ,
            cov=cov_stat,
            reasons=reasons,
            unitigs_from=unitigs_used_from,
            verdict=verdict,
        )
        safe_json_err(obj)

    if getattr(args, "diag_file", None):
        try:
            with open(args.diag_file, "wt") as fo:
                fo.write(summary + "\n")
        except Exception as e:
            safe_print_err(f"[judge] warn: cannot write diag_file: {e}")

    if getattr(args, "tsv_file", None):
        m = meta or {}
        k = _get(m, "k")
        smer = _get(m, "smer")
        a_arc = _get(m, "a")
        wxc = _get(m, "weak_x")
        uz = _get(m, "unzip")
        ecflag = _get(m, "ec")
        smer_pct = _get(m, "smer_singletons")
        kmer_pct = _get(m, "kmer_singletons")
        row = [
            _get(m, "name"),
            str(k),
            str(smer),
            str(a_arc),
            str(wxc),
            str(uz),
            str(ecflag),
            str(ust["n"]),
            str(ust["total"]),
            verdict,
            str(smer_pct),
            str(kmer_pct),
            tangle_s,
            circ_s,
            str(final_nodes),
            str(final_branch),
            str(utg_nodes),
            str(utg_branch),
            reasons_flat,  # LAST
        ]
        try:
            with open(args.tsv_file, "a") as fo:
                fo.write("\t".join(row) + "\n")
        except Exception as e:
            safe_print_err(f"[judge] warn: cannot write tsv_file: {e}")


# ---------------- main ----------------
def main():
    ignore_sigpipe()
    ap = argparse.ArgumentParser()
    ap.add_argument("--version", action="store_true", help="print version and exit")

    # GFA-only interface
    ap.add_argument(
        "--gfa-final", help="final unitig graph (.utg.final.gfa[.gz]) for decision"
    )
    ap.add_argument(
        "--gfa-utg", help="raw unitig graph (.utg.gfa[.gz]) for diagnostics"
    )

    # thresholds
    ap.add_argument("--mt-size-est", type=int, default=0)
    ap.add_argument("--max-unitigs", type=int, default=80)
    ap.add_argument("--min-total-bp", type=int, default=40000)
    ap.add_argument("--tangle-max", type=float, default=0.06)
    ap.add_argument("--min-n50", type=int, default=8000)
    ap.add_argument("--min-maxlen", type=int, default=12000)
    ap.add_argument("--circ-min-ovlp", type=int, default=1000)
    ap.add_argument("--circ-min-ident", type=float, default=0.95)

    # coverage options (optional)
    ap.add_argument("--reads")
    ap.add_argument("--sample-n", type=int, default=2000)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--cov-bin", type=int, default=100)
    ap.add_argument("--cov-cv-max", type=float, default=0.60)
    ap.add_argument("--cov-mean-min", type=float, default=3.0)

    # strictness / relaxers (optional)
    ap.add_argument("--require-circular", action="store_true")
    ap.add_argument("--min-contig-len-pass", type=int, default=0)
    ap.add_argument("--min-total-bp-pass", type=int, default=0)

    # outputs
    ap.add_argument("--diag-file", help="write SUMMARY to this file")
    ap.add_argument("--tsv-file", help="append TSV row (reasons LAST)")
    ap.add_argument("--quiet-diag", action="store_true")
    ap.add_argument("--trace", action="store_true")

    args = ap.parse_args()
    if args.version:
        print(f"polap-py-oatk-core-judge-assembly.py {VERSION}")
        return 0

    try:
        reasons = []
        # read both graphs (diagnostic & decision)
        gst_utg = (
            gfa_degrees(args.gfa_utg)
            if args.gfa_utg
            else dict(nodes=0, nbranch=0, tips=0)
        )
        gst_final = (
            gfa_degrees(args.gfa_final)
            if args.gfa_final
            else dict(nodes=0, nbranch=0, tips=0)
        )

        # unitigs: prefer final-gfa extraction; else utg
        unitigs_src = None
        unitigs_fa = None
        tmpdir = None
        if args.gfa_final and file_ok(args.gfa_final):
            tmpdir = tempfile.mkdtemp(prefix="judge_fa_")
            fa1 = os.path.join(tmpdir, "unitigs.final.fa")
            if extract_fasta_from_gfa(args.gfa_final, fa1):
                unitigs_src = "final"
                unitigs_fa = fa1
            else:
                reasons.append(f"final_gfa_has_no_sequences:{args.gfa_final}")
        if unitigs_fa is None and args.gfa_utg and file_ok(args.gfa_utg):
            if tmpdir is None:
                tmpdir = tempfile.mkdtemp(prefix="judge_fa_")
            fa2 = os.path.join(tmpdir, "unitigs.utg.fa")
            if extract_fasta_from_gfa(args.gfa_utg, fa2):
                unitigs_src = "utg"
                unitigs_fa = fa2
            else:
                reasons.append(f"utg_gfa_has_no_sequences:{args.gfa_utg}")

        if unitigs_fa is None:
            reasons.append("no_unitigs_source")
            ust = dict(n=0, total=0, max=0, N50=0, names=[])
            emit(
                "CONTINUE",
                ust,
                gst_final,
                gst_utg,
                None,
                float("inf"),
                None,
                reasons,
                "none",
                args,
                None,
            )
            if tmpdir:
                shutil.rmtree(tmpdir, ignore_errors=True)
            return 0

        # stats
        ust = fa_stats(unitigs_fa)
        circ = best_circular_candidate(
            unitigs_fa, args.circ_min_ovlp, args.circ_min_ident
        )

        # dynamic thresholds from mt-size-est
        min_total = args.min_total_bp
        min_maxlen = args.min_maxlen
        min_n50 = args.min_n50
        if args.mt_size_est and args.mt_size_est > 0:
            est = args.mt_size_est
            min_total = max(min_total, int(0.5 * est))
            min_maxlen = max(min_maxlen, int(0.35 * est))
            min_n50 = max(min_n50, int(0.20 * est))

        # choose which graph to judge on (final preferred)
        use_final = args.gfa_final and file_ok(args.gfa_final)
        gst_used = gst_final if use_final else gst_utg
        nodes = max(gst_used.get("nodes", 0), 1)
        tangle_final = (
            (gst_used.get("nbranch", 0) / nodes)
            if gst_used.get("nodes", 0) > 0
            else float("inf")
        )
        if not use_final:
            reasons.append("no_final_gfa_using_utg")

        # structural gating
        if ust["n"] == 0:
            reasons.append("no_unitigs")
        if ust["n"] > args.max_unitigs:
            reasons.append(f"too_many_unitigs({ust['n']}>{args.max_unitigs})")
        if ust["total"] < min_total:
            reasons.append(f"total_bp({ust['total']})<min({min_total})")
        if ust["max"] < min_maxlen:
            reasons.append(f"maxlen({ust['max']})<min({min_maxlen})")
        if ust["N50"] < min_n50:
            reasons.append(f"N50({ust['N50']})<min({min_n50})")
        if gst_used.get("nodes", 0) > 0 and tangle_final > args.tangle_max:
            reasons.append(f"tangle({tangle_final:.3f})>max({args.tangle_max})")

        # optional coverage evenness
        cov_stat = None
        if args.reads and have_tool("minimap2") and ust["n"] > 0:
            try:
                with tempfile.TemporaryDirectory() as td2:
                    sample_fq = os.path.join(td2, "sample.fastq")
                    paf = os.path.join(td2, "map.paf")
                    lens = {}
                    with open_gz(unitigs_fa) as fh:
                        name = None
                        L = 0
                        for ln in fh:
                            if ln.startswith(">"):
                                if name is not None:
                                    lens[name] = L
                                name = ln[1:].strip().split()[0]
                                L = 0
                            else:
                                L += len(ln.strip())
                        if name is not None:
                            lens[name] = L
                    sample_reads(args.reads, sample_fq, args.sample_n)
                    run_minimap2(unitigs_fa, sample_fq, args.threads, paf)
                    mean, stdev, cv = paf_cov_evenness(paf, lens, args.cov_bin)
                    cov_stat = dict(
                        mean=round(mean, 3),
                        stdev=round(stdev, 3),
                        cv=(round(cv, 3) if math.isfinite(cv) else float("inf")),
                    )
                    if mean < args.cov_mean_min:
                        reasons.append(f"cov_mean({mean:.2f})<min({args.cov_mean_min})")
                    if math.isfinite(cv) and cv > args.cov_cv_max:
                        reasons.append(f"cov_cv({cv:.3f})>max({args.cov_cv_max})")
            except subprocess.CalledProcessError as e:
                reasons.append(f"minimap2_failed:{e.returncode}")
            except Exception as e:
                reasons.append(f"cov_eval_error:{type(e).__name__}:{e}")

        # verdict
        verdict = "PASS" if not reasons else "CONTINUE"
        if verdict == "PASS" and args.require_circular and circ is None:
            verdict = "CONTINUE"
            reasons.append("no_circular_candidate")
        if (
            verdict == "CONTINUE"
            and args.min_contig_len_pass > 0
            and ust["max"] >= args.min_contig_len_pass
        ):
            verdict = "PASS"
        if (
            verdict == "CONTINUE"
            and args.min_total_bp_pass > 0
            and ust["total"] >= args.min_total_bp_pass
        ):
            verdict = "PASS"

        # minimal meta (optional; rows still print if missing)
        meta = dict(
            name=os.path.basename(
                os.path.dirname(args.gfa_final or args.gfa_utg or ".")
            ),
            k=".",
            smer=".",
            a=".",
            weak_x=".",
            unzip=".",
            ec=".",
            smer_singletons=".",
            kmer_singletons=".",
        )

        emit(
            verdict,
            ust,
            gst_final,
            gst_utg,
            circ,
            tangle_final,
            cov_stat,
            reasons,
            (unitigs_src or "none"),
            args,
            meta,
        )

        if tmpdir:
            shutil.rmtree(tmpdir, ignore_errors=True)
        return 0

    except Exception as e:
        if args.trace:
            traceback.print_exc(file=sys.stderr)
        print("CONTINUE")
        safe_print_err(
            f"SUMMARY v{VERSION} verdict=CONTINUE n=0 total=0 max=0 n50=0 tangle_final=NA circ=none from=none "
            f"final_nodes=0 final_branch=0 utg_nodes=0 utg_branch=0 reasons=judge_exception:{type(e).__name__}:{e}"
        )
        return 0


if __name__ == "__main__":
    sys.exit(main())
