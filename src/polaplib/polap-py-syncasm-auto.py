#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
polap-py-syncasm-auto.py
# version: v0.3.0 (2025-09-18)

Automated driver for running syncasm with parameter sweeps.

Key features:
  • Correct syncasm CLI: -o <prefix>, -t (fixed 8 per job), --weak-cross, --unzip-round, positional reads
  • Creates output folders (syncasm itself does not)
  • HPC via seqtk hpc (writes <out>/reads.hpc.fastq)
      - Default: DO NOT scale k/s
      - Use --hpc-scale to scale k/s by the HPC factor while roughly preserving (k - s + 1); clamp s<=31
  • Unified -c/--c (one value = override; multiple values = sweep list)
  • Auto-c sweep for mixed data (default): seed draft -> minimap2 -> depth histogram valley -> c-band
  • Parallel execution (GNU parallel preferred; ThreadPool fallback)
  • Skip detection via <run>/syncasm.asm.utg.gfa
  • Per-job logs: job.out / job.err
  • Summarization: count unitigs and total bp from utg.gfa

Examples
--------
# Simple run (no HPC, auto-c, default k/s)
python polap-py-syncasm-auto.py --out run_auto --reads reads.fastq -t 8

# HPC reads but keep your k/s
python polap-py-syncasm-auto.py --out run_auto --reads reads.fastq --hpc --k 121 --s 23 -t 16

# HPC reads AND scale k/s (optional)
python polap-py-syncasm-auto.py --out run_auto --reads reads.fastq --hpc --hpc-scale --k 121 --s 23 -t 16

# Force c or sweep specific c list
python polap-py-syncasm-auto.py --out run_auto --reads reads.fastq --k 121 --s 23 -c 30 -t 8
python polap-py-syncasm-auto.py --out run_auto --reads reads.fastq --k 121 --s 23 -c 26 30 34 --ascending-c -t 8
"""

import argparse, gzip, os, statistics, subprocess, sys, shutil, json
from pathlib import Path
from typing import List, Dict, Tuple, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed

# ─────────────────────────────
# Configs / mappings

TECH_ERROR = {
    # per-base error rates used only for quick q=(1-e)^k survival modeling
    "hifi": 0.01,  # PacBio HiFi/CCS (tune as needed)
    "duplex": 0.01,  # ONT Duplex consensus (ballpark)
    "sup": 0.03,  # ONT R10.4 SUP polished
    "ont": 0.07,  # vanilla ONT (pre-R10.4) / conservative
}

# ─────────────────────────────
# helpers


def okexe(name: str) -> bool:
    return shutil.which(name) is not None


def run_cmd(cmd: List[str], quiet=True) -> int:
    return subprocess.run(
        cmd,
        stdout=(subprocess.DEVNULL if quiet else None),
        stderr=(subprocess.DEVNULL if quiet else None),
    ).returncode


def readlen_sample(path: Path, limit: int = 20000) -> Tuple[int, float]:
    lens = []
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt", encoding="utf-8", errors="ignore") as fh:
        i = 0
        while True:
            h = fh.readline()
            if not h:
                break
            s = fh.readline()
            fh.readline()
            fh.readline()
            lens.append(len(s.rstrip()))
            i += 1
            if i >= limit:
                break
    if not lens:
        return 0, 0.0
    return len(lens), float(statistics.mean(lens))


def count_reads(path: Path) -> int:
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt", encoding="utf-8", errors="ignore") as fh:
        nlines = sum(1 for _ in fh)
    return nlines // 4


def survival_prob_exact_k(e: float, k: int) -> float:
    return (1.0 - e) ** k


def suggest_c_band_from_support(exp_support: float, floor: int = 10) -> List[int]:
    """Return a compact c-list around expected per-edge support."""
    if exp_support <= 0:
        return [max(3, floor)]
    low = max(floor, int(exp_support * 0.6))
    high = max(low + 1, int(exp_support * 1.6))
    span = high - low
    if span <= 2:
        return list(range(low, high + 1))
    step = max(1, span // 3)
    return list(range(low, high + 1, step))[:4]


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


# ─────────────────────────────
# GFA summary


def summarize_utg_gfa(utg: Path) -> Tuple[Optional[int], Optional[int]]:
    """Return (unitigs, total_bp) from GFA S-lines (sequence length or LN tag)."""
    if not utg.exists() or utg.stat().st_size == 0:
        return None, None
    n = 0
    L = 0
    try:
        with utg.open("r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if not line or line[0] != "S":
                    continue
                n += 1
                cols = line.rstrip("\n").split("\t")
                seq = cols[2] if len(cols) >= 3 else "*"
                if seq not in ("", "*"):
                    L += len(seq)
                else:
                    for t in cols[3:]:
                        if t.startswith("LN:i:"):
                            try:
                                L += int(t.split(":")[-1])
                            except Exception:
                                pass
                            break
        return n, L
    except Exception:
        return None, None


# ─────────────────────────────
# Mixed-data auto-c (seed draft -> depth valley)


def auto_c_for_mixed(
    reads: Path, k_max: int, tech: str, outroot: Path, min_floor: int = 10
) -> List[int]:
    """Seed a rough draft -> depth histogram -> find valley -> propose a small c sweep."""
    if not (okexe("minimap2") and okexe("samtools")):
        # Fallback: assume ~100x coverage
        e = TECH_ERROR[tech]
        q = survival_prob_exact_k(e, k_max)
        return suggest_c_band_from_support(100.0 * q, floor=min_floor)

    seed = outroot / "_seed"
    ensure_dir(seed)
    # quick & permissive seed (just to get targets for mapping)
    subprocess.run(
        [
            "syncasm",
            "-k",
            "111",
            "-s",
            "23",
            "-c",
            "10",
            "-a",
            "0.028",
            "--weak-cross",
            "0.010",
            "--unzip-round",
            "0",
            "-t",
            "8",
            "-o",
            str(seed / "syncasm.asm"),
            str(reads),
        ],
        check=False,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )

    # choose a reference for mapping
    ref = seed / "syncasm.asm.ctg.fa"
    if not ref.exists() or ref.stat().st_size == 0:
        # (Optional) could derive fasta from utg.gfa S sequences; if missing, fallback
        ref = seed / "syncasm.asm.utg.gfa"
    if not ref.exists() or ref.stat().st_size == 0:
        e = TECH_ERROR[tech]
        q = survival_prob_exact_k(e, k_max)
        return suggest_c_band_from_support(100.0 * q, floor=min_floor)

    # map & depth
    bam = seed / "seed.bam"
    p1 = subprocess.Popen(
        ["minimap2", "-t", "8", "-ax", "map-ont", str(ref), str(reads)],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=False,
    )
    p2 = subprocess.Popen(
        ["samtools", "sort", "-@", "2", "-o", str(bam)],
        stdin=p1.stdout,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    p1.stdout.close()
    p2.wait()
    if p2.returncode != 0:
        e = TECH_ERROR[tech]
        q = survival_prob_exact_k(e, k_max)
        return suggest_c_band_from_support(100.0 * q, floor=min_floor)
    subprocess.run(
        ["samtools", "index", str(bam)],
        check=False,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )

    # depth -> histogram (clip tail), smooth
    depths = []
    with subprocess.Popen(
        ["samtools", "depth", str(bam)], stdout=subprocess.PIPE, text=True
    ) as P:
        for line in P.stdout:
            try:
                d = int(line.rsplit("\t", 1)[-1])
                depths.append(d)
            except:  # noqa: E722
                pass
    if not depths:
        e = TECH_ERROR[tech]
        q = survival_prob_exact_k(e, k_max)
        return suggest_c_band_from_support(100.0 * q, floor=min_floor)

    depths.sort()
    clip = depths[int(len(depths) * 0.995)]
    bins = [0] * (clip + 1)
    for d in depths:
        if d <= clip:
            bins[d] += 1

    w = 7
    half = w // 2
    sm = [0.0] * len(bins)
    if len(bins) <= w:
        sm = [float(x) for x in bins]
    else:
        s = sum(bins[:w])
        sm[half] = s / float(w)
        for i in range(half + 1, len(bins) - half):
            s += bins[i + half] - bins[i - half - 1]
            sm[i] = s / float(w)
        for i in range(half):
            sm[i] = float(bins[i])
            sm[-1 - i] = float(bins[-1 - i])

    peaks = []
    for i in range(1, len(sm) - 1):
        if sm[i] >= sm[i - 1] and sm[i] >= sm[i + 1]:
            peaks.append((sm[i], i))
    peaks.sort(reverse=True)
    if len(peaks) < 2:
        valley = depths[int(0.75 * len(depths))]
    else:
        i1, i2 = peaks[0][1], peaks[1][1]
        lo, hi = sorted((i1, i2))
        valley = min(range(lo, hi + 1), key=lambda j: sm[j])

    # convert valley coverage to support, then c-band
    e = TECH_ERROR[tech]
    q = survival_prob_exact_k(e, k_max)
    exp_support = valley * q
    return suggest_c_band_from_support(exp_support, floor=min_floor)


# ─────────────────────────────
# main


def main():
    ap = argparse.ArgumentParser(description="Auto -c sweep and parallel syncasm grid.")
    ap.add_argument("--out", required=True, help="Output root directory")
    ap.add_argument("--reads", required=True, help="Reads (FASTQ/FASTA)")
    ap.add_argument(
        "-t",
        "--threads",
        type=int,
        default=os.cpu_count() or 8,
        help="Number of parallel jobs [default: CPU count]",
    )
    # sequencing technology -> error model
    ap.add_argument(
        "--tech",
        choices=list(TECH_ERROR.keys()),
        default="ont",
        help="Sequencing tech for survival model [ont]",
    )
    # k/s lists
    ap.add_argument("--k", type=int, nargs="+", default=[121], help="k values [121]")
    ap.add_argument("--s", type=int, nargs="+", default=[23], help="s values [23]")
    # c control (unified)
    ap.add_argument(
        "--c",
        type=int,
        nargs="+",
        help="One or more c values; one -> override, many -> sweep list",
    )
    ap.add_argument(
        "--ascending-c",
        action="store_true",
        help="Sweep c in ascending order (default descending)",
    )
    # bridge knobs
    ap.add_argument(
        "--a",
        type=float,
        nargs="+",
        default=[0.028, 0.030, 0.032, 0.034, 0.036],
        help="a (arc coverage) values",
    )
    ap.add_argument(
        "--wx",
        type=float,
        nargs="+",
        default=[0.010, 0.015, 0.020],
        help="--weak-cross values",
    )
    ap.add_argument(
        "--u", type=int, nargs="+", default=[0, 2], help="--unzip-round values"
    )
    # HPC
    ap.add_argument(
        "--hpc", action="store_true", help="Use seqtk hpc -> <out>/reads.hpc.fastq"
    )
    ap.add_argument(
        "--hpc-scale",
        action="store_true",
        help="Scale k/s by the HPC factor (default: keep user k/s)",
    )
    # misc
    ap.add_argument("--jobfile", type=Path, help="Write GNU parallel jobfile")
    ap.add_argument("--report-json", type=Path, help="Write JSON report with results")
    ap.add_argument(
        "--verbose", "-v", action="count", default=0, help="Increase verbosity"
    )
    args = ap.parse_args()

    outroot = Path(args.out)
    ensure_dir(outroot)
    reads = Path(args.reads)

    # (1) Optional HPC conversion (plain FASTQ)
    if args.hpc:
        if not okexe("seqtk"):
            print("ERROR: seqtk not found for --hpc", file=sys.stderr)
            sys.exit(2)
        hpc_reads = outroot / "reads.hpc.fastq"
        print(f"# [hpc] compressing {reads} -> {hpc_reads}")
        with open(hpc_reads, "w") as out:
            subprocess.run(["seqtk", "hpc", str(reads)], stdout=out, check=True)
        _, m_raw = readlen_sample(reads, limit=20000)
        _, m_hpc = readlen_sample(hpc_reads, limit=20000)
        if m_raw > 0 and m_hpc > 0:
            f = m_raw / m_hpc
            if args.hpc_scale:
                # preserve span = (k - s + 1); clamp s<=31 (syncasm limit)
                span = args.k[0] - args.s[0] + 1
                k2 = max(31, int(round(args.k[0] / f)))
                s2 = max(5, min(31, k2 - span + 1))
                args.k = [k2]
                args.s = [s2]
            print(f"# [hpc] factor≈{f:.3f}; k={args.k[0]}, s={args.s[0]}")
        reads = hpc_reads

    # (2) Build c_values (unified --c)
    if args.c:
        if len(args.c) == 1:
            c_values = [int(args.c[0])]
            print(f"# [c] override={c_values}")
        else:
            c_values = sorted(set(int(x) for x in args.c))
            print(f"# [c] list(before-order)={c_values}")
    else:
        # default mixed-data estimation
        c_values = auto_c_for_mixed(
            reads, k_max=max(args.k), tech=args.tech, outroot=outroot, min_floor=10
        )
        print(f"# [auto-c] sweep(before-order)={c_values}")

    # ordering: default descending, --ascending-c flips it
    c_values = sorted(c_values) if args.ascending_c else sorted(c_values, reverse=True)
    print(f"# [c-order] {'ascending' if args.ascending_c else 'descending'} {c_values}")

    # (3) Build job grid (ensure each run folder exists)
    jobs: List[Tuple[Path, List[str]]] = []
    for k in args.k:
        for s in args.s:
            for c in c_values:
                for a in args.a:
                    for wx in args.wx:
                        for u in args.u:
                            name = f"exp_a{int(round(a*1000)):03d}_k{k}_s{s}_c{c}_wx{int(round(wx*1000)):03d}_u{u}"
                            od = outroot / name
                            ensure_dir(od)
                            prefix = str(od / "syncasm.asm")
                            cmd = [
                                "syncasm",
                                "-k",
                                str(k),
                                "-s",
                                str(s),
                                "-c",
                                str(c),
                                "-a",
                                f"{a:.3f}",
                                "--weak-cross",
                                f"{wx:.3f}",
                                "--unzip-round",
                                str(u),
                                "-t",
                                "8",  # per-job threads fixed at 8
                                "-o",
                                prefix,
                                str(reads),
                            ]
                            jobs.append((od, cmd))
    print(f"# [grid] jobs={len(jobs)}")

    # (4) Optional jobfile for GNU parallel
    if args.jobfile:
        with args.jobfile.open("w") as fh:
            for od, cmd in jobs:
                utg = str(od / "syncasm.asm.utg.gfa")
                line = (
                    f'OUT="{od}"; UTG="{utg}"; [[ -s "$UTG" ]] && echo "[skip] $OUT" && exit 0; '
                    f'mkdir -p "$OUT"; '
                    + " ".join(cmd)
                    + ' >"$OUT/job.out" 2>"$OUT/job.err"\n'
                )
                fh.write(line)
        print(f"[ok] wrote jobfile: {args.jobfile}")

    # (5) Execute in parallel
    if okexe("parallel"):
        print(f"# [run] GNU parallel -j {args.threads}")
        p = subprocess.Popen(
            [
                "parallel",
                "-j",
                str(args.threads),
                "--eta",
                "--lb",
                "--halt",
                "now,fail=1",
                "bash -lc {}",
            ],
            stdin=subprocess.PIPE,
            text=True,
        )
        for od, cmd in jobs:
            utg = od / "syncasm.asm.utg.gfa"
            line = (
                f'OUT="{od}"; UTG="{utg}"; [[ -s "$UTG" ]] && echo "[skip] $OUT" && exit 0; '
                f'mkdir -p "$OUT"; '
                + " ".join(cmd)
                + ' >"$OUT/job.out" 2>"$OUT/job.err"'
            )
            p.stdin.write(line + "\n")
        p.stdin.close()
        p.wait()
    else:
        print(f"# [run] Python ThreadPool (workers={args.threads})")

        def worker(task: Tuple[Path, List[str]]):
            od, cmd = task
            ensure_dir(od)
            utg = od / "syncasm.asm.utg.gfa"
            if utg.exists() and utg.stat().st_size > 0:
                return "[skip]", str(od)
            rc = run_cmd(cmd, quiet=(args.verbose < 2))
            return ("[ok]" if rc == 0 else "[err]"), str(od)

        with ThreadPoolExecutor(max_workers=args.threads) as ex:
            futs = [ex.submit(worker, t) for t in jobs]
            for fu in as_completed(futs):
                tag, od = fu.result()
                print(tag, od)

    # (6) Summarize results
    rows = []
    for od, _ in jobs:
        utg = od / "syncasm.asm.utg.gfa"
        n, L = summarize_utg_gfa(utg)
        rows.append(
            {
                "name": od.name,
                "unitigs": str(n) if n is not None else ".",
                "total_bp": str(L) if L is not None else ".",
            }
        )
    # Simple display (top 10)
    print("\n# Top results (name, unitigs, total_bp)")
    for r in rows[:10]:
        print(f"{r['name']}\tunitigs={r['unitigs']}\ttotal_bp={r['total_bp']}")

    if args.report_json:
        args.report_json.write_text(json.dumps({"rows": rows}, indent=2))


if __name__ == "__main__":
    main()
