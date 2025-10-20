#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
polap-py-simulate-wgs-reads.py

Simulate WGS reads from a mixed plant sample for:
 - PacBio HiFi         (--platform hifi)
 - ONT simplex         (--platform ont)
 - ONT high-quality    (--platform onthq)  # R10.4.1 + Kit14 SUP/duplex-like
 - Illumina paired-end (--platform illumina)

Compartments:
 - Nuclear genome (random unless --nuc-fasta provided)
 - Mitochondrial genome (use --mt-fasta if provided; else random of --mt-size)
 - Plastid genome (use --pt-fasta if provided; else random of --pt-size)

Outputs:
  Long-read platforms (hifi/ont/onthq):
    <prefix>.fastq.gz
    <prefix>.nuclear.fastq.gz
    <prefix>.mito.fastq.gz
    <prefix>.plastid.fastq.gz

  Illumina (paired-end, FR):
    <prefix>.R1.fastq.gz, <prefix>.R2.fastq.gz (combined)
    <prefix>.nuclear.R1.fastq.gz,   <prefix>.nuclear.R2.fastq.gz
    <prefix>.mito.R1.fastq.gz,      <prefix>.mito.R2.fastq.gz
    <prefix>.plastid.R1.fastq.gz,   <prefix>.plastid.R2.fastq.gz

Optional truth (all platforms):
  --truth-tsv FILE  → per-read (per-end for Illumina) TSV with:
    read_id, platform, compartment, circular, ref_len, start, frag_len,
    pre_len, post_len, n_sub, n_ins, n_del, n_hp_events

Platform presets:
  --platform hifi   -> Gaussian length model, low indel/sub rates, Q≈40
  --platform ont    -> Log-normal lengths, indel-biased + homopolymer stutter, Q≈20
  --platform onthq  -> Log-normal lengths, strongly reduced errors/stutter (R10.4.1 + Kit14 SUP/duplex-like), Q≈30
  --platform illumina -> Fixed read length (R1/R2), FR pairs, sub-dominated errors with mild 3'-cycle ramp, Q≈30

Circularity:
  By default mt/pt are treated as circular and nuclear as linear. You can override with:
    --mt-circular/--no-mt-circular
    --pt-circular/--no-pt-circular
    --nuc-circular/--no-nuc-circular

References (qualitative behavior and presets):
  Illumina: Schirmer et al., ISME J. 2016; Minoche et al., Genome Biol. 2011; Nakamura et al., NAR 2011.
  PacBio HiFi: Wenger et al., Nat. Biotechnol. 2019.
  ONT indel/homopolymer patterns: Jain et al., Nat. Biotechnol. 2016/2018; De Coster et al., Bioinformatics 2018.
"""

import sys, gzip, random, argparse, math, csv
from typing import List, Tuple, Callable, Optional, Dict

DNA = "ACGT"
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

# ----------------------------- I/O helpers -----------------------------


def read_fasta(path: str) -> str:
    seqs, cur = [], []
    with open(path, "r", encoding="utf-8") as fh:
        for ln in fh:
            ln = ln.rstrip()
            if not ln:
                continue
            if ln.startswith(">"):
                if cur:
                    seqs.append("".join(cur).upper())
                    cur = []
            else:
                cur.append(ln)
        if cur:
            seqs.append("".join(cur).upper())
    if not seqs:
        raise ValueError(f"No FASTA sequence found in {path}")
    return "".join(seqs) if len(seqs) > 1 else seqs[0]


def random_dna(n: int, seed: Optional[int] = None) -> str:
    rng = random.Random(seed) if seed is not None else random
    return "".join(rng.choice(DNA) for _ in range(n))


def write_fastq_gz(
    path: str, entries: List[Tuple[str, str, Optional[str]]], default_qchar: str
):
    if not entries:
        return
    with gzip.open(path, "wt") as gz:
        for name, seq, q in entries:
            if q is None:
                q = default_qchar * len(seq)
            gz.write(f"@{name}\n{seq}\n+\n{q}\n")


def parse_bp_with_suffix(x: str) -> int:
    s = x.strip().lower()
    try:
        if s.endswith("k"):
            return int(float(s[:-1]) * 1_000)
        if s.endswith("m"):
            return int(float(s[:-1]) * 1_000_000)
        if s.endswith("g"):
            return int(float(s[:-1]) * 1_000_000_000)
        return int(float(s))
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid bp value: {x}")


def rc(seq: str) -> str:
    return seq.translate(COMP)[::-1]


# ----------------------------- Length models -----------------------------


def draw_length_gauss(
    rng: random.Random, mean: int, sd: int, min_len: int, max_len: int
) -> int:
    while True:
        L = int(rng.gauss(mean, sd))
        if min_len <= L <= max_len:
            return L


def _lognorm_mu_sigma_from_mean_sd(mean: float, sd: float) -> Tuple[float, float]:
    sigma2 = math.log(1.0 + (sd * sd) / (mean * mean))
    sigma = math.sqrt(sigma2)
    mu = math.log(mean) - 0.5 * sigma2
    return mu, sigma


def draw_length_lognorm(
    rng: random.Random, mean: int, sd: int, min_len: int, max_len: int
) -> int:
    mu, sigma = _lognorm_mu_sigma_from_mean_sd(float(mean), float(sd))
    while True:
        L = int(rng.lognormvariate(mu, sigma))
        if min_len <= L <= max_len:
            return L


# ----------------------------- Error models (HiFi/ONT) -----------------------------
# All mutators return (mut_seq, counts_dict)


def mutate_hifi(
    seq: str, sub_rate: float, ins_rate: float, del_rate: float, rng: random.Random
) -> Tuple[str, Dict[str, int]]:
    counts = {"n_sub": 0, "n_ins": 0, "n_del": 0, "n_hp": 0}
    out, i, L = [], 0, len(seq)
    while i < L:
        r = rng.random()
        if r < del_rate:
            counts["n_del"] += 1
            i += 1
            continue
        if r < del_rate + ins_rate:
            out.append(rng.choice(DNA))
            counts["n_ins"] += 1
            if rng.random() < 0.3:
                out.append(rng.choice(DNA))
                counts["n_ins"] += 1
        if r < del_rate + ins_rate + sub_rate:
            b = seq[i]
            out.append(rng.choice([x for x in DNA if x != b]))
            counts["n_sub"] += 1
        else:
            out.append(seq[i])
        i += 1
    return "".join(out), counts


def mutate_ont(
    seq: str,
    sub_rate: float,
    ins_rate: float,
    del_rate: float,
    hp_rate: float,
    rng: random.Random,
) -> Tuple[str, Dict[str, int]]:
    counts = {"n_sub": 0, "n_ins": 0, "n_del": 0, "n_hp": 0}
    out, i, L = [], 0, len(seq)
    while i < L:
        b = seq[i]
        j = i + 1
        while j < L and seq[j] == b:
            j += 1
        run_len = j - i
        stutter_prob = 0.0
        if run_len >= 3:
            stutter_prob = min(1.0, hp_rate * (run_len - 2) / 3.0)
        if rng.random() < stutter_prob:
            counts["n_hp"] += 1
            if rng.random() < 0.7:
                d = 1 + (rng.random() < 0.5) + (rng.random() < 0.3)  # 1..3
                k = max(1, run_len - min(d, run_len - 1))
            else:
                d = 1 + (rng.random() < 0.5) + (rng.random() < 0.3)
                k = run_len + d
            out.extend(b * k)
        else:
            for _ in range(run_len):
                r = rng.random()
                if r < del_rate:
                    counts["n_del"] += 1
                    continue
                if r < del_rate + ins_rate:
                    out.append(rng.choice(DNA))
                    counts["n_ins"] += 1
                    if rng.random() < 0.5:
                        out.append(rng.choice(DNA))
                        counts["n_ins"] += 1
                if r < del_rate + ins_rate + sub_rate:
                    out.append(rng.choice([x for x in DNA if x != b]))
                    counts["n_sub"] += 1
                else:
                    out.append(b)
        i = j
    return "".join(out), counts


# ----------------------------- Illumina paired-end model -----------------------------
# Returns ((name, seq, qual), counts_dict)


def illumina_mutate_read(
    template: str,
    sub0: float,
    sub_ramp: float,
    ins_rate: float,
    del_rate: float,
    rng: random.Random,
    qconst_override: Optional[int],
) -> Tuple[Tuple[str, str, Optional[str]], Dict[str, int]]:
    counts = {"n_sub": 0, "n_ins": 0, "n_del": 0, "n_hp": 0}
    L = len(template)
    out = []
    quals = []
    denom = max(1, L - 1)
    for c, base in enumerate(template):
        p_sub = min(0.25, max(0.0, sub0 + sub_ramp * (c / denom)))
        r = rng.random()
        if r < del_rate:
            counts["n_del"] += 1
            continue
        if r < del_rate + ins_rate:
            out.append(rng.choice(DNA))
            counts["n_ins"] += 1
        if r < del_rate + ins_rate + p_sub:
            out.append(rng.choice([x for x in DNA if x != base]))
            counts["n_sub"] += 1
            q = 2 if qconst_override is None else qconst_override
        else:
            out.append(base)
            if qconst_override is None:
                q = int(max(2, min(45, round(-10.0 * math.log10(max(1e-6, p_sub))))))
            else:
                q = qconst_override
        quals.append(chr(33 + q))
    seq = "".join(out)
    qual = "".join(quals[: len(seq)])
    if len(qual) < len(seq) and qual:
        qual = qual + qual[-1] * (len(seq) - len(qual))
    elif not qual:
        qual = chr(33 + (qconst_override if qconst_override is not None else 30)) * len(
            seq
        )
    return ("", seq, qual), counts


def illumina_pair_from_fragment(
    frag_seq: str,
    read_len: int,
    sub0: float,
    sub_ramp: float,
    ins_rate: float,
    del_rate: float,
    rng: random.Random,
    qconst_override: Optional[int],
) -> Tuple[
    Tuple[str, str, Optional[str]],
    Tuple[str, str, Optional[str]],
    Dict[str, int],
    Dict[str, int],
]:
    # R1 = leftmost read_len; R2 = rightmost read_len on RC
    if len(frag_seq) < 2 * read_len:
        r1_t = frag_seq[: min(read_len, len(frag_seq))]
        r2_t = rc(frag_seq[-min(read_len, len(frag_seq)) :])
    else:
        r1_t = frag_seq[:read_len]
        r2_t = rc(frag_seq[-read_len:])
    r1, c1 = illumina_mutate_read(
        r1_t, sub0, sub_ramp, ins_rate, del_rate, rng, qconst_override
    )
    r2, c2 = illumina_mutate_read(
        r2_t, sub0, sub_ramp, ins_rate, del_rate, rng, qconst_override
    )
    return r1, r2, c1, c2


# ----------------------------- Core helpers -----------------------------


def slice_chunk(ref: str, st: int, rl: int, circular: bool) -> str:
    n = len(ref)
    if rl <= 0 or n == 0:
        return ""
    if circular:
        st %= n
        if st + rl <= n:
            return ref[st : st + rl]
        else:
            k = (st + rl) - n
            return ref[st:] + ref[:k]
    else:
        st = max(0, min(st, max(0, n - 1)))
        en = min(n, st + rl)
        return ref[st:en]


# ----------------------------- Simulation (with truth logging) -----------------------------


def simulate_long_reads_from_ref(
    ref: str,
    total_bases_target: int,
    draw_length: Callable[[random.Random, int, int, int, int], int],
    draw_len_args: Tuple[int, int, int, int],
    mutate_fn: Callable[[str, random.Random], Tuple[str, Dict[str, int]]],
    rng: random.Random,
    name_prefix: str,
    circular: bool,
    platform: str,
    truth_rows: Optional[List[Dict[str, object]]],
    ref_len_for_truth: int,
) -> List[Tuple[str, str, Optional[str]]]:
    Lref = len(ref)
    if Lref == 0 or total_bases_target <= 0:
        return []
    out: List[Tuple[str, str, Optional[str]]] = []
    acc = 0
    n = 0
    mean, sd, min_len, max_len = draw_len_args
    while acc < total_bases_target:
        rl = draw_length(rng, mean, sd, min_len, max_len)
        st = rng.randrange(Lref)
        chunk = slice_chunk(ref, st, rl, circular=circular)
        pre_len = len(chunk)
        mut, counts = mutate_fn(chunk, rng)
        post_len = len(mut)
        name = f"{name_prefix}_{n:06d}"
        out.append((name, mut, None))
        acc += post_len
        if truth_rows is not None:
            truth_rows.append(
                {
                    "read_id": name,
                    "platform": platform,
                    "compartment": name_prefix.split("_")[0],
                    "circular": int(circular),
                    "ref_len": ref_len_for_truth,
                    "start": st,
                    "frag_len": rl,
                    "pre_len": pre_len,
                    "post_len": post_len,
                    "n_sub": counts["n_sub"],
                    "n_ins": counts["n_ins"],
                    "n_del": counts["n_del"],
                    "n_hp_events": counts["n_hp"],
                }
            )
        n += 1
    return out


def simulate_illumina_pairs_from_ref(
    ref: str,
    total_bases_target: int,
    read_len: int,
    frag_mean: int,
    frag_sd: int,
    rng: random.Random,
    name_prefix: str,
    sub0: float,
    sub_ramp: float,
    ins_rate: float,
    del_rate: float,
    qconst: Optional[int],
    circular: bool,
    truth_rows: Optional[List[Dict[str, object]]],
    ref_len_for_truth: int,
    platform: str,
) -> Tuple[List[Tuple[str, str, Optional[str]]], List[Tuple[str, str, Optional[str]]]]:
    Lref = len(ref)
    if Lref == 0 or total_bases_target <= 0:
        return [], []
    r1_out: List[Tuple[str, str, Optional[str]]] = []
    r2_out: List[Tuple[str, str, Optional[str]]] = []
    acc = 0
    n = 0
    min_frag = max(2 * read_len, 50)

    def draw_frag() -> int:
        for _ in range(1000):
            L = int(rng.gauss(frag_mean, frag_sd))
            if L >= min_frag:
                return L
        return min_frag

    while acc < total_bases_target:
        frag_len = draw_frag()
        if not circular and frag_len > Lref:
            frag_len = Lref
        st = rng.randrange(Lref)
        frag = slice_chunk(ref, st, frag_len, circular=circular)
        r1, r2, c1, c2 = illumina_pair_from_fragment(
            frag,
            read_len,
            sub0,
            sub_ramp,
            ins_rate,
            del_rate,
            rng,
            qconst_override=qconst,
        )
        name = f"{name_prefix}_{n:06d}"
        r1_out.append((f"{name}/1", r1[1], r1[2]))
        r2_out.append((f"{name}/2", r2[1], r2[2]))
        acc += len(r1[1]) + len(r2[1])
        if truth_rows is not None:
            # one truth row per end
            truth_rows.append(
                {
                    "read_id": f"{name}/1",
                    "platform": platform,
                    "compartment": name_prefix.split("_")[0],
                    "circular": int(circular),
                    "ref_len": ref_len_for_truth,
                    "start": st,
                    "frag_len": frag_len,
                    "pre_len": min(
                        read_len, frag_len
                    ),  # template length before per-base events
                    "post_len": len(r1[1]),
                    "n_sub": c1["n_sub"],
                    "n_ins": c1["n_ins"],
                    "n_del": c1["n_del"],
                    "n_hp_events": 0,
                }
            )
            truth_rows.append(
                {
                    "read_id": f"{name}/2",
                    "platform": platform,
                    "compartment": name_prefix.split("_")[0],
                    "circular": int(circular),
                    "ref_len": ref_len_for_truth,
                    "start": (
                        (st + max(0, frag_len - read_len)) % Lref
                        if circular
                        else max(0, min(Lref - 1, st + max(0, frag_len - read_len)))
                    ),
                    "frag_len": frag_len,
                    "pre_len": min(read_len, frag_len),
                    "post_len": len(r2[1]),
                    "n_sub": c2["n_sub"],
                    "n_ins": c2["n_ins"],
                    "n_del": c2["n_del"],
                    "n_hp_events": 0,
                }
            )
        n += 1
    return r1_out, r2_out


# ----------------------------- CLI -----------------------------


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description="Simulate WGS reads (HiFi, ONT, ONT-HQ, or Illumina PE) for nuclear + mito + plastid.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Platform preset
    p.add_argument(
        "--platform",
        choices=["hifi", "ont", "onthq", "illumina"],
        default="hifi",
        help="Sequencing platform to emulate.",
    )

    # References (either provide FASTA or sizes)
    p.add_argument(
        "--nuc-fasta", type=str, default=None, help="Nuclear reference FASTA (optional)"
    )
    p.add_argument(
        "--mt-fasta", type=str, default=None, help="Mitochondrial FASTA (optional)"
    )
    p.add_argument(
        "--pt-fasta", type=str, default=None, help="Plastid FASTA (optional)"
    )

    p.add_argument(
        "--nuc-size",
        type=parse_bp_with_suffix,
        default=parse_bp_with_suffix("10m"),
        help="Nuclear genome size if no --nuc-fasta",
    )
    p.add_argument(
        "--mt-size",
        type=parse_bp_with_suffix,
        default=parse_bp_with_suffix("500k"),
        help="Mito genome size if no --mt-fasta",
    )
    p.add_argument(
        "--pt-size",
        type=parse_bp_with_suffix,
        default=parse_bp_with_suffix("150k"),
        help="Plastid genome size if no --pt-fasta",
    )

    # Depths (fold coverage)
    p.add_argument("--nuc-depth", type=float, default=5.0, help="Nuclear coverage")
    p.add_argument("--mt-depth", type=float, default=50.0, help="Mito coverage")
    p.add_argument("--pt-depth", type=float, default=250.0, help="Plastid coverage")

    # Explicit circular controls
    p.add_argument(
        "--nuc-circular",
        dest="nuc_circ",
        action="store_true",
        help="Treat nuclear as circular",
    )
    p.add_argument("--no-nuc-circular", dest="nuc_circ", action="store_false")
    p.set_defaults(nuc_circ=False)

    p.add_argument(
        "--mt-circular",
        dest="mt_circ",
        action="store_true",
        help="Treat mt as circular",
    )
    p.add_argument("--no-mt-circular", dest="mt_circ", action="store_false")
    p.set_defaults(mt_circ=True)

    p.add_argument(
        "--pt-circular",
        dest="pt_circ",
        action="store_true",
        help="Treat pt as circular",
    )
    p.add_argument("--no-pt-circular", dest="pt_circ", action="store_false")
    p.set_defaults(pt_circ=True)

    # Length model (for long-read platforms; Illumina ignores)
    p.add_argument(
        "--len-model",
        choices=["gauss", "lognorm"],
        default=None,
        help="Length model; default based on --platform.",
    )
    p.add_argument("--mean", type=int, default=None, help="Mean read length (bp)")
    p.add_argument("--sd", type=int, default=None, help="Stddev read length (bp)")
    p.add_argument("--min", type=int, default=None, help="Min read length (bp)")
    p.add_argument("--max", type=int, default=None, help="Max read length (bp)")

    # Error model (explicit values override platform defaults)
    p.add_argument("--sub", type=float, default=None, help="Substitution rate per base")
    p.add_argument("--ins", type=float, default=None, help="Insertion rate per base")
    p.add_argument("--dele", type=float, default=None, help="Deletion rate per base")
    p.add_argument(
        "--hp-rate",
        type=float,
        default=None,
        help="Homopolymer stutter multiplier (ONT only; ignored for HiFi/Illumina)",
    )

    # Illumina-specific options
    p.add_argument(
        "--read-len", type=int, default=150, help="Illumina read length per end (bp)"
    )
    p.add_argument(
        "--frag-mean", type=int, default=350, help="Illumina fragment mean length (bp)"
    )
    p.add_argument(
        "--frag-sd", type=int, default=50, help="Illumina fragment stddev (bp)"
    )
    p.add_argument(
        "--cycle-ramp",
        type=float,
        default=0.0015,
        help="Per-cycle increase in substitution probability from 5'→ ' (Illumina)",
    )
    p.add_argument(
        "--pe-orientation",
        choices=["FR"],
        default="FR",
        help="Paired-end orientation (currently FR only).",
    )

    # Quality (constant Phred per base; Illumina may derive per-cycle if None)
    p.add_argument(
        "--qconst",
        type=int,
        default=None,
        help="Constant Phred quality to emit (0–45). For Illumina, if omitted, derive from per-cycle error.",
    )

    # Truth sidecar
    p.add_argument(
        "--truth-tsv",
        type=str,
        default=None,
        help="Optional path to write per-read truth table (TSV).",
    )

    # Output
    p.add_argument(
        "-o", "--out-prefix", type=str, default="wgs.sim", help="Output prefix."
    )
    p.add_argument("--seed", type=int, default=43, help="Random seed")

    args = p.parse_args(argv)

    rng = random.Random(args.seed)

    # ------------------ Platform presets & defaults ------------------
    if args.platform == "hifi":
        len_model = args.len_model if args.len_model is not None else "gauss"
        mean = args.mean if args.mean is not None else 10_000
        sd = args.sd if args.sd is not None else 1_500
        min_len = args.min if args.min is not None else 4_000
        max_len = args.max if args.max is not None else 20_000
        sub_rate = args.sub if args.sub is not None else 0.0020
        ins_rate = args.ins if args.ins is not None else 0.0003
        del_rate = args.dele if args.dele is not None else 0.0003
        hp_rate = 0.0
        qconst = args.qconst if args.qconst is not None else 40
        mutate_wrap = lambda s, r: mutate_hifi(s, sub_rate, ins_rate, del_rate, r)
        draw_len_fn = draw_length_gauss if len_model == "gauss" else draw_length_lognorm

    elif args.platform == "ont":
        len_model = args.len_model if args.len_model is not None else "lognorm"
        mean = args.mean if args.mean is not None else 15_000
        sd = args.sd if args.sd is not None else 5_000
        min_len = args.min if args.min is not None else 1_000
        max_len = args.max if args.max is not None else 150_000
        sub_rate = args.sub if args.sub is not None else 0.004
        ins_rate = args.ins if args.ins is not None else 0.008
        del_rate = args.dele if args.dele is not None else 0.008
        hp_rate = args.hp_rate if args.hp_rate is not None else 0.20
        qconst = args.qconst if args.qconst is not None else 20
        mutate_wrap = lambda s, r: mutate_ont(
            s, sub_rate, ins_rate, del_rate, hp_rate, r
        )
        draw_len_fn = (
            draw_length_lognorm if len_model == "lognorm" else draw_length_gauss
        )

    elif args.platform == "onthq":
        len_model = args.len_model if args.len_model is not None else "lognorm"
        mean = args.mean if args.mean is not None else 18_000
        sd = args.sd if args.sd is not None else 6_000
        min_len = args.min if args.min is not None else 1_000
        max_len = args.max if args.max is not None else 200_000
        sub_rate = args.sub if args.sub is not None else 0.0010
        ins_rate = args.ins if args.ins is not None else 0.0015
        del_rate = args.dele if args.dele is not None else 0.0015
        hp_rate = args.hp_rate if args.hp_rate is not None else 0.05
        qconst = args.qconst if args.qconst is not None else 30
        mutate_wrap = lambda s, r: mutate_ont(
            s, sub_rate, ins_rate, del_rate, hp_rate, r
        )
        draw_len_fn = (
            draw_length_lognorm if len_model == "lognorm" else draw_length_gauss
        )

    else:  # illumina
        read_len = int(args.read_len)
        frag_mean = int(args.frag_mean)
        frag_sd = int(args.frag_sd)
        sub0 = args.sub if args.sub is not None else 0.001
        sub_ramp = args.cycle_ramp
        ins_rate = args.ins if args.ins is not None else 1e-5
        del_rate = args.dele if args.dele is not None else 1e-5
        qconst = args.qconst
        if read_len <= 0:
            raise SystemExit("Illumina --read-len must be > 0")
        if frag_mean <= 0 or frag_sd < 0:
            raise SystemExit("Illumina --frag-mean/--frag-sd invalid")
        if args.pe_orientation != "FR":
            raise SystemExit("Only FR paired-end orientation is supported currently.")

    # Bounds & sanity (for long-read paths)
    if args.platform in ("hifi", "ont", "onthq"):
        if qconst < 0:
            qconst = 0
        if qconst > 45:
            qconst = 45
        default_qchar = chr(33 + qconst)
        if min_len <= 0 or max_len <= 0 or min_len > max_len:
            raise SystemExit("Invalid min/max read length.")
        if mean <= 0 or sd <= 0:
            raise SystemExit("Invalid mean/sd read length.")
        if not (0 <= sub_rate < 1 and 0 <= ins_rate < 1 and 0 <= del_rate < 1):
            raise SystemExit("Error rates must be in [0,1).")

    # Load or synthesize references
    if args.nuc_fasta:
        nuc_ref = read_fasta(args.nuc_fasta)
        nuc_size = len(nuc_ref)
    else:
        nuc_size = int(args.nuc_size)
        nuc_ref = random_dna(nuc_size, seed=args.seed + 11)

    if args.mt_fasta:
        mt_ref = read_fasta(args.mt_fasta)
        mt_size = len(mt_ref)
    else:
        mt_size = int(args.mt_size)
        mt_ref = random_dna(mt_size, seed=args.seed + 22)

    if args.pt_fasta:
        pt_ref = read_fasta(args.pt_fasta)
        pt_size = len(pt_ref)
    else:
        pt_size = int(args.pt_size)
        pt_ref = random_dna(pt_size, seed=args.seed + 33)

    # Target bases by coverage
    nuc_target = int(round(nuc_size * args.nuc_depth))
    mt_target = int(round(mt_size * args.mt_depth))
    pt_target = int(round(pt_size * args.pt_depth))

    sys.stderr.write(
        f"[sim] platform={args.platform}"
        + (
            f" len_model={args.len_model or ('gauss' if args.platform=='hifi' else 'lognorm')}"
            if args.platform != "illumina"
            else f" read_len={args.read_len} frag_mean={args.frag_mean} frag_sd={args.frag_sd}"
        )
        + f" qconst={args.qconst if args.qconst is not None else 'derived'}\n"
        f"[sim] circularity: nuc={args.nuc_circ} mt={args.mt_circ} pt={args.pt_circ}\n"
        f"[sim] nuclear: {nuc_size} bp @ {args.nuc_depth:.2f}× → ~{nuc_target} bp\n"
        f"[sim] mito   : {mt_size}  bp @ {args.mt_depth:.2f}× → ~{mt_target} bp\n"
        f"[sim] plastid: {pt_size}  bp @ {args.pt_depth:.2f}× → ~{pt_target} bp\n"
    )

    # Initialize truth collector
    truth_rows: Optional[List[Dict[str, object]]] = [] if args.truth_tsv else None

    out = args.out_prefix

    if args.platform in ("hifi", "ont", "onthq"):
        draw_args = (mean, sd, min_len, max_len)

        nuc_reads = simulate_long_reads_from_ref(
            nuc_ref,
            nuc_target,
            draw_len_fn,
            draw_args,
            mutate_wrap,
            rng,
            "nuc",
            circular=args.nuc_circ,
            platform=args.platform,
            truth_rows=truth_rows,
            ref_len_for_truth=nuc_size,
        )
        mt_reads = simulate_long_reads_from_ref(
            mt_ref,
            mt_target,
            draw_len_fn,
            draw_args,
            mutate_wrap,
            rng,
            "mt",
            circular=args.mt_circ,
            platform=args.platform,
            truth_rows=truth_rows,
            ref_len_for_truth=mt_size,
        )
        pt_reads = simulate_long_reads_from_ref(
            pt_ref,
            pt_target,
            draw_len_fn,
            draw_args,
            mutate_wrap,
            rng,
            "pt",
            circular=args.pt_circ,
            platform=args.platform,
            truth_rows=truth_rows,
            ref_len_for_truth=pt_size,
        )

        all_reads = nuc_reads + mt_reads + pt_reads
        rng.shuffle(all_reads)

        default_qchar = chr(33 + (args.qconst if args.qconst is not None else 30))
        write_fastq_gz(f"{out}.fastq.gz", all_reads, default_qchar)
        write_fastq_gz(f"{out}.nuclear.fastq.gz", nuc_reads, default_qchar)
        write_fastq_gz(f"{out}.mito.fastq.gz", mt_reads, default_qchar)
        write_fastq_gz(f"{out}.plastid.fastq.gz", pt_reads, default_qchar)

        total_bp = sum(len(seq) for _, seq, _ in all_reads)
        sys.stderr.write(
            f"[sim] wrote combined ~{total_bp} bp in {len(all_reads)} reads\n"
        )
        sys.stderr.write(
            f"[sim] outputs:\n"
            f"  {out}.fastq.gz (combined)\n"
            f"  {out}.nuclear.fastq.gz\n"
            f"  {out}.mito.fastq.gz\n"
            f"  {out}.plastid.fastq.gz\n"
        )

    else:  # Illumina paired-end
        nuc_r1, nuc_r2 = simulate_illumina_pairs_from_ref(
            nuc_ref,
            nuc_target,
            read_len,
            frag_mean,
            frag_sd,
            rng,
            "nuc",
            sub0,
            args.cycle_ramp,
            ins_rate,
            del_rate,
            args.qconst,
            circular=args.nuc_circ,
            truth_rows=truth_rows,
            ref_len_for_truth=nuc_size,
            platform=args.platform,
        )
        mt_r1, mt_r2 = simulate_illumina_pairs_from_ref(
            mt_ref,
            mt_target,
            read_len,
            frag_mean,
            frag_sd,
            rng,
            "mt",
            sub0,
            args.cycle_ramp,
            ins_rate,
            del_rate,
            args.qconst,
            circular=args.mt_circ,
            truth_rows=truth_rows,
            ref_len_for_truth=mt_size,
            platform=args.platform,
        )
        pt_r1, pt_r2 = simulate_illumina_pairs_from_ref(
            pt_ref,
            pt_target,
            read_len,
            frag_mean,
            frag_sd,
            rng,
            "pt",
            sub0,
            args.cycle_ramp,
            ins_rate,
            del_rate,
            args.qconst,
            circular=args.pt_circ,
            truth_rows=truth_rows,
            ref_len_for_truth=pt_size,
            platform=args.platform,
        )

        all_r1 = nuc_r1 + mt_r1 + pt_r1
        all_r2 = nuc_r2 + mt_r2 + pt_r2
        idx = list(range(len(all_r1)))
        rng.shuffle(idx)
        all_r1 = [all_r1[i] for i in idx]
        all_r2 = [all_r2[i] for i in idx]

        default_qchar = chr(33 + (args.qconst if args.qconst is not None else 30))
        write_fastq_gz(f"{out}.R1.fastq.gz", all_r1, default_qchar)
        write_fastq_gz(f"{out}.R2.fastq.gz", all_r2, default_qchar)
        write_fastq_gz(f"{out}.nuclear.R1.fastq.gz", nuc_r1, default_qchar)
        write_fastq_gz(f"{out}.nuclear.R2.fastq.gz", nuc_r2, default_qchar)
        write_fastq_gz(f"{out}.mito.R1.fastq.gz", mt_r1, default_qchar)
        write_fastq_gz(f"{out}.mito.R2.fastq.gz", mt_r2, default_qchar)
        write_fastq_gz(f"{out}.plastid.R1.fastq.gz", pt_r1, default_qchar)
        write_fastq_gz(f"{out}.plastid.R2.fastq.gz", pt_r2, default_qchar)

        total_bp = sum(len(seq) for _, seq, _ in all_r1) + sum(
            len(seq) for _, seq, _ in all_r2
        )
        sys.stderr.write(
            f"[sim] wrote combined ~{total_bp} bp across {len(all_r1)} pairs\n"
        )
        sys.stderr.write(
            f"[sim] outputs:\n"
            f"  {out}.R1.fastq.gz, {out}.R2.fastq.gz (combined)\n"
            f"  {out}.nuclear.R1.fastq.gz, {out}.nuclear.R2.fastq.gz\n"
            f"  {out}.mito.R1.fastq.gz, {out}.mito.R2.fastq.gz\n"
            f"  {out}.plastid.R1.fastq.gz, {out}.plastid.R2.fastq.gz\n"
        )

    # Write truth if requested
    if args.truth_tsv and truth_rows is not None:
        with open(args.truth_tsv, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(
                [
                    "read_id",
                    "platform",
                    "compartment",
                    "circular",
                    "ref_len",
                    "start",
                    "frag_len",
                    "pre_len",
                    "post_len",
                    "n_sub",
                    "n_ins",
                    "n_del",
                    "n_hp_events",
                ]
            )
            for row in truth_rows:
                w.writerow(
                    [
                        row["read_id"],
                        row["platform"],
                        row["compartment"],
                        row["circular"],
                        row["ref_len"],
                        row["start"],
                        row["frag_len"],
                        row["pre_len"],
                        row["post_len"],
                        row["n_sub"],
                        row["n_ins"],
                        row["n_del"],
                        row["n_hp_events"],
                    ]
                )
        sys.stderr.write(f"[sim] truth: {args.truth_tsv} ({len(truth_rows)} rows)\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
