#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
polap-py-make-nuclear-from-busco-aa.py v0.0.2

Build a synthetic nuclear genome by reverse-translating BUSCO amino-acid sequences
into CDS, optionally adding introns/UTRs, and stitching genes with intergenic
spacers until a target genome size is reached.

BUSCO dataset handling:
  - If --busco-faa is not given, the script will attempt to use BUSCO to download
    the lineage (default: viridiplantae_odb12) via: `busco --download <lineage>`.
  - If the `busco` command is missing and you're in a non-base conda env
    (CONDA_DEFAULT_ENV != 'base'), the script attempts: `conda install -y busco`.
  - If `busco` is still unavailable or conda env is not activated at all (no CONDA_DEFAULT_ENV),
    the script exits with an error and instructions.

Outputs (always placed under --out-dir):
  genome.fa
  genome.gff3
  genome.bed
  genome.truth.tsv

References:
  - BUSCO datasets: Simão et al., Bioinformatics 2015; Waterhouse et al., MBE 2018.
  - Reverse translation & codon usage basics: Sharp & Li, NAR 1987.
"""

import sys, os, argparse, random, gzip, csv, shutil, subprocess
from typing import Dict, List, Tuple, Optional

DNA = "ACGT"
COMP = str.maketrans("ACGTacgtnN", "TGCAtgcanN")

# Plant-like codon usage (heuristic). Override via --codon-usage CSV.
GENERIC_PLANT_CODON_USAGE = {
    "A": {"GCT": 0.25, "GCC": 0.28, "GCA": 0.22, "GCG": 0.25},
    "R": {"CGT": 0.20, "CGC": 0.22, "CGA": 0.18, "CGG": 0.18, "AGA": 0.11, "AGG": 0.11},
    "N": {"AAT": 0.48, "AAC": 0.52},
    "D": {"GAT": 0.54, "GAC": 0.46},
    "C": {"TGT": 0.45, "TGC": 0.55},
    "Q": {"CAA": 0.58, "CAG": 0.42},
    "E": {"GAA": 0.59, "GAG": 0.41},
    "G": {"GGT": 0.34, "GGC": 0.26, "GGA": 0.25, "GGG": 0.15},
    "H": {"CAT": 0.55, "CAC": 0.45},
    "I": {"ATT": 0.48, "ATC": 0.39, "ATA": 0.13},
    "L": {"TTA": 0.08, "TTG": 0.13, "CTT": 0.25, "CTC": 0.20, "CTA": 0.13, "CTG": 0.21},
    "K": {"AAA": 0.62, "AAG": 0.38},
    "M": {"ATG": 1.0},
    "F": {"TTT": 0.55, "TTC": 0.45},
    "P": {"CCT": 0.30, "CCC": 0.23, "CCA": 0.29, "CCG": 0.18},
    "S": {"TCT": 0.22, "TCC": 0.22, "TCA": 0.18, "TCG": 0.10, "AGT": 0.14, "AGC": 0.14},
    "T": {"ACT": 0.28, "ACC": 0.33, "ACA": 0.23, "ACG": 0.16},
    "W": {"TGG": 1.0},
    "Y": {"TAT": 0.55, "TAC": 0.45},
    "V": {"GTT": 0.35, "GTC": 0.27, "GTA": 0.17, "GTG": 0.21},
    "*": {"TAA": 0.45, "TAG": 0.10, "TGA": 0.45},
}


def read_fasta(path: str) -> List[Tuple[str, str]]:
    op = gzip.open if path.endswith(".gz") else open
    seqs = []
    with op(path, "rt", encoding="utf-8") as fh:
        name = None
        cur = []
        for ln in fh:
            ln = ln.rstrip()
            if not ln:
                continue
            if ln.startswith(">"):
                if name is not None:
                    seqs.append((name, "".join(cur).upper()))
                name = ln[1:].split()[0]
                cur = []
            else:
                cur.append(ln)
        if name is not None:
            seqs.append((name, "".join(cur).upper()))
    return seqs


def load_codon_usage(path: Optional[str]) -> Dict[str, Dict[str, float]]:
    if not path:
        return GENERIC_PLANT_CODON_USAGE
    tbl = {}
    with open(path, "r", encoding="utf-8") as fh:
        rd = csv.DictReader(fh)
        # columns: aa,codon,freq  (freq ~ sum to 1 per aa)
        for row in rd:
            aa = row["aa"].strip().upper()
            codon = row["codon"].strip().upper()
            freq = float(row["freq"])
            tbl.setdefault(aa, {})[codon] = freq
    # normalize
    for aa, m in tbl.items():
        s = sum(m.values())
        if s > 0:
            for k in list(m.keys()):
                m[k] = m[k] / s
    if "*" not in tbl:
        tbl["*"] = GENERIC_PLANT_CODON_USAGE["*"].copy()
    return tbl


def rc(seq: str) -> str:
    return seq.translate(COMP)[::-1]


def sample_codon(rng: random.Random, codon_freqs: Dict[str, float]) -> str:
    r = rng.random()
    acc = 0.0
    for cod, w in codon_freqs.items():
        acc += w
        if r <= acc:
            return cod
    return list(codon_freqs.keys())[-1]


def reverse_translate_protein(
    aa_seq: str,
    rng: random.Random,
    codon_usage: Dict[str, Dict[str, float]],
    start_codon: str = "ATG",
    stop_choice: Optional[str] = None,
) -> str:
    cds = [start_codon]
    for i, aa in enumerate(aa_seq):
        if aa == "*":
            break
        aaU = aa.upper()
        if aaU == "M" and i == 0:
            continue
        pool = (
            codon_usage.get(aaU)
            or GENERIC_PLANT_CODON_USAGE.get(aaU)
            or GENERIC_PLANT_CODON_USAGE["G"]
        )
        cds.append(sample_codon(rng, pool))
    cds.append(stop_choice if stop_choice else sample_codon(rng, codon_usage["*"]))
    return "".join(cds)


def make_intron(rng: random.Random, mean: int, sd: int) -> str:
    # Truncated normal, min 50 bp; GT-AG intron
    for _ in range(1000):
        L = int(rng.gauss(mean, sd))
        if L >= 50:
            break
    else:
        L = max(50, mean)
    mid = "".join(rng.choice(DNA) for _ in range(max(0, L - 4)))
    return "GT" + mid + "AG"


def inject_introns_between_codons(
    cds: str, rng: random.Random, prob: float, mean: int, sd: int
) -> Tuple[str, List[Tuple[int, int]]]:
    assert len(cds) % 3 == 0
    parts = []
    introns: List[Tuple[int, int]] = []
    pos = 0
    for i in range(0, len(cds), 3):
        codon = cds[i : i + 3]
        parts.append(codon)
        pos += 3
        if i + 3 < len(cds) and rng.random() < prob:
            intr = make_intron(rng, mean, sd)
            intr_start = pos
            parts.append(intr)
            pos += len(intr)
            introns.append((intr_start, pos))
    return "".join(parts), introns


def random_dna(rng: random.Random, n: int) -> str:
    return "".join(rng.choice(DNA) for _ in range(n))


def write_fasta(path: str, records: List[Tuple[str, str]]):
    with open(path, "w", encoding="utf-8") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i : i + 80] + "\n")


# ---------- BUSCO handling ----------


def which(cmd: str) -> Optional[str]:
    return shutil.which(cmd)


def ensure_busco_available_or_install() -> None:
    """
    Ensure 'busco' is available in PATH.
    If missing and we're in a non-base conda env, attempt 'conda install -y busco'.
    If still missing OR conda env not activated (no CONDA_DEFAULT_ENV), raise SystemExit.
    """
    if which("busco"):
        return

    conda_env = os.environ.get("CONDA_DEFAULT_ENV", None)
    if not conda_env:
        sys.stderr.write(
            "[ERR] 'busco' not found and no conda environment is active (CONDA_DEFAULT_ENV is unset).\n"
        )
        sys.stderr.write(
            "      Please activate a non-base conda environment and install busco, e.g.:\n"
        )
        sys.stderr.write(
            "        conda create -n polap-bio -y python=3.11 && conda activate polap-bio && conda install -y busco\n"
        )
        raise SystemExit(1)

    if conda_env == "base":
        sys.stderr.write("[ERR] 'busco' not found and conda environment is 'base'.\n")
        sys.stderr.write(
            "      Please activate a non-base environment and install busco, e.g.:\n"
        )
        sys.stderr.write(
            "        conda create -n polap-bio -y python=3.11 && conda activate polap-bio && conda install -y busco\n"
        )
        raise SystemExit(1)

    # We are in a non-base env; try to install with conda
    if not which("conda"):
        sys.stderr.write(
            "[ERR] 'busco' not found and 'conda' is not on PATH; cannot auto-install.\n"
        )
        raise SystemExit(1)

    sys.stderr.write(
        f"[info] 'busco' not found; attempting 'conda install -y busco' in env '{conda_env}'...\n"
    )
    try:
        # Use a conservative install call; advanced channel config is left to user's conda setup
        res = subprocess.run(
            ["conda", "install", "-y", "busco"],
            check=False,
            stdout=sys.stderr,
            stderr=sys.stderr,
        )
        if res.returncode != 0:
            sys.stderr.write(
                "[ERR] 'conda install -y busco' failed. Please install BUSCO manually.\n"
            )
            raise SystemExit(1)
    except FileNotFoundError:
        sys.stderr.write(
            "[ERR] conda executable not found; cannot auto-install BUSCO.\n"
        )
        raise SystemExit(1)

    if not which("busco"):
        sys.stderr.write(
            "[ERR] 'busco' still not found after conda install. Please verify your conda channels and try again.\n"
        )
        raise SystemExit(1)


def get_or_download_busco_faa(lineage: str, cwd: str) -> str:
    """
    Return path to the BUSCO amino-acid FASTA for a lineage, downloading if needed.
    By default BUSCO stores into ./busco_downloads under the current working directory.
    """
    rel_path = os.path.join("busco_downloads", "lineages", lineage, "refseq_db.faa.gz")
    abs_path = os.path.join(cwd, rel_path)
    if os.path.exists(abs_path):
        return abs_path

    # Need to download
    ensure_busco_available_or_install()
    sys.stderr.write(
        f"[busco] downloading lineage '{lineage}' into ./busco_downloads ...\n"
    )
    res = subprocess.run(
        ["busco", "--download", lineage], cwd=cwd, stdout=sys.stderr, stderr=sys.stderr
    )
    if res.returncode != 0:
        sys.stderr.write("[ERR] 'busco --download {lineage}' failed.\n")
        raise SystemExit(1)

    if not os.path.exists(abs_path):
        sys.stderr.write(
            f"[ERR] Expected BUSCO file not found after download: {abs_path}\n"
        )
        sys.stderr.write(
            "      Run 'busco --list-datasets' to see available datasets.\n"
        )
        raise SystemExit(1)
    return abs_path


# ---------- main pipeline ----------


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Build a synthetic nuclear genome from BUSCO protein sequences.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument(
        "--busco-faa",
        default=None,
        help="Path to BUSCO amino-acid FASTA (e.g., busco_downloads/lineages/viridiplantae_odb12/refseq_db.faa.gz). "
        "If omitted, the script will attempt to download the lineage via 'busco --download <lineage>'.",
    )
    ap.add_argument(
        "--busco-lineage",
        default="viridiplantae_odb12",
        help="BUSCO lineage name to download when --busco-faa is not provided (default: viridiplantae_odb12). "
        "List datasets with 'busco --list-datasets'.",
    )
    ap.add_argument(
        "-o",
        "--out-dir",
        default="nuclear_from_busco_out",
        help="Output directory (will be created if missing). Files: genome.fa/gff3/bed/truth.tsv",
    )
    ap.add_argument(
        "--genome-size",
        type=int,
        default=10_000_000,
        help="Target genome size in bp (approximate; padded with intergenic DNA)",
    )
    ap.add_argument(
        "--n-chr", type=int, default=1, help="Number of chromosomes/scaffolds"
    )
    ap.add_argument(
        "--codon-usage",
        type=str,
        default=None,
        help="CSV with columns: aa,codon,freq (optional; overrides built-in)",
    )
    ap.add_argument("--utr5", type=int, default=60, help="5' UTR length (bp)")
    ap.add_argument("--utr3", type=int, default=120, help="3' UTR length (bp)")
    ap.add_argument(
        "--intron-prob",
        type=float,
        default=0.35,
        help="Probability of intron insertion at each codon junction",
    )
    ap.add_argument("--intron-mean", type=int, default=250, help="Mean intron length")
    ap.add_argument("--intron-sd", type=int, default=80, help="SD of intron length")
    ap.add_argument(
        "--spacer-mean", type=int, default=5000, help="Mean intergenic spacer length"
    )
    ap.add_argument(
        "--spacer-sd", type=int, default=2000, help="SD of intergenic spacer length"
    )
    ap.add_argument("--seed", type=int, default=17, help="Random seed")
    args = ap.parse_args(argv)

    # Ensure output dir exists
    os.makedirs(args.out_dir, exist_ok=True)

    # Resolve BUSCO proteins path
    if args.busco_faa:
        busco_faa_path = args.busco_faa
        if not os.path.exists(busco_faa_path):
            sys.stderr.write(f"[ERR] --busco-faa not found: {busco_faa_path}\n")
            raise SystemExit(1)
    else:
        # Use or download the lineage into the current working dir
        busco_faa_path = get_or_download_busco_faa(args.busco_lineage, cwd=os.getcwd())

    rng = random.Random(args.seed)
    codon_usage = load_codon_usage(args.codon_usage)

    # Read proteins
    prots = [(name, seq) for name, seq in read_fasta(busco_faa_path) if len(seq) >= 30]
    if not prots:
        sys.stderr.write(f"[ERR] No proteins >=30 aa found in {busco_faa_path}\n")
        return 2
    rng.shuffle(prots)

    # Build gene sequences + annotations
    genes = []
    for name, aa in prots:
        cds = reverse_translate_protein(aa, rng, codon_usage, start_codon="ATG")
        if args.intron_prob > 0.0:
            tx, introns = inject_introns_between_codons(
                cds, rng, args.intron_prob, args.intron_mean, args.intron_sd
            )
        else:
            tx, introns = cds, []
        utr5 = random_dna(rng, args.utr5) if args.utr5 > 0 else ""
        utr3 = random_dna(rng, args.utr3) if args.utr3 > 0 else ""
        gene_seq = utr5 + tx + utr3

        # exon/intron blocks in gene-local coords; UTRs are exons too
        exons = []
        cur = 0
        if args.utr5 > 0:
            exons.append((0, args.utr5))
        cur = args.utr5
        prev = 0
        for s, e in introns:
            exons.append((cur + prev, cur + s))  # exon before intron
            exons.append(("INTRON", cur + s, cur + e))  # intron block
            prev = e
        exons.append((cur + prev, cur + len(tx)))  # last exon
        if args.utr3 > 0:
            exons.append((cur + len(tx), cur + len(tx) + args.utr3))

        genes.append(
            {
                "name": name,
                "aa_len": len(aa),
                "cds_len": len(cds),
                "seq": gene_seq,
                "exon_intron": exons,
            }
        )

    # Lay out across chromosomes with spacers
    chrom_seqs = [""] * args.n_chr
    gff_rows = []
    bed_rows = []
    truth_rows = []

    total_len = 0
    chrom_idx = 0

    def draw_spacer() -> int:
        for _ in range(1000):
            L = int(rng.gauss(args.spacer_mean, args.spacer_sd))
            if L >= 0:
                return L
        return max(0, args.spacer_mean)

    for g in genes:
        strand = rng.choice(["+", "-"])
        seq = g["seq"] if strand == "+" else rc(g["seq"])

        spacer_len = 0 if len(chrom_seqs[chrom_idx]) == 0 else draw_spacer()
        spacer_seq = random_dna(rng, spacer_len)
        start = len(chrom_seqs[chrom_idx]) + spacer_len
        chrom_seqs[chrom_idx] += spacer_seq + seq
        end = len(chrom_seqs[chrom_idx])

        gene_id = f"gene_{g['name']}"
        chr_name = f"chr{chrom_idx+1}"

        # gene feature
        gff_rows.append(
            [
                chr_name,
                "polap-synth",
                "gene",
                start + 1,
                end,
                ".",
                strand,
                ".",
                f"ID={gene_id};Name={g['name']}",
            ]
        )
        bed_rows.append([chr_name, str(start), str(end), g["name"], "0", strand])

        gene_len = len(seq)

        def flip(a: int, b: int) -> Tuple[int, int]:
            return (gene_len - b, gene_len - a)

        for item in g["exon_intron"]:
            if isinstance(item, tuple) and len(item) == 2:
                a, b = item
                if strand == "+":
                    A = start + a
                    B = start + b
                else:
                    A_local, B_local = flip(a, b)
                    A = start + A_local
                    B = start + B_local
                gff_rows.append(
                    [
                        chr_name,
                        "polap-synth",
                        "exon",
                        A + 1,
                        B,
                        ".",
                        strand,
                        ".",
                        f"Parent={gene_id}",
                    ]
                )
            elif isinstance(item, tuple) and len(item) == 3 and item[0] == "INTRON":
                _, a, b = item
                if strand == "+":
                    A = start + a
                    B = start + b
                else:
                    A_local, B_local = flip(a, b)
                    A = start + A_local
                    B = start + B_local
                gff_rows.append(
                    [
                        chr_name,
                        "polap-synth",
                        "intron",
                        A + 1,
                        B,
                        ".",
                        strand,
                        ".",
                        f"Parent={gene_id}",
                    ]
                )

        truth_rows.append(
            [
                g["name"],
                chr_name,
                start,
                end,
                strand,
                g["aa_len"],
                g["cds_len"],
                len(seq),
            ]
        )

        total_len += spacer_len + len(seq)
        chrom_idx = (chrom_idx + 1) % args.n_chr
        if total_len >= args.genome_size:
            break

    # Pad to target
    remaining = max(0, args.genome_size - sum(len(s) for s in chrom_seqs))
    i = 0
    while remaining > 0 and i < args.n_chr:
        pad = min(remaining, draw_spacer())
        chrom_seqs[i] += random_dna(rng, pad)
        remaining -= pad
        i += 1

    # Write outputs
    fa_path = os.path.join(args.out_dir, "genome.fa")
    gff_path = os.path.join(args.out_dir, "genome.gff3")
    bed_path = os.path.join(args.out_dir, "genome.bed")
    tsv_path = os.path.join(args.out_dir, "genome.truth.tsv")

    # fasta
    with open(fa_path, "w", encoding="utf-8") as fh:
        for i, seq in enumerate(chrom_seqs):
            name = f"chr{i+1}"
            fh.write(f">{name}\n")
            for j in range(0, len(seq), 80):
                fh.write(seq[j : j + 80] + "\n")

    # gff3
    with open(gff_path, "w", encoding="utf-8") as gff:
        gff.write("##gff-version 3\n")
        for row in gff_rows:
            gff.write("\t".join([str(x) for x in row]) + "\n")

    # bed
    with open(bed_path, "w", encoding="utf-8") as bed:
        for row in bed_rows:
            bed.write("\t".join([str(x) for x in row]) + "\n")

    # truth
    with open(tsv_path, "w", encoding="utf-8") as fh:
        fh.write("gene\tchr\tstart\tend\tstrand\taa_len\tcds_len\tgene_len\n")
        for r in truth_rows:
            fh.write("\t".join(map(str, r)) + "\n")

    total_bp = sum(len(s) for s in chrom_seqs)
    sys.stderr.write(f"[nuclear] wrote {fa_path} ({total_bp} bp, {args.n_chr} chr)\n")
    sys.stderr.write(f"[nuclear] annotations: {gff_path} / {bed_path} / {tsv_path}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
