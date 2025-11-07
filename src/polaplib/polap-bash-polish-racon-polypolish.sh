#!/usr/bin/env bash
# Version: v0.8.0
# Name: polap-bash-polish-racon-polypolish.sh
#
# MTPT-aware organelle polishing with:
#   1) Competitive read assignment (target + other organelle)
#   2) Long-read polish (racon xN; Medaka optional, OFF by default)
#   3) cp↔mt transfer detection → mask + reports
#   4) Short-read polish (Polypolish) only in mask-complement regions
#   5) Coverage diagnostics (all ONT vs assigned ONT) with MTPT shading
#
# Polypolish aligner choices (for all-per-read mapping): bowtie2 (default) | bwa | minimap2
#   bowtie2 : --very-sensitive -a
#   bwa-mem2: -a
#   minimap2: -x sr -a
#
# Helpers expected in: $POLAPLIB_DIR/scripts/
#   polap_py_paf_filter.py
#   polap_py_mask_from_paf.py
#   polap_py_bed_complement.py
#   polap_py_paf_to_hits_tsv.py
#   polap_py_bed_add_len.py
#   polap_py_assign_ont_ids.py
#   polap_py_assign_sr_ids.py
#   polap_py_fastq_filter_by_ids.py
#   polap_py_bam_coverage_bins.py
#   polap_r_covplot_mtpt_overlay.R
#   polap_r_mask_qc.R
#
set -Eeuo pipefail

# ---------- Paths ----------
POLAPLIB_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd -P)"
SCRIPTS_DIR="${POLAPLIB_DIR}/scripts"

PAF_FILTER="${SCRIPTS_DIR}/polap_py_paf_filter.py"
MASK_FROM_PAF="${SCRIPTS_DIR}/polap_py_mask_from_paf.py"
BED_COMPLEMENT="${SCRIPTS_DIR}/polap_py_bed_complement.py"
PAF_TO_HITS="${SCRIPTS_DIR}/polap_py_paf_to_hits_tsv.py"
BED_ADD_LEN="${SCRIPTS_DIR}/polap_py_bed_add_len.py"
ASSIGN_ONT_IDS="${SCRIPTS_DIR}/polap_py_assign_ont_ids.py"
ASSIGN_SR_IDS="${SCRIPTS_DIR}/polap_py_assign_sr_ids.py"
FASTQ_FILTER_BY_IDS="${SCRIPTS_DIR}/polap_py_fastq_filter_by_ids.py"
COV_BINS="${SCRIPTS_DIR}/polap_py_bam_coverage_bins.py"
COV_PLOT="${SCRIPTS_DIR}/polap_r_covplot_mtpt_overlay.R"
MASK_QC_R="${SCRIPTS_DIR}/polap_r_mask_qc.R"

# ---------- Defaults ----------
TARGET="mt" # mt | cp
FASTA=""
OTHER=""
ONT=""
SR1=""
SR2=""
OUTDIR=""
THREADS=16

# LR polish
RACON_ROUNDS=2
ONT_PRESET="map-ont"
PAF_MIN_IDENT=0.91
PAF_MIN_ALEN=3000
USE_MEDAKA=0 # default OFF
MEDAKA_MODEL="r941_min_sup_g507"

# Mask detection (cp→mt for mt target; mt→cp for cp target)
MASK_POLICY="auto" # auto:on for mt, off for cp
MASK_MIN_IDENT=0.85
MASK_MIN_LEN=150
MASK_PAD=1000
RUN_MASK_QC=1

# Coverage diagnostics
COV_DISABLE=0
COV_BIN=200
COV_MIN_MAPQ=0

# Assignment thresholds
ASSIGN_MAPQ=30

# Polypolish aligner choice
PP_ALIGNER="bowtie2" # bowtie2 | bwa | minimap2

# ---------- Logging ----------
ts() { date "+%Y-%m-%d %H:%M:%S"; }
log() { echo "[$(ts)] $*" >&2; }
die() {
	echo "[$(ts)] ERROR: $*" >&2
	exit 2
}

usage() {
	cat <<EOF
Usage:
  $(basename "$0") --target mt|cp --fasta asm.fa --other other.fa \\
    [--ont ont.fq.gz] [--sr1 R1.fq.gz --sr2 R2.fq.gz] \\
    --outdir OUT [--threads 16] \\
    [--racon-rounds 2] [--ont-preset map-ont] [--use-medaka 0|1] [--medaka-model MODEL] \\
    [--paf-min-ident 0.91] [--paf-min-alen 3000] \\
    [--mask auto|on|off] [--mask-min-ident 0.85] [--mask-min-len 150] [--mask-pad 1000] [--mask-qc 0|1] \\
    [--cov-disable 0|1] [--cov-bin 200] [--cov-min-mapq 0] \\
    [--assign-mapq 30] \\
    [--pp-aligner bowtie2|bwa|minimap2]

Outputs:
  polishing:  <outdir>/polished.fa
  masks:      <outdir>/mask/transfer.mask.bed + .hits.tsv + .merged.tsv (+ mask_qc.tsv/.pdf)
  coverage:   <outdir>/coverage/coverage_overlay.pdf + coverage_summary.tsv
EOF
}

# ---------- Parse CLI ----------
while (($#)); do
	case "$1" in
	--target)
		TARGET="${2:?}"
		shift 2
		;;
	--fasta)
		FASTA="${2:?}"
		shift 2
		;;
	--other)
		OTHER="${2:-}"
		shift 2
		;;
	--ont)
		ONT="${2:-}"
		shift 2
		;;
	--sr1)
		SR1="${2:-}"
		shift 2
		;;
	--sr2)
		SR2="${2:-}"
		shift 2
		;;
	--outdir)
		OUTDIR="${2:?}"
		shift 2
		;;
	--threads)
		THREADS="${2:?}"
		shift 2
		;;

	--racon-rounds)
		RACON_ROUNDS="${2:?}"
		shift 2
		;;
	--ont-preset)
		ONT_PRESET="${2:?}"
		shift 2
		;;
	--paf-min-ident)
		PAF_MIN_IDENT="${2:?}"
		shift 2
		;;
	--paf-min-alen)
		PAF_MIN_ALEN="${2:?}"
		shift 2
		;;
	--use-medaka)
		USE_MEDAKA="${2:?}"
		shift 2
		;;
	--medaka-model)
		MEDAKA_MODEL="${2:?}"
		shift 2
		;;

	--mask)
		MASK_POLICY="${2:?}"
		shift 2
		;;
	--mask-min-ident)
		MASK_MIN_IDENT="${2:?}"
		shift 2
		;;
	--mask-min-len)
		MASK_MIN_LEN="${2:?}"
		shift 2
		;;
	--mask-pad)
		MASK_PAD="${2:?}"
		shift 2
		;;
	--mask-qc)
		RUN_MASK_QC="${2:?}"
		shift 2
		;;

	--cov-disable)
		COV_DISABLE="${2:?}"
		shift 2
		;;
	--cov-bin)
		COV_BIN="${2:?}"
		shift 2
		;;
	--cov-min-mapq)
		COV_MIN_MAPQ="${2:?}"
		shift 2
		;;

	--assign-mapq)
		ASSIGN_MAPQ="${2:?}"
		shift 2
		;;
	--pp-aligner)
		PP_ALIGNER="${2:?}"
		shift 2
		;;

	-h | --help)
		usage
		exit 0
		;;
	*) die "Unknown option: $1" ;;
	esac
done

[[ -n "$FASTA" && -n "$OUTDIR" ]] || {
	usage
	die "Required: --fasta and --outdir"
}
[[ "$TARGET" == "mt" || "$TARGET" == "cp" ]] || die "--target must be mt or cp"
[[ -s "$FASTA" ]] || die "FASTA not found: $FASTA"
[[ -n "$OTHER" ]] && [[ ! -s "$OTHER" ]] && die "Other FASTA not found: $OTHER"

mkdir -p "$OUTDIR"/{stage0_assign,stage1_racon,mask,stage3_polypolish,coverage,final,logs}

# ---------- deps ----------
need() { command -v "$1" >/dev/null 2>&1 || die "Missing dependency: $1"; }
need samtools
need minimap2
need racon
need polypolish
need python3
# optional aligners for Polypolish
case "$PP_ALIGNER" in
bowtie2) need bowtie2 ;;
bwa) need bwa-mem2 ;;
minimap2) ;; # already required
*) die "--pp-aligner must be bowtie2|bwa|minimap2" ;;
esac
# helpers present?
for f in "$PAF_FILTER" "$MASK_FROM_PAF" "$BED_COMPLEMENT" "$PAF_TO_HITS" "$BED_ADD_LEN" \
	"$ASSIGN_ONT_IDS" "$ASSIGN_SR_IDS" "$FASTQ_FILTER_BY_IDS" "$COV_BINS" "$COV_PLOT"; do
	[[ -s "$f" ]] || die "Helper missing: $f"
done
# pysam check
python3 - <<'PY' || die "Python module 'pysam' not found (pip install pysam)"
import importlib, sys; importlib.import_module("pysam"); sys.exit(0)
PY

HAVE_SEQKIT=0
command -v seqkit >/dev/null 2>&1 && HAVE_SEQKIT=1
if [[ "$USE_MEDAKA" -eq 1 ]] && ! command -v medaka_consensus >/dev/null 2>&1; then
	log "Medaka requested but not found; proceeding without Medaka"
	USE_MEDAKA=0
fi

filter_fastq_by_ids() {
	local fq="$1" ids="$2" out="$3"
	if [[ "$HAVE_SEQKIT" -eq 1 ]]; then
		seqkit grep -f "$ids" "$fq" -o "$out" 2>>"$OUTDIR/logs/seqkit.err"
	else
		python3 "$FASTQ_FILTER_BY_IDS" --fastq "$fq" --ids "$ids" --out "$out"
	fi
}

# ---------- Combined reference ----------
COMBINED="$OUTDIR/stage0_assign/combined.fa"
TGT_PREFIX="${TARGET}|"
OTH_PREFIX="$([[ "$TARGET" == "mt" ]] && echo "cp|" || echo "mt|")"
{
	awk -vP="$TGT_PREFIX" '/^>/{sub(/^>/,">"P)}1' "$FASTA"
	if [[ -n "$OTHER" ]]; then awk -vP="$OTH_PREFIX" '/^>/{sub(/^>/,">"P)}1' "$OTHER"; fi
} >"$COMBINED"
samtools faidx "$COMBINED" >/dev/null

# ---------- Stage 0: assignment ----------
ASSIGN_DIR="$OUTDIR/stage0_assign"
LR_ASSIGNED="$ASSIGN_DIR/ont.target.fq.gz"
SR1_ASSIGNED="$ASSIGN_DIR/sr.target.R1.fq.gz"
SR2_ASSIGNED="$ASSIGN_DIR/sr.target.R2.fq.gz"
HAVE_LR=0
HAVE_SR=0

if [[ -n "$ONT" ]]; then
	log "Assign ONT → target (MAPQ≥$ASSIGN_MAPQ)"
	minimap2 -ax "$ONT_PRESET" -t "$THREADS" "$COMBINED" "$ONT" |
		samtools view -h -F 0x900 - |
		python3 "$ASSIGN_ONT_IDS" --sam - --target-prefix "$TGT_PREFIX" --mapq "$ASSIGN_MAPQ" \
			>"$ASSIGN_DIR/ont.target.ids"
	n_ont=$(wc -l <"$ASSIGN_DIR/ont.target.ids" || echo 0)
	[[ "$n_ont" -gt 0 ]] || die "No ONT reads assigned; adjust --assign-mapq or inputs"
	filter_fastq_by_ids "$ONT" "$ASSIGN_DIR/ont.target.ids" "$LR_ASSIGNED"
	HAVE_LR=1
else
	log "No ONT provided; LR polish & coverage will be skipped"
fi

if [[ -n "$SR1" && -n "$SR2" ]]; then
	log "Assign Illumina → target (both ends MAPQ≥$ASSIGN_MAPQ)"
	mkdir -p "$ASSIGN_DIR/bt2idx"
	bowtie2-build "$COMBINED" "$ASSIGN_DIR/bt2idx/comb" >"$OUTDIR/logs/assign_bt2_build.out" 2>"$OUTDIR/logs/assign_bt2_build.err"
	bowtie2 --very-sensitive -x "$ASSIGN_DIR/bt2idx/comb" -1 "$SR1" -2 "$SR2" -p "$THREADS" 2>"$OUTDIR/logs/assign_bt2.err" |
		samtools view -bh - |
		python3 "$ASSIGN_SR_IDS" --bam - --target-prefix "$TGT_PREFIX" --mapq "$ASSIGN_MAPQ" \
			>"$ASSIGN_DIR/sr.target.ids"
	n_sr=$(wc -l <"$ASSIGN_DIR/sr.target.ids" || echo 0)
	if [[ "$n_sr" -gt 0 ]]; then
		filter_fastq_by_ids "$SR1" "$ASSIGN_DIR/sr.target.ids" "$SR1_ASSIGNED"
		filter_fastq_by_ids "$SR2" "$ASSIGN_DIR/sr.target.ids" "$SR2_ASSIGNED"
		HAVE_SR=1
	else
		log "WARNING: no short-read pairs confidently assigned; SR polish will be skipped."
	fi
else
	log "No SR provided; SR polish will be skipped"
fi

((HAVE_LR == 1 || HAVE_SR == 1)) || die "No reads assigned; nothing to polish"

# ---------- Stage 1: Long-read polish ----------
CUR="$FASTA"
if ((HAVE_LR == 1)); then
	log "Stage 1: racon x${RACON_ROUNDS} (preset=${ONT_PRESET}; NO mask at LR stage)"
	for r in $(seq 1 "$RACON_ROUNDS"); do
		RAW="$OUTDIR/stage1_racon/r${r}.raw.paf"
		FLT="$OUTDIR/stage1_racon/r${r}.flt.paf"
		OUTFA="$OUTDIR/stage1_racon/polished.r${r}.fa"
		minimap2 -x "$ONT_PRESET" -t "$THREADS" -c --secondary=no "$CUR" "$LR_ASSIGNED" >"$RAW" 2>"$OUTDIR/logs/r${r}.minimap2.err"
		python3 "$PAF_FILTER" --min-ident "$PAF_MIN_IDENT" --min-alen "$PAF_MIN_ALEN" --in "$RAW" --out "$FLT" 2>>"$OUTDIR/logs/r${r}.paf_filter.err"
		[[ -s "$FLT" ]] || {
			log "WARN: filtered PAF empty; fallback to RAW"
			cp -f "$RAW" "$FLT"
		}
		racon -t "$THREADS" "$LR_ASSIGNED" "$FLT" "$CUR" >"$OUTFA" 2>"$OUTDIR/logs/r${r}.racon.err"
		CUR="$OUTFA"
	done
	if [[ "$USE_MEDAKA" -eq 1 ]] && command -v medaka_consensus >/dev/null 2>&1; then
		log "Stage 1b: medaka (model=$MEDAKA_MODEL)"
		MED_DIR="$OUTDIR/stage1_racon/medaka"
		mkdir -p "$MED_DIR"
		minimap2 -ax "$ONT_PRESET" -t "$THREADS" "$CUR" "$LR_ASSIGNED" | samtools sort -@ "$THREADS" -o "$MED_DIR/reads.sorted.bam" -
		samtools index -@ "$THREADS" "$MED_DIR/reads.sorted.bam"
		medaka_consensus -i "$LR_ASSIGNED" -d "$CUR" -o "$MED_DIR" -m "$MEDAKA_MODEL" -t "$THREADS" \
			1>"$OUTDIR/logs/medaka.out" 2>"$OUTDIR/logs/medaka.err" || log "Medaka failed; keep racon result"
		[[ -s "$MED_DIR/consensus.fasta" ]] && CUR="$MED_DIR/consensus.fasta"
	fi
fi

# ---------- Stage 2: Mask detection + reports ----------
if [[ "$MASK_POLICY" == "auto" ]]; then
	if [[ "$TARGET" == "mt" ]]; then MASK_POLICY="on"; else MASK_POLICY="off"; fi
fi
MASK_DIR="$OUTDIR/mask"
mkdir -p "$MASK_DIR"
MASK_BED="$MASK_DIR/transfer.mask.bed"
ALLOW_BED="$MASK_DIR/allow.bed"
A2A_PAF="$MASK_DIR/a2a.raw.paf"
A2A_FILT="$MASK_DIR/a2a.filt.paf"
MASK_HITS_TSV="$MASK_DIR/transfer.mask.hits.tsv"
MASK_MERGED_TSV="$MASK_DIR/transfer.mask.merged.tsv"

if [[ "$MASK_POLICY" == "off" || -z "$OTHER" || ! -s "$OTHER" ]]; then
	log "Stage 2: mask OFF or --other missing; ALLOW = whole assembly"
	python3 "$BED_COMPLEMENT" --fasta "$CUR" --bed /dev/null --out-bed "$ALLOW_BED" >/dev/null 2>&1
	: >"$MASK_BED"
else
	log "Stage 2: detect $([[ "$TARGET" == "mt" ]] && echo 'cp→mt' || echo 'mt→cp')) transfers"
	minimap2 -x asm10 -t "$THREADS" -c --secondary=no "$CUR" "$OTHER" >"$A2A_PAF" 2>"$OUTDIR/logs/a2a.minimap2.err"
	python3 "$PAF_FILTER" --min-ident "$MASK_MIN_IDENT" --min-alen "$MASK_MIN_LEN" --in "$A2A_PAF" --out "$A2A_FILT" 2>"$OUTDIR/logs/a2a.filter.err"
	python3 "$MASK_FROM_PAF" --paf "$A2A_FILT" --fasta "$CUR" --pad "$MASK_PAD" --out-bed "$MASK_BED" 2>"$OUTDIR/logs/mask_from_paf.err"
	python3 "$BED_COMPLEMENT" --fasta "$CUR" --bed "$MASK_BED" --out-bed "$ALLOW_BED" 2>"$OUTDIR/logs/bed_complement.err"
	python3 "$PAF_TO_HITS" --paf "$A2A_FILT" --out "$MASK_HITS_TSV" 2>"$OUTDIR/logs/paf_to_hits.err"
	python3 "$BED_ADD_LEN" --fasta "$CUR" --bed "$MASK_BED" --out "$MASK_MERGED_TSV" 2>"$OUTDIR/logs/bed_add_len.err"
	if [[ "$RUN_MASK_QC" -eq 1 ]] && command -v Rscript >/dev/null 2>&1 && [[ -s "$MASK_QC_R" ]]; then
		Rscript --vanilla "$MASK_QC_R" --fasta "$CUR" --bed "$MASK_BED" \
			--outpdf "$MASK_DIR/mask_qc.pdf" --outtsv "$MASK_DIR/mask_qc.tsv" >/dev/null 2>&1 || true
	fi
fi

# ---------- Stage 2.5: Coverage diagnostics ----------
if ((COV_DISABLE == 0)) && [[ -n "$ONT" ]] && ((HAVE_LR == 1)); then
	log "Stage 2.5: coverage (all ONT vs assigned ONT), bin=${COV_BIN}, minMAPQ=${COV_MIN_MAPQ}"
	COV_DIR="$OUTDIR/coverage"
	mkdir -p "$COV_DIR"

	minimap2 -ax "$ONT_PRESET" -t "$THREADS" "$CUR" "$ONT" |
		samtools sort -@ "$THREADS" -o "$COV_DIR/allONT.sorted.bam" -

	samtools index -@ "$THREADS" "$COV_DIR/allONT.sorted.bam"

	minimap2 -ax "$ONT_PRESET" -t "$THREADS" "$CUR" "$ASSIGN_DIR/ont.target.fq.gz" |
		samtools sort -@ "$THREADS" -o "$COV_DIR/assignedONT.sorted.bam" -

	samtools index -@ "$THREADS" "$COV_DIR/assignedONT.sorted.bam"

	python3 "$COV_BINS" --fasta "$CUR" --bam "$COV_DIR/allONT.sorted.bam" \
		--bin "$COV_BIN" --min-mapq "$COV_MIN_MAPQ" --source allONT \
		--out-tsv "$COV_DIR/coverage.allONT.bin${COV_BIN}.tsv"

	python3 "$COV_BINS" --fasta "$CUR" --bam "$COV_DIR/assignedONT.sorted.bam" \
		--bin "$COV_BIN" --min-mapq "$COV_MIN_MAPQ" --source assignedONT \
		--out-tsv "$COV_DIR/coverage.assignedONT.bin${COV_BIN}.tsv"

	if command -v Rscript >/dev/null 2>&1 && [[ -s "$COV_PLOT" ]]; then
		Rscript --vanilla "$COV_PLOT" \
			--fasta "$CUR" \
			--cov-all "$COV_DIR/coverage.allONT.bin${COV_BIN}.tsv" \
			--cov-assigned "$COV_DIR/coverage.assignedONT.bin${COV_BIN}.tsv" \
			--mask "$MASK_BED" \
			--bin "$COV_BIN" \
			--out-pdf "$COV_DIR/coverage_overlay.pdf" \
			--out-tsv "$COV_DIR/coverage_summary.tsv" \
			--title "$(basename "$OUTDIR"): coverage (all ONT vs assigned ONT)"
	fi
else
	log "Stage 2.5: coverage skipped"
fi

# ---------- Stage 3: SR Polypolish (choose aligner) ----------
FINAL="$CUR"
if [[ -n "${SR1_ASSIGNED:-}" && -s "$SR1_ASSIGNED" && -n "${SR2_ASSIGNED:-}" && -s "$SR2_ASSIGNED" ]]; then
	log "Stage 3: Polypolish using --pp-aligner=$PP_ALIGNER (mask enforced)"
	PP_DIR="$OUTDIR/stage3_polypolish"
	mkdir -p "$PP_DIR/bt2idx"

	case "$PP_ALIGNER" in
	bowtie2)
		bowtie2-build "$CUR" "$PP_DIR/bt2idx/asm" >"$OUTDIR/logs/pp_build.out" 2>"$OUTDIR/logs/pp_build.err"
		bowtie2 --very-sensitive -a -x "$PP_DIR/bt2idx/asm" -1 "$SR1_ASSIGNED" -2 "$SR2_ASSIGNED" -p "$THREADS" 2>"$OUTDIR/logs/pp_bt2.err" |
			samtools view -bh - |
			samtools view -bh -L "$ALLOW_BED" - |
			samtools sort -n -@ "$THREADS" -o "$PP_DIR/short.name.bam" - 2>"$OUTDIR/logs/pp_sortn.err"
		;;
	bwa)
		bwa-mem2 index "$CUR" >"$OUTDIR/logs/pp_bwa_index.out" 2>"$OUTDIR/logs/pp_bwa_index.err"
		bwa-mem2 mem -t "$THREADS" -a "$CUR" "$SR1_ASSIGNED" "$SR2_ASSIGNED" 2>"$OUTDIR/logs/pp_bwa.err" |
			samtools view -bh - |
			samtools view -bh -L "$ALLOW_BED" - |
			samtools sort -n -@ "$THREADS" -o "$PP_DIR/short.name.bam" - 2>"$OUTDIR/logs/pp_sortn.err"
		;;
	minimap2)
		minimap2 -x sr -a -t "$THREADS" "$CUR" "$SR1_ASSIGNED" "$SR2_ASSIGNED" 2>"$OUTDIR/logs/pp_mmsr.err" |
			samtools view -bh - |
			samtools view -bh -L "$ALLOW_BED" - |
			samtools sort -n -@ "$THREADS" -o "$PP_DIR/short.name.bam" - 2>"$OUTDIR/logs/pp_sortn.err"
		;;
	esac

	POLY_OUT="$PP_DIR/polished.fa"
	samtools view "$PP_DIR/short.name.bam" | polypolish polish "$CUR" /dev/stdin >"$POLY_OUT" 2>"$OUTDIR/logs/polypolish.err" || true
	[[ -s "$POLY_OUT" ]] && FINAL="$POLY_OUT"
else
	log "Stage 3: SR polish skipped (no confidently assigned pairs)"
fi

cp -f "$FINAL" "$OUTDIR/polished.fa"
log "DONE: $OUTDIR/polished.fa"
log "Summary: TARGET=$TARGET  LR=$([[ -n "$ONT" ]] && echo yes || echo no)  SR=$([[ -n "$SR1" && -n "$SR2" ]] && echo yes || echo no)"
log "Mask: policy=$MASK_POLICY  thresholds: ident>=$MASK_MIN_IDENT len>=$MASK_MIN_LEN pad=${MASK_PAD}bp"
log "Coverage: bin=$COV_BIN minMAPQ=$COV_MIN_MAPQ status=$([[ $COV_DISABLE -eq 0 && -n $ONT ]] && echo on || echo off)"
log "Polypolish aligner: $PP_ALIGNER"
