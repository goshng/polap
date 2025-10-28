#!/usr/bin/env bash
# Version: v0.5.1
# Name: polap-bash-polish-longreads-allreads-qc.sh
# Function: Perform two rounds of Racon polishing using ALL reads (filtered PAF by identity and aligned length),
#            then assess polishing quality (depth stats, violin plot, alignment QV via samtools stats).
#
# Usage:
#   polap-bash-polish-longreads-allreads-qc.sh \
#     --reads all_reads.fq.gz \
#     --fasta organelle.fa \
#     --outdir outdir \
#     [--platform ont|hifi] [--threads 32] \
#     [--rounds 2] [--min-ident 0.91] [--min-alen 3000]
#
# Dependencies:
#   minimap2, racon, samtools, python3
# Optional:
#   Rscript with ggplot2 for violin plots
#
# Helpers:
#   scripts/polap-py-filter-paf-by-iden-alen.py
#   scripts/polap-py-depth-summary.py
#   scripts/polap-py-mapping-summary-from-stats.py
#   scripts/polap-py-merge-kv-tables.py
#   scripts/polap-r-coverage-violin.R

set -Eeuo pipefail

# ------------------ Default parameters ------------------
READS=""
FASTA=""
OUT=""
PLAT="ont" # ont | hifi
THREADS=32
ROUNDS=2
MIN_IDENT=0.91
MIN_ALEN=3000

# ------------------ Logging ------------------
VERB=1
log() {
	local lvl="$1"
	shift
	((VERB >= lvl)) && printf '[%(%F %T)T] %s\n' -1 "$*" >&2 || true
}
die() {
	log 0 "ERR: $*"
	exit 2
}

# ------------------ CLI ------------------
print_help() {
	cat <<EOF
Usage: $0 --reads R.fq.gz --fasta A.fa --outdir out/ [options]
  --platform {ont|hifi}   Minimap2 preset (default: ont)
  --threads N             Threads (default: 32)
  --rounds N              Racon rounds (default: 2)
  --min-ident F           PAF filter minimum identity (default: 0.91)
  --min-alen N            PAF filter minimum aligned length (default: 3000)
  -v                      Increase verbosity (can repeat)
  -h|--help               Show this help
EOF
}

while (("$#")); do
	case "$1" in
	--reads)
		READS="${2:?}"
		shift 2
		;;
	--fasta)
		FASTA="${2:?}"
		shift 2
		;;
	--outdir)
		OUT="${2:?}"
		shift 2
		;;
	--platform)
		PLAT="${2:?}"
		shift 2
		;;
	--threads)
		THREADS="${2:?}"
		shift 2
		;;
	--rounds)
		ROUNDS="${2:?}"
		shift 2
		;;
	--min-ident)
		MIN_IDENT="${2:?}"
		shift 2
		;;
	--min-alen)
		MIN_ALEN="${2:?}"
		shift 2
		;;
	-v)
		VERB=$((VERB + 1))
		shift
		;;
	-h | --help)
		print_help
		exit 0
		;;
	*) die "Unknown option: $1" ;;
	esac
done

[[ -n "$READS" && -n "$FASTA" && -n "$OUT" ]] || {
	print_help
	die "Required: --reads, --fasta, --outdir"
}

# ------------------ Dependency checks ------------------
for tool in minimap2 racon samtools python3; do
	command -v "$tool" >/dev/null 2>&1 || die "Missing dependency: $tool"
done
HAVE_R=0
command -v Rscript >/dev/null 2>&1 && HAVE_R=1

# ------------------ Paths ------------------
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
S="${ROOT_DIR}/scripts"

PY_FILT="${S}/polap-py-filter-paf-by-iden-alen.py"
PY_DEPTH="${S}/polap-py-depth-summary.py"
PY_MAPSUM="${S}/polap-py-mapping-summary-from-stats.py"
PY_MERGE="${S}/polap-py-merge-kv-tables.py"
R_VIOL="${S}/polap-r-coverage-violin.R"

for f in "$PY_FILT" "$PY_DEPTH" "$PY_MAPSUM" "$PY_MERGE"; do
	[[ -s "$f" ]] || die "Missing helper: $f"
done

mkdir -p "$OUT"

# ------------------ Platform preset ------------------
case "$PLAT" in
ont) MM_PRESET="map-ont" ;;
hifi) MM_PRESET="map-hifi" ;;
*) die "--platform must be ont or hifi" ;;
esac

# --- Rounds of Racon polishing (primary-only overlaps) ---
CUR="$FASTA"
for r in $(seq 1 "$ROUNDS"); do
	RAW="$OUT/r${r}.raw.paf"
	FLT="$OUT/r${r}.flt.paf"
	OUTFA="$OUT/polished.r${r}.fa"

	log 1 "[Round $r] Mapping all reads -> PAF (primary-only emit)"
	minimap2 -x "$MM_PRESET" -t "$THREADS" -c --secondary=no "$CUR" "$READS" \
		>"$RAW" 2>"$RAW.err"

	log 1 "[Round $r] Filtering PAF (ident>=$MIN_IDENT, alen>=$MIN_ALEN, primary-only)"
	python3 "$PY_FILT" --paf "$RAW" --min-ident "$MIN_IDENT" --min-alen "$MIN_ALEN" --out "$FLT" \
		2>"$FLT.err"

	log 1 "[Round $r] Racon polishing"
	racon -t "$THREADS" "$READS" "$FLT" "$CUR" >"$OUTFA" 2>"$OUTFA.err"

	CUR="$OUTFA"
done
POL="$CUR"

log 1 "[Polish] Final polished FASTA: $POL"

# ------------------ 2. QC: mapping, depth, stats, violin, QV ------------------
BAM="$OUT/reads_to_polished.bam"
log 1 "[QC] Mapping reads -> BAM (primary-only)"
minimap2 -x "$MM_PRESET" -a -t "$THREADS" --secondary=no "$POL" "$READS" \
	2>"$OUT/mapping.err" |
	samtools view -F 0x900 -u |
	samtools sort -@ "$THREADS" -o "$BAM" -
samtools index "$BAM"

log 1 "[QC] Computing per-base depth"
samtools depth "$BAM" >"$OUT/depth.tsv" 2>"$OUT/depth.err"

log 1 "[QC] Summarizing depth"
python3 "$PY_DEPTH" --depth "$OUT/depth.tsv" --out "$OUT/depth.summary.tsv" \
	2>"$OUT/depth.summary.err"

log 1 "[QC] Generating samtools stats"
samtools stats "$BAM" >"$OUT/mapping.stats" 2>"$OUT/mapping.stats.err"
python3 "$PY_MAPSUM" --stats "$OUT/mapping.stats" --out "$OUT/mapping.summary.tsv" \
	2>"$OUT/mapping.summary.err"

log 1 "[QC] Merging summary tables"
python3 "$PY_MERGE" --tables "$OUT/depth.summary.tsv" "$OUT/mapping.summary.tsv" \
	--out "$OUT/polish_summary.tsv" 2>"$OUT/polish_summary.err"

if ((HAVE_R)) && [[ -s "$R_VIOL" ]]; then
	log 1 "[QC] Generating violin plot"
	Rscript "$R_VIOL" "$OUT/depth.tsv" "$OUT/coverage_violin.png" "Depth after polishing" \
		>"$OUT/violin.log" 2>&1
else
	log 1 "[QC] Skipping violin plot (Rscript or script missing)"
fi

log 1 "[DONE] Final polished FASTA: $POL"
log 1 "       Summary table: $OUT/polish_summary.tsv"
log 2 "       Depth TSV: $OUT/depth.tsv"
log 2 "       BAM: $BAM"
