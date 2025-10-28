#!/usr/bin/env bash
# Version: v0.2.0
# Name: polap-bash-compare-polish-longreads.sh
# Purpose: Compare two polished FASTA assemblies using the SAME long-read dataset.
#          Computes mapping-based QV (via samtools stats), identity, depth stats,
#          and renders violin plots for each; outputs a side-by-side comparison table.
#
# Usage:
#   polap-bash-compare-polish-longreads.sh \
#     --reads all_reads.fq.gz \
#     --fastaA polishedA.fa \
#     --fastaB polishedB.fa \
#     --outdir out_compare \
#     [--platform ont|hifi] [--threads 32] [-v ...]
#
# Outputs (in out_compare/):
#   A/reads_to_asm.{bam,bai}, A/depth.tsv, A/mapping.stats, A/coverage_violin.png
#   A/depth.summary.tsv, A/mapping.summary.tsv, A/summary.tsv
#   B/... (same as A)
#   compare.tsv  (metric  A  B  delta(B-A)  better)
#
# Notes:
#   - Primary-only QC mapping (--secondary=no) to avoid identity inflation from secondary/supplementary hits.
#   - ASCII-only logs; all stderr is redirected to *.err or *.log files.

set -Eeuo pipefail

# ---------------- Defaults ----------------
READS=""
FASTA_A=""
FASTA_B=""
OUT=""
PLAT="ont" # ont | hifi
THREADS=32

# ---------------- Logging ----------------
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

# ---------------- CLI ----------------
print_help() {
	cat <<EOF
Usage: $0 --reads R.fq.gz --fastaA A.fa --fastaB B.fa --outdir out/ [options]
  --platform {ont|hifi}   Minimap2 preset for QC mapping (default: ont)
  --threads N             Threads (default: 32)
  -v                      Increase verbosity (repeatable)
  -h|--help               Show help
EOF
}

while (("$#")); do
	case "$1" in
	--reads)
		READS="${2:?}"
		shift 2
		;;
	--fastaA)
		FASTA_A="${2:?}"
		shift 2
		;;
	--fastaB)
		FASTA_B="${2:?}"
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

[[ -n "$READS" && -n "$FASTA_A" && -n "$FASTA_B" && -n "$OUT" ]] || {
	print_help
	die "Required: --reads --fastaA --fastaB --outdir"
}

# ---------------- Deps ----------------
for t in minimap2 samtools python3; do
	command -v "$t" >/dev/null 2>&1 || die "Missing dependency: $t"
done
HAVE_R=0
command -v Rscript >/dev/null 2>&1 && HAVE_R=1

# ---------------- Paths ----------------
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
S="${ROOT_DIR}/scripts"

PY_DEPTH="${S}/polap-py-depth-summary.py"
PY_MAPSUM="${S}/polap-py-mapping-summary-from-stats.py"
PY_MERGE="${S}/polap-py-merge-kv-tables.py"
PY_COMPARE="${S}/polap-py-compare-polish-metrics.py"
R_VIOL="${S}/polap-r-coverage-violin.R"

for f in "$PY_DEPTH" "$PY_MAPSUM" "$PY_MERGE" "$PY_COMPARE"; do
	[[ -s "$f" ]] || die "Missing helper: $f"
done

mkdir -p "$OUT"/{A,B}

# ---------------- Preset ----------------
case "$PLAT" in
ont) MM_PRESET="map-ont" ;;
hifi) MM_PRESET="map-hifi" ;;
*) die "--platform must be ont or hifi" ;;
esac

# ---------------- Worker ----------------
assess_one() {
	local tag="$1" fasta="$2" dout="$3"
	mkdir -p "$dout"

	log 1 "[$tag] QC mapping primary-only"
	minimap2 -x "$MM_PRESET" -a -t "$THREADS" --secondary=no "$fasta" "$READS" \
		2>"$dout/qc.mapping.err" |
		samtools view -F 0x900 -u 2>"$dout/qc.samtools.view.err" |
		samtools sort -@ "$THREADS" -o "$dout/reads_to_asm.bam" - 2>"$dout/qc.samtools.sort.err"
	samtools index "$dout/reads_to_asm.bam" 2>"$dout/qc.samtools.index.err"

	log 1 "[$tag] Depth"
	samtools depth "$dout/reads_to_asm.bam" >"$dout/depth.tsv" 2>"$dout/depth.err"

	log 1 "[$tag] Depth summary"
	python3 "$PY_DEPTH" --depth "$dout/depth.tsv" --out "$dout/depth.summary.tsv" \
		2>"$dout/depth.summary.err"

	log 1 "[$tag] samtools stats -> QV and identity"
	samtools stats "$dout/reads_to_asm.bam" >"$dout/mapping.stats" 2>"$dout/mapping.stats.err"
	python3 "$PY_MAPSUM" --stats "$dout/mapping.stats" --out "$dout/mapping.summary.tsv" \
		2>"$dout/mapping.summary.err"

	log 1 "[$tag] Merge summaries"
	python3 "$PY_MERGE" --tables "$dout/depth.summary.tsv" "$dout/mapping.summary.tsv" \
		--out "$dout/summary.tsv" 2>"$dout/summary.err"

	if ((HAVE_R)) && [[ -s "$R_VIOL" ]]; then
		log 1 "[$tag] Violin plot"
		Rscript "$R_VIOL" "$dout/depth.tsv" "$dout/coverage_violin.png" "Depth ($tag)" \
			>"$dout/violin.log" 2>&1
	fi
}

# ---------------- Run ----------------
assess_one "A" "$FASTA_A" "$OUT/A"
assess_one "B" "$FASTA_B" "$OUT/B"

# ---------------- Compare ----------------
python3 "$PY_COMPARE" "$OUT/A/summary.tsv" "$OUT/B/summary.tsv" "$OUT/compare.tsv" \
	2>"$OUT/compare.err"

log 1 "[DONE] Compare table: $OUT/compare.tsv"
log 2 "A summary: $OUT/A/summary.tsv"
log 2 "B summary: $OUT/B/summary.tsv"
