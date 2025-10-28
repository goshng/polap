#!/usr/bin/env bash
# Version: v0.2.0
# Name: polap-bash-polish-fmlrc2.sh
# Purpose: Polish an assembly using short-read data with fmlrc2, then run QC:
#          map SR to BEFORE and AFTER assemblies, compute depth + mapping stats,
#          make violin plots, and produce a side-by-side comparison report + figure.
#
# Usage:
#   polap-bash-polish-fmlrc2.sh \
#     --fasta draft.fa \
#     --sr1 R1.fastq.gz \
#     --sr2 R2.fastq.gz \
#     --outdir out_fmlrc2 \
#     [--threads 32] [--prefix sample] [-v ...]
#
# Outputs (in outdir/):
#   comp_msbwt.npy
#   segments.fmlrc2.fa
#   <prefix>.polished.fmlrc2.fa
#   qc/BEFORE/... , qc/AFTER/...
#   qc/compare.tsv
#   qc/report.md
#   qc/coverage_before.png, qc/coverage_after.png, qc/coverage_compare.png (if ImageMagick present)
#
# Dependencies:
#   ropebwt2, fmlrc2-convert, fmlrc2, minimap2, samtools, python3
# Optional:
#   Rscript + ggplot2 (for violin plots)
#   ImageMagick "montage" (to make side-by-side figure)
#
# Helper scripts expected in ./scripts/ :
#   polap-py-mapping-summary-from-stats.py
#   polap-py-depth-summary.py
#   polap-py-merge-kv-tables.py
#   polap-py-compare-polish-metrics.py
#   polap-r-coverage-violin.R

set -Eeuo pipefail

# ---------------- defaults ----------------
FASTA=""
SR1=""
SR2=""
OUT=""
THREADS=32
PREFIX="polap"
VERB=1

# ---------------- logging ----------------
log() {
	local lvl="$1"
	shift
	((VERB >= lvl)) && printf '[fmlrc2] %s\n' "$*" >&2 || true
}
die() {
	printf '[ERROR] %s\n' "$*" >&2
	exit 2
}

# ---------------- cli ----------------
print_help() {
	cat <<EOF
Usage: $0 --fasta draft.fa --sr1 R1.fq.gz --sr2 R2.fq.gz --outdir out/ [options]
  --threads N     Default 32
  --prefix NAME   Default polap
  -v              Increase verbosity (repeatable)
  -h|--help       Show help
EOF
}

while (("$#")); do
	case "$1" in
	--fasta)
		FASTA="${2:?}"
		shift 2
		;;
	--sr1)
		SR1="${2:?}"
		shift 2
		;;
	--sr2)
		SR2="${2:?}"
		shift 2
		;;
	--outdir)
		OUT="${2:?}"
		shift 2
		;;
	--threads)
		THREADS="${2:?}"
		shift 2
		;;
	--prefix)
		PREFIX="${2:?}"
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

[[ -n "$FASTA" && -n "$SR1" && -n "$SR2" && -n "$OUT" ]] || {
	print_help
	die "Required: --fasta --sr1 --sr2 --outdir"
}
[[ -s "$FASTA" ]] || die "FASTA not found: $FASTA"
[[ -s "$SR1" ]] || die "SR1 not found: $SR1"
[[ -s "$SR2" ]] || die "SR2 not found: $SR2"

for t in ropebwt2 fmlrc2-convert fmlrc2 minimap2 samtools python3; do
	command -v "$t" >/dev/null 2>&1 || die "Missing dependency: $t"
done
HAVE_R=0
command -v Rscript >/dev/null 2>&1 && HAVE_R=1
HAVE_MONTAGE=0
command -v montage >/dev/null 2>&1 && HAVE_MONTAGE=1

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
S="${ROOT_DIR}/scripts"

# helpers
PY_MAPSUM="${S}/polap-py-mapping-summary-from-stats.py"
PY_DEPTH="${S}/polap-py-depth-summary.py"
PY_MERGE="${S}/polap-py-merge-kv-tables.py"
PY_COMPARE="${S}/polap-py-compare-polish-metrics.py"
R_VIOL="${S}/polap-r-coverage-violin.R"

for f in "$PY_MAPSUM" "$PY_DEPTH" "$PY_MERGE" "$PY_COMPARE"; do
	[[ -s "$f" ]] || die "Missing helper: $f"
done
if ((HAVE_R)) && [[ ! -s "$R_VIOL" ]]; then
	log 1 "Rscript is present but ${R_VIOL} is missing; violin plots will be skipped"
	HAVE_R=0
fi

mkdir -p "$OUT" "$OUT/qc"

# ---------------- build MSBWT from SR ----------------
log 1 "Building MSBWT from short reads (threads=$THREADS)"
# Extract read sequences (line 2 of FASTQ) from R1 and R2 and stream into ropebwt2 then fmlrc2-convert.
{
	if [[ "$SR1" =~ \.gz$ ]]; then gzip -cd "$SR1"; else cat "$SR1"; fi
	echo
	if [[ "$SR2" =~ \.gz$ ]]; then gzip -cd "$SR2"; else cat "$SR2"; fi
} |
	awk 'NR%4==2' |
	ropebwt2 -LR -t "$THREADS" 1>"$OUT/ropebwt2.out" 2>"$OUT/ropebwt2.err" |
	fmlrc2-convert "$OUT/comp_msbwt.npy" 1>"$OUT/fmlrc2-convert.out" 2>"$OUT/fmlrc2-convert.err"

[[ -s "$OUT/comp_msbwt.npy" ]] || die "MSBWT build failed: $OUT/comp_msbwt.npy missing or empty"

# ---------------- run fmlrc2 polish ----------------
POL_RAW_OUT="$OUT/segments.fmlrc2.fa"
FINAL="$OUT/${PREFIX}.polished.fmlrc2.fa"
log 1 "Running fmlrc2 polish"
fmlrc2 -t "$THREADS" "$OUT/comp_msbwt.npy" "$FASTA" "$POL_RAW_OUT" \
	1>"$OUT/fmlrc2.out" 2>"$OUT/fmlrc2.err"

[[ -s "$POL_RAW_OUT" ]] || die "fmlrc2 polish failed: $POL_RAW_OUT missing or empty"
cp -f "$POL_RAW_OUT" "$FINAL"
log 1 "Polished FASTA: $FINAL"

# ---------------- QC functions ----------------
qc_one() {
	local label="$1" fasta="$2" dir="$3"
	mkdir -p "$dir"
	log 1 "[QC:$label] Mapping paired SR (primary-only)"
	# minimap2 short-read preset; primary-only by filtering flag (0x900)
	# stderr to .err files to keep console quiet
	minimap2 -x sr -a -t "$THREADS" "$fasta" \
		<(if [[ "$SR1" =~ \.gz$ ]]; then gzip -cd "$SR1"; else cat "$SR1"; fi) \
		<(if [[ "$SR2" =~ \.gz$ ]]; then gzip -cd "$SR2"; else cat "$SR2"; fi) \
		2>"$dir/map.err" |
		samtools view -F 0x900 -u 2>"$dir/view.err" |
		samtools sort -@ "$THREADS" -o "$dir/reads_to_asm.bam" - 2>"$dir/sort.err"
	samtools index "$dir/reads_to_asm.bam" 2>"$dir/index.err"

	log 1 "[QC:$label] Depth and stats"
	samtools depth "$dir/reads_to_asm.bam" >"$dir/depth.tsv" 2>"$dir/depth.err"
	samtools stats "$dir/reads_to_asm.bam" >"$dir/mapping.stats" 2>"$dir/mapping.stats.err"

	python3 "$PY_DEPTH" --depth "$dir/depth.tsv" --out "$dir/depth.summary.tsv" 2>"$dir/depth.summary.err"
	python3 "$PY_MAPSUM" --stats "$dir/mapping.stats" --out "$dir/mapping.summary.tsv" 2>"$dir/mapping.summary.err"
	python3 "$PY_MERGE" --tables "$dir/depth.summary.tsv" "$dir/mapping.summary.tsv" \
		--out "$dir/summary.tsv" 2>"$dir/summary.err"

	if ((HAVE_R)); then
		Rscript "$R_VIOL" "$dir/depth.tsv" "$dir/coverage_${label}.png" "Depth ($label)" \
			>"$dir/violin.log" 2>&1 || true
	fi
}

# ---------------- QC BEFORE/AFTER ----------------
log 1 "QC BEFORE polishing"
qc_one "before" "$FASTA" "$OUT/qc/BEFORE"

log 1 "QC AFTER polishing"
qc_one "after" "$FINAL" "$OUT/qc/AFTER"

# ---------------- Compare and make report ----------------
python3 "$PY_COMPARE" "$OUT/qc/BEFORE/summary.tsv" "$OUT/qc/AFTER/summary.tsv" "$OUT/qc/compare.tsv" \
	2>"$OUT/qc/compare.err"

# Optional combined figure if ImageMagick is available
if ((HAVE_R)); then
	if ((HAVE_MONTAGE)) && [[ -s "$OUT/qc/BEFORE/coverage_before.png" && -s "$OUT/qc/AFTER/coverage_after.png" ]]; then
		montage "$OUT/qc/BEFORE/coverage_before.png" "$OUT/qc/AFTER/coverage_after.png" \
			-tile 2x1 -geometry +10+10 "$OUT/qc/coverage_compare.png" 2>"$OUT/qc/montage.err" || true
	fi
fi

# Write a simple Markdown report
{
	echo "# POLAP fmlrc2 polishing report"
	echo
	echo "## Inputs"
	echo "- FASTA (before): \`$FASTA\`"
	echo "- FASTA (after):  \`$FINAL\`"
	echo "- SR1: \`$SR1\`, SR2: \`$SR2\`"
	echo "- Threads: $THREADS"
	echo
	echo "## Summary metrics (before vs after)"
	echo
	echo '```'
	column -t -s$'\t' "$OUT/qc/compare.tsv" || cat "$OUT/qc/compare.tsv"
	echo '```'
	echo
	if [[ -s "$OUT/qc/coverage_compare.png" ]]; then
		echo "## Coverage violin (before vs after)"
		echo "![coverage compare](coverage_compare.png)"
	else
		if [[ -s "$OUT/qc/BEFORE/coverage_before.png" ]]; then
			echo "## Coverage violin (before)"
			echo "![before](BEFORE/coverage_before.png)"
		fi
		if [[ -s "$OUT/qc/AFTER/coverage_after.png" ]]; then
			echo "## Coverage violin (after)"
			echo "![after](AFTER/coverage_after.png)"
		fi
	fi
	echo
	echo "## Files"
	echo "- Comparison table: \`qc/compare.tsv\`"
	echo "- Before summary:   \`qc/BEFORE/summary.tsv\`"
	echo "- After summary:    \`qc/AFTER/summary.tsv\`"
	echo "- Mapping stats:    \`qc/BEFORE/mapping.stats\`, \`qc/AFTER/mapping.stats\`"
	echo "- Depth TSV:        \`qc/BEFORE/depth.tsv\`, \`qc/AFTER/depth.tsv\`"
} >"$OUT/qc/report.md"

log 1 "DONE"
log 1 "Polished FASTA: $FINAL"
log 1 "QC compare:     $OUT/qc/compare.tsv"
log 1 "Report:         $OUT/qc/report.md"
[[ -s "$OUT/qc/coverage_compare.png" ]] && log 1 "Figure:         $OUT/qc/coverage_compare.png"
