#!/usr/bin/env bash
# Version: v0.2.0
# Name: polap-bash-compare-assemblies.sh
# Purpose: Fairly compare TWO assemblies (A vs B) using the SAME input data:
#          - ONT/HiFi long reads (--reads-ont)
#          - and/or short reads (--reads-sr1/--reads-sr2)
#          For SR, also run Merqury per FASTA and merge QV into summaries.
#
# Outputs:
#   out_cmp/
#     ONT/ A/... B/... compare_ONT.tsv          (if --reads-ont provided)
#     SR/  A/... B/... compare_SR.tsv           (if SR provided)
#     SR/  merqury/ A_qv.tsv B_qv.tsv compare_MERQURY.tsv  (if Merqury available)
#     compare_merged.tsv   (single merged summary across available modalities)
#     report.md
#
# Dependencies: minimap2, samtools, python3
# Optional: Rscript + ggplot2 (violins), meryl + merqury.sh (for Merqury QV)
# Helpers expected in ./scripts/ :
#   polap-py-depth-summary.py
#   polap-py-mapping-summary-from-stats.py
#   polap-py-merge-kv-tables.py
#   polap-py-compare-polish-metrics.py
#   polap-py-extract-merqury-qv.py
#   polap-py-stack-compare-tables.py
#   polap-r-coverage-violin.R

set -Eeuo pipefail

# ---------------- Defaults ----------------
FASTA_A=""
FASTA_B=""
OUT=""
READS_ONT=""
READS_SR1=""
READS_SR2=""
ONT_PRESET="map-ont" # or map-hifi for HiFi
THREADS=32
VERB=1

# ---------------- Logging ----------------
log() {
	local lvl="$1"
	shift
	((VERB >= lvl)) && printf '[compare] %s\n' "$*" >&2 || true
}
die() {
	printf '[ERROR] %s\n' "$*" >&2
	exit 2
}

# ---------------- CLI ----------------
print_help() {
	cat <<EOF
Usage: $0 --fastaA A.fa --fastaB B.fa --outdir out/ [options]
  --reads-long FILE        ONT/HiFi long reads (FASTQ[.gz])
  --reads-sr1 FILE         Short-read R1 (FASTQ[.gz])
  --reads-sr2 FILE         Short-read R2 (FASTQ[.gz])
  --ont-preset PRESET      map-ont (default) or map-hifi
  --threads N              Default 32
  -v                       Increase verbosity (repeatable)
  -h|--help                Show help
EOF
}

while (("$#")); do
	case "$1" in
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
	--reads-ont | --reads-long)
		READS_ONT="${2:?}"
		shift 2
		;;
	--reads-sr1)
		READS_SR1="${2:?}"
		shift 2
		;;
	--reads-sr2)
		READS_SR2="${2:?}"
		shift 2
		;;
	--ont-preset)
		ONT_PRESET="${2:?}"
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

[[ -n "$FASTA_A" && -n "$FASTA_B" && -n "$OUT" ]] || {
	print_help
	die "Required: --fastaA --fastaB --outdir"
}
[[ -s "$FASTA_A" ]] || die "FASTA A not found: $FASTA_A"
[[ -s "$FASTA_B" ]] || die "FASTA B not found: $FASTA_B"
if [[ -z "$READS_ONT" && -z "$READS_SR1" ]]; then
	die "Provide at least one modality: --reads-ont or --reads-sr1/--reads-sr2"
fi
if [[ -n "$READS_SR1" && -z "$READS_SR2" ]]; then
	die "--reads-sr2 is required when --reads-sr1 is given"
fi

for t in minimap2 samtools python3; do
	command -v "$t" >/dev/null 2>&1 || die "Missing dependency: $t"
done
HAVE_R=0
command -v Rscript >/dev/null 2>&1 && HAVE_R=1
HAVE_MERYL=0
command -v meryl >/dev/null 2>&1 && HAVE_MERYL=1
HAVE_MERQURY=0
command -v merqury.sh >/dev/null 2>&1 && HAVE_MERQURY=1

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
S="${ROOT_DIR}/scripts"

PY_DEPTH="${S}/polap-py-depth-summary.py"
PY_MAPSUM="${S}/polap-py-mapping-summary-from-stats.py"
PY_MERGE="${S}/polap-py-merge-kv-tables.py"
PY_COMPARE="${S}/polap-py-compare-polish-metrics.py"
PY_MERQ="${S}/polap-py-extract-merqury-qv.py"
PY_STACK="${S}/polap-py-stack-compare-tables.py"
R_VIOL="${S}/polap-r-coverage-violin.R"

for f in "$PY_DEPTH" "$PY_MAPSUM" "$PY_MERGE" "$PY_COMPARE" "$PY_MERQ" "$PY_STACK"; do
	[[ -s "$f" ]] || die "Missing helper: $f"
done
if ((HAVE_R)) && [[ ! -s "$R_VIOL" ]]; then
	log 1 "Rscript present but ${R_VIOL} is missing; violin plots will be skipped"
	HAVE_R=0
fi

mkdir -p "$OUT"

# ---------------- Worker: assess one assembly with one modality ----------------
# Args: tag(A|B) mode(ONT|SR) fasta outdir
assess_one() {
	local tag="$1" mode="$2" fasta="$3" dout="$4"
	mkdir -p "$dout"

	if [[ "$mode" == "ONT" ]]; then
		log 1 "[$tag][$mode] Mapping (primary-only) with preset ${ONT_PRESET}"
		minimap2 -x "$ONT_PRESET" -a -t "$THREADS" --secondary=no "$fasta" "$READS_ONT" \
			2>"$dout/map.err" |
			samtools view -F 0x900 -u 2>"$dout/view.err" |
			samtools sort -@ "$THREADS" -o "$dout/reads_to_asm.bam" - 2>"$dout/sort.err"
	else
		log 1 "[$tag][$mode] Mapping paired SR (primary-only)"
		minimap2 -x sr -a -t "$THREADS" "$fasta" \
			<([[ "$READS_SR1" =~ \.gz$ ]] && gzip -cd "$READS_SR1" || cat "$READS_SR1") \
			<([[ "$READS_SR2" =~ \.gz$ ]] && gzip -cd "$READS_SR2" || cat "$READS_SR2") \
			2>"$dout/map.err" |
			samtools view -F 0x900 -u 2>"$dout/view.err" |
			samtools sort -@ "$THREADS" -o "$dout/reads_to_asm.bam" - 2>"$dout/sort.err"
	fi
	samtools index "$dout/reads_to_asm.bam" 2>"$dout/index.err"

	log 1 "[$tag][$mode] Depth and stats"
	samtools depth "$dout/reads_to_asm.bam" >"$dout/depth.tsv" 2>"$dout/depth.err"
	samtools stats "$dout/reads_to_asm.bam" >"$dout/mapping.stats" 2>"$dout/mapping.stats.err"

	python3 "$PY_DEPTH" --depth "$dout/depth.tsv" --out "$dout/depth.summary.tsv" 2>"$dout/depth.summary.err"
	python3 "$PY_MAPSUM" --stats "$dout/mapping.stats" --out "$dout/mapping.summary.tsv" 2>"$dout/mapping.summary.err"
	python3 "$PY_MERGE" --tables "$dout/depth.summary.tsv" "$dout/mapping.summary.tsv" \
		--out "$dout/summary.tsv" 2>"$dout/summary.err"

	if ((HAVE_R)); then
		Rscript "$R_VIOL" "$dout/depth.tsv" "$dout/coverage_${tag}_${mode}.png" "Depth (${tag}, ${mode})" \
			>"$dout/violin.log" 2>&1 || true
	fi
}

REPORT="${OUT}/report.md"
: >"$REPORT"
echo "# POLAP comparison report" >>"$REPORT"
echo "" >>"$REPORT"
echo "Assemblies:" >>"$REPORT"
echo "- A: \`$FASTA_A\`" >>"$REPORT"
echo "- B: \`$FASTA_B\`" >>"$REPORT"
echo "" >>"$REPORT"
echo "Threads: $THREADS" >>"$REPORT"
echo "" >>"$REPORT"

# ---------------- ONT modality ----------------
HAVE_ONT=0
if [[ -n "$READS_ONT" ]]; then
	HAVE_ONT=1
	log 1 "[ONT] Comparing with ONT/HiFi reads"
	mkdir -p "$OUT/ONT/A" "$OUT/ONT/B"

	assess_one "A" "ONT" "$FASTA_A" "$OUT/ONT/A"
	assess_one "B" "ONT" "$FASTA_B" "$OUT/ONT/B"

	python3 "$PY_COMPARE" "$OUT/ONT/A/summary.tsv" "$OUT/ONT/B/summary.tsv" "$OUT/ONT/compare_ONT.tsv" \
		2>"$OUT/ONT/compare.err"

	echo "## ONT/HiFi evaluation (preset: $ONT_PRESET)" >>"$REPORT"
	echo "" >>"$REPORT"
	echo '```' >>"$REPORT"
	column -t -s$'\t' "$OUT/ONT/compare_ONT.tsv" || cat "$OUT/ONT/compare_ONT.tsv"
	echo '```' >>"$REPORT"
	echo "" >>"$REPORT"
fi

# ---------------- SR modality ----------------
HAVE_SR=0
if [[ -n "$READS_SR1" ]]; then
	HAVE_SR=1
	log 1 "[SR] Comparing with short reads"
	mkdir -p "$OUT/SR/A" "$OUT/SR/B"

	assess_one "A" "SR" "$FASTA_A" "$OUT/SR/A"
	assess_one "B" "SR" "$FASTA_B" "$OUT/SR/B"

	python3 "$PY_COMPARE" "$OUT/SR/A/summary.tsv" "$OUT/SR/B/summary.tsv" "$OUT/SR/compare_SR.tsv" \
		2>"$OUT/SR/compare.err"

	echo "## Short-read evaluation" >>"$REPORT"
	echo "" >>"$REPORT"
	echo '```' >>"$REPORT"
	column -t -s$'\t' "$OUT/SR/compare_SR.tsv" || cat "$OUT/SR/compare_SR.tsv"
	echo '```' >>"$REPORT"
	echo "" >>"$REPORT"

	# -------- Merqury block (if tools available) --------
	if ((HAVE_MERYL)) && ((HAVE_MERQURY)); then
		log 1 "[Merqury] Building k-mer db from SR and computing QV per assembly"
		mkdir -p "$OUT/SR/merqury"
		# SR trusted k-mers
		if [[ ! -s "$OUT/SR/merqury/reads.meryl" ]]; then
			meryl count k=21 "$READS_SR1" "$READS_SR2" output "$OUT/SR/merqury/reads.meryl" \
				>"$OUT/SR/merqury/reads.meryl.out" 2>"$OUT/SR/merqury/reads.meryl.err"
		fi
		# Assembly k-mers
		for TAG in A B; do
			FASTA_VAR="$FASTA_A"
			[[ "$TAG" == "B" ]] && FASTA_VAR="$FASTA_B"
			if [[ ! -s "$OUT/SR/merqury/${TAG}.meryl" ]]; then
				meryl count k=21 "$FASTA_VAR" output "$OUT/SR/merqury/${TAG}.meryl" \
					>"$OUT/SR/merqury/${TAG}.meryl.out" 2>"$OUT/SR/merqury/${TAG}.meryl.err"
			fi
			merqury.sh "$OUT/SR/merqury/reads.meryl" "$OUT/SR/merqury/${TAG}.meryl" "$OUT/SR/merqury/${TAG}" \
				>"$OUT/SR/merqury/${TAG}.merqury.out" 2>"$OUT/SR/merqury/${TAG}.merqury.err" || true
			# Extract QV from merqury output into kv table
			python3 "$PY_MERQ" --merq-out "$OUT/SR/merqury/${TAG}" --out "$OUT/SR/merqury/${TAG}_qv.tsv" \
				2>"$OUT/SR/merqury/${TAG}_qv.err"
			# Merge Merqury QV into assembly summary.tsv (so compare_SR also reflects it if needed)
			python3 "$PY_MERGE" --tables "$OUT/SR/${TAG}/summary.tsv" "$OUT/SR/merqury/${TAG}_qv.tsv" \
				--out "$OUT/SR/${TAG}/summary.tsv" 2>"$OUT/SR/${TAG}/summary.merge.err"
		done
		# Produce a tiny compare table just for Merqury_QV
		python3 "$PY_COMPARE" "$OUT/SR/merqury/A_qv.tsv" "$OUT/SR/merqury/B_qv.tsv" "$OUT/SR/compare_MERQURY.tsv" \
			2>"$OUT/SR/compare_MERQURY.err"

		echo "### Merqury QV (reference-free, k=21)" >>"$REPORT"
		echo "" >>"$REPORT"
		echo '```' >>"$REPORT"
		column -t -s$'\t' "$OUT/SR/compare_MERQURY.tsv" || cat "$OUT/SR/compare_MERQURY.tsv"
		echo '```' >>"$REPORT"
		echo "" >>"$REPORT"
	else
		log 1 "[Merqury] Skipped (meryl or merqury.sh not found)"
	fi
fi

# ---------------- Single merged summary ----------------
# Stack available compare tables (label, file) into one compare_merged.tsv
STACK_ARGS=()
if ((HAVE_ONT)); then STACK_ARGS+=("ONT" "$OUT/ONT/compare_ONT.tsv"); fi
if ((HAVE_SR)); then STACK_ARGS+=("SR" "$OUT/SR/compare_SR.tsv"); fi
if ((HAVE_SR)) && ((HAVE_MERYL)) && ((HAVE_MERQURY)); then
	STACK_ARGS+=("MERQURY" "$OUT/SR/compare_MERQURY.tsv")
fi

if ((${#STACK_ARGS[@]} > 0)); then
	python3 "$PY_STACK" "${STACK_ARGS[@]}" --out "$OUT/compare_merged.tsv" \
		2>"$OUT/compare_merged.err"
	echo "## Merged summary" >>"$REPORT"
	echo "" >>"$REPORT"
	echo '```' >>"$REPORT"
	column -t -s$'\t' "$OUT/compare_merged.tsv" || cat "$OUT/compare_merged.tsv"
	echo '```' >>"$REPORT"
	echo "" >>"$REPORT"
fi

log 1 "[DONE] See report: $REPORT"
