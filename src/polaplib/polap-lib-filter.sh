################################################################################
# This file is part of polap.
#
# polap is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# polap is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# polap. If not, see <https://www.gnu.org/licenses/>.
################################################################################

################################################################################
# Convert numbers between different units.
################################################################################

################################################################################
# Ensure that the current script is sourced only once
source "${_POLAPLIB_DIR}/run-polap-function-include.sh"
_POLAP_INCLUDE_=$(_polap_include "${BASH_SOURCE[0]}")
set +u
if [[ -n "${!_POLAP_INCLUDE_}" ]]; then
	set -u
	return 0
fi
set -u
declare "$_POLAP_INCLUDE_=1"
#
################################################################################

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
	echo "[ERROR] This script must be sourced, not executed: use 'source $BASH_SOURCE'" >&2
	return 1 2>/dev/null || exit 1
fi
: "${_POLAP_DEBUG:=0}"
: "${_POLAP_RELEASE:=0}"

################################################################################
# TIPPo's GraphAligner-based filetring
################################################################################
_polap_lib_filter-reads-by-reference() {
	local outdir=""
	local fastq=""
	local platform="${_arg_data_type:-pacbio-hifi}"
	# local platform="pacbio-hifi"
	local reference=""

	# Argument parsing
	while [[ $# -gt 0 ]]; do
		case "$1" in
		-o)
			outdir="$2"
			shift 2
			;;
		-l)
			fastq="$2"
			shift 2
			;;
		--platform)
			platform="$2"
			shift 2
			;;
		--reference)
			reference="$2"
			shift 2
			;;
		-h | --help)
			echo "Usage: _polap_lib_filter-reads-by-reference -o OUTDIR -l FASTQ --reference GFA"
			return 0
			;;
		*)
			echo "[ERROR] Unknown option: $1"
			return 1
			;;
		esac
	done

	# Validate required arguments
	if [[ -z "$outdir" || -z "$fastq" || -z "$reference" ]]; then
		echo "[ERROR] Missing required arguments."
		echo "Usage: _polap_lib_filter-reads-by-reference -o OUTDIR -l FASTQ --reference GFA"
		return 1
	fi

	# Example command (replace with actual implementation)
	echo "[INFO] Filtering reads..."
	echo "  Output directory: $outdir"
	echo "  Input FASTQ:       $fastq"
	echo "  Reference GFA:     $reference"

	# use local variable

	local _arg_long_reads="$fastq"
	local _arg_reference="$reference"
	local _arg_outdir="$outdir"
	local _arg_data_type="$platform"

	_polap_filter-reads-by-reference

	# Do the actual work here (placeholder)
	echo "[INFO] Done."
}

# Usage:
#   _polap_filter-reads-by-reference \
#     --platform ont-q20 \
#     --id 0.86 --clip 250 --minlen 1000
#
# Platforms:
#   hifi     : PacBio HiFi (tight thresholds)
#   ont-q20  : ONT Q20+/SUP (moderate thresholds)
#   ont-raw  : Older/rough ONT (looser thresholds)

_polap_filter-reads-by-reference() {
	# === Inputs (positional already set in your wrapper) ===
	local MITO_READS="${_arg_long_reads}"  # reads.fastq(.gz)
	local CHLORO_GRAPH="${_arg_reference}" # chloroplast.gfa
	local OUTPUT_PREFIX="${_arg_outdir}/kmer"
	local THREADS="${_arg_threads}"

	# === Tunables (with sensible platform defaults) ===
	local PLATFORM="${_arg_data_type:-pacbio-hifi}"
	local ID_THRESH=""   # identity in [0,1]
	local CLIP_THRESH="" # max clipped bases tolerated
	local MIN_ALN_LEN="" # min aligned length (bp) to count as pt-like

	# Parse simple flags (optional)
	while [[ $# -gt 0 ]]; do
		case "$1" in
		--platform)
			PLATFORM="$2"
			shift 2
			;;
		--id)
			ID_THRESH="$2"
			shift 2
			;;
		--clip)
			CLIP_THRESH="$2"
			shift 2
			;;
		--minlen)
			MIN_ALN_LEN="$2"
			shift 2
			;;
		*) break ;;
		esac
	done

	# Defaults per platform
	case "$PLATFORM" in
	pacbio-hifi)
		: "${ID_THRESH:=0.95}"
		: "${CLIP_THRESH:=100}"
		: "${MIN_ALN_LEN:=800}" # HiFi: shorter good matches OK
		;;
	nano-hq)
		: "${ID_THRESH:=0.84}" # Q20/SUP ~98–99% raw ⇒ filtered identity ~0.86–0.9 is safe for pt removal
		: "${CLIP_THRESH:=250}"
		: "${MIN_ALN_LEN:=1000}"
		;;
	nano-raw)
		: "${ID_THRESH:=0.80}" # older ONT basecalls: be looser to avoid false negatives
		: "${CLIP_THRESH:=400}"
		: "${MIN_ALN_LEN:=2000}"
		;;
	*)
		echo "Unknown --platform '${PLATFORM}'. Use: hifi | ont-q20 | ont-raw" >&2
		return 2
		;;
	esac

	_polap_log0 "platform: ${PLATFORM}"
	_polap_log0 "  min. identity: ${ID_THRESH}"
	_polap_log0 "  min. clipping: ${CLIP_THRESH}"
	_polap_log0 "  min. alignment length: ${MIN_ALN_LEN}"

	# === Working paths ===
	mkdir -p "${OUTPUT_PREFIX}"
	local GAF="${OUTPUT_PREFIX}/aligned.gaf"
	local REJECT_IDS="${OUTPUT_PREFIX}/rejected.txt"
	local RETAIN_IDS="${OUTPUT_PREFIX}/retained.txt"
	local SUMMARY="${OUTPUT_PREFIX}/ref-summary.tsv"
	local FASTQ_REMOVED="${OUTPUT_PREFIX}/ref-chloroplast_like.fastq.gz"
	local FASTQ_FILTERED="${OUTPUT_PREFIX}/ref-filtered.fastq.gz"
	local RATIO="${OUTPUT_PREFIX}/chloroplast-ratio.txt"

	_polap_lib_conda-ensure_conda_env polap-graphaligner || exit 1

	echo "Aligning reads to graph with GraphAligner (${PLATFORM})..."
	GraphAligner \
		-g "$CHLORO_GRAPH" \
		-f "$MITO_READS" \
		-a "$GAF" \
		-t "$THREADS" \
		--precise-clipping 0.90 \
		-x vg

	conda deactivate

	# === Parse and filter (identity, clipping, min aligned length) ===
	# Assumes your script accepts: GAF reject.txt retain.txt summary.tsv id clip minlen
	echo "Parsing GAF with pysam helper..."
	python "${_POLAPLIB_DIR}/polap-py-parse-graphaligner-gaf.py" \
		"$GAF" "$REJECT_IDS" "$RETAIN_IDS" "$SUMMARY" "$ID_THRESH" "$CLIP_THRESH" "$MIN_ALN_LEN"

	# === Split FASTQ (gz outputs) ===
	_polap_log0 "Splitting FASTQ (gzipped outputs)..."
	# Removed (pt-like)
	seqkit grep -f "$REJECT_IDS" "$MITO_READS" | gzip -c >"$FASTQ_REMOVED"
	# Retained (mt-enriched)
	seqkit grep -f "$RETAIN_IDS" "$MITO_READS" | gzip -c >"$FASTQ_FILTERED"

	# === Compute pt-like ratio: removed / (removed + kept) ===
	local bases_removed bases_kept total
	bases_removed=$(seqkit stats -Ta "$FASTQ_REMOVED" | awk 'NR==2{print $5+0}')
	bases_kept=$(seqkit stats -Ta "$FASTQ_FILTERED" | awk 'NR==2{print $5+0}')
	total=$((bases_removed + bases_kept))
	if ((total > 0)); then
		echo "$(echo "scale=6; 100 * $bases_removed / $total" | bc)" >"$RATIO"
	else
		echo "0" >"$RATIO"
	fi

	_polap_log0 "Filtering complete (${PLATFORM})"
	_polap_log0 "  - Kept (mt-enriched): $FASTQ_FILTERED"
	_polap_log0 "  - Removed (pt-like):  $FASTQ_REMOVED"
	_polap_log0 "  - Summary:            $SUMMARY"
	_polap_log0 "  - pt-like % (bases):  $(cat "$RATIO")"
}

v1_polap_filter-reads-by-reference() {

	# === Inputs ===
	local MITO_READS="${_arg_long_reads}"     # e.g., mito_reads.fastq.gz
	local CHLORO_GRAPH="${_arg_reference}"    # e.g., chloroplast.gfa
	local OUTPUT_PREFIX="${_arg_outdir}/kmer" # e.g., filtered_output
	local THREADS="${_arg_threads}"
	local ID_THRESH="${5:-0.95}"
	local CLIP_THRESH="${6:-100}"

	# === Working paths ===
	# WORKDIR=$(mktemp -d)
	mkdir -p "${OUTPUT_PREFIX}"
	local SAM="${OUTPUT_PREFIX}/aligned.gaf"
	local REJECT_IDS="${OUTPUT_PREFIX}/rejected.txt"
	local RETAIN_IDS="${OUTPUT_PREFIX}/retained.txt"
	local SUMMARY="${OUTPUT_PREFIX}/ref-summary.tsv"
	local PLOT="${OUTPUT_PREFIX}/ref-summary.pdf"
	local FASTQ_REMOVED="${OUTPUT_PREFIX}/ref-chloroplast_like.fastq"
	local FASTQ_FILTERED="${OUTPUT_PREFIX}/ref-filtered.fastq"
	local RATIO="${OUTPUT_PREFIX}/chloroplast-ratio.txt"

	_polap_lib_conda-ensure_conda_env polap-graphaligner || exit 1

	# === Step 1: Align reads with GraphAligner
	echo "Aligning with GraphAligner..."
	GraphAligner \
		-g "$CHLORO_GRAPH" \
		-f "$MITO_READS" \
		-a "$SAM" \
		-t "$THREADS" \
		--precise-clipping 0.9 \
		-x vg

	conda deactivate

	# --seeds-mxm-length 19 \
	# --bandwidth 15 \

	# === Step 2: Parse and filter with Python
	echo "Parsing SAM with pysam..."
	python "${_POLAPLIB_DIR}"/polap-py-parse-graphaligner-gaf.py \
		"$SAM" "$REJECT_IDS" "$RETAIN_IDS" "$SUMMARY" "$ID_THRESH" "$CLIP_THRESH"

	# === Step 3: Extract FASTQ subsets
	echo "Splitting FASTQ..."
	seqkit grep -f "$REJECT_IDS" "$MITO_READS" >"$FASTQ_REMOVED"
	seqkit grep -f "$RETAIN_IDS" "$MITO_READS" >"$FASTQ_FILTERED"

	local bases_removed=$(seqkit stats -Ta "$FASTQ_REMOVED" | awk 'NR==2 {print $5}')
	local bases_filtered=$(seqkit stats -Ta "$FASTQ_FILTERED" | awk 'NR==2 {print $5}')
	echo $(echo "scale=0; 100 * $bases_removed/$bases_filtered" | bc) >"$RATIO"

	# seqkit grep -n -f "$REJECT_IDS" "$MITO_READS" | gzip -c >"$FASTQ_REMOVED"
	# seqkit grep -n -f "$RETAIN_IDS" "$MITO_READS" | gzip -c >"$FASTQ_FILTERED"

	# === Step 4: Plot with R
	# echo "📊 Generating plots..."
	# Rscript --vanialla "${_POLAPLIB_DIR}"/polap-r-plot-graphaligner-summary.R \
	# 	"$SUMMARY" "$PLOT"

	# === Summary
	_polap_log0 "Finished filtering mitochondrial reads:"
	_polap_log0 "  - Filtered reads: $FASTQ_FILTERED"
	_polap_log0 "  - Removed reads: $FASTQ_REMOVED"
	# _polap_log0 "- Summary table: $SUMMARY"
	# _polap_log0 "- Summary plot:  $PLOT"

}

_polap_lib_filter-reads-by-rmkc() {
	local outdir=""
	local fastq=""

	# Argument parsing
	while [[ $# -gt 0 ]]; do
		case "$1" in
		-o)
			outdir="$2"
			shift 2
			;;
		-l)
			fastq="$2"
			shift 2
			;;
		-h | --help)
			echo "Usage: _polap_lib_filter-reads-by-rmkc -o OUTDIR -l FASTQ"
			return 0
			;;
		*)
			echo "[ERROR] Unknown option: $1"
			return 1
			;;
		esac
	done

	# Validate required arguments
	if [[ -z "$outdir" || -z "$fastq" ]]; then
		echo "[ERROR] Missing required arguments."
		echo "Usage: _polap_lib_filter-reads-by-reference -o OUTDIR -l FASTQ --reference GFA"
		return 1
	fi

	# Example command (replace with actual implementation)
	echo "[INFO] Filtering reads..."
	echo "  Output directory: $outdir"
	echo "  Input FASTQ:       $fastq"

	local READS="${fastq}"

	# === Defaults ===
	MEM="100M"
	MINLEN=0
	VERBOSE=false
	LOGFILE=""
	K=21

	# Input variables
	# READS="${_arg_long_reads}"

	# PREFIX="${_arg_outdir}/kmer/rmkc"

	mkdir -p "${outdir}/kmer"
	PREFIX="${outdir}/kmer/rmkc"

	if [[ -z "$READS" || -z "$PREFIX" || -z "$K" ]]; then
		_polap_log0 "No read, no K, or no prefix"
		return
	fi

	# Output files
	JF="$PREFIX.jf"
	KMER_DUMP="$PREFIX.kmer_counts.txt"
	RMKC_TSV="$PREFIX.read_rmkc.tsv"
	MKC_OUT="$PREFIX.mkc.txt"
	CLEANED="$PREFIX.cleaned.fastq.gz"
	DISCARDED="$PREFIX.cleaned.discarded.fastq.gz"

	# === Step 1: Jellyfish counting
	_polap_log1 "[1/5] Counting k-mers with jellyfish..."
	jellyfish count -m "$K" -s "$MEM" -t 8 -C "$READS" -o "$JF"

	# === Step 2: Dump counts
	_polap_log1 "[2/5] Dumping k-mers..."
	jellyfish dump -c "$JF" >"$KMER_DUMP"

	# === Step 3: Compute rmkc
	_polap_log1 "[3/5] Computing per-read rmkc..."
	python "${_POLAPLIB_DIR}"/polap-py-compute-read-rmkc.py \
		--reads "$READS" \
		--kmer-counts "$KMER_DUMP" \
		--k "$K" \
		--output "$RMKC_TSV"

	# === Step 4: Compute mkc and thresholds
	_polap_log1 "[4/5] Computing mkc and thresholds..."
	Rscript -e "
x <- read.table('$RMKC_TSV', header=TRUE);
mkc <- median(x\$rmkc);
write(mkc, file='$MKC_OUT');
cat(sprintf('mkc = %.2f\\n', mkc));"

	local MKC=$(cat "$MKC_OUT")
	local LKC=$(awk -v mkc="$MKC" 'BEGIN { printf "%.2f", 0.3 * mkc }')
	local HKC=$(awk -v mkc="$MKC" 'BEGIN { printf "%.2f", 5.0 * mkc }')

	_polap_log0 "mkc = $MKC | lkc = $LKC | hkc = $HKC"

	# === Step 5: Filter nuclear outliers
	_polap_log1 "[5/5] Filtering out nuclear-like reads..."
	python "${_POLAPLIB_DIR}"/polap-py-remove-nuclear-outlier-reads.py \
		--reads "$READS" \
		--kmer-counts "$KMER_DUMP" \
		--k "$K" \
		--mkc "$MKC" \
		--minlen "$MINLEN" \
		--output "$PREFIX.cleaned"

	_polap_log1 "Final output:"
	_polap_log1 "  - Filtered FASTQ:    $CLEANED"
	_polap_log1 "  - Discarded FASTQ:   $DISCARDED"
	_polap_log1 "  - mkc:               $MKC_OUT"
	_polap_log1 "  - rmkc table:        $RMKC_TSV"

	# 	_polap_log1 "[6/6] Rendering HTML report..."
	# 	Rscript -e "rmarkdown::render(
	#   'rmkc_summary_report.Rmd',
	#   params = list(
	#     rmkc_file = '$RMKC_TSV',
	#     mkc_file = '$MKC_OUT',
	#     prefix = '$PREFIX'
	#   ),
	#   output_file = '${PREFIX}.rmkc.summary.html'
	# )"
	# 	_polap_log1 "Report created: ${PREFIX}.rmkc.summary.html"

}
