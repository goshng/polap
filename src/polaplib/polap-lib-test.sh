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

_polap_lib_test-polap-common-variable() {
	local outdir=""
	local fastq=""
	local reference=""
	local _arg_inum=""
	local _arg_jnum=""

	# Argument parsing
	while [[ $# -gt 0 ]]; do
		case "$1" in
		-o)
			outdir="$2"
			shift 2
			;;
		-i)
			_arg_inum="$2"
			shift 2
			;;
		-j)
			_arg_jnum="$2"
			shift 2
			;;
		-l)
			fastq="$2"
			shift 2
			;;
		--reference)
			reference="$2"
			shift 2
			;;
		-h | --help)
			echo "Usage: _polap_lib_test-reads-by-reference -o OUTDIR -l FASTQ --reference GFA"
			return 0
			;;
		*)
			echo "[ERROR] Unknown option: $1"
			return 1
			;;
		esac
	done

	local _arg_outdir="$outdir"
	local _arg_long_reads="$fastq"
	local _arg_reference="$reference"

	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	_polap_log0 "$FUNCNAME: _polap_var_oga_mt_fasta: $_polap_var_oga_mt_fasta"
}
################################################################################
# TIPPo's GraphAligner-based filetring
################################################################################
_polap_lib_test-variable-scope() {
	local outdir=""
	local fastq=""
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
		--reference)
			reference="$2"
			shift 2
			;;
		-h | --help)
			echo "Usage: _polap_lib_test-reads-by-reference -o OUTDIR -l FASTQ --reference GFA"
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
		echo "Usage: _polap_lib_test-reads-by-reference -o OUTDIR -l FASTQ --reference GFA"
		return 1
	fi

	# Example command (replace with actual implementation)
	echo "[INFO] Filtering reads..."
	echo "  Output directory: $outdir"
	echo "  Input FASTQ:       $fastq"
	echo "  Reference GFA:     $reference"

	local _arg_long_reads="$fastq"
	local _arg_reference="$reference"
	local _arg_outdir="outdir"

	_polap_test-level2

	# Do the actual work here (placeholder)
	echo "[INFO] Done."
}

_polap_test-level2() {

	# === Inputs ===
	MITO_READS="${_arg_long_reads}"     # e.g., mito_reads.fastq.gz
	CHLORO_GRAPH="${_arg_reference}"    # e.g., chloroplast.gfa
	OUTPUT_PREFIX="${_arg_outdir}/kmer" # e.g., filtered_output
	THREADS="${_arg_threads}"
	ID_THRESH="${5:-0.95}"
	CLIP_THRESH="${6:-100}"

	# === Working paths ===
	# WORKDIR=$(mktemp -d)
	mkdir -p "${OUTPUT_PREFIX}"
	SAM="${OUTPUT_PREFIX}/aligned.gaf"
	REJECT_IDS="${OUTPUT_PREFIX}/rejected.txt"
	RETAIN_IDS="${OUTPUT_PREFIX}/retained.txt"
	SUMMARY="${OUTPUT_PREFIX}/ref-summary.tsv"
	PLOT="${OUTPUT_PREFIX}/ref-summary.pdf"
	FASTQ_REMOVED="${OUTPUT_PREFIX}/ref-chloroplast_like.fastq"
	FASTQ_FILTERED="${OUTPUT_PREFIX}/ref-filtered.fastq"

	_polap_log0 GraphAligner \
		-g "$CHLORO_GRAPH" \
		-f "$MITO_READS" \
		-a "$SAM" \
		-t "$THREADS" \
		--precise-clipping 0.9 \
		-x vg

	return

	_polap_lib_conda-ensure_conda_env polap-graphaligner || exit 1

	# === Step 1: Align reads with GraphAligner
	echo "📌 Aligning with GraphAligner..."
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
	echo "🧠 Parsing SAM with pysam..."
	python "${_POLAPLIB_DIR}"/polap-py-parse-graphaligner-gaf.py \
		"$SAM" "$REJECT_IDS" "$RETAIN_IDS" "$SUMMARY" "$ID_THRESH" "$CLIP_THRESH"

	# === Step 3: Extract FASTQ subsets
	echo "✂️ Splitting FASTQ..."
	seqkit grep -f "$REJECT_IDS" "$MITO_READS" >"$FASTQ_REMOVED"
	seqkit grep -f "$RETAIN_IDS" "$MITO_READS" >"$FASTQ_FILTERED"
	# seqkit grep -n -f "$REJECT_IDS" "$MITO_READS" | gzip -c >"$FASTQ_REMOVED"
	# seqkit grep -n -f "$RETAIN_IDS" "$MITO_READS" | gzip -c >"$FASTQ_FILTERED"

	# === Step 4: Plot with R
	# echo "📊 Generating plots..."
	# Rscript --vanialla "${_POLAPLIB_DIR}"/polap-r-plot-graphaligner-summary.R \
	# 	"$SUMMARY" "$PLOT"

	# === Summary
	_polap_log0 "✅ Finished filtering mitochondrial reads:"
	_polap_log0 "- Filtered reads: $FASTQ_FILTERED"
	_polap_log0 "- Removed reads: $FASTQ_REMOVED"
	# _polap_log0 "- Summary table: $SUMMARY"
	# _polap_log0 "- Summary plot:  $PLOT"

}
