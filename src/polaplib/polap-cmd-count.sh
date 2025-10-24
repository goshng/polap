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

function _run_polap_count {
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

	# Print help message if requested
	help_message=$(
		cat <<'EOF'
Name:
  count - k-mers

Synopsis:
  polap count [ -acdfhklLnNrtvV19 ] [ name ... ]

Description:
  polap count uses JellyFish to count k-mers.

Options:
  -l FASTQ
    Long-read data

  -x-min FLOAT
    Left boundary of the X-axis in the k-mers distribution

  -x-max FLOAT
    Right boundary of the X-axis in the k-mers distribution

  --pacbio-hifi, --pacbio-raw, --nano-raw, --nano-hq
    Sequence data type

Examples:
  Count k-mers:
    polap count --pacbio-hifi -l o/reads.fq

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (1998–2018)

Author:
  Sang Chul Choi
EOF
	)

	# Display help message
	if [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]]; then
		local manfile=$(_polap_lib_man-convert_help_message "$help_message" "${_arg_menu[0]}")
		man "$manfile" >&3
		rm -f "$manfile"
		return
	fi

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		_polap_count_view_kmer

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
		exit $EXIT_SUCCESS
	fi

	if [[ "${_arg_menu[1]}" == "filter" ]]; then
		_polap_count_filter_kmer
	fi

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

_polap_count_kmer_jellyfish() {

	kmerdir="${_arg_outdir}"/kmer

	# === Defaults ===
	K="${_arg_knum}"
	K_HIFI=21
	K_ERR=17
	K_ONT_RAW=17
	JFSIZE=1G
	THREADS="${_arg_threads}"
	PREFIX=""
	PEAK_OVERRIDE=""
	if [[ "${_arg_coverage_is}" == "on" ]]; then
		PEAK_OVERRIDE="${_arg_coverage}"
	fi
	X_MIN="${_arg_x_min}"
	X_MAX="${_arg_x_max}"

	# === Parse args ===
	READS="${_arg_long_reads}"
	TECH="${_arg_data_type}"

	# === Auto-config based on tech ===
	case "$TECH" in
	pacbio-hifi)
		K=$K_HIFI
		;;
	pacbio-raw | nano-hq)
		K=$K_ERR
		JFSIZE=500M
		;;
	nano-raw)
		K=$K_ONT_RAW
		JFSIZE=500M
		;;
	*)
		K="${_arg_knum}"
		JFSIZE=500M
		;;
	esac

	if [[ "${_arg_knum}" != "1" ]]; then
		K="${_arg_knum}"
	fi

	[[ -z "$PREFIX" ]] && PREFIX="k$K"

	# === Output folders ===
	mkdir -p "${kmerdir}"
	KMER_TABLE="${kmerdir}/k${K}_all_reads_kmers.txt"
	HISTO="${kmerdir}/k${K}_reads.histo"
	PDF="${kmerdir}/k${K}_kmer_histogram.pdf"
	PEAKS_TSV="${kmerdir}/k${K}_peak_thresholds.tsv"
	FILTERED="${kmerdir}/k${K}_filtered.fastq.gz"
	DISCARDED="${kmerdir}/k${K}_discarded.fastq.gz"
	LABELS_TSV="${kmerdir}/read_labels.tsv"

	if [[ ! -s "$HISTO" ]]; then
		_polap_log1 "Count k-mers with Jellyfish (k=$K) ..."
		jellyfish count -m $K -s $JFSIZE -t $THREADS -C "$READS" -o "${kmerdir}/reads.jf"
		jellyfish dump -c "${kmerdir}/reads.jf" >"$KMER_TABLE"
		jellyfish histo "${kmerdir}/reads.jf" >"$HISTO"
	fi

	_polap_log1 "Detect peaks in the k-mers distribution..."
	Rscript "${_POLAPLIB_DIR}"/polap-r-plot-kmer-histogram-peaks-valleys.R "$HISTO" "$PDF" "$X_MIN" "$X_MAX" >"$PEAKS_TSV" 2>/dev/null

	# === Use override or compute ===
	if [[ -n "$PEAK_OVERRIDE" ]]; then
		_polap_log0 "⚙️ Skipping Jellyfish and using manual depth threshold: $PEAK_OVERRIDE"
		DEPTH_THRESHOLD="$PEAK_OVERRIDE"
	else
		ORGANELLE_VALLEY=$(awk '$1 == 4 {print $5}' "$PEAKS_TSV")
		PLASTID_VALLEY=$(awk '$1 == 5 {print $5}' "$PEAKS_TSV")
		NUCLEAR_PEAK=$(awk '$1 == 1 {print $5}' "$PEAKS_TSV")
		MITOCHONDRIAL_PEAK=$(awk '$1 == 2 {print $5}' "$PEAKS_TSV")
		PLASTID_PEAK=$(awk '$1 == 3 {print $5}' "$PEAKS_TSV")
		if [[ -z "$ORGANELLE_VALLEY" ]]; then
			DEPTH_THRESHOLD=100
			echo "⚠️ Could not detect peaks. Fallback depth threshold = $DEPTH_THRESHOLD"
		else
			DEPTH_THRESHOLD="$ORGANELLE_VALLEY"
		fi

		echo "✅ Nuclear peak ≈ ${NUCLEAR_PEAK}x"
		echo "✅ Mitochondrial peak ≈ ${MITOCHONDRIAL_PEAK}x"
		echo "✅ Plastid peak ≈ ${PLASTID_PEAK}x"
		echo "✅ Organelle valley ≈ ${ORGANELLE_VALLEY}x"
		echo "✅ Plastid valley ≈ ${PLASTID_VALLEY}x"
		echo "➡️  Using depth threshold: ${DEPTH_THRESHOLD}x"
		_polap_log0 "output: ${PDF}"
	fi

}

_polap_count_view_kmer() {

	# === Defaults ===
	local kmerdir
	local K
	local K_HIFI
	local K_ERR
	local K_ONT_RAW
	local JFSIZE
	local THREADS
	local PREFIX
	local PEAK_OVERRIDE
	local X_MIN
	local X_MAX

	# === Parse args ===
	local READS
	local TECH
	local PREFIX

	# === Output folders ===
	local KMER_TABLE
	local HISTO
	local PDF
	local PEAKS_TSV
	local FILTERED
	local DISCARDED
	local LABELS_TSV

	_polap_count_kmer_jellyfish

}

_polap_count_filter_kmer() {

	# === Defaults ===
	local kmerdir
	local K
	local K_HIFI
	local K_ERR
	local K_ONT_RAW
	local JFSIZE
	local THREADS
	local PREFIX
	local PEAK_OVERRIDE
	local X_MIN
	local X_MAX

	# === Parse args ===
	local READS
	local TECH
	local PREFIX

	# === Output folders ===
	local KMER_TABLE
	local HISTO
	local PDF
	local PEAKS_TSV
	local FILTERED
	local DISCARDED
	local LABELS_TSV

	_polap_count_kmer_jellyfish

	_polap_log1 "Filter reads based on mean k-mer depth..."

	REPORT_TSV="${kmerdir}/${PREFIX}_depth_report.tsv"

	# Step 1: Estimate mean k-mer depth
	_polap_log1 "Step 1: Estimate mean k-mer depths of all the reads: $READS"
	python "${_POLAPLIB_DIR}"/polap-py-estimate-kmer-depth.py \
		"$KMER_TABLE" "$READS" "$K" "$REPORT_TSV"

	# Step 2: Filter reads by depth threshold
	_polap_log1 "Step 2: Filter reads by the depth threshold: ${DEPTH_THRESHOLD}x"
	python "${_POLAPLIB_DIR}"/polap-py-filter-by-kmer-depth.py \
		"$READS" "$REPORT_TSV" "$DEPTH_THRESHOLD" "${kmerdir}/$PREFIX"

	_polap_log1 "Step 3: Generate read_labels.tsv for PCA/UMAP..."
	(
		zcat "$FILTERED" | awk 'NR % 4 == 1 {gsub(/^@/, "", $1); print $1 "\tplastid_or_mito"}'
		zcat "$DISCARDED" | awk 'NR % 4 == 1 {gsub(/^@/, "", $1); print $1 "\tnuclear"}'
	) >"$LABELS_TSV"

	_polap_log1 "📁 Output files in: ${kmerdir}/"
	_polap_log0 " - Organelle reads     -> $FILTERED"
	_polap_log1 " - Nuclear reads       -> $DISCARDED"
	_polap_log1 " - Read labels (TSV)   -> $LABELS_TSV"
	_polap_log1 " - Depth report        -> ${REPORT_TSV}"
	_polap_log1 " - Discarded IDs       -> ${kmerdir}/${PREFIX}_discarded.ids.txt"
	_polap_log1 " - k-mer histogram     -> $PDF"
}
