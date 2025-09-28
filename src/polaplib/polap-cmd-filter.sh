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

function _run_polap_filter {
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
		cat <<HEREDOC
Template for an external shell script

Arguments:
  -i ${_arg_inum}: source Flye (usually whole-genome) assembly number

Inputs:
  ${_polap_var_ga_annotation_all}

Outputs:
  ${_polap_var_mtcontigname}

See:
  run-polap-select-contigs-by-table-1.R for the description of --select-contig option

Example:
$(basename $0) ${_arg_menu[0]} [-i|--inum <arg>] [-j|--jnum <arg>] [--select-contig <number>]
HEREDOC
	)

	help_message=$(
		cat <<'EOF'
Name:
  count - k-mers

Synopsis:
  polap count [ -acdfhklLnNrtvV19 ] [ name ... ]

Description:
  polap count uses JellyFish to count k-mers.
  Option 1. we could do the blast to select pt reads.

Options:
  -l FASTQ
    Long-read data

  -x-min FLOAT
    Left boundary of the X-axis in the k-mers distribution

  -x-max FLOAT
    Right boundary of the X-axis in the k-mers distribution

  --pacbio-hifi, --pacbio-raw, --nano-raw, --nano-hq
    Sequence data type

  -k INT
    k-mer size, k=21 for --pacbio-hifi, k=17 for others,
    this option overrides k givent by sequencing data type

  -c INT
    Depth to filter out nuclear reads

  -y-min FLOAT
    Upper boundary of the Y-axis in PCA, UMAP, and the like

  -y-max FLOAT
    Lower boundary of the Y-axis in PCA, UMAP, and the like

Examples:
  Count k-mers:
    polap count --pacbio-hifi -l o/reads.fq
    polap count --nano-raw -l SRR7153095.fastq -o Eucalyptus_pauciflora-5

  Filter by k-mers depth:
    polap filter depth --pacbio-hifi -l o/reads.fq -c 50

  Filter by k-mers composition:
    polap filter pca -l o/reads.fq -k 2

  Filter by k-mers umap:
    polap filter umap -l o/reads.fq --outfile o/reads.organelle.fq
      --x-min 1 --x-max 5 --y-min 1 --y-max 10 -k 2

  Filter by k-mers dbscan:
    polap filter dbscan -l o/reads.fq -k 2

  Flye assembly:
    polap total-length-long -l o/reads.organelle.fq
    polap find-genome-size-for-pacbio -l o/reads.organelle.fq
    polap reduce-data -l o/reads.organelle.fq
    polap flye1 --pacbio-hifi

  Assemble ptDNA:
    polap blast-pt -l o/reads.organelle.fq -o o2

  Filter by ptDNA depth:
    polap filter depth --pacbio-hifi -l o/reads.fq -c 50

  Filter by reads:
  	cut -f5 o/kmer/selected.tsv |
		seqtk subseq o/reads.fq - \
			>o/reads_selected.fq

  Subtract ptDNA reads:
    polap filter reference -l o/reads.fq --reference ref-plastid.gfa
	    -> o/ref-filtered.fastq.gz

  Subtract ptDNA HiFi reads:
    polap filter reference hifi -l o/reads.fq --reference ref-plastid.gfa
	    -> o/ref-filtered.fastq.gz

  Select ptDNA gfa:
    use bandage to filter by depth
    automatically select ptDNA

  Remove nuclear reads or outliers based on k-mers composition:
    polap filter rmkc -l o/reads.fq

  Assemble mtDNA:

  Assembly using PCA:
   cut -f5 Eucalyptus_pauciflora-5/kmer/SRR7153095.k2.pca.selected.tsv | seqtk subseq SRR7153095.fastq - > SRR7153095.k2.pca.selected.fastq
   flye --nano-raw SRR7153095.k2.pca.selected.fastq -t 56 --asm-coverage 30 -g 3m --out-dir Eucalyptus_pauciflora-5/0

  

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

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
		exit $EXIT_SUCCESS
	fi

	# count and this by-kmer functions should be merged
	if [[ "${_arg_menu[1]}" == "depth" ]]; then
		_polap_filter-reads-by-kmer
	fi

	# by-pca, by-dbscan should be merged
	if [[ "${_arg_menu[1]}" == "pca" ]]; then
		_polap_filter-reads-by-pca
	fi

	if [[ "${_arg_menu[1]}" == "short" ]]; then
		_polap_filter-reads-by-pca-short
	fi

	if [[ "${_arg_menu[1]}" == "umap" ]]; then
		_polap_filter-reads-by-umap
	fi

	if [[ "${_arg_menu[1]}" == "dbscan" ]]; then
		_polap_filter-reads-by-dbscan
	fi

	if [[ "${_arg_menu[1]}" == "reference" ]]; then
		if [[ "${_arg_menu[2]}" == "hifi" ]]; then
			_polap_filter-reads-by-reference
		elif [[ "${_arg_menu[2]}" == "ont" ]]; then
			_polap_filter-ont-reads-by-reference
		fi
	fi

	if [[ "${_arg_menu[1]}" == "rmkc" ]]; then
		_polap_filter-reads-by-rmkc
	fi

	if [[ "${_arg_menu[1]}" == "pt" ]]; then
		_polap_lib_conda-ensure_conda_env polap || exit 1
		bash "${_POLAPLIB_DIR}/polap-bash-remove-ptdna-reads.sh" \
			-r "${_arg_long_reads}" \
			-p "${_arg_pt_ref}" \
			-o "${_arg_outdir}/pt-filter" \
			--min-ident 0.70 \
			--min-qcov 0.40 \
			--min-mapq 0
		conda deactivate
	fi

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

x_polap_filter-reads-by-kmer() {

	local kmerdir="${_arg_outdir}"/kmer

	# === Defaults ===
	local K="${_arg_knum}"
	local K_HIFI=21
	local K_ERR=17
	local K_ONT_RAW=13
	local JFSIZE=1G
	local THREADS="${_arg_threads}"
	local PREFIX=""
	local PEAK_OVERRIDE=""
	local X_MIN="${_arg_x_min}"
	local X_MAX="${_arg_x_max}"

	# === Parse args ===
	local READS="${_arg_long_reads}"
	local TECH="${_arg_data_type}"

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

	[[ -z "$PREFIX" ]] && PREFIX="k${K}"

	# === Output folders ===
	mkdir -p "${kmerdir}/figures"
	KMER_TABLE="${kmerdir}/k${K}_all_reads_kmers.txt"
	HISTO="${kmerdir}/k${K}_reads.histo"
	PDF="${kmerdir}/figures/k${K}_kmer_histogram.pdf"
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

	if [[ ! -s "$PEAKS_TSV" ]]; then
		_polap_log1 "Detect peaks in the k-mers distribution..."
		Rscript "${_POLAPLIB_DIR}"/polap-r-plot-kmer-histogram-peaks-valleys.R "$HISTO" "$PDF" "$X_MIN" "$X_MAX" >"$PEAKS_TSV" 2>/dev/null
	fi

	# === Use override or compute ===
	if [[ -n "$PEAK_OVERRIDE" ]]; then
		echo "⚙️ Skipping Jellyfish and using manual depth threshold: $PEAK_OVERRIDE"
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
	_polap_log0 " - Organelle reads     → $FILTERED"
	_polap_log1 " - Nuclear reads       → $DISCARDED"
	_polap_log1 " - Read labels (TSV)   → $LABELS_TSV"
	_polap_log1 " - Depth report        → ${REPORT_TSV}"
	_polap_log1 " - Discarded IDs       → ${kmerdir}/${PREFIX}_discarded.ids.txt"
	_polap_log1 " - k-mer histogram     → $PDF"
}

_polap_filter-reads-by-pca() {
	local kmerdir="${_arg_outdir}"/kmer

	local INFILE="${_arg_long_reads}"
	local KMERLEN="${_arg_knum}"
	if [[ "${KMERLEN}" == "1" ]]; then
		KMERLEN="2"
	fi

	mkdir -p "${kmerdir}"

	local LABELS_TSV="${kmerdir}/read_labels.tsv"

	# === Setup ===
	local BASENAME=$(basename "$INFILE")
	# local PREFIX="${kmerdir}/${BASENAME%.*}.k$KMERLEN"
	local PREFIX="${kmerdir}/k$KMERLEN"

	# === Run R script
	Rscript "${_POLAPLIB_DIR}/polap-r-kmer-dimensionality-reduction.R" \
		"$INFILE" "$PREFIX" "$KMERLEN" "$LABELS_TSV" 2>"${_polap_output_dest}"
}

_polap_filter-reads-by-pca-short() {
	local kmerdir="${_arg_outdir}"/kmer-short

	local INFILE="${_arg_long_reads}"
	local KMERLEN="${_arg_knum}"
	if [[ "${KMERLEN}" == "1" ]]; then
		KMERLEN="2"
	fi

	mkdir -p "${kmerdir}"

	local LABELS_TSV="${kmerdir}/read_labels.tsv"

	# === Setup ===
	local BASENAME=$(basename "$INFILE")
	# local PREFIX="${kmerdir}/${BASENAME%.*}.k$KMERLEN"
	local PREFIX="${kmerdir}/k$KMERLEN"

	# === Run R script
	Rscript "${_POLAPLIB_DIR}/polap-r-kmer-dimensionality-reduction.R" \
		"$INFILE" "$PREFIX" "$KMERLEN" "$LABELS_TSV" 2>"${_polap_output_dest}"
}

# automatic selection of reads using pca/umap with dbscan
_polap_filter-reads-by-dbscan() {
	local kmerdir="${_arg_outdir}"/kmer

	local INFILE="${_arg_long_reads}"
	local KMERLEN="${_arg_knum}"
	if [[ "${KMERLEN}" == "1" ]]; then
		KMERLEN="2"
	fi

	local LABELS_TSV="${kmerdir}/read_labels.tsv"

	# === Setup ===
	local BASENAME=$(basename "$INFILE")
	local PREFIX="${kmerdir}/${BASENAME%.*}.k$KMERLEN"

	# === Run R script
	Rscript "${_POLAPLIB_DIR}/polap-r-kmer-dimensionality-reduction-dbscan.R" \
		"$INFILE" "$PREFIX" "$KMERLEN" "$LABELS_TSV"
}

#
_polap_filter-reads-by-umap() {
	local kmerdir="${_arg_outdir}"/kmer

	local INFILE="${_arg_long_reads}"
	local KMERLEN="${_arg_knum}"
	if [[ "${KMERLEN}" == "1" ]]; then
		KMERLEN="2"
	fi

	local LABELS_TSV="${kmerdir}/read_labels.tsv"

	# === Setup ===
	local BASENAME=$(basename "$INFILE")
	local PREFIX="${kmerdir}/${BASENAME%.*}.k$KMERLEN"
	local input_tsv="${PREFIX}.umap.tsv"
	local output_tsv="${PREFIX}.umap.out.tsv"

	Rscript "${_POLAPLIB_DIR}/polap-r-filter-embedding-by-range.R" \
		"${input_tsv}" \
		"${output_tsv}" \
		--invert \
		--x-min "${_arg_x_min}" \
		--x-max "${_arg_x_max}" \
		--y-min "${_arg_y_min}" \
		--y-max "${_arg_y_max}"

	cut -f5 "${output_tsv}" |
		seqtk subseq "${_arg_long_reads}" - \
			>"${_arg_outfile}"

	# flye --pacbio-hifi selected_reads.fq -t 56 --asm-coverage 30 -g 1000000 --out-dir pt
}

_polap_filter-ont-reads-by-reference() {

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
	FASTQ_REMOVED="${OUTPUT_PREFIX}/ref-chloroplast_like.fastq.gz"
	FASTQ_FILTERED="${OUTPUT_PREFIX}/ref-filtered.fastq.gz"

	# === Step 1: Align reads with GraphAligner
	echo "📌 Aligning ONT reads with GraphAligner..."
	# GraphAligner \
	# 	-g "$CHLORO_GRAPH" \
	# 	-f "$MITO_READS" \
	# 	-a "$SAM" \
	# 	-t "$THREADS" \
	# 	--precise-clipping 0.9 \
	# 	-x vg
	#
	GraphAligner \
		-g "$CHLORO_GRAPH" \
		-f "$MITO_READS" \
		-a "$SAM" \
		-t "$THREADS" \
		-x dbg
	# --seeds-minimizer-length 15 \
	# --seeds-minimizer-windowsize 50 \
	# --bandwidth 50 \
	# --min-alignment-score 2000 \

	# --seeds-mxm-length 19 \
	# --bandwidth 15 \

	# === Step 2: Parse and filter with Python
	echo "🧠 Parsing SAM with pysam..."
	python "${_POLAPLIB_DIR}"/polap-py-parse-graphaligner-gaf.py \
		"$SAM" "$REJECT_IDS" "$RETAIN_IDS" "$SUMMARY" "$ID_THRESH" "$CLIP_THRESH"

	# === Step 3: Extract FASTQ subsets
	echo "✂️ Splitting FASTQ..."
	seqkit grep -f "$REJECT_IDS" "$MITO_READS" | gzip -c >"$FASTQ_REMOVED"
	seqkit grep -f "$RETAIN_IDS" "$MITO_READS" | gzip -c >"$FASTQ_FILTERED"

	# === Step 4: Plot with R
	echo "📊 Generating plots..."
	Rscript --vanialla "${_POLAPLIB_DIR}"/polap-r-plot-graphaligner-summary.R \
		"$SUMMARY" "$PLOT"

	# === Summary
	_polap_log0 "✅ Finished filtering mitochondrial reads:"
	_polap_log0 "- Filtered reads: $FASTQ_FILTERED"
	_polap_log0 "- Removed reads: $FASTQ_REMOVED"
	_polap_log0 "- Summary table: $SUMMARY"
	_polap_log0 "- Summary plot:  $PLOT"

}

# NOTE: use _polap_lib_filter-reads-by-rmkc at polap-lib-filter.sh
_polap_filter-reads-by-rmkc() {
	local READS="${1:-${_arg_long_reads}}"

	_polap_log0 "[Deprecated function: use _polap_lib_filter-reads-by-rmkc at polap-lib-filter.sh]"

	# === Defaults ===
	MEM="100M"
	MINLEN=0
	VERBOSE=false
	LOGFILE=""
	K=21

	# Input variables
	# READS="${_arg_long_reads}"

	PREFIX="${_arg_outdir}/kmer/rmkc"

	mkdir -p "${_arg_outdir}/kmer"

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
cat(sprintf('✅ mkc = %.2f\\n', mkc));"

	local MKC=$(cat "$MKC_OUT")
	local LKC=$(awk -v mkc="$MKC" 'BEGIN { printf "%.2f", 0.3 * mkc }')
	local HKC=$(awk -v mkc="$MKC" 'BEGIN { printf "%.2f", 5.0 * mkc }')

	_polap_log0 "ℹ️  mkc = $MKC | lkc = $LKC | hkc = $HKC"

	# === Step 5: Filter nuclear outliers
	_polap_log1 "[5/5] Filtering out nuclear-like reads..."
	python "${_POLAPLIB_DIR}"/polap-py-remove-nuclear-outlier-reads.py \
		--reads "$READS" \
		--kmer-counts "$KMER_DUMP" \
		--k "$K" \
		--mkc "$MKC" \
		--minlen "$MINLEN" \
		--output "$PREFIX.cleaned"

	_polap_log1 "✅ Final output:"
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
	# 	_polap_log1 "✅ Report created: ${PREFIX}.rmkc.summary.html"

}
