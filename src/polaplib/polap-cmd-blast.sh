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

_polap_lib_blast-fastq_to_gfa() {
	local fastq_file="$1"

	# Create temp files for S and P sections
	local tmp_s
	local tmp_p
	tmp_s=$(mktemp)
	tmp_p=$(mktemp)

	# Stream through FASTQ with seqkit and write S and P entries separately
	seqkit fx2tab -i -q "$fastq_file" | nl -n rz -w 8 |
		while IFS=$'\t' read -r lineno id seq qual; do
			echo -e "S\t$id\t$seq\tdp:i:1" >>"$tmp_s"
			echo -e "P\tcontig_$lineno\t${id}+\t*" >>"$tmp_p"
		done

	# Output final GFA
	echo -e "H\tVN:Z:1.0"
	cat "$tmp_s"
	cat "$tmp_p"

	# Clean up
	rm -f "$tmp_s" "$tmp_p"
}

_polap_lib_blast-annotate() {
	local _outdir="${1}"

	local _contigger_edges_gfa="${_outdir}/30-contigger/graph_final.gfa"
	local _contigger_edges_fasta="${_outdir}/30-contigger/graph_final.fasta"
	local _mtcontigname="${_outdir}/mt.contig.name"
	local _ga_annotation_all="${_outdir}/assembly_info_organelle_annotation_count-all.txt"

	polap_annotate "${_contigger_edges_gfa}" "${_ga_annotation_all}"
}

blast_fasta_to_gfa() {
	local fasta_file="$1"

	echo -e "H\tVN:Z:1.0"

	local ids=()

	# Print S lines and store sequence IDs
	while IFS=$'\t' read -r id seq; do
		echo -e "S\t$id\t$seq\tdp:i:32"
		ids+=("$id")
	done < <(seqkit fx2tab "$fasta_file")

	# Print P lines using stored IDs
	local i=1
	for id in "${ids[@]}"; do
		echo -e "P\tcontig_$i\t${id}+\t*"
		((i++))
	done
}

makeblastdb_nosplit() {
	local fasta="$1"
	local out="${2:-${fasta%.fa}}"

	if [[ ! -f "$fasta" ]]; then
		echo "Error: File not found: $fasta" >&2
		return 1
	fi

	local size_bytes
	size_bytes=$(stat --format=%s "$fasta")

	# Add 20% buffer
	local safe_bytes
	safe_bytes=$(awk -v s="$size_bytes" 'BEGIN {printf "%.0f", s * 1.2}')

	local max_sz
	if [[ "$safe_bytes" -ge 1000000000 ]]; then
		max_sz=$(awk -v s="$safe_bytes" 'BEGIN {printf "%.0fGB", s/1e9}')
	else
		max_sz=$(awk -v s="$safe_bytes" 'BEGIN {printf "%.0fMB", s/1e6}')
	fi

	# echo "Building BLAST DB:"
	# echo "  Input FASTA:     $fasta"
	# echo "  Output prefix:   $out"
	# echo "  Input size:      $size_bytes bytes"
	# echo "  max_file_sz:     $max_sz"
	# echo

	makeblastdb -dbtype nucl -in "$fasta" -out "$out" -max_file_sz "$max_sz"
}

polap_lib_blast() {
	local _ga="${1}"

	local _ann="${_ga}"
	# mkdir -p "${_ann}"
	local _ann_CONTIGDB="${_ann}/l"

	local MTAA="${_POLAPLIB_DIR}"/polap-mt.1.c70.3.faa
	local PTAA="${_POLAPLIB_DIR}"/polap-pt.2.c70.3.faa
	# local _ann_CONTIGDB="${_ann}/contig"
	local _ann_CONTIGFILE="${_ann}/contig.fasta"
	local _ann_CONTIGNAME="${_ann}/contig.name"
	local _ann_MTAABED="${_ann}/mtaa.bed"
	local _ann_MTAABLASTBED="${_ann}/mtaa.blast.bed"
	local _ann_MTAABLAST="${_ann}/mtaa.blast"
	local _ann_MTGENECOUNT_TMP="${_ann}/mt.gene.count.tmp"
	local _ann_PTAABED="${_ann}/ptaa.bed"
	local _ann_PTAABLASTBED="${_ann}/ptaa.blast.bed"
	local _ann_PTAABLAST="${_ann}/ptaa.blast"
	local _ann_PTGENECOUNT_TMP="${_ann}/pt.gene.count.tmp"

	tblastn -query ${MTAA} \
		-db ${_ann_CONTIGDB} \
		-out ${_ann_MTAABLAST} \
		-evalue 1e-30 \
		-outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle salltitles' \
		-num_threads ${_arg_threads}

	tblastn -query ${PTAA} \
		-db ${_ann_CONTIGDB} \
		-out ${_ann_PTAABLAST} \
		-evalue 1e-30 \
		-outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle salltitles' \
		-num_threads ${_arg_threads}

	cut -f2 ${_ann_MTAABLAST} | sort | uniq >"${_ga}"/mt.names
	cut -f2 ${_ann_PTAABLAST} | sort | uniq >"${_ga}"/pt.names
}

polap_lib_blast-mt() {
	local _ga="${1}"

	local _ann="${_ga}"
	# mkdir -p "${_ann}"
	local _ann_CONTIGDB="${_ann}/l"

	local MTAA="${_POLAPLIB_DIR}"/polap-mt.1.c70.3.faa
	local PTAA="${_POLAPLIB_DIR}"/polap-pt.2.c70.3.faa
	# local _ann_CONTIGDB="${_ann}/contig"
	local _ann_CONTIGFILE="${_ann}/contig.fasta"
	local _ann_CONTIGNAME="${_ann}/contig.name"
	local _ann_MTAABED="${_ann}/mtaa.bed"
	local _ann_MTAABLASTBED="${_ann}/mtaa.blast.bed"
	local _ann_MTAABLAST="${_ann}/mtaa.blast"
	local _ann_MTGENECOUNT_TMP="${_ann}/mt.gene.count.tmp"
	local _ann_PTAABED="${_ann}/ptaa.bed"
	local _ann_PTAABLASTBED="${_ann}/ptaa.blast.bed"
	local _ann_PTAABLAST="${_ann}/ptaa.blast"
	local _ann_PTGENECOUNT_TMP="${_ann}/pt.gene.count.tmp"

	tblastn -query ${MTAA} \
		-db ${_ann_CONTIGDB} \
		-out ${_ann_MTAABLAST} \
		-evalue 1e-30 \
		-outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle salltitles' \
		-num_threads ${_arg_threads}

	cut -f2 ${_ann_MTAABLAST} | sort | uniq >"${_ga}"/mt.names
}

polap_lib_blast-pt() {
	local _ga="${1}"

	local _ann="${_ga}"
	# mkdir -p "${_ann}"
	local _ann_CONTIGDB="${_ann}/l"

	local MTAA="${_POLAPLIB_DIR}"/polap-mt.1.c70.3.faa
	local PTAA="${_POLAPLIB_DIR}"/polap-pt.2.c70.3.faa
	# local _ann_CONTIGDB="${_ann}/contig"
	local _ann_CONTIGFILE="${_ann}/contig.fasta"
	local _ann_CONTIGNAME="${_ann}/contig.name"
	local _ann_MTAABED="${_ann}/mtaa.bed"
	local _ann_MTAABLASTBED="${_ann}/mtaa.blast.bed"
	local _ann_MTAABLAST="${_ann}/mtaa.blast"
	local _ann_MTGENECOUNT_TMP="${_ann}/mt.gene.count.tmp"
	local _ann_PTAABED="${_ann}/ptaa.bed"
	local _ann_PTAABLASTBED="${_ann}/ptaa.blast.bed"
	local _ann_PTAABLAST="${_ann}/ptaa.blast"
	local _ann_PTGENECOUNT_TMP="${_ann}/pt.gene.count.tmp"

	tblastn -query ${PTAA} \
		-db ${_ann_CONTIGDB} \
		-out ${_ann_PTAABLAST} \
		-evalue 1e-30 \
		-outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle salltitles' \
		-num_threads ${_arg_threads}

	cut -f2 ${_ann_PTAABLAST} | sort | uniq >"${_ga}"/pt.names
}

function _run_polap_blast {
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
Select reads using homology search

Arguments:
  -l ${_arg_long_reads}

Inputs:
  long-read data

Outputs:
  ${_arg_long_reads}-mt.fq.gz
  ${_arg_long_reads}-pt.fq.gz

Example:
$(basename $0) ${_arg_menu[0]} -l l.fq
$(basename $0) ${_arg_menu[0]} -l SRR15206231.fastq -o Solanum_tuberosum
HEREDOC
	)

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

	local long_reads_name="${_arg_long_reads%%.*}"

	#
	local blastdir="${_arg_outdir}"/blast
	mkdir -p "${blastdir}"

	if [[ ! -s "${blastdir}/l.fa" ]]; then
		seqtk seq -A "${_arg_long_reads}" >"${blastdir}"/l.fa
	fi

	makeblastdb -dbtype nucl -in "${blastdir}"/l.fa -out "${blastdir}"/l
	# makeblastdb -dbtype nucl -in "${blastdir}"/l.fa -out "${blastdir}"/l -max_file_sz 4
	# makeblastdb_nosplit "${blastdir}"/l.fa

	polap_lib_blast "${blastdir}"
	seqtk subseq "${_arg_long_reads}" "${blastdir}"/mt.names | gzip >"${_arg_long_reads}"_mt.fq.gz
	seqtk subseq "${_arg_long_reads}" "${blastdir}"/pt.names | gzip >"${_arg_long_reads}"_pt.fq.gz
	flye "${_arg_flye_data_type}" "${long_reads_name}"_mt.fq.gz --out-dir "${blastdir}"/mt
	flye "${_arg_flye_data_type}" "${long_reads_name}"_pt.fq.gz --out-dir "${blastdir}"/pt

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_polap_blast-mt {
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
Select reads using homology search

Arguments:
  -l ${_arg_long_reads}

Inputs:
  long-read data

Outputs:
  ${_arg_long_reads}-mt.fq.gz
  ${_arg_long_reads}-pt.fq.gz

Example:
$(basename $0) ${_arg_menu[0]} -l l.fq
$(basename $0) ${_arg_menu[0]} -l SRR15206231.fastq -o Solanum_tuberosum
HEREDOC
	)

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

	local long_reads_name="${_arg_long_reads%%.*}"

	#
	local blastdir="${_arg_outdir}"/blast
	mkdir -p "${blastdir}"

	if [[ ! -s "${blastdir}/l.fa" ]]; then
		seqtk seq -A "${_arg_long_reads}" >"${blastdir}"/l.fa
	fi

	makeblastdb -dbtype nucl -in "${blastdir}"/l.fa -out "${blastdir}"/l
	# makeblastdb_nosplit "${blastdir}"/l.fa

	polap_lib_blast-mt "${blastdir}"
	seqtk subseq "${_arg_long_reads}" "${blastdir}"/mt.names | gzip >"${_arg_long_reads}"_mt.fq.gz
	flye "${_arg_flye_data_type}" "${long_reads_name}"_mt.fq.gz --out-dir "${blastdir}"/mt

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}
function _run_polap_blast-pt {
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
Select reads using homology search

Arguments:
  -l ${_arg_long_reads}

Inputs:
  long-read data

Outputs:
  ${_arg_long_reads}-mt.fq.gz
  ${_arg_long_reads}-pt.fq.gz

Example:
$(basename $0) ${_arg_menu[0]} -l l.fq
$(basename $0) ${_arg_menu[0]} -l SRR15206231.fastq -o Solanum_tuberosum
HEREDOC
	)

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

	local long_reads_name="${_arg_long_reads%.*}"

	local blastdir="${_arg_outdir}"/blast
	mkdir -p "${blastdir}/30-contigger"

	if [[ "${_arg_redo}" == "on" ]]; then
		rm -rf "${blastdir}"
	fi

	local _contigger_edges_gfa="${blastdir}/30-contigger/graph_final.gfa"
	if [[ ! -s "${_contigger_edges_gfa}" ]]; then
		_polap_lib_blast-fastq_to_gfa "${_arg_long_reads}" >"${_contigger_edges_gfa}"
	fi

	local _pt_table="${blastdir}/pt-contig-annotation-depth-table.txt"
	local _mt_table="${blastdir}/contig-annotation-depth-table.txt"
	local _og_table="${blastdir}/assembly_info_organelle_annotation_count.txt"
	if [[ -s "${_pt_table}" ]]; then
		echo "found: ${_pt_table}"
	else
		_polap_lib_blast-annotate "${blastdir}"
	fi

	rm -rf "${blastdir}"/pt
	local total_bases_pt_reads=$(awk 'NR > 1 { sum += $2 } END { print sum }' "${_pt_table}")
	local total_bases_pt_reads_gb=$(_polap_lib_unit-convert_bp "${total_bases_pt_reads}")
	if ((total_bases_pt_reads < 1000000)); then
		echo "[INFO] PT long-reads are not enough: ${total_bases_pt_reads_gb} < 1 Mb"
		return 1
	else
		echo "[INFO] PT long-reads: ${total_bases_pt_reads_gb}"
	fi

	# genome size: 150k x 2
	awk '{print $1}' "${_pt_table}" | tail -n +2 >"${blastdir}"/pt.names
	seqtk subseq "${_arg_long_reads}" "${blastdir}"/pt.names | gzip >"${long_reads_name}"_pt.fq.gz

	command time -v flye "${_arg_flye_data_type}" "${long_reads_name}"_pt.fq.gz -t ${_arg_threads} --asm-coverage 30 -g 300000 --out-dir "${blastdir}"/pt 2>"${_arg_outdir}/timing-flye-blast-pt.txt"

	# genome size: 500k x 2
	awk '{print $1}' "${_mt_table}" | tail -n +2 >"${blastdir}"/mt.names
	seqtk subseq "${_arg_long_reads}" "${blastdir}"/mt.names | gzip >"${long_reads_name}"_mt.fq.gz
	command time -v flye "${_arg_flye_data_type}" "${long_reads_name}"_mt.fq.gz -t ${_arg_threads} --asm-coverage 30 -g 1000000 --out-dir "${blastdir}"/mt 2>"${_arg_outdir}/timing-flye-blast-mt.txt"

	awk '{print $1}' "${_og_table}" | tail -n +2 >"${blastdir}"/og.names
	# seqtk subseq "${_arg_long_reads}" "${blastdir}"/og.names | gzip >"${long_reads_name}"_og.fq.gz
	# flye "${_arg_flye_data_type}" "${long_reads_name}"_og.fq.gz -t ${_arg_threads} -g 300000 --out-dir "${blastdir}"/og

	# if [[ ! -s "${blastdir}/l.fa" ]]; then
	# 	seqtk seq -A "${_arg_long_reads}" >"${blastdir}"/l.fa
	# fi
	# makeblastdb -dbtype nucl -in "${blastdir}"/l.fa -out "${blastdir}"/l
	# polap_lib_blast-pt "${blastdir}"
	# seqtk subseq "${_arg_long_reads}" "${blastdir}"/pt.names | gzip >"${_arg_long_reads}"_pt.fq.gz
	# flye "${_arg_flye_data_type}" "${long_reads_name}"_pt.fq.gz --out-dir "${blastdir}"/pt

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}
