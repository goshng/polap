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
# Annotate contigs not edge sequences.
# See Also:
# polaplib/polap-cmd-annotate.sh
#
# TODO: delete this script in the future.
#
# The edge sequence version in polaplib/polap-cmd-annotate.sh replaces
# this script.
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

function _run_polap_blast-genome-contig { # v0.2.6: BLAST contig sequences on MT and PT genes
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"
	local MTAA="${_POLAPLIB_DIR}"/polap-mt.1.c70.3.faa
	local PTAA="${_POLAPLIB_DIR}"/polap-pt.2.c70.3.faa

	help_message=$(
		cat <<HEREDOC
# Blasts the genome assembly of a Flye run againt the plant organelle genes.
#
# Arguments:
#   -i ${_arg_inum}: a Flye genome assembly number
# Inputs:
#   ${_polap_var_ga_contigger_contigs_stats}
# Outputs:
#   ${_polap_var_ann_MTGENECOUNT}
#   ${_polap_var_ann_PTGENECOUNT}
Example: $(basename "$0") ${_arg_menu[0]} [-i|--inum <arg>]
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		# Display the BLAST genome output.

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		exit $EXIT_SUCCESS
	fi

	_polap_log0 "blasting contigs with mitochondrial and plastid genes on the assembly number ${_arg_inum} ..."

	# Checks the output files earlier than the input.
	if [[ -s "${_polap_var_ann_MTGENECOUNT}" ]] &&
		[[ -s "${_polap_var_ann_PTGENECOUNT}" ]] &&
		[ "${_arg_redo}" = "off" ]; then
		_polap_log0 "  found1: ${_polap_var_ann_MTGENECOUNT}"
		_polap_log0 "  found2: ${_polap_var_ann_PTGENECOUNT}"
		_polap_log0 "  so skipping the blast genome ..."
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return
	fi

	if [ ! -s "${_polap_var_ga_contigger_contigs_stats}" ]; then
		_polap_log0 "ERROR: no contig_stats.txt file: ${_polap_var_ga_contigger_contigs_stats}"
		exit $EXIT_ERROR
	fi

	_polap_log2 "  deleting and recreating the folder: ${_polap_var_ann}"
	if [ -d "${_polap_var_ann}" ]; then
		_polap_log3 rm -rf "${_polap_var_ann}"
		rm -rf "${_polap_var_ann}"
	fi
	_polap_log3 mkdir -p "${_polap_var_ann}"
	mkdir -p "${_polap_var_ann}"

	_polap_log1 "  input1: ${_polap_var_ga_contigger_contigs_stats}"
	_polap_log1 "  input2: ${_polap_var_ga_contigger_contigs_fasta}"

	#src/run-polap-select.R o/30-contigger/contigs_stats.txt o/50-annotation/contig.name
	_polap_log2 "  contig sequence names in file: ${_polap_var_ann_CONTIGNAME}"
	grep -v "#" "${_polap_var_ga_contigger_contigs_stats}" |
		cut -f 1 >"${_polap_var_ann_CONTIGNAME}"

	# seqkit grep --threads ${_arg_threads} -f "$CONTIGNAME" \
	# 	"$assembly_contigs_fasta" \
	# 	-o "${_polap_var_ann}"/contig.fasta \
	# 	>/dev/null 2>&1
	_polap_log2 "  contig sequence file: ${_polap_var_ann_CONTIGFILE}"
	cp "${_polap_var_ga_contigger_contigs_fasta}" "${_polap_var_ann_CONTIGFILE}"

	_polap_log2 "  making BLASTDB of the contig sequences: ${_polap_var_ann_CONTIGDB}"
	_polap_log3 makeblastdb -dbtype nucl \
		-in "${_polap_var_ga_contigger_contigs_fasta}" \
		-out "${_polap_var_ann_CONTIGDB}" \
		>${_polap_output_dest} 2>&1
	makeblastdb -dbtype nucl \
		-in "${_polap_var_ga_contigger_contigs_fasta}" \
		-out "${_polap_var_ann_CONTIGDB}" \
		>${_polap_output_dest} 2>&1

	# Mitochondrial gene annotation and counts
	_polap_log2 "  BLAST of the mitochondrial proteins against ${_polap_var_ann_CONTIGDB}"
	_polap_log3 "    executing the tblastn ... be patient!"
	tblastn -query "${MTAA}" \
		-db "${_polap_var_ann_CONTIGDB}" \
		-out "${_polap_var_ann_MTAABLAST}" \
		-evalue 1e-30 \
		-outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle salltitles' \
		-num_threads "${_arg_threads}" \
		>${_polap_output_dest} 2>&1

	Rscript "${_POLAPLIB_DIR}"/polap-r-genes.R \
		"${_polap_var_ann_MTAABLAST}" \
		"${_polap_var_ann_MTAABLASTBED}" \
		>${_polap_output_dest} 2>&1

	sort -k1,1 -k2,2n \
		"${_polap_var_ann_MTAABLASTBED}" \
		>"${_polap_var_ann_MTAABLAST}".sorted.bed

	mkdir "${_polap_var_ann_MTAABED}"

	_polap_log2 "  counting mitochondrial genes in the contigs ..."
	if [ "$_POLAP_DEBUG" -eq 1 ]; then set +x; fi
	while IFS= read -r contig; do
		grep -w "${contig}" "${_polap_var_ann_MTAABLAST}".sorted.bed \
			>"${_polap_var_ann_MTAABED}/${contig}".bed
		bedtools merge -i "${_polap_var_ann_MTAABED}/${contig}".bed \
			>"${_polap_var_ann_MTAABED}/${contig}".bed.txt
		printf "%s\t%d\n" "${contig}" \
			$(wc -l <"${_polap_var_ann_MTAABED}/${contig}".bed.txt)
	done <"${_polap_var_ann_CONTIGNAME}" |
		sort -k2 -rn >"${_polap_var_ann_MTGENECOUNT}"
	if [ "$_POLAP_DEBUG" -eq 1 ]; then set -x; fi

	_polap_log2 "  compressing the BLAST results of mitochondrial gene annotation"
	tar zcf "${_polap_var_ann_MTAABED}".tar.gz "${_polap_var_ann_MTAABED}"
	rm -rf "${_polap_var_ann_MTAABED}"

	_polap_log1 "  output1: ${_polap_var_ann_MTGENECOUNT}"

	# Plastid gene annotation and counts
	_polap_log2 "  BLAST of the plastid proteins against ${_polap_var_ann_CONTIGDB}"
	_polap_log3 "    executing the tblastn ... be patient!"
	tblastn -query "${PTAA}" \
		-db "${_polap_var_ann_CONTIGDB}" \
		-out "${_polap_var_ann_PTAABLAST}" \
		-evalue 1e-30 \
		-outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle salltitles' \
		-num_threads "${_arg_threads}" \
		>${_polap_output_dest} 2>&1

	Rscript "${_POLAPLIB_DIR}"/polap-r-genes.R \
		"${_polap_var_ann_PTAABLAST}" \
		"${_polap_var_ann_PTAABLASTBED}" \
		>${_polap_output_dest} 2>&1

	sort -k1,1 -k2,2n \
		"${_polap_var_ann_PTAABLASTBED}" \
		>"${_polap_var_ann_PTAABLAST}".sorted.bed

	mkdir "${_polap_var_ann_PTAABED}"

	_polap_log2 "  counting plastid genes in the contigs ..."
	if [ "$_POLAP_DEBUG" -eq 1 ]; then set +x; fi
	while IFS= read -r contig; do
		grep -w "${contig}" "${_polap_var_ann_PTAABLAST}".sorted.bed \
			>"${_polap_var_ann_PTAABED}/${contig}".bed
		bedtools merge -i "${_polap_var_ann_PTAABED}/${contig}".bed \
			>"${_polap_var_ann_PTAABED}/${contig}".bed.txt
		printf "%s\t%d\n" "${contig}" \
			$(wc -l <"${_polap_var_ann_PTAABED}/${contig}".bed.txt)
	done <"${_polap_var_ann_CONTIGNAME}" |
		sort -k2 -rn >"${_polap_var_ann_PTGENECOUNT}"
	if [ "$_POLAP_DEBUG" -eq 1 ]; then set -x; fi

	_polap_log2 "  compressing the BLAST results of plastid gene annotation"
	tar zcf "${_polap_var_ann_PTAABED}".tar.gz "${_polap_var_ann_PTAABED}"
	rm -rf "${_polap_var_ann_PTAABED}"

	_polap_log1 "  output2: ${_polap_var_ann_PTGENECOUNT}"

	_polap_log1 "NEXT (for testing purpose only): $(basename "$0") count-gene --test"
	_polap_log1 "NEXT: $(basename $0) count-gene -o ${_arg_outdir} [-i ${_arg_inum}]"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
}

function _run_polap_count-gene-contig { # v0.2.6: count MT and PT genes
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	help_message=$(
		cat <<HEREDOC
# Counts genes annotated on a genome assembly.
#
# Arguments:
#   -i ${_arg_inum}: a Flye genome assembly number
# Inputs:
#   ${_polap_var_ga}/30-contigger/contigs_stats.txt
#   ${_polap_var_ann_MTGENECOUNT}
#   ${_polap_var_ann_PTGENECOUNT}
# Outputs:
#   ${_polap_var_ga}/mt.contig.name-1
#   ${_polap_var_ga}/mt.contig.name-2
#   ${_polap_var_ga}/assembly_info_organelle_annotation_count.txt
#   ${_polap_var_ga}/assembly_info_organelle_annotation_count-all.txt
#   ${_polap_var_ga}/contig-annotation-table.txt 
#   ${_polap_var_ga}/contig-annotation-depth-table.txt 
# Usage:
#   assembly graph: ${_polap_var_ga}"/30-contigger/graph_final.gfa
# 	column -t ${_polap_var_ga}/assembly_info_organelle_annotation_count-all.txt | less -S
#   column -t ${_polap_var_ga}/assembly_info_organelle_annotation_count.txt | less -S
#   column -t ${_polap_var_ga}/contig-annotation-table.txt | less -S
# Example file: ${_polap_var_ga}/mt.contig.name-1
# 	edge_1
#   edge_2
#   edge_3
# edit ${_polap_var_ga}/mt.contig.name-1 for mtDNA contig candidates
# edit ${_polap_var_ga}/mt.contig.name-<destination flye number> for mtDNA contig candidates
Example: $(basename "$0") ${_arg_menu[0]} [-i|--inum <arg>]
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		# Display the BLAST genome output.

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		exit $EXIT_SUCCESS
	fi

	_polap_log0 "counting mitochondrial and plastid genes on the assembly number ${_arg_inum} ..."

	_polap_log0 "  input1: ${_polap_var_ga}/30-contigger/contigs_stats.txt"
	_polap_log0 "  input2: ${_polap_var_ann_MTGENECOUNT}"
	_polap_log0 "  input3: ${_polap_var_ann_PTGENECOUNT}"

	Rscript "${_POLAPLIB_DIR}"/polap-r-mtcontig-contig.R \
		"${_polap_var_ga}" \
		"${_polap_var_ga}"/mt.contig.name \
		"${_polap_var_ga_annotation}" \
		"${_polap_var_ga_annotation_all}" \
		--contigger \
		>${_polap_output_dest} 2>&1

	_polap_log0 "  output1: ${_polap_var_ga_annotation}"
	_polap_log0 "  output2: ${_polap_var_ga_annotation_all}"
	_polap_log0 "  output3: ${_polap_var_ga_annotation_table}"
	_polap_log0 "  output4: ${_polap_var_ga_annotation_depth_table}"

	if [[ ${_arg_test} == "on" ]]; then
		_polap_log0 "creating ${_polap_var_ga}/mt.contig.name-1 for testing purpose ..."
		_polap_log0 "  you would have to edit it for a real data-set."
		echo edge_1 >"${_polap_var_ga}"/mt.contig.name-1
		echo edge_2 >>"${_polap_var_ga}"/mt.contig.name-1
		echo edge_3 >>"${_polap_var_ga}"/mt.contig.name-1
	fi

	ANUMNEXT=$((_arg_inum + 1))
	_polap_log1 NEXT: $(basename "$0") select-reads -o "${_arg_outdir}" [-i ${_arg_inum}] [-j $ANUMNEXT]
	_polap_log1 NEXT: $(basename "$0") assemble2 -o "${_arg_outdir}" [-i ${_arg_inum}] [-j $ANUMNEXT]

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
}

################################################################################
# Annotates the genome assembly.
################################################################################
function _run_polap_annotate-contig { # annotate v0.2.6 for contigs_stats.txt
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	help_message=$(
		cat <<HEREDOC
# Annotates the genome assembly using contig sequences not edge ones.
#
# Arguments:
#   -i ${_arg_inum}: a Flye genome assembly number
# Inputs:
#   ${_arg_inum}
# Outputs:
#   ${_polap_var_ga_annotation_all}
#   ${_polap_var_ga_annotation}
#   ${_polap_var_ga_annotation_table}
# View:
#   all   -> ${_polap_var_ga_annotation_all}
#   mt    -> ${_polap_var_ga_annotation}
#   table -> ${_polap_var_ga_annotation_table}
Example: ${_polap_command_string} ${_arg_menu[0]} -i ${_arg_inum}
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then

		if [ ! -s "${_polap_var_ga_annotation_all}" ]; then
			_polap_log0 "ERROR: no such file: ${_polap_var_ga_annotation_all}"
		fi

		if [ ! -s "${_polap_var_ga_annotation}" ]; then
			_polap_log0 "ERROR: no such file: ${_polap_var_ga_annotation}"
		fi

		if [ ! -s "${_polap_var_ga_annotation_table}" ]; then
			_polap_log0 "ERROR: no such file: ${_polap_var_ga_annotation_table}"
		fi

		case "${_arg_menu[2]}" in
		all)
			_polap_log0_column "${_polap_var_ga_annotation_all}"
			;;
		mt)
			_polap_log0_column "${_polap_var_ga_annotation}"
			;;
		table | *)
			if [[ "${_arg_markdown}" == "off" ]]; then
				_polap_log0_column "${_polap_var_ga_annotation_table}"
			else
				csvtk space2tab "${_polap_var_ga_annotation_table}" | csvtk tab2csv | csvtk csv2md -a l,r,r,r,r,r >&3
			fi
			;;
		esac

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		exit $EXIT_SUCCESS
	fi

	_polap_log0 "annotating contigs with mitochondrial and plastid genes on the assembly number ${_arg_inum} ..."

	if [ -s "${_polap_var_ga_annotation_all}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log0 "  found: ${_polap_var_ga_annotation_all}, so skipping the annotation ..."
	else
		_run_polap_blast-genome-contig
		_run_polap_count-gene-contig
	fi

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
}
