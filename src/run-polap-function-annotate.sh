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
source "$script_dir/run-polap-function-include.sh"
_POLAP_INCLUDE_=$(_polap_include "${BASH_SOURCE[0]}")
set +u
[[ -n "${!_POLAP_INCLUDE_}" ]] && return 0
set -u
declare "$_POLAP_INCLUDE_=1"
#
################################################################################

################################################################################
#
################################################################################
function _run_polap_edges-stats() { # create an edge version of contigs_stats.txt
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "$script_dir/polap-variables-common.sh"

	# Print help message if requested
	help_message=$(
		cat <<HEREDOC
# Create an edge version of contigs_stats.txt: edges_stats.txt
#
# Arguments:
#   -i $INUM: source Flye (usually whole-genome) assembly number
# Inputs:
#   ${_polap_var_contigger_gfa}
# Outputs:
#   ${_polap_var_contigger_edges_stats}
Example: $0 ${_arg_menu[0]} -i ${INUM} [-o ${ODIR}]
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		if [[ -s "${_polap_var_contigger_edges_stats}" ]]; then
			_polap_log0_column "${_polap_var_contigger_edges_stats}"
		else
			_polap_log0 "No such file: ${_polap_var_contigger_edges_stats}"
		fi

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
		return
	fi

	mkdir -p "${_polap_var_ga_mtcontigs}"

	_polap_log0 "creating edges_stats.txt from graph_final.gfa ..."
	_polap_log1 "  creating GFA without sequence data: ${_polap_var_ga_gfa_all}"
	_polap_log2 "    input: ${_polap_var_contigger_gfa}"
	_polap_log2 "    output: ${_polap_var_ga_gfa_all}"

	if [[ ! -s "${_polap_var_contigger_gfa}" ]]; then
		_polap_log0 "ERROR: no such file: ${_polap_var_contigger_gfa}"
		return $RETURN_FAIL
	fi

	if [ -s "${_polap_var_ga_gfa_all}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log1 "    found: ${_polap_var_ga_gfa_all}, so skipping ..."
	else
		_polap_log3_pipe "gfatools view \
		  -S ${_polap_var_contigger_gfa} \
		  >${_polap_var_ga_gfa_all} \
		  2>$_polap_output_dest"
	fi

	_polap_log1 "  extracting sequence part of GFA: ${_polap_var_ga_gfa_seq_part}"
	_polap_log2 "    input: ${_polap_var_ga_gfa_all}"
	_polap_log2 "    output: ${_polap_var_ga_gfa_seq_part}"
	if [ -s "${_polap_var_ga_gfa_seq_part}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log1 "    found: ${_polap_var_ga_gfa_seq_part}, so skipping ..."
	else
		_polap_log3_pipe "grep ^S ${_polap_var_ga_gfa_all} >${_polap_var_ga_gfa_seq_part}"
	fi

	# Filter edges in GFA using depths.
	_polap_log1 "  filtering GFA sequence part using depth range"
	_polap_log2 "    input1: ${_polap_var_ga_gfa_seq_part}"
	_polap_log2 "    output1: ${_polap_var_contigger_edges_stats}"
	if [ -s "${_polap_var_contigger_edges_stats}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log1 "    found: ${_polap_var_contigger_edges_stats}, so skipping ..."
	else
		_polap_log3_pipe "Rscript $script_dir/run-polap-r-edges-stats.R \
		--gfa ${_polap_var_ga_gfa_seq_part} \
		--out ${_polap_var_contigger_edges_stats} \
		2>$_polap_output_dest"
	fi

	if [[ -s "${_polap_var_contigger_edges_stats}" ]]; then
		_polap_log2_column "${_polap_var_contigger_edges_stats}"
	else
		_polap_log2 "ERROR: no such file: ${_polap_var_contigger_edges_stats}"
	fi

	_polap_log1 "NEXT: $0 blast-genome -o $ODIR [-i ${INUM}]"
	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
# creates a depth distribution
################################################################################
function _run_polap_depth-distribution() { # creates a depth distribution
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "$script_dir/polap-variables-common.sh"

	# Print help message if requested
	help_message=$(
		cat <<HEREDOC
# Draw depth distributions of an annotation table.
#
# FIXME: Distribution by the cutoff of MT gene groups Copy average - 3 x SD
# FIXME: depth or copy? 
# FIXME: what is it for? The purpose or intended use of this function?
#
# Arguments:
#   -i $INUM: source Flye (usually whole-genome) assembly number
# Inputs:
#   ${_polap_var_annotation_table}
# Outputs:
#   ${_polap_var_annotation_table}.pdf
Example: $(basename $0) ${_arg_menu[0]}
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		_polap_log0_file "${_polap_var_annotation_table}".pdf
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
		exit $EXIT_SUCCESS
	fi

	Rscript "$script_dir"/run-polap-r-depth-distribution.R \
		-t "${_polap_var_annotation_table}" \
		-o "${_polap_var_annotation_table}".pdf \
		2>"$_polap_output_dest"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
# Annotates edge sequences of a flye genome assembly.
################################################################################
function _run_polap_annotate() { # annotate edge sequences in edges_stats.txt
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "$script_dir/polap-variables-common.sh"

	help_message=$(
		cat <<HEREDOC
# Annotate a Flye genome assembly.
#
# Arguments:
#   -i $INUM: a Flye genome assembly number
# Inputs:
#   $INUM
# Outputs:
#   ${_polap_var_ga_annotation_all}
#   ${_polap_var_ga_annotation}
#   ${_polap_var_ga_annotation_depth_table}
#   ${_polap_var_ga_annotation_table}
#   ${_polap_var_ga_pt_annotation_depth_table}
# View:
#   all
#   mt
#   table
#   pt-table
#   no-depth
Example: $0 ${_arg_menu[0]} -i ${INUM}
Example: $0 ${_arg_menu[0]} view depth
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		case "${_arg_menu[2]}" in
		all)
			_polap_log0_column "${_polap_var_ga_annotation_all}"
			;;
		table)
			_polap_log0_column "${_polap_var_ga_annotation_depth_table}"
			;;
		no-depth)
			_polap_log0_column "${_polap_var_ga_annotation_table}"
			;;
		contig)
			_polap_log0_column "${_polap_var_ga_annotation_table_contig}"
			;;
		mt)
			_polap_log0_column "${_polap_var_ga_annotation}"
			;;
		pt-table)
			_polap_log0_column "${_polap_var_ga_pt_annotation_depth_table}"
			;;
		*)
			_polap_log0 "menu3: all, table, no-depth, mt, pt-table"
			;;
		esac

		# _polap_log3_cmd touch "${ODIR}"
		# _polap_log3_cmd ln -f -s "${_polap_var_contigger_gfa#*/}" "${ODIR}"/wga.gfa

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	_polap_log0 "annotating edge sequence contigs with mitochondrial and plastid genes on the assembly number ${INUM} ..."

	if [ -s "${_polap_var_ga_annotation_all}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log0 "  found: ${_polap_var_ga_annotation_all}, so skipping the annotation ..."
	elif [ -s "${_polap_var_ga_annotation_all}" ] && [ "${_arg_redo}" = "on" ]; then
		_polap_log0 "  found: ${_polap_var_ga_annotation_all}, with --redo option, so backup the annotation ..."
		_polap_log3_cmd cp -p ${_polap_var_ga_annotation_all} ${_polap_var_ga_annotation_all_backup}
		_run_polap_blast-genome
		_run_polap_count-gene
	else
		_run_polap_blast-genome
		_run_polap_count-gene
	fi

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
# Counts genes annotated on a genome assembly.
#
# Arguments:
#   -i $INUM: a Flye genome assembly number
# Inputs:
#   ${_polap_var_ga}/30-contigger/edges_stats.txt
#   $MTGENECOUNT
#   $PTGENECOUNT
# Outputs:
#   ${_polap_var_ga}/mt.contig.name-1
#   ${_polap_var_ga}/mt.contig.name-2
#   ${_polap_var_ga}/assembly_info_organelle_annotation_count.txt
#   ${_polap_var_ga}/assembly_info_organelle_annotation_count-all.txt
#   ${_polap_var_ga}/contig-annotation-table.txt
################################################################################
function _run_polap_count-gene() { # count MT and PT genes using edges_stats.txt
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "$script_dir/polap-variables-common.sh"

	help_message=$(
		cat <<HEREDOC
# Count genes annotated on a genome assembly.
#
# Arguments:
#   -i $INUM: a Flye genome assembly number
# Inputs:
#   ${_polap_var_ga}/30-contigger/edges_stats.txt
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
# View:
#   all
#   short
#   table
#   depth
Example: $(basename "$0") ${_arg_menu[0]} -i ${INUM}
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		# Display the BLAST genome output.
		#
		case "${_arg_menu[2]}" in
		all)
			if [[ -s "${_polap_var_ga_annotation_all}" ]]; then
				_polap_log0_column "${_polap_var_ga_annotation_all}"
			else
				_polap_log0 "No such file: ${_polap_var_ga_annotation_all}"
			fi
			;;
		short)
			if [[ -s "${_polap_var_ga_annotation}" ]]; then
				_polap_log0_column "${_polap_var_ga_annotation}"
			else
				_polap_log0 "No such file: ${_polap_var_ga_annotation}"
			fi
			;;
		table)
			if [[ -s "${_polap_var_ga_annotation_table}" ]]; then
				_polap_log0_column "${_polap_var_ga_annotation_table}"
			else
				_polap_log0 "No such file: ${_polap_var_ga_annotation_table}"
			fi
			;;
		depth)
			if [[ -s "${_polap_var_ga_annotation_depth_table}" ]]; then
				_polap_log0_column "${_polap_var_ga_annotation_depth_table}"
			else
				_polap_log0 "No such file: ${_polap_var_ga_annotation_depth_table}"
			fi
			;;
		*)
			_polap_log0 "all, short, table, depth"
			;;
		esac

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
		exit $EXIT_SUCCESS
	fi

	_polap_log0 "counting mitochondrial and plastid genes on the assembly number ${INUM} ..."

	_polap_log1 "  input1: ${_polap_var_ga}/30-contigger/edges_stats.txt"
	_polap_log1 "  input2: ${_polap_var_ann_MTGENECOUNT}"
	_polap_log1 "  input3: ${_polap_var_ann_PTGENECOUNT}"
	_polap_log1 "  output1: ${_polap_var_ga_annotation}"
	_polap_log1 "  output2: ${_polap_var_ga_annotation_all}"
	_polap_log1 "  output3: ${_polap_var_ga_annotation_table}"
	_polap_log1 "  output4: ${_polap_var_ga_annotation_depth_table}"

	# Checks the output files earlier than the input.
	if [[ -s "${_polap_var_ga_annotation_all}" ]] &&
		[[ "${_arg_redo}" = "off" ]]; then
		_polap_log0 "  found1: ${_polap_var_ga_annotation_all}"
		_polap_log0 "  so skipping the blast genome ..."
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
		return
	fi

	_polap_log1 "  creating annotation tables"
	_polap_log2 "    input1: ${_polap_var_contigger_edges_stats}"
	_polap_log2 "    input2: ${_polap_var_ann_MTGENECOUNT}"
	_polap_log2 "    input3: ${_polap_var_ann_PTGENECOUNT}"
	_polap_log2 "    output1: ${_polap_var_ga_annotation}"
	_polap_log2 "    output2: ${_polap_var_ga_annotation_all}"
	_polap_log2 "    output3: ${_polap_var_ga_annotation_table}"
	_polap_log2 "    output4: ${_polap_var_ga_annotation_depth_table}"
	_polap_log2 "    output5: ${_polap_var_ga_pt_annotation_depth_table}"
	local _command1="Rscript $script_dir/run-polap-r-mtcontig.R \
		--flyeout-edges-stats ${_polap_var_contigger_edges_stats} \
    --mt-gene-count ${_polap_var_ann_MTGENECOUNT} \
    --pt-gene-count ${_polap_var_ann_PTGENECOUNT} \
		--out-annotation ${_polap_var_ga_annotation} \
		--out-annotation-all ${_polap_var_ga_annotation_all} \
		--out-annotation-table ${_polap_var_ga_annotation_table} \
		--out-annotation-depth-table ${_polap_var_ga_annotation_depth_table} \
		--out-pt-annotation-depth-table ${_polap_var_ga_pt_annotation_depth_table} \
		--contigger"
	if [[ "${_arg_plastid}" = "on" ]]; then
		_command1+=" \
        --plastid"
	fi
	_command1+=" \
				2>${_polap_output_dest}"
	_polap_log3_pipe "${_command1}"

	if [[ -s "${_polap_var_ga_annotation_all}" ]]; then
		_polap_log3_column "${_polap_var_ga_annotation_all}"
	fi

	# if [[ ${_arg_test} == "on" ]]; then
	# 	_polap_log0 "creating ${_polap_var_ga}/mt.contig.name-1 for testing purpose ..."
	# 	_polap_log0 "  you would have to edit it for a real data-set."
	# 	echo edge_1 >"${_polap_var_ga}"/mt.contig.name-1
	# 	echo edge_2 >>"${_polap_var_ga}"/mt.contig.name-1
	# 	echo edge_3 >>"${_polap_var_ga}"/mt.contig.name-1
	# fi

	ANUMNEXT=$((INUM + 1))
	_polap_log1 NEXT: $0 seeds -o "$ODIR" -i $INUM -j $ANUMNEXT

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
# Blasts the genome assembly of a Flye run.
# Arguments:
#   -i $INUM: a Flye genome assembly number
# Inputs:
#   graph_final.fasta not contigs.fasta
#
#   $assembly_contigs_fasta
# Outputs:
#   $CONTIGNAME
#   $MTGENECOUNT
#   $PTGENECOUNT
#   $ADIR/mtaa.blast
#   $ADIR/mtaa.blast.bed
#   $MTAABLAST.sorted.bed
#   $ADIR/mtaa.bed
#   $ADIR/ptaa.blast
#   $ADIR/ptaa.blast.bed
#   $PTAABLAST.sorted.bed
#   $ADIR/ptaa.bed
################################################################################
function _run_polap_blast-genome() { # BLAST edge sequences on MT and PT genes
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "$script_dir/polap-variables-common.sh"
	local MTAA="$script_dir"/polap-mt.1.c70.3.faa
	local PTAA="$script_dir"/polap-pt.2.c70.3.faa

	help_message=$(
		cat <<HEREDOC
# NCBI BLAST the genome assembly of a Flye run againt the plant organelle genes.
#
# Arguments:
#   -i $INUM: a Flye genome assembly number
# Inputs:
#   ${_polap_var_contigger_edges_fasta}
#   ${_polap_var_contigger_edges_stats}
# Outputs:
#   ${_polap_var_ann_MTGENECOUNT}
#   ${_polap_var_ann_PTGENECOUNT}
Example: $(basename "$0") ${_arg_menu[0]} -i ${INUM}
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		# Display the BLAST genome output.
		if [[ -s "${_polap_var_ann_MTGENECOUNT}" ]]; then
			_polap_log0_head "${_polap_var_ann_MTGENECOUNT}"
		else
			_polap_log0 "No such file: ${_polap_var_ann_MTGENECOUNT}"
		fi

		if [[ -s "${_polap_var_ann_PTGENECOUNT}" ]]; then
			_polap_log0_head "${_polap_var_ann_PTGENECOUNT}"
		else
			_polap_log0 "No such file: ${_polap_var_ann_PTGENECOUNT}"
		fi

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
		return
	fi

	_polap_log0 "blasting edge contigs with mitochondrial and plastid genes on the assembly number ${INUM} ..."
	_polap_log1 "  input1: ${_polap_var_contigger_edges_stats}"
	_polap_log1 "  input2: ${_polap_var_contigger_edges_fasta}"

	# Checks the output files earlier than the input.
	if [[ -s "${_polap_var_ann_MTGENECOUNT}" ]] &&
		[[ -s "${_polap_var_ann_PTGENECOUNT}" ]] &&
		[[ "${_arg_redo}" = "off" ]]; then
		_polap_log0 "  found1: ${_polap_var_ann_MTGENECOUNT}"
		_polap_log0 "  found2: ${_polap_var_ann_PTGENECOUNT}"
		_polap_log0 "  so skipping the blast genome ..."
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
		return
	fi

	if [ ! -s "${_polap_var_contigger_edges_stats}" ]; then
		_polap_log0 "ERROR: no edges_stats.txt file: ${_polap_var_contigger_edges_stats}"
		exit $EXIT_ERROR
	fi

	_polap_log2 "  deleting and recreating the folder: ${_polap_var_ann}"
	_polap_log3_cmd rm -rf "${_polap_var_ann}"
	_polap_log3_cmd mkdir -p "${_polap_var_ann}"

	#src/run-polap-select.R o/30-contigger/edges_stats.txt o/50-annotation/contig.name
	_polap_log2 "  contig sequence names in file: ${_polap_var_ann_CONTIGNAME}"
	_polap_log2 "    input: ${_polap_var_contigger_edges_stats}"
	_polap_log2 "    output: ${_polap_var_ann_CONTIGFILE}"
	_polap_log3_pipe "grep -v '#' ${_polap_var_contigger_edges_stats} |
		cut -f 1 >${_polap_var_ann_CONTIGNAME}"

	_polap_log2 "  copying edge contig sequence file: ${_polap_var_ann_CONTIGFILE}"
	_polap_log3_cmd cp "${_polap_var_contigger_edges_fasta}" "${_polap_var_ann_CONTIGFILE}"

	_polap_log2 "  making BLASTDB of the contig sequences: ${_polap_var_ann_CONTIGDB}"
	_polap_log2 "    input: ${_polap_var_ann_CONTIGFILE}"
	_polap_log2 "    output: ${_polap_var_ann_CONTIGDB}"
	_polap_log3_pipe "makeblastdb -dbtype nucl \
		-in ${_polap_var_ann_CONTIGFILE} \
		-out ${_polap_var_ann_CONTIGDB} \
		2>${_polap_output_dest}"

	# Mitochondrial gene annotation and counts
	_polap_log2 "  BLAST of the mitochondrial proteins against ${_polap_var_ann_CONTIGDB}"
	_polap_log2 "    input query: ${MTAA}"
	_polap_log2 "    input BLAST DB: ${_polap_var_ann_CONTIGDB}"
	_polap_log2 "    evalue cutoff: 1e-30"
	_polap_log2 "    number of CPUs: ${_arg_threads}"
	_polap_log3 "    executing the tblastn ... be patient!"
	_polap_log3_pipe "tblastn -query ${MTAA} \
		-db ${_polap_var_ann_CONTIGDB} \
		-out ${_polap_var_ann_MTAABLAST} \
		-evalue 1e-30 \
		-outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle salltitles' \
		-num_threads ${_arg_threads} \
		>${_polap_output_dest} 2>&1"

	_polap_log2 "  converting BLAST result in BED format for removing redundancy"
	_polap_log2 "    input1: ${_polap_var_ann_MTAABLAST}"
	_polap_log2 "    output1: ${_polap_var_ann_MTAABLASTBED}"
	_polap_log3_pipe "Rscript $script_dir/run-polap-r-genes.R \
		${_polap_var_ann_MTAABLAST} \
		${_polap_var_ann_MTAABLASTBED} \
		>${_polap_output_dest} 2>&1"

	_polap_log2 "  sorting BED format BLAST result"
	_polap_log2 "    input1: ${_polap_var_ann_MTAABLASTBED}"
	_polap_log2 "    output1: ${_polap_var_ann_MTAABLAST}.sorted.bed"
	_polap_log3_pipe "sort -k1,1 -k2,2n \
		${_polap_var_ann_MTAABLASTBED} \
		>${_polap_var_ann_MTAABLAST}.sorted.bed"

	_polap_log2 "  creating folder: ${_polap_var_ann_MTAABED}"
	_polap_log3_cmd mkdir "${_polap_var_ann_MTAABED}"

	_polap_log2 "  counting mitochondrial genes in the contigs ..."
	_polap_log2 "    input1: ${_polap_var_ann_MTAABLAST}.sorted.bed"
	_polap_log2 "    output1: ${_polap_var_ann_MTGENECOUNT}"
	_polap_log3_cmd rm -f "${_polap_var_ann_MTGENECOUNT_TMP}"
	while IFS= read -r contig; do
		_polap_log3_pipe "grep -w ${contig} ${_polap_var_ann_MTAABLAST}.sorted.bed \
			>${_polap_var_ann_MTAABED}/${contig}.bed"
		_polap_log3_pipe "bedtools merge -i ${_polap_var_ann_MTAABED}/${contig}.bed \
			>${_polap_var_ann_MTAABED}/${contig}.bed.txt"
		_polap_log3 commnad: printf \"%s\\t%d\\n\" "${contig}" $(wc -l <"${_polap_var_ann_MTAABED}/${contig}".bed.txt) ">>${_polap_var_ann_MTGENECOUNT_TMP}"
		printf "%s\t%d\n" ${contig} $(wc -l <${_polap_var_ann_MTAABED}/${contig}.bed.txt) >>"${_polap_var_ann_MTGENECOUNT_TMP}"
	done <"${_polap_var_ann_CONTIGNAME}"
	_polap_log3_pipe "sort -k2 -rn ${_polap_var_ann_MTGENECOUNT_TMP} >${_polap_var_ann_MTGENECOUNT}"

	_polap_log2 "  compressing the BLAST results of mitochondrial gene annotation"
	_polap_log2 "    input: ${_polap_var_ann_MTAABED}"
	_polap_log2 "    output: ${_polap_var_ann_MTAABED}.tar.gz"
	_polap_log3_cmd tar zcf "${_polap_var_ann_MTAABED}".tar.gz "${_polap_var_ann_MTAABED}"
	_polap_log2 "  deleting folder: ${_polap_var_ann_MTAABLASTBED}"
	_polap_log3_cmd rm -rf "${_polap_var_ann_MTAABED}"

	_polap_log1 "  output1: ${_polap_var_ann_MTGENECOUNT}"

	# Plastid gene annotation and counts
	_polap_log2 "  BLAST of the plastid proteins against ${_polap_var_ann_CONTIGDB}"
	_polap_log2 "    input query: ${PTAA}"
	_polap_log2 "    input BLAST DB: ${_polap_var_ann_CONTIGDB}"
	_polap_log2 "    evalue cutoff: 1e-30"
	_polap_log2 "    number of CPUs: ${_arg_threads}"
	_polap_log3 "    executing the tblastn ... be patient!"
	_polap_log3_pipe "tblastn -query ${PTAA} \
		-db ${_polap_var_ann_CONTIGDB} \
		-out ${_polap_var_ann_PTAABLAST} \
		-evalue 1e-30 \
		-outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle salltitles' \
		-num_threads ${_arg_threads} \
		>${_polap_output_dest} 2>&1"

	_polap_log2 "  converting BLAST result in BED format for removing redundancy"
	_polap_log2 "    input1: ${_polap_var_ann_PTAABLAST}"
	_polap_log2 "    output1: ${_polap_var_ann_PTAABLASTBED}"
	_polap_log3_pipe "Rscript $script_dir/run-polap-r-genes.R \
		${_polap_var_ann_PTAABLAST} \
		${_polap_var_ann_PTAABLASTBED} \
		>${_polap_output_dest} 2>&1"

	_polap_log2 "  sorting BED format BLAST result"
	_polap_log2 "    input1: ${_polap_var_ann_PTAABLASTBED}"
	_polap_log2 "    output1: ${_polap_var_ann_PTAABLAST}.sorted.bed"
	_polap_log3_pipe "sort -k1,1 -k2,2n \
		${_polap_var_ann_PTAABLASTBED} \
		>${_polap_var_ann_PTAABLAST}.sorted.bed"

	_polap_log2 "  creating folder: ${_polap_var_ann_PTAABED}"
	_polap_log3_cmd mkdir "${_polap_var_ann_PTAABED}"

	_polap_log2 "  counting mitochondrial genes in the contigs ..."
	_polap_log2 "    input1: ${_polap_var_ann_PTAABLAST}.sorted.bed"
	_polap_log2 "    output1: ${_polap_var_ann_PTGENECOUNT}"
	_polap_log3_cmd rm -f "${_polap_var_ann_PTGENECOUNT_TMP}"
	while IFS= read -r contig; do
		_polap_log3_pipe "grep -w ${contig} ${_polap_var_ann_PTAABLAST}.sorted.bed \
			>${_polap_var_ann_PTAABED}/${contig}.bed"
		_polap_log3_pipe "bedtools merge -i ${_polap_var_ann_PTAABED}/${contig}.bed \
			>${_polap_var_ann_PTAABED}/${contig}.bed.txt"
		_polap_log3 commnad: printf \"%s\\t%d\\n\" "${contig}" $(wc -l <"${_polap_var_ann_PTAABED}/${contig}".bed.txt) ">>${_polap_var_ann_PTGENECOUNT_TMP}"
		printf "%s\t%d\n" ${contig} $(wc -l <${_polap_var_ann_PTAABED}/${contig}.bed.txt) >>"${_polap_var_ann_PTGENECOUNT_TMP}"
	done <"${_polap_var_ann_CONTIGNAME}"
	_polap_log3_pipe "sort -k2 -rn ${_polap_var_ann_PTGENECOUNT_TMP} >${_polap_var_ann_PTGENECOUNT}"

	_polap_log2 "  compressing the BLAST results of mitochondrial gene annotation"
	_polap_log2 "    input: ${_polap_var_ann_PTAABED}"
	_polap_log2 "    output: ${_polap_var_ann_PTAABED}.tar.gz"
	_polap_log3_cmd tar zcf "${_polap_var_ann_PTAABED}".tar.gz "${_polap_var_ann_PTAABED}"
	_polap_log2 "  deleting folder: ${_polap_var_ann_PTAABLASTBED}"
	_polap_log3_cmd rm -rf "${_polap_var_ann_PTAABED}"

	_polap_log1 "  output2: ${_polap_var_ann_PTGENECOUNT}"

	# Wait for a while before deleting these plastid part
	# Original Plastid version
	# Plastid gene annotation and counts
	# _polap_log2 "  BLAST of the plastid proteins against ${_polap_var_ann_CONTIGDB}"
	# _polap_log3 "    executing the tblastn ... be patient!"
	# tblastn -query "${PTAA}" \
	# 	-db "${_polap_var_ann_CONTIGDB}" \
	# 	-out "${_polap_var_ann_PTAABLAST}" \
	# 	-evalue 1e-30 \
	# 	-outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle salltitles' \
	# 	-num_threads "${_arg_threads}" \
	# 	>${_polap_output_dest} 2>&1
	#
	# Rscript "$script_dir"/run-polap-r-genes.R \
	# 	"${_polap_var_ann_PTAABLAST}" \
	# 	"${_polap_var_ann_PTAABLASTBED}" \
	# 	>${_polap_output_dest} 2>&1
	#
	# sort -k1,1 -k2,2n \
	# 	"${_polap_var_ann_PTAABLASTBED}" \
	# 	>"${_polap_var_ann_PTAABLAST}".sorted.bed
	#
	# mkdir "${_polap_var_ann_PTAABED}"
	#
	# _polap_log2 "  counting plastid genes in the contigs ..."
	# if [ "$DEBUG" -eq 1 ]; then set +x; fi
	# while IFS= read -r contig; do
	# 	grep -w "${contig}" "${_polap_var_ann_PTAABLAST}".sorted.bed \
	# 		>"${_polap_var_ann_PTAABED}/${contig}".bed
	# 	bedtools merge -i "${_polap_var_ann_PTAABED}/${contig}".bed \
	# 		>"${_polap_var_ann_PTAABED}/${contig}".bed.txt
	# 	printf "%s\t%d\n" "${contig}" \
	# 		$(wc -l <"${_polap_var_ann_PTAABED}/${contig}".bed.txt)
	# done <"${_polap_var_ann_CONTIGNAME}" |
	# 	sort -k2 -rn >"${_polap_var_ann_PTGENECOUNT}"
	# if [ "$DEBUG" -eq 1 ]; then set -x; fi
	#
	# _polap_log2 "  compressing the BLAST results of plastid gene annotation"
	# tar zcf "${_polap_var_ann_PTAABED}".tar.gz "${_polap_var_ann_PTAABED}"
	# rm -rf "${_polap_var_ann_PTAABED}"
	#
	# _polap_log1 "  output2: ${_polap_var_ann_PTGENECOUNT}"

	_polap_log1 "NEXT (for testing purpose only): $0 count-gene --test"
	_polap_log1 "NEXT: $0 count-gene -o $ODIR [-i $INUM]"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}
