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

################################################################################
# Create an edge version of contigs_stats.txt.
#
# Flye provides contigs_stats where you find depths and mulitiplicity values
# for contig sequences. Because graph_final.gfa uses edge sequences rather than
# contig sequences, it would be preferable to have such depths and mulitiplicity
# for edge sequences. Depth values are extracted from graph_final.gfa genome
# assembly graph, but mulitiplicity values are not provided by Flye.
# We divide depth values by the median of all the depth values to have an edge
# version of the mulitiplicity values. See the following reference.
#
# Reference:
# https://github.com/mikolmogorov/Flye/issues/732#issuecomment-2444058108
################################################################################
function _run_polap_edges-stats { # create an edge version of contigs_stats.txt
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	# Print help message if requested
	help_message=$(
		cat <<HEREDOC
# Create an edge version of contigs_stats.txt: edges_stats.txt
#
# Arguments:
#   -i ${_arg_inum}: source Flye (usually whole-genome) assembly number
# Inputs:
#   ${_polap_var_ga_contigger_edges_gfa}
#   or
#   ${_polap_var_ga_assembly_graph_gfa} (not implemented yet!)
# Outputs:
#   ${_polap_var_ga_contigger_edges_stats}
#   or
#   ${_polap_var_ga_assembly_graph_edges_stats} (not implemented yet!)
Example: $(basename "$0") ${_arg_menu[0]} -i ${_arg_inum}
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		if [[ -s "${_polap_var_ga_contigger_edges_stats}" ]]; then
			_polap_log0_column "${_polap_var_ga_contigger_edges_stats}"
		else
			_polap_log0 "No such file: ${_polap_var_ga_contigger_edges_stats}"
		fi

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
		return
	fi

	# mkdir -p "${_polap_var_ga_mtcontigs}"

	_polap_log0 "creating edges_stats.txt from graph_final.gfa ..."
	_polap_log1 "  creating GFA without sequence data: ${_polap_var_ga_gfa_all}"
	_polap_log2 "    input: ${_polap_var_ga_contigger_edges_gfa}"
	_polap_log2 "    output: ${_polap_var_ga_gfa_all}"

	if [[ ! -s "${_polap_var_ga_contigger_edges_gfa}" ]]; then
		_polap_log0 "ERROR: no such file: ${_polap_var_ga_contigger_edges_gfa}"
		return $RETURN_FAIL
	fi

	if [ -s "${_polap_var_ga_gfa_all}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log1 "    found: ${_polap_var_ga_gfa_all}, so skipping ..."
	else
		_polap_log3_pipe "gfatools view \
		  -S ${_polap_var_ga_contigger_edges_gfa} \
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
	_polap_log2 "    output1: ${_polap_var_ga_contigger_edges_stats}"
	if [ -s "${_polap_var_ga_contigger_edges_stats}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log1 "    found: ${_polap_var_ga_contigger_edges_stats}, so skipping ..."
	else
		_polap_log3_pipe "Rscript ${_POLAPLIB_DIR}/run-polap-r-edges-stats.R \
		--gfa ${_polap_var_ga_gfa_seq_part} \
		--out ${_polap_var_ga_contigger_edges_stats} \
		2>$_polap_output_dest"
	fi

	if [[ -s "${_polap_var_ga_contigger_edges_stats}" ]]; then
		_polap_log2_column "${_polap_var_ga_contigger_edges_stats}"
	else
		_polap_log2 "ERROR: no such file: ${_polap_var_ga_contigger_edges_stats}"
	fi

	_polap_log1 "NEXT: $(basename "$0") blast-genome -o ${_arg_outdir} [-i ${_arg_inum}]"
	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
# creates a depth distribution
################################################################################
function _run_polap_depth-distribution { # creates a depth distribution
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	# Print help message if requested
	help_message=$(
		cat <<HEREDOC
# Draw depth distributions of an annotation table.
#
# Arguments:
#   -i ${_arg_inum}: source Flye (usually whole-genome) assembly number
# Inputs:
#   ${_polap_var_ga_annotation_all}
# Outputs:
#   ${_polap_var_ga_annotation_all}.pdf
Example: $(basename $0) ${_arg_menu[0]}
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		_polap_log0_file "${_polap_var_ga_annotation_all}".pdf
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	Rscript "${_POLAPLIB_DIR}"/run-polap-r-depth-distribution.R \
		-t "${_polap_var_ga_annotation_all}" \
		-o "${_polap_var_ga_annotation_all}".pdf \
		2>"$_polap_output_dest"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

# source
# from run-polap-function-annotate.sh
# from function _run_polap_edges-stats
#
# input1: a gfa file
# output: edges_stats.txt
#
# FIXME: the skipping makes trouble in resuming
polap_edges-stats() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	local _contigger_edges_gfa="$1"
	local _contigger_edges_stats="$2"
	local _contigger=$(dirname "$1")

	_polap_log1 "  create edges_stats.txt from graph_final.gfa ..."
	_polap_log2 "    input: ${_contigger_edges_gfa}"

	if [[ ! -s "${_contigger_edges_gfa}" ]]; then
		_polap_log0 "ERROR: no such file: ${_contigger_edges_gfa}"
		return $RETURN_FAIL
	fi

	local _gfa_all="${_contigger}/3-gfa.all.gfa"
	local _gfa_seq_part="${_contigger}/3-gfa.seq.part.tsv"

	_polap_log3_pipe "gfatools view \
		  -S ${_contigger_edges_gfa} \
		  >${_gfa_all} \
		  2>$_polap_output_dest"

	_polap_log2 "  extracting sequence part of GFA: ${_gfa_seq_part}"
	_polap_log2 "    input: ${_gfa_all}"
	_polap_log2 "    output: ${_gfa_seq_part}"
	_polap_log3_pipe "grep ^S ${_gfa_all} >${_gfa_seq_part}"

	# Filter edges in GFA using depths.
	_polap_log2 "  filtering GFA sequence part using depth range"
	_polap_log2 "    input1: ${_gfa_seq_part}"
	_polap_log2 "    output1: ${_contigger_edges_stats}"
	_polap_log3_pipe "Rscript ${_POLAPLIB_DIR}/run-polap-r-edges-stats.R \
		--gfa ${_gfa_seq_part} \
		--out ${_contigger_edges_stats} \
		2>$_polap_output_dest"

	if [[ -s "${_contigger_edges_stats}" ]]; then
		_polap_log3_column "${_contigger_edges_stats}"
	else
		_polap_log2 "ERROR: no such file: ${_contigger_edges_stats}"
	fi

	# Disable debugging if previously enabled
	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

# used by disassemble menu
polap_annotate() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	local _contigger=$(dirname "$1")
	local _ga=$(dirname "${_contigger}")
	local _ann="${_ga}/50-annotation"
	local _contigger_edges_gfa="$1"
	local _annotation_all="$2"

	_polap_log1 "  annotate edge sequence contigs with mitochondrial and plastid genes ..."

	if [ -s "${_annotation_all}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log1 "  found: ${_annotation_all}, so skipping the annotation ..."
	else
		local _contigger_edges_stats="${_contigger}/edges_stats.txt"
		polap_edges-stats "${_contigger_edges_gfa}" "${_contigger_edges_stats}"

		local _ann_MTGENECOUNT="${_ann}/mt.gene.count"
		local _ann_PTGENECOUNT="${_ann}/pt.gene.count"
		polap_blast-genome "${_contigger_edges_gfa}" \
			"${_contigger_edges_stats}" \
			"${_ann_MTGENECOUNT}" \
			"${_ann_PTGENECOUNT}"

		polap_count-gene "${_contigger_edges_stats}" \
			"${_ann_MTGENECOUNT}" \
			"${_ann_PTGENECOUNT}" \
			"${_ga_annotation_all}"
	fi

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}
################################################################################
# Annotates edge sequences of a flye genome assembly.
################################################################################
function _run_polap_annotate { # annotate edge sequences in edges_stats.txt
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	help_message=$(
		cat <<HEREDOC
# Annotate a Flye genome assembly with organelle genes.
#
# Arguments:
#   -i ${_arg_inum}: a Flye genome assembly number
#   --contigger (default: on)
#   --no-contigger (not implemented yet!)
# Inputs:
#   ${_arg_inum}
# Outputs:
#   ${_polap_var_ga_annotation_all}
#   ${_polap_var_ga_annotation}
#   ${_polap_var_ga_annotation_depth_table}
#   ${_polap_var_ga_annotation_table}
#   ${_polap_var_ga_pt_annotation_depth_table}
# View:
#   all      -> ${_polap_var_ga_annotation_all}
#   mt       -> ${_polap_var_ga_annotation}
#   table    -> ${_polap_var_ga_annotation_depth_table}
#   no-depth -> ${_polap_var_ga_annotation_table}
#   seed     -> ${_polap_var_ga_annotation_depth_table_seed_target}
#   pt-table -> ${_polap_var_ga_pt_annotation_depth_table}
Example: $(basename "$0") ${_arg_menu[0]} -i ${_arg_inum}
Example: $(basename "$0") ${_arg_menu[0]} view table
Example: $(basename "$0") ${_arg_menu[0]} view seed
Example: $(basename "$0") ${_arg_menu[0]} view seed -o Spirodela_polyrhiza/o2 --table-format docx --outfile sp.docx
HEREDOC
	)

	table_message=$(
		cat <<HEREDOC
--------------------------------------------------------------------------------
Depth: Read depths of a contig.
Copy: Copy number of the contig.
MT: Number of mitochondrial genes.
PT: Number of plastid genes.
Seed: Denoted as A if the contig is a seed for an organelle-genome assembly and
  has more MT genes than PT ones annotated. Denoted as G if the contig is a seed
  but has no gene annotations. Denoted as X if it is not selected as a seed.
HEREDOC
	)

	table_docx_message=$(
		cat <<HEREDOC

Depth: Read depths of a contig.

Copy: Copy number of the contig.

MT: Number of mitochondrial genes.

PT: Number of plastid genes.

Seed: Denoted as A if the contig is a seed for an organelle-genome assembly and
  has more MT genes than PT ones annotated. Denoted as G if the contig is a seed
  but has no gene annotations. Denoted as X if it is not selected as a seed.
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
		seed)
			if [[ "${_arg_table_format}" == "tsv" ]]; then
				_polap_log0 "$(basename $0) seeds annotation -i ${_arg_inum} -j ${_arg_jnum}"
				_polap_log0_column "${_polap_var_ga_annotation_depth_table_seed_target}"
				_polap_echo0 "${table_message}"
			elif [[ "${_arg_table_format}" == "markdown" ]]; then
				csvtk space2tab "${_polap_var_ga_annotation_depth_table_seed_target}" |
					sed 's/Length/Length (bp)/' |
					csvtk tab2csv | csvtk csv2md -a l,r,r,r,r,r,r >&3
				_polap_echo0 "${table_message}"
			elif [[ "${_arg_table_format}" == "docx" ]]; then
				_polap_log0 "outfile: ${_arg_outfile}"
				csvtk space2tab "${_polap_var_ga_annotation_depth_table_seed_target}" |
					awk -F'\t' 'BEGIN {OFS="\t"} NR==1 {print $0; next} { $2 = sprintf("%'\''d", $2); print }' |
					sed 's/Length/Length (bp)/' |
					pandoc -f tsv -t markdown |
					{
						cat
						echo ""
						echo "${table_docx_message}"
					} |
					pandoc -f markdown -t docx -o "${_arg_outfile}"
			fi
			;;
		no-depth)
			if [[ "${_arg_markdown}" == "off" ]]; then
				_polap_log0_column "${_polap_var_ga_annotation_table}"
			else
				csvtk space2tab "${_polap_var_ga_annotation_table}" | csvtk tab2csv | csvtk csv2md >&3
			fi
			;;
		contig)
			_polap_log0_column "${_polap_var_ga_annotation_table_contig}"
			;;
		mt)
			_polap_log0_column "${_polap_var_ga_annotation}"
			;;
		pt-table)
			if [[ "${_arg_markdown}" == "off" ]]; then
				_polap_log0_column "${_polap_var_ga_pt_annotation_depth_table}"
			else
				csvtk space2tab "${_polap_var_ga_pt_annotation_depth_table}" | csvtk tab2csv | csvtk csv2md >&3
			fi
			;;
		table | *)
			if [[ "${_arg_markdown}" == "off" ]]; then
				_polap_log0_column "${_polap_var_ga_annotation_depth_table}"
			else
				csvtk space2tab "${_polap_var_ga_annotation_depth_table}" | csvtk tab2csv | csvtk csv2md -a l,r,r,r,r,r,r >&3
			fi
			;;
		esac

		# _polap_log3_cmd touch "${_arg_outdir}"
		# _polap_log3_cmd ln -f -s "${_polap_var_ga_contigger_edges_gfa#*/}" "${_arg_outdir}"/wga.gfa

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	_polap_log0 "annotating edge sequence contigs with mitochondrial and plastid genes on the assembly number ${_arg_inum} ..."

	if [ -s "${_polap_var_ga_annotation_all}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log0 "  found: ${_polap_var_ga_annotation_all}, so skipping the annotation ..."
	elif [ -s "${_polap_var_ga_annotation_all}" ] && [ "${_arg_redo}" = "on" ]; then
		_polap_log0 "  found: ${_polap_var_ga_annotation_all}, with --redo option, so backup the annotation ..."
		_polap_log3_cmd cp -p ${_polap_var_ga_annotation_all} ${_polap_var_ga_annotation_all_backup}
		_run_polap_edges-stats
		_run_polap_blast-genome
		_run_polap_count-gene
	else
		_run_polap_edges-stats
		_run_polap_blast-genome
		_run_polap_count-gene
	fi

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

polap_count-gene() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	local _contigger=$(dirname "$1")
	local _ga=$(dirname "${_contigger}")
	local _ga_contigger_edges_stats="$1"
	local _ann_MTGENECOUNT=$2
	local _ann_PTGENECOUNT=$3
	local _ga_annotation_all="$4"

	# local _ga_annotation_all="${_ga}/assembly_info_organelle_annotation_count-all.txt"
	local _ga_annotation="${_ga}/assembly_info_organelle_annotation_count.txt"
	local _ga_annotation_depth_table="${_ga}/contig-annotation-depth-table.txt"
	local _ga_annotation_table="${_ga}/contig-annotation-table.txt"
	local _ga_pt_annotation_depth_table="${_ga}/pt-contig-annotation-depth-table.txt"

	_polap_log1 "  count mitochondrial and plastid genes ..."

	# Checks the output files earlier than the input.
	if [[ -s "${_ga_annotation_all}" ]] &&
		[[ "${_arg_redo}" = "off" ]]; then
		_polap_log0 "  found1: ${_ga_annotation_all}, so skipping the blast genome ..."
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	_polap_log1 "  creating annotation tables"
	_polap_log2 "    input1: ${_ga_contigger_edges_stats}"
	_polap_log2 "    input2: ${_ann_MTGENECOUNT}"
	_polap_log2 "    input3: ${_ann_PTGENECOUNT}"
	_polap_log2 "    output1: ${_ga_annotation}"
	_polap_log2 "    output2: ${_ga_annotation_all}"
	_polap_log2 "    output3: ${_ga_annotation_table}"
	_polap_log2 "    output4: ${_ga_annotation_depth_table}"
	_polap_log2 "    output5: ${_ga_pt_annotation_depth_table}"
	local _command1="Rscript ${_POLAPLIB_DIR}/run-polap-r-mtcontig.R \
		--flyeout-edges-stats ${_ga_contigger_edges_stats} \
    --mt-gene-count ${_ann_MTGENECOUNT} \
    --pt-gene-count ${_ann_PTGENECOUNT} \
		--out-annotation ${_ga_annotation} \
		--out-annotation-all ${_ga_annotation_all} \
		--out-annotation-table ${_ga_annotation_table} \
		--out-annotation-depth-table ${_ga_annotation_depth_table} \
		--out-pt-annotation-depth-table ${_ga_pt_annotation_depth_table} \
		--contigger"
	if [[ "${_arg_plastid}" = "on" ]]; then
		_command1+=" \
      --plastid"
	fi
	_command1+=" \
      2>${_polap_output_dest}"
	_polap_log3_pipe "${_command1}"

	if [[ -s "${_ga_annotation_all}" ]]; then
		_polap_log3_column "${_ga_annotation_all}"
	fi

	# Disable debugging if previously enabled
	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
# Counts genes annotated on a genome assembly.
#
# Arguments:
#   -i ${_arg_inum}: a Flye genome assembly number
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
function _run_polap_count-gene { # count MT and PT genes using edges_stats.txt
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	help_message=$(
		cat <<HEREDOC
# Count genes annotated on a genome assembly.
#
# Arguments:
#   -i ${_arg_inum}: a Flye genome assembly number
# Inputs:
#   ${_polap_var_ga}/30-contigger/edges_stats.txt
#   ${_polap_var_ann_MTGENECOUNT}
#   ${_polap_var_ann_PTGENECOUNT}
# Outputs:
#   ${_polap_var_ga_annotation_all}"
Example: $(basename "$0") ${_arg_menu[0]} -i ${_arg_inum}
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	_polap_log0 "counting mitochondrial and plastid genes on the assembly number ${_arg_inum} ..."

	# Checks the output files earlier than the input.
	if [[ -s "${_polap_var_ga_annotation_all}" ]] &&
		[[ "${_arg_redo}" = "off" ]]; then
		_polap_log0 "  found1: ${_polap_var_ga_annotation_all}, so skipping the blast genome ..."
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	_polap_log1 "  creating annotation tables"
	_polap_log2 "    input1: ${_polap_var_ga_contigger_edges_stats}"
	_polap_log2 "    input2: ${_polap_var_ann_MTGENECOUNT}"
	_polap_log2 "    input3: ${_polap_var_ann_PTGENECOUNT}"
	_polap_log2 "    output1: ${_polap_var_ga_annotation}"
	_polap_log2 "    output2: ${_polap_var_ga_annotation_all}"
	_polap_log2 "    output3: ${_polap_var_ga_annotation_table}"
	_polap_log2 "    output4: ${_polap_var_ga_annotation_depth_table}"
	_polap_log2 "    output5: ${_polap_var_ga_pt_annotation_depth_table}"
	local _command1="Rscript ${_POLAPLIB_DIR}/run-polap-r-mtcontig.R \
		--flyeout-edges-stats ${_polap_var_ga_contigger_edges_stats} \
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

	ANUMNEXT=$((_arg_inum + 1))
	_polap_log1 NEXT: $(basename "$0") seeds -o "${_arg_outdir}" -i ${_arg_inum} -j $ANUMNEXT

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

polap_blast-genome() { # BLAST edge sequences on MT and PT genes
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"
	local MTAA="${_POLAPLIB_DIR}"/polap-mt.1.c70.3.faa
	local PTAA="${_POLAPLIB_DIR}"/polap-pt.2.c70.3.faa

	local _contigger=$(dirname "$1")
	local _ga=$(dirname "${_contigger}")
	local _contigger_edges_gfa="$1"
	local _contigger_edges_stats="$2"
	local _ann_MTGENECOUNT="$3"
	local _ann_PTGENECOUNT="$4"
	local _contigger_edges_fasta="${_contigger}/graph_final.fasta"

	local _ann="${_ga}/50-annotation"
	local _ann_CONTIGDB="${_ann}/contig"
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

	_polap_log1 "  NCBI BLAST edge contigs with mitochondrial and plastid genes ..."
	_polap_log2 "  input1: ${_contigger_edges_stats}"
	_polap_log2 "  input2: ${_contigger_edges_fasta}"

	# Checks the output files earlier than the input.
	if [[ -s "${_ann_MTGENECOUNT}" ]] &&
		[[ -s "${_ann_PTGENECOUNT}" ]] &&
		[[ "${_arg_redo}" = "off" ]]; then
		_polap_log2 "  found1: ${_ann_MTGENECOUNT}"
		_polap_log2 "  found2: ${_ann_PTGENECOUNT}"
		_polap_log2 "  so skipping the blast genome ..."
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	if [ ! -s "${_contigger_edges_stats}" ]; then
		_polap_log0 "ERROR: no edges_stats.txt file: ${_contigger_edges_stats}"
		exit $EXIT_ERROR
	fi

	_polap_log2 "  deleting and recreating the folder: ${_ann}"
	_polap_log3_cmd rm -rf "${_ann}"
	_polap_log3_cmd mkdir -p "${_ann}"

	#src/run-polap-select.R o/30-contigger/edges_stats.txt o/50-annotation/contig.name
	_polap_log2 "  contig sequence names in file: ${_ann_CONTIGNAME}"
	_polap_log2 "    input: ${_contigger_edges_stats}"
	_polap_log2 "    output: ${_ann_CONTIGFILE}"
	_polap_log3_pipe "grep -v '#' ${_contigger_edges_stats} |
		cut -f 1 >${_ann_CONTIGNAME}"

	_polap_log2 "  preparing edge contig sequence file: ${_ann_CONTIGFILE}"
	if [[ -s "${_contigger_edges_fasta}" ]]; then
		_polap_log3 "    copying edge contig sequence file"
		_polap_log3_cmd cp "${_contigger_edges_fasta}" "${_ann_CONTIGFILE}"
	else
		if [[ -s "${_contigger_edges_gfa}" ]]; then
			_polap_log3 "    extracting edge contig sequence from the gfa"
			_polap_log3_pipe "gfatools gfa2fa \
		    ${_contigger_edges_gfa} \
		    >${_contigger_edges_fasta}"
			_polap_log3_cmd cp "${_contigger_edges_fasta}" "${_ann_CONTIGFILE}"
		else
			_polap_log0 "ERROR: no such file: ${_contigger_edges_gfa}"
			return $RETURN_FAIL
		fi
	fi

	_polap_log2 "  making BLASTDB of the contig sequences: ${_ann_CONTIGDB}"
	_polap_log2 "    input: ${_ann_CONTIGFILE}"
	_polap_log2 "    output: ${_ann_CONTIGDB}"
	_polap_log3_pipe "makeblastdb -dbtype nucl \
		-in ${_ann_CONTIGFILE} \
		-out ${_ann_CONTIGDB} \
		2>${_polap_output_dest}"

	# Mitochondrial gene annotation and counts
	_polap_log2 "  BLAST of the mitochondrial proteins against ${_ann_CONTIGDB}"
	_polap_log2 "    input query: ${MTAA}"
	_polap_log2 "    input BLAST DB: ${_ann_CONTIGDB}"
	_polap_log2 "    evalue cutoff: 1e-30"
	_polap_log2 "    number of CPUs: ${_arg_threads}"
	_polap_log3 "    executing the tblastn ... be patient!"
	_polap_log3_pipe "tblastn -query ${MTAA} \
		-db ${_ann_CONTIGDB} \
		-out ${_ann_MTAABLAST} \
		-evalue 1e-30 \
		-outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle salltitles' \
		-num_threads ${_arg_threads} \
		>${_polap_output_dest} 2>&1"

	_polap_log2 "  converting BLAST result in BED format for removing redundancy"
	_polap_log2 "    input1: ${_ann_MTAABLAST}"
	_polap_log2 "    output1: ${_ann_MTAABLASTBED}"
	_polap_log3_pipe "Rscript ${_POLAPLIB_DIR}/run-polap-r-genes.R \
		${_ann_MTAABLAST} \
		${_ann_MTAABLASTBED} \
		>${_polap_output_dest} 2>&1"

	if [[ -s "${_ann_MTAABLASTBED}" ]]; then
		_polap_log2 "  sorting BED format BLAST result"
		_polap_log2 "    input1: ${_ann_MTAABLASTBED}"
		_polap_log2 "    output1: ${_ann_MTAABLAST}.sorted.bed"
		_polap_log3_pipe "sort -k1,1 -k2,2n \
		${_ann_MTAABLASTBED} \
		>${_ann_MTAABLAST}.sorted.bed"
	else
		>"${_ann_MTAABLAST}.sorted.bed"
	fi

	_polap_log2 "  creating folder: ${_ann_MTAABED}"
	_polap_log3_cmd mkdir "${_ann_MTAABED}"

	_polap_log2 "  counting mitochondrial genes in the contigs ..."
	_polap_log2 "    input1: ${_ann_MTAABLAST}.sorted.bed"
	_polap_log2 "    output1: ${_ann_MTGENECOUNT}"
	_polap_log3_cmd rm -f "${_ann_MTGENECOUNT_TMP}"
	while IFS= read -r contig; do
		_polap_log3_pipe "grep -w ${contig} ${_ann_MTAABLAST}.sorted.bed \
			>${_ann_MTAABED}/${contig}.bed"
		_polap_log3_pipe "bedtools merge -i ${_ann_MTAABED}/${contig}.bed \
			>${_ann_MTAABED}/${contig}.bed.txt"
		_polap_log3 commnad: printf \"%s\\t%d\\n\" "${contig}" $(wc -l <"${_ann_MTAABED}/${contig}".bed.txt) ">>${_ann_MTGENECOUNT_TMP}"
		printf "%s\t%d\n" ${contig} $(wc -l <${_ann_MTAABED}/${contig}.bed.txt) >>"${_ann_MTGENECOUNT_TMP}"
	done <"${_ann_CONTIGNAME}"
	_polap_log3_pipe "sort -k2 -rn ${_ann_MTGENECOUNT_TMP} >${_ann_MTGENECOUNT}"

	_polap_log2 "  compressing the BLAST results of mitochondrial gene annotation"
	_polap_log2 "    input: ${_ann_MTAABED}"
	_polap_log2 "    output: ${_ann_MTAABED}.tar.gz"
	_polap_log3_cmd tar zcf "${_ann_MTAABED}".tar.gz "${_ann_MTAABED}"
	_polap_log2 "  deleting folder: ${_ann_MTAABLASTBED}"
	_polap_log3_cmd rm -rf "${_ann_MTAABED}"

	_polap_log1 "  output1: ${_ann_MTGENECOUNT}"

	# Plastid gene annotation and counts
	_polap_log2 "  BLAST of the plastid proteins against ${_ann_CONTIGDB}"
	_polap_log2 "    input query: ${PTAA}"
	_polap_log2 "    input BLAST DB: ${_ann_CONTIGDB}"
	_polap_log2 "    evalue cutoff: 1e-30"
	_polap_log2 "    number of CPUs: ${_arg_threads}"
	_polap_log3 "    executing the tblastn ... be patient!"
	_polap_log3_pipe "tblastn -query ${PTAA} \
		-db ${_ann_CONTIGDB} \
		-out ${_ann_PTAABLAST} \
		-evalue 1e-30 \
		-outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle salltitles' \
		-num_threads ${_arg_threads} \
		>${_polap_output_dest} 2>&1"

	_polap_log2 "  converting BLAST result in BED format for removing redundancy"
	_polap_log2 "    input1: ${_ann_PTAABLAST}"
	_polap_log2 "    output1: ${_ann_PTAABLASTBED}"
	_polap_log3_pipe "Rscript ${_POLAPLIB_DIR}/run-polap-r-genes.R \
		${_ann_PTAABLAST} \
		${_ann_PTAABLASTBED} \
		>${_polap_output_dest} 2>&1"

	if [[ -s "${_ann_PTAABLASTBED}" ]]; then
		_polap_log2 "  sorting BED format BLAST result"
		_polap_log2 "    input1: ${_ann_PTAABLASTBED}"
		_polap_log2 "    output1: ${_ann_PTAABLAST}.sorted.bed"
		_polap_log3_pipe "sort -k1,1 -k2,2n \
		${_ann_PTAABLASTBED} \
		>${_ann_PTAABLAST}.sorted.bed"
	else
		>"${_ann_PTAABLAST}.sorted.bed"
	fi

	_polap_log2 "  creating folder: ${_ann_PTAABED}"
	_polap_log3_cmd mkdir "${_ann_PTAABED}"

	_polap_log2 "  counting mitochondrial genes in the contigs ..."
	_polap_log2 "    input1: ${_ann_PTAABLAST}.sorted.bed"
	_polap_log2 "    output1: ${_ann_PTGENECOUNT}"
	_polap_log3_cmd rm -f "${_ann_PTGENECOUNT_TMP}"
	while IFS= read -r contig; do
		_polap_log3_pipe "grep -w ${contig} ${_ann_PTAABLAST}.sorted.bed \
			>${_ann_PTAABED}/${contig}.bed"
		_polap_log3_pipe "bedtools merge -i ${_ann_PTAABED}/${contig}.bed \
			>${_ann_PTAABED}/${contig}.bed.txt"
		_polap_log3 commnad: printf \"%s\\t%d\\n\" "${contig}" $(wc -l <"${_ann_PTAABED}/${contig}".bed.txt) ">>${_ann_PTGENECOUNT_TMP}"
		printf "%s\t%d\n" ${contig} $(wc -l <${_ann_PTAABED}/${contig}.bed.txt) >>"${_ann_PTGENECOUNT_TMP}"
	done <"${_ann_CONTIGNAME}"
	_polap_log3_pipe "sort -k2 -rn ${_ann_PTGENECOUNT_TMP} >${_ann_PTGENECOUNT}"

	_polap_log2 "  compressing the BLAST results of mitochondrial gene annotation"
	_polap_log2 "    input: ${_ann_PTAABED}"
	_polap_log2 "    output: ${_ann_PTAABED}.tar.gz"
	_polap_log3_cmd tar zcf "${_ann_PTAABED}".tar.gz "${_ann_PTAABED}"
	_polap_log2 "  deleting folder: ${_ann_PTAABLASTBED}"
	_polap_log3_cmd rm -rf "${_ann_PTAABED}"

	_polap_log1 "  output2: ${_ann_PTGENECOUNT}"

	# Disable debugging if previously enabled
	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
# Blasts the genome assembly of a Flye run.
# Arguments:
#   -i ${_arg_inum}: a Flye genome assembly number
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
function _run_polap_blast-genome { # BLAST edge sequences on MT and PT genes
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
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
# NCBI BLAST the genome assembly of a Flye run againt the plant organelle genes.
#
# Arguments:
#   -i ${_arg_inum}: a Flye genome assembly number
# Inputs:
#   ${_polap_var_ga_contigger_edges_fasta}
#   ${_polap_var_ga_contigger_edges_gfa} if no such file: ${_polap_var_ga_contigger_edges_fasta}
#   ${_polap_var_ga_contigger_edges_stats}
# Outputs:
#   ${_polap_var_ann_MTGENECOUNT}
#   ${_polap_var_ann_PTGENECOUNT}
Example: $(basename "$0") ${_arg_menu[0]} -i ${_arg_inum}
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

	_polap_log0 "blasting edge contigs with mitochondrial and plastid genes on the assembly number ${_arg_inum} ..."
	_polap_log1 "  input1: ${_polap_var_ga_contigger_edges_stats}"
	_polap_log1 "  input2: ${_polap_var_ga_contigger_edges_fasta}"

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
	fi

	if [ ! -s "${_polap_var_ga_contigger_edges_stats}" ]; then
		_polap_log0 "ERROR: no edges_stats.txt file: ${_polap_var_ga_contigger_edges_stats}"
		exit $EXIT_ERROR
	fi

	_polap_log2 "  deleting and recreating the folder: ${_polap_var_ann}"
	_polap_log3_cmd rm -rf "${_polap_var_ann}"
	_polap_log3_cmd mkdir -p "${_polap_var_ann}"

	#src/run-polap-select.R o/30-contigger/edges_stats.txt o/50-annotation/contig.name
	_polap_log2 "  contig sequence names in file: ${_polap_var_ann_CONTIGNAME}"
	_polap_log2 "    input: ${_polap_var_ga_contigger_edges_stats}"
	_polap_log2 "    output: ${_polap_var_ann_CONTIGFILE}"
	_polap_log3_pipe "grep -v '#' ${_polap_var_ga_contigger_edges_stats} |
		cut -f 1 >${_polap_var_ann_CONTIGNAME}"

	_polap_log2 "  preparing edge contig sequence file: ${_polap_var_ann_CONTIGFILE}"
	if [[ -s "${_polap_var_ga_contigger_edges_fasta}" ]]; then
		_polap_log3 "    copying edge contig sequence file"
		_polap_log3_cmd cp "${_polap_var_ga_contigger_edges_fasta}" "${_polap_var_ann_CONTIGFILE}"
	else
		if [[ -s "${_polap_var_ga_contigger_edges_gfa}" ]]; then
			_polap_log3 "    extracting edge contig sequence from the gfa"
			_polap_log3_pipe "gfatools gfa2fa \
		    ${_polap_var_ga_contigger_edges_gfa} \
		    >${_polap_var_ga_contigger_edges_fasta}"
			_polap_log3_cmd cp "${_polap_var_ga_contigger_edges_fasta}" "${_polap_var_ann_CONTIGFILE}"
		else
			_polap_log0 "ERROR: no such file: ${_polap_var_ga_contigger_edges_gfa}"
			return $RETURN_FAIL
		fi
	fi

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
	_polap_log3_pipe "Rscript ${_POLAPLIB_DIR}/run-polap-r-genes.R \
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
	_polap_log3_pipe "Rscript ${_POLAPLIB_DIR}/run-polap-r-genes.R \
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
	# Rscript "${_POLAPLIB_DIR}"/run-polap-r-genes.R \
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

	_polap_log1 "NEXT (for testing purpose only): $(basename "$0") count-gene --test"
	_polap_log1 "NEXT: $(basename "$0") count-gene -o ${_arg_outdir} [-i ${_arg_inum}]"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}
