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
# Subcommands for mtDNA or ptDNA sequences.
# 1. compare two sequences by a pairwise alignment.
# 2. fetch mtDNA or ptDNA sequences from the NCBI.
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

function _run_polap_mtdna-edges {
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Set variables
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	# Help message
	local help_message=$(
		cat <<HEREDOC
# Selects mtDNA sequences from a GFA.
#
# Arguments:
#   -j ${_arg_jnum}: organelle assembly number
#
# Inputs:
#   ${_polap_var_oga_assembly_graph_gfa}
#
# Outputs:
#   ${_polap_var_oga_mt_fasta}
#
Example: $(basename $0) ${_arg_menu[0]} -j <arg>
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	# Create necessary directories
	rm -rf "${_polap_var_mtdna}"
	mkdir -p "${_polap_var_mtdna}"

	# Prompt the user for input
	read -p "Enter edges using Bandage (e.g., edge_265, edge_520, edge_425): " edges

	# Convert the input string into new lines (replace ", " with newlines)
	local formatted_edges=$(echo "$edges" | tr ', ' '\n' | sed '/^ *$/d')
	echo "${formatted_edges}" >"${_polap_var_mtdna_names}"

	_polap_log1_file "Input GFA: ${_polap_var_oga_assembly_graph_gfa}"
	# _polap_log1_file "Input Annotation Table: ${_polap_var_ga_annotation_all}"
	#
	_polap_log3_pipe "gfatools view \
		-l @${_polap_var_mtdna_names} \
		${_polap_var_oga_assembly_graph_gfa} \
		2>$_polap_output_dest \
		>${_polap_var_mtdna_sub_gfa}"

	_polap_log0 "subgfa: ${_polap_var_mtdna_sub_gfa}"

	# Convert GFA to FASTA
	_polap_log3_pipe "gfatools gfa2fa ${_polap_var_mtdna_sub_gfa} \
		>${_polap_var_mtdna_fasta} \
		2>$_polap_output_dest"

	cp ${_polap_var_mtdna_fasta} ${_polap_var_oga_mt_fasta}

	_polap_log0 "FASTA of the input GFA: ${_polap_var_oga_mt_fasta}"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
# Selects mtDNA sequences from a GFA.
#
# annotation: uses 30-contigger/graph_final.gfa
#
# extraction: uses assembly_graph.gfa
################################################################################
function _run_polap_select-mtdna {
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Set variables
	# CHECK: local function
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	# local FDIR="${_arg_outdir}/${_arg_inum}"
	# local _polap_var_mtdna="$FDIR/mtdna"
	#
	# # File paths
	# local _polap_var_assembly_graph_gfa="$FDIR/assembly_graph.gfa"
	# local _polap_var_ga_annotation_all="$FDIR/assembly_info_organelle_annotation_count-all.txt"
	# local _polap_var_mt_fasta="$FDIR/mt.0.fasta"
	# local _polap_var_mt_edges="$FDIR/mt.0.edges"
	# local _polap_var_mtdna_1_gfa_all="${_polap_var_mtdna}/1-gfa.all.gfa"
	# local _polap_var_mtdna_gfa_links="${_polap_var_mtdna}/1-gfa.links.tsv"
	# local _polap_var_mtdna_gfa_links_edges="${_polap_var_mtdna}/1-gfa.links.edges.txt"
	# local _polap_var_mtdna_gfa_links_circular_path="${_polap_var_mtdna}/2-gfa.links.circular.path.txt"
	# local _polap_var_mtdna_circular_path="${_polap_var_mtdna}/3-circular.path.txt"
	# local _polap_var_mtdna_edge_fasta="${_polap_var_mtdna}/5-edge.fasta"
	# local _polap_var_mtdna_fasta="${_polap_var_mtdna}/4-gfa.fasta"

	# Help message
	local help_message=$(
		cat <<HEREDOC
# Selects mtDNA sequences from a GFA.
#
# Arguments:
#   -j ${_arg_jnum}: organelle assembly number
#
# Inputs:
#   ${_polap_var_oga_assembly_graph_gfa}
#
# Outputs:
#   ${_polap_var_oga_mt_fasta}
#
Example: $(basename $0) ${_arg_menu[0]} -j <arg>
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	if [[ "${_arg_menu[1]}" = "bandage" ]]; then
		# Prompt the user for input
		read -p "Enter edges using Bandage (e.g., edge_265, edge_520, edge_425): " edges

		# Convert the input string into new lines (replace ", " with newlines)
		local formatted_edges=$(echo "$edges" | tr ', ' '\n' | sed '/^ *$/d')
		echo "${formatted_edges}" >"${_polap_var_mtdna_names}"

		return $RETURN_SUCCESS
	fi

	# Check for required files
	# check_file_existence "${_polap_var_assembly_graph_gfa}"
	# check_file_existence "${_polap_var_ga_annotation_all}"

	# Create necessary directories
	# _polap_log3_cmd rm -rf "${_polap_var_mtdna}"
	# _polap_log3_cmd mkdir -p "${_polap_var_mtdna}"

	_polap_log1_file "Input GFA: ${_polap_var_oga_assembly_graph_gfa}"
	# _polap_log1_file "Input Annotation Table: ${_polap_var_ga_annotation_all}"
	#
	_polap_log3_pipe "gfatools view \
		-l @${_polap_var_mtdna_names} \
		${_polap_var_oga_assembly_graph_gfa} \
		2>$_polap_output_dest \
		>${_polap_var_mtdna_sub_gfa}"

	_polap_log0 "subgfa: ${_polap_var_mtdna_sub_gfa}"

	# Convert GFA to FASTA
	_polap_log3_pipe "gfatools gfa2fa ${_polap_var_mtdna_sub_gfa} \
		>${_polap_var_mtdna_fasta} \
		2>$_polap_output_dest"

	_polap_log2_file "FASTA of the input GFA: ${_polap_var_mtdna_fasta}"

	return

	# Step 1: List connected components (GFA processing)
	_polap_log3_pipe "gfatools view -S ${_polap_var_assembly_graph_gfa} \
		>${_polap_var_mtdna_1_gfa_all} \
		2>$_polap_output_dest"

	_polap_log2_file "GFA content: ${_polap_var_mtdna_1_gfa_all}"

	# Extract links from GFA and store them
	grep "^L" "${_polap_var_mtdna_1_gfa_all}" >"${_polap_var_mtdna_gfa_links}"

	local _polap_var_number_links_gfa=$(wc -l <"${_polap_var_mtdna_gfa_links}")
	if [ "${_polap_var_number_links_gfa}" -eq 0 ]; then
		echoerr "LOG: No circular sequences are found."
		return
	fi

	_polap_log2_file "GFA Links: ${_polap_var_mtdna_gfa_links}"

	# Process the GFA links
	"${_POLAPLIB_DIR}"/run-polap-r-select-mtdna-1-nx-gfa-links.R \
		"${_polap_var_mtdna_gfa_links}" \
		"${_polap_var_mtdna_gfa_links_edges}" \
		2>"$_polap_output_dest"

	_polap_log2_file "GFA Links Edges: ${_polap_var_mtdna_gfa_links_edges}"

	# Step 2: Find circular paths using Python script
	python "${_POLAPLIB_DIR}"/run-polap-py-select-mtdna-2-nx-simple-cycles.py \
		"${_polap_var_mtdna_gfa_links_edges}" \
		"${_polap_var_mtdna_gfa_links_circular_path}" \
		2>"$_polap_output_dest"
	# python "${_POLAPLIB_DIR}"/run-polap-py-select-mtdna-2-nx-find-circular-path.py \
	# 	"${_polap_var_mtdna_gfa_links_edges}" \
	# 	"${_polap_var_mtdna_gfa_links_circular_path}" \
	# 	2>"$_polap_output_dest"

	_polap_log2_file "Circular Path: ${_polap_var_mtdna_gfa_links_circular_path}"

	# return
	#
	# Process circular path
	# process_circular_path "${_polap_var_mtdna_gfa_links_circular_path}" "${_polap_var_mtdna_circular_path}"
	# _polap_log2_file "Processed Circular Path: ${_polap_var_mtdna_circular_path}"

	# Step 3: Concatenate sequences from circular path
	for _polap_file in "${_polap_var_mtdna_gfa_links_circular_path}"_*.tsv; do
		# Extract the number from the filename using parameter expansion
		local index=$(echo "${_polap_file}" | sed -E 's/.*_(.*)\.tsv/\1/')
		# _polap_log0 "${index}"
		concatenate_sequences "${_polap_file}" "${_polap_var_mtdna_fasta}" "${_polap_var_mtdna_edge_fasta}-${index}"
		finalize_mtdna_fasta "${_polap_var_mtdna_edge_fasta}-${index}" "${_polap_var_mt_fasta}-${index}.fasta"
		_polap_log2_file "Final mtDNA FASTA: ${_polap_var_mt_fasta}-${index}.fasta"
		cp "${_polap_file}" "${_polap_var_mt_edges}-${index}.edges"
		seqkit stats "${_polap_var_mt_fasta}-${index}.fasta" 1>&2
	done

	largest_file=$(ls -S "${_polap_var_mt_fasta}"-*.fasta | head -n 1)
	local index=$(echo "${largest_file}" | sed -E 's/.*-(.*)\.fasta/\1/')
	cp "${largest_file}" "${_polap_var_mt_fasta}"
	cp "${_polap_var_mt_edges}-${index}.edges" "${_polap_var_mt_edges}"

	# concatenate_sequences "${_polap_var_mtdna_circular_path}" "${_polap_var_mtdna_fasta}" "${_polap_var_mtdna_edge_fasta}"

	# Finalize the mtDNA FASTA
	# finalize_mtdna_fasta "${_polap_var_mtdna_edge_fasta}" "${_polap_var_mt_fasta}"
	# _polap_log2_file "Final mtDNA FASTA: ${_polap_var_mt_fasta}"

	# cp "${_polap_var_mtdna_circular_path}" "${_polap_var_mt_edges}"

	_polap_log1_file "Output mtDNA FASTA: ${_polap_var_mt_fasta}"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
# Before 2024-09-27
################################################################################
function _run_polap_select-mtdna-org {
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Set variables
	# CHECK: local function
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	# local FDIR="${_arg_outdir}/${_arg_inum}"
	# local _polap_var_mtdna="$FDIR/mtdna"
	#
	# # File paths
	# local _polap_var_assembly_graph_gfa="$FDIR/assembly_graph.gfa"
	# local _polap_var_ga_annotation_all="$FDIR/assembly_info_organelle_annotation_count-all.txt"
	# local _polap_var_mt_fasta="$FDIR/mt.0.fasta"
	# local _polap_var_mt_edges="$FDIR/mt.0.edges"
	# local _polap_var_mtdna_1_gfa_all="${_polap_var_mtdna}/1-gfa.all.gfa"
	# local _polap_var_mtdna_gfa_links="${_polap_var_mtdna}/1-gfa.links.tsv"
	# local _polap_var_mtdna_gfa_links_edges="${_polap_var_mtdna}/1-gfa.links.edges.txt"
	# local _polap_var_mtdna_gfa_links_circular_path="${_polap_var_mtdna}/2-gfa.links.circular.path.txt"
	# local _polap_var_mtdna_circular_path="${_polap_var_mtdna}/3-circular.path.txt"
	# local _polap_var_mtdna_edge_fasta="${_polap_var_mtdna}/5-edge.fasta"
	# local _polap_var_mtdna_fasta="${_polap_var_mtdna}/4-gfa.fasta"

	# Help message
	local help_message=$(
		cat <<HEREDOC
# Selects mtDNA sequences from a GFA.
#
# Arguments:
#   -i ${_arg_inum}: organelle assembly number
#
# Inputs:
#   ${_polap_var_assembly_graph_gfa}
#   ${_polap_var_ga_annotation_all}
#
# Outputs:
#   ${_polap_var_mt_fasta}
#
Example: $(basename $0) ${_arg_menu[0]} [-i|--inum <arg>]
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && echo "${help_message}" >&2 && exit $EXIT_SUCCESS

	# Check for required files
	check_file_existence "${_polap_var_assembly_graph_gfa}"
	check_file_existence "${_polap_var_ga_annotation_all}"

	# Create necessary directories
	rm -rf "${_polap_var_mtdna}"
	mkdir -p "${_polap_var_mtdna}"

	_polap_log1_file "Input GFA: ${_polap_var_assembly_graph_gfa}"
	_polap_log1_file "Input Annotation Table: ${_polap_var_ga_annotation_all}"

	# Convert GFA to FASTA
	gfatools gfa2fa "${_polap_var_assembly_graph_gfa}" \
		>"${_polap_var_mtdna_fasta}" \
		2>"$_polap_output_dest"

	_polap_log2_file "FASTA of the input GFA: ${_polap_var_mtdna_fasta}"

	# Step 1: List connected components (GFA processing)
	gfatools view -S "${_polap_var_assembly_graph_gfa}" \
		>"${_polap_var_mtdna_1_gfa_all}" \
		2>"$_polap_output_dest"

	_polap_log2_file "GFA content: ${_polap_var_mtdna_1_gfa_all}"

	# Extract links from GFA and store them
	grep "^L" "${_polap_var_mtdna_1_gfa_all}" >"${_polap_var_mtdna_gfa_links}"

	local _polap_var_number_links_gfa=$(wc -l <"${_polap_var_mtdna_gfa_links}")
	if [ "${_polap_var_number_links_gfa}" -eq 0 ]; then
		echoerr "LOG: No circular sequences are found."
		return
	fi

	_polap_log2_file "GFA Links: ${_polap_var_mtdna_gfa_links}"

	# Process the GFA links
	"${_POLAPLIB_DIR}"/run-polap-r-select-mtdna-1-nx-gfa-links.R \
		"${_polap_var_mtdna_gfa_links}" \
		"${_polap_var_mtdna_gfa_links_edges}" \
		2>"$_polap_output_dest"

	_polap_log2_file "GFA Links Edges: ${_polap_var_mtdna_gfa_links_edges}"

	# Step 2: Find circular paths using Python script
	# python "${_POLAPLIB_DIR}"/run-polap-py-select-mtdna-2-nx-simple-cycles.py \
	# 	"${_polap_var_mtdna_gfa_links_edges}" \
	# 	"${_polap_var_mtdna_gfa_links_circular_path}" \
	# 	2>"$_polap_output_dest"
	python "${_POLAPLIB_DIR}"/run-polap-py-select-mtdna-2-nx-find-circular-path.py \
		"${_polap_var_mtdna_gfa_links_edges}" \
		"${_polap_var_mtdna_gfa_links_circular_path}" \
		2>"$_polap_output_dest"

	_polap_log2_file "Circular Path: ${_polap_var_mtdna_gfa_links_circular_path}"

	# Process circular path
	process_circular_path "${_polap_var_mtdna_gfa_links_circular_path}" "${_polap_var_mtdna_circular_path}"
	_polap_log2_file "Processed Circular Path: ${_polap_var_mtdna_circular_path}"

	# Step 3: Concatenate sequences from circular path
	concatenate_sequences "${_polap_var_mtdna_circular_path}" "${_polap_var_mtdna_fasta}" "${_polap_var_mtdna_edge_fasta}"

	# Finalize the mtDNA FASTA
	finalize_mtdna_fasta "${_polap_var_mtdna_edge_fasta}" "${_polap_var_mt_fasta}"
	_polap_log2_file "Final mtDNA FASTA: ${_polap_var_mt_fasta}"

	cp "${_polap_var_mtdna_circular_path}" "${_polap_var_mt_edges}"

	_polap_log1_file "Output mtDNA FASTA: ${_polap_var_mt_fasta}"
	seqkit stats "${_polap_var_mt_fasta}" 1>&2

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

# Helper function to process the circular path and format the edge data
process_circular_path() {
	local input_file=$1
	local output_file=$2

	tail -n +2 "${input_file}" |
		while IFS=$'\t' read -r col1 col2; do
			local number1=$(echo "$col1" | grep -o '^[0-9]*')
			local sign1=$(echo "$col1" | grep -o '[+-]$')
			local new_col1="edge_${number1}\t${sign1}"
			echo -e "$new_col1"
		done >"${output_file}"
}

# Helper function to concatenate sequences based on circular path
concatenate_sequences() {
	local circular_path_file=$1
	local fasta_file=$2
	local output_fasta=$3

	# Create an empty file for the concatenated sequence
	>"${output_fasta}"

	# Loop through the circular path IDs and concatenate sequences
	while read -r id strand; do
		echoerr "Processing ID: $id with strand: $strand"
		if [[ "$strand" == "+" ]]; then
			seqkit grep -p "$id" "${fasta_file}" | seqkit seq -t dna -v >>"${output_fasta}"
		elif [[ "$strand" == "-" ]]; then
			seqkit grep -p "$id" "${fasta_file}" | seqkit seq -t dna -v -r -p >>"${output_fasta}"
		fi
	done <"${circular_path_file}"
}

# Helper function to finalize mtDNA FASTA by concatenating the sequences into one entry
finalize_mtdna_fasta() {
	local input_fasta=$1
	local output_fasta=$2

	echo ">concatenated_sequence" >"${output_fasta}"
	seqkit fx2tab "${input_fasta}" | cut -f2 | tr -d '\n' >>"${output_fasta}"
}