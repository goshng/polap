################################################################################
# Select contig seeds from MT and PT
################################################################################
function _run_polap_select-contigs-by-organelle-gene-group() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge 3 ] && _polap_output_dest="/dev/stderr"

	# Set variables
	local FDIR="$ODIR/$INUM"
	local MTCONTIGNAME="$FDIR/mt.contig.name-$JNUM"
	local _polap_var_assembly_graph_final_gfa="$FDIR/30-contigger/graph_final.gfa"
	local _polap_var_annotation_table="$FDIR/assembly_info_organelle_annotation_count-all.txt"
	local _polap_var_mtcontigs="$FDIR/$JNUM/mtcontigs"
	local _polap_var_mt_contig_name="$FDIR/$JNUM/mt.contig.name"
	local _polap_var_mtcontig_annotated="${_polap_var_mtcontigs}/1-mtcontig.annotated.txt"
	local _polap_var_mtcontigs_stats="${_polap_var_mtcontigs}/1-mtcontig.stats.txt"

	# Help message
	local help_message=$(
		cat <<HEREDOC
# Selects contigs from organelle contigs.
#
# 1. Collect all edges with MT/PT genes.
# 2. Filter out edges with too long sequences.
#
# Arguments:
#   -i $INUM: source Flye (usually organelle-genome) assembly number
#   -j $JNUM: destination Flye organelle assembly number
#
# Inputs:
#   ${_polap_var_assembly_graph_final_gfa}
#   ${_polap_var_annotation_table}
#
# Outputs:
#   ${_polap_var_mtcontigs}
#   $MTCONTIGNAME
Example: $(basename $0) ${_arg_menu[0]} [-i|--inum <arg>] [-j|--jnum <arg>]
HEREDOC
	)

	_polap_log0 "1"
	# Print help message if needed
	[[ ${_arg_menu[1]} == "help" ]] && _polap_log1 "$help_message" && exit $EXIT_SUCCESS
	# [[ ${_arg_menu[1]} == "help" ]] && echoerr "$help_message" && exit $EXIT_SUCCESS

	# Check for required files
	check_file_existence "${_polap_var_assembly_graph_final_gfa}"
	check_file_existence "${_polap_var_annotation_table}"

	# Create the directory for mtcontigs
	rm -rf "${_polap_var_mtcontigs}"
	mkdir -p "${_polap_var_mtcontigs}"

	# Log the important files used in the process
	_polap_log2_file "input1: ${_polap_var_assembly_graph_final_gfa}"
	_polap_log2_file "input2: ${_polap_var_annotation_table}"

	# Step 1: Select contigs based on gene density
	"$WDIR"/run-polap-select-contigs-by-organelle-gene-group-1-start.R \
		"${_polap_var_assembly_graph_final_gfa}" \
		"${_polap_var_annotation_table}" \
		"${_polap_var_mtcontig_annotated}" \
		"${_polap_var_mtcontigs_stats}"

	_polap_log2_file "${_polap_var_mtcontig_annotated}"

	# Save the first column (contig names) to the output file
	cut -f1 "${_polap_var_mtcontig_annotated}" >"${MTCONTIGNAME}"
	_polap_log1_file "output: ${MTCONTIGNAME}"

	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}
