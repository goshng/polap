function _run_polap_select-contigs-by-gene-density() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x

	# Set variables
	local FDIR="$ODIR/$INUM"
	local MTCONTIGNAME="$FDIR/mt.contig.name-$JNUM"
	local _polap_var_mtcontigs="$FDIR/$JNUM/mtcontigs"

	local _polap_var_assembly_graph_final_fasta="$FDIR/30-contigger/graph_final.fasta"
	local _polap_var_assembly_graph_final_gfa="$FDIR/30-contigger/graph_final.gfa"
	local _polap_var_annotation_table="$FDIR/assembly_info_organelle_annotation_count-all.txt"

	local _polap_var_mtcontig_annotated="${_polap_var_mtcontigs}/1-mtcontig.annotated.txt"
	local _polap_var_mtcontigs_stats="${_polap_var_mtcontigs}/1-mtcontig.stats.txt"

	# Help message
	local help_message=$(
		cat <<HEREDOC
# Selects contigs for an organelle-genome assembly based on gene density.
#
# 1. Collect all edges with gene density above a threshold.
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
#   $MTCONTIGNAME
Example: $(basename $0) ${_arg_menu[0]} [-i|--inum <arg>] [-j|--jnum <arg>]
HEREDOC
	)

	# Print help message if needed
	[[ ${_arg_menu[1]} == "help" ]] && echoerr "$help_message" && exit $EXIT_SUCCESS

	# Check for required files
	check_file_existence "${_polap_var_assembly_graph_final_gfa}"
	check_file_existence "${_polap_var_annotation_table}"

	# Create the directory for mtcontigs
	rm -rf "${_polap_var_mtcontigs}"
	mkdir -p "${_polap_var_mtcontigs}"

	# Log the important files used in the process
	_polap_log1_file "input1: ${_polap_var_assembly_graph_final_gfa}"
	_polap_log1_file "input2: ${_polap_var_annotation_table}"

	# Step 1: Select contigs based on gene density
	"$WDIR"/run-polap-select-contigs-by-gene-density-1-start.R \
		"${_polap_var_assembly_graph_final_gfa}" \
		"${_polap_var_annotation_table}" \
		"${_polap_var_mtcontig_annotated}" \
		"${_polap_var_mtcontigs_stats}"

	_polap_log2_flie "${_polap_var_mtcontig_annotated}"

	# Save the first column (contig names) to the output file
	cut -f1 "${_polap_var_mtcontig_annotated}" >"${MTCONTIGNAME}"
	_polap_log1_file "output: ${MTCONTIGNAME}"

	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}
