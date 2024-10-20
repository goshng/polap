################################################################################
# Selects contigs for an organelle-genome assembly.
#
# 1. We could select mitochondrial- or plastid-derived contigs using a contig annotation table.
# 2. We determine the range of sequencing depths for those candidate contigs: mean +/- sd \* 3.
# 3. For a given gfa of a genome assembly graph, subset the graph for selecting graph elements in the range.
# 4. Determine connected components in the subset.
# 5. Choose connected components with candidate edges.
#
# We need to read GFA files to manipulate.
# We need to determine connected components.
################################################################################
function _run_polap_select-contigs-by-graph() {
	if [ "$DEBUG" -eq 1 ]; then set -x; fi

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge 3 ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	local FDIR="$ODIR"/$INUM
	local MTCONTIGNAME="$FDIR"/mt.contig.name-"$JNUM"
	local _polap_var_mtcontigs="$FDIR"/"$JNUM"/mtcontigs
	local _polap_var_assembly_graph_final_fasta="${FDIR}/30-contigger/graph_final.fasta"
	local _polap_var_assembly_graph_final_gfa="${FDIR}/30-contigger/graph_final.gfa"
	local _polap_var_annotation_table="${FDIR}/assembly_info_organelle_annotation_count-all.txt"
	local _polap_var_mtcontig_annotated="${_polap_var_mtcontigs}/1-mtcontig.annotated.txt"
	local _polap_var_mtcontigs_stats="${_polap_var_mtcontigs}/1-mtcontig.stats.txt"

	local _polap_var_gfa_all="${_polap_var_mtcontigs}/2-gfa.all.gfa"
	local _polap_var_gfa_seq_part="${_polap_var_mtcontigs}/2-gfa.seq.part.tsv"
	local _polap_var_gfa_seq_filtered="${_polap_var_mtcontigs}/2-gfa.seq.filtered.txt"
	local _polap_var_gfa_seq_filtered_range="${_polap_var_mtcontigs}/2-gfa.seq.filtered.range.txt"
	local _polap_var_gfa_seq_filtered_edge="${_polap_var_mtcontigs}/2-gfa.seq.filtered.edge.txt"
	local _polap_var_gfa_filtered="${_polap_var_mtcontigs}/2-gfa.filtered.gfa"
	local _polap_var_gfa_links="${_polap_var_mtcontigs}/3-gfa.links.tsv"
	local _polap_var_gfa_links_number="${_polap_var_mtcontigs}/3-gfa.links.number.txt"
	local _polap_var_gfa_links_order="${_polap_var_mtcontigs}/3-gfa.links.order.txt"
	local _polap_var_gfa_links_contig="${_polap_var_mtcontigs}/3-gfa.links.contig.txt"
	local _polap_var_gfa_links_contig_na="${_polap_var_mtcontigs}/3-gfa.links.contig.na.txt"
	local _polap_var_gfa_links_seed="${_polap_var_mtcontigs}/4-gfa.links.seed.txt"
	local _polap_var_gfa_links_mtcontig="${_polap_var_mtcontigs}/5-gfa.links.mtcontig.txt"

	# Print help message if requested
	help_message=$(
		cat <<HEREDOC
# Selects contigs using three features.
#
# To identify seed contigs of mitochondrial origin, 
# a whole-genome assembly is evaluated for three criteria: 
# 1) the presence of mitochondrial or plastid genes, 
# 2) the number of read coverage, and
# 3) the connectivity of contigs in the genome assembly graph. 
#
# 1. We could select mitochondrial- or plastid-derived contigs using a contig annotation table.
# 2. We determine the range of sequencing depths for those candidate contigs: mean +/- sd \* 3.
# 3. For a given gfa of a genome assembly graph, subset the graph for selecting graph elements in the range.
# 4. Determine connected components in the subset.
# 5. Choose connected components with candidate edges.
#
# Arguments:
#   -i $INUM: source Flye (usually whole-genome) assembly number
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

	[[ ${_arg_menu[1]} == "help" ]] && echo "${help_message}" >&2 && exit $EXIT_SUCCESS

	# Check for required files
	check_file_existence "${_polap_var_assembly_graph_final_gfa}"
	check_file_existence "${_polap_var_annotation_table}"

	# Clean and create working directory
	rm -rf "${_polap_var_mtcontigs}"
	mkdir -p "${_polap_var_mtcontigs}"

	_polap_log1_file "input1: ${_polap_var_assembly_graph_final_gfa}"
	_polap_log1_file "input2: ${_polap_var_annotation_table}"

	# Run R script to select mitochondrial contigs
	"$WDIR"/run-polap-select-contigs-1-start.R \
		"${_polap_var_assembly_graph_final_gfa}" \
		"${_polap_var_annotation_table}" \
		"${_polap_var_mtcontig_annotated}" \
		"${_polap_var_mtcontigs_stats}" \
		2>"$_polap_output_dest"

	_polap_log2_file "${_polap_var_mtcontig_annotated}"
	_polap_log2_file "${_polap_var_mtcontigs_stats}"

	# Handle case with single starting contig
	local mtcontig_count=$(wc -l <"${_polap_var_mtcontig_annotated}")
	if [ "${mtcontig_count}" -eq 1 ]; then
		cut -f1 "${_polap_var_mtcontig_annotated}" >"${MTCONTIGNAME}"
		echoerr "LOG: Single starting contig - done."
		_polap_log2_file "${MTCONTIGNAME}"
		return
	fi

	# Extract sequences and filter GFA data
	gfatools view -S "${_polap_var_assembly_graph_final_gfa}" \
		>"${_polap_var_gfa_all}" \
		2>"$_polap_output_dest"
	_polap_log2_file "${_polap_var_gfa_all}"

	gfatools view -S "${_polap_var_assembly_graph_final_gfa}" 2>/dev/null |
		grep "^S" >"${_polap_var_gfa_seq_part}"
	_polap_log2_file "${_polap_var_gfa_seq_part}"

	# Run R script to filter sequences
	"$WDIR"/run-polap-select-contigs-2-gfa-filter.R \
		"${_polap_var_gfa_seq_part}" \
		"${_polap_var_mtcontig_annotated}" \
		"${_polap_var_gfa_seq_filtered}" \
		"${_polap_var_gfa_seq_filtered_range}" \
		2>"$_polap_output_dest"

	_polap_log2_file "${_polap_var_gfa_seq_filtered}"

	# Filter GFA based on filtered sequences
	cut -f1 "${_polap_var_gfa_seq_filtered}" >"${_polap_var_gfa_seq_filtered_edge}"
	gfatools view -S \
		-l @"${_polap_var_gfa_seq_filtered_edge}" \
		"${_polap_var_assembly_graph_final_gfa}" 2>/dev/null \
		>"${_polap_var_gfa_filtered}"

	_polap_log2_file "${_polap_var_gfa_filtered}"

	# Generate GFA links
	grep "^L" "${_polap_var_gfa_filtered}" | cut -f2,4 >"${_polap_var_gfa_links}"
	_polap_log2_file "${_polap_var_gfa_links}"

	# Run R script to analyze GFA links
	"$WDIR"/run-polap-select-contigs-3-gfa-links.R \
		"${_polap_var_mtcontig_annotated}" \
		"${_polap_var_gfa_links}" \
		"${_polap_var_gfa_links_number}" \
		"${_polap_var_gfa_links_order}" \
		"${_polap_var_gfa_links_contig}" \
		"${_polap_var_gfa_links_contig_na}" \
		2>"$_polap_output_dest"

	_polap_log2_file "${_polap_var_gfa_links_number}"
	_polap_log2_file "${_polap_var_gfa_links_order}"
	_polap_log2_file "${_polap_var_gfa_links_contig}"
	_polap_log2_file "${_polap_var_gfa_links_contig_na}"

	# Find connected components using Python script
	python "$WDIR"/run-polap-select-contigs-4-find-connected-components.py \
		"${_polap_var_gfa_links_number}" \
		"${_polap_var_gfa_links_contig}" \
		"${_polap_var_gfa_links_seed}" \
		2>"$_polap_output_dest"

	_polap_log2_file "${_polap_var_gfa_links_seed}"

	# Choose final mitochondrial contigs
	"$WDIR"/run-polap-select-contigs-5-gfa-mtcontig.R \
		"${_polap_var_gfa_links_seed}" \
		"${_polap_var_gfa_links_order}" \
		"${_polap_var_gfa_links_mtcontig}" \
		2>"$_polap_output_dest"

	_polap_log2_file "${_polap_var_gfa_links_mtcontig}"

	cat "${_polap_var_gfa_links_mtcontig}" "${_polap_var_gfa_links_contig_na}" >"${MTCONTIGNAME}"
	_polap_log1_file "output: ${MTCONTIGNAME}"

	if [ "$DEBUG" -eq 1 ]; then set +x; fi
}
