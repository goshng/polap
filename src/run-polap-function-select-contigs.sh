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
[[ -n "${!_POLAP_INCLUDE_}" ]] && return 0
declare "$_POLAP_INCLUDE_=1"
#
################################################################################

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
function _run_polap_select-contigs-by() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "$script_dir/polap-variables-mtcontig.sh"

	# Print help message if requested
	help_message=$(
		cat <<HEREDOC
# Selects contigs using three features using total lengths of contigs.
#
# Use the total length in the cumulative distribution of contig lengths
# to select the lower and upper bounds of contig lengths.
#
# To identify seed contigs of mitochondrial origin, 
# a whole-genome assembly is evaluated for three criteria: 
# 1) the presence of mitochondrial or plastid genes, 
# 2) the number of read coverage, and
# 3) the connectivity of contigs in the genome assembly graph. 
#
# 1. We could select mitochondrial- or plastid-derived contigs using a contig annotation table.
# 2. We determine the range of sequencing depths for those candidate contigs: mean +/- sd \* 3.
#   2.1 Construct the cumulative distribution of contig lengths.
#   2.2 Given L1=3Mb, determine the lower bound of the contig length.
#   2.3 Given L2=300 kb, determine the upper bound of the contig length.
# 3. For a given gfa of a genome assembly graph, subset the graph for selecting graph elements in the range.
# 4. Determine connected components in the subset.
#   Choose connected components with candidate edges.
#
# More options:
# 1. find connected components without depth filtering
#
# --select-contig Types
# 1. gene density (or manual depth range)
# 2. gene density and depth range with cumulative copy number distribution
# 3. gene density and depth range with mixture distribution of copy numbers
# 4. gene depth-mixture density (or manual depth range)
# 5. gene depth-mixture density and depth range with cumulative copy number distribution 
# 6. gene depth-mixture density and depth range with mixture distribution of copy numbers
#
# candidate contigs preselection
# 1. gene density only
# 2. gene depth-mixture density
#
# depth range inference
# 1. depth range with cumulative copy number distribution 
# 2. depth range with mixture distribution of copy numbers
# 3. manual depth range
#
# Steps:
# 1. preselect candidate contigs:
#   1-preselection.by.gene.density.txt
#   1-depth.range.by.gene.density.txt
#   1-preselection.by.depth.mixture.txt
#   1-depth.range.by.depth.mixture.txt
#   1-mixfit.txt
# 2. infer depth range:
#   2-depth.range.by.cdf.copy.number.txt
#   2-depth.range.by.manual.selection.txt
# 3. filter GFA using the depth range from step 2
# 4. prepare the graph for finding connected components
# 5. find connected components by the depth-filtered contigs
# 6. convert the depth-filtered contigs in edge with numbers
# 7. (type 3, 4, 5) rearranging the output files: .table.tsv
#
# Arguments:
#   -i $INUM: source Flye (usually whole-genome) assembly number
#
# Arguments:
#   -i $INUM: source Flye (usually whole-genome) assembly number
#   -j $JNUM: destination Flye organelle assembly number
#   --select-contig [1|2|3]
# Inputs:
#   ${_polap_var_assembly_graph_final_gfa}
#   ${_polap_var_annotation_table}
#   ${ODIR}/0/1-mtcontig.depth.stats.txt for manual depth range in --select-contig 1
# Outputs:
#   $MTCONTIGNAME
#   "${_polap_var_mtcontig_annotated}"
# See:
#   manual: creates a custom depth range file.
#   delete: deletes the custom depth range file.
# See:
#   run-polap-select-contigs-by-1-table.R for the description of --select-contig option
Example: $(basename $0) ${_arg_menu[0]} [-i|--inum <arg>] [-j|--jnum <arg>] [--select-contig <number>]
Example: $(basename $0) ${_arg_menu[0]} -o PRJNA914763 -i 0 -j 5 --select-contig 5
Example: $(basename $0) scb -o PRJNA914763 -j 1
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && _polap_echo0 "${help_message}" && exit $EXIT_SUCCESS

	if [[ "${_arg_menu[1]}" == "manual" ]]; then
		printf "depth_lower_bound\tdepth_upper_bound\n" >"${_polap_var_manual_depth_range}"
		# Prompt the user to enter the first number
		_polap_log0 "Select the desired depth values for your seed contigs."
		_polap_log0 "  1. First, calculate an estimated average value accordingly."
		_polap_log0 "  2. Choosing a lower bound that is one-third of the average value can be a suitable option."
		_polap_log0 "  3. Selecting an appropriate upper limit can be achieved by setting it to three times the average."
		read -p "  Enter the lower bound of a depth range: " _num1
		read -p "  Enter the upper bound of a depth range: " _num2
		printf "%d\t%d\n" "${_num1}" "${_num2}" >>"${_polap_var_manual_depth_range}"
		_polap_log0 "  manual depth range: ${_num1} ~ ${_num2}"
		_polap_log3_cat "${_polap_var_manual_depth_range}"
		return
	fi

	if [[ "${_arg_menu[1]}" == "delete" ]]; then
		_polap_log0 "We remove the manual depth range configuration file from the output folder."
		_polap_log3_cmd rm -f "${_polap_var_manual_depth_range}"
		return
	fi

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		_polap_log0 ""
		_polap_log0 "View: ${JNUM}"
		if [[ -s "${_polap_var_mtcontig_table}" ]]; then
			_polap_log0_cat "${_polap_var_mtcontig_table}"
		fi
		if [[ -s "${_polap_var_depth_range_by_gene_density}" ]]; then
			_polap_log0_cat "${_polap_var_depth_range_by_gene_density}"
		fi
		if [[ -s "${_polap_var_depth_range_by_depth_mixture}" ]]; then
			_polap_log0_cat "${_polap_var_depth_range_by_depth_mixture}"
		fi
		if [[ -s "${_polap_var_mtcontig_depth_range}" ]]; then
			_polap_log0_cat "${_polap_var_mtcontig_depth_range}"
		fi
		if [[ -s "${_polap_var_mixfit}" ]]; then
			_polap_log0_cat "${_polap_var_mixfit}"
		fi
		if [[ -s "${MTCONTIGNAME}" ]]; then
			wc -l "${MTCONTIGNAME}" >&2
			_polap_log0_cat "${MTCONTIGNAME}"
			_polap_log0 "for Bandage nodes:"
			paste -sd ',' "${MTCONTIGNAME}" >&2
		fi
		_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		return
	fi

	_polap_log0 "selecting seed contigs using $(echo $FUNCNAME | sed s/_run_polap_//)"
	_polap_log0 "  ${INUM} -> ${JNUM} with type ${_arg_select_contig} ..."

	# Check for required files
	check_file_existence "${_polap_var_assembly_graph_final_gfa}"
	check_file_existence "${_polap_var_annotation_table}"
	_polap_log0 "  input1: ${_polap_var_assembly_graph_final_gfa}"
	_polap_log0 "  input2: ${_polap_var_annotation_table}"

	if [[ -s "${_polap_var_manual_copy_range}" ]]; then
		_polap_log0 "  input3: ${_polap_var_manual_copy_range}"

		Rscript "$script_dir"/run-polap-r-select-contigs-by-0-copy2depth.R \
			-t "${_polap_var_annotation_table}" \
			-c "${_polap_var_manual_copy_range}" \
			-o "${_polap_var_manual_depth_range}" \
			2>"$_polap_output_dest"

		numbers=$(grep -Eo '[0-9]+' "${_polap_var_manual_copy_range}")
		array=($numbers)
		_polap_log0 "  manual copy range: ${array[0]} ~ ${array[1]}"
		_polap_log3_cat "${_polap_var_manual_depth_range}"
	else
		_polap_log0 "  input3 not found: ${_polap_var_manual_copy_range}"
	fi

	# Clean and create working directory
	_polap_log1 "  delete and create the base mtcontigs folder:${_polap_var_mtcontigs}"
	_polap_log3_cmd rm -rf "${_polap_var_mtcontigs}"
	_polap_log3_cmd mkdir -p "${_polap_var_mtcontigs}"

	# Step 1: Determine the depth range using the cumulative length distribution.
	_polap_log1 "step 1: select contigs based on organelle gene annotation"
	_polap_log2 "  run-polap-select-contigs-by-1-annotation.R"
	_polap_log2 "    gene density for mtDNA: 10"
	_polap_log2 "    gene copy number range from mix-r: filter(m1$mu[1] < Copy, Copy < m1$mu[3] + m1$sd[3])"
	_polap_log2 "    select-contig type: ${_arg_select_contig}"
	_polap_log2 "    input: ${_polap_var_annotation_table}"
	_polap_log2 "    output-base: ${_polap_var_mtcontig_base}"
	case "${_arg_select_contig}" in
	1 | 2 | 3)
		Rscript "$script_dir"/run-polap-r-select-contigs-by-1-preselection.R \
			-t "${_polap_var_annotation_table}" \
			--out "${_polap_var_preselection_by_gene_density}" \
			--depth "${_polap_var_depth_range_by_gene_density}" \
			-c -d 10 \
			2>"$_polap_output_dest"
		local _preselection="${_polap_var_preselection_by_gene_density}"
		local _depth_range="${_polap_var_depth_range_by_gene_density}"
		;;
	4 | 5 | 6)
		Rscript "$script_dir"/run-polap-r-select-contigs-by-1-preselection.R \
			-t "${_polap_var_annotation_table}" \
			--out "${_polap_var_preselection_by_depth_mixture}" \
			--depth "${_polap_var_depth_range_by_depth_mixture}" \
			--mixfit "${_polap_var_mixfit}" \
			-c -d 10 \
			-r 2 \
			2>"$_polap_output_dest"
		local _preselection="${_polap_var_preselection_by_depth_mixture}"
		local _depth_range="${_polap_var_depth_range_by_depth_mixture}"
		;;
	*)
		_polap_log0 "ERROR: invalid case: ${_arg_select_contig}"
		;;
	esac

	if [[ "${_arg_select_contig}" = "3" ]]; then
		Rscript "$script_dir"/run-polap-r-select-contigs-by-1-preselection.R \
			-t "${_polap_var_annotation_table}" \
			--out "${_polap_var_preselection_by_depth_mixture}" \
			--depth "${_polap_var_depth_range_by_depth_mixture}" \
			--mixfit "${_polap_var_mixfit}" \
			-c -d 10 \
			-r 2 \
			2>"$_polap_output_dest"
		local _depth_range="${_polap_var_depth_range_by_depth_mixture}"
	fi

	if [[ -s "${_preselection}" ]]; then
		_polap_log3_file "${_preselection}"
		_polap_log3_cat "${_preselection}"
	fi
	if [[ -s "${_depth_range}" ]]; then
		_polap_log3_file "${_depth_range}"
		_polap_log3_cat "${_depth_range}"
	else
		_polap_log0 "  no such file: ${_depth_range}"
	fi
	if [[ -s "${_polap_var_mixfit}" ]]; then
		_polap_log3_file "${_polap_var_mixfit}"
		_polap_log3_cat "${_polap_var_mixfit}"
	fi

	# For all types, stop processing here if there are no seed contigs available.
	# This scenario would be extremely uncommon, occurring only under exceptional
	# circumstances.
	if [ ! -s "${_preselection}" ] || [ ! -s "${_depth_range}" ]; then
		>"${MTCONTIGNAME}"
		_polap_log0 "  output: ${MTCONTIGNAME} -> empty"
		_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		[ "$DEBUG" -eq 1 ] && set +x
		return
	fi

	# Stop processing for the single-seed contig scenario.
	# This scenario, although uncommon, is still a plausible occurrence.
	local mtcontig_count=$(wc -l <"${_preselection}")
	if [ "${mtcontig_count}" -eq 1 ]; then
		cut -f1 "${_preselection}" >"${MTCONTIGNAME}"
		_polap_log0 "  output: ${MTCONTIGNAME} -> a single contig"
		_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		[ "$DEBUG" -eq 1 ] && set +x
		return
	fi

	# Check if the user has provided a depth range file.
	if [[ -s "${_polap_var_manual_depth_range}" ]]; then
		_polap_log0 "  found: ${_polap_var_manual_depth_range}"
		# _polap_log3_cmd cp "${_polap_var_manual_depth_range}" "${_polap_var_mtcontig_depth_stats}"
	else
		_polap_log0 "  not found: ${_polap_var_manual_depth_range} -> no manual depth range provided"
	fi

	# Stop processing immediately if no manual depth range is specified by the user.
	# --select-contig 1 or 4
	case "${_arg_select_contig}" in
	1 | 4)
		if [[ -s "${_preselection}" ]] &&
			[[ ! -s "${_polap_var_manual_depth_range}" ]]; then

			cut -f1 "${_preselection}" |
				sort | uniq >"${MTCONTIGNAME}"
			_polap_log0 "  output1: ${MTCONTIGNAME}"
			_polap_log3_file "${MTCONTIGNAME}"
			_polap_log3_cat "${MTCONTIGNAME}"

			_polap_log2 "step 7: rearranging the output files: .table.tsv"
			_polap_log2 "    input1: ${_polap_var_annotation_table}"
			_polap_log2 "    input2: ${MTCONTIGNAME}"
			_polap_log2 "    output2: ${_polap_var_mtcontig_table}"
			Rscript "$script_dir"/run-polap-r-select-contigs-by-7-table.R \
				-t "${_polap_var_annotation_table}" \
				-m "${MTCONTIGNAME}" \
				-o "${_polap_var_mtcontig_table}" \
				2>"$_polap_output_dest"

			_polap_log0 "  output2: ${_polap_var_mtcontig_table}"

			_polap_log2 "Function end (${_arg_select_contig}): $(echo $FUNCNAME | sed s/_run_polap_//)"
			[ "$DEBUG" -eq 1 ] && set +x
			return
		fi
		;;
	*)
		_polap_log0 "  depth range is determined in --select-contig ${_arg_select_contig}"
		;;
	esac

	# Step 2. determine the depth range of organelle contigs.
	_polap_log1 "Step 2: determine the depth range of organelle contigs."
	_polap_log2 "  run-polap-select-contigs-by-2-determine-depth-range.R"
	_polap_log2 "    cumulative length cutoff: 3e+6"
	_polap_log2 "    select-contig type: ${_arg_select_contig}"
	_polap_log2 "    input: ${_polap_var_annotation_table}"
	_polap_log2 "    output-base: ${_polap_var_mtcontig_base}"
	_polap_log2 "    output1: ${_polap_var_mtcontig_depth_range}"
	case "${_arg_select_contig}" in
	1 | 4)
		# We utilize the depth range specified by the user.
		check_file_existence "${_polap_var_manual_depth_range}"
		_polap_log3_cmd cp "${_polap_var_manual_depth_range}" "${_polap_var_mtcontig_depth_range}"
		;;
	2 | 5)
		Rscript "$script_dir"/run-polap-r-select-contigs-by-2-determine-depth-range.R \
			-t "${_polap_var_annotation_table}" \
			-o "${_polap_var_depth_range_by_cdf_copy_number}" \
			2>"$_polap_output_dest"
		check_file_existence "${_polap_var_depth_range_by_cdf_copy_number}"
		_polap_log3_cmd cp "${_polap_var_depth_range_by_cdf_copy_number}" "${_polap_var_mtcontig_depth_range}"
		;;
	3 | 6)
		# We utilize either the depth range provided by the Mix-R or the annotated statistics obtained from Step 1.
		check_file_existence "${_polap_var_depth_range_by_depth_mixture}"
		_polap_log3_cmd cp "${_polap_var_depth_range_by_depth_mixture}" "${_polap_var_mtcontig_depth_range}"
		;;
	*)
		_polap_log0 "ERROR: invalid case: ${_arg_select_contig}"
		;;
	esac

	if [[ -s "${_polap_var_mtcontig_depth_range}" ]]; then
		_polap_log3_file "${_polap_var_mtcontig_depth_range}"
		_polap_log3_cat "${_polap_var_mtcontig_depth_range}"
	fi

	# Continue for the select contig types 3, 4, or 5
	_polap_log1 "step 3: filtering GFA using the depth range from step 2"
	_polap_log2 "  creating GFA without sequence data: ${_polap_var_gfa_all}"
	_polap_log2 "    input: ${_polap_var_assembly_graph_final_gfa}"
	_polap_log2 "    output: ${_polap_var_gfa_all}"
	gfatools view \
		-S "${_polap_var_assembly_graph_final_gfa}" \
		>"${_polap_var_gfa_all}" \
		2>"$_polap_output_dest"

	_polap_log2 "  extracting sequence part of GFA: ${_polap_var_gfa_seq_part}"
	_polap_log2 "    input: ${_polap_var_gfa_all}"
	_polap_log2 "    output: ${_polap_var_gfa_seq_part}"
	grep "^S" "${_polap_var_gfa_all}" >"${_polap_var_gfa_seq_part}"

	# Filter edges in GFA using depths.
	_polap_log2 "  filtering GFA sequence part using depth range"
	_polap_log2 "    input1: ${_polap_var_gfa_seq_part}"
	_polap_log2 "    input2: ${_polap_var_mtcontig_depth_range}"
	_polap_log2 "    output1: ${_polap_var_gfa_seq_filtered}"
	Rscript "$script_dir"/run-polap-r-select-contigs-by-3-subset-by-depth.R \
		--gfa "${_polap_var_gfa_seq_part}" \
		--depth "${_polap_var_mtcontig_depth_range}" \
		--out "${_polap_var_gfa_seq_filtered}" \
		2>"$_polap_output_dest"

	if [[ -s "${_polap_var_gfa_seq_filtered}" ]]; then
		_polap_log3_head "${_polap_var_gfa_seq_filtered}"
	fi

	# Recreate GFA based on filtered edge sequences.
	_polap_log1 "step 4: preparing the graph for finding connected components"
	_polap_log2 "  subsetting GFA using the depth-filtered GFA sequence part: ${_polap_var_gfa_filtered}"
	_polap_log2 "    input: ${_polap_var_gfa_seq_filtered}"
	_polap_log2 "    output: ${_polap_var_gfa_filtered}"
	cut -f1 "${_polap_var_gfa_seq_filtered}" >"${_polap_var_gfa_seq_filtered_edge}"
	gfatools view -S \
		-l @"${_polap_var_gfa_seq_filtered_edge}" \
		"${_polap_var_assembly_graph_final_gfa}" \
		2>"$_polap_output_dest" \
		>"${_polap_var_gfa_filtered}"

	if [[ -s "${_polap_var_gfa_filtered}" ]]; then
		_polap_log3_head "${_polap_var_gfa_filtered}"
	fi

	# Prepare links for finding connected components.
	_polap_log2 "  preparing links for finding connected components: ${_polap_var_gfa_links}"
	_polap_log2 "    input: ${_polap_var_gfa_filtered}"
	_polap_log2 "    output: ${_polap_var_links_tsv}"
	grep "^L" "${_polap_var_gfa_filtered}" | cut -f2,4 >"${_polap_var_links_tsv}"

	# No links case
	if [[ -s "${_polap_var_links_tsv}" ]]; then
		_polap_log3_head "${_polap_var_links_tsv}"
	else
		# no links -> just current selection
		cut -f2 "${_polap_var_gfa_filtered}" |
			sort | uniq >"${MTCONTIGNAME}"
		_polap_log0 "  output1: ${MTCONTIGNAME}"
		_polap_log3_file "${MTCONTIGNAME}"
		_polap_log3_cat "${MTCONTIGNAME}"

		_polap_log2 "step 7: rearranging the output files: .table.tsv"
		_polap_log2 "    input1: ${_polap_var_annotation_table}"
		_polap_log2 "    input2: ${MTCONTIGNAME}"
		_polap_log2 "    output2: ${_polap_var_mtcontig_table}"
		Rscript "$script_dir"/run-polap-r-select-contigs-by-7-table.R \
			-t "${_polap_var_annotation_table}" \
			-m "${MTCONTIGNAME}" \
			-o "${_polap_var_mtcontig_table}" \
			2>"$_polap_output_dest"

		_polap_log0 "  output2: ${_polap_var_mtcontig_table}"

		_polap_log2 "Function end (${_arg_select_contig}): $(echo $FUNCNAME | sed s/_run_polap_//)"
		[ "$DEBUG" -eq 1 ] && set +x
		return
	fi

	_polap_log3_file "${_polap_var_links_tsv}"
	_polap_log3_file "${_preselection}"

	# Run R script to analyze GFA links
	_polap_log2 "  preparing for finding connected components"
	_polap_log2 "    input1: ${_preselection}"
	_polap_log2 "    input2: ${_polap_var_links_tsv}"
	_polap_log2 "    output-base: ${_polap_var_links}"
	_polap_log2 "    output1: ${_polap_var_links_number}"
	_polap_log2 "    output2: ${_polap_var_links_order}"
	_polap_log2 "    output3: ${_polap_var_links_contig}"
	_polap_log2 "    output4: ${_polap_var_links_contig_na}"
	Rscript "$script_dir"/run-polap-r-select-contigs-by-4-prepare-for-connected-components.R \
		--edge "${_preselection}" \
		--gfa "${_polap_var_links_tsv}" \
		--out "${_polap_var_links}" \
		2>"$_polap_output_dest"

	if [[ -s "${_polap_var_links_number}" ]]; then
		_polap_log3_head "${_polap_var_links_number}"
	fi
	if [[ -s "${_polap_var_links_order}" ]]; then
		_polap_log3_head "${_polap_var_links_order}"
	fi
	if [[ -s "${_polap_var_links_contig}" ]]; then
		_polap_log3_head "${_polap_var_links_contig}"
	fi
	if [[ -s "${_polap_var_links_contig_na}" ]]; then
		_polap_log3_head "${_polap_var_links_contig_na}"
	fi

	# Find connected components using Python script
	_polap_log1 "step 5:"
	_polap_log2 "  finding connected components by the depth-filtered contigs"
	_polap_log2 "    input1: ${_polap_var_links_number}"
	_polap_log2 "    input2: ${_polap_var_links_contig}"
	_polap_log2 "    output1: ${_polap_var_links_seed}"
	python "$script_dir"/run-polap-select-contigs-by-5-find-connected-components.py \
		"${_polap_var_links_number}" \
		"${_polap_var_links_contig}" \
		"${_polap_var_links_seed}" \
		2>"$_polap_output_dest"

	if [[ -s "${_polap_var_links_seed}" ]]; then
		_polap_log3_head "${_polap_var_links_seed}"
	fi

	# Choose final mitochondrial contigs
	_polap_log1 "step 6:"
	_polap_log2 "  converting the depth-filtered contigs in edge with numbers"
	_polap_log2 "    input1: ${_polap_var_links_seed}"
	_polap_log2 "    input2: ${_polap_var_links_order}"
	_polap_log2 "    output: ${_polap_var_links_mtcontig}"
	Rscript "$script_dir"/run-polap-r-select-contigs-by-6-to-links-mtcontig.R \
		--seed "${_polap_var_links_seed}" \
		--order "${_polap_var_links_order}" \
		--out "${_polap_var_links_mtcontig}" \
		2>"$_polap_output_dest"

	if [[ -s "${_polap_var_links_mtcontig}" ]]; then
		_polap_log3_head "${_polap_var_links_mtcontig}"
	fi

	_polap_log2 "  concatenating the depth-filtered edges and NA edges: ${MTCONTIGNAME}"
	_polap_log2 "    input1: ${_polap_var_links_mtcontig}"
	_polap_log2 "    input2: ${_polap_var_links_contig_na}"
	_polap_log2 "    output1: ${MTCONTIGNAME}"
	cat "${_polap_var_links_mtcontig}" "${_polap_var_links_contig_na}" |
		sort | uniq >"${MTCONTIGNAME}"

	if [[ -s "${_polap_var_links_contig_na}" ]]; then
		_polap_log3_head "${_polap_var_links_contig_na}"
	fi

	_polap_log0 "  output1: ${MTCONTIGNAME}"

	if [[ -s "${MTCONTIGNAME}" ]]; then
		_polap_log3_cat "${MTCONTIGNAME}"
	fi

	_polap_log2 "step 7: (type 3, 4, 5) rearranging the output files: .table.tsv"
	_polap_log2 "    input1: ${_polap_var_annotation_table}"
	_polap_log2 "    input2: ${MTCONTIGNAME}"
	_polap_log2 "    output2: ${_polap_var_mtcontig_table}"
	Rscript "$script_dir"/run-polap-r-select-contigs-by-7-table.R \
		-t "${_polap_var_annotation_table}" \
		-m "${MTCONTIGNAME}" \
		-o "${_polap_var_mtcontig_table}" \
		2>"$_polap_output_dest"

	_polap_log0 "  output2: ${_polap_var_mtcontig_table}"

	if [[ -s "${_polap_var_mtcontig_table}" ]]; then
		_polap_log3_cat "${_polap_var_mtcontig_table}"
	fi

	_polap_log2 "Function end (${_arg_select_contig}): $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}

################################################################################
# Temporary
# scb: JNUM -> --select-contig
################################################################################
function _run_polap_scb() {
	# temporary
	_arg_select_contig="${JNUM}"
	_run_polap_select-contigs-by
}

################################################################################
# Select seed contigs using multiple methods
################################################################################
function _run_polap_select-contigs() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	help_message=$(
		cat <<HEREDOC
# Select seed contigs using multiple methods.
#
# Arguments:
#   -o $ODIR
#   -i $INUM: source Flye (usually whole-genome) assembly number
# Inputs:
#   $ODIR/$INUM
# Outputs:
#   $ODIR/$INUM/1
#   $ODIR/$INUM/2
#   $ODIR/$INUM/3
#   $ODIR/$INUM/4
#   $ODIR/$INUM/5
Example: $(basename $0) ${_arg_menu[0]} [-o $ODIR] [-i <number>]
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && _polap_echo0 "${help_message}" && exit $EXIT_SUCCESS
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		if [[ "${_arg_menu[2]}" == "more" ]]; then
			for i in "${_arg_select_contig_numbers[@]}"; do
				local MTCONTIGNAME="${ODIR}/${INUM}/mt.contig.name-${i}"

				if [ -s "$MTCONTIGNAME" ]; then
					_arg_select_contig="${i}"
					JNUM="${i}"
					_run_polap_select-contigs-by
				fi
			done
		else
			wc "${ODIR}/${INUM}"/mt.contig.name-* >&2
			exit $EXIT_SUCCESS
		fi

		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		exit $EXIT_SUCCESS
	fi

	_polap_log0 "selecting contigs from the genome assembly number ${INUM} ..."

	# Loop over numbers from 1 to 5
	for i in "${_arg_select_contig_numbers[@]}"; do
		# Call the function corresponding to the current number (index is i-1)
		if [[ "${INUM}" -ne 0 ]]; then
			_polap_log0 "ASSERT: only -i 0 is implemented"
			exit $EXIT_FAIL
		fi
		local MTCONTIGNAME="${ODIR}/${INUM}/mt.contig.name-${i}"

		if [ -e "$MTCONTIGNAME" ] && [ "${_arg_redo}" = "off" ]; then
			_polap_log0 "  found: ${MTCONTIGNAME}, so skipping the select-contig step for ${MTCONTIGNAME}"
		else
			_arg_select_contig="${i}"
			JNUM="${i}"
			_run_polap_select-contigs-by
		fi

		# check the mt.contig.name-1
		if [ -s "$MTCONTIGNAME" ]; then
			_polap_log1_file "${MTCONTIGNAME}"
		else
			_polap_log0 "  $MTCONTIGNAME is empty; choose seed contigs by yourself."
		fi

	done

	if [ "${_arg_verbose}" -ge 0 ]; then
		wc "${ODIR}/${INUM}"/mt.contig.name-* >&2
	fi

	_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}

function _run_polap_sc() {
	_run_polap_select-contigs
}
