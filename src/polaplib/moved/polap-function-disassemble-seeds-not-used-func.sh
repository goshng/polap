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
# This script has functions for ptDNA contig seed selection.
################################################################################

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  echo "[ERROR] This script must be sourced, not executed: use 'source $BASH_SOURCE'" >&2
  return 1 2>/dev/null || exit 1
fi
: "${_POLAP_DEBUG:=0}"
: "${_POLAP_RELEASE:=0}"

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
# Finalize the mt.contig.name.
################################################################################
# TODO
# Called by:
# nothing
_polap_disassemble_seeds_final-mtcontig() {
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	local _ga_annotation_all="$1"
	local _ga_annotation_depth_table="$2"
	local _mtcontigs_7mtcontigname="$3"
	local _mtcontig_table="$4"
	local _mtcontigs_annotation_table_seed="$5"

	# local _is_auto=$1 not using it anymore? for --length 1000000 option below

	_polap_log1 "  step 8-1: filtering by max length 1 Mb of contig seeds"
	_polap_log1 "    1 Mb length filtering for auto"
	_polap_log2 "    input1: ${_ga_annotation_all}"
	_polap_log2 "    input2: ${_mtcontigs_7mtcontigname}"
	_polap_log2 "    output1: ${_mtcontig_table}"
	local _command1="Rscript ${_POLAPLIB_DIR}/polap-r-final-filter-mtcontig.R \
		-t ${_ga_annotation_all} \
		-m ${_mtcontigs_7mtcontigname} \
    -o ${_mtcontig_table} \
    --length 1000000 \
		2>$_polap_output_dest"
	_polap_log3_pipe "${_command1}"

	_polap_disassemble_seeds_final-seeds-mtcontig \
		"${_ga_annotation_all}" \
		"${_ga_annotation_depth_table}" \
		"${_mtcontigs_7mtcontigname}" \
		"${_mtcontigs_annotation_table_seed}"
}

# TODO
# Called by:
# _polap_disassemble_seeds_final-mtcontig() {
_polap_disassemble_seeds_final-seeds-mtcontig() {
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	local _ga_annotation_all="$1"
	local _ga_annotation_depth_table="$2"
	local _mt_contig_name="$3"
	local _annotation_table_seed_target="$4"

	_polap_log1 "  annotation table with contig seed marks"
	_polap_log2 "    input1: ${_ga_annotation_all}"
	_polap_log2 "    input2: ${_ga_annotation_depth_table}"
	_polap_log2 "    input3: ${_mt_contig_name}"
	_polap_log2 "    output1: ${_annotation_table_seed_target}"
	local _command1="Rscript ${_POLAPLIB_DIR}/polap-r-final-seed-mtcontig.R \
		-t ${_ga_annotation_all} \
		-a ${_ga_annotation_depth_table} \
		-m ${_mt_contig_name} \
    -o ${_annotation_table_seed_target} \
		2>$_polap_output_dest"
	_polap_log3_pipe "${_command1}"

}

################################################################################
# Prepare data for connected components
################################################################################
_polap_disassemble_seeds_prepare-cc() {
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	local _mtcontigs=$(dirname "$1")
	local _mtcontigs_preselection="$1"
	local _mtcontigs_gfa_seq_filtered="$2"
	local _mtcontigs_links_number="$3"
	local _mtcontigs_links_order="$4"
	local _mtcontigs_links_contig="$5"
	local _mtcontigs_links_contig_na="$6"
	local _mtcontigname="$7"

	local _mtcontigs_gfa_seq_filtered_edge="${_mtcontigs}/4-gfa.seq.depthfiltered.edge.txt"
	local _mtcontigs_gfa_depthfiltered_gfa="${_mtcontigs}/4-gfa.depthfiltered.gfa"
	local _mtcontigs_links="${_mtcontigs}/4-gfa.links"
	local _mtcontigs_links_tsv="${_mtcontigs}/4-gfa.links.tsv"

	_polap_log2 "    step 4-1: subsetting GFA using the depth-filtered GFA sequence part: ${_mtcontigs_gfa_depthfiltered_gfa}"
	_polap_log2 "      input1: ${_mtcontigs_gfa_seq_filtered}"
	_polap_log2 "      output1: ${_mtcontigs_gfa_seq_filtered_edge}"
	_polap_log2 "      output2: ${_mtcontigs_gfa_depthfiltered_gfa}"

	_polap_log3_pipe "cut -f1 \
    ${_mtcontigs_gfa_seq_filtered} \
    >${_mtcontigs_gfa_seq_filtered_edge}"

	_polap_log3_pipe "gfatools view -S \
		-l @${_mtcontigs_gfa_seq_filtered_edge} \
		${_ga_contigger_edges_gfa} \
		2>$_polap_output_dest \
		>${_mtcontigs_gfa_depthfiltered_gfa}"

	if [[ -s "${_mtcontigs_gfa_depthfiltered_gfa}" ]]; then
		_polap_log2_head "${_mtcontigs_gfa_depthfiltered_gfa}"
	fi

	# Prepare links for finding connected components.
	_polap_log2 "    step 4-2: preparing links for finding connected components"
	_polap_log2 "      input: ${_mtcontigs_gfa_depthfiltered_gfa}"
	_polap_log2 "      output: ${_mtcontigs_links_tsv}"
	_polap_log3_pipe "grep ^L ${_mtcontigs_gfa_depthfiltered_gfa} |\
    cut -f2,4 >${_mtcontigs_links_tsv}"

	# No links case
	if [[ -s "${_mtcontigs_links_tsv}" ]]; then
		_polap_log2_head "${_mtcontigs_links_tsv}"
	else
		# no links -> just current selection
		_polap_log3_pipe "cut -f2 ${_mtcontigs_gfa_depthfiltered_gfa} |
			sort | uniq >${_mtcontigname}"
		_polap_log0 "  output1: ${_mtcontigname}"
		_polap_log2_cat "${_mtcontigname}"

		# if [[ "${_arg_plastid}" = "off" ]]; then
		# 	_polap_disassemble_seeds_final-mtcontig \
		# 		"${_ga_annotation_all}" \
		# 		"${_ga_annotation_depth_table}" \
		# 		"${_mtcontigs_7mtcontigname}" \
		# 		"${_mtcontig_table}" \
		# 		"${_mtcontigs_annotation_table_seed}"
		# else
		# 	_polap_disassemble_seeds_final-mtcontig \
		# 		"${_ga_annotation_all}" \
		# 		"${_ga_pt_annotation_depth_table}" \
		# 		"${_mtcontigs_7mtcontigname}" \
		# 		"${_mtcontig_table}" \
		# 		"${_mtcontigs_annotation_table_seed}"
		# fi

		_polap_log2 "Function end (${_arg_select_contig}): $(echo $FUNCNAME | sed s/_run_polap_//)"
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
	fi

	local _preselection="${_mtcontigs_preselection}"

	# Run R script to analyze GFA links
	_polap_log2 "    step 4-3: preparing for finding connected components"
	_polap_log2 "      input1: ${_preselection}"
	_polap_log2 "      input2: ${_mtcontigs_links_tsv}"
	_polap_log2 "      output-base: ${_mtcontigs_links}"
	_polap_log2 "      output1: ${_mtcontigs_links_number}"
	_polap_log2 "      output2: ${_mtcontigs_links_order}"
	_polap_log2 "      output3: ${_mtcontigs_links_contig}"
	_polap_log2 "      output4: ${_mtcontigs_links_contig_na}"
	_polap_log3_pipe "Rscript ${_POLAPLIB_DIR}/polap-r-prepare-cc.R \
		--edge ${_preselection} \
		--gfa ${_mtcontigs_links_tsv} \
		--out ${_mtcontigs_links} \
		--unit-annotation-edge \
		2>$_polap_output_dest"

	if [[ -s "${_mtcontigs_links_number}" ]]; then
		_polap_log2_head "${_mtcontigs_links_number}"
	fi
	if [[ -s "${_mtcontigs_links_order}" ]]; then
		_polap_log2_head "${_mtcontigs_links_order}"
	fi
	if [[ -s "${_mtcontigs_links_contig}" ]]; then
		_polap_log2_head "${_mtcontigs_links_contig}"
	fi
	if [[ -s "${_mtcontigs_links_contig_na}" ]]; then
		_polap_log2_head "${_mtcontigs_links_contig_na}"
	else
		_polap_log2 "No NA edge: ${_mtcontigs_links_contig_na}"
	fi
}
