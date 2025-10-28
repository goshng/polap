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
# Functions for subcommand template ...
# Describe what they are and what they do.
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

function _run_polap_seed-plastid-fixed-form {
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

	help_message=$(
		cat <<'EOF'
Name:
  polap seed-plastid-fixed-form - select plastid seed of LSC-IR-SSC or one ring

Synopsis:
  polap seed-plastid-fixed-form [options]

Description:
  polap seed-plastid uses a gfa to select plastid contigs.
  Annotate the gfa with organelle genes to select one starting plastid contig.
  Select contigs that are connected by tha starting plastid contig.

Options:
  -o DIR
    Polap output directory

  -i STR
    Polap source assembly name

  -j STR
    Polap destination assembly name

Examples:
  Get organelle genome sequences:
    polap annotate-read -l l.fq

TODO:
  Dev.

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	# Display help message
	if [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]]; then
		local manfile=$(_polap_lib_man-convert_help_message "$help_message" "${_arg_menu[0]}")
		man "$manfile" >&3
		rm -f "$manfile"
		return
	fi

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
		exit $EXIT_SUCCESS
	fi

	# reannotate

	# _polap_log0 "${_polap_var_mtcontigname}"
	# _polap_log0 "${_polap_var_ga_pt_annotation_depth_table}"
	# _polap_log0 "${_polap_var_wga_contigger_edges_gfa}"
	# _polap_log0 "${_polap_var_ga_contigger_edges_gfa}"
	# _polap_log0 "${_polap_var_oga_contigger_edges_gfa}"

	local start_contig=$(awk 'NR==2 {print $1}' "${_polap_var_ga_pt_annotation_depth_table}")

	python "${_POLAPLIB_DIR}"/polap-py-find-plastid-gfa.py \
		--seed "${start_contig}" \
		--mtcontig "${_polap_var_mtcontigname}" \
		--fasta "${_polap_var_ga_contigger}/ptdna.fa" \
		"${_polap_var_ga_contigger_edges_gfa}"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_polap_seed-plastid {
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

	help_message=$(
		cat <<'EOF'
Name:
  polap seed-plastid - select plastid seed

Synopsis:
  polap seed-plastid [options]

Description:
  polap seed-plastid uses a gfa to select plastid contigs.
  Annotate the gfa with organelle genes to select one starting plastid contig.
  Select contigs that are connected by tha starting plastid contig.

Options:
  -o DIR
    Polap output directory

  -i STR
    Polap source assembly name

  -j STR
    Polap destination assembly name

Examples:
  Get organelle genome sequences:
    polap annotate-read -l l.fq

TODO:
  Dev.

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	# Display help message
	if [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]]; then
		local manfile=$(_polap_lib_man-convert_help_message "$help_message" "${_arg_menu[0]}")
		man "$manfile" >&3
		rm -f "$manfile"
		return
	fi

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
		exit $EXIT_SUCCESS
	fi

	# reannotate

	# _polap_log0 "${_polap_var_mtcontigname}"
	# _polap_log0 "${_polap_var_ga_pt_annotation_depth_table}"
	# _polap_log0 "${_polap_var_wga_contigger_edges_gfa}"
	# _polap_log0 "${_polap_var_ga_contigger_edges_gfa}"
	# _polap_log0 "${_polap_var_oga_contigger_edges_gfa}"

	local start_contig=$(awk 'NR==2 {print $1}' "${_polap_var_ga_pt_annotation_depth_table}")

	python "${_POLAPLIB_DIR}"/polap-py-find-mito-gfa.py \
		--seed "${start_contig}" \
		--mtcontig "${_polap_var_mtcontigname}" \
		"${_polap_var_ga_contigger_edges_gfa}"

	if [[ ! -s "${_polap_var_mtcontigname}" ]]; then
		echo "${start_contig}" >>"${_polap_var_mtcontigname}"
	fi

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_polap_seed-mito {
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

	help_message=$(
		cat <<'EOF'
Name:
  polap seed-mito - select mitochondrial seed

Synopsis:
  polap seed-plastid [options]

Description:
  polap seed-mito uses a gfa to select mitochondrial contigs.
  Annotate the gfa with organelle genes to select starting mitochondrial contigs.
  Select contigs that are connected by tha starting plastid contigs.

Options:
  -o DIR
    Polap output directory

  -i STR
    Polap source assembly name

  -j STR
    Polap destination assembly name

Examples:
  Get organelle genome sequences:
    polap annotate-read -l l.fq

TODO:
  Dev.

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	# Display help message
	if [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]]; then
		local manfile=$(_polap_lib_man-convert_help_message "$help_message" "${_arg_menu[0]}")
		man "$manfile" >&3
		rm -f "$manfile"
		return
	fi

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
		exit $EXIT_SUCCESS
	fi

	# _polap_log0 "${_polap_var_mtcontigname}"
	# _polap_log0 "${_polap_var_ga_pt_annotation_depth_table}"
	# _polap_log0 "${_polap_var_wga_contigger_edges_gfa}"
	# _polap_log0 "${_polap_var_ga_contigger_edges_gfa}"
	# _polap_log0 "${_polap_var_oga_contigger_edges_gfa}"

	# local start_contig=$(awk 'NR>1 {print $1}' "${_polap_var_ga_annotation_depth_table}")
	# _polap_log0 "${start_contig}"

	mkdir -p "${_polap_var_ga_contig}"
	local mtcontigname="${_polap_var_ga_contig}/mt.contig.name.txt"
	local mtcontignameall="${_polap_var_ga_contig}/mt.contig.name.all.txt"
	local mtcontignameset="${_polap_var_ga_contig}/mt.contig.name.set.txt"
	rm -f "${mtcontignameall}"
	while read -r start_contig; do
		rm -f "${mtcontigname}"
		_polap_log3_cmd python "${_POLAPLIB_DIR}"/polap-py-find-mito-gfa.py \
			--seed "${start_contig}" \
			--mtcontig "${mtcontigname}" \
			"${_polap_var_ga_contigger_edges_gfa}"
		if [[ -s "${mtcontigname}" ]]; then
			cat "${mtcontigname}" >>"${mtcontignameall}"
		else
			echo "${start_contig}" >>"${mtcontignameall}"
		fi
	done < <(awk 'NR>1 {print $1}' "${_polap_var_ga_annotation_depth_table}")
	if [[ -s "${mtcontignameall}" ]]; then
		sort "${mtcontignameall}" | uniq >"${mtcontignameset}"
	else
		_polap_log0 "No mitochondrial contigs found in the gfa file."
		return 1
	fi
	_polap_log0 "mtcontigname set: ${mtcontignameset}"

	# pt
	# local ptcontigname="${_polap_var_oga_contig}/pt.contig.name.txt"
	# local ptcontignameall="${_polap_var_oga_contig}/pt.contig.name.all.txt"
	# >"${ptcontignameall}"
	# while read -r start_contig; do
	# 	python "${_POLAPLIB_DIR}"/polap-py-find-mito-gfa.py \
	# 		--seed "${start_contig}" \
	# 		--mtcontig "${ptcontigname}" \
	# 		"${_polap_var_ga_contigger_edges_gfa}"
	# 	cat "${ptcontigname}" >>"${ptcontignameall}"
	# done < <(awk 'NR==2 {print $1}' "${_polap_var_ga_pt_annotation_depth_table}")
	# sort "${ptcontignameall}" | uniq >"${_polap_var_ptcontigname}"
	# _polap_log0 "ptcontigname: ${_polap_var_ptcontigname}"
	#
	# # mt - pt contig name
	# grep -Fxv -f "${_polap_var_ptcontigname}" "${mtcontignameset}" >"${_polap_var_mtcontigname}"
	#
	# _polap_log0 "mtcontigname (mtcontignameset - ptcontigname): ${_polap_var_mtcontigname}"

	cp -p "${mtcontignameset}" "${_polap_var_mtcontigname}"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}
