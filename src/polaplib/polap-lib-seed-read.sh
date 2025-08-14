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
# Convert numbers between different units.
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

_polap_lib_seed-plastid() {
	# we can use all polap_var_ variables.
	# They are determined by output, i, and j.
	source "${_POLAPLIB_DIR}/polap-variables-option.sh"
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	local start_contig=$(awk 'NR==2 {print $1}' "${_polap_var_ga_pt_annotation_depth_table}")
	# local start_contig=$(awk 'NR==2 {print $1}' "${_polap_var_ga_annotation_depth_table}")

	python "${_POLAPLIB_DIR}"/polap-py-find-mito-gfa.py \
		--seed "${start_contig}" \
		--mtcontig "${_polap_var_mtcontigname}" \
		"${_polap_var_ga_contigger_edges_gfa}"

	if [[ ! -s "${_polap_var_mtcontigname}" ]]; then
		echo "${start_contig}" >>"${_polap_var_mtcontigname}"
	fi
}

################################################################################
# Function to convert base pairs to the highest appropriate unit
# Example usage
# bp=31846726397
# convert_bp $bp
################################################################################
_polap_lib_seed-mito() {

	# we can use all polap_var_ variables.
	# They are determined by output, i, and j.
	source "${_POLAPLIB_DIR}/polap-variables-option.sh"
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

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
	_polap_log1 "mtcontigname set: ${mtcontignameset}"

	cp -p "${mtcontignameset}" "${_polap_var_mtcontigname}"
}

_polap_lib_seed-mito-2() {

	# we can use all polap_var_ variables.
	# They are determined by output, i, and j.
	source "${_POLAPLIB_DIR}/polap-variables-option.sh"
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	local _arg_seeds_scheme="5,6"
	_run_polap_seeds
	# local _polap_var_mtcontigname="${_polap_var_ga}/mt.contig.name-${_arg_jnum}"
	if [[ -s "${_polap_var_ga}/mt.contig.name-1" ]]; then
		_polap_log0_cat "${_polap_var_ga}/mt.contig.name-1"
	fi

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
	done < <(awk '{print $1}' "${_polap_var_ga}/mt.contig.name-1")
	# done < <(awk 'NR>1 {print $1}' "${_polap_var_ga_annotation_depth_table}")

	if [[ -s "${mtcontignameall}" ]]; then
		sort "${mtcontignameall}" | uniq >"${mtcontignameset}"
	else
		_polap_log0 "No mitochondrial contigs found in the gfa file."
		return 1
	fi
	_polap_log1 "mtcontigname set: ${mtcontignameset}"

	cp -p "${mtcontignameset}" "${_polap_var_mtcontigname}"
}
