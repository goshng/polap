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

################################################################################
# input1: gfa
# output: output folder
################################################################################
function _polap_lib_pt-extract-dna {
	local gfa="${1}"
	local outdir="${2}"

	# avoid polap.log removal
	local ptdir="${2}/ptdna"
	local mtcontigname="${ptdir}/mt.contig.name"
	local _ptdna_extracted_fasta="${ptdir}/circular_path_1_concatenated.fa"
	local _ptdna_fasta="${outdir}/pt.0.fa"

	# annotate the gfa for contigs with PT genes
	#
	if [[ -d "${ptdir}" ]]; then
		_polap_log3_cmd rm -rf "${ptdir}"
	fi
	_polap_log3_cmd mkdir -p "${ptdir}/30-contigger"
	_polap_log3_cmd cp "${gfa}" "${ptdir}/30-contigger/graph_final.gfa"
	local _contigger_edges_gfa="${ptdir}/30-contigger/graph_final.gfa"
	local _ga_annotation_all="${ptdir}/assembly_info_organelle_annotation_count-all.txt"
	polap_annotate "${_contigger_edges_gfa}" "${_ga_annotation_all}"

	# prepare mt.contig.name
	# local _ga_pt_annotation_depth_table="${ptdir}/pt-contig-annotation-depth-table.txt"
	local _ga_pt_annotation_depth_table="${ptdir}/contig-annotation-depth-table.txt"
	awk 'NR==2 {print $1}' "${_ga_pt_annotation_depth_table}" >"${mtcontigname}"

	# extract ptdna
	_polap_log3_pipe "python \
          ${_POLAPLIB_DIR}/polap-py-find-plastid-gfa2fasta.py \
		        --gfa ${gfa} \
		        --seed ${mtcontigname} \
		        --out ${ptdir} \
		        2>$_polap_output_dest"

	if [[ -s "${ptdir}/circular_path_1_concatenated.fa" ]]; then
		cp "${_ptdna_extracted_fasta}" "${_ptdna_fasta}"
		_polap_log0 "output: ${_ptdna_fasta}"

		# NOTE: we could use this instead.
		#
		local mtcontigname2="${ptdir}/mt.contig.name2.txt"
		local start_contig=$(awk 'NR==2 {print $1}' "${_ga_pt_annotation_depth_table}")

		python "${_POLAPLIB_DIR}"/polap-py-find-plastid-gfa.py \
			--seed "${start_contig}" \
			--mtcontig "${mtcontigname2}" \
			--fasta "${outdir}/pt2.0.fa" \
			"${gfa}"
	else
		_polap_log0 "output: no ptDNA"
	fi
}
