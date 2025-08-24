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
# Function to convert base pairs to the highest appropriate unit
# input1: gfa
# output: output folder
################################################################################
function _polap_lib_mt-extract-dna {
	local gfa="${1}"
	local outdir="${2}"

	# avoid polap.log removal
	local mtdir="${2}/mtdna"
	local mtcontigname="${mtdir}/mt.contig.name"
	local _mtdna_extracted_fasta="${mtdir}/circular_path_1_concatenated.fa"
	local _mtdna_fasta="${outdir}/mt.0.fa"

	# annotate the gfa for contigs with mt genes
	#
	if [[ -d "${mtdir}" ]]; then
		_polap_log3_cmd rm -rf "${mtdir}"
	fi
	_polap_log3_cmd mkdir -p "${mtdir}/30-contigger"
	_polap_log3_cmd cp "${gfa}" "${mtdir}/30-contigger/graph_final.gfa"
	local _contigger_edges_gfa="${mtdir}/30-contigger/graph_final.gfa"
	local _ga_annotation_all="${mtdir}/assembly_info_organelle_annotation_count-all.txt"
	_polap_lib_annotate \
		-o "${outdir}" \
		-i "mtdna"
	# polap_annotate "${_contigger_edges_gfa}" "${_ga_annotation_all}"

	# prepare mt.contig.name
	local _ga_mt_annotation_depth_table="${mtdir}/contig-annotation-depth-table.txt"
	awk 'NR==2 {print $1}' "${_ga_mt_annotation_depth_table}" >"${mtcontigname}"

	_polap_log0_dev "TODO: implement mtDNA sequence extraction to output a FASTA file: ${_mtdna_fasta}"

	# assemble2 based on the mtcontigname

	# extract mtdna
	# _polap_log3_pipe "python \
	#          ${_POLAPLIB_DIR}/polap-py-find-plastid-gfa2fasta.py \
	# 	        --gfa ${gfa} \
	# 	        --seed ${mtcontigname} \
	# 	        --out ${mtdir} \
	# 	        2>$_polap_output_dest"
	#
	# if [[ -s "${mtdir}/circular_path_1_concatenated.fa" ]]; then
	# 	cp "${_mtdna_extracted_fasta}" "${_mtdna_fasta}"
	# 	_polap_log0 "output: ${_mtdna_fasta}"
	# else
	# 	_polap_log0 "output: no mtDNA"
	# fi
}
