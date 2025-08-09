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
# Example usage
# bp=31846726397
# convert_bp $bp
################################################################################
_polap_lib_assemble-rate() {

	_polap_log0 "Assemble using seed contigs by adjusting data sampling"

	# we can use all polap_var_ variables.
	# They are determined by output, i, and j.
	source "${_POLAPLIB_DIR}/polap-variables-option.sh"
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	if ! _polap_gfatools-gfa2fasta; then
		_polap_error_message $?
	fi

	check_file_existence "${_polap_var_mtcontigname}"
	check_file_existence "${_polap_var_ga_contigger_edges_fasta}"

	_run_polap_map-reads

	_polap_lib_oga-estimate-read-sampling-rate

	# remove NUMT/NUPT using rkmerrc
	# if [[ "${_arg_data_type}" == "pacbio-hifi" ]] && [[ "${_arg_plastid}" == "off" ]]; then
	if [[ "${_arg_data_type}" == "pacbio-hifi" ]]; then
		_polap_log0 "  use TIPPo's rmkc method to remove unusual PacBio HiFi reads before assembly"
		local _pread_sel="ptgaul-reads"
		local index=$(<"${_polap_var_oga_contig}/index.txt")
		local fq="${_polap_var_oga_seeds}/${_pread_sel}/${index}.fq"
		gunzip "${fq}.gz"
		local PREFIX="${_arg_outdir}/kmer/rmkc"
		local CLEANED="$PREFIX.cleaned.fastq.gz"
		_polap_log0 "rmkc on ${fq}"
		_polap_filter-reads-by-rmkc "${fq}"
		_polap_log0 "rmkc produces ${CLEANED}"
		_polap_lib_oga-flye-select-reads "${CLEANED}"
	else
		_polap_lib_oga-flye-select-reads
	fi

	return 0
}

# assemble after adjusting the read data downsampling rate so that
# the subsampling rate is between 0.1 and 0.5.
#
_polap_lib_assemble-omega() {

	_polap_log0 "Assemble using seed contigs by adjusting omega"

	# we can use all polap_var_ variables.
	# They are determined by output, i, and j.
	source "${_POLAPLIB_DIR}/polap-variables-option.sh"
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	if ! _polap_gfatools-gfa2fasta; then
		_polap_error_message $?
	fi

	check_file_existence "${_polap_var_mtcontigname}"
	check_file_existence "${_polap_var_ga_contigger_edges_fasta}"

	_run_polap_map-reads

	_polap_lib_oga-estimate-omega

	# remove NUMT/NUPT using rkmerrc
	# if [[ "${_arg_data_type}" == "pacbio-hifi" ]] && [[ "${_arg_plastid}" == "off" ]]; then
	if [[ "${_arg_data_type}" == "pacbio-hifi" ]]; then
		_polap_log0 "  use TIPPo's rmkc method to remove unusual PacBio HiFi reads before assembly"
		local _pread_sel="ptgaul-reads"
		local index=$(<"${_polap_var_oga_contig}/index.txt")
		local fq="${_polap_var_oga_seeds}/${_pread_sel}/${index}.fq"
		gunzip "${fq}.gz"
		local PREFIX="${_arg_outdir}/kmer/rmkc"
		local CLEANED="$PREFIX.cleaned.fastq.gz"
		_polap_log0 "rmkc on ${fq}"
		_polap_filter-reads-by-rmkc "${fq}"
		_polap_log0 "rmkc produces ${CLEANED}"
		_polap_lib_oga-flye-select-reads "${CLEANED}"
	else
		# _polap_lib_oga-flye-select-reads
		_arg_redo="on"
		_run_polap_assemble2
	fi

	return 0
}
