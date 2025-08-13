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

# Estimate the genome size from a short-read dataset
#
# input1: a short-read fastq data
# input2: an output folder _outdir
# input3: redo default is on
# output: the genome size
# _outdir_genome_size="${_outdir}/short_expected_genome_size.txt"
_polap_lib_wga-find-genome-size() {
	local _short_read1="$1"
	local _outdir="$2"
	local _redo="${3:-on}"
	local _rstatus
	local _outdir_genome_size
	local _outdir_jellyfish_out
	local _outdir_jellyfish_out_histo
	local unzipped_file
	local _EXPECTED_GENOME_SIZE
	local _expected_genome_size_bp

	# Set paths for bioproject data
	# source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'
	# source "${_POLAPLIB_DIR}/run-polap-function-utilities.sh"

	_polap_log2 "    estimate the genome size using the short-read data ${_short_read1} ..."

	if [[ -d "${_outdir}" ]]; then
		_polap_log3 "    output folder: ${_outdir}"
	else
		_polap_log3 mkdir -p "${_outdir}"
		mkdir -p "${_outdir}"
	fi
	check_file_existence "${_short_read1}"

	_polap_log3 "    input1: ${_short_read1}"
	_outdir_genome_size="${_outdir}/short_expected_genome_size.txt"
	_outdir_jellyfish_out="${_outdir}/jellyfish_out"
	_outdir_jellyfish_out_histo="${_outdir}/jellyfish_out.histo"

	# See https://bioinformatics.uconn.edu/genome-size-estimation-tutorial/
	if [ -s "${_outdir_genome_size}" ] && [ "${_redo}" = "off" ]; then
		_polap_log3 "    found: ${_outdir_genome_size}, so skipping the genome size estimation ..."
		_polap_log3_file "${_outdir_genome_size}"
	else
		if [ -s "${_outdir_jellyfish_out}" ] && [ "${_redo}" = "off" ]; then
			_polap_log3 "    found: ${_outdir_jellyfish_out}, so skipping the JellyFish counting ..."
			_polap_log3_file "${_outdir_jellyfish_out}"
		else
			unzipped_file=$(_polap_lib_file-gunzip "${_short_read1}")
			# local comp_size=$(stat -c%s "${_short_read1}")
			_rstatus="$?"
			if [[ "$_rstatus" -eq 0 ]]; then
				_polap_log3 "    unzipped file: $unzipped_file"
				_short_read1="$unzipped_file"
			fi
			# File sizes in bytes
			# local orig_size=$(stat -c%s "${unzipped_file}")

			# Compute compression ratio: original / compressed
			# local ratio=$(echo "scale=2; $orig_size / $comp_size" | bc)
			# echo "${ratio}" >"${_outdir}/ratio-l-fq-gz.txt"

			if [[ -s "${_short_read1}" ]]; then
				_polap_log3_pipe "command time -v jellyfish count \
					-t ${_arg_threads} \
					-C -m 21 \
					-s ${_arg_jellyfish_s} \
					-o ${_outdir_jellyfish_out} \
          ${_short_read1} 2>${_outdir}/timing-jellyfish.txt"
			else
				die "ASSERT: we must have at least one short-read fastq file."
			fi
			check_file_existence "${_outdir_jellyfish_out}"
		fi

		check_file_existence "${_outdir_jellyfish_out}"
		_polap_log3_pipe "jellyfish histo \
			-t ${_arg_threads} \
			-o ${_outdir_jellyfish_out_histo} \
      ${_outdir_jellyfish_out} 2>${_outdir}/timing-jellyfish-histo.txt"
		_polap_log3_file "${_outdir_jellyfish_out_histo}"

		if [[ "${_arg_data_type}" == "pacbio-hifi" ]]; then
			_polap_log3_pipe "genomescope2 \
					-i ${_outdir_jellyfish_out_histo} \
    		  -k 21 \
					-o ${_outdir_jellyfish_out}.dir \
					>${_outdir_genome_size}.genomescope2.txt"
			# Check the exit status
			if [ $? -ne 0 ]; then
				# Take action if needed, e.g., logging, sending a notification, etc.
				_polap_log0 "ERROR: disk space might be full"
				die "ERROR: not being to estimate the genome size using short-read data: short-read may be too small."
			fi
			cut -d: -f6 "${_outdir_genome_size}.genomescope2.txt" |
				tail -n 1 >"${_outdir_genome_size}"
		else
			_polap_log3_pipe "Rscript --vanilla ${_POLAPLIB_DIR}/polap-r-jellyfish.R \
					${_outdir_jellyfish_out_histo} \
					${_outdir_genome_size}"
			# Check the exit status
			if [ $? -ne 0 ]; then
				# Take action if needed, e.g., logging, sending a notification, etc.
				_polap_log0 "ERROR: disk space might be full"
				die "ERROR: not being to estimate the genome size using short-read data: short-read may be too small."
			fi
		fi
	fi
	_polap_log3 "    output: ${_outdir_genome_size}"
	_polap_log3_cat "${_outdir_genome_size}"

	local _EXPECTED_GENOME_SIZE=$(<"${_outdir_genome_size}")
	local _EXPECTED_GENOME_SIZE=${_EXPECTED_GENOME_SIZE%.*}
	local _expected_genome_size_bp=$(_polap_utility_convert_bp ${_EXPECTED_GENOME_SIZE})
	_polap_log3 "    expected genome size using short-read data (bases): ${_expected_genome_size_bp}"

	_polap_log3_cmd rm -f "${_outdir_jellyfish_out}"
	_polap_log3_cmd rm -f "${_outdir_jellyfish_out_histo}"
	_polap_log3_cmd rm -rf "${_outdir_jellyfish_out}.dir"

	return 0
}
