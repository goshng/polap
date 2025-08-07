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
# Do something about FASTQ files.
#
# Function:
# function _run_polap_fastq {
#
# TODO: document the help message.
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

function _run_polap_fastq {
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

	# Print help message if requested
	help_message=$(
		cat <<HEREDOC
Deal with FASTQ sequence files.

Menu:
view
sample
subsample
subsample2

Example: 
$(basename $0) ${_arg_menu[0]} view infile.fq
$(basename $0) ${_arg_menu[0]} sample infile.fq outfile.fq <rate:0.1>
$(basename $0) ${_arg_menu[0]} subsample infile.fq outfile.fq.gz -c 10 -o o
$(basename $0) ${_arg_menu[0]} subsample2 s1.fq s2.fq s1.fq.gz s2.fq.gz -c 10 -o o
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then

		local _infile="${_arg_menu[2]}"

		seqkit stats "${_infile}" >&3

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
	fi

	# sample infile.fq outfile.fq.gz 0.1
	if [[ "${_arg_menu[1]}" == "sample" ]]; then

		local _infile="${_arg_menu[2]}"
		local _outfile="${_arg_menu[3]}"
		local _rate="${_arg_menu[4]}"
		_polap_lib_random-get
		local _seed=${_polap_var_random_number}
		rm -f "${_outfile}"

		seqkit sample \
			-p "${_rate}" \
			-s "${_seed}" \
			"${_infile}" \
			-o ${_outfile} 2>${_polap_output_dest}

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
	fi

	# subsample infile.fq outfile.fq.gz -c 10 --species
	if [[ "${_arg_menu[1]}" == "subsample" ]]; then

		local _infile="${_arg_menu[2]}"
		local _outfile="${_arg_menu[3]}"
		local _rate
		local _seed

		_polap_log0 "subsample the long-read data using a given target coverage: ${_arg_coverage}x"
		if [[ ! -s "${_arg_outdir}/l.fq.txt" ]] || [[ "${_arg_redo}" == "on" ]]; then
			_polap_lib_fastq-total-length-of "${_infile}" "${_arg_outdir}/l.fq.txt"
		fi
		local _l=$(<"${_arg_outdir}/l.fq.txt")

		local _v
		if [[ -z "${_arg_genomesize}" ]]; then
			if [[ -s "${_arg_outdir}/ncbi-species-genome-size.txt" ]]; then
				_v=$(<"${_arg_outdir}/ncbi-species-genome-size.txt")
			else
				_v=$(_polap_lib_ncbi-query-genome-size "${_arg_species}")
			fi
		else
			_v="${_arg_genomesize}"
		fi
		echo "${_v}" >"${_arg_outdir}/ncbi-species-genome-size.txt"

		local _coverage_long=$(echo "scale=5; ${_l} / ${_v}" | bc)
		local _rate=$(echo "scale=5; ${_arg_coverage} / ${_coverage_long}" | bc)

		_polap_lib_random-get
		_seed=${_polap_var_random_number}

		_polap_log1 "  arg1: ${_arg_species}"
		_polap_log1 "  input1: ${_infile}"
		_polap_log1 "  output1: ${_outfile}"
		_polap_log1 "  random seed: ${_seed}"
		_polap_log1 "  long-read: ${_l} (bp)"
		_polap_log1 "  genome size: ${_v} (bp)"
		_polap_log1 "  long-read coverage: ${_coverage_long}x"
		_polap_log1 "  target coverage: ${_arg_coverage}x"
		_polap_log1 "  sampling rate: ${_rate}"

		local result=$(echo "$_rate < 1" | bc)

		if [ "$result" -eq 1 ]; then
			# echo "The rate value is less than 1"
			if [[ "${_arg_dry}" == "off" ]]; then
				rm -f "${_outfile}"
				seqkit sample \
					-p "${_rate}" \
					-s "${_seed}" \
					"${_infile}" \
					-o ${_outfile} 2>${_polap_output_dest}
				# gzip "${_outfile}"
			fi

		else
			# echo "The value is not less than 1"
			_polap_log1 "  sampling rate is not less than 1: ${_rate}"
			_polap_log1 "  no subsampling of input: ${_infile}"
			_polap_log0 "  no subsampling: ${_outfile}"
			rm -f "${_outfile}"
			# ln -s "$(realpath ${_infile})" "$(realpath -m ${_outfile})"
			_polap_lib_make_relative_symlink "${_infile}" "${_outfile}"
		fi
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
	fi

	# subsample2 infile1.fq infile2.fq outfile.fq.gz -c 10 --species
	if [[ "${_arg_menu[1]}" == "subsample2" ]]; then

		local _infile1="${_arg_menu[2]}"
		local _infile2="${_arg_menu[3]}"
		local _outfile1="${_arg_menu[4]}"
		local _outfile2="${_arg_menu[5]}"
		local _rate
		local _seed

		_polap_log0 "subsample the paired-read data using a given target coverage: ${_arg_coverage}x"
		if [[ ! -s "${_arg_outdir}/s1.fq.txt" ]]; then
			_polap_lib_fastq-total-length-of "${_infile1}" "${_arg_outdir}/s1.fq.txt"
		fi
		if [[ ! -s "${_arg_outdir}/s2.fq.txt" ]]; then
			_polap_lib_fastq-total-length-of "${_infile2}" "${_arg_outdir}/s2.fq.txt"
		fi
		local _s1=$(<"${_arg_outdir}/s1.fq.txt")
		local _s2=$(<"${_arg_outdir}/s2.fq.txt")
		local _s=$((_s1 + _s2))
		if [[ "${_s1}" == "${_s2}" ]]; then
			_polap_log1 "Two of the pair are the same in the number of reads."

		else
			_polap_log0 "  short-read1: ${_s1} (bp)"
			_polap_log0 "  short-read2: ${_s2} (bp)"
			_polap_log0 "ERROR: two of the pair are different in the number of reads."
		fi

		local _v
		if [[ -z "${_arg_genomesize}" ]]; then
			if [[ -s "${_arg_outdir}/ncbi-species-genome-size.txt" ]]; then
				_v=$(<"${_arg_outdir}/ncbi-species-genome-size.txt")
			else
				_v=$(_polap_lib_ncbi-query-genome-size "${_arg_species}")
			fi
		else
			_v="${_arg_genomesize}"
		fi
		echo "${_v}" >"${_arg_outdir}/ncbi-species-genome-size.txt"

		local _coverage_short=$(echo "scale=5; ${_s} / ${_v}" | bc)
		local _rate=$(echo "scale=5; ${_arg_coverage} / ${_coverage_short}" | bc)

		_polap_lib_random-get
		local _seed=${_polap_var_random_number}

		_polap_log1 "  arg1: ${_arg_species}"
		_polap_log1 "  input1: ${_infile1}"
		_polap_log1 "  input2: ${_infile2}"
		_polap_log1 "  output1: ${_outfile1}"
		_polap_log1 "  output2: ${_outfile2}"
		_polap_log1 "  random seed: ${_seed}"
		_polap_log1 "  short-read1: ${_s1} (bp)"
		_polap_log1 "  short-read2: ${_s2} (bp)"
		_polap_log1 "  short-read: ${_s} (bp)"
		_polap_log1 "  short-read coverage: ${_coverage_short}x"
		_polap_log1 "  genome size: ${_v} (bp)"
		_polap_log1 "  target coverage: ${_arg_coverage}x"
		_polap_log1 "  sampling rate: ${_rate}"

		local result=$(echo "$_rate < 1" | bc)

		if [ "$result" -eq 1 ]; then
			# echo "The rate value is less than 1"

			# Example:
			# seqtk sample -s100 read1.fq 0.1 >sub1.fq
			# seqtk sample -s100 read2.fq 0.1 >sub2.fq

			rm -f "${_outfile1}"
			rm -f "${_outfile2}"
			if [[ "${_arg_dry}" == "off" ]]; then
				seqtk sample \
					-s"${_seed}" \
					"${_infile1}" \
					"${_rate}" \
					>"${_outfile1}"

				seqtk sample \
					-s"${_seed}" \
					"${_infile2}" \
					"${_rate}" \
					>"${_outfile2}"

				# gzip "${_outfile1}"
				# gzip "${_outfile2}"
			fi
		else
			# echo "The value is not less than 1"
			_polap_log1 "  sampling rate is not less than 1: ${_rate}"
			_polap_log1 "  no subsampling of input: ${_infile1}"
			_polap_log1 "  no subsampling of input: ${_infile2}"
			_polap_log0 "  no subsampling: ${_outfile1}"
			_polap_log0 "  no subsampling: ${_outfile2}"
			rm -f "${_outfile1}"
			rm -f "${_outfile2}"
			_polap_lib_make_relative_symlink "${_infile1}" "${_outfile1}"
			_polap_lib_make_relative_symlink "${_infile2}" "${_outfile2}"
			# ln -s "$(realpath ${_infile1})" "$(realpath -m ${_outfile1})"
			# ln -s "$(realpath ${_infile2})" "$(realpath -m ${_outfile2})"
		fi

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
	fi

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_polap_fastq-inspect {
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	# Print help message if requested
	help_message=$(
		cat <<HEREDOC
Determine the fastq data type.

Example: 
$(basename $0) ${_arg_menu[0]} infile.fq outfile.txt
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then

		local _infile="${_arg_menu[2]}"

		seqkit stats "${_infile}" >&3

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
	fi

	local infile="${_arg_menu[1]}"
	local outfile="${_arg_menu[2]}"

	_polap_lib_fastq-pacbio-type "${infile}" "${outfile}"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_polap_fastq-sample-to {
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
  polap fastq-sample-to - sample down to a fixed size

Synopsis:
  polap fastq-sample-to [option]

Description:
  polap fastq-sample-to samples the long-read or short-read data down to
  a fixed size.

Options:
  -l FASTQ
    reads data file

  -g INT
    the fixed size

  --outfile FILE
    output file

Examples:
  Sample data:
    polap fastq-sample-to -l l.fq -g 30m --outfile o.fq

TODO:
  Dev.

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (1998–2018)

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

	if [[ "${_arg_long_reads_is}" == "on" ]]; then
		if [[ -n "${_arg_genomesize}" ]]; then
			_polap_lib_fastq-sample-to "${_arg_long_reads}" "${_arg_outfile}" "${_arg_genomesize}"
		else
			_polap_log0 "use -g to set the max size of the output fastq file"
		fi
	elif [[ "${_arg_long_reads_is}" == "off" ]]; then
		_polap_log0 "use -l to set the input fastq file"
	else
		_polap_log0 "infile: no such file or empty file: ${_arg_long_reads}"
	fi

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}
