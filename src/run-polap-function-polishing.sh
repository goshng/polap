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
# Prepares the polishing using FMLRC.
# Arguments:
#   -a s1.fq
#   -b s2.fq
# Inputs:
#   s1.fq
#   s2.fq
# Outputs:
#   $$ODIR/msbwt
################################################################################
function _run_polap_prepare-polishing() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	_polap_set-variables-short-read
	source "$script_dir/polap-variables-base.sh" # '.' means 'source'

	help_message=$(
		cat <<HEREDOC
# Prepares the polishing using FMLRC.
# Arguments:
#   -a $SR1
#   -b $SR2
#   or
#   --bioproject ${_arg_bioproject}
# Inputs:
#   $SR1
#   $SR2
# Outputs:
#   $ODIR/msbwt/comp_msbwt.npy
# Precondition:
#   get-bioproject --bioproject ${_arg_bioproject}
Example: $(basename "$0") ${_arg_menu[0]} [-a|--short-read1 <arg>] [-b|--short-read2 <arg>]
Example: $(basename "$0") ${_arg_menu[0]} -o $ODIR
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && _polap_echo0 "${help_message}" && exit $EXIT_SUCCESS

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		ls -l "${ODIR}/msbwt" >&2
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		exit $EXIT_SUCCESS
	fi

	if [ -s "${_polap_var_base_msbwt_tar_gz}" ]; then
		_polap_log1_file "${_polap_var_base_msbwt_tar_gz}"
		if [[ -s "${_polap_var_base_msbwt}" ]]; then
			_polap_log1_file "${_polap_var_base_msbwt}"
			_polap_log1 "  skipping the short-read polishing preparation."
		else
			tar -zxf "${_polap_var_base_msbwt_tar_gz}" -C "${ODIR}"
		fi
	elif [[ -s "${_polap_var_base_msbwt}" ]]; then
		_polap_log1_file "${_polap_var_base_msbwt}"
		_polap_log1 "  skipping the short-read polishing preparation."
	else

		source $HOME/miniconda3/bin/activate polap-fmlrc

		if ! run_check2; then
			echoerr "ERROR: change your conda environment to polap-fmlrc."
			echoerr "INFO: (base) $ conda env create -f src/environment-fmlrc.yaml"
			echoerr "INFO: (base) $ conda activate polap-fmlrc"
			exit $EXIT_ERROR
		fi

		check_file_existence "${SR1}"
		check_file_existence "${SR2}"

		_polap_log1 "excuting ropebwt2 and msbwt on the short reads ... be patient!"
		if [[ $SR1 = *.fastq || $SR1 = *.fq ]]; then
			cat "$SR1" "$SR2" |
				awk 'NR % 4 == 2' | sort | tr NT TN |
				ropebwt2 -LR 2>"${_polap_output_dest}" |
				tr NT TN |
				msbwt convert "$ODIR"/msbwt \
					>/dev/null 2>&1
		elif [[ $SR1 = *.fq.gz ]] || [[ $SR1 = *.fastq.gz ]]; then
			zcat "$SR1" "$SR2" |
				awk 'NR % 4 == 2' | sort | tr NT TN |
				ropebwt2 -LR 2>"${_polap_output_dest}" |
				tr NT TN |
				msbwt convert "$ODIR"/msbwt \
					>/dev/null 2>&1
		fi
		conda deactivate
	fi

	_polap_log1 "NEXT: $(basename $0) polish [-p mt.0.fasta] [-f mt.1.fa]"

	_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}

################################################################################
# Polishes using FMLRC.
# Arguments:
#   -p mt.0.fasta
#   -f mt.1.fa
# Inputs:
#   $ODIR/msbwt/comp_msbwt.npy
#   $PA
# Outputs:
#   $FA
################################################################################
function _run_polap_polish() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	help_message=$(
		cat <<HEREDOC
# Polishes using FMLRC.
# Arguments:
#   -p $PA: a long-read draft genome assembly
#   -f $FA: a final genome assembly sequence name
# Inputs:
#   $ODIR/msbwt/comp_msbwt.npy
#   $PA
# Outputs:
#   $FA
Example: $(basename "$0") ${_arg_menu[0]} [-p|--unpolished-fasta <arg>] [-f|--final-assembly <arg>]
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && _polap_echo0 "${help_message}" && exit $EXIT_SUCCESS

	source $HOME/miniconda3/bin/activate polap-fmlrc

	if ! run_check2; then
		_polap_log0 "ERROR: change your conda environment to polap-fmlrc."
		_polap_log0 "INFO: (base) $ conda env create -f src/environment-fmlrc.yaml"
		_polap_log0 "INFO: (base) $ conda activate polap-fmlrc"
		exit $EXIT_ERROR
	fi

	if [[ ! -s "$ODIR/msbwt/comp_msbwt.npy" ]]; then
		_polap_log0 "ERROR: no msbwt at $ODIR/msbwt/comp_msbwt.npy"
		_polap_log0 "HINT: $0 prepare-polishing [-a s1.fq] [-b s2.fq]"
		exit $EXIT_ERROR
	fi

	_polap_log1 "INFO: executing fmlrc on the draft sequence $PA ... be patient!"
	if [[ -s "${PA}" ]]; then
		fmlrc -p "${NT}" "$ODIR"/msbwt/comp_msbwt.npy "${PA}" "${FA}" >/dev/null 2>&1
	else
		_polap_log0 "ERROR: no unpolished fasta file: [$PA]"
		exit $EXIT_ERROR
	fi

	conda deactivate

	_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}