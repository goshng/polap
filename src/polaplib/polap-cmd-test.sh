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
# Subcommands for testing polap codes.
#
# Example:
# polap test <menu>
################################################################################

################################################################################
# Ensure that the current script is sourced only once
source "${_POLAPLIB_DIR}/run-polap-function-include.sh"
_POLAP_INCLUDE_=$(_polap_include "${BASH_SOURCE[0]}")
set +u
[[ -n "${!_POLAP_INCLUDE_}" ]] && return 0
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

function _run_polap_test {
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'
	source "${_POLAPLIB_DIR}/polap-package-common.sh"   # '.' means 'source'

	# Print help message if requested
	help_message=$(
		cat <<HEREDOC
Test codes.

Arguments:
  -i ${_arg_inum}: source Flye (usually whole-genome) assembly number

Inputs:
  ${_polap_var_ga_annotation_all}

Outputs:
  ${_polap_var_mtcontigname}

See:
  run-polap-select-contigs-by-table-1.R for the description of --select-contig option

Example:
$(basename $0) ${_arg_menu[0]} [-i|--inum <arg>] [-j|--jnum <arg>] [--select-contig <number>]
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
		exit $EXIT_SUCCESS
	fi

	if [[ "${_arg_menu[1]}" == "gfa2fasta" ]]; then
		if _polap_gfatools-gfa2fasta; then
			_polap_log0 "success: $?"
		else
			_polap_error_message $?
		fi
	fi

	if [[ "${_arg_menu[1]}" == "archive" ]]; then
		_polap_log0 "_arg_archive=${_arg_archive}"
		_polap_log0 "_ppack_var_outdir=${_ppack_var_outdir}"
		_polap_log0 "_polap_var_ga=${_polap_var_ga}"
		_polap_log0 "_ppack_var_ga=${_ppack_var_ga}"
		_polap_log0 "_ppack_var_mt_edges=${_ppack_var_mt_edges}"
	fi

	if [[ "${_arg_menu[1]}" == "gunzip-short-read" ]]; then
		if [[ -s "${_arg_short_read1}" ]]; then
			_arg_short_read1=$(_polap_gunzip-fastq "${_arg_short_read1}")
		else
			_arg_short_read1=""
		fi
		if [[ -s "${_arg_short_read2}" ]]; then
			_arg_short_read2=$(_polap_gunzip-fastq "${_arg_short_read2}")
		else
			_arg_short_read2=""
		fi

		_polap_log0 "short-read1: ${_arg_short_read1}"
		_polap_log0 "short-read2: ${_arg_short_read2}"
	fi

	# All log to the polap.log
	# -q: no log to screen
	# -v: log0 and log1 to screen
	# -v -v: log0, log1 and log2 to screen
	# -v -v -v: log0, log1, log2, and log3 to screen
	if [[ "${_arg_menu[1]}" == "log" ]]; then
		_polap_log0 "log0"
		_polap_log1 "log1"
		_polap_log2 "log2"
		_polap_log3 "log3"

		_polap_log0 "log0 to >test.log" >test.log
		_polap_log1 "log1 to >test.log" >>test.log
		_polap_log2 "log2 to >test.log" >>test.log
		_polap_log3 "log3 to >test.log" >>test.log
	fi

	if [[ "${_arg_menu[1]}" == "config" ]]; then
		_polap_create_config config.cfg
		# Set some configuration values.
		_polap_set_config_value config.cfg "name" "Alice Lee"
		_polap_set_config_value config.cfg "username" "alice"
		_polap_set_config_value config.cfg "password" "secret123"
		_polap_set_config_value config.cfg "port" "8080"

		# Read configuration values.
		_polap_log0 "Name: $(_polap_get_config_value config.cfg name)"
		_polap_log0 "Username: $(_polap_get_config_value config.cfg username)"
		_polap_log0 "Password: $(_polap_get_config_value config.cfg password)"
		_polap_log0 "Port: $(_polap_get_config_value config.cfg port)"
	fi

	# src/polap.sh test disassemble-r -v -v -v -o Juncus_roemerianus
	if [[ "${_arg_menu[1]}" == "disassemble-r" ]]; then
		local _disassemble_dir="${_arg_outdir}/disassemble"
		_arg_disassemble_i=1
		_disassemble_i="${_disassemble_dir}/${_arg_disassemble_i}"
		_disassemble_i_stage="${_disassemble_i}/1"
		_polap_log3_pipe "Rscript --vanilla ${_POLAPLIB_DIR}/polap-r-disassemble.R \
      --table ${_disassemble_i_stage}/summary1.txt \
      --out ${_disassemble_i_stage}/summary1-ordered.txt \
      --plot ${_disassemble_i_stage}/summary1-ordered.pdf \
      2>${_polap_output_dest}"
		# 2>&1"
		_rstatus=$?
		if [[ "${_rstatus}" -ne 0 ]]; then
			_polap_log0 "error in Rscript: ${_rstatus}"
		else
			_polap_log0 "okay"
		fi
	fi

	# src/polap.sh test compare2ptdna -v -v -v
	if [[ "${_arg_menu[1]}" == "compare2ptdna" ]]; then
		_arg_disassemble_c="ptdna-Juncus_roemerianus-ptgaul.fa"
		_ptdna="Juncus_roemerianus/disassemble/1/1/1/52-mtdna/ptdna.1.fa"
		_ptdir="Juncus_roemerianus/disassemble/1/1/1/52-mtdna/b"
		_var_mtdna="Juncus_roemerianus/disassemble/1/1/1/52-mtdna"
		_polap_log3_pipe "python ${_POLAPLIB_DIR}/polap-py-compare2ptdna.py \
    		      --seq1 ${_arg_disassemble_c} \
	    	      --seq2 ${_ptdna} \
		          --out ${_ptdir} \
		          2>$_polap_output_dest"
		pident=$(<"${_var_mtdna}/b/pident.txt")
		_polap_log0 "pident: ${pident}"
	fi

	if [[ "${_arg_menu[1]}" == "system" ]]; then
		system_genus_species >&3
	fi

	if [[ "${_arg_menu[1]}" == "increment-label" ]]; then
		local inum_next=$(_polap_lib_string-increment_label "${_arg_inum}")
		_polap_log0 "j: ${_arg_inum}"
		_polap_log0 "next i: ${inum_next}"
	fi

	if [[ "${_arg_menu[1]}" == "unit" ]]; then
		local _a=1.54m
		_polap_log0 "_a before: ${_a}"
		local _a=$(_polap_utility_convert_unit_to_bp "${_a}")
		_polap_log0 "_a converted: ${_a}"
	fi

	if [[ "${_arg_menu[1]}" == "link" ]]; then
		local infile="${_arg_menu[2]}"
		local outfile="${_arg_menu[3]}"
		_polap_lib_filepath-smart_ln_s2 "${infile}" "${outfile}"
	fi

	if [[ "${_arg_menu[1]}" == "var" ]]; then
		local annotatedir="annotate-read-mt"
		local resolved_fastq="resolved.fastq"
		_polap_lib_test-variable-scope \
			-o "${annotatedir}" \
			-l "${resolved_fastq}" \
			--reference "${annotatedir}"/pt.2.gfa
		_polap_log0 "top-level outdir: ${_arg_outdir}"
		_polap_log0 "top-level lfq: ${_arg_long_reads}"
	fi

	if [[ "${_arg_menu[1]}" == "polap-common-var" ]]; then
		local annotatedir="annotate-read-mt"
		local resolved_fastq="resolved.fastq"
		_polap_lib_test-polap-common-variable \
			-o "${annotatedir}" \
			-l "${resolved_fastq}" \
			-i 5 -j 7 \
			--reference "${annotatedir}"/pt.2.gfa
		_polap_log0 "top-level outdir: ${_arg_outdir}"
		_polap_log0 "top-level lfq: ${_arg_long_reads}"
		_polap_log0 "top-level i: ${_arg_inum}"
		_polap_log0 "top-level j: ${_arg_jnum}"
	fi

	if [[ "${_arg_menu[1]}" == "yaml" ]]; then
		_polap_log0 "nuc size: ${_arg_sim_nuc_size}"
		_polap_log0 "depth nuc: ${_arg_sim_depth_nuc}"
		_polap_log0 "depth mt: ${_arg_sim_depth_mt}"
		_polap_log0 "depth pt: ${_arg_sim_depth_pt}"
		_polap_log0 "a-b-c: ${PCFG_A_B_C}"
		_polap_log0 "x-y-z: ${PCFG-Y}"
		compgen -v | grep '^PCFG_'
	fi

	if [[ "${_arg_menu[1]}" == "step" ]]; then
		local _include="${_arg_steps_include}"
		local _exclude="${_arg_steps_exclude}" # Optional range or list of steps to exclude
		local _step_array=()

		_step_array=($(_polap_parse_steps "${_include}" "${_exclude}"))
		_rstatus="$?"
		if [[ "${_rstatus}" -ne 0 ]]; then
			_polap_log0 "ERROR: Error parsing steps."
			return "${_POLAP_ERR_CMD_OPTION_STEPS}"
		fi

		if _polap_contains_step 1 "${_step_array[@]}"; then
			_polap_log0 "step 1"
		fi

		if _polap_contains_step 2 "${_step_array[@]}"; then
			_polap_log0 "step 2"
		fi

		if _polap_contains_step 3 "${_step_array[@]}"; then
			_polap_log0 "step 3"
		fi

		if _polap_contains_step 4 "${_step_array[@]}"; then
			_polap_log0 "step 4"
		fi

		if _polap_contains_step 5 "${_step_array[@]}"; then
			_polap_log0 "step 5"
		fi
	fi

	if [[ "${_arg_menu[1]}" == "math" ]]; then
		local v=$(_polap_lib_math-percentify 0.87345)
		_polap_log0 "0.87345 -> $v"
		local v=$(_polap_lib_math-percentify 0.87345 3 floor) # 87.345  (floor to 3)
		_polap_log0 "0.87345 -> $v"
		_polap_lib_math-percentify 1e-3 1 ceil # 0.1     (ceil to 1)
		# _polap_lib_math-percentify -0.1            # ERROR with caller file:line (exits 1)
		# _polap_lib_math-percentify 1.2 # ERROR (out of [0,1])
	fi

	if [[ "${_arg_menu[1]}" == "fastq-estimate" ]]; then
		local v=$(_polap_lib_fastq-estimate-bases ${_arg_long_reads})
		local v_str=$(_polap_lib_unit-convert_bp ${v})
		_polap_log0 "${_arg_long_reads} size: $v ($v_str)"
	fi

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}
