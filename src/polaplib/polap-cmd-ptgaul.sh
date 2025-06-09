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

function _run_polap_ptgaul {
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
# Template for an external shell script
#
# Arguments:
#   -i ${_arg_inum}: source Flye (usually whole-genome) assembly number
#
# Inputs:
#   ${_polap_var_ga_annotation_all}
#
# Outputs:
#   ${_polap_var_mtcontigname}
#
# See:
#   run-polap-select-contigs-by-table-1.R for the description of --select-contig option
Example: $(basename $0) ${_arg_menu[0]} [-i|--inum <arg>] [-j|--jnum <arg>] [--select-contig <number>]
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

	# 2025-06-04: This used to be in the disassemble module.
	# I wish to have extraction and polishing in a separate menu.
	# This might be in other module like polap-cmd-gfa because we extract ptDNA from a gfa file.
	# Because we use the ptGAUL's gfa or Flye's gfa, it is fine to have this here.
	if [[ "${_arg_menu[1]}" == "extract-ptdna" ]]; then

		_contigger_edges_gfa="${_arg_outdir}/ptgaul/flye_cpONT/assembly_graph.gfa"
		_outdir="${_arg_outdir}/ptgaul/flye_cpONT/ptdna"

		_mtcontigname="${_outdir}/mt.contig.name"
		_arg_unpolished_fasta="${_outdir}/circular_path_1_concatenated.fa"
		_arg_final_assembly="${_outdir}/pt.1.fa"

		if [[ "${_arg_menu[2]}" == "outfile" ]]; then
			_polap_log3_cmd mkdir -p "${_outdir}"
			gfatools stat -l "${_contigger_edges_gfa}" \
				>"${_outdir}/graph_final.txt" \
				2>"$_polap_output_dest"
			_summary_gfa_number_segments=$(grep "Number of segments:" "${_outdir}/graph_final.txt" | awk '{print $4}')

			if [[ "${_summary_gfa_number_segments}" -eq 3 ]]; then
				_polap_log1 "ptGAUL assembly's segments: ${_summary_gfa_number_segments}"
				_polap_log1 "you may use ptgaul menu:"
				_polap_log1 "  ptgaul 1"
				_polap_log1 "  ptgaul 1 1"
				_polap_log1 "  ptgaul 1 3"
				_polap_log1 "  ptgaul 2"
				_polap_log1 "  ptgaul 3"
				# mt.contig.name
				_polap_log1 "create file: ${_mtcontigname}"
				local string="edge_1,edge_2,edge_3"
				echo "$string" | tr ',' '\n' >"${_mtcontigname}"
				# extract ptDNA from ptGAUL's result
				_polap_log3_pipe "python \
          ${_POLAPLIB_DIR}/run-polap-py-find-plastid-gfa2fasta.py \
		        --gfa ${_contigger_edges_gfa} \
		        --seed ${_mtcontigname} \
		        --out ${_outdir} \
		        2>$_polap_output_dest"

			elif [[ "${_summary_gfa_number_segments}" -eq 1 ]]; then
				_polap_log1 "  ptGAUL assembly's segments: ${_summary_gfa_number_segments}"
				local _final_assembly_fasta="${_arg_outdir}/ptgaul/ptGAUL_final_assembly/final_assembly.fasta"
				cp "${_final_assembly_fasta}" "${_arg_unpolished_fasta}"
			else
				# extract three unique edges from the gfa for a ptDNA candidate
				#
				# L	edge_3	+	edge_4	+	0M	RC:i:41
				# L	edge_3	+	edge_4	-	0M	RC:i:40
				# L	edge_3	-	edge_5	-	0M	RC:i:42
				# L	edge_3	-	edge_5	+	0M	RC:i:41
				# P	contig_5	edge_4-,edge_3-,edge_5+,edge_3+,edge_4+,edge_3-	*
				# P	contig_1	edge_1+	*
				# P	contig_2	edge_2+	*
				#
				if [[ ! -s "${_mtcontigname}" ]]; then
					bash ${_POLAPLIB_DIR}/run-polap-sh-extract-three-edges-of-ptdna.sh \
						"${_contigger_edges_gfa}" \
						"${_mtcontigname}"
				fi

				if [[ -s "${_mtcontigname}" ]]; then
					# extract ptDNA from ptGAUL's result
					_polap_log3_pipe "python \
          ${_POLAPLIB_DIR}/run-polap-py-find-plastid-gfa2fasta.py \
		        --gfa ${_contigger_edges_gfa} \
		        --seed ${_mtcontigname} \
		        --out ${_outdir} \
		        2>$_polap_output_dest"

				else
					_polap_log0 "ptGAUL assembly's segment counts: ${_summary_gfa_number_segments}"
					_polap_log0 "you may not use ptgaul menu at the moment."
					_polap_log0 "check: ${_contigger_edges_gfa}"
					_polap_log0 "create: ${_mtcontigname} with e.g., edge_1 that is originated from the ptDNA."
					_polap_log0 "execute: polap ptgaul extract-ptdna"
					_polap_log0 "or execute: polap-data-cflye run extract-ptdna-ptgaul <outdir>"
				fi
			fi
		fi

		# create _mtcontigname for manual extraction
		if [[ "${_arg_menu[2]}" == "1" ]]; then
			if [[ "${_arg_menu[3]}" == "1" ]]; then
				_polap_log0 "create file: ${_mtcontigname}"
				local string="edge_1"
				echo "$string" | tr ',' '\n' >"${_mtcontigname}"
			elif [[ "${_arg_menu[3]}" == "3" ]]; then
				_polap_log0 "create file: ${_mtcontigname}"
				local string="edge_1,edge_2,edge_3"
				echo "$string" | tr ',' '\n' >"${_mtcontigname}"
			fi
		fi

		# extract
		if [[ "${_arg_menu[2]}" == "2" ]]; then
			_polap_log0 "extract ptDNA: ${_arg_unpolished_fasta}"
			_polap_log3_pipe "python \
          ${_POLAPLIB_DIR}/run-polap-py-find-plastid-gfa2fasta.py \
		        --gfa ${_contigger_edges_gfa} \
		        --seed ${_mtcontigname} \
		        --out ${_outdir} \
		        2>$_polap_output_dest"
		fi

		# Disable debugging if previously enabled
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
	fi

	# menu: ptgaul
	# ptgaul
	if [[ "${_arg_menu[1]}" == "polish-ptdna" ]]; then

		_contigger_edges_gfa="${_arg_outdir}/ptgaul/flye_cpONT/assembly_graph.gfa"
		_outdir="${_arg_outdir}/ptgaul/flye_cpONT/ptdna"

		_mtcontigname="${_outdir}/mt.contig.name"
		_arg_unpolished_fasta="${_outdir}/circular_path_1_concatenated.fa"
		_arg_final_assembly="${_outdir}/pt.1.fa"

		if [[ "${_arg_menu[2]}" == "outfile" ]]; then
			_summary_gfa_number_segments=$(grep "Number of segments:" "${_outdir}/graph_final.txt" | awk '{print $4}')

			if [[ "${_summary_gfa_number_segments}" -eq 3 ]]; then
				# polish one of the ptDNA sequences
				_polap_log1 "  polishing ptDNA: ${_arg_unpolished_fasta}"
				_run_polap_polish
				exit_code=$?
				# If a function or command is killed due to an out-of-memory (OOM) event,
				# it usually exits with signal 137 (SIGKILL).
				if [ $exit_code -eq 137 ]; then
					die "polishing was killed due to OOM: an out-of-memory event"
				fi
				if [[ -s "${_arg_final_assembly}" ]]; then
					_polap_log1 "  ptGAUL polished assembly: ${_arg_final_assembly}"
				else
					die "ERROR: no such ptGAUL polished assembly: ${_arg_final_assembly}"
				fi
			elif [[ "${_summary_gfa_number_segments}" -eq 1 ]]; then
				_polap_log1 "  ptGAUL assembly's segments: ${_summary_gfa_number_segments}"
				_arg_unpolished_fasta="${_arg_outdir}/ptgaul/ptGAUL_final_assembly/final_assembly.fasta"
				_polap_log1 "  ptGAUL assembly: ${_arg_unpolished_fasta}"
				_polap_log1 "  polishing ptDNA: ${_arg_unpolished_fasta}"
				_run_polap_polish
				exit_code=$?
				# If a function or command is killed due to an out-of-memory (OOM) event,
				# it usually exits with signal 137 (SIGKILL).
				if [ $exit_code -eq 137 ]; then
					die "polishing was killed due to OOM: an out-of-memory event"
				fi
				if [[ -s "${_arg_final_assembly}" ]]; then
					_polap_log1 "  ptGAUL polished assembly: ${_arg_final_assembly}"
				else
					die "ERROR: no such ptGAUL polished assembly: ${_arg_final_assembly}"
				fi
			else
				if [[ -s "${_mtcontigname}" || -s "${_arg_unpolished_fasta}" ]]; then
					# polish one of the ptDNA sequences
					_polap_log0 "polishing ptDNA: ${_arg_unpolished_fasta}"
					_run_polap_polish
					exit_code=$?
					# If a function or command is killed due to an out-of-memory (OOM) event,
					# it usually exits with signal 137 (SIGKILL).
					if [ $exit_code -eq 137 ]; then
						die "polishing was killed due to OOM: an out-of-memory event"
					fi
					if [[ -s "${_arg_final_assembly}" ]]; then
						_polap_log1 "  ptGAUL polished assembly: ${_arg_final_assembly}"
					else
						die "ERROR: no such ptGAUL polished assembly: ${_arg_final_assembly}"
					fi
				else
					_polap_log0 "ptGAUL assembly's segment counts: ${_summary_gfa_number_segments}"
					_polap_log0 "you may not use ptgaul menu at the moment."
					_polap_log0 "check: ${_contigger_edges_gfa}"
					_polap_log0 "create: ${_mtcontigname} with e.g., edge_1 that is originated from the ptDNA."
					_polap_log0 "then execute:disassemble ptgaul 2"
					_polap_log0 "polap disassemble ptgaul 2 -o ${_arg_outdir}"
					_polap_log0 "polap disassemble ptgaul 3 -o ${_arg_outdir}"
				fi
			fi
		fi

		# polishing
		if [[ "${_arg_menu[2]}" == "3" ]]; then
			_polap_log0 "polishing ptDNA: ${_arg_unpolished_fasta}"
			_run_polap_polish
			_polap_log0 "ptGAUL polished assembly: ${_arg_final_assembly}"
		fi

		if [[ "${_arg_menu[2]}" == "4" ]]; then
			for ((i = 1; i <= 4; i++)); do
				_arg_unpolished_fasta="${_outdir}/circular_path_${i}_concatenated.fa"
				_arg_final_assembly="${_outdir}/pt.${i}.fa"
				_polap_log0 "polishing ptDNA: ${_arg_unpolished_fasta}"
				_run_polap_polish
				_polap_log0 "ptGAUL polished assembly: ${_arg_final_assembly}"
			done
		fi

		# Disable debugging if previously enabled
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
	fi

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}
