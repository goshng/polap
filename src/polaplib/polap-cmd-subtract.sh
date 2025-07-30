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

function _run_polap_subtract {
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
Subtract reads mapped on a ptDNA sequence from the long-read data used
for assembling the ptDNA.

Arguments:
  -l ${_arg_long_reads}
  -p ptdna.fa

Inputs:
  ${_arg_long_reads}: the input long-read fastq data
  ptdna.fa: 4 DNA sequences from a plastid DNA sequence, 
    ptDNA has LSC, SSC, and two IR regions, from which we have 4 possible linear sequences.
    We could say that we have two if reverse complemenaty sequences are the same.
    Anyway, we map reads on these ptDNA sequences to remove them from the input long-read.
    In order to map reads on the region of the beginning and ending of each sequence,
    we duplicate each sequence and concatenate the duplicates before mapping
    reads on the ptDNA fasta sequences.

Outputs:
  a new l.fq without reads mapped on the ptDNA.

See Also:
  disassemble/infer-1/2/3/52-mtdna
  the 2nd stage has the index that is selected for the 3rd stage.
  Use the index to locate 52-mtdna and then
  cat *_concatenated.fa >ptdna4.fa
  do not polish it before using it because the assembly is from the long-read.

Example:
$(basename $0) ${_arg_menu[0]} -l l.fq -p ptdna.fa
$(basename $0) ${_arg_menu[0]} -l SRR7153095.fastq -p ptdna4.fa
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

	# map the long-read data on the ptdna
	if [[ "${_arg_menu[1]}" == "1" ]]; then

		_polap_log0 "mapping long-read data on the seed contigs ..."

		# _polap_log2_cat "${_polap_var_oga_contig}/contig_total_length.txt"
		# local _contig_length_bp=$(_polap_utility_convert_bp ${CONTIG_LENGTH})
		# _polap_log1 "    organelle genome size based on the seed contig selection: ${_contig_length_bp}"
		mkdir -p "${_arg_outdir}"/subtract

		_polap_log1 "  mapping long-read data on the seed contigs using minimap2 ..."
		_polap_log2 "    input1: ${_arg_unpolished_fasta}"
		_polap_log2 "    input2: ${_arg_long_reads}"
		_polap_log2 "    output: ${_arg_outdir}/subtract/contig.paf"
		if [[ -s "${_arg_outdir}"/subtract/contig.paf ]] && [[ "${_arg_redo}" = "off" ]]; then
			_polap_log1 "  found: ${_arg_outdir}/subtract/contig.paf, skipping the minimap2 mapping step ..."
		else
			_polap_log3_pipe "minimap2 -cx \
        ${_arg_minimap2_data_type} \
        ${_arg_unpolished_fasta} \
        ${_arg_long_reads} \
        -t ${_arg_threads} \
        -o ${_arg_outdir}/subtract/contig.paf \
        >${_polap_output_dest} 2>&1"
		fi

		_polap_log1 "  converting PAF to TAB ..."
		_polap_log2 "    input1: ${_arg_outdir}/subtract/contig.paf"
		_polap_log2 "    output: ${_arg_outdir}/subtract/contig.tab"
		_polap_log3_cmd bash "${_POLAPLIB_DIR}/polap-bash-minimap2-paf2tab.sh" "${_arg_min_read_length}" "${_arg_outdir}/subtract/contig.paf" "${_arg_outdir}/subtract/contig.tab"

	fi

	# grep '^@' input.fastq | paste - - - - | grep -Po '^@\K\S+'
	if [[ "${_arg_menu[1]}" == "2x" ]]; then

		local all_names=$(mktemp)
		local keep_names=$(mktemp)

		grep '^@' ${_arg_long_reads} | paste - - - - | grep -Po '^@\K\S+' >"$all_names"

		grep -vFf "$remove_names" "$all_names" >"$keep_names"
		seqtk subseq "$input" "$keep_names" >"$output"

		rm -f "$all_names" "$keep_names"
	fi

	if [[ "${_arg_menu[1]}" == "2" ]]; then

		# mkdir -p ${_polap_var_oga_reads}/${_pread_sel}/${i}
		# mkdir -p ${_polap_var_oga_seeds}/${_pread_sel}/${i}
		# mkdir -p ${_polap_var_oga_subsample}/${_pread_sel}/${i}

		local _pread_sel="ptgaul-reads"
		mkdir -p ${_arg_outdir}/subtract/${_pread_sel}

		_polap_log3_pipe "Rscript --vanilla ${_POLAPLIB_DIR}/polap-r-minimap2readname.R \
		  -t ${_arg_outdir}/subtract/contig.tab \
		  --out ${_arg_outdir}/subtract/${_pread_sel} \
      -w ${_arg_single_min} \
		  -r ${_arg_pair_min} \
		  -x ${_arg_bridge_min} \
      --create-ptgaul \
		  >${_polap_output_dest} 2>&1"

		# local all_names=$(mktemp)
		local all_names="${_arg_outdir}/subtract/${_pread_sel}/ptgaul.all.names"
		local keep_names="${_arg_outdir}/subtract/${_pread_sel}/ptgaul.keep.names"
		local remove_names="${_arg_outdir}/subtract/${_pread_sel}/ptgaul.names"

		# grep '^@' ${_arg_long_reads} | paste - - - - | grep -Po '^@\K\S+' >"$all_names"
		awk 'NR % 4 == 1 {sub(/^@/, "", $0); print $1}' "$_arg_long_reads" >"$all_names"

		# grep -vxFf ptgaul.names ptgaul.all.names > ptgaul.keep.names
		# -x: match entire line.
		# -v: invert.
		# -F: fixed string.
		# -f: pattern file.
		grep -vxFf "$remove_names" "$all_names" >"$keep_names"

		_polap_log3_pipe "seqtk subseq \
		    ${_arg_long_reads} \
		    ${keep_names} |\
		    gzip >${_arg_outdir}/subtract/${_pread_sel}/0.fq.gz"
	fi

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}
