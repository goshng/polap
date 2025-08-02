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
# Functions for subcommand template ...
# Describe what they are and what they do.
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

function _run_polap_fasta {
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
  polap fasta - process FASTA files

Synopsis:
  polap fasta [options]

Description:
  polap fasta processes FASTA files for various purposes, such as
	sequence extraction, filtering, and processing.

Options:
  -l FASTA
    reads data file

Examples:
  Plot sequence length distribution:
		polap fasta length -l sequences.fasta

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

	if [[ "${_arg_menu[1]}" == "length" ]]; then
		local _arg_fasta="${_arg_long_reads}"
		if [[ -z "${_arg_fasta}" ]]; then
			echo "[ERROR] FASTA file is required." >&2
			return 1
		fi

		if [[ ! -f "${_arg_fasta}" ]]; then
			echo "[ERROR] FASTA file '${_arg_fasta}' does not exist." >&2
			return 1
		fi

		_polap_log3 "Processing FASTA file: ${_arg_fasta}"
		local seq_count=$(grep -c '^>' "${_arg_fasta}")
		_polap_log0 "Number of sequences in ${_arg_fasta}: ${seq_count}"

		plot_seq_length_dist_pdf "${_arg_fasta}" "${_arg_outfile:-read_length_distribution.pdf}"
		return 0
	fi

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

plot_seq_length_dist_pdf() {
  local input="$1"
  local output_pdf="${2:-read_length_distribution.pdf}"
  local tmp_len
  tmp_len=$(mktemp)

  # Extract sequence lengths
  seqkit fx2tab -n -l "$input" | cut -f2 > "$tmp_len"

  # Generate PDF plot using R
  Rscript - <<EOF
x <- scan("$tmp_len")
pdf("$output_pdf", width=8, height=6)
hist(x, breaks=100, main="Read Length Distribution", xlab="Sequence Length (bp)", col="skyblue", border="white")
dev.off()
EOF

  echo "PDF saved to $output_pdf"
  rm -f "$tmp_len"
}
