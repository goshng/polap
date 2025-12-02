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

function _run_polap_convert {
	# Enable debugging if _POLAP_DEBUG is set
	[ "${_POLAP_DEBUG:-0}" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose:-0}" -ge "${_polap_var_function_verbose:-2}" ] && _polap_output_dest="/dev/stderr"

	# Common variables & helpers
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	local polap_cmd="${FUNCNAME##*_}"
	local help_message=$(
		cat <<EOF
Name:
  polap ${polap_cmd} - convert between common graph/sequence formats

Synopsis:
  polap ${polap_cmd} <in2out> <input> <output> [options]

Description:
  Conversions are delegated to well-tested command-line tools.

  Supported menu:
    - gfa2fasta : Use gfatools gfa2fa; optional ID subset via seqkit/seqtk.
    - himt2prefix: Extract himt HTML to a report out directory

Positional arguments:
  in2out
    A conversion key (e.g., gfa2fasta).

  input
    Input path (file, file.gz, or '-' for stdin).

  output
    Output path (file, file.gz, or '-' for stdout).

Options:
  --ids LIST_OR_FILE
    Comma/space list or a file of IDs (one per line). '+/-' suffix stripped.
  --line-width N
    FASTA to width N (0 = single line). Uses seqkit or seqtk if present.

Examples:
  polap ${polap_cmd} gfa2fasta graph.gfa out.fa

  polap ${polap_cmd} gfa2fasta graph.gfa.gz out.fa.gz

  polap ${polap_cmd} gfa2fasta - out.fa --ids ids.txt

  polap ${polap_cmd} himt2prefix himt-in.html out-prefix

Notes:
  - Implementation uses: gfatools gfa2fa; seqkit grep/seq or seqtk subseq/seq for filtering/rewrap.

  - Unknown '*' sequences in S-lines are not invented (left out by gfatools).
EOF
	)

	_polap_lib_help-maybe-show3 "convert" help_message || return 0

	# Reserved "view"
	if [[ "${_arg_menu[1]:-}" == "view" ]]; then
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		[ "${_POLAP_DEBUG:-0}" -eq 1 ] && set +x
		return 0
	fi

	# Ensure env
	_polap_lib_conda-ensure_conda_env polap || return 1

	# Parse "convert in2out: ${_arg_menu[1]}"
	local conv="${_arg_menu[1]:-}"
	local inpath="${_arg_menu[2]:-}"
	local outpath="${_arg_menu[3]:-}"

	# Global-option variables (may be unset)
	local ids="${_arg_ids:-}"
	local line_width="${_arg_line_width:-}"
	local stable="${_arg_stable:-off}"

	if [[ -z "$conv" || -z "$inpath" || -z "$outpath" ]]; then
		_polap_log0 "ERROR: usage: polap ${polap_cmd} <in2out> <input> <output> [--ids ...]"
		echo "$help_message" >&2
		conda deactivate || true
		[ "${_POLAP_DEBUG:-0}" -eq 1 ] && set +x
		return 2
	fi

	_polap_log1 "convert: $conv"
	_polap_log1 "input  : $inpath"
	_polap_log1 "output : $outpath"
	_polap_log2 "ids    : ${ids:-<none>}"

	case "$conv" in
	gfa2fasta | gfa2fa)
		# Exactly this name, as requested
		local runner="${_POLAPLIB_DIR}/polap-bash-convert-gfa2fasta.sh"
		if [[ ! -s "$runner" ]]; then
			_polap_log0 "ERROR: Not found: $runner"
			conda deactivate || true
			[ "${_POLAP_DEBUG:-0}" -eq 1 ] && set +x
			return 3
		fi

		# Assemble arguments for the converter
		local -a conv_args=()
		[[ -n "$ids" ]] && conv_args+=(--ids "${_arg_ids}")
		# [[ -n "$line_width" ]] && conv_args+=(--line-width "$line_width")
		# [[ "$stable" == "on" ]] && conv_args+=(--stable)

		_polap_log2 "dispatch: bash \"$runner\" ${conv_args[*]} \"$inpath\" \"$outpath\""
		if ! bash "$runner" "${conv_args[@]}" "$inpath" "$outpath" 2>$_polap_output_dest; then
			_polap_log0 "ERROR: polap-bash-convert-gfa2fasta.sh failed."
			conda deactivate || true
			[ "${_POLAP_DEBUG:-0}" -eq 1 ] && set +x
			return 4
		fi
		;;

	himt2prefix)
		# Exactly this name, as requested
		# local runner="${_POLAPLIB_DIR}/polap-bash-parse-himt-mito-report.sh"
		local runner="${_POLAPLIB_DIR}/polap-bash-assess-mito.sh"
		if [[ ! -s "$runner" ]]; then
			_polap_log0 "ERROR: Not found: $runner"
			conda deactivate || true
			[ "${_POLAP_DEBUG:-0}" -eq 1 ] && set +x
			return 3
		fi

		mkdir -p $(dirname "$outpath")
		# if ! bash "$runner" -p "$outpath" "$inpath" 2>$_polap_output_dest; then
		if [[ -s "$inpath" ]]; then
			if ! bash "$runner" -i "$inpath" -o "$outpath" 2>$_polap_output_dest; then
				_polap_log0 "ERROR: $runner failed."
				_polap_log0 "inpath: $inpath"
				_polap_log0 "outpath: $outpath"
				_polap_log0 "bash $runner -i $inpath -o $outpath"
				conda deactivate || true
				[ "${_POLAP_DEBUG:-0}" -eq 1 ] && set +x
				return 4
			fi
		else
			_polap_log0 "Skipping because of no input HiMT HTML: $inpath"
		fi
		;;
	*)
		_polap_log0 "ERROR: unknown conversion key: ${conv}"
		echo "$help_message" >&2
		conda deactivate || true
		[ "${_POLAP_DEBUG:-0}" -eq 1 ] && set +x
		return 2
		;;
	esac

	conda deactivate || true

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	[ "${_POLAP_DEBUG:-0}" -eq 1 ] && set +x
	return 0
}
