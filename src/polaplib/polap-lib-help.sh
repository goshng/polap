#!/usr/bin/env bash
# polaplib/polap-lib-help.sh
# Version: v0.1.0
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Common helper to show a help_message as a man page if _brg_help=on

set -euo pipefail

# _polap_lib_help-maybe-show "bolap_cmd" help_message
# --------------------------------------------------
# Args:
#   $1: bolap_cmd or command name (used for man title)
#   $2: variable name that holds the heredoc help text (e.g. help_message)
#
# Behavior:
#   - If _brg_help == on, converts the help message to a temporary man file
#     via _bolap_lib_man-convert_help_message and opens it with `man`.
#   - Always cleans up the temp file.
#   - Safe under `set -euo pipefail` and uses declare -n for the reference.
#
# Usage:
# local bolap_cmd="${FUNCNAME%%_*}"
# 	help_message=$(
# 		cat <<EOF
# Name:
#   bolap $bolap_cmd
# EOF
# 	)
#
# for bolap
# _polap_lib_help-maybe-show "$bolap_cmd" help_message || return 0
#
# for polap
# _polap_lib_help-maybe-show3 "$polap_cmd" help_message || return 0

_polap_lib_help-maybe-show() {
	local cmd="${1:?}"
	local msgvar="${2:?}"

	# nothing to do unless help requested
	if [[ "${_brg_help:-off}" == "off" && "${_arg_help:-off}" == "off" ]]; then
		return 0
	fi

	# reference to the variable that holds the help_message string
	declare -n ref="$msgvar"

	local manfile
	manfile=$(_bolap_lib_man-convert_help_message "$ref" "$cmd")
	man "$manfile"
	rm -f "$manfile"
	return 1
}

_polap_lib_help-maybe-show3() {
	local cmd="${1:?}"
	local msgvar="${2:?}"

	# nothing to do unless help requested
	if [[ "${_brg_help:-off}" == "off" && "${_arg_help:-off}" == "off" ]]; then
		return 0
	fi

	# reference to the variable that holds the help_message string
	declare -n ref="$msgvar"

	local manfile
	manfile=$(_bolap_lib_man-convert_help_message "$ref" "$cmd")
	man "$manfile" >&3
	rm -f "$manfile"
	return 1
}
