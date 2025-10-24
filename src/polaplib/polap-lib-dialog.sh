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

#!/usr/bin/env bash
# polaplib/polap-lib-dialog.sh
#
# Version: v0.1.0
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Purpose:
#   Lightweight, portable TTY dialog utilities for polap/bolap scripts.
#   Provides yes/no confirmation compatible with `set -euo pipefail`.
#
# Usage:
#   source "${POLAPLIB_DIR}/polap-lib-dialog.sh"
#   if _polap_lib_dialog-yes-no "Proceed with deletion?"; then
#       echo "yes"
#   else
#       echo "no"
#   fi
#
# Behavior:
#   • If stdin is a terminal, prompt the user and read the answer.
#   • Accept "y", "yes" (case-insensitive) -> return 0.
#   • Anything else (including Enter or EOF) -> return 1.
#   • In non-interactive mode (no TTY):
#       - If POLAP_ASSUME_YES=1 -> return 0.
#       - Else POLAP_ASSUME_NO=1 -> return 1.
#       - Otherwise: log a warning and return 1 (safe default: NO).
#

set -euo pipefail

# stderr logger
_polap_echo_stderr() { printf '%s\n' "$*" >&2; }

# --------------------------------------------------------------------------
# Function: _polap_lib_dialog-yes-no
# Ask a yes/no question. Returns 0 for "yes", 1 for "no".
# --------------------------------------------------------------------------
_polap_lib_dialog-yes-no() {
	local question="${1:-Proceed?}"
	local default="${2:-y}" # default: 'y' = yes, 'n' = no
	local ans=""

	# non-interactive guard
	if [[ ! -t 0 ]]; then
		if [[ "${POLAP_ASSUME_YES:-0}" == "1" ]]; then
			_polap_echo_stderr "[INFO] Non-interactive; assuming YES for: ${question}"
			return 0
		elif [[ "${POLAP_ASSUME_NO:-0}" == "1" ]]; then
			_polap_echo_stderr "[INFO] Non-interactive; assuming NO for: ${question}"
			return 1
		else
			_polap_echo_stderr "[WARN] Non-interactive; defaulting to NO for: ${question}"
			return 1
		fi
	fi

	# build prompt label
	local prompt
	if [[ "$default" =~ ^[Yy]$ ]]; then
		prompt="[Y/n]"
	else
		prompt="[y/N]"
	fi

	while true; do
		printf "%s %s " "$question" "$prompt" >/dev/tty
		IFS= read -r ans </dev/tty || ans=""
		ans="${ans,,}" # lowercase
		case "$ans" in
		y | yes) return 0 ;;
		n | no) return 1 ;;
		"") # empty (use default)
			if [[ "$default" =~ ^[Yy]$ ]]; then
				return 0
			else
				return 1
			fi
			;;
		*) _polap_echo_stderr "[INFO] Please answer 'y' or 'n'." ;;
		esac
	done
}
