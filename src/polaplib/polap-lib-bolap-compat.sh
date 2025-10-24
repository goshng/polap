#!/usr/bin/env bash
# polap-lib-bolap-compat.sh
# Version: v0.2.0
# SPDX-License-Identifier: GPL-3.0-or-later
set -euo pipefail

# Rewrite legacy invocations into new flag style.
# Usage: bolap_compat_rewrite_menu _brg_menu
#
# Supported legacy forms (head ∈ {list, find, search}):
#   1)  bolap list data                 -> list --query data --where any
#   2)  bolap list data start           -> list --query data --where start
#   3)  bolap list '^read' end -x 1     -> list --query '^read' --where end -x 1
#   4)  bolap find foo ...              -> find --query foo --where any ...
#   5)  bolap search bar start ...      -> search --query bar --where start ...
#
# Notes:
# - If second positional matches {any,start,end}, it is treated as --where.
# - Only rewrites when the 2nd token is a non-option (does not start with '-').
# - Leaves already-flagged invocations untouched.
bolap_compat_rewrite_menu() {
	local __arrname="$1"
	# Safe copy by indirection (array may be empty)
	local -a __tmp=()
	eval "__tmp=( \"\${${__arrname}[@]:-}\" )"

	[[ ${#__tmp[@]} -eq 0 ]] && return 0

	local head="${__tmp[0]}"
	case "$head" in
	list | find | search)
		# already flag style? (token 2 starts with '-') -> do nothing
		if [[ ${#__tmp[@]} -ge 2 && "${__tmp[1]}" != -* ]]; then
			local q="${__tmp[1]}"
			local where="any"
			local rest_start=2

			if [[ ${#__tmp[@]} -ge 3 && "${__tmp[2]}" != -* ]]; then
				case "${__tmp[2]}" in
				any | start | end)
					where="${__tmp[2]}"
					rest_start=3
					;;
				*) : ;; # treat as rest (unknown positional)
				esac
			fi

			# Rebuild: head --query q --where where [rest...]
			local -a rebuilt=("$head" --query "$q" --where "$where")
			if [[ ${#__tmp[@]} -gt $rest_start ]]; then
				rebuilt+=("${__tmp[@]:$rest_start}")
			fi

			# Write back
			eval "$__arrname=()"
			eval "$__arrname+=( \"\${rebuilt[@]}\" )"
		fi
		;;
	*) : ;; # other heads untouched
	esac
}
