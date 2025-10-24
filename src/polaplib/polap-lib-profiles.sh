#!/usr/bin/env bash
# polaplib/polap-lib-profiles.sh
# Version: v0.1.0
# SPDX-License-Identifier: GPL-3.0-or-later
# Helpers to ensure ~/.polap/profiles exists and is initialized.

set -euo pipefail

# Print to stderr
_polap_echo_err() { printf '%s\n' "$*" >&2; }

# Create ~/.polap/profiles skeleton
_polap_init_profiles() {
	local base="${1:?}" prof="${2:?}"
	mkdir -p "${prof}"

	# Default profile content (edit to taste)
	local def="${prof}/default.env"
	if [[ ! -s "$def" ]]; then
		cat >"$def" <<'ENV'
# ~/.polap/profiles/default.env
# Version: v0.1.0
# Customize these for your site; sourced by polap/bolap at startup (if you want).
# Examples:
# export _LOCAL_HOST="thorne"
# export _MEDIA_DIR="/media/h2/sra"
# export _MEDIA1_DIR="/media/h1/sra"
# export _MEDIA2_DIR="/media/h2/sra"
# export _POLAP_RELEASE=0
# export _POLAP_DEBUG=0
ENV
	fi

	# Point 'current' to default if not present
	if [[ ! -e "${prof}/current" ]]; then
		ln -sfn "default.env" "${prof}/current"
	fi

	# Optional README
	if [[ ! -s "${base}/README.txt" ]]; then
		cat >"${base}/README.txt" <<'TXT'
This directory holds your Polap profiles.

- profiles/default.env  : default environment for polap/bolap (bash syntax)
- profiles/current      : symlink to the active profile (sourced by your wrapper if you wire it)
To add a new profile, copy default.env -> mylab.env and retarget 'current' to it.
TXT
	fi

	_polap_echo_err "[INFO] Created skeleton: ${prof} (default.env, current)"
}

# Ensure ~/.polap/profiles; prompt if interactive
_polap_ensure_profiles_dir() {
	local base="${HOME}/.polap"
	local prof="${base}/profiles"

	# Already present? nothing to do.
	[[ -d "$prof" ]] && return 0

	# Non-interactive: create only if asked via env
	if [[ ! -t 0 && "${POLAP_AUTO_INIT:-0}" != "1" ]]; then
		_polap_echo_err "[WARN] Missing ${prof}. Non-interactive; set POLAP_AUTO_INIT=1 to auto-create."
		return 0
	fi

	# Interactive: ask user
	_polap_echo_err "[INFO] Missing ${prof}."
	local ans=""
	if [[ -t 0 ]]; then
		printf "Create %s now? [Y/n] " "$prof" >/dev/tty
		IFS= read -r ans </dev/tty || ans=""
	fi
	case "${ans,,}" in
	"" | y | yes)
		_polap_init_profiles "$base" "$prof"
		;;
	n | no)
		_polap_echo_err "[INFO] Skipping profile creation (you can run again later)."
		;;
	*)
		_polap_echo_err "[INFO] Unrecognized answer; skipping."
		;;
	esac
}
