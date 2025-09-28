#!/usr/bin/env bash
# polap-bash-apply-pcfg.sh
# Apply PCFG_* environment variables to _arg_* variables used by parse_commandline.
# Usage:
#   source polap-bash-apply-pcfg.sh
#   polap_apply_pcfg            # uses default prefix PCFG
#   polap_apply_pcfg CUSTOM     # uses CUSTOM_* instead of PCFG_*

set -euo pipefail

polap_apply_pcfg() {
	local prefix="${1:-PCFG}" # env var prefix: e.g., PCFG_
	local prefix_="${prefix}_"
	local var

	# Map of _arg names that should trigger a companion _is flag when sourced from env
	# e.g., PCFG_LONG_READS -> _arg_long_reads and _arg_long_reads_is=on
	# Extend this list if you add more *_is companions.
	local -a with_is=(
		long_reads
		short_read1
		short_read2
		coverage
		disassemble_b
		disassemble_c
		disassemble_p
		disassemble_q
		disassemble_n
		genomesize_b
		jellyfish_s
		log
		archive
		select_read_range
		report_x
	)

	# Normalize/guard lookups even under `set -u`
	_pcfg_get() {
		local envname="$1"
		# indirection-safe check
		if [[ -n "${!envname-}" ]]; then
			printf '%s' "${!envname}"
			return 0
		fi
		return 1
	}

	# bool/flag normalizer: on/off for many flags that your parser sets
	_norm_bool() {
		case "${1,,}" in
		1 | true | yes | on) echo "on" ;;
		0 | false | no | off) echo "off" ;;
		*) echo "$1" ;;
		esac
	}

	# helper: check if a name is in with_is array
	_has_is() {
		local name="$1" x
		for x in "${with_is[@]}"; do [[ "$x" == "$name" ]] && return 0; done
		return 1
	}

	# Iterate over all exported variables that start with PREFIX_
	# Note: compgen lists variable NAMES; we then read values via indirect expansion.
	while IFS= read -r var; do
		# var is like PCFG_OUTDIR; compute the _arg_ name
		local key="${var#${prefix_}}" # OUTDIR
		local key_lc
		key_lc="$(tr '[:upper:]' '[:lower:]' <<<"$key")" # outdir
		local arg_name="_arg_${key_lc}"                  # _arg_outdir
		local val
		val="$(_pcfg_get "$var")" || continue

		# Apply special transforms **before** assignment when needed
		case "$key_lc" in
		outdir | anotherdir)
			# strip trailing slash
			val="${val%/}"
			;;
		species)
			# underscores -> space; remove slashes
			val="${val//_/ }"
			val="${val//\//}"
			;;
		random_seed)
			# if <= 0, use RANDOM
			if [[ "$val" =~ ^-?[0-9]+$ ]] && ((val <= 0)); then
				val="$RANDOM"
			fi
			;;
		esac

		# For known boolean/flag args (detected by your parser usage), normalize on/off
		case "$key_lc" in
		reduction_reads | markdown | flye | contigger | polish | dry | all_annotate | log_stderr | polap_reads | \
			bridge_same_strand | coverage_check | clock | timing | nucleotide | plastid | animal | noncoding | yes | \
			resume | circularize | disassemble_align_reference | disassemble_simple_polishing | blast | directional | \
			test | redo)
			val="$(_norm_bool "$val")"
			;;
		esac

		# Assign to the target _arg_* variable (using indirect declare -g for globals)
		# shellcheck disable=SC2086
		declare -g "$arg_name=$val"

		# Handle companions or derived fields
		case "$key_lc" in
		outdir)
			# set _arg_logdir when _arg_outdir is set
			declare -g _arg_logdir="${val}/log"
			;;
		esac

		# If this name is in the "with_is" list, also set its *_is=on
		if _has_is "$key_lc"; then
			# e.g., _arg_long_reads_is
			declare -g "${arg_name}_is=on"
		fi
	done < <(compgen -v | awk -v pfx="${prefix_}" 'index($0, pfx)==1 {print $0}')
}

polap_args_dump_yaml() {
	local out=""
	while (($#)); do
		case "$1" in
		--out)
			out="${2:?}"
			shift 2
			;;
		*)
			echo "Unknown option: $1" >&2
			return 2
			;;
		esac
	done

	mapfile -t __vars < <(compgen -v | awk '/^_arg/ {print}' | sort)

	if [[ -n "$out" ]]; then : >"$out"; fi
	__emit() { if [[ -n "$out" ]]; then printf '%s\n' "$*" >>"$out"; else printf '%s\n' "$*"; fi; }

	local v key key_lc val low
	for v in "${__vars[@]}"; do
		key="${v#_arg_}"
		key_lc="$(tr '[:upper:]' '[:lower:]' <<<"$key")"
		val="${!v-}"
		val="${val#"${val%%[![:space:]]*}"}"
		val="${val%"${val##*[![:space:]]}"}"
		low="${val,,}"

		# Multiline -> block scalar
		if [[ "$val" == *$'\n'* ]]; then
			__emit "${key_lc}: |"
			while IFS= read -r line; do __emit "  ${line}"; done <<<"$val"
			continue
		fi

		# Numbers -> plain
		if [[ "$val" =~ ^-?[0-9]+([.][0-9]+)?$ ]]; then
			__emit "${key_lc}: ${val}"
		# Booleans -> true/false
		elif [[ "$low" =~ ^(on|true|yes|1)$ ]]; then
			__emit "${key_lc}: true"
		elif [[ "$low" =~ ^(off|false|no|0)$ ]]; then
			__emit "${key_lc}: false"
		else
			# Everything else (including 1,2,3,5) stays literal string
			val="${val//\'/\'\'}"
			__emit "${key_lc}: '${val}'"
		fi
	done
}

# # polap_args_dump_yaml: dump all _arg* variables to YAML
# # Usage:
# #   polap_args_dump_yaml --out config.yaml
# #   polap_args_dump_yaml > config.yaml
# polap_args_dump_yaml() {
# 	local out=""
# 	while (($#)); do
# 		case "$1" in
# 		--out)
# 			out="${2:?}"
# 			shift 2
# 			;;
# 		-h | --help)
# 			cat <<EOF
# Usage: polap_args_dump_yaml [--out FILE]
#
# Scans environment for all variables named _arg* and writes them as YAML.
# Strips "_arg_" prefix to form keys. Values are rendered as:
#   - true/false for on/off/yes/no/0/1
#   - numbers without quotes
#   - comma-separated lists -> YAML sequences
#   - multiline strings -> block scalars
#   - other strings quoted safely
# EOF
# 			return 0
# 			;;
# 		*)
# 			echo "Unknown option: $1" >&2
# 			return 2
# 			;;
# 		esac
# 	done
#
# 	mapfile -t __vars < <(compgen -v | awk '/^_arg/ {print}' | sort)
#
# 	if [[ -n "$out" ]]; then : >"$out"; fi
# 	__emit() { if [[ -n "$out" ]]; then printf '%s\n' "$*" >>"$out"; else printf '%s\n' "$*"; fi; }
#
# 	local v key key_lc val low
# 	for v in "${__vars[@]}"; do
# 		key="${v#_arg_}" # strip prefix
# 		key_lc="$(tr '[:upper:]' '[:lower:]' <<<"$key")"
# 		val="${!v-}" # safe expansion under set -u
# 		val="${val#"${val%%[![:space:]]*}"}"
# 		val="${val%"${val##*[![:space:]]}"}"
# 		low="${val,,}"
#
# 		if [[ "$val" == *$'\n'* ]]; then
# 			__emit "${key_lc}: |"
# 			while IFS= read -r line; do __emit "  ${line}"; done <<<"$val"
# 			continue
# 		fi
#
# 		if [[ "$val" == *,* ]]; then
# 			__emit "${key_lc}:"
# 			IFS=',' read -r -a __arr <<<"$val"
# 			for item in "${__arr[@]}"; do
# 				item="${item#"${item%%[![:space:]]*}"}"
# 				item="${item%"${item##*[![:space:]]}"}"
# 				ilow="${item,,}"
# 				if [[ "$item" =~ ^-?[0-9]+([.][0-9]+)?$ ]]; then
# 					__emit "  - $item"
# 				elif [[ "$ilow" =~ ^(on|true|yes|1)$ ]]; then
# 					__emit "  - true"
# 				elif [[ "$ilow" =~ ^(off|false|no|0)$ ]]; then
# 					__emit "  - false"
# 				else
# 					item="${item//\'/\'\'}"
# 					__emit "  - '${item}'"
# 				fi
# 			done
# 			continue
# 		fi
#
# 		if [[ "$val" =~ ^-?[0-9]+([.][0-9]+)?$ ]]; then
# 			__emit "${key_lc}: ${val}"
# 		elif [[ "$low" =~ ^(on|true|yes|1)$ ]]; then
# 			__emit "${key_lc}: true"
# 		elif [[ "$low" =~ ^(off|false|no|0)$ ]]; then
# 			__emit "${key_lc}: false"
# 		else
# 			val="${val//\'/\'\'}"
# 			__emit "${key_lc}: '${val}'"
# 		fi
# 	done
# }
#
# # If this file is sourced directly for a quick test, try applying current env:
# if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
# 	echo "This script is meant to be sourced. Example:"
# 	echo "  export PCFG_OUTDIR=o"
# 	echo "  export PCFG_SPECIES='Arabidopsis_thaliana'"
# 	echo "  source ./polap-bash-apply-pcfg.sh && polap_apply_pcfg"
# fi
