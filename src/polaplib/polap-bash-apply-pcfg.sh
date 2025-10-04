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

	# ---------- helpers ----------
	_pcfg_get() { # safe env lookup under set -u
		local envname="$1"
		[[ -n "${!envname-}" ]] || return 1
		printf '%s' "${!envname}"
	}

	_trim() {
		local s="$*"
		s="${s#"${s%%[![:space:]]*}"}"
		printf '%s' "${s%"${s##*[![:space:]]}"}"
	}

	# bool normalizer: on/off
	_norm_bool() {
		case "$(_trim "${1,,}")" in
		1 | true | t | yes | y | on) echo "on" ;;
		0 | false | f | no | n | off) echo "off" ;;
		*) echo "$1" ;;
		esac
	}

	# test membership in a bash array: _in_array <needle> <array...>
	_in_array() {
		local needle="$1"
		shift
		local x
		for x in "$@"; do [[ "$x" == "$needle" ]] && return 0; done
		return 1
	}

	# looks-like-bool? (respects assume_bool_01 for 0/1)
	_looks_bool() {
		local v="$(_trim "${1,,}")"
		case "$v" in
		true | t | false | f | yes | y | no | n | on | off) return 0 ;;
		0 | 1) [[ "$assume_bool_01" == "on" ]] && return 0 || return 1 ;;
		*) return 1 ;;
		esac
	}

	_has_is() {
		local name="$1" x
		for x in "${with_is[@]}"; do [[ "$x" == "$name" ]] && return 0; done
		return 1
	}

	# Keys that should also get a companion *_is=on when present
	local -a with_is=(
		long_reads short_read1 short_read2 coverage
		disassemble_b disassemble_c disassemble_p disassemble_q disassemble_n
		genomesize_b jellyfish_s log archive select_read_range report_x
	)

	# Numeric-leaning keys: never auto-boolify even if the value is 0/1
	# (add here if you see misclassification)
	local -a numeric_prefer=(coverage genomesize_b jellyfish_s random_seed)

	# Control: treat bare 0/1 as booleans?
	#   export PCFG_BOOLIFY_01=on    # opt in
	#   (default: off)
	local assume_bool_01
	assume_bool_01="$(_pcfg_get "${prefix_}BOOLIFY_01" && echo "${!REPLY}" || echo "off")"
	assume_bool_01="$(_norm_bool "$assume_bool_01")"

	# ---------- main ----------
	while IFS= read -r var; do
		local key="${var#${prefix_}}" # e.g., OUTDIR
		local key_lc
		key_lc="$(tr '[:upper:]' '[:lower:]' <<<"$key")"
		local arg_name="_arg_${key_lc}"
		local val
		val="$(_pcfg_get "$var")" || continue

		# Pre-assignment transforms
		case "$key_lc" in
		outdir | anotherdir) val="${val%/}" ;;
		species)
			val="${val//_/ }"
			val="${val//\//}"
			;;
		random_seed)
			if [[ "$val" =~ ^-?[0-9]+$ ]] && ((val <= 0)); then val="$RANDOM"; fi
			;;
		esac

		# If the value *looks* boolean (true/false/on/off/yes/no and optionally 0/1),
		# and the key isn't in the numeric-prefer list, normalize to on/off.
		if _looks_bool "$val" && ! _in_array "$key_lc" "${numeric_prefer[@]}"; then
			val="$(_norm_bool "$val")"
		fi

		# (Optional) explicit allowlist of known boolean args — keeps past behavior
		case "$key_lc" in
		reduction_reads | markdown | flye | contigger | polish | dry | all_annotate | log_stderr | polap_reads | \
			bridge_same_strand | coverage_check | clock | timing | nucleotide | plastid | animal | noncoding | yes | \
			resume | circularize | disassemble_align_reference | disassemble_simple_polishing | blast | directional | \
			test | redo)
			val="$(_norm_bool "$val")"
			;;
		esac

		# Assign
		declare -g "$arg_name=$val"

		# Deriveds
		case "$key_lc" in
		outdir) declare -g _arg_logdir="${val}/log" ;;
		esac

		# *_is companion
		if _has_is "$key_lc"; then
			declare -g "${arg_name}_is=on"
		fi
	done < <(compgen -v | awk -v pfx="${prefix_}" 'index($0, pfx)==1 {print $0}')
}

# polap_apply_pcfg() {
# 	local prefix="${1:-PCFG}" # env var prefix: e.g., PCFG_
# 	local prefix_="${prefix}_"
# 	local var
#
# 	# Map of _arg names that should trigger a companion _is flag when sourced from env
# 	# e.g., PCFG_LONG_READS -> _arg_long_reads and _arg_long_reads_is=on
# 	# Extend this list if you add more *_is companions.
# 	local -a with_is=(
# 		long_reads
# 		short_read1
# 		short_read2
# 		coverage
# 		disassemble_b
# 		disassemble_c
# 		disassemble_p
# 		disassemble_q
# 		disassemble_n
# 		genomesize_b
# 		jellyfish_s
# 		log
# 		archive
# 		select_read_range
# 		report_x
# 	)
#
# 	# Normalize/guard lookups even under `set -u`
# 	_pcfg_get() {
# 		local envname="$1"
# 		# indirection-safe check
# 		if [[ -n "${!envname-}" ]]; then
# 			printf '%s' "${!envname}"
# 			return 0
# 		fi
# 		return 1
# 	}
#
# 	# bool/flag normalizer: on/off for many flags that your parser sets
# 	_norm_bool() {
# 		case "${1,,}" in
# 		1 | true | yes | on) echo "on" ;;
# 		0 | false | no | off) echo "off" ;;
# 		*) echo "$1" ;;
# 		esac
# 	}
#
# 	# helper: check if a name is in with_is array
# 	_has_is() {
# 		local name="$1" x
# 		for x in "${with_is[@]}"; do [[ "$x" == "$name" ]] && return 0; done
# 		return 1
# 	}
#
# 	# Iterate over all exported variables that start with PREFIX_
# 	# Note: compgen lists variable NAMES; we then read values via indirect expansion.
# 	while IFS= read -r var; do
# 		# var is like PCFG_OUTDIR; compute the _arg_ name
# 		local key="${var#${prefix_}}" # OUTDIR
# 		local key_lc
# 		key_lc="$(tr '[:upper:]' '[:lower:]' <<<"$key")" # outdir
# 		local arg_name="_arg_${key_lc}"                  # _arg_outdir
# 		local val
# 		val="$(_pcfg_get "$var")" || continue
#
# 		# Apply special transforms **before** assignment when needed
# 		case "$key_lc" in
# 		outdir | anotherdir)
# 			# strip trailing slash
# 			val="${val%/}"
# 			;;
# 		species)
# 			# underscores -> space; remove slashes
# 			val="${val//_/ }"
# 			val="${val//\//}"
# 			;;
# 		random_seed)
# 			# if <= 0, use RANDOM
# 			if [[ "$val" =~ ^-?[0-9]+$ ]] && ((val <= 0)); then
# 				val="$RANDOM"
# 			fi
# 			;;
# 		esac
#
# 		# For known boolean/flag args (detected by your parser usage), normalize on/off
# 		case "$key_lc" in
# 		reduction_reads | markdown | flye | contigger | polish | dry | all_annotate | log_stderr | polap_reads | \
# 			bridge_same_strand | coverage_check | clock | timing | nucleotide | plastid | animal | noncoding | yes | \
# 			resume | circularize | disassemble_align_reference | disassemble_simple_polishing | blast | directional | \
# 			test | redo)
# 			val="$(_norm_bool "$val")"
# 			;;
# 		esac
#
# 		# Assign to the target _arg_* variable (using indirect declare -g for globals)
# 		# shellcheck disable=SC2086
# 		declare -g "$arg_name=$val"
#
# 		# Handle companions or derived fields
# 		case "$key_lc" in
# 		outdir)
# 			# set _arg_logdir when _arg_outdir is set
# 			declare -g _arg_logdir="${val}/log"
# 			;;
# 		esac
#
# 		# If this name is in the "with_is" list, also set its *_is=on
# 		if _has_is "$key_lc"; then
# 			# e.g., _arg_long_reads_is
# 			declare -g "${arg_name}_is=on"
# 		fi
# 	done < <(compgen -v | awk -v pfx="${prefix_}" 'index($0, pfx)==1 {print $0}')
# }

# polap_args_dump_yaml() {
# 	local out=""
# 	while (($#)); do
# 		case "$1" in
# 		--out)
# 			out="${2:?}"
# 			shift 2
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
# 		key="${v#_arg_}"
# 		key_lc="$(tr '[:upper:]' '[:lower:]' <<<"$key")"
# 		val="${!v-}"
# 		val="${val#"${val%%[![:space:]]*}"}"
# 		val="${val%"${val##*[![:space:]]}"}"
# 		low="${val,,}"
#
# 		# Multiline -> block scalar
# 		if [[ "$val" == *$'\n'* ]]; then
# 			__emit "${key_lc}: |"
# 			while IFS= read -r line; do __emit "  ${line}"; done <<<"$val"
# 			continue
# 		fi
#
# 		# Numbers -> plain
# 		if [[ "$val" =~ ^-?[0-9]+([.][0-9]+)?$ ]]; then
# 			__emit "${key_lc}: ${val}"
# 		# Booleans -> true/false
# 		elif [[ "$low" =~ ^(on|true|yes|1)$ ]]; then
# 			__emit "${key_lc}: true"
# 		elif [[ "$low" =~ ^(off|false|no|0)$ ]]; then
# 			__emit "${key_lc}: false"
# 		else
# 			# Everything else (including 1,2,3,5) stays literal string
# 			val="${val//\'/\'\'}"
# 			__emit "${key_lc}: '${val}'"
# 		fi
# 	done
# }

polap_args_dump_yaml() {
	# --- default skip list (match lowercased _arg_* suffixes) ---
	# Add more keys here to always exclude by default.
	local -a _default_skip=(
		plastid plastid_is animal noncoding use_oatk verbose help
		dry log_is log_stderr coverage_is nucleotide markdown
		flye directional blast reduction_reads contigger polish
		all_annotate ploap_reads bridge_same_strand resume clock
		yes timing circularize redo test debug downsample_is disassemble_p_is
		disassemble_n_is disassemble_c_is
		disassemble_align_reference
		disassemble_simple_polishing
		genomesize_b_is
		steps_is stages_is
		jellyfish_s_is
		data_type
		flye_data_type
		minimap2_data_type
		oatk_no_read_ec
		oatk_no_hpc
		parallel
		log
		logdir
	)

	local out=""
	local -a _skip_keys=()        # will be populated below
	local _use_default_skip="yes" # can be turned off by --no-default-skip

	# ---- options ----
	while (($#)); do
		case "$1" in
		--out)
			out="${2:?}"
			shift 2
			;;
		--skip)
			_skip_keys+=("$(tr '[:upper:]' '[:lower:]' <<<"${2:?}")")
			shift 2
			;;
		--skip-csv)
			IFS=, read -r -a __tmp <<<"${2:?}"
			for __k in "${__tmp[@]}"; do
				__k="$(tr '[:upper:]' '[:lower:]' <<<"${__k}")"
				[[ -n "$__k" ]] && _skip_keys+=("$__k")
			done
			unset __tmp __k
			shift 2
			;;
		--no-default-skip)
			_use_default_skip="no"
			shift
			;;
		*)
			echo "Unknown option: $1" >&2
			return 2
			;;
		esac
	done

	# Merge default skip list unless disabled
	if [[ "$_use_default_skip" == "yes" ]]; then
		_skip_keys+=("${_default_skip[@]}")
	fi

	# Deduplicate skip list
	if ((${#_skip_keys[@]})); then
		mapfile -t _skip_keys < <(printf '%s\n' "${_skip_keys[@]}" | awk 'NF' | tr '[:upper:]' '[:lower:]' | sort -u)
	fi

	# helper: membership test
	__in_skip() {
		local needle="$1" k
		for k in "${_skip_keys[@]}"; do
			[[ "$k" == "$needle" ]] && return 0
		done
		return 1
	}

	mapfile -t __vars < <(compgen -v | awk '/^_arg/ {print}' | sort)

	# output init
	if [[ -n "$out" ]]; then : >"$out"; fi
	__emit() { if [[ -n "$out" ]]; then printf '%s\n' "$*" >>"$out"; else printf '%s\n' "$*"; fi; }

	local v key key_lc val low line
	for v in "${__vars[@]}"; do
		key="${v#_arg_}"
		key_lc="$(tr '[:upper:]' '[:lower:]' <<<"$key")"

		# Skip requested (or default) keys
		if __in_skip "$key_lc"; then
			continue
		fi

		val="${!v-}"
		# trim
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
			# Everything else as quoted string (escape single quotes)
			val="${val//\'/\'\'}"
			__emit "${key_lc}: '${val}'"
		fi
	done
}
