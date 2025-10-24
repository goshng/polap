#!/usr/bin/env bash
################################################################################
# polap-bash-conda-sync.sh
# Version: v0.8.0 (2025-10-18)
# License: GPL-3.0+
#
# Reproducibly manage all polap-* Conda environments **without activating them**.
#
# Commands
#   export     -> write envs/<ENV>/environment.yml  +  locks/<ENV>-<platform>.txt
#   recreate   -> (re)create envs from YAML (portable) or from explicit lock (exact)
#   update     -> apply YAML to existing envs (optional --prune)
#   verify     -> probe tools inside env (via `conda run`, no activate)
#
# Key Features
#   • Default channels: conda-forge, bioconda  (override with --channels)
#   • No activation needed for any operation
#   • Skipping controls for bulk runs:
#       --skip-env NAME     (repeatable)
#       --skip-file FILE    (one NAME per line; '#' comments allowed)
#       --skip-regex REGEX  (exclude by regex)
#   • Recreate helpers:
#       --replace           remove env first, then create
#       --strict-lock       create from locks/<ENV>-<platform>.txt (exact)
#       --skip-existing     silently skip envs that already exist
#
# Environment Variables
#   POLAP_CONDA_SYNC_ROOT   Base folder for envs/ and locks/ (default: $PWD)
#   POLAP_CONDA_PLATFORM    Lockfile platform triplet (default: linux-64)
#   USE_MAMBA               auto|yes|no (default: auto)
#
# Examples
#   bash polap-bash-conda-sync.sh export
#   bash polap-bash-conda-sync.sh recreate --env polap-fmlrc2
#   bash polap-bash-conda-sync.sh recreate --env polap-fmlrc2 --replace --strict-lock
#   bash polap-bash-conda-sync.sh recreate --skip-env polap-dev --skip-existing
#   bash polap-bash-conda-sync.sh update   --env polap-graphaligner --prune
#   bash polap-bash-conda-sync.sh verify   --env polap-polish --tools "python samtools blastn R"
################################################################################
set -euo pipefail
IFS=$'\n\t'

# ---- Config (env-overridable) -----------------------------------------------
ROOT="${POLAP_CONDA_SYNC_ROOT:-$PWD}"
ENV_DIR="${ROOT}/envs"
LOCK_DIR="${ROOT}/locks"
PLATFORM="${POLAP_CONDA_PLATFORM:-linux-64}"
USE_MAMBA="${USE_MAMBA:-auto}" # auto|yes|no
DEFAULT_CHANNELS=("conda-forge" "bioconda")

mkdir -p "$ENV_DIR" "$LOCK_DIR"

# ---- Utilities ---------------------------------------------------------------
_log() { printf '[conda-sync] %s\n' "$*" >&2; }
_die() {
	_log "ERROR: $*"
	exit 2
}

_boot_conda() {
	if ! command -v conda >/dev/null 2>&1; then
		for d in "$HOME/miniconda3" "$HOME/mambaforge" "$HOME/anaconda3" "/opt/conda"; do
			[[ -f "$d/etc/profile.d/conda.sh" ]] && . "$d/etc/profile.d/conda.sh" && break
		done
	fi
	command -v conda >/dev/null 2>&1 || _die "conda not found in PATH"
	# conda hook can reference unset vars; relax nounset temporarily
	set +u
	eval "$(conda shell.bash hook)"
	set -u
}

_solver() {
	case "$USE_MAMBA" in
	yes | auto) command -v mamba >/dev/null 2>&1 && {
		echo mamba
		return
	} ;;
	esac
	echo conda
}

# list polap-* env names known to conda
_list_envs_from_conda() {
	conda info --envs 2>/dev/null |
		awk '{gsub(/\*/,""); if (NF && $1 ~ /^polap(-|$)/) print $1}' || true
}

_env_exists() {
	local env="$1"
	conda info --envs 2>/dev/null | awk '{gsub(/\*/,""); if(NF) print $1}' | grep -qx "$env"
}

_usage() {
	cat <<EOF
Usage:
  $(basename "$0") export   [--env ENV ...] [--full] [--channels conda-forge,bioconda]
  $(basename "$0") recreate [--env ENV ...] [--strict-lock] [--replace] [--skip-existing]
                            [--skip-env NAME ...] [--skip-file FILE] [--skip-regex REGEX]
                            [--channels conda-forge,bioconda]
  $(basename "$0") update   [--env ENV ...] [--prune]
                            [--skip-env NAME ...] [--skip-file FILE] [--skip-regex REGEX]
  $(basename "$0") verify   [--env ENV ...] [--tools "python samtools blastn R"]
                            [--skip-env NAME ...] [--skip-file FILE] [--skip-regex REGEX]

Options:
  --env ENV            Operate on specific env(s). If omitted, target all conda envs named ^polap
  --channels LIST      Comma-separated channels list (default: conda-forge,bioconda)
  --full               (export) write full YAML (no --from-history)
  --strict-lock        (recreate) use locks/<ENV>-<platform>.txt (exact, platform-specific)
  --replace            (recreate) remove env first, then create fresh
  --skip-existing      (recreate) skip envs that already exist
  --prune              (update) remove packages not listed in YAML
  --tools LIST         (verify) tools to probe (default: "python samtools nucmer blastn R")

Skip filters (apply to recreate/update/verify target selection):
  --skip-env NAME      Skip a specific env (repeatable)
  --skip-file FILE     File with env names to skip (one per line; '#' comments ok)
  --skip-regex REGEX   Skip envs whose name matches the regex

Environment:
  POLAP_CONDA_SYNC_ROOT   Base folder for envs/ and locks/ (default: current directory)
  POLAP_CONDA_PLATFORM    Lockfile platform triplet (default: linux-64)
  USE_MAMBA               auto|yes|no (default: auto)
EOF
}

_yaml_ensure_name_and_channels() {
	local yaml="$1" env="$2"
	shift 2 || true
	local -a channels=("$@")

	# Ensure name:
	if ! grep -qE '^name:' "$yaml"; then
		{
			echo "name: ${env}"
			cat "$yaml"
		} >"${yaml}.tmp" && mv "${yaml}.tmp" "$yaml"
	else
		sed -i -E "s/^name:.*/name: ${env}/" "$yaml"
	fi

	# Ensure channels:
	if ! grep -qE '^channels:' "$yaml"; then
		{
			echo "channels:"
			for c in "${channels[@]}"; do echo "  - ${c}"; done
			echo ""
			cat "$yaml"
		} >"${yaml}.tmp" && mv "${yaml}.tmp" "$yaml"
	fi
}

# ---- Parse CLI ---------------------------------------------------------------
cmd="${1:-}"
shift || true
[[ -z "$cmd" ]] && {
	_usage
	exit 1
}

_boot_conda
solver="$(_solver)"

declare -a targets=()
declare -a channels=("${DEFAULT_CHANNELS[@]}")
declare -a skip_envs=()
skip_file=""
skip_regex=""
full_export="false"
strict_lock="false"
replace="false"
skip_existing="false"
prune="false"
tools="python samtools nucmer blastn R"

while (($#)); do
	case "$1" in
	--env)
		targets+=("${2:?}")
		shift 2
		;;
	--channels)
		IFS=',' read -r -a channels <<<"${2:?}"
		shift 2
		;;
	--full)
		full_export="true"
		shift
		;;
	--strict-lock)
		strict_lock="true"
		shift
		;;
	--replace)
		replace="true"
		shift
		;;
	--skip-existing)
		skip_existing="true"
		shift
		;;
	--prune)
		prune="true"
		shift
		;;
	--tools)
		tools="${2:?}"
		shift 2
		;;
	--skip-env)
		skip_envs+=("${2:?}")
		shift 2
		;;
	--skip-file)
		skip_file="${2:?}"
		shift 2
		;;
	--skip-regex)
		skip_regex="${2:?}"
		shift 2
		;;
	-h | --help)
		_usage
		exit 0
		;;
	*) break ;;
	esac
done

# Default target set = all conda envs whose name starts with "polap"
if [[ "${#targets[@]}" -eq 0 ]]; then
	mapfile -t targets < <(_list_envs_from_conda)
fi
[[ "${#targets[@]}" -eq 0 ]] && _die "No polap-* environments found. Use --env to specify."

# Load skip list from file (if any)
if [[ -n "$skip_file" ]]; then
	while IFS= read -r line || [[ -n "$line" ]]; do
		line="${line%%#*}"                      # strip comments
		line="${line#"${line%%[![:space:]]*}"}" # trim leading space
		line="${line%"${line##*[![:space:]]}"}" # trim trailing space
		[[ -n "$line" ]] && skip_envs+=("$line")
	done <"$skip_file"
fi

# Build a set for fast exact-name skips
declare -A _skip_set=()
for s in "${skip_envs[@]}"; do _skip_set["$s"]=1; done

# Filter targets by exact-name and optional regex
declare -a _filtered=()
for e in "${targets[@]}"; do
	if [[ -n "${_skip_set[$e]:-}" ]]; then
		_log "Skip (exact): $e"
		continue
	fi
	if [[ -n "$skip_regex" ]] && [[ "$e" =~ $skip_regex ]]; then
		_log "Skip (regex): $e"
		continue
	fi
	_filtered+=("$e")
done
targets=("${_filtered[@]}")

[[ "${#targets[@]}" -eq 0 ]] && _die "All envs were skipped by the filters."

# ---- Commands ---------------------------------------------------------------
case "$cmd" in
export)
	for env in "${targets[@]}"; do
		[[ -z "$env" ]] && continue
		_log "Exporting ${env}"

		mkdir -p "${ENV_DIR}/${env}"

		# portable YAML (no activation)
		if [[ "$full_export" == "true" ]]; then
			conda env export -n "$env" >"${ENV_DIR}/${env}/environment.yml"
		else
			conda env export -n "$env" --from-history >"${ENV_DIR}/${env}/environment.yml"
		fi
		_yaml_ensure_name_and_channels "${ENV_DIR}/${env}/environment.yml" "$env" "${channels[@]}"

		# exact lock (no activation)
		conda list -n "$env" --explicit --md5 >"${LOCK_DIR}/${env}-${PLATFORM}.txt"

		_log "  -> ${ENV_DIR}/${env}/environment.yml"
		_log "  -> ${LOCK_DIR}/${env}-${PLATFORM}.txt"
	done
	;;

recreate)
	for env in "${targets[@]}"; do
		[[ -z "$env" ]] && continue
		local_yaml="${ENV_DIR}/${env}/environment.yml"
		local_lock="${LOCK_DIR}/${env}-${PLATFORM}.txt"

		# Skip if env already exists (when requested)
		if [[ "$skip_existing" == "true" ]] && _env_exists "$env"; then
			_log "Env '$env' already exists — skipping."
			continue
		fi

		# Replace (remove) existing env first if requested
		if [[ "$replace" == "true" ]]; then
			conda remove -n "$env" --all -y >/dev/null 2>&1 || true
		fi

		# If env still exists and we weren't told to replace, bail early with message
		if _env_exists "$env"; then
			_die "Env '$env' already exists. Use --replace or --skip-existing."
		fi

		if [[ "$strict_lock" == "true" ]]; then
			[[ -f "$local_lock" ]] || _die "Missing lock: $local_lock"
			_log "Recreate (lock): $env"
			"$solver" create -n "$env" -y --file "$local_lock"
		else
			[[ -f "$local_yaml" ]] || _die "Missing YAML: $local_yaml"
			_yaml_ensure_name_and_channels "$local_yaml" "$env" "${channels[@]}"
			_log "Recreate (yaml): $env"
			"$solver" env create -n "$env" -f "$local_yaml"
		fi

		_log "Recreated: $env"
	done
	;;

update)
	# (fall-through permitted if user passed 'recreate' then wants update for remaining targets;
	# in practice, bash doesn't do switch fall-through—above ';|' is no-op stylistically. Keep distinct.)
	if [[ "$cmd" == "update" ]]; then :; fi
	if [[ "$cmd" == "update" ]]; then
		:
	fi
	if [[ "$cmd" == "update" ]]; then
		for env in "${targets[@]}"; do
			[[ -z "$env" ]] && continue
			local_yaml="${ENV_DIR}/${env}/environment.yml"
			[[ -f "$local_yaml" ]] || _die "Missing YAML: $local_yaml"
			_yaml_ensure_name_and_channels "$local_yaml" "$env" "${channels[@]}"

			_log "Update (yaml): $env"
			if [[ "$prune" == "true" ]]; then
				conda env update -n "$env" -f "$local_yaml" --prune
			else
				conda env update -n "$env" -f "$local_yaml"
			fi
			_log "Updated: $env"
		done
	fi
	;;

verify)
	for env in "${targets[@]}"; do
		_log "Verify $env"
		# probe tools using conda run (never activate)
		for t in $tools; do
			if conda run -n "$env" bash -lc "command -v $t >/dev/null 2>&1"; then
				p="$(conda run -n "$env" bash -lc "command -v $t")"
				printf "  ✔ %-12s -> %s\n" "$t" "$p"
			else
				printf "  ✘ %-12s (missing)\n" "$t"
			fi
		done
	done
	;;

*)
	_usage
	exit 1
	;;
esac

exit 0
