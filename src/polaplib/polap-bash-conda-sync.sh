#!/usr/bin/env bash
################################################################################
# polap-bash-conda-sync.sh
# Version: v0.7.0 (2025-10-18)
# License: GPL-3.0+
#
# Reproducibly manage all polap-* Conda environments **without activating them**.
#
# Features
#   • export   -> envs/<ENV>/environment.yml  +  locks/<ENV>-<platform>.txt
#   • recreate -> from YAML (portable) or from explicit lock (exact, platform-pinned)
#   • update   -> apply YAML to existing envs (optional --prune)
#   • verify   -> probe tools inside env (via `conda run`, no activate)
#
# Defaults
#   • channels:   conda-forge, bioconda
#   • platform:   linux-64     (override: export POLAP_CONDA_PLATFORM=linux-aarch64)
#   • solver:     mamba if present or USE_MAMBA=yes/auto, else conda
#
# Examples
#   ./polap-bash-conda-sync.sh export
#   ./polap-bash-conda-sync.sh recreate --env polap-fmlrc2
#   ./polap-bash-conda-sync.sh recreate --env polap-fmlrc2 --strict-lock --replace
#   ./polap-bash-conda-sync.sh update   --env polap-graphaligner --prune
#   ./polap-bash-conda-sync.sh verify   --env polap-polish --tools "python samtools blastn R"
#
# Notes
#   • Uses `conda env export -n <ENV>` and `conda list -n <ENV> --explicit --md5`
#     so broken activate.d hooks cannot break exports.
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
	# no `set -u` crash: conda hook may reference unset vars
	set +u
	eval "$(conda shell.bash hook)"
	set -u
}

_solver() {
	case "$USE_MAMBA" in
	yes) command -v mamba >/dev/null 2>&1 && {
		echo mamba
		return
	} ;;
	auto) command -v mamba >/dev/null 2>&1 && {
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

_usage() {
	cat <<EOF
Usage:
  $(basename "$0") export   [--env ENV ...] [--full]
  $(basename "$0") recreate [--env ENV ...] [--strict-lock] [--replace]
  $(basename "$0") update   [--env ENV ...] [--prune]
  $(basename "$0") verify   [--env ENV ...] [--tools "python samtools blastn R"]

Options:
  --env ENV          Operate on specific env(s). If omitted, target all conda envs named ^polap
  --full             Export full YAML (no --from-history)
  --strict-lock      Recreate from explicit lock (exact, platform-specific)
  --replace          Remove env first (fresh rebuild)
  --prune            Remove packages not listed in YAML during update
  --tools "<list>"   Tools to probe in verify (default: "python samtools nucmer blastn R")

Environment variables:
  POLAP_CONDA_SYNC_ROOT   Base folder for envs/ and locks/ (default: \$PWD)
  POLAP_CONDA_PLATFORM    Lockfile platform triplet (default: linux-64)
  USE_MAMBA               auto|yes|no (default: auto)

Defaults:
  channels: conda-forge, bioconda
EOF
}

# ensure YAML has `name: <env>` and includes `channels:` block if missing
_yaml_ensure_name_and_channels() {
	local yaml="$1" env="$2"
	shift 2 || true
	local -a channels=("$@")

	# ensure name:
	if ! grep -qE '^name:' "$yaml"; then
		# prepend name at top
		{
			echo "name: ${env}"
			cat "$yaml"
		} >"${yaml}.tmp" && mv "${yaml}.tmp" "$yaml"
	else
		# replace any existing name with target env
		# (keep idempotent if already same)
		sed -i -E "s/^name:.*/name: ${env}/" "$yaml"
	fi

	# ensure channels:
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
full_export="false"
strict_lock="false"
replace="false"
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
		;; # optional override
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
	--prune)
		prune="true"
		shift
		;;
	--tools)
		tools="${2:?}"
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

		_log "  → ${ENV_DIR}/${env}/environment.yml"
		_log "  → ${LOCK_DIR}/${env}-${PLATFORM}.txt"
	done
	;;

recreate)
	for env in "${targets[@]}"; do
		[[ -z "$env" ]] && continue
		local_yaml="${ENV_DIR}/${env}/environment.yml"
		local_lock="${LOCK_DIR}/${env}-${PLATFORM}.txt"

		if [[ "$replace" == "true" ]]; then
			conda remove -n "$env" --all -y >/dev/null 2>&1 || true
		fi

		if [[ "$strict_lock" == "true" ]]; then
			[[ -f "$local_lock" ]] || _die "Missing lock: $local_lock"
			_log "Recreate (lock): $env"
			"$solver" create -n "$env" -y --file "$local_lock"
		else
			[[ -f "$local_yaml" ]] || _die "Missing YAML: $local_yaml"
			# Ensure name matches env before create (some solvers read from YAML)
			_yaml_ensure_name_and_channels "$local_yaml" "$env" "${channels[@]}"
			_log "Recreate (yaml): $env"
			"$solver" env create -n "$env" -f "$local_yaml"
		fi

		_log "Recreated: $env"
	done
	;;

update)
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
	;;

verify)
	for env in "${targets[@]}"; do
		_log "Verify $env"
		# probe tools using conda run (never activate)
		for t in $tools; do
			if conda run -n "$env" bash -lc "command -v $t >/dev/null 2>&1"; then
				# show resolved path
				p="$(conda run -n "$env" bash -lc "command -v $t")"
				printf "  ✔ %-10s -> %s\n" "$t" "$p"
			else
				printf "  ✘ %-10s (missing)\n" "$t"
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
