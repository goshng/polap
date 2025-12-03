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
# Polap's test scripts such as polap-data-cflye or polap-data-v2.sh call
# other tools, which are installed in their own conda environments.
# Before executing such tools, we make sure that such execution happens
# in its conda environment by 'conda activate'.
# Activating a conda environment need a little more than conda deactivate
# using the following function:
# _polap_lib_conda-ensure_conda_env
#
# Example:
# _polap_lib_conda-ensure_conda_env polap-fmlrc
# <your execution>
# conda deactivate
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

_polap_lib_conda-ensure_conda_env() {
	local env_name="$1"

	if [[ -z "$env_name" ]]; then
		[[ "$_POLAP_DEBUG" == "1" ]] && echo "[ERROR] No conda environment name provided"
		return 1
	fi

	# Try to find conda.sh
	local conda_sh=""
	if command -v conda >/dev/null 2>&1; then
		conda_sh="$(conda info --base)/etc/profile.d/conda.sh"
	elif [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
		conda_sh="$HOME/miniconda3/etc/profile.d/conda.sh"
	elif [[ -f "$HOME/etc/profile.d/conda.sh" ]]; then
		conda_sh="$HOME/etc/profile.d/conda.sh"
	fi

	if [[ -z "$conda_sh" || ! -f "$conda_sh" ]]; then
		[[ "$_POLAP_DEBUG" == "1" ]] && echo "[ERROR] Cannot find conda.sh"
		return 1
	fi

	# Source conda.sh to enable 'conda' command
	# shellcheck disable=SC1090
	source "$conda_sh"

	# Activate the environment
	if [[ "${CONDA_DEFAULT_ENV:-}" != "$env_name" ]]; then
		[[ "$_POLAP_DEBUG" == "1" ]] && echo "[INFO] Activating conda environment '$env_name'..."
		set +u
		conda activate "$env_name"
		set -u
	fi

	# Verify activation
	if [[ "${CONDA_DEFAULT_ENV:-}" != "$env_name" ]]; then
		[[ "$_POLAP_DEBUG" == "1" ]] && echo "[ERROR] Failed to enter conda environment '$env_name'"
		return 1
	fi

	return 0
}

# Robust, generic Conda env installer (safe with: set -euo pipefail)
# Usage:
#   ensure_conda_env ENV_NAME [PKG ...]
#   ensure_conda_env ENV_NAME [PKG ...] --python 3.11 -y --channel conda-forge --channel bioconda --replace --no-mamba
#
# Flags:
#   -y | --yes              Non-interactive (skip prompt). Also honored: opt_y_flag=true
#   --python VER            Pin Python (e.g., 3.11)
#   --channel NAME          Add a Conda channel (repeatable). Defaults: conda-forge, bioconda
#   --no-default-channels   Don’t prepend default channels
#   --replace               If env exists, remove it and re-create (instead of install/upgrade)
#   --no-mamba              Don’t use mamba even if available
#   --quiet                 Less solver output
# Notes:
#   • Prefers mamba for speed if available (for create/install), uses conda to activate.
#   • If the env exists and --replace is NOT given, this will “install/upgrade” the requested packages into it.
#   • Works without being in “base”; no stacking: it switches directly.

# --- helpers ---------------------------------------------------------------

_polap_lib_conda-bootstrap() {
	if ! command -v conda >/dev/null 2>&1; then
		# Try common install locations
		for d in "$HOME/miniconda3" "$HOME/mambaforge" "$HOME/anaconda3" "/opt/conda"; do
			if [[ -f "$d/etc/profile.d/conda.sh" ]]; then
				# shellcheck disable=SC1090
				source "$d/etc/profile.d/conda.sh"
				break
			fi
		done
	fi
	command -v conda >/dev/null 2>&1 || {
		echo "ERROR: conda not found." >&2
		return 1
	}
	# Enable `conda activate` in non-interactive shells
	eval "$(conda shell.bash hook)" || {
		echo "ERROR: Failed to init conda hook." >&2
		return 1
	}
}

_polap_lib_conda-current-env() {
	if [[ -n "${CONDA_PREFIX-}" ]]; then
		basename -- "$CONDA_PREFIX"
	elif [[ -n "${CONDA_DEFAULT_ENV-}" ]]; then
		printf '%s' "$CONDA_DEFAULT_ENV"
	else printf ''; fi
}

_polap_lib_conda-env-exists() {
	local env="$1"
	conda info --envs 2>/dev/null |
		awk '{gsub(/\*/,""); if (NF) print $1}' |
		grep -Fxq -- "$env"
}

# --- main ------------------------------------------------------------------

_polap_lib_conda-create-env() {
	local env_name
	local python_ver=""
	local auto_yes="true"
	local replace="false"
	local quiet="false"
	local use_mamba="false"
	local use_default_channels="true"
	local -a channels=() pkgs=() create_args=() install_args=() solver_cmd=()

	# honor external boolean if present
	if [[ "${opt_y_flag-}" == "true" ]]; then auto_yes="true"; fi

	# parse args
	if [[ $# -lt 1 ]]; then
		echo "Usage: ensure_conda_env ENV_NAME [PKG ...] [flags]" >&2
		return 2
	fi
	env_name="$1"
	shift || true
	while (($#)); do
		case "$1" in
		--python)
			python_ver="${2-}"
			shift 2
			;;
		--no)
			auto_yes="false"
			shift
			;;
		-y | --yes)
			auto_yes="true"
			shift
			;;
		--replace)
			replace="true"
			shift
			;;
		--channel)
			channels+=("${2-}")
			shift 2
			;;
		--no-default-channels)
			use_default_channels="false"
			shift
			;;
		--no-mamba)
			use_mamba="false"
			shift
			;;
		--quiet)
			quiet="true"
			shift
			;;
		--)
			shift
			break
			;;
		-*)
			echo "ERROR: Unknown flag: $1" >&2
			return 2
			;;
		*)
			pkgs+=("$1")
			shift
			;;
		esac
	done
	# any tail as pkgs
	while (($#)); do
		pkgs+=("$1")
		shift
	done

	# sensible defaults
	if [[ "${#channels[@]}" -eq 0 && "$use_default_channels" == "true" ]]; then
		channels=(conda-forge bioconda)
	fi
	if [[ -n "$python_ver" ]]; then
		pkgs=("python=${python_ver}" "${pkgs[@]}")
	fi

	# prompt unless auto-yes
	if [[ "$auto_yes" != "true" ]]; then
		local reply
		IFS= read -r -p "Create/update Conda env '${env_name}' with packages: ${pkgs[*]:-(none)} ? (y/N): " reply || true
		reply="${reply,,}"
		if [[ "$reply" != "y" && "$reply" != "yes" ]]; then
			echo "Canceled."
			return 0
		fi
	fi

	# conda bootstrap
	_polap_lib_conda-bootstrap || return 1

	# choose solver
	if [[ "$use_mamba" == "true" ]] && command -v mamba >/dev/null 2>&1; then
		solver_cmd=(mamba)
	else
		solver_cmd=(conda)
	fi

	# channel flags
	local -a chan_flags=()
	if [[ "$use_default_channels" == "false" ]]; then
		chan_flags+=(--override-channels)
	fi
	for c in "${channels[@]}"; do
		chan_flags+=(-c "$c")
	done

	# quiet?
	local -a quiet_flag=()
	if [[ "$quiet" == "true" ]]; then quiet_flag+=(--quiet); fi

	# create or update
	if _polap_lib_conda-env-exists "$env_name"; then
		if [[ "$replace" == "true" ]]; then
			echo "Removing existing env '$env_name'..."
			conda remove -n "$env_name" --all -y || return 1
			echo "Creating env '$env_name'..."
			"${solver_cmd[@]}" "${quiet_flag[@]}" create -y -n "$env_name" "${chan_flags[@]}" "${pkgs[@]}" || return 1
		else
			echo "Env '$env_name' exists; installing/upgrading requested packages..."
			if [[ "${#pkgs[@]}" -gt 0 ]]; then
				"${solver_cmd[@]}" "${quiet_flag[@]}" install -y -n "$env_name" "${chan_flags[@]}" "${pkgs[@]}" || return 1
			else
				echo "No packages specified; leaving packages as-is."
			fi
		fi
	else
		echo "Creating env '$env_name'..."
		"${solver_cmd[@]}" "${quiet_flag[@]}" create -y -n "$env_name" "${chan_flags[@]}" "${pkgs[@]}" || return 1
	fi

	# activate & report
	conda activate "$env_name" || {
		echo "ERROR: failed to activate '$env_name'." >&2
		return 1
	}
	echo "Activated env: $(_polap_lib_conda-current-env)"

	# quick PATH sanity for a few common tools if present in spec
	# local -a probe=(python R samtools minimap2 blastn)
	local -a probe=(python)
	local t
	for t in "${probe[@]}"; do
		if command -v "$t" >/dev/null 2>&1; then
			printf "  ✔ %s -> %s\n" "$t" "$(command -v "$t")"
		fi
	done
}

# -------------------- Examples --------------------
# pkgs=(python=3.11 samtools mummer4 blast minimap2 biopython pysam r-base r-data.table r-ggplot2 r-gridextra)
# _polap_lib_conda-create-env myenv pip numpy pandas
# _polap_lib_conda-create-env polap-evo -y --channel conda-forge --channel bioconda "${pkgs[@]}"
# _polap_lib_conda-create-env myenv -y --python 3.11 --channel conda-forge pip numpy pandas
# _polap_lib_conda-create-env scratch --replace -y --no-default-channels --channel bioconda fastp bwa
