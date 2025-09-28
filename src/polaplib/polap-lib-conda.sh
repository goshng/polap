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
		conda activate "$env_name"
	fi

	# Verify activation
	if [[ "${CONDA_DEFAULT_ENV:-}" != "$env_name" ]]; then
		[[ "$_POLAP_DEBUG" == "1" ]] && echo "[ERROR] Failed to enter conda environment '$env_name'"
		return 1
	fi

	return 0
}

# _polap_lib_conda-ensure_conda_env() {
# 	local env_name="$1"
#
# 	if [[ -z "$env_name" ]]; then
# 		if [[ "$_POLAP_DEBUG" == "1" ]]; then
# 			echo "[ERROR] No conda environment name provided"
# 		fi
# 		return 1
# 	fi
#
# 	# initialize conda and activate the base conda environment
# 	source "$(conda info --base)/etc/profile.d/conda.sh"
#
# 	# activate the conda environment provided by the input
# 	if [[ "${CONDA_DEFAULT_ENV:-}" != "$env_name" ]]; then
# 		if [[ "$_POLAP_DEBUG" == "1" ]]; then
# 			echo "[INFO] Activating conda environment '$env_name'..."
# 		fi
# 		conda activate "$env_name"
# 	fi
#
# 	# check the activation of the conda environment
# 	if [[ "${CONDA_DEFAULT_ENV:-}" != "$env_name" ]]; then
# 		if [[ "$_POLAP_DEBUG" == "1" ]]; then
# 			echo "[ERROR] Failed to enter conda environment '$env_name'"
# 		fi
# 		return 1
# 	fi
# }
