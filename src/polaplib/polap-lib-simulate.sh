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

_polap_lib_simulate-ont-with-reference() {
	############################################
	# defaults
	############################################
	SIM_MT_REF=""
	SIM_PT_REF=""
	SIM_NUC_REF=""
	SIM_NUC_SIZE=1000000

	SIM_DEPTH_NUC="10x"
	SIM_DEPTH_MT="50x"
	SIM_DEPTH_PT="500x"

	SIM_OUT="sim_ont"
	SIM_THREADS=8

	# ONT tech preset
	SIM_TECH="r10.4" # r9, r9.4, r9.4.1, r10, r10.4, kit14, sup, sup_r10, dorado, duplex

	# Badread knobs (will be overridden by set_ont_tech)
	SIM_LEN_MEAN=15000

	local SIM_BIN="${_POLAPLIB_DIR}/polap-bash-sim-ont.sh"
	ARGS=()

	# pass-through arguments to the simulator, also capture --out if present
	while [[ $# -gt 0 ]]; do
		case "$1" in
		--out)
			SIM_OUT="$2"
			ARGS+=("$1" "$2")
			shift 2
			;;
		*)
			ARGS+=("$1")
			shift
			;;
		esac
	done

	mkdir -p "$SIM_OUT"
	# pushd "$SIM_OUT" >/dev/null

	# Run simulation
	echo "[sim] running: $SIM_BIN ${ARGS[*]}"
	# echo bash "$SIM_BIN" "${ARGS[@]}" >&2
	bash "$SIM_BIN" "${ARGS[@]}"

	# popd >/dev/null
}
