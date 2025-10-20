#!/usr/bin/env bash
set -euo pipefail

if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
	echo "[ERROR] This script must be executed, not sourced: use 'bash $BASH_SOURCE'" >&2
	return 1 2>/dev/null || exit 1
fi

: "${_POLAP_DEBUG:=0}"
export _POLAP_DEBUG
: "${_POLAP_RELEASE:=0}"
export _POLAP_RELEASE

_local_host="thorne"
_media_dir="/media/h2/sra"
_media1_dir="/media/h1/sra"
_media2_dir="/media/h2/sra"

# Local site defaults (override via env if needed)
: "${_LOCAL_HOST:=thorne}"
: "${_MEDIA_DIR:=/media/h2/sra}"

# Determine base
_polap_script_bin_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" || {
	echo "Couldn't determine the script's running directory, bailing out" >&2
	exit 2
}
_POLAPLIB_DIR="${_polap_script_bin_dir}/polaplib"

if [[ -r "${_POLAPLIB_DIR}/polap-lib-profiles.sh" ]]; then
	source "${_POLAPLIB_DIR}/polap-lib-profiles.sh"
	_polap_ensure_profiles_dir
fi

# ── NEW: trace-load flag
for arg in "$@"; do
	case "$arg" in
	--trace-load)
		export POLAP_AUTOLOAD_TRACE=1
		shift
		;;
	esac
done

source "${_POLAPLIB_DIR}/polap-lib-version.sh"
source "${_POLAPLIB_DIR}/bolap-parsing.sh"
_polap_output_dest="/dev/null"

# ── NEW: Use the bolap autoloader
source "${_POLAPLIB_DIR}/polap-bash-bolap-autoload-min.sh"
: "${_bolap_type:=read}"
bolap_autoload_min "${_POLAPLIB_DIR}" "${_bolap_type}"

source "${_POLAPLIB_DIR}/polap-lib-command.sh"
source "${_POLAPLIB_DIR}/polap-lib-bolap-compat.sh"
source "${_POLAPLIB_DIR}/polap-lib-seqkit.sh"
source "${_POLAPLIB_DIR}/bolap-lib-dataset.sh"
source "${_POLAPLIB_DIR}/polap-lib-dialog.sh"

# MAIN
if [ $# -eq 0 ]; then
	print_help
	exit 0
fi

_run_bolap
# Dispatch
# if declare -f "_run_bolap_${_brg_menu[0]}" >/dev/null 2>&1; then
# 	_run_bolap_${_brg_menu[0]}
# else
# 	_run_bolap "${_brg_menu[0]}"
# fi
