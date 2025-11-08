# FILE: polap-bash-step4-gain-loss.sh
#!/usr/bin/env bash
# Version: v0.2.1
# SPDX-License-Identifier: GPL-3.0-or-later
set -euo pipefail
IFS=$'\n\t'
: "${_POLAPLIB_DIR:=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

# ------------------------------------------------------------------------------
# enable failsafe
# re-exec in bash if invoked via /bin/sh
# strict + trap propagation
# ------------------------------------------------------------------------------
had_u=0
case $- in *u*) had_u=1 ;; esac
set +u
if [ -z "${BASH_VERSION+x}" ]; then exec /usr/bin/env bash "$0" "$@"; fi
[ "$had_u" -eq 1 ] && set -u
set -Eeuo pipefail
set -o errtrace
set -o functrace
source "${_POLAPLIB_DIR}/polap-lib-failsafe.sh"
polap_enable_failsafe

usage() {
	cat <<EOF
polap-bash-step4-gain-loss v${VERSION}
Map presence/absence of MTPT clusters on PT tree using parsimony (Fitch/Dollo-like)
and/or Mk(ER) likelihood. Supports --resume (skip if up-to-date and params unchanged).

USAGE:
  \$(basename "\$0") [-b <analysis_base>] [-t <treefile>]
                     [--model mk|fitch|dollo|all] [--resume]
Inputs:
  <base>/clusters/cluster_map.tsv
  <base>/pt_tree/concat/iqtree.treefile (default if -t not given)
Outputs:
  <base>/gain_loss/{presence.tsv, events.tsv, timeline.tsv}
EOF
}
VERSION="0.2.1"

ANALYSIS_BASE="man/analysis"
TREE=""
MODEL="all"
RESUME=0

while [[ $# -gt 0 ]]; do
	case "$1" in
	-b)
		ANALYSIS_BASE="$2"
		shift 2
		;;
	-t)
		TREE="$2"
		shift 2
		;;
	--model)
		MODEL="$2"
		shift 2
		;;
	--resume)
		RESUME=1
		shift
		;;
	-h | --help)
		usage
		exit 0
		;;
	*)
		echo "Unknown: $1"
		usage
		exit 1
		;;
	esac
done

need() { command -v "$1" >/dev/null 2>&1 || {
	echo "Missing $1" >&2
	exit 127
}; }
uptodate() { # uptodate <out...> -- <in...>
	local outs=() ins=() sep=0 a m min_out="" max_in=0
	for a in "$@"; do
		[[ "$a" == "--" ]] && {
			sep=1
			continue
		}
		((sep == 0)) && outs+=("$a") || ins+=("$a")
	done
	for a in "${outs[@]}"; do [[ -e "$a" ]] || return 1; done
	for a in "${outs[@]}"; do
		m=$(stat -c %Y "$a" 2>/dev/null || echo 0)
		[[ -z "$min_out" || $m -lt $min_out ]] && min_out=$m
	done
	for a in "${ins[@]}"; do
		[[ -e "$a" ]] || continue
		m=$(stat -c %Y "$a" 2>/dev/null || echo 0)
		((m > max_in)) && max_in=$m
	done
	[[ -n "$min_out" && $min_out -ge $max_in ]]
}

need python3
need Rscript

BASE="${ANALYSIS_BASE}"
[[ -z "${TREE}" ]] && TREE="${BASE}/pt_tree/concat/iqtree.treefile"
CLMAP="${BASE}/clusters/cluster_map.tsv"
OUTD="${BASE}/gain_loss"
mkdir -p "${OUTD}" "${BASE}/logs"
LOG="${BASE}/logs/step4_gainloss.log"

PRES="${OUTD}/presence.tsv"
EVT="${OUTD}/events.tsv"
TML="${OUTD}/timeline.tsv"
PARAMS="${OUTD}/.params"
PARAM_STR="v=${VERSION}|model=${MODEL}|tree=$(basename "${TREE}")"

params_match() { [[ -s "$1" ]] && [[ "$(cat "$1")" == "$2" ]]; }

# Resume check: both presence.tsv and events.tsv should be newer than inputs & params
if ((RESUME == 1)) && uptodate "${PRES}" "${EVT}" -- "${CLMAP}" "${TREE}" && params_match "${PARAMS}" "${PARAM_STR}"; then
	echo "[step4] up-to-date (params match) → skip" | tee "${LOG}"
	exit 0
fi

echo "Create presence/absence matrix ..." | tee "${LOG}"
# python3 "${_POLAPLIB_DIR}/scripts/create_presence_matrix.py" \
# 	"${CLMAP}" \
#   "${PRES}" \
# 	--tree "${TREE}" --strip-tree-suffix "-0" 2>&1 | tee -a "${LOG}"

# FILE: polap-bash-mtpt-step4-gain-loss.sh (orchestrator already updated to new name)
python3 "${_POLAPLIB_DIR}/scripts/create_presence_matrix.py" \
	"${CLMAP}" \
	"${PRES}" \
	--tree "${TREE}" \
	--strip-tree-suffix "-0" \
	--strip-presence-suffix "-0" 2>&1 | tee -a "${LOG}"

echo "Run ASR mapping (${MODEL}) ..." | tee -a "${LOG}"
Rscript "${_POLAPLIB_DIR}/scripts/mk_mapping.R" "${TREE}" "${PRES}" "${OUTD}" "${MODEL}" | tee -a "${LOG}"

echo "${PARAM_STR}" >"${PARAMS}"
echo "Done. See ${OUTD}/"
