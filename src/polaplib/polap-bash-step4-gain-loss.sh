# FILE: polap-bash-step4-gain-loss.sh
#!/usr/bin/env bash
# Version: v0.2.0
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Map presence/absence of MTPT clusters on the PT tree using parsimony (Fitch/Dollo)
# and/or Mk likelihood.

set -euo pipefail
IFS=$'\n\t'

# Resolve Polap library base (so helpers are looked up relative to this script)
: "${_POLAPLIB_DIR:=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

usage() {
	cat <<EOF
polap-bash-step4-gain-loss v${VERSION}
Map presence/absence of MTPT clusters on PT tree using parsimony (Fitch/Dollo) and/or Mk likelihood.

USAGE:
  \$(basename "\$0") [-b <analysis_base>] [-t <treefile>] [--model mk|fitch|dollo|all]

Inputs:
  <base>/clusters/cluster_map.tsv     (from Step 2)
  <base>/pt_tree/concat/iqtree.treefile  (from Step 3, default)

Outputs:
  <base>/gain_loss/{presence.tsv, events.tsv, timeline.tsv}
EOF
}

VERSION="0.2.0"

# Defaults
ANALYSIS_BASE="man/analysis"
TREE=""
MODEL="all"
OUTPDF=""

# Parse CLI
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
	--outpdf)
		OUTPDF="$2"
		shift 2
		;;
	-h | --help)
		usage
		exit 0
		;;
	*)
		echo "Unknown: $1" >&2
		usage
		exit 1
		;;
	esac
done

BASE="${ANALYSIS_BASE}"
[[ -z "${TREE}" ]] && TREE="${BASE}/pt_tree/concat/iqtree.treefile"
CLMAP="${BASE}/clusters/cluster_map.tsv"
OUTD="${BASE}/gain_loss"
mkdir -p "${OUTD}" "${BASE}/logs"
LOG="${BASE}/logs/step4_gainloss.log"

need() { command -v "$1" >/dev/null 2>&1 || {
	echo "Missing $1" >&2
	exit 127
}; }
need python3
need Rscript

echo "Create presence/absence matrix ..." | tee "${LOG}"
python3 "${_POLAPLIB_DIR}/scripts/create_presence_matrix.py" "${CLMAP}" "${OUTD}/presence.tsv"

echo "Run ASR mapping (${MODEL}) ..." | tee -a "${LOG}"
Rscript "${_POLAPLIB_DIR}/scripts/mk_mapping.R" "${TREE}" "${OUTD}/presence.tsv" "${OUTD}" "${MODEL}" | tee -a "${LOG}"

LOG="${ANALYSIS_BASE}/logs/plot_pt_gainloss.log"
PRES="${OUTD}/presence.tsv"
TIP_SIZE=2.5
LAMBDA=1.0
CHRONOS=0
MODEL=fitch

Rscript "${_POLAPLIB_DIR}/scripts/plot_mtpt_timescaled.R" \
	--tree "${TREE}" \
	--presence "${PRES}" \
	--out "${OUTPDF}" \
	--model "${MODEL}" \
	--tip-size "${TIP_SIZE}" \
	$( ((CHRONOS == 1)) && echo --chronos) \
	--lambda "${LAMBDA}" 2>&1 | tee "${LOG}"

echo "Wrote ${OUTPDF}"

echo "Done. See ${OUTD}/"
