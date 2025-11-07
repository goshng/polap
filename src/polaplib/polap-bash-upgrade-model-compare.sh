# FILE: polap-bash-upgrade-model-compare.sh
#!/usr/bin/env bash
set -euo pipefail
VERSION="0.1.0"

usage() {
	cat <<EOF
polap-bash-upgrade-model-compare v${VERSION}
Compare Dollo-like vs Mk models for MTPT cluster presence using AICc.
Also estimates expected #gains/#losses via stochastic mapping.

USAGE:
  $(basename "$0") [-b <analysis_base>] [-t <treefile>] [-p <presence.tsv>]
                   [--nsim 200] [--qratio-thr 0.1]

Defaults:
  -b man/analysis
  -t <b>/cp_tree/concat/iqtree.treefile
  -p <b>/gain_loss/presence.tsv
Outputs:
  <b>/gain_loss/model_compare.tsv
EOF
}

ANALYSIS_BASE="man/analysis"
TREE=""
PRES=""
NSIM=200
QRATIO_THR=0.1

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
	-p)
		PRES="$2"
		shift 2
		;;
	--nsim)
		NSIM="$2"
		shift 2
		;;
	--qratio-thr)
		QRATIO_THR="$2"
		shift 2
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

[[ -z "${TREE}" ]] && TREE="${ANALYSIS_BASE}/cp_tree/concat/iqtree.treefile"
[[ -z "${PRES}" ]] && PRES="${ANALYSIS_BASE}/gain_loss/presence.tsv"

mkdir -p "${ANALYSIS_BASE}/gain_loss" "${ANALYSIS_BASE}/logs"
LOG="${ANALYSIS_BASE}/logs/upgrade_model_compare.log"

need() { command -v "$1" >/dev/null 2>&1 || {
	echo "Missing $1"
	exit 1
}; }
need Rscript

echo "polap-bash-upgrade-model-compare v${VERSION}" | tee "${LOG}"
echo "Tree: ${TREE}" | tee -a "${LOG}"
echo "Presence: ${PRES}" | tee -a "${LOG}"

Rscript scripts/mk_model_compare.R "${TREE}" "${PRES}" "${ANALYSIS_BASE}/gain_loss/model_compare.tsv" "${NSIM}" "${QRATIO_THR}" | tee -a "${LOG}"

echo "Output -> ${ANALYSIS_BASE}/gain_loss/model_compare.tsv" | tee -a "${LOG}"
