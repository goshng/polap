# FILE: polap-bash-upgrade-bootstrap-asr.sh
#!/usr/bin/env bash
set -euo pipefail
VERSION="0.1.0"

usage() {
	cat <<EOF
polap-bash-upgrade-bootstrap-asr v${VERSION}
Propagate topology uncertainty using IQ-TREE bootstrap trees to quantify
uncertainty on expected #gains/#losses per MTPT cluster.

USAGE:
  $(basename "$0") [-b <analysis_base>] [-p <presence.tsv>] [-u <ufboot.trees>]
                   [--B 200] [--nsim 10] [--model auto|ER|ARD]

Defaults:
  -b man/analysis
  -p <b>/gain_loss/presence.tsv
  -u <b>/cp_tree/concat/iqtree.ufboot
  --model auto uses best model per cluster from model_compare.tsv (falls back to ER)
Outputs:
  <b>/gain_loss/boot_asr_summary.tsv
EOF
}

ANALYSIS_BASE="man/analysis"
PRES=""
UFBOOT=""
B=200
NSIM=10
MODEL="auto"

while [[ $# -gt 0 ]]; do
	case "$1" in
	-b)
		ANALYSIS_BASE="$2"
		shift 2
		;;
	-p)
		PRES="$2"
		shift 2
		;;
	-u)
		UFBOOT="$2"
		shift 2
		;;
	--B)
		B="$2"
		shift 2
		;;
	--nsim)
		NSIM="$2"
		shift 2
		;;
	--model)
		MODEL="$2"
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

[[ -z "${PRES}" ]] && PRES="${ANALYSIS_BASE}/gain_loss/presence.tsv"
[[ -z "${UFBOOT}" ]] && UFBOOT="${ANALYSIS_BASE}/cp_tree/concat/iqtree.ufboot"

mkdir -p "${ANALYSIS_BASE}/gain_loss" "${ANALYSIS_BASE}/logs"
LOG="${ANALYSIS_BASE}/logs/upgrade_bootstrap_asr.log"

need() { command -v "$1" >/dev/null 2>&1 || {
	echo "Missing $1"
	exit 1
}; }
need Rscript

MC="${ANALYSIS_BASE}/gain_loss/model_compare.tsv"
if [[ ! -f "${MC}" && "${MODEL}" == "auto" ]]; then
	echo "[WARN] model_compare.tsv not found; using ER for all clusters." | tee "${LOG}"
fi

Rscript scripts/mk_bootstrap_asr.R "${PRES}" "${UFBOOT}" "${ANALYSIS_BASE}/gain_loss/boot_asr_summary.tsv" "${B}" "${NSIM}" "${MODEL}" "${MC}" | tee -a "${LOG}"

echo "Output -> ${ANALYSIS_BASE}/gain_loss/boot_asr_summary.tsv" | tee -a "${LOG}"
