# FILE: polap-bash-extras-erosion-spectrum.sh
#!/usr/bin/env bash
set -euo pipefail
VERSION="0.1.0"

usage() {
	cat <<EOF
polap-bash-extras-erosion-spectrum v${VERSION}
Quantile regression of MTPT length vs divergence (1-PID) per species or clade.

USAGE:
  $(basename "$0") [-b <analysis_base>] [-o <outdir>] [--clades <species_to_clade.tsv>]
                   [--taus 0.5,0.9] [--boots 500]

Inputs:
  <b>/mtpt_calls/<species>/mtpt.tsv
Outputs:
  <o>/erosion_slopes.tsv
  <o>/erosion_plot.pdf
EOF
}

ANALYSIS_BASE="man/analysis"
OUTDIR=""
CLADES=""
TAUS="0.5,0.9"
BOOTS=500

while [[ $# -gt 0 ]]; do
	case "$1" in
	-b)
		ANALYSIS_BASE="$2"
		shift 2
		;;
	-o)
		OUTDIR="$2"
		shift 2
		;;
	--clades)
		CLADES="$2"
		shift 2
		;;
	--taus)
		TAUS="$2"
		shift 2
		;;
	--boots)
		BOOTS="$2"
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

[[ -z "${OUTDIR}" ]] && OUTDIR="${ANALYSIS_BASE}/results"
mkdir -p "${OUTDIR}" "${ANALYSIS_BASE}/logs"
LOG="${ANALYSIS_BASE}/logs/extras_erosion.log"

need() { command -v "$1" >/dev/null 2>&1 || {
	echo "Missing $1"
	exit 1
}; }
need Rscript

Rscript scripts/quantile_slope.R "${ANALYSIS_BASE}/mtpt_calls" "${OUTDIR}" "${CLADES}" "${TAUS}" "${BOOTS}" | tee "${LOG}"

echo "Results -> ${OUTDIR}/erosion_slopes.tsv and ${OUTDIR}/erosion_plot.pdf"
