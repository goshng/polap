# FILE: polap-bash-upgrade-ir-standardize.sh
#!/usr/bin/env bash
set -euo pipefail
VERSION="0.1.0"

usage() {
	cat <<EOF
polap-bash-upgrade-ir-standardize v${VERSION}
Detect IRs in cp.fa via self-alignment (minimap2), drop one IR (keep IRa),
and rotate genome so rbcL starts at pos 1. Emits standardized cp FASTA.

USAGE:
  $(basename "$0") [-b <analysis_base>] -H <plastid_hmms.hmm> [--min-ir 5000]

Inputs:
  <b>/assemblies/<species>/cp.fa
Outputs:
  <b>/cp_tree/standardized/<species>/cp.std.fa
  <b>/cp_tree/standardized/manifest.tsv (appended)
EOF
}

ANALYSIS_BASE="man/analysis"
HMM=""
MIN_IR=5000
while [[ $# -gt 0 ]]; do
	case "$1" in
	-b)
		ANALYSIS_BASE="$2"
		shift 2
		;;
	-H)
		HMM="$2"
		shift 2
		;;
	--min-ir)
		MIN_IR="$2"
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
[[ -z "${HMM}" ]] && {
	usage
	exit 1
}

need() { command -v "$1" >/dev/null 2>&1 || {
	echo "Missing $1"
	exit 1
}; }
need minimap2
need hmmsearch
need python

outbase="${ANALYSIS_BASE}/cp_tree/standardized"
mkdir -p "${outbase}" "${ANALYSIS_BASE}/logs"
MAN="${outbase}/manifest.tsv"
[[ ! -f "${MAN}" ]] && echo -e "species\tir_len\tkeep_ir_start\tkeep_ir_end\trbcL_pos" >"${MAN}"

while IFS= read -r -d '' cp; do
	sp=$(basename "$(dirname "$cp")")
	spdir="${outbase}/${sp}"
	mkdir -p "${spdir}"
	python scripts/cp_ir_standardize.py "${cp}" "${HMM}" "${spdir}/cp.std.fa" "${MIN_IR}" >"${spdir}/standardize.log" 2>&1
	if [[ -f "${spdir}/cp.std.fa" ]]; then
		# pull summary
		awk -F'\t' 'BEGIN{irec=""} $1=="IR_LEN"{irlen=$2} $1=="KEEP_IR_START"{kst=$2} $1=="KEEP_IR_END"{ked=$2} $1=="RBCL_POS"{rb=$2}
                END{print "'"${sp}"'", irlen, kst, ked, rb}' OFS="\t" "${spdir}/standardize.log" >>"${MAN}" || true
		echo "[OK] ${sp} standardized -> ${spdir}/cp.std.fa"
	else
		echo "[WARN] Failed to standardize ${sp}; see ${spdir}/standardize.log"
	fi
done < <(find "${ANALYSIS_BASE}/assemblies" -type f -name "cp.fa" -print0)
