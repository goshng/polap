# FILE: polap-bash-step3-pt-phylogeny.sh
#!/usr/bin/env bash
# Version: v0.2.0
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Build a plastid (PT) species tree via:
#   MODE=coding (default): Oatk hmmannot → per-gene MAFFT (L-INS-i) → back-translate → mask → IQ-TREE2
#   MODE=mauve: progressiveMauve core LCBs → concat → IQ-TREE2
#
# Inputs (required):
#   <ANALYSIS_BASE>/assemblies/*/pt.fa
#
# Outputs:
#   <ANALYSIS_BASE>/pt_tree/{genes,aln,concat}/...
#   <ANALYSIS_BASE>/pt_tree/concat/iqtree.treefile
#
# Usage:
#   polap-bash-step3-pt-phylogeny.sh [-b <analysis_base>] [-t <threads>] -H <plastid_fam.fam> [--mode coding|mauve]
#                                    [--min-sites 100000] [--gap-frac 0.7]

set -Eeuo pipefail
IFS=$'\n\t'

# Resolve library dir for helpers (robust to sourcing/symlinks)
: "${_POLAPLIB_DIR:=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

# Optional call-stack logger
if [[ -r "${_POLAPLIB_DIR}/polap-lib-logcallstack.sh" ]]; then
	# shellcheck disable=SC1091
	source "${_POLAPLIB_DIR}/polap-lib-logcallstack.sh"
	polap_trap_enable || true
fi

VERSION="0.2.0"

usage() {
	cat <<EOF
polap-bash-step3-pt-phylogeny v${VERSION}
Build a PT species tree via:
  coding (default): Oatk hmmannot → MAFFT (L-INS-i) → back-translate → mask → IQ-TREE2
  mauve: progressiveMauve core LCBs → concat → IQ-TREE2

USAGE:
  \$(basename "\$0") [-b <analysis_base>] [-t <threads>] -H <plastid_fam.fam> [--mode coding|mauve]
                     [--min-sites 100000] [--gap-frac 0.7]

Inputs:
  <base>/assemblies/*/pt.fa

Outputs:
  <base>/pt_tree/concat/iqtree.treefile (+ support files)
EOF
}

# Defaults
ANALYSIS_BASE="man/analysis"
THREADS=8
MODE="coding"
HMM_LIB=""
MIN_SITES=100000
GAP_FRAC=0.7

# Parse CLI
while [[ $# -gt 0 ]]; do
	case "$1" in
	-b)
		ANALYSIS_BASE="$2"
		shift 2
		;;
	-t)
		THREADS="$2"
		shift 2
		;;
	-H)
		HMM_LIB="$2"
		shift 2
		;;
	--mode)
		MODE="$2"
		shift 2
		;;
	--min-sites)
		MIN_SITES="$2"
		shift 2
		;;
	--gap-frac)
		GAP_FRAC="$2"
		shift 2
		;;
	-h | --help)
		usage
		exit 0
		;;
	*)
		echo "Unknown option: $1" >&2
		usage
		exit 1
		;;
	esac
done

# Require HMM library only for coding mode
if [[ -z "${HMM_LIB}" && "${MODE}" == "coding" ]]; then
	usage
	exit 1
fi

need() { command -v "$1" >/dev/null 2>&1 || {
	echo "Missing $1" >&2
	exit 127
}; }

# Tool checks
if [[ "${MODE}" == "coding" ]]; then
	# Either hmmannot or hmm_annotation; and nhmmscan MUST exist (crash if missing)
	if command -v hmmannot >/dev/null 2>&1; then
		OATK_HMMANNOT="hmmannot"
	elif command -v hmm_annotation >/dev/null 2>&1; then
		OATK_HMMANNOT="hmm_annotation"
	else
		echo "Missing hmmannot/hmm_annotation" >&2
		exit 127
	fi
	need nhmmscan
	need mafft
	need python3
	need iqtree2
else
	need progressiveMauve
	need python3
	need iqtree2
fi

BASE="${ANALYSIS_BASE}"
PTDIR="${BASE}/pt_tree"
mkdir -p "${PTDIR}/genes" "${PTDIR}/aln" "${PTDIR}/concat" "${BASE}/logs"
LOG="${BASE}/logs/step3_pt.log"

# Collect pt.fa paths
declare -a PTS=()
while IFS= read -r -d '' f; do PTS+=("$f"); done < <(find "${BASE}/assemblies" -type f -name "pt.fa" -print0)
((${#PTS[@]})) || {
	echo "No pt.fa found under ${BASE}/assemblies" | tee "${LOG}"
	exit 2
}

SCDIR="${_POLAPLIB_DIR}/scripts"

if [[ "${MODE}" == "coding" ]]; then
	echo "[coding] annotate PT genes with Oatk hmmannot, then per-gene align/back-translate ..." | tee "${LOG}"

	ANN_PY="${SCDIR}/annotate_pt_genes_hmmannot.py"
	ALN_BT_PY="${SCDIR}/align_and_backtranslate.py"
	MASK_PY="${SCDIR}/mask_alignment_by_gap.py"
	PART_PY="${SCDIR}/partition_builder.py"
	for req in "${ANN_PY}" "${ALN_BT_PY}" "${MASK_PY}" "${PART_PY}"; do
		[[ -s "${req}" ]] || {
			echo "Missing ${req}" | tee -a "${LOG}"
			exit 3
		}
	done

	# Per-species: run hmmannot first (produces annotation table), then parse/extract with Python
	for pt in "${PTS[@]}"; do
		sp=$(basename "$(dirname "$pt")")
		out="${PTDIR}/genes/${sp}"
		mkdir -p "${out}"

		ann_out="${out}/pt_vs_hmms.annot.txt"
		echo "  [${sp}] hmmannot -> ${ann_out}" | tee -a "${LOG}"
		"${OATK_HMMANNOT}" -t "${THREADS}" --nhmmscan "$(command -v nhmmscan)" \
			-o "${ann_out}" "${HMM_LIB}" "${pt}" 2>>"${LOG}"

		echo "  [${sp}] extract CDS/protein per gene ..." | tee -a "${LOG}"
		python3 "${ANN_PY}" \
			--pt-fasta "${pt}" \
			--hmmannot "${ann_out}" \
			--outdir "${out}" \
			--species "${sp}"
	done

	echo "[coding] pool proteins by gene, align (MAFFT L-INS-i), back-translate, mask >${GAP_FRAC} gaps ..." | tee -a "${LOG}"
	python3 "${ALN_BT_PY}" "${PTDIR}/genes" "${PTDIR}/aln" "${THREADS}"
	python3 "${MASK_PY}" "${PTDIR}/aln" "${GAP_FRAC}"

	echo "[coding] concatenate, build partitions, IQ-TREE 2 ..." | tee -a "${LOG}"
	python3 "${PART_PY}" "${PTDIR}/aln" "${PTDIR}/concat/concat.fna" "${PTDIR}/concat/partitions.nex"
	# pushd "${PTDIR}/concat" >/dev/null
	iqtree2 -s "${PTDIR}/concat/concat.fna" -spp "${PTDIR}/concat"/partitions.nex -m MFP+MERGE -B 1000 --alrt 1000 -T AUTO -redo | tee -a "${LOG}"
	# popd >/dev/null

elif [[ "${MODE}" == "mauve" ]]; then
	echo "[mauve] progressiveMauve on PT set ..." | tee "${LOG}"
	XMFA="${PTDIR}/aln/progressiveMauve.xmfa"
	progressiveMauve --output="${XMFA}" "${PTS[@]}"

	XMFA_PY="${SCDIR}/xmfa_concat_core_lcbs.py"
	[[ -s "${XMFA_PY}" ]] || {
		echo "Missing ${XMFA_PY}" | tee -a "${LOG}"
		exit 5
	}

	python3 "${XMFA_PY}" "${XMFA}" "${MIN_SITES}" "${GAP_FRAC}" \
		"${PTDIR}/concat/concat.fna" "${PTDIR}/concat/partitions.nex"
	# pushd "${PTDIR}/concat" >/dev/null
	iqtree2 -s "${PTDIR}/concat/concat.fna" -spp "${PTDIR}/concat"/partitions.nex -m MFP+MERGE -B 1000 --alrt 1000 -T AUTO -redo | tee -a "${LOG}"
	# iqtree2 -s concat.fna -spp partitions.nex -m MFP+MERGE -B 1000 --alrt 1000 -T AUTO | tee -a "${LOG}"
	# popd >/dev/null

else
	echo "Unknown --mode ${MODE}" | tee -a "${LOG}"
	exit 1
fi

echo "Tree at ${PTDIR}/concat/iqtree.treefile"
