#!/usr/bin/env bash
# FILE: polap-bash-step1-mtpt-detector.sh
# Detect PT->MT tracts (MTPTs), extract, annotate plastid genes with Oatk's hmmannot/hmm_annotation,
# and optionally score ONT boundary support + coverage overlay.
# Version: v0.1.2
# SPDX-License-Identifier: GPL-3.0-or-later
set -Eeuo pipefail
IFS=$'\n\t'
: "${_POLAPLIB_DIR:=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${_POLAPLIB_DIR}/polap-lib-logcallstack.sh"
polap_trap_enable

# --- Resolve _POLAPLIB_DIR to the base directory of this script ---
_POLAPLIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
	cat <<EOF
polap-bash-step1-mtpt-detector.sh v0.1.2
Detect PT->MT tracts (MTPTs), extract, annotate plastid genes with Oatk's hmmannot/hmm_annotation,
and optionally score ONT boundary support + coverage overlay.

USAGE:
  $(basename "$0") -s <species> -H <plastid_fam.fam> [-b <analysis_base>] [-t <threads>]
                   [--min-len 150] [--min-pid 0.85] [--merge-dist 0]
                   [--ont-bam allONT.bam --mt-assigned-bam mtAssigned.bam
                    --anchor 200 --mapq 20]

Inputs (expected):
  man/analysis/assemblies/<species>/mt.fa
  man/analysis/assemblies/<species>/pt.fa

Outputs:
  man/analysis/mtpt_calls/<species>/{mtpt.bed,mtpt.tsv,fasta/*.fa}
  optional coverage_overlay.pdf

Notes:
  - All helper scripts are invoked from: \${_POLAPLIB_DIR}/scripts
  - Plastid annotation uses Oatk's 'hmmannot' or 'hmm_annotation' and requires 'nhmmscan' in PATH.
    Example DB: OatkDB/v20230921/embryophyta_pltd.fam
EOF
}

ANALYSIS_BASE="man/analysis"
THREADS=4
MIN_LEN=150
MIN_PID=0.85
MERGE_DIST=0
ONT_BAM=""
MTASSIGN_BAM=""
ANCHOR=200
MAPQ=20
SPECIES=""
HMM_LIB="" # e.g., OatkDB/v20230921/embryophyta_pltd.fam

while [[ $# -gt 0 ]]; do
	case "$1" in
	-s)
		SPECIES="$2"
		shift 2
		;;
	-H)
		HMM_LIB="$2"
		shift 2
		;;
	-b)
		ANALYSIS_BASE="$2"
		shift 2
		;;
	-t)
		THREADS="$2"
		shift 2
		;;
	--min-len)
		MIN_LEN="$2"
		shift 2
		;;
	--min-pid)
		MIN_PID="$2"
		shift 2
		;;
	--merge-dist)
		MERGE_DIST="$2"
		shift 2
		;;
	--ont-bam)
		ONT_BAM="$2"
		shift 2
		;;
	--mt-assigned-bam)
		MTASSIGN_BAM="$2"
		shift 2
		;;
	--anchor)
		ANCHOR="$2"
		shift 2
		;;
	--mapq)
		MAPQ="$2"
		shift 2
		;;
	-h | --help)
		usage
		exit 0
		;;
	*)
		echo "Unknown arg $1"
		usage
		exit 1
		;;
	esac
done

[[ -z "${SPECIES}" || -z "${HMM_LIB}" ]] && {
	usage
	exit 1
}

# Paths
BASE="${ANALYSIS_BASE}"
ASM_DIR="${BASE}/assemblies/${SPECIES}"
MT="${ASM_DIR}/mt.fa"
PT="${ASM_DIR}/pt.fa"
OUT_DIR="${BASE}/mtpt_calls/${SPECIES}"
WORK="${BASE}/work/${SPECIES}"
LOG_DIR="${BASE}/logs"
mkdir -p "${OUT_DIR}/fasta" "${WORK}" "${LOG_DIR}"
LOG="${LOG_DIR}/step1_${SPECIES}.log"

die() {
	echo "[ERR] $*" | tee -a "${LOG}"
	exit 1
}
need() { command -v "$1" >/dev/null 2>&1 || die "Missing required tool '$1'."; }

echo "polap-bash-step1-mtpt-detector.sh v0.1.2" | tee "${LOG}"
echo "Species=${SPECIES}; HMM_LIB=$(realpath -m "${HMM_LIB}")" | tee -a "${LOG}"
echo "_POLAPLIB_DIR=${_POLAPLIB_DIR}" | tee -a "${LOG}"

# Required tools
need minimap2
need bedtools
need seqkit
need python
need samtools

# Oatk annotators (one of them must exist)
if command -v hmmannot >/dev/null 2>&1; then
	OATK_HMMANNOT="hmmannot"
elif command -v hmm_annotation >/dev/null 2>&1; then
	OATK_HMMANNOT="hmm_annotation"
else
	die "neither 'hmmannot' nor 'hmm_annotation' found in PATH."
fi

# nhmmscan is mandatory (crash if absent) — requirement #3
need nhmmscan

[[ -f "${MT}" ]] || die "Missing ${MT}"
[[ -f "${PT}" ]] || die "Missing ${PT}"

PAF="${WORK}/pt_vs_mt.paf"
BED_HITS="${WORK}/ptmt_hits.bed"
TRACTS_BED="${OUT_DIR}/mtpt.bed"
TRACTS_TSV="${OUT_DIR}/mtpt.tsv"
ALL_MTPT_FASTA="${WORK}/all_mtpts.fa"
HMMANN_OUT="${WORK}/mtpts_vs_plastid.annot.txt"
BOUNDARY_TSV="${WORK}/boundary_support.tsv"

echo "[1/7] minimap2 PT->MT ..." | tee -a "${LOG}"
minimap2 -x asm10 -t "${THREADS}" "${MT}" "${PT}" >"${PAF}" 2>>"${LOG}"

echo "[2/7] filter PAF -> BED (len>=${MIN_LEN}, pid>=${MIN_PID}) ..." | tee -a "${LOG}"
python "${_POLAPLIB_DIR}/scripts/paf_to_bed_with_pid.py" "${PAF}" "${MIN_LEN}" "${MIN_PID}" >"${BED_HITS}"

echo "[3/7] merge hits on MT axis into tracts (merge_dist=${MERGE_DIST}) ..." | tee -a "${LOG}"
python "${_POLAPLIB_DIR}/scripts/merge_mt_hits_to_tracts.py" \
	"${BED_HITS}" "${MERGE_DIST}" "${SPECIES}" >"${TRACTS_TSV}"

# also write BED from TSV
awk 'BEGIN{OFS="\t"} NR>1{print $2,$3,$4,$1, int($7*1000), $6}' "${TRACTS_TSV}" >"${TRACTS_BED}"

echo "[4/7] extract tract FASTAs (from MT) ..." | tee -a "${LOG}"
python "${_POLAPLIB_DIR}/scripts/extract_fasta_regions.py" "${MT}" "${TRACTS_BED}" "${OUT_DIR}/fasta" "${ALL_MTPT_FASTA}"

echo "[5/7] Oatk hmmannot over pooled MTPT FASTA (requires nhmmscan) ..." | tee -a "${LOG}"
"${OATK_HMMANNOT}" -t "${THREADS}" --nhmmscan "$(command -v nhmmscan)" \
	-o "${HMMANN_OUT}" "${HMM_LIB}" "${ALL_MTPT_FASTA}" 2>>"${LOG}"

echo "[6/7] name tracts by ordered plastid gene content (from hmmannot) ..." | tee -a "${LOG}"
python "${_POLAPLIB_DIR}/scripts/name_mtpts_from_hmmannot.py" "${TRACTS_TSV}" "${HMMANN_OUT}" >"${TRACTS_TSV}.tmp"
mv "${TRACTS_TSV}.tmp" "${TRACTS_TSV}"

# Optional: boundary support + coverage overlay if BAMs provided (mosdepth optional)
if [[ -n "${ONT_BAM}" ]]; then
	echo "[7/7] boundary support with ONT (anchor=${ANCHOR}, MAPQ>=${MAPQ}) ..." | tee -a "${LOG}"
	python "${_POLAPLIB_DIR}/scripts/boundary_support_counter.py" "${MT}" "${TRACTS_BED}" "${ONT_BAM}" "${ANCHOR}" "${MAPQ}" >"${BOUNDARY_TSV}" || true
	# merge into TSV if counts exist
	if [[ -s "${BOUNDARY_TSV}" ]]; then
		python "${_POLAPLIB_DIR}/scripts/merge_boundary_into_mtpt_tsv.py" "${TRACTS_TSV}" "${BOUNDARY_TSV}" >"${TRACTS_TSV}.tmp2"
		mv "${TRACTS_TSV}.tmp2" "${TRACTS_TSV}"
	fi
	if [[ -n "${MTASSIGN_BAM}" ]]; then
		if command -v mosdepth >/dev/null 2>&1; then
			echo "[opt] coverage overlay (all ONT vs MT-assigned) ..." | tee -a "${LOG}"
			mosdepth -t "${THREADS}" -n --by 200 "${WORK}/allONT" "${ONT_BAM}" 2>>"${LOG}"
			mosdepth -t "${THREADS}" -n --by 200 "${WORK}/mtAssigned" "${MTASSIGN_BAM}" 2>>"${LOG}"
			Rscript "${_POLAPLIB_DIR}/scripts/coverage_overlay.R" \
				"${WORK}/allONT.regions.bed.gz" "${WORK}/mtAssigned.regions.bed.gz" \
				"${TRACTS_BED}" "${OUT_DIR}/coverage_overlay.pdf"
		else
			echo "[opt] mosdepth not found; skipping coverage overlay." | tee -a "${LOG}"
		fi
	fi
fi

echo "Done. Outputs at ${OUT_DIR}/" | tee -a "${LOG}"
