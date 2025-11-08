#!/usr/bin/env bash
# FILE: polap-bash-mtpt-step3-pt-phylogeny.sh
# Version: v0.2.1
# SPDX-License-Identifier: GPL-3.0-or-later
set -Eeuo pipefail
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
polap-bash-step3-pt-phylogeny v0.2.1
Build a plastid (PT) tree via coding genes (default) or progressiveMauve core LCBs.
Supports --resume (skip if outputs newer than inputs and params unchanged).

USAGE:
  \$(basename "\$0") [-b <base>] [-t <threads>] -H <plastid.fam> [--mode coding|mauve]
                     [--min-sites 100000] [--gap-frac 0.7] [--resume]
EOF
}

# Defaults
ANALYSIS_BASE="man/analysis"
THREADS=8
MODE="coding"
HMM_LIB=""
MIN_SITES=100000
GAP_FRAC=0.7
RESUME=0

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
	--resume)
		RESUME=1
		shift
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

[[ -z "${HMM_LIB}" && "${MODE}" == "coding" ]] && {
	usage
	exit 1
}

need() {
	command -v "$1" >/dev/null 2>&1 || {
		echo "Missing $1" >&2
		exit 127
	}
}

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

if [[ "${MODE}" == "coding" ]]; then
	need hmmannot || need hmm_annotation
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
CONCAT="${PTDIR}/concat/concat.fna"
PART="${PTDIR}/concat/partitions.nex"
TREE="${PTDIR}/concat/iqtree.treefile"
PARAMS="${PTDIR}/concat/.params"

# Collect pt.fa inputs
mapfile -t PTS < <(find "${BASE}/assemblies" -type f -name "pt.fa" -print0 | xargs -0 -I{} echo {} | sort)
((${#PTS[@]})) || {
	echo "No pt.fa found under ${BASE}/assemblies" | tee "${LOG}"
	exit 2
}

PARAM_STR="v=0.2.1|mode=${MODE}|gap=${GAP_FRAC}|minsites=${MIN_SITES}|threads=${THREADS}|hmm=$(basename "${HMM_LIB:-NA}")"

params_match() { [[ -s "$1" ]] && [[ "$(cat "$1")" == "$2" ]]; }

# Determine inputs set (affects resume)
IN_LIST=("${PTS[@]}")
[[ "${MODE}" == "coding" ]] && IN_LIST+=("${HMM_LIB}")

if ((RESUME == 1)) && uptodate "${TREE}" "${CONCAT}" "${PART}" -- "${IN_LIST[@]}" && params_match "${PARAMS}" "${PARAM_STR}"; then
	echo "[step3] up-to-date (params match) → skip" | tee "${LOG}"
	exit 0
fi

SCDIR="${_POLAPLIB_DIR}/scripts"

if [[ "${MODE}" == "coding" ]]; then
	echo "[coding] annotate PT genes (hmmannot), extract CDS/proteins → align/back-translate" | tee "${LOG}"
	ANN_PY="${SCDIR}/annotate_pt_genes_hmmannot.py"
	ALN_BT_PY="${SCDIR}/align_and_backtranslate.py"
	MASK_PY="${SCDIR}/mask_alignment_by_gap.py"
	PART_PY="${SCDIR}/partition_builder.py"
	[[ -s "${ANN_PY}" && -s "${ALN_BT_PY}" && -s "${MASK_PY}" && -s "${PART_PY}" ]] || {
		echo "Missing helper scripts in ${SCDIR}" | tee -a "${LOG}"
		exit 4
	}

	# --- run per PT: Oatk hmmannot -> parse & extract CDS/AA per gene ---
	# Resolve nhmmscan path (must exist; tool check should have run above)
	NHMMSCAN_BIN="$(command -v nhmmscan)"
	if [[ -z "${NHMMSCAN_BIN}" ]]; then
		echo "[ERROR] nhmmscan not found on PATH." | tee -a "${LOG}"
		exit 127
	fi
	# Prefer hmmannot; fall back to hmm_annotation
	if command -v hmmannot >/dev/null 2>&1; then
		OATK_HMMANNOT="hmmannot"
	elif command -v hmm_annotation >/dev/null 2>&1; then
		OATK_HMMANNOT="hmm_annotation"
	else
		echo "[ERROR] Neither hmmannot nor hmm_annotation is available." | tee -a "${LOG}"
		exit 127
	fi

	for pt in "${PTS[@]}"; do
		sp="$(basename "$(dirname "$pt")")"
		out="${PTDIR}/genes/${sp}"
		mkdir -p "${out}"

		ann_tbl="${out}/pt_vs_hmms.annotations.txt"

		# Reuse existing annotation table if up-to-date; else (re)run hmmannot
		if [[ -s "${ann_tbl}" && "${ann_tbl}" -nt "${pt}" && "${ann_tbl}" -nt "${HMM_LIB}" ]]; then
			echo "[coding:${sp}] reuse ${ann_tbl}" | tee -a "${LOG}"
		else
			echo "[coding:${sp}] Oatk ${OATK_HMMANNOT} → ${ann_tbl}" | tee -a "${LOG}"
			# Oatk accepts FASTA as query; pass nhmmscan explicitly
			"${OATK_HMMANNOT}" -t "${THREADS}" \
				--nhmmscan "${NHMMSCAN_BIN}" \
				-o "${ann_tbl}" \
				"${HMM_LIB}" "${pt}" 2>>"${LOG}"
		fi

		# Parse table and emit per-gene {gene}.prot.faa / {gene}.cds.fna for this species
		echo "[coding:${sp}] parse & extract CDS/proteins → ${out}" | tee -a "${LOG}"
		python3 "${ANN_PY}" \
			--pt-fasta "${pt}" \
			--hmmannot "${ann_tbl}" \
			--outdir "${out}" \
			--species "${sp}" 2>>"${LOG}"
	done

	# Pool proteins per gene across species → align (MAFFT L-INS-i) → back-translate → mask → partitions
	echo "[coding] pool & align; back-translate; mask; partition" | tee -a "${LOG}"
	python3 "${ALN_BT_PY}" "${PTDIR}/genes" "${PTDIR}/aln" "${THREADS}" 2>>"${LOG}"
	python3 "${MASK_PY}" "${PTDIR}/aln" "${GAP_FRAC}" 2>>"${LOG}"
	python3 "${PART_PY}" "${PTDIR}/aln" "${CONCAT}" "${PART}" 2>>"${LOG}"

	# pushd "${PTDIR}/concat" >/dev/null
	# iqtree2 -redo -s "$(basename "${CONCAT}")" -spp "$(basename "${PART}")" -m MFP+MERGE -B 1000 --alrt 1000 -T AUTO | tee -a "${LOG}"
	# popd >/dev/null

elif [[ "${MODE}" == "mauve" ]]; then
	echo "[mauve] progressiveMauve on pt.fa set → core LCB concat" | tee "${LOG}"
	XMFA="${PTDIR}/aln/progressiveMauve.xmfa"
	progressiveMauve --output="${XMFA}" "${PTS[@]}"
	XMFA_PY="${SCDIR}/xmfa_concat_core_lcbs.py"
	[[ -s "${XMFA_PY}" ]] || {
		echo "Missing ${XMFA_PY}" | tee -a "${LOG}"
		exit 5
	}
	python3 "${XMFA_PY}" "${XMFA}" "${MIN_SITES}" "${GAP_FRAC}" "${CONCAT}" "${PART}"

	# pushd "${PTDIR}/concat" >/dev/null
	# iqtree2 -s "$(basename "${CONCAT}")" -spp "$(basename "${PART}")" -m MFP+MERGE -B 1000 --alrt 1000 -T AUTO | tee -a "${LOG}"
	# popd >/dev/null

else
	echo "Unknown --mode ${MODE}" | tee -a "${LOG}"
	exit 1
fi

# Write all IQ-TREE artifacts under ${PTDIR}/concat with prefix "iqtree"
IQ_PREFIX="${PTDIR}/concat/iqtree"

# Ensure output dir exists (it is, but keep this robust)
mkdir -p "${PTDIR}/concat"

# Run IQ-TREE without changing directories.
# Use full paths for -s and -spp, and -pre to place outputs as ${PTDIR}/concat/iqtree.*
iqtree2 -redo \
	-s "${CONCAT}" \
	-spp "${PART}" \
	-pre "${IQ_PREFIX}" \
	-m MFP+MERGE \
	-B 1000 --alrt 1000 -T AUTO \
	2>&1 | tee -a "${LOG}"

echo "${PARAM_STR}" >"${PARAMS}"
echo "Tree at ${TREE}"
