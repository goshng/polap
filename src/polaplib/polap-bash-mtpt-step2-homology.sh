# FILE: polap-bash-mtpt-step2-homology.sh
#!/usr/bin/env bash
# Version: v0.2.0
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
polap-bash-mtpt-step2-homology v0.2.0
Pool MTPT tracts across species, cluster with MMseqs2 linclust, build cluster_map.tsv,
and align each cluster. Supports --resume (skip if up-to-date and params unchanged).

USAGE:
  \$(basename "\$0") [-b man/analysis] [-t 8] [--min-id 0.85] [--min-cov 0.75] [--resume]
EOF
}

ANALYSIS_BASE="man/analysis"
THREADS=8
MIN_ID="0.85"
MIN_COV="0.75"
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
	--min-id)
		MIN_ID="$2"
		shift 2
		;;
	--min-cov)
		MIN_COV="$2"
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
need mmseqs
need python3
need mafft

uptodate() { # uptodate <out1> [out2 ...] -- <in1> [in2 ...]
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

BASE="${ANALYSIS_BASE}"
CLDIR="${BASE}/clusters"
LOG="${BASE}/logs/step2_homology.log"
mkdir -p "${CLDIR}" "${BASE}/logs"

# Inputs: all tract FASTAs produced in Step1
mapfile -t TRACTS < <(find "${BASE}/mtpt_calls" -type f -name "*.fa" | sort)
((${#TRACTS[@]})) || {
	echo "No MTPT FASTAs found under ${BASE}/mtpt_calls/" | tee "${LOG}"
	exit 3
}

ALLFA="${CLDIR}/all_mtpts.fa"
TSV="${CLDIR}/linclust.tsv"
CMAP="${CLDIR}/cluster_map.tsv"
PARAMS="${CLDIR}/.params"
ALIGN_OK="${CLDIR}/_alignments.done.ok"

# Param fingerprint
PARAM_STR="v=0.2.0|min_id=${MIN_ID}|min_cov=${MIN_COV}|threads=${THREADS}"

params_match() { [[ -s "$1" ]] && [[ "$(cat "$1")" == "$2" ]]; }

if ((RESUME == 1)) && uptodate "${CMAP}" "${ALIGN_OK}" -- "${TRACTS[@]}" && params_match "${PARAMS}" "${PARAM_STR}"; then
	echo "[step2] up-to-date (params match) → skip" | tee "${LOG}"
	exit 0
fi

echo "[1/5] Pool MTPT FASTAs → ${ALLFA}" | tee "${LOG}"
cat "${TRACTS[@]}" >"${ALLFA}"

# [2/5] MMseqs2 linclust (id=${MIN_ID}, cov=${MIN_COV})  — refresh run
echo "[2/5] MMseqs2 linclust (id=${MIN_ID}, cov=${MIN_COV})" | tee -a "${LOG}"

# Remove old outputs if present (safe: ignores if absent)
[[ -f "${CLDIR}/clu.dbtype" ]] && mmseqs rmdb "${CLDIR}/clu" 2>>"${LOG}" || true
[[ -f "${CLDIR}/db.dbtype" ]] && mmseqs rmdb "${CLDIR}/db" 2>>"${LOG}" || true

# Fresh temporary workspace
rm -rf "${CLDIR}/tmp" 2>>"${LOG}" || true
mkdir -p "${CLDIR}/tmp"

# Rebuild input DB and re-cluster
mmseqs createdb "${ALLFA}" "${CLDIR}/db" 2>>"${LOG}"
mmseqs linclust "${CLDIR}/db" "${CLDIR}/clu" "${CLDIR}/tmp" \
	--min-seq-id "${MIN_ID}" -c "${MIN_COV}" --cov-mode 0 --threads "${THREADS}" 2>>"${LOG}"

# Export cluster pairs (rep \t member)
mmseqs createtsv "${CLDIR}/db" "${CLDIR}/db" "${CLDIR}/clu" "${TSV}" 2>>"${LOG}"

echo "[3/5] Build cluster_map.tsv" | tee -a "${LOG}"
python3 "${_POLAPLIB_DIR}/scripts/mmseqs2_tsv_to_cluster_map.py" "${TSV}" "${CMAP}"

echo "[4/5] Per-cluster alignments (MAFFT L-INS-i)" | tee -a "${LOG}"
python3 "${_POLAPLIB_DIR}/scripts/align_clusters_mafft.py" "${ALLFA}" "${CMAP}" "${CLDIR}" "${THREADS}"
date >"${ALIGN_OK}"

echo "${PARAM_STR}" >"${PARAMS}"
echo "[5/5] Done: ${CMAP}" | tee -a "${LOG}"
