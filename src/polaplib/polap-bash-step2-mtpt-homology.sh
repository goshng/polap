#!/usr/bin/env bash
# FILE: polap-bash-step2-mtpt-homology.sh
# Version: v0.1.1
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Pool MTPT FASTAs, cluster with MMseqs2 linclust, align each cluster (MAFFT L-INS-i),
# and compute alignment stats. Optionally build HMMs per cluster.
#
# Inputs:
#   man/analysis/mtpt_calls/*/fasta/*.fa
#
# Outputs:
#   man/analysis/clusters/cluster_map.tsv
#   man/analysis/clusters/<CID>/align.fa
#   man/analysis/clusters/<CID>/align.stats.tsv
#
# Usage:
#   polap-bash-step2-mtpt-homology.sh [-b <analysis_base>] [-t <threads>] [--id 0.85] [--cov 0.7] [--hmmbuild]

set -Eeuo pipefail

# Resolve library dir (so helpers & scripts are looked up relative to here)
: "${_POLAPLIB_DIR:=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

# Optional: if you have your logging/traps library, source it (safe if absent)
if [[ -r "${_POLAPLIB_DIR}/polap-lib-logcallstack.sh" ]]; then
	# shellcheck disable=SC1091
	source "${_POLAPLIB_DIR}/polap-lib-logcallstack.sh"
	polap_trap_enable || true
fi

VERSION="0.1.1"

usage() {
	cat <<EOF
polap-bash-step2-mtpt-homology v${VERSION}
Pool MTPT FASTAs, cluster with MMseqs2, align each cluster (MAFFT L-INS-i), compute alignment stats.

USAGE:
  $(basename "$0") [-b <analysis_base>] [-t <threads>] [--id 0.85] [--cov 0.7] [--hmmbuild]

Inputs:
  man/analysis/mtpt_calls/*/fasta/*.fa

Outputs:
  man/analysis/clusters/cluster_map.tsv
  man/analysis/clusters/<CID>/align.fa, align.stats.tsv
EOF
}

# Defaults
ANALYSIS_BASE="man/analysis"
THREADS=8
ID=0.85
COV=0.7
HMMBUILD=0

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
	--id)
		ID="$2"
		shift 2
		;;
	--cov)
		COV="$2"
		shift 2
		;;
	--hmmbuild)
		HMMBUILD=1
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

# Require tools
need() { command -v "$1" >/dev/null 2>&1 || {
	echo "ERROR: Missing dependency: $1" >&2
	exit 127
}; }
need mmseqs
need mafft
need python3

# Layout
BASE="${ANALYSIS_BASE}"
CLDIR="${BASE}/clusters"
LOGDIR="${BASE}/logs"
TMP="${CLDIR}/tmp"
POOLED="${CLDIR}/all_mtpts.fa"
RAWTSV="${CLDIR}/cluster_raw.tsv"
MAPTSV="${CLDIR}/cluster_map.tsv"

rm -rf "${CLDIR}"
mkdir -p "${CLDIR}" "${LOGDIR}" "${TMP}"

LOG="${LOGDIR}/step2_clusters.log"

echo "[$(date '+%F %T')] Pooling MTPT FASTAs ..." | tee "${LOG}"
# Be strict but predictable about file order for reproducibility.
# Use -print0 to handle odd filenames safely.
find "${BASE}/mtpt_calls" -type f -name 'MTPT_*.fa' -print0 |
	sort -z |
	xargs -0 cat >"${POOLED}"

# MMseqs2 pipeline
echo "[$(date '+%F %T')] MMseqs2 linclust (id=${ID}, cov=${COV}) ..." | tee -a "${LOG}"
# mmseqs createdb "${POOLED}" "${TMP}/db"
mmseqs createdb "${POOLED}" "${TMP}/db" >/dev/null
# mmseqs linclust "${TMP}/db" "${TMP}/clu" "${TMP}/work" \
# 	--min-seq-id "${ID}" -c "${COV}" --cov-mode 1 --threads "${THREADS}"
mmseqs linclust "${TMP}/db" "${TMP}/clu" "${TMP}/work" \
	--min-seq-id "${ID}" -c "${COV}" --cov-mode 1 --threads "${THREADS}" >/dev/null
# mmseqs createtsv "${TMP}/db" "${TMP}/db" "${TMP}/clu" "${RAWTSV}"
mmseqs createtsv "${TMP}/db" "${TMP}/db" "${TMP}/clu" "${RAWTSV}" >/dev/null

# Ensure helper script exists (multi-line => real file under scripts/)
SCDIR="${_POLAPLIB_DIR}/scripts"
CLMAP_PY="${SCDIR}/cluster_map_from_mmseqs.py"
mkdir -p "${SCDIR}"
if [[ ! -s "${CLMAP_PY}" ]]; then
	cat >"${CLMAP_PY}" <<'PY'
#!/usr/bin/env python3
# cluster_map_from_mmseqs.py
# Version: v0.1.1
# Convert "rep<TAB>member" MMseqs TSV to "MTPT_id<TAB>CID" map with CID = C0001, C0002, ...
import sys, csv
if len(sys.argv) != 3:
    sys.stderr.write("Usage: cluster_map_from_mmseqs.py <raw.tsv> <out.tsv>\n")
    sys.exit(2)
raw, out = sys.argv[1], sys.argv[2]
clmap = {}
with open(raw, 'r', encoding='utf-8') as fh:
    for line in fh:
        line=line.strip()
        if not line: continue
        rep, mem = line.split("\t", 1)
        clmap.setdefault(rep, []).append(mem)
with open(out, 'w', encoding='utf-8', newline='') as o:
    w = csv.writer(o, delimiter="\t")
    w.writerow(["MTPT_id", "CID"])
    for idx, (rep, mems) in enumerate(clmap.items(), start=1):
        cid = f"C{idx:04d}"
        for m in mems:
            # MMseqs createtsv yields sequence IDs as they appeared; keep first token, strip leading >
            mtpt = m.split()[0]
            if mtpt.startswith(">"): mtpt = mtpt[1:]
            if mtpt.endswith(".fa"): mtpt = mtpt[:-3]
            w.writerow([mtpt, cid])
PY
	chmod +x "${CLMAP_PY}"
fi

echo "[$(date '+%F %T')] Build cluster map ..." | tee -a "${LOG}"
python3 "${CLMAP_PY}" "${RAWTSV}" "${MAPTSV}"

# Emit per-cluster alignments + stats (expects existing helper)
EMIT_PY="${SCDIR}/emit_cluster_alignments.py"
if [[ ! -s "${EMIT_PY}" ]]; then
	echo "ERROR: Missing ${EMIT_PY}. Please provide this helper." | tee -a "${LOG}"
	exit 2
fi

echo "[$(date '+%F %T')] Per-cluster alignment + stats (MAFFT L-INS-i) ..." | tee -a "${LOG}"
python3 "${EMIT_PY}" "${POOLED}" "${MAPTSV}" "${CLDIR}" "${THREADS}"

# Optional HMMs per cluster (placeholder loop)
if [[ "${HMMBUILD}" -eq 1 ]]; then
	need hmmbuild
	echo "[$(date '+%F %T')] Optional hmmbuild per cluster ..." | tee -a "${LOG}"
	# Example scaffold: iterate unique CIDs from map
	while IFS=$'\t' read -r _mtpt cid; do
		[[ "$cid" == "CID" ]] && continue
		# e.g., hmmbuild "${CLDIR}/${cid}/${cid}.hmm" "${CLDIR}/${cid}/align.fa"
		:
	done < <(cut -f2 "${MAPTSV}" | sort -u | awk 'NR==1{print "CID"; next} {print "X\t"$0}') # keeps header guard
fi

echo "[$(date '+%F %T')] Done. Outputs under: ${CLDIR}/" | tee -a "${LOG}"
