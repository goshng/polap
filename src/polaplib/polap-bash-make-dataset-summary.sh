#!/usr/bin/env bash
################################################################################
# polap-bash-make-dataset-summary-v0.1.0.sh
#
# Version : v0.1.0
# Author  : Sang Chul Choi  (POLAP project)
# Created : 2025-10-05
# Purpose : Summarize per-species ONT (or other long-read) datasets using
#           seqkit ≥ 2.10, which provides the AvgQual column.
# License : GPL-3.0+
#
# Changelog:
#   v0.1.0 – initial defensive implementation
################################################################################
set -euo pipefail
IFS=$'\n\t'

################################################################################
# Configuration defaults
################################################################################
readonly SCRIPT_VERSION="v0.1.0"
readonly REQUIRED_SEQKIT_VERSION="2.10"

: "${_POLAPLIB_DIR:?POLAP library directory must be set}"

SET="some" # all|some|test|each|file:<list>|<species>
ROOT="."
TIER="t5"
INUM="0"
OUT="md/dataset-summary.tsv"
PARALLEL=1
THREADS="$(nproc || echo 4)"
CSV_FOR_SETS="${_polap_data_csv:-polap-bash-report.csv}"

TMPDIR="$(mktemp -d -t polap-dsetsum.XXXXXX)"
trap 'rm -rf "${TMPDIR}"' EXIT INT TERM

################################################################################
# Defensive utilities
################################################################################
_die() {
	echo "[ERR] $*" >&2
	exit 1
}
_warn() { echo "[WARN] $*" >&2; }
_info() { echo "[INFO] $*" >&2; }

_check_cmd() {
	command -v "$1" >/dev/null 2>&1 || _die "Required command not found: $1"
}

version_ge() { # $1>= $2 ?
	[[ "$(printf '%s\n' "$2" "$1" | sort -V | head -n1)" == "$2" ]]
}

################################################################################
# Help message
################################################################################
usage() {
	cat <<EOF
polap-bash-make-dataset-summary-v0.1.0.sh  —  Dataset summary table builder

Usage:
  $(basename "$0") [options]

Options:
  --set SET        Species set: all|some|test|each|file:<list>|<species>
  --root DIR       Root directory (default: .)
  --t TIER         t-folder name (default: t5)
  --inum N         index folder (default: 0)
  -o, --out FILE   Output TSV (default: md/dataset-summary.tsv)
  -j, --threads N  Number of parallel workers (default: nproc)
  --no-parallel    Disable parallel mode
  -h, --help       Show this help and exit
  --version        Print script version

Output columns:
  species  total_bases  read_count  mean_length  N50  mean_Q  min_len  max_len

Requires:
  seqkit ≥ ${REQUIRED_SEQKIT_VERSION}

EOF
}

################################################################################
# Parse arguments
################################################################################
while (($#)); do
	case "$1" in
	--set)
		SET="${2:?}"
		shift 2
		;;
	--root)
		ROOT="${2:?}"
		shift 2
		;;
	--t)
		TIER="${2:?}"
		shift 2
		;;
	--inum)
		INUM="${2:?}"
		shift 2
		;;
	-o | --out)
		OUT="${2:?}"
		shift 2
		;;
	-j | --threads)
		THREADS="${2:?}"
		shift 2
		;;
	--no-parallel)
		PARALLEL=0
		shift
		;;
	-h | --help)
		usage
		exit 0
		;;
	--version)
		echo "${SCRIPT_VERSION}"
		exit 0
		;;
	*) _die "Unknown option: $1" ;;
	esac
done

################################################################################
# Pre-flight checks
################################################################################
_check_cmd seqkit
_check_cmd awk
_check_cmd csvtk || _warn "csvtk not found; --set all will fail if used."

SEQKIT_VER="$(seqkit version 2>/dev/null | awk '/Version/ {print $2}' | tr -d 'v')"
version_ge "${SEQKIT_VER:-0}" "${REQUIRED_SEQKIT_VERSION}" ||
	_warn "seqkit ${SEQKIT_VER:-0} < ${REQUIRED_SEQKIT_VERSION}; AvgQual may be missing."

mkdir -p "$(dirname "$OUT")"

################################################################################
# Functions
################################################################################
get_species_list() {
	case "$SET" in
	all | each)
		[[ -s "$CSV_FOR_SETS" ]] ||
			_die "Species CSV missing: $CSV_FOR_SETS"
		csvtk -t cut -f species "$CSV_FOR_SETS" | tail -n +2
		;;
	some)
		printf '%s\n' Aegilops_umbellulata Anthoceros_agrestis Eucalyptus_pauciflora
		;;
	test)
		printf '%s\n' Aegilops_umbellulata Anthoceros_agrestis
		;;
	file:*)
		local f="${1#file:}"
		[[ -s "$f" ]] || _die "Species list not found: $f"
		cat "$f"
		;;
	*)
		echo "$SET"
		;;
	esac
}

find_fastq_for_species() {
	local sp="$1" base="${ROOT}/${sp}/${TIER}/${INUM}"
	local fq
	for fq in \
		"${base}/summary-data/l.fq.gz" \
		"${base}/summary-data/l.fq"; do
		[[ -s "$fq" ]] && {
			echo "$fq"
			return 0
		}
	done
	find "$base" -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" \) -print -quit 2>/dev/null || true
}

summarize_one() {
	local sp="$1"
	local fq
	fq="$(find_fastq_for_species "$sp" || true)"
	if [[ -z "$fq" ]]; then
		echo -e "${sp}\tNA\tNA\tNA\tNA\tNA\tNA\tNA"
		return 0
	fi

	if ! seqkit stats -T -a "$fq" 2>"${TMPDIR}/seqkit.err" |
		awk -v sp="$sp" -f "${_POLAPLIB_DIR}/scripts/parse-seqkit-stats-v0.1.0.awk" \
			>>"${OUT}.tmp"; then
		_warn "seqkit failed for ${sp}: $(<"${TMPDIR}/seqkit.err")"
		echo -e "${sp}\tERR\tERR\tERR\tERR\tERR\tERR\tERR" >>"${OUT}.tmp"
	fi
}

################################################################################
# Main execution
################################################################################
echo -e "species\ttotal_bases\tread_count\tmean_length\tN50\tmean_Q\tmin_len\tmax_len" >"${OUT}.tmp"

if [[ $PARALLEL -eq 1 && $(command -v parallel) ]]; then
	get_species_list | parallel -j"${THREADS}" summarize_one {}
else
	get_species_list | while read -r sp; do summarize_one "$sp"; done
fi

mv "${OUT}.tmp" "${OUT}"
_info "Wrote: ${OUT}"
