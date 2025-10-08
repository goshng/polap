#!/usr/bin/env bash
################################################################################
# polap-bash-report-v0.1.0.sh
#
# Version : v0.1.0
# Author  : Sang Chul Choi
# Date    : 2025-10-05
# License : GPL-3.0+
#
# Purpose :
#   Orchestrate report generation for POLAP:
#     - Build tables (md/table-*.tsv, md/*.md)
#     - Build figures (Bandage PNGs + R plots)
#     - Build manuscript (PDF/HTML via Makefile/Quarto)
#
# Usage   :
#   polap-bash-report.sh tables   [--set SET] [--csv CSV] [--cache]
#   polap-bash-report.sh figures  [--set SET] [--csv CSV] [--out DIR]
#                                 [--type pt|mt] [--no-bandage]
#                                 [--cache] [--parallel] [--view]
#   polap-bash-report.sh manuscript [--type pdf|html] [--out DIR] [--view]
#   polap-bash-report.sh --version | --help
#
# Notes   :
#   - All dependent scripts live under "\$(_POLAPLIB_DIR)/scripts/".
#   - _POLAPLIB_DIR is set to the base directory of this script automatically.
#   - Requires: bash ≥ 4, awk, sed, grep, coreutils, csvtk, Rscript, python3, Bandage, make
################################################################################
set -euo pipefail
IFS=$'\n\t'

# Resolve _POLAPLIB_DIR to this script's base folder
# shellcheck disable=SC2164
_POLAPLIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# delete
# export _POLAPLIB_DIR

SCRIPT_VERSION="v0.1.0"

# Defaults
SET_NAME="some"                                 # all|some|test|each|file:<path>|<species>
IN_CSV="${_POLAPLIB_DIR}/polap-bash-report.csv" # must include 'species' if used
OUTDIR="${_POLAPLIB_DIR}/man/v0.5.4/figures"
TYPE="pt"         # pt|mt (used by figures/manuscript)
BANDAGE="bandage" # bandage|no-bandage
CACHE=0
PARALLEL=0
VIEW=0

# Defensive helpers
die() {
	echo "[ERR] $*" >&2
	exit 1
}
warn() { echo "[WARN] $*" >&2; }
info() { echo "[INFO] $*" >&2; }
have() { command -v "$1" >/dev/null 2>&1; }

usage() {
	# sed -n '1,80p' "${BASH_SOURCE[0]}" | sed -n '1,/^################################################################################/p' | sed -n '1,/^set -euo pipefail/p' | sed -n '1,40p'

	cat <<'EOF'
polap-bash-report.sh tables   [--set SET] [--csv CSV] [--cache]
polap-bash-report.sh figures  [--set SET] [--csv CSV] [--out DIR]
                              [--type pt|mt] [--no-bandage]
                              [--cache] [--parallel] [--view]
polap-bash-report.sh manuscript [--type pdf|html] [--out DIR] [--view]

Options:
  --set SET        Species set: all|some|test|each|file:<path>|<species>
  --csv FILE       Config CSV (must include 'species' column for all/each)
  --out DIR        Figures/manuscript output directory
  --type TYPE      pt|mt (figures/manuscript)
  --no-bandage     Skip Bandage PNG rendering; create placeholders instead
  --cache          Skip steps if outputs newer than inputs
  --parallel       Use parallel/xargs -P for per-species operations
  --view           Open resulting artifact (xdg-open/open)
  --help           Show this help and exit
  --version        Show version and exit
EOF
}

# Parse top-level command
CMD="${1:-}"
shift || true
case "${CMD}" in
tables | figures | manuscript | "") : ;;
--help | -h)
	usage
	exit 0
	;;
--version)
	echo "${SCRIPT_VERSION}"
	exit 0
	;;
*) die "Unknown subcommand: ${CMD}" ;;
esac

# Parse flags
while (($#)); do
	case "$1" in
	--set)
		SET_NAME="${2:?}"
		shift 2
		;;
	--csv)
		IN_CSV="${2:?}"
		shift 2
		;;
	--out)
		OUTDIR="${2:?}"
		shift 2
		;;
	--type)
		TYPE="${2:?}"
		shift 2
		;;
	--no-bandage)
		BANDAGE="no-bandage"
		shift
		;;
	--cache)
		CACHE=1
		shift
		;;
	--parallel)
		PARALLEL=1
		shift
		;;
	--view)
		VIEW=1
		shift
		;;
	--help | -h)
		usage
		exit 0
		;;
	--version)
		echo "${SCRIPT_VERSION}"
		exit 0
		;;
	*) die "Unknown option: $1" ;;
	esac
done

# Checks
have awk || die "awk not found"
have sed || die "sed not found"
have grep || die "grep not found"
have python3 || die "python3 not found"
have Rscript || die "Rscript not found"
have csvtk || warn "csvtk not found; 'all/each' sets will not work"
have make || warn "make not found; manuscript building may fail"
# Bandage is optional unless BANDAGE=bandage
[[ "${BANDAGE}" == "no-bandage" ]] || have polap.sh || warn "polap.sh not in PATH; Bandage PNGs may fail (we will fallback/placeholder)."

# Utility: species list
get_species_list() {
	local set="$1"
	case "$set" in
	all | each)
		[[ -s "${IN_CSV}" ]] || die "CSV missing: ${IN_CSV}"
		csvtk -t cut -f species "${IN_CSV}" | tail -n +2
		;;
	some)
		printf '%s\n' Aegilops_umbellulata Anthoceros_agrestis Eucalyptus_pauciflora
		;;
	test)
		printf '%s\n' Aegilops_umbellulata Anthoceros_agrestis
		;;
	file:*)
		local path="${set#file:}"
		[[ -s "${path}" ]] || die "List not found: ${path}"
		cat "${path}"
		;;
	*)
		# treat as a single species literal
		printf '%s\n' "${set}"
		;;
	esac
}

# need_build: return 0 if we should (re)build tgt (cache miss), 1 to skip
need_build() {
	local tgt="$1"
	shift || true
	[[ ${CACHE} -eq 0 ]] && return 0
	[[ ! -e "$tgt" ]] && return 0
	local dep
	for dep in "$@"; do
		[[ -e "$dep" && "$dep" -nt "$tgt" ]] && return 0
	done
	return 1
}

# Subcommands
run_tables() {
	mkdir -p "${OUTDIR}"
	local out_tsv="${OUTDIR}/table-${SET_NAME}-0.tsv"

	local csv_file="${IN_CSV}"
	[[ -s "$csv_file" ]] || csv_file="${_POLAPLIB_DIR}/${csv_file}"

	if [[ ! -f "$csv_file" ]]; then
		echo "[ERROR] CSV file not found: $csv_file" >&2
		return 0
	fi

	if need_build "${out_tsv}" "${IN_CSV}"; then
		python3 "${_POLAPLIB_DIR}/scripts/build_benchmark_tables.py" \
			--csv "${csv_file}" \
			--out "${OUTDIR}" \
			--species-set "${SET_NAME}"
	fi
	echo "${out_tsv}"
}

render_png_for_species() {
	local sp="$1"
	bash "${_POLAPLIB_DIR}/scripts/gfa_to_png-v0.1.0.sh" "$sp" "$TYPE" "$BANDAGE" || true
}

run_figures() {
	mkdir -p "${OUTDIR}"
	# 1) ensure assembly PNGs exist (Bandage or placeholders)
	if [[ ${PARALLEL} -eq 1 ]] && have parallel; then
		get_species_list "${SET_NAME}" | parallel -j"$(nproc)" render_png_for_species {}
	else
		get_species_list "${SET_NAME}" | xargs -r -n1 -P"$(nproc)" -I{} bash -c 'render_png_for_species "$@"' _ {}
	fi

	# 2) benchmark plots from TSV
	local tsv="${_POLAPLIB_DIR}/md/table-${SET_NAME}-0.tsv"
	[[ -s "${tsv}" ]] || tsv="$(run_tables)"
	Rscript --vanilla "${_POLAPLIB_DIR}/scripts/make_benchmark_plots-v0.1.0.R" \
		--table "${tsv}" \
		--out "${OUTDIR}" \
		--type "time"
	Rscript --vanilla "${_POLAPLIB_DIR}/scripts/make_benchmark_plots-v0.1.0.R" \
		--table "${tsv}" \
		--out "${OUTDIR}" \
		--type "memory"

	# 3) manifest (best-effort)
	bash "${_POLAPLIB_DIR}/scripts/build_manifest-v0.1.0.sh" \
		--set "${SET_NAME}" \
		--type "${TYPE}" \
		--out "${OUTDIR}" || true

	# Optional view
	if [[ ${VIEW} -eq 1 ]]; then
		local f="${OUTDIR}/benchmark-time.pdf"
		{ command -v xdg-open >/dev/null && xdg-open "${f}"; } || { command -v open >/dev/null && open "${f}"; } || true
	fi
}

run_manuscript() {
	# Delegate to Makefile living in ${OUTDIR}/.. (same convention as earlier)
	local fmt="${TYPE}" # pdf|html
	local workdir
	workdir="$(cd "${OUTDIR}/.." && pwd)"
	(cd "${workdir}" && make "${SET_NAME}" FORMAT="${fmt}") || warn "make manuscript failed"
	local artifact="${workdir}/manuscript-${SET_NAME}.${fmt}"
	if [[ ${VIEW} -eq 1 && -s "${artifact}" ]]; then
		{ command -v xdg-open >/dev/null && xdg-open "${artifact}"; } || { command -v open >/dev/null && open "${artifact}"; } || true
	fi
	echo "${artifact}"
}

# Main
case "${CMD}" in
tables) run_tables ;;
figures) run_figures ;;
manuscript) run_manuscript ;;
*)
	usage
	exit 2
	;;
esac
