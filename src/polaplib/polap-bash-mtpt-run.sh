#!/usr/bin/env bash
# FILE: polap-bash-mtpt-run.sh
# Version: v0.2.0
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Orchestrate MTPT pipeline:
#   Step1 (per species)  -> Step2 homology (pooled)
#   Step3 PT phylogeny   -> Step4 gain/loss mapping
#
# Features:
#   - Lockfiles per step (safe in --parallel)
#   - --resume: skip when outputs are newer than inputs
#   - Failsafe wiring: polap-lib-failsafe.sh / run-* (optional)
#
# Usage:
#   polap-bash-mtpt-run.sh -H OatkDB/v20230921/embryophyta_pltd.fam \
#     [-b man/analysis] [-t 8] [--parallel 4] [--mode coding|mauve] \
#     [--gap-frac 0.7] [--min-sites 100000] [--resume] \
#     [--skip-step2] [--skip-step3] [--skip-step4]
#
set -Eeuo pipefail
IFS=$'\n\t'

: "${_POLAPLIB_DIR:=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
VERSION="0.2.0"

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
polap-bash-mtpt-run v${VERSION}
Run MTPT pipeline end-to-end with locks + resume.

USAGE:
  \$(basename "\$0") -H [Oatk_plastid.fam] [-b man/analysis] [-t 8] [--parallel 4]
                     [--mode coding|mauve] [--gap-frac 0.7] [--min-sites 100000]
                     [--resume] [--skip-step2] [--skip-step3] [--skip-step4]

Notes:
  - Step 1 is per species (assemblies/<sp>/{mt.fa,pt.fa})
  - Steps 2–4 are pooled across species
  - Oatk_plastid.fam (default): OatkDB/v20230921/embryophyta_pltd.fam
EOF
}

source "${_POLAPLIB_DIR}/polap-lib-log.sh"
_polap_log0 "mtpt-run: starting (log=$POLAP_LOG_FILE, verbose=$POLAP_VERBOSE)"

# ---------- CLI ----------
ANALYSIS_BASE="man/analysis"
THREADS=8
PARALLEL=1
HMM_LIB="OatkDB/v20230921/embryophyta_pltd.fam"
MODE="coding"
GAP_FRAC="0.7"
MIN_SITES="100000"
RESUME=0
SKIP2=0
SKIP3=0
SKIP4=0

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
	--parallel)
		PARALLEL="$2"
		shift 2
		;;
	--mode)
		MODE="$2"
		shift 2
		;;
	--gap-frac)
		GAP_FRAC="$2"
		shift 2
		;;
	--min-sites)
		MIN_SITES="$2"
		shift 2
		;;
	--resume)
		RESUME=1
		shift
		;;
	--skip-step2)
		SKIP2=1
		shift
		;;
	--skip-step3)
		SKIP3=1
		shift
		;;
	--skip-step4)
		SKIP4=1
		shift
		;;
	-h | --help)
		usage
		exit 0
		;;
	*)
		echo "Unknown: $1" >&2
		usage
		exit 1
		;;
	esac
done

[[ -z "${HMM_LIB}" ]] && {
	echo "ERROR: -H <Oatk .fam> is required." >&2
	usage
	exit 2
}

# ---------- Optional libraries (w/ graceful fallbacks) ----------
# These are your libraries; if absent, we provide built-ins.
# if [[ -r "${_POLAPLIB_DIR}/polap-lib-failsafe.sh" ]]; then
# 	# shellcheck disable=SC1091
# 	source "${_POLAPLIB_DIR}/polap-lib-failsafe.sh"
# fi

[[ -r "${_POLAPLIB_DIR}/polap-lib-run-simple.sh" ]] && source "${_POLAPLIB_DIR}/polap-lib-run-simple.sh" || true
[[ -r "${_POLAPLIB_DIR}/polap-lib-run-argv.sh" ]] && source "${_POLAPLIB_DIR}/polap-lib-run-argv.sh" || true

# ---------- Minimal shims if lib funcs are absent ----------
# polap_log() { printf '%s %s\n' "[$(date +'%F %T')]" "$*" >&2; }
polap_log() {
	_polap_log0 "$*"
}

need() { command -v "$1" >/dev/null 2>&1 || {
	echo "Missing $1" >&2
	exit 127
}; }

# Lock acquire/release (prefer your lib if provided)
_lock_mkdir() {
	local l="$1"
	mkdir "$l" 2>/dev/null && {
		echo $$ >"$l/pid"
		return 0
	}
	return 1
}

_lock_release() {
	local l="$1"
	[[ -d "$l" ]] && rm -rf "$l"
}

polap_lock_acquire_wrapper() {
	local lock="$1" wait_s="${2:-0}"
	if declare -F polap_lock_acquire >/dev/null; then
		polap_lock_acquire "$@"
	else
		local t0=$(date +%s)
		while ! _lock_mkdir "$lock"; do
			((wait_s == 0)) && return 1
			sleep 1
			(($(date +%s) - t0 >= wait_s)) && return 2
		done
		return 0
	fi
}
polap_lock_release_wrapper() {
	local lock="$1"
	if declare -F polap_lock_release >/dev/null; then
		polap_lock_release "$lock"
	else
		_lock_release "$lock"
	fi
}

# Run wrapper (prefer your run helpers)
polap_run_wrapper() {
	# Usage: polap_run_wrapper <tag> <cmd ...argv...>
	local tag="$1"
	shift

	if declare -F polap_run_argv >/dev/null 2>&1; then
		# safest: no eval; pass argv as-is
		polap_run_argv --tag "$tag" -- "$@"

	elif declare -F polap_run_simple >/dev/null 2>&1; then
		# fallback: string runner (only when you need shell features like pipes)
		# join argv into a preview string, but still execute the argv via eval string
		local cmd
		cmd="$(printf '%q ' "$@")"
		cmd="${cmd# }"
		polap_run_simple --tag "$tag" -- "$cmd"

	else
		# last resort: run directly, but at least log with tag
		_polap_log0 "[RUN:${tag}] $(printf '%q ' "$@")"
		"$@"
	fi
}

# Up-to-date check:
#   uptodate <out1> [out2 ...] -- <in1> [in2 ...]
# returns 0 if all outputs exist and oldest(out) >= newest(in); else 1
uptodate() {
	local outs=() ins=() seen_sep=0
	for a in "$@"; do
		if [[ "$a" == "--" ]]; then
			seen_sep=1
			continue
		fi
		((seen_sep == 0)) && outs+=("$a") || ins+=("$a")
	done
	# All outputs must exist
	for o in "${outs[@]}"; do [[ -e "$o" ]] || return 1; done
	# Compute min(out mtime) and max(in mtime)
	local min_out=
	local max_in=0 m
	for o in "${outs[@]}"; do
		m=$(stat -c %Y "$o" 2>/dev/null || echo 0)
		if [[ -z "$min_out" || "$m" -lt "$min_out" ]]; then min_out="$m"; fi
	done
	for i in "${ins[@]}"; do
		[[ -e "$i" ]] || continue
		m=$(stat -c %Y "$i" 2>/dev/null || echo 0)
		((m > max_in)) && max_in="$m"
	done
	[[ -n "$min_out" && "$min_out" -ge "$max_in" ]]
}

# ---------- Paths ----------
BASE="${ANALYSIS_BASE}"
ASM="${BASE}/assemblies"
LOCKS="${BASE}/locks"
LOGD="${BASE}/logs"
mkdir -p "${LOCKS}" "${LOGD}"

# ---------- Tool sanity ----------
need bash
need python3
need Rscript

# ---------- Detect species (must have both mt.fa and pt.fa) ----------
mapfile -t SPECIES < <(find "${ASM}" -mindepth 1 -maxdepth 1 -type d -print0 |
	xargs -0 -I{} bash -lc 'sp=$(basename "{}"); [[ -s "{}/mt.fa" && -s "{}/pt.fa" ]] && echo "$sp"' | sort)
((${#SPECIES[@]})) || {
	echo "No species with both mt.fa and pt.fa under ${ASM}" | tee "${LOGD}/run_mtpt.log"
	exit 3
}

# ---------- Step 1: per-species detector ----------
run_step1_one() {
	local sp="$1"
	local lock="${LOCKS}/step1_${sp}.lock"
	local out1="${BASE}/mtpt_calls/${sp}/mtpt.tsv"
	local out2="${BASE}/mtpt_calls/${sp}/fasta" # dir
	local in1="${ASM}/${sp}/mt.fa"
	local in2="${ASM}/${sp}/pt.fa"
	local in3="${HMM_LIB}"

	# Resume: skip if up-to-date
	if ((RESUME == 1)) && uptodate "${out1}" "${out2}" -- "${in1}" "${in2}" "${in3}"; then
		_polap_log0 "[Step1:${sp}] up-to-date -> skip"
		return 0
	fi

	# Lock
	if ! polap_lock_acquire_wrapper "${lock}" 1; then
		_polap_log0 "[Step1:${sp}] lock busy, skipping this attempt"
		return 0
	fi

	# Double-check after acquiring lock (avoid races)
	if ((RESUME == 1)) && uptodate "${out1}" "${out2}" -- "${in1}" "${in2}" "${in3}"; then
		_polap_log0 "[Step1:${sp}] up-to-date after lock -> skip"
		polap_lock_release_wrapper "${lock}"
		return 0
	fi

	# Run
	_polap_log0 "[Step1:${sp}] running detector"

	if ! polap_run_wrapper "step1:${sp}" \
		bash "${_POLAPLIB_DIR}/polap-bash-mtpt-step1-mtpt-detector.sh" \
		-s "${sp}" -H "${HMM_LIB}" -b "${BASE}" -t "${THREADS}" \
		--min-len 150 --min-pid 0.85; then
		_polap_log0 "[Step1:${sp}] FAILED"
		touch "${LOGD}/step1_${sp}.failed"
		polap_lock_release_wrapper "${lock}"
		return 1
	fi

	touch "${LOGD}/step1_${sp}.ok"
	polap_lock_release_wrapper "${lock}"
	_polap_log0 "[Step1:${sp}] done"
}

_polap_log0 "[Step1] MTPT detector per species (n=${#SPECIES[@]}; P=${PARALLEL})"
export -f run_step1_one polap_log polap_run_wrapper polap_lock_acquire_wrapper polap_lock_release_wrapper uptodate
export _POLAPLIB_DIR BASE ASM HMM_LIB THREADS RESUME LOCKS LOGD

if ((PARALLEL > 1)); then
	printf "%s\n" "${SPECIES[@]}" | xargs -P "${PARALLEL}" -I{} bash -lc 'run_step1_one "$@"' _ {}
else
	for sp in "${SPECIES[@]}"; do run_step1_one "${sp}"; done
fi

# ---------- Step 2: homology (pooled) ----------
if ((SKIP2 == 0)); then
	_polap_log0 "[Step2] homology clustering check"
	CLDIR="${BASE}/clusters"
	mkdir -p "${CLDIR}"
	local_out="${CLDIR}/cluster_map.tsv"
	# Inputs: all tract FASTAs
	mapfile -t tract_fastas < <(find "${BASE}/mtpt_calls" -type f -name "*.fa" | sort)
	if ((RESUME == 1)) && uptodate "${local_out}" -- "${tract_fastas[@]}"; then
		_polap_log0 "[Step2] up-to-date -> skip"
	else
		local_lock="${LOCKS}/step2.lock"
		if polap_lock_acquire_wrapper "${local_lock}" 10; then
			# re-evaluate after lock
			if ((RESUME == 1)) && uptodate "${local_out}" -- "${tract_fastas[@]}"; then
				_polap_log0 "[Step2] up-to-date after lock -> skip"
			else
				_polap_log0 "[Step2] running"
				if ! polap_run_wrapper "step2" \
					bash "${_POLAPLIB_DIR}/polap-bash-mtpt-step2-homology.sh" \
					-b "${BASE}" -t "${THREADS}"; then
					_polap_log0 "[Step2] FAILED"
					touch "${LOGD}/step2.failed"
					polap_lock_release_wrapper "${local_lock}"
					exit 1
				fi
				touch "${LOGD}/step2.ok"
			fi
			polap_lock_release_wrapper "${local_lock}"
		else
			_polap_log0 "[Step2] lock busy, skipping"
		fi
	fi
else
	_polap_log0 "[Step2] skipped via --skip-step2"
fi

# ---------- Step 3: PT phylogeny (pooled) ----------
if ((SKIP3 == 0)); then
	_polap_log0 "[Step3] PT phylogeny check"
	PTDIR="${BASE}/pt_tree/concat"
	mkdir -p "${PTDIR}"
	local_out_tree="${PTDIR}/iqtree.treefile"
	# Inputs: all pt.fa (+ HMM_LIB if coding mode)
	mapfile -t pt_fastas < <(find "${ASM}" -type f -name "pt.fa" | sort)
	in_list=("${pt_fastas[@]}")
	if [[ "${MODE}" == "coding" ]]; then in_list+=("${HMM_LIB}"); fi

	if ((RESUME == 1)) && uptodate "${local_out_tree}" -- "${in_list[@]}"; then
		_polap_log0 "[Step3] up-to-date -> skip"
	else
		local_lock="${LOCKS}/step3.lock"
		if polap_lock_acquire_wrapper "${local_lock}" 10; then
			if ((RESUME == 1)) && uptodate "${local_out_tree}" -- "${in_list[@]}"; then
				_polap_log0 "[Step3] up-to-date after lock -> skip"
			else
				_polap_log0 "[Step3] running (mode=${MODE})"
				if ! polap_run_wrapper "step3" \
					bash "${_POLAPLIB_DIR}/polap-bash-mtpt-step3-pt-phylogeny.sh" \
					-b "${BASE}" -t "${THREADS}" -H "${HMM_LIB}" --mode "${MODE}" \
					--gap-frac "${GAP_FRAC}" --min-sites "${MIN_SITES}"; then
					_polap_log0 "[Step3] FAILED"
					touch "${LOGD}/step3.failed"
					polap_lock_release_wrapper "${local_lock}"
					exit 1
				fi
				touch "${LOGD}/step3.ok"
			fi
			polap_lock_release_wrapper "${local_lock}"
		else
			_polap_log0 "[Step3] lock busy, skipping"
		fi
	fi
else
	_polap_log0 "[Step3] skipped via --skip-step3"
fi

# ---------- Step 4: Gain/Loss mapping (pooled) ----------
if ((SKIP4 == 0)); then
	_polap_log0 "[Step4] gain/loss check"
	GLDIR="${BASE}/gain_loss"
	mkdir -p "${GLDIR}"
	local_out_ev="${GLDIR}/events.tsv"
	local_in_clmap="${BASE}/clusters/cluster_map.tsv"
	local_in_tree="${BASE}/pt_tree/concat/iqtree.treefile"

	if ((RESUME == 1)) && uptodate "${local_out_ev}" -- "${local_in_clmap}" "${local_in_tree}"; then
		_polap_log0 "[Step4] up-to-date -> skip"
	else
		local_lock="${LOCKS}/step4.lock"
		if polap_lock_acquire_wrapper "${local_lock}" 10; then
			if ((RESUME == 1)) && uptodate "${local_out_ev}" -- "${local_in_clmap}" "${local_in_tree}"; then
				_polap_log0 "[Step4] up-to-date after lock -> skip"
			else
				_polap_log0 "[Step4] running"
				if ! polap_run_wrapper "step4" \
					bash "${_POLAPLIB_DIR}/polap-bash-mtpt-step4-gain-loss.sh" \
					-b "${BASE}" --model all; then
					_polap_log0 "[Step4] FAILED"
					touch "${LOGD}/step4.failed"
					polap_lock_release_wrapper "${local_lock}"
					exit 1
				fi
				touch "${LOGD}/step4.ok"
			fi
			polap_lock_release_wrapper "${local_lock}"
		else
			_polap_log0 "[Step4] lock busy, skipping"
		fi
	fi
else
	_polap_log0 "[Step4] skipped via --skip-step4"
fi

_polap_log0 "Done. See results/ and logs/."
