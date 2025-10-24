#!/usr/bin/env bash
set -euo pipefail

# polap-launch-autotune.sh
# Load a flat YAML config -> apply CLI overrides (CLI wins) -> call polap-bash-autotune-oatk.sh.
#
# Addressing (one of):
#   --config-path FILE.yaml
#   --config-dir DIR --preset NAME         # resolves DIR/NAME.yaml (default DIR=~/.polap/profiles)
#
# Any flat YAML key can be overridden on CLI by replacing '_' with '-'.
# Canonical output locator is out_prefix (or --out-prefix); outdir = dirname(out_prefix).
#
# Anchors policy (WGS):
#   Use only --mt-anchors / --pt-anchors (never a generic --anchors). Selection is label-aware.
#
# Verbosity:
#   -v / --verbose    (repeatable; -vv enables bash xtrace with file:func:line PS4)

# ---------- library paths ----------
POLAPLIB_DIR="${_POLAPLIB_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)}"
LOAD_PY="${POLAPLIB_DIR}/polap-py-config-load.py"
AUTOTUNE_SH="${POLAPLIB_DIR}/polap-bash-autotune-oatk.sh"
LOG_LIB="${POLAPLIB_DIR}/polap-bash-oatk-log.sh"

# ---------- source logging ----------
if [[ -r "$LOG_LIB" ]]; then
	# shellcheck disable=SC1090
	source "$LOG_LIB"
else
	# tiny fallback if log lib missing
	OATK_LOG_LEVEL=1
	oatk_log() { ((OATK_LOG_LEVEL >= ${1:-1})) && shift && printf '[%s:%s:%s] %s\n' "${BASH_SOURCE[1]##*/}" "${FUNCNAME[1]:-main}" "${BASH_LINENO[0]:-0}" "$*" >&2; }
	oatk_run() {
		oatk_log 1 "exec: $*"
		"$@"
		local st=$?
		oatk_log 1 "exit($st): $*"
		return $st
	}
	oatk_log_set_level() { OATK_LOG_LEVEL="${1:-1}"; }
	oatk_enable_trace_if_needed() { :; }
	oatk_filter_overrides() { printf '%s\n' "$@"; }
	oatk_provenance_line() { printf '[provenance] config=%s overrides=%s\n' "${1:-?}" "${*:2:-'(none)'}"; }
	oatk_provenance_line1() { printf '[provenance] config=%s overrides=%s\n' "${1:-?}" "${*:2:-'(none)'}"; }
	oatk_provenance_append() { :; }
fi

# ---------- parse addressing + verbosity + collect all args ----------
CFG_PATH="" # new: --config-path (was --path)
CONFIG_DIR=""
PRESET=""

LOG_COUNT=0
ARGS=()
while [[ $# -gt 0 ]]; do
	case "$1" in
	--config-path)
		CFG_PATH="${2}"
		shift 2
		;;
	--config-dir)
		CONFIG_DIR="${2}"
		shift 2
		;;
	--preset)
		PRESET="${2}"
		shift 2
		;;
	-v | --verbose)
		LOG_COUNT=$((LOG_COUNT + 1))
		shift
		;;
	*)
		ARGS+=("$1")
		shift
		;;
	esac
done

# set logging level & trace
oatk_log_set_level "$LOG_COUNT"
oatk_enable_trace_if_needed

# ---------- minimal deps ----------
need() { command -v "$1" >/dev/null 2>&1 || {
	oatk_log 1 "$1 missing"
	exit 127
}; }
need python3
[[ -f "$LOAD_PY" ]] || {
	oatk_log 1 "not found: $LOAD_PY"
	exit 127
}

# default config dir if not provided
if [[ -z "$CFG_PATH" && -z "$CONFIG_DIR" ]]; then
	CONFIG_DIR="${HOME}/.polap/profiles"
	oatk_log 2 "default CONFIG_DIR=${CONFIG_DIR}"
fi

# ---------- load YAML -> PCFG_* ----------
RESOLVED_CFG=""
if [[ -n "$CFG_PATH" ]]; then
	oatk_log 1 "loading YAML --config-path ${CFG_PATH}"
	eval "$(
		python3 "$LOAD_PY" --path "$CFG_PATH" --format env --prefix PCFG 2>/dev/null
	)" || {
		oatk_log 1 "failed to load config: $CFG_PATH"
		exit 4
	}
	RESOLVED_CFG="$CFG_PATH"
else
	[[ -n "$PRESET" ]] || {
		oatk_log 1 "need --preset when --config-path is not given"
		exit 2
	}
	if [[ -z "$CONFIG_DIR" ]]; then CONFIG_DIR="${HOME}/.polap/profiles"; fi
	oatk_log 1 "loading YAML --config-dir ${CONFIG_DIR} --preset ${PRESET}"
	eval "$(
		python3 "$LOAD_PY" --config-dir "$CONFIG_DIR" --preset "$PRESET" --format env --prefix PCFG 2>/dev/null
	)" || {
		oatk_log 1 "failed to load config: ${CONFIG_DIR}/${PRESET}.yaml"
		exit 4
	}
	RESOLVED_CFG="${CONFIG_DIR}/${PRESET}.yaml"
fi

# ---------- apply ALL CLI overrides (CLI wins) ----------
# any --foo-bar becomes PCFG_FOO_BAR; plus boolean toggles without values
apply_overrides() {
	local j=0
	while [[ $j -lt ${#ARGS[@]} ]]; do
		local tok="${ARGS[$j]}"
		case "$tok" in
		# toggles (no value; normalize to true/false if the script uses those keys)
		--wgs-mode)
			PCFG_WGS_MODE="true"
			j=$((j + 1))
			continue
			;;
		--no-wgs-mode)
			PCFG_WGS_MODE="false"
			j=$((j + 1))
			continue
			;;
		--hpc)
			PCFG_HPC="true"
			j=$((j + 1))
			continue
			;;
		--no-hpc)
			PCFG_HPC="false"
			j=$((j + 1))
			continue
			;;
		--hpc-qv)
			PCFG_HPC_QV="true"
			j=$((j + 1))
			continue
			;;
		--no-hpc-qv)
			PCFG_HPC_QV="false"
			j=$((j + 1))
			continue
			;;
		--ppr-hpc)
			PCFG_PPR_HPC="true"
			j=$((j + 1))
			continue
			;;
		--no-ppr-hpc)
			PCFG_PPR_HPC="false"
			j=$((j + 1))
			continue
			;;
		--edge-norm)
			PCFG_EDGE_NORM="true"
			j=$((j + 1))
			continue
			;;
		--no-edge-norm)
			PCFG_EDGE_NORM="false"
			j=$((j + 1))
			continue
			;;
		--emit-fastq)
			PCFG_EMIT_FASTQ="true"
			j=$((j + 1))
			continue
			;;
		--no-emit-fastq)
			PCFG_EMIT_FASTQ="false"
			j=$((j + 1))
			continue
			;;
		--keep-intermediate)
			PCFG_KEEP_INTERMEDIATE="true"
			j=$((j + 1))
			continue
			;;
		--no-keep-intermediate)
			PCFG_KEEP_INTERMEDIATE="false"
			j=$((j + 1))
			continue
			;;
		--redo)
			PCFG_REDO="true"
			j=$((j + 1))
			continue
			;;
		--no-redo)
			PCFG_REDO="false"
			j=$((j + 1))
			continue
			;;

		# K/V mapping
		--*)
			local key="${tok#--}"
			key="${key//-/_}"
			local val="${ARGS[$((j + 1))]:-}"
			[[ -z "$val" || "$val" == --* ]] && {
				oatk_log 1 "missing value for $tok"
				exit 2
			}
			case "$tok" in
			--mt-anchors) PCFG_MT_ANCHORS="$val" ;;
			--pt-anchors) PCFG_PT_ANCHORS="$val" ;;
			--out-prefix) PCFG_OUT_PREFIX="$val" ;;
			--label) PCFG_LABEL="$val" ;;
			*)
				local varname="PCFG_${key^^}"
				printf -v "$varname" "%s" "$val"
				;;
			esac
			j=$((j + 2))
			;;
		*)
			j=$((j + 1))
			;;
		esac
	done
}
apply_overrides

# ---------- normalize core fields ----------
LABEL="${PCFG_LABEL:-mt}"
[[ -n "${PCFG_READS:-}" ]] || {
	oatk_log 1 "reads missing (YAML or --reads)"
	exit 3
}
[[ -n "${PCFG_HMM_DB:-}" ]] || {
	oatk_log 1 "hmm_db missing (YAML or --hmm-db)"
	exit 3
}
[[ -n "${PCFG_OUT_PREFIX:-}" ]] || {
	oatk_log 1 "out_prefix missing (YAML or --out-prefix)"
	exit 3
}

OUTDIR="$(dirname -- "${PCFG_OUT_PREFIX}")"
mkdir -p "$OUTDIR"

# anchors must be label-aware in WGS mode
PCFG_WGS_MODE="${PCFG_WGS_MODE:-true}"
if [[ "$PCFG_WGS_MODE" == "true" ]]; then
	if [[ "$LABEL" == "mt" ]]; then
		[[ -n "${PCFG_MT_ANCHORS:-}" ]] || {
			oatk_log 1 "wgs_mode=true but mt_anchors missing (--mt-anchors)"
			exit 3
		}
	else
		[[ -n "${PCFG_PT_ANCHORS:-}" ]] || {
			oatk_log 1 "wgs_mode=true but pt_anchors missing (--pt-anchors)"
			exit 3
		}
	fi
fi

# threads/k/s/hpc fallbacks
PCFG_THREADS="${PCFG_THREADS:-16}"
PCFG_HPC="${PCFG_HPC:-false}"
PCFG_K="${PCFG_K:-$([[ "$PCFG_HPC" == "true" ]] && echo 41 || echo 121)}"
PCFG_S="${PCFG_S:-$([[ "$PCFG_HPC" == "true" ]] && echo 21 || echo 27)}"

# derive platform preset
case "${PCFG_PRESET,,}" in
hifi | ont) TUNE_PRESET="${PCFG_PRESET,,}" ;;
*) TUNE_PRESET=$([[ "$PCFG_HPC" == "true" ]] && echo "ont" || echo "hifi") ;;
esac

# optional knobs
PCFG_FINAL_IF="${PCFG_FINAL_IF:-}"
PCFG_C_RADIUS="${PCFG_C_RADIUS:-0,10,20}"
PCFG_MIN_SHARED="${PCFG_MIN_SHARED:-5}"
PCFG_JACCARD_MIN="${PCFG_JACCARD_MIN:-$([[ "$TUNE_PRESET" == "ont" ]] && echo 0.0075 || echo 0.015)}"
PCFG_TOPK_NEI="${PCFG_TOPK_NEI:-$([[ "$TUNE_PRESET" == "ont" ]] && echo 40 || echo 50)}"
PCFG_STEPS="${PCFG_STEPS:-$([[ "$TUNE_PRESET" == "ont" ]] && echo 1 || echo 2)}"
PCFG_MAX_RUNS="${PCFG_MAX_RUNS:-12}"

# ---------- summary ----------
oatk_log 1 "config=${RESOLVED_CFG}"
oatk_log 1 "preset=${PCFG_PRESET:-} platform=${TUNE_PRESET} label=${LABEL}"
oatk_log 1 "wgs_mode=${PCFG_WGS_MODE} reads=${PCFG_READS}"
[[ "$PCFG_WGS_MODE" == "true" && "$LABEL" == "mt" ]] && oatk_log 1 "mt_anchors=${PCFG_MT_ANCHORS}"
[[ "$PCFG_WGS_MODE" == "true" && "$LABEL" == "pt" ]] && oatk_log 1 "pt_anchors=${PCFG_PT_ANCHORS}"
oatk_log 1 "hmm_db=${PCFG_HMM_DB} out_prefix=${PCFG_OUT_PREFIX} outdir=${OUTDIR}"
oatk_log 1 "threads=${PCFG_THREADS} k=${PCFG_K} s=${PCFG_S} hpc=${PCFG_HPC}"
[[ -n "${PCFG_FINAL_IF}" ]] && oatk_log 1 "final_if=${PCFG_FINAL_IF}"
oatk_log 2 "c_radius=${PCFG_C_RADIUS} min_shared=${PCFG_MIN_SHARED} jaccard_min=${PCFG_JACCARD_MIN} topk_nei=${PCFG_TOPK_NEI} steps=${PCFG_STEPS} max_runs=${PCFG_MAX_RUNS}"

# provenance append-only
filtered_overrides=$(oatk_filter_overrides "${ARGS[@]}")
oatk_provenance_line "$RESOLVED_CFG" $filtered_overrides
oatk_provenance_append "${OUTDIR}/provenance.autotune.txt" "$RESOLVED_CFG" $filtered_overrides

# ---------- build and run child command (log full cmd) ----------
CMD=(bash "$AUTOTUNE_SH")

# mode
if [[ "$PCFG_WGS_MODE" == "true" ]]; then
	CMD+=(--wgs-mode)
else
	CMD+=(--no-wgs-mode)
fi

# core
CMD+=(--reads "$PCFG_READS")
CMD+=(--label "$LABEL")
CMD+=(--hmm-db "$PCFG_HMM_DB")
CMD+=(--out-prefix "$PCFG_OUT_PREFIX")
CMD+=(--preset "$TUNE_PRESET")
CMD+=(--threads "$PCFG_THREADS")
CMD+=(--k "$PCFG_K" --s "$PCFG_S")
CMD+=(--c-radius "$PCFG_C_RADIUS")
CMD+=(--max-runs "$PCFG_MAX_RUNS")

# anchors (label-aware; no generic --anchors)
if [[ "$PCFG_WGS_MODE" == "true" ]]; then
	CMD+=(--mt-anchors "$PCFG_MT_ANCHORS")
	CMD+=(--pt-anchors "$PCFG_PT_ANCHORS")
fi

# optionals
[[ -n "${PCFG_FINAL_IF}" ]] && CMD+=(--final-if "$PCFG_FINAL_IF")
[[ -n "${PCFG_MIN_SHARED}" ]] && CMD+=(--min-shared "$PCFG_MIN_SHARED")
[[ -n "${PCFG_JACCARD_MIN}" ]] && CMD+=(--jaccard-min "$PCFG_JACCARD_MIN")
[[ -n "${PCFG_TOPK_NEI}" ]] && CMD+=(--topk-nei "$PCFG_TOPK_NEI")
[[ -n "${PCFG_STEPS}" ]] && CMD+=(--steps "$PCFG_STEPS")
[[ -n "${PCFG_NUC_CUT_LOG10:-}" ]] && CMD+=(--nuc-cut-log10 "$PCFG_NUC_CUT_LOG10")
[[ -n "${PCFG_X_SLOPE:-}" ]] && CMD+=(--x-slope "$PCFG_X_SLOPE")

# Append -v LOG_COUNT times
if [[ -n "${LOG_COUNT:-}" && "$LOG_COUNT" -gt 0 ]]; then
	for ((i = 0; i < LOG_COUNT; i++)); do
		CMD+=(-v)
	done
fi

# log + run child
oatk_run "${CMD[@]}"
status=$?

oatk_log 1 "autotune finished with status=$status"
exit "$status"
