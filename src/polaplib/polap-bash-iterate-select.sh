#!/usr/bin/env bash
set -euo pipefail

# polap-bash-iterate-select.sh
# Iteratively recruit organelle reads (mt/pt) until gain ratio < epsilon.
#
# Addressing (one of):
#   --config-path FILE.yaml
#   --config-dir DIR --preset NAME     # resolves DIR/NAME.yaml (default: ~/.polap/profiles)
#
# Any flat YAML key can be overridden by CLI (--foo-bar -> YAML key foo_bar).
# Required keys (YAML or CLI): reads, out_prefix, label, (mt_anchors|pt_anchors for first round in WGS)
#
# Loop:
#   A_i = anchors for iteration i (union from previous)
#   select-mt (with A_i) -> C_i (PPR ids)
#   gain = C_i - A_i
#   A_{i+1} = A_i ∪ C_i
#   stop if |gain|/|A_{i+1}| < epsilon  or i >= max_iters
#
# Writes:
#   <outdir>/iter<i>/{anchors.txt,candidates.txt,gain.txt,union.txt}
#   <outdir>/iter-log.csv       (iter,|A_i|,|C_i|,|gain|,ratio)
#   provenance.iterate.txt
#
# Verbosity:
#   -v / --verbose  (repeatable; -vv enables xtrace with file:function:line PS4)

# ---------- paths ----------
POLAPLIB_DIR="${_POLAPLIB_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)}"
LOAD_PY="${POLAPLIB_DIR}/polap-py-config-load.py"
SELECT_MT_SH="${POLAPLIB_DIR}/polap-bash-select-mt.sh"
LOG_LIB="${POLAPLIB_DIR}/polap-bash-oatk-log.sh"

# ---------- logging ----------
if [[ -r "$LOG_LIB" ]]; then
	# shellcheck disable=SC1090
	source "$LOG_LIB"
else
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

die() {
	oatk_log 1 "$*"
	exit 2
}
need() { command -v "$1" >/dev/null 2>&1 || {
	oatk_log 1 "$1 missing"
	exit 127
}; }

# ---------- CLI parse ----------
CFG_PATH=""
CONFIG_DIR=""
PRESET=""
LOG_COUNT=0
# loop knobs
EPSILON="0.01"
MAX_ITERS=8

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
	--epsilon)
		EPSILON="${2}"
		shift 2
		;;
	--max-iters)
		MAX_ITERS="${2}"
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

oatk_log_set_level "$LOG_COUNT"
oatk_enable_trace_if_needed

# build -v flags for cascading to child scripts
VERBOSE_FLAGS=()
if [[ -n "${LOG_COUNT:-}" && "$LOG_COUNT" -gt 0 ]]; then
	for ((i = 0; i < LOG_COUNT; i++)); do VERBOSE_FLAGS+=(-v); done
fi

# ---------- deps ----------
need python3
[[ -f "$LOAD_PY" ]] || die "not found: $LOAD_PY"
[[ -f "$SELECT_MT_SH" ]] || die "not found: $SELECT_MT_SH"

# ---------- load YAML -> ICFG_* ----------
RESOLVED_CFG=""
if [[ -n "$CFG_PATH" ]]; then
	oatk_log 1 "loading YAML --config-path ${CFG_PATH}"
	eval "$(
		python3 "$LOAD_PY" --path "$CFG_PATH" --format env --prefix ICFG 2>/dev/null
	)" || die "failed to load config: $CFG_PATH"
	RESOLVED_CFG="$CFG_PATH"
else
	[[ -n "$PRESET" ]] || die "need --preset when --config-path is not given"
	[[ -n "$CONFIG_DIR" ]] || CONFIG_DIR="${HOME}/.polap/profiles"
	oatk_log 1 "loading YAML --config-dir ${CONFIG_DIR} --preset ${PRESET}"
	eval "$(
		python3 "$LOAD_PY" --config-dir "$CONFIG_DIR" --preset "$PRESET" --format env --prefix ICFG 2>/dev/null
	)" || die "failed to load config: ${CONFIG_DIR}/${PRESET}.yaml"
	RESOLVED_CFG="${CONFIG_DIR}/${PRESET}.yaml"
fi

# ---------- apply overrides (CLI wins) ----------
apply_overrides() {
	local j=0
	while [[ $j -lt ${#ARGS[@]} ]]; do
		local tok="${ARGS[$j]}"
		case "$tok" in
		-v | --verbose)
			j=$((j + 1))
			continue
			;;
		--config-path | --config-dir | --preset | --epsilon | --max-iters)
			j=$((j + 2))
			continue
			;;
		--*)
			local key="${tok#--}"
			key="${key//-/_}"
			local val="${ARGS[$((j + 1))]:-}"
			[[ -z "$val" || "$val" == --* ]] && die "missing value for $tok"
			case "$tok" in
			--out-prefix) ICFG_OUT_PREFIX="$val" ;;
			--label) ICFG_LABEL="$val" ;;
			--reads) ICFG_READS="$val" ;;
			--mt-anchors) ICFG_MT_ANCHORS="$val" ;;
			--pt-anchors) ICFG_PT_ANCHORS="$val" ;;
			*)
				local varname="ICFG_${key^^}"
				printf -v "$varname" "%s" "$val"
				;;
			esac
			j=$((j + 2))
			;;
		*) j=$((j + 1)) ;;
		esac
	done
}
apply_overrides

# ---------- normalize core ----------
LABEL="${ICFG_LABEL:-mt}"
[[ "$LABEL" == "mt" || "$LABEL" == "pt" ]] || die "label must be mt|pt"
[[ -n "${ICFG_READS:-}" ]] || die "reads missing (YAML or --reads)"
[[ -n "${ICFG_OUT_PREFIX:-}" ]] || die "out_prefix missing (YAML or --out-prefix)"

OUTDIR="$(dirname -- "${ICFG_OUT_PREFIX}")"
mkdir -p "$OUTDIR"

# anchors for iteration 0
if [[ "$LABEL" == "mt" ]]; then
	[[ -n "${ICFG_MT_ANCHORS:-}" && -s "${ICFG_MT_ANCHORS}" ]] || die "mt_anchors (--mt-anchors) missing/empty"
	A0="${ICFG_MT_ANCHORS}"
else
	[[ -n "${ICFG_PT_ANCHORS:-}" && -s "${ICFG_PT_ANCHORS}" ]] || die "pt_anchors (--pt-anchors) missing/empty"
	A0="${ICFG_PT_ANCHORS}"
fi

# threads/preset defaults (pass down to select-mt)
ICFG_THREADS="${ICFG_THREADS:-16}"
ICFG_PRESET="${ICFG_PRESET:-hifi}"

# ---------- helpers ----------
uniq_ids() { awk 'NF{print $1}' "$1" | tr -d '\r' | LC_ALL=C sort -u; }
union_ids() { LC_ALL=C sort -u "$1" "$2"; }
diff_ids() { LC_ALL=C comm -23 "$1" "$2"; }
count_ids() { LC_ALL=C awk 'END{print NR}' "$1"; }

# provenance
filtered_overrides=$(oatk_filter_overrides "${ARGS[@]}")
oatk_provenance_line "$RESOLVED_CFG" $filtered_overrides
oatk_provenance_append "${OUTDIR}/provenance.iterate.txt" "$RESOLVED_CFG" $filtered_overrides

# summary CSV
ILOG="${OUTDIR}/${LABEL}.iter-log.csv"
if [[ ! -s "$ILOG" ]]; then
	echo "iter,|A_i|,|C_i|,|gain|,ratio" >"$ILOG"
fi

# ---------- iterate ----------
A_CUR="${OUTDIR}/${LABEL}.iter0.anchors.txt"
uniq_ids "$A0" >"$A_CUR"
iter=0

while [[ $iter -lt ${MAX_ITERS} ]]; do
	ITERDIR="${OUTDIR}/iter${iter}"
	mkdir -p "$ITERDIR"

	# snapshot A_i
	cp -f "$A_CUR" "${ITERDIR}/anchors.txt"

	# run select-mt with A_i as anchors (label-aware)
	RECR_PREFIX="${OUTDIR}/${LABEL}.wgs.iter${iter}"
	CMD=(bash "$SELECT_MT_SH" select-mt
		--config-dir "${CONFIG_DIR:-}"
		--preset "${ICFG_PRESET}"
		--label "$LABEL"
		--reads "${ICFG_READS}"
		--out "$RECR_PREFIX"
		--threads "${ICFG_THREADS}"
		--ensemble majority
		--emit-fastq
		"${VERBOSE_FLAGS[@]}"
	)
	if [[ "$LABEL" == "mt" ]]; then
		CMD+=(--mt-anchors "${A_CUR}")
	else
		CMD+=(--pt-anchors "${A_CUR}")
	fi

	oatk_run "${CMD[@]}" || true

	# candidates C_i: prefer PPR ids; fallback to xcover band ids
	C_RAW="${ITERDIR}/candidates.txt"
	if [[ -s "${RECR_PREFIX}.ppr.${LABEL}.ids" ]]; then
		uniq_ids "${RECR_PREFIX}.ppr.${LABEL}.ids" >"$C_RAW"
	elif [[ -s "${RECR_PREFIX}.xcover.${LABEL}.ids" ]]; then
		uniq_ids "${RECR_PREFIX}.xcover.${LABEL}.ids" >"$C_RAW"
	else
		: >"$C_RAW"
	fi

	# gain = C_i - A_i
	GAIN="${ITERDIR}/gain.txt"
	diff_ids "$C_RAW" "$A_CUR" >"$GAIN"

	# A_{i+1} = A_i ∪ C_i
	A_NEXT="${OUTDIR}/${LABEL}.iter$((iter + 1)).anchors.txt"
	union_ids "$A_CUR" "$C_RAW" >"$A_NEXT"

	# ratios
	nA=$(count_ids "$A_CUR")
	nC=$(count_ids "$C_RAW")
	nG=$(count_ids "$GAIN")
	nNext=$(count_ids "$A_NEXT")
	ratio="0"
	if [[ "$nNext" -gt 0 ]]; then
		ratio=$(echo "scale=6; $nG / $nNext" | bc -l)
	fi
	printf "%d,%d,%d,%d,%.6f\n" "$iter" "$nA" "$nC" "$nG" "$ratio" >>"$ILOG"
	oatk_log 1 "[iter $iter] |A|=$nA |C|=$nC |gain|=$nG ratio=$ratio"

	# decide stop
	ok=$(echo "$ratio < $EPSILON" | bc -l || echo 0)
	if [[ "$ok" -eq 1 ]]; then
		oatk_log 1 "[stop] ratio=$ratio < epsilon=$EPSILON"
		break
	fi

	# next
	cp -f "$C_RAW" "${ITERDIR}/candidates.txt"
	A_CUR="$A_NEXT"
	iter=$((iter + 1))
done

# final
oatk_log 1 "[done] A_final=${A_CUR}"
oatk_log 1 "[log]  ${ILOG}"
exit 0
