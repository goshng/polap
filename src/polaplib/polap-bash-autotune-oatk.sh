#!/usr/bin/env bash
set -euo pipefail

# polap-bash-autotune-oatk.sh
# YAML-config + CLI overrides (CLI wins). Orchestrates:
#   (optional) WGS recruitment via select-mt  ->
#   tune-c (coverage valley)                  ->
#   seed sweep (k/s, -c grid) + scoring       ->
#   optional final pass if --final-if holds   ->
#   summary CSV and provenance.
#
# Addressing (one of):
#   --config-path FILE.yaml
#   --config-dir DIR --preset NAME         # resolves DIR/NAME.yaml (default DIR=~/.polap/profiles)
#
# Any flat YAML key can be overridden on CLI: replace '_' with '-'.
# Canonical output locator is out_prefix (or --out-prefix); outdir = dirname(out_prefix).
#
# Booleans (toggles without values):
#   --wgs-mode/--no-wgs-mode, --hpc/--no-hpc, --hpc-qv/--no-hpc-qv, --ppr-hpc/--no-ppr-hpc,
#   --edge-norm/--no-edge-norm, --emit-fastq/--no-emit-fastq, --keep-intermediate/--no-keep-intermediate,
#   --redo/--no-redo, --plot/--no-plot
#
# Anchors (label-aware; never use --anchors):
#   --label mt|pt
#   --mt-anchors PATH   (required if --label mt and wgs_mode=true)
#   --pt-anchors PATH   (required if --label pt and wgs_mode=true)

# ---------- paths ----------
POLAPLIB_DIR="${_POLAPLIB_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)}"
LOAD_PY="${POLAPLIB_DIR}/polap-py-config-load.py"
SELECT_MT_SH="${POLAPLIB_DIR}/polap-bash-select-mt.sh"
TUNE_C_SH="${POLAPLIB_DIR}/polap-bash-oatk-tune-c.sh"
SCORE_SH="${POLAPLIB_DIR}/polap-bash-oatk-score-seeds.sh"
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

# ---------- helpers ----------
_norm_bool() {
	local v="${1:-}"
	case "${v,,}" in
	true | 1 | yes | y | on) echo "true" ;;
	false | 0 | no | n | off | "") echo "false" ;;
	*) echo "$v" ;;
	esac
}

# ---------- parse addressing + verbosity + collect overrides ----------
CFG_PATH=""
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

oatk_log_set_level "$LOG_COUNT"
oatk_enable_trace_if_needed

# ---------- deps ----------
need() { command -v "$1" >/dev/null 2>&1 || {
	oatk_log 1 "$1 missing"
	exit 127
}; }
need python3
need bash
need syncasm
need minimap2
need samtools
need seqtk
[[ -f "$LOAD_PY" && -f "$SELECT_MT_SH" && -f "$TUNE_C_SH" && -f "$SCORE_SH" ]] || {
	oatk_log 1 "ensure helper scripts exist in ${POLAPLIB_DIR}"
	exit 127
}

# default config dir
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
apply_overrides() {
	local j=0
	while [[ $j -lt ${#ARGS[@]} ]]; do
		local tok="${ARGS[$j]}"
		case "$tok" in
		# boolean toggles (no value)
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
		--plot)
			PCFG_PLOT="true"
			j=$((j + 1))
			continue
			;;
		--no-plot)
			PCFG_PLOT="false"
			j=$((j + 1))
			continue
			;;

		# K/V: --foo-bar -> PCFG_FOO_BAR
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

# anchors (label-aware) in WGS mode
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

# compute defaults
PCFG_THREADS="${PCFG_THREADS:-16}"
PCFG_HPC="${PCFG_HPC:-false}"
PCFG_K="${PCFG_K:-$([[ "$PCFG_HPC" == "true" ]] && echo 41 || echo 121)}"
PCFG_S="${PCFG_S:-$([[ "$PCFG_HPC" == "true" ]] && echo 21 || echo 27)}"

case "${PCFG_PRESET,,}" in
hifi | ont) TUNE_PRESET="${PCFG_PRESET,,}" ;;
*) TUNE_PRESET=$([[ "$PCFG_HPC" == "true" ]] && echo "ont" || echo "hifi") ;;
esac

# optionals
PCFG_FINAL_IF="${PCFG_FINAL_IF:-}"
PCFG_C_RADIUS="${PCFG_C_RADIUS:-0,10,20}"
PCFG_MIN_SHARED="${PCFG_MIN_SHARED:-5}"
PCFG_JACCARD_MIN="${PCFG_JACCARD_MIN:-$([[ "$TUNE_PRESET" == "ont" ]] && echo 0.0075 || echo 0.015)}"
PCFG_TOPK_NEI="${PCFG_TOPK_NEI:-$([[ "$TUNE_PRESET" == "ont" ]] && echo 40 || echo 50)}"
PCFG_STEPS="${PCFG_STEPS:-$([[ "$TUNE_PRESET" == "ont" ]] && echo 1 || echo 2)}"
PCFG_MAX_RUNS="${PCFG_MAX_RUNS:-12}"
PCFG_PLOT="$(_norm_bool "${PCFG_PLOT:-false}")"

# ---------- summary + provenance ----------
oatk_log 1 "config=${RESOLVED_CFG}"
oatk_log 1 "preset=${PCFG_PRESET:-} platform=${TUNE_PRESET} label=${LABEL}"
oatk_log 1 "wgs_mode=${PCFG_WGS_MODE} reads=${PCFG_READS}"
[[ "$PCFG_WGS_MODE" == "true" && "$LABEL" == "mt" ]] && oatk_log 1 "mt_anchors=${PCFG_MT_ANCHORS}"
[[ "$PCFG_WGS_MODE" == "true" && "$LABEL" == "pt" ]] && oatk_log 1 "pt_anchors=${PCFG_PT_ANCHORS}"
oatk_log 1 "hmm_db=${PCFG_HMM_DB} out_prefix=${PCFG_OUT_PREFIX} outdir=${OUTDIR}"
oatk_log 1 "threads=${PCFG_THREADS} k=${PCFG_K} s=${PCFG_S} hpc=${PCFG_HPC}"
[[ -n "${PCFG_FINAL_IF}" ]] && oatk_log 1 "final_if=${PCFG_FINAL_IF}"
oatk_log 2 "c_radius=${PCFG_C_RADIUS} min_shared=${PCFG_MIN_SHARED} jaccard_min=${PCFG_JACCARD_MIN} topk_nei=${PCFG_TOPK_NEI} steps=${PCFG_STEPS} max_runs=${PCFG_MAX_RUNS}"

filtered_overrides=$(oatk_filter_overrides "${ARGS[@]}")
oatk_provenance_line "$RESOLVED_CFG" $filtered_overrides
oatk_provenance_append "${OUTDIR}/provenance.autotune.txt" "$RESOLVED_CFG" $filtered_overrides

# Build verbosity flags (-v repeated LOG_COUNT times)
VERBOSE_FLAGS=()
if [[ -n "${LOG_COUNT:-}" && "$LOG_COUNT" -gt 0 ]]; then
	for ((i = 0; i < LOG_COUNT; i++)); do
		VERBOSE_FLAGS+=(-v)
	done
fi

# ---------- Step 0: WGS recruitment (if needed) ----------
READS_IN="${PCFG_READS}"
if [[ "$PCFG_WGS_MODE" == "true" ]]; then
	RECR_PREFIX="${OUTDIR}/${LABEL}.wgs"
	oatk_log 1 "[wgs] recruiting $LABEL via select-mt"
	if [[ "$LABEL" == "pt" ]]; then
		oatk_run bash "$SELECT_MT_SH" select-mt \
			--config-dir "${CONFIG_DIR:-}" --preset "${PCFG_PRESET:-$TUNE_PRESET}" \
			--label pt \
			--pt-anchors "${PCFG_PT_ANCHORS}" \
			--mt-anchors "${PCFG_MT_ANCHORS}" \
			--reads "$PCFG_READS" \
			--out-prefix "$RECR_PREFIX" \
			--threads "$PCFG_THREADS" \
			--ensemble majority \
			"${VERBOSE_FLAGS[@]}" \
			--emit-fastq || true
		if [[ -s "${RECR_PREFIX}.pt.fastq.gz" ]]; then
			READS_IN="${RECR_PREFIX}.pt.fastq"
			gzip -dc "${RECR_PREFIX}.pt.fastq.gz" >"$READS_IN"
		elif [[ -s "${RECR_PREFIX}.ensemble.pt.ids" ]]; then
			awk 'NF{print $1}' "${RECR_PREFIX}.ensemble.pt.ids" | sort -u >"${RECR_PREFIX}.pt.names"
			seqtk subseq "$PCFG_READS" "${RECR_PREFIX}.pt.names" >"${RECR_PREFIX}.pt.fastq"
			READS_IN="${RECR_PREFIX}.pt.fastq"
		else
			oatk_log 1 "WGS recruitment failed for pt"
			exit 3
		fi
	else

		# oatk_run bash "$SELECT_MT_SH" select-mt \
		# 	--config-dir "${CONFIG_DIR:-}" \
		# 	--preset "${PCFG_PRESET:-$TUNE_PRESET}" \
		# 	--label mt \
		# 	--pt-anchors "${PCFG_PT_ANCHORS}" \
		# 	--mt-anchors "${PCFG_MT_ANCHORS}" \
		# 	--reads "$PCFG_READS" \
		# 	--out-prefix "$RECR_PREFIX" \
		# 	--threads "$PCFG_THREADS" \
		# 	--ensemble majority \
		# 	"${VERBOSE_FLAGS[@]}" \
		# 	--emit-fastq || true

		# Build verbosity flags, if you use -v cascaded:
		CMD=(bash "$SELECT_MT_SH" select-mt
			--config-dir "${CONFIG_DIR:-}"
			--preset "${PCFG_PRESET:-$TUNE_PRESET}"
			--label mt
			--pt-anchors "${PCFG_PT_ANCHORS}"
			--mt-anchors "${PCFG_MT_ANCHORS}"
			--reads "$PCFG_READS"
			--out-prefix "$RECR_PREFIX"
			--threads "$PCFG_THREADS"
			--ensemble majority
			--emit-fastq
			"${VERBOSE_FLAGS[@]}"
		)

		# Log + run (no truncation; arguments preserved)
		oatk_run "${CMD[@]}" || true

		if [[ -s "${RECR_PREFIX}.mt.fastq.gz" ]]; then
			READS_IN="${RECR_PREFIX}.mt.fastq"
			gzip -dc "${RECR_PREFIX}.mt.fastq.gz" >"$READS_IN"
		elif [[ -s "${RECR_PREFIX}.ensemble.mt.ids" ]]; then
			awk 'NF{print $1}' "${RECR_PREFIX}.ensemble.mt.ids" | sort -u >"${RECR_PREFIX}.mt.names"
			seqtk subseq "$PCFG_READS" "${RECR_PREFIX}.mt.names" >"${RECR_PREFIX}.mt.fastq"
			READS_IN="${RECR_PREFIX}.mt.fastq"
		else
			oatk_log 1 "WGS recruitment failed for mt"
			exit 3
		fi
	fi
	oatk_log 1 "[wgs] recruited -> $READS_IN"
fi

# ---------- Step 1: Tune -c ----------
TUNE_PREFIX="${OUTDIR}/${LABEL}.tuneC"
TUNE_CMD=(bash "$TUNE_C_SH"
	--config-dir "${CONFIG_DIR:-}" --preset "${PCFG_PRESET:-$TUNE_PRESET}"
	--reads "$READS_IN"
	--out-prefix "$TUNE_PREFIX"
	--threads "$PCFG_THREADS"
	--k "$PCFG_K" --s "$PCFG_S"
)
[[ "$PCFG_PLOT" == "true" ]] && TUNE_CMD+=(--plot)
oatk_run "${TUNE_CMD[@]}" || true

SUG_LINE=$(tail -n 1 "${TUNE_PREFIX}.suggest.txt" 2>/dev/null || echo "")
SUG_C=$(echo "$SUG_LINE" | awk '{for(i=1;i<=NF;i++) if($i=="-c") print $(i+1)}' 2>/dev/null)
[[ -z "$SUG_C" ]] && SUG_C=50
oatk_log 1 "[tune] suggested -c ≈ $SUG_C"

# ---------- Step 2: Build c-grid ----------
declare -a CAND_C=() GRID=()
declare -A SEEN
IFS=',' read -r -a TOKS <<<"${PCFG_C_RADIUS}"
for tok in "${TOKS[@]}"; do
	tok="$(echo "$tok" | tr -d ' ')"
	if [[ "$tok" =~ %$ ]]; then
		pct="${tok%\%}"
		pct="${pct//+/}"
		pct="${pct//-/}"
		if [[ "$tok" == -*% ]]; then
			off=$(awk -v c="$SUG_C" -v p="$pct" 'BEGIN{printf "%.0f", c*(1.0 - p/100.0)}')
		else
			off=$(awk -v c="$SUG_C" -v p="$pct" 'BEGIN{printf "%.0f", c*(1.0 + p/100.0)}')
		fi
		[[ -n "${off:-}" ]] && CAND_C+=("$off")
	else
		off=$((SUG_C + ${tok#+}))
		CAND_C+=("$off")
	fi
done
for c in "${CAND_C[@]}"; do
	[[ "$c" -lt 1 ]] && c=1
	if [[ -z "${SEEN[$c]:-}" ]]; then
		SEEN[$c]=1
		GRID+=("$c")
	fi
done
[[ "${#GRID[@]}" -eq 0 ]] && GRID=("$SUG_C")
oatk_log 1 "[grid] c-values: ${GRID[*]}"

# ---------- Step 3: Sweep seeds & score ----------
SUMMARY="${OUTDIR}/${LABEL}.summary.csv"
echo "phase,k,s,label,c,opts,out_prefix,score_csv,score" >"$SUMMARY"

runs=0
for C in "${GRID[@]}"; do
	RUN_SEED="${OUTDIR}/${LABEL}.k${PCFG_K}.c${C}.seed"
	oatk_log 1 "[seed] assembling k=${PCFG_K} c=$C"
	oatk_run syncasm -k "$PCFG_K" -s "$PCFG_S" $([[ "$TUNE_PRESET" == "ont" ]] && echo --hpc || :) \
		-t "$PCFG_THREADS" -o "${RUN_SEED}.asm" -c "$C" -a 0.05 \
		--unzip-round 0 --max-bubble 0 --max-tip 0 --weak-cross 0 --no-read-ec \
		"$READS_IN" || true

	GFA=""
	for g in "${RUN_SEED}.asm.utg.final.gfa" "${RUN_SEED}.asm.utg.gfa" "${RUN_SEED}.asm.utg.clean.gfa"; do
		[[ -s "$g" ]] && {
			GFA="$g"
			break
		}
	done

	if [[ -z "$GFA" ]]; then
		echo "seed,${PCFG_K},${PCFG_S},${LABEL},${C},\"seed\",${RUN_SEED},NA,0.000000" >>"$SUMMARY"
	else
		SEEDS_FA="${RUN_SEED}.seeds.fa"
		awk -F'\t' '$1=="S"{print ">"$2"\n"$3}' "$GFA" >"$SEEDS_FA"

		SCORE_CMD=(bash "$SCORE_SH"
			--config-dir "${CONFIG_DIR:-}" --preset "${PCFG_PRESET:-$TUNE_PRESET}"
			--seeds "$SEEDS_FA"
			--reads "$READS_IN"
			--hmm-db "$PCFG_HMM_DB"
			--out-prefix "${RUN_SEED}"
			--threads "$PCFG_THREADS"
		)
		[[ -n "${PCFG_NUC_CUT_LOG10:-}" ]] && SCORE_CMD+=(--nuc-cut-log10 "$PCFG_NUC_CUT_LOG10")
		[[ -n "${PCFG_X_SLOPE:-}" ]] && SCORE_CMD+=(--x-slope "$PCFG_X_SLOPE")
		[[ -n "${PCFG_MIN_SHARED:-}" ]] && SCORE_CMD+=(--min-shared "$PCFG_MIN_SHARED")
		[[ -n "${PCFG_JACCARD_MIN:-}" ]] && SCORE_CMD+=(--jaccard-min "$PCFG_JACCARD_MIN")
		[[ -n "${PCFG_TOPK_NEI:-}" ]] && SCORE_CMD+=(--topk-nei "$PCFG_TOPK_NEI")
		[[ -n "${PCFG_STEPS:-}" ]] && SCORE_CMD+=(--steps "$PCFG_STEPS")
		oatk_run "${SCORE_CMD[@]}"

		SCV="${RUN_SEED}.seed_score.csv"
		SCORE=$(tail -n 1 "$SCV" 2>/dev/null | awk -F',' '{print $NF}')
		[[ -z "$SCORE" ]] && SCORE="0.000000"
		echo "seed,${PCFG_K},${PCFG_S},${LABEL},${C},\"seed\",${RUN_SEED},${SCV},${SCORE}" >>"$SUMMARY"
	fi

	runs=$((runs + 1))
	[[ $runs -ge ${PCFG_MAX_RUNS} ]] && break
done

# ---------- Step 4: Final decision ----------
BEST_SEED=$(awk -F',' 'NR>1 && $1=="seed"{print $0}' "$SUMMARY" | sort -t',' -k9,9gr | head -n1 || true)
[[ -n "$BEST_SEED" ]] || {
	oatk_log 1 "no successful seed builds"
	exit 4
}

best_c=$(echo "$BEST_SEED" | awk -F',' '{print $5}')
seed_pref=$(echo "$BEST_SEED" | awk -F',' '{print $7}')
oatk_log 1 "[best-seed] c=$best_c out=$seed_pref"

RUN_FINAL=1
if [[ -n "${PCFG_FINAL_IF:-}" ]]; then
	seed_csv="${seed_pref}.seed_score.csv"
	if [[ -s "$seed_csv" ]]; then
		read -r genes_score breadth depth n50 len_total x_score score < <(awk -F',' 'NR==1{next} END{printf "%s %s %s %s %s %s %s", $6,$7,$8,$9,$10,$12,$13}' "$seed_csv")
		expr="${PCFG_FINAL_IF}"
		expr=${expr//genes_score/$genes_score}
		expr=${expr//breadth/$breadth}
		expr=${expr//depth/$depth}
		expr=${expr//n50/$n50}
		expr=${expr//len_total/$len_total}
		expr=${expr//x_score/$x_score}
		expr=${expr//score/$score}
		ok=$(echo "$expr" | bc -l 2>/dev/null || echo 1)
		[[ "$ok" == "1" ]] && RUN_FINAL=1 || RUN_FINAL=0
		oatk_log 1 "[final-if] '$PCFG_FINAL_IF' => $ok (1=run,0=skip)"
	fi
fi

if [[ $RUN_FINAL -eq 1 ]]; then
	RUN_FIN="${seed_pref}.final"
	oatk_log 1 "[final] assembling k=${PCFG_K} c=$best_c (relaxed)…"
	oatk_run syncasm -k "$PCFG_K" -s "$PCFG_S" $([[ "$TUNE_PRESET" == "ont" ]] && echo --hpc || :) \
		-t "$PCFG_THREADS" -o "${RUN_FIN}.asm" -c "$best_c" -a 0.05 \
		--unzip-round 1 --max-bubble 1 --max-tip 1 --weak-cross 1 --no-read-ec \
		"$READS_IN" || true

	FGFA=""
	for g in "${RUN_FIN}.asm.utg.final.gfa" "${RUN_FIN}.asm.utg.gfa" "${RUN_FIN}.asm.utg.clean.gfa"; do
		[[ -s "$g" ]] && {
			FGFA="$g"
			break
		}
	done
	if [[ -n "$FGFA" ]]; then
		FFA="${RUN_FIN}.seeds.fa"
		awk -F'\t' '$1=="S"{print ">"$2"\n"$3}' "$FGFA" >"$FFA"

		SCOREF_CMD=(bash "$SCORE_SH"
			--config-dir "${CONFIG_DIR:-}" --preset "${PCFG_PRESET:-$TUNE_PRESET}"
			--seeds "$FFA" --reads "$READS_IN" --hmm-db "$PCFG_HMM_DB"
			--out-prefix "${RUN_FIN}" --threads "$PCFG_THREADS"
		)
		[[ -n "${PCFG_NUC_CUT_LOG10:-}" ]] && SCOREF_CMD+=(--nuc-cut-log10 "$PCFG_NUC_CUT_LOG10")
		[[ -n "${PCFG_X_SLOPE:-}" ]] && SCOREF_CMD+=(--x-slope "$PCFG_X_SLOPE")
		[[ -n "${PCFG_MIN_SHARED:-}" ]] && SCOREF_CMD+=(--min-shared "$PCFG_MIN_SHARED")
		[[ -n "${PCFG_JACCARD_MIN:-}" ]] && SCOREF_CMD+=(--jaccard-min "$PCFG_JACCARD_MIN")
		[[ -n "${PCFG_TOPK_NEI:-}" ]] && SCOREF_CMD+=(--topk-nei "$PCFG_TOPK_NEI")
		[[ -n "${PCFG_STEPS:-}" ]] && SCOREF_CMD+=(--steps "$PCFG_STEPS")
		oatk_run "${SCOREF_CMD[@]}"

		FSC="${RUN_FIN}.seed_score.csv"
		FSCORE=$(tail -n 1 "$FSC" 2>/dev/null | awk -F',' '{print $NF}')
		[[ -z "$FSCORE" ]] && FSCORE="0.000000"
		echo "final,${PCFG_K},${PCFG_S},${LABEL},${best_c},\"final\",${RUN_FIN},${FSC},${FSCORE}" >>"${OUTDIR}/${LABEL}.summary.csv"
	else
		echo "final,${PCFG_K},${PCFG_S},${LABEL},${best_c},\"final\",${RUN_FIN},NA,0.000000" >>"${OUTDIR}/${LABEL}.summary.csv"
	fi
else
	echo "final,${PCFG_K},${PCFG_S},${LABEL},${best_c},\"skipped\",${seed_pref}.final,NA,0.000000" >>"${OUTDIR}/${LABEL}.summary.csv"
fi

oatk_log 1 "[done] summary: ${OUTDIR}/${LABEL}.summary.csv"
exit 0
