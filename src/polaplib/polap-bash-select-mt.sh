#!/usr/bin/env bash
set -euo pipefail

# polap-bash-select-mt.sh
# YAML-config + CLI overrides (CLI wins). Select mt/pt reads via:
#   syncfilter quickview (coverage X) → cut-auto (bands) → PPR (syncmer connectivity) → ensemble → optional FASTQs
#
# Subcommands:
#   select-mt         : run selection pipeline
#   select-mt-report  : render HTML report (after you generate GFAs with Flye)
#
# Addressing (one of):
#   --config-path FILE.yaml
#   --config-dir DIR --preset NAME        # resolves DIR/NAME.yaml (default DIR=~/.polap/profiles)
#
# Any YAML key can be overridden by CLI (replace '_' → '-'), e.g.:
#   --reads, --out-prefix, --threads, --k-qv, --s-qv, --hpc-qv, --x-method, --x-tail, --tail-mt, --tail-pt,
#   --x-win, --band-max, --sigma-floor, --ppr-k, --ppr-s, --ppr-hpc, --max-occ, --min-shared, --jaccard-min,
#   --edge-norm, --topk-nei, --steps, --ppr-alpha, --ppr-iter, --score-th, --nuc-cut-log10, --x-slope, --ensemble,
#   --emit-fastq, --keep-intermediate, --redo, --sf-verbose, --label, --mt-anchors, --pt-anchors
#
# Verbosity:
#   -v / --verbose  (repeatable; -vv enables bash xtrace with file:function:line PS4)
#
# Requires:
#   syncfilter (>= v0.7.0 with --emit-read-syncmers), seqtk, python3
#   ${_POLAPLIB_DIR}/polap-bash-oatk-log.sh
#   ${_POLAPLIB_DIR}/polap-py-config-load.py
#   ${_POLAPLIB_DIR}/polap-py-syncfilter-cut-auto.py
#   ${_POLAPLIB_DIR}/polap-py-syncmer-connectivity-select-mt.py
#   ${_POLAPLIB_DIR}/polap-py-ensemble-mt.py
#   ${_POLAPLIB_DIR}/polap-py-ensemble-report.py  (select-mt-report)

# ---------- paths ----------
POLAPLIB_DIR="${_POLAPLIB_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)}"
LOG_LIB="${POLAPLIB_DIR}/polap-bash-oatk-log.sh"
LOAD_PY="${POLAPLIB_DIR}/polap-py-config-load.py"
CUT_AUTO="${POLAPLIB_DIR}/polap-py-syncfilter-cut-auto.py"
PPR_SEL="${POLAPLIB_DIR}/polap-py-syncmer-connectivity-select-mt.py"
ENSEMBLE_PY="${POLAPLIB_DIR}/polap-py-ensemble-mt.py"
REPORT_PY="${POLAPLIB_DIR}/polap-py-ensemble-report.py"

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
need() { command -v "$1" >/dev/null 2>&1 || die "$1 missing"; }

_norm_bool() {
	local v="${1:-}"
	case "${v,,}" in true | 1 | yes | y | on) echo "true" ;;
	false | 0 | no | n | off | "") echo "false" ;;
	*) echo "$v" ;;
	esac
}

_usage() {
	cat <<'EOF'
Usage:
  polap-bash-select-mt.sh select-mt [--config-path FILE.yaml | --config-dir DIR --preset NAME] [overrides...]
  polap-bash-select-mt.sh select-mt-report [--config-path FILE.yaml | --config-dir DIR --preset NAME] [--out PREFIX] [--bandage-size WxH]

Notes:
  • Config is flat YAML; any key can be overridden on the CLI by replacing '_' with '-'.
  • Canonical output locus is 'out_prefix' (YAML or --out-prefix). Script derives dirs from it.
EOF
}

# ---------------- YAML loader + CLI overrides ----------------
_load_config_env() {
	local CFG_PATH="" CONFIG_DIR="" PRESET="" RESOLVED_CFG=""
	local LOG_COUNT=0
	local ARGV=("$@") i=0

	# capture verbosity flags for this stage
	while [[ $i -lt ${#ARGV[@]} ]]; do
		case "${ARGV[$i]}" in
		-v | --verbose)
			LOG_COUNT=$((LOG_COUNT + 1))
			i=$((i + 1))
			;;
		--config-path)
			CFG_PATH="${ARGV[$((i + 1))]:-}"
			i=$((i + 2))
			;;
		--config-dir)
			CONFIG_DIR="${ARGV[$((i + 1))]:-}"
			i=$((i + 2))
			;;
		--preset)
			PRESET="${ARGV[$((i + 1))]:-}"
			i=$((i + 2))
			;;
		*) i=$((i + 1)) ;;
		esac
	done
	oatk_log_set_level "$LOG_COUNT"
	oatk_enable_trace_if_needed

	[[ -z "$CFG_PATH" && -z "$CONFIG_DIR" ]] && CONFIG_DIR="${HOME}/.polap/profiles"

	need python3
	[[ -f "$LOAD_PY" ]] || die "not found: $LOAD_PY"

	if [[ -n "$CFG_PATH" ]]; then
		oatk_log 1 "loading YAML --config-path ${CFG_PATH}"
		eval "$(python3 "$LOAD_PY" --path "$CFG_PATH" --format env --prefix SELCFG 2>/dev/null)" ||
			die "failed to load config: $CFG_PATH"
		RESOLVED_CFG="$CFG_PATH"
	else
		[[ -n "$PRESET" ]] || die "need --preset when --config-path is not given"
		oatk_log 1 "loading YAML --config-dir ${CONFIG_DIR} --preset ${PRESET}"
		eval "$(python3 "$LOAD_PY" --config-dir "$CONFIG_DIR" --preset "$PRESET" --format env --prefix SELCFG 2>/dev/null)" ||
			die "failed to load config: ${CONFIG_DIR}/${PRESET}.yaml"
		RESOLVED_CFG="${CONFIG_DIR}/${PRESET}.yaml"
	fi
	SELCFG_RESOLVED_CFG="$RESOLVED_CFG"
}

_apply_overrides() {
	local ARGS=("$@") j=0
	while [[ $j -lt ${#ARGS[@]} ]]; do
		local tok="${ARGS[$j]}"
		case "$tok" in
		-v | --verbose)
			j=$((j + 1))
			continue
			;;

		# addressing (already handled)
		--config-path | --config-dir | --preset)
			j=$((j + 2))
			continue
			;;

		# boolean toggles (no value)
		--wgs-mode)
			SELCFG_WGS_MODE="true"
			j=$((j + 1))
			continue
			;;
		--no-wgs-mode)
			SELCFG_WGS_MODE="false"
			j=$((j + 1))
			continue
			;;
		--hpc-qv)
			SELCFG_HPC_QV="true"
			j=$((j + 1))
			continue
			;;
		--no-hpc-qv)
			SELCFG_HPC_QV="false"
			j=$((j + 1))
			continue
			;;
		--ppr-hpc)
			SELCFG_PPR_HPC="true"
			j=$((j + 1))
			continue
			;;
		--no-ppr-hpc)
			SELCFG_PPR_HPC="false"
			j=$((j + 1))
			continue
			;;
		--edge-norm)
			SELCFG_EDGE_NORM="true"
			j=$((j + 1))
			continue
			;;
		--no-edge-norm)
			SELCFG_EDGE_NORM="false"
			j=$((j + 1))
			continue
			;;
		--emit-fastq)
			SELCFG_EMIT_FASTQ="true"
			j=$((j + 1))
			continue
			;;
		--no-emit-fastq)
			SELCFG_EMIT_FASTQ="false"
			j=$((j + 1))
			continue
			;;
		--keep-intermediate)
			SELCFG_KEEP_INTERMEDIATE="true"
			j=$((j + 1))
			continue
			;;
		--no-keep-intermediate)
			SELCFG_KEEP_INTERMEDIATE="false"
			j=$((j + 1))
			continue
			;;
		--redo)
			SELCFG_REDO="true"
			j=$((j + 1))
			continue
			;;
		--no-redo)
			SELCFG_REDO="false"
			j=$((j + 1))
			continue
			;;

		# K/V flags (generic)
		--*)
			local key="${tok#--}"
			key="${key//-/_}"
			local val="${ARGS[$((j + 1))]:-}"
			[[ -z "$val" || "$val" == --* ]] && die "missing value for $tok"
			case "$tok" in
			--mt-anchors) SELCFG_MT_ANCHORS="$val" ;;
			--pt-anchors) SELCFG_PT_ANCHORS="$val" ;;
			--out-prefix) SELCFG_OUT_PREFIX="$val" ;;
			--label) SELCFG_LABEL="$val" ;;
			*)
				local varname="SELCFG_${key^^}"
				printf -v "$varname" "%s" "$val"
				;;
			esac
			j=$((j + 2))
			;;
		*) j=$((j + 1)) ;;
		esac
	done
}

# ---------------- core: select-mt ----------------
_select_mt() {
	_load_config_env "$@"
	_apply_overrides "$@"

	need syncfilter
	need seqtk
	[[ -f "$CUT_AUTO" && -f "$PPR_SEL" && -f "$ENSEMBLE_PY" ]] || die "missing python tools in $POLAPLIB_DIR"

	local OUT="${SELCFG_OUT_PREFIX:-}"
	[[ -n "$OUT" ]] || die "out_prefix missing (YAML key: out_prefix or --out-prefix)"
	local OUTDIR
	OUTDIR="$(dirname -- "$OUT")"
	mkdir -p "$OUTDIR" "${OUTDIR}/tmp"

	# provenance
	local filtered_overrides
	filtered_overrides=$(oatk_filter_overrides "$@")
	oatk_provenance_line "${SELCFG_RESOLVED_CFG:-?}" $filtered_overrides
	oatk_provenance_append "${OUTDIR}/provenance.select-mt.txt" "${SELCFG_RESOLVED_CFG:-?}" $filtered_overrides

	# Core inputs
	local READS="${SELCFG_READS:-}"
	[[ -n "$READS" && -s "$READS" ]] || die "reads missing/empty"
	local THREADS="${SELCFG_THREADS:-8}"
	local PRESET="${SELCFG_PRESET:-hifi}"
	local LABEL="${SELCFG_LABEL:-mt}"
	local REDO="$(_norm_bool "${SELCFG_REDO:-false}")"
	[[ "$REDO" == "true" && -d "$OUTDIR" ]] && rm -rf "$OUTDIR"
	mkdir -p "$OUTDIR" "${OUTDIR}/tmp"

	# Quickview params
	local K_QV="${SELCFG_K_QV:-$([[ "$PRESET" == "ont" ]] && echo 41 || echo 121)}"
	local S_QV="${SELCFG_S_QV:-$([[ "$PRESET" == "ont" ]] && echo 17 || echo 27)}"
	local _hpc_qv_default=$([[ "$PRESET" == "ont" ]] && echo true || echo false)
	local HPC_QV="$(_norm_bool "${SELCFG_HPC_QV:-_hpc_qv_default}")"
	[[ "$HPC_QV" == "_hpc_qv_default" ]] && HPC_QV="$(_norm_bool "${_hpc_qv_default}")"
	local HPC_QV_FLAG=""
	[[ "$HPC_QV" == "true" ]] && HPC_QV_FLAG="--hpc"
	local SFV="${SELCFG_SF_VERBOSE:-0}"

	# Coverage banding
	local XMETHOD="${SELCFG_X_METHOD:-hybrid}"
	local XTAIL="${SELCFG_X_TAIL:-0.05}"
	local TAIL_MT="${SELCFG_TAIL_MT:-}"
	local TAIL_PT="${SELCFG_TAIL_PT:-}"
	local XWIN="${SELCFG_X_WIN:-0.20}"
	local BAND_MAX="${SELCFG_BAND_MAX:-2.0}"
	local SIGMA_FLOOR="${SELCFG_SIGMA_FLOOR:-0.03}"

	# Anchors
	local MT_ANCH="${SELCFG_MT_ANCHORS:-}"
	local PT_ANCH="${SELCFG_PT_ANCHORS:-}"

	# PPR params
	local PPR_K="${SELCFG_PPR_K:-$K_QV}"
	local PPR_S="${SELCFG_PPR_S:-$S_QV}"
	local _ppr_hpc_default=$([[ "$PRESET" == "ont" ]] && echo true || echo false)
	local PPR_HPC="$(_norm_bool "${SELCFG_PPR_HPC:-$_ppr_hpc_default}")"
	local PPR_HPC_FLAG=""
	[[ "$PPR_HPC" == "true" ]] && PPR_HPC_FLAG="--hpc"
	local MAX_OCC="${SELCFG_MAX_OCC:-200}"
	local MIN_SHARED="${SELCFG_MIN_SHARED:-$([[ "$PRESET" == "ont" ]] && echo 5 || echo 4)}"
	local JAC_MIN="${SELCFG_JACCARD_MIN:-$([[ "$PRESET" == "ont" ]] && echo 0.0075 || echo 0.015)}"
	local EDGE_NORM="$(_norm_bool "${SELCFG_EDGE_NORM:-true}")"
	local EDGE_NORM_FLAG=""
	[[ "$EDGE_NORM" == "true" ]] && EDGE_NORM_FLAG="--edge-norm"
	local TOPK="${SELCFG_TOPK_NEI:-$([[ "$PRESET" == "ont" ]] && echo 40 || echo 50)}"
	local STEPS="${SELCFG_STEPS:-$([[ "$PRESET" == "ont" ]] && echo 1 || echo 2)}"
	local PPR_ALPHA="${SELCFG_PPR_ALPHA:-0.85}"
	local PPR_ITER="${SELCFG_PPR_ITER:-30}"
	local SCORE_TH="${SELCFG_SCORE_TH:-0.0}"

	# X-prior
	local NUC_CUT_LOG10="${SELCFG_NUC_CUT_LOG10:-}"
	local X_SLOPE="${SELCFG_X_SLOPE:-}"

	# Ensemble / outputs
	local ENSEMBLE="${SELCFG_ENSEMBLE:-majority}"
	local EMIT_FASTQ="$(_norm_bool "${SELCFG_EMIT_FASTQ:-true}")"
	local KEEP_INT="$(_norm_bool "${SELCFG_KEEP_INTERMEDIATE:-false}")"

	# 1) Quickview
	oatk_log 1 "[step] quickview (coverage X)"
	oatk_run syncfilter --mode quickview $HPC_QV_FLAG -k "$K_QV" -s "$S_QV" -t "$THREADS" -v "$SFV" -o "$OUT" "$READS"
	[[ -f "${OUT}.syncfilter.tsv" ]] && mv -f "${OUT}.syncfilter.tsv" "${OUT}.quickview.tsv"

	# 2) Coverage banding (anchors → bands)
	oatk_log 1 "[step] coverage banding"
	oatk_run python3 "$CUT_AUTO" \
		-i "${OUT}.quickview.tsv" \
		--method "$XMETHOD" --tail "$XTAIL" \
		${TAIL_MT:+--taili-mt "$TAIL_MT"} \
		${TAIL_PT:+--tail-pt "$TAIL_PT"} \
		--window-width "$XWIN" --band-max "$BAND_MAX" --sigma-floor "$SIGMA_FLOOR" \
		${MT_ANCH:+--mt-anchors "$MT_ANCH"} ${PT_ANCH:+--pt-anchors "$PT_ANCH"} \
		--cuts-out "${OUT}.xcover.cuts.tsv" \
		--png "${OUT}.xcover.x_hist.png" \
		--tsv-out "${OUT}.xcover.tsv"

	# helper: run PPR + ensemble for a specific label (mt or pt)
	_ppr_for_label() {
		local L="$1"
		local IDS="${OUT}.xcover.${L}.ids"
		local NAMES="${OUT}.xcover.${L}.names"
		local FASTQ="${OUT}.xcover.${L}.fastq"
		local ANCH="${MT_ANCH}"
		[[ "$L" == "pt" ]] && ANCH="${PT_ANCH}"
		[[ -s "$IDS" ]] || {
			oatk_log 1 "[warn] no ${L}-band IDs found"
			return 1
		}

		oatk_log 1 "[step] subset ${L}-band reads"
		awk 'NF{print $1}' "$IDS" | sort -u >"$NAMES"
		oatk_run seqtk subseq "$READS" "$NAMES" >"$FASTQ"

		oatk_log 1 "[step] dump read → syncmers (${L})"
		oatk_run syncfilter --mode quickview $PPR_HPC_FLAG -k "$PPR_K" -s "$PPR_S" -t "$THREADS" -v "$SFV" \
			--emit-read-syncmers "${OUT}.xcover.${L}" \
			-o "${OUT}.xcover.${L}.pprprep" \
			"$FASTQ"

		local SM_DUMP="${OUT}.xcover.${L}.read2sm.tsv"
		[[ -s "$SM_DUMP" ]] || die "missing $SM_DUMP"

		oatk_log 1 "[step] PPR (${L}) with sm-dump"
		local PPR_CMD=(python3 "$PPR_SEL"
			--reads "$FASTQ"
			--sm-dump "$SM_DUMP"
			--mt "$ANCH"
			--label "$L"
			--max-occ "$MAX_OCC" --min-shared "$MIN_SHARED" --jaccard-min "$JAC_MIN"
			--topk-nei "$TOPK" --steps "$STEPS"
			--ppr-alpha "$PPR_ALPHA" --ppr-iter "$PPR_ITER"
			--tmpdir "${OUTDIR}/tmp"
			--score-th "$SCORE_TH"
			$([[ "$EDGE_NORM" == "true" ]] && echo --edge-norm)
			-o "${OUT}.ppr"
		)
		if [[ -n "$NUC_CUT_LOG10" ]]; then
			PPR_CMD+=(--x-prior --x-tsv "${OUT}.quickview.tsv" --nuc-cut-log10 "$NUC_CUT_LOG10" --x-slope "$X_SLOPE")
		fi
		oatk_run "${PPR_CMD[@]}"

		oatk_log 1 "[step] ensemble ${L} ($ENSEMBLE)"
		oatk_run python3 "$ENSEMBLE_PY" \
			--xlda "${OUT}.xcover.${L}.ids" \
			--ppr "${OUT}.ppr.${L}.ids" \
			--mode "$ENSEMBLE" \
			--label "$L" \
			-o "${OUT}.ensemble"

		if [[ "$EMIT_FASTQ" == "true" && -s "${OUT}.ensemble.${L}.ids" ]]; then
			oatk_log 1 "[emit] ${L} FASTQ"
			awk 'NF{print $1}' "${OUT}.ensemble.${L}.ids" | sort -u >"${OUT}.ensemble.${L}.names"
			oatk_run seqtk subseq "$READS" "${OUT}.ensemble.${L}.names" | gzip -c >"${OUT}.${L}.fastq.gz"
		fi

		[[ "$KEEP_INT" == "false" ]] && rm -f "$NAMES" "$FASTQ" || true
		return 0
	}

	# Always try mt; coverage-only fallback if none
	_ppr_for_label "mt" || {
		oatk_log 1 "[warn] coverage-only mt"
		cp -f "${OUT}.xcover.mt.ids" "${OUT}.ensemble.mt.ids" 2>/dev/null || : >"${OUT}.ensemble.mt.ids"
		: >"${OUT}.ensemble.mt.nuclear.ids"
		: >"${OUT}.ensemble.mt.report.tsv"
	}

	# Run pt only if anchors provided & IDs exist
	if [[ -n "${PT_ANCH:-}" && -s "${OUT}.xcover.pt.ids" ]]; then
		_ppr_for_label "pt" || {
			oatk_log 1 "[warn] coverage-only pt"
			cp -f "${OUT}.xcover.pt.ids" "${OUT}.ensemble.pt.ids" 2>/dev/null || : >"${OUT}.ensemble.pt.ids"
			: >"${OUT}.ensemble.pt.nuclear.ids"
			: >"${OUT}.ensemble.pt.report.tsv"
		}
	else
		: >"${OUT}.ensemble.pt.ids"
		: >"${OUT}.ensemble.pt.nuclear.ids"
		: >"${OUT}.ensemble.pt.report.tsv"
	fi

	oatk_log 1 "[done] mt: ${OUT}.ensemble.mt.ids  pt: ${OUT}.ensemble.pt.ids"
}

# ---------------- REPORT (select-mt-report) ----------------
_select_mt_report() {
	local cfg_path="" cfg_dir="" preset="" BANDAGE_SIZE="3000x3000" OUT=""
	local LOG_COUNT=0
	local ARGV=("$@") i=0
	while [[ $i -lt ${#ARGV[@]} ]]; do
		case "${ARGV[$i]}" in
		-v | --verbose)
			LOG_COUNT=$((LOG_COUNT + 1))
			i=$((i + 1))
			;;
		--config-path)
			cfg_path="${ARGV[$((i + 1))]:-}"
			i=$((i + 2))
			;;
		--config-dir)
			cfg_dir="${ARGV[$((i + 1))]:-}"
			i=$((i + 2))
			;;
		--preset)
			preset="${ARGV[$((i + 1))]:-}"
			i=$((i + 2))
			;;
		--bandage-size)
			BANDAGE_SIZE="${ARGV[$((i + 1))]:-}"
			i=$((i + 2))
			;;
		--out)
			OUT="${ARGV[$((i + 1))]:-}"
			i=$((i + 2))
			;;
		*) i=$((i + 1)) ;;
		esac
	done
	oatk_log_set_level "$LOG_COUNT"
	oatk_enable_trace_if_needed

	_load_config_env "$@"
	[[ -z "$OUT" ]] && OUT="${SELCFG_OUT_PREFIX:-}"
	[[ -n "$OUT" ]] || die "need --out or YAML out_prefix"
	[[ -f "$REPORT_PY" ]] || die "not found: $REPORT_PY"

	_find_gfa_twins() {
		local out="$1" outdir
		outdir="$(dirname -- "$out")"
		[[ "$outdir" == "." ]] && outdir="$PWD"
		local cand_mt=("${outdir}/mt1/assembly_graph.gfa" "${out}.asm.utg.final.gfa")
		local cand_pt=("${outdir}/pt1/assembly_graph.gfa" "${out}.asm.utg.final.gfa")
		shopt -s nullglob
		local any=("${outdir}"/*.gfa)
		shopt -u nullglob
		MT_GFA_PICK=""
		PT_GFA_PICK=""
		for f in "${cand_mt[@]}"; do [[ -s "$f" ]] && {
			MT_GFA_PICK="$f"
			break
		}; done
		for f in "${cand_pt[@]}"; do [[ -s "$f" ]] && {
			PT_GFA_PICK="$f"
			break
		}; done
		OVERVIEW_GFA_PICK="${MT_GFA_PICK:-$PT_GFA_PICK}"
		[[ -z "$OVERVIEW_GFA_PICK" && ${#any[@]} -gt 0 ]] && OVERVIEW_GFA_PICK="${any[0]}"
	}

	_bandage_image() {
		local gfa="$1" png="$2" size="$3"
		local h=${size#*x}
		Bandage image "$gfa" "$png" --height "$h" || {
			oatk_log 1 "[bandage] failed: $png"
			rm -f "$png" || true
		}
	}

	_find_gfa_twins "$OUT"
	if command -v Bandage >/dev/null 2>&1; then
		[[ -n "${OVERVIEW_GFA_PICK:-}" && -r "$OVERVIEW_GFA_PICK" ]] && _bandage_image "$OVERVIEW_GFA_PICK" "${OUT}.bandage.overview.png" "$BANDAGE_SIZE"
		[[ -s "${OUT}.ensemble.mt.ids" && -n "${MT_GFA_PICK:-}" && -r "$MT_GFA_PICK" ]] && _bandage_image "$MT_GFA_PICK" "${OUT}.bandage.mt.png" "$BANDAGE_SIZE"
		[[ -s "${OUT}.ensemble.pt.ids" && -n "${PT_GFA_PICK:-}" && -r "$PT_GFA_PICK" ]] && _bandage_image "$PT_GFA_PICK" "${OUT}.bandage.pt.png" "$BANDAGE_SIZE"
	else
		oatk_log 1 "[bandage] Bandage not found; skipping PNGs"
	fi

	oatk_run python3 "$REPORT_PY" \
		--prefix "$OUT" \
		--quickview "${OUT}.quickview.tsv" \
		--xlda-prefix "${OUT}.xcover" \
		--ppr-prefix "${OUT}.ppr.mt" \
		--ensemble-prefix "${OUT}.ensemble.mt" \
		${PT_GFA_PICK:+--ppr-pt-prefix "${OUT}.ppr.pt"} \
		${PT_GFA_PICK:+--ensemble-pt-prefix "${OUT}.ensemble.pt"} \
		--title "POLAP organelle selection (coverage + PPR + Bandage)" \
		-o "${OUT}.report.html"

	oatk_log 1 "[report] ${OUT}.report.html"
}

# ---------------- dispatch ----------------
sub="${1:-}"
case "$sub" in
select-mt)
	shift
	_select_mt "$@"
	;;
select-mt-report)
	shift
	_select_mt_report "$@"
	;;
"" | -h | --help) _usage ;;
*)
	oatk_log 1 "unknown subcommand: $sub"
	_usage
	exit 2
	;;
esac
