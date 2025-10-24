#!/usr/bin/env bash
set -euo pipefail

# polap-bash-oatk-tune-c.sh
# Probe SyncAsm with -c 1 and suggest a coverage cutoff (-c) from the valley
# between nuclear/organelle peaks (log10 coverage histogram).
#
# Addressing (one of):
#   --config-path FILE.yaml
#   --config-dir DIR --preset NAME   # resolves DIR/NAME.yaml (default: ~/.polap/profiles)
#
# Any flat YAML key may be overridden on the CLI by replacing '_' with '-'.
# Canonical output locator: out_prefix (YAML or --out-prefix); outdir = dirname(out_prefix).
#
# Common keys for this tool:
#   reads, out_prefix, preset, threads, k, s, hpc, plot
#
# Booleans (no value): --hpc/--no-hpc, --plot/--no-plot
# Verbose: -v/--verbose (repeatable). -vv enables xtrace with file:function:line PS4.
#
# Outputs:
#   <out_prefix>.suggest.txt
#   <out_prefix>.hist.tsv
#   <out_prefix>.probe.<gfa>
#   <out_prefix>.png               (if plot=true)
#   <outdir>/provenance.tune-c.txt

# ---------- paths ----------
POLAPLIB_DIR="${_POLAPLIB_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)}"
LOAD_PY="${POLAPLIB_DIR}/polap-py-config-load.py"
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
need python3
need syncasm
[[ -f "$LOAD_PY" ]] || die "not found: $LOAD_PY"

# ---------- load YAML -> TCFG_* ----------
RESOLVED_CFG=""
if [[ -n "$CFG_PATH" ]]; then
	oatk_log 1 "loading YAML --config-path ${CFG_PATH}"
	eval "$(
		python3 "$LOAD_PY" --path "$CFG_PATH" --format env --prefix TCFG 2>/dev/null
	)" || die "failed to load config: $CFG_PATH"
	RESOLVED_CFG="$CFG_PATH"
else
	[[ -n "$PRESET" ]] || die "need --preset when --config-path is not given"
	[[ -n "$CONFIG_DIR" ]] || CONFIG_DIR="${HOME}/.polap/profiles"
	oatk_log 1 "loading YAML --config-dir ${CONFIG_DIR} --preset ${PRESET}"
	eval "$(
		python3 "$LOAD_PY" --config-dir "$CONFIG_DIR" --preset "$PRESET" --format env --prefix TCFG 2>/dev/null
	)" || die "failed to load config: ${CONFIG_DIR}/${PRESET}.yaml"
	RESOLVED_CFG="${CONFIG_DIR}/${PRESET}.yaml"
fi

# ---------- apply ALL CLI overrides (CLI wins) ----------
apply_overrides() {
	local j=0
	while [[ $j -lt ${#ARGS[@]} ]]; do
		local tok="${ARGS[$j]}"
		case "$tok" in
		# boolean toggles (no value)
		--hpc)
			TCFG_HPC="true"
			j=$((j + 1))
			continue
			;;
		--no-hpc)
			TCFG_HPC="false"
			j=$((j + 1))
			continue
			;;
		--plot)
			TCFG_PLOT="true"
			j=$((j + 1))
			continue
			;;
		--no-plot)
			TCFG_PLOT="false"
			j=$((j + 1))
			continue
			;;

		# K/V mappings
		--*)
			local key="${tok#--}"
			key="${key//-/_}"
			local val="${ARGS[$((j + 1))]:-}"
			[[ -z "$val" || "$val" == --* ]] && die "missing value for $tok"
			case "$tok" in
			--out-prefix) TCFG_OUT_PREFIX="$val" ;;
			*)
				local varname="TCFG_${key^^}"
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

# ---------- normalize core fields ----------
[[ -n "${TCFG_READS:-}" ]] || die "reads is missing (YAML or --reads)"
[[ -n "${TCFG_OUT_PREFIX:-}" ]] || die "out_prefix is missing (YAML or --out-prefix)"
OUTDIR="$(dirname -- "${TCFG_OUT_PREFIX}")"
mkdir -p "$OUTDIR"

# preset defaults for k/s/hpc (allow override from YAML/CLI)
TCFG_PRESET="${TCFG_PRESET:-hifi}"
TCFG_THREADS="${TCFG_THREADS:-16}"

if [[ -z "${TCFG_K:-}" || -z "${TCFG_S:-}" ]]; then
	if [[ "${TCFG_PRESET,,}" == "ont" ]]; then
		TCFG_K="${TCFG_K:-41}"
		TCFG_S="${TCFG_S:-21}"
	else
		TCFG_K="${TCFG_K:-121}"
		TCFG_S="${TCFG_S:-27}"
	fi
fi

if [[ -z "${TCFG_HPC:-}" ]]; then
	if [[ "${TCFG_PRESET,,}" == "ont" ]]; then TCFG_HPC="true"; else TCFG_HPC="false"; fi
else
	TCFG_HPC="$(_norm_bool "${TCFG_HPC}")"
fi
HPC_FLAG=""
[[ "$TCFG_HPC" == "true" ]] && HPC_FLAG="--hpc"

TCFG_PLOT="$(_norm_bool "${TCFG_PLOT:-false}")"
PLOT_ON=0
[[ "$TCFG_PLOT" == "true" ]] && PLOT_ON=1

# ---------- provenance ----------
filtered_overrides=$(oatk_filter_overrides "${ARGS[@]}")
oatk_provenance_line "$RESOLVED_CFG" $filtered_overrides
oatk_provenance_append "${OUTDIR}/provenance.tune-c.txt" "$RESOLVED_CFG" $filtered_overrides

# ---------- Step 1: probe syncasm at -c 1 ----------
PROBE_PREFIX="${TCFG_OUT_PREFIX}.probe"
oatk_log 1 "[probe] syncasm -k ${TCFG_K} -s ${TCFG_S} $([[ "$TCFG_HPC" == "true" ]] && echo --hpc || :) -t ${TCFG_THREADS} -c 1"
oatk_run syncasm -k "$TCFG_K" -s "$TCFG_S" $([[ "$TCFG_HPC" == "true" ]] && echo --hpc || :) \
	-t "$TCFG_THREADS" -o "$PROBE_PREFIX" -c 1 -a 0.05 \
	--unzip-round 0 --max-bubble 0 --max-tip 0 --weak-cross 0 --no-read-ec \
	"${TCFG_READS}" || true

GFA=""
for g in "${PROBE_PREFIX}.utg.final.gfa" "${PROBE_PREFIX}.utg.gfa" "${PROBE_PREFIX}.utg.clean.gfa"; do
	[[ -s "$g" ]] && {
		GFA="$g"
		break
	}
done
[[ -n "$GFA" ]] || {
	oatk_log 1 "probe produced no GFA"
	echo "suggest: -c 50  (no GFA)" >"${TCFG_OUT_PREFIX}.suggest.txt"
	exit 3
}

cp -f "$GFA" "${TCFG_OUT_PREFIX}.probe.$(basename "$GFA")"

# ---------- Step 2: analyze coverage & suggest c ----------
HIST_TSV="${TCFG_OUT_PREFIX}.hist.tsv"
SUG_TXT="${TCFG_OUT_PREFIX}.suggest.txt"
PNG="${TCFG_OUT_PREFIX}.png"

oatk_run python3 - "$GFA" "$HIST_TSV" "$SUG_TXT" "$PLOT_ON" "$PNG" <<'PY'
import sys, math

gfa, hist_tsv, sug_txt, do_plot, png = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), sys.argv[5]

def parse_cov_from_gfa(path):
    cov = []
    with open(path, 'r', encoding='utf-8') as fh:
        for ln in fh:
            if not ln.startswith('S\t'):
                continue
            sp = ln.rstrip('\n').split('\t')
            if len(sp) < 3:
                continue
            seq = sp[2]
            L = float(len(seq))
            sc = None
            kc = None
            for t in sp[3:]:
                if t.startswith("SC:i:"):
                    try:
                        sc = float(t[5:])
                        break
                    except:
                        pass
            if sc is None:
                for t in sp[3:]:
                    if t.startswith("KC:i:"):
                        try:
                            kc = float(t[5:])
                            break
                        except:
                            pass
            if sc is not None and sc >= 0:
                cov.append(sc)
            elif kc is not None and L > 0:
                cov.append(kc / L)
    return cov

def smooth(y, rounds=2):
    x = y[:]
    for _ in range(rounds):
        if len(x) < 3: break
        x = [x[0]] + [0.25*x[i-1] + 0.5*x[i] + 0.25*x[i+1] for i in range(1, len(x)-1)] + [x[-1]]
    return x

def find_valley(logx, nbins=200):
    if not logx:
        return (None, None, (None,None), [], [])
    lo, hi = min(logx), max(logx)
    if hi <= lo:
        hi = lo + 1e-6
    step = (hi - lo) / nbins
    edges = [lo + i*step for i in range(nbins+1)]
    ctrs  = [(edges[i] + edges[i+1]) * 0.5 for i in range(nbins)]
    hist = [0]*nbins
    for v in logx:
        j = int((v - lo) / step)
        if j < 0: j = 0
        if j >= nbins: j = nbins - 1
        hist[j] += 1
    sh = smooth(hist, 2)

    # left mode
    m1 = max(range(nbins), key=lambda i: sh[i]) if nbins>0 else None
    if m1 is None:
        return (None, None, (None,None), ctrs, hist)

    # second mode at least ~10% bins away
    m2 = None; best = -1
    for i in range(nbins):
        if i == m1: continue
        if abs(i - m1) >= nbins//10 and sh[i] > best:
            best = sh[i]; m2 = i

    if m2 is None:
        # fallback valley to the right of first mode
        vpos = None
        for i in range(m1+1, nbins-1):
            if sh[i] <= sh[i-1] and sh[i] <= sh[i+1]:
                vpos = i; break
        if vpos is None:
            vpos = min(nbins-1, m1 + nbins//8)
        cut_log = ctrs[vpos]; cut_lin = 10**cut_log
        return (cut_log, cut_lin, (edges[vpos], edges[vpos+1]), ctrs, hist)

    # ensure m1<m2
    if m2 < m1: m1, m2 = m2, m1

    vpos = m1 + 1; bestv = sh[vpos]; besti = vpos
    for i in range(m1+1, m2):
        if sh[i] < bestv:
            bestv = sh[i]; besti = i
    cut_log = ctrs[besti]; cut_lin = 10**cut_log
    return (cut_log, cut_lin, (edges[besti], edges[besti+1]), ctrs, hist)

cov = parse_cov_from_gfa(gfa)
if not cov:
    open(hist_tsv, "w").write("")
    open(sug_txt, "w").write("suggest: -c 50  (no coverage derived)\n")
    print("NO_COVERAGE")
    sys.exit(0)

logx = [math.log10(max(1e-9, c)) for c in cov]
cut_log, cut_lin, (blo, bhi), ctrs, hist = find_valley(logx, nbins=200)

with open(hist_tsv, "w") as fh:
    fh.write("bin_center_log10\tcount\n")
    for c,h in zip(ctrs,hist):
        fh.write(f"{c:.6f}\t{h}\n")

suggest = int(round(cut_lin)) if (cut_lin and cut_lin>0) else 50
with open(sug_txt, "w") as fh:
    if cut_log is None:
        fh.write(f"suggest: -c {suggest}  (no valley; fallback)\n")
    else:
        fh.write(f"suggest: -c {suggest}  (log10 valley≈{cut_log:.3f}; bin={blo:.2f}–{bhi:.2f})\n")

if do_plot:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        if len(ctrs) > 1:
            bw = (ctrs[1] - ctrs[0])
        else:
            bw = 0.01
        plt.figure(figsize=(7.0,4.0))
        plt.bar(ctrs, hist, width=bw, alpha=0.6)
        if cut_log is not None:
            plt.axvline(cut_log, linestyle="--", linewidth=2)
        plt.xlabel("log10(segment coverage)")
        plt.ylabel("count (segments)")
        plt.tight_layout(); plt.savefig(png, dpi=150); plt.close()
    except Exception:
        pass
PY

if grep -q "no coverage" "${SUG_TXT}" 2>/dev/null; then
	oatk_log 1 "[warn] no coverage derived; wrote fallback suggestion"
fi

oatk_log 1 "[tune] wrote ${SUG_TXT}"
oatk_log 1 "[tune] wrote ${HIST_TSV}"
[[ $PLOT_ON -eq 1 ]] && oatk_log 1 "[tune] wrote ${PNG}"

exit 0
