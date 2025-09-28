#!/usr/bin/env bash
set -euo pipefail

# polap-bash-oatk-score-seeds.sh
# YAML-config + CLI overrides (CLI wins). Scores seed contigs by:
#   minimap2 mapping → coverage (breadth, depth, N50, len_total) →
#   syncfilter quickview on seeds (X proxy) →
#   HMM gene content (hmmannot or hmmsearch) →
#   composite score → append CSV row.
#
# Addressing (one of):
#   --config-path FILE.yaml
#   --config-dir DIR --preset NAME        # resolves DIR/NAME.yaml (default: ~/.polap/profiles)
#
# Any flat YAML key can be overridden on CLI by replacing '_' with '-'.
# Canonical output locator: out_prefix (YAML or --out-prefix); outdir = dirname(out_prefix).
#
# Required (YAML or CLI):
#   seeds, reads, hmm_db, out_prefix, preset, threads
#
# Optional knobs:
#   k, s                                  (only used for X-proxy quickview)
#   nuc_cut_log10, x_slope               (X-prior for x_score bending)
#   min_shared, jaccard_min, topk_nei, steps   (tracked in CSV for provenance)
#
# Verbosity:
#   -v / --verbose   (repeatable; -vv enables xtrace with file:function:line PS4)

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
need minimap2
need samtools
need syncfilter
[[ -f "$LOAD_PY" ]] || die "not found: $LOAD_PY"

# ---------- load YAML → SCFG_* ----------
RESOLVED_CFG=""
if [[ -n "$CFG_PATH" ]]; then
	oatk_log 1 "loading YAML --config-path ${CFG_PATH}"
	eval "$(
		python3 "$LOAD_PY" --path "$CFG_PATH" --format env --prefix SCFG 2>/dev/null
	)" || die "failed to load config: $CFG_PATH"
	RESOLVED_CFG="$CFG_PATH"
else
	[[ -n "$PRESET" ]] || die "need --preset when --config-path is not given"
	[[ -n "$CONFIG_DIR" ]] || CONFIG_DIR="${HOME}/.polap/profiles"
	oatk_log 1 "loading YAML --config-dir ${CONFIG_DIR} --preset ${PRESET}"
	eval "$(
		python3 "$LOAD_PY" --config-dir "$CONFIG_DIR" --preset "$PRESET" --format env --prefix SCFG 2>/dev/null
	)" || die "failed to load config: ${CONFIG_DIR}/${PRESET}.yaml"
	RESOLVED_CFG="${CONFIG_DIR}/${PRESET}.yaml"
fi

# ---------- apply ALL CLI overrides (CLI wins) ----------
apply_overrides() {
	local j=0
	while [[ $j -lt ${#ARGS[@]} ]]; do
		local tok="${ARGS[$j]}"
		case "$tok" in
		-v | --verbose)
			j=$((j + 1))
			continue
			;;
		--config-path | --config-dir | --preset)
			j=$((j + 2))
			continue
			;;
		--*) # --foo-bar => SCFG_FOO_BAR
			local key="${tok#--}"
			key="${key//-/_}"
			local val="${ARGS[$((j + 1))]:-}"
			[[ -z "$val" || "$val" == --* ]] && die "missing value for $tok"
			case "$tok" in
			--out-prefix) SCFG_OUT_PREFIX="$val" ;;
			*)
				local varname="SCFG_${key^^}"
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
SCFG_PRESET="${SCFG_PRESET:-hifi}"
SCFG_THREADS="${SCFG_THREADS:-16}"

[[ -n "${SCFG_SEEDS:-}" ]] || die "seeds missing (YAML or --seeds)"
[[ -n "${SCFG_READS:-}" ]] || die "reads missing (YAML or --reads)"
[[ -n "${SCFG_HMM_DB:-}" ]] || die "hmm_db missing (YAML or --hmm-db)"
[[ -n "${SCFG_OUT_PREFIX:-}" ]] || die "out_prefix missing (YAML or --out-prefix)"

OUTDIR="$(dirname -- "${SCFG_OUT_PREFIX}")"
mkdir -p "$OUTDIR"

# k/s defaults by preset (only for X-proxy quickview)
if [[ -z "${SCFG_K:-}" || -z "${SCFG_S:-}" ]]; then
	if [[ "${SCFG_PRESET,,}" == "ont" ]]; then
		SCFG_K="${SCFG_K:-41}"
		SCFG_S="${SCFG_S:-21}"
	else
		SCFG_K="${SCFG_K:-121}"
		SCFG_S="${SCFG_S:-27}"
	fi
fi

# optional knobs
SCFG_NUC_CUT_LOG10="${SCFG_NUC_CUT_LOG10:-}"
SCFG_X_SLOPE="${SCFG_X_SLOPE:-}"
SCFG_MIN_SHARED="${SCFG_MIN_SHARED:-}"
SCFG_JACCARD_MIN="${SCFG_JACCARD_MIN:-}"
SCFG_TOPK_NEI="${SCFG_TOPK_NEI:-}"
SCFG_STEPS="${SCFG_STEPS:-}"

# ---------- provenance ----------
filtered_overrides=$(oatk_filter_overrides "${ARGS[@]}")
oatk_provenance_line "$RESOLVED_CFG" $filtered_overrides
oatk_provenance_append "${OUTDIR}/provenance.score-seeds.txt" "$RESOLVED_CFG" $filtered_overrides

# ---------- 1) Mapping & coverage ----------

# Ensure BAM target is defined (set -u safe)
BAM="${SCFG_OUT_PREFIX}.bam"

# Choose minimap2 preset as an ARRAY (no word-splitting surprises)
MMIDX_OPTS=(-x map-hifi)
if [[ "${SCFG_PRESET,,}" == "ont" ]]; then
	MMIDX_OPTS=(-x map-ont)
fi

# Build commands as arrays (force SAM with -a so samtools sees a header)
MM_CMD=(minimap2 "${MMIDX_OPTS[@]}" -a -t "${SCFG_THREADS}" "${SCFG_SEEDS}" "${SCFG_READS}")
SORT_CMD=(samtools sort -@ "${SCFG_THREADS}" -O BAM -o "$BAM" -)

# Log + run the pipeline with pipefail, so failures propagate
oatk_log 1 "exec: $(oatk_cmdstr "${MM_CMD[@]}") | $(oatk_cmdstr "${SORT_CMD[@]}")"
set -o pipefail
"${MM_CMD[@]}" | "${SORT_CMD[@]}"
st=$?
set +o pipefail
oatk_log 1 "exit(${st}): $(oatk_cmdstr "${MM_CMD[@]}") | $(oatk_cmdstr "${SORT_CMD[@]}")"

# Index BAM
oatk_run samtools index -@ "${SCFG_THREADS}" "$BAM"

DEPTH_TSV="${SCFG_OUT_PREFIX}.depth.tsv"
oatk_log 1 "[depth] samtools depth -aa"
oatk_run samtools depth -aa "$BAM" >"$DEPTH_TSV"

COV_JSON="$(dirname -- "${SCFG_OUT_PREFIX}")/._tmp_cov_${RANDOM}.json"
oatk_run python3 - "${SCFG_SEEDS}" "$DEPTH_TSV" "$COV_JSON" <<'PY'
import sys, json
from collections import defaultdict

fa, depth_tsv, outj = sys.argv[1], sys.argv[2], sys.argv[3]

def fa_lengths(path):
    lens = {}
    with open(path, 'r') as fh:
        name=None; L=0
        for ln in fh:
            if ln.startswith('>'):
                if name is not None: lens[name]=L
                name = ln[1:].strip().split()[0]; L=0
            else:
                L += len(ln.strip())
        if name is not None: lens[name]=L
    return lens

lens = fa_lengths(fa)
total = sum(lens.values()) or 1
covered = 0
sum_depth = 0

with open(depth_tsv, 'r') as fh:
    for ln in fh:
        sp = ln.rstrip('\n').split('\t')
        if len(sp) < 3: continue
        d = int(sp[2])
        if d>0: covered += 1
        sum_depth += d

breadth = covered / total
depth_mean = sum_depth / total

# N50
lengths_sorted = sorted(lens.values(), reverse=True)
half = total/2
acc=0; n50=0
for L in lengths_sorted:
    acc+=L
    if acc>=half:
        n50=L; break

with open(outj, 'w') as fh:
    json.dump({
        "breadth_weighted": breadth,
        "depth_mean": depth_mean,
        "len_total": total,
        "n50": n50
    }, fh)
PY

# ---------- 2) X proxy via syncfilter quickview on seeds ----------
PSEUDO_FQ="$(dirname -- "${SCFG_OUT_PREFIX}")/._tmp_seeds_${RANDOM}.fq"
XTSV="${SCFG_OUT_PREFIX}.xcover.tsv"

oatk_run python3 - "${SCFG_SEEDS}" "$PSEUDO_FQ" <<'PY'
import sys
fa, fq = sys.argv[1], sys.argv[2]
with open(fa) as fi, open(fq,'w') as fo:
    name=None; seq=[]
    def flush():
        if name is None: return
        s=''.join(seq); q='~'*len(s)
        fo.write(f"@{name}\n{s}\n+\n{q}\n")
    for ln in fi:
        if ln.startswith('>'):
            flush(); name = ln[1:].strip().split()[0]; seq=[]
        else:
            seq.append(ln.strip())
    flush()
PY

HPC_FLAG=""
[[ "${SCFG_PRESET,,}" == "ont" ]] && HPC_FLAG="--hpc"
oatk_log 1 "[xcover] syncfilter quickview k=${SCFG_K} s=${SCFG_S} ${HPC_FLAG}"
oatk_run syncfilter --mode quickview $HPC_FLAG -k "${SCFG_K}" -s "${SCFG_S}" -t "${SCFG_THREADS}" -o "$(dirname -- "$PSEUDO_FQ")/x" "$PSEUDO_FQ"
mv "$(dirname -- "$PSEUDO_FQ")/x.syncfilter.tsv" "$XTSV"

XJSON="$(dirname -- "${SCFG_OUT_PREFIX}")/._tmp_x_${RANDOM}.json"
oatk_run python3 - "$XTSV" "$XJSON" "${SCFG_NUC_CUT_LOG10:-}" "${SCFG_X_SLOPE:-}" <<'PY'
import sys, json, math
tsv, outj, cut, slope = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
xs=[]
with open(tsv) as fh:
    hdr=None
    for ln in fh:
        sp=ln.rstrip('\n').split('\t')
        if sp[0]=='read_id': hdr=sp; continue
        try:
            med=float(sp[3])
            if med>0: xs.append(math.log10(med))
        except: pass
xmu = (sum(xs)/len(xs)) if xs else 0.0
def sigmoid(z): 
    try: import math; return 1.0/(1.0+math.exp(-z))
    except Exception: return 0.5
def bend(xmu, cut, slope):
    try:
        c=float(cut); m=float(slope)
        return sigmoid(m*(xmu - c))
    except:
        # fallback center ~1.0 log10
        return sigmoid(3.0*(xmu - 1.0))
xscore = bend(xmu, cut, slope)
with open(outj,'w') as fh:
    json.dump({"x_mu_log10": xmu, "x_score": xscore}, fh)
PY

# ---------- 3) HMM gene content ----------
HMM_TSV="${SCFG_OUT_PREFIX}.hmm.tsv"
if command -v hmmannot >/dev/null 2>&1; then
	oatk_log 1 "[hmm] hmmannot"
	oatk_run hmmannot --db "${SCFG_HMM_DB}" --seq "${SCFG_SEEDS}" --out "$HMM_TSV" || true
elif command -v hmmsearch >/dev/null 2>&1; then
	oatk_log 1 "[hmm] hmmsearch fallback"
	oatk_run hmmsearch --noali --tblout "${SCFG_OUT_PREFIX}.hmm.tbl" "${SCFG_HMM_DB}" "${SCFG_SEEDS}" >/dev/null 2>&1 || true
	awk '!/^#/ && NF>=18 {print $1"\t"$3"\t"$18}' "${SCFG_OUT_PREFIX}.hmm.tbl" >"$HMM_TSV" || true
else
	oatk_log 1 "[warn] neither hmmannot nor hmmsearch available; gene score=0"
	: >"$HMM_TSV"
fi

GENE_JSON="$(dirname -- "${SCFG_OUT_PREFIX}")/._tmp_gene_${RANDOM}.json"
oatk_run python3 - "$HMM_TSV" "$GENE_JSON" <<'PY'
import sys, json
tsv, outj = sys.argv[1], sys.argv[2]
genes=set()
try:
    with open(tsv) as fh:
        for ln in fh:
            ln=ln.strip()
            if not ln or ln.startswith('#'): continue
            genes.add(ln.split()[0])
    found=len(genes)
except:
    found=0
d={"genes_found": found, "genes_total": 0, "genes_score": (1.0 if found>0 else 0.0)}
with open(outj,'w') as fh:
    json.dump(d, fh)
PY

# ---------- 4) Compose CSV row ----------
CSV="${SCFG_OUT_PREFIX}.seed_score.csv"
if [[ ! -s "$CSV" ]]; then
	echo "out,seeds_fa,reads_fq,genes_found,genes_total,genes_score,breadth_weighted,depth_mean,n50,len_total,x_mu_log10,x_score,score,k,s,hpc,preset,threads,nuc_cut_log10,x_slope,min_shared,jaccard_min,topk_nei,steps" >"$CSV"
fi

oatk_run python3 - "${SCFG_OUT_PREFIX}" "${SCFG_SEEDS}" "${SCFG_READS}" "$COV_JSON" "$XJSON" "$GENE_JSON" "$CSV" "${SCFG_PRESET}" "${SCFG_THREADS}" "${SCFG_K}" "${SCFG_S}" "${SCFG_NUC_CUT_LOG10:-}" "${SCFG_X_SLOPE:-}" "${SCFG_MIN_SHARED:-}" "${SCFG_JACCARD_MIN:-}" "${SCFG_TOPK_NEI:-}" "${SCFG_STEPS:-}" <<'PY'
import sys, json, math
out, seeds, reads, covj, xj, gj, csv, preset, threads, k, s, cut, slope, min_sh, jmin, topk, steps = sys.argv[1:]

def jload(p, dflt):
    try:
        with open(p) as fh: return json.load(fh)
    except Exception: return dflt

cov = jload(covj, {"breadth_weighted":0.0,"depth_mean":0.0,"len_total":0,"n50":0})
xv  = jload(xj,  {"x_mu_log10":0.0,"x_score":0.0})
gen = jload(gj,  {"genes_found":0,"genes_total":0,"genes_score":0.0})

breadth = float(cov.get("breadth_weighted",0.0))
depth   = float(cov.get("depth_mean",0.0))
len_tot = int(cov.get("len_total",0))
n50     = int(cov.get("n50",0))
xmu     = float(xv.get("x_mu_log10",0.0))
xscore  = float(xv.get("x_score",0.0))
gf      = int(gen.get("genes_found",0))
gt      = int(gen.get("genes_total",0))
gscore  = float(gen.get("genes_score",0.0))

# composite score (weights can be tuned)
score = 0.50*gscore + 0.30*breadth + 0.15*xscore + 0.05*min(1.0,(n50/100000.0 if 100000.0>0 else 0.0))
hpc = ("true" if preset=="ont" else "false")

def fmt(x, nd=6):
    try: return f"{float(x):.{nd}f}"
    except Exception: return str(x)

row = [
    out, seeds, reads,
    str(gf), str(gt), fmt(gscore),
    fmt(breadth), fmt(depth), str(n50), str(len_tot),
    fmt(xmu), fmt(xscore), fmt(score),
    str(k), str(s), hpc, preset, str(threads),
    (cut if cut else ""), (slope if slope else ""),
    (min_sh if min_sh else ""), (jmin if jmin else ""),
    (topk if topk else ""), (steps if steps else "")
]
with open(csv,"a",encoding="utf-8") as fh:
    fh.write(",".join(row) + "\n")
PY

# cleanup temps
rm -f "$COV_JSON" "$XJSON" "$GENE_JSON" "$PSEUDO_FQ" 2>/dev/null || true

oatk_log 1 "[score] wrote ${CSV}"
exit 0
