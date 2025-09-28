#!/usr/bin/env bash
set -euo pipefail

POLAP_ORGSEL_VERSION="0.9.0"

# ---- locate library dir (scripts live here) ----
_POLAPLIB_DIR="${_POLAPLIB_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)}"

# ---- defaults ----
PRESET=""
STAGE="graph-early" # quickview | graph-early | graph-final
OUTDIR="out"
THREADS=8
SIMULATE=false
REAL_READS=""
REDO=false

TAIL=0.10
WIN=0.20
CUT_METHOD="hybrid"
KEEP_NUCLEAR=false
PT_IDS="" # --pt pt.id.all.txt
MT_IDS="" # --mt mt.id.all.txt

# sim defaults
NUC_SIZE="10m"
NUC_COV=5
MT_SIZE="500k"
MT_COV=50
PT_SIZE="150k"
PT_COV=250

HIFI_RL_MEAN=10000
HIFI_RL_SD=1500
HIFI_SUB=0.002
HIFI_INS=0.0003
HIFI_DEL=0.0003
ONT_RL_MEAN=10000
ONT_RL_SD=3000
ONT_SUB=0.01
ONT_INS=0.01
ONT_DEL=0.01

# k/s defaults
K_HIFI_QV=121
S_HIFI_QV=27
HPC_HIFI=""
K_HIFI_GE=121
S_HIFI_GE=27
K_ONT_QV=41
S_ONT_QV=17
HPC_ONT="--hpc"
K_ONT_GE=41
S_ONT_GE=21

usage() {
	cat <<EOF
polap-bash-organelle-select.sh v${POLAP_ORGSEL_VERSION}

Usage:
  $(basename "$0") --preset {hifi|ont} [--simulate | --reads FILE] --stage {quickview|graph-early|graph-final} [options]

Options:
  --out DIR         Output dir (default: out)
  --threads INT     Threads (default: 8)
  --redo            Remove output dir first
  --tail FLOAT      Right-tail for cut (default: 0.10)
  --win FLOAT       Window width for cut-auto (log10; default: 0.20)
  --method STR      Cut method (mad|window|std|valley|hybrid; default: hybrid)
  --keep-nuclear    Use permissive syncasm cleanup
  --polaplib DIR    Override library directory

Anchors (recommended):
  --pt FILE         File of plastid-anchored read IDs (one per line)
  --mt FILE         File of mitochondrial-anchored read IDs (one per line)

Meta:
  --version         Print version and exit
  -h | --help       Show help and exit
EOF
	exit 1
}

# ---- parse args ----
while [[ $# -gt 0 ]]; do
	case "$1" in
	--preset)
		PRESET="$2"
		shift 2
		;;
	--stage)
		STAGE="$2"
		shift 2
		;;
	--simulate)
		SIMULATE=true
		shift
		;;
	--reads)
		REAL_READS="$2"
		shift 2
		;;
	--out)
		OUTDIR="$2"
		shift 2
		;;
	--threads)
		THREADS="$2"
		shift 2
		;;
	--redo)
		REDO=true
		shift
		;;
	--tail)
		TAIL="$2"
		shift 2
		;;
	--win)
		WIN="$2"
		shift 2
		;;
	--method)
		CUT_METHOD="$2"
		shift 2
		;;
	--keep-nuclear)
		KEEP_NUCLEAR=true
		shift
		;;
	--polaplib)
		_POLAPLIB_DIR="$2"
		shift 2
		;;
	--pt)
		PT_IDS="$2"
		shift 2
		;;
	--mt)
		MT_IDS="$2"
		shift 2
		;;
	--version)
		echo "${POLAP_ORGSEL_VERSION}"
		exit 0
		;;
	-h | --help) usage ;;
	*)
		echo "[error] unknown option: $1"
		usage
		;;
	esac
done

# ---- tools ----
CUT_AUTO="${_POLAPLIB_DIR}/polap-py-syncfilter-cut-auto.py"
SIM_HIFI="${_POLAPLIB_DIR}/simulate_wgs_hifi.py"
SIM_ONT="${_POLAPLIB_DIR}/simulate_wgs_ont.py"

[[ -z "$PRESET" ]] && {
	echo "[error] --preset required"
	usage
}
if ! $SIMULATE && [[ -z "$REAL_READS" ]]; then
	echo "[error] need --simulate or --reads"
	usage
fi
case "$STAGE" in quickview | graph-early | graph-final) : ;; *)
	echo "[error] bad --stage"
	exit 2
	;;
esac

need() { command -v "$1" >/dev/null 2>&1 || {
	echo "[error] missing $1"
	exit 127
}; }
need syncfilter
need syncasm
need python3
need seqtk
[[ -f "$CUT_AUTO" ]] || {
	echo "[error] not found: $CUT_AUTO"
	exit 127
}

# ---- setup ----
if [[ "$REDO" == true && -d "$OUTDIR" ]]; then rm -rf "$OUTDIR"; fi
mkdir -p "$OUTDIR"
PREFIX="$OUTDIR/sample"
ASM_PREFIX="$OUTDIR/asm"
READS_OUT="$OUTDIR/wgs_mix.fq" # plain FASTQ working input

# ---- k/s presets ----
if [[ "$PRESET" == "hifi" ]]; then
	K_QV=$K_HIFI_QV
	S_QV=$S_HIFI_QV
	HPC_FLAG="$HPC_HIFI"
	K_GE=$K_HIFI_GE
	S_GE=$S_HIFI_GE
	SIM_SCRIPT="$SIM_HIFI"
	SIM_ARGS=(--rl-mean "$HIFI_RL_MEAN" --rl-sd "$HIFI_RL_SD" --sub "$HIFI_SUB" --ins "$HIFI_INS" --dele "$HIFI_DEL")
else
	K_QV=$K_ONT_QV
	S_QV=$S_ONT_QV
	HPC_FLAG="$HPC_ONT"
	K_GE=$K_ONT_GE
	S_GE=$S_ONT_GE
	SIM_SCRIPT="$SIM_ONT"
	SIM_ARGS=(--rl-mean "$ONT_RL_MEAN" --rl-sd "$ONT_RL_SD" --sub "$ONT_SUB" --ins "$ONT_INS" --dele "$ONT_DEL")
fi

# ---- syncasm cleanup policy ----
if $KEEP_NUCLEAR; then
	SYNCASM_OPTS_EARLY=(-c 1 -a 0.05 --unzip-round 0 --max-bubble 0 --max-tip 0 --weak-cross 0 --no-read-ec)
	SYNCASM_OPTS_FINAL=(-c 1 -a 0.05 --unzip-round 0 --max-bubble 0 --max-tip 0 --weak-cross 0 --no-read-ec)
else
	SYNCASM_OPTS_EARLY=()
	SYNCASM_OPTS_FINAL=()
fi

# ---- acquire plain FASTQ input ----
if $SIMULATE; then
	echo "[step] simulate WGS ($PRESET) -> plain FASTQ"
	python3 "$SIM_SCRIPT" \
		--nuclear-size "$NUC_SIZE" --nuclear-cov "$NUC_COV" \
		--mito-size "$MT_SIZE" --mito-cov "$MT_COV" \
		--plastid-size "$PT_SIZE" --plastid-cov "$PT_COV" \
		"${SIM_ARGS[@]}" --emit-refs \
		-o "$READS_OUT" \
		--summary "$OUTDIR/wgs_mix.summary.tsv" || true
	if [[ ! -s "$READS_OUT" && -s "$OUTDIR/wgs_mix.fq.gz" ]]; then
		gzip -dc "$OUTDIR/wgs_mix.fq.gz" >"$READS_OUT"
	fi
	[[ -s "$READS_OUT" ]] || {
		echo "[error] no reads at $READS_OUT"
		exit 3
	}
else
	echo "[step] real reads -> ensure plain FASTQ"
	if [[ "$REAL_READS" == *.gz ]]; then gzip -dc "$REAL_READS" >"$READS_OUT"; else ln -sf "$(realpath "$REAL_READS")" "$READS_OUT"; fi
fi

run_cut_auto() {
	local TSV_IN="$1"
	local TAG="$2"
	local PNG="${PREFIX}.${TAG}.png"
	local CUTS="${PREFIX}.${TAG}.cuts.tsv"
	local TSV_OUT="${PREFIX}.${TAG}.reclass.tsv" # class rewritten to nuclear/pt/mt
	local EXTRA=()
	[[ -n "$PT_IDS" ]] && EXTRA+=(--pt "$PT_IDS")
	[[ -n "$MT_IDS" ]] && EXTRA+=(--mt "$MT_IDS")
	python3 "$CUT_AUTO" -i "$TSV_IN" \
		--method "$CUT_METHOD" --tail "$TAIL" --window-width "$WIN" \
		--cuts-out "$CUTS" --png "$PNG" --tsv-out "$TSV_OUT" "${EXTRA[@]}"
	echo "$CUTS"
}

emit_bins_from_ids() {
	local READS="$1"
	local PFX="$2"
	# for t in nuclear pt mt; do
	for t in pt mt; do
		if [[ -s "${PFX}.${t}.ids" ]]; then
			awk 'NF{print $1}' "${PFX}.${t}.ids" | tr -d '\r' | sort -u >"${PFX}.${t}.names"
			seqtk subseq "$READS" "${PFX}.${t}.names" | gzip -c >"${PFX}.${t}.fastq.gz"
			echo "[emit] ${PFX}.${t}.fastq.gz"
		fi
	done
}

# ---- Stage 1: quickview (plain FASTQ) ----
echo "[cmd] syncfilter --mode quickview $HPC_FLAG -k $K_QV -s $S_QV -t $THREADS -o $PREFIX $READS_OUT"
syncfilter --mode quickview $HPC_FLAG -k "$K_QV" -s "$S_QV" -t "$THREADS" -o "$PREFIX" "$READS_OUT"
mv -f "${PREFIX}.syncfilter.tsv" "${PREFIX}.quickview.tsv"

# use anchors (pt/mt) to learn bands; also writes ids files
Q_CUTS=$(run_cut_auto "${PREFIX}.quickview.tsv" "quickview")

if [[ "$STAGE" == "quickview" ]]; then
	emit_bins_from_ids "$READS_OUT" "$PREFIX.$STAGE"
	echo "[done] quickview"
	exit 0
fi

# ---- Stage 2: graph-early ----
echo "[step] syncasm raw graph"
syncasm -k "$K_GE" -s "$S_GE" -t "$THREADS" "${SYNCASM_OPTS_EARLY[@]}" -o "$ASM_PREFIX" "$READS_OUT"
RAW_GFA="${ASM_PREFIX}.utg.gfa"
[[ -f "$RAW_GFA" ]] || RAW_GFA="${ASM_PREFIX}.utg.raw.gfa"
[[ -f "$RAW_GFA" ]] || {
	echo "[error] raw GFA not found"
	exit 4
}

syncfilter --mode graph-early --gfa "$RAW_GFA" $HPC_FLAG -k "$K_GE" -s "$S_GE" -t "$THREADS" -o "$PREFIX" "$READS_OUT"
mv -f "${PREFIX}.syncfilter.tsv" "${PREFIX}.graph.early.tsv"
GE_CUTS=$(run_cut_auto "${PREFIX}.graph.early.tsv" "graph.early")

if [[ "$STAGE" == "graph-early" ]]; then
	emit_bins_from_ids "$READS_OUT" "$PREFIX.$STAGE"
	echo "[done] graph-early"
	exit 0
fi

# ---- Stage 3: graph-final ----
FINAL_GFA="${ASM_PREFIX}.utg.final.gfa"
[[ -f "$FINAL_GFA" ]] || { [[ -f "${ASM_PREFIX}.utg.clean.gfa" ]] && FINAL_GFA="${ASM_PREFIX}.utg.clean.gfa"; }
if [[ ! -f "$FINAL_GFA" ]]; then
	syncasm -k "$K_GE" -s "$S_GE" -t "$THREADS" "${SYNCASM_OPTS_FINAL[@]}" -o "$ASM_PREFIX" "$READS_OUT"
	FINAL_GFA="${ASM_PREFIX}.utg.final.gfa"
	[[ -f "$FINAL_GFA" ]] || { [[ -f "${ASM_PREFIX}.utg.clean.gfa" ]] && FINAL_GFA="${ASM_PREFIX}.utg.clean.gfa"; }
fi
[[ -f "$FINAL_GFA" ]] || {
	echo "[error] final GFA not found"
	exit 5
}

syncfilter --mode graph-final --gfa "$FINAL_GFA" $HPC_FLAG -k "$K_GE" -s "$S_GE" -t "$THREADS" -o "$PREFIX" "$READS_OUT"
mv -f "${PREFIX}.syncfilter.tsv" "${PREFIX}.graph.final.tsv"
GF_CUTS=$(run_cut_auto "${PREFIX}.graph.final.tsv" "graph.final")

emit_bins_from_ids "$READS_OUT" "$PREFIX"
emit_bins_from_ids "$READS_OUT" "$PREFIX.$STAGE"
echo "[done] graph-final"
