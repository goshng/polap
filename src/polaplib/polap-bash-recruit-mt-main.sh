#!/usr/bin/env bash
# polap-bash-recruit-mt-main.sh  v0.4.0
# Driver for permissive mitochondrial read recruitment:
#   1) (optional) plastid two-isomer builder
#   2) plastid subtraction (competitive vs PT isomers)
#   3) nuclear subtraction (vs --nuc-ref)
#   4) iterative k-mer signature recruit (5% gain stop)
#   5) iterative read-overlap recruit (5% gain stop)
#   6) (optional) quick coverage check
#
# Logs & TSVs are written into --outdir (single rolling log file).
# All external scripts are expected in ${_POLAPLIB_DIR}.
# Required tools depend on the called sub-scripts (minimap2, samtools, sourmash/mash, etc).

set -euo pipefail

VERSION="v0.4.0"

# ──────────────────────────────────────────────────────────────────────────────
# Helpers (no ugly quoting)
#   _run: echo + optionally exec a *command with args*, safely; honors --dry
#   _sh : echo + optionally exec a *shell line* (for pipes/redirs); honors --dry
# ──────────────────────────────────────────────────────────────────────────────

_DR=0   # dry-run flag (0=execute, 1=print only)
_LOG="" # path to unified log file (set after parsing)

_log0() { printf "[%s] %s\n" "$(date +'%F %T')" "$*" | tee -a "$_LOG" >&2; }
_log1() { [[ $VERBOSE -gt 0 ]] && printf "[%s] %s\n" "$(date +'%F %T')" "$*" | tee -a "$_LOG" >&2; }

_run() {
	# usage: _run <cmd> [arg...]
	if [[ ${_DR} -eq 1 ]]; then
		printf "[DRY] %q" "$1" | tee -a "$_LOG"
		shift
		for a in "$@"; do printf " %q" "$a" | tee -a "$_LOG"; done
		printf "\n" | tee -a "$_LOG"
		return 0
	fi
	# print and exec (args are already split)
	printf "[RUN] %q" "$1" | tee -a "$_LOG"
	shift
	for a in "$@"; do printf " %q" "$a" | tee -a "$_LOG"; done
	printf "\n" | tee -a "$_LOG"
	"$@"
}

_sh() {
	# usage: _sh "full shell line with pipes/redirs"
	local line="$1"
	if [[ ${_DR} -eq 1 ]]; then
		printf "[DRY] %s\n" "$line" | tee -a "$_LOG"
		return 0
	fi
	printf "[RUN] %s\n" "$line" | tee -a "$_LOG"
	bash -c "$line"
}

# ──────────────────────────────────────────────────────────────────────────────
# Defaults & CLI
# ──────────────────────────────────────────────────────────────────────────────

_POLAPLIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
: "${_POLAPLIB_DIR:?need _POLAPLIB_DIR in env (path to polap scripts)}"

PT_TWOFORMS="${_POLAPLIB_DIR}/polap-bash-pt-twoforms-from-IR-SSC.sh"
RUN_COMP="${_POLAPLIB_DIR}/polap-bash-recruit-competitive.sh"
RUN_KSIG="${_POLAPLIB_DIR}/polap-bash-recruit-kmer-signature.sh"
RUN_OVLP="${_POLAPLIB_DIR}/polap-bash-recruit-read-overlap.sh"
RUN_COVQ="${_POLAPLIB_DIR}/polap-bash-cov-quickcheck.sh"

# inputs
READS=""   # input FASTQ(.gz)/FASTA(.gz)
PT_REF=""  # plastid reference (single FASTA)
NUC_REF="" # nuclear decoy/reference (FASTA) — small panel ok
MT_SEED="" # current mito seed/contigs (FASTA)

OUTDIR="mt-recruit"
THREADS=16
VERBOSE=0
QUIET=0
MODE=1                             # 0=safe, 1=balanced, 2=aggressive
DO_PT_ISOMERS=1                    # try to build two plastid forms when possible
DO_COV=0                           # run quick coverage check
ID_MIN="" QCOV_MIN="" TSPAN_MIN="" # tuned by MODE
DRY=0

print_help() {
	cat <<EOF
polap-bash-recruit-mt-main.sh ${VERSION}

USAGE
  bash polap-bash-recruit-mt-main.sh -r reads.fq.gz -m mito_seed.fa -p plastid.fa -n nuc.fa -o outdir [options]

REQUIRED
  -r, --reads PATH         Input reads (FASTQ/FASTA(.gz))
  -m, --mito PATH          Mito seed contigs (FASTA) — can be tiny/fragmented
  -p, --plastid PATH       Plastid reference (FASTA) for subtraction
  -n, --nuc-ref PATH       Nuclear reference/decoy (FASTA) for subtraction

OPTIONS
  -o, --outdir DIR         Output directory [${OUTDIR}]
  -t, --threads INT        Threads [${THREADS}]
  --mode 0|1|2             0=safe, 1=balanced, 2=aggressive [${MODE}]
  --no-pt-iso              Skip plastid two-isomer builder
  --cov-check              Run quick coverage check at the end
  --id-min F               Override identity minimum (competitive step)
  --qcov-min F             Override query coverage minimum
  --tspan-min INT          Override minimal on-target span (bp)
  --dry                    Print commands only (no execution)
  -v, --verbose            Verbose logging
  --quiet                  Silence info logs
  -h, --help               Show help

NOTES
  * Single rolling log: <outdir>/pipeline.log
  * Step TSVs are written side-by-side in <outdir>.
  * External scripts are called via _run/_sh to avoid quoting issues.
EOF
}

# Parse
while [[ $# -gt 0 ]]; do
	case "$1" in
	-r | --reads)
		READS="$2"
		shift 2
		;;
	-m | --mito)
		MT_SEED="$2"
		shift 2
		;;
	-p | --plastid)
		PT_REF="$2"
		shift 2
		;;
	-n | --nuc-ref)
		NUC_REF="$2"
		shift 2
		;;
	-o | --outdir)
		OUTDIR="$2"
		shift 2
		;;
	-t | --threads)
		THREADS="${2}"
		shift 2
		;;
	--mode)
		MODE="${2}"
		shift 2
		;;
	--no-pt-iso)
		DO_PT_ISOMERS=0
		shift
		;;
	--cov-check)
		DO_COV=1
		shift
		;;
	--id-min)
		ID_MIN="$2"
		shift 2
		;;
	--qcov-min)
		QCOV_MIN="$2"
		shift 2
		;;
	--tspan-min)
		TSPAN_MIN="$2"
		shift 2
		;;
	--dry)
		DRY=1
		shift
		;;
	-v | --verbose)
		VERBOSE=$((VERBOSE + 1))
		shift
		;;
	--quiet)
		QUIET=1
		VERBOSE=0
		shift
		;;
	-h | --help)
		print_help
		exit 0
		;;
	--version)
		echo "$VERSION"
		exit 0
		;;
	*)
		echo "[ERROR] unknown arg: $1" >&2
		exit 2
		;;
	esac
done

# sanity
[[ -z "$READS" || -z "$MT_SEED" || -z "$PT_REF" || -z "$NUC_REF" ]] && {
	print_help
	exit 2
}
mkdir -p "$OUTDIR"
_LOG="${OUTDIR}/pipeline.log"
_DR=$DRY

# tune thresholds by MODE (unless overridden)
case "$MODE" in
0) # safe
	: "${ID_MIN:=0.90}"
	: "${QCOV_MIN:=0.70}"
	: "${TSPAN_MIN:=1500}"
	;;
1) # balanced
	: "${ID_MIN:=0.85}"
	: "${QCOV_MIN:=0.60}"
	: "${TSPAN_MIN:=1200}"
	;;
2) # aggressive intergenic
	: "${ID_MIN:=0.80}"
	: "${QCOV_MIN:=0.50}"
	: "${TSPAN_MIN:=800}"
	;;
*)
	_log0 "[ERROR] invalid --mode $MODE (0|1|2)"
	exit 2
	;;
esac

_log0 "polap-bash-recruit-mt-main.sh ${VERSION}"
_log1 "reads=$READS"
_log1 "mt_seed=$MT_SEED"
_log1 "plastid=$PT_REF"
_log1 "nuc_ref=$NUC_REF"
_log1 "mode=$MODE id_min=$ID_MIN qcov_min=$QCOV_MIN tspan_min=$TSPAN_MIN"
_log1 "threads=$THREADS dry=$DRY"

# ──────────────────────────────────────────────────────────────────────────────
# Step 1: Plastid two-isomer builder (optional)
#   Produces two rotated/isomerized plastid forms (PT1/PT2) where possible.
#   If IR/SSC inference fails, the helper falls back to a doubled linear PT.
#   Outputs: <outdir>/pt_form1.fa, <outdir>/pt_form2.fa, step1.tsv
# ──────────────────────────────────────────────────────────────────────────────

PT1="${OUTDIR}/pt_form1.fa"
PT2="${OUTDIR}/pt_form2.fa"
if [[ $DO_PT_ISOMERS -eq 1 ]]; then
	_log0 "[1/6] plastid two-isomer builder"
	_run bash "$PT_TWOFORMS" \
		--input "$PT_REF" \
		--out1 "$PT1" \
		--out2 "$PT2" \
		--threads "$THREADS" \
		-o "$OUTDIR/step1.tsv"
else
	_log0 "[1/6] plastid two-isomer builder (skipped)"
	_run cp -f "$PT_REF" "$PT1"
	_run cp -f "$PT_REF" "$PT2"
	_sh "printf 'step\tpt_forms\tcomment\n1\tfallback\tno-pt-iso\n' > ${OUTDIR}/step1.tsv"
fi

# ──────────────────────────────────────────────────────────────────────────────
# Step 2: Plastid subtraction (competitive vs both PT isomers)
#   Removes obvious plastid reads. Writes:
#     step2.tsv, reads.step2.keep.fq(.gz)  (retained after subtraction)
# ──────────────────────────────────────────────────────────────────────────────

KEEP_PT_SUB="${OUTDIR}/reads.step2.keep.fq.gz"
_log0 "[2/6] plastid subtraction (competitive vs two PT forms)"
_run bash "$RUN_COMP" \
	--reads "$READS" \
	--mito /dev/null \
	--pt1 "$PT1" \
	--pt2 "$PT2" \
	--nuc /dev/null \
	--threads "$THREADS" \
	--id-min "$ID_MIN" \
	--qcov-min "$QCOV_MIN" \
	--tspan-min "$TSPAN_MIN" \
	--mode "$MODE" \
	-o "$OUTDIR/step2.tsv" \
	--emit-keep "$KEEP_PT_SUB"

# ──────────────────────────────────────────────────────────────────────────────
# Step 3: Nuclear subtraction (competitive vs nuc decoy/ref)
#   Removes obvious nuclear reads (BUSCO set or small decoy panel).
#   Input  : KEEP_PT_SUB
#   Output : reads.step3.keep.fq.gz, step3.tsv
# ──────────────────────────────────────────────────────────────────────────────

# step 3: nuclear subtraction
if [[ -n "${NUC_REF:-}" ]]; then
	_log1 "[3/6] nuclear subtraction using ${NUC_REF}"
	_run bash "$RUN_NUC_SUB" \
		--reads "$ROUND_DIR/step2.recruited.fq.gz" \
		--nuc-ref "$NUC_REF" \
		-o "$ROUND_DIR/step3"
else
	_log1 "[3/6] nuclear subtraction skipped (no --nuc-ref provided)"
	ln -s "$ROUND_DIR/step2.recruited.fq.gz" \
		"$ROUND_DIR/step3.recruited.fq.gz"
fi

KEEP_NUC_SUB="${OUTDIR}/reads.step3.keep.fq.gz"
_log0 "[3/6] nuclear subtraction (competitive vs nuc decoy/ref)"
_run bash "$RUN_COMP" \
	--reads "$KEEP_PT_SUB" \
	--mito /dev/null \
	--pt1 /dev/null \
	--pt2 /dev/null \
	--nuc "$NUC_REF" \
	--threads "$THREADS" \
	--id-min "$ID_MIN" \
	--qcov-min "$QCOV_MIN" \
	--tspan-min "$TSPAN_MIN" \
	--mode "$MODE" \
	-o "$OUTDIR/step3.tsv" \
	--emit-keep "$KEEP_NUC_SUB"

# ──────────────────────────────────────────────────────────────────────────────
# Step 4: Iterative k-mer signature recruitment (5% gain stop)
#   Seeds: MT_SEED; reads: KEEP_NUC_SUB
#   Updates the recruited set in-place; emits per-round TSV; final:
#     reads.step4.kmer.fq.gz, step4.kmer.tsv
# ──────────────────────────────────────────────────────────────────────────────

KMER_OUT="${OUTDIR}/reads.step4.kmer.fq.gz"
_log0 "[4/6] iterative k-mer signature recruitment (5% incremental stop)"
_run bash "$RUN_KSIG" \
	--reads "$KEEP_NUC_SUB" \
	--seed "$MT_SEED" \
	--threads "$THREADS" \
	--mode "$MODE" \
	-o "$OUTDIR/step4.kmer.tsv" \
	--emit "$KMER_OUT"

# ──────────────────────────────────────────────────────────────────────────────
# Step 5: Iterative read-overlap recruitment (5% gain stop)
#   Starts from KMER_OUT and grows by read-to-read overlaps
#   Final: reads.step5.ovlp.fq.gz, step5.ovlp.tsv
# ──────────────────────────────────────────────────────────────────────────────

OVLP_OUT="${OUTDIR}/reads.step5.ovlp.fq.gz"
_log0 "[5/6] iterative read-overlap recruitment (5% incremental stop)"
_run bash "$RUN_OVLP" \
	--reads "$KMER_OUT" \
	--threads "$THREADS" \
	--mode "$MODE" \
	-o "$OUTDIR/step5.ovlp.tsv" \
	--emit "$OVLP_OUT"

# ──────────────────────────────────────────────────────────────────────────────
# Step 6: (Optional) Quick coverage check on current seed or a provisional assembly
#   Uses RUN_COVQ; you can pass a contig set with --mito (MT_SEED) or
#   replace with your latest draft assembly.
# ──────────────────────────────────────────────────────────────────────────────

if [[ $DO_COV -eq 1 ]]; then
	_log0 "[6/6] quick coverage check (reads vs mito seed)"
	_run bash "$RUN_COVQ" \
		--reads "$OVLP_OUT" \
		--contig "$MT_SEED" \
		--threads "$THREADS" \
		-o "$OUTDIR/step6.cov.tsv"
else
	_log0 "[6/6] quick coverage check (skipped)"
fi

_log0 "DONE  ${VERSION}"
