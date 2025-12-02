#!/usr/bin/env bash
# File: polap-bash-polish-racon-polypolish.sh
# Version: v0.9.2
# MTPT-aware organelle polishing with failsafe step wrappers.
set -Eeuo pipefail

# ---------- Locate polaplib and source run lib ----------
POLAPLIB_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd -P)"
RUNLIB="${POLAPLIB_DIR}/polap-lib-run-command.sh"
[[ -s "$RUNLIB" ]] || {
	echo "Missing: $RUNLIB" >&2
	exit 2
}
# shellcheck source=polap-lib-run-command.sh
source "$RUNLIB"

# polap_run_wrapper() {
# 	local name="$1"
# 	shift
# 	local stepdir="$OUTDIR/.polap-run/$name"
# 	local steplog="$stepdir/step.log"
# 	polap__ensure_dir "$stepdir"
#
# 	if ((POLAP_FORCE == 0)) && [[ -e "$stepdir/.ok" ]]; then
# 		polap_log1 "SKIP (ok): $name"
# 		return 0
# 	fi
#
# 	printf '%q ' "$@" >"$stepdir/cmd.txt"
# 	: >"$stepdir/.running"
# 	polap_log1 "RUN : $name"
#
# 	if ((POLAP_DRYRUN == 1)); then
# 		rm -f "$stepdir/.running"
# 		return 0
# 	fi
#
# 	local rc=0
# 	(
# 		set -Eeuo pipefail
# 		"$@"
# 	) >>"$steplog" 2>&1 || rc=$?
#
# 	if ((rc != 0)); then
# 		rm -f "$stepdir/.running"
# 		polap_log0 "FAIL: $name (rc=$rc) see $steplog"
# 		return "$rc"
# 	fi
#
# 	touch "$stepdir/.ok"
# 	rm -f "$stepdir/.running"
# 	polap_log1 "DONE: $name"
# }

# ---------- Helpers & step script paths ----------
SCRIPTS_DIR="${POLAPLIB_DIR}/scripts"

PAF_FILTER="${SCRIPTS_DIR}/polap_py_paf_filter.py"
MASK_FROM_PAF="${SCRIPTS_DIR}/polap_py_mask_from_paf.py"
BED_COMPLEMENT="${SCRIPTS_DIR}/polap_py_bed_complement.py"
PAF_TO_HITS="${SCRIPTS_DIR}/polap_py_paf_to_hits_tsv.py"
BED_ADD_LEN="${SCRIPTS_DIR}/polap_py_bed_add_len.py"
ASSIGN_ONT_IDS="${SCRIPTS_DIR}/polap_py_assign_ont_ids.py"
ASSIGN_SR_IDS="${SCRIPTS_DIR}/polap_py_assign_sr_ids.py"
FASTQ_FILTER_BY_IDS="${SCRIPTS_DIR}/polap_py_fastq_filter_by_ids.py"
COV_BINS="${SCRIPTS_DIR}/polap_py_bam_coverage_bins.py"
COV_PLOT="${SCRIPTS_DIR}/polap_r_covplot_mtpt_overlay.R"
MASK_QC_R="${SCRIPTS_DIR}/polap_r_mask_qc.R"
MAKE_COMBINED="${SCRIPTS_DIR}/polap_py_make_combined_fasta.py"

STEP1_ASSIGN_ONT="${POLAPLIB_DIR}/polap-bash-polish-raconpoly-step1-assign-ont.sh"
# STEP2_FILTER_FQ="${POLAPLIB_DIR}/polap-bash-polish-raconpoly-step2-filter-fastq-by-ids.sh"
STEP3_ASSIGN_SR="${POLAPLIB_DIR}/polap-bash-polish-raconpoly-step3-assign-sr.sh"
STEP4_RACON_RND="${POLAPLIB_DIR}/polap-bash-polish-raconpoly-step4-racon-round.sh"
STEP5_MEDAKA="${POLAPLIB_DIR}/polap-bash-polish-raconpoly-step5-medaka.sh"
STEP6_A2A_MASK="${POLAPLIB_DIR}/polap-bash-polish-raconpoly-step6-a2a-mask.sh"
STEP7_COVERAGE="${POLAPLIB_DIR}/polap-bash-polish-raconpoly-step7-coverage.sh"
STEP8_POLYPOL="${POLAPLIB_DIR}/polap-bash-polish-raconpoly-step8-polypolish.sh"

# ---------- Defaults ----------
TARGET="mt" # mt | cp
FASTA=""
OTHER=""
ONT=""
SR1=""
SR2=""
OUTDIR=""
THREADS=16
RACON_ROUNDS=2
ONT_PRESET="map-ont"
PAF_MIN_IDENT=0.91
PAF_MIN_ALEN=3000
USE_MEDAKA=0
MEDAKA_MODEL="r941_min_sup_g507"
MASK_POLICY="auto"
MASK_MIN_IDENT=0.85
MASK_MIN_LEN=150
MASK_PAD=1000
RUN_MASK_QC=1
COV_DISABLE=0
COV_BIN=200
COV_MIN_MAPQ=0
ASSIGN_MAPQ=30
PP_ALIGNER="bowtie2" # bowtie2 | bwa | minimap2
: "${PP_MMSR_SECONDARIES:=0}"
: "${POLAP_FORCE:=0}"
: "${POLAP_DRYRUN:=0}"

usage() {
	cat <<EOF
Usage:
  $(basename "$0") --target mt|cp --fasta asm.fa --other other.fa \\
    [--ont ont.fq.gz] [--sr1 R1.fq.gz --sr2 R2.fq.gz] \\
    --outdir OUT [--threads 16] \\
    [--racon-rounds 2] [--ont-preset map-ont] [--use-medaka 0|1] [--medaka-model MODEL] \\
    [--paf-min-ident 0.91] [--paf-min-alen 3000] \\
    [--mask auto|on|off] [--mask-min-ident 0.85] [--mask-min-len 150] [--mask-pad 1000] [--mask-qc 0|1] \\
    [--cov-disable 0|1] [--cov-bin 200] [--cov-min-mapq 0] \\
    [--assign-mapq 30] \\
    [--pp-aligner bowtie2|bwa|minimap2] [--pp-mm2-secondary 0|1] \\
    [--force 0|1] [--dry-run 0|1]
EOF
}

# ---------- Parse CLI ----------
while (($#)); do
	case "$1" in
	--target)
		TARGET="${2:?}"
		shift 2
		;;
	--fasta)
		FASTA="${2:?}"
		shift 2
		;;
	--other)
		OTHER="${2:-}"
		shift 2
		;;
	--ont)
		ONT="${2:-}"
		shift 2
		;;
	--sr1)
		SR1="${2:-}"
		shift 2
		;;
	--sr2)
		SR2="${2:-}"
		shift 2
		;;
	--outdir)
		OUTDIR="${2:?}"
		shift 2
		;;
	--threads)
		THREADS="${2:?}"
		shift 2
		;;
	--racon-rounds)
		RACON_ROUNDS="${2:?}"
		shift 2
		;;
	--ont-preset)
		ONT_PRESET="${2:?}"
		shift 2
		;;
	--paf-min-ident)
		PAF_MIN_IDENT="${2:?}"
		shift 2
		;;
	--paf-min-alen)
		PAF_MIN_ALEN="${2:?}"
		shift 2
		;;
	--use-medaka)
		USE_MEDAKA="${2:?}"
		shift 2
		;;
	--medaka-model)
		MEDAKA_MODEL="${2:?}"
		shift 2
		;;
	--mask)
		MASK_POLICY="${2:?}"
		shift 2
		;;
	--mask-min-ident)
		MASK_MIN_IDENT="${2:?}"
		shift 2
		;;
	--mask-min-len)
		MASK_MIN_LEN="${2:?}"
		shift 2
		;;
	--mask-pad)
		MASK_PAD="${2:?}"
		shift 2
		;;
	--mask-qc)
		RUN_MASK_QC="${2:?}"
		shift 2
		;;
	--cov-disable)
		COV_DISABLE="${2:?}"
		shift 2
		;;
	--cov-bin)
		COV_BIN="${2:?}"
		shift 2
		;;
	--cov-min-mapq)
		COV_MIN_MAPQ="${2:?}"
		shift 2
		;;
	--assign-mapq)
		ASSIGN_MAPQ="${2:?}"
		shift 2
		;;
	--pp-aligner)
		PP_ALIGNER="${2:?}"
		shift 2
		;;
	--pp-mm2-secondary)
		PP_MMSR_SECONDARIES="${2:?}"
		shift 2
		;;
	--force)
		POLAP_FORCE="${2:?}"
		shift 2
		;;
	--dry-run)
		POLAP_DRYRUN="${2:?}"
		shift 2
		;;
	-h | --help)
		usage
		exit 0
		;;
	*) polap_die "Unknown option: $1" ;;
	esac
done

[[ -n "$FASTA" && -n "$OUTDIR" ]] || {
	usage
	polap_die "Required: --fasta and --outdir"
}

[[ "$TARGET" == "mt" || "$TARGET" == "cp" ]] || polap_die "--target must be mt or cp"

[[ -s "$FASTA" ]] || polap_die "FASTA not found: $FASTA"

[[ -n "$OTHER" ]] && [[ ! -s "$OTHER" ]] && polap_die "Other FASTA not found: $OTHER"

# ---------- Prepare dirs, logging, lock ----------
install -d -- "$OUTDIR"/{stage0_assign,stage1_racon,mask,stage3_polypolish,coverage,final,logs}
LOG_FILE="$OUTDIR/logs/main.log"
LOCKDIR="$OUTDIR/.lock"
polap_lock_acquire_wrapper "$LOCKDIR" || polap_die "Another process is using $OUTDIR (lock busy)."

# ---------- deps ----------
polap_require samtools minimap2 racon polypolish python3

case "$PP_ALIGNER" in
bowtie2) polap_require bowtie2 ;;
bwa) polap_require bwa-mem2 ;;
minimap2) : ;;
*) polap_die "--pp-aligner must be bowtie2|bwa|minimap2" ;;
esac

for f in "$PAF_FILTER" "$MASK_FROM_PAF" "$BED_COMPLEMENT" "$PAF_TO_HITS" "$BED_ADD_LEN" \
	"$ASSIGN_ONT_IDS" "$ASSIGN_SR_IDS" "$FASTQ_FILTER_BY_IDS" "$COV_BINS" "$MAKE_COMBINED"; do
	if [[ ! -s "$f" ]]; then
		polap_die "Helper missing: $f"
	fi
done

python3 - <<'PY' || polap_die "Python module 'pysam' not found (pip install pysam)"
import importlib, sys; importlib.import_module("pysam"); sys.exit(0)
PY

# ---------- Combined reference (Python; no AWK) ----------
ASSIGN_DIR="$OUTDIR/stage0_assign"
COMBINED="$ASSIGN_DIR/combined.fa"
TGT_PREFIX="${TARGET}|"
OTH_PREFIX="$([[ "$TARGET" == "mt" ]] && echo "cp|" || echo "mt|")"

POLAP_RUN_BACKEND="direct"
if ! polap_run_wrapper "s0_combined_fa" \
	python3 "$MAKE_COMBINED" \
	--target-fasta "$FASTA" --target-prefix "$TGT_PREFIX" \
	${OTHER:+--other-fasta "$OTHER"} ${OTHER:+--other-prefix "$OTH_PREFIX"} \
	--out "$COMBINED"; then
	polap_log0 "[s0_combined_fa] FAILED"
	polap_lock_release_wrapper "$LOCKDIR"
	exit 1
fi

POLAP_RUN_BACKEND=""
if ! polap_run_wrapper "s0_combined_fai" \
	samtools faidx "$COMBINED"; then
	polap_log0 "[s0_combined_fai] FAILED"
	polap_lock_release_wrapper "$LOCKDIR"
	exit 1
fi

# ---------- Stage 0: assignment ----------
LR_ASSIGNED="$ASSIGN_DIR/ont.target.fq.gz"
SR1_ASSIGNED="$ASSIGN_DIR/sr.target.R1.fq.gz"
SR2_ASSIGNED="$ASSIGN_DIR/sr.target.R2.fq.gz"
HAVE_LR=0
HAVE_SR=0

# ONT
IDS_ONT="$ASSIGN_DIR/ont.target.ids"
if [[ -n "$ONT" ]]; then

	POLAP_RUN_BACKEND="direct"
	if ! polap_run_wrapper "s0_assign_ont_ids" \
		bash "$STEP1_ASSIGN_ONT" -r "$COMBINED" -q "$ONT" -x "$TGT_PREFIX" \
		-m "$ASSIGN_MAPQ" -p "$ONT_PRESET" -t "$THREADS" -o "$IDS_ONT"; then
		polap_log0 "[s0_assign_ont_ids] FAILED"
		polap_lock_release_wrapper "$LOCKDIR"
		exit 1
	fi

	POLAP_RUN_BACKEND=""
	if ! [[ -s "$IDS_ONT" ]]; then
		polap_log0 "No ONT reads assigned; adjust --assign-mapq or inputs"
		polap_lock_release_wrapper "$LOCKDIR"
		exit 1
	fi

	if ! polap_run_wrapper "s0_filter_ont_to_target" \
		seqkit grep -f "$IDS_ONT" "$ONT" -o "$LR_ASSIGNED"; then
		polap_log0 "[s0_filter_ont_to_target] FAILED"
		polap_lock_release_wrapper "$LOCKDIR"
		exit 1
	fi

	# if ! polap_run_wrapper "s0_filter_ont_to_target" \
	# 	bash "$STEP2_FILTER_FQ" -i "$ONT" -l "$IDS_ONT" -o "$LR_ASSIGNED"; then
	# 	polap_log0 "[s0_filter_ont_to_target] FAILED"
	# 	polap_lock_release_wrapper "$LOCKDIR"
	# 	exit 1
	# fi

	HAVE_LR=1
else
	polap_log1 "No ONT provided; LR polish & coverage will be skipped"
fi

# SR pairs
IDS_SR="$ASSIGN_DIR/sr.target.ids"
if [[ -n "$SR1" && -n "$SR2" ]]; then
	POLAP_RUN_BACKEND="direct"
	if ! polap_run_wrapper "s0_assign_sr_ids" \
		bash "$STEP3_ASSIGN_SR" -r "$COMBINED" -x "$TGT_PREFIX" -m "$ASSIGN_MAPQ" -t "$THREADS" \
		-1 "$SR1" -2 "$SR2" -o "$IDS_SR"; then
		polap_log0 "[s0_assign_sr_ids] FAILED"
		polap_lock_release_wrapper "$LOCKDIR"
		exit 1
	fi

	POLAP_RUN_BACKEND=""
	if [[ -s "$IDS_SR" ]]; then
		# Filter R1 reads by assigned IDs
		if ! polap_run_wrapper "s0_filter_sr_R1" \
			seqkit grep -f "$IDS_SR" "$SR1" -o "$SR1_ASSIGNED"; then
			polap_log0 "[s0_filter_sr_R1] FAILED"
			polap_lock_release_wrapper "$LOCKDIR"
			exit 1
		fi

		# Filter R2 reads by assigned IDs
		if ! polap_run_wrapper "s0_filter_sr_R2" \
			seqkit grep -f "$IDS_SR" "$SR2" -o "$SR2_ASSIGNED"; then
			polap_log0 "[s0_filter_sr_R2] FAILED"
			polap_lock_release_wrapper "$LOCKDIR"
			exit 1
		fi

		HAVE_SR=1
	else
		polap_log1 "WARNING: no short-read pairs confidently assigned; SR polish will be skipped."
	fi

	# if [[ -s "$IDS_SR" ]]; then
	# 	if ! polap_run_wrapper "s0_filter_sr_R1" bash "$STEP2_FILTER_FQ" -i "$SR1" -l "$IDS_SR" -o "$SR1_ASSIGNED"; then
	# 		polap_log0 "[s0_filter_sr_R1] FAILED"
	# 		polap_lock_release_wrapper "$LOCKDIR"
	# 		exit 1
	# 	fi
	# 	if ! polap_run_wrapper "s0_filter_sr_R2" bash "$STEP2_FILTER_FQ" -i "$SR2" -l "$IDS_SR" -o "$SR2_ASSIGNED"; then
	# 		polap_log0 "[s0_filter_sr_R2] FAILED"
	# 		polap_lock_release_wrapper "$LOCKDIR"
	# 		exit 1
	# 	fi
	# 	HAVE_SR=1
	# else
	# 	polap_log1 "WARNING: no short-read pairs confidently assigned; SR polish will be skipped."
	# fi

else
	polap_log1 "No SR provided; SR polish will be skipped"
fi

((HAVE_LR == 1 || HAVE_SR == 1)) || {
	polap_log0 "No reads assigned; nothing to polish"
	polap_lock_release_wrapper "$LOCKDIR"
	exit 1
}

# ---------- Stage 1: Racon xN, optional Medaka ----------
CUR="$FASTA"
if ((HAVE_LR == 1)); then
	polap_log1 "Stage 1: racon x${RACON_ROUNDS} (preset=${ONT_PRESET})"
	for r in $(seq 1 "$RACON_ROUNDS"); do
		OUTFA="$OUTDIR/stage1_racon/polished.r${r}.fa"
		if ! polap_run_wrapper "s1_r${r}_racon_round" \
			bash "$STEP4_RACON_RND" -R "$LR_ASSIGNED" -A "$CUR" -p "$ONT_PRESET" -t "$THREADS" \
			-i "$PAF_MIN_IDENT" -a "$PAF_MIN_ALEN" -n "$r" -O "$OUTFA" -D "$OUTDIR/stage1_racon"; then
			polap_log0 "[s1_r${r}_racon_round] FAILED"
			polap_lock_release_wrapper "$LOCKDIR"
			exit 1
		fi

		# Index the new consensus so any downstream tool (R coverage plot) can read .fai
		if ! polap_run_wrapper "s1_r${r}_faidx" \
			samtools faidx "$OUTFA"; then
			polap_log0 "[s1_r${r}_faidx] FAILED"
			polap_lock_release_wrapper "$LOCKDIR"
			exit 1
		fi

		CUR="$OUTFA"
	done

	# We do not use medaka because we do not know the exact sequencing platform.
	#
	# if [[ "$USE_MEDAKA" -eq 1 ]] && command -v medaka_consensus >/dev/null 2>&1; then
	# 	MED_DIR="$OUTDIR/stage1_racon/medaka"
	# 	if ! polap_run_wrapper "s1_medaka" \
	# 		bash "$STEP5_MEDAKA" -R "$LR_ASSIGNED" -A "$CUR" -m "$MEDAKA_MODEL" -t "$THREADS" -o "$MED_DIR"; then polap_log0 "[s1_medaka] FAILED (continue with racon)"; fi
	# 	[[ -s "$MED_DIR/consensus.fasta" ]] && CUR="$MED_DIR/consensus.fasta"
	# else
	# 	polap_log1 "Medaka OFF or not found"
	# fi
fi

# ---------- Stage 2: Mask detect + reports ----------

# 3 cases but only 2 cases.
# 1. No mask
# 2. Detect MTPT candidates
# 3. or we could use detected MTPT.
#
[[ "$MASK_POLICY" == "auto" ]] && { if [[ "$TARGET" == "mt" ]]; then MASK_POLICY="on"; else MASK_POLICY="off"; fi; }
MASK_DIR="$OUTDIR/mask"
install -d -- "$MASK_DIR"
MASK_BED="$MASK_DIR/transfer.mask.bed"
ALLOW_BED="$MASK_DIR/allow.bed"

if [[ "$MASK_POLICY" == "off" || -z "$OTHER" || ! -s "$OTHER" ]]; then

	polap_log1 "Stage 2: mask OFF or --other missing; ALLOW = whole assembly"
	if ! polap_run_wrapper "s2_allow_whole" \
		python3 "$BED_COMPLEMENT" --fasta "$CUR" --bed /dev/null --out-bed "$ALLOW_BED"; then
		polap_log0 "[s2_allow_whole] FAILED"
		polap_lock_release_wrapper "$LOCKDIR"
		exit 1
	fi
	: >"$MASK_BED"

else

	if ! polap_run_wrapper "s2_mask_detect" \
		bash "$STEP6_A2A_MASK" -A "$CUR" -O "$OTHER" -t "$THREADS" \
		-i "$MASK_MIN_IDENT" -l "$MASK_MIN_LEN" -p "$MASK_PAD" -d "$MASK_DIR" -Q "$RUN_MASK_QC"; then
		polap_log0 "[s2_mask_detect] FAILED"
		polap_lock_release_wrapper "$LOCKDIR"
		exit 1
	fi

fi

# ---------- Stage 2.5: Coverage ----------
if ((COV_DISABLE == 0)) && [[ -n "$ONT" ]] && ((HAVE_LR == 1)); then

	# Ensure the assembly used for coverage/overlay has a .fai
	if ! polap_run_wrapper "s25_cur_faidx" \
		samtools faidx "$CUR"; then
		polap_log0 "[s25_cur_faidx] FAILED"
		polap_lock_release_wrapper "$LOCKDIR"
		exit 1
	fi

	if ! polap_run_wrapper "s25_coverage" \
		bash "$STEP7_COVERAGE" -A "$CUR" -Q "$ONT" -E "$ASSIGN_DIR/ont.target.fq.gz" \
		-t "$THREADS" -b "$COV_BIN" -q "$COV_MIN_MAPQ" -m "$MASK_BED" -o "$OUTDIR/coverage" -R "$COV_PLOT"; then
		polap_log0 "[s25_coverage] FAILED"
		polap_lock_release_wrapper "$LOCKDIR"
		exit 1
	fi
else
	polap_log1 "Stage 2.5: coverage skipped"
fi

# ---------- Stage 3: Polypolish ----------
FINAL="$CUR"
PP_DIR="$OUTDIR/stage3_polypolish"
install -d -- "$PP_DIR"
if [[ -s "${SR1_ASSIGNED:-/dev/null}" && -s "${SR2_ASSIGNED:-/dev/null}" ]]; then
	if ! polap_run_wrapper "s3_polypolish" \
		bash "$STEP8_POLYPOL" -a "$PP_ALIGNER" -S "$PP_MMSR_SECONDARIES" \
		-A "$CUR" -1 "$SR1_ASSIGNED" -2 "$SR2_ASSIGNED" -B "$ALLOW_BED" \
		-t "$THREADS" -d "$PP_DIR" -o "$PP_DIR/polished.fa" -b "$PP_DIR/short.name.bam"; then polap_log0 "[s3_polypolish] FAILED (continue)"; else
		[[ -s "$PP_DIR/polished.fa" ]] && FINAL="$PP_DIR/polished.fa"
	fi
else
	polap_log1 "Stage 3: SR polish skipped (no confidently assigned pairs)"
fi

# ---------- Finalize ----------
if ! polap_run_wrapper "s4_final_copy" cp -f "$FINAL" "$OUTDIR/polished.fa"; then
	polap_log0 "[s4_final_copy] FAILED"
	polap_lock_release_wrapper "$LOCKDIR"
	exit 1
fi

polap_log1 "DONE: $OUTDIR/polished.fa"
polap_log1 "Summary: TARGET=$TARGET LR=$([[ -n "$ONT" ]] && echo yes || echo no) SR=$([[ -n "$SR1" && -n "$SR2" ]] && echo yes || echo no)"
polap_log1 "Mask: policy=$MASK_POLICY thresholds: ident>=$MASK_MIN_IDENT len>=$MASK_MIN_LEN pad=${MASK_PAD}bp"
polap_log1 "Coverage: bin=$COV_BIN minMAPQ=$COV_MIN_MAPQ status=$([[ $COV_DISABLE -eq 0 && -n $ONT ]] && echo on || echo off)"
polap_log1 "Polypolish aligner: $PP_ALIGNER (mm2 secondaries=${PP_MMSR_SECONDARIES})"

polap_lock_release_wrapper "$LOCKDIR"
