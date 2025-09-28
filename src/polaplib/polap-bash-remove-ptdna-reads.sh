#!/usr/bin/env bash
# polap-bash-remove-ptdna-reads.sh
# Version: v0.5.0
#
# Purpose:
#   1) Remove plastid-origin reads (PT) from ONT long reads using a
#      data-driven identity cutoff (or --identity-min override).
#   2) [--nuclear-deplete] On PT-depleted reads, run one all-vs-all
#      minimap2 and remove nuclear-origin reads based on overlapness
#      thresholds guided by BUSCO nuclear labels (miniprot).
#   3) [--miniprot-sample] Do BUSCO labeling on a length-stratified
#      sample (fast) and apply thresholds to the full set.
#
# Inputs:
#   -r, --reads            FASTQ(.gz) reads            [required]
#   -p, --pt-ref           plastid reference FASTA     [required]
#
# PT removal knobs:
#   --pt-origin            PT-origin ids (optional, helps identity cut)
#   --nuc-origin           Nuclear-origin ids (optional, helps identity cut)
#   --alen-min             aligned length guard for PT selection [1000]
#   --identity-min         override auto identity cutoff (e.g., 0.93)
#   --fpr                  nuclear FPR target for identity [0.01]
#   --tpr                  PT TPR target for identity [0.95]
#
# Nuclear deplete (overlapness):
#   --nuclear-deplete      run all-vs-all + overlapness nuclear filtering
#   --busco                BUSCO proteins FAA (required with --nuclear-deplete)
#   --minimap2-ava-opts    extra options for all-vs-all minimap2
#   --min-olen             graph edge min aligned length [1200]
#   --min-ident            graph edge min identity       [0.84]
#   --w-floor              graph edge min weight         [0.12]
#   --fpr-nuc              FPR target for nuclear tail   [0.01]
#
# Miniprot sampling (faster BUSCO labeling):
#   --miniprot-sample      enable length-stratified sampling
#   --miniprot-sample-n    total reads to sample         [25000]
#   --miniprot-sample-seed RNG seed                      [42]
#   --miniprot-sample-bins length cutpoints (bp)         [2000,5000,10000]
#
# General:
#   -o, --outdir           output directory              [ptfilter.out]
#   -t, --threads          threads                       [16]
#   --dry                  log only; no execution
#   -v, --verbose          verbose notes
#   --quiet                quiet
#   --profiling            write per-step timings to log/profile.tsv
#
# Outputs:
#   out/reads.nonpt.fq.gz                 (PT removed)
#   out/reads.mt_enriched.fq.gz           (if --nuclear-deplete)
#   out/pt.ids, out/pt_thresh.vars, out/pt_thresh.diag.tsv
#   out/01_allvsall/allvsall.paf, out/overlapness.tsv, out/nuc.ids
#
set -euo pipefail

# ───────────── Logging & profiling ─────────────
VERBOSE=1
QUIET=0
DRY=0
LOG_FILE=""
PROF=0
PROF_FILE=""

note() {
	[[ $QUIET -eq 1 ]] && return 0
	local ts src ln
	ts="$(date +'%F %T')"
	src="${BASH_SOURCE[1]##*/}"
	ln="${BASH_LINENO[0]}"
	printf "[%s][%s:%s] %s\n" "$ts" "$src" "$ln" "$*" |
		{ if [[ -n "$LOG_FILE" ]]; then tee -a "$LOG_FILE"; else cat; fi; } >&2
}
_prof_t0=""
_prof_step=""
prof_start() {
	((PROF)) || return 0
	_prof_step="$1"
	_prof_t0=$(date +%s)
}
prof_end() {
	((PROF)) || return 0
	local t1
	t1=$(date +%s)
	local dt=$((t1 - _prof_t0))
	echo -e "${_prof_step}\t${dt}" >>"$PROF_FILE"
}

# ───────────── CLI ─────────────
READS=""
PT_REF=""
OUTDIR="ptfilter.out"
THREADS=16
PT_ORIGIN=""
NUC_ORIGIN=""
ALEN_MIN=1000
IDENTITY_MIN=""
FPR=0.01
TPR=0.95

NUCLEAR_DEPLETE=0
BUSCO=""
AVA_OPTS=""
MIN_OLEN=1200
MIN_IDENT=0.84
W_FLOOR=0.12
FPR_NUC=0.01

MP_SAMPLE=0
MP_SAMPLE_N=25000
MP_SAMPLE_SEED=42
MP_SAMPLE_BINS="2000,5000,10000"

usage() {
	cat <<EOF
Usage:
  $0 -r reads.fq.gz -p plastid.fa -o outdir [options]
  Try: $0 -h
EOF
}

while [[ $# -gt 0 ]]; do
	case "$1" in
	-r | --reads)
		READS="$2"
		shift 2
		;;
	-p | --pt-ref)
		PT_REF="$2"
		shift 2
		;;
	-o | --outdir)
		OUTDIR="$2"
		shift 2
		;;
	-t | --threads)
		THREADS="$2"
		shift 2
		;;
	--pt-origin)
		PT_ORIGIN="$2"
		shift 2
		;;
	--nuc-origin)
		NUC_ORIGIN="$2"
		shift 2
		;;
	--alen-min)
		ALEN_MIN="$2"
		shift 2
		;;
	--identity-min)
		IDENTITY_MIN="$2"
		shift 2
		;;
	--fpr)
		FPR="$2"
		shift 2
		;;
	--tpr)
		TPR="$2"
		shift 2
		;;
	--nuclear-deplete)
		NUCLEAR_DEPLETE=1
		shift
		;;
	--busco)
		BUSCO="$2"
		shift 2
		;;
	--minimap2-ava-opts)
		AVA_OPTS="$2"
		shift 2
		;;
	--min-olen)
		MIN_OLEN="$2"
		shift 2
		;;
	--min-ident)
		MIN_IDENT="$2"
		shift 2
		;;
	--w-floor)
		W_FLOOR="$2"
		shift 2
		;;
	--fpr-nuc)
		FPR_NUC="$2"
		shift 2
		;;
	--miniprot-sample)
		MP_SAMPLE=1
		shift
		;;
	--miniprot-sample-n)
		MP_SAMPLE_N="$2"
		shift 2
		;;
	--miniprot-sample-seed)
		MP_SAMPLE_SEED="$2"
		shift 2
		;;
	--miniprot-sample-bins)
		MP_SAMPLE_BINS="$2"
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
		shift
		;;
	--profiling)
		PROF=1
		shift
		;;
	-h | --help)
		usage
		exit 0
		;;
	*)
		note "ERR unknown arg: $1"
		usage
		exit 2
		;;
	esac
done

[[ -s "$READS" ]] || {
	note "ERR missing --reads"
	exit 2
}
[[ -s "$PT_REF" ]] || {
	note "ERR missing --pt-ref"
	exit 2
}
mkdir -p "$OUTDIR"/{log,panel,map,01_allvsall}
LOG_FILE="$OUTDIR/log/remove_ptdna.log"
PROF_FILE="$OUTDIR/log/profile.tsv"
((PROF)) && echo -e "step\tseconds" >"$PROF_FILE"
((DRY)) && note "--dry on (no execution)"

# ───────────── Helpers & tool paths ─────────────
_POLAPLIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HELP_REMOVE="${_POLAPLIB_DIR}/remove-ptdna-reads"
HELP_FAST="${_POLAPLIB_DIR}/fast-mtseed-ont"

PY_PT_THRESH="${HELP_REMOVE}/polap-py-pt-ident-threshold.py"
PY_OVERLAP="${HELP_FAST}/polap-py-overlapness-from-paf.py"
PY_NUC_THRESH="${HELP_FAST}/polap-py-threshold-from-nuclear.py"
PT_ISOFORM_SH="${_POLAPLIB_DIR}/polap-bash-pt-isoform.sh"

need() { command -v "$1" >/dev/null 2>&1 || {
	note "ERR tool not in PATH: $1"
	exit 127
}; }
need minimap2
need seqkit
need awk
need sort
need comm
need python3
need samtools
need miniprot
for f in "$PY_PT_THRESH" "$PY_OVERLAP" "$PY_NUC_THRESH" "$PT_ISOFORM_SH"; do
	[[ -s "$f" ]] || {
		note "ERR missing helper: $f"
		exit 127
	}
done

# ───────────── 1) Build PT isoforms & map reads to doubled A/B ─────────────
PANEL_DIR="$OUTDIR/panel"
mkdir -p "$PANEL_DIR"
prof_start "pt_isoform"
note "Build plastid isoforms via $PT_ISOFORM_SH"
((!DRY)) && bash "$PT_ISOFORM_SH" -r "$PT_REF" -o "$PANEL_DIR" -t "$THREADS"
prof_end

ISO_A="$PANEL_DIR/pt_isomerA.fa"
ISO_B="$PANEL_DIR/pt_isomerB.fa"
[[ -s "$ISO_A" ]] || {
	note "ERR missing pt_isomerA.fa"
	exit 3
}

DBL_A="$PANEL_DIR/pt_isomerA.double.fa"
DBL_B="$PANEL_DIR/pt_isomerB.double.fa"
build_double() {
	local in="$1" out="$2"
	python - "$in" "$out" <<'PY'
import sys
inp,outp=sys.argv[1],sys.argv[2]
name=None; seq=[]
with open(inp) as f:
  for ln in f:
    if ln.startswith('>'):
      if name is None: name=ln[1:].strip().split()[0]
      else: break
    else: seq.append(ln.strip())
S=''.join(seq)
with open(outp,'w') as w:
  w.write('>pt.double\n'); w.write(S+S+'\n')
PY
}
prof_start "pt_double"
((!DRY)) && build_double "$ISO_A" "$DBL_A"
((!DRY)) && { [[ -s "$ISO_B" ]] && build_double "$ISO_B" "$DBL_B" || true; }
prof_end

PAF_A="$OUTDIR/map/formA.paf"
PAF_B="$OUTDIR/map/formB.paf"
prof_start "map_ptA"
note "Map reads → isomerA.double > $PAF_A"
((!DRY)) && minimap2 -x map-ont --secondary=yes -N 50 -t "$THREADS" \
	"$DBL_A" "$READS" >"$PAF_A"
prof_end
if [[ -s "$DBL_B" ]]; then
	prof_start "map_ptB"
	note "Map reads → isomerB.double > $PAF_B"
	((!DRY)) && minimap2 -x map-ont --secondary=yes -N 50 -t "$THREADS" \
		"$DBL_B" "$READS" >"$PAF_B"
	prof_end
fi

# ───────────── 2) Choose identity cutoff (or override) & emit pt.ids ─────────────
PT_VARS="$OUTDIR/pt_thresh.vars"
PT_DIAG="$OUTDIR/pt_thresh.diag.tsv"
PT_IDS="$OUTDIR/pt.ids"
[[ -n "$PT_ORIGIN" && -s "$PT_ORIGIN" ]] || note "WARN no --pt-origin; PT labels empty"
[[ -n "$NUC_ORIGIN" && -s "$NUC_ORIGIN" ]] || note "WARN no --nuc-origin; nuclear labels empty"

if [[ -n "$IDENTITY_MIN" ]]; then
	# Manual override: use AWK on PAFs with ident>=IDENTITY_MIN and alen>=ALEN_MIN
	note "Using user-provided --identity-min=$IDENTITY_MIN"
	prof_start "emit_pt_ids_override"
	if ((!DRY)); then
		awk -v ID="$IDENTITY_MIN" -v AL="$ALEN_MIN" 'BEGIN{FS=OFS="\t"}
      NF>=12{ id=($11>0?$10/$11:0); if(id>=ID && $11+0>=AL) print $1 }' \
			"$PAF_A" | sort -u >"$OUTDIR/ptA.ids"
		if [[ -s "$PAF_B" ]]; then
			awk -v ID="$IDENTITY_MIN" -v AL="$ALEN_MIN" 'BEGIN{FS=OFS="\t"}
        NF>=12{ id=($11>0?$10/$11:0); if(id>=ID && $11+0>=AL) print $1 }' \
				"$PAF_B" | sort -u >"$OUTDIR/ptB.ids"
			sort -u "$OUTDIR/ptA.ids" "$OUTDIR/ptB.ids" >"$PT_IDS"
		else
			mv "$OUTDIR/ptA.ids" "$PT_IDS"
		fi
		printf "ident_min=%s\nalen_min=%d\n" "$IDENTITY_MIN" "$ALEN_MIN" >"$PT_VARS"
		: >"$PT_DIAG"
	fi
	prof_end
else
	# Auto: use helper (mapq_min always 0; we removed the user knob)
	prof_start "pt_ident_threshold"
	((!DRY)) && python "${HELP_REMOVE}/polap-py-pt-ident-threshold.py" \
		--paf "$PAF_A" ${DBL_B:+ "$PAF_B"} \
		--pt-ids "${PT_ORIGIN:-/dev/null}" \
		--nuc-ids "${NUC_ORIGIN:-/dev/null}" \
		--alen-min "$ALEN_MIN" \
		--fpr "$FPR" --tpr "$TPR" \
		--diag "$PT_DIAG" \
		--emit-pt-ids "$PT_IDS" \
		>"$PT_VARS"
	prof_end
fi

# Parse vars (no eval)
IDENT_MIN_OUT="$FPR"
if [[ -s "$PT_VARS" ]]; then
	while IFS== read -r kv; do
		case "$kv" in
		ident_min=*) IDENT_MIN_OUT="${kv#ident_min=}" ;;
		alen_min=*) ALEN_MIN="${kv#alen_min=}" ;;
		esac
	done <"$PT_VARS"
fi
note "PT thresholds: ident_min=${IDENT_MIN_OUT}  alen_min=${ALEN_MIN}"

# ───────────── 3) Remove PT reads from full dataset ─────────────
prof_start "write_nonPT"
if ((!DRY)); then
	seqkit fx2tab -n "$READS" | sort -u >"$OUTDIR/map/all.ids"
	sort -u "$PT_IDS" >"$OUTDIR/map/pt.ids.sorted"
	comm -23 "$OUTDIR/map/all.ids" "$OUTDIR/map/pt.ids.sorted" >"$OUTDIR/keep.nonpt.ids"
	note "Write reads.nonpt.fq.gz"
	seqkit grep -f "$OUTDIR/keep.nonpt.ids" "$READS" -o "$OUTDIR/reads.nonpt.fq.gz"
fi
prof_end

# ───────────── 4) Optional nuclear deplete via overlapness ─────────────
if ((NUCLEAR_DEPLETE)); then
	[[ -s "$BUSCO" ]] || {
		note "ERR --busco required with --nuclear-deplete"
		exit 4
	}

	R1="$OUTDIR/reads.nonpt.fq.gz"
	[[ -s "$R1" ]] || {
		note "ERR missing reads.nonpt.fq.gz"
		exit 4
	}

	# 4A all-vs-all
	AVA_DIR="$OUTDIR/01_allvsall"
	mkdir -p "$AVA_DIR"
	PAF_ALL="$AVA_DIR/allvsall.paf"
	prof_start "allvsall"
	note "All-vs-all on non-PT reads → $PAF_ALL"
	((!DRY)) && minimap2 -x ava-ont -t "$THREADS" \
		-k19 -w7 --secondary=yes -N 30 --mask-level 0.60 \
		--min-occ-floor=10 -I4g $AVA_OPTS \
		"$R1" "$R1" >"$PAF_ALL"
	prof_end

	# 4B overlapness table
	prof_start "overlapness"
	note "Compute overlapness.tsv (min_olen=$MIN_OLEN, min_ident=$MIN_IDENT, w_floor=$W_FLOOR)"
	((!DRY)) && python "$PY_OVERLAP" "$PAF_ALL" \
		--min_olen "$MIN_OLEN" --min_ident "$MIN_IDENT" --w_floor "$W_FLOOR" \
		>"$OUTDIR/overlapness.tsv"
	prof_end

	# 4C BUSCO nuclear labels (sample or full)
	prof_start "miniprot_label"
	if ((MP_SAMPLE)); then
		note "Miniprot sampling mode ON (N=${MP_SAMPLE_N}, bins=${MP_SAMPLE_BINS})"
		# Build stratified sample ids for R1
		IFS=',' read -r B1 B2 B3 <<<"$MP_SAMPLE_BINS"
		seqkit fx2tab -nl "$R1" |
			awk -v b1="$B1" -v b2="$B2" -v b3="$B3" 'NR>1{
        len=$2;
        if(len<b1) b="b1"; else if(len<b2) b="b2"; else if(len<b3) b="b3"; else b="b4";
        print $1"\t"b
      }' >"$OUTDIR/reads.lenbin.tsv"
		# split evenly across bins
		per=$((MP_SAMPLE_N / 4))
		: >"$OUTDIR/sample.ids"
		for bb in b1 b2 b3 b4; do
			awk -v bb="$bb" '$2==bb{print $1}' "$OUTDIR/reads.lenbin.tsv" |
				seqkit sample -s "$MP_SAMPLE_SEED" -n "$per" -w 0 \
					>>"$OUTDIR/sample.ids" || true
		done
		sort -u "$OUTDIR/sample.ids" -o "$OUTDIR/sample.ids"
		seqkit grep -f "$OUTDIR/sample.ids" "$R1" -o "$OUTDIR/reads.sample.fq.gz"

		# miniprot on sample only → nuc.sample.ids
		miniprot -d "$OUTDIR/nonpt.sample.mpi" "$OUTDIR/reads.sample.fq.gz"
		miniprot -t "$THREADS" -S -N 3 --outc 0.4 \
			"$OUTDIR/nonpt.sample.mpi" "$BUSCO" >"$OUTDIR/nonpt.sample.busco.paf"
		awk 'BEGIN{FS=OFS="\t"} NF>=12 && $11+0>=150 {print $1}' \
			"$OUTDIR/nonpt.sample.busco.paf" | sort -u >"$OUTDIR/nuc.ids"
	else
		# full non-PT set
		miniprot -d "$OUTDIR/nonpt.mpi" "$R1"
		miniprot -t "$THREADS" -S -N 3 --outc 0.4 \
			"$OUTDIR/nonpt.mpi" "$BUSCO" >"$OUTDIR/nonpt.busco.paf"
		awk 'BEGIN{FS=OFS="\t"} NF>=12 && $11+0>=150 {print $1}' \
			"$OUTDIR/nonpt.busco.paf" | sort -u >"$OUTDIR/nuc.ids"
	fi
	prof_end

	# 4D thresholds from nuclear overlapness
	prof_start "overlap_thresh"
	note "Threshold from nuclear (FPR=$FPR_NUC)"
	((!DRY)) && python "$PY_NUC_THRESH" "$OUTDIR/overlapness.tsv" "$OUTDIR/nuc.ids" \
		--mode fpr --fpr "$FPR_NUC" \
		--diag "$OUTDIR/threshold_from_nuclear.tsv" >"$OUTDIR/nuc_thresh.vars"
	WDEG_MIN=0
	DEG_MIN=0
	if [[ -s "$OUTDIR/nuc_thresh.vars" ]]; then
		while IFS== read -r kv; do
			case "$kv" in
			wdeg_min=*) WDEG_MIN="${kv#wdeg_min=}" ;;
			deg_min=*) DEG_MIN="${kv#deg_min=}" ;;
			esac
		done <"$OUTDIR/nuc_thresh.vars"
	fi
	note "Overlapness thresholds: wdeg_min=$WDEG_MIN  deg_min=$DEG_MIN"
	prof_end

	# 4E write mt_enriched reads
	prof_start "write_mt_enriched"
	if ((!DRY)); then
		awk -v W="$WDEG_MIN" -v D="$DEG_MIN" 'BEGIN{FS=OFS="\t"}
      NR>1 && ($3+0)>=W && ($2+0)>=D {print $1}' "$OUTDIR/overlapness.tsv" |
			sort -u >"$OUTDIR/organelle.ids"

		comm -23 "$OUTDIR/organelle.ids" "$OUTDIR/nuc.ids" >"$OUTDIR/mt_enriched.ids"
		note "Write reads.mt_enriched.fq.gz"
		seqkit grep -f "$OUTDIR/mt_enriched.ids" "$R1" -o "$OUTDIR/reads.mt_enriched.fq.gz"
	fi
	prof_end
else
	note "--nuclear-deplete not set; stopping after PT removal."
fi

note "DONE remove-ptDNA (and optional nuclear-deplete)"
