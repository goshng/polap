#!/usr/bin/env bash
# polap-bash-fast-mtseed-ont.sh v0.5.1
# Version: v0.5.1 (PT removal + fast mtseed; steps start at 1)
#
# v0.5.1 highlights
#  - Step 2 switches to **edge-first hybrid**:
#      shortlist×ALL (targets sharded) -> per-shard **edge parts** (no PAF),
#      combine to a single **edges_loose.tsv.gz**, then **refilter** by a single
#      edge-weight threshold (eweight ladder) to produce **overlapness.tsv**.
#  - New CLI defaults (lowercase):
#      --hybrid-n-shards, --shortlist-frac, --shortlist-step-frac,
#      --shortlist-max-frac, --shortlist-target,
#      --eweight, --eweight-step, --eweight-minimum
#
set -euo pipefail

# options to adjust to have a smaller number of mito reads
#
# ../SRRxxx -> all of the reads
# mtseed/keep.nonpt.ids -> number of nuclear and mito reads
# mtseed/05-round/select_ids.txt -> number of mito reads
# then, miniasm plus minimap2

# ────────────────────────────────────────────────────────────────────
# Paths & environment
# ────────────────────────────────────────────────────────────────────
_POLAPLIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${_POLAPLIB_DIR}/polap-lib-conda.sh" || true
source "${_POLAPLIB_DIR}/polap-lib-math.sh" || true

# Helpers (fast-mtseed)
HELP="${_POLAPLIB_DIR}/fast-mtseed-ont"
SCAN_QCR="${HELP}/polap-r-scan-autotune-qc.R"
AUTOSCANPY="${HELP}/polap-py-autotune-scan.py"
OLPY="${HELP}/polap-py-overlapness-from-paf.py"
TOPPY="${HELP}/polap-py-topfrac-by-wdeg.py"
FILTPY="${HELP}/polap-py-filter-paf-by-ids.py"
CDSPY="${HELP}/polap-py-cds-coverage-from-paf.py"
SEEDREADSPY="${HELP}/polap-py-seed-mapped-reads.py"
STOPPY="${HELP}/polap-py-stop-delta.py"
AWK_CONS="${HELP}/polap-awk-filter-conservative.awk"
RECRUIT_SH="${HELP}/polap-bash-recruit-count.sh"
PAF2MREADS_SH="${HELP}/polap-bash-paf2mreads.sh"
R_PICK="${HELP}/polap-r-pick-high-recruit.R"
THRESH_NUC_PY="${HELP}/polap-py-threshold-from-nuclear.py"
PY_PT_THRESH="${HELP}/polap-py-pt-ident-threshold.py"
PT_ISOFORM_SH="${_POLAPLIB_DIR}/polap-bash-pt-isoform.sh"
EDGE_DUMP="${HELP}/polap-py-edge-dump.py"
REFILTER_EDGES="${HELP}/polap-py-refilter-edges-to-overlapness.py"
TOPKPY="${HELP}/paf-topk-per-q.py"
EDGPY="${HELP}/paf-to-edges-stream.py"
EDGECOMP="${HELP}/polap-py-edges-components.py"

# map-ont	Align noisy long reads of ~10% error rate to a reference genome. This is the default mode.
# ava-ont	Oxford Nanopore all-vs-all overlap mapping (-k15 -Xw5 -e0 -m100 -r2k).
# MINIASM_MAP_OPTS="--secondary=yes -N 30 --mask-level 0.60 --min-occ-floor 10"
# MINIASM_MAP_OPTS="-k19 -w10 -N 5 --mask-level 0.95 --min-occ-floor 50 -m 200"
MINIASM_MAP_OPTS="-k25 -w10 -N 1 -m 200"

# ────────────────────────────────────────────────────────────────────
# Logging & profiling
# ────────────────────────────────────────────────────────────────────
VERBOSE=1
DRY=0
LOG_FD="/dev/stderr"
TEE_ACTIVE=0
PROF=0
PROF_FILE=""

_note_base() {
	local lvl=$1
	shift
	local msg="$*"
	local src="${BASH_SOURCE[2]##*/}"
	local ln="${BASH_LINENO[1]}"
	local ts
	# ts="$(date +'%F %T')"
	ts="$(date +'%T')"
	local line="$ts [${src}:${ln}] $msg"
	if [[ "${VERBOSE:-0}" -gt "$lvl" ]]; then
		if [[ "${_arg_log_stderr:-off}" == "off" ]]; then
			printf "%s\n" "$line" >&3 || printf "%s\n" "$line" >&2
		else
			printf "%s\n" "$line" >&2
		fi
	fi
	printf "%s\n" "$line" >/dev/null
}
note0() { _note_base 0 "$@"; }
note1() { _note_base 1 "$@"; }
note2() { _note_base 2 "$@"; }
note3() { _note_base 3 "$@"; }

cmd() {
	if ((DRY)); then
		local p="$1"
		shift
		printf "[DRY] %q" "$p" | tee -a "$LOG_FD"
		for a in "$@"; do printf " %q" "$a" | tee -a "$LOG_FD"; done
		printf "\n" | tee -a "$LOG_FD"
		return 0
	fi
	local p="$1"
	shift
	printf "[RUN] %q" "$p" | tee -a "$LOG_FD"
	for a in "$@"; do printf " %q" "$a" | tee -a "$LOG_FD"; done
	printf "\n" | tee -a "$LOG_FD"
	"$p" "$@"
}
prof_start() { ((PROF)) && PROF_T0=$(date +%s) && PROF_STEP="$1" || true; }
prof_end() {
	((PROF)) || return 0
	local t1=$(date +%s)
	echo -e "${PROF_STEP}\t$((t1 - PROF_T0))" >>"$PROF_FILE"
}

require_tools() {
	local missing=() t
	for t in "$@"; do command -v "$t" >/dev/null 2>&1 || {
		echo "[ERR] missing tool: $t" >&2
		missing+=("$t")
	}; done
	if ((${#missing[@]} > 0)); then
		echo "[ERR] ${#missing[@]} tool(s) missing: ${missing[*]}" >&2
		echo "[HINT] activate env / load modules" >&2
		exit 127
	fi
	return 0
}

# ---------- helpers for integer-scaled thresholds (milli = 1/1000) ----------
to_milli() {
	local x="$1" int frac
	x="${x# +}"
	case "$x" in
	"")
		int=0
		frac=000
		;;
	*.*)
		int="${x%%.*}"
		frac="${x#*.}000"
		frac="${frac:0:3}"
		;;
	*)
		int="$x"
		frac=000
		;;
	esac
	[[ -z "$int" ]] && int=0
	printf '%d\n' "$((int * 1000 + 10#$frac))"
}

milli_to_str() {
	local m="$1" sgn=""
	if ((m < 0)); then
		sgn="-"
		m=$((-m))
	fi
	printf '%s%d.%03d\n' "$sgn" $((m / 1000)) $((m % 1000))
}

# Version: v0.1.1
# exists FILE
# Exit 2 if FILE does not exist or is empty.
# Version: v0.2.0
# exists FILE
# Check if a file exists and is non-empty.
# If missing or empty, print an error message including caller info and exit 2.

exists() {
	local file="$1"

	if [[ ! -s "$file" ]]; then
		# --- caller info ---
		local caller_func="${FUNCNAME[1]:-MAIN}"
		local caller_file="$(basename "${BASH_SOURCE[1]:-unknown}")"
		local caller_line="${BASH_LINENO[0]:-?}"

		printf "ERR: missing or empty file: %s\n" "$file" >&2
		printf "     (called from %s() in %s:%s)\n" \
			"$caller_func" "$caller_file" "$caller_line" >&2
		exit 2
	fi
}

# Version: v0.1.0
# have CMD
# Check if a command exists in PATH.
# Returns 0 if found, otherwise prints error and returns 2.
#
# examples:
# have seqkit || return 2
# have bc || return 2
have() {
	local cmd="$1"
	if ! command -v "$cmd" >/dev/null 2>&1; then
		echo "ERR: $cmd not found" >&2
		return 2
	fi
}

count_bases() {
	local seqkit_stats_ta_tsv="${1}.seqkit.stats.ta.tsv"
	if [[ ! -s "$seqkit_stats_ta_tsv" ]]; then
		seqkit stats -Ta "$1" -o "$seqkit_stats_ta_tsv"
	fi
	awk 'NR==2 {print $5}' "$seqkit_stats_ta_tsv"
}

# Sample either 10% of reads or ~1 Gbp of sequence, whichever is smaller.
# Version: v0.1.0
# sample_min_10pct_or_1gb IN_FASTQ OUT_FASTQ [SEED]
# Subsample to the smaller of 10% or 1 GB of bases, without using awk.
# Requires: seqkit, bc
sample_min_10pct_or_1gb() {
	local in="$1" out="$2" seed="${3:-13}"

	if [[ -z "$in" || -z "$out" ]]; then
		echo "usage: sample_min_10pct_or_1gb IN OUT [SEED]" >&2
		return 2
	fi
	have seqkit || return 2
	have bc || return 2

	# total bases from seqkit stats (TSV): header line then data line; col 5 = sum_len
	local total_bases="$(count_bases "$in")"
	# r1g = min(1.0, 1e9 / total_bases), scale to 6 decimals
	local r1g="$(bc -l <<<"scale=6; r=(10^9)/$total_bases; if (r>1) r=1; r")"
	# ratio = min(0.10, r1g) using bc only (no awk)
	local ratio="$(bc -l <<<"scale=6; a=0.10; b=$r1g; if (a<b) a else b")"

	note2 "total_bases: $total_bases"
	note2 "r1g: $r1g"
	note2 "ratio: $ratio"
	note1 seqkit sample -s "$seed" -p "$ratio" "$in" -o "$out"
	seqkit sample -s "$seed" -p "$ratio" "$in" -o "$out"
}

# ────────────────────────────────────────────────────────────────────
# CLI & defaults
# ────────────────────────────────────────────────────────────────────
reads=""
pt_ref=""
outdir="fast_mtseed_ont_out"
threads=32
pt_origin=""
nuc_origin=""
alen_min=3000
identity_min=""
fpr=0.01
tpr=0.95

# kept for Step 4 assembly branch
overlap_mode="scan+stream+edges"
topk_per_q=0
topk_edges_per_node=0
keep_scan_paf=0
gzip_edges=1
top_frac="0.20"
min_olen=1200
min_ident="0.84"
w_floor="0.12" # not used in edge-refilter ladder
assembler="miniasm"
mcds=""
nprot="busco_downloads/lineages/viridiplantae_odb12/refseq_db.faa.gz"
nuc_ids_opt=""
rounds_max=1
delta_stop="0.05"
seed=13
mm2_mode="${mm2_mode:-0}" # hybrid (edge-first)
map_mode="${map_mode:-0}"
use_parallel=0
do_shuffle=0
do_polap=1

steps=""

# HYBRID edge-first defaults (lowercase)
ava_mode="${ava_mode:-hybrid}"           # hybrid (edge-first)
hybrid_n_shards="${hybrid_n_shards:-8}"  # target shards
shortlist_frac="${shortlist_frac:-0.50}" # start fraction by length
shortlist_step_frac="${shortlist_step_frac:-0.10}"
shortlist_max_frac="${shortlist_max_frac:-0.80}"
shortlist_target="${shortlist_target:-0.35}" # coverage target (fraction of R1)

# single tuning knob for refilter -> eweight (edge weight floor)
eweight="${eweight:-0.80}"
eweight_step="${eweight_step:-0.10}"
eweight_minimum="${eweight_minimum:-0.01}"

usage() {
	cat <<EOF
polap-bash-fast-mtseed-ont.sh v0.5.1
Usage: $0 -r reads.fq.gz -p pt_ref.fa -o outdir [options]

# Hybrid (edge-first) knobs:
  --ava-mode hybrid
  --hybrid-n-shards INT        (default ${hybrid_n_shards})
  --shortlist-frac FLOAT       (default ${shortlist_frac})
  --shortlist-step-frac FLOAT  (default ${shortlist_step_frac})
  --shortlist-max-frac FLOAT   (default ${shortlist_max_frac})
  --shortlist-target FLOAT     (default ${shortlist_target})
  --eweight FLOAT              (initial w0, default ${eweight})
  --eweight-step FLOAT         (default ${eweight_step})
  --eweight-minimum FLOAT      (default ${eweight_minimum})
EOF
}

# Command-line processing
while [[ $# -gt 0 ]]; do
	case "$1" in
	-r | --reads)
		reads="$2"
		shift 2
		;;
	-p | --pt-ref)
		pt_ref="$2"
		shift 2
		;;
	-o | --outdir)
		outdir="$2"
		shift 2
		;;
	-t | --threads)
		threads="$2"
		TOTAL_CORES="$2"
		shift 2
		;;
	--ava-mode)
		ava_mode="$2"
		shift 2
		;;
	--hybrid-n-shards)
		hybrid_n_shards="$2"
		shift 2
		;;
	--shortlist-frac)
		shortlist_frac="$2"
		shift 2
		;;
	--shortlist-step-frac)
		shortlist_step_frac="$2"
		shift 2
		;;
	--shortlist-max-frac)
		shortlist_max_frac="$2"
		shift 2
		;;
	--shortlist-target)
		shortlist_target="$2"
		shift 2
		;;
	--eweight)
		eweight="$2"
		shift 2
		;;
	--eweight-step)
		eweight_step="$2"
		shift 2
		;;
	--eweight-minimum)
		eweight_minimum="$2"
		shift 2
		;;
	--nuc-ids)
		nuc_ids_opt="$2"
		shift 2
		;;
	--mm2-mode)
		mm2_mode="$2"
		shift 2
		;;
	--pt-origin)
		pt_origin="$2"
		shift 2
		;;
	--mt-origin)
		mt_origin="$2"
		shift 2
		;;
	--nuc-origin)
		nuc_origin="$2"
		shift 2
		;;
	--alen-min)
		alen_min="$2"
		shift 2
		;;
	--identity-min)
		identity_min="$2"
		shift 2
		;;
	--fpr)
		fpr="$2"
		shift 2
		;;
	--tpr)
		tpr="$2"
		shift 2
		;;
	--topk-per-q)
		topk_per_q="$2"
		shift 2
		;;
	--topk-edges-per-node)
		topk_edges_per_node="$2"
		shift 2
		;;
	--keep-scan-paf)
		keep_scan_paf=1
		shift
		;;
	--no-do-polap)
		do_polap=0
		shift
		;;
	--do-polap)
		do_polap=1
		shift
		;;
	--do-shuffle)
		do_shuffle=1
		shift
		;;
	--use-parallel)
		use_parallel=1
		shift
		;;
	--no-gzip-edges)
		gzip_edges=0
		shift
		;;
	--top-frac)
		top_frac="$2"
		shift 2
		;;
	--min-olen)
		min_olen="$2"
		shift 2
		;;
	--min-ident)
		min_ident="$2"
		shift 2
		;;
	--w-floor)
		w_floor="$2"
		shift 2
		;;
	-m | --mt-cds)
		mcds="$2"
		shift 2
		;;
	-n | --nuc-prot)
		nprot="$2"
		shift 2
		;;
	--rounds-max)
		rounds_max="$2"
		shift 2
		;;
	--delta-stop)
		delta_stop="$2"
		shift 2
		;;
	--seed)
		seed="$2"
		shift 2
		;;
	--step)
		steps="$2"
		shift 2
		;;
	-v | --verbose)
		VERBOSE=$((VERBOSE + 1))
		shift
		;;
	--quiet | -q)
		VERBOSE=0
		shift
		;;
	--dry)
		DRY=1
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
		echo "ERR unknown arg: $1" >&2
		usage
		exit 2
		;;
	esac
done

# Step gating (robust)
_step_set="" # e.g. " 1 2 3 "
_add_steps() {
	local spec="$1"
	_step_set=""

	# Normalize separators: commas as split points, convert Unicode dashes to ASCII
	# and strip spaces around tokens.
	local raw="$spec"
	# replace unicode dashes with '-'
	raw="${raw//–/-}"
	raw="${raw//—/-}"
	raw="${raw//−/-}"
	# split on commas into array
	IFS=',' read -r -a arr <<<"$raw"

	local tok a b i
	for tok in "${arr[@]}"; do
		# trim spaces in token
		tok="${tok//[[:space:]]/}"
		[[ -z "$tok" ]] && continue

		if [[ "$tok" =~ ^([0-9]+)-([0-9]+)$ ]]; then
			a=${BASH_REMATCH[1]}
			b=${BASH_REMATCH[2]}
			# allow reversed ranges like 3-1
			if ((a <= b)); then
				for ((i = a; i <= b; i++)); do _step_set+=" $i "; done
			else
				for ((i = b; i <= a; i++)); do _step_set+=" $i "; done
			fi
		elif [[ "$tok" =~ ^[0-9]+$ ]]; then
			_step_set+=" $tok "
		else
			# ignore invalid token but log once at verbose level
			note2 "Ignoring malformed --step token: '$tok' (from '$spec')"
		fi
	done

	# optional: collapse duplicates
	if [[ -n "$_step_set" ]]; then
		# sort unique numeric tokens
		local uniq
		uniq=$(tr ' ' '\n' <<<"$_step_set" | awk 'NF' | sort -n | uniq | xargs)
		_step_set=" $uniq "
	fi
	return 0
}

_should_run() {
	local s="$1"

	# no --step provided -> run everything
	if [[ -z "${steps:-}" ]]; then
		return 0
	fi

	# materialize once (or re-materialize if you want to allow updates)
	if [[ -z "${_step_set:-}" ]]; then
		_add_steps "$steps" || return 1
		note2 "Step gating: steps='${steps}' -> _step_set='${_step_set}'"
	fi

	# explicit return codes (safe under set -e)
	if [[ " ${_step_set} " == *" $s "* ]]; then
		return 0
	else
		return 1
	fi
}

if [[ -n "${steps:-}" ]]; then
	_add_steps "$steps"
	note2 "Parsed --step: steps='${steps}' -> _step_set='${_step_set}'"
fi

note2 "Check: tools"
require_tools minimap2 samtools seqkit python Rscript awk sort comm cut gzip
note2 "Check: scripts"
for s in "$EDGE_DUMP" "$REFILTER_EDGES" "$OLPY" "$TOPPY" "$FILTPY" "$CDSPY" \
	"$SEEDREADSPY" "$STOPPY" "$AWK_CONS" "$RECRUIT_SH" "$PAF2MREADS_SH" \
	"$R_PICK" "$THRESH_NUC_PY" "$PT_ISOFORM_SH"; do

	if [[ ! -s "$s" ]]; then
		note0 "[ERR] missing helper: $s"
		exit 1
	fi
done

if [[ -z "$reads" || ! -s "$reads" ]]; then
	note0 "ERR missing --reads"
	exit 2
fi

if [[ -z "$pt_ref" || ! -s "$pt_ref" ]]; then
	note0 "ERR missing --pt-ref"
	exit 2
fi

mkdir -p "$outdir" "$outdir/log"
LOG="${outdir}/pipeline.log"
LOG_FD="$LOG"
if ((!DRY)); then
	exec > >(tee -a "$LOG") 2>&1
	TEE_ACTIVE=1
fi

if ((PROF)); then
	PROF_FILE="$outdir/log/profile.tsv"
	echo -e "step\tseconds" >"$PROF_FILE"
fi

note0 "pipeline: outdir=$outdir assembler=$assembler threads=$threads"

# ────────────────────────────────────────────────────────────────────
# Step 0: mm2 options
# ────────────────────────────────────────────────────────────────────
# Presets for ONT shard mapping with minimap2
# You will call:
#   minimap2 -x ava-ont $INDEX_OPTS -d "$idx" "$shard"
#   minimap2 -x ava-ont -t "$threads" $MAP_OPTS "$idx" "$infile"
#
# Inputs (override before calling):
#   TOTAL_CORES   (default 48)   # total CPU cores on node
#   TOTAL_RAM_GB  (default 128)  # total RAM (for clamping -K if needed)
#   mm2_mode      (0|1|2)        # 0=lean, 1=balanced, 2=sensitive
#   HPC           (0|1)          # 1 => add -H at index time (HPC seeds)
#
# Outputs set by this function:
#   hybrid_n_shards   # GNU parallel -j (concurrency)
#   threads           # minimap2 -t per job
#   INDEX_OPTS        # index-only options (includes seeding: -H/-k/-w and index chunking -I)
#   MAP_OPTS          # map-only options (NO seeding: reporting/filters/batch size)

TOTAL_CORES="${threads:-48}"
TOTAL_RAM_GB="${TOTAL_RAM_GB:-128}"
mm2_mode="${mm2_mode:-1}"
HPC="${HPC:-0}"

choose_mm2_preset() {
	local hflag=""
	((HPC == 1)) && hflag="-H"

	case "$mm2_mode" in
	0) # LEAN: lowest RAM/job, highest concurrency
		hybrid_n_shards=8
		threads=$((TOTAL_CORES / hybrid_n_shards))
		((threads < 4)) && threads=4
		# index-only: embed seeding & index chunking
		INDEX_OPTS="$hflag -k 27 -w 11 -I 1g"
		# map-only: reporting/filters/batch; NO -k/-w/-H here
		MAP_OPTS="--secondary=yes -N 10 --mask-level 0.70 --min-occ-floor 20 -K 256m"
		;;

	1) # BALANCED: good default
		hybrid_n_shards=6
		threads=$((TOTAL_CORES / hybrid_n_shards))
		((threads < 6)) && threads=6
		INDEX_OPTS="$hflag -k 23 -w 9 -I 2g"
		MAP_OPTS="--secondary=yes -N 30 --mask-level 0.60 --min-occ-floor 10 -K 512m"
		;;

	2) # SENSITIVE: heavier
		hybrid_n_shards=4
		threads=$((TOTAL_CORES / hybrid_n_shards))
		((threads < 8)) && threads=8
		INDEX_OPTS="$hflag -k 19 -w 7 -I 2g"
		MAP_OPTS="--secondary=yes -N 50 --mask-level 0.55 --min-occ-floor 8 -K 512m"
		;;

	*)
		echo "WARN: unknown mm2_mode='$mm2_mode' -> using balanced (1)" >&2
		mm2_mode=1
		hybrid_n_shards=6
		threads=$((TOTAL_CORES / hybrid_n_shards))
		((threads < 6)) && threads=6
		INDEX_OPTS="$hflag -k 17 -w 5 -I 2g"
		MAP_OPTS="--secondary=yes -N 30 --mask-level 0.60 --min-occ-floor 10 -K 512m"
		;;
	esac

	# Do not oversubscribe cores
	if ((hybrid_n_shards * threads > TOTAL_CORES)); then
		threads=$((TOTAL_CORES / hybrid_n_shards))
		((threads < 4)) && threads=4
	fi

	# If node RAM is tight, clamp map batch size
	if ((TOTAL_RAM_GB < 64)); then
		MAP_OPTS="${MAP_OPTS/-K 512m/-K 256m}"
		INDEX_OPTS="${INDEX_OPTS/-I 2g/-I 1g}"
	fi
}

ST0_check() {
	local outdir="$1"

	echo 1 >&2

	exists "$nprot"
}

# ST1: Build plastid (PT) isoforms and doubled single-record references
# Usage:
#   ST1_pt_isoform "$outdir"
# Requires (from caller env):
#   PT_ISOFORM_SH : script that produces pt_isomerA.fa (and optionally pt_isomerB.fa)
#   pt_ref        : input plastid reference (FA/GBK/etc. per PT_ISOFORM_SH)
#   threads       : cores to use
#   DRY           : 0/1 (skip execution when DRY==1)
#   note0, prof_start, prof_end : logging/timing helpers
ST1_pt_isoform() {
	local outdir="$1"

	note2 "1a) Build PT isoforms A/B"
	mkdir -p "$PANEL_DIR"
	prof_start "pt_isoform"
	cmd bash "$PT_ISOFORM_SH" -r "$pt_ref" -o "$PANEL_DIR" -t "$threads"
	prof_end

	local ISO_A="$PANEL_DIR/pt_isomerA.fa"
	local ISO_B="$PANEL_DIR/pt_isomerB.fa"

	if [[ ! -s "$ISO_A" ]]; then
		echo "ERR missing pt_isomerA.fa at $ISO_A" >&2
		exit 3
	fi

	note2 "1b) Double isoforms (single-record)"
	dbld() {
		python - "$1" "$2" <<'PY'
import sys
inp,outp=sys.argv[1],sys.argv[2]
name=None; seq=[]
with open(inp) as f:
    for ln in f:
        if ln.startswith('>'):
            if name is None:
                name=ln[1:].strip().split()[0]
            else:
                break
        else:
            seq.append(ln.strip())
S=''.join(seq)
with open(outp,'w') as w:
    w.write('>pt.double\n')
    w.write(S+S+'\n')
PY
	}

	local DBL_A="$PANEL_DIR/pt_isomerA.double.fa"
	local DBL_B="$PANEL_DIR/pt_isomerB.double.fa"

	prof_start "pt_double"
	((!DRY)) && dbld "$ISO_A" "$DBL_A"
	((!DRY)) && { [[ -s "$ISO_B" ]] && dbld "$ISO_B" "$DBL_B" || true; }
	prof_end

	# Export for downstream steps if you like
	# export PANEL_DIR ISO_A ISO_B DBL_A DBL_B
}

# ST2: Map reads to doubled plastid isoforms (A/B) and remove PT reads
# Usage:
#   ST2_pt_map "$outdir" "$reads"
# Requires (from env / prior stages):
#   DBL_A, DBL_B, threads, DRY, note0/note1, prof_start/prof_end, cmd
#   identity_min (opt), alen_min (req for identity_min path), PY_PT_THRESH (else path)
#   pt_origin (opt), mt_origin (opt), fpr (opt), tpr (opt)
ST2_pt_map() {
	local outdir="$1"
	local reads="$2"

	note2 "1c) Map reads -> doubled A/B"
	local MAP_DIR="$outdir/02-map"
	mkdir -p "$MAP_DIR"
	local PAF_A="$MAP_DIR/formA.paf"
	local PAF_B="$MAP_DIR/formB.paf"
	local DBL_A="$PANEL_DIR/pt_isomerA.double.fa"
	local DBL_B="$PANEL_DIR/pt_isomerB.double.fa"

	prof_start "map_ptA"
	cmd minimap2 -x map-ont --secondary=yes -N 50 -t "${TOTAL_CORES}" "$DBL_A" "$reads" >"$PAF_A"
	prof_end

	if [[ -s "${DBL_B:-}" ]]; then
		prof_start "map_ptB"
		cmd minimap2 -x map-ont --secondary=yes -N 50 -t "${TOTAL_CORES}" "$DBL_B" "$reads" >"$PAF_B"
		prof_end
	else
		PAF_B=""
	fi

	note2 "1d) MT-guided identity cutoff & emit pt.ids"
	local PT_VARS="$outdir/pt_thresh.vars"
	local PT_DIAG="$outdir/pt_thresh.diag.tsv"
	local PT_IDS="$outdir/pt.ids"

	# skip because identity_min is null.
	if [[ -n "${identity_min:-}" ]]; then
		note2 "using --identity-min=${identity_min} ; alen_min=${alen_min}"
		if ((!DRY)); then
			awk -v ID="$identity_min" -v AL="${alen_min:-0}" 'BEGIN{FS=OFS="\t"} NF>=12{ id=($11>0?$10/$11:0); if(id>=ID && $11+0>=AL) print $1 }' "$PAF_A" | sort -u >"$outdir/ptA.ids"
			if [[ -s "$PAF_B" ]]; then
				awk -v ID="$identity_min" -v AL="${alen_min:-0}" 'BEGIN{FS=OFS="\t"} NF>=12{ id=($11>0?$10/$11:0); if(id>=ID && $11+0>=AL) print $1 }' "$PAF_B" | sort -u >"$outdir/ptB.ids"
				sort -u "$outdir/ptA.ids" "$outdir/ptB.ids" >"$PT_IDS"
			else
				mv "$outdir/ptA.ids" "$PT_IDS"
			fi
			printf "ident_min=%s\nalen_min=%d\n" "$identity_min" "${alen_min:-0}" >"$PT_VARS"
			: >"$PT_DIAG"
		fi
	else
		note1 "using $PY_PT_THRESH -> file: $PT_VARS and PT IDs: $PT_IDS"
		prof_start "pt_ident_threshold"
		if ((!DRY)); then
			note1 python "$PY_PT_THRESH" \
				--paf "$PAF_A" ${PAF_B:+ "$PAF_B"} \
				--pt-ids "${pt_origin:-/dev/null}" --mt-ids "${mt_origin:-/dev/null}" \
				--alen-min "${alen_min:-0}" --fpr "${fpr:-0.10}" --tpr "${tpr:-0.90}" \
				--diag "$PT_DIAG" --emit-pt-ids "$PT_IDS" ">$PT_VARS"
			python "$PY_PT_THRESH" \
				--paf "$PAF_A" ${PAF_B:+ "$PAF_B"} \
				--pt-ids "${pt_origin:-/dev/null}" --mt-ids "${mt_origin:-/dev/null}" \
				--alen-min "${alen_min:-0}" --fpr "${fpr:-0.10}" --tpr "${tpr:-0.90}" \
				--diag "$PT_DIAG" --emit-pt-ids "$PT_IDS" >"$PT_VARS"
		fi
		prof_end
	fi

	note2 "1e) Reading $PT_VARS to remove PT reads -> R1"

	# delete this: not used identity_min and alen_min values later
	#
	# local ident_min_out="${identity_min:-}"
	# if [[ -s "$PT_VARS" ]]; then
	# 	# robust parsing: read key=value pairs safely
	# 	while IFS='=' read -r key val; do
	# 		case "$key" in
	# 		ident_min) ident_min_out="$val" ;;
	# 		alen_min) alen_min="$val" ;;
	# 		esac
	# 	done <"$PT_VARS"
	# fi

	if ((!DRY)); then
		seqkit fx2tab -ni "$reads" | sort -u >"$MAP_DIR/all.ids"
		sort -u "$PT_IDS" >"$MAP_DIR/pt.ids.sorted"
		comm -23 "$MAP_DIR/all.ids" "$MAP_DIR/pt.ids.sorted" >"$outdir/keep.nonpt.ids"
		seqkit grep -f "$outdir/keep.nonpt.ids" "$reads" -o "$outdir/reads.nonpt.fq.gz"
	fi

	R1="$outdir/reads.nonpt.fq.gz"
	seqkit stats -Ta "$R1" >"$R1".seqkit.stats.ta.tsv
	# export R1
}

# ────────────────────────────────────────────────────────────────────
# Shared helpers for Step 2
# ────────────────────────────────────────────────────────────────────

# coverage = (#nonzero rows in overlapness) / (#reads in R1)
_ovl_cov() {
	local ovl="$1"
	local reads="$2"
	local nz tot
	nz=$(awk 'NR>1{c++} END{print (c+0)}' "$ovl" 2>/dev/null || echo 0)

	tot=0

	local seqkit_stats_ta="${reads}.seqkit.stats.Ta.tsv"
	if [[ ! -s "$seqkit_stats_ta" ]]; then
		seqkit stats -Ta "$reads" -o "$seqkit_stats_ta"
	fi
	tot=$(cat "$seqkit_stats_ta" 2>/dev/null | awk 'NR==2{print $4+0}')

	awk -v a="$nz" -v b="$tot" 'BEGIN{ if (b>0) printf "%.6f", a/b; else print 0 }'
}

# shortlist top fraction by length -> <outdir>/R1.short.fq.gz
_build_shortlist_len() {
	local reads="$1"
	local outfile="$2"
	local frac="${3:-0.5}"
	seqkit fx2tab -nil "$reads" |
		sort -k2,2nr |
		awk -v f="$frac" \
			'{ids[NR-1]=$1} END{lim=int((NR-1)*f); for(i=0;i<lim;i++) print ids[i]}' \
			>"$outfile"
}

# shard target once (ALL reads)
_shard_targets_once() {
	local reads="$1"
	local shard_dir="$2"
	local B="${3:-4}"

	# Pre-flight checks
	[[ -n "$reads" && -s "$reads" ]] || {
		note0 "ERR: reads not set or empty: '$reads'"
		return 1
	}
	[[ -n "$shard_dir" ]] || {
		note0 "ERR: shard_dir not set"
		return 1
	}

	# If shard_dir already has files, consider sharding done
	if [[ -d "$shard_dir" && -n "$(ls -A "$shard_dir" 2>/dev/null)" ]]; then
		note1 "reusing existing shards at $shard_dir"
		return 0
	fi

	# Create dir
	if ! mkdir -p "$shard_dir"; then
		note0 "ERR: cannot create shard_dir '$shard_dir'"
		return 1
	fi

	# Run the split; capture exit status explicitly
	note0 "hybrid: sharding targets into $B blocks at $shard_dir"
	if ! seqkit split2 -p "$B" -O "$shard_dir" "$reads"; then
		note0 "ERR: seqkit split2 failed"
		return 1
	fi

	# Verify shards exist
	if [[ -z "$(ls -A "$shard_dir" 2>/dev/null)" ]]; then
		note0 "ERR: no shards found after split"
		return 1
	fi

	return 0
}

# MM2_THOROUGH_OPTS="-k 17 -w 5 --secondary=yes -N 40 --mask-level 0.55 --min-occ-floor 8 -I 4g -K 2g"

# 2b) map shortlist×each shard -> per-shard edge parts
# _map_edges_parts <indir> <infile> <outdir>
# Map shortlist (infile) × each shard in indir and dump edges to outdir.
# Produces outdir/edges_part_XX.tsv.gz for each shard XX.
_map_edges_parts() {
	# ---- Inputs ----
	local indir="$1"
	local infile="$2"
	local outdir="$3"

	# ---- Pre-flight checks ----
	if [[ -z "$indir" || -z "$infile" || -z "$outdir" ]]; then
		note0 "ERR[_map_edges_parts]: missing args (indir='$indir' infile='$infile' outdir='$outdir')"
		return 1
	fi
	if [[ ! -s "$infile" ]]; then
		note0 "ERR[_map_edges_parts]: infile not found or empty: '$infile'"
		return 1
	fi
	if [[ ! -d "$indir" ]]; then
		note0 "ERR[_map_edges_parts]: indir not a directory: '$indir'"
		return 1
	fi
	# ensure there is at least one shard
	if ! ls "$indir"/*part_*.* >/dev/null 2>&1; then
		note0 "ERR[_map_edges_parts]: no shard files found in '$indir'"
		return 1
	fi
	if [[ -z "${hybrid_n_shards:-}" ]]; then
		note1 "WARN[_map_edges_parts]: hybrid_n_shards is unset; defaulting to 4"
		hybrid_n_shards=4
	fi
	if ! mkdir -p "$outdir"; then
		note0 "ERR[_map_edges_parts]: cannot create outdir '$outdir'"
		return 1
	fi

	# ---- Clean old parts ----
	rm -f "$outdir/"*.tsv.gz 2>/dev/null || true

	# ---- Map per shard -> edges_part_XX.tsv.gz ----
	note2 "mapping shortlist × shards (thorough) ->  edge parts"
	if [[ "${use_parallel}" -eq 1 ]]; then

		if command -v parallel >/dev/null 2>&1; then

			parallel -j 4 --halt now,fail=1 --line-buffer \
				'minimap2 -x ava-ont -t '"$threads"' '"$MM2_THOROUGH_OPTS"' "'"$infile"'" {} \
       | python "'"$EDGE_DUMP"'" --paf - --out "'"$PARTS_DIR"'/edges_part_{#}.tsv.gz" \
       --min_olen 0 --min_ident 0 --w_floor 0 --dedup' \
				::: "$indir"/*.part_*.fq*

		fi

	else

		# build index per shard (if needed), then map with the .mmi
		# target = shard index (.mmi), query = infile
		# keep memory low: streaming dump (NO in-Python dedup), then disk dedup per shard
		local j j2 shard idx
		local INDEX_OPTS="${INDEX_OPTS:--x ava-ont}" # e.g., "-x ava-ont -H -k 17 -w 5" for HPC seeds
		local INDEX_I="${INDEX_I:-2g}"               # index build chunk size
		local SORT_MEM="${SORT_MEM:-2G}"             # memory cap for external sort
		local SORT_TMP="${SORT_TMP:-$outdir/tmp}"    # temp dir for sort spills
		mkdir -p "$SORT_TMP"

		for j in $(seq 1 "$hybrid_n_shards"); do
			j2=$(printf "%03d" "$j")
			shard=$(ls "$indir"/*.part_${j2}.* 2>/dev/null | head -n1)
			if [[ -z "$shard" ]]; then
				note0 "ERR[_map_edges_parts]: missing shard for index $j2 in '$indir'"
				return 1
			fi

			# choose index path next to the shard (or move to a dedicated mmi dir if you prefer)
			case "$shard" in
			*.gz) idx="${shard%.gz}.mmi" ;;
			*) idx="${shard}.mmi" ;;
			esac

			# (re)build the index if missing or older than the shard
			# if [[ ! -s "$idx" || "$idx" -ot "$shard" ]]; then
			rm -f "$idx"
			note2 "hybrid: indexing shard $j2 ->  $(basename "$idx")"
			if ! minimap2 -x ava-ont $INDEX_OPTS -d "$idx" "$shard"; then
				note0 "ERR[_map_edges_parts]: index build failed for shard $j2 ($shard)"
				return 1
			fi
			# fi

			# paths per shard
			RAW="$outdir/edges_raw_${j2}.tsv.gz"
			PART="$outdir/edges_part_${j2}.tsv.gz"

			# map with index, stream-dump edges (NO dedup, NO gates) — low RAM
			note2 "hybrid: shard $j2/$hybrid_n_shards (map with index -> raw edges)"
			if ! minimap2 -x ava-ont -t "$TOTAL_CORES" $MAP_OPTS \
				"$idx" "$infile" |
				python "$EDGE_DUMP" --paf - --out "$RAW" \
					--min_olen 0 --min_ident 0 --w_floor 0; then
				note0 "ERR[_map_edges_parts]: mapping or raw edge dump failed for shard $j2"
				return 1
			fi

			# disk dedup per shard (bounded by sort -S): keep one max-weight edge per unordered pair
			# input cols assumed: u v alen ident weight
			note2 "hybrid: shard $j2 dedup on disk (bounded RAM)"
			if ! zcat "$RAW" |
				awk 'BEGIN{FS=OFS="\t"} NF>=5{u=$1;v=$2; if(u>v){t=u;u=v;v=t} print u,v,$3,$4,$5}' |
				LC_ALL=C sort -S "$SORT_MEM" -T "$SORT_TMP" -k1,1 -k2,2 -k5,5gr |
				awk 'BEGIN{FS=OFS="\t"} !seen[$1 FS $2]++' |
				gzip -1 >"$PART"; then
				note0 "ERR[_map_edges_parts]: disk dedup failed for shard $j2"
				return 1
			fi
		done

	fi

	# ---- Post-check results ----
	if ! ls "$outdir/"*.tsv.gz >/dev/null 2>&1; then
		note0 "ERR[_map_edges_parts]: no edge parts produced in '$outdir'"
		return 1
	fi

	# ---- Success ----
	return 0
}

# 2c) combine to edges_loose.tsv.gz
_combine_edges() {
	local indir="$1"
	local outfile="$2"

	note2 "hybrid: combine edge parts ->  $outfile"

	# sanity checks
	if [[ -z "$indir" || -z "$outfile" ]]; then
		note0 "ERR[_combine_edges]: missing args (indir='$indir' outfile='$outfile')"
		return 1
	fi
	if ! ls "$indir"/*.tsv.gz >/dev/null 2>&1 2>/dev/null; then
		note0 "ERR[_combine_edges]: no edge parts (*.tsv.gz) found in '$indir'"
		return 1
	fi

	# combine
	if ! zcat "$indir"/*.tsv.gz | gzip -1 >"$outfile"; then
		note0 "ERR[_combine_edges]: combine failed"
		return 1
	fi

	# check result
	if [[ ! -s "$outfile" ]]; then
		note0 "ERR[_combine_edges]: output '$oroundutfile' is empty"
		return 1
	fi

	return 0
}

# Expect these to be set:
#   EDGES_COMBINED="$ST3/edges_loose.tsv.gz"
#   OVL_STRICT="$ST3/overlapness.tsv"
#   REFILTER_EDGES="${HELP}/polap-py-refilter-edges-to-overlapness.py"
#   R1 (reads), SHORT_DIR, SHARD_DIR, PARTS_DIR
#   shortlist_frac, shortlist_step_frac, shortlist_max_frac, shortlist_target
#   eweight (start), eweight_step, eweight_minimum
# Also expect helpers:
#   _ovl_cov  (prints coverage as fraction 0..1)
#   _map_edges_parts  (maps shortlist x shards -> PARTS_DIR/edges_part_XX.tsv.gz)
#   _combine_edges    (zcat parts -> EDGES_COMBINED)

# 2d) refilter with eweight ladder (inline loop) until coverage >= target or floor reached
_refilter_once_eweight_loop() {
	# returns 0 on success (target met), 1 if floor reached without target, 2 on hard error
	local ew cov cov_pct tgt_pct
	local best_cov="0"
	local best_ew=""

	# sanity checks
	if [[ ! -s "$EDGES_COMBINED" ]]; then
		note0 "ERR: edges archive missing: $EDGES_COMBINED"
		return 2
	fi
	if awk -v s="$eweight_step" 'BEGIN{exit !(s>0)}'; then :; else
		note0 "ERR: edge weight step size must be > 0 (got '$eweight_step')"
		return 2
	fi

	ew=$(awk -v x="$eweight" 'BEGIN{printf "%.3f", x}')
	while :; do
		note0 "hybrid: refilter edges -> overlapness (eweight=$ew)"
		if ! python "$REFILTER_EDGES" \
			--edges "$EDGES_COMBINED" \
			--w_floor "$ew" \
			--out "$OVL_STRICT"; then
			note0 "ERR: refilter failed at eweight=$ew"
			return 2
		fi

		cov=$(_ovl_cov "$OVL_STRICT" "$R1")
		cov="${cov:-0}"
		cov_pct=$(awk -v c="$cov" 'BEGIN{printf "%.2f", 100*c}')

		local tgt_pct=$(_polap_lib_math-percentify $shortlist_target)
		# tgt_pct=$(awk -v t="$shortlist_target" 'BEGIN{printf "%.2f", 100*t}')
		# tgt_pct=$(printf "%.2f" "$(bc -l <<<"$shortlist_target * 100")")

		note1 "coverage=${cov_pct}% target=${tgt_pct}% (eweight=$ew)"

		# track best attempt
		if awk -v a="$cov" -v b="$best_cov" 'BEGIN{exit !(a>b)}'; then
			best_cov="$cov"
			best_ew="$ew"
		fi

		# success?
		if awk -v c="$cov" -v tgt="$shortlist_target" 'BEGIN{exit !(c>=tgt)}'; then
			eweight="$ew"
			return 0
		fi

		# descend eweight; stop at floor
		ew=$(awk -v x="$ew" -v s="$eweight_step" -v f="$eweight_minimum" 'BEGIN{v=x-s; if(v<f) v=f; printf "%.3f", v}')
		if awk -v x="$ew" -v f="$eweight_minimum" 'BEGIN{exit !(x<=f)}'; then
			# floor hit; keep best attempt (regenerate if needed)
			if [[ -n "$best_ew" && "$best_ew" != "$ew" ]]; then
				note0 "hybrid: floor reached; restoring best attempt at eweight=$best_ew (best_cov=$(awk -v c="$best_cov" "BEGIN{printf \"%.2f\",100*c}")%)"
				if ! python "$REFILTER_EDGES" \
					--edges "$EDGES_COMBINED" \
					--w_floor "$best_ew" \
					--out "$OVL_STRICT"; then
					note0 "ERR: failed to restore best attempt at eweight=$best_ew"
					return 2
				fi
				eweight="$best_ew"
			fi
			return 1
		fi
	done
}

#!/usr/bin/env bash
# 03-partvspart: randomized sharding + intra-shard mapping -> edges.tsv
# Deps: minimap2, seqkit, gzip, python (for EDGE_DUMP only)
# Presets: always use -x ava-ont per your snippet
#
# Usage:
#   ST3_run READS OUTDIR TOTAL_CORES [TARGET_GBP] [SEED]
#   Example:
#     ST3_run reads.fq.gz out 16 1 13
#
# Env knobs (optional):
#   INDEX_OPTS=""      # extra flags for 'minimap2 -d'
#   MAP_OPTS=""        # extra flags for mapping
#   EDGE_DUMP="polap-py-edge-dump.py"   # path to your edge dumper (PAF -> SV)
#
# Output:
#   OUTDIR/03-partvspart/
#     shuffled/reads.shuf.fq.gz
#     shards/shard_000.fq.gz ... shard_{B-1}.fq.gz
#     edges/shard_000.edges.tsv.gz ... shard_{B-1}.edges.tsv.gz
#     edges/all.edges.tsv.gz

TARGET_GBP="${TARGET_GBP:-1}" # target ~ Gbp per shard (default 1 Gbp)
ST3_part_vs_part() {
	local OUTDIR="$1"
	local READS="$2"

	# Use R1 not reads to compute the overlapness
	# local READS="${reads}"

	local SEED="${SEED:-13}"

	: "${INDEX_OPTS:=}"
	: "${MAP_OPTS:=}"

	local STDIR="${OUTDIR}/03-allvsall"

	local SHARD_DIR="${STDIR}/01-shards"
	local EDGE_DIR="${STDIR}/02-edges"

	mkdir -p "$SHARD_DIR" "$EDGE_DIR"

	# --- helpers ---
	_zcat() { case "$1" in *.gz) gzip -cd -- "$1" ;; *) cat -- "$1" ;; esac }

	# Version: v0.1.2
	# calc_shards TOTAL_BASES TARGET_GBP
	# Computes number of shards (B) from total bases and target gigabases.
	# Ensures at least one shard and prints the result.
	calc_shards() {
		local TOTAL_BASES="$1"
		local TARGET_GBP="$2"
		local TARGET_BASES B

		# 1 gigabase = 10^9 bases
		local GIGA=1000000000

		TARGET_BASES=$((TARGET_GBP * GIGA))
		B=$(((TOTAL_BASES + TARGET_BASES - 1) / TARGET_BASES))

		if [[ "$B" -lt 1 ]]; then
			B=1
		fi

		echo "$B"
	}

	local TOTAL_BASES="$(count_bases "$READS")"
	note1 "[ST3] Count total bases in: $READS -> $TOTAL_BASES"
	local TARGET_GBP=1
	local B=$(calc_shards "$TOTAL_BASES" "$TARGET_GBP")
	note1 "[ST3] Number of shards of $TARGET_GBP Gb: $B"

	note1 "[ST3] Splitting shuffled reads into $B shards ->  $SHARD_DIR"
	# seqkit split2 creates reads.shuf.part_001.gz, ...
	rm -rf "$SHARD_DIR"
	mkdir -p "$SHARD_DIR"
	seqkit split2 -p "$B" -O "$SHARD_DIR" "$READS"
	local i=1
	for p in "$SHARD_DIR"/*.fq*; do
		printf -v idx "%03d" "$i"
		mv -f "$p" "${SHARD_DIR}/shard_${idx}.fq.gz"
		((i++))
	done

	# Parallelized per-shard index + map -> edges.tsv.gz using GNU parallel
	# Requires: GNU parallel. Falls back to serial if not available.
	# Uses env vars: B, SHARD_DIR, EDGE_DIR, INDEX_OPTS, MAP_OPTS, TOTAL_CORES, EDGE_DUMP, note1
	# Worker for a single shard index (000-padded) like "001"
	_st3_process_shard() {
		local idx="$1"
		local shard="${SHARD_DIR}/shard_${idx}.fq.gz"
		local idxpath="${EDGE_DIR}/shard_${idx}.mmi"
		local infile="${SHARD_DIR}/shard_${idx}.fq.gz"
		local raw="${EDGE_DIR}/shard_${idx}.edges.tsv.gz"
		local log="${EDGE_DIR}/shard_${idx}.minimap2.log"

		if ! minimap2 -x ava-ont ${INDEX_OPTS:-} -d "$idxpath" "$shard" 2>"$log"; then
			return 0
		fi

		# minimap2 -> dump:
		# edge_u edge_v alen ident weight
		# edge_u, edge_v: read IDs
		# alen: PAF column 11
		# ident: col10/col11 of PAF or number of match / alignment length
		# weight: ident * (alen / min(qlen,tlen))
		#
		# where
		# qlen: col01 of PAF
		# tlen: col06 of PAF
		if ! minimap2 -x ava-ont -t "${MM2_TOTAL_CORES:-8}" ${MAP_OPTS:-} \
			"$idxpath" "$infile" 2>>"$log" |
			python "$EDGE_DUMP" --paf - --out - --min_olen 0 --min_ident 0 --w_floor 0 |
			gzip -c >"$raw"; then
			return 0
		fi
	}

	if [[ "${use_parallel:-0}" -eq 1 ]] && command -v parallel >/dev/null 2>&1; then
		# Export the worker for GNU parallel’s subshells
		MM2_TOTAL_CORES="${threads}"
		export -f _st3_process_shard
		export SHARD_DIR EDGE_DIR INDEX_OPTS MAP_OPTS MM2_TOTAL_CORES EDGE_DUMP

		note1 "[ST3] Per-shard indexing + mapping with minimap2 (-x ava-ont) [parallel]"
		# run with GNU parallel
		seq -f "%03g" 1 "$B" |
			parallel -j "${PARALLEL_JOBS:-4}" --no-notice _st3_process_shard {}
	else
		note1 "[ST3] Parallel disabled or not available; running serially"
		for ((i = 1; i <= B; i++)); do
			printf -v idx "%03d" "$i"
			_st3_process_shard "$idx"
		done
	fi

	# 4) concat per-shard TSV -> one big edges.tsv.gz
	local ALL="${EDGE_DIR}/all.edges.tsv.gz"
	note1 "[ST3] Concatenating shard edges -> $ALL"
	# Using zcat (works on .gz files)
	ls "$EDGE_DIR"/shard_*.edges.tsv.gz | sort -V |
		xargs zcat |
		gzip -c >"$ALL"

	local EDGES_COMBINED="$STDIR/edges_loose.tsv.gz"
	# ln -sf "02-edges/all.edges.tsv.gz" "$EDGES_COMBINED"
	mv "${EDGE_DIR}/all.edges.tsv.gz" "$EDGES_COMBINED"
	ln -sf "../edges_loose.tsv.gz" "${EDGE_DIR}/all.edges.tsv.gz"
	note1 "[ST3] combined: $EDGES_COMBINED"
}

ST3_overlapness() {
	local OUTDIR="$1"
	local READS="$2"

	local STDIR="${OUTDIR}/03-allvsall"
	local EDGES_COMBINED="$STDIR/edges_loose.tsv.gz"
	local OVL_STRICT="$STDIR/overlapness.tsv"

	note1 "[ST3] input: Read-by-Read overlapness: $EDGES_COMBINED"
	note1 "[ST3] output: Overlapness per read: $OVL_STRICT"
	# Check the input files
	exists "$EDGES_COMBINED"

	note1 "ovelapness 35%: $STDIR/edges_loose.tsv.35.txt"
	zcat "${EDGES_COMBINED}" | cut -f5 |
		Rscript -e 'x <- scan(file="stdin", quiet=TRUE); q <- quantile(x, probs=0.35, type=7); cat(q, "\n")' \
			>"$STDIR/edges_loose.tsv.35.txt"

	note1 "2e) loop over edge weights to get $OVL_STRICT from $EDGES_COMBINED"

	# One-pass eweight thresholding loop (uses to_milli/milli_to_str and calls REFILTER_EDGES each step)
	# Assumes:
	#   EDGES_COMBINED="$ST3/edges_loose.tsv.gz"
	#   OVL_STRICT="$ST3/overlapness.tsv"
	#   REFILTER_EDGES="${HELP}/polap-py-refilter-edges-to-overlapness.py"
	#   _ovl_cov prints coverage fraction (0..1)
	# Knobs (defaults if unset):
	eweight_m=$(to_milli "${eweight}")      # start at 0.900
	step_m=$(to_milli "${eweight_step}")    # decrement 0.050
	min_m=$(to_milli "${eweight_minimum}")  # floor 0.010
	tgt_m=$(to_milli "${shortlist_target}") # coverage target (0.35 -> 350)

	# Sanity checks
	if ((step_m <= 0)); then
		note0 "ERR: eweight_step must be > 0 (got '${eweight_step:-unset}')"
		exit 2
	fi

	local ew_str=$(<"$ST3/edges_loose.tsv.35.txt")
	note1 python "$REFILTER_EDGES" \
		--edges "$EDGES_COMBINED" \
		--w_floor "$ew_str" \
		--out "$OVL_STRICT"
	python "$REFILTER_EDGES" \
		--edges "$EDGES_COMBINED" \
		--w_floor "$ew_str" \
		--out "$OVL_STRICT"

	exists "$OVL_STRICT"
	note1 "_ovl_cov $OVL_STRICT $READS"
	cov_str=$(_ovl_cov "$OVL_STRICT" "$READS")

	cov_str="${cov_str:-0}"
	cov_m=$(to_milli "$cov_str")

	local cov_pct=$(_polap_lib_math-percentify $cov_str)
	local tgt_pct=$(_polap_lib_math-percentify $shortlist_target)
	# cov_pct=$(awk -v c="$cov_str" 'BEGIN{printf "%.2f", 100*c}')
	# tgt_pct=$(awk -v t="${shortlist_target:-0}" 'BEGIN{printf "%.2f", 100*t}')
	note1 "coverage=${cov_pct}% ($cov_str) target=${tgt_pct}% (eweight=$ew_str)"

	note1 "[ST3] 04-qc"
	local QC_DIR="$STDIR/04-qc"
	mkdir -p "$QC_DIR"
	local OVL_VARS="$QC_DIR/overlap_qc.vars"
	note1 Rscript --vanilla "$SCAN_QCR" --input "$OVL_STRICT" --outdir "$QC_DIR" --target 0.80 ">$OVL_VARS"
	Rscript --vanilla "$SCAN_QCR" --input "$OVL_STRICT" --outdir "$QC_DIR" --target 0.80 >"$OVL_VARS" || true
}

# 2. miniprot -t "${threads:-8}" …
#
# This runs miniprot, a protein-to-genome aligner.
# 	•	-t "${threads:-8}" -> use the number of threads in $threads, or default to 8 if $threads is unset.
# 	•	-S -> output SAM-like secondary alignment tags (depends on version; in miniprot it means “spliced alignment” mode).
# 	•	-N 3 -> keep at most 3 secondary alignments per query.
# 	•	--outc 0.4 -> minimum coverage of the query protein required to report an alignment (40%).
# 	•	"$RDIR/nonpt.sample.mpi" -> the genome index built by miniprot -d.
# 	•	"$nprot" -> the input protein database (e.g. BUSCO protein set).
#
# 3. >"$RDIR/nonpt.sample.busco.paf"
# 	•	Redirect stdout into a file named nonpt.sample.busco.paf under $RDIR.
# 	•	This file contains the alignment results in PAF format (Pairwise mApping Format), because miniprot outputs PAF by default.
ST4_busco() {
	local OUTDIR="$1"
	local READS="$2"

	local ST3="$OUTDIR/03-allvsall"
	local BDIR="$OUTDIR/04-busco"

	if ! [[ -n "${nprot:-}" && -s "${nprot}" ]]; then
		note0 "No busco reference protein sequence dataset"
	fi

	note1 "BUSCO QC on a 10% length-stratified subsample -> $BDIR/nuc.ids.sample"
	local SAMP_FQ="$BDIR/reads.nonpt.sample.fq.gz"
	if ((!DRY)); then
		note2 "$READS -> $SAMP_FQ min(10% or 1 Gb)"
		sample_min_10pct_or_1gb "$READS" "$SAMP_FQ" "${seed:-13}"

		note2 "[MP] $SAMP_FQ -> $BDIR/nonpt.sample.mpi"
		miniprot -d "$BDIR/nonpt.sample.mpi" "$SAMP_FQ"

		note2 "[MP] nonpt.sample.mpi vs BUSCO -> nonpt.sample.busco.paf"
		note2 miniprot -t "${TOTAL_CORES}" -S -N 3 --outc 0.4 \
			"$BDIR/nonpt.sample.mpi" "$nprot" \
			">$BDIR/nonpt.sample.busco.paf"
		miniprot -t "${TOTAL_CORES}" -S -N 3 --outc 0.4 \
			"$BDIR/nonpt.sample.mpi" "$nprot" \
			>"$BDIR/nonpt.sample.busco.paf"

		# awk 'BEGIN{FS=OFS="\t"} NF>=12 && $11+0>=150 {print $6}' \
		# 	"$BDIR/nonpt.sample.busco.paf" |
		# 	sort -u >"$BDIR/nuc.ids.sample"

		local NUC_SAMPLE="$BDIR/nuc.ids.sample"
		note2 "[PAF] $BDIR/nonpt.sample.busco.paf -> $NUC_SAMPLE"
		while IFS=$'\t' read -r -a f; do
			((${#f[@]} >= 12)) || continue
			[[ ${f[10]} =~ ^[0-9]+$ ]] && ((f[10] >= 150)) || continue # col11
			printf '%s\n' "${f[5]}"                                    # col6
		done <"$BDIR/nonpt.sample.busco.paf" | LC_ALL=C sort -u >"$NUC_SAMPLE"
	fi

}

# ST4: BUSCO-guided QC + selection round
# Usage:
#   ST4_round "$outdir" ["$st3_dir"]
# Notes:
#   Expects env/tools/vars already defined in caller:
#     R1, nprot, seed, threads, DRY, nuc_ids_opt, THRESH_NUC_PY, TOPPY, note1
#   Defaults ST3 to "$outdir/03-allvsall" if not provided.
ST5_round() {
	local outdir="$1"

	local ST3="$outdir/03-allvsall"
	local RDIR="$outdir/05-round"

	local OTSV="$ST3/overlapness.tsv"
	local SELECT_IDS="$RDIR/select_ids.txt"

	note1 "3b) selection (prefer BUSCO sample labels) -> $SELECT_IDS"
	# [[ -s "$RDIR/overlapness.tsv" ]] || cp -f "$OTSV_GLOBAL" "$RDIR/overlapness.tsv"

	local NUC_SAMPLE="$BDIR/nuc.ids.sample"
	local NUC_FOR_THRESH="${nuc_ids_opt:-"$NUC_SAMPLE"}"

	if [[ -s "$NUC_FOR_THRESH" ]]; then
		# THRESH_NUC_PY="${HELP}/polap-py-threshold-from-nuclear.py"
		note1 "nuclear-guided selection using $NUC_FOR_THRESH"
		note1 python "$THRESH_NUC_PY" "$OTSV" "$NUC_FOR_THRESH" \
			--mode fpr --fpr 0.20 \
			--diag "$RDIR/threshold_from_nuclear.tsv" \
			">$RDIR/nuc_thresh.vars"
		python "$THRESH_NUC_PY" "$OTSV" "$NUC_FOR_THRESH" \
			--mode fpr --fpr 0.20 \
			--diag "$RDIR/threshold_from_nuclear.tsv" \
			>"$RDIR/nuc_thresh.vars"

		# Get the minimum overlapness values.
		local wdeg_min=0 deg_min=0 kv
		while IFS== read -r kv; do
			case "$kv" in
			wdeg_min=*) wdeg_min="${kv#wdeg_min=}" ;;
			deg_min=*) deg_min="${kv#deg_min=}" ;;
			esac
		done <"$RDIR/nuc_thresh.vars"

		# mtseed/keep.nonpt.ids -> number of nuclear and mito reads
		# mtseed/05-round/select_ids.txt -> number of mito reads
		awk -v W="$wdeg_min" -v D="$deg_min" \
			'BEGIN{FS=OFS="\t"} NR>1 && ($3+0)>W && ($2+0)>D {print $1}' \
			"$OTSV" | sort -u >"$RDIR/organelle.ids"

		comm -23 <(sort -u "$RDIR/organelle.ids") <(sort -u "$NUC_FOR_THRESH") >"$SELECT_IDS"
	else
		note1 "No BUSCO; use $top_frac"
		note1 "select top-frac by wdegree ($top_frac)"
		((!DRY)) && python "$TOPPY" "$OTSV" --top_frac "$top_frac" >"$SELECT_IDS"
	fi
}

# ST6: Assembly with miniasm (or raven fallback)
# Usage:
#   ST6_miniasm "$outdir" "$assembler"
# Notes:
#   Expects env/tools already defined:
#     R1, threads, DRY, assembler, seqkit, minimap2, miniasm, raven, etc.
#   Paths used inside:
#     $RDIR/   : round dir from ST4
#     $ADIR/   : assembly dir under $outdir
#     $SELECT_IDS : ids selected in ST4_round

ST6_miniasm() {
	local OUTDIR="$1"
	local READS="$2"
	local assembler="${3:-miniasm}"

	local RDIR="$OUTDIR/05-round"
	local ADIR="$OUTDIR/06-miniasm"

	local SELECT_IDS="$RDIR/select_ids.txt"
	local TOPFA="$ADIR/top_reads.fa.gz"
	local TOPFQ="$ADIR/top_reads.fq.gz"
	local GFA SEEDS_ROUND
	local SEEDS_CUR

	if [[ "$assembler" == "miniasm" ]]; then

		note1 "4a) extract selected reads -> $TOPFA"
		seqkit fq2fa "$READS" | seqkit grep -f "$SELECT_IDS" -o "$TOPFA"
		seqkit grep -f "$SELECT_IDS" "$READS" -o "$TOPFQ"

		GFA="$ADIR/miniasm.gfa"

		# recompute overlaps among selected reads
		local SELPAF="$ADIR/selected_allvsall.paf.gz"
		note1 "4b) recompute overlaps among selected reads -> $SELPAF"
		minimap2 -x ava-ont -t "$TOTAL_CORES" \
			"${MINIASM_MAP_OPTS}" \
			"$TOPFA" "$TOPFA" | gzip -1 >"$SELPAF"

		miniasm -f "$TOPFA" <(zcat -f "$SELPAF") >"$GFA"

		note1 "4c) extract contigs from GFA -> seeds"
		SEEDS_ROUND="$ADIR/m_seeds_raw.fa"
		awk '/^S/{print ">"$2"\n"$3}' "$GFA" >"$SEEDS_ROUND" || :

	else
		note1 "4a) raven : extract selected reads (fq) -> $TOPFQ"
		seqkit grep -f "$SELECT_IDS" "$READS" -o "$TOPFQ"

		SEEDS_ROUND="$ADIR/raven_contigs.fasta"
		note1 "raven --disable-polishing -> $SEEDS_ROUND"
		raven --threads "$threads" --disable-polishing -o "$SEEDS_ROUND" "$TOPFQ"
	fi
}

# ────────────────────────────────────────────────────────────────────
# Main
# ────────────────────────────────────────────────────────────────────

# ────────────────────────────────────────────────────────────────────
# Step 0: minimap2 setting
# ────────────────────────────────────────────────────────────────────
# echo "std1"
# echo "std2" 1>&2
# exists "$nprot"
note1 "Step0. choose presets for minimap2"
choose_mm2_preset
note1 "Step0. checking"
ST0_check "${outdir}"
note1 "Step0. checked"

# ────────────────────────────────────────────────────────────────────
# Step 1: PT removal (unchanged)
# ────────────────────────────────────────────────────────────────────
R1="$reads"
PANEL_DIR="$outdir/01-panel"
if _should_run 1; then
	note1 "Step1. prepare ptDNA sequences"
	ST1_pt_isoform "${outdir}"
fi

# ────────────────────────────────────────────────────────────────────
# Step 2: PT removal (unchanged)
# ────────────────────────────────────────────────────────────────────

note1 "R1 before ST2_pt_map: $R1 compare $R1.seqkit.stats.ta.tsv not .txt"
if _should_run 2; then
	note1 "Step2. map reads on ptDNA"
	ST2_pt_map "${outdir}" "$R1"
else
	R1="$outdir/reads.nonpt.fq.gz"
fi
note1 "R1 after ST2_pt_map: $R1 compare $R1.seqkit.stats.ta.tsv not .txt"

if [[ ! -s "$R1" ]]; then
	note0 "ERR missing R1 ($R1)"
	exit 3
fi

# ────────────────────────────────────────────────────────────────────
# Step 3: minimap2 mapping all vs all
# ────────────────────────────────────────────────────────────────────
ST3="$outdir/03-allvsall"
if _should_run 3; then
	mkdir -p "$ST3"
	note1 "Step3. map reads on themselves - part vs part"
	ST3_part_vs_part "$outdir" "$R1"
	ST3_overlapness "$outdir" "$R1"
fi

# ────────────────────────────────────────────────────────────────────
# Step 4: BUSCO QC (10% subsample), then selection (uses $RDIR/overlapness.tsv)
# ────────────────────────────────────────────────────────────────────
BDIR="$outdir/04-busco"
if _should_run 4; then
	note1 "Step4. Detect nuclear genes"
	mkdir -p "$BDIR"
	ST4_busco "$outdir" "$R1"
fi

# ────────────────────────────────────────────────────────────────────
# Step 5: selection (uses $RDIR/overlapness.tsv)
# ────────────────────────────────────────────────────────────────────
RDIR="$outdir/05-round"
if _should_run 5; then
	note1 "Step5. Get overlapness"
	mkdir -p "$RDIR"
	ST5_round "$outdir"
fi

# ────────────────────────────────────────────────────────────────────
# Step 6: assemble selected reads (miniasm|raven)
# ────────────────────────────────────────────────────────────────────
ADIR="$outdir/06-miniasm"
if _should_run 6; then
	note1 "Step6. assemble selected reads using miniasm"
	mkdir -p "$ADIR"
	ST6_miniasm "$outdir" "$R1"
fi

# ────────────────────────────────────────────────────────────────────
# Step 7: finisher (polap readassemble) using miniasm GFA
# ────────────────────────────────────────────────────────────────────
FDIR="$outdir/07-flye"
if _should_run 7; then
	note1 "Step7. assemble mtDNA using selected contigs"
	contigger_dir="${FDIR}/30-contigger"
	mkdir -p "${contigger_dir}"

	note2 "7a) convert miniasm gfa for flye run -> ${contigger_dir}/graph_final.gfa"
	sed 's/LN:i/dp:i/' "${ADIR}/miniasm.gfa" >"${contigger_dir}/graph_final.gfa"

	# note2 "7b) polap readassemble using miniasm seeds ..."
	# if ((do_polap == 1)); then
	# 	bash "${_POLAPLIB_DIR}/../polap.sh" readassemble annotated -o "${outdir}" -i "${FDIR_NAME}" -l "${reads}"
	# fi
fi
