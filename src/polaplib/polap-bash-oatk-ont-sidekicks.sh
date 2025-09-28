#!/usr/bin/env bash
# polap-bash-oatk-ont-sidekicks.sh
# Unified ONT preprocessing wrapper for Oatk (two-stage run):
#   Stage-0 (optional prefilters): trim, scrub, lenfilt, HMM mt-bait/drop plastid, duplex-only, rare-mask, correction
#   Stage-1: (HPC default ON) → qc,assemble,summary  → backbone
#   Bait   : backbone-based k-mer bait (bbduk|meryl) → ONT mt-subset (small & fast)
#   Stage-2: recover with lift (RLE) or polish (racon/medaka) on subset → annotate → pathfinder → summary
#
# Requirements (checked on demand):
#   Always: seqkit, awk, grep, polap-bash-oatk-ont.sh
#   Stage-0:
#     --trim              : porechop_abi
#     --scrub             : minimap2, yacrd
#     --lenfilt           : filtlong (preferred) or seqkit
#     --mt-bait|--drop-plastid : hmmannot (or hmm_annotation), nhmmscan
#     --duplex-only       : none (awk/seqkit)
#     --rare-mask         : kmc (and kmc_tools)
#     --correct-canu      : canu
#     --correct-ratatosk  : ratatosk (and short reads)
#     HPC (default ON)    : seqtk
#   Bait: bbduk.sh (default) or meryl,meryl-lookup for --bait meryl
#   Stage-2 speed: minimap2 chunking is handled in polap-bash-oatk-ont.sh via --map-chunks

set -euo pipefail

#############################################
# Logging & helpers
#############################################
: "${POLAP_LOG_LEVEL:=1}" # 0=quiet, 1=info, 2=verbose
log() {
	local lvl=$1
	shift
	[[ $POLAP_LOG_LEVEL -ge $lvl ]] && echo "$@" >&2
}
die() {
	echo "[ERR]" "$@" >&2
	exit 1
}
need() { command -v "$1" >/dev/null 2>&1 || die "missing dependency: $1"; }

#############################################
# Defaults
#############################################
READS="" # input ONT reads (fa/fq[.gz])
OUT=""   # output folder
THREADS=32

# Oatk DB & clade (for HMM bait prefiltering)
OATKDB="" # must contain v20230921
CLADE="magnoliopsida"
HMMBIN="hmmannot" # or hmm_annotation
NHMMSCAN="nhmmscan"

# Side-steps toggles (all optional)
DO_TRIM=0    # porechop_abi
DO_SCRUB=0   # yacrd chimera scrub
DO_LENFILT=0 # length/Q filter
LEN_MIN=5000
KEEP_PCT=70       # for filtlong
DO_MT_BAIT=0      # HMM-based mt bait (prefilter)
DO_DROP_PLASTID=0 # HMM-based plastid removal
DO_DUPLEX_ONLY=0  # grep header tag pattern
DUPLEX_TAG="duplex"
DO_RARE_MASK=0 # KMC rare k-mer masking (heavy)
RARE_K=31
RARE_MAX=1   # <=1 occurrence → rare
DO_CORRECT=0 # 1=Canu, 2=Ratatosk
SHORT_READS=""

# HPC (assembly-time compression)
HPC_ENABLE=1 # default ON (assembly uses HPC reads)
HPC_OUT=""   # default OUT/reads.hpc.fa

# Two-stage bait to speed mapping
BAIT_METHOD="bbduk" # bbduk|meryl (backbone-based)
BBDUK_BIN="bbduk.sh"
MERYL_K=21

# Mapping acceleration to polish/lift stage
MAP_CHUNKS=4 # split mapping reads into N chunks (0=disabled)
MM2_EXTRA="" # extra minimap2 flags (e.g., --cs=short)

# Downstream ONT runner (step-wise)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ONT_SCRIPT="${ONT_SCRIPT:-"$SCRIPT_DIR/polap-bash-oatk-ont.sh"}"

# Stage-1 assemble knobs (backbone)
K_LIST_A="251,151,121,91"
SMER=31
AARC=0.25
WEAKX=0.20
UNZIP=6
NO_EC=0 # pass --no-read-ec if 1

# Stage-2 recovery
USE_LIFT=0 # --lift to use RLE lifter, otherwise polish (racon/medaka)
RACON_ROUNDS=2
MEDAKA_MODEL="r104_e81_sup_g615"
NO_MEDAKA=0

# Misc
KEEP_INTER=1
DRYRUN=0

usage() {
	cat <<EOF
Usage:
  $(basename "$0") --reads READS.fq.gz --out OUTDIR --oatkdb /path/OatkDB [options]

Required:
  --reads FILE            ONT reads (fa/fq[.gz])
  --out   DIR             output folder
  --oatkdb DIR            OatkDB root (must have v20230921)

General:
  --threads INT           (default ${THREADS})
  -v|--verbose            verbose logs
  -q|--quiet              quiet logs
  --dry-run               print plan, do nothing

Stage-0 Prefilters (optional):
  --trim                  adapter trim (porechop_abi)
  --scrub                 chimera scrub (minimap2 ava-ont + yacrd)
  --lenfilt               length filter (filtlong preferred; else seqkit)
  --len-min INT           min length (default ${LEN_MIN})
  --keep-pct INT          filtlong keep percent (default ${KEEP_PCT})
  --mt-bait               HMM bait: keep reads hitting mito fam (pre-assembly)
  --drop-plastid          HMM drop: remove reads hitting plastid fam
  --duplex-only           keep only reads with header tag (default '${DUPLEX_TAG}'; change via --duplex-tag)
  --duplex-tag STR        header pattern to select duplex reads
  --rare-mask             rare k-mer masking (KMC) [heavy]
  --rare-k INT            K for rare masking (default ${RARE_K})
  --rare-max INT          max count to mark rare (default ${RARE_MAX})
  --correct-canu          Canu correction-only (pre-consensus)
  --correct-ratatosk      Ratatosk correction (needs --short-reads)
  --short-reads FILE      short reads for Ratatosk

HPC (assembly-time):
  --hpc                   ENABLE HPC (seqtk hpc) [default]
  --no-hpc                disable HPC
  --hpc-out PATH          HPC FASTA path (default OUT/reads.hpc.fa)

Stage-1 Assemble (backbone):
  --k-list "251,151,121,91"
  --smer INT              syncasm -s (<=31; default ${SMER})
  --a FLOAT               min arc coverage (default ${AARC})
  --weak-cross FLOAT      (default ${WEAKX})
  --unzip-round INT       (default ${UNZIP})
  --no-ec                 pass --no-read-ec

Backbone-based bait (between stages):
  --bait bbduk|meryl      method for fast mt-subset (default bbduk)
  --meryl-k INT           meryl k (default ${MERYL_K})

Stage-2 Recovery + fast mapping:
  --lift                  use RLE lifter (Option B) instead of polish (Option A)
  --racon-rounds INT      racon rounds (default ${RACON_ROUNDS})
  --medaka-model STR      medaka model (default ${MEDAKA_MODEL})
  --no-medaka             skip medaka
  --map-chunks INT        split mapping reads into INT chunks (default ${MAP_CHUNKS})
  --mm2-extra "FLAGS"     extra minimap2 flags (default "")
EOF
}

#############################################
# Parse CLI
#############################################
while [[ $# -gt 0 ]]; do
	case "$1" in
	--reads)
		READS="$2"
		shift 2
		;;
	--out)
		OUT="$2"
		shift 2
		;;
	--oatkdb)
		OATKDB="$2"
		shift 2
		;;
	--threads)
		THREADS="$2"
		shift 2
		;;
	--clade)
		CLADE="$2"
		shift 2
		;;
	-v | --verbose)
		POLAP_LOG_LEVEL=2
		shift
		;;
	-q | --quiet)
		POLAP_LOG_LEVEL=0
		shift
		;;
	--dry-run)
		DRYRUN=1
		shift
		;;

	--trim)
		DO_TRIM=1
		shift
		;;
	--scrub)
		DO_SCRUB=1
		shift
		;;
	--lenfilt)
		DO_LENFILT=1
		shift
		;;
	--len-min)
		LEN_MIN="$2"
		shift 2
		;;
	--keep-pct)
		KEEP_PCT="$2"
		shift 2
		;;
	--mt-bait)
		DO_MT_BAIT=1
		shift
		;;
	--drop-plastid)
		DO_DROP_PLASTID=1
		shift
		;;
	--duplex-only)
		DO_DUPLEX_ONLY=1
		shift
		;;
	--duplex-tag)
		DUPLEX_TAG="$2"
		shift 2
		;;
	--rare-mask)
		DO_RARE_MASK=1
		shift
		;;
	--rare-k)
		RARE_K="$2"
		shift 2
		;;
	--rare-max)
		RARE_MAX="$2"
		shift 2
		;;
	--correct-canu)
		DO_CORRECT=1
		shift
		;;
	--correct-ratatosk)
		DO_CORRECT=2
		shift
		;;
	--short-reads)
		SHORT_READS="$2"
		shift 2
		;;

	--hpc)
		HPC_ENABLE=1
		shift
		;;
	--no-hpc)
		HPC_ENABLE=0
		shift
		;;
	--hpc-out)
		HPC_OUT="$2"
		shift 2
		;;

	--ont-script)
		ONT_SCRIPT="$2"
		shift 2
		;;

	--k-list)
		K_LIST_A="$2"
		shift 2
		;;
	--smer)
		SMER="$2"
		shift 2
		;;
	--a)
		AARC="$2"
		shift 2
		;;
	--weak-cross)
		WEAKX="$2"
		shift 2
		;;
	--unzip-round)
		UNZIP="$2"
		shift 2
		;;
	--no-ec)
		NO_EC=1
		shift
		;;

	--bait)
		BAIT_METHOD="$2"
		shift 2
		;;
	--meryl-k)
		MERYL_K="$2"
		shift 2
		;;

	--lift)
		USE_LIFT=1
		shift
		;;
	--racon-rounds)
		RACON_ROUNDS="$2"
		shift 2
		;;
	--medaka-model)
		MEDAKA_MODEL="$2"
		shift 2
		;;
	--no-medaka)
		NO_MEDAKA=1
		shift
		;;

	--map-chunks)
		MAP_CHUNKS="$2"
		shift 2
		;;
	--mm2-extra)
		MM2_EXTRA="$2"
		shift 2
		;;

	-h | --help)
		usage
		exit 0
		;;
	*)
		echo "[ERR] unknown arg: $1" >&2
		usage
		exit 1
		;;
	esac
done
[[ -z "$READS" || -z "$OUT" || -z "$OATKDB" ]] && {
	usage
	exit 1
}

#############################################
# Prep
#############################################
READS="$(readlink -f "$READS")"
OUT="$(readlink -f "$OUT")"
OATKDB="$(readlink -f "$OATKDB")"
mkdir -p "$OUT"
cd "$OUT"
log 1 "[info] OUT=$OUT threads=$THREADS"

# fams for HMM-based prefilters
FAM_M="$OATKDB/v20230921/${CLADE}_mito.fam"
FAM_P="$OATKDB/v20230921/${CLADE}_pltd.fam"
[[ -s "$FAM_M" ]] || die "mito fam not found: $FAM_M"
[[ $DO_DROP_PLASTID -eq 0 || -s "$FAM_P" ]] || die "plastid fam not found: $FAM_P"

# dependency checks (on demand)
need seqkit
[[ $DO_TRIM -eq 1 ]] && need porechop_abi
[[ $DO_SCRUB -eq 1 ]] && {
	need minimap2
	need yacrd
}
[[ $DO_LENFILT -eq 1 ]] && { command -v filtlong >/dev/null 2>&1 || need seqkit; }
[[ $DO_MT_BAIT -eq 1 || $DO_DROP_PLASTID -eq 1 ]] && {
	command -v "$HMMBIN" >/dev/null 2>&1 || die "missing $HMMBIN"
	need "$NHMMSCAN"
}
[[ $DO_RARE_MASK -eq 1 ]] && {
	need kmc
	command -v kmc_tools >/dev/null 2>&1 || need kmc_tools
}
[[ $DO_CORRECT -eq 1 ]] && need canu
[[ $DO_CORRECT -eq 2 ]] && {
	need ratatosk || need Ratatosk || true
	[[ -n "$SHORT_READS" ]] || die "--short-reads required for Ratatosk"
}
[[ "$BAIT_METHOD" == "bbduk" ]] && need "$BBDUK_BIN" || { [[ "$BAIT_METHOD" == "meryl" ]] && {
	need meryl
	need meryl-lookup
}; }
[[ $HPC_ENABLE -eq 1 ]] && need seqtk
# need "$ONT_SCRIPT"

[[ $DRYRUN -eq 1 ]] && {
	log 1 "[dry-run] would run prefilters on $READS and call $ONT_SCRIPT"
	exit 0
}

#############################################
# Stage-0: prefilters
#############################################
CUR="$READS"

# 0.1 Trim
if [[ $DO_TRIM -eq 1 ]]; then
	log 1 "[trim] porechop_abi"
	porechop_abi -i "$CUR" -o reads.trim.fq -t "$THREADS" || die "porechop_abi failed"
	CUR="$OUT/reads.trim.fq"
fi

# 0.2 Chimera scrub
if [[ $DO_SCRUB -eq 1 ]]; then
	log 1 "[scrub] minimap2 ava-ont + yacrd"
	minimap2 -x ava-ont -g 500 -t "$THREADS" "$CUR" "$CUR" >overlaps.paf
	yacrd -i overlaps.paf -o yacrd.tsv -c 4 -n -p 0.4
	awk '$3=="chimera"{next} $3=="unchanged"{print $1}' yacrd.tsv >keep.ids
	seqkit grep -f keep.ids "$CUR" -o reads.scrub.fq
	CUR="$OUT/reads.scrub.fq"
fi

# 0.3 Length/Q filter
if [[ $DO_LENFILT -eq 1 ]]; then
	if command -v filtlong >/dev/null 2>&1; then
		log 1 "[lenfilt] filtlong --min_length $LEN_MIN --keep_percent $KEEP_PCT"
		filtlong --min_length "$LEN_MIN" --keep_percent "$KEEP_PCT" "$CUR" >reads.long.fq
		CUR="$OUT/reads.long.fq"
	else
		log 1 "[lenfilt] seqkit seq -m $LEN_MIN"
		seqkit seq -m "$LEN_MIN" "$CUR" -o reads.len.fq
		CUR="$OUT/reads.len.fq"
	fi
fi

# 0.4 HMM mt-bait (prefilter)
if [[ $DO_MT_BAIT -eq 1 ]]; then
	log 1 "[mt-bait pre] $HMMBIN vs mito fam"
	"$HMMBIN" -t "$THREADS" --nhmmscan "$NHMMSCAN" -o reads.mito.txt "$FAM_M" "$CUR"
	awk -v FS='[[:space:]]+' '$1!~/^#/{print $3}' reads.mito.txt | sort -u >mt.ids
	seqkit grep -f mt.ids "$CUR" -o reads.mt.fq
	CUR="$OUT/reads.mt.fq"
fi

# 0.5 Plastid down-select (remove plastid-like reads)
if [[ $DO_DROP_PLASTID -eq 1 ]]; then
	log 1 "[plastid pre] dropping reads hitting plastid fam"
	"$HMMBIN" -t "$THREADS" --nhmmscan "$NHMMSCAN" -o reads.pt.txt "$FAM_P" "$CUR"
	awk -v FS='[[:space:]]+' '$1!~/^#/{print $3}' reads.pt.txt | sort -u >pt.ids
	seqkit grep -v -f pt.ids "$CUR" -o reads.mt.noPT.fq
	CUR="$OUT/reads.mt.noPT.fq"
fi

# 0.6 Duplex-only
if [[ $DO_DUPLEX_ONLY -eq 1 ]]; then
	log 1 "[duplex] selecting reads with tag: $DUPLEX_TAG"
	awk -v tag="$DUPLEX_TAG" 'BEGIN{OFS="\n"} NR%4==1{keep=index($0,tag)>0;h=$0} NR%4==2{s=$0} NR%4==3{p=$0} NR%4==0{q=$0; if(keep){print h,s,p,q}}' "$CUR" >reads.duplex.fq
	CUR="$OUT/reads.duplex.fq"
fi

# 0.7 Rare k-mer masking (heavy; placeholder pass-through)
if [[ $DO_RARE_MASK -eq 1 ]]; then
	log 1 "[rare-k] KMC on k=${RARE_K} ; rare <= ${RARE_MAX}"
	mkdir -p kmc_tmp
	kmc -k"${RARE_K}" -t"$THREADS" -ci1 -cs100000000 "$CUR" kmc.db kmc_tmp/
	kmc_tools transform kmc.db dump -ci1 -cx${RARE_MAX} rare_k.txt
	log 1 "[rare-k] (placeholder) pass-through reads → reads.mask.fq"
	cp "$CUR" reads.mask.fq
	CUR="$OUT/reads.mask.fq"
fi

# 0.8 Pre-correction (Canu or Ratatosk)
if [[ $DO_CORRECT -eq 1 ]]; then
	log 1 "[correct] Canu correction-only"
	mkdir -p corr
	canu -correct -p ontcorr -d corr genomeSize=500k -nanopore-raw "$CUR" 2>corr/canu.log || log 1 "[warn] canu correction failed?"
	if ls corr/ontcorr.correctedReads.fasta.gz >/dev/null 2>&1; then
		CUR="$OUT/corr/ontcorr.correctedReads.fasta.gz"
	fi
elif [[ $DO_CORRECT -eq 2 ]]; then
	log 1 "[correct] Ratatosk"
	need ratatosk || need Ratatosk || true
	[[ -n "$SHORT_READS" ]] || die "--short-reads required for Ratatosk"
	ratatosk -s "$SHORT_READS" -l "$CUR" -o reads.rtk.fq -t "$THREADS" || log 1 "[warn] ratatosk failed?"
	CUR="$OUT/reads.rtk.fq"
fi

# 0.9 HPC for assembly
HPC_FLAG=(--hpc)
if [[ $HPC_ENABLE -eq 1 ]]; then
	need seqtk
	HPC_FILE="${HPC_OUT:-$OUT/reads.hpc.fa}"
	log 1 "[hpc] seqtk hpc → $HPC_FILE"
	seqtk hpc "$CUR" >"$HPC_FILE"
	CUR_ASM="$HPC_FILE"
	HPC_FLAG=(--hpc --hpc-out "$HPC_FILE")
else
	log 1 "[hpc] disabled; assemble on prefiltered reads"
	CUR_ASM="$CUR"
	HPC_FLAG=(--no-hpc)
fi

#############################################
# Stage-1: assemble (backbone only)
#############################################
log 1 "[stage1] assemble backbone (k=${K_LIST_A}; smer=${SMER}; -a=${AARC})"
bash "$ONT_SCRIPT" \
	--reads "$CUR_ASM" \
	--out "$OUT/ont_stage1" \
	--oatkdb "$OATKDB" --clade "$CLADE" \
	--threads "$THREADS" \
	"${HPC_FLAG[@]}" \
	--k-list "$K_LIST_A" --smer "$SMER" --a "$AARC" --weak-cross "$WEAKX" --unzip-round "$UNZIP" \
	$([[ $NO_EC -eq 1 ]] && echo --no-read-ec) \
	-s qc,assemble,summary -v

BACKBONE="$OUT/ont_stage1/k1/unitigs.fa"
[[ -s "$BACKBONE" ]] || die "backbone unitigs not found: $BACKBONE"

#############################################
# Bait: backbone-based k-mer bait → subset
#############################################
mkdir -p bait
SUBSET="bait/ont.mt.fq"

if [[ "$BAIT_METHOD" == "bbduk" ]]; then
	need "$BBDUK_BIN"
	log 1 "[bait] bbduk: ref=$BACKBONE → $SUBSET"
	"$BBDUK_BIN" in="$CUR" outm="$SUBSET" outu=bait/nonmt.fq ref="$BACKBONE" k=31 hdist=1 threads="$THREADS"
else
	need meryl
	need meryl-lookup
	log 1 "[bait] meryl k=$MERYL_K"
	meryl k="$MERYL_K" count "$BACKBONE" output bait/seeds.k"$MERYL_K".meryl
	meryl-lookup -existence -sequence "$CUR" bait/seeds.k"$MERYL_K".meryl |
		awk '/^>/{print substr($0,2)}' >bait/mt.ids
	seqkit grep -f bait/mt.ids "$CUR" -o "$SUBSET"
fi

log 1 "[bait] subset size: $( (wc -c <"$SUBSET") 2>/dev/null || echo 0) bytes"

#############################################
# Stage-2: lift OR polish → annotate → PF → summary
#############################################
if [[ $USE_LIFT -eq 1 ]]; then
	# Lifter (RLE) + 1× polish after lift
	log 1 "[stage2] LIFT (RLE) + 1× polish (fast), then annotate & PF"
	bash "$ONT_SCRIPT" \
		--reads "$CUR_ASM" \
		--reads-map "$SUBSET" \
		--map-chunks "$MAP_CHUNKS" --mm2-extra "$MM2_EXTRA" \
		--out "$OUT/ont_stage2" \
		--oatkdb "$OATKDB" --clade "$CLADE" \
		--threads "$THREADS" \
		--lift --polish-after-lift \
		--racon-rounds 1 \
		$([[ $NO_MEDAKA -eq 1 ]] && echo --no-medaka || echo --medaka-model "$MEDAKA_MODEL") \
		-s lift,polish,annotate,pathfinder,summary -v
else
	# Classic polish on subset
	log 1 "[stage2] POLISH (racon/medaka) on subset, then annotate & PF"
	bash "$ONT_SCRIPT" \
		--reads "$CUR_ASM" \
		--reads-map "$SUBSET" \
		--map-chunks "$MAP_CHUNKS" --mm2-extra "$MM2_EXTRA" \
		--out "$OUT/ont_stage2" \
		--oatkdb "$OATKDB" --clade "$CLADE" \
		--threads "$THREADS" \
		--racon-rounds "$RACON_ROUNDS" \
		$([[ $NO_MEDAKA -eq 1 ]] && echo --no-medaka || echo --medaka-model "$MEDAKA_MODEL") \
		-s polish,annotate,pathfinder,summary -v
fi

#############################################
# Side-report
#############################################
log 1 "[EMIT] sidekicks report"
{
	echo -e "key\tvalue"
	echo -e "reads_in\t$READS"
	echo -e "out\t$OUT"
	echo -e "hpc\t$([[ $HPC_ENABLE -eq 1 ]] && echo ON || echo OFF)"
	echo -e "trim\t$([[ $DO_TRIM -eq 1 ]] && echo ON || echo OFF)"
	echo -e "scrub\t$([[ $DO_SCRUB -eq 1 ]] && echo ON || echo OFF)"
	echo -e "len_min\t$([[ $DO_LENFILT -eq 1 ]] && echo $LEN_MIN || echo .)"
	echo -e "keep_pct\t$([[ $DO_LENFILT -eq 1 && $(
		command -v filtlong >/dev/null
		echo $?
	) -eq 0 ]] && echo $KEEP_PCT || echo .)"
	echo -e "mt_bait_prefilter\t$([[ $DO_MT_BAIT -eq 1 ]] && echo ON || echo OFF)"
	echo -e "drop_plastid_prefilter\t$([[ $DO_DROP_PLASTID -eq 1 ]] && echo ON || echo OFF)"
	echo -e "duplex_only\t$([[ $DO_DUPLEX_ONLY -eq 1 ]] && echo $DUPLEX_TAG || echo OFF)"
	echo -e "rare_mask\t$([[ $DO_RARE_MASK -eq 1 ]] && echo K=$RARE_K,max=$RARE_MAX || echo OFF)"
	echo -e "corrector\t$([[ $DO_CORRECT -eq 1 ]] && echo CANU || ([[ $DO_CORRECT -eq 2 ]] && echo RATATOSK || echo OFF))"
	echo -e "bait_method\t$BAIT_METHOD"
	echo -e "map_chunks\t$MAP_CHUNKS"
	echo -e "mm2_extra\t${MM2_EXTRA:-.}"
	echo -e "ont_script\t$ONT_SCRIPT"
	echo -e "k_list_stage1\t$K_LIST_A"
	echo -e "racon_rounds\t$RACON_ROUNDS"
	echo -e "medaka_model\t$([[ $NO_MEDAKA -eq 1 ]] && echo OFF || echo $MEDAKA_MODEL)"
} >"$OUT/sidekicks.report.tsv"

log 1 "[done] Two-stage ONT run completed."
log 1 "  Stage-1 backbone : $OUT/ont_stage1/k1/unitigs.fa"
log 1 "  Stage-2 outputs  : $OUT/ont_stage2 (polished/lifted, annotated, pathfinder, summary)"
