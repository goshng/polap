#!/usr/bin/env bash
# polap-bash-oatk-ont-s1.sh
# Stage-1 ONT prefilters for organelle assembly:
#  trim → scrub → lenfilt → HMM (mt-bait / pt-drop) → duplex → rare-k (kmc) → correction
#  NEW: filtlong ref-bias + add-back of non-mapping reads to guard junctions
set -euo pipefail

: "${POLAP_LOG_LEVEL:=1}" # 0=quiet,1=info,2=debug
log() {
	local l=$1
	shift
	[[ $POLAP_LOG_LEVEL -ge $l ]] && echo "[INFO]" "$@" >&2
}
die() {
	echo "[ERR]" "$@" >&2
	exit 1
}
need() { command -v "$1" >/dev/null 2>&1 || die "missing dependency: $1"; }

# ───────── CLI defaults ─────────
READS=""
OUT=""
OATKDB=""
CLADE="magnoliopsida"
THREADS=32
HMMBIN="hmmannot"
NHMMSCAN="nhmmscan"

DO_TRIM=0
DO_SCRUB=0
DO_LENFILT=0
LEN_MIN=5000
KEEP_PCT=70
DO_MT_BAIT=0
DO_DROP_PLASTID=0
DO_DUPLEX_ONLY=0
DUPLEX_TAG="duplex"
DO_RARE_MASK=0
RARE_K=31
RARE_MAX=1
DO_CORRECT=0
SHORT_READS=""

# filtlong ref-bias + add-back
FILT_REF=""
FILT_TB=""
FILT_KP=""
FILT_MINLEN=3000
FILT_LW=1.0
FILT_QW=1.0
FILT_ADDBACK=5

usage() {
	cat <<EOF
Usage:
  $(basename "$0") --reads READS.fq.gz --out OUTDIR --oatkdb OATKDB [options]

Required:
  --reads FILE            ONT reads (fa/fq[.gz])
  --out   DIR             output directory
  --oatkdb DIR            OatkDB root (must contain v20230921/)

General:
  --threads INT           (default ${THREADS})
  --clade NAME            fam prefix (default ${CLADE})
  --hmm-bin PATH          hmmannot/hmm_annotation (default ${HMMBIN})
  --nhmmscan PATH         nhmmscan (default ${NHMMSCAN})
  -v|--verbose | -q|--quiet

Stage-1 filters (all optional):
  --trim                  porechop_abi
  --scrub                 minimap2 ava-ont + yacrd
  --lenfilt               filtlong (preferred) or seqkit
  --len-min INT           (default ${LEN_MIN})
  --keep-pct INT          kept percent (default ${KEEP_PCT})
  --mt-bait               keep reads hitting mito fam (HMM)
  --drop-plastid          drop reads hitting plastid fam (HMM)
  --duplex-only           keep headers containing tag (default '${DUPLEX_TAG}')
  --duplex-tag STR
  --rare-mask             KMC rare-k masking (placeholder)
  --rare-k INT            (default ${RARE_K})
  --rare-max INT          (default ${RARE_MAX})
  --correct-canu          Canu correction-only
  --correct-ratatosk      Ratatosk (needs --short-reads)
  --short-reads FILE

filtlong ref-bias + add-back:
  --filtlong-ref PATH
  --filtlong-target-bases N   (mutually exclusive with --filtlong-keep-percent)
  --filtlong-keep-percent P
  --filtlong-min-length N     (default ${FILT_MINLEN})
  --filtlong-length-weight W  (default ${FILT_LW})
  --filtlong-meanq-weight W   (default ${FILT_QW})
  --filtlong-addback-percent P  (default ${FILT_ADDBACK})
EOF
}

# ───────── parse CLI ─────────
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
	--hmm-bin)
		HMMBIN="$2"
		shift 2
		;;
	--nhmmscan)
		NHMMSCAN="$2"
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

	--filtlong-ref)
		FILT_REF="$2"
		shift 2
		;;
	--filtlong-target-bases)
		FILT_TB="$2"
		shift 2
		;;
	--filtlong-keep-percent)
		FILT_KP="$2"
		shift 2
		;;
	--filtlong-min-length)
		FILT_MINLEN="$2"
		shift 2
		;;
	--filtlong-length-weight)
		FILT_LW="$2"
		shift 2
		;;
	--filtlong-meanq-weight)
		FILT_QW="$2"
		shift 2
		;;
	--filtlong-addback-percent)
		FILT_ADDBACK="$2"
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
if [[ -n "$FILT_REF" && -z "$FILT_TB" && -z "$FILT_KP" ]]; then
	die "with --filtlong-ref, specify one of --filtlong-target-bases or --filtlong-keep-percent"
fi

# ───────── prep & deps ─────────
READS="$(readlink -f "$READS")"
OUT="$(readlink -f "$OUT")"
OATKDB="$(readlink -f "$OATKDB")"
mkdir -p "$OUT/pre"
cd "$OUT/pre"

need seqkit
[[ $DO_TRIM -eq 1 ]] && need porechop_abi
[[ $DO_SCRUB -eq 1 ]] && {
	need minimap2
	need yacrd
}
[[ $DO_LENFILT -eq 1 ]] && { command -v filtlong >/dev/null 2>&1 || need seqkit; }
FAM_M="$OATKDB/v20230921/${CLADE}_mito.fam"
FAM_P="$OATKDB/v20230921/${CLADE}_pltd.fam"
[[ $DO_MT_BAIT -eq 0 || -s "$FAM_M" ]] || die "mito fam not found: $FAM_M"
[[ $DO_DROP_PLASTID -eq 0 || -s "$FAM_P" ]] || die "plastid fam not found: $FAM_P"
[[ $DO_MT_BAIT -eq 0 && $DO_DROP_PLASTID -eq 0 ]] || {
	command -v "$HMMBIN" >/dev/null 2>&1 || die "missing $HMMBIN"
	need "$NHMMSCAN"
}
if [[ -n "$FILT_REF" ]]; then
	need filtlong
	need minimap2
	FILT_REF="$(readlink -f "$FILT_REF")"
fi
[[ $DO_RARE_MASK -eq 1 ]] && {
	need kmc
	need kmc_tools
}

# ───────── Stage-1 run ─────────
CUR="$READS"

# (1) trim
if [[ $DO_TRIM -eq 1 ]]; then
	log 1 "[trim] porechop_abi"
	porechop_abi -i "$CUR" -o reads.trim.fq -t "$THREADS"
	CUR="$PWD/reads.trim.fq"
fi

# (2) scrub (ava-ont + yacrd)
if [[ $DO_SCRUB -eq 1 ]]; then
	log 1 "[scrub] minimap2 ava-ont + yacrd"
	minimap2 -x ava-ont -g 500 -t "$THREADS" "$CUR" "$CUR" >overlaps.paf
	yacrd -i overlaps.paf -o yacrd.tsv -c 4 -n -p 0.4
	awk '$3=="chimera"{next} $3=="unchanged"{print $1}' yacrd.tsv >keep.ids
	seqkit grep -f keep.ids "$CUR" -o reads.scrub.fq
	CUR="$PWD/reads.scrub.fq"
fi

# (3) length filter
if [[ $DO_LENFILT -eq 1 ]]; then
	if command -v filtlong >/dev/null 2>&1; then
		log 1 "[lenfilt] filtlong --min_length $LEN_MIN --keep_percent $KEEP_PCT"
		filtlong --min_length "$LEN_MIN" --keep_percent "$KEEP_PCT" "$CUR" >reads.long.fq
		CUR="$PWD/reads.long.fq"
	else
		log 1 "[lenfilt] seqkit seq -m $LEN_MIN"
		seqkit seq -m "$LEN_MIN" "$CUR" -o reads.len.fq
		CUR="$PWD/reads.len.fq"
	fi
fi

# (4) HMM pre-bait / pt-drop
if [[ $DO_MT_BAIT -eq 1 ]]; then
	log 1 "[mt-bait] $HMMBIN vs mito fam"
	"$HMMBIN" -t "$THREADS" --nhmmscan "$NHMMSCAN" -o reads.mito.txt "$FAM_M" "$CUR"
	awk -v FS='[[:space:]]+' '$1!~/^#/{print $3}' reads.mito.txt | LC_ALL=C sort -u >mt.ids
	seqkit grep -f mt.ids "$CUR" -o reads.mt.fq
	CUR="$PWD/reads.mt.fq"
fi
if [[ $DO_DROP_PLASTID -eq 1 ]]; then
	log 1 "[pt-drop] drop plastid-like reads"
	"$HMMBIN" -t "$THREADS" --nhmmscan "$NHMMSCAN" -o reads.pt.txt "$FAM_P" "$CUR"
	awk -v FS='[[:space:]]+' '$1!~/^#/{print $3}' reads.pt.txt | LC_ALL=C sort -u >pt.ids
	seqkit grep -v -f pt.ids "$CUR" -o reads.mt.noPT.fq
	CUR="$PWD/reads.mt.noPT.fq"
fi

# (5) duplex-only
if [[ $DO_DUPLEX_ONLY -eq 1 ]]; then
	log 1 "[duplex] keep header tag: $DUPLEX_TAG"
	awk -v tag="$DUPLEX_TAG" 'BEGIN{OFS="\n"}
      NR%4==1{keep=index($0,tag)>0;h=$0}
      NR%4==2{s=$0}
      NR%4==3{p=$0}
      NR%4==0{q=$0; if(keep){print h,s,p,q}}' "$CUR" >reads.duplex.fq
	CUR="$PWD/reads.duplex.fq"
fi

# (6) rare-k masking (placeholder)
if [[ $DO_RARE_MASK -eq 1 ]]; then
	log 1 "[rare-k] KMC k=$RARE_K, <=$RARE_MAX (placeholder)"
	mkdir -p kmc_tmp
	kmc -k"$RARE_K" -t"$THREADS" -ci1 -cs100000000 "$CUR" kmc.db kmc_tmp/
	kmc_tools transform kmc.db dump -ci1 -cx"$RARE_MAX" rare_k.txt
	cp "$CUR" reads.mask.fq
	CUR="$PWD/reads.mask.fq"
fi

# (7) correction
if [[ $DO_CORRECT -eq 1 ]]; then
	log 1 "[correct] Canu correction-only"
	mkdir -p corr
	canu -correct -p ontcorr -d corr genomeSize=500k -nanopore-raw "$CUR" 2>corr/canu.log || log 1 "[warn] canu failed?"
	if ls corr/ontcorr.correctedReads.fasta.gz >/dev/null 2>&1; then CUR="$PWD/corr/ontcorr.correctedReads.fasta.gz"; fi
elif [[ $DO_CORRECT -eq 2 ]]; then
	log 1 "[correct] Ratatosk"
	need ratatosk || need Ratatosk || true
	[[ -n "$SHORT_READS" ]] || die "--short-reads required for Ratatosk"
	ratatosk -s "$SHORT_READS" -l "$CUR" -o reads.rtk.fq -t "$THREADS" || log 1 "[warn] ratatosk failed?"
	CUR="$PWD/reads.rtk.fq"
fi

# (8) filtlong ref-bias + add-back
FINAL="$PWD/reads.pre.fq"
if [[ -n "$FILT_REF" ]]; then
	need filtlong
	need minimap2
	log 1 "[filtlong] ref=$FILT_REF min_len=$FILT_MINLEN lw=$FILT_LW qw=$FILT_QW"
	if [[ -n "$FILT_TB" ]]; then
		filtlong --ref "$FILT_REF" --min_length "$FILT_MINLEN" \
			--target_bases "$FILT_TB" --length_weight "$FILT_LW" --mean_q_weight "$FILT_QW" \
			"$CUR" >pre.mtbiased.fq
	else
		filtlong --ref "$FILT_REF" --min_length "$FILT_MINLEN" \
			--keep_percent "$FILT_KP" --length_weight "$FILT_LW" --mean_q_weight "$FILT_QW" \
			"$CUR" >pre.mtbiased.fq
	fi
	log 1 "[filtlong] add-back: collect non-mapping reads"
	minimap2 -x map-ont -t "$THREADS" "$FILT_REF" "$CUR" >pre.map.paf
	awk '{if($1!~/^@/ && NF>=12 && $12+0>=1) print $1}' pre.map.paf | LC_ALL=C sort -u >mapped.ids
	seqkit seq -n "$CUR" | LC_ALL=C sort -u >all.ids
	comm -23 all.ids mapped.ids >nomap.ids
	if [[ -s nomap.ids && "$FILT_ADDBACK" -gt 0 ]]; then
		seqkit grep -f nomap.ids "$CUR" -o nomap.fq
		filtlong --keep_percent "$FILT_ADDBACK" nomap.fq >nomap.best.fq
		cat pre.mtbiased.fq nomap.best.fq >"$FINAL"
	else
		cp -f pre.mtbiased.fq "$FINAL"
	fi
else
	cp -f "$CUR" "$FINAL"
fi

# report
{
	echo -e "key\tvalue"
	echo -e "reads_in\t$READS"
	echo -e "reads_pre\t$FINAL"
	echo -e "trim\t$([[ $DO_TRIM -eq 1 ]] && echo ON || echo OFF)"
	echo -e "scrub\t$([[ $DO_SCRUB -eq 1 ]] && echo ON || echo OFF)"
	echo -e "len_min\t$([[ $DO_LENFILT -eq 1 ]] && echo $LEN_MIN || echo .)"
	if [[ $DO_LENFILT -eq 1 ]] && command -v filtlong >/dev/null 2>&1; then
		echo -e "keep_pct\t${KEEP_PCT}"
	else
		echo -e "keep_pct\t."
	fi
	echo -e "mt_bait_prefilter\t$([[ $DO_MT_BAIT -eq 1 ]] && echo ON || echo OFF)"
	echo -e "drop_plastid_prefilter\t$([[ $DO_DROP_PLASTID -eq 1 ]] && echo ON || echo OFF)"
	echo -e "duplex_only\t$([[ $DO_DUPLEX_ONLY -eq 1 ]] && echo $DUPLEX_TAG || echo OFF)"
	echo -e "rare_mask\t$([[ $DO_RARE_MASK -eq 1 ]] && echo K=$RARE_K,max=$RARE_MAX || echo OFF)"
	echo -e "corrector\t$([[ $DO_CORRECT -eq 1 ]] && echo CANU || ([[ $DO_CORRECT -eq 2 ]] && echo RATATOSK || echo OFF))"
	echo -e "filtlong_ref\t${FILT_REF:-.}"
	echo -e "filtlong_tb\t${FILT_TB:-.}"
	echo -e "filtlong_kp\t${FILT_KP:-.}"
	echo -e "filtlong_minlen\t${FILT_MINLEN}"
	echo -e "filtlong_lw\t${FILT_LW}"
	echo -e "filtlong_qw\t${FILT_QW}"
	echo -e "filtlong_addback\t${FILT_ADDBACK}%"
} >"$OUT/pre/stage1.report.tsv"

log 1 "[Stage-1 done] Preprocessed reads: $FINAL"
