#!/usr/bin/env bash
# polap-bash-oatk-core-ont.sh v0.0.5
# Core ONT/HIFI k-ladder assembly using syncasm + judge per k.
# - Captures syncasm stderr and parses singleton % (smer/kmer) into qc/summary.k*.txt
# - Extracts unitigs from GFA (prefers .utg.final.gfa, falls back to .utg.gfa)
# - Calls judge with --gfa-final/--gfa-utg; records SUMMARY/JSON/TSV
# - Appends rich rows to qc/assemble.out even if NO PASS

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

: "${POLAP_LOG_LEVEL:=1}" # 0=quiet,1=info,2=debug
log() {
	local l=$1
	shift
	[[ $POLAP_LOG_LEVEL -ge $l ]] && echo "[INFO]" "$@" >&2
}
dbg() {
	local l=$1
	shift
	[[ $POLAP_LOG_LEVEL -ge $l ]] && echo "[DBG]" "$@" >&2
}
die() {
	echo "[ERR]" "$@" >&2
	exit 1
}
need() { command -v "$1" >/dev/null 2>&1 || die "missing dependency: $1"; }
TIME_BIN="$(command -v /usr/bin/time || command -v gtime || echo time)"

# ---------------- CLI defaults ----------------
READS=""
OUT=""
THREADS=32
K_LIST="401,251,151,121"
SMER=31
NO_READ_EC=0
AARC=""
MAX_BUB="100000"
MAX_TIP=""
WEAK_X="0.20"
UNZIP_ROUND="4"
COV=""
MT_SIZE_EST=0

# Judge thresholds (tunable)
J_MAX_UNITIGS=80
J_MIN_TOTAL_BP=40000
J_TANGLE_MAX=0.06
J_MIN_N50=8000
J_MIN_MAXLEN=12000
J_CIRC_MIN_OVLP=1000
J_CIRC_MIN_IDENT=0.95

usage() {
	cat <<EOF
polap-bash-oatk-core-ont.sh v0.0.5
Usage:
  $(basename "$0") --reads READS.fq.gz --out OUTDIR [options]

Required:
  --reads FILE                 reads (fa/fq[.gz])
  --out   DIR                  output directory

General:
  --threads INT                (default ${THREADS})
  -v|--verbose                 verbose logs
  -q|--quiet                   quiet logs

syncasm ladder:
  --k-list "401,251,151,121"   k ladder (first PASS wins)
  --smer INT                   -s (default ${SMER})
  --no-read-ec                 pass --no-read-ec
  --a FLOAT                    -a min arc coverage (optional)
  --max-bubble INT             --max-bubble (default ${MAX_BUB})
  --max-tip INT                --max-tip (optional)
  --weak-cross FLOAT           --weak-cross (default ${WEAK_X})
  --unzip-round INT            --unzip-round (default ${UNZIP_ROUND})
  --c INT                      -c min k-mer coverage (optional)
  --mt-size-est INT            estimated organelle size (tightens judge)

Judge thresholds:
  --j-max-unitigs INT          (default ${J_MAX_UNITIGS})
  --j-min-total-bp INT         (default ${J_MIN_TOTAL_BP})
  --j-tangle-max FLOAT         (default ${J_TANGLE_MAX})
  --j-min-n50 INT              (default ${J_MIN_N50})
  --j-min-maxlen INT           (default ${J_MIN_MAXLEN})
  --j-circ-min-ovlp INT        (default ${J_CIRC_MIN_OVLP})
  --j-circ-min-ident FLOAT     (default ${J_CIRC_MIN_IDENT})
EOF
}

# ---------------- parse args ----------------
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
	--threads)
		THREADS="$2"
		shift 2
		;;
	--k-list)
		K_LIST="$2"
		shift 2
		;;
	--smer)
		SMER="$2"
		shift 2
		;;
	--no-read-ec)
		NO_READ_EC=1
		shift
		;;
	--a)
		AARC="$2"
		shift 2
		;;
	--max-bubble)
		MAX_BUB="$2"
		shift 2
		;;
	--max-tip)
		MAX_TIP="$2"
		shift 2
		;;
	--weak-cross)
		WEAK_X="$2"
		shift 2
		;;
	--unzip-round)
		UNZIP_ROUND="$2"
		shift 2
		;;
	--c)
		COV="$2"
		shift 2
		;;
	--mt-size-est)
		MT_SIZE_EST="$2"
		shift 2
		;;
	--j-max-unitigs)
		J_MAX_UNITIGS="$2"
		shift 2
		;;
	--j-min-total-bp)
		J_MIN_TOTAL_BP="$2"
		shift 2
		;;
	--j-tangle-max)
		J_TANGLE_MAX="$2"
		shift 2
		;;
	--j-min-n50)
		J_MIN_N50="$2"
		shift 2
		;;
	--j-min-maxlen)
		J_MIN_MAXLEN="$2"
		shift 2
		;;
	--j-circ-min-ovlp)
		J_CIRC_MIN_OVLP="$2"
		shift 2
		;;
	--j-circ-min-ident)
		J_CIRC_MIN_IDENT="$2"
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
	-h | --help)
		usage
		exit 0
		;;
	*) die "unknown arg: $1" ;;
	esac
done
[[ -z "$READS" || -z "$OUT" ]] && {
	usage
	exit 1
}

READS="$(readlink -f "$READS")"
OUT="$(readlink -f "$OUT")"
mkdir -p "$OUT"/{qc,k1}
cd "$OUT"

# ---------------- deps ----------------
need syncasm
need python3
[[ -s "$SCRIPT_DIR/polap-py-oatk-core-judge-assembly.py" ]] || die "judge script missing: $SCRIPT_DIR/polap-py-oatk-core-judge-assembly.py"

# ---------------- helpers ----------------

# parse singleton % from syncasm stderr
# returns "smer_pct<TAB>kmer_pct" or ".\t." if missing
parse_sr_singletons() {
	local err="$1"
	local smer kmer
	smer=$(awk '/sr_db_stat] number uniqe smer:/ { if (match($0, /\(([0-9.]+)%\)/, m)) s=m[1] } END{ if (s=="") print "."; else printf "%.3f", s }' "$err" 2>/dev/null)
	kmer=$(awk '/sr_db_stat] number uniqe kmer:/ { if (match($0, /\(([0-9.]+)%\)/, m)) k=m[1] } END{ if (k=="") print "."; else printf "%.3f", k }' "$err" 2>/dev/null)
	printf "%s\t%s\n" "$smer" "$kmer"
}

extract_unitigs_from_gfa() { # $1=gfa $2=out.fa
	awk '$1=="S"{print ">"$2"\n"$3}' "$1" >"$2"
	[[ -s "$2" ]] || return 1
}

# judge (new GFA-only interface). arg: outd=k_try_kXXX
judge_stop() {
	local outd="$1"
	local gfa_final="$outd/syncasm.asm.utg.final.gfa"
	local gfa_utg="$outd/syncasm.asm.utg.gfa"
	local diag="qc/judge.$(basename "$outd").summary.txt"
	local tsv="qc/judge.all.tsv"
	local json="qc/judge.$(basename "$outd").json"

	local jargs=()
	[[ -s "$gfa_final" ]] && jargs+=(--gfa-final "$gfa_final")
	[[ -s "$gfa_utg" ]] && jargs+=(--gfa-utg "$gfa_utg")

	jargs+=(--threads "$THREADS"
		--max-unitigs "$J_MAX_UNITIGS"
		--min-total-bp "$J_MIN_TOTAL_BP"
		--min-maxlen "$J_MIN_MAXLEN"
		--min-n50 "$J_MIN_N50"
		--tangle-max "$J_TANGLE_MAX"
		--circ-min-ovlp "$J_CIRC_MIN_OVLP"
		--circ-min-ident "$J_CIRC_MIN_IDENT")
	((MT_SIZE_EST > 0)) && jargs+=(--mt-size-est "$MT_SIZE_EST")

	# optional coverage evenness (comment out to skip)
	# [[ -n "$READS" ]] && jargs+=( --reads "$READS" --sample-n 2000 --cov-bin 100 --cov-cv-max 0.60 --cov-mean-min 3.0 )

	jargs+=(--diag-file "$diag" --tsv-file "$tsv")

	local verdict
	verdict="$(python3 "$SCRIPT_DIR/polap-py-oatk-core-judge-assembly.py" "${jargs[@]}" 2>"$json")"

	# log human SUMMARY if available
	if [[ -s "$diag" ]]; then
		local summary
		summary="$(tail -n1 "$diag" || true)"
		[[ -n "$summary" ]] && log 1 "[judge] ${summary#SUMMARY }"
	fi
	echo "$verdict"
}

# run one syncasm attempt. returns a TSV line:
# ufa<TAB>gfa<TAB>smer_singletons<TAB>kmer_singletons  (or return 2 on failure)
run_syncasm_try() {
	local k="$1" outd="k_try_k${k}"
	mkdir -p "$outd"

	local cmd=(syncasm -k "$k" -s "$SMER" -t "$THREADS" -o "$outd/syncasm.asm")
	((NO_READ_EC == 1)) && cmd+=(--no-read-ec)
	[[ -n "$AARC" ]] && cmd+=(-a "$AARC")
	[[ -n "$MAX_BUB" ]] && cmd+=(--max-bubble "$MAX_BUB")
	[[ -n "$MAX_TIP" ]] && cmd+=(--max-tip "$MAX_TIP")
	[[ -n "$WEAK_X" ]] && cmd+=(--weak-cross "$WEAK_X")
	[[ -n "$UNZIP_ROUND" ]] && cmd+=(--unzip-round "$UNZIP_ROUND")
	[[ -n "$COV" ]] && cmd+=(-c "$COV")

	log 1 "[assemble] k=$k → ${cmd[*]} ${READS##*/}"

	# capture stdout & stderr
	$TIME_BIN -f "TIME(s) %E  RSS(KB) %M" \
		"${cmd[@]}" "$READS" \
		>"$outd/run.stdout" \
		2>"$outd/run.stderr" || true

	# small stderr summary + singleton %
	grep -E 'sr_db_stat] number uniqe (smer|kmer)|singletons|peak_hom|peak_het' \
		"$outd/run.stderr" >"qc/summary.k${k}.txt" || true

	local SMER_SING="." KMER_SING="."
	if [[ -s "$outd/run.stderr" ]]; then
		local parsed
		parsed="$(parse_sr_singletons "$outd/run.stderr" 2>/dev/null || echo ".\t.")"
		IFS=$'\t' read -r SMER_SING KMER_SING <<<"$parsed"
		[[ -z "${SMER_SING:-}" ]] && SMER_SING="."
		[[ -z "${KMER_SING:-}" ]] && KMER_SING="."
	fi
	log 1 "[syncasm] k=${k} smer_singletons=${SMER_SING}% kmer_singletons=${KMER_SING}%"

	# discover outputs; prefer final gfa for unitigs
	local gfaF="$outd/syncasm.asm.utg.final.gfa"
	local gfaI="$outd/syncasm.asm.utg.gfa"
	local gfa="$outd/syncasm.asm.gfa"
	local ufa="" gfapath=""

	if [[ -s "$gfaF" ]]; then
		gfapath="$gfaF"
	elif [[ -s "$gfaI" ]]; then
		gfapath="$gfaI"
	elif [[ -s "$gfa" ]]; then
		gfapath="$gfa"
	fi

	if [[ -n "$gfapath" ]]; then
		extract_unitigs_from_gfa "$gfapath" "$outd/unitigs.fa" || true
		[[ -s "$outd/unitigs.fa" ]] && ufa="$outd/unitigs.fa"
	fi

	if [[ -n "$ufa" && -s "$ufa" && -n "$gfapath" && -s "$gfapath" ]]; then
		printf "%s\t%s\t%s\t%s\n" "$ufa" "$gfapath" "$SMER_SING" "$KMER_SING"
	else
		return 2
	fi
}

# ---------------- init QC ----------------
SUMMARY="qc/assemble.out"
: >"$SUMMARY"
log 1 "[stage] polap-bash-oatk-core-ont.sh v0.0.5"
log 1 "[ladder] k-list=${K_LIST}; smer=${SMER}; threads=${THREADS}; EC=$((NO_READ_EC ? 0 : 1))"

# ---------------- run ladder ----------------
K_USED=""
UFA_FINAL=""
best_key=""
best_row=""

IFS=',' read -r -a KS <<<"$K_LIST"
for k in "${KS[@]}"; do
	k="${k// /}"
	[[ -z "$k" ]] && continue

	if out_pair="$(run_syncasm_try "$k")"; then
		IFS=$'\t' read -r ufa gfa SMER_SING KMER_SING <<<"$out_pair"

		# quick FA stats
		n_unitigs="$(awk '/^>/{c++} END{print c+0}' "$ufa")"
		TBP="$(awk '/^>/{if(len){t+=len}len=0;next}{len+=length}END{t+=len;print t+0}' "$ufa")"

		# judge
		log 1 "before: judge"
		verdict="$(judge_stop "k_try_k${k}")"
		log 1 "after: judge"

		# log & append row (include singleton %)
		log 1 "k=${k} unitigs=${n_unitigs} total_bp=${TBP} judge=${verdict} smer_singletons=${SMER_SING}% kmer_singletons=${KMER_SING}%"
		printf "k=%s\tunitigs=%s\ttotal_bp=%s\tjudge=%s\tsmer_singletons=%s%%\tkmer_singletons=%s%%\n" \
			"$k" "$n_unitigs" "$TBP" "$verdict" "$SMER_SING" "$KMER_SING" >>"$SUMMARY"

		# keep best by (total_bp, then max_len via quick awk)
		mx_len="$(awk 'BEGIN{mx=0} /^>/{next} {L+=length} END{if(L>mx)mx=L; print mx}' "$ufa" 2>/dev/null || echo 0)"
		key=$(printf "%012d:%012d" "$TBP" "$mx_len")
		[[ -z "$best_key" || "$key" > "$best_key" ]] && {
			best_key="$key"
			best_row="k=${k} unitigs=${n_unitigs} total_bp=${TBP} max=${mx_len}"
		}

		if [[ "$verdict" == "PASS" ]]; then
			K_USED="$k"
			UFA_FINAL="$ufa"
			ln -sf "../${UFA_FINAL}" k1/unitigs.fa 2>/dev/null || cp -f "${UFA_FINAL}" k1/unitigs.fa
			echo "$K_USED" >k1/USED_K.txt
			log 1 "[done] PASS at k=$K_USED → k1/unitigs.fa"
			log 1 "[emit] QC table: $SUMMARY"
			exit 0
		fi
	else
		log 1 "k=${k} no_outputs"
		printf "k=%s\tno_outputs\n" "$k" >>"$SUMMARY"
	fi
done

# ---------------- NO PASS report ----------------
log 1 "[ladder] no PASS, but assemblies exist. See $SUMMARY"
if [[ -n "$best_row" ]]; then
	log 1 "[best] ${best_row}"
fi
exit 2
