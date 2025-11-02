#!/usr/bin/env bash
# polap-bash-standardize-plastome.sh
# Version: v0.2.0
# Rotate a plastome so that the chosen anchor's 5' base sits at position 1 on '+'
# Dependencies: minimap2, seqkit (for reverse-complement + circular restart), python3
# References: minimap2 [Li 2018], seqkit docs/manpage

set -euo pipefail
IFS=$'\n\t'

usage() {
	cat <<'USAGE'
Usage:
  polap-bash-standardize-plastome.sh \
    -i plastome.fa \
    -a anchors.fa \
    -o outdir \
    [--seqid <id>] \
    [--minimap2-preset asm5] \
    [--prefer-plus] \
    [--report-all-hits]

Notes:
  * anchors.fa may contain multiple candidates (psbA, trnH-GUG, rbcL, matK, ycf1, ndhF).
  * If plastome.fa has multiple sequences, --seqid selects which to standardize.
  * Output: outdir/standardized.fasta (anchor at pos1, '+'), outdir/summary.tsv
USAGE
}

INP_FA=""
ANCHORS=""
OUT=""
SEQID=""
MM_PRESET="asm5"
PREFER_PLUS=0
REPORT_ALL=0

while [[ $# -gt 0 ]]; do
	case "$1" in
	-i)
		INP_FA="$2"
		shift 2
		;;
	-a)
		ANCHORS="$2"
		shift 2
		;;
	-o)
		OUT="$2"
		shift 2
		;;
	--seqid)
		SEQID="$2"
		shift 2
		;;
	--minimap2-preset)
		MM_PRESET="$2"
		shift 2
		;;
	--prefer-plus)
		PREFER_PLUS=1
		shift 1
		;;
	--report-all-hits)
		REPORT_ALL=1
		shift 1
		;;
	-h | --help)
		usage
		exit 0
		;;
	*)
		echo "[ERR] Unknown arg: $1"
		usage
		exit 1
		;;
	esac
done

[[ -z "$INP_FA" || -z "$ANCHORS" || -z "$OUT" ]] && {
	usage
	exit 1
}
mkdir -p "$OUT"

# 1) Map anchors -> plastome (PAF)
paf="$OUT/hits.paf"
minimap2 -x "$MM_PRESET" -c --secondary=no "$INP_FA" "$ANCHORS" >"$paf" 2>"$OUT/minimap2.log"

# 2) Pick the best anchor hit
PICKER="$(dirname "$0")/scripts/pick_anchor_from_paf.py"
python3 "$PICKER" \
	--paf "$paf" \
	--prefer-plus $PREFER_PLUS \
	--report-all $REPORT_ALL \
	--out-tsv "$OUT/hits.parsed.tsv" >"$OUT/anchor.choice.tsv"

read -r QNAME TNAME TSTART TEND STRAND SCORE < <(tail -n1 "$OUT/anchor.choice.tsv" || true)
[[ -z "${TSTART:-}" ]] && {
	echo "[ERR] No valid anchor hit found"
	exit 2
}

# 3) Determine working seqid and sequence length
if [[ -n "$SEQID" ]]; then
	WORK_ID="$SEQID"
else
	WORK_ID="$TNAME"
fi

LEN=$(seqkit fx2tab -l "$INP_FA" | awk -v id="$WORK_ID" '$1==id{print $2; found=1} END{if(!found)exit 1}')
# If the header was not a simple single token, fallback: take first record
if [[ -z "${LEN:-}" ]]; then
	WORK_ID=$(seqkit seq -n "$INP_FA" | head -n1)
	LEN=$(seqkit fx2tab -l "$INP_FA" | head -n1 | awk '{print $2}')
fi

# 4) If the anchor aligned on '-', reverse-complement the whole plastome
tmp_rc="$OUT/tmp.rc.fasta"
tmp_in="$OUT/tmp.in.fasta"
tmp_work="$OUT/tmp.work.fasta"
mkdir -p "$OUT/tmp"

# Extract just the target sequence if multi-FASTA
seqkit grep -n -p "$WORK_ID" "$INP_FA" >"$tmp_in"

ORI="+"
Astart="$TSTART"
Aend="$TEND"

if [[ "$STRAND" == "-" ]]; then
	# reverse-complement
	seqkit seq -r -p "$tmp_in" >"$tmp_rc"
	mv "$tmp_rc" "$tmp_work"
	ORI="-"
	# new anchor 5' position after RC becomes: L - Aend + 1
	POS=$((LEN - Aend + 1))
else
	cp "$tmp_in" "$tmp_work"
	POS="$Astart"
fi

# 5) Rotate circular start so that POS -> 1
# seqkit restart uses 1-based positions for circular genomes
OUT_FA="$OUT/standardized.fasta"
seqkit restart -i "$POS" "$tmp_work" -o "$OUT_FA"

# 6) Write summary
SUM="$OUT/summary.tsv"
{
	echo -e "input\tseqid\tlength\tanchor\tpos_before\tstrand_before\torientation\trotation_start\tpos_after"
	echo -e "$(basename "$INP_FA")\t$WORK_ID\t$LEN\t$QNAME\t$TSTART\t$STRAND\t$ORI\t$POS\t1"
} >"$SUM"

echo "[OK] Standardized FASTA: $OUT_FA"
echo "[OK] Summary: $SUM"
