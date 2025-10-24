#!/usr/bin/env bash
# polap-bash-pt-two-isoforms.sh  v0.0.1
# Make two plastid isoforms (SSC forward / SSC inverted) from a circular plastome.
# Requires: nucmer, show-coords (MUMmer4), seqkit, awk, grep, paste, mktemp

set -euo pipefail

usage() {
	cat <<'EOF'
Usage:
  polap-bash-pt-two-isoforms.sh -i plastome.fa -o outprefix [options]

Options:
  -i FILE   Input plastome FASTA (single circular sequence; header first record used)
  -o PREFIX Output prefix (default: ptiso)
  -t INT    Threads for nucmer (default: 8)
  --min-ir INT   Minimum IR length bp to consider (default: 8000)
  --pad INT      Duplicate this many bp at the end for boundary-safe mapping (default: 5000)
  -v        Verbose
  -h        Help

Outputs (all FASTA):
  PREFIX.formA.fa                  LSC–IRa–SSC–IRb
  PREFIX.formB.fa                  LSC–IRa–revcomp(SSC)–IRb
  PREFIX.formA.pad.fa              same as A, but ends padded by --pad bp
  PREFIX.formB.pad.fa              same as B, but ends padded by --pad bp
EOF
}

# ---- defaults
IN=""
OUT="ptiso"
THREADS=8
MIN_IR=8000
PAD=5000
VERB=0

# ---- parse
while [[ $# -gt 0 ]]; do
	case "$1" in
	-i)
		IN="$2"
		shift 2
		;;
	-o)
		OUT="$2"
		shift 2
		;;
	-t)
		THREADS="$2"
		shift 2
		;;
	--min-ir)
		MIN_IR="$2"
		shift 2
		;;
	--pad)
		PAD="$2"
		shift 2
		;;
	-v)
		VERB=1
		shift
		;;
	-h | --help)
		usage
		exit 0
		;;
	*)
		echo "[ERR] unknown arg: $1"
		usage
		exit 2
		;;
	esac
done

[[ -s "$IN" ]] || {
	echo "[ERR] -i plastome FASTA required"
	exit 2
}

# ---- check tools
for x in nucmer show-coords seqkit awk grep sed paste; do
	command -v "$x" >/dev/null 2>&1 || {
		echo "[ERR] missing tool: $x"
		exit 2
	}
done

# ---- get the first (and ideally only) sequence ID and length
SEQID=$(seqkit seq -ni "$IN" | head -n1)
[[ -n "$SEQID" ]] || {
	echo "[ERR] no sequence ID found in $IN"
	exit 2
}
LEN=$(seqkit fx2tab -l "$IN" | awk '{print $2}' | head -n1)
[[ "$LEN" =~ ^[0-9]+$ ]] || {
	echo "[ERR] failed to get length"
	exit 2
}

((VERB)) && echo "[info] seq: $SEQID len=$LEN"

# ---- self-alignment with nucmer
TMPD=$(mktemp -d)
trap 'rm -rf "$TMPD"' EXIT

NUCPFX="$TMPD/self"
nucmer --maxmatch -t "$THREADS" -p "$NUCPFX" "$IN" "$IN" >/dev/null 2>&1
show-coords -rclT "$NUCPFX".delta >"$NUCPFX".coords

# coords format: start1 end1 start2 end2 len1 len2 %id ...
# We want reverse-complement matches (inverted) that are long; heuristic:
#  - orientation is opposite (start2>end2 OR start1>end1) -> inverted
#  - take top two longest non-overlapping hits as IRa/IRb
awk -v MIN="$MIN_IR" '
  NR>5 {
    s1=$1; e1=$2; s2=$3; e2=$4; l1=$5; l2=$6; pid=$7
    inv = ( (s1<e1 && s2>e2) || (s1>e1 && s2<e2) )
    len = (l1<l2?l1:l2)
    if (inv && len>=MIN) {
      if (s1>e1){t=s1;s1=e1;e1=t}  # normalize coords so s1<e1 etc.
      if (s2>e2){t=s2;s2=e2;e2=t}
      print s1, e1, s2, e2, len, pid
    }
  }
' "$NUCPFX".coords | sort -k5,5nr >"$NUCPFX".inv.tsv

[[ -s "$NUCPFX".inv.tsv ]] || {
	echo "[ERR] no inverted repeats >= ${MIN_IR}bp found"
	exit 2
}

# Pick two longest non-overlapping on ref coordinates as IRs
# Simple greedy: take first row as IRa; then pick next row whose (s1..e1) doesn’t overlap IRa.
read -r A1 A2 B1 B2 ALEN APID < <(head -n1 "$NUCPFX".inv.tsv)
IR_A1=$A1
IR_A2=$A2
IR_B1=""
IR_B2=""
while read -r s1 e1 s2 e2 l pid; do
	# non-overlap on reference axis
	if ( (e1 <IR_A1) || (s1 >IR_A2)); then
		IR_B1=$s1
		IR_B2=$e1
		BLEN=$l
		BPID=$pid
		break
	fi
done < <(tail -n +2 "$NUCPFX".inv.tsv)

[[ -n "$IR_B1" ]] || {
	echo "[ERR] could not find a second non-overlapping IR match"
	exit 2
}

# sort IRs so IR_A1 < IR_A2 < IR_B1 < IR_B2 along linearized coordinates
if ((IR_B1 < IR_A1)); then
	# swap to keep A before B
	t1=$IR_A1
	t2=$IR_A2
	IR_A1=$IR_B1
	IR_A2=$IR_B2
	IR_B1=$t1
	IR_B2=$t2
fi

((VERB)) && echo "[info] IRa=${IR_A1}-${IR_A2}  IRb=${IR_B1}-${IR_B2} (lenA≈$ALEN lenB≈$BLEN)"

# For a typical plastome: IRa < IRb; SSC = (IR_A2+1 .. IR_B1-1); LSC = wrap (IR_B2+1 .. LEN) + (1 .. IR_A1-1)
# Guard indices to 1..LEN
wrap() { # clamp to [1..LEN]
	local v=$1
	if ((v < 1)); then v=$((v + LEN)); fi
	if ((v > LEN)); then v=$((v - LEN)); fi
	echo "$v"
}

SSC_S=$((IR_A2 + 1))
SSC_E=$((IR_B1 - 1))
if ((SSC_S > SSC_E)); then
	# circle wrap (unlikely if IR order is correct)
	echo "[ERR] inferred SSC wraps; IR detection odd. Aborting."
	exit 2
fi

LSC1_S=$((IR_B2 + 1))
LSC1_E=$LEN
LSC2_S=1
LSC2_E=$((IR_A1 - 1))

# extract slices with seqkit subseq --range
# helper to extract [S..E], allowing empty if S>E (skip)
grab_range() {
	local f=$1 s=$2 e=$3
	if ((s <= e)); then
		seqkit subseq --range ${s}:${e} "$f"
	fi
}

# Build Form A: LSC (wrap) + IRa + SSC + IRb
# IR intervals themselves we take exactly as detected (IR_A1..IR_A2), (IR_B1..IR_B2)
# 1) LSC = [IR_B2+1..LEN] + [1..IR_A1-1]
grab_range "$IN" "$LSC1_S" "$LSC1_E" >"$TMPD/_LSC1.fa"
grab_range "$IN" "$LSC2_S" "$LSC2_E" >"$TMPD/_LSC2.fa"
# 2) IRa
seqkit subseq --range ${IR_A1}:${IR_A2} "$IN" >"$TMPD/_IRa.fa"
# 3) SSC
seqkit subseq --range ${SSC_S}:${SSC_E} "$IN" >"$TMPD/_SSC.fa"
# 4) IRb
seqkit subseq --range ${IR_B1}:${IR_B2} "$IN" >"$TMPD/_IRb.fa"

# Concatenate pieces into one record; keep original header but annotate
# Use seqkit to standardize single-record sequences then paste sequences (drop headers) carefully.
seqkit seq -w 0 "$TMPD/_LSC1.fa" | awk 'NR==2{print}' >"$TMPD/LSC.txt"
seqkit seq -w 0 "$TMPD/_LSC2.fa" | awk 'NR==2{print}' >>"$TMPD/LSC.txt"
seqkit seq -w 0 "$TMPD/_IRa.fa" | awk 'NR==2{print}' >"$TMPD/IRa.txt"
seqkit seq -w 0 "$TMPD/_SSC.fa" | awk 'NR==2{print}' >"$TMPD/SSC.txt"
seqkit seq -w 0 "$TMPD/_IRb.fa" | awk 'NR==2{print}' >"$TMPD/IRb.txt"

cat >"${OUT}.formA.fa" <<EOF
>${SEQID}|formA LSC-IRa-SSC-IRb IRa=${IR_A1}-${IR_A2} IRb=${IR_B1}-${IR_B2} SSC=${SSC_S}-${SSC_E}
$(cat "$TMPD/LSC.txt" "$TMPD/IRa.txt" "$TMPD/SSC.txt" "$TMPD/IRb.txt")
EOF

# Build Form B: LSC–IRa–revcomp(SSC)–IRb
seqkit seq -r -p "$TMPD/_SSC.fa" >"$TMPD/_SSCrc.fa"
seqkit seq -w 0 "$TMPD/_SSCrc.fa" | awk 'NR==2{print}' >"$TMPD/SSCrc.txt"
cat >"${OUT}.formB.fa" <<EOF
>${SEQID}|formB LSC-IRa-SSCrc-IRb IRa=${IR_A1}-${IR_A2} IRb=${IR_B1}-${IR_B2} SSC=${SSC_S}-${SSC_E}(rc)
$(cat "$TMPD/LSC.txt" "$TMPD/IRa.txt" "$TMPD/SSCrc.txt" "$TMPD/IRb.txt")
EOF

# Make boundary-padded versions (duplicate first PAD bp at end)
pad_one() {
	local infa=$1 outfa=$2 pad=$3
	local headseq
	headseq=$(seqkit subseq --range 1:${pad} "$infa" | seqkit seq -w 0 | awk 'NR==2{print}')
	{
		read -r hdr < <(grep "^>" "$infa")
		seq=$(seqkit seq -w 0 "$infa" | awk 'NR==2{print}')
		echo "${hdr}|pad=${pad}"
		echo "${seq}${headseq}" | fold -w 60
	} >"$outfa"
}
pad_one "${OUT}.formA.fa" "${OUT}.formA.pad.fa" "$PAD"
pad_one "${OUT}.formB.fa" "${OUT}.formB.pad.fa" "$PAD"

((VERB)) && {
	echo "[done] Wrote:"
	printf "  %s\n" "${OUT}.formA.fa" "${OUT}.formB.fa" "${OUT}.formA.pad.fa" "${OUT}.formB.pad.fa"
}
