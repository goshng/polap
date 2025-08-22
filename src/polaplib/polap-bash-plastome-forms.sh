#!/usr/bin/env bash
set -euo pipefail

# plastome_forms.sh
# Build four plastome forms with canonical printing:
#   1: LSC-IR-SSC-IR′
#   2: LSC-IR-SSC′-IR′
#   3: LSC′-IR-SSC′-IR′
#   4: LSC′-IR-SSC-IR′
#
# Key features:
# - IR detection via MUMmer3 `palindrome` (exact inverted repeats), fallback to MUMmer4 `nucmer` self-vs-self (reverse hits only).
# - IRa chosen so the neighboring arcs satisfy LSC–IRa–SSC (LSC longer), with optional rRNA preference (IRa adjacent to SSC contains rRNA).
# - IRb is taken from genome coordinates (not rc(IRa)); then both IR copies are trimmed to the same core length to preserve N.
# - Final printed forms always end with IR′ = rc(IRa), matching common reference orientation.
# - Outputs: DIR/1.fa..4.fa, boundary tables, and an optional thorough check report.
#
# Usage:
#   plastome_forms.sh [options] ptDNA.fa
#
# Options:
#   -o, --out DIR            Output directory (default: out)
#   -c, --check              Ensure one output equals the input exactly (by rotation), and write boundary files
#       --check-report FILE  Write a consolidated check report (implies --check)
#       --min-ir INT         Minimum IR length (default: $MIN_IR_MIN or 8000)
#       --min-id FLOAT       Minimum identity (0–1) for nucmer fallback (default: $MIN_ID or 0.90)
#       --rrna-ref FILE      FASTA of plastid 16S/23S (optional, prefers IRa with rRNA adjacent to SSC)
#       --rrna-min-len INT   Min aligned length for rRNA hits (default: $RRNA_MINLEN or 600)
#       --rrna-min-id FLOAT  Min identity (0–1) for rRNA hits with blastn (default: $RRNA_MINID or 0.85)

# ---- ensure bash ----
if [ -z "${BASH_VERSION:-}" ]; then
	echo "Please run with bash, not sh." >&2
	exit 2
fi

# ---- defaults / env fallbacks ----
OUTDIR="out"
CHECK=0
CHECK_REPORT=""
MIN_IR=${MIN_IR_MIN:-8000}
MIN_ID=${MIN_ID:-0.90}
RRNA_REF=""
RRNA_MINLEN=${RRNA_MINLEN:-600}
RRNA_MINID=${RRNA_MINID:-0.85}

print_help() {
	sed -n '1,120p' "$0" | sed 's/^# \{0,1\}//'
	exit 0
}

# ---- parse args ----
args=()
while [[ $# -gt 0 ]]; do
	case "$1" in
	-h | --help) print_help ;;
	-o | --out)
		OUTDIR="$2"
		shift 2
		;;
	-c | --check)
		CHECK=1
		shift
		;;
	--check-report)
		CHECK_REPORT="$2"
		CHECK=1
		shift 2
		;;
	--min-ir)
		MIN_IR="$2"
		shift 2
		;;
	--min-id)
		MIN_ID="$2"
		shift 2
		;;
	--rrna-ref)
		RRNA_REF="$2"
		shift 2
		;;
	--rrna-min-len)
		RRNA_MINLEN="$2"
		shift 2
		;;
	--rrna-min-id)
		RRNA_MINID="$2"
		shift 2
		;;
	--)
		shift
		break
		;;
	-*)
		echo "Unknown option: $1" >&2
		exit 1
		;;
	*)
		args+=("$1")
		shift
		;;
	esac
done
set -- "${args[@]}"
FA=${1:?Input plastome FASTA required}
mkdir -p "$OUTDIR"

# ---- tool checks ----
need() { command -v "$1" >/dev/null 2>&1; }
if ! need seqkit; then
	echo "Missing: seqkit" >&2
	exit 1
fi
PAL_OK=0
if need palindrome; then PAL_OK=1; fi
NUC_OK=0
if need nucmer && need show-coords && need delta-filter; then NUC_OK=1; fi
if ((PAL_OK == 0 && NUC_OK == 0)); then
	echo "Need MUMmer3 (palindrome) or MUMmer4 (nucmer, show-coords, delta-filter)." >&2
	exit 1
fi

work=$(mktemp -d -p . plastome_forms.XXXXXX)
trap 'rm -rf "$work"' EXIT

# ---- helpers ----
minmax() { [[ $1 -le $2 ]] && echo "$1 $2" || echo "$2 $1"; }

rotate_str() { # $1=S, $2=pos1based
	local S k n
	S=$1
	k=$2
	n=${#S}
	((n == 0)) && {
		echo ""
		return
	}
	((k = (k - 1) % n, k < 0 ? k += n : 0))
	echo "${S:k}${S:0:k}"
}

slice_circ() { # $1=S, $2=a, $3=b (1-based, half-open [a,b))
	local S a b n len start S2
	S=$1
	a=$2
	b=$3
	n=${#S}
	if ((b > a)); then len=$((b - a)); else len=$((n - (a - b))); fi
	start=$((a - 1))
	S2="${S}${S}"
	echo "${S2:start:len}"
}

advance_pos() { # $1=pos1based $2=delta $3=n -> new pos1based
	awk -v p="$1" -v d="$2" -v n="$3" 'BEGIN{print ((p-1+d)%n)+1}'
}

slice_circ_len() { # $1=S $2=start1based $3=len
	local S a len n start S2
	S=$1
	a=$2
	len=$3
	n=${#S}
	start=$((a - 1))
	S2="${S}${S}"
	echo "${S2:start:len}"
}

rc_of() { printf ">x\n%s\n" "$1" | seqkit seq -t DNA -w 0 -rp | tail -n +2; }

circular_find_offset() { # $1=A $2=B  -> offset(1-based) or 0
	local A=$1 B=$2 n=${#B}
	[[ ${#A} -ne $n || $n -eq 0 ]] && {
		echo 0
		return
	}
	printf "%s" "$A" >"$work/needle.txt"
	printf "%s%s" "$B" "$B" >"$work/hay.txt"
	awk '
    BEGIN{
      while ((getline l < ARGV[1]) > 0) A = A l; close(ARGV[1])
      while ((getline l < ARGV[2]) > 0) D = D l; close(ARGV[2])
      n = ARGV[3]+0
      pos = index(D, A)
      if (pos>=1 && pos<=n) print pos; else print 0
    }
  ' "$work/needle.txt" "$work/hay.txt" "$n"
}

arc_len() { # [a,b) on circle length n
	local a=$1 b=$2 n=$3
	if ((b > a)); then echo $((b - a)); else echo $((n - (a - b))); fi
}

# ---- load input ----
seqkit seq -w 0 "$FA" >"$work/in.fa"
[[ $(grep -c '^>' "$work/in.fa") -eq 1 ]] || {
	echo "Input must be exactly one FASTA record." >&2
	exit 1
}
name=$(head -1 "$work/in.fa" | sed 's/^>//; s/ .*//')
INSEQ=$(tail -n +2 "$work/in.fa")
N=${#INSEQ}
[[ $N -gt 0 ]] || {
	echo "Empty sequence." >&2
	exit 1
}

# ---- IR detection ----
IRa_start=0
IRa_end=0
IRb_start=0
IRb_end=0

if ((PAL_OK)); then
	printf ">%s\n%s\n" "$name" "$INSEQ" >"$work/p.fa"
	# exact inverted repeats; we pick the longest >= MIN_IR
	palindrome -l 20 -b -c "$work/p.fa" >"$work/pal.out" || true
	best=$(awk -v MINL="$MIN_IR" '
    $1 ~ /^[0-9]+$/ && NF>=6 { L=($5<$6?$5:$6); if (L>=MINL && L>best){best=L; line=$0} }
    END{if(line!="")print line}
  ' "$work/pal.out" || true)
	if [[ -n "$best" ]]; then
		read s1 e1 s2 e2 l1 l2 <<<"$best"
		read IRa_start IRa_end <<<"$(minmax "$s1" "$e1")"
		read IRb_start IRb_end <<<"$(minmax "$s2" "$e2")"
	fi
fi

if ((IRa_start == 0)); then
	# fallback to nucmer self-vs-self, reverse hits only
	printf ">%s\n%s\n" "$name" "$INSEQ" >"$work/self.fa"
	nucmer --maxmatch -l 100 -c 500 -p "$work/self" "$work/self.fa" "$work/self.fa" >/dev/null
	# parse raw to preserve strand; reverse hits have S2 > E2
	show-coords -HT "$work/self.delta" >"$work/coords.txt"
	best=$(awk -v MINL="$MIN_IR" -v MINIDPCT="$(awk -v x="$MIN_ID" 'BEGIN{print 100*x}')" '
    $1 ~ /^[0-9]+$/ && NF>=7 {
      s1=$1; e1=$2; s2=$3; e2=$4; L1=$5; L2=$6; ID=$7+0;
      L=(L1<L2?L1:L2);
      # reverse
      if (L>=MINL && ID>=MINIDPCT && s2>e2) {
        score=L*(ID/100.0);
        if (score>best){best=score; line=$0}
      }
    }
    END{if(line!="")print line}
  ' "$work/coords.txt" || true)
	[[ -n "$best" ]] || {
		echo "No IR candidate (min-ir=$MIN_IR, min-id=$MIN_ID)." >&2
		exit 1
	}
	read s1 e1 s2 e2 l1 l2 idy <<<"$best"
	read IRa_start IRa_end <<<"$(minmax "$s1" "$e1")"
	read IRb_start IRb_end <<<"$(minmax "$e2" "$s2")" # reverse arm normalized
fi

# ---- Decide IRa as the copy that yields LSC–IRa–SSC; prefer rRNA if provided ----
# Candidate A: treat (IRa_start,IRa_end) as IRa
lenIR_A=$((IRa_end - IRa_start + 1))
qS1_A=$((((IRb_start - IRa_start + N) % N) + 1))
qE1_A=$((((IRb_end - IRa_start + N) % N) + 1))
after_A=$(arc_len $((lenIR_A + 1)) "$qS1_A" "$N")
before_A=$(arc_len $((qE1_A + 1)) 1 "$N")

# Candidate B: swap roles
lenIR_B=$((IRb_end - IRb_start + 1))
qS1_B=$((((IRa_start - IRb_start + N) % N) + 1))
qE1_B=$((((IRa_end - IRb_start + N) % N) + 1))
after_B=$(arc_len $((lenIR_B + 1)) "$qS1_B" "$N")
before_B=$(arc_len $((qE1_B + 1)) 1 "$N")

geom_pref_A=$((before_A > after_A ? 1 : 0))
geom_pref_B=$((before_B > after_B ? 1 : 0))

rrna_hits_A=0
rrna_hits_B=0
if [[ -n "$RRNA_REF" && -s "$RRNA_REF" ]]; then
	# extract the two IRs from original genome
	GA=$(rotate_str "$INSEQ" "$IRa_start")
	IRseqA="${GA:0:lenIR_A}"
	GB=$(rotate_str "$INSEQ" "$IRb_start")
	IRseqB="${GB:0:lenIR_B}"
	printf ">IR_A\n%s\n" "$IRseqA" >"$work/IR_A.fa"
	printf ">IR_B\n%s\n" "$IRseqB" >"$work/IR_B.fa"
	if command -v blastn >/dev/null 2>&1; then
		blastn -subject "$work/IR_A.fa" -query "$RRNA_REF" -outfmt '6 pident length' >"$work/rrnaA.tsv" || true
		blastn -subject "$work/IR_B.fa" -query "$RRNA_REF" -outfmt '6 pident length' >"$work/rrnaB.tsv" || true
		rrna_hits_A=$(awk -v L="$RRNA_MINLEN" -v P="$(awk -v x="$RRNA_MINID" 'BEGIN{print 100*x}')" '$1+0>=P && $2+0>=L{c++} END{print c+0}' "$work/rrnaA.tsv")
		rrna_hits_B=$(awk -v L="$RRNA_MINLEN" -v P="$(awk -v x="$RRNA_MINID" 'BEGIN{print 100*x}')" '$1+0>=P && $2+0>=L{c++} END{print c+0}' "$work/rrnaB.tsv")
	elif command -v minimap2 >/dev/null 2>&1; then
		minimap2 -x sr -N 10 -t 1 "$work/IR_A.fa" "$RRNA_REF" 2>/dev/null | awk -v L="$RRNA_MINLEN" '$10+0>=L{c++} END{print c+0}' >"$work/rrnaA.count"
		minimap2 -x sr -N 10 -t 1 "$work/IR_B.fa" "$RRNA_REF" 2>/dev/null | awk -v L="$RRNA_MINLEN" '$10+0>=L{c++} END{print c+0}' >"$work/rrnaB.count"
		rrna_hits_A=$(cat "$work/rrnaA.count")
		rrna_hits_B=$(cat "$work/rrnaB.count")
	fi
fi

choose_A=0
if ((rrna_hits_A > 0 || rrna_hits_B > 0)); then
	if ((rrna_hits_A > 0 && rrna_hits_B == 0 && geom_pref_A == 1)); then
		choose_A=1
	elif ((rrna_hits_B > 0 && rrna_hits_A == 0 && geom_pref_B == 1)); then
		choose_A=0
	else
		diffA=$((before_A - after_A))
		diffB=$((before_B - after_B))
		if ((diffA >= diffB)); then choose_A=1; else choose_A=0; fi
	fi
else
	if ((geom_pref_A == 1 && geom_pref_B == 0)); then
		choose_A=1
	elif ((geom_pref_B == 1 && geom_pref_A == 0)); then
		choose_A=0
	else
		diffA=$((before_A - after_A))
		diffB=$((before_B - after_B))
		if ((diffA >= diffB)); then choose_A=1; else choose_A=0; fi
	fi
fi

if ((choose_A)); then
	: # keep IRa_start/IRa_end as IRa; IRb_* as IRb
else
	tmpS=$IRa_start
	tmpE=$IRa_end
	IRa_start=$IRb_start
	IRa_end=$IRb_end
	IRb_start=$tmpS
	IRb_end=$tmpE
	lenIR_A=$lenIR_B
	qS1_A=$qS1_B
	qE1_A=$qE1_B
	after_A=$after_B
	before_A=$before_B
fi

# ---- Rotate to IRa and equalize IR lengths (common core) ----
G=$(rotate_str "$INSEQ" "$IRa_start") # IRa now at [1..]
N=${#G}
Ia_raw=$((IRa_end - IRa_start + 1))

# partner IRb in rotated frame
qS1=$((((IRb_start - IRa_start + N) % N) + 1))
qE1=$((((IRb_end - IRa_start + N) % N) + 1))

Ib_raw=$(arc_len "$qS1" "$qE1" "$N")
# common core length (preserve total N after printing two equal IRs)
lenIR=$((Ia_raw < Ib_raw ? Ia_raw : Ib_raw))
IRa="${G:0:lenIR}"
IRb_trunc="$(slice_circ_len "$G" "$qS1" "$lenIR")"

# arcs around the trimmed IRs
qE1_trim=$(advance_pos "$qS1" $((lenIR - 1)) "$N")
seg_after="$(slice_circ "$G" $((lenIR + 1)) "$qS1")" # IRa_end -> IRb_start  (candidate SSC)
seg_before="$(slice_circ "$G" $((qE1_trim + 1)) 1)"  # IRb_end -> IRa_start  (candidate LSC)

if ((${#seg_before} >= ${#seg_after})); then
	LSC="$seg_before"
	SSC="$seg_after"
else
	LSC="$seg_after"
	SSC="$seg_before"
fi

# ---- Build forms with primed last IR: IR′ = rc(IRa) ----
LSCp="$(rc_of "$LSC")"
SSCp="$(rc_of "$SSC")"
IRp="$(rc_of "$IRa")" # the printed tail block

FORM1="${LSC}${IRa}${SSC}${IRp}"   # LSC-IR-SSC-IR′
FORM2="${LSC}${IRa}${SSCp}${IRp}"  # LSC-IR-SSC′-IR′
FORM3="${LSCp}${IRa}${SSCp}${IRp}" # LSC′-IR-SSC′-IR′
FORM4="${LSCp}${IRa}${SSC}${IRp}"  # LSC′-IR-SSC-IR′

# lengths for printed form (two equal IRs by construction)
L=${#LSC}
Ia=${#IRa}
S=${#SSC}
Ib=${#IRp}
if ((L + Ia + S + Ib != N)); then
	echo "[WARN] Printed form length $(($L + $Ia + $S + $Ib)) != genome length $N; check IR detection/trim." >&2
fi

# ---- write outputs ----
printf ">%s|LSC-IR-SSC-IRp\n" "$name" >"$OUTDIR/1.fa"
echo "$FORM1" | fold -w 80 >>"$OUTDIR/1.fa"
printf ">%s|LSC-IR-SSCp-IRp\n" "$name" >"$OUTDIR/2.fa"
echo "$FORM2" | fold -w 80 >>"$OUTDIR/2.fa"
printf ">%s|LSCp-IR-SSCp-IRp\n" "$name" >"$OUTDIR/3.fa"
echo "$FORM3" | fold -w 80 >>"$OUTDIR/3.fa"
printf ">%s|LSCp-IR-SSC-IRp\n" "$name" >"$OUTDIR/4.fa"
echo "$FORM4" | fold -w 80 >>"$OUTDIR/4.fa"

echo "[OK] Wrote:"
echo "  $OUTDIR/1.fa"
echo "  $OUTDIR/2.fa"
echo "  $OUTDIR/3.fa"
echo "  $OUTDIR/4.fa"
echo "[INFO] N=$N; LSC=$L; IR=$Ia; SSC=$S (raw IRs: Ia=$Ia_raw, Ib=$Ib_raw; trimmed to $lenIR)"

# ---- thorough check / report ----
if ((CHECK)); then
	[[ -z "$CHECK_REPORT" ]] && CHECK_REPORT="${OUTDIR%/}/pt.check.report.txt"
	report="$CHECK_REPORT"
	: >"$report"
	say() { echo "$@" | tee -a "$report" >&2; }

	# Which outputs are rotations of input?
	declare -a matched_forms=() matched_offsets=()
	for idx in 1 2 3 4; do
		var="FORM$idx"
		s="${!var-}"
		[[ -z "$s" ]] && continue
		off=$(circular_find_offset "$INSEQ" "$s")
		((off > 0)) && {
			matched_forms+=("$idx")
			matched_offsets+=("$off")
		}
	done

	chosen_idx=0
	chosen_off=0
	if ((${#matched_forms[@]})); then
		# prefer FORM1 if it matches
		for i in "${!matched_forms[@]}"; do
			if [[ "${matched_forms[$i]}" == "1" ]]; then
				chosen_idx=1
				chosen_off="${matched_offsets[$i]}"
				break
			fi
		done
		if ((!chosen_idx)); then
			chosen_idx="${matched_forms[0]}"
			chosen_off="${matched_offsets[0]}"
		fi
		var="FORM$chosen_idx"
		s="${!var-}"
		s=$(rotate_str "$s" "$chosen_off")
		printf -v "$var" '%s' "$s"
		say "[CHECK] Matched input to ${OUTDIR}/${chosen_idx}.fa (rotation offset ${chosen_off})."
		if ((${#matched_forms[@]} > 1)); then
			other=()
			for i in "${!matched_forms[@]}"; do
				[[ "${matched_forms[$i]}" == "$chosen_idx" ]] && continue
				other+=("${matched_forms[$i]}@off=${matched_offsets[$i]}")
			done
			say "[CHECK] Other matching forms: ${other[*]}"
		fi
	else
		FORM1="$INSEQ"
		chosen_idx=1
		chosen_off=1
		say "[CHECK] No form was a rotation of the input; wrote input exactly to ${OUTDIR}/1.fa."
	fi

	# Canonical rotation check (use printed canonical FORM1)
	CANON="$FORM1"
	off2=$(circular_find_offset "$INSEQ" "$CANON")
	if ((off2 > 0)); then
		say "[CHECK] Input is a rotation of canonical LSC-IR-SSC-IR′ (offset ${off2})."
	else
		say "[WARN] Input is NOT a rotation of canonical LSC-IR-SSC-IR′."
	fi

	# IR symmetry sanity (genomic partner vs rc(IRa))
	IRb_vs_rc_mism=0
	IRa_rc="$(rc_of "$IRa")"
	if [[ ${#IRa_rc} -eq ${#IRb_trunc} ]]; then
		IRb_vs_rc_mism=$(awk -v A="$IRa_rc" -v B="$IRb_trunc" 'BEGIN{m=0; n=length(A); for(i=1;i<=n;i++) if(substr(A,i,1)!=substr(B,i,1)) m++; print m}')
		say "[CHECK] Genomic IRb_trunc vs rc(IRa): mismatches=${IRb_vs_rc_mism} / ${#IRa_rc} bp."
	else
		say "[WARN] IRb_trunc and rc(IRa) have different lengths (${#IRb_trunc} vs ${#IRa_rc})."
	fi

	# Boundaries for printed form (two equal IRs)
	{
		echo -e "segment\tstart\tend\tlength"
		echo -e "LSC\t1\t${L}\t${L}"
		echo -e "IR\t$((L + 1))\t$((L + Ia))\t${Ia}"
		echo -e "SSC\t$((L + Ia + 1))\t$((L + Ia + S))\t${S}"
		echo -e "IR′\t$((L + Ia + S + 1))\t$((L + Ia + S + Ia))\t${Ia}"
	} >"${OUTDIR%/}/boundaries.tsv"

	# Junction 20-mers from the matched output string
	case "$chosen_idx" in
	1) s="${FORM1}" ;; 2) s="${FORM2}" ;; 3) s="${FORM3}" ;; 4) s="${FORM4}" ;; *) s="${FORM1}" ;;
	esac
	K=20
	S2="${s}${s}"
	p1=$L
	p2=$((L + Ia))
	p3=$((L + Ia + S))
	p4=$((L + Ia + S + Ia))
	b1L=${S2:$((p1 - K)):K}
	b1R=${S2:$p1:K}
	b2L=${S2:$((p2 - K)):K}
	b2R=${S2:$p2:K}
	b3L=${S2:$((p3 - K)):K}
	b3R=${S2:$p3:K}
	b4L=${S2:$((p4 - K)):K}
	b4R=${S2:$p4:K}
	{
		echo "boundary	left_${K}mer	right_${K}mer"
		echo -e "LSC|IR\t$b1L\t$b1R"
		echo -e "IR|SSC\t$b2L\t$b2R"
		echo -e "SSC|IR′\t$b3L\t$b3R"
		echo -e "IR′|LSC\t$b4L\t$b4R"
	} >"${OUTDIR%/}/boundaries.context.txt"

	{
		echo "=== SUMMARY ==="
		echo "Genome length (N): $N"
		echo "LSC=$L, IR_print=$Ia, SSC=$S (raw IRs: Ia=$Ia_raw, Ib=$Ib_raw; trimmed=$lenIR)"
		echo "Chosen match: ${chosen_idx}.fa (offset=$chosen_off)"
		echo "Canonical-rotation offset: $off2"
		echo ""
		echo "=== BOUNDARIES ==="
		cat "${OUTDIR%/}/boundaries.tsv"
		echo ""
		echo "=== JUNCTION 20-mers ==="
		cat "${OUTDIR%/}/boundaries.context.txt"
	} >>"$CHECK_REPORT"

	say "[CHECK] Report written: $CHECK_REPORT"
fi
