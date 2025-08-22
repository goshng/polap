#!/usr/bin/env bash
################################################################################
# This file is part of polap.
#
# polap is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# polap is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# polap. If not, see <https://www.gnu.org/licenses/>.
################################################################################

set -euo pipefail

# align_ptdna.sh
#
# Quiet by default: prints ONLY the final percent identity to stdout.
#   - default format: "99.932%"
#   - with --raw:      "99.932"
# Use --verbose to see progress/info/warnings/errors on stderr.
#
# Pipeline (plastome-specific):
#   1) Build ONE reference plastome from input1 (GFA with 1 or ≥3 segments, or FASTA).
#   2) From input2, enumerate alignable candidates:
#        - GFA (3 seg): two isomers
#            A = LSC + IR + SSC + rc(IR)
#            B = LSC + IR + rc(SSC) + rc(IR)
#          Include reverse-complements → up to 4 candidates total.
#        - GFA (1 seg) or FASTA: forward + RC → 2 candidates.
#   3) For each candidate, BLAST candidate→reference:
#        - Compute REFERENCE coverage using ONLY PLUS-strand HSPs
#          (union on subject coordinates).
#        - Keep candidates with coverage ≥ --min-cov (default 0.90).
#        - Rank kept candidates:
#            primary: plus_cov (higher is better)
#            tie #1 : plus_aligned_bases × plus_weighted_identity (higher)
#            tie #2 : plus_best_bitscore (higher)
#   4) Rotate chosen candidate so REF position 1 aligns to CANDIDATE position 1
#      using BLAST(ref→candidate) to map coordinates.
#   5) MAFFT global alignment (ref, rotated-cand).
#   6) Print global percent identity over UNGAPPED columns.
#
# Requirements: gfatools, seqkit, blastn (and makeblastdb), mafft
#
# Usage:
#   align_ptdna.sh <in1.gfa|fa[.gz]> <in2.gfa|fa[.gz]> \
#       [--out PREFIX] [--threads N] [--min-ir N] [--max-ir N] [--min-cov F] [--raw] [--verbose]

# ---------------------------- CLI & logging ----------------------------
die() {
	echo "[ERROR] $*" >&2
	echo 0
	exit 1
}
log() { ((VERBOSE)) && echo "[INFO] $*" >&2 || true; }
warn() { ((VERBOSE)) && echo "[WARN] $*" >&2 || true; }

if (($# < 2)); then
	echo "Usage: $0 <in1> <in2> [--out PREFIX] [--threads N] [--min-ir N] [--max-ir N] [--min-cov F] [--raw] [--verbose]" >&2
	exit 1
fi

in1="$1"
shift
in2="$1"
shift

prefix="ptalign"
threads="$(nproc 2>/dev/null || echo 4)"
min_ir=8000
max_ir=80000
min_cov="0.90"
VERBOSE=0
RAW=0

while [[ $# -gt 0 ]]; do
	case "$1" in
	--out)
		prefix="${2:-ptalign}"
		shift 2
		;;
	--threads)
		threads="${2:-4}"
		shift 2
		;;
	--min-ir)
		min_ir="${2:-8000}"
		shift 2
		;;
	--max-ir)
		max_ir="${2:-80000}"
		shift 2
		;;
	--min-cov)
		min_cov="${2:-0.90}"
		shift 2
		;;
	--verbose)
		VERBOSE=1
		shift
		;;
	--raw)
		RAW=1
		shift
		;;
	*) die "Unknown option: $1" ;;
	esac
done

need() { command -v "$1" >/dev/null 2>&1 || die "Missing required tool: $1"; }
need gfatools
need seqkit
need blastn
need makeblastdb
need mafft

# run commands silently unless VERBOSE
run_cmd() {
	if ((VERBOSE)); then
		"$@"
	else
		"$@" >/dev/null 2>&1
	fi
}

# ---------------------------- temp dirs ----------------------------
tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT
canddir="$prefix.cand"
mkdir -p "$canddir"

# ---------------------------- helpers ----------------------------
detect_type() { # gfa|fasta|unknown
	local f="$1" line=""
	if [[ "$f" =~ \.gz$ ]]; then
		line="$(zcat -- "$f" 2>/dev/null | sed -n '/^[^#]/p' | head -n1 || true)"
	else
		line="$(sed -n '/^[^#]/p' "$f" | head -n1 || true)"
	fi
	case "$line" in
	H*VN:Z:* | S$'\t'* | L$'\t'* | P$'\t'*) echo gfa ;;
	'>'* | '@'*) echo fasta ;;
	*) echo unknown ;;
	esac
}

get_len_first() { seqkit fx2tab -n -l "$1" 2>/dev/null | awk 'NR==1{print $2; exit}'; }

rc_fa() { seqkit seq -r -p -w 0 "$1" 2>/dev/null; }

rotate_fa_start() { # rotate single-record FASTA so new 1-based start is $2
	local fa="$1" start="${2:-1}" out="$3"
	local hdr seq len a b
	hdr="$(grep '^>' "$fa" | head -n1)"
	seq="$(grep -v '^>' "$fa" | tr -d '\n')"
	len="${#seq}"
	((len > 0)) || die "Empty sequence in $fa"
	[[ "$start" =~ ^[0-9]+$ ]] || start=1
	if ((start <= 1 || start > len)); then
		cp "$fa" "$out"
		return 0
	fi
	a="${seq:start-1}"
	b="${seq:0:start-1}"
	{
		echo "$hdr"
		printf "%s%s\n" "$a" "$b"
	} | seqkit seq -w 0 2>/dev/null >"$out"
}

# ---------------------------- reference build (in1) ----------------------------
gfa_to_plastome_one() {
	local gfa="$1" outfa="$2"
	local sfa="$tmpdir/$(basename "$gfa").S.fa"
	if [[ "$gfa" =~ \.gz$ ]]; then
		zcat -- "$gfa" 2>/dev/null | gfatools gfa2fa - 2>/dev/null >"$sfa"
	else gfatools gfa2fa "$gfa" 2>/dev/null >"$sfa"; fi
	local tab="$tmpdir/$(basename "$gfa").tab"
	seqkit fx2tab -n -l "$sfa" 2>/dev/null >"$tab" || die "Failed to read segments from $gfa"
	local n
	n="$(wc -l <"$tab" | tr -d ' ')"
	: "${n:=0}"

	get_by_id() { seqkit grep -nr -p "^${1}\$" "$sfa" 2>/dev/null | seqkit seq -w 0 2>/dev/null; }

	if [[ "$n" == "1" ]]; then
		local id len
		read -r id len < <(awk '{print $1" "$2}' "$tab")
		[[ "$len" =~ ^[0-9]+$ ]] || die "Non-numeric length in GFA S: '$len'"
		((len >= 100000 && len <= 250000)) || die "Single segment $len not in [100k,250k] ($(basename "$gfa"))"
		get_by_id "$id" | sed "1s/^>.*/>ptdna_${id}/" >"$outfa"
		return 0
	fi

	((n >= 3)) || die "GFA needs 1 or ≥3 S-segments (found $n) in $(basename "$gfa")"
	local top="$tmpdir/$(basename "$gfa").top"
	sort -k2,2nr "$tab" | head -n 6 >"$top"

	awk -vMINIR="$min_ir" -vMAXIR="$max_ir" '
    {ids[NR]=$1; lens[NR]=$2} END{
      m=NR
      for(i=1;i<=m;i++)for(j=i+1;j<=m;j++)for(k=j+1;k<=m;k++){
        L1=lens[i]+0; L2=lens[j]+0; L3=lens[k]+0
        I1=ids[i]; I2=ids[j]; I3=ids[k]
        if(L1<L2){t=L1;L1=L2;L2=t; u=I1;I1=I2;I2=u}
        if(L2<L3){t=L2;L2=L3;L3=t; u=I2;I2=I3;I3=u}
        if(L1<L2){t=L1;L1=L2;L2=t; u=I1;I1=I2;I2=u}
        LSC=L1; IR=L2; SSC=L3; total=LSC+SSC+2*IR
        if(total>=100000 && total<=250000 && IR>=MINIR && IR<=MAXIR) printf "%s\t%s\t%s\n", I1,I2,I3
      }
    }' "$top" | head -n1 >"$tmpdir/candidate.ref.txt"

	[[ -s "$tmpdir/candidate.ref.txt" ]] || die "No plausible (LSC,IR,SSC) triplet in $(basename "$gfa")"
	local LSC_ID IR_ID SSC_ID
	read -r LSC_ID IR_ID SSC_ID <"$tmpdir/candidate.ref.txt"

	local LSC IR SSC IRrc
	LSC="$(get_by_id "$LSC_ID" | sed '1d' | tr -d '\n')"
	IR="$(get_by_id "$IR_ID" | sed '1d' | tr -d '\n')"
	SSC="$(get_by_id "$SSC_ID" | sed '1d' | tr -d '\n')"
	printf ">tmp\n%s\n" "$IR" | seqkit seq -r -p -w 0 2>/dev/null | sed '1d' | tr -d '\n' >"$tmpdir/IRrc.ref.txt"
	IRrc="$(cat "$tmpdir/IRrc.ref.txt")"

	{
		echo ">ptdna_${LSC_ID}_${IR_ID}_${SSC_ID}"
		printf "%s%s%s%s\n" "$LSC" "$IR" "$SSC" "$IRrc"
	} |
		seqkit seq -w 0 2>/dev/null >"$outfa"

	local len
	len="$(get_len_first "$outfa")"
	: "${len:=0}"
	[[ "$len" =~ ^[0-9]+$ ]] || die "Could not parse reconstructed length: '$len'"
	((len >= 100000 && len <= 250000)) || die "Reconstructed length $len not in [100k,250k] ($(basename "$gfa"))"
}

fasta_to_plastome_one() {
	local fa="$1" outfa="$2"
	local pick id len
	pick="$(seqkit fx2tab -n -l "$fa" 2>/dev/null | sort -k2,2nr | awk 'NR==1{print $1" "$2}')"
	[[ -n "$pick" ]] || die "No records in $fa"
	read -r id len <<<"$pick"
	[[ "$len" =~ ^[0-9]+$ ]] || die "Non-numeric FASTA length: '$len'"
	((len >= 100000 && len <= 250000)) || die "Top FASTA record $len not in [100k,250k] ($(basename "$fa"))"
	seqkit grep -nr -p "^${id}\$" "$fa" 2>/dev/null | seqkit seq -w 0 2>/dev/null | sed "1s/^>.*/>ptdna_${id}/" >"$outfa"
}

# ---------------------------- candidate build (in2) ----------------------------
gfa_to_all_candidates() {
	local gfa="$1" outprefix="$2"
	local sfa="$tmpdir/$(basename "$gfa").S2.fa"
	if [[ "$gfa" =~ \.gz$ ]]; then
		zcat -- "$gfa" 2>/dev/null | gfatools gfa2fa - 2>/dev/null >"$sfa"
	else gfatools gfa2fa "$gfa" 2>/dev/null >"$sfa"; fi

	local tab="$tmpdir/$(basename "$gfa").tab2"
	seqkit fx2tab -n -l "$sfa" 2>/dev/null >"$tab"
	local n
	n="$(wc -l <"$tab" | tr -d ' ')"
	: "${n:=0}"

	get_by_id_fa() { seqkit grep -nr -p "^${1}\$" "$sfa" 2>/dev/null | seqkit seq -w 0 2>/dev/null; }

	if [[ "$n" == "1" ]]; then
		local id len
		read -r id len < <(awk '{print $1" "$2}' "$tab")
		[[ "$len" =~ ^[0-9]+$ ]] || die "Non-numeric length in GFA S: '$len'"
		((len >= 100000 && len <= 250000)) || die "Single segment $len not in [100k,250k] ($(basename "$gfa"))"
		get_by_id_fa "$id" | sed "1s/^>.*/>${outprefix}_A/" >"$canddir/${outprefix}.A.fw.fa"
		rc_fa "$canddir/${outprefix}.A.fw.fa" >"$canddir/${outprefix}.A.rc.fa"
		return 0
	fi

	((n >= 3)) || die "GFA needs 1 or ≥3 S-segments (found $n) in $(basename "$gfa")"
	local top="$tmpdir/$(basename "$gfa").top2"
	sort -k2,2nr "$tab" | head -n 6 >"$top"

	awk -vMINIR="$min_ir" -vMAXIR="$max_ir" '
    {ids[NR]=$1; lens[NR]=$2} END{
      m=NR
      for(i=1;i<=m;i++)for(j=i+1;j<=m;j++)for(k=j+1;k<=m;k++){
        L1=lens[i]+0; L2=lens[j]+0; L3=lens[k]+0
        I1=ids[i]; I2=ids[j]; I3=ids[k]
        if(L1<L2){t=L1;L1=L2;L2=t; u=I1;I1=I2;I2=u}
        if(L2<L3){t=L2;L2=L3;L3=t; u=I2;I2=I3;I3=u}
        if(L1<L2){t=L1;L1=L2;L2=t; u=I1;I1=I2;I2=u}
        LSC=L1; IR=L2; SSC=L3; total=LSC+SSC+2*IR
        if(total>=100000 && total<=250000 && IR>=MINIR && IR<=MAXIR) printf "%s\t%s\t%s\n", I1,I2,I3
      }
    }' "$top" | head -n1 >"$tmpdir/cand3.txt"

	[[ -s "$tmpdir/cand3.txt" ]] || die "No plausible (LSC,IR,SSC) triplet in $(basename "$gfa")"
	local LSC_ID IR_ID SSC_ID
	read -r LSC_ID IR_ID SSC_ID <"$tmpdir/cand3.txt"

	local LSC IR SSC IRrc SSCrc
	LSC="$(get_by_id_fa "$LSC_ID" | sed '1d' | tr -d '\n')"
	IR="$(get_by_id_fa "$IR_ID" | sed '1d' | tr -d '\n')"
	SSC="$(get_by_id_fa "$SSC_ID" | sed '1d' | tr -d '\n')"
	printf ">tmp\n%s\n" "$IR" | seqkit seq -r -p -w 0 2>/dev/null | sed '1d' | tr -d '\n' >"$tmpdir/IRrc2.txt"
	IRrc="$(cat "$tmpdir/IRrc2.txt")"
	printf ">tmp\n%s\n" "$SSC" | seqkit seq -r -p -w 0 2>/dev/null | sed '1d' | tr -d '\n' >"$tmpdir/SSCrc2.txt"
	SSCrc="$(cat "$tmpdir/SSCrc2.txt")"

	# A: LSC + IR + SSC + rc(IR)
	{
		echo ">${outprefix}_A"
		printf "%s%s%s%s\n" "$LSC" "$IR" "$SSC" "$IRrc"
	} | seqkit seq -w 0 2>/dev/null >"$canddir/${outprefix}.A.fw.fa"
	rc_fa "$canddir/${outprefix}.A.fw.fa" >"$canddir/${outprefix}.A.rc.fa"

	# B: LSC + IR + rc(SSC) + rc(IR)  (SSC-flipped isomer)
	{
		echo ">${outprefix}_B"
		printf "%s%s%s%s\n" "$LSC" "$IR" "$SSCrc" "$IRrc"
	} | seqkit seq -w 0 2>/dev/null >"$canddir/${outprefix}.B.fw.fa"
	rc_fa "$canddir/${outprefix}.B.fw.fa" >"$canddir/${outprefix}.B.rc.fa"
}

fasta_to_all_candidates() {
	local fa="$1" outprefix="$2"
	local tmpfa="$tmpdir/${outprefix}.long.fa"
	local pick id len
	pick="$(seqkit fx2tab -n -l "$fa" 2>/dev/null | sort -k2,2nr | awk 'NR==1{print $1" "$2}')"
	[[ -n "$pick" ]] || die "No records in $fa"
	read -r id len <<<"$pick"
	[[ "$len" =~ ^[0-9]+$ ]] || die "Non-numeric FASTA length: '$len'"
	((len >= 100000 && len <= 250000)) || die "Top FASTA record $len not in [100k,250k] ($(basename "$fa"))"
	seqkit grep -nr -p "^${id}\$" "$fa" 2>/dev/null | seqkit seq -w 0 2>/dev/null | sed "1s/^>.*/>${outprefix}_A/" >"$tmpfa"
	cp "$tmpfa" "$canddir/${outprefix}.A.fw.fa"
	rc_fa "$tmpfa" >"$canddir/${outprefix}.A.rc.fa"
}

# ---------------------------- build reference & candidates ----------------------------
log "Selecting reference plastome from input1..."
case "$(detect_type "$in1")" in
gfa) gfa_to_plastome_one "$in1" "$prefix.ref.fa" ;;
fasta) fasta_to_plastome_one "$in1" "$prefix.ref.fa" ;;
*) die "Unknown type for $in1" ;;
esac

log "Enumerating candidates from input2..."
case "$(detect_type "$in2")" in
gfa) gfa_to_all_candidates "$in2" "cand" ;;
fasta) fasta_to_all_candidates "$in2" "cand" ;;
*) die "Unknown type for $in2" ;;
esac
ls -1 "$canddir"/*.fa >/dev/null 2>&1 || die "No candidates produced from input2"

# ---------------------------- BLAST selection ----------------------------
log "Scoring candidates via BLAST against reference (ranking by plus-strand coverage ≥ ${min_cov})..."
run_cmd makeblastdb -in "$prefix.ref.fa" -dbtype nucl

# header (kept for debugging file; not printed)
echo -e "candidate\tplus_cov\tplus_len\tplus_wid\tplus_bestbits\tall_aln" >"$prefix.blast.choice.tsv"

for cand in "$canddir"/*.fa; do
	out="$tmpdir/$(basename "$cand").choice.tsv"
	run_cmd blastn -task megablast -query "$cand" -db "$prefix.ref.fa" \
		-num_threads "$threads" \
		-outfmt '6 qstart qend sstart send sstrand length pident qlen slen bitscore' \
		-out "$out" || true

	if [[ ! -s "$out" ]]; then
		echo -e "$(basename "$cand")\t0.000\t0\t0.000\t0\t0" >>"$prefix.blast.choice.tsv"
		continue
	fi

	all_aln="$(awk '{s+=$6} END{print s+0}' "$out")"

	plus_file="$tmpdir/plus.$(basename "$cand").tsv"
	awk '$5=="plus"' "$out" >"$plus_file" || true
	plus_cov="0.000"
	plus_len=0
	plus_pid=0
	plus_best=0
	if [[ -s "$plus_file" ]]; then
		read -r plus_len plus_pid plus_best <<<"$(awk 'BEGIN{s=0;w=0;best=0}{s+=$6; w+=$6*$7; if($10>best)best=$10}END{printf "%.0f %.6f %.0f\n", s, (s>0?w/s:0), best}' "$plus_file")"
		slen="$(awk 'NR==1{print $9; exit}' "$plus_file")"
		: "${slen:=0}"
		[[ "$slen" =~ ^[0-9]+$ ]] || slen=0
		ivals="$tmpdir/ivals.$(basename "$cand").bed"
		awk '{s=($3<$4)?$3:$4; e=($3<$4)?$4:$3; print s"\t"e}' "$plus_file" | sort -k1,1n -k2,2n >"$ivals"
		merged="$tmpdir/merged.$(basename "$cand").bed"
		awk 'NR==1{S=$1;E=$2;next}{if($1<=E){if($2>E)E=$2}else{print S"\t"E;S=$1;E=$2}}END{if(NR>0)print S"\t"E}' "$ivals" >"$merged" || true
		covbp=0
		[[ -s "$merged" ]] && covbp="$(awk '{c+=($2-$1+1)} END{print c+0}' "$merged")"
		[[ "$slen" -gt 0 ]] && plus_cov="$(awk -v c="$covbp" -v s="$slen" 'BEGIN{printf "%.3f", (s>0?c/s:0)}')"
	fi

	echo -e "$(basename "$cand")\t${plus_cov}\t${plus_len}\t${plus_pid}\t${plus_best}\t${all_aln}" >>"$prefix.blast.choice.tsv"
done

# filter & choose best
mapfile -t pass_list < <(awk -v t="$min_cov" -F'\t' 'NR>1 && ($2+0.0)>=t {print $1}' "$prefix.blast.choice.tsv")
[[ "${#pass_list[@]}" -gt 0 ]] || die "No candidate reached plus-strand reference coverage ≥ ${min_cov}."

best_cand=""
best_cov="-1"
best_primary=0
best_bits=0
for name in "${pass_list[@]}"; do
	row="$(awk -v n="$name" -F'\t' '$1==n{print; exit}' "$prefix.blast.choice.tsv")"
	pcov="$(awk -F'\t' '{print $2}' <<<"$row")"
	plen="$(awk -F'\t' '{print $3}' <<<"$row")"
	ppid="$(awk -F'\t' '{print $4}' <<<"$row")"
	pbits="$(awk -F'\t' '{print $5}' <<<"$row")"
	pprimary="$(awk -v L="$plen" -v P="$ppid" 'BEGIN{printf "%.6f", L*P}')"

	if [[ -z "${best_cand:-}" ]]; then
		best_cand="$name"
		best_cov="$pcov"
		best_primary="$pprimary"
		best_bits="$pbits"
	else
		# prefer higher coverage; then higher L*PID; then higher bitscore
		if awk -v bc="$best_cov" -v pc="$pcov" -v bp="$best_primary" -v pp="$pprimary" -v bb="$best_bits" -v pb="$pbits" 'BEGIN{
          if ((pc - bc) > 1e-9) exit 0;
          if ((bc - pc) > 1e-9) exit 1;
          if ((pp - bp) > 1e-9) exit 0;
          if ((bp - pp) > 1e-9) exit 1;
          if ((pb - bb) > 1e-9) exit 0;
          exit 1
        }'; then
			best_cand="$name"
			best_cov="$pcov"
			best_primary="$pprimary"
			best_bits="$pbits"
		fi
	fi
done

cp "$canddir/$best_cand" "$prefix.best.raw.fa"

# ---------------------------- rotation (ref→candidate) ----------------------------
log "Mapping ref position 1 to chosen candidate for rotation..."
run_cmd blastn -task megablast -query "$prefix.ref.fa" -subject "$prefix.best.raw.fa" \
	-num_threads "$threads" \
	-outfmt '6 qstart qend sstart send sstrand length pident qlen slen bitscore' \
	-out "$prefix.blast.rotate.tsv" || true

if [[ ! -s "$prefix.blast.rotate.tsv" ]]; then
	cp "$prefix.best.raw.fa" "$prefix.best.rotated.fa"
else
	rot_start="$(
		awk '
      function abs(x){return x<0?-x:x}
      BEGIN{bestCoverLen=-1; bestNearDist=1e18; selS=1; haveCover=0}
      {
        q1=$1+0; q2=$2+0; s1=$3+0; s2=$4+0; strand=$5; L=$6+0; qlen=$9+0; slen=$10+0
        qlo=(q1<q2?q1:q2); qhi=(q1>q2?q1:q2);
        if (qlo<=1 && 1<=qhi){
          off = 1 - q1;
          if (strand=="plus"){ s = s1 + off } else { s = s2 - off }
          if (s<1) s=1; if (s>slen) s=slen;
          if (L>bestCoverLen){bestCoverLen=L; selS=s; haveCover=1}
        } else {
          d = (abs(1-qlo) < abs(1-qhi)) ? abs(1-qlo) : abs(1-qhi);
          if (d < bestNearDist && haveCover==0){
            if (abs(1-qlo) < abs(1-qhi)){
              if (strand=="plus"){ s = s1 } else { s = s2 }
            } else {
              if (strand=="plus"){ s = s2 } else { s = s1 }
            }
            if (s<1) s=1; if (s>slen) s=slen;
            bestNearDist=d; selS=s;
          }
        }
      }
      END{print (selS>0?selS:1)}
    ' "$prefix.blast.rotate.tsv"
	)"
	: "${rot_start:=1}"
	rotate_fa_start "$prefix.best.raw.fa" "$rot_start" "$prefix.best.rotated.fa"
fi

# ---------------------------- MAFFT (quiet) ----------------------------
log "Running MAFFT..."
# per request: only --thread and --quiet
if ((VERBOSE)); then
	mafft --thread "$threads" --quiet \
		<(cat "$prefix.ref.fa" "$prefix.best.rotated.fa") >"$prefix.mafft.aln.fasta"
else
	mafft --thread "$threads" --quiet \
		<(cat "$prefix.ref.fa" "$prefix.best.rotated.fa") >"$prefix.mafft.aln.fasta" 2>/dev/null
fi

# ---------------------------- Global % identity (UNGAPPED) ----------------------------
if [[ -s "$prefix.mafft.aln.fasta" ]]; then
	seqkit fx2tab "$prefix.mafft.aln.fasta" 2>/dev/null >"$tmpdir/aln.fx2tab" || true
	nseq="$(wc -l <"$tmpdir/aln.fx2tab" | tr -d ' ')"
	if ((nseq != 2)); then
		warn "Alignment should contain exactly 2 sequences, got $nseq"
		printf "%s\n" "NA"
		exit 1
	fi
	# compute PID and print only the number (with % unless --raw)
	awk -F '\t' -v RAW="$RAW" '
    NR==1 { s1=$2; next }
    NR==2 { s2=$2; }
    END {
      len1=length(s1); len2=length(s2);
      if (len1 != len2) { print "NA"; exit 2 }
      for (i=1; i<=len1; i++) {
        a = substr(s1,i,1); b = substr(s2,i,1);
        A = toupper(a); B = toupper(b);
        if (a != "-" && b != "-") { ungap++; if (A == B) match_ungap++; }
      }
      pid = (ungap>0)?(100.0*match_ungap/ungap):0.0;
      if (RAW==1) { printf("%.3f\n", pid); }
      else        { printf("%.3f%%\n", pid); }
    }
  ' "$tmpdir/aln.fx2tab"
else
	printf "%s\n" "NA"
	exit 1
fi
