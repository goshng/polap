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

################################################################################
# Convert numbers between different units.
################################################################################

################################################################################
# Ensure that the current script is sourced only once
source "${_POLAPLIB_DIR}/run-polap-function-include.sh"
_POLAP_INCLUDE_=$(_polap_include "${BASH_SOURCE[0]}")
set +u
if [[ -n "${!_POLAP_INCLUDE_}" ]]; then
	set -u
	return 0
fi
set -u
declare "$_POLAP_INCLUDE_=1"
#
################################################################################

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
	echo "[ERROR] This script must be sourced, not executed: use 'source $BASH_SOURCE'" >&2
	return 1 2>/dev/null || exit 1
fi
: "${_POLAP_DEBUG:=0}"
: "${_POLAP_RELEASE:=0}"

################################################################################
# Function to convert base pairs to the highest appropriate unit
# Example usage
# bp=31846726397
# convert_bp $bp
################################################################################
# FIXME: rewrite the gfa2fa with seed contig
# 1. find the path for the plastid genome sequence
#   1.1 find the longest contig with PT genes
#   1.2 find the path starting at the contig
# 2. stitch the genome sequence of the path from the gfa: LSC-IR-SSC-IR'
# 3. create 1 from ref and 4 from the query

# input: gfa vs. gfa
# input: fasta vs. gfa
# input: fasta vs. fasta
# input: gfa vs. fasta
#
# reference can be fasta or gfa: rearrange so that LSC starts first.
# if query is a fasta, then warnings are necessary because rotation is possible,
# but rearranging LSC, IR, SSC may not biologically make sense.
# However, we could still compute the percent identity anyway.
#
# reference: gfa or fasta -> make LSC start
# query: gfa or fasta -> 4 sequences with LSC starting

function _polap_lib_salign-pt {
	_polap_log0 "Align two ptDNAs"

	if (($# < 2)); then
		echo "Usage: $0 <in1.gfa|fa[.gz]> <in2.gfa|fa[.gz]> [--out PREFIX|DIR/PREFIX] [--threads N] [--min-ir N] [--max-ir N] [--min-cov F] [--raw] [--verbose] [--tmp-in-out]" >&2
		exit 1
	fi

	in1="$1"
	shift
	in2="$1"
	shift

	# Defaults
	prefix="ptalign"
	threads="$(nproc 2>/dev/null || echo 4)"
	min_ir=8000
	max_ir=80000
	min_cov="0.90"
	VERBOSE=0
	RAW=0
	TMP_IN_OUT=0

	# Where AWK helpers live
	SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/scripts"

	# Parse options
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
		--tmp-in-out)
			TMP_IN_OUT=1
			shift
			;;
		*) die "Unknown option: $1" ;;
		esac
	done
}

# alignpt.sh
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
#          Include reverse-complements -> up to 4 candidates total.
#        - GFA (1 seg) or FASTA: forward + RC -> 2 candidates.
#   3) For each candidate, BLAST candidate -> eference:
#        - Compute REFERENCE coverage using ONLY PLUS-strand HSPs
#          (union on subject coordinates).
#        - Keep candidates with coverage ≥ --min-cov (default 0.90).
#        - Rank kept candidates:
#            primary: plus_cov (higher is better)
#            tie #1 : plus_aligned_bases × plus_weighted_identity (higher)
#            tie #2 : plus_best_bitscore (higher)
#   4) Rotate chosen candidate so REF position 1 aligns to CANDIDATE position 1
#      using BLAST(ref -> andidate) to map coordinates.
#   5) MAFFT global alignment (ref, rotated-cand).
#   6) Print global percent identity over UNGAPPED columns.
#
# Requirements: gfatools, seqkit, blastn (and makeblastdb), mafft
#
# Usage:
#   alignpt.sh <in1.gfa|fa[.gz]> <in2.gfa|fa[.gz]> #       [--out PREFIX] [--threads N] [--min-ir N] [--max-ir N] [--min-cov F] [--raw] [--verbose]

# alignpt.sh
#
# Quiet by default: prints ONLY the final percent identity to stdout.
#   - default format: "99.932%"
#   - with --raw:      "99.932"
# Use --verbose to see progress/info/warnings/errors on stderr.
#
# Requirements: gfatools, seqkit, blastn (makeblastdb), mafft
#
# Notable options:
#   --out OUTSPEC     OUTSPEC can include a path, e.g. out1/pt -> makes out1/, prefix=pt, outputs out1/pt.*
#   --tmp-in-out      Use OUTDIR/tmp as tmpdir and DO NOT delete it after run (kept for inspection)

die() {
	echo "[ERROR] $*" >&2
	echo 0
	exit 1
}
log() { ((VERBOSE)) && echo "[INFO] $*" >&2 || true; }
warn() { ((VERBOSE)) && echo "[WARN] $*" >&2 || true; }

if (($# < 2)); then
	echo "Usage: $0 <in1.gfa|fa[.gz]> <in2.gfa|fa[.gz]> [--out PREFIX|DIR/PREFIX] [--threads N] [--min-ir N] [--max-ir N] [--min-cov F] [--raw] [--verbose] [--tmp-in-out]" >&2
	exit 1
fi

in1="$1"
shift
in2="$1"
shift

# Defaults
prefix="ptalign"
threads="$(nproc 2>/dev/null || echo 4)"
min_ir=8000
max_ir=80000
min_cov="0.90"
VERBOSE=0
RAW=0
TMP_IN_OUT=0

# Where AWK helpers live
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/scripts"

# Parse options
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
	--tmp-in-out)
		TMP_IN_OUT=1
		shift
		;;
	*) die "Unknown option: $1" ;;
	esac
done

# Tools
need() { command -v "$1" >/dev/null 2>&1 || die "Missing required tool: $1"; }
need gfatools
need seqkit
need blastn
need makeblastdb
need mafft
[[ -s "$SCRIPTDIR/select_triplet.awk" ]] || die "Missing $SCRIPTDIR/select_triplet.awk"
[[ -s "$SCRIPTDIR/merge_intervals.awk" ]] || die "Missing $SCRIPTDIR/merge_intervals.awk"
[[ -s "$SCRIPTDIR/rotation_anchor.awk" ]] || die "Missing $SCRIPTDIR/rotation_anchor.awk"
[[ -s "$SCRIPTDIR/compute_pid.awk" ]] || die "Missing $SCRIPTDIR/compute_pid.awk"
[[ -s "$SCRIPTDIR/plus_stats.awk" ]] || die "Missing $SCRIPTDIR/plus_stats.awk"

# Handle --out possibly containing a path
outdir=""
if [[ "$prefix" == */* ]]; then
	outdir="$(dirname -- "$prefix")"
	prefix="$(basename -- "$prefix")"
	mkdir -p -- "$outdir"
fi
# All product files use workprefix (which includes outdir if set)
workprefix="${outdir:+$outdir/}$prefix"

# Temp dir setup
cleanup() { rm -rf "$tmpdir"; }
if ((TMP_IN_OUT)); then
	# Must have an outdir to place tmp in
	[[ -n "${outdir:-}" ]] || die "--tmp-in-out requires --out DIR/PREFIX (so we know OUTDIR)"
	tmpdir="$outdir/tmp"
	mkdir -p "$tmpdir"
	# Do NOT clean up
else
	tmpdir="$(mktemp -d)"
	trap cleanup EXIT
fi

# Candidate dir lives next to workprefix
canddir="${workprefix}.cand"
mkdir -p "$canddir"

# Runner
run_cmd() {
	if ((VERBOSE)); then
		"$@"
	else
		"$@" >/dev/null 2>&1
	fi
}

# ---------- helpers ----------
detect_type() {
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

get_len_first() {
	seqkit fx2tab -n -l "$1" 2>/dev/null | awk 'NR==1{print $2; exit}'
}

rc_fa() {
	seqkit seq -r -p -w 0 "$1" 2>/dev/null
}

rotate_fa_start() {
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

# ---------- reference from input1 ----------
gfa_to_plastome_one() {
	local gfa="$1" outfa="$2"
	local sfa="$tmpdir/$(basename "$gfa").S.fa"

	if [[ "$gfa" =~ \.gz$ ]]; then
		zcat -- "$gfa" 2>/dev/null | gfatools gfa2fa - 2>/dev/null >"$sfa"
	else
		gfatools gfa2fa "$gfa" 2>/dev/null >"$sfa"
	fi

	local tab="$tmpdir/$(basename "$gfa").tab"

	seqkit fx2tab -n -l "$sfa" 2>/dev/null >"$tab" || die "Failed to read segments from $gfa"

	local n
	n="$(wc -l <"$tab" | tr -d ' ')"
	: "${n:=0}"

	get_by_id() {
		seqkit grep -nr -p "^${1}\$" "$sfa" 2>/dev/null | seqkit seq -w 0 2>/dev/null
	}

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

	# FIXME: LSC-IR-SSC selection is incorrect.
	# Use the gfa graph structure not the order of the lengths of the segments.
	awk -vMINIR="$min_ir" -vMAXIR="$max_ir" -f "$SCRIPTDIR/select_triplet.awk" "$top" | head -n1 >"$tmpdir/candidate.ref.txt"
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

# ---------- candidates from input2 ----------
gfa_to_all_candidates() {
	local gfa="$1" outprefix="$2"
	local sfa="$tmpdir/$(basename "$gfa").S2.fa"
	if [[ "$gfa" =~ \.gz$ ]]; then
		zcat -- "$gfa" 2>/dev/null | gfatools gfa2fa - 2>/dev/null >"$sfa"
	else
		gfatools gfa2fa "$gfa" 2>/dev/null >"$sfa"
	fi

	local tab="$tmpdir/$(basename "$gfa").tab2"
	seqkit fx2tab -n -l "$sfa" 2>/dev/null >"$tab"

	local n
	n="$(wc -l <"$tab" | tr -d ' ')"
	: "${n:=0}"
	get_by_id_fa() {
		seqkit grep -nr -p "^${1}\$" "$sfa" 2>/dev/null | seqkit seq -w 0 2>/dev/null
	}

	if [[ "$n" == "1" ]]; then
		local id len
		read -r id len < <(awk '{print $1" "$2}' "$tab")
		[[ "$len" =~ ^[0-9]+$ ]] || die "Non-numeric length in GFA S: '$len'"
		((len >= 100000 && len <= 250000)) || die "Single segment $len not in [100k,250k] ($(basename "$gfa"))"
		get_by_id_fa "$id" | sed "1s/^>.*/>${outprefix}_A/" >"${canddir}/${outprefix}.A.fw.fa"
		rc_fa "${canddir}/${outprefix}.A.fw.fa" >"${canddir}/${outprefix}.A.rc.fa"
		return 0
	fi

	((n >= 3)) || die "GFA needs 1 or ≥3 S-segments (found $n) in $(basename "$gfa")"

	local top="$tmpdir/$(basename "$gfa").top2"
	sort -k2,2nr "$tab" | head -n 6 >"$top"

	# FIXME: LSC-IR-SSC selection is incorrect.
	# Use the gfa graph structure not the order of the lengths of the segments.
	awk -vMINIR="$min_ir" -vMAXIR="$max_ir" -f "$SCRIPTDIR/select_triplet.awk" "$top" | head -n1 >"$tmpdir/cand3.txt"
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

	{
		echo ">${outprefix}_A"
		printf "%s%s%s%s\n" "$LSC" "$IR" "$SSC" "$IRrc"
	} | seqkit seq -w 0 2>/dev/null >"${canddir}/${outprefix}.A.fw.fa"
	rc_fa "${canddir}/${outprefix}.A.fw.fa" >"${canddir}/${outprefix}.A.rc.fa"

	{
		echo ">${outprefix}_B"
		printf "%s%s%s%s\n" "$LSC" "$IR" "$SSCrc" "$IRrc"
	} | seqkit seq -w 0 2>/dev/null >"${canddir}/${outprefix}.B.fw.fa"
	rc_fa "${canddir}/${outprefix}.B.fw.fa" >"${canddir}/${outprefix}.B.rc.fa"
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
	cp "$tmpfa" "${canddir}/${outprefix}.A.fw.fa"
	rc_fa "$tmpfa" >"${canddir}/${outprefix}.A.rc.fa"
}

################################################################################
# --- MAIN ---
################################################################################

# ---------- Build reference & candidates ----------
case "$(detect_type "$in1")" in
gfa)
	gfa_to_plastome_one "$in1" "${workprefix}.ref.fa"
	;;
fasta)
	fasta_to_plastome_one "$in1" "${workprefix}.ref.fa"
	;;
*)
	die "Unknown type for $in1"
	;;
esac

case "$(detect_type "$in2")" in
gfa) gfa_to_all_candidates "$in2" "cand" ;;
fasta) fasta_to_all_candidates "$in2" "cand" ;;
*) die "Unknown type for $in2" ;;
esac
ls -1 "$canddir"/*.fa >/dev/null 2>&1 || die "No candidates produced from input2"

# FIXME: rewrite the gfa2fa with seed contig
# 1. find the path for the plastid genome sequence
#   1.1 find the longest contig with PT genes
#   1.2 find the path starting at the contig
# 2. stitch the genome sequence of the path from the gfa: LSC-IR-SSC-IR'
# 3. create 1 from ref and 4 from the query

# input: gfa vs. gfa
# input: fasta vs. gfa
# input: fasta vs. fasta
# input: gfa vs. fasta
#
# reference can be fasta or gfa: rearrange so that LSC starts first.
# if query is a fasta, then warnings are necessary because rotation is possible,
# but rearranging LSC, IR, SSC may not biologically make sense.
# However, we could still compute the percent identity anyway.
#
# reference: gfa or fasta -> make LSC start
# query: gfa or fasta -> 4 sequences with LSC starting

# ---------- BLAST candidate -> eference, rank by plus-only ref coverage ----------
run_cmd makeblastdb -in "${workprefix}.ref.fa" -dbtype nucl
echo -e "candidate\tplus_cov\tplus_len\tplus_wid\tplus_bestbits\tall_aln" >"${workprefix}.blast.choice.tsv"

for cand in "$canddir"/*.fa; do
	out="$tmpdir/$(basename "$cand").choice.tsv"
	echo "run blast: $cand"
	echo blastn -task megablast -query "$cand" -db "${workprefix}.ref.fa" \
		-num_threads "$threads" \
		-outfmt '6 qstart qend sstart send sstrand length pident qlen slen bitscore' \
		>"$out.sh" || true
	blastn -task megablast -query "$cand" -db "${workprefix}.ref.fa" \
		-num_threads "$threads" \
		-outfmt '6 qstart qend sstart send sstrand length pident qlen slen bitscore' \
		>"$out" || true

	if [[ ! -s "$out" ]]; then
		echo -e "$(basename "$cand")\t0.000\t0\t0.000\t0\t0" >>"${workprefix}.blast.choice.tsv"
		continue
	fi

	all_aln="$(awk -f "$SCRIPTDIR/plus_stats.awk" -v MODE=all "$out")"

	plus_file="$tmpdir/plus.$(basename "$cand").tsv"
	awk '$5=="plus"' "$out" >"$plus_file" || true
	plus_cov="0.000"
	plus_len=0
	plus_pid=0
	plus_best=0
	if [[ -s "$plus_file" ]]; then
		read -r plus_len plus_pid plus_best <<<"$(awk -f "$SCRIPTDIR/plus_stats.awk" -v MODE=plus "$plus_file")"
		slen="$(awk 'NR==1{print $9; exit}' "$plus_file")"
		: "${slen:=0}"
		[[ "$slen" =~ ^[0-9]+$ ]] || slen=0
		ivals="$tmpdir/ivals.$(basename "$cand").bed"
		awk '{s=($3<$4)?$3:$4; e=($3<$4)?$4:$3; print s"\t"e}' "$plus_file" | sort -k1,1n -k2,2n >"$ivals"
		merged="$tmpdir/merged.$(basename "$cand").bed"
		awk -f "$SCRIPTDIR/merge_intervals.awk" "$ivals" >"$merged" || true
		covbp=0
		[[ -s "$merged" ]] && covbp="$(awk '{c+=($2-$1+1)} END{print c+0}' "$merged")"
		[[ "$slen" -gt 0 ]] && plus_cov="$(awk -v c="$covbp" -v s="$slen" 'BEGIN{printf "%.3f", (s>0?c/s:0)}')"
	fi

	echo -e "$(basename "$cand")\t${plus_cov}\t${plus_len}\t${plus_pid}\t${plus_best}\t${all_aln}" >>"${workprefix}.blast.choice.tsv"
done

# Choose best candidate
mapfile -t pass_list < <(awk -v t="$min_cov" -F'\t' 'NR>1 && ($2+0.0)>=t {print $1}' "${workprefix}.blast.choice.tsv")
[[ "${#pass_list[@]}" -gt 0 ]] || die "No candidate reached plus-strand reference coverage ≥ ${min_cov}."

best_cand=""
best_cov="-1"
best_primary=0
best_bits=0
for name in "${pass_list[@]}"; do
	IFS=$'\t' read -r _ pcov plen ppid pbits _ <<<"$(awk -v n="$name" -F'\t' '$1==n{print $0; exit}' "${workprefix}.blast.choice.tsv")"
	pprimary="$(awk -v L="$plen" -v P="$ppid" 'BEGIN{printf "%.6f", L*P}')"
	if [[ -z "${best_cand:-}" ]]; then
		best_cand="$name"
		best_cov="$pcov"
		best_primary="$pprimary"
		best_bits="$pbits"
	else
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

cp "$canddir/$best_cand" "${workprefix}.best.raw.fa"

# ---------- Rotate chosen candidate so ref pos1 aligns to cand pos1 ----------
blastn -task megablast -query "${workprefix}.ref.fa" -subject "${workprefix}.best.raw.fa" \
	-num_threads "$threads" \
	-outfmt '6 qstart qend sstart send sstrand length pident qlen slen bitscore' \
	>"${workprefix}.blast.rotate.tsv" || true

if [[ ! -s "${workprefix}.blast.rotate.tsv" ]]; then
	cp "${workprefix}.best.raw.fa" "${workprefix}.best.rotated.fa"
else
	rot_start="$(awk -f "$SCRIPTDIR/rotation_anchor.awk" "${workprefix}.blast.rotate.tsv")"
	: "${rot_start:=1}"
	rotate_fa_start "${workprefix}.best.raw.fa" "$rot_start" "${workprefix}.best.rotated.fa"
fi

# ---------- MAFFT (quiet) ----------
if ((VERBOSE)); then
	mafft --thread "$threads" --quiet \
		<(cat "${workprefix}.ref.fa" "${workprefix}.best.rotated.fa") >"${workprefix}.mafft.aln.fasta"
else
	mafft --thread "$threads" --quiet \
		<(cat "${workprefix}.ref.fa" "${workprefix}.best.rotated.fa") >"${workprefix}.mafft.aln.fasta" 2>/dev/null
fi

# ---------- Global % identity (UNGAPPED) ----------
if [[ -s "${workprefix}.mafft.aln.fasta" ]]; then
	seqkit fx2tab "${workprefix}.mafft.aln.fasta" 2>/dev/null >"$tmpdir/aln.fx2tab" || true
	nseq="$(wc -l <"$tmpdir/aln.fx2tab" | tr -d ' ')"
	if ((nseq != 2)); then
		warn "Alignment should contain exactly 2 sequences, got $nseq"
		printf "%s\n" "NA"
		exit 1
	fi
	if ((RAW)); then
		awk -f "$SCRIPTDIR/compute_pid.awk" -v RAW=1 "$tmpdir/aln.fx2tab"
	else
		awk -f "$SCRIPTDIR/compute_pid.awk" -v RAW=0 "$tmpdir/aln.fx2tab"
	fi
else
	printf "%s\n" "NA"
	exit 1
fi
