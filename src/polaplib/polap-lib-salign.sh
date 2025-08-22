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

warn() { ((VERBOSE)) && echo "[WARN] $*" >&2 || true; }

# Temp dir setup
cleanup() {
	rm -rf "$tmpdir"
}

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

need() {
	command -v "$1" >/dev/null 2>&1 || die "Missing required tool: $1"
}

get_by_id() {
	seqkit grep -nr -p "^${1}\$" "$sfa" 2>/dev/null | seqkit seq -w 0 2>/dev/null
}

# 1. find the path for the plastid genome sequence
#   1.1 find the longest contig with PT genes
#   1.2 find the path starting at the contig
# 2. stitch the genome sequence of the path from the gfa: LSC-IR-SSC-IR'
# 3. create 1 from ref and 4 from the query
gfa_to_plastome_one() {
	local gfa="$1"
	local outfa="$2"
	local outdir="$3"

	# 1.1 find the longest contig with PT genes
	_polap_lib_annotate-prepare \
		-o "${outdir}" \
		-i "pt0" \
		--gfa "${gfa}"

	_polap_lib_annotate \
		-o "${outdir}" \
		-i "pt0"

	local start_seed=$(awk 'NR==2 {print $1}' "${outdir}/pt0/pt-contig-annotation-depth-table.txt")
	_polap_log0 "reference seed: ${start_seed}"

	# 1.2 find the path starting at the contig
	_polap_log3_cmdout python "${_POLAPLIB_DIR}/polap-py-find-plastid-path-of-gfa.py" \
		--gfa "${gfa}" \
		--seed-edge "${start_seed}" \
		--out "${outdir}/pt0/plastid-path.txt"

	# 2. stitch the genome sequence of the path from the gfa: LSC-IR-SSC-IR'
	local path=$(sed -n "1p" "${outdir}/pt0/plastid-path.txt")

	if [[ -n "${path}" ]]; then
		python "${_POLAPLIB_DIR}/polap-py-stitch-path-of-gfa.py" \
			--gfa "${gfa}" \
			--path "${path}" \
			--circular-path \
			--out "${outfa}"
	fi
}

# 1. find the path for the plastid genome sequence
#   1.1 find the longest contig with PT genes
#   1.2 find the path starting at the contig
# 2. stitch the genome sequence of the path from the gfa: LSC-IR-SSC-IR'
# 3. create 1 from ref and 4 from the query
gfa_to_all_candidates() {
	local gfa="$1"
	local outfa="$2"
	local outdir="$3"

	# 1.1 find the longest contig with PT genes
	_polap_lib_annotate-prepare \
		-o "${outdir}" \
		-i "pt0" \
		--gfa "${gfa}"

	_polap_lib_annotate \
		-o "${outdir}" \
		-i "pt0"

	local start_seed=$(awk 'NR==2 {print $1}' "${outdir}/pt0/pt-contig-annotation-depth-table.txt")
	_polap_log0 "start seed: ${start_seed}"

	# 1.2 find the path starting at the contig
	python "${_POLAPLIB_DIR}/polap-py-find-plastid-path-of-gfa.py" \
		--gfa "${gfa}" \
		--seed-edge "${start_seed}" \
		--out "${outdir}/pt0/plastid-path.txt"

	if [[ -s "${outdir}/pt0/plastid-path.txt" ]]; then
		# 2. stitch the genome sequence of the path from the gfa: LSC-IR-SSC-IR'
		local i
		for i in {1..4}; do
			local path=$(sed -n "${i}p" "${outdir}/pt0/plastid-path.txt")
			if [[ -n "${path}" ]]; then
				python "${_POLAPLIB_DIR}/polap-py-stitch-path-of-gfa.py" \
					--gfa "${gfa}" \
					--path "${path}" \
					--circular-path \
					--out "${outdir}/${i}.fa"
			fi
		done

		# if there are only two
		# find LSC, IR, SSC
		if [[ ! -s "${outdir}/3.fa" ]]; then
			cp "${outdir}/1.fa" "${outdir}/../1.fa"

			_polap_log0 "  searching 1-segment ptDNA for LSC-IR-SSC-IR' using MUMmer ..."
			local IR_MIN_ID=${IR_MIN_ID:-95}
			local IR_MIN_LEN=${IR_MIN_LEN:-10000}
			_polap_log3_cmdout bash "${_POLAPLIB_DIR}/polap-bash-plastome-forms.sh" \
				-o "${outdir}" \
				--check-report "${outdir}/pt.check.report.txt" \
				"${outdir}/../1.fa"

		fi

	else
		_polap_log0 "No candidate ptDNA in GFA: $gfa"
	fi
}

find_irs_mummer() {
	local fa=$1 min_id=${2:-95} min_len=${3:-10000}
	local tmp
	tmp=$(mktemp -d)
	# self-align
	nucmer --maxmatch --nosimplify -p "$tmp/self" "$fa" "$fa" >/dev/null 2>&1
	# filter by identity and length
	delta-filter -i "$min_id" -l "$min_len" "$tmp/self.delta" >"$tmp/flt.delta"
	# tabular coords (headers in first line; -T keeps tabs)
	show-coords -rclTH "$tmp/flt.delta" >"$tmp/coords.tsv"
	# pick the longest reverse-strand self-hit (non-diagonal)
	awk -f "$SCRIPTDIR/coords_pick_ir.awk" "$tmp/coords.tsv"
	local rc=$?
	rm -rf "$tmp"
	return $rc
}

function _polap_lib_salign-pt {
	# _polap_log0 "Align two ptDNAs"

	if (($# < 2)); then
		echo "Usage: $0 <in1.gfa|fa[.gz]> <in2.gfa|fa[.gz]> [--out PREFIX|DIR/PREFIX] [--threads N] [--min-ir N] [--max-ir N] [--min-cov F] [--raw] [--verbose] [--tmp-in-out]" >&2
		exit 1
	fi

	in1="$1"
	shift
	in2="$1"
	shift

	# Defaults
	local prefix="ptalign"
	local threads="$(nproc 2>/dev/null || echo 4)"
	local min_ir=8000
	local max_ir=80000
	local min_cov="0.90"
	local VERBOSE=0
	local RAW=0
	local TMP_IN_OUT=0

	# Where AWK helpers live
	SCRIPTDIR="${_POLAPLIB_DIR}/scripts"

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

	out_percent="${prefix}.percent-identity.txt"

	warn() { ((VERBOSE)) && echo "[WARN] $*" >&2 || true; }

	# Temp dir setup
	cleanup() {
		rm -rf "$tmpdir"
	}

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

	need() {
		command -v "$1" >/dev/null 2>&1 || die "Missing required tool: $1"
	}

	#######################################################################
	# MAIN

	# _polap_log0 "in1: $in1"
	# _polap_log0 "in2: $in2"
	# _polap_log0 "Check tools ..."
	# Tools
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
	local outdir=""
	if [[ "$prefix" == */* ]]; then
		outdir="$(dirname -- "$prefix")"
		prefix="$(basename -- "$prefix")"
		mkdir -p -- "$outdir"
	fi
	# All product files use workprefix (which includes outdir if set)
	local workprefix="${outdir:+$outdir/}$prefix"

	local tmpdir
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
	local canddir="${workprefix}.cand"
	mkdir -p "$canddir"

	# ---------- Build reference & candidates ----------
	case "$(detect_type "$in1")" in
	gfa)
		gfa_to_plastome_one "$in1" "${workprefix}.ref.fa" "${outdir}"
		;;
	fasta)
		fasta_to_plastome_one "$in1" "${workprefix}.ref.fa" "${outdir}"
		;;
	*)
		die "Unknown type for $in1"
		;;
	esac

	case "$(detect_type "$in2")" in
	gfa)
		gfa_to_all_candidates "$in2" "${workprefix}.cand" "${canddir}"
		;;
	fasta)
		fasta_to_all_candidates "$in2" "${workprefix}.cand" "${canddir}"
		;;
	*)
		die "Unknown type for $in2"
		;;
	esac
	if ls -1 "$canddir"/*.fa >/dev/null 2>&1; then
		_polap_log1 "candidates produced from input2"
	else
		_polap_log0 "No candidates produced from input2"
		return
	fi

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
		-outfmt '6 qstart qend sstart send sstrand length pident qlen slen bitscore' \
		>"${workprefix}.blast.rotate.tsv" || true

	if [[ ! -s "${workprefix}.blast.rotate.tsv" ]]; then
		cp "${workprefix}.best.raw.fa" "${workprefix}.best.rotated.fa"
	else
		rot_start="$(awk -f "$SCRIPTDIR/rotation_anchor.awk" "${workprefix}.blast.rotate.tsv")"
		: "${rot_start:=1}"
		echo "${rot_start}" >"${workprefix}.best.rot_start.txt"
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

	local

	# ---------- Global % identity (UNGAPPED) ----------
	if [[ -s "${workprefix}.mafft.aln.fasta" ]]; then
		seqkit fx2tab "${workprefix}.mafft.aln.fasta" 2>/dev/null >"$tmpdir/aln.fx2tab" || true
		nseq="$(wc -l <"$tmpdir/aln.fx2tab" | tr -d ' ')"
		if ((nseq != 2)); then
			warn "Alignment should contain exactly 2 sequences, got $nseq"
			printf "%s\n" "NA" >"${out_percent}"
		fi
		if ((RAW)); then
			awk -f "$SCRIPTDIR/compute_pid.awk" -v RAW=1 "$tmpdir/aln.fx2tab" >"${out_percent}"
		else
			awk -f "$SCRIPTDIR/compute_pid.awk" -v RAW=0 "$tmpdir/aln.fx2tab" >"${out_percent}"
		fi
	else
		printf "%s\n" "NA" >"${out_percent}"
	fi

}

function _polap_lib_salign-check-pt {
	# _polap_log0 "Align two ptDNAs"

	if (($# < 2)); then
		echo "Usage: $0 <in1.gfa|fa[.gz]> [--out PREFIX|DIR/PREFIX] [--threads N] [--min-ir N] [--max-ir N] [--min-cov F] [--raw] [--verbose] [--tmp-in-out]" >&2
		exit 1
	fi

	in1="$1"
	shift

	# Defaults
	local prefix="ptalign"
	local threads="$(nproc 2>/dev/null || echo 4)"
	local min_ir=8000
	local max_ir=80000
	local min_cov="0.90"
	local VERBOSE=0
	local RAW=0
	local TMP_IN_OUT=0

	# Where AWK helpers live
	SCRIPTDIR="${_POLAPLIB_DIR}/scripts"

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

	out_percent="${prefix}.percent-identity.txt"

	#######################################################################
	# MAIN

	# _polap_log0 "in1: $in1"
	# _polap_log0 "in2: $in2"
	# _polap_log0 "Check tools ..."
	# Tools
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
	local outdir=""
	if [[ "$prefix" == */* ]]; then
		outdir="$(dirname -- "$prefix")"
		prefix="$(basename -- "$prefix")"
		mkdir -p -- "$outdir"
	fi
	# All product files use workprefix (which includes outdir if set)
	local workprefix="${outdir:+$outdir/}$prefix"

	local tmpdir
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
	local canddir="${workprefix}.cand"
	mkdir -p "$canddir"

	# ---------- Build reference & candidates ----------
	case "$(detect_type "$in1")" in
	gfa)
		gfa_to_plastome_one "$in1" "${workprefix}.ref.fa" "${outdir}"
		;;
	fasta)
		fasta_to_plastome_one "$in1" "${workprefix}.ref.fa" "${outdir}"
		;;
	*)
		die "Unknown type for $in1"
		;;
	esac
}
