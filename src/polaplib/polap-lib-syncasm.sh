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
# polap-lib-syncasm.sh  (v0.0.1)
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

# ──────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ──────────────────────────────────────────────────────────────────────────────

# check required executables are in PATH
_polap__need_exec() {
	local miss=0 x
	for x in "$@"; do
		if ! command -v "$x" >/dev/null 2>&1; then
			_polap_log0 "[ERROR] missing executable: $x"
			miss=1
		fi
	done
	return "$miss"
}

# expand "A", "A,B", "A,B,C", or "A,B,C,D,..." into a whitespace list (int or float)
_polap__expand_series() {
	local spec="$1" is_float="${2:-0}"
	spec="${spec//[[:space:]]/}"
	IFS=',' read -r -a toks <<<"$spec"
	local n="${#toks[@]}"
	((n == 0)) && {
		echo ""
		return 1
	}

	_seq_int() {
		local A="$1" B="$2" cnt="$3"
		local out=() step num i
		if ((cnt > 0)); then
			((cnt < 2)) && cnt=2
			if ((A <= B)); then
				step=$(awk -v a="$A" -v b="$B" -v c="$cnt" 'BEGIN{printf "%.10f",(b-a)/(c-1)}')
				for ((i = 0; i < cnt; i++)); do
					num=$(awk -v a="$A" -v s="$step" -v i="$i" 'BEGIN{printf "%.0f",a+i*s}')
					out+=("$num")
				done
			else
				step=$(awk -v a="$A" -v b="$B" -v c="$cnt" 'BEGIN{printf "%.10f",(a-b)/(c-1)}')
				for ((i = 0; i < cnt; i++)); do
					num=$(awk -v a="$A" -v s="$step" -v i="$i" 'BEGIN{printf "%.0f",a-i*s}')
					out+=("$num")
				done
			fi
			echo "${out[@]}"
			return 0
		else
			_polap__expand_series "${A},${B},5" 0
			return $?
		fi
	}

	_seq_float() {
		local A="$1" B="$2" cnt="$3" out=() i val
		((cnt < 2)) && cnt=2
		if awk -v a="$A" -v b="$B" 'BEGIN{exit (a<=b)?0:1}'; then
			for ((i = 0; i < cnt; i++)); do
				val=$(
					python - "$A" "$B" "$cnt" "$i" <<'PY'
import sys
A=float(sys.argv[1]); B=float(sys.argv[2]); C=int(sys.argv[3]); i=int(sys.argv[4])
if C<2: C=2
print(f"{A + (B-A)*i/(C-1):.3f}")
PY
				)
				out+=("$val")
			done
		else
			for ((i = 0; i < cnt; i++)); do
				val=$(
					python - "$A" "$B" "$cnt" "$i" <<'PY'
import sys
A=float(sys.argv[1]); B=float(sys.argv[2]); C=int(sys.argv[3]); i=int(sys.argv[4])
if C<2: C=2
print(f"{A - (A-B)*i/(C-1):.3f}")
PY
				)
				out+=("$val")
			done
		fi
		echo "${out[@]}"
	}

	if ((n == 1)); then
		echo "${toks[0]}"
	elif ((n == 2)); then
		((is_float)) && _seq_float "${toks[0]}" "${toks[1]}" 5 || _seq_int "${toks[0]}" "${toks[1]}" 0
	elif ((n == 3)); then
		((is_float)) && _seq_float "${toks[0]}" "${toks[1]}" "${toks[2]}" || _seq_int "${toks[0]}" "${toks[1]}" "${toks[2]}"
	else
		echo "${toks[@]}"
	fi
}

# v0.0.2  — generate k values on the “...1” ladder with safe set -u handling
_polap__expand_series_k1() {
	local spec="$1"
	spec="${spec//[[:space:]]/}"
	local A B C
	IFS=',' read -r A B C <<<"$spec"

	# default count when only A,B given
	local want="${C:-5}"

	# Single token: snap to nearest ...1 (tie -> down), return one value
	if [[ -z "$B" ]]; then
		[[ "$A" =~ ^[0-9]+$ ]] || {
			echo ""
			return 1
		}
		if ((A % 10 == 1)); then
			echo "$A"
		else
			local a1=$((A - ((A - 1) % 10)))
			echo "$a1"
		fi
		return 0
	fi

	# Two/three tokens
	[[ "$A" =~ ^[0-9]+$ && "$B" =~ ^[0-9]+$ ]] || {
		echo ""
		return 1
	}
	[[ "$want" =~ ^[0-9]+$ ]] || {
		echo ""
		return 1
	}
	((want >= 2)) || want=2

	local up=1 lo hi step
	if ((A <= B)); then
		up=1
		lo=$(((A % 10 == 1) ? A : A + ((11 - (A % 10)) % 10)))
		hi=$(((B % 10 == 1) ? B : B - ((B - 1) % 10)))
		step=10
		((lo > hi)) && {
			echo ""
			return 1
		}
	else
		up=0
		lo=$(((B % 10 == 1) ? B : B + ((11 - (B % 10)) % 10)))
		hi=$(((A % 10 == 1) ? A : A - ((A - 1) % 10)))
		step=-10
		((lo > hi)) && {
			echo ""
			return 1
		}
	fi

	# Build grid on the “…1” ladder
	local grid=() x
	if ((up)); then
		for ((x = lo; x <= hi; x += 10)); do grid+=("$x"); done
	else
		for ((x = hi; x >= lo; x -= 10)); do grid+=("$x"); done
	fi
	local n=${#grid[@]}
	((n > 0)) || {
		echo ""
		return 1
	}

	# If grid is short, return all
	if ((want >= n)); then
		echo "${grid[@]}"
		return 0
	fi

	# Evenly subsample 'want' points including endpoints
	local out_idxs=() i idx
	for ((i = 0; i < want; i++)); do
		idx=$(((i * (n - 1) + (want - 1)) / (want - 1)))
		((idx < 0)) && idx=0
		((idx >= n)) && idx=$((n - 1))
		# Avoid duplicates from integer rounding
		if ((i > 0 && idx == out_idxs[i - 1])); then
			idx=$((idx + (up ? 1 : -1)))
			((idx < 0)) && idx=0
			((idx >= n)) && idx=$((n - 1))
		fi
		out_idxs+=("$idx")
	done

	# Emit unique values; set -u-safe existence test
	declare -A seen=()
	local uniq=() key val
	for idx in "${out_idxs[@]}"; do
		val="${grid[idx]}"
		key="k$val"
		if [[ -v 'seen[$key]' ]]; then
			continue
		fi
		seen[$key]=1
		uniq+=("$val")
	done

	echo "${uniq[@]}"
}

# If you still want the generic dispatcher, add a k-mode shim:
#   _polap__expand_series_k "251,101,5" → calls _polap__expand_series_k1
_polap__expand_series_k() {
	_polap__expand_series_k1 "$1"
}

# validate constraints and normalize values; emit only valid members
_polap__validate_list() {
	local type="$1"
	shift
	local ok=() x
	for x in "$@"; do
		case "$type" in
		c) [[ "$x" =~ ^[0-9]+$ ]] && ((x >= 1 && x <= 1000)) && ok+=("$x") || (: && :) ;;
		k) [[ "$x" =~ ^[0-9]+$ ]] && ((x >= 51 && x <= 8001)) && [[ "$x" =~ 1$ ]] && ok+=("$x") || (: && :) ;;
		s) [[ "$x" =~ ^[0-9]+$ ]] && ((x >= 11 && x <= 31 && x % 2 == 1)) && ok+=("$x") || (: && :) ;;
		a) if [[ "$x" =~ ^0\.[0-9]+$ ]] && awk -v v="$x" 'BEGIN{exit (v>0&&v<1)?0:1}'; then
			ok+=("$(awk -v v="$x" 'BEGIN{printf "%.3f", v+0}')")
		fi ;;
		esac
	done
	echo "${ok[@]}"
}

# GFA helpers
_polap__gfa_total_bp() {
	local gfa="$1"
	# awk: sum S segment lengths (sequence or LN tag)
	awk '$1=="S"{ if($3!="*") len+=length($3); else for(i=4;i<=NF;i++) if($i~/^LN:i:/){split($i,a,":");len+=a[3]} } END{print len+0}' "$gfa"
}

_polap__gfa_stats() {
	local gfa="$1"
	[[ -s "$gfa" ]] || {
		echo "0 0"
		return 0
	}
	# awk: count S and sum lengths
	awk 'BEGIN{n=0;L=0} $1=="S"{n++; if($3!="*") L+=length($3); else for(i=4;i<=NF;i++) if($i~/^LN:i:/){split($i,a,":");L+=a[3]+0;break}} END{print L,n}' "$gfa"
}

# parse exp_* name into an associative hash (flexible digits)
_polap__parse_run_name() {
	local name="$1"
	declare -gA SYNCASM_PARSE
	SYNCASM_PARSE=()
	[[ "$name" =~ _k([0-9]+) ]] && SYNCASM_PARSE[k]="${BASH_REMATCH[1]}"
	[[ "$name" =~ _s([0-9]+) ]] && SYNCASM_PARSE[s]="${BASH_REMATCH[1]}"
	[[ "$name" =~ _c([0-9]+) ]] && SYNCASM_PARSE[c]="${BASH_REMATCH[1]}"
	[[ "$name" =~ _u([0-9]+) ]] && SYNCASM_PARSE[u]="${BASH_REMATCH[1]}"
	[[ "$name" =~ _a([0-9]+) ]] && SYNCASM_PARSE[a]="0.${BASH_REMATCH[1]}"
	[[ "$name" =~ _wx([0-9]+) ]] && SYNCASM_PARSE[wx]="0.${BASH_REMATCH[1]}"
}

# # float normalization & tags
# _polap__fmt_float3() {
# 	local v="$1"
# 	[[ "$v" =~ ^0\.[0-9]+$ ]] && awk -v x="$v" 'BEGIN{exit (x>0&&x<1)?0:1}' && printf '0.%03d\n' "$((10#${v#0.}))" || {
# 		echo ""
# 		return 1
# 	}
# }
#
# v0.0.2
_polap__fmt_float3() {
	local v="$1"
	# accept 0.x, .x, 0.xx, 0.xxx (no leading/trailing spaces)
	[[ "$v" =~ ^0?\.[0-9]+$ ]] || {
		echo ""
		return 1
	}
	awk -v x="$v" 'BEGIN{
    if (x>0 && x<1) { printf "%.3f\n", x; exit 0 }
    exit 1
  }' || {
		echo ""
		return 1
	}

}
_polap__tag_from_float3() {
	local v="${1#0.}"
	v="${v//[^0-9]/}"
	printf '%03d' "$((10#$v))"
}
_polap__pad2() { printf '%02d' "$((10#$1))"; }
_polap__pad3() { printf '%03d' "$((10#$1))"; }
_polap__pad4() { printf '%04d' "$((10#$1))"; }

normalize_name() {
	local raw="$1" c a k s wx u
	[[ "$raw" =~ _c([0-9]+) ]] && c="$(_polap__pad3 "${BASH_REMATCH[1]}")"
	[[ "$raw" =~ _a([0-9]+) ]] && a="$(_polap__pad3 "${BASH_REMATCH[1]}")"
	[[ "$raw" =~ _k([0-9]+) ]] && k="$(_polap__pad4 "${BASH_REMATCH[1]}")"
	[[ "$raw" =~ _s([0-9]+) ]] && s="$(_polap__pad2 "${BASH_REMATCH[1]}")"
	[[ "$raw" =~ _wx([0-9]+) ]] && wx="$(_polap__pad3 "${BASH_REMATCH[1]}")"
	[[ "$raw" =~ _u([0-9]+) ]] && u="$(_polap__pad2 "${BASH_REMATCH[1]}")"
	echo "exp_c${c}_a${a}_k${k}_s${s}_wx${wx}_u${u}"
}

# _polap_lib_subset_raw_by_hpc_mapping v0.2.0
# Map HPC reads to HPC contigs, then subset RAW reads by mapped IDs.
# Adds --min-aln-len (PAF col11 alignment block length threshold).
_polap_lib_syncasm-subset_raw_by_hpc_mapping() {
	local _ver="v0.2.0"
	local _verbose=0 _quiet=0

	# inputs
	local raw=""       # -r/--raw RAW reads (fastq[.gz])  [required]
	local contigs=""   # -c/--contigs HPC contigs (fasta) [required]
	local hpc_reads="" # --hpc-reads precomputed HPC reads (optional)

	# knobs
	local outprefix="subset_out" # -o/--out
	local threads=16             # -t/--threads
	local preset="map-ont"       # -x/--preset
	local kmer=14                # -k/--kmer
	local min_id="0.80"          # --min-id (matches/alnlen)
	local min_qcov="0.80"        # --min-qcov ((qend-qstart)/qlen)
	local min_aln_len=3000       # --min-aln-len (PAF col11); 0 = off

	# parse CLI
	while [[ $# -gt 0 ]]; do
		case "$1" in
		--version)
			echo "_polap_lib_subset_raw_by_hpc_mapping ${_ver}"
			return 0
			;;
		-v | --verbose)
			_verbose=$((_verbose + 1))
			shift
			;;
		--quiet)
			_quiet=1
			_verbose=0
			shift
			;;
		-r | --raw)
			raw="$2"
			shift 2
			;;
		-c | --contigs)
			contigs="$2"
			shift 2
			;;
		--hpc-reads)
			hpc_reads="$2"
			shift 2
			;;
		-o | --out)
			outprefix="$2"
			shift 2
			;;
		-t | --threads)
			threads="$2"
			shift 2
			;;
		-x | --preset)
			preset="$2"
			shift 2
			;;
		-k | --kmer)
			kmer="$2"
			shift 2
			;;
		--min-id)
			min_id="$2"
			shift 2
			;;
		--min-qcov)
			min_qcov="$2"
			shift 2
			;;
		--min-aln-len)
			min_aln_len="$2"
			shift 2
			;;
		-h | --help)
			cat <<EOF
_polap_lib_subset_raw_by_hpc_mapping ${_ver}
Map HPC reads -> HPC contigs, then subset RAW reads by mapped IDs.

Required:
  -r, --raw FILE       RAW reads (FASTQ[.gz])
  -c, --contigs FILE   HPC contigs (FASTA)

Optional:
      --hpc-reads FILE  HPC reads (if omitted, generate via seqtk hpc)
  -o, --out PREFIX      Output prefix [subset_out]
  -t, --threads INT     Threads [16]
  -x, --preset PRESET   minimap2 preset [map-ont]
  -k, --kmer INT        minimap2 -k [14]
      --min-id FLOAT    identity cutoff (matches/alnlen) [0.80]
      --min-qcov FLOAT  read coverage cutoff (qcov=(qend-qstart)/qlen) [0.50]
      --min-aln-len INT minimum alignment block length (PAF col11) [0 = off]
  -v, --verbose         Verbose logs
      --quiet           Suppress info logs
      --version         Print version

Outputs (using PREFIX):
  PREFIX.paf
  PREFIX.mapped.ids
  PREFIX.raw.mapped.fastq.gz
  PREFIX.raw.unmapped.fastq.gz
EOF
			return 0
			;;
		--)
			shift
			break
			;;
		*)
			_polap_log0 "[ERROR] unknown option: $1"
			return 2
			;;
		esac
	done

	# info logger (only when verbose and not quiet)
	_info() { [[ $_verbose -gt 0 && $_quiet -eq 0 ]] && _polap_log1 "$*"; }

	# sanity
	[[ -s "$raw" ]] || {
		_polap_log0 "[ERROR] missing --raw $raw"
		return 2
	}
	[[ -s "$contigs" ]] || {
		_polap_log0 "[ERROR] missing --contigs $contigs"
		return 2
	}
	local exe
	for exe in minimap2 seqtk awk sort comm; do
		command -v "$exe" >/dev/null 2>&1 || {
			_polap_log0 "[ERROR] missing executable: $exe"
			return 2
		}
	done

	# paths
	local paf="${outprefix}.paf"
	local mapped_ids="${outprefix}.mapped.ids"
	local raw_mapped="${outprefix}.raw.mapped.fastq.gz"
	local raw_unmapped="${outprefix}.raw.unmapped.fastq.gz"
	local tmp_all_ids="${outprefix}.all.ids"
	local hpc="$hpc_reads"

	# HPC reads
	if [[ -n "$hpc" && -s "$hpc" ]]; then
		_info "[HPC] using provided HPC reads: $hpc"
	else
		hpc="${outprefix}.reads.hpc.fastq"
		_info "[HPC] generating HPC reads → $hpc"
		seqtk hpc "$raw" >"$hpc" || {
			_polap_log0 "[ERROR] seqtk hpc failed"
			return 2
		}
	fi

	# map
	_info "[map] minimap2 -x ${preset} -k ${kmer} (HPC reads → HPC contigs)"
	minimap2 -t "$threads" -x "$preset" -k "$kmer" --secondary=yes -c \
		"$contigs" "$hpc" >"$paf" || {
		_polap_log0 "[ERROR] minimap2 failed"
		return 2
	}

	# filter PAF → IDs (identity, qcov, alnlen)
	_info "[map] filtering PAF: min_id=${min_id}, min_qcov=${min_qcov}, min_aln_len=${min_aln_len}"
	awk -v ID="$min_id" -v QC="$min_qcov" -v ML="$min_aln_len" '
    # PAF: 1=qname 2=qlen 3=qstart 4=qend 6=tname 8=tstart 9=tend 10=nmatch 11=alnlen
    ($10>0 && $11>0 && $2>0){
      id   = $10/$11
      qcov = ($4-$3)/$2
      ok   = (id>=ID && qcov>=QC)
      if (ML>0) ok = ok && ($11>=ML)
      if (ok) print $1
    }' "$paf" | sort -u >"$mapped_ids"

	local nmap
	nmap=$(wc -l <"$mapped_ids" 2>/dev/null || echo 0)
	_info "[map] mapped HPC read IDs: $nmap"

	# subset RAW reads by IDs
	_info "[raw] extracting mapped RAW reads →  $raw_mapped"
	seqtk subseq "$raw" "$mapped_ids" | gzip -c >"$raw_mapped" ||
		{
			_polap_log0 "[ERROR] seqtk subseq (mapped) failed"
			return 2
		}

	_info "[raw] computing unmapped RAW reads →  $raw_unmapped"
	if [[ "$raw" =~ \.gz$ ]]; then
		zcat -- "$raw" | awk 'NR%4==1{print substr($1,2)}' | sort -u >"$tmp_all_ids"
	else
		awk 'NR%4==1{print substr($1,2)}' "$raw" | sort -u >"$tmp_all_ids"
	fi
	comm -23 "$tmp_all_ids" "$mapped_ids" >"${outprefix}.unmapped.ids"
	seqtk subseq "$raw" "${outprefix}.unmapped.ids" | gzip -c >"$raw_unmapped" ||
		{
			_polap_log0 "[ERROR] seqtk subseq (unmapped) failed"
			return 2
		}

	local n_all n_unmap
	n_all=$(wc -l <"$tmp_all_ids" 2>/dev/null || echo 0)
	n_unmap=$(wc -l <"${outprefix}.unmapped.ids" 2>/dev/null || echo 0)
	_info "[done] total RAW: $n_all | mapped: $nmap | unmapped: $n_unmap"
	_info "[out] $paf"
	_info "[out] $mapped_ids"
	_info "[out] $raw_mapped"
	_info "[out] $raw_unmapped"
	return 0
}

# _polap_lib_pt_iterate v0.1.0
# Run the ptDNA iterative assembly/annotation loop (your Stage-7 block) with a clean CLI.
_polap_lib_syncsam-pt_iterate() {
	local _ver="v0.1.0"
	local _verbose=0 _quiet=0

	# user options
	local annotatedir=""                      # -o/--outdir
	local type="annotate-read"                # -t/--type  (used in final symlinks)
	local long_reads="${_arg_long_reads}"     # -l/--long-reads (input to _polap_lib_assemble-rate)
	local single_min="${_arg_single_min}"     # --single-min  (pass-through to _polap_lib_assemble-rate -w)
	local data_type="${_arg_data_type}"       # --data-type   ("pacbio-raw" triggers omega fallback)
	local resolved_fastq="${_arg_long_reads}" # --resolved-fastq (for omega fallback)
	local outdir_final=""                     # --final-outdir (target for final symlinks); default: parent of annotatedir

	# parse CLI
	while [[ $# -gt 0 ]]; do
		case "$1" in
		--version)
			echo "_polap_lib_pt_iterate ${_ver}"
			return 0
			;;
		-v | --verbose)
			_verbose=$((_verbose + 1))
			shift
			;;
		--quiet)
			_quiet=1
			_verbose=0
			shift
			;;
		-o | --outdir)
			annotatedir="$2"
			shift 2
			;;
		-t | --type)
			type="$2"
			shift 2
			;;
		-l | --long-reads)
			long_reads="$2"
			shift 2
			;;
		--single-min)
			single_min="$2"
			shift 2
			;;
		--data-type)
			data_type="$2"
			shift 2
			;;
		--resolved-fastq)
			resolved_fastq="$2"
			shift 2
			;;
		--final-outdir)
			outdir_final="$2"
			shift 2
			;;
		-h | --help)
			cat <<EOF
_polap_lib_pt_iterate ${_ver}
Run the ptDNA iterative seed/annotate loop and create stage artifacts & symlinks.

Required:
  -o, --outdir DIR            Annotate/working directory (annotatedir)
  -l, --long-reads FILE       Long reads (FASTQ[.gz]) for assemble-rate

Optional:
  -t, --type STR              Label for final symlinks [annotate-read]
      --single-min INT        Pass-through to assemble-rate (-w)
      --data-type STR         Data type (e.g., pacbio-raw triggers omega fallback)
      --resolved-fastq FILE   Reads for omega fallback (pacbio-raw only)
      --final-outdir DIR      Final symlink root [parent of annotatedir]
  -v, --verbose               Verbose logs
      --quiet                 Suppress info logs
      --version               Print version
EOF
			return 0
			;;
		--)
			shift
			break
			;;
		*)
			_polap_log0 "[ERROR] unknown option: $1"
			return 2
			;;
		esac
	done

	# logger
	_info() { [[ $_verbose -gt 0 && $_quiet -eq 0 ]] && _polap_log1 "$*"; }

	# sanity
	[[ -n "$annotatedir" ]] || {
		_polap_log0 "[ERROR] --outdir is required"
		return 2
	}
	[[ -n "$long_reads" ]] || {
		_polap_log0 "[ERROR] --long-reads is required"
		return 2
	}
	[[ -d "$annotatedir" ]] || mkdir -p "$annotatedir"

	# default final outdir = parent of annotatedir if unset
	if [[ -z "$outdir_final" ]]; then
		outdir_final="$(cd "$annotatedir/.." && pwd)"
	fi

	# --- Stage 0 inputs & visualization ---
	local ptdna0_gfa="${annotatedir}/pt/30-contigger/graph_final.gfa"
	if [[ ! -s "$ptdna0_gfa" ]]; then
		_polap_log0 "No ptDNA assembly stage 0: $ptdna0_gfa"
		return 0
	fi

	# Extract DNA, plot with bandage, and stage-0 symlinks
	_polap_lib_pt-extract-dna \
		"$ptdna0_gfa" \
		"${annotatedir}/pt/ptdna"

	_polap_lib_bandage \
		"$ptdna0_gfa" \
		"${annotatedir}/pt/30-contigger/graph_final.png"

	ln -sf "pt/ptdna/pt.0.fa" "${annotatedir}/pt.0.fa"
	ln -sf "pt/30-contigger/graph_final.gfa" "${annotatedir}/pt.0.gfa"
	ln -sf "pt/30-contigger/graph_final.png" "${annotatedir}/pt.0.png"

	# alias pt -> pt0 for the loop
	ln -sfn pt "${annotatedir}/pt0"

	# --- Iterative loop pt0 -> pt1 -> ... -> pt7 (max) ---
	local i j
	for ((i = 4; i < 8; i++)); do
		j=$((i + 1))

		# 1) annotate (for seeding); select connected PT components
		_info "[pt$i] annotate"
		_polap_lib_annotate -o "$annotatedir" -i "pt$i"

		# 2) seed plastid: produce pt$i/mt.contig.name-pt$j
		_info "[pt$i] seed plastid -> pt$j"
		_polap_lib_seed-plastid -o "$annotatedir" -i "pt$i" -j "pt$j"

		if [[ ! -s "${annotatedir}/pt$i/mt.contig.name-pt$j" ]]; then
			_polap_log0 "No ptDNA seed for pt$j"
			if [[ "$data_type" == "pacbio-raw" ]]; then
				_info "pacbio-raw: use input long reads with adjusted omega for final mtDNA assembly"
				_polap_lib_assemble-omega \
					-o "$annotatedir" \
					-l "$resolved_fastq" \
					-t pt \
					-i pt0 -j ptx
			fi
			return 0
		fi

		# 3) assemble-rate: build next stage pt$j
		_info "[pt$i] assemble-rate -> pt$j"
		_polap_lib_assemble-rate \
			-o "$annotatedir" \
			-l "$long_reads" \
			-w "$single_min" \
			-i "pt$i" -j "pt$j"

		# 4) stage pt$j artifacts
		ptdna0_gfa="${annotatedir}/pt$j/30-contigger/graph_final.gfa"
		if [[ -s "$ptdna0_gfa" ]]; then
			_polap_lib_pt-extract-dna "$ptdna0_gfa" "${annotatedir}/pt$j/ptdna"

			# NOTE: using assembly_graph.gfa/png for pt$j visualization (as in your code)
			_polap_lib_bandage \
				"${annotatedir}/pt$j/assembly_graph.gfa" \
				"${annotatedir}/pt$j/assembly_graph.png"

			ln -sf "pt$j/ptdna/pt.0.fa" "${annotatedir}/pt.$j.fa"
			ln -sf "pt$j/assembly_graph.gfa" "${annotatedir}/pt.$j.gfa"
			ln -sf "pt$j/assembly_graph.png" "${annotatedir}/pt.$j.png"

			if [[ -s "${annotatedir}/pt.$j.fa" ]]; then
				_info "ptDNA assembly reached stage ${j}: ${annotatedir}/pt.$j.gfa"
				break
			fi
		else
			_polap_log0 "No ptDNA assembly stage $j"
		fi
	done

	# --- Post-loop: annotate the final stage (unless loop exhausted) ---
	j=$((i + 1))
	if [[ "$i" == "7" ]]; then
		_info "No ptDNA assembly stage $j"
	else
		_polap_lib_annotate -o "$annotatedir" -i "pt$j"
	fi

	# column print (if present)
	_polap_log0_column "${annotatedir}/pt$j/pt-contig-annotation-depth-table.txt"

	# --- Final symlinks (type-scoped) ---
	mkdir -p "${outdir_final}"
	ln -sf "annotate-read-${type}/pt.0.fa" "${outdir_final}/${type}-pt.0.fa"
	ln -sf "annotate-read-${type}/pt.0.gfa" "${outdir_final}/${type}-pt.0.gfa"
	ln -sf "annotate-read-${type}/pt.0.png" "${outdir_final}/${type}-pt.0.png"

	ln -sf "annotate-read-${type}/pt.$j.fa" "${outdir_final}/${type}-pt.1.fa"
	ln -sf "annotate-read-${type}/pt.$j.gfa" "${outdir_final}/${type}-pt.1.gfa"
	ln -sf "annotate-read-${type}/pt.$j.png" "${outdir_final}/${type}-pt.1.png"

	return 0
}

# _polap_lib_mt_iterate v0.1.0
# Wraps the mtDNA iterative assemble/annotate loop and exports staged symlinks.
_polap_lib_syncsam-mt_iterate() {
	local _ver="v0.1.0"
	local _verbose=0 _quiet=0

	# -------- CLI --------
	local annotatedir=""                # -o/--outdir   (required)
	local long_reads=""                 # -l/--long-reads (required)
	local data_type="${_arg_data_type}" # --data-type   (pacbio-hifi | nano-raw | ...)
	local single_min="3000"             # --single-min  (pass to assemble-rate -w)
	local pt_ref=""                     # --pt-ref      (pt.0.gfa; if set & pacbio-hifi -> filter reads)
	local final_outdir=""               # --final-outdir (default: parent of annotatedir)

	while [[ $# -gt 0 ]]; do
		case "$1" in
		--version)
			echo "_polap_lib_mt_iterate ${_ver}"
			return 0
			;;
		-v | --verbose)
			_verbose=$((_verbose + 1))
			shift
			;;
		--quiet)
			_quiet=1
			_verbose=0
			shift
			;;
		-o | --outdir)
			annotatedir="$2"
			shift 2
			;;
		-l | --long-reads)
			long_reads="$2"
			shift 2
			;;
		--data-type)
			data_type="$2"
			shift 2
			;;
		--single-min)
			single_min="$2"
			shift 2
			;;
		--pt-ref)
			pt_ref="$2"
			shift 2
			;;
		--final-outdir)
			final_outdir="$2"
			shift 2
			;;
		-h | --help)
			cat <<EOF
_polap_lib_mt_iterate ${_ver}
Run mitochondrial iterative assemble/annotate loop and export staged symlinks.

Required:
  -o, --outdir DIR         Annotate/working directory (annotatedir)
  -l, --long-reads FILE    Long reads (FASTQ[.gz]) for mt assembly

Optional:
      --data-type STR      e.g., pacbio-hifi | nano-raw (controls assembler)
      --single-min INT     Pass-through to assemble-rate (-w)
      --pt-ref FILE        ptDNA GFA used to filter HiFi reads for mt (filter-reads-by-reference)
      --final-outdir DIR   Root for final symlinks [parent of --outdir]
  -v, --verbose            Verbose logs
      --quiet              Suppress info logs
      --version            Print version
EOF
			return 0
			;;
		--)
			shift
			break
			;;
		*)
			_polap_log0 "[ERROR] unknown option: $1"
			return 2
			;;
		esac
	done

	# logger
	_info() { [[ $_verbose -gt 0 && $_quiet -eq 0 ]] && _polap_log1 "$*"; }

	# -------- sanity --------
	[[ -n "$annotatedir" ]] || {
		_polap_log0 "[ERROR] --outdir is required"
		return 2
	}
	[[ -n "$long_reads" ]] || {
		_polap_log0 "[ERROR] --long-reads is required"
		return 2
	}
	[[ -d "$annotatedir" ]] || mkdir -p "$annotatedir"
	if [[ -z "$final_outdir" ]]; then
		final_outdir="$(cd "$annotatedir/.." && pwd)"
	fi

	# -------- stage 0: extract/plot mt0 if present --------
	local mt0_gfa="${annotatedir}/mt/assembly_graph.gfa"
	if [[ -s "$mt0_gfa" ]]; then
		_polap_lib_mt-extract-dna "$mt0_gfa" "${annotatedir}/mt/mtdna"
		_polap_lib_bandage "$mt0_gfa" "${annotatedir}/mt/assembly_graph.png"

		ln -sf "mt/mtdna/mt.0.fa" "${annotatedir}/mt.0.fa"
		ln -sf "mt/assembly_graph.gfa" "${annotatedir}/mt.0.gfa"
		ln -sf "mt/assembly_graph.png" "${annotatedir}/mt.0.png"
	else
		_polap_log0 "No mtDNA assembly stage 0: $mt0_gfa"
	fi

	# alias mt -> mt0
	ln -sfn mt "${annotatedir}/mt0"

	# -------- optional: HiFi read filtering using pt reference --------
	local resolved_fastq="$long_reads"
	if [[ "$data_type" == "pacbio-hifi" && -n "$pt_ref" ]]; then
		if [[ -s "$pt_ref" ]]; then
			_info "filtering HiFi reads by ptDNA reference: $pt_ref"
			_polap_lib_filter-reads-by-reference \
				-o "$annotatedir" \
				-l "$long_reads" \
				--reference "$pt_ref"
			resolved_fastq="${annotatedir}/kmer/ref-filtered.fastq"
		else
			_polap_log0 "pt reference missing for HiFi filtering: $pt_ref"
		fi
	fi

	# -------- iterative loop mt0..mt6 -> mt1..mt7 --------
	local i j
	for ((i = 0; i < 6; i++)); do
		j=$((i + 1))

		# annotate (for seeding; select connected mt components)
		_info "[mt$i] annotate"
		_polap_lib_annotate -o "$annotatedir" -i "mt$i"

		# seed mito
		_info "[mt$i] seed-mito -> mt$j"
		_polap_lib_seed-mito -o "$annotatedir" -i "mt$i" -j "mt$j"

		if [[ ! -s "${annotatedir}/mt$i/mt.contig.name-mt$j" ]]; then
			_polap_log0 "No mt seed for mt$j"
			break
		fi

		# assemble next stage
		if [[ "$data_type" == "pacbio-hifi" ]]; then
			_info "[mt$i] assemble-rate (HiFi filtered) -> mt$j"
			_polap_lib_assemble-rate \
				-o "$annotatedir" \
				-l "$resolved_fastq" \
				-w "$single_min" \
				-t mt \
				-i "mt$i" -j "mt$j"
		elif [[ "$data_type" == "nano-raw" ]]; then
			_info "[mt$i] assemble-rate (nano selected) -> mt$j"
			_polap_lib_assemble-rate \
				-o "$annotatedir" \
				-l "${resolved_fastq}" \
				-w "$single_min" \
				-t mt -i "mt$i" -j "mt$j"
		else
			# default path if data_type not recognized: try assemble-rate on provided long_reads
			_info "[mt$i] assemble-rate (default long reads) -> mt$j"
			_polap_lib_assemble-rate \
				-o "$annotatedir" \
				-l "$long_reads" \
				-w "$single_min" \
				-t mt -i "mt$i" -j "mt$j"
		fi

		# stage artifacts
		local mtj_gfa="${annotatedir}/mt$j/assembly_graph.gfa"
		if [[ -s "$mtj_gfa" ]]; then
			_polap_lib_mt-extract-dna "$mtj_gfa" "${annotatedir}/mt$j/mtdna"
			_polap_lib_bandage "$mtj_gfa" "${annotatedir}/mt$j/assembly_graph.png"

			ln -sf "mt$j/mtdna/mt.0.fa" "${annotatedir}/mt.$j.fa"
			ln -sf "mt$j/assembly_graph.gfa" "${annotatedir}/mt.$j.gfa"
			ln -sf "mt$j/assembly_graph.png" "${annotatedir}/mt.$j.png"

			_polap_log0 "mtDNA assembly: ${annotatedir}/mt.$j.gfa"
		else
			_polap_log0 "No mt assembly $j"
		fi
	done

	# -------- extra pass: annotate->seed and final assemble --------
	j=$((i + 1))
	_polap_lib_annotate -o "$annotatedir" -i "mt$i"
	_polap_lib_seed-mito -o "$annotatedir" -i "mt$i" -j "mt$j"

	if [[ "$data_type" == "pacbio-hifi" ]]; then
		_info "final assemble-rate with HiFi filtered reads"
		_polap_lib_assemble-rate \
			-o "$annotatedir" -l "$resolved_fastq" -w "$single_min" -i "mt$i" -j "mt$j"
	elif [[ "$data_type" == "nano-raw" ]]; then
		_info "final assemble-omega with nano RAW reads"
		_polap_lib_assemble-omega \
			-o "$annotatedir" -l "$long_reads" -i "mt$i" -j "mt$j"
	else
		_info "final assemble-rate (default long reads)"
		_polap_lib_assemble-rate \
			-o "$annotatedir" -l "$long_reads" -w "$single_min" -i "mt$i" -j "mt$j"
	fi

	_polap_lib_annotate -o "$annotatedir" -i "mt$j"
	_polap_log0_column "${annotatedir}/mt$j/contig-annotation-depth-table.txt"

	# plot & link final stage if present
	if [[ -s "${annotatedir}/mt$j/assembly_graph.gfa" ]]; then
		_polap_lib_bandage \
			"${annotatedir}/mt$j/assembly_graph.gfa" \
			"${annotatedir}/mt$j/assembly_graph.png"
		ln -sf "mt$j/assembly_graph.gfa" "${annotatedir}/mt.$j.gfa"
		ln -sf "mt$j/assembly_graph.png" "${annotatedir}/mt.$j.png"
		_polap_log0 "mtDNA assembly: ${annotatedir}/mt.$j.gfa"
	else
		_polap_log0 "No MT assembly $j"
	fi

	# -------- export final symlinks --------
	i=$((j - 1))
	mkdir -p "$final_outdir"
	ln -sf "annotate-read-mt/mt.$i.gfa" "${final_outdir}/mt.0.gfa"
	ln -sf "annotate-read-mt/mt.$i.png" "${final_outdir}/mt.0.png"
	ln -sf "annotate-read-mt/mt.$j.gfa" "${final_outdir}/mt.1.gfa"
	ln -sf "annotate-read-mt/mt.$j.png" "${final_outdir}/mt.1.png"

	return 0
}

# ──────────────────────────────────────────────────────────────────────────────
# _polap_lib_syncasm (monolith, oatk + GNU parallel, c→a→k→s)
# ──────────────────────────────────────────────────────────────────────────────
_polap_lib_syncasm() {
	# CLI defaults (renamed to _crg_)
	local _crg_verbose=0
	local _crg_version="off"
	local _crg_steps="1-7" # default; use --steps all for 1-16
	local _crg_long_reads=""
	local _crg_out_prefix="syncasm.asm"
	local _crg_threads=4 # GNU parallel -j
	local _crg_k="121"
	local _crg_s="23"
	local _crg_c="30"
	local _crg_a="0.030"
	local _crg_max_bubble=100000
	local _crg_max_tip=10000
	local _crg_weak_cross="0.010"
	local _crg_unzip_round=1
	local _crg_no_read_ec="false"
	local _crg_no_hpc="false"
	local _crg_oatkdb="${OATKDB:-$HOME/OatkDB}"
	local _crg_max_total_bp=3000000 # 3 Mb cutoff per-GFA
	local _crg_min_seq_len=10000    # keep contigs ≥ 10 kb
	local _crg_minimap2_min_identity=0.80
	local _crg_minimap2_min_qcov=0.50
	local _crg_read_min_qcov=0.50
	local _crg_plastid=false
	local _crg_parallel=false
	local _crg_redo=true

	# parse CLI
	while [[ $# -gt 0 ]]; do
		case "$1" in
		--version)
			_crg_version="on"
			shift
			;;
		-v | --verbose)
			_crg_verbose=$((_crg_verbose + 1))
			shift
			;;
		--steps)
			_crg_steps="$2"
			shift 2
			;;
		-l)
			_crg_long_reads="$2"
			shift 2
			;;
		-o)
			_crg_out_prefix="$2"
			shift 2
			;;
		-t)
			_crg_threads="$2"
			shift 2
			;;
		--k)
			_crg_k="$2"
			shift 2
			;;
		--s)
			_crg_s="$2"
			shift 2
			;;
		--c)
			_crg_c="$2"
			shift 2
			;;
		--a)
			_crg_a="$2"
			shift 2
			;;
		--max-bubble)
			_crg_max_bubble="$2"
			shift 2
			;;
		--max-tip)
			_crg_max_tip="$2"
			shift 2
			;;
		--weak-cross)
			_crg_weak_cross="$2"
			shift 2
			;;
		--unzip-round)
			_crg_unzip_round="$2"
			shift 2
			;;
		--no-read-ec)
			_crg_no_read_ec="true"
			shift
			;;
		--no-hpc)
			_crg_no_hpc="true"
			shift
			;;
		--plastid)
			_crg_plastid="true"
			shift
			;;
		--parallel)
			_crg_parallel="true"
			shift
			;;
		--oatkdb)
			_crg_oatkdb="$2"
			shift 2
			;;
		--max-total-bp)
			_crg_max_total_bp="$2"
			shift 2
			;;
		--min-seq-len)
			_crg_min_seq_len="$2"
			shift 2
			;;
		--minimap2-min-identity)
			_crg_minimap2_min_identity="$2"
			shift 2
			;;
		--minimap2-min-qcov)
			_crg_minimap2_min_qcov="$2"
			shift 2
			;;
		--read-min-qcov)
			_crg_read_min_qcov="$2"
			shift 2
			;;
		--)
			shift
			break
			;;
		*)
			((_crg_verbose)) && _polap_log1 "[INFO] ignoring unknown option: $1"
			shift
			;;
		esac
	done

	[[ "$_crg_version" == "on" ]] && {
		((_crg_verbose)) && _polap_log1 "polap-lib-syncasm v0.0.1"
		return 0
	}

	# sanity & execs
	[[ -z "$_crg_long_reads" || ! -s "$_crg_long_reads" ]] && {
		_polap_log0 "[ERROR] missing -l FASTQ(.gz) input"
		return 2
	}

	_polap_log0 "_crg_parallel=${_crg_parallel}"

	_polap__need_exec awk sed sort grep python Rscript bash oatk seqtk minimap2 samtools || return 2

	# expand + validate
	local K_LIST=($(_polap__expand_series_k "$_crg_k" 0))
	local S_LIST=($(_polap__expand_series "$_crg_s" 0))
	local C_LIST=($(_polap__expand_series "$_crg_c" 0))
	local A_LIST=($(_polap__expand_series "$_crg_a" 1))
	K_LIST=($(_polap__validate_list k "${K_LIST[@]}"))
	S_LIST=($(_polap__validate_list s "${S_LIST[@]}"))
	C_LIST=($(_polap__validate_list c "${C_LIST[@]}"))
	A_LIST=($(_polap__validate_list a "${A_LIST[@]}"))
	((${#K_LIST[@]} && ${#S_LIST[@]} && ${#C_LIST[@]} && ${#A_LIST[@]})) || {
		_polap_log0 "[ERROR] empty k/s/c/a after validation"
		return 2
	}
	((_crg_verbose)) && _polap_log1 "[INFO] k: ${K_LIST[*]}"
	((_crg_verbose)) && _polap_log1 "[INFO] s: ${S_LIST[*]}"
	((_crg_verbose)) && _polap_log1 "[INFO] c: ${C_LIST[*]}"
	((_crg_verbose)) && _polap_log1 "[INFO] a: ${A_LIST[*]}"

	local OATKDB_MT="${_crg_oatkdb}/v20230921/magnoliopsida_mito.fam"
	local OATKDB_PT="${_crg_oatkdb}/v20230921/magnoliopsida_pltd.fam"

	# steps
	local include_spec="$_crg_steps"
	[[ "$include_spec" == "all" ]] && include_spec="1-7"
	local _stage_array
	_stage_array=($(_polap_parse_steps "${include_spec}" ""))
	((_crg_verbose)) && _polap_log1 "[INFO] steps: ${_stage_array[*]}"

	# workspace
	local workdir
	_polap_log3_cmdout mkdir -p "${_crg_out_prefix}"
	workdir="${_crg_out_prefix}"
	# workdir="$(dirname -- "${_crg_out_prefix}")"
	# workdir == _arg_outdir
	[[ -z "$workdir" || "$workdir" == "." ]] && workdir="$(pwd)"
	_polap_log3_cmdout mkdir -p "$workdir"
	_polap_log1 "workdir: $workdir"

	# STEP 1: HPC (unless --no-hpc)
	local reads_in="${_crg_long_reads}"
	local reads_hpc="${workdir}/reads.hpc.fa"
	if _polap_contains_step x "${_stage_array[@]}"; then
		if [[ "$_crg_no_hpc" == "false" ]]; then
			((_crg_verbose)) && _polap_log1 "[STEP 1] seqtk hpc →  ${reads_hpc}"
			if [[ -s "${reads_hpc}" ]]; then
				_polap_log0 "[STEP 1] reuse HPC: $reads_hpc"
			else
				seqtk hpc "$reads_in" >"$reads_hpc"
			fi
		else
			((_crg_verbose)) && _polap_log1 "[STEP 1] skipping HPC"
			reads_hpc="$reads_in"
		fi
	fi

	# ───────────────────────────────────────────────────────────────
	# STEP 0: (optional) trim adapters + longest ~11Gb + HPC
	# Requires: porechop_abi, seqkit, filtlong
	# Controls (add to CLI parser if you want custom flags):
	#   _crg_step0_enable="true|false"
	#   _crg_step0_target_bases="11e9"
	#   _crg_step0_min_len="1000"
	#   _crg_step0_subset_n="100000"
	#   _crg_step0_seed="42"
	#   _crg_step0_threads="$((_crg_threads))"
	# ───────────────────────────────────────────────────────────────
	local _crg_step0_enable="${_crg_step0_enable:-true}"
	local _crg_step0_target_bases="${_crg_step0_target_bases:-11g}"
	local _crg_step0_min_len="${_crg_step0_min_len:-10000}"
	local _crg_step0_subset_n="${_crg_step0_subset_n:-100000}"
	local _crg_step0_seed="${_crg_step0_seed:-42}"
	local _crg_step0_threads="${_crg_step0_threads:-$((_arg_threads))}"

	# Output dir for step 0
	local step0_dir="${workdir}/trim11g"
	local step0_hpc="${step0_dir}/reads.hpc.fa"
	local trimmed_all="${step0_dir}/trimmed.fq"
	local trimmed_11g="${step0_dir}/trimmed.filtlong.fq"

	if _polap_contains_step 1 "${_stage_array[@]}"; then
		if [[ ! -s "${step0_hpc}" ]]; then
			((_crg_verbose)) && _polap_log1 "[STEP 0] porechop_abi + filtlong + HPC → ${step0_hpc}"

			# sanity: tools
			for _exe in seqkit filtlong; do
				command -v "$_exe" >/dev/null 2>&1 || {
					_polap_log0 "[ERROR] STEP 0: missing $_exe"
					return 2
				}
			done
			mkdir -p "${step0_dir}"

			_polap_lib_conda-ensure_conda_env polap-ont || exit 1
			porechop_abi -abi \
				-i "${reads_in}" \
				-o "${trimmed_all}" \
				--threads "${_crg_step0_threads}"
			conda deactivate

			# 0.3 keep longest reads up to target bases (optional min length)
			_polap_log1 "filtlong: ${_crg_step0_target_bases} of ${trimmed_all} -> ${trimmed_11g}"
			if ((_crg_step0_min_len > 0)); then
				filtlong --target_bases "${_crg_step0_target_bases}" --min_length "${_crg_step0_min_len}" \
					"${trimmed_all}" >"${trimmed_11g}" 2>/dev/null
			else
				filtlong --target_bases "${_crg_step0_target_bases}" \
					"${trimmed_all}" >"${trimmed_11g}" 2>/dev/null
			fi

			# 0.4 HPC (run AFTER trimming & downselect)
			seqtk hpc "${trimmed_11g}" >"${step0_hpc}"

		else
			((_crg_verbose)) && _polap_log1 "[STEP 0] disabled; skipping"
		fi
	fi

	reads_hpc="${step0_hpc}"

	# STEP 2: grid over c→a→k→s (GNU parallel; oatk direct; per-job -t 8)
	local grid_root="${workdir}/grid"
	if [[ "${_crg_redo}" == true ]]; then
		rm -rf "$grid_root"
	fi
	mkdir -p "$grid_root"
	if _polap_contains_step 2 "${_stage_array[@]}"; then
		((_crg_verbose)) && _polap_log1 "[STEP 2] oatk grid (parallel jobs = ${_crg_threads})"
		local cmdfile="${grid_root}/cmds.txt"
		: >"$cmdfile"
		local c a k s
		for c in "${C_LIST[@]}"; do
			for a in "${A_LIST[@]}"; do
				for k in "${K_LIST[@]}"; do
					for s in "${S_LIST[@]}"; do
						# normalize floats for CLI
						local a3 wx3 a_tag wx_tag c_tag k_tag s_tag u_tag
						a3="$(_polap__fmt_float3 "$a")" || {
							((_crg_verbose)) && _polap_log1 "[WARN] skip invalid -a '$a'"
							continue
						}
						wx3="$(_polap__fmt_float3 "$_crg_weak_cross")" || {
							((_crg_verbose)) && _polap_log1 "[WARN] skip invalid --weak-cross '$_crg_weak_cross'"
							continue
						}

						# _polap_log0 "a: $a"
						# _polap_log0 "a3: $a3"
						# _polap_log0 "wx: $_crg_weak_cross"
						# _polap_log0 "wx3: $wx3"
						#
						a_tag="$(_polap__tag_from_float3 "$a3")"
						wx_tag="$(_polap__tag_from_float3 "$wx3")"
						c_tag="$(_polap__pad3 "$c")"
						k_tag="$(_polap__pad4 "$k")"
						s_tag="$(_polap__pad2 "$s")"
						u_tag="$(_polap__pad2 "$_crg_unzip_round")"

						local run="exp_c${c_tag}_a${a_tag}_k${k_tag}_s${s_tag}_wx${wx_tag}_u${u_tag}"
						local rundir="${grid_root}/${run}"
						mkdir -p "$rundir"
						local prefix="${rundir}/syncasm.asm"

						# printf 'OUT=%q; mkdir -p "$OUT"; oatk -k %q -s %q -c %q -a %q --weak-cross %q --unzip-round %q %s -m %s -p %s -t 8 -o %q %q >"%s/job.out" 2>"%s/job.err"\n' \
						# 	"$rundir" "$k" "$s" "$c" "$a3" "$wx3" "$_crg_unzip_round" \
						# 	$([[ "$_crg_no_read_ec" == "true" ]] && echo "--no-read-ec" || echo "") \
						# 	"$OATKDB_MT" "$OATKDB_PT" \
						# 	"$prefix" "$reads_hpc" "$rundir" "$rundir" >>"$cmdfile"
						#
						# Assume these are set earlier
						# local _crg_no_read_ec="false"
						# local rundir k s c a3 wx3 _crg_unzip_round OATKDB_MT OATKDB_PT prefix reads_hpc grid_root cmdfile

						# 1) map boolean(s) to helpers
						local opt_no_read_ec=""
						[[ "$_crg_no_read_ec" == "true" ]] && opt_no_read_ec=1

						# 2) printf with parameter expansion for optional flag
						printf 'OUT=%q; mkdir -p "$OUT"; oatk -k %q -s %q -c %q -a %q --weak-cross %q --unzip-round %q %s -m %q -p %q -t 8 -o %q %q >"%s/job.out" 2>"%s/job.err"\n' \
							"$rundir" "$k" "$s" "$c" "$a3" "$wx3" "$_crg_unzip_round" \
							"${opt_no_read_ec:+--no-read-ec}" \
							"$OATKDB_MT" "$OATKDB_PT" \
							"$prefix" "$reads_hpc" "$rundir" "$rundir" >>"$cmdfile"

						((_crg_verbose)) && _polap_log1 "[INFO] job: ${run} (c=$c a=$a3 k=$k s=$s wx=$wx3 u=${_crg_unzip_round})"
					done
				done
			done
		done

		if [[ "${_crg_parallel}" == true ]]; then
			if command -v parallel >/dev/null 2>&1; then
				parallel -j "${_crg_threads}" --eta --lb --joblog "${grid_root}/parallel.log" "bash -lc {}" <"$cmdfile"
				# print failures after all jobs
				[[ -s "${grid_root}/parallel.log" ]] && awk -F'\t' 'NR>1 && $7!=0{print "# FAIL EXIT=" $7 "\tCMD=" $9}' "${grid_root}/parallel.log" | sed -n '1,120p'
			else
				((_crg_verbose)) && _polap_log1 "[INFO] GNU parallel not found; running serially"
				while IFS= read -r line; do bash -lc "$line"; done <"$cmdfile"
			fi
		else
			((_crg_verbose)) && _polap_log1 "[INFO] running serially"
			while IFS= read -r line; do bash -lc "$line"; done <"$cmdfile"
		fi
	fi

	# STEP 3: per-run summary (normalized name) for utg.final.gfa and utg.gfa
	if _polap_contains_step 3 "${_stage_array[@]}"; then
		local outtsv="${grid_root}/summary-utg.tsv"
		: >"$outtsv"
		echo -e "name\ttotal_bp_final\tn_unitigs_final\ttotal_bp_utg\tn_unitigs_utg" >"$outtsv"
		shopt -s nullglob
		local d raw name gfa_final gfa_utg L n Lu nu
		for d in "${grid_root}"/*/; do
			raw="$(basename "$d")"
			name="$(normalize_name "$raw")"
			gfa_final="${d}/syncasm.asm.utg.final.gfa"
			gfa_utg="${d}/syncasm.asm.utg.gfa"

			L="." n="." Lu="." nu="."
			# awk: count S and sum lengths (final)
			if [[ -s "$gfa_final" ]]; then
				read L n < <(_polap__gfa_stats "$gfa_final")
			fi
			# awk: count S and sum lengths (utg)
			if [[ -s "$gfa_utg" ]]; then
				read Lu nu < <(_polap__gfa_stats "$gfa_utg")
			fi
			printf "%s\t%s\t%s\t%s\t%s\n" "$name" "$L" "$n" "$Lu" "$nu" >>"$outtsv"
		done

		# sort by preferred metric: final total_bp desc, then final n asc (fallback to utg if final missing)
		local tmp="${outtsv}.tmp"
		awk 'BEGIN{FS=OFS="\t"} NR==1{print; next}{
			name=$1; bF=$2+0; nF=$3+0; bU=$4+0; nU=$5+0;
			kbp=(bF>0?bF:bU); kutg=(bF>0?nF:nU);
			print $0, kbp, kutg
		}' "$outtsv" | (
			read h
			echo "$h"
			sort -k6,6nr -k7,7n
		) | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,$5}' >"$tmp" && mv "$tmp" "$outtsv"

		((_crg_verbose)) && { head -n2 "$outtsv" | sed -n '2p' | awk '{print "[TOP1 summary]",$0}'; }
	fi

	# ────────────────────────────────────────────────────────────────
	# STEP 4: Merge/decide → ≥10kb → dedup → rename edge_#
	# Requirements: gfatools, seqkit
	# Outputs:
	#   ${final_dir}/final10kb.edge.fa  (numbered edge_1..N, description kept)
	# ────────────────────────────────────────────────────────────────
	local max_total_bp="${_crg_max_total_bp}" # 3 Mb cutoff per-GFA
	local min_len="${_crg_min_seq_len}"       # keep contigs ≥ 10 kb

	local final_dir="${workdir}/filtered"
	local merged_raw="${final_dir}/final10kb.merged.fa"
	local filtered_fa="${final_dir}/final10kb.filtered.fa"
	local dedup_fa="${final_dir}/final10kb.rmdup.fa"
	local final_fa="${final_dir}/final10kb.edge.fa" # ← renumbered output (edge_#)

	if _polap_contains_step 4 "${_stage_array[@]}"; then
		mkdir -p "$final_dir"
		: >"$merged_raw"

		shopt -s nullglob
		local gfa run total_bp
		for gfa in "${grid_root}"/*/syncasm.asm.utg.final.gfa; do
			[[ -s "$gfa" ]] || continue
			run="$(basename "$(dirname "$gfa")")"

			# 1) decide by total sequence length (sum_len) using gfatools + seqkit stats
			total_bp=$(
				gfatools gfa2fa "$gfa" 2>/dev/null |
					seqkit stats -Ta |
					awk 'NR==2{print $6}' # sum_len column
			)
			[[ -z "$total_bp" ]] && total_bp=0
			((total_bp < max_total_bp)) || continue
			((_crg_verbose)) && _polap_log1 "[STEP 4] keep: $gfa (sum_len=${total_bp})"

			# 2) convert S-lines to FASTA, preserve tags in description (LN/KC/SC run=…)
			#    >SEGID LN:i:.. KC:i:.. SC:f:.. run=RUN
			awk -v RUN="$run" '
      $1=="S"{
        name=$2; seq=$3; ln=""; kc=""; sc="";
        for(i=4;i<=NF;i++){
          if($i ~ /^LN:i:/) ln=$i;
          else if($i ~ /^KC:i:/) kc=$i;
          else if($i ~ /^SC:f:/) sc=$i;
        }
        if(seq!="*"){
          printf(">%s", name);
          if(ln!="") printf(" %s", ln);
          if(kc!="") printf(" %s", kc);
          if(sc!="") printf(" %s", sc);
          printf(" run=%s\n", RUN);
          print seq
        }
      }
    ' "$gfa" >>"$merged_raw"
		done

		# 3) filter by length (≥ 10 kb)
		_polap_log3_cmdout seqkit seq -m "$min_len" "$merged_raw" -o "$filtered_fa"

		# 4) deduplicate by sequence content (IDs ignored)
		_polap_log3_cmdout seqkit rmdup -s "$filtered_fa" -o "$dedup_fa"

		# 5) renumber to edge_1..N and keep the description (everything after the first space)
		#    If a header has no description, we write just ">edge_N".
		awk '
    /^>/ {
      line = substr($0,2)
      # find description (everything after first whitespace)
      sp = index(line," ")
      desc = (sp ? substr(line, sp+1) : "")
      ++i
      if (desc != "") printf(">edge_%d %s\n", i, desc)
      else            printf(">edge_%d\n", i)
      next
    }
    { print }
  ' "$dedup_fa" >"$final_fa"

		# verbose summary
		((_crg_verbose)) && {
			_polap_log1 "[STEP 4] merged:   $(grep -c '^>' "$merged_raw") entries"
			_polap_log1 "[STEP 4] ≥${min_len}bp: $(grep -c '^>' "$filtered_fa") entries"
			_polap_log1 "[STEP 4] unique:   $(grep -c '^>' "$dedup_fa") entries"
			_polap_log1 "[STEP 4] final:    $(grep -c '^>' "$final_fa") entries →  $final_fa"
		}
	fi

	_polap_lib_conda-ensure_conda_env polap || exit 1

	# ---------- STEP 5: spans without double counting (Option 2) ----------
	# Assumes Step 4 produced ${final_fa} or ${workdir}/filtered/final10kb.rmdup.fa

	local mt_cds="${_POLAPLIB_DIR}/polap-mt.1.c70.3.fna.hpc.fa"
	local pt_cds="${_POLAPLIB_DIR}/polap-pt.2.c70.3.fna.hpc.fa"
	local assembly_fa="${final_fa}" # or "${workdir}/filtered/final10kb.rmdup.fa"
	local outdir="${workdir}/annotation"
	local pt_paf="${outdir}/pt.paf"
	local mt_paf="${outdir}/mt.paf"
	local id_len_depth="${outdir}/id_len_depth.tsv"
	local pt_bed="${outdir}/pt.filtered.bed"
	local mt_bed="${outdir}/mt.filtered.bed"
	local pt_merged="${outdir}/pt.merged.bed"
	local mt_merged="${outdir}/mt.merged.bed"
	local pt_bp="${outdir}/pt.bp.tsv"
	local mt_bp="${outdir}/mt.bp.tsv"
	local master="${outdir}/assembly_info_organelle_annotation_count-all.txt"
	local mt_out="${outdir}/contig-annotation-depth-table.txt"
	local pt_out="${outdir}/pt-contig-annotation-depth-table.txt"
	if _polap_contains_step 5 "${_stage_array[@]}"; then
		mkdir -p "$outdir"

		# thresholds
		local MIN_ID="${_crg_minimap2_min_identity}"
		local MIN_QCOV="${_crg_minimap2_min_qcov}"
		# local MIN_ID="0.80"   # identity = nmatch/alnlen
		# local MIN_QCOV="0.50" # CDS (query) coverage = (qend-qstart)/qlen

		# sanity
		[[ -s "$assembly_fa" ]] || {
			_polap_log0 "[ERROR] STEP 5: assembly not found: $assembly_fa"
			return 2
		}

		# A) map CDS -> assembly (spliced; keep secondaries)
		if [[ -n "$pt_cds" && -s "$pt_cds" ]]; then
			minimap2 -t 16 -x splice -k11 --splice-flank=no --secondary=yes -N 50 -p 0.5 \
				-c "$assembly_fa" "$pt_cds" >"$pt_paf" 2>/dev/null
		else : >"$pt_paf"; fi

		if [[ -n "$mt_cds" && -s "$mt_cds" ]]; then
			minimap2 -t 16 -x splice -k11 --splice-flank=no --secondary=yes -N 50 -p 0.5 \
				-c "$assembly_fa" "$mt_cds" >"$mt_paf" 2>/dev/null
		else : >"$mt_paf"; fi

		# B) FASTA -> id,len,depth (python)
		python "${_POLAPLIB_DIR}/polap-py-fasta-len-depth.py" "$assembly_fa" "$id_len_depth" ||
			{
				_polap_log0 "[ERROR] STEP 5: python len/depth failed"
				return 2
			}

		# C) PAF -> BED (filter by ID/QCOV), sorted
		_paf_to_bed() {
			local paf="$1" bed="$2" idmin="$3" qcovmin="$4"
			awk -v ID="$idmin" -v QC="$qcovmin" '
      ($10>0 && $11>0 && $2>0){
        id=$10/$11; qcov=($4-$3)/$2;
        if (id>=ID && qcov>=QC){
          s=$8; e=$9; if (s>e){t=s;s=e;e=t}
          print $6"\t"s"\t"e
        }
      }' "$paf" | sort -k1,1 -k2,2n >"$bed"
		}
		_paf_to_bed "$pt_paf" "$pt_bed" "$MIN_ID" "$MIN_QCOV"
		_paf_to_bed "$mt_paf" "$mt_bed" "$MIN_ID" "$MIN_QCOV"

		# D) merge overlaps per contig
		if command -v bedtools >/dev/null 2>&1; then
			bedtools merge -i "$pt_bed" >"$pt_merged" 2>/dev/null || : >"$pt_merged"
			bedtools merge -i "$mt_bed" >"$mt_merged" 2>/dev/null || : >"$mt_merged"
		else
			_merge_bed() {
				awk 'BEGIN{OFS="\t"}
        NR==1{c=$1;s=$2;e=$3; next}
        { if($1==c && $2<=e){ if($3>e) e=$3 } else { print c,s,e; c=$1; s=$2; e=$3 } }
        END{ if(NR>0) print c,s,e }'
			}
			_merge_bed <"$pt_bed" >"$pt_merged"
			_merge_bed <"$mt_bed" >"$mt_merged"
		fi

		# E) sum merged spans per contig -> TSV
		# NEW: counts = number of merged intervals per contig (gene-like regions)
		local pt_n="${outdir}/pt.n.tsv"
		local mt_n="${outdir}/mt.n.tsv"
		awk 'BEGIN{OFS="\t"} {c[$1]++} END{for(k in c) print k,c[k]}' "$pt_merged" | sort -k1,1 >"$pt_n" || : >"$pt_n"
		awk 'BEGIN{OFS="\t"} {c[$1]++} END{for(k in c) print k,c[k]}' "$mt_merged" | sort -k1,1 >"$mt_n" || : >"$mt_n"

		# F) R (tidyverse) to write the 3 tables
		Rscript --vanilla "${_POLAPLIB_DIR}/polap-r-anno-span.R" \
			"$outdir" "$id_len_depth" "$pt_n" "$mt_n" \
			"$master" "$mt_out" "$pt_out" ||
			{
				_polap_log0 "[ERROR] STEP 5: R table generation failed"
				return 2
			}

		((_crg_verbose)) && {
			_polap_log1 "[STEP 5/opt2] wrote:"
			_polap_log1 "  - $master"
			_polap_log1 "  - $mt_out"
			_polap_log1 "  - $pt_out"
		}
	fi

	# STEP 6: seed-based assembly (skeleton)
	local pt_final_fa="${final_dir}/final10kb.edge.pt.fa" # ← renumbered output (edge_#)
	local mt_final_fa="${final_dir}/final10kb.edge.mt.fa" # ← renumbered output (edge_#)
	if _polap_contains_step 6 "${_stage_array[@]}"; then
		((_crg_verbose)) && _polap_log1 "[STEP 5] seed-based assembly kick"

		if [[ "${_crg_plastid}" == "true" ]]; then
			mkdir -p "${outdir}/pt/30-contigger"
			_polap_lib_lines-skip1 "${pt_out}" | cut -f1 >"${outdir}/pt/mt.contig.name-pt1"
			seqtk subseq "${final_fa}" "${outdir}/pt/mt.contig.name-pt1" >"${pt_final_fa}"
			# fasta -> gfa
			bash "${_POLAPLIB_DIR}/polap-bash-fa2gfa.sh" \
				-o "${outdir}/pt/30-contigger/graph_final.gfa" "$final_fa"
			final_fa="${pt_final_fa}"
		else
			mkdir -p "${outdir}/mt/30-contigger"
			_polap_lib_lines-skip1 "${mt_out}" | cut -f1 >"${outdir}/mt/mt.contig.name-mt1"
			seqtk subseq "${final_fa}" "${outdir}/mt/mt.contig.name-mt1" >"${mt_final_fa}"
			# fasta -> gfa
			bash "${_POLAPLIB_DIR}/polap-bash-fa2gfa.sh" \
				-o "${outdir}/mt/30-contigger/graph_final.gfa" "$final_fa"
			final_fa="${mt_final_fa}"
		fi

		#
		_polap_lib_syncasm-subset_raw_by_hpc_mapping -v \
			--raw "${_crg_long_reads}" \
			--contigs "${final_fa}" \
			--hpc-reads "${reads_hpc}" \
			--min-id "${_crg_minimap2_min_identity}" \
			--min-qcov "${_crg_read_min_qcov}" \
			--min-aln-len 3000 \
			--out "${outdir}"

		if [[ "${_crg_plastid}" == "true" ]]; then
			flye "${_arg_flye_data_type}" \
				"${outdir}.raw.mapped.fastq.gz" \
				-t "${_arg_threads}" \
				--out-dir "${outdir}"/pt \
				2>"${_polap_output_dest}"
		else
			flye "${_arg_flye_data_type}" \
				"${outdir}.raw.mapped.fastq.gz" \
				-t "${_arg_threads}" \
				--out-dir "${outdir}"/mt \
				2>"${_polap_output_dest}"
		fi

		# _polap_lib_assemble-rate \
		# 	-o "${outdir}" \
		# 	-l "${_crg_long_reads}" \
		# 	-w "${_arg_single_min}" \
		# 	-i pt0 -j pt1

	fi

	if _polap_contains_step 7 "${_stage_array[@]}"; then

		local annotatedir="${outdir}"
		local ptdna0_gfa="${annotatedir}/pt/30-contigger/graph_final.gfa"

		if [[ "${_crg_plastid}" == "true" ]]; then
			_polap_lib_syncsam-pt_iterate \
				--outdir "${annotatedir}" \
				--long-reads "${_arg_long_reads}"

		else
			_polap_lib_syncsam-mt_iterate \
				--outdir "${annotatedir}" \
				--long-reads "${_arg_long_reads}"
		fi

	fi

	if _polap_contains_step 0 "${_stage_array[@]}"; then

		local annotatedir="${outdir}"
		local ptdna0_gfa="${annotatedir}/pt/30-contigger/graph_final.gfa"

		_polap_log0_cmdout bash "${_POLAPLIB_DIR}/polap-bash-paf-diagnose.sh" \
			-i "${mt_paf}" \
			-o "${outdir}/mt_diag_primary" \
			--fixed-qcov 0.50 \
			--fixed-id 0.80 \
			--min-alnlen 200 \
			--primary-only

		_polap_log0_cmdout bash "${_POLAPLIB_DIR}/polap-bash-paf-diagnose.sh" \
			-i "${_crg_out_prefix}/annotation.paf" \
			-o "${_crg_out_prefix}/mt_read_diag_primary" \
			--fixed-qcov 0.50 \
			--fixed-id 0.80 \
			--min-alnlen 200 \
			--primary-only

	fi
	conda deactivate

	((_crg_verbose)) && _polap_log1 "[DONE] _polap_lib_syncasm finished."
	return 0
}
