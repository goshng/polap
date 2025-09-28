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

# Fetch organelle (mtDNA/ptDNA) sequences in FASTA using NCBI E-utilities.
# Requires: esearch, efetch  (Entrez Direct)
#
# Usage examples:
#   polap_fetch_organelle_fasta -s "Macadamia tetraphylla" -t mt -o o/mt.fa
#   polap_fetch_organelle_fasta -s "Arabidopsis thaliana" -t pt -o o/pt.fa --without-genome
#
# Returns:
#   0 on success (outfile exists and is non-empty)
#   1 on failure (no records written)
#   2 on usage error
#   3 if required tools are missing
_polap_lib_mtdna-fetch_organelle_fasta() {
	local species=""
	local otype="mt"
	local outfile=""
	local without_genome=0

	local usage
	usage() {
		cat <<'EOF' >&2
Usage:
  polap_fetch_organelle_fasta -s <species name> -t <mt|pt> -o <out.fa> [--without-genome]

Options:
  -s  Species scientific name (e.g., "Arabidopsis thaliana")
  -t  Organelle type: mt (mitochondrion) or pt (plastid/chloroplast)
  -o  Output FASTA file path
      --without-genome  Omit 'genome[Title]' in the query (fetch complete sequences not tagged as 'genome')
EOF
	}

	# parse args
	while (($#)); do
		case "$1" in
		-s)
			species="${2:?}"
			shift 2
			;;
		-t)
			otype="${2:?}"
			shift 2
			;;
		-o)
			outfile="${2:?}"
			shift 2
			;;
		--without-genome)
			without_genome=1
			shift
			;;
		-h | --help)
			usage
			return 0
			;;
		*)
			echo "Unknown option: $1" >&2
			usage
			return 2
			;;
		esac
	done

	[[ -n "$species" && -n "$otype" && -n "$outfile" ]] || {
		usage
		return 2
	}

	# tool checks
	command -v esearch >/dev/null 2>&1 || {
		echo "ERROR: missing tool: esearch" >&2
		return 3
	}
	command -v efetch >/dev/null 2>&1 || {
		echo "ERROR: missing tool: efetch" >&2
		return 3
	}

	# ensure output dir exists
	mkdir -p -- "$(dirname -- "$outfile")"

	# build query & fetch
	if [[ "$otype" == "mt" ]]; then
		local query
		if ((without_genome)); then
			query="(mitochondrion[Title] AND complete[Title]) AND ${species}[Organism]"
		else
			query="(mitochondrion[Title] AND complete[Title] AND genome[Title]) AND ${species}[Organism]"
		fi
		esearch -db nuccore -query "$query" | efetch -format fasta >"$outfile" || true

	elif [[ "$otype" == "pt" ]]; then
		# try chloroplast first, fallback to plastid
		local q1 q2
		if ((without_genome)); then
			q1="(chloroplast[Title] AND complete[Title]) AND ${species}[Organism]"
			q2="(plastid[Title] AND complete[Title]) AND ${species}[Organism]"
		else
			q1="(chloroplast[Title] AND complete[Title] AND genome[Title]) AND ${species}[Organism]"
			q2="(plastid[Title] AND complete[Title] AND genome[Title]) AND ${species}[Organism]"
		fi
		esearch -db nuccore -query "$q1" | efetch -format fasta >"$outfile" || true
		[[ -s "$outfile" ]] || esearch -db nuccore -query "$q2" | efetch -format fasta >"$outfile" || true

	else
		echo "ERROR: -t must be 'mt' or 'pt'." >&2
		return 2
	fi

	# final check
	if [[ ! -s "$outfile" ]]; then
		echo "No organelle FASTA retrieved for species='${species}', type='${otype}'" >&2
		return 1
	fi
	return 0
}

# Requires: esearch, efetch (Entrez Direct). Optional: seqkit (for min length filter).
# Returns non-zero on failure.

# 1) Fetch organelle FASTA (mt or pt), optionally omitting 'genome[Title]'.
polap_fetch_organelle_fasta() {
	local species="" otype="" outfile="" without_genome=0
	while (($#)); do
		case "$1" in
		-s | --species)
			species="${2:?}"
			shift 2
			;;
		-t | --type)
			otype="${2:?}"
			shift 2
			;; # mt|pt
		-o | --out)
			outfile="${2:?}"
			shift 2
			;;
		--without-genome)
			without_genome=1
			shift
			;;
		-h | --help)
			echo "Usage: polap_fetch_organelle_fasta -s 'Species name' -t mt|pt -o out.fa [--without-genome]" >&2
			return 0
			;;
		*)
			echo "Unknown option: $1" >&2
			return 2
			;;
		esac
	done

	[[ -n "$species" && -n "$otype" && -n "$outfile" ]] || {
		echo "Missing -s/-t/-o" >&2
		return 2
	}
	command -v esearch >/dev/null 2>&1 || {
		echo "Missing esearch" >&2
		return 3
	}
	command -v efetch >/dev/null 2>&1 || {
		echo "Missing efetch" >&2
		return 3
	}
	mkdir -p -- "$(dirname -- "$outfile")"

	if [[ "$otype" == "mt" ]]; then
		local q
		if ((without_genome)); then
			q="(mitochondrion[Title] AND complete[Title]) AND ${species}[Organism]"
		else
			q="(mitochondrion[Title] AND complete[Title] AND genome[Title]) AND ${species}[Organism]"
		fi
		esearch -db nuccore -query "$q" | efetch -format fasta >"$outfile" || true

	elif [[ "$otype" == "pt" ]]; then
		local q1 q2
		if ((without_genome)); then
			q1="(chloroplast[Title] AND complete[Title]) AND ${species}[Organism]"
			q2="(plastid[Title] AND complete[Title]) AND ${species}[Organism]"
		else
			q1="(chloroplast[Title] AND complete[Title] AND genome[Title]) AND ${species}[Organism]"
			q2="(plastid[Title] AND complete[Title] AND genome[Title]) AND ${species}[Organism]"
		fi
		esearch -db nuccore -query "$q1" | efetch -format fasta >"$outfile" || true
		[[ -s "$outfile" ]] || esearch -db nuccore -query "$q2" | efetch -format fasta >"$outfile" || true
	else
		echo "-t must be mt or pt" >&2
		return 2
	fi

	[[ -s "$outfile" ]] || {
		echo "No organelle FASTA for ${species} (${otype})" >&2
		return 1
	}
}

# 2) List sequences (id, length, description) and choose one representative.
# Strategies: longest (default) | first | match (requires --pattern REGEX).
polap_list_and_select_rep() {
	set -euo pipefail
	local fasta="" tsv="" rep_fa="" rep_id="" strategy="longest" pattern=""
	while (($#)); do
		case "$1" in
		--fasta)
			fasta="${2:?}"
			shift 2
			;;
		--tsv)
			tsv="${2:?}"
			shift 2
			;;
		--rep-fasta)
			rep_fa="${2:?}"
			shift 2
			;;
		--rep-id)
			rep_id="${2:?}"
			shift 2
			;;
		--strategy)
			strategy="${2:?}"
			shift 2
			;;
		--pattern)
			pattern="${2:-}"
			shift 2
			;;
		-h | --help)
			echo "Usage: polap_list_and_select_rep --fasta in.fa --tsv list.tsv --rep-fasta rep.fa --rep-id rep.id [--strategy longest|first|match --pattern REGEX]" >&2
			return 0
			;;
		*)
			echo "Unknown option: $1" >&2
			return 2
			;;
		esac
	done
	[[ -n "$fasta" && -n "$tsv" && -n "$rep_fa" && -n "$rep_id" ]] || {
		echo "Missing required args" >&2
		return 2
	}
	[[ -s "$fasta" ]] || {
		echo "FASTA empty: $fasta" >&2
		return 1
	}
	mkdir -p -- "$(dirname -- "$tsv")" "$(dirname -- "$rep_fa")" "$(dirname -- "$rep_id")"

	# id \t length \t description
	awk 'BEGIN{OFS="\t"}
    /^>/{
      if (NR>1) print id, len, desc;
      line=substr($0,2);
      split(line,a,/[\t ]+/); id=a[1];
      desc=line; sub(/^[^ \t]+[ \t]*/,"",desc);
      len=0; next
    }
    !/^>/{
      gsub(/[ \t\r]/,""); len+=length($0)
    }
    END{ if (NR>0) print id, len, desc }
  ' "$fasta" >"$tsv"

	local chosen_id=""
	case "$strategy" in
	longest) chosen_id="$(awk -F'\t' '($2+0)>m{m=$2;id=$1} END{print id}' "$tsv")" ;;
	first) chosen_id="$(head -n1 "$tsv" | cut -f1)" ;;
	match)
		if [[ -z "$pattern" ]]; then
			chosen_id="$(awk -F'\t' '($2+0)>m{m=$2;id=$1} END{print id}' "$tsv")"
		else
			chosen_id="$(
				awk -F'\t' -v pat="$pattern" 'BEGIN{IGNORECASE=1}
            $1 ~ pat || $3 ~ pat { if(($2+0)>m){m=$2; id=$1} }
            END{ if(id!="") print id }
          ' "$tsv"
			)"
			[[ -n "$chosen_id" ]] || chosen_id="$(awk -F'\t' '($2+0)>m{m=$2;id=$1} END{print id}' "$tsv")"
		fi
		;;
	*) chosen_id="$(awk -F'\t' '($2+0)>m{m=$2;id=$1} END{print id}' "$tsv")" ;;
	esac
	[[ -n "$chosen_id" ]] || {
		echo "Failed to choose representative" >&2
		return 1
	}
	printf '%s\n' "$chosen_id" >"$rep_id"

	# extract representative record by exact id (first token after '>')
	awk -v ID="$chosen_id" '
    BEGIN{p=0}
    /^>/{
      hdr=substr($0,2); split(hdr,a,/[\t ]+/); p = (a[1]==ID)
    }
    { if(p) print }
  ' "$fasta" >"$rep_fa"

	[[ -s "$rep_fa" ]] || {
		echo "Representative FASTA empty" >&2
		return 1
	}
}

# 3) Wrapper: fetch + list + select (one call)
# Example:
#   polap_fetch_select_organelle \
#     --species "Macadamia tetraphylla" --type mt \
#     --out-all o/mt.all.fa --out-tsv o/mt.list.tsv \
#     --out-rep-fa o/mt.rep.fa --out-rep-id o/mt.rep.id \
#     --strategy longest  [--without-genome]  [--minlen 100000]  [--pattern REGEX]
polap_fetch_select_organelle() {
	local species="" type="" out_all="" out_tsv="" out_rep_fa="" out_rep_id=""
	local without_genome=0 strategy="longest" pattern="" minlen=""
	while (($#)); do
		case "$1" in
		--species)
			species="${2:?}"
			shift 2
			;;
		--type)
			type="${2:?}"
			shift 2
			;; # mt|pt
		--out-all)
			out_all="${2:?}"
			shift 2
			;;
		--out-tsv)
			out_tsv="${2:?}"
			shift 2
			;;
		--out-rep-fa)
			out_rep_fa="${2:?}"
			shift 2
			;;
		--out-rep-id)
			out_rep_id="${2:?}"
			shift 2
			;;
		--strategy)
			strategy="${2:?}"
			shift 2
			;;
		--pattern)
			pattern="${2:-}"
			shift 2
			;;
		--without-genome)
			without_genome=1
			shift
			;;
		--minlen)
			minlen="${2:?}"
			shift 2
			;; # optional, needs seqkit
		-h | --help)
			cat <<EOF >&2
Usage:
  polap_fetch_select_organelle \\
    --species "Arabidopsis thaliana" --type mt|pt \\
    --out-all o/all.fa --out-tsv o/list.tsv \\
    --out-rep-fa o/rep.fa --out-rep-id o/rep.id \\
    [--strategy longest|first|match] [--pattern REGEX] [--without-genome] [--minlen N]
EOF
			return 0
			;;
		*)
			echo "Unknown option: $1" >&2
			return 2
			;;
		esac
	done

	[[ -n "$species" && -n "$type" && -n "$out_all" && -n "$out_tsv" && -n "$out_rep_fa" && -n "$out_rep_id" ]] ||
		{
			echo "Missing required args" >&2
			return 2
		}

	polap_fetch_organelle_fasta -s "$species" -t "$type" -o "$out_all" $([[ $without_genome -eq 1 ]] && echo --without-genome)
	[[ -s "$out_all" ]] || {
		echo "Fetch returned no sequences" >&2
		return 1
	}

	# Optional min length filter if seqkit available
	if [[ -n "$minlen" ]] && command -v seqkit >/dev/null 2>&1; then
		local tmp="${out_all}.min.$$"
		seqkit seq -m "$minlen" "$out_all" >"$tmp" && mv "$tmp" "$out_all"
		[[ -s "$out_all" ]] || {
			echo "No sequences passed --minlen $minlen" >&2
			return 1
		}
	fi

	polap_list_and_select_rep --fasta "$out_all" --tsv "$out_tsv" \
		--rep-fasta "$out_rep_fa" --rep-id "$out_rep_id" \
		--strategy "$strategy" $([[ "$strategy" == match && -n "$pattern" ]] && printf -- "--pattern %s" "$pattern")
}
