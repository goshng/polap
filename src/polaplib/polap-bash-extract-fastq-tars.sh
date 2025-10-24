#!/usr/bin/env bash
# polap-bash-extract-fastq-tars.sh
# Version: v0.2.0
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Extract every *.fastq.tar.gz in INDIR:
#   - 1 FASTQ  -> <stem>.fastq.gz
#   - 2 FASTQs -> <stem>_1.fastq.gz and <stem>_2.fastq.gz
#   - >2 or unusual names -> <stem>_<basename>.fastq.gz
#
# pv stages:
#   scan:<tar>  — listing (bytes through archive; ETA)
#   x:<member>  — copy gz member
#   gz:<member> — gzip plain FASTQ member
#
# pigz:
#   If --threads N is given and pigz is installed, scanning and extraction
#   use `pigz -dc -p N` for parallel inflate (faster than tar -z).
set -euo pipefail
IFS=$'\n\t'

usage() {
	cat <<'EOF'
Usage:
  polap-bash-extract-fastq-tars.sh [--indir DIR] [--outdir DIR]
                                   [--threads N]
                                   [--overwrite | --no-overwrite]
                                   [--dry-run | --no-dry-run]
                                   [--rm-archives | --move-archives DIR]
                                   [-h|--help]

Defaults:
  --threads 4
  --overwrite
  --dry-run

Notes:
  - If --outdir is not specified, outputs go to --indir.
  - Archives are NEVER removed/moved unless you use --rm-archives or --move-archives.
EOF
}

# ---- options (defaults) ----
INDIR="."
OUTDIR="."
THREADS=4   # default pigz threads
OVERWRITE=1 # default overwrite ON
DRYRUN=1    # default dry-run
RM_ARCHIVES=0
MOVE_ARCHIVES=""
# defaults
THREADS=4 # you already set this
LEVEL=6   # gzip/pigz compression level (1=fast .. 9=max)

# helpers (keep your existing have()/pigz_ok())
comp_ok() { ((THREADS > 0)) && command -v pigz >/dev/null 2>&1; }

while (($#)); do
	case "$1" in
	# in arg parser
	--level)
		LEVEL="${2:?}"
		shift 2
		;;
	--indir)
		INDIR="${2:?}"
		shift 2
		;;
	--outdir)
		OUTDIR="${2:?}"
		shift 2
		;;
	--threads)
		THREADS="${2:?}"
		shift 2
		;;
	--overwrite)
		OVERWRITE=1
		shift
		;;
	--no-overwrite)
		OVERWRITE=0
		shift
		;;
	--dry-run)
		DRYRUN=1
		shift
		;;
	--no-dry-run)
		DRYRUN=0
		shift
		;;
	--rm-archives)
		RM_ARCHIVES=1
		shift
		;;
	--move-archives)
		MOVE_ARCHIVES="${2:?}"
		shift 2
		;;
	-h | --help)
		usage
		exit 0
		;;
	*)
		echo "[ERROR] Unknown option: $1" >&2
		usage
		exit 2
		;;
	esac
done

# ---- dir normalization (prevents // in paths) ----
normdir() {
	local d="$1"
	[[ "$d" == "/" ]] && {
		printf "/"
		return
	}
	printf "%s" "${d%/}"
}
INDIR="$(normdir "$INDIR")"
if [[ "$OUTDIR" == "." ]]; then OUTDIR="$INDIR"; fi
OUTDIR="$(normdir "$OUTDIR")"
[[ -n "$MOVE_ARCHIVES" ]] && MOVE_ARCHIVES="$(normdir "$MOVE_ARCHIVES")"

mkdir -p "$OUTDIR"

have() { command -v "$1" >/dev/null 2>&1; }
pigz_ok() { ((THREADS > 0)) && have pigz; }
log() { printf '[extract] %s\n' "$*" >&2; }

# cross-platform file size helper (Linux/macOS)
file_size() {
	local f="$1"
	if stat --version >/dev/null 2>&1; then
		stat -c %s -- "$f"
	else
		stat -f%z -- "$f"
	fi
}

# Stream one member to gz output; recompress if needed (pv if present)
# Args: archive member out.gz size_bytes_or_empty
stream_member_to_gz() {
	local arc="$1" mem="$2" out="$3" sz="${4:-}"

	if [[ -e "$out" && $OVERWRITE -eq 0 ]]; then
		log "skip (exists): $out"
		return 0
	fi
	mkdir -p "$(dirname "$out")"

	local base
	base="${mem##*/}"

	# Build pv -s as an array (robust with our IFS)
	local -a pv_s_opt=()
	[[ -n "$sz" ]] && pv_s_opt=(-s "$sz")

	if [[ "$mem" =~ \.fastq\.gz$ ]]; then
		# Already gzipped member: copy out
		if ((DRYRUN)); then
			if have pv; then
				if pigz_ok; then
					echo "[DRYRUN] pigz -dc -p $THREADS '$arc' | tar -xOf - '$mem' | pv -N 'x:${base}' -petbr${sz:+ -s $sz} > '$out'"
				else
					echo "[DRYRUN] tar -xOzf '$arc' '$mem' | pv -N 'x:${base}' -petbr${sz:+ -s $sz} > '$out'"
				fi
			else
				if pigz_ok; then
					echo "[DRYRUN] pigz -dc -p $THREADS '$arc' | tar -xOf - '$mem' > '$out'"
				else
					echo "[DRYRUN] tar -xOzf '$arc' '$mem' > '$out'"
				fi
			fi
		else
			if have pv; then
				if pigz_ok; then
					pigz -dc -p "$THREADS" -- "$arc" | tar -xOf - "$mem" | pv -N "x:${base}" -petbr "${pv_s_opt[@]}" >"$out"
				else
					tar -xOzf "$arc" "$mem" | pv -N "x:${base}" -petbr "${pv_s_opt[@]}" >"$out"
				fi
			else
				if pigz_ok; then
					pigz -dc -p "$THREADS" -- "$arc" | tar -xOf - "$mem" >"$out"
				else
					tar -xOzf "$arc" "$mem" >"$out"
				fi
			fi
		fi
	else

		# Plain FASTQ member: compress on the fly (prefer pigz -c -p N)
		if ((DRYRUN)); then
			if have pv; then
				if comp_ok; then
					echo "[DRYRUN] pigz -dc -p $THREADS '$arc' | tar -xOf - '$mem' | pv -N 'gz:${base}' -petbr${sz:+ -s $sz} | pigz -c -p $THREADS -$LEVEL > '$out'"
				else
					echo "[DRYRUN] tar -xOzf '$arc' '$mem' | pv -N 'gz:${base}' -petbr${sz:+ -s $sz} | gzip -c -$LEVEL > '$out'"
				fi
			else
				if comp_ok; then
					echo "[DRYRUN] pigz -dc -p $THREADS '$arc' | tar -xOf - '$mem' | pigz -c -p $THREADS -$LEVEL > '$out'"
				else
					echo "[DRYRUN] tar -xOzf '$arc' '$mem' | gzip -c -$LEVEL > '$out'"
				fi
			fi
		else
			if have pv; then
				if comp_ok; then
					# shellcheck disable=SC2086
					pigz -dc -p "$THREADS" -- "$arc" | tar -xOf - "$mem" |
						pv -N "gz:${base}" -petbr "${pv_s_opt[@]}" |
						pigz -c -p "$THREADS" -"$LEVEL" >"$out"
				else
					# shellcheck disable=SC2086
					tar -xOzf "$arc" "$mem" |
						pv -N "gz:${base}" -petbr "${pv_s_opt[@]}" |
						gzip -c -"$LEVEL" >"$out"
				fi
			else
				if comp_ok; then
					pigz -dc -p "$THREADS" -- "$arc" | tar -xOf - "$mem" |
						pigz -c -p "$THREADS" -"$LEVEL" >"$out"
				else
					tar -xOzf "$arc" "$mem" | gzip -c -"$LEVEL" >"$out"
				fi
			fi
		fi

	fi
}

# Collect candidate archives from INDIR (non-recursive)
shopt -s nullglob
tars=("$INDIR"/*.fastq.tar.gz)
shopt -u nullglob

if ((${#tars[@]} == 0)); then
	log "no *.fastq.tar.gz files found under: $INDIR"
	exit 0
fi

for arc in "${tars[@]}"; do
	base_arc="$(basename "$arc")"
	stem="${base_arc%.fastq.tar.gz}"
	log "archive: $base_arc (stem: $stem)"

	# ---- LIST once with pv (and pigz if enabled): collect names + sizes ----
	declare -a mems=()
	declare -A SZMAP=()

	if have pv; then
		arc_bytes="$(file_size "$arc" || echo "")"
		if ((DRYRUN)); then
			if [[ -n "$arc_bytes" ]]; then
				if pigz_ok; then
					echo "[DRYRUN] pv -N 'scan:${base_arc}' -petbr -s $arc_bytes '$arc' | pigz -dc -p $THREADS | tar -tvf -"
				else
					echo "[DRYRUN] pv -N 'scan:${base_arc}' -petbr -s $arc_bytes '$arc' | tar -tzvf -"
				fi
			else
				if pigz_ok; then
					echo "[DRYRUN] pv -N 'scan:${base_arc}' -petbr '$arc' | pigz -dc -p $THREADS | tar -tvf -"
				else
					echo "[DRYRUN] pv -N 'scan:${base_arc}' -petbr '$arc' | tar -tzvf -"
				fi
			fi
		fi

		# Real listing is required once to know members
		if [[ -n "$arc_bytes" ]]; then
			if pigz_ok; then
				while IFS=$'\t' read -r sz nm; do
					[[ "$nm" =~ \.fastq(\.gz)?$ ]] || continue
					mems+=("$nm")
					SZMAP["$nm"]="$sz"
				done < <(
					pv -N "scan:${base_arc}" -petbr -s "$arc_bytes" "$arc" |
						pigz -dc -p "$THREADS" |
						tar -tvf - 2>/dev/null |
						awk '{
              name=$NF; sz="";
              for(i=1;i<=NF;i++) if($i~/^[0-9]{4}-[0-9]{2}-[0-9]{2}$/){sz=$(i-1);break}
              if(sz!="" && name!="") print sz "\t" name
            }'
				)
			else
				while IFS=$'\t' read -r sz nm; do
					[[ "$nm" =~ \.fastq(\.gz)?$ ]] || continue
					mems+=("$nm")
					SZMAP["$nm"]="$sz"
				done < <(
					pv -N "scan:${base_arc}" -petbr -s "$arc_bytes" "$arc" |
						tar -tzvf - 2>/dev/null |
						awk '{
              name=$NF; sz="";
              for(i=1;i<=NF;i++) if($i~/^[0-9]{4}-[0-9]{2}-[0-9]{2}$/){sz=$(i-1);break}
              if(sz!="" && name!="") print sz "\t" name
            }'
				)
			fi
		else
			if pigz_ok; then
				while IFS=$'\t' read -r sz nm; do
					[[ "$nm" =~ \.fastq(\.gz)?$ ]] || continue
					mems+=("$nm")
					SZMAP["$nm"]="$sz"
				done < <(
					pv -N "scan:${base_arc}" -petbr "$arc" |
						pigz -dc -p "$THREADS" |
						tar -tvf - 2>/dev/null |
						awk '{
              name=$NF; sz="";
              for(i=1;i<=NF;i++) if($i~/^[0-9]{4}-[0-9]{2}-[0-9]{2}$/){sz=$(i-1);break}
              if(sz!="" && name!="") print sz "\t" name
            }'
				)
			else
				while IFS=$'\t' read -r sz nm; do
					[[ "$nm" =~ \.fastq(\.gz)?$ ]] || continue
					mems+=("$nm")
					SZMAP["$nm"]="$sz"
				done < <(
					pv -N "scan:${base_arc}" -petbr "$arc" |
						tar -tzvf - 2>/dev/null |
						awk '{
              name=$NF; sz="";
              for(i=1;i<=NF;i++) if($i~/^[0-9]{4}-[0-9]{2}-[0-9]{2}$/){sz=$(i-1);break}
              if(sz!="" && name!="") print sz "\t" name
            }'
				)
			fi
		fi
	else
		# no pv at all
		if pigz_ok; then
			while IFS=$'\t' read -r sz nm; do
				[[ "$nm" =~ \.fastq(\.gz)?$ ]] || continue
				mems+=("$nm")
				SZMAP["$nm"]="$sz"
			done < <(
				pigz -dc -p "$THREADS" -- "$arc" |
					tar -tvf - 2>/dev/null |
					awk '{
            name=$NF; sz="";
            for(i=1;i<=NF;i++) if($i~/^[0-9]{4}-[0-9]{2}-[0-9]{2}$/){sz=$(i-1);break}
            if(sz!="" && name!="") print sz "\t" name
          }'
			)
		else
			while IFS=$'\t' read -r sz nm; do
				[[ "$nm" =~ \.fastq(\.gz)?$ ]] || continue
				mems+=("$nm")
				SZMAP["$nm"]="$sz"
			done < <(
				tar -tzvf "$arc" 2>/dev/null |
					awk '{
            name=$NF; sz="";
            for(i=1;i<=NF;i++) if($i~/^[0-9]{4}-[0-9]{2}-[0-9]{2}$/){sz=$(i-1);break}
            if(sz!="" && name!="") print sz "\t" name
          }'
			)
		fi
	fi

	if ((${#mems[@]} == 0)); then
		log "  no FASTQ members found -> skip"
		continue
	fi

	if ((${#mems[@]} == 1)); then
		out="$OUTDIR/${stem}.fastq.gz"
		log "  1 member -> $out"
		sz="${SZMAP[${mems[0]}]:-}"
		stream_member_to_gz "$arc" "${mems[0]}" "$out" "$sz"
	else
		# detect R1/R2
		r1=""
		r2=""
		for m in "${mems[@]}"; do
			base="$(basename "$m")"
			if [[ "$base" =~ _1\.fastq(\.gz)?$ ]]; then
				r1="$m"
			elif [[ "$base" =~ _2\.fastq(\.gz)?$ ]]; then
				r2="$m"
			fi
		done
		if [[ -n "$r1" || -n "$r2" ]]; then
			if [[ -n "$r1" ]]; then
				out1="$OUTDIR/${stem}_1.fastq.gz"
				log "  R1 -> $out1"
				stream_member_to_gz "$arc" "$r1" "$out1" "${SZMAP[$r1]:-}"
			fi
			if [[ -n "$r2" ]]; then
				out2="$OUTDIR/${stem}_2.fastq.gz"
				log "  R2 -> $out2"
				stream_member_to_gz "$arc" "$r2" "$out2" "${SZMAP[$r2]:-}"
			fi
			# extras, if any
			for m in "${mems[@]}"; do
				[[ "$m" == "$r1" || "$m" == "$r2" ]] && continue
				base="$(basename "$m")"
				safe="${base//[^A-Za-z0-9._-]/_}"
				out="$OUTDIR/${stem}_${safe%.fastq*}.fastq.gz"
				log "  extra -> $out"
				stream_member_to_gz "$arc" "$m" "$out" "${SZMAP[$m]:-}"
			done
		else
			# multi-members without _1/_2 -> map each to <stem>_<basename>.fastq.gz
			for m in "${mems[@]}"; do
				base="$(basename "$m")"
				safe="${base//[^A-Za-z0-9._-]/_}"
				out="$OUTDIR/${stem}_${safe%.fastq*}.fastq.gz"
				log "  member -> $out"
				stream_member_to_gz "$arc" "$m" "$out" "${SZMAP[$m]:-}"
			done
		fi
	fi

	# post-action on archive
	if ((DRYRUN == 0)); then
		if [[ -n "$MOVE_ARCHIVES" ]]; then
			mkdir -p "$MOVE_ARCHIVES"
			log "  move -> $MOVE_ARCHIVES/$base_arc"
			mv -n -- "$arc" "$MOVE_ARCHIVES"/
		elif ((RM_ARCHIVES)); then
			log "  remove -> $base_arc"
			rm -f -- "$arc"
		fi
	else
		if [[ -n "$MOVE_ARCHIVES" ]]; then
			echo "[DRYRUN] mv -n -- '$arc' '$MOVE_ARCHIVES/'"
		elif ((RM_ARCHIVES)); then
			echo "[DRYRUN] rm -f -- '$arc'"
		fi
	fi
done

log "done."
