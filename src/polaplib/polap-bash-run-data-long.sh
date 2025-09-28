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

#!/usr/bin/env bash
# polap-bash-run-data-long.sh
# Normalize a long-read dataset to: $outdir/tmp/l.fq
# Extension priority: .fastq, .fq, .fastq.gz, .fq.gz, .fastq.tar.gz, .fq.tar.gz, -10x.fastq.tar.gz
# Directory order: $PWD first, then --media, --media1, --media2, then remote (if provided).
set -euo pipefail

usage() {
	cat <<'EOF'
Usage:
  polap-bash-run-data-long.sh -l <SRA_ID> -o <outdir> [options]

Required:
  -l, --long-sra  SRA accession (e.g., SRR123456)
  -o, --outdir    Output directory (final file: outdir/tmp/l.fq)

Optional:
  --media         First media directory (searched after $PWD)
  --media1        Second media directory
  --media2        Third media directory
  -r, --remote    Remote host for ssh/scp (e.g., user@host)
  -n, --dry-run   Print actions without making changes
  -h, --help      Show this help

Exit codes:
  0  success
  2  usage error
  3  required tool missing
  20 fetch produced no files
  21 fetch failed
  22 normalize failed
EOF
}

# -------- parse args
long_sra=""
outdir=""
media_dir="/media/h2/sra"
media1_dir="/media/h1/sra"
media2_dir="/media/h2/sra"
remote_host="thorne"
dry_run="false"

while (("$#")); do
	case "$1" in
	-l | --long-sra)
		long_sra="${2:?}"
		shift 2
		;;
	-o | --outdir)
		outdir="${2:?}"
		shift 2
		;;
	--media)
		media_dir="${2:-}"
		shift 2
		;;
	--media1)
		media1_dir="${2:-}"
		shift 2
		;;
	--media2)
		media2_dir="${2:-}"
		shift 2
		;;
	-r | --remote)
		remote_host="${2:-}"
		shift 2
		;;
	-n | --dry-run)
		dry_run="true"
		shift
		;;
	-h | --help)
		usage
		exit 0
		;;
	--)
		shift
		break
		;;
	-*)
		echo "Error: Unknown option $1" >&2
		usage
		exit 2
		;;
	*)
		echo "Error: Unexpected positional $1" >&2
		usage
		exit 2
		;;
	esac
done

[[ -n "$long_sra" ]] || {
	echo "Error: -l|--long-sra is required." >&2
	usage
	exit 2
}
[[ -n "$outdir" ]] || {
	echo "Error: -o|--outdir is required." >&2
	usage
	exit 2
}

# tools (optional here; we’ll guard fetch calls anyway)
need_tool() { command -v "$1" >/dev/null 2>&1 || return 1; }

# -------- paths & helpers
tmpd="${outdir%/}/tmp"
mkdir -p "$tmpd"
norm="${tmpd}/l.fq"
work="${tmpd}/.${long_sra}.work.$$"
mkdir -p "$work"

log() { printf '[polap-long] %s\n' "$*" >&2; }
maybe() { if [[ "$dry_run" == "true" ]]; then log "(dry-run) $*"; else eval "$@"; fi; }

cleanup() {
	rm -rf -- "$work"
}
trap cleanup EXIT

# ---- OpenSSL mismatch guard ----------------------------------------------
# Runs a command; on failure, if stderr includes "OpenSSL version mismatch",
# warns and returns 100 (skip), otherwise returns the real non-zero code.
openssl_mismatch_guard() {
	local stderrf="${work}/.stderr.$RANDOM"
	set +e
	"$@" 2>"$stderrf"
	local rc=$?
	set -e
	if ((rc != 0)); then
		if grep -qi 'OpenSSL version mismatch' "$stderrf"; then
			log "WARN(OpenSSL mismatch): $* — skipping this step"
			rm -f "$stderrf"
			return 100
		fi
		cat "$stderrf" >&2
		rm -f "$stderrf"
		return $rc
	fi
	rm -f "$stderrf"
	return 0
}

# -------- normalization helpers
concat_to_norm() {
	# args: files... (mix of .fastq/.fq and gz)
	if [[ "$dry_run" == "true" ]]; then
		log "(dry-run) would concatenate ${#@} file(s) -> $norm"
		return 0
	fi
	: >"$norm"
	local f
	for f in "$@"; do
		if [[ "$f" == *.gz ]]; then
			gzip -dc -- "$f" >>"$norm"
		else
			cat -- "$f" >>"$norm"
		fi
	done
}

emit_from_one_file() {
	# copy/decompress into fixed filename $norm (final name is literally l.fq)
	local src="$1"

	# size of the source file (human-readable)
	if [[ -f "$src" ]]; then
		log "Source $src size: $(du -h "$src" | awk '{print $1}')"
	else
		log "WARN: source $src not found"
	fi

	case "$src" in
	*.fq | *.fastq)
		maybe "cat -- '$src' > '$norm'"
		;;
	*.fq.gz | *.fastq.gz)
		log "Decompressing gzip -> $norm"
		maybe "gzip -dc -- '$src' > '$norm'"
		;;
	*)
		return 1
		;;
	esac

	# size of the normalized file (human-readable)
	if [[ -f "$norm" ]]; then
		log "Created $norm size: $(du -h "$norm" | awk '{print $1}')"
	else
		log "WARN: $norm not created"
	fi
}

emit_from_tar() {
	local tarpath="$1"
	log "Extracting archive: $tarpath"
	if [[ "$dry_run" == "true" ]]; then
		log "(dry-run) would 'tar -tzf' and 'tar -zxf' into $work, then concat *.{fq,fastq}[.gz] -> $norm"
		return 0
	fi
	tar -tzf "$tarpath" >/dev/null || {
		log "Archive list failed: $tarpath"
		return 1
	}
	tar -zxf "$tarpath" -C "$work"
	mapfile -t files < <(find "$work" -type f \( -iname "*.fq" -o -iname "*.fastq" -o -iname "*.fq.gz" -o -iname "*.fastq.gz" \) | sort -V)
	((${#files[@]} > 0)) || {
		log "No FASTQ found inside archive: $tarpath"
		return 1
	}
	concat_to_norm "${files[@]}"
}

emit_from_tar() {
	# Extract an archive to $work and concat FASTQs to $norm.
	# Supports: .tar.gz (gzipped tar), .tar (plain), .zip. Logs sizes before/after.
	# Uses globals: dry_run, work, norm, log(), concat_to_norm()

	local apath="$1"
	[[ -f "$apath" ]] || {
		log "ERROR: archive not found: $apath"
		return 1
	}

	# Human-readable source size
	local src_size
	src_size="$(du -h -- "$apath" | awk '{print $1}')" || src_size="?"
	log "Extracting archive: $apath"
	log "Archive size: $src_size"

	# Identify container
	local mt
	mt="$(file -b --mime-type -- "$apath" 2>/dev/null || true)"
	# Common mimes: application/x-gzip, application/gzip, application/x-tar, application/zip
	[[ -z "$mt" ]] && mt="(unknown)"

	# Dry-run: try to list members anyway for better UX
	if [[ "${dry_run:-false}" == "true" ]]; then
		case "$mt" in
		application/gzip | application/x-gzip)
			if tar -tzf -- "$apath" >/dev/null 2>&1; then
				local n
				n="$(tar -tzf -- "$apath" |
					awk 'BEGIN{IGNORECASE=1} /\.f(ast)?q(\.gz)?$/ {c++} END{print c+0}')"
				log "(dry-run) would extract to: $work"
				log "(dry-run) would find ${n} FASTQ files and concat → $norm"
			else
				log "(dry-run) would try: tar -tzf (list) then tar -zxf (extract) → $work → $norm"
			fi
			return 0
			;;
		application/x-tar)
			if tar -tf -- "$apath" >/dev/null 2>&1; then
				local n
				n="$(tar -tf -- "$apath" |
					awk 'BEGIN{IGNORECASE=1} /\.f(ast)?q(\.gz)?$/ {c++} END{print c+0}')"
				log "(dry-run) would extract to: $work"
				log "(dry-run) would find ${n} FASTQ files and concat → $norm"
			else
				log "(dry-run) would try: tar -tf (list) then tar -xf (extract) → $work → $norm"
			fi
			return 0
			;;
		application/zip)
			if unzip -l -- "$apath" >/dev/null 2>&1; then
				local n
				n="$(unzip -l -- "$apath" 2>/dev/null |
					awk 'BEGIN{IGNORECASE=1} /\.f(ast)?q(\.gz)?$/ {c++} END{print c+0}')"
				log "(dry-run) would unzip to: $work ; find ${n} FASTQ files → concat → $norm"
			else
				log "(dry-run) would try: unzip -l then unzip → $work → $norm"
			fi
			return 0
			;;
		*)
			log "(dry-run) unknown mime-type: $mt ; would attempt tar -tzf → tar -tf → unzip -l."
			return 0
			;;
		esac
	fi

	# Real run: list → extract with fallbacks
	mkdir -p -- "$work"

	local listed=0
	case "$mt" in
	application/gzip | application/x-gzip)
		# Verify gzip integrity; if OK, try tar -tzf, else maybe it's a single .fastq.gz misnamed.
		log "tar -xzf $apath -C $work"
		tar -xzf "$apath" -C "$work"
		# if ! tar -tzf "$apath" 2>/dev/null 2>&1; then
		# 	log "ERROR: gzip integrity check failed (corrupt .gz?): $apath"
		# 	return 1
		# fi
		# if tar -tzf -- "$apath" >/dev/null 2>&1; then
		# 	listed=1
		# 	tar -zxf -- "$apath" -C "$work"
		# else
		# 	# Not a tarball inside gzip (e.g., just a single .fastq.gz, misnamed .tar.gz)
		# 	log "WARN: gzip OK but not a tarball; treating as a single gzipped FASTQ"
		# 	local fastq_out="$work/$(basename "${apath%.tar.gz}").fastq"
		# 	log "tar -xOzf $apath >$fastq_out"
		# 	log tar -xOzf "$apath" >"$fastq_out"
		# fi
		;;
	application/x-tar)
		if tar -tf -- "$apath" >/dev/null 2>&1; then
			listed=1
			tar -xf -- "$apath" -C "$work"
		else
			log "ERROR: tar list failed (plain tar): $apath"
			return 1
		fi
		;;
	application/zip)
		if unzip -l -- "$apath" >/dev/null 2>&1; then
			listed=1
			unzip -q -- "$apath" -d "$work"
		else
			log "ERROR: unzip list failed (zip): $apath"
			return 1
		fi
		;;
	*)
		# Fallback probe order if mime-type was unknown/misleading
		if tar -tzf -- "$apath" >/dev/null 2>&1; then
			listed=1
			tar -zxf -- "$apath" -C "$work"
		elif tar -tf -- "$apath" >/dev/null 2>&1; then
			listed=1
			tar -xf -- "$apath" -C "$work"
		elif unzip -l -- "$apath" >/dev/null 2>&1; then
			listed=1
			unzip -q -- "$apath" -d "$work"
		else
			log "ERROR: could not recognize or list archive format: $apath"
			return 1
		fi
		;;
	esac

	# Gather FASTQs (case-insensitive), then concat → $norm
	# This finds all FASTQ/FASTQ.GZ files under $work, sorts them naturally
	# (version-aware), and saves the list into the Bash array files.
	mapfile -t files < <(
		find "$work" -type f \( -iname "*.fq" -o -iname "*.fastq" -o -iname "*.fq.gz" -o -iname "*.fastq.gz" \) \
			-print0 | xargs -0 -I{} printf "%s\n" "{}" | sort -V
	)
	if ((${#files[@]} == 0)); then
		log "No FASTQ files found after extraction: $apath"
		return 1
	fi

	log "Found ${#files[@]} FASTQ files; concatenating → $norm"
	concat_to_norm "${files[@]}"

	if [[ -s "$norm" ]]; then
		local out_size
		out_size="$(du -h -- "$norm" | awk '{print $1}')" || out_size="?"
		log "Created $norm size: $out_size"
		return 0
	else
		log "❌ Normalized output missing or empty: $norm"
		return 1
	fi
}

# search functions (exact-name only, in priority order)
exts=(
	".fastq"
	".fq"
	".fastq.gz"
	".fq.gz"
	".fastq.tar.gz"
	".fq.tar.gz"
	"-10x.fastq.tar.gz" # NEW
)

try_exact_in_dir() {
	local dir="$1" p
	for ext in "${exts[@]}"; do
		p="${dir%/}/${long_sra}${ext}"
		if [[ -s "$p" ]]; then
			log "Found: $p"
			case "$p" in
			*.tar.gz) emit_from_tar "$p" ;;
			*) emit_from_one_file "$p" ;;
			esac
			return 0
		fi
	done
	return 1
}

try_exact_remote_dir() {
	local dir="$1" rp base local_copy
	[[ -z "$remote_host" ]] && return 1
	for ext in "${exts[@]}"; do
		rp="${dir%/}/${long_sra}${ext}"
		# Check existence on remote (guarded)
		set +e
		openssl_mismatch_guard ssh "$remote_host" "test -s '$rp'"
		local rc=$?
		set -e
		if ((rc == 0)); then
			log "Found remote: ${remote_host}:$rp"
			base="$(basename "$rp")"
			local_copy="${work}/${base}"
			# Copy (guarded)
			set +e
			openssl_mismatch_guard scp "$remote_host:$rp" "$local_copy"
			rc=$?
			set -e
			if ((rc == 0)); then
				case "$local_copy" in
				*.tar.gz) emit_from_tar "$local_copy" ;;
				*) emit_from_one_file "$local_copy" ;;
				esac
				return 0
			elif ((rc == 100)); then
				# OpenSSL mismatch during scp → warn already printed; skip remote
				return 1
			else
				# other scp error: keep searching other patterns/dirs
				continue
			fi
		elif ((rc == 100)); then
			# OpenSSL mismatch on ssh; skip remote entirely for this dir
			return 1
		fi
	done
	return 1
}

finalize_and_check() {
	if [[ -s "$norm" ]]; then
		log "OK -> $norm"
		printf '%s\n' "$norm"
		return 0
	fi
	log "❌ Normalized output missing or empty: $norm"
	return 22
}

# -------- NCBI fetch (prefetch → vdb-validate → fasterq-dump), guarded
fetch_from_ncbi() {
	log "Fetching from NCBI: $long_sra"
	if [[ "$dry_run" == "true" ]]; then
		log "(dry-run) prefetch '$long_sra' --quiet --max-size u"
		log "(dry-run) vdb-validate '$long_sra' --quiet 2>/dev/null"
		log "(dry-run) fasterq-dump '$long_sra' --quiet 2>/dev/null"
		return 0
	fi

	local rc

	# prefetch
	set +e
	openssl_mismatch_guard prefetch "$long_sra" --quiet --max-size u
	rc=$?
	set -e
	if ((rc == 100)); then
		# OpenSSL mismatch → skip fetch path
		return 20
	elif ((rc != 0)); then
		return 21
	fi

	# vdb-validate (non-fatal; guard for mismatch)
	set +e
	openssl_mismatch_guard vdb-validate "$long_sra" --quiet
	rc=$?
	set -e
	# rc==100 -> mismatch warned; continue anyway
	# rc!=0 and !=100 -> validation error, but keep going (as before)

	# fasterq-dump
	set +e
	openssl_mismatch_guard fasterq-dump "$long_sra" --quiet
	rc=$?
	set -e
	if ((rc == 100)); then
		return 20
	elif ((rc != 0)); then
		return 21
	fi

	# Collect outputs from CWD (SRA tools default)
	mapfile -t outs < <(find . -maxdepth 1 -type f \( \
		-name "${long_sra}.fastq" -o \
		-name "${long_sra}.fq" -o \
		-name "${long_sra}_1.fastq" -o -name "${long_sra}_2.fastq" -o \
		-name "${long_sra}_1.fq" -o -name "${long_sra}_2.fq" \
		\) | sed 's#^\./##' | sort -V)

	if ((${#outs[@]} == 0)); then
		log "No FASTQ produced by fasterq-dump."
		return 20
	fi

	concat_to_norm "${outs[@]}"
}

# ==================== MAIN FLOW ====================

# 1) PWD first
if try_exact_in_dir "$PWD"; then
	finalize_and_check
	exit $?
fi

# 2) media dirs (if provided)
for d in "${media_dir:-}" "${media1_dir:-}" "${media2_dir:-}"; do
	[[ -n "$d" ]] || continue
	if try_exact_in_dir "$d"; then
		finalize_and_check
		exit $?
	fi
done

# 3) remote (same dir order)
if [[ -n "$remote_host" ]]; then
	for d in "${media_dir:-}" "${media1_dir:-}" "${media2_dir:-}"; do
		[[ -n "$d" ]] || continue
		if try_exact_remote_dir "$d"; then
			finalize_and_check
			exit $?
		fi
	done
fi

# 4) fetch from NCBI (guarded)
# Run fetch in a temporary work subdir to avoid cluttering PWD
(
	cd "$work"
	fetch_from_ncbi
)
rc=$?
if ((rc == 0)); then
	finalize_and_check
	exit $?
elif ((rc == 20)); then
	log "Fetch skipped or produced no files (likely OpenSSL mismatch)."
	exit 21
else
	log "❌ NCBI fetch failed."
	exit 21
fi
