#!/usr/bin/env bash
################################################################################
# polap-bash-run-data-long.sh
# Version: v0.3.0
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Prepare a long-read FASTQ under <sfolder>/tmp/ as:
#   - l.fq.gz              (default)
#   - or l.fq              (with --as-fq / --as-fastq)
#
# Policy (mirrors short-read script):
#   • If the FASTQ exists on *this* machine ➜ create a symlink (gz mode).
#     - In --as-fq: decompress into a local l.fq (keep /media gz intact).
#   • Else if it exists on a remote ➜ scp to tmp (gz mode keeps gz; --as-fq gunzips).
#   • Else ➜ download from NCBI via `polap-ncbitools fetch sra <SRR>` (fallback: fasterq-dump).
#       - default gz mode: if .fastq, compress to .fq.gz; if .gz given, keep .gz.
#       - --as-fq: save as plain .fq
#
# Notes
#   • Only symlink when source is local *and* the requested suffix matches (gz mode).
#   • When only *.fastq (not .gz) exists locally and gz mode is requested:
#       - we avoid recompressing media; we link as l.fq and warn.
#   • Records SRA ID to <sfolder>/tmp/l.sra.txt
################################################################################
set -euo pipefail
IFS=$'\n\t'

usage() {
	cat <<'EOF'
Usage:
  polap-bash-run-data-long.sh -l SRRxxxxxx -o <sfolder>
                              [--as-fq | --as-fastq]
                              [--threads N]
                              [--remote user@host]
                              [--redo] [--cleanup]
                              [-h|--help]

Inputs:
  -l, --long-sra   SRA accession for the long reads (single), e.g., SRR10190639
  -o, --outdir     Species folder (we write to <outdir>/tmp/)

Options:
  --as-fq | --as-fastq  Write plain l.fq (decompress if needed). Default: l.fq.gz
  --threads N           Threads for pigz/gzip (default: 4)
  --remote HOST         Remote host alias (default: $_local_host if set)
  --redo                Re-create outputs even if they already exist
  --cleanup             Remove transient files fetched into tmp (e.g., .tmp files)
  -h, --help            Show this help

Outputs:
  <outdir>/tmp/l.fq.gz  (default)   or   <outdir>/tmp/l.fq  (with --as-fq)
  <outdir>/tmp/l.sra.txt  (contains the SRR)
EOF
}

# ---- defaults / environment ---------------------------------------------------
LONG_SRA=""
OUTDIR=""
AS_FQ=0
THREADS=4
REMOTE="${_local_host:-}" # optional override
REDO=0
CLEANUP=0

# site-local media hints (optional)
: "${_media_dir:=/media/h2/sra}"
: "${_media1_dir:=/media/h1/sra}"
: "${_media2_dir:=/media/h2/sra}"

# ---- args ---------------------------------------------------------------------
while (($#)); do
	case "$1" in
	-l | --long-sra)
		LONG_SRA="${2:?}"
		shift 2
		;;
	-o | --outdir)
		OUTDIR="${2:?}"
		shift 2
		;;
	--as-fq | --as-fastq)
		AS_FQ=1
		shift
		;;
	--threads)
		THREADS="${2:?}"
		shift 2
		;;
	--remote)
		REMOTE="${2:?}"
		shift 2
		;;
	--redo)
		REDO=1
		shift
		;;
	--cleanup)
		CLEANUP=1
		shift
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

[[ -n "$LONG_SRA" ]] || {
	echo "[ERROR] --long-sra required" >&2
	usage
	exit 2
}
[[ -n "$OUTDIR" ]] || {
	echo "[ERROR] --outdir required" >&2
	usage
	exit 2
}

TMP="$OUTDIR/tmp"
mkdir -p "$TMP"

# ---- helpers ------------------------------------------------------------------
have() { command -v "$1" >/dev/null 2>&1; }
pigz_dc() { if have pigz; then pigz -dc -p "$THREADS"; else gzip -dc; fi; }
pigz_c() { if have pigz; then pigz -c -p "$THREADS" -"${1:-6}"; else gzip -c -"${1:-6}"; fi; }
show() { printf '[long] %s\n' "$*" >&2; }

# locate first existing path among candidates
first_existing() {
	local p
	for p in "$@"; do
		[[ -s "$p" ]] && {
			printf '%s' "$p"
			return 0
		}
	done
	return 1
}

# remote test
remote_has() {
	local path="$1"
	[[ -n "$REMOTE" ]] || return 1
	ssh -o BatchMode=yes -o ConnectTimeout=5 "$REMOTE" "test -s '$path'" 2>/dev/null
}

# scp helper
pull_remote() {
	local rpath="$1" lpath="$2"
	mkdir -p "$(dirname "$lpath")"
	show "scp: $REMOTE:$rpath -> $lpath"
	scp "$REMOTE":"$rpath" "$lpath"
}

# choose output based on mode
target_path() {
	if ((AS_FQ)); then
		printf '%s' "$TMP/l.fq"
	else
		printf '%s' "$TMP/l.fq.gz"
	fi
}

# verify existence unless REDO
check_skip() {
	local t="$1"
	if ((REDO == 0)) && [[ -s "$t" ]]; then
		show "exists: $(basename "$t") (use --redo to rebuild)"
		return 0
	fi
	return 1
}

# write from gz/plain sources
write_from_gz() { # src.gz -> dst(.fq or .fq.gz)
	local src="$1" dst="$2"
	if [[ "$dst" == *.gz ]]; then
		ln -sf "$(basename $src)" "$dst"
	else
		show "gunzip -> $(basename "$dst")"
		pigz_dc <"$src" >"$dst"
	fi
}

write_from_fastq() { # src.fastq -> dst(.fq or .fq.gz)
	local src="$1" dst="$2"
	if [[ "$dst" == *.gz ]]; then
		show "gzip -> $(basename "$dst")"
		pigz_c 6 <"$src" >"$dst"
	else
		ln -sf "$(basename $src)" "$dst"
	fi
}

# polap-bash-normalize-downloaded.sh
# Version: v0.1.0
# GPL-3.0+
normalize_downloaded_into() {
	local base="$1" dst="$2"
	local cand_gz="$TMP/${base}.fastq.gz"
	local cand_fq="$TMP/${base}.fastq"

	if [[ -s "$cand_gz" ]]; then
		write_from_gz "$cand_gz" "$dst" || return 2
		if ((CLEANUP)); then rm -f -- "$cand_gz"; fi
		return 0
	elif [[ -s "$cand_fq" ]]; then
		write_from_fastq "$cand_fq" "$dst" || return 3
		if ((CLEANUP)); then rm -f -- "$cand_fq"; fi
		return 0
	else
		return 1
	fi
}

# ---- main ---------------------------------------------------------------------
# record SRA
printf '%s\n' "$LONG_SRA" >"$TMP/l.sra.txt"

TDST="$(target_path)"
check_skip "$TDST" && {
	show "done: $(basename "$TDST")"
	exit 0
}

base="${LONG_SRA}"
local_gz="$(first_existing "./${base}.fastq.gz" \
	"${_media_dir}/${base}.fastq.gz" \
	"${_media1_dir}/${base}.fastq.gz" \
	"${_media_dir}/${base}-10x.fastq.gz" \
	"${_media1_dir}/${base}-10x.fastq.gz" \
	"${_media2_dir}/${base}.fastq.gz")" || true
local_fq="$(first_existing "./${base}.fastq" \
	"${_media_dir}/${base}.fastq" \
	"${_media1_dir}/${base}.fastq" \
	"${_media_dir}/${base}-10x.fastq" \
	"${_media1_dir}/${base}-10x.fastq" \
	"${_media2_dir}/${base}.fastq")" || true

# 1) LOCAL
if [[ -n "$local_gz" ]]; then
	if ((AS_FQ)); then
		show "local -> $(basename "$TDST") (gunzip copy)"
		pigz_dc <"$local_gz" >"$TDST"
	else
		show "local -> symlink $(basename "$TDST")"
		ln -sf "$local_gz" "$TDST"
	fi
	show "done: $(basename "$TDST")"
	exit 0
fi

if [[ -n "$local_fq" ]]; then
	if ((AS_FQ)); then
		show "local -> symlink $(basename "$TDST")"
		ln -sf "$local_fq" "$TDST"
	else
		show "WARN: only *.fastq found locally; linking as uncompressed l.fq"
		ln -sf "$local_fq" "$TMP/l.fq"
	fi
	show "done: $(basename "$(target_path)")"
	exit 0
fi

# 2) REMOTE
if [[ -n "$REMOTE" ]]; then
	for r in "${_media_dir}/${base}.fastq.gz" \
		"${_media1_dir}/${base}.fastq.gz" \
		"${_media1_dir}/${base}-10x.fastq.gz" \
		"${_media2_dir}/${base}-10x.fastq.gz" \
		"${_media2_dir}/${base}.fastq.gz"; do
		if remote_has "$r"; then
			tmpgz="$TMP/.${base}.fastq.gz.tmp"
			pull_remote "$r" "$tmpgz"
			if ((AS_FQ)); then
				show "gunzip -> $(basename "$TDST")"
				pigz_dc <"$tmpgz" >"$TDST"
				((CLEANUP)) && rm -f -- "$tmpgz"
			else
				mv -f "$tmpgz" "$TDST"
			fi
			show "done: $(basename "$TDST")"
			exit 0
		fi
	done
	for r in "${_media_dir}/${base}.fastq" \
		"${_media1_dir}/${base}.fastq" \
		"${_media1_dir}/${base}-10x.fastq" \
		"${_media2_dir}/${base}-10x.fastq" \
		"${_media2_dir}/${base}.fastq"; do
		if remote_has "$r"; then
			tmpfq="$TMP/.${base}.fastq.tmp"
			pull_remote "$r" "$tmpfq"
			if ((AS_FQ)); then
				mv -f "$tmpfq" "$TDST"
			else
				show "gzip -> $(basename "$TDST")"
				pigz_c 6 <"$tmpfq" >"$TDST"
				((CLEANUP)) && rm -f -- "$tmpfq"
			fi
			show "done: $(basename "$TDST")"
			exit 0
		fi
	done
fi

# 3) NCBI: fallback download into TMP then normalize
show "NCBI fetch: ${LONG_SRA}"
if command -v _polap_lib_conda-ensure_conda_env >/dev/null 2>&1; then
	_polap_lib_conda-ensure_conda_env polap-ncbitools || true
fi
(
	cd "$TMP"
	if command -v polap-ncbitools >/dev/null 2>&1; then
		polap-ncbitools fetch sra "${LONG_SRA}"
	else
		if command -v fasterq-dump >/dev/null 2>&1; then
			fasterq-dump -e "$THREADS" -p "${LONG_SRA}"
			# compress if gz mode requested
			if ((AS_FQ == 0)); then
				if command -v pigz >/dev/null 2>&1; then
					[[ -s "${LONG_SRA}.fastq" ]] && pigz -p "$THREADS" "${LONG_SRA}.fastq" || true
				else
					[[ -s "${LONG_SRA}.fastq" ]] && gzip "${LONG_SRA}.fastq" || true
				fi
			fi
		else
			echo "[ERROR] neither polap-ncbitools nor fasterq-dump found" >&2
			exit 3
		fi
	fi
)

if normalize_downloaded_into "${LONG_SRA}" "$TDST"; then
	show "done: $(basename "$TDST")"
else
	echo "[ERROR] could not locate downloaded ${LONG_SRA}.fastq[.gz]" >&2
	exit 4
fi

exit 0
