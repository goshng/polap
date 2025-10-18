#!/usr/bin/env bash
################################################################################
# polap-bash-run-data-short.sh
# Version: v0.3.0
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Prepare paired short reads under <sfolder>/tmp/ as:
#   - s_1.fq.gz, s_2.fq.gz           (default)
#   - or s_1.fq,   s_2.fq            (with --as-fq / --as-fastq)
#
# Policy
#   • If the FASTQ exists on *this* machine ➜ create a symlink (default gz mode).
#     - In --as-fq mode, we *decompress into a local copy* (keep /media gz intact).
#   • Else if it exists on a remote host ➜ scp to tmp (gz mode keeps gz; --as-fq decompresses).
#   • Else ➜ download from NCBI via `polap-ncbitools fetch sra <SRR>`:
#       - default gz mode: keep as .gz if downloader returns .gz, otherwise compress
#       - --as-fq mode: save as .fq
#
# Notes
#   • Only symlink when the source is local *and* the requested suffix matches (gz mode).
#   • When source is local but only *.fastq (not .gz) exists and we are in gz mode:
#       - we will symlink to *.fastq as s_1.fq (and warn), to avoid copying/compressing local media.
#   • Records SRA ID to <sfolder>/tmp/s.sra.txt
################################################################################
set -euo pipefail
IFS=$'\n\t'

usage() {
	cat <<'EOF'
Usage:
  polap-bash-run-data-short.sh -s SRRxxxxxx -o <sfolder>
                               [--as-fq | --as-fastq]
                               [--threads N]
                               [--remote user@host]
                               [--redo] [--cleanup]
                               [-h|--help]

Inputs:
  -s, --short-sra  SRA accession for the short reads (paired), e.g., SRR10250248
  -o, --outdir     Species folder (we write to <outdir>/tmp/)
  --remote         Remote host alias (default: $ _local_host if set)
  --threads N      Threads for pigz/gzip (default: 4)

Outputs (default gz mode):
  <outdir>/tmp/s_1.fq.gz
  <outdir>/tmp/s_2.fq.gz (if mate 2 exists)

Outputs (--as-fq / --as-fastq):
  <outdir>/tmp/s_1.fq
  <outdir>/tmp/s_2.fq (if mate 2 exists)

Behavior:
  1) LOCAL: if /media or CWD has SRRxxxxxx_1.fastq.gz ➜ symlink to s_1.fq.gz
     (In --as-fq: decompress from gz into s_1.fq; keep /media gz intact.)
  2) REMOTE: if remote has SRRxxxxxx_1.fastq.gz ➜ scp; (in --as-fq: gunzip to .fq)
  3) NCBI: polap-ncbitools fetch sra SRRxxxxxx (into tmp) then normalize names.

Misc:
  --redo     Re-create outputs even if they already exist
  --cleanup  Remove any transient files fetched into tmp (e.g., *_1.fastq.gz.tmp)
EOF
}

# ---- defaults / environment ---------------------------------------------------
SHORT_SRA=""
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
	-s | --short-sra)
		SHORT_SRA="${2:?}"
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

[[ -n "$SHORT_SRA" ]] || {
	echo "[ERROR] --short-sra required" >&2
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
show() { printf '[short] %s\n' "$*" >&2; }

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
	scp -q "$REMOTE":"$rpath" "$lpath"
}

# choose outputs based on mode
target_paths() {
	local mate="$1" # 1 or 2
	if ((AS_FQ)); then
		printf '%s' "$TMP/s_${mate}.fq"
	else
		printf '%s' "$TMP/s_${mate}.fq.gz"
	fi
}

# verify existence unless REDO
check_skip() {
	local t1="$1"
	if ((REDO == 0)) && [[ -s "$t1" ]]; then
		show "exists: $(basename "$t1") (use --redo to rebuild)"
		return 0
	fi
	return 1
}

# decompress/copy pipeline writers
write_from_gz() { # src.gz -> dst(.fq or .fq.gz)
	local src="$1" dst="$2"
	if [[ "$dst" == *.gz ]]; then
		# just link/copy gz
		ln -sf "$src" "$dst"
	else
		# write uncompressed copy
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
		ln -sf "$src" "$dst"
	fi
}

# normalize downloaded names to our s_1/s_2
normalize_downloaded_into() {
	local base="$1" mate="$2" dst="$3"
	# search temp files under TMP for either gz or plain
	local cand_gz="$TMP/${base}_${mate}.fastq.gz"
	local cand_fq="$TMP/${base}_${mate}.fastq"
	if [[ -s "$cand_gz" ]]; then
		write_from_gz "$cand_gz" "$dst"
		((CLEANUP)) && rm -f -- "$cand_gz"
	elif [[ -s "$cand_fq" ]]; then
		write_from_fastq "$cand_fq" "$dst"
		((CLEANUP)) && rm -f -- "$cand_fq"
	else
		return 1
	fi
}

# ---- main: per-mate processing -----------------------------------------------
process_mate() {
	local mate="$1" # "1" or "2"
	local tdst
	tdst="$(target_paths "$mate")"
	check_skip "$tdst" && return 0

	local base="${SHORT_SRA}_${mate}"
	local local_gz local_fq

	# 1) LOCAL: try PWD then media dirs
	local_gz="$(first_existing "./${base}.fastq.gz" \
		"${_media_dir}/${base}.fastq.gz" \
		"${_media1_dir}/${base}.fastq.gz" \
		"${_media2_dir}/${base}.fastq.gz")" || true
	local_fq="$(first_existing "./${base}.fastq" \
		"${_media_dir}/${base}.fastq" \
		"${_media1_dir}/${base}.fastq" \
		"${_media2_dir}/${base}.fastq")" || true

	if [[ -n "$local_gz" ]]; then
		if ((AS_FQ)); then
			# local gz -> make an uncompressed copy
			show "local -> $(basename "$tdst") (gunzip copy)"
			pigz_dc <"$local_gz" >"$tdst"
		else
			# symlink to gz
			show "local -> symlink $(basename "$tdst")"
			ln -sf "$local_gz" "$tdst"
		fi
		return 0
	fi

	if [[ -n "$local_fq" ]]; then
		if ((AS_FQ)); then
			show "local -> symlink $(basename "$tdst")"
			ln -sf "$local_fq" "$tdst"
		else
			# prefer symlink: extension mismatch (warn)
			show "WARN: only *.fastq found locally; linking as uncompressed s_${mate}.fq"
			ln -sf "$local_fq" "$TMP/s_${mate}.fq"
			# if caller expected .gz, leave both s_${mate}.fq and (missing) .gz; or create .gz?
			# We keep policy: avoid copying/compressing local media.
		fi
		return 0
	fi

	# 2) REMOTE: try media dirs if remote specified/available
	if [[ -n "$REMOTE" ]]; then
		local r
		for r in "${_media_dir}/${base}.fastq.gz" \
			"${_media1_dir}/${base}.fastq.gz" \
			"${_media2_dir}/${base}.fastq.gz"; do
			if remote_has "$r"; then
				local tmpgz="$TMP/.${base}.fastq.gz.tmp"
				pull_remote "$r" "$tmpgz"
				if ((AS_FQ)); then
					show "gunzip -> $(basename "$tdst")"
					pigz_dc <"$tmpgz" >"$tdst"
					((CLEANUP)) && rm -f -- "$tmpgz"
				else
					mv -f "$tmpgz" "$tdst"
				fi
				return 0
			fi
		done
		for r in "${_media_dir}/${base}.fastq" \
			"${_media1_dir}/${base}.fastq" \
			"${_media2_dir}/${base}.fastq"; do
			if remote_has "$r"; then
				local tmpfq="$TMP/.${base}.fastq.tmp"
				pull_remote "$r" "$tmpfq"
				if ((AS_FQ)); then
					mv -f "$tmpfq" "$tdst"
				else
					show "gzip -> $(basename "$tdst")"
					pigz_c 6 <"$tmpfq" >"$tdst"
					((CLEANUP)) && rm -f -- "$tmpfq"
				fi
				return 0
			fi
		done
	fi

	# 3) NCBI: fallback download
	show "NCBI fetch: ${SHORT_SRA} (mate ${mate})"
	# ensure conda env if you use it (optional hook)
	if command -v _polap_lib_conda-ensure_conda_env >/dev/null 2>&1; then
		_polap_lib_conda-ensure_conda_env polap-ncbitools || true
	fi
	(
		cd "$TMP"
		# This wrapper is assumed to create ${SHORT_SRA}_1.fastq[.gz] and _2.* if paired.
		if command -v polap-ncbitools >/dev/null 2>&1; then
			polap-ncbitools fetch sra "${SHORT_SRA}"
		else
			# fallback: try fasterq-dump if present (no prefetch here)
			if command -v fasterq-dump >/dev/null 2>&1; then
				fasterq-dump -e "$THREADS" -p "${SHORT_SRA}"
				# gzip outputs to save space if default gz mode requested
				if ((AS_FQ == 0)); then
					[[ -s "${SHORT_SRA}_1.fastq" ]] && pigz -p "$THREADS" "${SHORT_SRA}_1.fastq" || true
					[[ -s "${SHORT_SRA}_2.fastq" ]] && pigz -p "$THREADS" "${SHORT_SRA}_2.fastq" || true
				fi
			else
				echo "[ERROR] neither polap-ncbitools nor fasterq-dump found" >&2
				exit 3
			fi
		fi
	)

	# normalize into tdst
	if normalize_downloaded_into "${SHORT_SRA}" "${mate}" "$tdst"; then
		:
	else
		echo "[ERROR] could not locate downloaded mate ${mate} for ${SHORT_SRA}" >&2
		exit 4
	fi
	return 0
}

# ---- run mates ---------------------------------------------------------------
# record SRA
printf '%s\n' "$SHORT_SRA" >"$TMP/s.sra.txt"

process_mate "1" || exit $?
# mate 2 is optional (some runs are single-end)
if process_mate "2"; then
	:
else
	show "mate 2 not found; continuing with single mate"
fi

show "done: $(basename "$(target_paths 1)") $([[ -e "$(target_paths 2)" ]] && printf ' %s' "$(basename "$(target_paths 2)")" || true)"
exit 0
