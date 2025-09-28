#!/usr/bin/env bash
# polap-bash-bolap-env-check.sh
# Validate bolap host/media config and optionally mirror from a remote host.
# Usage: polap-bash-bolap-env-check.sh --host thorne --media /media/h2/sra --mirror github --dry-run
set -euo pipefail

HOST="thorne"
MEDIA="/media/h2/sra"
MIRROR_DIR="github"
DRY=0
RELEASE="${_POLAP_RELEASE:-0}"

while [[ $# -gt 0 ]]; do
	case "$1" in
	--host)
		HOST="$2"
		shift 2
		;;
	--media)
		MEDIA="$2"
		shift 2
		;;
	--mirror)
		MIRROR_DIR="$2"
		shift 2
		;;
	--dry-run)
		DRY=1
		shift
		;;
	*)
		echo "Unknown arg: $1" >&2
		exit 2
		;;
	esac
done

# echo "[bolap/env] hostname=$(hostname)"
if [[ ! -d "$MEDIA" ]]; then
	echo "[bolap/env] WARN: media dir not present: $MEDIA" >&2
fi

if [[ "${RELEASE}" == "0" ]]; then
	if [[ "$HOST" != "$(hostname)" ]]; then
		# mirror ${PWD}/github from remote ${HOST}:${PWD}/github
		if [[ $DRY -eq 1 ]]; then
			echo "[bolap/env] DRY-RUN: rsync -aq ${HOST}:${PWD}/${MIRROR_DIR}/ ${PWD}/${MIRROR_DIR}/"
		else
			rsync -aq "${HOST}:${PWD}/${MIRROR_DIR}/" "${PWD}/${MIRROR_DIR}/"
		fi
	fi
else
	echo "[bolap/env] RELEASE=1: skipping rsync mirror"
fi
