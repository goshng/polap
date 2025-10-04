#!/usr/bin/env bash
# dialog_splice_copy_any.sh — Copy either a single file or an entire directory.
# Dependencies: dialog, pv; Optional: splice (from moreutils)
set -euo pipefail

TITLE="Splice Copy"
SRC=""
DST=""
USE_SPLICE=0

require_cmd() { command -v "$1" >/dev/null 2>&1 || {
	echo "Need '$1'"
	exit 1
}; }
d_msg() { dialog --title "$TITLE" --msgbox "$1" 10 70; }
d_yesno() { dialog --title "$TITLE" --yesno "$1" 10 70; }
d_input() { dialog --title "$TITLE" --inputbox "$1" 10 80 2>&1 >/dev/tty; }
d_fsel() { dialog --title "$TITLE" --fselect "$1" 20 80 2>&1 >/dev/tty; }
d_dsel() { dialog --title "$TITLE" --dselect "$1" 20 80 2>&1 >/dev/tty; }

cleanup() { [[ -p "${GAUGE_FIFO:-}" ]] && rm -f "$GAUGE_FIFO"; }
trap cleanup EXIT

emit_progress() {
	local p="$1"
	shift
	local msg="$*"
	{
		echo "XXX"
		echo "$p"
		printf "%s\n" "$msg"
		echo "XXX"
	} >"$GAUGE_FIFO"
}

toggle_splice() {
	if command -v splice >/dev/null 2>&1; then
		USE_SPLICE=$((1 - USE_SPLICE))
	else
		USE_SPLICE=0
		d_msg "'splice' not found (install with: sudo apt install moreutils)"
	fi
}

pick_src() {
	local hint="${SRC:-$PWD/}"
	SRC="$(d_fsel "$hint")" || return
	[[ -n "$SRC" ]] || return
	SRC="${SRC%/}"
	[[ -e "$SRC" ]] || d_msg "No such file or directory: $SRC"
}

pick_dst() {
	local hint="${DST:-$PWD/}"
	DST="$(d_input "Enter destination path or directory:\n\nHint: $hint")" || return
	[[ -n "$DST" ]] || return
	if [[ ! -e "$DST" ]]; then
		if d_yesno "Create path or directory?\n$DST"; then
			mkdir -p -- "$DST"
		fi
	fi
}

show_help() {
	d_msg "How it works:\n
• Choose SOURCE (file or directory) and DESTINATION.\n
• If SOURCE is a directory, it copies recursively.\n
• If SOURCE is a single file, it copies that file only.\n
• pv shows total progress; dialog displays a gauge.\n
• splice (optional) performs zero-copy streaming.\n
Install:\n  sudo apt install dialog pv moreutils"
}

copy_file() {
	local src="$1" dst="$2" sofar="$3" total="$4"
	mkdir -p "$(dirname "$dst")"

	local sz
	sz=$(stat -c '%s' "$src")
	if [[ "$USE_SPLICE" -eq 1 ]]; then
		pv -n "$src" 2> >(while read -r p; do
			local cur=$((sofar + (p * sz) / 100))
			local pct=$(((cur * 100) / (total == 0 ? 1 : total)))
			emit_progress "$pct" "Copying: $(basename "$src")"
		done) | splice >"$dst"
	else
		pv -n "$src" 2> >(while read -r p; do
			local cur=$((sofar + (p * sz) / 100))
			local pct=$(((cur * 100) / (total == 0 ? 1 : total)))
			emit_progress "$pct" "Copying: $(basename "$src")"
		done) >"$dst"
	fi
}

start_copy() {
	require_cmd dialog
	require_cmd pv
	command -v splice >/dev/null 2>&1 || true

	if [[ -z "$SRC" || -z "$DST" ]]; then
		d_msg "Please set both SOURCE and DESTINATION."
		return
	fi
	[[ -e "$SRC" ]] || {
		d_msg "Source does not exist."
		return
	}

	# Determine mode
	local mode="file"
	[[ -d "$SRC" ]] && mode="dir"

	# Prepare file list
	local -a FILES SIZES
	if [[ "$mode" == "file" ]]; then
		FILES=("$SRC")
		SIZES=($(stat -c '%s' "$SRC"))
	else
		mapfile -d '' FILES < <(cd "$SRC" && find . -type f -print0)
		for f in "${FILES[@]}"; do SIZES+=($(stat -c '%s' "$SRC/${f#./}")); done
	fi

	local total=0
	for s in "${SIZES[@]}"; do total=$((total + s)); done
	((total > 0)) || {
		d_msg "No files to copy."
		return
	}

	GAUGE_FIFO="$(mktemp -u)"
	mkfifo "$GAUGE_FIFO"
	(dialog --title "$TITLE" --gauge "Preparing…" 10 80 0 <"$GAUGE_FIFO") &
	local gauge_pid=$!
	emit_progress 0 "Starting copy…"

	local sofar=0
	for i in "${!FILES[@]}"; do
		local rel="${FILES[i]}"
		local abs_src abs_dst
		if [[ "$mode" == "file" ]]; then
			abs_src="$SRC"
			if [[ -d "$DST" ]]; then abs_dst="$DST/$(basename "$SRC")"; else abs_dst="$DST"; fi
		else
			abs_src="$SRC/${rel#./}"
			abs_dst="$DST/$(basename "$SRC")/${rel#./}"
		fi
		copy_file "$abs_src" "$abs_dst" "$sofar" "$total"
		sofar=$((sofar + SIZES[i]))
		emit_progress $((sofar * 100 / total)) "Copied: $(basename "$abs_src")"
	done

	emit_progress 100 "Done."
	sleep 0.3
	kill "$gauge_pid" 2>/dev/null || true
	rm -f "$GAUGE_FIFO"

	d_msg "Copy complete.\n\nFrom: $SRC\nTo:   $DST\nTotal bytes: $total"
}

main_menu() {
	local sstate=$([ "$USE_SPLICE" -eq 1 ] && echo ON || echo OFF)
	local tmp
	tmp="$(mktemp)"
	trap 'rm -f "$tmp"' RETURN
	dialog --title "$TITLE" --menu "Splice Copy Tool" 15 75 7 \
		1 "Pick SOURCE (file or dir)   [${SRC:-not set}]" \
		2 "Pick DESTINATION path       [${DST:-not set}]" \
		3 "Toggle splice (zero-copy):  [$sstate]" \
		4 "Start copy" \
		5 "Help" \
		6 "Quit" 2>"$tmp" || return
	case "$(cat "$tmp")" in
	1) pick_src ;;
	2) pick_dst ;;
	3) toggle_splice ;;
	4) start_copy ;;
	5) show_help ;;
	6 | "") exit 0 ;;
	esac
}

require_cmd dialog
require_cmd pv
command -v splice >/dev/null 2>&1 && USE_SPLICE=1 || USE_SPLICE=0

while true; do
	main_menu || exit 0
done
