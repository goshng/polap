#!/usr/bin/env bash
# polap-bash-patch-polish.sh
# Version: v0.1.0
# SPDX-License-Identifier: GPL-3.0-or-later
set -euo pipefail

F="${1:-polap-bash-polish.sh}"

[[ -f "$F" ]] || {
	echo "[ERROR] file not found: $F" >&2
	exit 2
}

# Backup
cp -p "$F" "${F}.bak"

# Replace any line that consists solely of optional spaces + ';*' with ';;'
awk '
  /^[[:space:]]*;\*[[:space:]]*$/ { sub(/;\*/, ";;"); patched=1 }
  { print }
  END {
    if (patched != 1) {
      # not fatal; maybe already fixed or multiple; we don’t fail here
      # (If you expect exactly one, change to: if (patched!=1) exit 3)
    }
  }
' "$F" >"${F}.tmp" && mv "${F}.tmp" "$F"

# Show the neighborhood where the error likely was (around your line 454)
nl -ba "$F" | sed -n '440,470p'

# Validate
echo "[INFO] bash -n check…"
bash -n "$F" && echo "[OK] syntax clean"

if command -v shellcheck >/dev/null 2>&1; then
	echo "[INFO] shellcheck…"
	shellcheck -x "$F" || true
fi

# Show any remaining literal ';*' just in case
if grep -n ';\*' "$F" >/dev/null; then
	echo "[WARN] Remaining ';*' patterns:"
	grep -n ';\*' "$F" || true
fi
