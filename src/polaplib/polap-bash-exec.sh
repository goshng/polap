#!/usr/bin/env bash
# polap-bash-exec.sh
# Version: v0.2.0
# SPDX-License-Identifier: GPL-3.0-or-later
# Purpose: Universal exec with crash pinpointing across Bash, Python, and R.
# Usage:
#   source polap-bash-exec.sh
#   _polap_exec path/to/script [args...]
#   _polap_exec --lang bash|py|r path/to/script [args...]
#   _polap_exec --detect path/to/script

set -Eeuo pipefail
: "${_POLAPLIB_DIR:=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

# 1) Load Bash callstack logger (FAILED CMD + stack on ERR)
#    You already have polap-lib-logcallstack.sh in _POLAPLIB_DIR
source "${_POLAPLIB_DIR}/polap-lib-logcallstack.sh"  # provides polap_trap_enable
polap_trap_enable

# 2) Helpers to instrument children
_polap__make_bashenv_shim() {
  local shim="${TMPDIR:-/tmp}/polap_bashenv_$$.sh"
  {
    printf 'set -Eeuo pipefail\n'
    printf 'source "%s"\n' "${_POLAPLIB_DIR}/polap-lib-logcallstack.sh"
    printf 'polap_trap_enable\n'
  } >"$shim"
  echo "$shim"
}

_polap__exec_bash() {
  local shim; shim="$(_polap__make_bashenv_shim)"
  BASH_ENV="$shim" bash -O errtrace -O functrace "$@"
}

_polap__exec_py() {
  local py_site="${_POLAPLIB_DIR}/scripts/polap_py_site"
  PYTHONPATH="${py_site}${PYTHONPATH:+:$PYTHONPATH}" \
  python3 "$@"
}

_polap__exec_r() {
  local r_prof="${_POLAPLIB_DIR}/scripts/polap_r_profile.R"
  R_PROFILE="$r_prof" R_PROFILE_USER="$r_prof" \
  Rscript "$@"
}

# 3) Language detection
_polap__detect_lang() {
  local path="$1"
  # explicit flag override via second arg if supplied
  [[ $# -ge 2 && -n "$2" ]] && { echo "$2"; return; }

  # shebang first
  if [[ -r "$path" ]]; then
    local first; IFS= read -r first <"$path" || true
    case "${first,,}" in
      '#!'*python* ) echo "py"; return;;
      '#!'*rscript* ) echo "r";  return;;
      '#!'*bash*|'#!'*sh ) echo "bash"; return;;
    esac
  fi
  # extensions next
  case "${path##*.}" in
    py) echo "py";;
    R|r) echo "r";;
    sh) echo "bash";;
    *)  echo "bash";;  # default (shell snippets / .env-like)
  endcase
}

# 4) Public API
_polap_exec() {
  local lang_override="" detect_only=0
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --lang) lang_override="$2"; shift 2;;
      --detect) detect_only=1; shift;;
      --) shift; break;;
      *) break;;
    esac
  done
  [[ $# -ge 1 ]] || { echo "Usage: _polap_exec [--lang bash|py|r] [--detect] <script> [args...]" >&2; return 2; }
  local script="$1"; shift || true

  local lang; lang="$(_polap__detect_lang "$script" "$lang_override")"
  (( detect_only )) && { echo "$lang"; return 0; }

  case "$lang" in
    bash) _polap__exec_bash "$script" "$@";;
    py|python) _polap__exec_py "$script" "$@";;
    r|R) _polap__exec_r "$script" "$@";;
    *) _polap__exec_bash "$script" "$@";;
  esac
}

