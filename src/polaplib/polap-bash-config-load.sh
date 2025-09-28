#!/usr/bin/env bash
# Load a config as PCFG_* variables using preset + config-dir.
# Required:
#   PCFG_PRESET=<name>
# Optional:
#   PCFG_CONFIG_DIR=<dir>   (default: $HOME/.polap/profiles)
#   PCFG_PREFIX=<prefix>    (default: PCFG)

set -euo pipefail

: "${PCFG_PRESET:?Set PCFG_PRESET to the YAML config before sourcing this file}"
: "${PCFG_CONFIG_DIR:=$HOME/.polap/profiles}"
: "${PCFG_PREFIX:=PCFG}"

# call python loader to emit shell assignments, then eval them
tmp_env="/tmp/.pcfg.$$.$RANDOM.env"
python3 "${_POLAPLIB_DIR}/polap-py-config-load.py" \
	--preset "$PCFG_PRESET" \
	--config-dir "$PCFG_CONFIG_DIR" \
	--format env \
	--prefix "$PCFG_PREFIX" >"$tmp_env"

# shellcheck disable=SC1090
source "$tmp_env"
rm -f "$tmp_env"
