#!/usr/bin/env bash
set -euo pipefail

# polap-bash-template.sh
# Template: load config into environment via Python, then use PCFG_* vars.

# Defaults: allow PCFG_CONFIG_DIR override; else ~/.polap/profiles
PCFG_CONFIG_DIR="${PCFG_CONFIG_DIR:-$HOME/.polap/profiles}"
POLAPLIB_DIR="${_POLAPLIB_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)}"

# Parse addressing
PRESET=""
PATH_YAML=""

ARGS=("$@")
i=0
while [[ $i -lt ${#ARGS[@]} ]]; do
	case "${ARGS[$i]}" in
	--preset)
		PRESET="${ARGS[$((i + 1))]:-}"
		i=$((i + 2))
		;;
	--config-dir)
		PCFG_CONFIG_DIR="${ARGS[$((i + 1))]:-}"
		i=$((i + 2))
		;;
	--path)
		PATH_YAML="${ARGS[$((i + 1))]:-}"
		i=$((i + 2))
		;;
	-h | --help)
		cat <<EOF
Usage:
  $(basename "$0") [--path FILE.yaml | --config-dir DIR --preset NAME]

Notes:
  - If --path absent, defaults to --config-dir \$HOME/.polap/profiles
  - Exposes PCFG_* environment variables for downstream commands
EOF
		exit 0
		;;
	*)
		echo "[error] unknown option: ${ARGS[$i]}" >&2
		exit 2
		;;
	esac
done

# Build call to polap-py-config-load.py to emit PCFG_* env lines
LOAD_PY="${POLAPLIB_DIR}/polap-py-config-load.py"
[[ -f "$LOAD_PY" ]] || {
	echo "[error] not found: $LOAD_PY" >&2
	exit 127
}

if [[ -n "$PATH_YAML" ]]; then
	ENV_LINES=$(python3 "$LOAD_PY" --path "$PATH_YAML" --format env --prefix "PCFG_")
else
	if [[ -z "$PRESET" ]]; then
		echo "[error] need --preset when --path is not given" >&2
		exit 2
	fi
	ENV_LINES=$(python3 "$LOAD_PY" --config-dir "$PCFG_CONFIG_DIR" --preset "$PRESET" --format env --prefix "PCFG_")
fi

# Export PCFG_* to the current shell
# shellcheck disable=SC1090
eval "$ENV_LINES"

# Demonstrate usage
echo "[config] preset=${PCFG_PRESET:-}"
echo "         reads=${PCFG_READS:-}"
echo "         threads=${PCFG_THREADS:-}"
echo "         k=${PCFG_K:-} s=${PCFG_S:-} hpc=${PCFG_HPC:-}"

# Now your Bash logic can use PCFG_* vars, e.g.:
# minimap2 -t "${PCFG_THREADS:-8}" ...
