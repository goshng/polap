#!/usr/bin/env bash
set -euo pipefail

# Usage: polap-bash-generate-species-codes.sh input.txt [output.txt]
# Reads species names with underscores and prepends a 2-letter + 2-digit code.
# Example:
#   Aristolochia_californica           -> Ac01 Aristolochia_californica
#   Brassica_napus                     -> Bn02 Brassica_napus
#   Cucumis_sativus_var_hardwickii     -> Cs03 Cucumis_sativus_var_hardwickii

if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
	echo "Usage: $0 input.txt [output.txt]" >&2
	exit 1
fi

input="$1"
output="${2:-/dev/stdout}"

awk '
BEGIN {
    FS = "_";              # split on underscore
    n  = 0;                # line counter for numbering
}
NF > 0 {                   # skip completely empty lines
    n++
    first  = substr($1, 1, 1)                    # first letter of the first word
    second = (NF >= 2 ? substr($2, 1, 1) : "x")  # first letter after first underscore (or x if missing)

    code = sprintf("%s%s%02d", first, second, n) # 2 letters + 2-digit zero-padded number
    print code, $0
}
' "$input" >"$output"
