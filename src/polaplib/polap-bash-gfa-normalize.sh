#!/usr/bin/env bash
# polap-bash-gfa-normalize.sh
# Version: v0.7.0
# Normalize a Flye GFA for gfatk:
#   - Ensure GFA1 header
#   - Renumber S IDs to 1..N (and emit idmap)
#   - On S: ensure ll:f (from dp:i/f if present; else 1.0)
#   - On L: coerce CIGAR to <int>M (sum M; '*' → '0M'), ensure ec:i (from endpoints' ll)
#   - Ensure reciprocal L edges exist (with same ec:i)
# Output: normalized GFA + idmap
set -euo pipefail
IFS=$'\n\t'

usage() {
	cat <<'USAGE'
Usage:
  polap-bash-gfa-normalize.sh --gfa assembly_graph.gfa --out norm.gfa --idmap idmap.tsv \
      [--ec-mode mean|min|max|const] [--ec-const 1] [--ec-scale 1.0] [--ec-round round|ceil|floor]
USAGE
}

gfa="" out="" idmap=""
ec_mode="mean" ec_const="1" ec_scale="1.0" ec_round="round"

while (($#)); do
	case "$1" in
	--gfa)
		gfa="$2"
		shift 2
		;;
	--out)
		out="$2"
		shift 2
		;;
	--idmap)
		idmap="$2"
		shift 2
		;;
	--ec-mode)
		ec_mode="$2"
		shift 2
		;;
	--ec-const)
		ec_const="$2"
		shift 2
		;;
	--ec-scale)
		ec_scale="$2"
		shift 2
		;;
	--ec-round)
		ec_round="$2"
		shift 2
		;;
	-h | --help)
		usage
		exit 0
		;;
	*)
		echo "ERR: unknown arg $1" >&2
		exit 2
		;;
	esac
done

[[ -s "$gfa" ]] || {
	echo "ERR: --gfa missing/empty" >&2
	exit 2
}
[[ -n "$out" && -n "$idmap" ]] || {
	echo "ERR: --out and --idmap required" >&2
	exit 2
}

command -v python3 >/dev/null 2>&1 || {
	echo "ERR: need python3" >&2
	exit 127
}

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
py="${here}/scripts/polap_py_gfa_normalize_ec_bidirectional.py"

python3 "$py" \
	--gfa "$gfa" \
	--out "$out" \
	--idmap "$idmap" \
	--ec-mode "$ec_mode" \
	--ec-const "$ec_const" \
	--ec-scale "$ec_scale" \
	--ec-round "$ec_round"

echo "[OK] norm: $out"
echo "[OK] idmap: $idmap"
