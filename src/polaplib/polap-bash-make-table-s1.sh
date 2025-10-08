#!/usr/bin/env bash
# polap-bash-make-table-s1.sh
# Version : v0.1.1
set -euo pipefail
IFS=$'\n\t'

# Base dir
_POLAPLIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export _POLAPLIB_DIR

OUT="md/tableS1-dataset-summary.tsv"
MARKDOWN=0
DEDUP="species"
SELECT="max-bases"
SUMMARY_SCOPES=""
SUMMARY_STAT="median"
MANIFESTS=()
MANIFEST_LIST=""
MANIFEST_GLOB=""
PLATFORM_MAP=""

usage() {
	cat <<EOF
Usage:
  $(basename "$0") --out FILE [--markdown]
                   [--dedup species|species-tier|none]
                   [--select max-bases|latest|first]
                   [--summary-scopes "platform,tier,platform-tier,overall"]
                   [--summary-stat median|mean|both]
                   [--platform-map FILE]
                   --manifest FILE ... | --manifest-list FILE | --manifest-glob "glob"
EOF
}

while (($#)); do
	case "$1" in
	--out)
		OUT="${2:?}"
		shift 2
		;;
	--markdown)
		MARKDOWN=1
		shift
		;;
	--dedup)
		DEDUP="${2:?}"
		shift 2
		;;
	--select)
		SELECT="${2:?}"
		shift 2
		;;
	--summary-scopes)
		SUMMARY_SCOPES="${2:?}"
		shift 2
		;;
	--summary-stat)
		SUMMARY_STAT="${2:?}"
		shift 2
		;;
	--platform-map)
		PLATFORM_MAP="${2:?}"
		shift 2
		;;
	--manifest)
		MANIFESTS+=("${2:?}")
		shift 2
		;;
	--manifest-list)
		MANIFEST_LIST="${2:?}"
		shift 2
		;;
	--manifest-glob)
		MANIFEST_GLOB="${2:?}"
		shift 2
		;;
	-h | --help)
		usage
		exit 0
		;;
	*)
		echo "[ERR] Unknown option: $1" >&2
		usage
		exit 2
		;;
	esac
done

mkdir -p "$(dirname "$OUT")"

# polaplib/polap-bash-make-table-s1.sh (only the list processing section)

# Build a list file for manifests
tmp_list="$(mktemp -t polap-manifests.XXXXXX.txt)"
trap 'rm -f "$tmp_list"' EXIT

# explicit
for m in "${MANIFESTS[@]:-}"; do
	[[ -n "${m:-}" ]] || continue
	printf '%s\n' "$m" >>"$tmp_list"
done
# from list file
if [[ -n "${MANIFEST_LIST}" ]]; then
	[[ -s "$MANIFEST_LIST" ]] || {
		echo "[ERR] manifest-list not found: $MANIFEST_LIST" >&2
		exit 2
	}
	cat "$MANIFEST_LIST" >>"$tmp_list"
fi
# from glob
if [[ -n "${MANIFEST_GLOB}" ]]; then
	shopt -s nullglob
	found=($MANIFEST_GLOB)
	shopt -u nullglob
	((${#found[@]})) || {
		echo "[ERR] glob matched no files: $MANIFEST_GLOB" >&2
		exit 2
	}
	printf '%s\n' "${found[@]}" >>"$tmp_list"
fi

# normalize: trim, drop blanks and comments, dedupe
norm_list="$(mktemp -t polap-manifests.norm.XXXXXX.txt)"
awk 'BEGIN{FS=OFS="\n"} {gsub(/^[ \t]+|[ \t]+$/,""); if($0!="" && $0!~/^#/){print $0}}' "$tmp_list" | sort -u >"$norm_list"

echo "[INFO] Manifests to merge:"
nl -ba "$norm_list" | sed 's/^/[INFO]   /'

# validate
while IFS= read -r m; do
	[[ -s "$m" ]] || {
		echo "[ERR] manifest not found: '$m'" >&2
		exit 2
	}
done <"$norm_list"

args=(--manifest-list "$norm_list" --out "$OUT" --dedup "$DEDUP" --select "$SELECT")
[[ -n "$PLATFORM_MAP" ]] && args+=(--platform-map "$PLATFORM_MAP")
[[ -n "$SUMMARY_SCOPES" ]] && args+=(--summary-scopes "$SUMMARY_SCOPES")
[[ -n "$SUMMARY_STAT" ]] && args+=(--summary-stat "$SUMMARY_STAT")
((MARKDOWN == 1)) && args+=(--markdown)

python3 "${_POLAPLIB_DIR}/scripts/make_table_dataset_summary_merge.py" "${args[@]}"
