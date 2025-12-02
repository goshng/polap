# delete
#!/usr/bin/env bash
# File: polap-bash-polish-raconpoly-step2-filter-fastq-by-ids.sh
# Version: v0.9.2
set -Eeuo pipefail
# Args: -i reads.fq[.gz] -l ids.txt -o out.fq.gz
while getopts ":i:l:o:" opt; do
	case "$opt" in
	i) FQ="$OPTARG" ;;
	l) IDS="$OPTARG" ;;
	o) OUT="$OPTARG" ;;
	*)
		echo "bad opt"
		exit 2
		;;
	esac
done
[[ -s "$FQ" && -s "$IDS" && -n "$OUT" ]] || {
	echo "args?"
	exit 2
}
if command -v seqkit >/dev/null 2>&1; then
	seqkit grep -f "$IDS" "$FQ" -o "$OUT"
else
	POLAPLIB_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd -P)"
	python3 "${POLAPLIB_DIR}/scripts/polap_py_fastq_filter_by_ids.py" \
		--fastq "$FQ" --ids "$IDS" --out "$OUT"
fi
