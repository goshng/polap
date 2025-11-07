# FILE: polap-bash-init-analysis
#!/usr/bin/env bash
set -euo pipefail
VERSION="0.1.0"

usage() {
	cat <<EOF
polap-bash-init-analysis v${VERSION}
Link/copy mt.fa and pt.fa for a species into man/analysis/assemblies/<species>/.

USAGE:
  $(basename "$0") -s <species> -m <mt.fa> -c <pt.fa> [-b <analysis_base>] [--copy]

OPTIONS:
  -s  Species name (no spaces)
  --mt  Path to mitochondrial FASTA (mt.fa)
  --pt  Path to plastid FASTA (pt.fa)
  -b  Analysis base (default: man/analysis)
  --copy  Copy instead of symlink

EOF
}

ANALYSIS_BASE="man/analysis"
COPY=0
species="" mtfa="" ptfa=""
while [[ $# -gt 0 ]]; do
	case "$1" in
	-s)
		species="$2"
		shift 2
		;;
	--mt)
		mtfa="$2"
		shift 2
		;;
	--pt)
		ptfa="$2"
		shift 2
		;;
	-b)
		ANALYSIS_BASE="$2"
		shift 2
		;;
	--copy)
		COPY=1
		shift
		;;
	-h | --help)
		usage
		exit 0
		;;
	*)
		echo "Unknown arg: $1"
		usage
		exit 1
		;;
	esac
done

[[ -z "${species}" || -z "${mtfa}" || -z "${ptfa}" ]] && {
	usage
	exit 1
}

mkdir -p "${ANALYSIS_BASE}/assemblies/${species}" "${ANALYSIS_BASE}/logs"
dest="${ANALYSIS_BASE}/assemblies/${species}"
log="${ANALYSIS_BASE}/logs/init_${species}.log"

echo "[$(date)] init ${species}" | tee "${log}"
if [[ ${COPY} -eq 1 ]]; then
	echo cp -f "${mtfa}" "${dest}/mt.fa"
	echo cp -f "${ptfa}" "${dest}/pt.fa"
	cp -f "${mtfa}" "${dest}/mt.fa"
	cp -f "${ptfa}" "${dest}/pt.fa"
else
	ln -sf "$(realpath "${mtfa}")" "${dest}/mt.fa"
	ln -sf "$(realpath "${ptfa}")" "${dest}/pt.fa"
fi
echo "Done. Files at ${dest}/." | tee -a "${log}"
