#!/bin/bash

# genome_size=$(<"${_arg_outdir}"/o/short_expected_genome_size.txt)
# command time -v PMAT autoMito -i ${_arg_outdir}/0/cns.fa -o ${_arg_outdir}/0/pmat-fc-${_arg_fc} -st ont -g ${genome_size} --task p1 --type mt -cpu 24 -m -fc ${_arg_fc} >${_arg_outdir}/0/stdout-pmat-nextdenovo.txt 2>${_arg_outdir}/0/pmat/timing-pmat-nextdenovo.txt

#!/bin/bash

# Check if required arguments are provided
if [[ -z "$1" || -z "$2" || -z "$3" ]]; then
	echo "Usage: $0 <outdir> <inum> <type> [fc]"
	echo
	echo "Arguments:"
	echo "  outdir   Output directory path (must exist)"
	echo "  inum     Input identifier (must not already exist under outdir)"
	echo "  type     Must be one of: all, mt, pt"
	echo "  fc       Optional. A number in [0, 1]. Default: 0"
	exit 1
fi

outdir="$1"
inum="$2"
type="${3:-all}"
fc="${4:-0}" # default value is 0

# Check if outdir exists
if [[ ! -d "$outdir" ]]; then
	echo "Error: Directory '$outdir' does not exist." >&2
	exit 1
fi

# Check if outdir/inum already exists
if [[ ! -d "$outdir/$inum" ]]; then
	echo "Error: Directory '$outdir/$inum' does not exist." >&2
	exit 1
fi

if [[ ! -s "$outdir/$inum/cns.fa" ]]; then
	echo "Error: file '$outdir/$inum/cns.fa' does not exist." >&2
	exit 1
fi

# Check if type is valid
if [[ ! "$type" =~ ^(all|mt|pt)$ ]]; then
	echo "Error: type ($type) must be one of: all, mt, pt" >&2
	exit 1
fi

# Check if fc is a valid number in [0, 1]
if ! [[ "$fc" =~ ^[0-9]*\.?[0-9]+$ ]] || (($(echo "$fc < 0 || $fc > 1" | bc -l))); then
	echo "Error: fc must be a number in the range [0, 1]" >&2
	exit 1
fi

if [[ "${fc}" -eq 0 ]]; then
	awk 'BEGIN {for (i = 0.1; i <= 0.9; i += 0.1) printf "%.1f\n", i}' | while read i; do
		bash src/polap-data-v2.sh pmat-nextdenovo2 ${outdir} ${inum} ${type} $i
	done
else
	bash src/polap-data-v2.sh pmat-nextdenovo2 ${outdir} ${inum} ${type} ${fc}
fi
