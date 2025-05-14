################################################################################
# This file is part of polap.
#
# polap is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# polap is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# polap. If not, see <https://www.gnu.org/licenses/>.
################################################################################

################################################################################
# Ensure that the current script is sourced only once
source "${_POLAPLIB_DIR}/run-polap-function-include.sh"
_POLAP_INCLUDE_=$(_polap_include "${BASH_SOURCE[0]}")
set +u
if [[ -n "${!_POLAP_INCLUDE_}" ]]; then
	set -u
	return 0
fi
set -u
declare "$_POLAP_INCLUDE_=1"
#
################################################################################

function _run_polap_nextdenovo {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

	# Print help message if requested
	help_message=$(
		cat <<HEREDOC
NextDenovo substitutes flye assembler.

This is the most high level function for nextdenovo command.

1. run nextDenovo using run.cfg
2. convert gfa to assembly stats and fasta
3. annotate the fasta and update the assembly stats with MT/PT gene count
4. seed selection and mapping
5. use the selected reads to apply to nextDenovo assembler for mtDNA

Applications:
We could use PacBio HiFi data.
We could compare the results using two assemblers: NextDenovo and Flye.
  We need to make Flye to use PacBio HiFi data.
We could directly compare PMAT with Polap.

TODO:
We need a working script of using NextDenovo to assemble.

Arguments:
  -i ${_arg_inum}: source Flye (usually whole-genome) assembly number

Inputs:
  ${_polap_var_ga_annotation_all}

Outputs:
  ${_polap_var_mtcontigname}

See:
  run-polap-select-contigs-by-table-1.R for the description of --select-contig option
Example: $(basename $0) ${_arg_menu[0]} [-i|--inum <arg>] [-j|--jnum <arg>] [--select-contig <number>]
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == h* || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
		exit $EXIT_SUCCESS
	fi

	echo "verbose level: ${_arg_verbose}" >&2
	echoall "command: $0"
	echoall "function: $FUNCNAME"
	echoerr "LOG: echoerr"
	echoall "LOG: echoall"

	echoerr "LOG: echoerr"
	verbose_echo 0 "Log level   - screen        polap.log file" 1>&2
	_polap_log0 "Log level 0 - nothing        minimal log - --quiet"
	_polap_log1 "Log level 1 - minimal        step info and main io files"
	_polap_log2 "Log level 2 - main io files  inside of function: file input/output --verbose"
	_polap_log3 "Log level 3 - files inside   all log or details of file contents --verbose --verbose"
	_polap_log0_file "log0.file: main assembly input/output"
	_polap_log1_file "log1.file: step main input/output"
	_polap_log2_file "log2.file: inside detail input/output"
	_polap_log3_file "log3.file: all input/output"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

# use this to compute the depths of contigs
#
# Input:
# reads: 02.cns_align/01.seed_cns.sh.work/seed_cns*/cns.fasta
# assembly: 03.ctg_graph/nd.asm.fasta
#
# Output:
# Contig         Length(bp)    AverageDepth
# ctg000000      98060         28.32
# ctg000010      155563        31.01
# ...
#
_polap_nextdenovo-compute_nd_coverage() {
	local asm="$1"           # Assembly FASTA (e.g., 03_ctg_graph/nd.asm.fasta)
	local reads="$2"         # Long-read file (e.g., reads.fq.gz)
	local outdir="$3"        # Output folder
	local threads="${4:-16}" # Optional: threads (default 16)

	if [[ ! -f "$asm" || ! -f "$reads" ]]; then
		echo "Usage: compute_nd_coverage <assembly.fasta> <reads.fq.gz> <output_dir> [threads]" >&2
		return 1
	fi

	mkdir -p "$outdir"

	echo "🧬 Mapping reads to assembly..."
	minimap2 -x map-ont -t "$threads" "$asm" "$reads" |
		samtools view -@ "$threads" -b -o "$outdir/mapped.bam"

	echo "📦 Sorting BAM..."
	samtools sort -@ "$threads" -o "$outdir/mapped.sorted.bam" "$outdir/mapped.bam"
	samtools index "$outdir/mapped.sorted.bam"

	echo "📊 Calculating depth..."
	samtools depth -a "$outdir/mapped.sorted.bam" >"$outdir/depth.txt"

	echo "📈 Averaging depth per contig..."
	awk '
    {
      sum[$1] += $3;
      count[$1]++;
    }
    END {
      printf "Contig\tLength(bp)\tAverageDepth\n";
      for (c in sum) {
        printf "%s\t%d\t%.2f\n", c, count[c], sum[c] / count[c];
      }
    }
    ' "$outdir/depth.txt" >"$outdir/average_depth.tsv"

	echo "✅ Done. Output: $outdir/average_depth.tsv"
}
