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
# Functions for subcommand template ...
# Describe what they are and what they do.
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

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
	echo "[ERROR] This script must be sourced, not executed: use 'source $BASH_SOURCE'" >&2
	return 1 2>/dev/null || exit 1
fi
: "${_POLAP_DEBUG:=0}"
: "${_POLAP_RELEASE:=0}"

function _run_polap_subcommand {
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

	local polap_cmd="${FUNCNAME##*_}"

	help_message=$(
		cat <<EOF
Name:
  polap - assemble plant organelle-genome sequences using ONT long-read data

Synopsis:
  polap [subcommand] [-l|--long-reads FASTQ]
                 [-a|--short-read1 FASTQ] [-b|--short-read2 FASTQ]
                 [-o|--outdir STR] [--help] [-v] [--version]
                 [-t|--threads INT] [--random-seed INT]

Description:
  polap assemble runs the whole-genome assembly using Flye and then
organelle-genome assembly using seed contig sequences.

To assemble mitochondrial DNA (mtDNA), follow a series of sequential steps:

  ${_polap_command_string} init -o <arg>

  ${_polap_command_string} summary-reads -a <arg> [-b <arg>]

  ${_polap_command_string} total-length-long -l <arg>

  ${_polap_command_string} find-genome-size -a <arg> [-b <arg>]

  ${_polap_command_string} reduce-data -l <arg> [-m <arg>]

  ${_polap_command_string} flye1 [-t <arg>]

  ${_polap_command_string} edges-stats -i <arg>

  ${_polap_command_string} annotate -i <arg>

  ${_polap_command_string} seeds [-i <arg>] -j <arg>

  ${_polap_command_string} map-reads [-i <arg>] -j <arg>

  ${_polap_command_string} test-reads [-i <arg>] -j <arg> -s <begin>,<end>,<count> [-c <arg>]

  ${_polap_command_string} select-reads [-i <arg>] -j <arg> -w <arg> [-c <arg>]

  ${_polap_command_string} flye2 [-i <arg>] -j <arg>

  init - inititalize a polap output folder

  summary-reads - summarize short-read data

  total-length-long - count the total length of the long-read data

  find-genome-size - estimate the genome size using short-read data

  reduce-data - reduce the long-read data for the whole-genome assembly
  
  flye1 - execute Flye for a whole-genome assembly
  
  edges-stats - prepare the assembly graph summary from Flye's graph_final.gfa
  
  map-reads - map long reads on the seed contigs using minimap2
  
  test-reads - iterate organelle genome assemblies on a range of -w option values
  
  select-reads - use a -w option value after test-reads
  
  flye2 - execute Flye for an organelle-genome assembly
  
  blast-genome - execute NCBI's BLAST on seed contigs using organelle amino acid sequences
  
  conut-genes - count genes on seed contigs after blast-genome
  
  flye-polishing - execute Flye's polishing stage
  
  make-menus - create empty files for easy command typing
  
  clean-menus - delete the empty menu-name files
  
  list - list commands
  
  get-bioproject - fetch NCBI's BioProject info
  
  bioproject-prepare - arrange BioProject info to determine which data to download
  
  get-bioproject-sra - fetch NCBI's SRA data 
  
  get-mtdna - fetch known organelle genome sequences

  simulate - simulate long-read data from known organelle genome sequences

Options:
  -l FASTQ
    reads data file
	
  -a FASTQ
    short-read data file 1

  -b FASTQ
    short-read data file 2

  -o OUTDIR
    output folder

  -m INT
    minimum read length [default: ${_arg_min_read_length}]

  -t INT
    the number of CPU cores [default: ${_arg_threads}]

  -c INT
    the coverage option [default: ${_arg_coverage}]

  -g <arg>
    computed by find-genome-size polap command or given by users

  --no-reduction-reads
    no data reduction in a whole-genome assembly

  --flye-asm-coverage INT
    Flye option --asm-coverage [default: 50]

  --random-seed INT
    a 5-digit number

  --redo
    do not use previously generated intermediate results

  --stopafter {data,flye1}
    stop after data or flye1 step

  -o, --outdir: output folder name (default: ${_arg_outdir})
    The option '-o' or '--outdir' specifies the output folder name, 
    with a default value of 'o'. The output folder typically contains input 
    files that are long-read and short-read data files. Input data files can 
    be specified using the options provided by -l, -a, and -b.

  -l, --long-reads: long-reads data file in fastq format (default: ${_arg_long_reads})
    The option '-l' or '--long-reads' specifies the location of a long-reads 
    data file in fastq format, with a default filename of 'l.fq'.

  -a, --short-read1: short-read fastq file 1 (default: ${_arg_short_read1})
    The option '-a' or '--short-read1' specifies the first short-read fastq 
    file to be used, with a default value of "s1.fq".

  -b, --short-read2: short-read fastq file 2 (default: ${_arg_short_read2})
    The option '-b' or '--short-read2' specifies a short-read fastq file 2, 
    with a default value of 's2.fq'. The second short-read data file, 
    if provided, is considered optional.
    (Note: not tested yet; -a & -b are required.)

  -i, --inum: previous output number of organelle-genome assembly (default: ${_arg_inum})
    The option '-i' or '--inum' specifies the previous output number 
    of an organelle-genome assembly, with a default value of '0'.
    The zero for this option specifies the whole-genome assembly.

  -j, --jnum: current output number of organelle-genome assembly (default: ${_arg_jnum})
    The option '-j' or '--jnum' allows users to specify the current output number 
    for an organelle-genome assembly, with a default value of '1'.

  -m, --min-read-length: minimum length of long reads (default: ${_arg_min_read_length})
    The option '-m' or '--min-read-length' specifies the minimum length of 
    long reads, with a default value of 3000. 

  -t, --threads: number of CPUs (default: maximum number of cores)
    The option '-t' or '--threads' specifies the number of CPU threads to 
    utilize, with a default value equal to the maximum number of available 
    cores in the computer you execute 'polap'.

  -c, --coverage: coverage for the organelle-genome assembly (default: ${_arg_coverage})
    The option '-c' or '--coverage' specifies the coverage percentage for the 
    organelle-genome assembly for controlling the data size.

  -w, --single-min: minimum mapped bases or PAF 11th column (default: ${_arg_single_min})
    This parameter ensures that the alignment level between a long-read and a
    seed contig is properly controlled. For plant mitochondrial DNAs, a DNA
    fragment size of approximately 3 kilobases appears to be more effective
    than the smaller 1-kilobase fragment. In the case of plastid DNAs, a
    fragment size of 1 kilobase (kb) might be more suitable, requiring an
    adjustment to the -m option accordingly.

  -g, --genomesize: expected genome size (default: estimated with a short-read dataset)
    Users can assemble an organelle genome when they have a genome size
    estimate. But, we require a short-read dataset to determine the genome size
    for the whole-genome assembly process. Polishing a long-read assembly
    necessitates the use of a short-read dataset.

  -p, --unpolished-fasta: polishing sequence in fasta format (default: ${_arg_unpolished_fasta})
    The option enables the polishing of sequences in a FASTA format, 
    with the default output file being named 'mt.0.fasta'. 

  -f, --final-assembly: final assembly in fasta format (default: ${_arg_final_assembly})
    The final assembly in FASTA format, with a default file name of 'mt.1.fa'. 

  --no-reduction-reads: reduction of long-read data before assemble1
    In the process of whole-genome assembly, we utilize a reduced amount of
    long-read data. By default, we reduce the size of a long-read dataset prior
    to performing a whole-genome assembly.
    Note: The size of coverage is set by --coverage (default: ${_arg_coverage})

  --no-coverage-check: coverage check before assemble2 step
    By default, in the process of assembling organelle genomes, we reduce
    the size of selected seed reads.
    Note: The size of coverage is set by --coverage (default: ${_arg_coverage})

  --yes: always yes for a question or deletes output completely (off by default)

  --redo: redo a POLAP pipeline (off by default)
    The command specifies that any previously generated intermediate results 
    should be disregarded and new calculations performed from scratch.

  -s, --select-read-range: start,end,number for the range of read selection (default: ${_arg_select_read_range})
    It specifies the values for ptGAUL read-selection minimum number of 
    bases or ratios. For the start and end values of a ratio, real numbers
		must 
    fall within the range of 0 to 1.
    Note: refer to the menu "test-reads" for help.

  --start-index: used by test-reads

  --end-index: used by test-reads

  --random-seed: 5-digit number (default automatically assigned)
    To ensure reproducibility, you can supply a random number seed 
    to facilitate sampling of reads.
    0 or negative for automatically assigned
    seqkit sample random seed; 11 used in seqkit sample.

  --flye-asm-coverage: Flye --asm-coverage (default: ${_arg_flye_asm_coverage})
    Flye --asm-coverage is a parameter used with the assembly coverage of
		Flye.

  --no-flye-asm-coverage: no use of Flye --asm-coverage
    The flag '--no-flye-asm-coverage' indicates that we use Flye option 
    neither --asm-coverage nor --genome-size in flye execution.
    This option is the same as --flye-asm-coverage set to 0.
    Note: not tested yet!

  --polap-reads: use intra- and inter-contig read selection (default: ${_arg_polap_reads})
    The default read selection is ptGAUL's approach.
    This option allows long reads that are mapped within a seed contig and
    between two contigs.

  --directional: (default: ${_arg_directional})
    Use only forward strands of seed contigs and reads mapped on the seeds.
    Seed contigs need to be specific to direction; edge_1+ not just edge_1.
    Each line of mt.contig.name file could have multiple edge name with plus or
    minus sign separated by a comma.
    Example: edge_1+, edge_2-, edge_7+, edge_2+

  --blast: (default: ${_arg_blast})

  --bridge-same-strand: (default: ${_arg_bridge_same_strand})
    When linking two inverted repeats, enabling this feature ensures that 
    the strands are equal for the two mapped IR contigs.
    Note: currently only plus strand is used.

  --log: log file (default: <output>/polap.log)
    The log file option allows users to specify a custom log file location, 
    with a default setting of '<output>/polap.log'.

  --clock: display the start and ending time (default: ${_arg_clock})
    The clock option allows users to display both the start and end times.

  --markdown: display the table in markdown format (default: ${_arg_markdown})

  disassemble options:
    Use help menu of the subcommand to see more example commands.

  --downsample : maximum genome coverage to downsample (default: ${_arg_downsample})
    The coverage for downsampling before assembling the plastid genome.
    The genome size for the coverage is computed using short-read data.
    Use option -i for a new downsampled set of data.
    The default is a recommended downsample depth relative to a genome
		size estimate.

  --disassemble-p: the percentile of the largest long-read (default: ${_arg_disassemble_p})
    The maximum percentage for the long-read data subsampling. Use this or
    option --disassemble-b for the maximum subsample size for a range of the
    subsampling.
    The default is a recommended value given the downsample value.

  --disassemble-n: the number of steps (default: ${_arg_disassemble_n})
    The number of iterations in the first stage of the subsampling-based assembly.
    The default is a recommended number of iterations.

  --disassemble-r: the number of replicates (default: ${_arg_disassemble_r})
    The number of iterations in the second and third stages.
    The default is a recommended number of iterations.

  --disassemble-q: the percentile of the largest short-read (default: ${_arg_disassemble_p})
    The maximum percentage for the short-read data subsampling.
    If not specified, we use the option value of --disassemble-p as the value
    for --disassemble-q.

  --disassemble-i: the index used for a separate plastid assemblies (default: ${_arg_disassemble_i})
    Use the same downsampled data to assemble ptDNA with a different set of
    options. Use option -i to have different downsampled data.

  --disassemble-a: the smallest base pairs for a subsampling range (default: ${_arg_disassemble_a})
    The smallest subsample size: default is 10 Mb. This number must be greater
    than the largest one using option --disassemble-b.
    The option --disassemble-p should be chosen so that the largest subsample
    size should be greater than the value of option --disassemble-a.

  --disassemble-b: the largest base pairs for a subsampling range (default: ${_arg_disassemble_b})
    The largest subsample size in base pairs. 

  --disassemble-m: the upper bound for a Flye assembly (default: ${_arg_disassemble_m})
    The upper bound for an initital preassembly size.
    The default for plastid genome assembly is 500 kb.
    We keep the Flye's preassembly size under this value by adjusting
    the read-coverage threshold.

  --disassemble-memory: the maximum memory in Gb (default: ${_arg_disassemble_memory})
    The memory maximum requirement for the plastid genome assembly.
    Note that this value does not guarantee the overrun of memory. 
    You can choose to allocate a smaller amount of memory than your computer's 
    total capacity if you are constrained by limited resources.
    You should adjust this based on the specifics of your plastid genome 
    assembly outcome.

  --disassemble-alpha: the starting Flye's disjointig coverage (default: ${_arg_disassemble_alpha})
    This option should remain unchanged initially.
    Use it if you really have a long-run analysis that you want to be shortened.

  --disassemble-delta: the move size of alpha (0.1 - 1.0) (default: ${_arg_disassemble_delta})
    Leave this option as it is, without any changes.

  --disassemble-s: subsample size for stage 2, skipping stage 1 (default: ${_arg_disassemble_s})
    If you are sure of the subsample size and the read-coverage threshold,
    set this option to bypass Stage 1.

  --disassemble-beta: subsample rate for stage 2, skipping stage 1 (default: ${_arg_disassemble_beta})
    Use either this option or --disassemble-s to set the subsample size.
    Setting both of them does not make sense.
    If you are sure of the subsample size and the read-coverage threshold,
    set this option to bypass Stage 1.

  --disassemble-align-reference: (default: ${_arg_disassemble_align_reference})
    Use this and --disassemble-c so that you can compare a known ptDNA and
    the subsampling-based assembly.

  --disassemble-c: a single reference sequence in FASTA (no default)
    If you want to compare a ptDNA sequence with a subsampling-based assembly,
    set this option to a FASTA file with a single sequence. It will be compared
    with the assembled ptDNA.

  --jellyfish-s: JellyFish's -s option (default: ${_arg_jellyfish_s})
    For disassemble subcommand, it is set to 2G if not set here.
    The -s 2G of JellyFish option claims memory of about 6 GB.
    If it is slow in assembly size estimate, increase this value to
    e.g., 5g or 10g. It would requires more memory.

  Experimental (not implemented yet!):

  --nano-raw (default), --nano-corr, --nano-hq:
  --pacbio-raw, --pacbio-corr, --pacbio-hifi:
    The Flye program requires a specific input data type.
    If one switch is activated, the other switches are automatically deactivated.
    Note: not tested yet!
  --timing: turn on timing and memory usage

  --species: Species scientific name (no default)
  --sra: SRA data (no default)

  -v, --verbose: use multiple times to increase the verbose level
  --version: Prints version
  -h: Prints polap global help
  --help: Prints menu help


Outputs:
  ${_polap_var_outdir_s1_fq_stats}

  ${_polap_var_outdir_s2_fq_stats}

  ${_polap_var_outdir_long_total_length}

  ${_polap_var_outdir_genome_size}

  ${_polap_var_outdir_nk_fq_gz}

  ${_polap_var_outdir_lk_fq_gz}

  ${_polap_var_wga_contigger_edges_gfa}

  ${_polap_var_ga_contigger_edges_stats}

Examples:
  Get organelle genome sequences:
    polap assemble -l l.fq -a s1.fq -b s2.fq -o outdir

See also:
  init, download, download-sra, download-bioproject, 
  download-ftp, download-genbank, download-refseq, 
  download-accession, download-species, 
  download-sra-fastq, download-sra-fasta, 
  test-reads, total-length-long, find-genome-size,
  reduce-data, flye1, edges-stats,
  assemble, assemble1, annotate, assemble2, flye-polishing, 
  make-menus, list, clean-menus, cleanup, init,
  summary-reads, total-length-long, find-genome-size, reduce-data, flye1
  blast-genome, count-gene, seeds,
  prepare-seeds, map-reads, test-reads, select-reads, flye2,
  flye-polishing, prepare-polishing, polish,
  version

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi

EOF
	)

	# Display help message
	_polap_lib_help-maybe-show3 "$polap_cmd" help_message || return 0

	_run_polap_assemble

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}
