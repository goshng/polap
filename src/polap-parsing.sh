#!/bin/bash

# Created by argbash-init v2.10.0
# ARG_OPTIONAL_SINGLE([long-reads],[l],[long-reads data file in fastq format],[l.fq])
# ARG_OPTIONAL_SINGLE([outdir],[o],[output folder name],[o])
# ARG_OPTIONAL_SINGLE([short-read1],[a],[short-read fastq file 1],[s1.fq])
# ARG_OPTIONAL_SINGLE([short-read2],[b],[short-read fastq file 2],[s2.fq])
# ARG_OPTIONAL_SINGLE([sra],[],[SRA data])
# ARG_OPTIONAL_SINGLE([unpolished-fasta],[p],[polishing sequence in fasta format],[mt.0.fasta])
# ARG_OPTIONAL_SINGLE([final-assembly],[f],[final assembly in fasta format],[mt.1.fa])
# ARG_OPTIONAL_SINGLE([min-read-length],[m],[minimum length of long reads],[3000])
# ARG_OPTIONAL_SINGLE([threads],[t],[number of CPUs],[$(cat /proc/cpuinfo | grep -c processor)])
# ARG_OPTIONAL_SINGLE([log],[],[log file],[polap.log])
# ARG_OPTIONAL_SINGLE([coverage],[c],[step4: coverage for the 2nd assembly],[30])
# ARG_OPTIONAL_SINGLE([pair-min],[r],[minimum mapped bases or PAF 11th column],[3000])
# ARG_OPTIONAL_SINGLE([bridge-min],[x],[minimum bridging read length or PAF 7th column],[3000])
# ARG_OPTIONAL_SINGLE([single-min],[w],[minimum mapped bases or PAF 11th column],[3000])
# ARG_OPTIONAL_SINGLE([inum],[i],[previous output number of organelle-genome assembly],[0])
# ARG_OPTIONAL_SINGLE([jnum],[j],[current output number of organelle-genome assembly],[1])
# ARG_OPTIONAL_SINGLE([genomesize],[g],[expected genome size])
# ARG_OPTIONAL_SINGLE([bioproject],[],[NCBI BioProject ID])
# ARG_OPTIONAL_SINGLE([species],[],[Species scientific name])
# ARG_OPTIONAL_SINGLE([accession],[],[NCBI accession ID])
# ARG_OPTIONAL_SINGLE([query],[],[query sequence for blastn])
# ARG_OPTIONAL_SINGLE([subject],[],[subject sequence for blastn])
# ARG_OPTIONAL_REPEATED([minimum],[M],[Pair, Bridge, and Single minimum],[])
# ARG_OPTIONAL_BOOLEAN([reduction-reads],[],[redo: no reduction of long-read data],[on])
# ARG_OPTIONAL_BOOLEAN([contigger],[],[step1: use flye's 40-polishing result])
# ARG_OPTIONAL_BOOLEAN([all-annotate],[],[step2: annotate all contigs])
# ARG_OPTIONAL_BOOLEAN([bridge-same-strand],[],[step4: use flye's edges not contigs])
# ARG_OPTIONAL_BOOLEAN([coverage-check],[],[step4: no coverage check for step 4])
# ARG_OPTIONAL_BOOLEAN([resume],[],[step1,step4: flye option resume])
# ARG_OPTIONAL_BOOLEAN([circularize],[u],[step4: circularize a contig])
# ARG_POSITIONAL_MULTI([menu],[Polap menu],[3],[assemble],[infile],[outfile])
# ARG_VERSION([echo $0 v0.2.6])
# ARG_VERBOSE([])
# ARG_DEFAULTS_POS([])
# ARG_HELP(['P'lant 'o'rganelle DNA 'l'ong-read 'a'ssembly 'p'ipeline.])
# DEFINE_SCRIPT_DIR([])
# ARGBASH_GO()
# needed because of Argbash --> m4_ignore([
### START OF CODE GENERATED BY Argbash v2.10.0 one line above ###
# Argbash is a bash code generator used to get arguments parsing right.
# Argbash is FREE SOFTWARE, see https://argbash.io for more info

die() {
	local _ret="${2:-1}"
	test "${_PRINT_HELP:-no}" = yes && print_help >&2
	echo "$1" >&2
	exit "${_ret}"
}

begins_with_short_option() {
	local first_option all_short_options='loabpfmtcrxwijgMuvhs'
	first_option="${1:0:1}"
	test "$all_short_options" = "${all_short_options/$first_option/}" && return 1 || return 0
}

# THE DEFAULTS INITIALIZATION - POSITIONALS
_positionals=()
_arg_menu=("assemble" "infile" "outfile")
# THE DEFAULTS INITIALIZATION - OPTIONALS
_arg_long_reads="l.fq"
_arg_long_reads_is="off"
_arg_outdir="o"
_arg_archive="a"
_arg_archive_is="off"
_arg_short_read1="s1.fq"
_arg_short_read1_is="off"
_arg_short_read2="s2.fq"
_arg_short_read2_is="off"
_arg_sra=
_arg_unpolished_fasta="mt.0.fasta"
_arg_final_assembly="mt.1.fa"
_arg_table_format="tsv"
_arg_outfile="o.tsv"
_arg_min_read_length="3000"
_arg_threads="$(cat /proc/cpuinfo | grep -c processor)"
_arg_log="polap.log"
_arg_log_is="off"
_arg_log_stderr="off"
_arg_coverage="50" # 2024-10-29 was change to 50
_arg_flye_asm_coverage="30"
_arg_single_min="3000"
_arg_pair_min="3000"
_arg_bridge_min="0"
_arg_max_seeds="35" # Lolium perenne: over 30 less than 35
_arg_inum="0"
_arg_jnum="1"
_arg_knum="1"
_arg_start_index=0
_arg_select_contig="1"
_arg_select_contig_numbers=(1 2 3 4 5 6)
_arg_select_read_range="3000,39000,7"
_arg_select_read_range_is="off"
_arg_report_x="5000,7000,9000,11000,13000,15000,17000"
_arg_report_x_is="off"
_arg_random_seed=
_arg_genomesize=
_arg_bioproject=
_arg_species=
_arg_accession=
_arg_query=
_arg_subject=
# Add options of on and off
_arg_markdown="off"
_arg_flye="on"
_arg_reduction_reads="on"
_arg_contigger="on"
_arg_all_annotate="off"
_arg_polap_reads="off"
_arg_bridge_same_strand="off"
_arg_coverage_check="on"
_arg_resume="off"
_arg_plastid="off"
_arg_clock="off"
_arg_yes="off"
_arg_circularize="off"
_arg_redo="off"
_arg_test="off"
_arg_verbose=1
_arg_help="off"
# flye options
_arg_flye_data_type="--nano-raw"

source "$script_dir/polap-git-hash-version.sh"
_polap_version=v0.3.8.1-"${_polap_git_hash_version}"
_polap_command_string=polap

print_help() {

	help_message=$(
		cat <<HEREDOC
POLAP - Plant organelle DNA long-read assembly pipeline.
version ${_polap_version}

Usage: ${_polap_command_string} [<menu> [<menu2> [<menu3>]]] [-o|--outdir <arg>]
      [-l|--long-reads <arg>] [-a|--short-read1 <arg>] [-b|--short-read2 <arg>]
      [-i|--inum <arg>] [-j|--jnum <arg>] [-w|--single-min <arg>]
      [-m|--min-read-length <arg>] [-t|--threads <arg>] [--test] [--log <arg>] 
      [--random-seed <arg>] [--version] [-h|--help]

Assemble mitochondrial DNA (mtDNA) in a single command (not tested yet):
  ${_polap_command_string} -l <arg> -a <arg> [-b <arg>]
  ${_polap_command_string} assemble -l <arg> -a <arg> [-b <arg>]
 
Perform a polishing of the mtDNA sequence utilizing the FMLRC protocol:
  ${_polap_command_string} prepare-polishing  -a <arg> [-b <arg>]
  ${_polap_command_string} polish -p <arg> -f <arg>

To assemble mitochondrial DNA (mtDNA), follow a two-step process involving
manual seed contig selection:
  ${_polap_command_string} assemble1 -l <arg> -a <arg> [-b <arg>] [-m <arg>]
  ${_polap_command_string} assemble2 -i <arg> -j <arg> [-w <arg>] [-c <arg>]

To assemble mitochondrial DNA (mtDNA), follow a three-step process
that utilizes semi-automatic seed contig selection:
  ${_polap_command_string} assemble1 -l <arg> -a <arg> [-b <arg>] [-m <arg>]
  ${_polap_command_string} seeds -i <arg> -j <arg>
  ${_polap_command_string} assemble2 -i <arg> -j <arg> [-w <arg>] [-c <arg>]

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

Others menus:
  ${_polap_command_string} blast-genome -i <arg>
  ${_polap_command_string} count-genes -i <arg>
  ${_polap_command_string} flye-polishing -j <arg>
  ${_polap_command_string} make-menus
  ${_polap_command_string} clean-menus
  ${_polap_command_string} list

BioProject menus:
  ${_polap_command_string} get-bioproject --bioproject <arg>
  ${_polap_command_string} bioproject-prepare -o <arg>
  ${_polap_command_string} get-bioproject-sra --sra <arg>
  ${_polap_command_string} get-mtdna --species <arg>

Other options:
      [-p|--unpolished-fasta <arg>] [-f|--final-assembly <arg>]
      [-c|--coverage <arg>] [--flye-asm-coverage <arg>]
      [--bioproject <arg>] [--species <arg>] [--accession <arg>]
      [--query <arg>] [--subject <arg>]
      [--no-reduction-reads] [--no-coverage-check]
      [--plastid]
      [--archive <arg>]
      [--sra <arg>] [-g|--genomesize <arg>]

menu: assemble, assemble1, annotate, assemble2, flye-polishing, 
      make-menus, list, clean-menus, cleanup, init,
      summary-reads, total-length-long, find-genome-size, reduce-data, flye1
      blast-genome, count-gene, seeds,
      prepare-seeds, map-reads, test-reads, select-reads, flye2,
      flye-polishing, prepare-polishing, polish,
      version

Options:
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
    cores.

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
    bases or ratios. For the start and end values of a ratio, real numbers must 
    fall within the range of 0 to 1.
    Note: refer to the menu "test-reads" for help.

  --start-index: used by test-reads

  --random-seed: 5-digit number (default automatically assigned)
    To ensure reproducibility, you can supply a random number seed 
    to facilitate sampling of reads.
    seqkit sample random seed; 11 used in seqkit sample.

  --flye-asm-coverage: Flye --asm-coverage (default: ${_arg_flye_asm_coverage})
    Flye --asm-coverage is a parameter used with the assembly coverage of Flye.

  --no-flye-asm-coverage: no use of Flye --asm-coverage
    The flag '--no-flye-asm-coverage' indicates that we use Flye option 
    neither --asm-coverage nor --genome-size in flye execution.
    This option is the same as --flye-asm-coverage set to 0.
    Note: not tested yet!

  --polap-reads: use intra- and inter-contig read selection (default: ${_arg_polap_reads})
    The default read selection is ptGAUL's approach.
    This option allows long reads that are mapped within a seed contig and
    between two contigs.

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

  Experimental (not tested yet!):
  --flye-nano-raw (default), --flye-nano-corr, --flye-nano-hq:
  --flye-pacbio-raw, --flye-pacbio-corr, --flye-pacbio-hifi:
    The Flye program requires a specific input data type.
    If one switch is activated, the other switches are automatically deactivated.
    Note: not tested yet!

  --species: Species scientific name (no default)
	--sra: SRA data (no default)

  -v, --verbose: use multiple times to increase the verbose level
  --version: Prints version
  -h: Prints polap global help
  --help: Prints menu help

Example:
git clone https://github.com/goshng/polap.git
cd polap/test
polap --test

Place your long-read and short-read files at a folder:
long-read file: l.fq
short-read file: s1.fq, s2.fq
Execute: polap init
HEREDOC
	)

	# Display help message
	echo "${help_message}"
}

parse_commandline() {
	source "$script_dir/polap-version.sh" # '.' means 'source'
	_positionals_count=0
	while test $# -gt 0; do
		_key="$1"
		case "$_key" in
		-l | --long-reads)
			test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
			_arg_long_reads="$2"
			_arg_long_reads_is="on"
			shift
			;;
		--long-reads=*)
			_arg_long_reads="${_key##--long-reads=}"
			_arg_long_reads_is="on"
			;;
		-l*)
			_arg_long_reads="${_key##-l}"
			_arg_long_reads_is="on"
			;;
		-o | --outdir)
			test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
			_arg_outdir="$2"
			_arg_outdir="${_arg_outdir%/}"
			shift
			;;
		--outdir=*)
			_arg_outdir="${_key##--outdir=}"
			_arg_outdir="${_arg_outdir%/}"
			;;
		-o*)
			_arg_outdir="${_key##-o}"
			_arg_outdir="${_arg_outdir%/}"
			;;
		-a | --short-read1)
			test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
			_arg_short_read1="$2"
			_arg_short_read1_is="on"
			shift
			;;
		--short-read1=*)
			_arg_short_read1="${_key##--short-read1=}"
			_arg_short_read1_is="on"
			;;
		-a*)
			_arg_short_read1="${_key##-a}"
			_arg_short_read1_is="on"
			;;
		-b | --short-read2)
			test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
			_arg_short_read2="$2"
			_arg_short_read2_is="on"
			shift
			;;
		--short-read2=*)
			_arg_short_read2="${_key##--short-read2=}"
			_arg_short_read2_is="on"
			;;
		-b*)
			_arg_short_read2="${_key##-b}"
			_arg_short_read2_is="on"
			;;
		--sra)
			test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
			_arg_sra="$2"
			shift
			;;
		--sra=*)
			_arg_sra="${_key##--sra=}"
			;;
		-p | --unpolished-fasta)
			test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
			_arg_unpolished_fasta="$2"
			shift
			;;
		--unpolished-fasta=*)
			_arg_unpolished_fasta="${_key##--unpolished-fasta=}"
			;;
		-p*)
			_arg_unpolished_fasta="${_key##-p}"
			;;
		-f | --final-assembly)
			test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
			_arg_final_assembly="$2"
			shift
			;;
		--final-assembly=*)
			_arg_final_assembly="${_key##--final-assembly=}"
			;;
		-f*)
			_arg_final_assembly="${_key##-f}"
			;;
		-m | --min-read-length)
			test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
			_arg_min_read_length="$2"
			shift
			;;
		--min-read-length=*)
			_arg_min_read_length="${_key##--min-read-length=}"
			;;
		-m*)
			_arg_min_read_length="${_key##-m}"
			;;
		-t | --threads)
			test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
			_arg_threads="$2"
			shift
			;;
		--threads=*)
			_arg_threads="${_key##--threads=}"
			;;
		-t*)
			_arg_threads="${_key##-t}"
			;;
		--flye-asm-coverage)
			test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
			_arg_flye_asm_coverage="$2"
			shift
			;;
		--flye-asm-coverage=*)
			_arg_flye_asm_coverage="${_key##--flye-asm-coverage=}"
			;;
		--no-flye-asm-coverage)
			_arg_flye_asm_coverage="0"
			shift
			;;
		-c | --coverage)
			test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
			_arg_coverage="$2"
			shift
			;;
		--coverage=*)
			_arg_coverage="${_key##--coverage=}"
			;;
		-c*)
			_arg_coverage="${_key##-c}"
			;;
		-r | --pair-min)
			test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
			_arg_pair_min="$2"
			shift
			;;
		--pair-min=*)
			_arg_pair_min="${_key##--pair-min=}"
			;;
		-r*)
			_arg_pair_min="${_key##-r}"
			;;
		-x | --bridge-min)
			test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
			_arg_bridge_min="$2"
			shift
			;;
		--bridge-min=*)
			_arg_bridge_min="${_key##--bridge-min=}"
			;;
		-x*)
			_arg_bridge_min="${_key##-x}"
			;;
		-w | --single-min)
			test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
			_arg_single_min="$2"
			shift
			;;
		--single-min=*)
			_arg_single_min="${_key##--single-min=}"
			;;
		-w*)
			_arg_single_min="${_key##-w}"
			;;
		--rw)
			test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
			_arg_pair_min="$2"
			_arg_single_min="$2"
			_arg_bridge_min="0"
			shift
			;;
		--rw=*)
			_arg_pair_min="${_key##--rw=}"
			_arg_single_min="${_key##--rw=}"
			_arg_bridge_min="0"
			;;
		--max-seeds)
			test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
			_arg_max_seeds="$2"
			shift
			;;
		--max-seeds=*)
			_arg_max_seeds="${_key##--max-seeds=}"
			;;
		-i | --inum)
			test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
			_arg_inum="$2"
			shift
			;;
		--inum=*)
			_arg_inum="${_key##--inum=}"
			;;
		-i*)
			_arg_inum="${_key##-i}"
			;;
		-j | --jnum)
			test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
			_arg_jnum="$2"
			shift
			;;
		--jnum=*)
			_arg_jnum="${_key##--jnum=}"
			;;
		-j*)
			_arg_jnum="${_key##-j}"
			;;
		-k | --knum)
			test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
			_arg_knum="$2"
			shift
			;;
		--knum=*)
			_arg_knum="${_key##--knum=}"
			;;
		-k*)
			_arg_knum="${_key##-k}"
			;;
		--start-index)
			test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
			_arg_start_index="$2"
			shift
			;;
		--start-index=*)
			_arg_start_index="${_key##--start-index=}"
			;;
		--select-contig)
			test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
			_arg_select_contig="$2"
			shift
			;;
		--select-contig=*)
			_arg_select_contig="${_key##--select-contig=}"
			;;
		--random-seed)
			test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
			_arg_random_seed="$2"
			shift
			;;
		--random-seed=*)
			_arg_random_seed="${_key##--random-seed=}"
			;;
		-g | --genomesize)
			test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
			_arg_genomesize="$2"
			shift
			;;
		--genomesize=*)
			_arg_genomesize="${_key##--genomesize=}"
			;;
		-g*)
			_arg_genomesize="${_key##-g}"
			;;
		--archive)
			test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
			_arg_archive="$2"
			_arg_archive_is="on"
			shift
			;;
		--archive=*)
			_arg_archive="${_key##--archive=}"
			_arg_archive_is="on"
			;;
		--bioproject)
			test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
			_arg_bioproject="$2"
			shift
			;;
		--bioproject=*)
			_arg_bioproject="${_key##--bioproject=}"
			;;
		--species)
			test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
			_arg_species="$2"
			shift
			;;
		--species=*)
			_arg_species="${_key##--species=}"
			;;
		--accession)
			test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
			_arg_accession="$2"
			shift
			;;
		--accession=*)
			_arg_accession="${_key##--accession=}"
			;;
		--query)
			test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
			_arg_query="$2"
			shift
			;;
		--query=*)
			_arg_query="${_key##--query=}"
			;;
		--subject)
			test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
			_arg_subject="$2"
			shift
			;;
		--subject=*)
			_arg_subject="${_key##--subject=}"
			;;
		-s | --select-read-range)
			test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
			_arg_select_read_range="$2"
			_arg_select_read_range_is="on"
			shift
			;;
		--report-x)
			test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
			_arg_report_x="$2"
			_arg_report_x_is="on"
			shift
			;;
		--no-reduction-reads | --reduction-reads)
			_arg_reduction_reads="on"
			test "${1:0:5}" = "--no-" && _arg_reduction_reads="off"
			;;
			# Add options of on and off
		--no-markdown | --markdown)
			_arg_markdown="on"
			test "${1:0:5}" = "--no-" && _arg_markdown="off"
			;;
		--no-flye | --flye)
			_arg_flye="on"
			test "${1:0:5}" = "--no-" && _arg_flye="off"
			;;
		--no-contigger | --contigger)
			_arg_contigger="on"
			test "${1:0:5}" = "--no-" && _arg_contigger="off"
			;;
		--no-all-annotate | --all-annotate)
			_arg_all_annotate="on"
			test "${1:0:5}" = "--no-" && _arg_all_annotate="off"
			;;
		--no-log-stderr | --log-stderr)
			_arg_log_stderr="on"
			test "${1:0:5}" = "--no-" && _arg_log_stderr="off"
			;;
		--no-polap-reads | --polap-reads)
			_arg_polap_reads="on"
			test "${1:0:5}" = "--no-" && _arg_polap_reads="off"
			;;
		--no-bridge-same-strand | --bridge-same-strand)
			_arg_bridge_same_strand="on"
			test "${1:0:5}" = "--no-" && _arg_bridge_same_strand="off"
			;;
		--no-coverage-check | --coverage-check)
			_arg_coverage_check="on"
			test "${1:0:5}" = "--no-" && _arg_coverage_check="off"
			;;
		--no-clock | --clock)
			_arg_clock="on"
			test "${1:0:5}" = "--no-" && _arg_clock="off"
			;;
		--no-plastid | --plastid)
			_arg_plastid="on"
			test "${1:0:5}" = "--no-" && _arg_plastid="off"
			;;
		--no-yes | --yes)
			_arg_yes="on"
			test "${1:0:5}" = "--no-" && _arg_yes="off"
			;;
		--no-resume | --resume)
			_arg_resume="on"
			test "${1:0:5}" = "--no-" && _arg_resume="off"
			;;
		-u | --no-circularize | --circularize)
			_arg_circularize="on"
			test "${1:0:5}" = "--no-" && _arg_circularize="off"
			;;
		-u*)
			_arg_circularize="on"
			_next="${_key##-u}"
			if test -n "$_next" -a "$_next" != "$_key"; then
				{ begins_with_short_option "$_next" && shift && set -- "-u" "-${_next}" "$@"; } || die "The short option '$_key' can't be decomposed to ${_key:0:2} and -${_key:2}, because ${_key:0:2} doesn't accept value and '-${_key:2:1}' doesn't correspond to a short option."
			fi
			;;
		--flye-pacbio-raw)
			_arg_flye_data_type="--pacbio-raw"
			;;
		--flye-pacbio-corr)
			_arg_flye_data_type="--pacbio-corr"
			;;
		--flye-pacbio-hifi)
			_arg_flye_data_type="--pacbio-hifi"
			;;
		--flye-nano-raw)
			_arg_flye_data_type="--nano-raw"
			;;
		--flye-nano-corr)
			_arg_flye_data_type="--nano-corr"
			;;
		--flye-nano-fq)
			_arg_flye_data_type="--nano-fq"
			;;
		--outfile)
			test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
			_arg_outfile="$2"
			shift
			;;
		--outfile=*)
			_arg_outfile="${_key##--outfile=}"
			;;
		--table-format)
			test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
			_arg_table_format="$2"
			shift
			;;
		--table-format=*)
			_arg_table_format="${_key##--table-format=}"
			;;
		--log)
			test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
			_arg_log="$2"
			_arg_log_is="on"
			shift
			;;
		--log=*)
			_arg_log="${_key##--log=}"
			_arg_log_is="on"
			;;
		--no-redo | --redo)
			_arg_redo="on"
			test "${1:0:5}" = "--no-" && _arg_redo="off"
			;;
		--no-test | --test)
			_arg_test="on"
			test "${1:0:5}" = "--no-" && _arg_test="off"
			if [[ "${_arg_test}" == "on" ]]; then
				_arg_plastid="on"
			fi
			;;
		--version)
			echo $0 ${_polap_version}
			exit 0
			;;
		-q | --quiet)
			_arg_verbose=0
			;;
		-v | --verbose)
			_arg_verbose=$((_arg_verbose + 1))
			;;
		--help)
			_arg_help="on"
			test "${1:0:5}" = "--no-" && _arg_test="off"
			;;
		-h)
			print_help
			exit 0
			;;
		*)
			_last_positional="$1"
			_positionals+=("$_last_positional")
			_positionals_count=$((_positionals_count + 1))
			;;
		esac
		shift
	done
}

handle_passed_args_count() {
	test "${_positionals_count}" -le 3 || _PRINT_HELP=yes die "FATAL ERROR: There were spurious positional arguments --- we expect between 0 and 3, but got ${_positionals_count} (the last one was: '${_last_positional}')." 1
}

assign_positional_args() {
	local _positional_name _shift_for=$1
	_positional_names="_arg_menu[0] _arg_menu[1] _arg_menu[2] "

	shift "$_shift_for"
	for _positional_name in ${_positional_names}; do
		test $# -gt 0 || break
		eval "$_positional_name=\${1}" || die "Error during argument parsing, possibly an Argbash bug." 1
		shift
	done
}

parse_commandline "$@"
handle_passed_args_count
set +u
assign_positional_args 1 "${_positionals[@]}"
set -u

# OTHER STUFF GENERATED BY Argbash
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" || {
	echo "Couldn't determine the script's running directory, which probably matters, bailing out" >&2
	exit 2
}

### END OF CODE GENERATED BY Argbash (sortof) ### ])
