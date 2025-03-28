#!/usr/bin/bash

# TODO:
# [ ] A manual for proper understanding and execution of tasks.
# [ ] Sall read from the CSV file

# How-To
#
# Edit "${script_dir}/polap-data-v2.csv"
# It looks for the CSV file in the following order.
# 1. The current folder where you execute this script
# 2. The folder in which this script resides
#
# Eucalyptus_pauciflora is set to the main result.
# If needed, this should be changed by editing Smain array below.
#
# Edit the following 3 variables if necessary:
# You may not need to edit them if you are not developing this script.
# _polap_cmd=${_polap_cmd}
# _media_dir=${_media_dir}
# _local_host=${_local_host}
# host=$(hostname)

Smain=(
	'Eucalyptus_pauciflora'
)

_arg_default_target_dir="$HOME/all/manuscript/polap-v0.4/"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" || {
	echo "Couldn't determine the script's running directory, which probably matters, bailing out" >&2
	exit 2
}
LIB_DIR="${script_dir}/polaplib"

source "${LIB_DIR}/polap-lib-timing.sh"

_log_echo() {
	echo "$(date '+%Y-%m-%d %H:%M:%S') [$subarg1] - $1" >>"${output_dir}/polap-data-v2.txt"
	echo "$1"
}

# Data28
Sall=(
	Anthoceros_agrestis
	Arabidopsis_thaliana
	Canavalia_ensiformis
	Cinchona_pubescens
	Codonopsis_lanceolata
	Cucumis_sativus_var_hardwickii
	Dioscorea_japonica
	Dunaliella_tertiolecta
	Eucalyptus_pauciflora
	Euonymus_alatus
	Gossypium_herbaceum
	Juncus_effusus
	Juncus_inflexus
	Juncus_roemerianus
	Juncus_validus
	Leiosporoceros_dussii
	Macadamia_jansenii
	Musa_acuminata_subsp_malaccensis
	Notothylas_orbicularis
	Ophrys_lutea
	Oryza_rufipogon
	Phaeomegaceros_chiloensis
	Populus_x_sibirica
	Prunus_mandshurica
	Solanum_lycopersicum
	Spirodela_polyrhiza
	Vaccinium_vitis-idaea
	Vitis_vinifera
)

S31=(
	Anthoceros_agrestis
	Arabidopsis_thaliana
	Brassica_carinata
	Canavalia_ensiformis
	Cinchona_pubescens
	Codonopsis_lanceolata
	Cucumis_sativus_var_hardwickii
	Dioscorea_japonica
	Dunaliella_tertiolecta
	Eucalyptus_pauciflora
	Euonymus_alatus
	Gossypium_herbaceum
	Juncus_effusus
	Juncus_inflexus
	Juncus_roemerianus
	Juncus_validus
	Leiosporoceros_dussii
	Macadamia_jansenii
	Musa_acuminata_subsp_malaccensis
	Notothylas_orbicularis
	Ophrys_lutea
	Oryza_rufipogon
	Phaeomegaceros_chiloensis
	Pisum_sativum
	Populus_x_sibirica
	Prunus_mandshurica
	Pterocarpus_santalinus
	Solanum_lycopersicum
	Spirodela_polyrhiza
	Vaccinium_vitis-idaea
	Vitis_vinifera
)

Snot=(
	Delphinium_montanum
)

Stable5=(
	x_y
	# 'Juncus_effusus'
	# 'Juncus_inflexus'
	# 'Juncus_roemerianus'
	# 'Juncus_validus'
	# 'Eucalyptus_pauciflora'
)

_local_host="thorne"

# Input parameter
subcmd1="${1:-help}"

_polap_subcmd=(
	'batch'
	'recover'
	'mkdir'
	'refs'
	'getorganelle'
	'ptgaul'
	'msbwt'
	'extract-ptdna-of-ptgaul'
	'coverage'
	'infer'
	'check'
	'compare'
	'best'
	'bandage1'
	'bandage2'
	'copy-figure'
	'downsample'
	'use-downsample'
	'mauve'
	'restart'
	'write-config'
	'wga'
	'archive'
	'get'
	'report'
	'table1'
	'table2'
	'suptable1'
	'supfigure1'
	'supfigure2'
	'clean'
	'menu'
	'help'
)

# Check if subcmd1 is an integer (0 or positive)
if [[ "$subcmd1" =~ ^[0-9]+$ ]]; then
	# Convert to an index and replace with the corresponding _polap_subcmd value
	index=$subcmd1
	if [[ $index -ge 0 && $index -lt ${#_polap_subcmd[@]} ]]; then
		subcmd1="${_polap_subcmd[$index]}"
	else
		echo "Error: Index $index is out of range (0-${#_polap_subcmd[@]})"
		exit 1
	fi

fi

if [[ -d "src" ]]; then
	_polap_cmd="src/polap.sh"
else
	_polap_cmd="polap"
fi
_polap_version="0.4.3.7"
_media1_dir="/media/h1/sra"
_media2_dir="/media/h2/sra"
_media_dir="/media/h2/sra"

help_message=$(
	cat <<HEREDOC
# Polap data analysis for subsampling-based plastid genome assembly

## Install this script with either bioconda or github:

### using bioconda

1. Follow the quick start available at 
   https://github.com/goshng/polap?tab=readme-ov-file#quick-start
   to install polap. The followings are the essential commands for the installation.
2. conda create -y --name polap polap=0.4.3.7.x
3. conda activate polap
4. install polap-fmlrc using polap-conda-environment-fmlrc.yaml
5. conda install goshng::cflye

### using github

1. Follow the quick start available at 
   https://github.com/goshng/polap?tab=readme-ov-file#using-github-source
   to use the source code. The followings are the essential commands for the installation.
2. git clone https://github.com/goshng/polap.git
3. ln -s polap/src

## Configuration

### Data file

It looks for a CSV file in the following order:
1. The current folder where you execute this script
2. The folder in which this script resides

The CSV file should have its first line serving as a header, followed by
additional lines containing the options for an analysis. Example:
species,_taxon,_folder,_long,_short,_host,_ptgaul_genomesize,_compare_n,_compare_p,_compare_r,_polish_n,_polish_p,_random_seed,_ssh,_memory,_downsample,_inum,_table1,_table2,_mainfigure,dummy,_status
Anthoceros_agrestis-0,Anthoceros_agrestis,Anthoceros_agrestis,SRR10190639,SRR10250248,siepel,160000,20,10,10,5,5,157,siepel,16,10,0,T,F,F,dummy,done

The following descriptions outline the column items present in the CSV file.
The first three columns can be considered redundant, admittedly.

- species: index, <genus>_<species>-<inum>
- _taxon: species name, <genus>_<species>
- _folder: <genus>_<species>
- _long: long-read NCBI SRA accession
- _short: short-read NCBI SRA accession
- _host: the remote comuputer host name echoed by hostname command
- _ptgaul_genomesize: the genome size used by ptGAUL assembly
- _compare_n: the number of steps in the first stage of polap's disassemble
- _compare_p: the maximum rate in the first stage of polap's disassemble
- _compare_r: the number of steps in the second and third stages
- _polish_n: NA or any number
- _polish_p: NA or any number
- _random_seed: seed for random number generation
- _ssh: the ssh name to connect to the remote
- _memory: the maximum memory to use in polap's disassemble
- _downsample: the downsample coverage
- _inum: analysis number e.g., 2, 0, or 3
- _table1: T or F indicating the use of the analysis in the main table
- _table2: T/F indicating reference-based ptGAUL run
- _mainfigure: T or F indicating the use of the analysis as the main result
- dummy: dummy always
- _status: comment any string

### Variables

At the beginning of the script, several settings require initial setup before proceeding.

- Smain: some of the results go to the main not supplementary materials.
- _arg_default_target_dir: folder where you write a manuscript using tables and figures.

## Help

The following displays a simple usage.
$0 help <subcommand>
$0 <subcommand>

## Subcommands

- system: to display the Linux CPU and memory
- batch <species_folder> [number]: To execute all subcommands.
- sra <species_folder>: scp data files to the remote 
- send <species_folder>: scp data files to the remote 
- recover <species_folder>
- mkdir <species_folder>
- mkdir-all: to create empty folders for the analyses 
- refs <species_folder>: get-ptdna-from-ncbi <species_folder>
- getorganelle <species_folder>
- ptgaul <species_folder>: ptGAUL analysis on the data
- msbwt <species_folder>: prepare short-read polishing
- extract-ptdna-of-ptgaul <species_folder>: extract ptDNA from ptGAUL's result
- coverage <species_folder> [number]: ptGAUL's result
- infer <species_folder> [number] [--disassemble-simple-polishing]: assemble the ptDNA without a reference ptDNA
- infer2 <species_folder> [number] [--disassemble-simple-polishing]: assemble the ptDNA without a reference ptDNA
- check <species_folder> [number] [index:0] [--disassemble-simple-polishing]: assemble ptDNA with subsampling by comparing it with the ptGAUL assembly
- compare <species_folder>: compare ptDNAs with the ptGAUL ptDNA
- archive <species_folder>: archive the results
- get: fetch or recover the data from the archive from local or remote
- report [species_folder]: report the results
- maintable1: to create the main table for all the assemblies
- supptable1: to create the 3-stage tables
- suppfigure1: to create the assembly figure for the first stage
- suppfigure3: to create the figure for all the assemblies
- wga <species_folder>: whole-genome assembly
- simple-polish <species_folder>: simple polish the ptDNA using the short-read data
- subsample-polish <species_folder>: subsampling polish the ptDNA using the short-read data
- mauve <species_folder>: subsampling polish the ptDNA using the short-read data
- clean: [DANGER] delete the species folder
- write-config src/polap-data-v2.csv: write config
#
#
downsample2infer <outdir> <inum>
infer <outdir> <inum> [polish|simple]
check <outdir> <inmu> [polish|simple]
#

# call the following menu to create tables and figures
#
bandage1 <species_folder> 2
bandage2 <species_folder> 2
copy-figure: copy all figures to the target directory
table1 2
table1 0
table2 2
table2 0

for i in 0 2; do
  bash src/polap-data-v2.sh table1 \$i
  bash src/polap-data-v2.sh table2 \$i
done

for i in 1 2 3-infer; do
  for j in 0 2; do
    for k in on off; do
      bash src/polap-data-v2.sh suptable1 \$j \$k \$i
    done
  done
done

suptable1 2 off 1
suptable1 2 off 2
suptable1 2 off 3-infer
suptable1 2 on 1
suptable1 2 on 2
suptable1 2 on 3-infer

for i in 0 2; do
  bash src/polap-data-v2.sh supfigure1 \$i no off
  bash src/polap-data-v2.sh supfigure2 \$i off
done

supfigure1 2 no off
supfigure1 2 no on

supfigure2 2 off <- for sup
supfigure2 2 on

copy-figure: copy all figures to the target directory


# Main Table 1 and others
maintable1
supptable1
suppfigure1
suppfigure3
copy-figures

# Example
bash src/polap-data-v2.sh get all 2 off
bash src/polap-data-v2.sh maintable1 all 2 infer-1 1
bash src/polap-data-v2.sh supptable1 all 2 infer-1 x
bash src/polap-data-v2.sh supptable1 Eucalyptus_pauciflora 2 infer-1 x
bash src/polap-data-v2.sh suppfigure1 all 2 infer-1 1 1 yes
bash src/polap-data-v2.sh suppfigure3 2 infer-1 yes
bash src/polap-data-v2.sh copy-figures

# files

HEREDOC
)

help_message_get=$(
	cat <<HEREDOC

  Archive files are named in either one of these two forms:
  Eucalyptus_pauciflora-a.tar.gz or Eucalyptus_pauciflora-a-0.tar.gz
  If you confirm by responding with 'yes', this action will replace
  the existing folder named Eucalyptus_pauciflora with the extracted folder.
  If you respond with 'add', it will replace only the <inum> folder in the folder
  named Eucalyptus_pauciflora with the extracted folder's <inum> folder.
  If <inum> is -1, it will look for Eucalyptus_pauciflora-a.tar.gz.
  Ottherwise, it will use Eucalyptus_pauciflora-a-<inum>.tar.gz.
  So, using -1 of <inum>, you will replace the Eucalyptus_pauciflora folder.

  If option confirm is off, then it will get your response. You can skip the
  confirmation step by setting it to yes or add.
HEREDOC
)

help_message_maintable1=$(
	cat <<HEREDOC

  The main table for the summary of all analyses is created. In the CSV data
  file, only those with table1 column being T and <inum> being given as option
  are included in the table.
HEREDOC
)

help_message_supptable1=$(
	cat <<HEREDOC

  The supplementary table for the three-stage run is created. In the CSV data
  file, only those with table1 column being T and <inum> being given as option
  are included in the table.
HEREDOC
)

help_message_report=$(
	cat <<HEREDOC

  Create a set of report tables using markdown formatting to present data.
  It creates only this kind of file: 0/disassemble/infer-1/3-infer/summary1x.md
HEREDOC
)

help_message_remote_batch=$(
	cat <<HEREDOC

  Execute the batch job in the remote in dev mode only.

  per_host: run only jobs assigned to a particular remote host in the CSV.
HEREDOC
)

help_message_batch=$(
	cat <<HEREDOC

  The main batch command does the followings.
  - get the data
  - GetOrganelle assembly
  - msbwt short-read polishing
  - fetch reference from NCBI
  - ptGAUL assembly
  - ptGAUL's polishing
  - subsampling-based assembly
  <>: 0, meaning archiving it with the file name of -a-0.tar.gz
HEREDOC
)

help_message_archive=$(
	cat <<HEREDOC

  Archive the result.
  <inum>: -1 is default, meaning archiving it with the file name of -a.tar.gz
  <inum>: 0, meaning archiving it with the file name of -a-0.tar.gz
HEREDOC
)

declare -A _folder
declare -A _taxon
declare -A _mainfigure
declare -A _memory
declare -A _inum
declare -A _table1
declare -A _table2
declare -A _dummy
declare -A _downsample
declare -A _status
declare -A _host
declare -A _long
declare -A _short
declare -A _ptgaul_genomesize
declare -A _compare_n
declare -A _compare_p
declare -A _compare_r
declare -A _random_seed
declare -A _ssh
declare -A _polish_n
declare -A _polish_p

set +u

# Read the config files
read-a-tsv-file-into-associative-arrays() {
	# Define input TSV file
	csv_file="${PWD}/polap-data-v2.csv"
	if [[ ! -s "${csv_file}" ]]; then
		csv_file="${script_dir}/polap-data-v2.csv"
	fi

	# Read the TSV file (skip header)
	while IFS=$',' read -r species taxon folder long short host ptgaul_genomesize compare_n compare_p compare_r polish_n polish_p random_seed ssh memory downsample inum table1 table2 mainfigure dummy status; do
		# Skip header line
		[[ "$species" == "species" ]] && continue

		# Store in associative arrays
		_folder["$species"]="$folder"
		_taxon["$species"]="$taxon"
		_long["$species"]="$long"
		_short["$species"]="$short"
		_host["$species"]="$host"
		_ptgaul_genomesize["$species"]="$ptgaul_genomesize"
		_compare_n["$species"]="$compare_n"
		_compare_p["$species"]="$compare_p"
		_compare_r["$species"]="$compare_r"
		_polish_n["$species"]="$polish_n"
		_polish_p["$species"]="$polish_p"
		_random_seed["$species"]="$random_seed"
		_ssh["$species"]="$ssh"
		_memory["$species"]="$memory"
		_downsample["$species"]="$downsample"
		_inum["$species"]="$inum"
		_table1["$species"]="$table1"
		_table2["$species"]="$table2"
		_mainfigure["$species"]="$mainfigure"
		_dummy["$species"]="$dummy"
		_status["$species"]="$status"
	done <"$csv_file"
}

read-a-tsv-file-into-associative-arrays

keys_array=($(for key in "${!_long[@]}"; do echo "$key"; done | sort))

_arg1=${1:-arg1}
# Check if the species folder is provided
if [[ "${_arg1}" == "arg1" ]]; then
	echo "Usage: $0 <subcommand> [species_folder] ..."
	echo "${help_message}"
	exit 1
fi

_arg2=${2:-arg2}
if [[ "${_arg2}" != "arg2" ]]; then
	_arg2="${2%/}"
fi

_arg3=${3:-arg3}
_arg4=${4:-arg4}
_arg5=${5:-arg5}
_arg6=${6:-arg6}
_arg7=${7:-arg7}
_arg8=${8:-arg8}

################################################################################
# Part of genus_species
#
batch_genus_species() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local _brg_ref="${3:-off}"
	local _brg_yes="${4:-off}"
	local _brg_redo="${5:-off}"

	local do_simple_polishing="on"
	local target_index="${output_dir}-${isuffix}"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local ssh_remote="${_ssh["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	# rm -rf "${output_dir}"
	# tar -zxf "${output_dir}-a.tar.gz"
	# mv "${output_dir}-a" "${output_dir}"

	echo host: $(hostname)
	echo ssh-remote: $ssh_remote
	echo output: $output_dir
	echo inum: $extracted_inum
	echo species: $species_name
	echo random seed: $random_seed
	local long_data="${_media_dir}/${long_sra}.fastq.tar.gz"
	local short_data="${_media_dir}/${short_sra}.fastq.tar.gz"

	# if [[ "${_local_host}" == "$(hostname)" ]]; then
	# 	if [[ "${_brg_remote}" == "off" ]]; then
	# 		if [[ -s "${long_data}" ]]; then
	# 			echo long: $(du -h ${long_data})
	# 			echo short: $(du -h ${short_data})
	# 		else
	# 			echo long: $(du -h "${_media1_dir}/${long_sra}.fastq")
	# 			echo short1: $(du -h "${_media1_dir}/${short_sra}_1.fastq")
	# 			echo short2: $(du -h "${_media1_dir}/${short_sra}_2.fastq")
	# 		fi
	# 		return 0
	# 	else
	# 		echo "remote run"
	# 	fi
	# fi

	if [[ -s "${long_sra}.fastq.tar.gz" ]] ||
		[[ -s "${long_sra}.fq.tar.gz" ]] ||
		[[ -s "${long_sra}.fq" ]] ||
		[[ -s "${long_sra}.fastq" ]]; then
		echo "found: long-read data file"
	else
		echo "ERROR: No long-read data file"
		return 1
	fi

	if [[ "${_brg_yes}" == "off" ]]; then
		read -p "Do you want to execute the batch procedure? (y/N): " confirm
		case "$confirm" in
		[yY] | [yY][eE][sS])
			_brg_yes="on"
			;;
		*)
			echo "Batch procedure is canceled."
			;;
		esac
	fi

	if [[ "${_brg_yes}" == "off" ]]; then
		return 0
	fi

	# main
	if [[ -s "${output_dir}/polap.log" ]]; then
		echo "Folder ${output_dir} is not empty"
	else
		echo "Folder ${output_dir} is empty"
		recover_genus_species "${output_dir}" "${isuffix}"
	fi

	mkdir -p "${output_dir_i}"

	if [[ -s "${long_sra}.fastq" ]] &&
		[[ -s "${short_sra}_1.fastq" ]] &&
		[[ -s "${short_sra}_2.fastq" ]]; then
		_log_echo "Found: sequencing data"
	else
		mkdir_genus_species "${output_dir}" "${isuffix}"
		if [[ -s "${long_sra}.fastq" ]] &&
			[[ -s "${short_sra}_1.fastq" ]] &&
			[[ -s "${short_sra}_2.fastq" ]]; then
			_log_echo "Success: sequencing data"
		else
			_log_echo "Fail: sequencing data"
			return 1
		fi
	fi

	if [[ -d "${output_dir}/getorganelle" ]] &&
		find "${output_dir}/getorganelle" -maxdepth 1 -type f -name 'embplant_pt.*.gfa' -size +0c | grep -q .; then
		_log_echo "Found: GetOrganelle assembled ptDNA"
	else
		getorganelle_genus_species "${output_dir}"
		if find "${output_dir}/getorganelle" -maxdepth 1 -type f -name 'embplant_pt.*.gfa' -size +0c | grep -q .; then
			_log_echo "Success: GetOrganelle assembled ptDNA"
		else
			_log_echo "Fail: GetOrganelle assembled ptDNA"
			return 1
		fi
	fi

	if [[ -s "${output_dir}/msbwt/comp_msbwt.npy" ]]; then
		_log_echo "Found: FMLRC msbwt"
	else
		msbwt_genus_species "${output_dir}"
		if [[ -s "${output_dir}/msbwt/comp_msbwt.npy" ]]; then
			_log_echo "Success: FMLRC msbwt"
		else
			_log_echo "Fail: FMLRC msbwt"
			return 1
		fi
	fi

	# skip if the reference is not known

	if [[ "${_brg_ref}" == "on" ]]; then
		if [[ -s "${output_dir}/ptdna-reference.fa" ]]; then
			_log_echo "Found: reference ptDNA"
		else
			get-ptdna-from-ncbi_genus_species "${output_dir}"
			if [[ -s "${output_dir}/ptdna-reference.fa" ]]; then
				_log_echo "Success: reference ptDNA"
			else
				_log_echo "Fail: reference ptDNA"
				return 1
			fi
		fi

		if [[ -s "${output_dir}/ptgaul/flye_cpONT/assembly_graph.gfa" ]]; then
			_log_echo "Found: ptGAUL assembly"
		else
			ptgaul_genus_species "${output_dir}"
			if [[ -s "${output_dir}/ptgaul/flye_cpONT/assembly_graph.gfa" ]]; then
				_log_echo "Success: ptGAUL assembly"
			else
				_log_echo "Fail: ptGAUL assembly"
				return 1
			fi
		fi

		if [[ -s "${output_dir}/ptdna-ptgaul.fa" ]]; then
			_log_echo "Found: ptGAUL polished genome"
		else
			extract-ptdna-of-ptgaul_genus_species "${output_dir}"
			copy-ptdna-of-ptgaul_genus_species "${output_dir}"
			if [[ -s "${output_dir}/ptdna-ptgaul.fa" ]]; then
				_log_echo "Success: ptGAUL polished genome"
			else
				_log_echo "Fail: ptGAUL polished genome"
				return 1
			fi
		fi
	fi

	# delete if you want to redo the analysis
	if [[ "${_brg_redo}" == "on" ]]; then
		_log_echo "INFO: delete ${output_dir_i} to reanalyze the subsampling part."
		rm -rf "${output_dir_i}"
	fi

	if [[ -s "${output_dir_i}/disassemble/infer-1/pt.subsample-polishing.1.fa" ]]; then
		_log_echo "Found: infer case"
	else
		infer_genus_species "${output_dir}" "${isuffix}"
		if [[ -s "${output_dir_i}/disassemble/infer-1/pt.subsample-polishing.1.fa" ]]; then
			_log_echo "Success: infer case"
		else
			_log_echo "Fail: infer case"
			return 1
		fi
	fi

	if [[ "${do_simple_polishing}" == "on" ]]; then

		if [[ -s "${output_dir_i}/disassemble/infer-1/pt.simple-polishing.1.fa" ]]; then
			_log_echo "Found: infer case - simple polishing"
		else
			infer_genus_species "${output_dir}" "${isuffix}" simple
			if [[ -s "${output_dir_i}/disassemble/infer-1/pt.simple-polishing.1.fa" ]]; then
				_log_echo "Success: infer case - simple polishing"
			else
				_log_echo "Fail: infer case - simple polishing"
				_log_echo "  potenial error: fmlrc not enough memory"
				return 1
			fi
		fi
	fi

	if [[ -s "${output_dir}/ptgaul/flye_cpONT/assembly_graph.gfa" ]]; then

		if [[ -s "${output_dir_i}/disassemble/infer-1/pt.subsample-polishing.reference.aligned.1.fa" ]]; then
			_log_echo "Found: check case"
		else
			check_genus_species "${output_dir}" "${isuffix}"
			if [[ -s "${output_dir_i}/disassemble/infer-1/pt.subsample-polishing.reference.aligned.1.fa" ]]; then
				_log_echo "Success: check case"
			else
				_log_echo "Fail: check case"
				return 1
			fi
		fi

		if [[ "${do_simple_polishing}" == "on" ]]; then
			if [[ -s "${output_dir_i}/disassemble/infer-1/pt.simple-polishing.reference.aligned.1.fa" ]]; then
				_log_echo "Found: check case - simple polishing"
			else
				check_genus_species "${output_dir}" "${isuffix}" simple
				if [[ -s "${output_dir_i}/disassemble/infer-1/pt.simple-polishing.reference.aligned.1.fa" ]]; then
					_log_echo "Success: check case - simple polishing"
				else
					_log_echo "Fail: check case - simple polishing"
					_log_echo "  potenial error: fmlrc not enough memory"
					return 1
				fi
			fi
		fi

		if printf '%s\n' "${Stable5[@]}" | grep -qx "${output_dir}"; then
			if [[ -s "${output_dir_i}/disassemble/compare-1/pt.simple-polishing.1.fa" ]]; then
				_log_echo "Found: compare case"
			else
				compare_genus_species "${output_dir}" "${isuffix}"
				if [[ -s "${output_dir_i}/disassemble/compare-1/pt.simple-polishing.1.fa" ]]; then
					_log_echo "Success: compare case"
				else
					_log_echo "Fail: compare case"
					return 1
				fi
			fi
		fi
	else
		echo "No such file: ${output_dir}/ptgaul/flye_cpONT/assembly_graph.gfa"
		echo "  so, skip comparing the ptGAUL and the subsampling-based assembly"
	fi

}

# Gets the short- and long-read data and the archive
# Batch run on them
# Send the result back
# Clean up the data
remote-batch_genus_species_for() {
	local output_dir="$1"
	local _brg_inum="${2:-0}"

	# Gets the datasets
	local target_index="${output_dir}-${_brg_inum}"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"

	if [[ ! -v _long["$target_index"] ]]; then
		echo "INFO: no element for key $target_index"
		return 1
	fi

	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local ssh_remote="${_ssh["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	local _switch_ref="on"
	if [[ "${_table2["$target_index"]}" == "T" ]]; then
		_switch_ref="off"
	fi

	# Copy the processed dataset if not local
	if [[ "${_local_host}" != "$(hostname)" ]]; then
		scp -p ${_local_host}:$PWD/${output_dir}-a.tar.gz .
	fi

	local long_data="${_media_dir}/${long_sra}.fastq.tar.gz"
	local short_data="${_media_dir}/${short_sra}.fastq.tar.gz"

	if [[ ! -s "${long_sra}.fastq" ]]; then
		if ssh ${_local_host} "test -f ${long_data}"; then
			scp ${_local_host}:${long_data} .
		else
			scp ${_local_host}:"${_media1_dir}/${long_sra}.fastq" .
		fi
	fi

	if [[ ! -s "${short_sra}_1.fastq" ]]; then
		if ssh ${_local_host} "test -f ${short_data}"; then
			scp ${_local_host}:${short_data} .
		else
			scp ${_local_host}:"${_media1_dir}/${short_sra}_1.fastq" .
			scp ${_local_host}:"${_media1_dir}/${short_sra}_2.fastq" .
		fi
	fi

	batch_genus_species ${output_dir} ${_brg_inum} "${_switch_ref}" on on
	archive_genus_species ${output_dir}

	if [[ "${_local_host}" != "$(hostname)" ]]; then
		scp -p ${output_dir}-a.tar.gz ${_local_host}:$PWD/${output_dir}-a-${_brg_inum}.tar.gz
	fi

	clean_genus_species "${output_dir}" on

	return
}

remote-batch_genus_species() {
	local output_dir="${1:-all}"
	local _brg_inum="${2:-0}"
	local _brg_per_host="${3:-off}"

	if [[ "${_local_host}" == "$(hostname)" ]]; then
		echo host: $(hostname)
		echo "Run at the remote."
		return
	fi

	if [[ "${output_dir}" == "all" ]]; then
		for _v1 in "${Sall[@]}"; do
			if [[ "${_brg_per_host}" == "on" ]]; then
				local _target_index="${_v1}-0"
				local _host_remote="${_host["$_target_index"]}"
				if [[ "${_host_remote}" == "$(hostname)" ]]; then
					remote-batch_genus_species_for "${_v1}" "${_brg_inum}"
				fi
			else
				remote-batch_genus_species_for "${_v1}" "${_brg_inum}"
			fi
		done
	else
		remote-batch_genus_species_for "${output_dir}" "${_brg_inum}"
	fi
}

local-batch_genus_species_for() {
	local output_dir="$1"
	local _brg_inum="${2:-0}"

	# Gets the datasets
	local target_index="${output_dir}-0"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local ssh_remote="${_ssh["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	local _switch_ref="on"
	if [[ "${_table2["$target_index"]}" == "T" ]]; then
		_switch_ref="off"
	fi

	# Copy the processed dataset if not local
	mkdir -p backup
	timestamp=$(date +"%Y%m%d%H%M%S") # Get the current date and time
	archive_genus_species ${output_dir}
	cp -p ${output_dir}-a.tar.gz backup/${output_dir}-a-${timestamp}.tar.gz

	local long_data="${_media_dir}/${long_sra}.fastq.tar.gz"
	local short_data="${_media_dir}/${short_sra}.fastq.tar.gz"

	if [[ ! -s "${long_sra}.fastq" ]]; then
		if [[ -s "${long_data}" ]]; then
			cp "${long_data}" .
		elif [[ -s "${_media1_dir}/${long_sra}.fastq" ]]; then
			cp "${_media1_dir}/${long_sra}.fastq" .
		else
			echo "fetching ${long_sra} from NCBI ... it takes time ... be patient!"
			${_polap_cmd} x-ncbi-fetch-sra --sra ${long_sra}
		fi
	fi

	if [[ ! -s "${short_sra}.fastq" ]]; then
		if [[ -s "${short_data}" ]]; then
			cp "${short_data}" .
		elif [[ -s "${_media1_dir}/${short_sra}_1.fastq" ]]; then
			cp "${_media1_dir}/${short_sra}_1.fastq" .
			cp "${_media1_dir}/${short_sra}_2.fastq" .
		else
			echo "fetching ${short_sra} from NCBI ... it takes time ... be patient!"
			${_polap_cmd} x-ncbi-fetch-sra --sra ${short_sra}
		fi
	fi

	if [[ -d "${long_sra}" ]]; then
		rm -rf "${long_sra}"
	fi
	if [[ -d "${short_sra}" ]]; then
		rm -rf "${short_sra}"
	fi

	batch_genus_species ${output_dir} 2 "${_switch_ref}" on on
	batch_genus_species ${output_dir} 0 "${_switch_ref}" on on

	return
}

local-batch_genus_species() {
	local output_dir="${1:-all}"
	local _brg_inum="${2:-0}"

	if [[ "${_local_host}" != "$(hostname)" ]]; then
		echo host: $(hostname)
		echo "Run at the local."
		return
	fi

	if [[ "${output_dir}" == "all" ]]; then
		for _v1 in "${Sall[@]}"; do
			local-batch_genus_species_for "${_v1}" "${_brg_inum}"
		done
	else
		local-batch_genus_species_for "${output_dir}" "${_brg_inum}"
	fi
}

send-data_genus_species() {
	local output_dir="$1"
	local inum="${2:-0}"
	local target_index="${output_dir}-${inum}"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local ssh_name="${_ssh["$target_index"]}"

	echo host: $(hostname)
	echo remote dir: ${ssh_name}:$PWD
	echo output: ${output_dir}
	echo target_index: ${target_index}
	echo species: $species_name
	echo long: $long_sra
	echo short: $short_sra

	local long_data="${_media_dir}/${long_sra}.fastq.tar.gz"
	local short_data="${_media_dir}/${short_sra}.fastq.tar.gz"

	if [[ "${_local_host}" == "$(hostname)" ]]; then
		if [[ -s "${long_data}" ]]; then
			echo long: $(du -h ${long_data})
			echo short: $(du -h ${short_data})
		else
			echo long: $(du -h "${_media1_dir}/${long_sra}.fastq")
			echo short1: $(du -h "${_media1_dir}/${short_sra}_1.fastq")
			echo short2: $(du -h "${_media1_dir}/${short_sra}_2.fastq")
		fi
	else
		echo "You can send the data from the central host: ${_local_host}"
		return 1
	fi

	if [[ "${_local_host}" == "$(hostname)" ]]; then

		ssh ${ssh_name} "mkdir -p $PWD/${output_dir}"

		if [[ -s "${long_data}" ]]; then
			scp "${long_data}" ${ssh_name}:$PWD/
		else
			scp \
				"${_media1_dir}/${long_sra}.fastq" \
				${ssh_name}:$PWD/
		fi

		if [[ -s "${short_data}" ]]; then
			scp "${short_data}" ${ssh_name}:$PWD/
		else
			scp \
				"${_media1_dir}/${short_sra}_1.fastq" \
				"${_media1_dir}/${short_sra}_2.fastq" \
				${ssh_name}:$PWD/
		fi
	else
		echo "ERROR: run at the local host."
	fi
}

# send the input data to the remote
send_genus_species() {
	local output_dir="$1"
	local inum="${2:-0}"
	local target_index="${output_dir}-${inum}"

	if [[ "${_local_host}" == "$(hostname)" ]]; then
		archive_genus_species "${output_dir}"
		scp "${output_dir}-a.tar.gz" ${_ssh["${target_index}"]}:$PWD/
	else
		echo "ERROR: run at the local host."
	fi

	send-data_genus_species "${output_dir}" "${inum}"
}

send-archive_genus_species() {
	local output_dir="$1"
	local inum="${2:-0}"
	local target_index="${output_dir}-${inum}"

	if [[ "${_local_host}" == "$(hostname)" ]]; then
		archive_genus_species "${output_dir}"
		scp "${output_dir}-a.tar.gz" ${_ssh["${target_index}"]}:$PWD/
	else
		echo "ERROR: run at the local host."
	fi
}

# extract the archive at the remote
recover_genus_species() {
	local output_dir="$1"

	if [[ -s "${output_dir}-a.tar.gz" ]]; then
		echo "Deleting ${output_dir} ..."
		rm -rf "${output_dir}"
		tar -zxf "${output_dir}-a.tar.gz"
		mv "${output_dir}-a" "${output_dir}"
		echo "  we have recreated ${output_dir}!"
		echo "Next: $0 mkdir ${output_dir}"
	else
		mkdir -p "${output_dir}"
	fi
}

# create input files
mkdir_genus_species() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local target_index="${output_dir}-${isuffix}"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local ssh_remote="${_ssh["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	local long_data="${_media_dir}/${long_sra}.fastq.tar.gz"
	local short_data="${_media_dir}/${short_sra}.fastq.tar.gz"

	echo "create ${output_dir} ..."
	mkdir -p "${output_dir}/timing"
	if [[ -s "${long_sra}.fastq" ]]; then
		echo "  found: long SRA: ${long_sra}.fastq"
	else
		if tar -zxf "$(basename ${long_data})"; then
			rm "$(basename ${long_data})"
			echo "Extraction: ${long_sra} deleted."
		fi
	fi
	if [[ -s "${short_sra}_1.fastq" ]] && [[ -s "${short_sra}_2.fastq" ]]; then
		echo "  found: short SRA1: ${short_sra}_1.fastq"
		echo "  found: short SRA1: ${short_sra}_2.fastq"
	else
		if tar -zxf "$(basename ${short_data})"; then
			rm "$(basename ${short_data})"
			echo "Extraction: ${short_sra} deleted."
		fi
	fi

	echo "Next: $0 refs ${output_dir} [number] to download reference ptDNAs from NCBI"
	echo "Next: $0 coverage ${output_dir} [number] to overview your data"
}

system_genus_species() {
	echo "Host: $(hostname)"
	echo "  CPU: lscpu"
	lscpu | grep Model | head -1
	echo "  Memory: free -h"
	free -h | grep Mem
}

sra_genus_species() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local target_index="${output_dir}-${isuffix}"

	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"

	echo "create ${output_dir} ..."
	mkdir -p "${output_dir}/timing"
	run-polap-ncbitools fetch sra "$long_sra"
	run-polap-ncbitools fetch sra "$short_sra"

	echo "Next: $0 refs ${output_dir} [number] to download reference ptDNAs from NCBI"
	echo "Next: $0 coverage ${output_dir} [number] to overview your data"
}

get-ptdna-from-ncbi_genus_species() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local target_index="${output_dir}-${isuffix}"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local ssh_remote="${_ssh["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	if [[ "${output_dir}" == "Juncus_inflexus" ]]; then
		species_name="Juncus effusus"
		echo "No ptDNA for ${output_dir}, so we use ${species_name}"
	fi

	if [[ -s "${output_dir}/ptdna-reference.fa" ]]; then
		echo "found: ptDNA reference: ${output_dir}/ptdna-reference.fa"
	else

		${_polap_cmd} get-mtdna \
			--plastid \
			--species "${species_name}" \
			-o ${output_dir}

		if [[ -s "${output_dir}/00-bioproject/2-mtdna.fasta" ]]; then
			echo "copy ${output_dir}/ptdna-reference.fa"
			cp -p "${output_dir}/00-bioproject/2-mtdna.fasta" \
				"${output_dir}/ptdna-reference.fa"
		else
			echo "No such file: ${output_dir}/ptdna-reference.fa"
			local _genus_name=$(echo ${species_name} | awk '{print $1}')
			echo "  trying to search NCBI plastid genomes for genus name only: ${_genus_name}"
			${_polap_cmd} get-mtdna \
				--plastid \
				--species "${_genus_name}" \
				-o ${output_dir}
			if [[ -s "${output_dir}/00-bioproject/2-mtdna.fasta" ]]; then
				echo "copy ${output_dir}/ptdna-reference.fa"
				cp -p "${output_dir}/00-bioproject/2-mtdna.fasta" \
					"${output_dir}/ptdna-reference.fa"
			else
				echo "  we could not find one even in the genus level."
				echo "No such file: ${output_dir}/ptdna-reference.fa"
			fi
		fi
	fi
}

copy-ptdna-of-ncbi-as-reference_genus_species() {
	local output_dir="$1"

	echo "copy ${output_dir}/ptdna-reference.fa"
	cp -p "${output_dir}/00-bioproject/2-mtdna.fasta" \
		"${output_dir}/ptdna-reference.fa"
}

clean_genus_species() {
	local output_dir="$1"
	local _brg_yes="${2:-off}"
	local target_index="${output_dir}-0"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local ssh_remote="${_ssh["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	local long_data="${_media_dir}/${long_sra}.fastq.tar.gz"
	local short_data="${_media_dir}/${short_sra}.fastq.tar.gz"

	if [[ "${_brg_yes}" == "off" ]]; then
		read -p "Do you want to delete it? (y/N): " confirm
		case "$confirm" in
		[yY] | [yY][eE][sS])
			_brg_yes="on"
			;;
		*)
			echo "Batch procedure is canceled."
			;;
		esac
	fi

	if [[ "${_brg_yes}" == "off" ]]; then
		return 0
	fi

	# Check if the directory exists
	if [[ -d "$output_dir" ]]; then
		echo "Directory '$output_dir' exists."

		echo "cleaning ${output_dir} ..."
		rm -rf "${output_dir}"
		rm -rf "${output_dir}-a"
		rm -f "${output_dir}"-a*.tar.gz
		rm -f "${long_sra}.fastq"
		rm -f "${short_sra}"_?.fastq
		rm -f "${long_sra}.fastq.tar.gz"
		rm -f "${short_sra}.fastq.tar.gz"

	else
		echo "Directory '$output_dir' does not exist."
		echo "  but, deleting its fastq.tar.gz files ... if any"
		rm -f "${long_sra}.fastq.tar.gz"
		rm -f "${short_sra}.fastq.tar.gz"
	fi

}

getorganelle_genus_species_for() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local target_index="${output_dir}-${isuffix}"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local ssh_remote="${_ssh["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	if [[ "${_local_host}" == "$(hostname)" ]]; then
		decompress-short_genus_species "${output_dir}"
	else
		echo "  run at the local host for the decompressing the data"
	fi
	echo "${short_sra}_1.fastq"

	if [[ -s "${short_sra}_1.fastq" ]] &&
		[[ -s "${short_sra}_2.fastq" ]]; then
		# Initialize Conda
		source $HOME/miniconda3/etc/profile.d/conda.sh
		conda activate getorganelle

		# Example
		# -1 Arabidopsis_simulated.1.fq.gz -2 Arabidopsis_simulated.2.fq.gz -t 1 -o Arabidopsis_simulated.plastome -F embplant_pt -R 10
		# default
		# -R 15 \
		mkdir -p "${output_dir}/timing"
		command time -v get_organelle_from_reads.py \
			-1 ${short_sra}_1.fastq \
			-2 ${short_sra}_2.fastq \
			-o ${output_dir}-getorganelle \
			-t 24 \
			-F embplant_pt \
			2>${output_dir}/timing/timing-getorganelle.txt

		rsync -aq \
			--max-size=5M \
			"${output_dir}-getorganelle/" \
			"${output_dir}/getorganelle/"
		rm -rf "${output_dir}-getorganelle"

		conda deactivate

		if [[ "${_local_host}" == "$(hostname)" ]]; then
			delete-short_genus_species "${output_dir}"
		else
			echo "  run at the local host for deleting the data"
		fi
	else
		echo "  no run getorganelle"
	fi
}

decompress-short_genus_species() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local target_index="${output_dir}-${isuffix}"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local ssh_remote="${_ssh["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	local short_data="${_media_dir}/${short_sra}.fastq.tar.gz"

	if [[ -s "${short_data}" ]]; then
		echo "decompressing the short-read data ..."
		tar zxf "${short_data}"
	else
		echo "copying the short-read data ..."
		cp -p "/media/h1/sra/${short_sra}_1.fastq" .
		cp -p "/media/h1/sra/${short_sra}_2.fastq" .
	fi

	if [[ -s "${short_sra}_1.fastq" ]] &&
		[[ -s "${short_sra}_2.fastq" ]]; then
		echo short: $short_sra
	else
		echo short: no such fastq file: $short_sra
		echo run $0 send-data-to at ${_local_host}
	fi
}

delete-short_genus_species() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local target_index="${output_dir}-${isuffix}"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local ssh_remote="${_ssh["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	if [[ -s "${short_sra}_1.fastq" ]] &&
		[[ -s "${short_sra}_2.fastq" ]]; then
		echo deleteing short-read data: $short_sra
		rm -f "${short_sra}_1.fastq"
		rm -f "${short_sra}_2.fastq"
	else
		echo short: no such fastq file: $short_sra
		echo run $0 send-data-to at ${_local_host}
	fi
}

getorganelle_genus_species() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local target_index="${output_dir}-${isuffix}"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local ssh_remote="${_ssh["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	if [[ "${output_dir}" == "default" ]]; then
		for _v1 in "${Sall[@]}"; do
			if [[ -d "${_v1}" ]]; then
				getorganelle_genus_species_for "${_v1}" "${isuffix}"
			else
				echo "No such folder: ${_v1}"
			fi
		done
	else
		getorganelle_genus_species_for "${output_dir}" "${isuffix}"
	fi
}

ptgaul_genus_species() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local target_index="${output_dir}-${isuffix}"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local ssh_remote="${_ssh["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"
	local extracted_ptgaul_genomesize="${_ptgaul_genomesize["$target_index"]}"

	mkdir -p "${output_dir}/timing"
	command time -v bash ${script_dir}/ptGAUL1.sh \
		-o ${output_dir}-ptgaul \
		-r ${output_dir}/ptdna-reference.fa \
		-g "${extracted_ptgaul_genomesize}" \
		-l ${long_sra}.fastq \
		-t 24 \
		2>${output_dir}/timing/timing-ptgaul.txt

	mv ${output_dir}-ptgaul/result_3000 ${output_dir}/ptgaul
	rm -rf "${output_dir}-ptgaul"

	msbwt_genus_species "${output_dir}"
}

msbwt_genus_species() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local target_index="${output_dir}-${isuffix}"

	local short_sra="${_short["$target_index"]}"

	mkdir -p "${output_dir}/timing"
	echo "This polishing preparation takes long ..."
	command time -v ${_polap_cmd} prepare-polishing \
		-a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
		-o ${output_dir} \
		2>${output_dir}/timing/timing-prepare-polishing.txt
}

# Extraction of ptDNA from the assembly: a first try
extract-ptdna-of-ptgaul_genus_species() {
	local output_dir="$1"

	mkdir -p "${output_dir}/timing"
	# extract ptGAUL result
	echo "extract ptDNA from the ptGAUL result with fmlrc polishing"
	command time -v ${_polap_cmd} disassemble ptgaul \
		-o ${output_dir} \
		2>${output_dir}/timing/timing-ptgaul-polishing.txt
	echo "use extract-ptdna-of-ptgaul2 <species_folder> if not working"
	copy-ptdna-of-ptgaul_genus_species
}

# Extraction of ptDNA from the assembly: another try
# not using it much
extract-ptdna-of-ptgaul2_genus_species() {
	local output_dir="$1"

	mkdir -p "${output_dir}/timing"
	# extract ptGAUL result
	echo "extract ptDNA from the ptGAUL result with fmlrc polishing"
	${_polap_cmd} disassemble ptgaul 2 \
		-o ${output_dir}
	command time -v ${_polap_cmd} disassemble ptgaul 3 \
		-v -v -v \
		-o ${output_dir} \
		2>${output_dir}/timing/timing-ptgaul-polishing.txt

	_outdir="${output_dir}/ptgaul/flye_cpONT/ptdna"
	_arg_final_assembly="${_outdir}/pt.1.fa"
	if [[ -s "${_arg_final_assembly}" ]]; then
		# copy ptGAUL result
		echo "copy ${output_dir}/ptdna-ptgaul.fa"
		cp -p ${_arg_final_assembly} ${output_dir}/ptdna-ptgaul.fa
	else
		echo "No such file: ${_arg_final_assembly}"
	fi
}

copy-ptdna-of-ptgaul_genus_species() {
	local output_dir="$1"

	# copy ptGAUL result
	echo "copy ${output_dir}/ptdna-ptgaul.fa"
	_outdir="${output_dir}/ptgaul/flye_cpONT/ptdna"
	_arg_final_assembly="${_outdir}/pt.1.fa"
	cp -pu ${_arg_final_assembly} ${output_dir}/ptdna-ptgaul.fa
}

# Case of the infer menu
# no --disassemble-c
coverage_genus_species() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local target_index="${output_dir}-${isuffix}"

	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	local extracted_downsample="${_downsample["$target_index"]}"

	local i=1
	mkdir -p "${output_dir_i}"
	echo "Analysis (coverage): ${output_dir} at ${isuffix}"
	echo "  downsample=${extracted_downsample}x"

	local _d_i="infer-$i"
	local _x_i="coverage-$i"
	local _stages="--stages-include 0"

	# NOTE: "${_stages}" is a bug.
	# use it without quotations.
	command time -v ${_polap_cmd} disassemble \
		${_stages} \
		--downsample ${extracted_downsample} \
		-i ${extracted_inum} \
		-o ${output_dir} \
		-l ${long_sra}.fastq \
		-a ${short_sra}_1.fastq \
		-b ${short_sra}_2.fastq \
		--disassemble-i "${_d_i}" \
		--random-seed "${random_seed}" \
		2>"${output_dir_i}/timing-${_x_i}.txt"

	echo "See ${output_dir_i}/{lx.txt,sx.txt} for the basic statistics of your data"
	echo "Next: $0 infer ${output_dir} [number] to assemble the ptDNA genome"
}

delete-disassemble_genus_species() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local target_index="${output_dir}-${isuffix}"

	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	local i=1
	echo "Deleting disassemble: ${output_dir_i}/disassemble"
	read -p "Do you want to delete it? (y/N): " confirm

	case "$confirm" in
	[yY] | [yY][eE][sS])
		rm -rf "${output_dir_i}/disassemble"
		;;
	*)
		echo "Deletion canceled."
		;;
	esac

	echo "See ${output_dir_i}/{lx.txt,sx.txt} for the basic statistics of your data"
	echo "Next: $0 infer ${output_dir} [number] to assemble the ptDNA genome"
}

# Case of the infer menu
# no --disassemble-c
#
# infer <outdir> <inum> [default|polishing|simple]
#
# default: all stages
# polish: stage 3 only
# simple: stage 3 only but simple polishing
infer_genus_species() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local simple_polishing="${3:-default}"
	local target_index="${output_dir}-${isuffix}"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local ssh_remote="${_ssh["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	local i=0
	local n
	local p
	IFS=':' read -r -a extracted_array_n <<<"${_compare_n["$target_index"]}"
	IFS=':' read -r -a extracted_array_p <<<"${_compare_p["$target_index"]}"
	local extracted_r="${_compare_r["$target_index"]}"
	local extracted_memory="${_memory["$target_index"]}"
	local extracted_downsample="${_downsample["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	mkdir -p "${output_dir_i}"

	for n in "${extracted_array_n[@]}"; do
		for p in "${extracted_array_p[@]}"; do
			i=$((i + 1))
			echo "Analysis (inference): ${output_dir} at ${isuffix}"
			echo "($i) n=$n, p=$p, r=${extracted_r} memory=${extracted_memory}G, downsample=${extracted_downsample}x"

			local _d_i="infer-$i"
			local _x_i="infer-$i"
			local _s_i="subsample-polish"
			local _stages="--stages-include 0-3"
			if [[ "${simple_polishing}" == "simple" ]]; then
				simple_polishing="--disassemble-simple-polishing"
				_stages="--stages-include 3"
				_s_i="simple-polish"
			elif [[ "${simple_polishing}" == "polish" ]]; then
				simple_polishing=""
				_stages="--stages-include 3"
				_s_i="subsample-polish-only"
			fi

			# NOTE: "${_stages}" is a bug.
			# use it without quotations.
			command time -v ${_polap_cmd} disassemble \
				${_stages} \
				--downsample ${extracted_downsample} \
				-i ${extracted_inum} \
				-o ${output_dir} \
				-l ${long_sra}.fastq \
				-a ${short_sra}_1.fastq \
				-b ${short_sra}_2.fastq \
				${simple_polishing} \
				--disassemble-i "${_d_i}" \
				--disassemble-n $n \
				--disassemble-p $p \
				--disassemble-r ${extracted_r} \
				--disassemble-memory ${extracted_memory} \
				--disassemble-alpha 1.0 \
				--random-seed "${random_seed}" \
				2>"${output_dir_i}/timing-${_x_i}-${_s_i}.txt"

			if [[ -d "${output_dir_i}/disassemble/${_d_i}/3" ]]; then
				rm -rf "${output_dir_i}/disassemble/${_d_i}/3-infer"
				mv "${output_dir_i}/disassemble/${_d_i}/3" \
					"${output_dir_i}/disassemble/${_d_i}/3-infer"
			fi
		done
	done
	echo "Next: $0 check ${output_dir} [number] to compare the assembly of the ptDNA genome with ptGAUL's"
}

downsample2infer_genus_species() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local simple_polishing="${3:-default}"
	local target_index="${output_dir}-${isuffix}"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local ssh_remote="${_ssh["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	local i=0
	local n
	local p
	IFS=':' read -r -a extracted_array_n <<<"${_compare_n["$target_index"]}"
	IFS=':' read -r -a extracted_array_p <<<"${_compare_p["$target_index"]}"
	local extracted_r="${_compare_r["$target_index"]}"
	local extracted_memory="${_memory["$target_index"]}"
	local extracted_downsample="${_downsample["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	mkdir -p "${output_dir_i}"

	for n in "${extracted_array_n[@]}"; do
		for p in "${extracted_array_p[@]}"; do
			i=$((i + 1))
			echo "Analysis (inference): ${output_dir} at ${isuffix}"
			echo "($i) n=$n, p=$p, r=${extracted_r} memory=${extracted_memory}G, downsample=${extracted_downsample}x"

			local _d_i="infer-$i"
			local _x_i="infer-$i"
			local _s_i="downsample"

			# NOTE: "${_stages}" is a bug.
			# use it without quotations.
			command time -v ${_polap_cmd} disassemble downsample \
				--downsample ${extracted_downsample} \
				--disassemble-i ${_d_i} \
				-i ${extracted_inum} \
				-o ${output_dir} \
				-l ${long_sra}.fastq \
				-a ${short_sra}_1.fastq \
				-b ${short_sra}_2.fastq \
				--random-seed "${random_seed}" \
				2>"${output_dir_i}/timing-${_x_i}-${_s_i}.txt"

		done
	done
	echo "Next: $0 check ${output_dir} [number] to compare the assembly of the ptDNA genome with ptGAUL's"
}

# what is this? fix the sampling size?
# fix: alpha and beta
infer2_genus_species() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local simple_polishing="${3:-default}"
	local target_index="${output_dir}-${isuffix}"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local ssh_remote="${_ssh["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	local i=0
	local n
	local p
	IFS=':' read -r -a extracted_array_n <<<"${_compare_n["$target_index"]}"
	IFS=':' read -r -a extracted_array_p <<<"${_compare_p["$target_index"]}"
	local extracted_r="${_compare_r["$target_index"]}"
	local extracted_memory="${_memory["$target_index"]}"
	local extracted_downsample="${_downsample["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	mkdir -p "${output_dir_i}"

	for n in "${extracted_array_n[@]}"; do
		for p in "${extracted_array_p[@]}"; do
			i=$((i + 2))
			echo "Analysis (inference): ${output_dir} at ${isuffix}"
			echo "($i) n=$n, p=$p, r=${extracted_r} memory=${extracted_memory}G, downsample=${extracted_downsample}x"

			local _d_i="infer-$i"
			local _x_i="infer-$i"
			local _s_i="subsample-polish"
			local _stages="--stages-include 0-3"
			if [[ "${simple_polishing}" != "default" ]]; then
				_stages="--stages-include 3"
				_s_i="simple-polish"
			fi

			# NOTE: "${_stages}" is a bug.
			# use it without quotations.
			command time -v ${_polap_cmd} disassemble \
				${_stages} \
				--downsample ${extracted_downsample} \
				-i ${extracted_inum} \
				-o ${output_dir} \
				-l ${long_sra}.fastq \
				-a ${short_sra}_1.fastq \
				-b ${short_sra}_2.fastq \
				${simple_polishing} \
				--disassemble-i "${_d_i}" \
				--disassemble-n $n \
				--disassemble-p $p \
				--disassemble-r ${extracted_r} \
				--disassemble-memory ${extracted_memory} \
				--disassemble-alpha 1.5 \
				--disassemble-beta 0.05 \
				--random-seed "${random_seed}" \
				2>"${output_dir_i}/timing-${_x_i}-${_s_i}.txt"

			if [[ -d "${output_dir_i}/disassemble/${_d_i}/3" ]]; then
				rm -rf "${output_dir_i}/disassemble/${_d_i}/3-infer"
				mv "${output_dir_i}/disassemble/${_d_i}/3" \
					"${output_dir_i}/disassemble/${_d_i}/3-infer"
			fi
		done
	done
	echo "Next: $0 check ${output_dir} [number] to compare the assembly of the ptDNA genome with ptGAUL's"
}

# Case of the check menu
# --disassemble-c
# --disassemble-align-reference
# --disassemble-simple-polishing
#
# check <outdir> <inum> [simple|polish]
#
check_genus_species() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local simple_polishing="${3:-default}"
	local icount="${3:-0}" # not used any more
	local target_index="${output_dir}-${isuffix}"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local ssh_remote="${_ssh["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	local i="${icount}"
	local n
	local p
	IFS=':' read -r -a extracted_array_n <<<"${_compare_n["$target_index"]}"
	IFS=':' read -r -a extracted_array_p <<<"${_compare_p["$target_index"]}"
	local extracted_r="${_compare_r["$target_index"]}"
	local extracted_memory="${_memory["$target_index"]}"
	local extracted_downsample="${_downsample["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	mkdir -p "${output_dir_i}"

	for n in "${extracted_array_n[@]}"; do
		for p in "${extracted_array_p[@]}"; do
			i=$((i + 1))
			echo "Analysis (check): ${output_dir} at ${isuffix}"
			echo "($i) n=$n, p=$p, memory=${extracted_memory}G, downsample=${extracted_downsample}x"

			local _d_i="infer-$i"
			local _x_i="check-$i"
			local _s_i="subsample-polish"
			local _stages="--stages-include 3"
			if [[ "${simple_polishing}" == "simple" ]]; then
				simple_polishing="--disassemble-simple-polishing"
				_s_i="simple-polish"
			elif [[ "${simple_polishing}" == "polish" ]]; then
				simple_polishing=""
				_stages="--stages-include 3"
				_s_i="subsample-polish-only"
			else
				simple_polishing=""
			fi

			# "${_stages}" \
			command time -v ${_polap_cmd} disassemble \
				${_stages} \
				--downsample ${extracted_downsample} \
				-i ${extracted_inum} \
				-o ${output_dir} \
				-l ${long_sra}.fastq \
				-a ${short_sra}_1.fastq \
				-b ${short_sra}_2.fastq \
				--disassemble-c ${output_dir}/ptdna-ptgaul.fa \
				--disassemble-align-reference \
				${simple_polishing} \
				--disassemble-i "${_d_i}" \
				--disassemble-n $n \
				--disassemble-p $p \
				--disassemble-r ${extracted_r} \
				--disassemble-memory ${extracted_memory} \
				--disassemble-alpha 1.0 \
				--random-seed "${random_seed}" \
				2>"${output_dir_i}/timing-${_x_i}-${_s_i}.txt"

			# compare the results

			if [[ -d "${output_dir_i}/disassemble/${_d_i}/3" ]]; then
				rm -rf "${output_dir_i}/disassemble/${_d_i}/3-check"
				mv "${output_dir_i}/disassemble/${_d_i}/3" \
					"${output_dir_i}/disassemble/${_d_i}/3-check"
			fi

			local mauve_dir="${output_dir_i}/mauve/${i}"
			local blast_dir="${output_dir_i}/blast/${i}"
			local mafft_dir="${output_dir_i}/mafft/${i}"
			mkdir -p "${mauve_dir}"
			mkdir -p "${blast_dir}"
			mkdir -p "${mafft_dir}"
			if [[ -s "${output_dir_i}/disassemble/${_d_i}/pt.subsample-polishing.reference.aligned.1.fa" ]]; then

				${_polap_cmd} mauve-mtdna -a "${output_dir}/ptdna-ptgaul.fa" \
					-b "${output_dir_i}/disassemble/${_d_i}/pt.subsample-polishing.reference.aligned.1.fa" \
					-o "${mauve_dir}" \
					>"${mauve_dir}/log.txt"
				echo "see ${mauve_dir}/log.txt"
				cat "${mauve_dir}/log.txt"

				${_polap_cmd} compare2ptdna -a "${output_dir}/ptdna-ptgaul.fa" \
					-b "${output_dir_i}/disassemble/${_d_i}/pt.subsample-polishing.reference.aligned.1.fa" \
					-o "${blast_dir}"
				echo "see ${blast_dir}/pident.txt"
				cat "${blast_dir}/pident.txt"

				${_polap_cmd} mafft-mtdna -a "${output_dir}/ptdna-ptgaul.fa" \
					-b "${output_dir_i}/disassemble/${_d_i}/pt.subsample-polishing.reference.aligned.1.fa" \
					-o "${mafft_dir}" \
					>"${mafft_dir}/log.txt"
				echo "see ${mafft_dir}/pident.txt"
				cat "${mafft_dir}/pident.txt"
			else
				echo "ERROR: no such file: ${output_dir_i}/disassemble/${_d_i}/pt.subsample-polishing.reference.aligned.1.fa"
			fi

		done
	done
}

# Case of the compare menu
# --disassemble-c
# --no-disassemble-align-reference
compare_genus_species() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local simple_polishing="${3:-default}"
	local target_index="${output_dir}-${isuffix}"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local ssh_remote="${_ssh["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	local i=0
	local n
	local p
	IFS=':' read -r -a extracted_array_n <<<"${_compare_n["$target_index"]}"
	IFS=':' read -r -a extracted_array_p <<<"${_compare_p["$target_index"]}"
	local extracted_r="${_compare_r["$target_index"]}"
	local extracted_memory="${_memory["$target_index"]}"
	local extracted_downsample="${_downsample["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	mkdir -p "${output_dir_i}"

	for n in "${extracted_array_n[@]}"; do
		for p in "${extracted_array_p[@]}"; do
			i=$((i + 1))
			echo "Analysis (check stage 1 or compare): ${output_dir} at ${isuffix}"
			echo "($i) n=$n, p=$p, memory=${extracted_memory}G, downsample=${extracted_downsample}x"

			local _d_i
			_d_i="compare-$i"
			local _stages="--stages-include 1-3" # to avoid another down-sampling

			command time -v ${_polap_cmd} disassemble \
				${_stages} \
				--downsample ${extracted_downsample} \
				-i ${extracted_inum} \
				-o ${output_dir} \
				-l ${long_sra}.fastq \
				-a ${short_sra}_1.fastq \
				-b ${short_sra}_2.fastq \
				--disassemble-c ${output_dir}/ptdna-ptgaul.fa \
				--disassemble-i "${_d_i}" \
				--disassemble-n $n \
				--disassemble-p $p \
				--disassemble-r ${extracted_r} \
				--disassemble-memory ${extracted_memory} \
				--disassemble-alpha 1.0 \
				--random-seed "${random_seed}" \
				2>${output_dir_i}/timing-${_d_i}.txt
		done
	done
}

best-genus_species() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local simple_polishing="${3:-default}"
	local target_index="${output_dir}-${isuffix}"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local ssh_remote="${_ssh["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	${_polap_cmd} disassemble best \
		-i ${extracted_inum} \
		-o ${output_dir} \
		--disassemble-i infer-1

}

bandage1_genus_species() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local simple_polishing="${3:-default}"
	local target_index="${output_dir}-${isuffix}"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local ssh_remote="${_ssh["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"
	local images=()
	local captions=()

	${_polap_cmd} disassemble bandage \
		-i ${extracted_inum} \
		-o ${output_dir} \
		--disassemble-i infer-1
	local _gfa_polap="${output_dir}/${isuffix}/disassemble/infer-1/pt.1.gfa"
	local _png_polap="${output_dir}/${isuffix}/disassemble/infer-1/pt.1.png"

	images+=(${_png_polap})
	captions+=("Polap")

	local _gfa_ptgaul="${output_dir}/ptgaul/flye_cpONT/assembly_graph.gfa"
	local _png_ptgaul="${output_dir}/ptgaul.png"
	${_polap_cmd} bandage png \
		${_gfa_ptgaul} \
		${_png_ptgaul}
	images+=(${_png_ptgaul})
	captions+=("ptGAUL")

	echo "ptGAUL ptDNA gfa: ${_gfa_ptgaul}"
	echo "ptGAUL ptDNA png: ${_png_ptgaul}"

	local _gfa_getorganelle=$(find "${output_dir}/getorganelle" -type f -name 'embplant_pt*.gfa' | head -n 1)
	local _png_getorganelle="${output_dir}/getorganelle.png"
	${_polap_cmd} bandage png \
		${_gfa_getorganelle} \
		${_png_getorganelle}
	images+=(${_png_getorganelle})
	captions+=("GetOrganelle")

	echo "GetOrganelle ptDNA gfa: ${_gfa_getorganelle}"
	echo "GetOrganelle ptDNA png: ${_png_getorganelle}"

	local _figure_md="figure1.md"

	# echo "Warning: no such file: ${_gfa}"
	# echo "FIX QT5 problem:"
	# problem solved
	# echo "export QT_QPA_PLATFORM=offscreen"
	# Output Markdown file
	#
	output="${_figure_md}"

	# Start writing to the markdown file
	cat <<EOF >"$output"
---
title: "${species_name}"
geometry: margin=1in
---

# Image Grid

| Three Organelle genome assemblies |  |  |
|-----------------|-----------------|-----------------|
EOF

	# Generate the image table with subcaptions
	count=0
	row="| "
	caption_row="| "

	for ((i = 0; i < ${#images[@]}; i++)); do
		row+="![${captions[i]}](figures/${images[i]}){width=100px} | "
		caption_row+="**${captions[i]}** | "
		((count++))

		# End row if 3 images are added
		if ((count % 3 == 0)); then
			echo "$row" >>"$output"
			echo "|-----------------|-----------------|-----------------|" >>"$output"
			echo "$caption_row" >>"$output"
			# echo "" >>"$output"
			row="| "
			caption_row="| "
		fi
	done
}

# figure1's all bandage assembly_graph
# figure2.md
bandage2_genus_species() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local simple_polishing="${3:-default}"
	local target_index="${output_dir}-${isuffix}"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local ssh_remote="${_ssh["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"
	local extracted_n="${_compare_n["$target_index"]}"
	local i
	local width=13
	local images=()
	local captions=()

	local _figure2_md="figure2.md"
	local k=0
	for ((i = 0; i < extracted_n; i++)); do

		local _gfa_infer="${output_dir}/${extracted_inum}/disassemble/infer-1/1/${i}/30-contigger/graph_final.gfa"
		local _png_infer="${output_dir}/${extracted_inum}/disassemble/infer-1/1/${i}/30-contigger/graph_final.png"
		if [[ -s "${_gfa_infer}" ]]; then
			echo "gfa file: ${_gfa_infer}"
			${_polap_cmd} bandage png \
				${_gfa_infer} \
				${_png_infer}
			printf "| ![polap %s](figures/%s){ width=%s%% } " "${i}" "${_png_infer}" "${width}" >>"${_figure2_md}"
			images+=("figures/${_png_infer}")
			captions+=(${i})
			((k++))
		else
			echo "no such file: ${_gfa_infer}"
		fi
	done

	printf "Figure. Polap inference - %s\n" "${species_name}" >>"${_figure2_md}"

	#!/bin/bash

	# Define an array of PNG files and their subcaptions
	# images=("fig1.png" "fig2.png" "fig3.png" "fig4.png" "fig5.png" "fig6.png" "fig7.png" "fig8.png")
	# captions=("Subcaption 1" "Subcaption 2" "Subcaption 3" "Subcaption 4" "Subcaption 5" "Subcaption 6" "Subcaption 7" "Subcaption 8")

	# Output Markdown file
	output="${_figure2_md}"

	# Start writing to the markdown file
	cat <<EOF >"$output"
---
title: "${species_name}"
geometry: margin=1in
---

# Image Grid

| Organelle genome assemblies |  |  |
|-----------------|-----------------|-----------------|
EOF

	# Generate the image table with subcaptions
	count=0
	row="| "
	caption_row="| "

	for ((i = 0; i < ${#images[@]}; i++)); do
		row+="![${captions[i]}](${images[i]}){width=100px} | "
		caption_row+="**${captions[i]}** | "
		((count++))

		# End row if 3 images are added
		if ((count % 3 == 0)); then
			echo "$row" >>"$output"
			echo "|-----------------|-----------------|-----------------|" >>"$output"
			echo "$caption_row" >>"$output"
			# echo "" >>"$output"
			row="| "
			caption_row="| "
		fi
	done

	# Handle the last incomplete row (if any)
	remaining=$((3 - count % 3))
	if ((remaining < 3)); then
		# Fill empty image columns
		for ((i = 0; i < remaining; i++)); do
			row+=" | "
			caption_row+=" | "
		done
		echo "$row" >>"$output"
		echo "|-----------------|-----------------|-----------------|" >>"$output"
		echo "$caption_row" >>"$output"
	fi

	echo "" >>"${output}"
	printf "Figure. Polap inference - %s\n" "${species_name}" >>"${output}"

}

copy-figures_genus_species() {
	local _arg_t_dir="${1:-"${_arg_default_target_dir}"/figures}"
	rsync -av --include='*/' --include='*.png' --exclude='*' ./ "${_arg_t_dir}/"
	rsync -av --include='*/' --include='*.pdf' --exclude='*' ./ "${_arg_t_dir}/"
	# rsync -av --include='*/' --include='*.png' --exclude='*' ./ ../../manuscript/polap-v0.4/figures/
	# rsync -av --include='*/' --include='*.pdf' --exclude='*' ./ ../../manuscript/polap-v0.4/figures/
}

# downsample <folder> <isuffix> [--dry]
# e.g., downsample Pisum_sativum 0 --dry
downsample_genus_species() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local _dry="${3:-default_value}"
	local target_index="${output_dir}-${isuffix}"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local ssh_remote="${_ssh["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"
	local extracted_downsample="${_downsample["$target_index"]}"

	_c=${extracted_downsample}
	if [[ "${_dry}" == "default_value" ]]; then
		_dry=""
	fi

	local _d_i
	_d_i="downsample"

	mkdir -p "${output_dir}"

	local long_data="${_media_dir}/${long_sra}.fastq.tar.gz"
	local short_data="${_media_dir}/${short_sra}.fastq.tar.gz"

	if [[ "${_local_host}" == "$(hostname)" ]]; then
		if [[ -s "${long_data}" ]]; then
			tar -zxf ${long_data}
			tar -zxf ${short_data}
		else
			ln -s "${_media1_dir}/${long_sra}.fastq"
			ln -s "${_media1_dir}/${short_sra}_1.fastq"
			ln -s "${_media1_dir}/${short_sra}_2.fastq"
		fi
	# else
	# 	echo "You can send the data from the central host: ${_local_host}"
	# 	return 1
	fi

	if [[ -s "$output_dir/short_expected_genome_size.txt" ]]; then
		echo "Found: $output_dir/short_expected_genome_size.txt"
	else
		${_polap_cmd} find-genome-size \
			-a ${short_sra}_1.fastq \
			-b ${short_sra}_2.fastq \
			-o "${output_dir}"
	fi
	local _genome_size=$(<"$output_dir/short_expected_genome_size.txt")
	echo "Genome size estimate: ${_genome_size}"

	${_polap_cmd} fastq subsample -v ${_dry} \
		"${long_sra}.fastq" \
		"${output_dir}/l${_c}x.fq" \
		-c "${_c}" \
		-o "${output_dir}" \
		--random-seed "${random_seed}" \
		--genomesize "${_genome_size}" -v \
		>"${output_dir}/l${_c}x.txt"

	${_polap_cmd} fastq subsample2 -v ${_dry} \
		"${short_sra}_1.fastq" \
		"${short_sra}_2.fastq" \
		"${output_dir}/s${_c}x_1.fq" \
		"${output_dir}/s${_c}x_2.fq" \
		-c "${_c}" \
		-o "${output_dir}" \
		--random-seed "${random_seed}" \
		--genomesize "${_genome_size}" -v \
		>"${output_dir}/s${_c}x.txt"

}

# use-downsample <folder> [inum]
# replace the input data with the downsampled data
use-downsample_genus_species() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local target_index="${output_dir}-${isuffix}"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local ssh_remote="${_ssh["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"
	local extracted_downsample="${_downsample["$target_index"]}"

	_c=${extracted_downsample}

	# mv "${long_sra}.fastq" "${output_dir}/l.fq"
	rm -f "${long_sra}.fastq"
	# gunzip -f "${output_dir}/l${_c}x.fq.gz"
	mv "${output_dir}/l${_c}x.fq" "${long_sra}.fastq"

	# mv "${short_sra}_1.fastq" "${output_dir}/s_1.fq"
	rm -f "${short_sra}_1.fastq"
	# gunzip -f "${output_dir}/s${_c}x_1.fq.gz"
	mv "${output_dir}/s${_c}x_1.fq" "${short_sra}_1.fastq"

	# mv "${short_sra}_2.fastq" "${output_dir}/s_2.fq"
	rm -f "${short_sra}_2.fastq"
	# gunzip -f "${output_dir}/s${_c}x_2.fq.gz"
	mv "${output_dir}/s${_c}x_2.fq" "${short_sra}_2.fastq"
}

mauve_genus_species() {
	local output_dir="${1:-default}"
	local isuffix="${2:-0}"
	local target_index="${output_dir}-${isuffix}"

	if [[ "${output_dir}" == "default" ]]; then

		# Extract and sort keys
		sorted_keys=($(for key in "${!_table1[@]}"; do echo "$key"; done | sort))

		# Iterate over sorted keys and check if value is "T"
		for key in "${sorted_keys[@]}"; do
			if [[ "${_table1[$key]}" == "F" ]]; then
				continue
			fi
			local _v1=${_folder[$key]}

			mauve_genus_species_for "${_v1}" "${isuffix}"
		done
	else
		mauve_genus_species_for "${output_dir}" "${isuffix}"
	fi
}

mauve_genus_species_for() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local target_index="${output_dir}-${isuffix}"

	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	local i=0
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	mkdir -p "${output_dir_i}"

	i=$((i + 1))
	echo "($i)"

	local mauve_dir="${output_dir_i}/mauve/${i}"
	local mafft_dir="${output_dir_i}/mafft/${i}"
	local blast_dir="${output_dir_i}/blast/${i}"
	mkdir -p "${mauve_dir}"
	mkdir -p "${mafft_dir}"
	mkdir -p "${blast_dir}"
	if [[ -s "${output_dir_i}/disassemble/infer-${i}/pt.subsample-polishing.reference.aligned.1.fa" ]]; then
		${_polap_cmd} mafft-mtdna -a "${output_dir}/ptdna-ptgaul.fa" \
			-b "${output_dir_i}/disassemble/infer-${i}/pt.subsample-polishing.reference.aligned.1.fa" \
			-o "${mafft_dir}" \
			>"${mafft_dir}/log.txt"
		echo "see ${mafft_dir}/pident.txt"
		cat "${mafft_dir}/pident.txt"

		# ${_polap_cmd} mauve-mtdna -a "${output_dir}/ptdna-ptgaul.fa" \
		# 	-b "${output_dir_i}/disassemble/infer-${i}/pt.subsample-polishing.reference.aligned.1.fa" \
		# 	-o "${mauve_dir}" \
		# 	>"${mauve_dir}/log.txt"
		# echo "see ${mauve_dir}/log.txt"
		# cat "${mauve_dir}/log.txt"
		#
		# ${_polap_cmd} compare2ptdna -a "${output_dir}/ptdna-ptgaul.fa" \
		# 	-b "${output_dir_i}/disassemble/infer-${i}/pt.subsample-polishing.reference.aligned.1.fa" \
		# 	-o "${blast_dir}"
		# echo "see ${blast_dir}/pident.txt"
		# cat "${blast_dir}/pident.txt"
	else
		echo "ERROR: no such file: ${output_dir_i}/disassemble/infer-${i}/pt.subsample-polishing.reference.aligned.1.fa"
	fi

}

restart_genus_species() {
	local output_dir="$1"

	local target_dir="${output_dir}-t"
	mkdir "${target_dir}"

	for i in \
		00-bioproject \
		0-bioproject \
		getorganelle \
		ptgaul \
		msbwt \
		timing \
		polap.log; do
		mv "${output_dir}/$i" "${target_dir}"
	done
	mv "${output_dir}"/*.fa "${target_dir}"

	mv "${output_dir}" "${output_dir}-backup"
	mv "${target_dir}" "${output_dir}"
}

write-config_genus_species() {
	local csv_file="$1"

	# Write the output to CSV
	{
		# Print the header
		echo "species,_long,_short,_host,_ptgaul_genomesize,_compare_n,_compare_p,_compare_r,_polish_n,_polish_p,_random_seed"

		# Loop through species and print their values
		local seed=101
		for key in "${!_long[@]}"; do
			if [[ -z "${_ptgaul_genomesize[$key]}" ]]; then
				_ptgaul_genomesize[$key]=160000
			fi
			_random_seed[$key]="${seed}"
			_polish_n[$key]="5"
			_polish_p[$key]="5"
			echo "$key,${_long[$key]},${_short[$key]},${_host[$key]},${_ptgaul_genomesize[$key]},${_compare_n[$key]},${_compare_p[$key]},${_compare_r[$key]},${_polish_n[$key]},${_polish_p[$key]},${_random_seed[$key]}"
			seed=$((seed + 2))
		done
	} >"$csv_file"
}

wga_genus_species() {
	local output_dir="$1"
	local species_name="$(echo $1 | sed 's/_/ /')"
	local long_sra="${_long["$1"]}"
	local short_sra="${_short["$1"]}"
	# copy_data

	rm -rf "${output_dir}/0"

	${_polap_cmd} assemble1 \
		-o ${output_dir} \
		-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
		--stopafter data

	rm -rf "${output_dir}/0"

	command time -v ${_polap_cmd} flye1 polishing \
		-o ${output_dir} \
		-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
		2>${output_dir}/timing-flye1.txt
}

archive_genus_species_for() {
	local output_dir="$1"
	# copy_data

	rm -rf "${output_dir}-a"
	rm -f "${output_dir}-a.tar.gz"
	${_polap_cmd} disassemble archive \
		--max-filesize 5M \
		-o ${output_dir}
}

archive_genus_species() {
	local output_dir="${1:-all}"

	if [[ "${output_dir}" == "all" ]]; then
		for _v1 in "${Sall[@]}"; do
			archive_genus_species_for "${_v1}"
		done
	else
		archive_genus_species_for "${output_dir}"
	fi
}

get_genus_species() {
	local output_dir="${1:-all}"
	local _brg_inum="${2:--1}"
	local _brg_confirm="${3:-off}"
	local _brg_download="${4:-on}"
	local target_index="${output_dir}-0"

	if [[ "${output_dir}" == "all" ]]; then
		for _v1 in "${Sall[@]}"; do
			get_genus_species_for "${_v1}" "${_brg_inum}" "${_brg_confirm}" "${_brg_download}"
		done
	else
		get_genus_species_for "${output_dir}" "${_brg_inum}" "${_brg_confirm}" "${_brg_download}"
	fi
}

get_genus_species_for() {
	local output_dir="${1:-default}"
	local _brg_inum="${2:--1}"
	local _brg_confirm="${3:-off}"
	local _brg_download="${4:-on}"
	local target_index="${output_dir}-0"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local ssh_remote="${_ssh["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	# if [[ "${_local_host}" == "$(hostname)" ]]; then

	# Define the directory name (can be passed as an argument)
	dir_to_check="${output_dir}"

	if [[ "${_brg_confirm}" == "off" ]]; then
		# Check if the directory exists
		if [[ -d "$dir_to_check" ]]; then
			echo "Directory '$dir_to_check' exists."
			# Ask for user confirmation before deleting
			read -p "Do you want to delete or add ${_brg_inum} to it? (a/N/y): " confirm
		else
			echo "Directory '$dir_to_check' does not exist."
			confirm="yes"
		fi
	else
		confirm="${_brg_confirm}"
	fi

	case "$confirm" in
	[yY] | [yY][eE][sS])

		local _atargz="${output_dir}-a.tar.gz"
		if [[ "${_brg_inum}" != "-1" ]]; then
			_atargz="${output_dir}-a-${_brg_inum}.tar.gz"
		fi

		if [[ "${_brg_download}" == "on" ]]; then
			rm -f "${_atargz}"
			scp -p "${ssh_remote}:$PWD/${_atargz}" .
		fi

		if [[ -s "${_atargz}" ]]; then
			tar zxf "${_atargz}"
			mv "$dir_to_check/0-bioproject" "${output_dir}-a/"
			mv "$dir_to_check/bioproject.txt" "${output_dir}-a/"
			rm -rf "$dir_to_check"
			echo "Directory '$dir_to_check' deleted."
			mv ${output_dir}-a ${output_dir}
			echo "${output_dir} is replaced by ${output_dir}-a."
		else
			echo "INFO: no such file: ${_atargz}, so we do nothing."
		fi
		;;
	[aA] | [aA][dD][dD])
		if [[ "${_brg_inum}" == "-1" ]]; then
			echo "ERRO: <inum> must be zero or positive for add operation."
		else
			echo "Replacing ${output_dir}/${_brg_inum} ..."
			local _atargz="${output_dir}-a-${_brg_inum}.tar.gz"
			if [[ "${_brg_download}" == "on" ]]; then
				echo "  fetching ${ssh_remote}:$PWD/${_atargz} ..."
				rm -f "${_atargz}"
				scp -p "${ssh_remote}:$PWD/${_atargz}" .
			fi

			if [[ -s "${_atargz}" ]]; then
				tar zxf "${_atargz}"
				rm -rf "${output_dir}/${_brg_inum}"
				mv ${output_dir}-a/${_brg_inum} ${output_dir}
				echo "  ${output_dir}/${_brg_inum} is replaced by ${output_dir}-a/${_brg_inum}"
				rm -rf "${output_dir}-a"
			else
				echo "INFO: no such file: ${_atargz}, so we skip adding for ${output_dir}/${_brg_inum}"
			fi
		fi

		;;
	*)
		echo "Deletion canceled."
		;;
	esac
}

# TODO
report_genus_species() {
	local _brg_outdir="${1:-all}"
	local _brg_inum="${2:-2}"
	local _brg_d_index="${3:-infer-1}"

	if [[ "${_brg_outdir}" == "all" ]]; then
		for _v1 in "${Sall[@]}"; do
			report_genus_species_for "${_v1}" "${_brg_inum}" "${_brg_d_index}"
		done
	else
		report_genus_species_for "$@"
	fi
}

# TODO
report_genus_species_for() {
	local _brg_outdir="${1}"
	local _brg_inum="${2:-2}"
	local _brg_d_index="${3:-infer-1}"
	local _key="${_brg_outdir}-${_brg_inum}"

	local i=0
	local n
	local p
	IFS=':' read -r -a extracted_array_n <<<"${_compare_n["$_key"]}"
	IFS=':' read -r -a extracted_array_p <<<"${_compare_p["$_key"]}"
	for n in "${extracted_array_n[@]}"; do
		for p in "${extracted_array_p[@]}"; do
			i=$((i + 1))
			k="3x"
			${_polap_cmd} disassemble report ${k} 3-infer \
				-o "${_brg_outdir}" \
				-i "${_brg_inum}" \
				--disassemble-i "${_brg_d_index}"

			# for k in {1..2}; do
			# 	${_polap_cmd} disassemble report ${k} infer \
			# 		-o ${output_dir} \
			# 		--disassemble-i infer-$i
			# 	${_polap_cmd} disassemble report ${k} \
			# 		-o ${output_dir} \
			# 		--disassemble-i compare-$i
			# done

		done
	done
}

################################################################################
# create tables and figures
# some helper functions first and then the main table and figure functions
#

# Example usage
# echo "Hours: $(convert_to_hours '4:05:17')"  # h:mm:ss
# echo "Hours: $(convert_to_hours '58:07.72')" # mm:ss.ss
# echo $(convert_to_hours "0:00.53")
# echo $(convert_to_hours "58:07.72")
convert_to_hours() {
	local time_string="$1"
	local hours=0

	if [[ "$time_string" =~ ^([0-9]+):([0-9]{2}):([0-9]{2})$ ]]; then
		# h:mm:ss format
		local h="${BASH_REMATCH[1]}"
		local m="${BASH_REMATCH[2]}"
		local s="${BASH_REMATCH[3]}"
		hours=$(bc <<<"scale=1; $h + $m / 60 + $s / 3600")
	elif [[ "$time_string" =~ ^([0-9]+):([0-9]{2})\.([0-9]{2})$ ]]; then
		# mm:ss.ss format
		local m="${BASH_REMATCH[1]}"
		local s="${BASH_REMATCH[2]}.${BASH_REMATCH[3]}"

		hours=$(bc <<<"scale=1; $m / 60 + $s / 3600")
	else
		echo "Invalid time format: $time_string" >&2
		return 1
	fi

	echo "$hours"
}

# Example usage
# result=$(parse_params "params.txt")
# read -r I P N R <<< "$result"
parse_params() {
	local file="$1" # Input file
	local I P N R   # Declare variables for the parameters

	# Read the file and extract the values
	while IFS=": " read -r key value; do
		case "$key" in
		"I") I="$value" ;;
		"P") P="$value" ;;
		"N") N="$value" ;;
		"R") R="$value" ;;
		esac
	done <"$file"

	# Print the variables (return as output)
	echo "$I $P $N $R"
}

# read -r memory1 time1 < <(parse_timing Eucalyptus_pauciflora 3)
parse_timing() {
	local _v1="${1}"
	local j="${2}"
	local _memory_gb
	local _total_hours
	local _timing_file=${_v1}/timing-${j}.txt
	local _params_txt=${_v1}/disassemble/${j}/params.txt
	# percent_identity=$(<"${_v1}/disassemble/${j}/c/coverage.txt")
	# ipn=$(parse_params "${params_txt}")
	# read -r I P N <<<"$ipn"

	if [[ -s "${_timing_file}" ]]; then
		local _time_wga=$(grep 'Elapsed' "${_timing_file}" | head -1)
		local _memory_wga=$(grep 'Maximum resident set size' "${_timing_file}" | head -1)
		# Extract the number in kilobytes
		local _memory_kbytes=$(echo "$_memory_wga" | grep -oE "[0-9]+")
		# Convert kilobytes to gigabytes
		# Extract the time portion using grep with regex
		local _memory_gb=$(echo "scale=2; $_memory_kbytes / 1048576" | bc)
		# time_only=$(echo "$_time_wga" | grep -oE "[0-9]+(:[0-9]{2}){1,2}")
		local time_only=$(grep "Elapsed (wall clock) time (h:mm:ss or m:ss):" "$_timing_file" | awk -F': ' '{print $2}')

		local _total_hours=$(convert_to_hours "${time_only}")
	else
		_memory_gb=0
		_total_hours=0
	fi

	echo "${_memory_gb} ${_total_hours}"
}

function _polap_utility_convert_bp {
	local bp=${1%.*}
	if ((bp >= 1000000000000)); then
		echo "$(bc <<<"scale=1; $bp/1000000000000") Tbp"
	elif ((bp >= 1000000000)); then
		echo "$(bc <<<"scale=1; $bp/1000000000") Gbp"
	elif ((bp >= 1000000)); then
		echo "$(bc <<<"scale=1; $bp/1000000") Mbp"
	elif ((bp >= 1000)); then
		echo "$(bc <<<"scale=1; $bp/1000") kbp"
	else
		echo "$bp bp"
	fi
}

# extract short-read coverage
get_short_read_coverage() {
	local file="$1"
	if [[ ! -f "$file" ]]; then
		echo "Error: File not found!"
		return 1
	fi

	local coverage
	# coverage=$(grep -oP 'short-read coverage:\s*\K\d+' "$file")
	coverage=$(grep -oP 'short-read coverage:\s*\K[\d.]+(?=x)' "$file" | tail -1)

	if [[ -n "$coverage" ]]; then
		coverage=$(echo "$coverage" | tr -d '[:space:]') # Remove spaces
		printf "%.1f\n" "$coverage"
	else
		echo "Error: Coverage information not found!"
		return 1
	fi
}

get_long_read_coverage() {
	local file="$1"
	if [[ ! -f "$file" ]]; then
		echo "Error: File not found!"
		return 1
	fi

	local coverage
	# coverage=$(grep -oP 'long-read coverage:\s*\K\d+' "$file")
	coverage=$(grep -oP 'long-read coverage:\s*\K[\d.]+(?=x)' "$file" | tail -1)

	if [[ -n "$coverage" ]]; then
		coverage=$(echo "$coverage" | tr -d '[:space:]') # Remove spaces
		printf "%.1f\n" "$coverage"
	else
		echo "Error: Coverage information not found!"
		return 1
	fi
}

get_target_read_coverage() {
	local file="$1"
	if [[ ! -f "$file" ]]; then
		echo "Error: File not found!"
		return 1
	fi

	local coverage
	# coverage=$(grep -oP 'target coverage:\s*\K\d+' "$file")
	coverage=$(grep -oP 'target coverage:\s*\K[\d.]+(?=x)' "$file" | tail -1)

	if [[ -n "$coverage" ]]; then
		coverage=$(echo "$coverage" | tr -d '[:space:]') # Remove spaces
		printf "%.1f\n" "$coverage"
	else
		echo "Error: Coverage information not found!"
		return 1
	fi
}

maintable1_genus_species_header() {
	local _arg_table="${1:-1}"

	local _items=(
		"Species"
		"Order"
		"Family"
		"I"
		"C"
		"L_SRA"
		"L_size"
		"L_cov"
		"S_SRA"
		"S_size"
		"S_cov"
		"G"
		"N"
		"P"
		"R"
		"Rate"
		"Size"
		"Alpha"
		"ptDNA"
		"Length1"
		"Length2"
		"Pident"
		"N1"
		"Mode"
		"SD"
		"N2"
		"M"
		"M_g"
		"M_p"
		"M_t1"
		"M_t2"
		"M_s"
		"M_f"
		"T"
	)

	if ((_arg_table == 2)); then
		_items=(
			"Species"
			"Order"
			"Family"
			"N"
			"C"
			"L_SRA"
			"L_size"
			"L_cov"
			"S_SRA"
			"S_size"
			"S_cov"
			"P"
			"R"
			"Rate"
			"Size"
			"Alpha"
			"Length2"
			"N1"
			"Mode"
			"SD"
			"N2"
			"M"
			"M_g"
			"M_f"
			"T"
		)
	fi

	printf "%s\t" "${_items[@]::${#_items[@]}-1}"
	printf "%s\n" "${_items[-1]}"
}

# Table 1's row for a given outdir and inum.
#
# arg1: outdir
# arg2: inum
# arg3: disassemble-i
# arg4: 1 or 2
# arg5: target directory to copy the result
#
# outdir-inum is used as the key as well.
#
# return:
#
maintable1_genus_species_for() {
	local _arg_outdir="${1:-all}"
	local _arg_inum="${2:-2}"
	local _arg_d_index="${3:-infer-1}"
	local _arg_table="${4:-1}"
	local _arg_t_dir="${5:-"${_arg_default_target_dir}"}"
	local _key="${_arg_outdir}-${_arg_inum}"

	local _v1="${_arg_outdir}"
	local _polap_log=${_v1}/polap.log

	# Taxon names: species, genus, family, and order
	local _species=${_taxon[$_key]}
	local _species="${_species//_/ }"
	local _genus=${_species%% *}
	local _family=$(grep "${_genus}" ${script_dir}/polaplib/taxonomy_output.tsv | cut -f 7 | head -1)
	local _order=$(grep "${_genus}" ${script_dir}/polaplib/taxonomy_output.tsv | cut -f 6 | head -1)

	# Long SRA
	local _l_sra=$(awk '/long-read/ {match($0, /([A-Z]RR[0-9]+)/, arr); print arr[0]}' "$_polap_log" | sort -u | grep -v '^$')
	if [[ -s "${_v1}/long_total_length.txt" ]]; then
		_l_sra_size=$(<"${_v1}/long_total_length.txt")
	else
		_l_sra_size=$(<"${_v1}/l.fq.txt")
	fi
	local _l_sra_size_gb=$(_polap_utility_convert_bp "${_l_sra_size}")

	# Short SRA
	local _s_sra=$(awk '/short-read1/ {match($0, /([A-Z]RR[0-9]+)/, arr); print arr[0]}' "$_polap_log" | sort -u | grep -v '^$')
	if [[ -s "${_v1}/short_total_length.txt" ]]; then
		local _s_sra_size=$(<"${_v1}/short_total_length.txt")
	else
		local _s_sra_size1=$(<"${_v1}/s1.fq.txt")
		local _s_sra_size2=$(<"${_v1}/s2.fq.txt")
		local _s_sra_size=$((_s_sra_size1 + _s_sra_size2))
	fi
	local _s_sra_size_gb=$(_polap_utility_convert_bp "${_s_sra_size}")

	# Known ptDNA NCBI accession
	local _known_mtdna=$(grep 'NCBI accession:' ${_polap_log} | cut -d: -f4 | tail -n 1 || true)
	_known_mtdna=${_known_mtdna:-'NA'}

	# Genome size estimates
	if [[ -s "${_v1}/short_expected_genome_size.txt" ]]; then
		_genome_size=$(<"${_v1}/short_expected_genome_size.txt")
	else
		echo "ERROR: no such file: ${_v1}/short_expected_genome_size.txt"
		exit 1
	fi

	if ((_arg_table == 1)); then
		local _ptdna_ptgaul=${_v1}/ptdna-ptgaul.fa
		local _ptdna_reference=${_v1}/ptdna-reference.fa
		if [[ -s "${_ptdna_ptgaul}" ]]; then
			local _seq_length_ptgaul=$(bioawk -c fastx 'NR==1 {print length($seq); exit}' "${_ptdna_ptgaul}")
		else
			_seq_length_ptgaul='NA'
		fi
	fi

	# read -r _memory_gb_flye1 _total_hours_flye1 < <(parse_timing "${_v1}" "infer12-4")

	local j=1
	local _i="${_arg_inum}"
	local _v1_inum="${_v1}/${_arg_inum}"
	local _disassemble_index="${_arg_d_index}"
	local _extracted_memory="${_memory["$_key"]}"
	local _summary1_ordered_txt="${_v1_inum}/disassemble/${_disassemble_index}/1/summary1-ordered.txt"

	local _ptdna_subsample="${_v1_inum}/disassemble/${_disassemble_index}/pt.subsample-polishing.1.fa"
	local _seq_length_subsample=0
	if [[ -s "${_ptdna_subsample}" ]]; then
		# Count the number of sequences (lines starting with '>')
		local _seq_count=$(grep -c "^>" "$_ptdna_subsample")

		# Ensure there is exactly one sequence
		if ((_seq_count != 1)); then
			echo "Error: FASTA file does not contain exactly one sequence: ${_ptdna_subsample}"
			exit 1
		fi

		# Compute the sequence length
		_seq_length_subsample=$(bioawk -c fastx 'NR==1 {print length($seq); exit}' "${_ptdna_subsample}")
	fi

	# Extract mode value
	# Extract SD value
	# Extract the first index value
	local _mode=$(grep "^#mode:" "$_summary1_ordered_txt" | awk '{print $2}')
	local _sd=$(grep "^#sd:" "$_summary1_ordered_txt" | awk '{print $2}')
	local _first_index=$(grep "^#index:" "$_summary1_ordered_txt" | awk 'NR==1 {print $2}')
	local _n1=$(grep "^#n:" "$_summary1_ordered_txt" | awk 'NR==1 {print $2}')
	local _output=$(awk -F'\t' 'NR==2 {print $1}' "${_summary1_ordered_txt}")
	read -r _second_line_index <<<"$_output"
	if ((_first_index != _second_line_index)); then
		echo "ERROR: #index: ${_first_index} vs. #2nd line index: ${_second_line_index}"
		exit 1
	fi

	local _params_txt=${_v1_inum}/disassemble/${_disassemble_index}/params.txt
	local _ipn=$(parse_params "${_params_txt}")
	local _I _P _N _R
	read -r _I _P _N _R <<<"$_ipn"

	# The determined subsampling rate and alpha
	local _summary2_ordered_txt=${_v1_inum}/disassemble/${_disassemble_index}/2/summary1-ordered.txt
	local _n2=$(grep "^#n:" "$_summary2_ordered_txt" | awk 'NR==1 {print $2}')
	local _output=$(awk -F'\t' 'NR==2 {print $1, $2, $4, $11}' "${_summary2_ordered_txt}")
	local _summary2_index
	local _summary2_size
	local _summary2_size_gb
	local _summary2_rate_rounded
	local _summary2_alpha
	local _summary2_alpha_formatted
	read -r _summary2_index _summary2_size _summary2_rate _summary2_alpha <<<"$_output"
	_summary2_size_gb=$(_polap_utility_convert_bp "${_summary2_size}")
	local _summary2_rate_decimal=$(printf "%.10f" "$_summary2_rate")
	_summary2_rate_rounded=$(echo "scale=4; $_summary2_rate_decimal / 1" | bc)
	_summary2_alpha_formatted=$(echo "scale=2; $_summary2_alpha / 1" | bc | awk '{printf "%.2f\n", $1}')

	# Timing
	if ((_arg_table == 1)); then
		read -r _memory_gb_getorganelle _total_hours_getorganelle < <(_polap_lib_timing-parse-timing "${_v1}/timing/timing-getorganelle.txt")
		if [[ -s "${_v1}/timing/timing-ptgaul.txt" ]]; then
			read -r _memory_gb_ptgaul _total_hours_ptgaul < <(_polap_lib_timing-parse-timing "${_v1}/timing/timing-ptgaul.txt")
			read -r _memory_gb_prepare_polishing _total_hours_prepare_polishing < <(_polap_lib_timing-parse-timing "${_v1}/timing/timing-prepare-polishing.txt")
			read -r _memory_gb_ptgaul_polishing _total_hours_ptgaul_polishing < <(_polap_lib_timing-parse-timing "${_v1}/timing/timing-ptgaul-polishing.txt")
			read -r _memory_gb_subsampling_polishing _total_hours_subsampling_polishing < <(_polap_lib_timing-parse-timing "${_v1_inum}/timing-check-${j}-subsample-polish.txt")
		else
			_memory_gb_ptgaul='NA'
			_memory_gb_prepare_polishing='NA'
			_memory_gb_ptgaul_polishing='NA'
			_memory_gb_subsampling_polishing='NA'
		fi
		read -r _memory_gb _total_hours < <(_polap_lib_timing-parse-timing "${_v1_inum}/timing-infer-${j}-subsample-polish.txt")
	elif ((_arg_table == 2)); then
		read -r _memory_gb_getorganelle _total_hours_getorganelle < <(_polap_lib_timing-parse-timing "${_v1}/timing/timing-getorganelle.txt")
		read -r _memory_gb_prepare_polishing _total_hours_prepare_polishing < <(_polap_lib_timing-parse-timing "${_v1}/timing/timing-prepare-polishing.txt")
		read -r _memory_gb_subsampling_polishing _total_hours_subsampling_polishing < <(_polap_lib_timing-parse-timing "${_v1_inum}/timing-check-${j}-subsample-polish.txt")
		read -r _memory_gb _total_hours < <(_polap_lib_timing-parse-timing "${_v1_inum}/timing-infer-${j}-subsample-polish.txt")

		#   read -r _memory_gb_getorganelle _total_hours_getorganelle < <(_polap_lib_timing-parse-timing "${_v1}/timing/timing-getorganelle.txt")
		# read -r _memory_gb_polishing _total_hours_polishing < <(_polap_lib_timing-parse-timing "${_v1}/timing/timing-prepare-polishing.txt")
		#   read -r _memory_gb_subsampling_polishing _total_hours_subsampling_polishing < <(_polap_lib_timing-parse-timing "${_v1_inum}/timing-infer-${j}-subsample-polish.txt")
		#   read -r _memory_gb _total_hours < <(_polap_lib_timing-parse-timing "${_v1_inum}/timing-infer-${j}-subsample-polish.txt")
	fi

	# Percent identity from pairwise sequence alignment
	local _mafft_pident="NA"
	if [[ -s "${_v1_inum}/mafft/${j}/pident.txt" ]]; then
		_mafft_pident="$(<${_v1_inum}/mafft/${j}/pident.txt)"
	fi

	# local _blast_pident="NA"
	# if [[ -s "${_v1_inum}/blast/${j}/pident.txt" ]]; then
	# 	_blast_pident="$(<${_v1_inum}/blast/${j}/pident.txt)"
	# else
	# 	_blast_pident=0
	# fi

	# local _mauve_lcb_coverage="NA"
	# if [[ -s "${_v1_inum}/mauve/${j}/log.txt" ]]; then
	# 	_mauve_lcb_coverage=$(awk '{print $2}' "${_v1_inum}/mauve/${j}/log.txt")
	# else
	# 	_mauve_lcb_coverage=0
	# fi

	local _pident="${_mafft_pident}"

	# sequencing data coverage
	local _short_coverage=$(get_short_read_coverage "${_v1_inum}/sx.txt")
	local _long_coverage=$(get_long_read_coverage "${_v1_inum}/lx.txt")
	local _target_coverage=$(get_target_read_coverage "${_v1_inum}/lx.txt")

	local _items
	if ((_arg_table == 1)); then
		local _items=(
			"_${_species}_"
			"${_order}"
			"${_family}"
			"${_i}"
			"${_target_coverage}"
			"${_l_sra}"
			"${_l_sra_size_gb}"
			"${_long_coverage}"
			"${_s_sra}"
			"${_s_sra_size_gb}"
			"${_short_coverage}"
			"${_genome_size}"
			"${_N}"
			"${_P}"
			"${_R}"
			"${_summary2_rate_rounded}"
			"${_summary2_size_gb}"
			"${_summary2_alpha_formatted}"
			"${_known_mtdna}"
			"${_seq_length_ptgaul}"
			"${_seq_length_subsample}"
			"${_pident}"
			"${_n1}"
			"${_mode}"
			"${_sd}"
			"${_n2}"
			"${_extracted_memory}"
			"${_memory_gb_getorganelle}"
			"${_memory_gb_ptgaul}"
			"${_memory_gb_prepare_polishing}"
			"${_memory_gb_ptgaul_polishing}"
			"${_memory_gb_subsampling_polishing}"
			"${_memory_gb}"
			"${_total_hours}"
		)
	elif ((_arg_table == 2)); then
		_items=(
			"_${_species}_"
			"${_order}"
			"${_family}"
			"${_i}"
			"${_target_coverage}"
			"${_l_sra}"
			"${_l_sra_size_gb}"
			"${_long_coverage}"
			"${_s_sra}"
			"${_s_sra_size_gb}"
			"${_short_coverage}"
			"${_P}"
			"${_R}"
			"${_summary2_rate_rounded}"
			"${_summary2_size_gb}"
			"${_summary2_alpha_formatted}"
			"${_seq_length_subsample}"
			"${_n1}"
			"${_mode}"
			"${_sd}"
			"${_n2}"
			"${_extracted_memory}"
			"${_memory_gb_getorganelle}"
			"${_memory_gb}"
			"${_total_hours}"
		)
	fi

	printf "%s\t" "${_items[@]::${#_items[@]}-1}" >>"${_table_tsv}"
	printf "%s\n" "${_items[-1]}" >>"${_table_tsv}"
}

# Table 1's row for a given outdir and inum.
#
# arg1: outdir
# arg2: inum
# arg3: disassemble-i
# arg4: 1 or 2
# arg5: target directory to copy the result
#
maintable1_genus_species() {
	local _arg_outdir="${1:-all}"
	local _arg_inum="${2:-2}"
	local _arg_d_index="${3:-infer-1}"
	local _arg_table="${4:-1}"
	local _arg_t_dir="${5:-"${_arg_default_target_dir}"}"
	# local _table_name=$(echo "${FUNCNAME[0]}" | cut -d'_' -f1)

	local _table_name="maintable1"
	if ((_arg_table == 2)); then
		_table_name="maintable2"
	fi
	local _table_tsv="${_table_name}-${_arg_inum}.tsv"

	maintable1_genus_species_header "${_arg_table}" >"${_table_tsv}"

	if [[ "${_arg_outdir}" == "all" ]]; then

		# Extract and sort keys
		sorted_keys=($(for key in "${!_table1[@]}"; do echo "$key"; done | sort))

		# Iterate over sorted keys and check if value is "T"
		for key in "${sorted_keys[@]}"; do
			if ((_arg_table == 1)); then
				if [[ "${_table1[$key]}" == "F" ]]; then
					continue
				fi
			elif ((_arg_table == 2)); then
				if [[ "${_table2[$key]}" == "F" ]]; then
					continue
				fi
			fi
			local _extracted_inum=${_inum[$key]}
			if [[ "${_extracted_inum}" != "${_arg_inum}" ]]; then
				continue
			fi
			local _v_folder=${_folder[$key]}
			local _v_folder_inum="${_v_folder}/${_arg_inum}"
			if [[ ! -d "${_v_folder_inum}" ]]; then
				continue
			fi

			maintable1_genus_species_for "${_v_folder}" \
				"${_arg_inum}" "${_arg_d_index}" "${_arg_table}" "${_arg_t_dir}" \
				>>"${_table_tsv}"
		done
	else
		maintable1_genus_species_for "$@" >>"${_table_tsv}"
	fi

	if ((_arg_table == 1)); then
		csvtk -t cut -f Species,C,N,P,R,Rate,Alpha,Length1,Length2,Pident,N1,Mode,SD,M,M_g,M_p,M_t1,M_t2,M_s,M_f,T \
			${_table_tsv} |
			csvtk -t rename -f 1-21 -n Species,C,N,P,R,Rate,Alpha,L1,L2,Pident,N1,Mode,SD,M,Mg,Mp,Mt1,Mt2,Ms,Mf,T |
			csvtk -t csv2md -a right -o ${_table_name}-${_arg_inum}-analysis.md

	elif ((_arg_table == 2)); then
		csvtk -t cut -f Species,C,N,P,R,Rate,Alpha,Length2,N1,Mode,SD,M,M_g,M_f,T \
			${_table_tsv} |
			csvtk -t rename -f 1-15 -n Species,C,N,P,R,Rate,Alpha,L2,N1,Mode,SD,M,Mg,Mf,T |
			csvtk -t csv2md -a right -o ${_table_name}-${_arg_inum}-analysis.md
	fi
	csvtk -t cut -f Species,Order,Family,L_SRA,L_size,L_cov,S_SRA,S_size,S_cov ${_table_tsv} |
		csvtk -t rename -f 1-9 -n Species,Order,Family,L_SRA,L_size,L_cov,S_SRA,S_size,S_cov |
		csvtk -t csv2md -a right -o ${_table_name}-${_arg_inum}-data.md

	echo "  copying to ${_arg_t_dir}"
	cp -p ${_table_name}-${_arg_inum}-analysis.md "${_arg_t_dir}"
	cp -p ${_table_name}-${_arg_inum}-data.md "${_arg_t_dir}"

	csvtk -t csv2md -a right ${_table_tsv} -o ${_table_name}-${_arg_inum}.md

	# Rscript src/polap-data-v2.R

	cat "${_table_name}-${_arg_inum}.md"
	echo "ouput: ${_table_name}-${_arg_inum}-analysis.md"
	echo "ouput: ${_table_name}-${_arg_inum}-data.md"
	echo "ouput: ${_table_name}-${_arg_inum}.md"
	echo "ouput: ${_table_tsv}"
}

supptable1_genus_species_for() {
	local _arg_outdir="${1}"
	local _arg_inum="${2:-2}"
	local _arg_d_index="${3:-infer-1}"
	local _arg_stage="${4:-1}"
	local _arg_table="${5:-1}"
	local _arg_t_dir="${6:-"${_arg_default_target_dir}"}"
	local _arg_one="${7:-0}"
	local _key="${_arg_outdir}-${_arg_inum}"

	local _table_name="supptable1"
	if ((_arg_table == 2)); then
		_table_name="supptable2"
	fi
	local _supptable_md="${_table_name}-${_arg_inum}-${_arg_stage}.md"
	if [[ "${_arg_one}" == "1" ]]; then
		_supptable_md="${_table_name}-${_arg_inum}-${_arg_stage}-${_arg_outdir}.md"
	fi

	local _i=${_arg_inum}
	local j="${_arg_stage}"

	local _species=${_taxon[$_key]}
	local _species="${_species//_/ }"
	local _v1=${_folder[$_key]}
	local _label_base="${_v1/_/-}"
	_label_base=$(echo "$_label_base" | awk '{print tolower($0)}')

	local _v1_inum="${_v1}/${_arg_inum}"
	local _params_txt=${_v1_inum}/disassemble/${_arg_d_index}/params.txt
	_ipn=$(parse_params "${_params_txt}")
	local _I _P _N _R
	read -r _I _P _N _R <<<"$_ipn"
	local _label="${_arg_inum}-${_label_base}"
	if [[ "${_arg_one}" == "1" ]]; then
		_label="main-${_arg_inum}-${_label_base}"
	fi

	case "${_arg_stage}" in
	1*)
		printf "Table: Plastid genome assemblies with the increasing subsample size upto the maximum subsampling rate of ${_P}%% from Stage 1 of Polap's subsampling-based analysis for the dataset of _${_species}_. {#tbl:supptable1-${_label}}\n\n" \
			>>"${_supptable_md}"
		;;
	2*)
		printf "Table: Plastid genome assemblies with a fixed subsample size and subsampling-rate from Stage 2 of Polap's subsampling-based analysis for the dataset of _${_species}_. {#tbl:supptable2-${_label}}\n\n" \
			>>"${_supptable_md}"
		;;
	3*)
		printf "Table: Plastid genome assemblies with the subsampling-based short-read polishing from Stage 3 of Polap's subsampling-based analysis for the dataset of _${_species}_. {#tbl:supptable3-${_label}}\n\n" \
			>>"${_supptable_md}"
		;;
	x)
		printf "Table: Three stages for subsampling-based plastid genome assemblies with the increasing subsample size upto the maximum subsampling rate of ${_P}%% from Stage 1 of the subsampling-based analysis for the dataset of _${_species}_. {#tbl:supptable1-${_label}}\n\n" \
			>>"${_supptable_md}"
		;;
	*)
		echo "ERROR: no such stage: ${_arg_stage}"
		exit 1
		;;
	esac

	case "${_arg_stage}" in
	1* | 2* | 3*)
		cat "${_v1_inum}/disassemble/${_arg_d_index}/${_arg_stage}/summary1.md" \
			>>"${_supptable_md}"
		;;
	2*)
		cat "${_v1_inum}/disassemble/${_arg_d_index}/${_arg_stage}/summary1.md" \
			>>"${_supptable_md}"
		;;
	3*)
		cat "${_v1_inum}/disassemble/${_arg_d_index}/${_arg_stage}/summary1.md" \
			>>"${_supptable_md}"
		;;
	x)
		#!/bin/bash

		# Output file
		output="${_supptable_md}"

		# Print the header for Stage 1
		echo "| Stage 1 |   I |   Rate | Alpha | Pmem |       G | Time |   N |      L |   C | Length |" >>$output
		echo "| -----: | --: | -----: | ----: | ---: | ------: | ---: | --: | -----: | --: | -----: |" >>$output

		# Append the content of 1.md with Stage 1 formatting
		awk 'NR>2 {print "|        " $0}' "${_v1_inum}/disassemble/${_arg_d_index}/1/summary1.md" >>$output

		# Print the separation line
		echo "|        |      |               |          |        |                |          |       |              |       |               |" >>$output

		# Print the header for Stage 2
		echo "| Stage 2 |   I |   Rate | Alpha | Pmem |       G | Time |   N |       L |   C | Length |" >>$output
		echo "|        |      |               |          |        |                |          |       |              |       |               |" >>$output

		# Append the content of 2.md with Stage 2 formatting
		awk 'NR>2 {print "|       " $0}' "${_v1_inum}/disassemble/${_arg_d_index}/2/summary1.md" >>$output

		# Print the separation line
		echo "|        |      |               |          |        |                |          |       |              |       |               |" >>$output

		# Print the header for Stage 3
		echo "| Stage 3 |   I |     Rate |   Size  | Seed   |  Tp |  Mp |  Ts |  Ms | Pident |    Length |" >>$output
		echo "|        |      |               |          |        |                |          |       |              |       |               |" >>$output

		# Append the content of 3.md with Stage 3 formatting
		awk 'NR>2 {print "|      " $0}' "${_v1_inum}/disassemble/${_arg_d_index}/3-infer/summary1x.md" >>$output

		# cat "${_v1_inum}/disassemble/${_arg_d_index}/1/summary1.md" \
		# 	>>"${_supptable_md}"
		# cat "${_v1_inum}/disassemble/${_arg_d_index}/2/summary1.md" \
		# 	>>"${_supptable_md}"
		# cat "${_v1_inum}/disassemble/${_arg_d_index}/3-infer/summary1.md" \
		# 	>>"${_supptable_md}"
		;;
	*)
		echo "ERROR: no such stage: ${_arg_stage}"
		exit 1
		;;
	esac

	# cat "${_v1_inum}/disassemble/${_arg_d_index}/${_arg_stage}/summary1.md" \
	# 	>>"${_supptable_md}"

	printf "\n" \
		>>"${_supptable_md}"

	case "${_arg_stage}" in
	1*)
		cat "${script_dir}"/polaplib/polap-data-v2-suptable1_footnote.tex \
			>>"${_supptable_md}"
		;;
	2*)
		cat "${script_dir}"/polaplib/polap-data-v2-suptable2_footnote.tex \
			>>"${_supptable_md}"
		;;
	3*)
		cat "${script_dir}"/polaplib/polap-data-v2-suptable3_footnote.tex \
			>>"${_supptable_md}"
		;;
	x)
		cat "${script_dir}"/polaplib/polap-data-v2-suptable1_footnote.tex \
			>>"${_supptable_md}"
		;;
	*)
		echo "ERROR: no such stage: ${_arg_stage}"
		exit 1
		;;
	esac

	# printf "\\\elandscape\n\n\\\newpage\n\n" \
	# 	>>"${_supptable_md}"

	printf "\n\n\\\newpage\n\n" \
		>>"${_supptable_md}"

}

supptable1_genus_species() {
	local _arg_outdir="${1:-all}"
	local _arg_inum="${2:-2}"
	local _arg_d_index="${3:-infer-1}"
	local _arg_stage="${4:-1}"
	local _arg_table="${5:-1}"
	local _arg_t_dir="${6:-"${_arg_default_target_dir}"}"

	local _table_name="supptable1"
	if ((_arg_table == 2)); then
		_table_name="supptable2"
	fi
	local _supptable_md="${_table_name}-${_arg_inum}-${_arg_stage}.md"
	if [[ "${_arg_outdir}" != "all" ]]; then
		_supptable_md="${_table_name}-${_arg_inum}-${_arg_stage}-${_arg_outdir}.md"
	fi

	rm -f "${_supptable_md}"

	# maintable1_genus_species_header "${_arg_table}" >"${_table_tsv}"

	if [[ "${_arg_outdir}" == "all" ]]; then

		# Extract and sort keys
		sorted_keys=($(for key in "${!_table1[@]}"; do echo "$key"; done | sort))

		# Iterate over sorted keys and check if value is "T"
		for key in "${sorted_keys[@]}"; do
			if ((_arg_table == 1)); then
				if [[ "${_table1[$key]}" == "F" ]]; then
					continue
				fi
			elif ((_arg_table == 2)); then
				if [[ "${_table2[$key]}" == "F" ]]; then
					continue
				fi
			fi
			local _extracted_inum=${_inum[$key]}
			if [[ "${_extracted_inum}" != "${_arg_inum}" ]]; then
				continue
			fi
			local _v_folder=${_folder[$key]}
			local _v_folder_inum="${_v_folder}/${_arg_inum}"
			if [[ ! -d "${_v_folder_inum}" ]]; then
				continue
			fi
			# skip the main figure
			if [[ " ${Smain[@]} " =~ " ${_v_folder} " ]]; then
				continue
			fi

			supptable1_genus_species_for "${_v_folder}" \
				"${_arg_inum}" "${_arg_d_index}" "${_arg_stage}" \
				"${_arg_table}" "${_arg_t_dir}" 0
		done
	else
		supptable1_genus_species_for "$@" 1
	fi

	echo "See ${_supptable_md}"
	cp -p ${_supptable_md} "${_arg_t_dir}"
}

suppfigure1_genus_species_for() {
	local _arg_outdir="${1:-all}"
	local _arg_inum="${2:-2}"
	local _arg_d_index="${3:-infer-1}"
	local _arg_stage="${4:-1}"
	local _arg_table="${5:-1}"
	local _arg_bandage="${6:-no}"
	local _arg_t_dir="${7:-"${_arg_default_target_dir}"}"
	local _arg_one="${8:-0}"
	local _key="${_arg_outdir}-${_arg_inum}"

	local _table_name="suppfigure1"
	if ((_arg_table == 2)); then
		_table_name="suppfigure2"
	fi
	local _supptable_md="${_table_name}-${_arg_inum}-${_arg_stage}.md"
	if [[ "${_arg_one}" == "1" ]]; then
		_supptable_md="${_table_name}-${_arg_inum}-${_arg_stage}-${_arg_outdir}.md"
	fi

	local _supptable_md="${_supptable_md}"

	local _species=${_taxon[$_key]}
	local _species="${_species//_/ }"
	local _v1=${_folder[$_key]}
	local _label_base="${_v1/_/-}"
	_label_base=$(echo "$_label_base" | awk '{print tolower($0)}')

	output_dir=${_v1}
	echo output_dir: ${output_dir}

	local j=1

	local _v1_inum="${_v1}/${_arg_inum}"
	local _params_txt=${_v1_inum}/disassemble/${_arg_d_index}/params.txt
	_ipn=$(parse_params "${_params_txt}")
	local _I _P _N _R
	read -r _I _P _N _R <<<"$_ipn"
	local _label="${_arg_inum}-${_label_base}"
	if [[ "${_arg_one}" == "1" ]]; then
		_label="main-${_arg_inum}-${_label_base}"
	fi

	printf "\\\newpage\n\n" \
		>>"${_supptable_md}"

	printf "\n\n" \
		>>"${_supptable_md}"

	#########################################################
	# all figures
	local extracted_n="${_compare_n["$_key"]}"
	local i
	local width=13
	local images=()
	local captions=()

	local k=0
	for ((i = 0; i < extracted_n; i++)); do

		local _gfa_infer="${_v1_inum}/disassemble/${_arg_d_index}/${_arg_stage}/${i}/30-contigger/graph_final.gfa"
		local _png_infer="${_v1_inum}/disassemble/${_arg_d_index}/${_arg_stage}/${i}/30-contigger/graph_final.png"
		if [[ -s "${_gfa_infer}" ]]; then
			echo "gfa file: ${_gfa_infer}"
			if [[ "${_arg_bandage}" == "yes" ]]; then
				${_polap_cmd} bandage png \
					${_gfa_infer} \
					${_png_infer}
			fi
			# printf "| ![polap %s](figures/%s){ width=%s%% } " "${i}" "${_png_infer}" "${width}" >>"${_supptable_md}"
			images+=("figures/${_png_infer}")
			captions+=(${i})
			((k++))
		else
			echo "no such file: ${_gfa_infer}"
		fi
	done

	# 	cat <<EOF >>"$_supptable_md"
	#
	# Table: Plastid genome assembly graphs generated from Stage 1 of the subsampling-based method for the dataset of ${_species}. The graphs were generated with Bandage software. Each number corresponds to the iteration index of Stage 1 ([@tbl:supptable1-${_label}]). {#fig:suppfigure1-${_label}}
	#
	# EOF

	# Start writing to the markdown file
	cat <<EOF >>"$_supptable_md"

| Organelle genome assemblies |  |  |
|-----------------|-----------------|-----------------|
EOF

	# Generate the image table with subcaptions
	count=0
	row="| "
	caption_row="| "

	for ((i = 0; i < ${#images[@]}; i++)); do
		row+="![${captions[i]}](${images[i]}){width=100px} | "
		caption_row+="**${captions[i]}** | "
		((count++))

		# End row if 3 images are added
		if ((count % 3 == 0)); then
			echo "$row" >>"$_supptable_md"
			echo "$caption_row" >>"$_supptable_md"
			echo "|-----------------|-----------------|-----------------|" >>"$_supptable_md"
			# echo "" >>"$_supptable_md"
			row="| "
			caption_row="| "
		fi
	done

	# Handle the last incomplete row (if any)
	remaining=$((3 - count % 3))
	if ((remaining < 3)); then
		# Fill empty image columns
		for ((i = 0; i < remaining; i++)); do
			row+=" | "
			caption_row+=" | "
		done
		echo "$row" >>"$_supptable_md"
		echo "$caption_row" >>"$_supptable_md"
		echo "|-----------------|-----------------|-----------------|" >>"$_supptable_md"
	fi

	cat <<EOF >>"$_supptable_md"

![Plastid genome assembly graphs generated from Stage 1 of the subsampling-based method for the dataset of _${_species}_. The graphs were generated with Bandage software. Each number corresponds to the iteration index of Stage 1 of @tbl:supptable1-${_label}.](empty.png){#fig:suppfigure1-${_label}} 

EOF

	echo "" >>"${_supptable_md}"

}

suppfigure1_genus_species() {
	local _arg_outdir="${1:-all}"
	local _arg_inum="${2:-2}"
	local _arg_d_index="${3:-infer-1}"
	local _arg_stage="${4:-1}"
	local _arg_table="${5:-1}"
	local _arg_bandage="${6:-no}"
	local _arg_t_dir="${7:-"${_arg_default_target_dir}"}"

	local _table_name="suppfigure1"
	if ((_arg_table == 2)); then
		_table_name="suppfigure2" # we had figure 2
	fi
	local _supptable_md="${_table_name}-${_arg_inum}-${_arg_stage}.md"
	if [[ "${_arg_outdir}" != "all" ]]; then
		_supptable_md="${_table_name}-${_arg_inum}-${_arg_stage}-${_arg_outdir}.md"
	fi

	rm -f "${_supptable_md}"

	# maintable1_genus_species_header "${_arg_table}" >"${_table_tsv}"

	if [[ "${_arg_outdir}" == "all" ]]; then

		# Extract and sort keys
		sorted_keys=($(for key in "${!_table1[@]}"; do echo "$key"; done | sort))

		# Iterate over sorted keys and check if value is "T"
		for key in "${sorted_keys[@]}"; do
			if ((_arg_table == 1)); then
				if [[ "${_table1[$key]}" == "F" ]]; then
					continue
				fi
			elif ((_arg_table == 2)); then
				if [[ "${_table2[$key]}" == "F" ]]; then
					continue
				fi
			fi
			local _extracted_inum=${_inum[$key]}
			if [[ "${_extracted_inum}" != "${_arg_inum}" ]]; then
				continue
			fi
			local _v_folder=${_folder[$key]}
			local _v_folder_inum="${_v_folder}/${_arg_inum}"
			if [[ ! -d "${_v_folder_inum}" ]]; then
				continue
			fi
			# skip the main figure
			if [[ " ${Smain[@]} " =~ " ${_v_folder} " ]]; then
				continue
			fi

			suppfigure1_genus_species_for "${_v_folder}" \
				"${_arg_inum}" "${_arg_d_index}" "${_arg_stage}" \
				"${_arg_table}" "${_arg_bandage}" "${_arg_t_dir}" 0
		done
	else
		suppfigure1_genus_species_for "$@" 1
	fi

	echo "See ${_supptable_md}"
	cp -p ${_supptable_md} "${_arg_t_dir}"
	echo "rsync figures as well"
}

# figures of Polap, ptGAUL, and GetOrganelle's assemblies

# supfigure2.md
# 2 no
# 2 yes to extract bandage graph png
suppfigure3_genus_species() {
	local _arg_inum="${1:-2}"
	local _arg_d_index="${2:-infer-1}"
	local _brg_mainfigure="${3:-off}"
	local _arg_bandage="${4:-no}"
	local _arg_t_dir="${5:-"${_arg_default_target_dir}"}"

	# local _table_name="suppfigure1"
	# if ((_arg_table == 2)); then
	# 	_table_name="suppfigure2"
	# fi
	# local _supptable_md="${_table_name}-${_arg_inum}-${_arg_stage}.md"

	local _tf_mainfigure="F"
	local _supfigure_file="suppfigure3-${_arg_inum}.md"
	if [[ "${_brg_mainfigure}" == "on" ]]; then
		_tf_mainfigure="T"
		_supfigure_file="mainfigure3-${_arg_inum}.md"
	fi

	# printf "# Supplementary Figures: Polap, ptGAUL, and GetOrganelle\n\n" \
	# 	>"${_supfigure_file}"
	rm -f "${_supfigure_file}"

	# Start writing to the markdown file
	cat <<EOF >>"${_supfigure_file}"

  | Species (Order) | Polap | ptGAUL | GetOrganelle |
|-----------------|-----------------|-----------------|-----------------|
EOF

	# Extract and sort keys
	sorted_keys=($(for key in "${!_table1[@]}"; do echo "$key"; done | sort))

	# Iterate over sorted keys and check if value is "T"
	for key in "${sorted_keys[@]}"; do

		if [[ "${_mainfigure[$key]}" != ${_tf_mainfigure} ]]; then
			continue
		fi
		if [[ "${_table1[$key]}" == "F" ]]; then
			continue
		fi

		# if [[ "${_table1[$key]}" == "F" ]]; then
		# 	continue
		# fi
		# echo "key: $key"

		local _v1=${_folder[$key]}
		local _species=${_taxon[$key]}
		local _species="${_species//_/ }"
		local _genus=${_species%% *}
		_label_base="${_v1/_/-}"
		_label_base=$(echo "$_label_base" | awk '{print tolower($0)}')

		# bandage graph
		extracted_inum=${_inum[$key]}
		if [[ "${extracted_inum}" != "${_arg_inum}" ]]; then
			continue
		fi

		local _order=$(grep "${_genus}" ${script_dir}/polaplib/taxonomy_output.tsv | cut -f 6 | head -1)
		output_dir=${_v1}
		# echo inum: ${extracted_inum}
		# echo output_dir: ${output_dir}
		# echo species: ${_species}
		# echo order: ${_order}

		if [[ "${_arg_bandage}" == "yes" ]]; then
			${_polap_cmd} disassemble bandage \
				-i ${extracted_inum} \
				-o ${output_dir} \
				--disassemble-i infer-1
		fi
		local _gfa_polap="${output_dir}/${extracted_inum}/disassemble/infer-1/pt.1.gfa"
		local _png_polap="${output_dir}/${extracted_inum}/disassemble/infer-1/pt.1.png"

		local images=()
		local captions=()
		images+=("${_species}")
		captions+=("${_order}")
		images+=(${_png_polap})
		captions+=("Polap")

		local _gfa_ptgaul="${output_dir}/ptgaul/flye_cpONT/assembly_graph.gfa"
		local _png_ptgaul="${output_dir}/ptgaul.png"
		if [[ "${_arg_bandage}" == "yes" ]]; then
			${_polap_cmd} bandage png \
				${_gfa_ptgaul} \
				${_png_ptgaul}
		fi
		images+=(${_png_ptgaul})
		captions+=("ptGAUL")

		local _gfa_getorganelle=$(find "${output_dir}/getorganelle" -type f -name 'embplant_pt*.gfa' | head -n 1)
		local _png_getorganelle="${output_dir}/getorganelle.png"
		if [[ "${_arg_bandage}" == "yes" ]]; then
			${_polap_cmd} bandage png \
				${_gfa_getorganelle} \
				${_png_getorganelle}
		fi
		images+=(${_png_getorganelle})
		captions+=("GetOrganelle")

		# echo "Warning: no such file: ${_gfa}"
		# echo "FIX QT5 problem:"
		# problem solved
		# echo "export QT_QPA_PLATFORM=offscreen"
		# Output Markdown file
		#
		output="${_supfigure_file}"

		# Generate the image table with subcaptions
		count=0
		row="| "
		caption_row="| "

		for ((i = 0; i < ${#images[@]}; i++)); do
			if [[ ${i} == "0" ]]; then
				row+="_${images[i]}_ | "
			else
				row+="![${captions[i]}](figures/${images[i]}){width=100px} | "
			fi
			caption_row+="**${captions[i]}** | "
			((count++))

			# End row if 3 images are added
			if ((count % 4 == 0)); then
				echo "$row" >>"$output"
				echo "$caption_row" >>"$output"
				echo "|-----------------|-----------------|-----------------|-----------------|" >>"$output"
				# echo "" >>"$output"
				row="| "
				caption_row="| "
			fi
		done

	done

	cp -p ${_supfigure_file} "${_arg_t_dir}"
	echo ${_supfigure_file}
}

function _run_polap_menu { # Interactive menu interface
	local species_key=0

	# Function to display the main menu
	show_main_menu() {
		echo "===================="
		echo "     Main Menu      "
		echo "===================="
		echo "1. BioProject Folder"
		echo "2. BioProject-species"
		echo "E. Exit"
		echo "--------------------------------------------------"
		local _species_name="${keys_array[$species_key]}"
		local _remote_name="${_ssh["$_species_name"]}"
		local _status_name="${_status["$_species_name"]}"
		echo "POLAP Species Folder ([$species_key] ${keys_array[$species_key]} [$_remote_name] [$_status_name])"
		echo "===================="
		echo -n "Enter your choice [2] (or just press ENTER): "
	}

	# Function to display the species operations submenu
	show_bioproject_menu() {
		echo "===================="
		echo "  BioProject"
		echo "===================="
		for i in "${!keys_array[@]}"; do
			local _species_index=${keys_array[$i]}
			echo "$i. ${keys_array[$i]} (${_ssh[$_species_index]}: ${_status[$_species_index]})"
		done
		echo "M. Back to Main Menu"
		echo "--------------------------------------------------"
		local _species_name="${keys_array[$species_key]}"
		local _remote_name="${_ssh["$_species_name"]}"
		local _status_name="${_status["$_species_name"]}"
		echo "POLAP Species Folder ([$species_key] ${keys_array[$species_key]} [$_remote_name] [$_status_name])"
		echo "===================="
		echo -n "Enter your choice: "
	}

	# Function to display the system operations submenu
	#
	# BioProject folder has species name appended
	#
	#
	# Download BioProject runinfo (creating folders for all species)
	# Download
	show_species_menu() {
		local _species_name="${keys_array[$species_key]}"
		local _remote_name="${_ssh["$_species_name"]}"
		local _species_folder="${keys_array[$species_key]}"
		local _status_name="${_status["$_species_name"]}"
		echo "===================="
		echo "  BioProject-species"
		echo "===================="
		echo "0. Batch"
		echo "1. Send data to the remote: ${_remote_name}"
		echo "2. Make a starting folder (3)"
		echo "3. Get ptDNA from NCBI"
		echo "4. GetOrganelle assembly with the short-read data"
		echo "5. ptGAUL assembly with the reference (6)"
		echo "6. Prepare the short-read polishing"
		echo "7. Extract the ptGAUL assembly (8)"
		echo "8. Copy the ptGAUL assembly"
		echo "9. Infer with POLAP"
		echo "10. Check with POLAP (needs step 5,7)"
		echo "11. (only for 5 datasets) Compare with POLAP (needs step 5,7)"
		echo "L. List the folder: ${_species_folder}"
		echo "S. Estimate the coverage ${_species_folder}"
		echo "A. Archive the ${_species_folder}"
		echo "G. Get the result from the remote"
		echo "C. Clean up the ${_species_folder}"
		echo "M. Back to Main Menu"
		echo "--------------------------------------------------"
		echo current host: "$(hostname)"
		echo rement host name: "${_host["$_species_name"]}"
		echo "POLAP Species Folder ([$species_key] ${keys_array[$species_key]} [$_remote_name] [$_status_name])"
		echo "===================="
		# echo -n "Enter your choice: "
	}

	# Function to navigate between menus
	navigate_menu() {
		local current_menu="$1"
		while true; do
			date +"%Y-%m-%d %H:%M:%S"
			case "$current_menu" in
			main)
				show_main_menu
				read -r choice
				if [[ -z "${choice}" ]]; then
					choice=2
				fi
				# Convert the input to lowercase for case-insensitive comparison
				local choice="${choice,,}"
				case $choice in
				1) current_menu="bioproject" ;;
				2) current_menu="species" ;;
				e | q)
					echo "Exiting..."
					exit 0
					;;
				*) echo "Invalid choice, please select a valid option." ;;
				esac
				;;
			bioproject)
				show_bioproject_menu
				read -r choice
				# Convert the input to lowercase for case-insensitive comparison
				local choice="${choice,,}"
				case $choice in
				m | e) current_menu="main" ;; # Go back to Main Menu
				A)
					echo "Case A detected"
					;;
				'' | *[!0-9]*)
					echo "Invalid number"
					;;
				*)
					current_menu="species"
					species_key=$choice
					;;
				esac
				;;
			species)
				show_species_menu
				read -r -p "Enter your chioces: " input
				for choice in $input; do
					case $choice in
					0)
						# src/polap.sh -o "${_arg_outdir}" get-mtdna file

						echo "$species_key: ${keys_array[$species_key]}"
						batch_genus_species "${keys_array[$species_key]}"
						# current_menu="bioproject"

						;;
					1)
						send-data-to_genus_species "${keys_array[$species_key]}"
						;;
					2)
						mkdir_genus_species "${keys_array[$species_key]}"
						;;
					3)
						get-ptdna-from-ncbi_genus_species "${keys_array[$species_key]}"
						;;
					4)
						getorganelle_genus_species \
							"${keys_array[$species_key]}"
						;;
					5)
						ptgaul_genus_species \
							"${keys_array[$species_key]}"
						;;
					6)
						msbwt_genus_species \
							"${keys_array[$species_key]}"
						;;
					7)
						extract-ptdna-of-ptgaul_genus_species \
							"${keys_array[$species_key]}"
						;;
					8)
						copy-ptdna-of-ptgaul_genus_species \
							"${keys_array[$species_key]}"
						;;
					9)
						infer_genus_species \
							"${keys_array[$species_key]}"
						infer_genus_species \
							"${keys_array[$species_key]}" \
							--disassemble-simple-polishing
						;;
					10)
						check_genus_species \
							"${keys_array[$species_key]}"
						check_genus_species \
							"${keys_array[$species_key]}" \
							--disassemble-simple-polishing
						;;
					11)
						compare_genus_species \
							"${keys_array[$species_key]}"
						;;
					l)
						ls -l \
							"${keys_array[$species_key]}"
						;;
					g)
						get_genus_species \
							"${keys_array[$species_key]}"
						;;
					a)
						archive_genus_species \
							"${keys_array[$species_key]}"
						;;
					s)
						downsample_genus_species \
							"${keys_array[$species_key]}" \
							5 --dry
						;;
					c)
						clean_genus_species \
							"${keys_array[$species_key]}"
						;;
					m | e) current_menu="main" ;; # Go back to Main Menu
					*) echo "Invalid choice, please select a valid option." ;;
					esac
				done
				;;
			*)
				echo "Unknown menu. Returning to main menu."
				current_menu="main"
				;;
			esac
			echo
		done
	}

	# Start at the main menu
	# navigate_menu "main"
	navigate_menu "bioproject"
}

if [[ "${subcmd1}" == "help" ]]; then
	subcmd1="${_arg2}"
	_arg2="arg2"
fi

# Main case statement
case "$subcmd1" in
'system' | \
	'copy-figures')
	${subcmd1}_genus_species
	;;
'refs')
	get-ptdna-from-ncbi_genus_species "${_arg2}"
	;;
'send-data' | \
	'send' | \
	'local-batch' | \
	'send-archive' | \
	'recover' | \
	'mkdir' | \
	'sra' | \
	'get-ptdna-from-ncbi' | \
	'copy-ptdna-of-ncbi-as-reference' | \
	'getorganelle' | \
	'ptgaul' | \
	'msbwt' | \  | \
	'extract-ptdna-of-ptgaul' | \
	'extract-ptdna-of-ptgaul2' | \
	'copy-ptdna-of-ptgaul' | \
	'polishing' | \
	'wga' | \
	'polish' | \
	'restart' | \
	'clean-infer' | \
	'write-config' | \
	'table1' | \
	'table2' | \
	'single-argument')
	${subcmd1}_genus_species "${_arg2}"
	;;
'supfigure2' | \
	'compare' | \
	'use-downsample' | \
	'best' | \
	'mauve' | \
	bandage* | \
	'coverage' | \
	'delete-disassemble' | \
	'two-arguments')
	${subcmd1}_genus_species "${_arg2}" "${_arg3}"
	;;
'infer' | \
	'infer2' | \
	'infer3only' | \
	'downsample2infer' | \
	'suptable1' | \
	'supfigure1' | \
	'check')
	${subcmd1}_genus_species "${_arg2}" "${_arg3}" "${_arg4}"
	# ${subcmd1}_genus_species "${_arg2}" --disassemble-simple-polishing
	;;
'downsample')
	${subcmd1}_genus_species "${_arg2}" "${_arg3}" "${_arg4}"
	;;
'clean')
	if [[ "${_arg2}" == arg2 ]]; then
		echo "Help: ${subcmd1} <outdir> [confirm:off|on]"
		echo "  polap-data-v2.sh ${subcmd1} Arabidopsis_thaliana"
		echo "  polap-data-v2.sh ${subcmd1} Arabidopsis_thaliana on"
		exit 0
	fi
	[[ "${_arg3}" == arg3 ]] && _arg3=""
	[[ "${_arg4}" == arg4 ]] && _arg4=""
	${subcmd1}_genus_species "${_arg2}" "${_arg3}"
	;;
'archive')
	if [[ "${_arg2}" == arg2 ]]; then
		echo "Help: ${subcmd1} <outdir> <inum:-1|N>"
		echo "  polap-data-v2.sh ${subcmd1} Arabidopsis_thaliana"
		echo "  polap-data-v2.sh ${subcmd1} Arabidopsis_thaliana 0"
		echo "${help_message_archive}"
		exit 0
	fi
	[[ "${_arg3}" == arg3 ]] && _arg3=""
	${subcmd1}_genus_species "${_arg2}" "${_arg3}"
	;;
'batch')
	if [[ "${_arg2}" == arg2 ]]; then
		echo "Help: ${subcmd1} <outdir> <inum:N> <ref:off|on> <confirm:off|on> [redo:off|on]"
		echo "  polap-data-v2.sh ${subcmd1} Arabidopsis_thaliana 2 on off off"
		echo "  polap-data-v2.sh ${subcmd1} Arabidopsis_thaliana 2 on on off"
		echo "  polap-data-v2.sh ${subcmd1} Arabidopsis_thaliana 2 off on on"
		echo "${help_message_batch}"
		exit 0
	fi
	[[ "${_arg3}" == arg3 ]] && _arg3=""
	[[ "${_arg4}" == arg4 ]] && _arg4=""
	[[ "${_arg5}" == arg5 ]] && _arg5=""
	[[ "${_arg6}" == arg6 ]] && _arg6=""
	${subcmd1}_genus_species "${_arg2}" "${_arg3}" "${_arg4}" "${_arg5}" "${_arg6}"
	;;
'remote-batch')
	if [[ "${_arg2}" == arg2 ]]; then
		echo "Help: ${subcmd1} <outdir|all> <inum:N> <per_host:off|on>"
		echo "  polap-data-v2.sh ${subcmd1} Arabidopsis_thaliana 0 off"
		echo "  polap-data-v2.sh ${subcmd1} all 0"
		echo "  polap-data-v2.sh ${subcmd1} all 2 on"
		echo "${help_message_remote_batch}"
		exit 0
	fi
	[[ "${_arg3}" == arg3 ]] && _arg3=""
	[[ "${_arg4}" == arg4 ]] && _arg4=""
	${subcmd1}_genus_species "${_arg2}" "${_arg3}" "${_arg4}"
	;;
'get')
	if [[ "${_arg2}" == arg2 ]]; then
		echo "Help: ${subcmd1} <outdir|all> <inum:N> <confirm-yes:off|yes|add> <download:on|off>"
		echo "  polap-data-v2.sh ${subcmd1} Arabidopsis_thaliana -1 off off"
		echo "  polap-data-v2.sh ${subcmd1} Arabidopsis_thaliana 0 off"
		echo "  polap-data-v2.sh ${subcmd1} all 2 add off"
		echo "${help_message_get}"
		exit 0
	fi
	[[ "${_arg3}" == arg3 ]] && _arg3=""
	[[ "${_arg4}" == arg4 ]] && _arg4=""
	[[ "${_arg5}" == arg5 ]] && _arg5=""
	${subcmd1}_genus_species "${_arg2}" "${_arg3}" "${_arg4}" "${_arg5}"
	;;
'report')
	if [[ "${_arg2}" == arg2 ]]; then
		echo "Help: ${subcmd1} <outdir|all> <inum:N> <disassemble index:infer-1>"
		echo "  polap-data-v2.sh ${subcmd1} all 2 infer-1"
		echo "  polap-data-v2.sh ${subcmd1} Arabidopsis_thaliana 2 infer-1"
		echo "${help_message_report}"
		exit 0
	fi
	[[ "${_arg3}" == arg3 ]] && _arg3=""
	[[ "${_arg4}" == arg4 ]] && _arg4=""
	${subcmd1}_genus_species "${_arg2}" "${_arg3}" "${_arg4}"
	;;
'maintable1')
	if [[ "${_arg2}" == arg2 ]]; then
		echo "Help: ${subcmd1} <outdir|all> <inum:N> <disassemble index:infer-1> <table:1|2> <targe dir>"
		echo "  polap-data-v2.sh maintable1 all 2 infer-1"
		echo "  polap-data-v2.sh maintable1 all 0 infer-1"
		echo "  polap-data-v2.sh maintable1 Eucalyptus_pauciflora 2 infer-1"
		echo "${help_message_maintable1}"
		exit 0
	fi
	[[ "${_arg3}" == arg3 ]] && _arg3=""
	[[ "${_arg4}" == arg4 ]] && _arg4=""
	[[ "${_arg5}" == arg5 ]] && _arg5=""
	[[ "${_arg6}" == arg6 ]] && _arg6=""
	${subcmd1}_genus_species "${_arg2}" "${_arg3}" "${_arg4}" "${_arg5}" "${_arg6}"
	;;
'supptable1')
	if [[ "${_arg2}" == arg2 ]]; then
		echo "Help: ${subcmd1} <outdir|all> <inum:N> <disassemble index:infer-1> <stage:x|1|2|3-infer> <table:1|2> <targe dir>"
		echo "  polap-data-v2.sh supptable1 all 2 infer-1 x"
		echo "  polap-data-v2.sh supptable1 Eucalyptus_pauciflora 2 infer-1 x"
		echo "  polap-data-v2.sh supptable1 Eucalyptus_pauciflora 0 infer-1 x"
		echo "  polap-data-v2.sh supptable1 Eucalyptus_pauciflora 8 infer-1 x"
		echo "${help_message_supptable1}"
		exit 0
	fi
	[[ "${_arg3}" == arg3 ]] && _arg3=""
	[[ "${_arg4}" == arg4 ]] && _arg4=""
	[[ "${_arg5}" == arg5 ]] && _arg5=""
	[[ "${_arg6}" == arg6 ]] && _arg6=""
	[[ "${_arg7}" == arg7 ]] && _arg7=""
	${subcmd1}_genus_species "${_arg2}" "${_arg3}" "${_arg4}" "${_arg5}" "${_arg6}" "${_arg7}"
	;;
'suppfigure1')
	if [[ "${_arg2}" == arg2 ]]; then
		echo "Help: ${subcmd1} <outdir> <inum> <disassemble index> <stage> <table:1|2> <bandage:no|yes> <targe dir>"
		echo "  polap-data-v2.sh suppfigure1 all 2 infer-1 1 1 yes"
		echo "  polap-data-v2.sh suppfigure1 Eucalyptus_pauciflora 2 infer-1 1 1 yes"
		exit 0
	fi
	[[ "${_arg3}" == arg3 ]] && _arg3=""
	[[ "${_arg4}" == arg4 ]] && _arg4=""
	[[ "${_arg5}" == arg5 ]] && _arg5=""
	[[ "${_arg6}" == arg6 ]] && _arg6=""
	[[ "${_arg7}" == arg7 ]] && _arg7=""
	[[ "${_arg8}" == arg8 ]] && _arg8=""
	${subcmd1}_genus_species "${_arg2}" "${_arg3}" "${_arg4}" "${_arg5}" "${_arg6}" "${_arg7}" "${_arg8}"
	;;
'suppfigure3')
	if [[ "${_arg2}" == arg2 ]]; then
		echo "Help: ${subcmd1} <inum> <disassemble index> <main:off|on> <bandage:no|yes> <targe dir>"
		echo "  polap-data-v2.sh suppfigure3 2 infer-1 on yes"
		echo "  polap-data-v2.sh suppfigure3 2 infer-1 off yes"
		exit 0
	fi
	[[ "${_arg3}" == arg3 ]] && _arg3=""
	[[ "${_arg4}" == arg4 ]] && _arg4=""
	[[ "${_arg5}" == arg5 ]] && _arg5=""
	[[ "${_arg6}" == arg6 ]] && _arg6=""
	${subcmd1}_genus_species "${_arg2}" "${_arg3}" "${_arg4}" "${_arg5}" "${_arg6}"
	;;
"menu")
	_run_polap_menu
	;;
"install-conda")
	mkdir -p ~/miniconda3
	wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
	bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
	rm ~/miniconda3/miniconda.sh
	echo After installing Miniconda3, close and reopen your terminal application.
	;;
"setup-conda")
	source ~/miniconda3/bin/activate
	conda config --add channels bioconda
	conda config --add channels conda-forge
	conda config --set channel_priority strict
	;;
"install-polap" | "install")
	conda create --name polap bioconda::polap
	;;
"download-polap")
	wget https://github.com/goshng/polap/archive/refs/tags/${_polap_version}.zip
	unzip -o -q ${_polap_version}.zip
	cd polap-${_polap_version}
	;;
"clean")
	rm -f ${_polap_version}.zip
	rm -rf polap-${_polap_version}
	;;
"install-fmlrc")
	wget https://github.com/goshng/polap/archive/refs/tags/${_polap_version}.zip
	unzip -o -q ${_polap_version}.zip
	cd polap-${_polap_version}
	conda env create -f src/polap-conda-environment-fmlrc.yaml
	;;
"patch-polap" | "patch")
	wget https://github.com/goshng/polap/archive/refs/tags/${_polap_version}.zip
	unzip -o -q ${_polap_version}.zip
	cd polap-${_polap_version}/src
	bash polap-build.sh >../build.sh
	cd ..
	PREFIX="$HOME/miniconda3/envs/polap" bash build.sh
	;;
"test-polap")
	cd polap-${_polap_version}
	cd test
	source ~/miniconda3/bin/activate polap
	polap assemble --test
	;;
"uninstall")
	# source ~/miniconda3/bin/deactivate
	conda remove -n polap --all
	conda remove -n polap-fmlrc --all
	;;
"install-getorganelle")
	conda create --name getorganelle bioconda::getorganelle
	conda activate getorganelle
	get_organelle_config.py --add embplant_pt,embplant_mt
	;;
"mkdir-all")
	for i in "${Sall[@]}"; do
		echo "creating folder $i ..."
		mkdir ${i}
	done
	;;
"rm")
	read -p "Do you want to delete all? (y/N): " confirm
	case "$confirm" in
	[yY] | [yY][eE][sS])
		for i in "${Sall[@]}"; do
			echo "deleting folder $i ..."
			rm -rf ${i}
		done
		;;
	*)
		echo "Deleting all is canceled."
		;;
	esac
	;;
*)
	echo "Usage: $0 <subcommand> [species_folder]"
	echo "${help_message}"
	echo "subcommand '$subcmd1' is not recognized."
	exit 1
	;;
esac
