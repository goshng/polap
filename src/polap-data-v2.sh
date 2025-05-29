#!/usr/bin/env bash
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

# polap and polap-data-cflye are release-versions with _POLAP_RELEASE to be 1.
: "${_POLAP_DEBUG:=0}"
export _POLAP_DEBUG
: "${_POLAP_RELEASE:=0}"
export _POLAP_RELEASE

# Data directories: we download data from the NCBI SRA database
# unless they exist in the following folders.
_local_host="thorne"
_media1_dir="/media/h1/sra"
_media2_dir="/media/h2/sra"
_media_dir="/media/h2/sra"

# TODO:
# [ ] A manual for proper understanding and execution of tasks.
# [ ] assembly graphs for all plastid genomes
# [ ] species list and summary table: monocots, eudicots, mosses, liverwort, ...
# [ ] genome assembly results from other tools
# [ ] consider what we should do if no assembly is built in Stage 1 or 2.
#
# 2025-05-29
# [x] Sall read from the CSV file
# [x] cflye conda env polap-cflye execution
# [x] benchmarking: GetOrganelle, ptGAUL, PMAT, TIPPo, oatk

# NOTE: we might consider how we change the options if no assembly is built.
#
# TODO: we need to consider much more than just adding this, we wil consider
# Stage 1: -n => -sn, -a -sa -sap, -p => -sb -sap
# Stage 2: -r => -rn, -ra -rap, -rb -rbp
# Stage 3: -l -ln, -la -lap, -lb -lbp
# CSV: add sa, sb, and others more
# This is another idea that might need to consider:
# we might use the second candidate if Stage 2 produces no assembly.
# Then, we need to move or copy the check of candidate front of Stage 2
# not at the end of Stage 1.
# In Stage 1, we could do something else if it does not produce assembly.
# Although it may not be obvious, we might increase the max subsampling rate.
# Or, we could increase the max memory limit. Something we already mentioned
# in the manuscript. To do such thing, we need to move around the candidate
# checking part.
# How-To
#
# Edit "${_POLAPLIB_DIR}/polap-data-v2.csv"
# It looks for the CSV file in the following order.
# 1. The current folder where you execute this script
# 2. The folder in which this script resides
#
# Eucalyptus_pauciflora is set to the main result.
# If needed, this should be changed by editing Smain array below.
#
# _polap_cmd=${_polap_cmd}
# _media_dir=${_media_dir}
# _local_host=${_local_host}
# host=$(hostname)
#
# Development
# 1. Add a function named test_genus_species
# 2. Add its caller in a case statement: search for 'test'
# 3. Add a function named test_genus_species_for to iterate over all datasets.

# Use man folder for a release-version
# otherwise use a custom folder.
_brg_default_target_dir="$HOME/all/manuscript/polap-v0.4/figures"
_brg_default_target_dir="man/v0.4/figures"
if [[ -d "src" || -d "polap" ]]; then
  _brg_default_target_dir="$HOME/all/manuscript/polap-v0.4/figures"
else
  if [[ -d "man" ]]; then
    _brg_default_target_dir="man/v0.4/figures"
  fi
  mkdir -p "${_brg_default_target_dir}"
fi

help_message=$(
  cat <<HEREDOC
usage: polap-data-cflye [-h] [-y] [-c CSV] [--version] COMMAND [-h] ...

polap-data-cflye is a tool for data analysis of subsampling-based plastid genome assembly

options:
  -h, --help          Show this help message and exit.
  -y                  Enable -y flag to say YES to any question.
  -c <arg>            Set value for -c option (default: off)
  -t <arg>            Set value for -t option (default: t1)
  -m <arg>            Set value for -m option figure folder (default: ${_brg_default_target_dir})
  -e <ame>            Call <name>_genus_species function and exit.
  -v, --verbose       [Not Yet] Can be used multiple times. Once for detailed output,
                      twice for INFO logging, thrice for _POLAP_DEBUG logging, four
                      times for TRACE logging.
  --version           Show the conda version number and exit.

commands:
  The following commands are available for the script.

  COMMAND
    install            Install a list of tools to some conda environments.
    setup              Setup installed tools.
    update             Update tools.
    list (search)      List tools.
    run                Run an executable in a conda environment.
    remove (uninstall) Remove a list of tools.
    benchmark          Benchmark GetOrganelle, ptGAUL, PMAT, TIPPo, and Oatk.
    get                Get results.
    archive            Archive results.
    man                Generate reports.
    help               Print help message for commands and others.

    batch              Batch run per data. (not implemented yet!)
    build (assemble)   Build plastid genomes.
    compare            Compare the performance of two tools.
    clean              Remove unnecessary results.
    download           Download data.
    delete             Delete data.
    print-help-all     List all help messages.
HEREDOC
)

# Tesing groups of datasets
declare -a Sall

Smain=(
  'Eucalyptus_pauciflora'
)

# Data28
Stest=(
  Euonymus_alatus
  Vaccinium_vitis_idaea
  Vitis_vinifera
  Test_species
  Taxon_genus
)

Ssome=(
  Anthoceros_agrestis
)

S28=(
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
  Vaccinium_vitis_idaea
  Vitis_vinifera
)

# Data28
S7=(
  Anthoceros_agrestis
  Arabidopsis_thaliana
  Canavalia_ensiformis
  Eucalyptus_pauciflora
  Spirodela_polyrhiza
  Vaccinium_vitis-idaea
  Vitis_vinifera
)

# OTHER STUFF GENERATED BY Argbash
_polap_script_bin_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" || {
  echo "Couldn't determine the script's running directory, which probably matters, bailing out" >&2
  exit 2
}
_POLAPLIB_DIR="${_polap_script_bin_dir}/polaplib"

# Target version is 0.4.4
source "${_POLAPLIB_DIR}/polap-git-hash-version.sh"
_polap_version=v0.4.3.7-"${_polap_git_hash_version}"
if [ -z "${_polap_version+x}" ]; then
  _polap_version="0.4.3.7"
fi

source "${_POLAPLIB_DIR}/polap-lib-conda.sh"
source "${_POLAPLIB_DIR}/polap-lib-timing.sh"
source "${_POLAPLIB_DIR}/polap-lib-unit.sh"
source "${_POLAPLIB_DIR}/polap-lib-array.sh"
source "${_POLAPLIB_DIR}/polap-lib-number.sh"
source "${_POLAPLIB_DIR}/polap-lib-data.sh"
source "${_POLAPLIB_DIR}/polap-lib-file.sh"
source "${_POLAPLIB_DIR}/polap-lib-process.sh"
source "${_POLAPLIB_DIR}/polap-lib-extract.sh"
source <(echo 'export PATH="$PWD/bin:$PATH"')
# source <(echo 'export QT_QPA_PLATFORM=offscreen')
source <(echo 'export QT_QPA_PLATFORM=minimal')
# echo "export QT_QPA_PLATFORM=offscreen"

: "${_brg_outdir:=.}"
_polap_data_csv="$(basename "$0" .sh).csv"
_polap_data_data="$(basename "$0" .sh).data"
_polap_data_txt="$(basename "$0" .sh).txt"

_log_echo() {
  if [[ -d "${_brg_outdir}" ]]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') [$subcmd1] - $1" >>"${_brg_outdir}/${_polap_data_txt}"
  else
    echo "$(date '+%Y-%m-%d %H:%M:%S') [$subcmd1] - $1" >>"./${_polap_data_txt}"
  fi
  echo "$(date '+%Y-%m-%d %H:%M:%S') [$subcmd1] - $1"
}

setup-csv_genus_species() {
  local _brg_data="${1:-${_polap_data_csv}}"

  local _data=$(
    cat <<HEREDOC
species,long,short
Anthoceros_agrestis,SRR10190639,SRR10250248
Arabidopsis_thaliana,ERR2173373,ERR2173372
Canavalia_ensiformis,SRR18714551,SRR18714547
Cinchona_pubescens,SRR20784020,SRR20784021
Codonopsis_lanceolata,SRR11585869,SRR11585868
Cucumis_sativus_var_hardwickii,SRR28091980,SRR28091977
Dioscorea_japonica,SRR16108312,SRR16108386
Dunaliella_tertiolecta,SRR22857204,SRR22857205
Eucalyptus_pauciflora,SRR7153095,SRR7161123
Euonymus_alatus,SRR16133411,SRR16122871
Gossypium_herbaceum,SRR17224483,SRR17211914
Juncus_effusus,SRR14298760,SRR14298746
Juncus_inflexus,SRR14298751,SRR14298745
Juncus_roemerianus,SRR21976090,SRR21976092
Juncus_validus,SRR21976089,SRR21976091
Leiosporoceros_dussii,SRR25387688,SRR25387689
Macadamia_jansenii,SRR11191910,SRR11191912
Musa_acuminata_subsp_malaccensis,ERR5455028,ERR3606950
Notothylas_orbicularis,SRR25405055,SRR25405056
Ophrys_lutea,ERR5167480,ERR5303530
Oryza_rufipogon,SRR12104676,SRR12102351
Phaeomegaceros_chiloensis,SRR25430413,SRR25430414
Populus_x_sibirica,SRR15146668,SRR12963707
Prunus_mandshurica,ERR4656977,ERR4656978
Solanum_lycopersicum,SRR11073833,SRR11205474
Spirodela_polyrhiza,SRR11472010,SRR11472009
Vaccinium_vitis-idaea,SRR25468450,SRR25477290
Vitis_vinifera,SRR26163227,SRR26163231
HEREDOC
  )

  local confirm="no"
  if [[ -s "${_brg_data}" ]]; then
    read -p "Do you want to replace ${_brg_data}? (y/N): " confirm
  fi

  if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
    # Write the output to CSV
    {
      # Print the header
      echo species,taxon,long,short,host,random,down,memory,p,n,r,alpha,delta,ptgaul,pmat,tippo,oatk,dummy,status

      local seed=101
      # Skip the header and read each line
      while IFS=',' read -r species long short; do
        # Skip header
        if [[ "$species" == "species" ]]; then
          continue
        fi
        echo "$species,$species,$long,$short,hostname,$seed,10,16,10,10,5,1.0,0.25,160000,0.1,ont,30,dummy,done"
        seed=$((seed + 2))
      done <<<"$_data"
    } >"${_brg_data}"
    _log_echo "A new CSV: ${_brg_data}"
  else
    echo "Canceled: A new CSV: ${_brg_data}"
  fi
}

############################################################
# main command arguments used before a subcommand
#
# Default values for options
opt_c_arg="off"
opt_t_arg="t1"
opt_m_arg="off"
opt_y_flag=false
opt_e_arg=""

print_help() {
  echo "${help_message}"
}

print_version() {
  echo "polap-data-cflye ${_polap_version}"
}

# Parse options
while [[ "${1-}" == -* ]]; do
  case "$1" in
  -c)
    shift
    if [[ -z "${1-}" || "${1-}" == -* ]]; then
      echo "Error: -c requires an argument"
      exit 1
    fi
    opt_c_arg="$1"
    ;;
  -t)
    shift
    if [[ -z "${1-}" || "${1-}" == -* ]]; then
      echo "Error: -t requires an argument"
      exit 1
    fi
    opt_t_arg="$1"
    ;;
  -m)
    shift
    if [[ -z "${1-}" || "${1-}" == -* ]]; then
      echo "Error: -t requires an argument"
      exit 1
    fi
    opt_m_arg="$1"
    ;;
  -y)
    opt_y_flag=true
    ;;
  -e)
    shift
    if [[ -z "${1-}" || "${1-}" == -* ]]; then
      echo "Error: -e requires an argument"
      exit 1
    fi
    opt_e_arg="$1"
    func_name="${opt_e_arg}_genus_species"
    if declare -f "$func_name" >/dev/null; then
      "$func_name"
      exit 0
    else
      echo "Error: function '$func_name' not found"
      exit 1
    fi
    ;;
  --version)
    print_version
    exit 0
    ;;
  -h | --help)
    print_help
    exit 0
    ;;
  --) # End of options
    shift
    break
    ;;
  -*)
    echo "Unknown option: $1"
    exit 1
    ;;
  esac
  shift || break
done

if [[ "${opt_m_arg}" != "off" ]]; then
  _brg_default_target_dir="${opt_m_arg}"
fi

if [[ "${opt_c_arg}" != "off" ]]; then
  csv_file="${opt_c_arg}"
fi

# Input parameter
subcmd1="${1:-help}"

if [[ "${subcmd1}" == "search" ]]; then
  subcmd1="list"
fi

# TODO: depending on _POLAP_RELEASE, we could use polap
# but then polap needs to be the current conda env, which is not.
_polap_cmd="${_polap_script_bin_dir}/polap.sh"
if [[ "${_POLAP_RELEASE}" == "1" ]]; then
  _polap_cmd="${_polap_script_bin_dir}/polap"
else
  _polap_cmd="${_polap_script_bin_dir}/polap.sh"
fi

help_message_requirement=$(
  cat <<HEREDOC
Polap data analysis of subsampling-based plastid genome assembly

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
species,taxon,long,short,host,random,down,memory,p,n,r,alpha,delta,ptgaul,pmat,tippo,oatk,dummy,status

The following descriptions outline the column items present in the CSV file.

  - species: index, <genus>_<species>-<inum>
  - taxon: species name, <genus>_<species>
  - long: long-read NCBI SRA accession
  - short: short-read NCBI SRA accession
  - host: the remote comuputer host name echoed by hostname command
  - random: seed for random number generation
  - down: the downsample coverage
  - memory: the maximum memory to use in polap's disassemble
  - p: the maximum rate in the first stage of polap's disassemble
  - n: the number of steps in the first stage of polap's disassemble
  - r: the number of steps in the second and third stages
  - alpha: NA or any number
  - delta: NA or any number
  - ptgaul: the genome size used by ptGAUL assembly
  - pmat: -fc option e.g., 0.1
  - tippo: e.g., ont onthq
  - oatk: -c option e.g., 30
  - dummy: dummy always
  - status: comment any string

### Variables

At the beginning of the script, several settings require initial setup before proceeding.

- Smain: some of the results go to the main not supplementary materials.
- _brg_default_target_dir: folder where you write a manuscript using tables and figures.

HEREDOC
)

help_message_development=$(
  cat <<HEREDOC
For development only

List of subcommands:
  local-edge-polap
  remote-get-data
  install-tippo,
  download-pmat, install-pmat, pmat
  download-man, install-man, man
HEREDOC
)

help_message_man=$(
  cat <<HEREDOC

  man table-benchmark test 0
  man table-benchmark test 0 benchmark-memory
  man table-benchmark test 0 benchmark-polap
  man table-benchmark test 0 benchmark-time
  man table-benchmark test 0 data
  man table-data test 0
  man table-polap-disassemble Test_species
  man figure-alpha Eucalyptus_pauciflora
  man figure-delta Eucalyptus_pauciflora
HEREDOC
)

help_message_extract_ptgaul_ptdna=$(
  cat <<HEREDOC

  Extract ptDNA sequence from ptGAUL result and save it in a FASTA file.

check <outdir>/t1/0/ptgaul/flye_cpONT/assembly_graph.gfa using Bandage
extract path_sequence.fasta
at <outdir>/t1/0/ptgaul/flye_cpONT/ptdna: ln -s path_sequence.fasta circular_path_1_concatenated.fa
conda activate polap-fmlrc
fmlrc -p 1 <outdir>/t1/0/msbwt/comp_msbwt.npy <outdir>/t1/0/ptgaul/flye_cpONT/ptdna/circular_path_1_concatenated.fa <outdir>/t1/0/ptgaul/flye_cpONT/ptdna/pt.1.fa

e.g., fmlrc -p 1 Cinchona_pubescens/t1/0/msbwt/comp_msbwt.npy Cinchona_pubescens/t1/0/ptgaul/flye_cpONT/ptdna/circular_path_1_concatenated.fa Cinchona_pubescens/t1/0/ptgaul/flye_cpONT/ptdna/pt.1.fa

cp <outdir>/t1/0/ptgaul/flye_cpONT/ptdna/pt.1.fa <outdir>/t1/0/ptdna-ptgaul.fa
conda deactivate
HEREDOC
)

help_message_sample_csv=$(
  cat <<HEREDOC

  Copy all or parts of polap-data-v2.csv.

  1: Spirodela_polyrhiza
  2: Spirodela_polyrhiza and Eucalyptus_pauciflora
  all: all
  force: delete the CSV file.
  inum: used with each option. 
    select all with -0 key but replacing the key with this inum
HEREDOC
)

##### INSERT_HELP_HERE #####
help_message_man_figure_delta=$(
  cat <<HEREDOC

  menu title
HEREDOC
)

help_message_man_figure_alpha=$(
  cat <<HEREDOC

  menu title
HEREDOC
)

help_message_man_table_data=$(
  cat <<HEREDOC

  menu title
HEREDOC
)

help_message_man_figure_polap=$(
  cat <<HEREDOC

  menu title
HEREDOC
)

help_message_man_figure_sheet_pmat=$(
  cat <<HEREDOC

  menu title
HEREDOC
)

help_message_man_figure_sheet_tippo=$(
  cat <<HEREDOC

  menu title
HEREDOC
)

help_message_man_figure_sheet_oatk=$(
  cat <<HEREDOC

  menu title
HEREDOC
)

help_message_man_figure_sheet_polap=$(
  cat <<HEREDOC

  menu title
HEREDOC
)

help_message_man_figure_sheet=$(
  cat <<HEREDOC

  menu title
HEREDOC
)

help_message_man_table_polap_disassemble=$(
  cat <<HEREDOC

  menu title
HEREDOC
)

help_message_man_table_benchmark=$(
  cat <<HEREDOC

  menu title
HEREDOC
)

help_message_kickstart=$(
  cat <<HEREDOC

  Log in a terminal of a Linux computer:
  mkdir -p all/polap/cflye1
  cd all/polap/cflye1
  git clone https://github.com/goshng/polap.git
  bash polap/src/polap-data-v2.sh install conda
  bash polap/src/polap-data-v2.sh setup conda
  bash polap/src/polap-data-v2.sh -y install all
  bash polap/src/polap-data-v2.sh setup polap
  bash polap/src/polap-data-v2.sh setup pmat
  Log out and back into the terminal:
  source ~/miniconda3/bin/activate
  cd all/polap/cflye1
  p -y mkdir-all
  p -y benchmark <outdir>
HEREDOC
)

help_message_setup_csv=$(
  cat <<HEREDOC

  Create a csv config file.
  csv: a default csv config file.
HEREDOC
)

help_message_test=$(
  cat <<HEREDOC

  Test the code.
HEREDOC
)

_v1_help_message_get=$(
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

help_message_benchmark=$(
  cat <<HEREDOC

  The main batch command does the followings.
  - get the data
  - GetOrganelle assembly
  - msbwt short-read polishing
  - fetch reference from NCBI
  - ptGAUL assembly
  - ptGAUL's polishing
  - NextDenovo long-read data error correction
  - Oatk assembly with the error-corrected long-read data
  - TIPPo assembly with the error-corrected long-read data
  - PMAT assembly with the error-corrected long-read data
  - Polap subsampling-based assembly

  See Also:
  benchmark-command
HEREDOC
)

help_message_benchmark_command=$(
  cat <<HEREDOC

  The main batch command does the following polap subcommands:
  data-long <outdir> [inum]
  data-short <outdir> [inum]
  run getorganelle <outdir> [inum]
  run msbwt <outdir> [inum]
  download ptdna <outdir> [inum]
  run ptgaul <outdir> [inum]
  run extract-ptdna-ptgau <outdir> [inum]
  run estimate-genomesize <outdir> [inum]
  run nextdenovo-polish <outdir> [inum]
  run oatk <outdir> [inum]
  run tippo <outdir> [inum]
  run pmat <outdir> [inum]
  run polap-disassemble <outdir> <inum> default
  run polap-disassemble <outdir> <inum> simple
  run polap-disassemble-check <outdir> <inum> default
  run polap-disassemble-check <outdir> <inum> simple
  run polap-disassemble-compare <outdir> <inum>
HEREDOC
)

help_message_archive=$(
  cat <<HEREDOC

  Archive the result.
  <inum>: -1 is default, meaning archiving it with the file name of -a.tar.gz
  <inum>: 0, meaning archiving it with the file name of -a-0.tar.gz
HEREDOC
)

############################################################
# CSV setting for each analysis
#
declare -A _taxon
declare -A _long
declare -A _short
declare -A _host
declare -A _random_seed
declare -A _downsample
declare -A _memory
declare -A _compare_p
declare -A _compare_n
declare -A _compare_r
declare -A _disassemble_alpha
declare -A _disassemble_delta
declare -A _ptgaul_genomesize
declare -A _bench_pmat
declare -A _bench_tippo
declare -A _bench_oatk
declare -A _dummy
declare -A _status

set +u

# Read the config files
read-a-tsv-file-into-associative-arrays() {
  # Define input TSV file
  if [ -z "${csv_file+x}" ]; then
    csv_file="${PWD}/${_polap_data_csv}"
  fi
  if [[ ! -s "${csv_file}" ]]; then
    csv_file="${_POLAPLIB_DIR}/${_polap_data_csv}"
  fi

  # Read the TSV file (skip header)
  #
  # v0.4.3.7.5
  # while IFS=$',' read -r species taxon folder long short host ptgaul_genomesize compare_n compare_p compare_r disassemble_alpha disassemble_delta disassemble_beta1 disassemble_beta2 disassemble_type random_seed ssh memory downsample inum table1 table2 mainfigure dummy status; do
  #
  # v0.4.3.7.5 2025-05-24 14:45
  # while IFS=$',' read -r species taxon folder long short host ssh ptgaul_genomesize disassemble_type inum disassemble_index random_seed downsample memory compare_p compare_n compare_r disassemble_alpha disassemble_delta disassemble_a disassemble_b table1 table2 mainfigure dummy status; do
  while IFS=$',' read -r species taxon long short host random_seed downsample memory compare_p compare_n compare_r disassemble_alpha disassemble_delta ptgaul_genomesize bench_pmat bench_tippo bench_oatk dummy status; do
    # Skip header line
    [[ "$species" == "species" ]] && continue
    [[ "$species" == \#* ]] && continue

    # Store in associative arrays
    if [[ -z "${species:-}" ]]; then
      continue
    fi
    _taxon["$species"]="$taxon"
    _long["$species"]="$long"
    _short["$species"]="$short"
    _host["$species"]="$host"
    _random_seed["$species"]="$random_seed"
    _downsample["$species"]="$downsample"
    _memory["$species"]="$memory"
    _compare_p["$species"]="$compare_p"
    _compare_n["$species"]="$compare_n"
    _compare_r["$species"]="$compare_r"
    _disassemble_alpha["$species"]="$disassemble_alpha"
    _disassemble_delta["$species"]="$disassemble_delta"
    _ptgaul_genomesize["$species"]="$ptgaul_genomesize"
    _bench_pmat["$species"]="$bench_pmat"
    _bench_tippo["$species"]="$bench_tippo"
    _bench_oatk["$species"]="$bench_oatk"
    _dummy["$species"]="$dummy"
    _status["$species"]="$status"
  done <"$csv_file"

  # Create Sall with all species folder names
  mapfile -t Sall < <(
    for key in "${!_long[@]}"; do
      echo "${key%%-*}"
    done | sort -u
  )
}

read-a-tsv-file-into-associative-arrays

# Create all keys
keys_array=($(for key in "${!_long[@]}"; do echo "$key"; done | sort))
Skeys=("${keys_array[@]}")

# Command arguments
_arg1=${1:-arg1}
# Check if the species folder is provided
if [[ "${_arg1}" == "arg1" ]]; then
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
_arg9=${9:-arg9}
_arg10=${10:-arg10}

################################################################################
# Part of genus_species
#
test_genus_species_for() {
  local _brg_outdir="${1}"
  local _brg_inum="${2:-0}"
  local key="${_brg_outdir}-${_brg_inum}"

  local long_sra="${_long["$key"]}"
  local short_sra="${_short["$key"]}"

  echo "Key: $key"
  if [[ "${_POLAP_RELEASE}" == "1" ]]; then
    echo "  long_sra: ${long_sra}"
    echo "  short_sra: ${short_sra}"
  fi
}

test_genus_species() {
  local _brg_outdir="${1:-all}"
  local _brg_inum="${2:-0}"

  if [[ "${_brg_outdir}" == "all" ]]; then
    for _v1 in "${Sall[@]}"; do
      test_genus_species_for "${_v1}" "${_brg_inum}"
    done
  elif [[ "${_brg_outdir}" == "each" ]]; then
    for key in "${Skeys[@]}"; do
      _brg_outdir="${key%-*}"
      _brg_inum="${key##*-}"
      test_genus_species_for "${_brg_outdir}" "${_brg_inum}"
    done
  else
    test_genus_species_for "$@"
  fi
}

##### INSERT_FUNCTION_HERE #####
benchmark-command_genus_species_for() {
  local _brg_outdir="$1"
  local _brg_inum="${2:-0}"
  local _brg_inum_0=0

  local full_name="${FUNCNAME[0]}"
  local middle_part="${full_name#man-}"
  local run_title="${middle_part%%_*}"

  local _brg_ref="${3:-off}"
  local _brg_confirm="${4:-off}"
  local _brg_redo="${5:-off}"
  local _brg_getorganelle="${6:-off}"
  local _brg_random="${7:-off|on}"

  local target_index="${_brg_outdir}-${_brg_inum}"
  local species_name="$(echo ${_brg_outdir} | sed 's/_/ /')"
  local long_sra="${_long["$target_index"]}"
  local short_sra="${_short["$target_index"]}"
  local random_seed="${_random_seed["$target_index"]}"
  local _brg_outdir_t="${_brg_outdir}/${opt_t_arg}"
  local _brg_outdir_i="${_brg_outdir_t}/${_brg_inum}"
  local _brg_outdir_0="${_brg_outdir_t}/0"

  if [[ "${opt_y_flag}" == false ]]; then
    echo "We will execute benchmarking using GetOrganelle, ptGAUL, PMAT, TIPPo, Oatk on ${_brg_outdir}/${opt_t_arg}/0 ..."
    read -p "Do you want to execute assemble on ${_brg_outdir}/${opt_t_arg}/${_brg_inum}? (yes/no): " confirm
  else
    confirm="yes"
  fi

  if [[ "${confirm}" != "yes" ]]; then
    echo "assemble processing with benchmarking analysis is canceled."
    return
  fi

  _log_echo "START: polap assembling using cflye with banchmarking analysis"

  # Prepare the input data
  polap-analysis-data_genus_species "${_brg_outdir}"

  # main
  mkdir -p "${_brg_outdir_i}"

  # Execute summary-data
  local _run_title="summary-data"
  local _brg_rundir="${_brg_outdir_0}/${_run_title}"
  if [[ -s "${_brg_rundir}"/l.fq.stats ]]; then
    _log_echo "Found: the data summary"
  else
    run-summary-data_genus_species "${_brg_outdir}" "${_brg_inum}"
    if [[ -s "${_brg_rundir}"/l.fq.stats ]]; then
      _log_echo "Success: the data summary"
    else
      _log_echo "Fail: the data summary"
      return 1
    fi
  fi

  # Execute GetOrganelle
  local _run_title="getorganelle"
  local _brg_rundir="${_brg_outdir_0}/${_run_title}"
  if [[ -d "${_brg_rundir}" ]] &&
    find "${_brg_rundir}" -maxdepth 1 -type f -name 'embplant_pt.*.gfa' -size +0c | grep -q .; then
    _log_echo "Found: GetOrganelle assembled ptDNA"
  else
    run-getorganelle_genus_species "${_brg_outdir}" "${_brg_inum}"
    if find "${_brg_rundir}" -maxdepth 1 -type f -name 'embplant_pt.*.gfa' -size +0c | grep -q .; then
      _log_echo "Success: GetOrganelle assembled ptDNA"
    else
      _log_echo "Fail: GetOrganelle assembled ptDNA"
      return 1
    fi
  fi

  # Execute FMLRC's msbwt preparation
  local _run_title="msbwt"
  local _brg_rundir="${_brg_outdir_0}/${_run_title}"
  if [[ -s "${_brg_rundir}/comp_msbwt.npy" ]]; then
    _log_echo "Found: FMLRC msbwt"
  else
    run-msbwt_genus_species "${_brg_outdir}" "${_brg_inum}"
    if [[ -s "${_brg_rundir}/comp_msbwt.npy" ]]; then
      _log_echo "Success: FMLRC msbwt"
    else
      _log_echo "Fail: FMLRC msbwt"
      return 1
    fi
  fi

  # Execute NCBI edirect download ptDNA
  local _run_title="ncbi-ptdna"
  local _brg_rundir="${_brg_outdir_0}/${_run_title}"
  if [[ -s "${_brg_rundir}/ptdna-reference.fa" ]]; then
    _log_echo "Found: reference ptDNA"
  else
    download-ptdna_genus_species "${_brg_outdir}"
    if [[ -s "${_brg_rundir}/ptdna-reference.fa" ]]; then
      _log_echo "Success: reference ptDNA"
    else
      _log_echo "Fail: reference ptDNA"
      return 1
    fi
  fi

  # Execute ptGAUL
  local _run_title="ptgaul"
  local _brg_rundir="${_brg_outdir_0}/${_run_title}"
  if [[ -s "${_brg_rundir}/flye_cpONT/assembly_graph.gfa" ]]; then
    _log_echo "Found: ptGAUL assembly"
  else
    run-ptgaul_genus_species "${_brg_outdir}" "${_brg_inum}"
    if [[ -s "${_brg_rundir}/flye_cpONT/assembly_graph.gfa" ]]; then
      _log_echo "Success: ptGAUL assembly"
    else
      _log_echo "Fail: ptGAUL assembly"
      return 1
    fi
  fi

  # Extract ptDNA from ptGAUL
  local _run_title="ptgaul"
  local _brg_rundir="${_brg_outdir_0}/${_run_title}"
  if [[ -s "${_brg_outdir_0}/ptdna-ptgaul.fa" ]]; then
    _log_echo "Found: ptGAUL polished genome"
  else
    run-extract-ptdna-ptgaul_genus_species "${_brg_outdir}" "${_brg_inum_0}"
    if [[ -s "${_brg_outdir_0}/ptdna-ptgaul.fa" ]]; then
      _log_echo "Success: ptGAUL polished genome"
    else
      _log_echo "Fail: ptGAUL polished genome"
      _log_echo "${help_message_extract_ptgaul_ptdna}"
      return 1
    fi
  fi

  # run-estimate-genomesize
  local _run_title="estimate-genomesize"
  local _brg_rundir="${_brg_outdir_0}/${_run_title}"
  if [[ -s "${_brg_rundir}/short_expected_genome_size.txt" ]]; then
    _log_echo "Found: the genome size estimate"
  else
    run-estimate-genomesize_genus_species "${_brg_outdir}" "${_brg_inum_0}"
    if [[ -s "${_brg_rundir}/short_expected_genome_size.txt" ]]; then
      _log_echo "Success: the genome size estimate"
    else
      _log_echo "Fail: the genome size estimate"
      return 1
    fi
  fi

  # run-nextdenovo-polish
  local _run_title="nextdenovo-polish"
  local _brg_rundir="${_brg_outdir_0}/${_run_title}"
  if [[ -s "${_brg_outdir_0}/cns.fa" ]]; then
    _log_echo "Found: the long-read NextDenovo polished data"
  else
    run-nextdenovo-polish_genus_species "${_brg_outdir}" "${_brg_inum_0}"
    if [[ -s "${_brg_outdir_0}/cns.fa" ]]; then
      _log_echo "Success: the long-read NextDenovo polished data"
    else
      _log_echo "Fail: the long-read NextDenovo polished data"
      return 1
    fi
  fi

  # run-oatk
  local _run_title="oatk-nextdenovo"
  local _brg_rundir="${_brg_outdir_0}/${_run_title}"
  if [[ -s "${_brg_rundir}/oatk-nextdenovo-30.utg.gfa" ]]; then
    _log_echo "Found: the oatk result"
  else
    run-oatk_genus_species "${_brg_outdir}" "${_brg_inum_0}"
    if [[ -s "${_brg_rundir}/oatk-nextdenovo-30.utg.gfa" ]]; then
      _log_echo "Success: the oatk result"
    else
      _log_echo "Fail: the oatk result"
      return 1
    fi
  fi

  # run-tippo
  local _run_title="tippo-nextdenovo"
  local _brg_rundir="${_brg_outdir_0}/${_run_title}"
  if [[ -s "${_brg_rundir}/ont/log_cns.fa.tiara.out" ]]; then
    _log_echo "Found: the tippo result"
  else
    run-tippo_genus_species "${_brg_outdir}" "${_brg_inum_0}"
    if [[ -s "${_brg_rundir}/ont/log_cns.fa.tiara.out" ]]; then
      _log_echo "Success: the tippo result"
    else
      _log_echo "Fail: the tippo result"
      return 1
    fi
  fi

  # run-pmat
  local _run_title="pmat-nextdenovo"
  local _brg_rundir="${_brg_outdir_0}/${_run_title}"
  # NOTE: PMAT run first try to run withouth time limit.
  # If it fails for some reason, then it crashes. If the run time is over 12 hours,
  # then we would use PMATContigGraph.txt to skip the run and just report the run time.
  local _check_file_for_finish="${_brg_rundir}/0.1/assembly_result/PMATContigGraph.txt"
  local _check_file_for_finish="${_brg_rundir}/0.1/PMAT_mt_blastn.txt"
  if [[ -s "${_check_file_for_finish}" ]]; then
    _log_echo "Found: the pmat result"
  else
    run-pmat_genus_species "${_brg_outdir}" "${_brg_inum_0}"
    if [[ -s "${_brg_rundir}/0.1/PMAT_mt_blastn.txt" ]]; then
      _log_echo "Success: the pmat result"
    else
      _log_echo "Fail: the pmat result"
      return 1
    fi
  fi

  # run-polap
  local _run_title="polap-disassemble"
  local _brg_rundir="${_brg_outdir_i}/${_run_title}"
  if [[ -s "${_brg_rundir}/${_brg_inum}/disassemble/infer-1/pt.subsample-polishing.1.fa" ]]; then
    _log_echo "Found: infer case"
  else
    run-polap-disassemble_genus_species "${_brg_outdir}" "${_brg_inum}" default "${_brg_random}"
    run-polap-disassemble_genus_species "${_brg_outdir}" "${_brg_inum}" simple "${_brg_random}"
    if [[ -s "${_brg_rundir}/${_brg_inum}/disassemble/infer-1/pt.subsample-polishing.1.fa" ]]; then
      _log_echo "Success: infer case"
    else
      _log_echo "Fail: infer case"
      return 1
    fi
  fi

  if [[ -s "${_brg_outdir_0}/ptdna-ptgaul.fa" ]]; then
    if [[ "${_brg_outdir_0}" != "${_brg_outdir_i}" ]]; then
      cp -p "${_brg_outdir_0}/ptdna-ptgaul.fa" "${_brg_outdir_i}"
    fi
    local _run_title="polap-disassemble"
    local _brg_rundir="${_brg_outdir_i}/${_run_title}"
    if [[ -s "${_brg_rundir}/${_brg_inum}/disassemble/infer-1/pt.subsample-polishing.reference.aligned.1.fa" ]]; then
      _log_echo "Found: check case"
    else
      run-polap-disassemble-check_genus_species "${_brg_outdir}" "${_brg_inum}" default "${_brg_random}"
      run-polap-disassemble-check_genus_species "${_brg_outdir}" "${_brg_inum}" simple "${_brg_random}"
      run-polap-disassemble-compare_genus_species "${_brg_outdir}" "${_brg_inum}"
      if [[ -s "${_brg_rundir}/${_brg_inum}/disassemble/infer-1/pt.subsample-polishing.reference.aligned.1.fa" ]]; then
        _log_echo "Success: check case"
      else
        _log_echo "Fail: check case"
        return 1
      fi
    fi
  fi

  _log_echo "END: polap assembling using cflye with banchmarking analysis"
}

benchmark-command_genus_species() {
  local _brg_outdir="${1:-all}"
  local _brg_inum="${2:-0}"

  if [[ "${_brg_outdir}" == "all" ]]; then
    for _v1 in "${Sall[@]}"; do
      benchmark-command_genus_species_for "${_v1}" "${@:2}"
    done
  elif [[ "${_brg_outdir}" == "each" ]]; then
    for _v1 in "${Sall[@]}"; do
      benchmark-command_genus_species_for "${_v1}" "${@:2}"
    done
  else
    benchmark-command_genus_species_for "$@"
  fi
}

benchmark_genus_species_for() {
  local _brg_outdir="${1}"

  _brg_outdir="${_brg_outdir%/}"
  # batch-cmd_genus_species "${_brg_outdir}" 0
  for i in 0 1 2 4; do
    benchmark-command_genus_species "${_brg_outdir}" $i
  done
}

benchmark_genus_species() {
  local args=("$@")

  for item in "$@"; do
    benchmark_genus_species_for "${item}"
  done
}

################################################################################
# Figures and Tables
#
man-figure-delta_genus_species() {
  local _brg_outdir="${1:-Eucalyptus_pauciflora}"
  local _brg_inum="${2:-0}"
  local _brg_t_dir="${3:-"${_brg_default_target_dir}"}"

  # Folders
  local _brg_outdir_t="${_brg_outdir}/${opt_t_arg}"
  local _brg_outdir_i="${_brg_outdir_t}/${_brg_inum}"

  local _suppfigure_file="delta.pdf"

  # for i in {21..27}; do p get Eucalyptus_pauciflora $i add off; done

  local -A number2delta
  number2delta["21"]="0.25"
  number2delta["22"]="0.50"
  number2delta["23"]="0.75"
  number2delta["24"]="1.25"
  number2delta["25"]="1.50"
  number2delta["26"]="1.75"
  number2delta["27"]="2.00"

  _polap_lib_conda-ensure_conda_env polap

  rm ?.??.tsv
  for i in {21..27}; do
    csvtk -t cut -f index,alpha "${_brg_outdir_i}"/polap-disassemble/$i/disassemble/infer-1/1/summary1.txt >${number2delta[$i]}.tsv
  done

  Rscript ${_POLAPLIB_DIR}/run-polap-r-data-v2-alpha0.R ?.??.tsv -l delta -o "${_suppfigure_file}"

  conda deactivate

  if [[ -d "${_brg_t_dir}" ]]; then
    cp -p ${_suppfigure_file} "${_brg_t_dir}"
    echo copy ${_suppfigure_file} to ${_brg_t_dir}
  else
    echo "Error: no such target dir: ${_brg_t_dir}"
  fi
}

man-figure-alpha_genus_species() {
  local _brg_outdir="${1:-Eucalyptus_pauciflora}"
  local _brg_inum="${2:-0}"
  local _brg_t_dir="${3:-"${_brg_default_target_dir}"}"

  # Folders
  local _brg_outdir_t="${_brg_outdir}/${opt_t_arg}"
  local _brg_outdir_i="${_brg_outdir_t}/${_brg_inum}"

  local _suppfigure_file="alpha0.pdf"

  local -A number2alpha0
  number2alpha0["11"]="0.00"
  number2alpha0["12"]="1.00"
  number2alpha0["13"]="2.00"
  number2alpha0["14"]="3.00"
  number2alpha0["15"]="4.00"
  number2alpha0["16"]="5.00"
  number2alpha0["17"]="6.00"
  number2alpha0["18"]="7.00"

  _polap_lib_conda-ensure_conda_env polap

  rm ?.??.tsv
  for i in {11..18}; do
    csvtk -t cut -f index,alpha "${_brg_outdir_i}"/polap-disassemble/$i/disassemble/infer-1/1/summary1.txt >${number2alpha0[$i]}.tsv
  done

  Rscript ${_POLAPLIB_DIR}/run-polap-r-data-v2-alpha0.R ?.??.tsv -l alpha0 -o "${_suppfigure_file}"

  conda deactivate

  if [[ -d "${_brg_t_dir}" ]]; then
    cp -p ${_suppfigure_file} "${_brg_t_dir}"
    echo copy ${_suppfigure_file} to ${_brg_t_dir}
  else
    echo "Error: no such target dir: ${_brg_t_dir}"
  fi
}

man-figure-polap_genus_species_for() {
  local _brg_outdir="${1}"
  local _brg_inum="${2:-0}"

}

man-figure-polap_genus_species() {
  local _brg_outdir="${1:-all}"
  local _brg_inum="${2:-0}"

  if [[ "${_brg_outdir}" == "all" ]]; then
    for _v1 in "${Sall[@]}"; do
      man-figure-polap_genus_species_for "${_v1}" "${@:2}"
    done
  elif [[ "${_brg_outdir}" == "each" ]]; then
    for _v1 in "${Sall[@]}"; do
      man-figure-polap_genus_species_for "${_v1}" "${@:2}"
    done
  else
    man-figure-polap_genus_species_for "$@"
  fi
}

man-figure-sheet-pmat_genus_species_for() {
  local _brg_outdir="${1}"
  local _brg_inum="${2:-0}"
  local _brg_bandage="${3:-off}"
  local _brg_csv="${4:-sheet_pmat.csv}"
  local _run_title="pmat-nextdenovo"

  local _base_figure="."

  local _key="${_brg_outdir}-${_brg_inum}"

  local _species=${_taxon[$_key]}
  local _species="${_species//_/ }"
  local _genus=${_species%% *}
  local _order=$(grep "${_genus}" ${_POLAPLIB_DIR}/taxonomy_output.tsv | cut -f 6 | head -1)
  local _label_base="${_brg_outdir/_/-}"

  # Result folders
  local _brg_outdir_i="${_brg_outdir}/${opt_t_arg}/${_brg_inum}"
  local _brg_rundir="${_brg_outdir_i}/${_run_title}"

  local _brg_fc="0.1,1.0"
  local fc_list

  _polap_lib_array-csv_to_array "${_brg_fc}" fc_list

  for _fc in "${fc_list[@]}"; do

    # local formatted_fc=$(printf "%02d" "${_fc}")
    local _gfa_infer="${_brg_rundir}/${_fc}/gfa_result/PMAT_pt_master.gfa"
    local _png_infer="${_brg_rundir}/${_fc}/gfa_result/PMAT_pt_master.png"

    if [[ -s "${_gfa_infer}" ]]; then
      # echo "gfa file: ${_gfa_infer}"
      if [[ "${_brg_bandage}" == "on" ]]; then
        ${_polap_cmd} bandage png \
          ${_gfa_infer} \
          ${_png_infer}
      fi
      images+=("figures/${_png_infer}")
      captions+=("${i}, ${_disjointig_min_coverage}x, ${_memory_gb_cflye} GB")
      ((k++))

      if [[ "$_POLAP_DEBUG" == "0" ]]; then
        printf "%s,%s,%s,%s\n" "${_run_title}" "${_species}" "${_fc}" "${_base_figure}/${_png_infer}" >>"${_brg_csv}"
      else
        printf "%s,%s,%s,%s\n" "${_run_title}" "${_species}" "${_fc}" "${_base_figure}/${_png_infer}"
      fi
    else
      printf "%s,%s,%s,%s\n" "${_run_title}" "${_species}" "${_fc}" "${_base_figure}/empty.png" >>"${_brg_csv}"
      echo "no such file: ${_gfa_infer}"
    fi
  done
}

man-figure-sheet-pmat_genus_species() {
  local _brg_outdir="${1:-all}"
  local _brg_inum="${2:-0}"
  local _brg_bandage="${3:-off}"
  local _brg_csv="${4:-sheet_pmat.csv}"

  if [[ "${_brg_outdir}" == "all" ]]; then
    for _v1 in "${Sall[@]}"; do
      man-figure-sheet-pmat_genus_species_for "${_v1}" "${@:2}"
    done
  elif [[ "${_brg_outdir}" == "each" ]]; then
    for _v1 in "${Sall[@]}"; do
      man-figure-sheet-pmat_genus_species_for "${_v1}" "${@:2}"
    done
  else
    man-figure-sheet-pmat_genus_species_for "$@"
  fi
}

man-figure-sheet-tippo_genus_species_for() {
  local _brg_outdir="${1}"
  local _brg_inum="${2:-0}"
  local _brg_bandage="${3:-off}"
  local _brg_csv="${4:-sheet_tippo.csv}"
  local _run_title="tippo-nextdenovo"

  local _base_figure="."

  local _key="${_brg_outdir}-${_brg_inum}"

  #
  local _species=${_taxon[$_key]}
  local _species="${_species//_/ }"
  local _genus=${_species%% *}
  local _order=$(grep "${_genus}" ${_POLAPLIB_DIR}/taxonomy_output.tsv | cut -f 6 | head -1)
  local _label_base="${_brg_outdir/_/-}"

  # Result folders
  local _brg_outdir_i="${_brg_outdir}/${opt_t_arg}/${_brg_inum}"
  local _brg_rundir="${_brg_outdir_i}/${_run_title}"

  local _brg_fc="hifi,clr,ont,onthq"
  local fc_list

  _polap_lib_array-csv_to_array "${_brg_fc}" fc_list

  for _fc in "${fc_list[@]}"; do
    # local formatted_fc=$(printf "%02d" "${_fc}")
    local _gfa_infer="${_brg_rundir}/${_fc}/cns.fa.chloroplast.fasta.filter.800.round1.fasta.chloroplast.flye/assembly_graph.gfa"
    local _png_infer="${_brg_rundir}/${_fc}/cns.fa.chloroplast.fasta.filter.800.round1.fasta.chloroplast.flye/assembly_graph.png"

    if [[ -s "${_gfa_infer}" ]]; then
      # echo "gfa file: ${_gfa_infer}"
      if [[ "${_brg_bandage}" == "on" ]]; then
        ${_polap_cmd} bandage png \
          ${_gfa_infer} \
          ${_png_infer}
      fi
      images+=("figures/${_png_infer}")
      captions+=("${i}, ${_disjointig_min_coverage}x, ${_memory_gb_cflye} GB")
      ((k++))

      if [[ "$_POLAP_DEBUG" == "0" ]]; then
        printf "%s,%s,%s,%s\n" "${_run_title}" "${_species}" "${_fc}" "${_base_figure}/${_png_infer}" >>"${_brg_csv}"
      else
        printf "%s,%s,%s,%s\n" "${_run_title}" "${_species}" "${_fc}" "${_base_figure}/${_png_infer}"
      fi
    else
      printf "%s,%s,%s,%s\n" "${_run_title}" "${_species}" "${_fc}" "${_base_figure}/empty.png" >>"${_brg_csv}"
      echo "no such file: ${_gfa_infer}"
    fi
  done

}

man-figure-sheet-tippo_genus_species() {
  local _brg_outdir="${1:-all}"
  local _brg_inum="${2:-0}"
  local _brg_bandage="${3:-off}"
  local _brg_csv="${4:-sheet_tippo.csv}"
  local _run_title="tippo-nextdenovo"

  if [[ "${_brg_outdir}" == "all" ]]; then
    for _v1 in "${Sall[@]}"; do
      man-figure-sheet-tippo_genus_species_for "${_v1}" "${@:2}"
    done
  elif [[ "${_brg_outdir}" == "each" ]]; then
    for _v1 in "${Sall[@]}"; do
      man-figure-sheet-tippo_genus_species_for "${_v1}" "${@:2}"
    done
  else
    man-figure-sheet-tippo_genus_species_for "$@"
  fi
}

man-figure-sheet-oatk_genus_species_for() {
  local _brg_outdir="${1}"
  local _brg_inum="${2:-0}"
  local _brg_bandage="${3:-off}"
  local _brg_csv="${4:-sheet_oatk.csv}"
  local _run_title="oatk-nextdenovo"

  local _base_figure="."
  local _key="${_brg_outdir}-${_brg_inum}"

  local _species=${_taxon[$_key]}
  local _species="${_species//_/ }"
  local _genus=${_species%% *}
  local _order=$(grep "${_genus}" ${_POLAPLIB_DIR}/taxonomy_output.tsv | cut -f 6 | head -1)
  local _label_base="${_brg_outdir/_/-}"

  # Result folders
  local _brg_outdir_i="${_brg_outdir}/${opt_t_arg}/${_brg_inum}"
  local _brg_rundir="${_brg_outdir_i}/oatk-nextdenovo"

  local _brg_fc="1,10,20,30"
  local fc_list

  _polap_lib_array-csv_to_array "${_brg_fc}" fc_list

  for _fc in "${fc_list[@]}"; do

    local formatted_fc=$(printf "%02d" "${_fc}")
    local _gfa_infer="${_brg_rundir}/oatk-nextdenovo-${formatted_fc}.pltd.gfa"
    local _png_infer="${_brg_rundir}/oatk-nextdenovo-${formatted_fc}.pltd.png"

    if [[ -s "${_gfa_infer}" ]]; then
      # echo "gfa file: ${_gfa_infer}"
      if [[ "${_brg_bandage}" == "on" ]]; then
        ${_polap_cmd} bandage png \
          ${_gfa_infer} \
          ${_png_infer}
      fi
      images+=("figures/${_png_infer}")
      captions+=("${i}, ${_disjointig_min_coverage}x, ${_memory_gb_cflye} GB")
      ((k++))

      if [[ "$_POLAP_DEBUG" == "0" ]]; then
        printf "%s,%s,%s,%s\n" "${_run_title}" "${_species}" "${_fc}" "${_base_figure}/${_png_infer}" >>"${_brg_csv}"
      else
        printf "%s,%s,%s,%s\n" "${_run_title}" "${_species}" "${_fc}" "${_base_figure}/${_png_infer}"
      fi
    else
      printf "%s,%s,%s,%s\n" "${_run_title}" "${_species}" "${_fc}" "${_base_figure}/empty.png" >>"${_brg_csv}"
      echo "no such file: ${_gfa_infer}"
    fi

  done

  # debug local variables
  # for debugging: Inline local printing local var
  # if [[ "$_POLAP_DEBUG" == "1" ]]; then
  #   while IFS= read -r line; do
  #     if [[ $line =~ ^declare\ --\ ([^=]+)= ]]; then
  #       var="${BASH_REMATCH[1]}"
  #       printf "%s=%q\n" "$var" "${!var}"
  #     fi
  #   done < <(local -p 2>/dev/null)
  #   return
  # fi
}

man-figure-sheet-oatk_genus_species() {
  local _brg_outdir="${1:-all}"
  local _brg_inum="${2:-0}"
  local _brg_bandage="${3:-off}"
  local _brg_csv="${4:-sheet_oatk.csv}"

  if [[ "${_brg_outdir}" == "all" ]]; then
    for _v1 in "${Sall[@]}"; do
      man-figure-sheet-oatk_genus_species_for "${_v1}" "${@:2}"
    done
  elif [[ "${_brg_outdir}" == "each" ]]; then
    for _v1 in "${Sall[@]}"; do
      man-figure-sheet-oatk_genus_species_for "${_v1}" "${@:2}"
    done
  else
    man-figure-sheet-oatk_genus_species_for "$@"
  fi
}

# All figures in a sheet from polap disassemble
man-figure-sheet-polap_genus_species_for() {
  local _brg_outdir="${1}"
  local _brg_inum="${2:-0}"
  local _brg_bandage="${3:-off}"
  local _brg_csv="${4:-sheet_polap.csv}"
  local _run_title="polap"

  local _base_figure="."

  local _key="${_brg_outdir}-${_brg_inum}"
  local _csv="sheet-polap.csv"

  local _species=${_taxon[$_key]}
  local _species="${_species//_/ }"
  local _genus=${_species%% *}
  local _order=$(grep "${_genus}" ${_POLAPLIB_DIR}/taxonomy_output.tsv | cut -f 6 | head -1)
  local _label_base="${_brg_outdir/_/-}"

  local extracted_n="${_compare_n["$_key"]}"

  local _brg_d_index="infer-1"
  local _brg_stage="1"
  local _brg_bandage="yes"
  local _v1_inum="${_brg_outdir}/${opt_t_arg}/${_brg_inum}/polap-disassemble/0"
  local k=0
  for ((i = 0; i < extracted_n; i++)); do
    local _gfa_infer="${_v1_inum}/disassemble/${_brg_d_index}/${_brg_stage}/${i}/30-contigger/graph_final.gfa"
    local _png_infer="${_v1_inum}/disassemble/${_brg_d_index}/${_brg_stage}/${i}/30-contigger/graph_final.png"
    local _timing_cflye="${_v1_inum}/disassemble/${_brg_d_index}/${_brg_stage}/${i}/timing-cflye.txt"
    read -r _memory_gb_cflye _total_hours_cflye < <(_polap_lib_timing-parse-timing "${_timing_cflye}")
    local _disjointig_min_coverage=$(grep "Command being timed" "${_timing_cflye}" | awk -F '--disjointig-min-coverage ' '{print $2}' | awk '{print $1}')
    if [[ -s "${_gfa_infer}" ]]; then
      # echo "gfa file: ${_gfa_infer}"
      if [[ "${_brg_bandage}" == "yes" ]]; then
        ${_polap_cmd} bandage png \
          ${_gfa_infer} \
          ${_png_infer}
      fi
      # printf "| ![polap %s](figures/%s){ width=%s%% } " "${i}" "${_png_infer}" "${width}" >>"${_supptable_md}"
      images+=("figures/${_png_infer}")
      captions+=("${i}, ${_disjointig_min_coverage}x, ${_memory_gb_cflye} GB")
      ((k++))

      printf "polap,%s,%d,%s\n" "${_species}" "$i" "figures/${_png_infer}" >>"${_brg_csv}"

      if [[ "$_POLAP_DEBUG" == "0" ]]; then
        printf "%s,%s,%s,%s\n" "${_run_title}" "${_species}" "${i}" "${_base_figure}/${_png_infer}" >>"${_brg_csv}"
      else
        printf "%s,%s,%s,%s\n" "${_run_title}" "${_species}" "${i}" "${_base_figure}/${_png_infer}"
      fi
    else
      printf "%s,%s,%s,%s\n" "${_run_title}" "${_species}" "${i}" "${_base_figure}/empty.png" >>"${_brg_csv}"
      echo "no such file: ${_gfa_infer}"
    fi
  done

  # debug local variables
  # for debugging: Inline local printing local var
  if [[ "$_POLAP_DEBUG" == "1" ]]; then
    while IFS= read -r line; do
      if [[ $line =~ ^declare\ --\ ([^=]+)= ]]; then
        var="${BASH_REMATCH[1]}"
        printf "%s=%q\n" "$var" "${!var}"
      fi
    done < <(local -p 2>/dev/null)
    return
  fi
}

man-figure-sheet-polap_genus_species() {
  local _brg_outdir="${1:-all}"
  local _brg_inum="${2:-0}"
  local _brg_bandage="${3:-off}"
  local _brg_csv="${4:-sheet_polap.csv}"

  # rm -f "${_brg_csv}"

  if [[ "${_brg_outdir}" == "all" ]]; then
    for _v1 in "${Sall[@]}"; do
      man-figure-sheet-polap_genus_species_for "${_v1}" "${@:2}"
    done
  elif [[ "${_brg_outdir}" == "each" ]]; then
    for _v1 in "${Sall[@]}"; do
      man-figure-sheet-polap_genus_species_for "${_v1}" "${@:2}"
    done
  else
    man-figure-sheet-polap_genus_species_for "$@"
  fi
}

man-figure-sheet_genus_species_for() {
  local _brg_outdir="${1}"
  local _brg_inum="${2:-0}"
  local _brg_bandage="${3:-off}"
  local _brg_csv="${4:-sheet_benchmark.csv}"

  local _base_figure="."

  local _key="${_brg_outdir}-${_brg_inum}"

  # Species name
  local _species=${_taxon[$_key]}
  local _species="${_species//_/ }"
  local _genus=${_species%% *}
  local _order=$(grep "${_genus}" ${_POLAPLIB_DIR}/taxonomy_output.tsv | cut -f 6 | head -1)

  # Result folders
  local _brg_outdir_i="${_brg_outdir}/${opt_t_arg}/${_brg_inum}"

  # GetOrganelle
  local _run_title="getorganelle"
  local _brg_rundir="${_brg_outdir_i}/${_run_title}"
  local _gfa_infer="${_brg_rundir}/embplant_pt.K115.complete.graph1.selected_graph.gfa"
  local _png_infer="${_brg_rundir}/embplant_pt.K115.complete.graph1.selected_graph.png"
  if [[ -s "${_gfa_infer}" ]]; then
    if [[ "${_brg_bandage}" == "on" ]]; then
      ${_polap_cmd} bandage png ${_gfa_infer} ${_png_infer}
    fi
    printf "%s,%s,%s,%s\n" "${_run_title}" "${_species}" "GetOrganelle" "${_base_figure}/${_png_infer}" >>"${_brg_csv}"
  else
    echo "no such file: ${_gfa_infer}"
  fi

  # ptGAUL
  local _run_title="ptgaul"
  local _brg_rundir="${_brg_outdir_i}/${_run_title}"
  local _gfa_infer="${_brg_rundir}/flye_cpONT/assembly_graph.gfa"
  local _png_infer="${_brg_rundir}/flye_cpONT/assembly_graph.png"
  if [[ -s "${_gfa_infer}" ]]; then
    if [[ "${_brg_bandage}" == "on" ]]; then
      ${_polap_cmd} bandage png ${_gfa_infer} ${_png_infer}
    fi
    printf "%s,%s,%s,%s\n" "${_run_title}" "${_species}" "ptGAUL" "${_base_figure}/${_png_infer}" >>"${_brg_csv}"
  else
    echo "no such file: ${_gfa_infer}"
  fi

  # Polap

  # PMAT
  local _run_title="pmat-nextdenovo"
  local _brg_rundir="${_brg_outdir_i}/${_run_title}"
  local _fc="${_bench_pmat[$_key]}"
  local _gfa_infer="${_brg_rundir}/${_fc}/gfa_result/PMAT_pt_master.gfa"
  local _png_infer="${_brg_rundir}/${_fc}/gfa_result/PMAT_pt_master.png"
  if [[ -s "${_gfa_infer}" ]]; then
    if [[ "${_brg_bandage}" == "on" ]]; then
      ${_polap_cmd} bandage png ${_gfa_infer} ${_png_infer}
    fi
    printf "%s,%s,%s,%s\n" "${_run_title}" "${_species}" "PMAT" "${_base_figure}/${_png_infer}" >>"${_brg_csv}"
  else
    printf "%s,%s,%s,%s\n" "${_run_title}" "${_species}" "PMAT" "${_base_figure}/empty.png" >>"${_brg_csv}"
  fi

  # TIPPo
  local _run_title="tippo-nextdenovo"
  local _brg_rundir="${_brg_outdir_i}/${_run_title}"
  local _fc="${_bench_tippo[$_key]}"
  local _gfa_infer="${_brg_rundir}/${_fc}/cns.fa.chloroplast.fasta.filter.800.round1.fasta.chloroplast.flye/assembly_graph.gfa"
  local _png_infer="${_brg_rundir}/${_fc}/cns.fa.chloroplast.fasta.filter.800.round1.fasta.chloroplast.flye/assembly_graph.png"
  if [[ -s "${_gfa_infer}" ]]; then
    if [[ "${_brg_bandage}" == "on" ]]; then
      ${_polap_cmd} bandage png ${_gfa_infer} ${_png_infer}
    fi
    printf "%s,%s,%s,%s\n" "${_run_title}" "${_species}" "TIPPo" "${_base_figure}/${_png_infer}" >>"${_brg_csv}"
  else
    printf "%s,%s,%s,%s\n" "${_run_title}" "${_species}" "TIPPo" "${_base_figure}/empty.png" >>"${_brg_csv}"
  fi

  # oatk
  local _run_title="oatk-nextdenovo"
  local _brg_rundir="${_brg_outdir_i}/${_run_title}"
  local _fc="${_bench_oatk[$_key]}"
  local formatted_fc=$(printf "%02d" "${_fc}")
  local _gfa_infer="${_brg_rundir}/oatk-nextdenovo-${formatted_fc}.pltd.gfa"
  local _png_infer="${_brg_rundir}/oatk-nextdenovo-${formatted_fc}.pltd.png"
  if [[ -s "${_gfa_infer}" ]]; then
    if [[ "${_brg_bandage}" == "on" ]]; then
      ${_polap_cmd} bandage png ${_gfa_infer} ${_png_infer}
    fi
    printf "%s,%s,%s,%s\n" "${_run_title}" "${_species}" "Oatk" "${_base_figure}/${_png_infer}" >>"${_brg_csv}"
  else
    printf "%s,%s,%s,%s\n" "${_run_title}" "${_species}" "Oatk" "${_base_figure}/empty.png" >>"${_brg_csv}"
  fi

}

# All gfa figures from polap and other tools
man-figure-sheet_genus_species() {
  local _brg_outdir="${1:-all}"
  local _brg_inum="${2:-0}"
  local _brg_bandage="${3:-off}"
  local _brg_csv="${4:-sheet_benchmark.csv}"

  if [[ "${_brg_outdir}" == "all" ]]; then
    for _v1 in "${Sall[@]}"; do
      man-figure-sheet_genus_species_for "${_v1}" "${@:2}"
    done
  elif [[ "${_brg_outdir}" == "each" ]]; then
    for _v1 in "${Sall[@]}"; do
      man-figure-sheet_genus_species_for "${_v1}" "${@:2}"
    done
  else
    man-figure-sheet_genus_species_for "$@"
  fi
}

man_genus_species() {
  local first_arg="$1"
  local remaining_args=("${@:2}")

  man-${first_arg}_genus_species "${remaining_args[@]}"
}

# copy figures from the analysis to the manuscript folder.
copy-figures_genus_species() {
  local _brg_t_dir_figures="${1:-"${_brg_default_target_dir}"}"
  rsync -qav --include='*/' --include='*.png' --exclude='*' ./ "${_brg_t_dir_figures}/"
  rsync -qav --include='*/' --include='*.pdf' --exclude='*' ./ "${_brg_t_dir_figures}/"
  # rsync -av --include='*/' --include='*.png' --exclude='*' ./ ../../manuscript/polap-v0.4/figures/
  # rsync -av --include='*/' --include='*.pdf' --exclude='*' ./ ../../manuscript/polap-v0.4/figures/
}

archive_genus_species_for() {
  local _brg_outdir="$1"
  local _brg_inum="${2:--1}"

  # if [[ -d "src" ]]; then
  # 	_polap_cmd="src/polap.sh"
  # else
  # 	_polap_cmd="polap"
  # fi

  rm -rf "${_brg_outdir}-a"
  rm -f "${_brg_outdir}-a.tar.gz"
  ${_polap_cmd} disassemble archive \
    --max-filesize 5M \
    -o ${_brg_outdir}
  if [[ "${_brg_inum}" != "-1" ]]; then
    mv "${_brg_outdir}-a.tar.gz" "${_brg_outdir}-a-${_brg_inum}.tar.gz"
    _log_echo "  creating ${_brg_outdir}-a-${_brg_inum}.tar.gz ..."
    echo "  creating ${_brg_outdir}-a-${_brg_inum}.tar.gz ..."
  else
    _log_echo "  creating ${_brg_outdir}-a.tar.gz ..."
    echo "  creating ${_brg_outdir}-a.tar.gz ..."
  fi

  rm -rf "${_brg_outdir}-a"
}

archive_genus_species() {
  local _brg_outdir="${1:-all}"
  local _brg_inum="${2:--1}"

  if [[ "${_brg_outdir}" == "all" ]]; then
    for _v1 in "${Sall[@]}"; do
      archive_genus_species_for "${_v1}" "${_brg_inum}"
    done
  else
    archive_genus_species_for "$@"
  fi
}

# kill nextdenovo jobs.
nextdenovo-polish-cleanup_genus_species() {
  local _brg_outdir="${1}"

  if [[ ! -d "${_brg_outdir}" ]]; then
    echo "Error: no such outdir: ${_brg_outdir}"
    return
  fi

  rm -rf "${_brg_outdir}"-nextdenovo
  rm -rf "${_brg_outdir}"-pmat
  pkill -f execute_nextdenovo_"${_brg_outdir}".sh
  pkill -f nextDenovo
  pkill -f minimap2-nd
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
#
# result=$(parse_params "params.txt")
# read -r I P N R <<< "$result"
#
# local _ipn=$(parse_params "${_params_txt}")
# local _I _P _N _R _A _B _M _D _Alpha _Memory
# read -r _I _P _N _R _A _B _M _D _Alpha _Memory <<<"$_ipn"
parse_params() {
  local file="$1" # Input file
  local I=-1
  local P=-1
  local N=-1
  local R=-1
  local A=-1
  local B=-1
  local M=-1
  local D=-1
  local Alpha=-1
  local Memory=-1 # Declare variables for the parameters

  # Read the file and extract the values
  while IFS=": " read -r key value; do
    case "$key" in
    "I") I="$value" ;;
    "P") P="$value" ;;
    "N") N="$value" ;;
    "R") R="$value" ;;
    "A") A="$value" ;;
    "B") B="$value" ;;
    "M") M="$value" ;;
    "D") D="$value" ;;
    "Alpha") Alpha="$value" ;;
    "Memory") Memory="$value" ;;
    esac
  done <"$file"

  # Print the variables (return as output)
  echo "$I $P $N $R $A $B $M $D $Alpha $Memory"
}

# functions
extract_and_replace_suffix() {
  local input_file="$1"
  local replacement_suffix="$2"
  local output_file="${3:-filtered_output.csv}"

  {
    # Print the header
    head -n 1 "$input_file"

    # Replace -0 with given suffix in first column only
    awk -F',' -v OFS=',' -v suffix="$replacement_suffix" '
      NR > 1 && $1 ~ /-0$/ {
        sub(/-0$/, suffix, $1)
        print
      }
    ' "$input_file"
  } >"$output_file"
}

############################################################
# Table for benchmarking
#
# a copy of maintable1 to focus on the memory and time used
# by GetOrganelle, ptGAUL, PMAT, and Polap's disassemble menu
# to construct plastid genomes.
#
# Table 1. Benchmarking by count of success cases
# Figure 1. ptDNA analyses in a figure and graphs
#
# Supplementary materials:
#
# Table S1. data table
# Table S2. banchmarking table comparing with other tools
# Table S3. table for polap disassemble 5 % - n=r=5
# Table S4. table for polap disassemble 10 % - n=r=5
# Table S5. table for polap disassemble 1 % - n=r=5
# Table S6. table for polap disassemble 10 % - n=50, r=10 (only Eucalyptus_pauciflora)
#
# Others:
# Table S. polap disassemble stage tables for all the datasets
# Figure 1. method description
# Figures S1 type. comparison of polap, pmat, tippo, and oatk gfa figures
# Figures S2 type. polap gfa figures
# Figures S3 type. pmat gfa figures
# Figures S4 type. tippo gfa figures
# Figures S5 type. oatk gfa figures
# Figure S. alpha
# Figure S. delta
#
man-table-data_genus_species() {
  local _brg_outdir="${1:-all}"
  local _brg_inum="${2:-0}"

  man-table-benchmark_genus_species "${_brg_outdir}" "${_brg_inum}" data
}

man-table-benchmark_genus_species_header() {
  local _brg_table="${1:-1}"

  local _items=(
    _species
    _species_italic
    _genus
    _order
    _family
    _target_coverage
    _l_sra
    _l_sra_size_gb
    _long_coverage
    _s_sra
    _s_sra_size_gb
    _short_coverage
    _genome_size
    _memory_gb_genomesize
    _total_hours_genomesize
    _memory_gb_getorganelle
    _total_hours_getorganelle
    _memory_gb_msbwt
    _total_hours_msbwt
    _memory_gb_ptgaul
    _total_hours_ptgaul
    _memory_gb_nextdenovo_polish
    _total_hours_nextdenovo_polish
    _memory_gb_oatk_nextdenovo
    _total_hours_oatk_nextdenovo
    _memory_gb_tippo_nextdenovo
    _total_hours_tippo_nextdenovo
    _memory_gb_pmat_nextdenovo
    _total_hours_pmat_nextdenovo
    _memory_gb_polap_disassemble_default
    _total_hours_polap_disassemble_default
    _memory_gb_polap_disassemble_simple
    _total_hours_polap_disassemble_simple
    _memory_gb_polap_disassemble_check_default
    _total_hours_polap_disassemble_check_default
    _memory_gb_polap_disassemble_check_simple
    _total_hours_polap_disassemble_check_simple
    _known_mtdna
    _seq_length_ptgaul
    _num_seq_ptgaul
    _seq_length_subsample
    _num_seq_subsample
    _D
    _P
    _N
    _R
    _Alpha
    _Memory
    _extracted_memory
    _n1
    _mode1
    _sd1
    _index1
    _n2
    _index2
    _summary2_rate_rounded
    _summary2_size_gb
    _summary2_alpha_formatted
    _mafft_pident
  )

  # local _items=(
  # 	species
  # 	genus
  # 	order
  # 	family
  # 	target_coverage
  # 	l_sra
  # 	l_sra_size_gb
  # 	long_coverage
  # 	s_sra
  # 	s_sra_size_gb
  # 	short_coverage
  # 	genome_size
  # 	memory_gb_genomesize
  # 	total_hours_genomesize
  # 	memory_gb_getorganelle
  # 	total_hours_getorganelle
  # 	memory_gb_msbwt
  # 	total_hours_msbwt
  # 	memory_gb_ptgaul
  # 	total_hours_ptgaul
  # 	memory_gb_nextdenovo_polish
  # 	total_hours_nextdenovo_polish
  # 	memory_gb_oatk_nextdenovo
  # 	total_hours_oatk_nextdenovo
  # 	memory_gb_tippo_nextdenovo
  # 	total_hours_tippo_nextdenovo
  # 	memory_gb_pmat_nextdenovo
  # 	total_hours_pmat_nextdenovo
  # 	memory_gb_polap_disassemble_default
  # 	total_hours_polap_disassemble_default
  # 	memory_gb_polap_disassemble_simple
  # 	total_hours_polap_disassemble_simple
  # 	memory_gb_polap_disassemble_check_default
  # 	total_hours_polap_disassemble_check_default
  # 	memory_gb_polap_disassemble_check_simple
  # 	total_hours_polap_disassemble_check_simple
  # 	known_mtdna
  # 	seq_length_ptgaul
  # 	num_seq_ptgaul
  # 	seq_length_subsample
  # 	num_seq_subsample
  # 	D
  # 	P
  # 	N
  # 	R
  # 	Alpha
  # 	Memory
  # 	extracted_memory
  # 	n1
  # 	mode1
  # 	sd1
  # 	index1
  # 	n2
  # 	index2
  # 	summary2_rate_rounded
  # 	summary2_size_gb
  # 	summary2_alpha_formatted
  # 	mafft_pident
  # )

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
man-table-benchmark_genus_species_for() {
  local _brg_outdir="${1}"
  local _brg_inum="${2:-0}"
  local _brg_type="${3:-benchmark}"

  local _brg_d_index="${3:-infer-1}"
  local _brg_table="${4:-1}"
  local _brg_t_dir="${5:-"${_brg_default_target_dir}"}"
  local _key="${_brg_outdir}-${_brg_inum}"

  # Used the upper bound of the memory
  local _extracted_memory="${_memory["$_key"]}"
  local _v1="${_brg_outdir}"
  local _polap_log=${_v1}/polap.log

  # Folders
  local _brg_outdir_t="${_brg_outdir}/${opt_t_arg}"
  local _brg_outdir_i="${_brg_outdir_t}/${_brg_inum}"
  local _brg_outdir_0="${_brg_outdir_t}/0"

  # Taxon names: species, genus, family, and order
  local _species=${_taxon[$_key]}
  local _species="${_species//_/ }"
  local _species_italic="_${_species}_"
  local _genus=${_species%% *}
  local _order=$(grep "${_genus}" ${_POLAPLIB_DIR}/taxonomy_output.tsv | cut -f 6 | head -1)
  local _family=$(grep "${_genus}" ${_POLAPLIB_DIR}/taxonomy_output.tsv | cut -f 7 | head -1)

  # Extract long-read SRA ID
  local _l_fq_stats="${_brg_outdir_0}/summary-data/l.fq.stats"
  local _l_sra=$(_polap_lib_extract-seqkit_stats "${_l_fq_stats}")
  local _l_sra_size=$(_polap_lib_extract-seqkit_stats_sum_len "${_l_fq_stats}")
  local _l_sra_size_gb=$(_polap_lib_unit-convert_bp "${_l_sra_size}")

  # Extract short-read SRA ID
  local _s1_fq_stats="${_brg_outdir_0}/summary-data/s1.fq.stats"
  local _s1_sra=$(_polap_lib_extract-seqkit_stats "${_s1_fq_stats}")
  local _s_sra="${_s1_sra%%_*}"
  local _s1_sra_size=$(_polap_lib_extract-seqkit_stats_sum_len "${_s1_fq_stats}")
  local _s1_sra_size_gb=$(_polap_lib_unit-convert_bp "${_s1_sra_size}")
  local _s2_fq_stats="${_brg_outdir_0}/summary-data/s2.fq.stats"
  local _s2_sra=$(_polap_lib_extract-seqkit_stats "${_s2_fq_stats}")
  local _s2_sra_size=$(_polap_lib_extract-seqkit_stats_sum_len "${_s2_fq_stats}")
  local _s2_sra_size_gb=$(_polap_lib_unit-convert_bp "${_s2_sra_size}")
  local _s_sra_size=$((_s1_sra_size + _s2_sra_size))
  local _s_sra_size_gb=$(_polap_lib_unit-convert_bp "${_s_sra_size}")

  # Extract the genome size estimate
  local _genome_size_dir="${_brg_outdir_0}/estimate-genomesize"
  local _genome_size=$(_polap_lib_extract-genome_size "${_genome_size_dir}")
  local _genome_size_gb=$(_polap_lib_unit-convert_bp "${_genome_size}")

  # sequencing data coverage
  local _long_coverage=$(echo "scale=2; ${_l_sra_size} / ${_genome_size}" | bc)
  local _short_coverage=$(echo "scale=2; ${_s_sra_size} / ${_genome_size}" | bc)

  # sequencing data coverage
  local _target_coverage=$(_polap_lib_extract-target_read_coverage "${_brg_outdir_i}/polap-disassemble/${_brg_inum}/lx.txt")

  # Extract the timing for the corresponing coverage x case
  local _timing_dir="${_brg_outdir_0}"

  # Extract the time and memory for genome size estimate
  _polap_lib_extract-time_memory_timing_summary_file \
    "${_timing_dir}/summary-estimate-genomesize.txt"
  local _memory_gb_genomesize=$_memory_gb
  local _total_hours_genomesize=$_total_hours

  # Extract the time and memory for GetOrganelle
  _polap_lib_extract-time_memory_timing_summary_file \
    "${_timing_dir}/summary-getorganelle.txt"
  local _memory_gb_getorganelle=$_memory_gb
  local _total_hours_getorganelle=$_total_hours

  # Extract the time and memory for FMLRC's msbwt
  _polap_lib_extract-time_memory_timing_summary_file \
    "${_timing_dir}/summary-msbwt.txt"
  local _memory_gb_msbwt=$_memory_gb
  local _total_hours_msbwt=$_total_hours

  # Extract the time and memory for ptgaul
  _polap_lib_extract-time_memory_timing_summary_file \
    "${_timing_dir}/summary-ptgaul.txt"
  local _memory_gb_ptgaul=$_memory_gb
  local _total_hours_ptgaul=$_total_hours

  # Extract the time and memory for nextdenovo-polish
  _polap_lib_extract-time_memory_timing_summary_file \
    "${_timing_dir}/summary-nextdenovo-polish.txt"
  local _memory_gb_nextdenovo_polish=$_memory_gb
  local _total_hours_nextdenovo_polish=$_total_hours

  # Extract the time and memory for oatk-nextdenovo
  _polap_lib_extract-time_memory_timing_summary_file \
    "${_timing_dir}/summary-oatk-nextdenovo.txt"
  local _memory_gb_oatk_nextdenovo=$_memory_gb
  local _total_hours_oatk_nextdenovo=$_total_hours

  # Extract the time and memory for tippo-nextdenovo
  _polap_lib_extract-time_memory_timing_summary_file \
    "${_timing_dir}/summary-tippo-nextdenovo.txt"
  local _memory_gb_tippo_nextdenovo=$_memory_gb
  local _total_hours_tippo_nextdenovo=$_total_hours

  # Extract the time and memory for pmat-nextdenovo
  _polap_lib_extract-time_memory_timing_summary_file \
    "${_timing_dir}/summary-pmat-nextdenovo.txt"
  local _memory_gb_pmat_nextdenovo=$_memory_gb
  local _total_hours_pmat_nextdenovo=$_total_hours

  local _timing_dir="${_brg_outdir_i}"

  # Extract the time and memory for polap-disassemble-default
  _polap_lib_extract-time_memory_timing_summary_file \
    "${_timing_dir}/summary-polap-disassemble-default.txt"
  local _memory_gb_polap_disassemble_default=$_memory_gb
  local _total_hours_polap_disassemble_default=$_total_hours

  # Extract the time and memory for polap-disassemble-simple
  _polap_lib_extract-time_memory_timing_summary_file \
    "${_timing_dir}/summary-polap-disassemble-simple.txt"
  local _memory_gb_polap_disassemble_simple=$_memory_gb
  local _total_hours_polap_disassemble_simple=$_total_hours

  # Extract the time and memory for polap-disassemble-check-default
  if [[ -s "${_timing_dir}/summary-polap-disassemble-check-default.txt" ]]; then
    _polap_lib_extract-time_memory_timing_summary_file \
      "${_timing_dir}/summary-polap-disassemble-check-default.txt"
  else
    _memory_gb=NA
    _total_hours=NA
  fi
  local _memory_gb_polap_disassemble_check_default=$_memory_gb
  local _total_hours_polap_disassemble_check_default=$_total_hours

  # Extract the time and memory for polap-disassemble-check-simple
  if [[ -s "${_timing_dir}/summary-polap-disassemble-check-simple.txt" ]]; then
    _polap_lib_extract-time_memory_timing_summary_file \
      "${_timing_dir}/summary-polap-disassemble-check-simple.txt"
  else
    _memory_gb=NA
    _total_hours=NA
  fi
  local _memory_gb_polap_disassemble_check_simple=$_memory_gb
  local _total_hours_polap_disassemble_check_simple=$_total_hours

  # Known ptDNA NCBI accession
  local _known_mtdna=$(<"${_brg_outdir_0}/ncbi-ptdna/00-bioproject/2-mtdna.accession")
  _known_mtdna=${_known_mtdna:-'NA'}

  local _ptdna_ptgaul="${_brg_outdir_0}/ptdna-ptgaul.fa"
  local _seq_length_ptgaul=$(_polap_lib_extract-fasta_seqlen "${_ptdna_ptgaul}")
  local _num_seq_ptgaul=$(_polap_lib_extract-fasta_numseq "${_ptdna_ptgaul}")

  # NOTE: use tools like seqkit not command tools such as awk
  # then, delete these two lines.
  local _seq_length_xyz=$(_polap_lib_extract-fasta_seqlen "${_brg_outdir_i}/xyz.fa")
  local _num_seq_xyz=$(_polap_lib_extract-fasta_numseq "${_brg_outdir_i}/xyz.fa")

  local _ptdna_subsample="${_brg_outdir_i}/polap-disassemble/${_brg_inum}/disassemble/infer-1/pt.subsample-polishing.1.fa"
  local _seq_length_subsample=$(_polap_lib_extract-fasta_seqlen "${_ptdna_subsample}")
  local _num_seq_subsample=$(_polap_lib_extract-fasta_numseq "${_ptdna_subsample}")

  # NOTE: simplify this part.
  # Files to parse in polap-disassemble result.
  local _polap_disassemble_dir="${_brg_outdir_i}/polap-disassemble/${_brg_inum}/disassemble/infer-1"
  local _params_txt="${_polap_disassemble_dir}/params.txt"
  local _summary1_ordered_txt="${_polap_disassemble_dir}/1/summary1-ordered.txt"
  local _summary2_ordered_txt="${_polap_disassemble_dir}/2/summary1-ordered.txt"

  # Extract polap-disassemble parameters
  local _ipn=$(parse_params "${_params_txt}")
  local _I _P _N _R _A _B _M _D _Alpha _Memory
  read -r _I _P _N _R _A _B _M _D _Alpha _Memory <<<"$_ipn"

  # Extract mode value
  # Extract SD value
  # Extract the first index value
  local _n1=$(grep "^#n:" "$_summary1_ordered_txt" | awk 'NR==1 {print $2}')
  local _mode1=$(grep "^#mode:" "$_summary1_ordered_txt" | awk '{print $2}')
  local _sd1=$(grep "^#sd:" "$_summary1_ordered_txt" | awk '{print $2}')
  local _index1=$(grep "^#index:" "$_summary1_ordered_txt" | awk 'NR==1 {print $2}')

  # NOTE: so which one is selected for the next stage?
  # problem: the first sorted and the selected one are different.
  # use either one of the two.
  # Then, we should not check the following unnecessary step of checking.
  # Use #index the first one and do not sort them.
  # Or, sort them and use the first one. Problem is we are not sure which index
  # is used.
  local _output=$(awk -F'\t' 'NR==2 {print $1}' "${_summary1_ordered_txt}")
  read -r _second_line_index <<<"$_output"
  if ((_index1 != _second_line_index)); then
    echo "ERROR: #index: ${_index1} vs. #2nd line index: ${_second_line_index}" >&2
    echo "See ${_summary1_ordered_txt}" >&2
    echo "------------------------------" >&2
    cat "${_summary1_ordered_txt}" >&2
    echo "------------------------------> skip it!" >&2
    # exit 1
  fi

  # The determined subsampling rate and alpha
  local _n2=$(grep "^#n:" "$_summary2_ordered_txt" | awk 'NR==1 {print $2}')
  local _output=$(awk -F'\t' 'NR==2 {print $1, $2, $4, $11}' "${_summary2_ordered_txt}")
  local _index2
  local _summary2_rate
  local _summary2_size
  local _summary2_alpha
  read -r _index2 _summary2_size _summary2_rate _summary2_alpha <<<"$_output"
  local _summary2_rate_decimal=$(printf "%.10f" "$_summary2_rate")
  local _summary2_rate_rounded=$(echo "scale=4; $_summary2_rate_decimal / 1" | bc)
  local _summary2_size_gb=$(_polap_lib_unit-convert_bp "${_summary2_size}")
  local _summary2_alpha_formatted=$(echo "scale=2; $_summary2_alpha / 1" | bc | awk '{printf "%.2f\n", $1}')

  # Percent identity from pairwise sequence alignment
  local _mafft_pindent_txt="${_brg_outdir_i}/polap-disassemble/mafft/1/pident.txt"
  local _mafft_pident="NA"
  if [[ -s "${_mafft_pindent_txt}" ]]; then
    _mafft_pident=$(<"${_mafft_pindent_txt}")
  fi

  # Check any errors
  if [[ -s "${_ptdna_subsample}" ]]; then
    # Ensure there is exactly one sequence
    if ((_num_seq_subsample != 1)); then
      echo "Error: FASTA file does not contain exactly one sequence: ${_ptdna_subsample}"
      exit 1
    fi
  else
    echo "Error: no such file: ${_ptdna_subsample}" >&2
    exit 1
  fi

  # debug local variables
  # for debugging: Inline local printing local var
  # while IFS= read -r line; do
  # 	if [[ $line =~ ^declare\ --\ ([^=]+)= ]]; then
  # 		var="${BASH_REMATCH[1]}"
  # 		printf "%s=%q\n" "$var" "${!var}"
  # 	fi
  # done < <(local -p 2>/dev/null)
  # return

  local _items=(
    "${_species}"
    "${_species_italic}"
    "${_genus}"
    "${_order}"
    "${_family}"
    "${_target_coverage}"
    "${_l_sra}"
    "${_l_sra_size_gb}"
    "${_long_coverage}"
    "${_s_sra}"
    "${_s_sra_size_gb}"
    "${_short_coverage}"
    "${_genome_size}"
    "${_memory_gb_genomesize}"
    "${_total_hours_genomesize}"
    "${_memory_gb_getorganelle}"
    "${_total_hours_getorganelle}"
    "${_memory_gb_msbwt}"
    "${_total_hours_msbwt}"
    "${_memory_gb_ptgaul}"
    "${_total_hours_ptgaul}"
    "${_memory_gb_nextdenovo_polish}"
    "${_total_hours_nextdenovo_polish}"
    "${_memory_gb_oatk_nextdenovo}"
    "${_total_hours_oatk_nextdenovo}"
    "${_memory_gb_tippo_nextdenovo}"
    "${_total_hours_tippo_nextdenovo}"
    "${_memory_gb_pmat_nextdenovo}"
    "${_total_hours_pmat_nextdenovo}"
    "${_memory_gb_polap_disassemble_default}"
    "${_total_hours_polap_disassemble_default}"
    "${_memory_gb_polap_disassemble_simple}"
    "${_total_hours_polap_disassemble_simple}"
    "${_memory_gb_polap_disassemble_check_default}"
    "${_total_hours_polap_disassemble_check_default}"
    "${_memory_gb_polap_disassemble_check_simple}"
    "${_total_hours_polap_disassemble_check_simple}"
    "${_known_mtdna}"
    "${_seq_length_ptgaul}"
    "${_num_seq_ptgaul}"
    "${_seq_length_subsample}"
    "${_num_seq_subsample}"
    "${_D}"
    "${_P}"
    "${_N}"
    "${_R}"
    "${_Alpha}"
    "${_Memory}"
    "${_extracted_memory}"
    "${_n1}"
    "${_mode1}"
    "${_sd1}"
    "${_index1}"
    "${_n2}"
    "${_index2}"
    "${_summary2_rate_rounded}"
    "${_summary2_size_gb}"
    "${_summary2_alpha_formatted}"
    "${_mafft_pident}"
  )

  printf "%s\t" "${_items[@]::${#_items[@]}-1}"
  printf "%s\n" "${_items[-1]}"

  # local _items=(
  # 	_species
  # 	_species_italic
  # 	_genus
  # 	_order
  # 	_family
  # 	_target_coverage
  # 	_l_sra
  # 	_l_sra_size_gb
  # 	_long_coverage
  # 	_s_sra
  # 	_s_sra_size_gb
  # 	_short_coverage
  # 	_genome_size
  # 	_memory_gb_genomesize
  # 	_total_hours_genomesize
  # 	_memory_gb_getorganelle
  # 	_total_hours_getorganelle
  # 	_memory_gb_msbwt
  # 	_total_hours_msbwt
  # 	_memory_gb_ptgaul
  # 	_total_hours_ptgaul
  # 	_memory_gb_nextdenovo_polish
  # 	_total_hours_nextdenovo_polish
  # 	_memory_gb_oatk_nextdenovo
  # 	_total_hours_oatk_nextdenovo
  # 	_memory_gb_tippo_nextdenovo
  # 	_total_hours_tippo_nextdenovo
  # 	_memory_gb_pmat_nextdenovo
  # 	_total_hours_pmat_nextdenovo
  # 	_memory_gb_polap_disassemble_default
  # 	_total_hours_polap_disassemble_default
  # 	_memory_gb_polap_disassemble_simple
  # 	_total_hours_polap_disassemble_simple
  # 	_memory_gb_polap_disassemble_check_default
  # 	_total_hours_polap_disassemble_check_default
  # 	_memory_gb_polap_disassemble_check_simple
  # 	_total_hours_polap_disassemble_check_simple
  # 	_known_mtdna
  # 	_seq_length_ptgaul
  # 	_num_seq_ptgaul
  # 	_seq_length_subsample
  # 	_num_seq_subsample
  # 	_D
  # 	_P
  # 	_N
  # 	_R
  # 	_Alpha
  # 	_Memory
  # 	_extracted_memory
  # 	_n1
  # 	_mode1
  # 	_sd1
  # 	_index1
  # 	_n2
  # 	_index2
  # 	_summary2_rate_rounded
  # 	_summary2_size_gb
  # 	_summary2_alpha_formatted
  # 	_mafft_pident
  # )
  #
  # printf "%s\t" "${!_items[@]::${#_items[@]}-1}"
  # printf "%s\n" "${!_items[-1]}"
}

# Table 1's row for a given outdir and inum.
#
# arg1: outdir
# arg2: inum
# arg3: disassemble-i
# arg4: 1 or 2
# arg5: target directory to copy the result
#
man-table-benchmark_genus_species() {
  local _brg_outdir="${1:-all}"
  local _brg_inum="${2:-0}"
  local _brg_type="${3:-benchmark}"

  # NOTE: delete these
  local _brg_d_index="${3:-infer-1}"
  local _brg_table="${4:-1}"
  local _brg_format="${5:-1}"
  local _brg_t_dir="${6:-"${_brg_default_target_dir}"}"

  # Set the run title
  local full_name="${FUNCNAME[0]}"
  local middle_part="${full_name#man-}"
  local run_title="${middle_part%%_*}"
  local run_title="${run_title%-*}"
  local table_md="${run_title}-${_brg_type}-${_brg_outdir}-${_brg_inum}.md"
  local table_tsv="${run_title}-${_brg_type}-${_brg_outdir}-${_brg_inum}.tsv"

  if [[ "$_POLAP_DEBUG" == "0" ]]; then
    man-table-benchmark_genus_species_header >"${table_tsv}"
  else
    man-table-benchmark_genus_species_header
  fi

  if [[ "${_brg_outdir}" == "all" ]]; then
    for key in "${Sall[@]}"; do
      man-table-benchmark_genus_species_for "$key" "${@:2}" >>"${table_tsv}"
      # echo man-table-benchmark_genus_species_for "$key" "${@:2}" ">>${table_tsv}"
    done
  elif [[ "${_brg_outdir}" == "test" ]]; then
    for key in "${Stest[@]}"; do
      man-table-benchmark_genus_species_for "$key" "${@:2}" >>"${table_tsv}"
      # echo man-table-benchmark_genus_species_for "$key" "${@:2}" ">>${table_tsv}"
    done
  elif [[ "${_brg_outdir}" == "some" ]]; then
    for key in "${Ssome[@]}"; do
      man-table-benchmark_genus_species_for "$key" "${@:2}" >>"${table_tsv}"
      # echo man-table-benchmark_genus_species_for "$key" "${@:2}" ">>${table_tsv}"
    done
  elif [[ "${_brg_outdir}" == "each" ]]; then
    for key in "${Skeys[@]}"; do
      man-table-benchmark_genus_species_for "$key" "${@:2}" >>"${table_tsv}"
    done
  else
    if [[ "$_POLAP_DEBUG" == "0" ]]; then
      man-table-benchmark_genus_species_for "$@" >>"${table_tsv}"
    else
      man-table-benchmark_genus_species_for "$@"
    fi
  fi

  _polap_lib_conda-ensure_conda_env polap || exit 1

  # Table S1. data table
  if [[ "${_brg_type}" == "data" ]]; then
    csvtk -t cut -f _species_italic,_order,_family,_l_sra,_l_sra_size_gb,_long_coverage,_s_sra,_s_sra_size_gb,_short_coverage "${table_tsv}" |
      csvtk -t rename -f 1-9 -n Species,Order,Family,"Long SRA","Long Size","Long Coverage","Short SRA","Short Size","Short Coverage" |
      csvtk -t csv2md -a right -o "${table_md}"
    # Table S2. banchmarking table comparing with other tools
  elif [[ "${_brg_type}" == "benchmark-memory" ]]; then
    csvtk -t cut -f _species_italic,_long_coverage,_short_coverage,_memory_gb_getorganelle,_memory_gb_ptgaul,_memory_gb_nextdenovo_polish,_memory_gb_pmat_nextdenovo,_memory_gb_tippo_nextdenovo,_memory_gb_oatk_nextdenovo,_memory_gb_polap_disassemble_default "${table_tsv}" |
      csvtk -t rename -f 1-10 -n Species,"Long Coverage","Short Coverage",GetOrganelle,ptGAUL,NextDenovo,PMAT,TIPPo,Oatk,Polap |
      csvtk -t csv2md -a right -o "${table_md}"
    # Table S2. banchmarking table comparing with other tools
  elif [[ "${_brg_type}" == "benchmark-time" ]]; then
    csvtk -t cut -f _species_italic,_long_coverage,_short_coverage,_total_hours_getorganelle,_total_hours_ptgaul,_total_hours_nextdenovo_polish,_total_hours_pmat_nextdenovo,_total_hours_tippo_nextdenovo,_total_hours_oatk_nextdenovo,_total_hours_polap_disassemble_default "${table_tsv}" |
      csvtk -t rename -f 1-10 -n Species,"Long Coverage","Short Coverage",GetOrganelle,ptGAUL,NextDenovo,PMAT,TIPPo,Oatk,Polap |
      csvtk -t csv2md -a right -o "${table_md}"
    # Table S3. table for polap disassemble 5 % - n=r=5
  elif [[ "${_brg_type}" == "benchmark-polap" ]]; then
    csvtk -t cut -f _species_italic,_long_coverage,_short_coverage,_target_coverage,_P,_N,_R,_Alpha,_seq_length_ptgaul,_seq_length_subsample,_mafft_pident "${table_tsv}" |
      csvtk -t rename -f 1-11 -n Species,"Long Coverage","Short Coverage","Downsample Depth",P,N,R,Alpha,"Length ptGAUL","Length Polap","Percent identity" |
      csvtk -t csv2md -a right -o "${table_md}"
  fi

  conda deactivate

  echo "type: ${_brg_type}"
  echo "${table_tsv}"
  echo "${table_md}"
  echo cp -p "${table_md}" "${_brg_t_dir}"
}

man-table-polap-disassemble_genus_species_for() {
  local _brg_outdir="${1}"
  local _brg_inum="${2:-0}"
  local _key="${_brg_outdir}-${_brg_inum}"

  # Set the run title
  local full_name="${FUNCNAME[0]}"
  local middle_part="${full_name#man-}"
  local _run_title="${middle_part%%_*}"
  local table_md="${_run_title}-${_brg_outdir}.md"
  local table_tsv="${_run_title}.tsv"

  # local _brg_d_index="${3:-infer-1}"
  # local _brg_stage="${4:-x}"
  # local _brg_table="${5:-1}"
  # local _brg_t_dir="${6:-"${_brg_default_target_dir}"}"
  # local _brg_one="${7:-0}"

  # local _i=${_brg_inum}
  # local j="${_brg_stage}"

  local _species=${_taxon[$_key]}
  local _species="${_species//_/ }"

  # local _label_base="${_v1/_/-}"
  # _label_base=$(echo "$_label_base" | awk '{print tolower($0)}')

  # local _v1_inum="${_v1}/${_brg_inum}"
  # local _params_txt=${_v1_inum}/disassemble/${_brg_d_index}/params.txt
  # local _ipn=$(parse_params "${_params_txt}")
  # local _I _P _N _R _A _B _M _D _Alpha _Memory
  # read -r _I _P _N _R _A _B _M _D _Alpha _Memory <<<"$_ipn"
  # local _label="${_brg_inum}-${_label_base}"
  # if [[ "${_brg_one}" == "1" ]]; then
  #   _label="main-${_brg_inum}-${_label_base}"
  # fi

  # printf "Table: Three stages of subsampling-based plastid genome assembly for the _${_species}_ dataset. The configuration includes an increasing subsample size up to a maximum subsampling rate of ${_P}%%, a step size of ${_N} in Stage 1, ${_R} replicates in Stages 2 and 3, and a maximum memory limit of ${_Memory} GB. {#tbl:supptable1-${_label}}\n\n"

  # Result folders
  local _brg_outdir_i="${_brg_outdir}/${opt_t_arg}/${_brg_inum}"
  local _brg_rundir="${_brg_outdir_i}/polap-disassemble"
  local _brg_rundir_i="${_brg_rundir}/0"
  local x1="${_brg_rundir_i}/disassemble/infer-1/1/summary1.txt"
  local x2="${_brg_rundir_i}/disassemble/infer-1/2/summary1.txt"
  local x3="${_brg_rundir_i}/disassemble/infer-1/3-infer/summary1.txt"

  # Initialize Conda and activate polap environment
  source "$(conda info --base)/etc/profile.d/conda.sh"
  if [[ "${CONDA_DEFAULT_ENV:-}" != "polap" ]]; then
    if [[ "$_POLAP_DEBUG" == "1" ]]; then
      echo "[INFO] Activating conda environment 'polap'..."
    fi
    conda activate polap
  fi

  # Check the activation of the polap conda environment
  if [[ "${CONDA_DEFAULT_ENV:-}" != "polap" ]]; then
    echo "[ERROR] Failed to enter conda environment 'polap'"
    return
  fi

  # Stage Rate Alpha D N L C Length Memory Time
  {
    csvtk -t cut -f index,long_rate_sample,alpha,draft_assembly_size,gfa_number_segments,gfa_total_segment_length,num_circular_paths,peak_ram_size_gb,time,length ${x1} |
      csvtk -t rename -f 1-10 -n Index,Rate,Alpha,DraftL,N,L,C,Memory,Time,Length |
      sed 's/-1/NA/g' |
      csvtk -t round -n 2 -f Rate |
      csvtk -t round -n 2 -f Alpha |
      csvtk -t mutate2 -n Stage -e "'1'" --at 1

    csvtk -t cut -f index,long_rate_sample,alpha,draft_assembly_size,gfa_number_segments,gfa_total_segment_length,num_circular_paths,peak_ram_size_gb,time,length ${x2} |
      csvtk -t rename -f 1-10 -n Index,Rate,Alpha,DraftL,N,L,C,Memory,Time,Length |
      sed 's/-1/NA/g' |
      csvtk -t round -n 2 -f Rate |
      csvtk -t round -n 2 -f Alpha |
      csvtk -t mutate2 -n Stage -e "'2'" --at 1 |
      csvtk del-header

    csvtk -t cut -f index,short_rate_sample,pident,short_total,short_total,short_total,short_total,memory_prepare,time_prepare,length ${x3} |
      csvtk -t rename -f 1-10 -n Index,Rate,Alpha,DraftL,N,L,C,Memory,Time,Length |
      sed 's/-1/NA/g' |
      csvtk -t round -n 2 -f Rate |
      csvtk -t round -n 2 -f Alpha |
      csvtk -t mutate2 -n Stage -e "'3'" --at 1 |
      csvtk del-header
  } | csvtk -t csv2md -a right -o ${table_md}

  conda deactivate

  echo "${table_md}"

  # printf "\n" \
  #   >>"${_supptable_md}"
  #
  # cat "${_POLAPLIB_DIR}"/polap-data-v2-suptable1_footnote.tex \
  #   >>"${_supptable_md}"
  #
  # printf "\n\n\\\newpage\n\n" \
  #   >>"${_supptable_md}"

}

man-table-polap-disassemble_genus_species() {
  local _brg_outdir="${1:-all}"
  local _brg_inum="${2:-0}"

  # Set the run title
  local full_name="${FUNCNAME[0]}"
  local middle_part="${full_name#man-}"
  local _run_title="${middle_part%%_*}"
  local table_md="${_run_title}.md"

  local _brg_d_index="${3:-infer-1}"
  local _brg_stage="${4:-x}"
  local _brg_table="${5:-1}"
  local _brg_t_dir="${6:-"${_brg_default_target_dir}"}"

  rm -f "${table_md}"

  if [[ "${_brg_outdir}" == "all" ]]; then
    for _v1 in "${Sall[@]}"; do
      man-table-polap-disassemble_genus_species_for "${_v1}" "${@:2}"
    done
  else
    man-table-polap-disassemble_genus_species_for "$@"
  fi
}

################################################################################
# main cases
#
# if [[ "${subcmd1}" == "help" ]]; then
#   if [[ "${_arg2}" == "arg2" ]]; then
#     subcmd1="help"
#   else
#     subcmd1="${_arg2}"
#     _arg2="arg2"
#   fi
# fi

if [[ "$_POLAP_DEBUG" == "1" ]]; then
  echo "option -y: ${opt_y_flag}"
  echo "option -c: ${opt_c_arg}"
  echo "option -t: ${opt_t_arg}"
  for item in "$@"; do
    echo "A: $item"
  done
fi

all_args=("$@")               # Save all arguments to an array
cmd_args=("${all_args[@]:1}") # Slice from index 1 onward

# Call common case first
common_handled=1
# _polap_lib_data-execute-common-subcommand "$subcmd1" "${_arg2}" "${_arg3}" "$opt_y_flag"
_polap_lib_data-execute-common-subcommand "$subcmd1" "$opt_y_flag" cmd_args
common_handled=$?

# Main case statement
case "$subcmd1" in
system)
  ${subcmd1}_genus_species
  ;;
copy-figures)
  ${subcmd1}_genus_species
  ;;
coverage)
  ${subcmd1}_genus_species "${_arg2}" "${_arg3}"
  ;;
extract-ptgaul-ptdna)
  if [[ "${_arg2}" == arg2 ]]; then
    echo "Help: ${subcmd1} <outdir>"
    echo "  polap-data-v2.sh ${subcmd1} Arabidopsis_thaliana"
    echo "${help_message_extract_ptgaul_ptdna}"
    exit 0
  fi
  ${subcmd1}_genus_species "${_arg2}"
  ;;
sample-csv)
  if [[ "${_arg2}" == arg2 ]]; then
    echo "Help: ${subcmd1} <csv:polap-data-v2.csv> [number:1|2|7|test|all|each] [force:off|on] [inum:0|N]"
    echo "  $(basename $0) ${subcmd1} 1.csv 1"
    echo "${help_message_sample_csv}"
    exit 0
  fi
  [[ "${_arg3}" == arg3 ]] && _arg3="1"
  [[ "${_arg4}" == arg4 ]] && _arg4="off"
  [[ "${_arg5}" == arg5 ]] && _arg5="0"

  if [[ -s "${_arg2}" ]] && [[ "${_arg4}" == "off" ]]; then
    echo "ERROR: you already have ${_arg2}"
    echo "  delete it if you want to create a new one."
  else
    if [[ "${_arg3}" == "1" ]]; then
      head -1 ${_POLAPLIB_DIR}/polap-data-v2.csv >"${_arg2}"
      grep Spirodela_polyrhiza ${_POLAPLIB_DIR}/polap-data-v2.csv >>"${_arg2}"
    elif [[ "${_arg3}" == "2" ]]; then
      head -1 ${_POLAPLIB_DIR}/polap-data-v2.csv >"${_arg2}"
      grep Spirodela_polyrhiza ${_POLAPLIB_DIR}/polap-data-v2.csv >>"${_arg2}"
      grep Eucalyptus_pauciflora ${_POLAPLIB_DIR}/polap-data-v2.csv >>"${_arg2}"
    elif [[ "${_arg3}" == "test" ]]; then
      head -1 ${_POLAPLIB_DIR}/polap-data-v2.csv >"${_arg2}"
      grep test ${_POLAPLIB_DIR}/polap-data-v2.csv >>"${_arg2}"
    elif [[ "${_arg3}" == "7" ]]; then
      head -1 ${_POLAPLIB_DIR}/polap-data-v2.csv >"${_arg2}"
      for item in $(printf "%s\n" "${S7[@]}" | sort); do
        grep $item ${_POLAPLIB_DIR}/polap-data-v2.csv |
          grep -v test >>"${_arg2}"
      done
    elif [[ "${_arg3}" == "each" ]]; then
      extract_and_replace_suffix "${_POLAPLIB_DIR}"/polap-data-v2.csv "-${_arg5}" "${_arg2}"
    elif [[ "${_arg3}" == "all" ]]; then
      grep -v test ${_POLAPLIB_DIR}/polap-data-v2.csv >"${_arg2}"
    fi
    echo "create CSV: ${_arg2}"
  fi
  ;;
archive)
  if [[ "${_arg2}" == arg2 ]]; then
    echo "Help: ${subcmd1} <outdir|all> <inum:-1|N>"
    echo "  polap-data-v2.sh ${subcmd1} all"
    echo "  polap-data-v2.sh ${subcmd1} Arabidopsis_thaliana"
    echo "  polap-data-v2.sh ${subcmd1} Arabidopsis_thaliana 0"
    echo "${help_message_archive}"
    exit 0
  fi
  [[ "${_arg3}" == arg3 ]] && _arg3=""
  ${subcmd1}_genus_species "${_arg2}" "${_arg3}"
  ;;
benchmark-command)
  if [[ "${_arg2}" == arg2 ]]; then
    echo "Help: ${subcmd1} <outdir> [inum:0|N]"
    echo "  $(basename $0) ${subcmd1} Arabidopsis_thaliana"
    _subcmd1_clean="${subcmd1//-/_}"
    declare -n ref="help_message_${_subcmd1_clean}"
    echo "$ref"
    exit 0
  fi
  [[ "${_arg3}" == arg3 ]] && _arg3=""
  ${subcmd1}_genus_species "${_arg2}" "${_arg3}"
  ;;
man)
  if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
    echo "Help: ${subcmd1} <command> ..."
    echo "  $(basename $0) ${subcmd1} table-benchmark"
    _subcmd1_clean="${subcmd1//-/_}"
    declare -n ref="help_message_${_subcmd1_clean}"
    echo "$ref"
    exit 0
  fi
  ${subcmd1}_genus_species "${cmd_args[@]}"
  ;;
setup-csv)
  if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
    echo "Help: ${subcmd1} <${_polap_data_csv}>"
    _subcmd1_clean="${subcmd1//-/_}"
    declare -n ref="help_message_${_subcmd1_clean}"
    echo "$ref"
    exit 0
  fi
  ${subcmd1}_genus_species "${cmd_args[@]}"
  ;;
man-figure-delta)
  if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
    echo "Help: ${subcmd1} <outdir>"
    echo "  ${0} ${subcmd1} Eucalyptus_pauciflora"
    _subcmd1_clean="${subcmd1//-/_}"
    declare -n ref="help_message_${_subcmd1_clean}"
    echo "$ref"
    exit 0
  fi
  ${subcmd1}_genus_species "${cmd_args[@]}"
  ;;
man-figure-alpha)
  if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
    echo "Help: ${subcmd1} <outdir>"
    echo "  ${0} ${subcmd1} Eucalyptus_pauciflora"
    _subcmd1_clean="${subcmd1//-/_}"
    declare -n ref="help_message_${_subcmd1_clean}"
    echo "$ref"
    exit 0
  fi
  ${subcmd1}_genus_species "${cmd_args[@]}"
  ;;
man-table-data)
  if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
    echo "Help: ${subcmd1} <all|outdir> [inum:0|N]"
    echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
    _subcmd1_clean="${subcmd1//-/_}"
    declare -n ref="help_message_${_subcmd1_clean}"
    echo "$ref"
    exit 0
  fi
  ${subcmd1}_genus_species "${cmd_args[@]}"
  ;;
man-figure-polap)
  if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
    echo "Help: ${subcmd1} <outdir> [inum:0|N]"
    echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
    _subcmd1_clean="${subcmd1//-/_}"
    declare -n ref="help_message_${_subcmd1_clean}"
    echo "$ref"
    exit 0
  fi
  ${subcmd1}_genus_species "${cmd_args[@]}"
  ;;
man-figure-sheet-pmat)
  if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
    echo "Help: ${subcmd1} <outdir> [inum:0|N]"
    echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
    _subcmd1_clean="${subcmd1//-/_}"
    declare -n ref="help_message_${_subcmd1_clean}"
    echo "$ref"
    exit 0
  fi
  ${subcmd1}_genus_species "${cmd_args[@]}"
  ;;
man-figure-sheet-tippo)
  if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
    echo "Help: ${subcmd1} <outdir> [inum:0|N]"
    echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
    _subcmd1_clean="${subcmd1//-/_}"
    declare -n ref="help_message_${_subcmd1_clean}"
    echo "$ref"
    exit 0
  fi
  ${subcmd1}_genus_species "${cmd_args[@]}"
  ;;
man-figure-sheet-oatk)
  if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
    echo "Help: ${subcmd1} <outdir> [inum:0|N]"
    echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
    _subcmd1_clean="${subcmd1//-/_}"
    declare -n ref="help_message_${_subcmd1_clean}"
    echo "$ref"
    exit 0
  fi
  ${subcmd1}_genus_species "${cmd_args[@]}"
  ;;
man-figure-sheet-polap)
  if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
    echo "Help: ${subcmd1} <outdir> [inum:0|N]"
    echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
    _subcmd1_clean="${subcmd1//-/_}"
    declare -n ref="help_message_${_subcmd1_clean}"
    echo "$ref"
    exit 0
  fi
  ${subcmd1}_genus_species "${cmd_args[@]}"
  ;;
man-figure-sheet)
  if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
    echo "Help: ${subcmd1} <outdir> [inum:0|N]"
    echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
    _subcmd1_clean="${subcmd1//-/_}"
    declare -n ref="help_message_${_subcmd1_clean}"
    echo "$ref"
    exit 0
  fi
  ${subcmd1}_genus_species "${cmd_args[@]}"
  ;;
man-table-polap-disassemble)
  if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
    echo "Help: ${subcmd1} <analysis:t1|string>"
    echo "  $(basename $0) ${subcmd1}"
    _subcmd1_clean="${subcmd1//-/_}"
    declare -n ref="help_message_${_subcmd1_clean}"
    echo "$ref"
    exit 0
  fi
  ${subcmd1}_genus_species "${cmd_args[@]}"
  ;;
man-table-benchmark)
  if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
    echo "Help: ${subcmd1} <all|outdir> [inum:0|N] [benchmark|data]"
    echo "  $(basename $0) ${subcmd1}"
    _subcmd1_clean="${subcmd1//-/_}"
    declare -n ref="help_message_${_subcmd1_clean}"
    echo "$ref"
    exit 0
  fi
  ${subcmd1}_genus_species "${cmd_args[@]}"
  ;;
benchmark)
  if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
    echo "Help: ${subcmd1} <outdir> ..."
    echo "  $(basename $0) ${subcmd1}"
    _subcmd1_clean="${subcmd1//-/_}"
    declare -n ref="help_message_${_subcmd1_clean}"
    echo "$ref"
    exit 0
  fi
  ${subcmd1}_genus_species "${cmd_args[@]}"
  ;;
nextdenovo-polish-cleanup)
  if [[ -z "${_arg2}" || "${_arg2}" == arg2 ]]; then
    echo "Help: ${subcmd1} <outdir>"
    echo "  ${0##*/} ${subcmd1} Arabidopsis_thaliana"
    _subcmd1_clean="${subcmd1//-/_}"
    declare -n ref="help_message_${_subcmd1_clean}"
    echo "$ref"
    exit 0
  fi
  ${subcmd1}_genus_species "${_arg2}"
  ;;
test)
  if [[ "${_arg2}" == arg2 ]]; then
    echo "Help: ${subcmd1} <outdir|all> [inum:0|N]"
    echo "  polap-data-v2.sh ${subcmd1} all"
    echo "  polap-data-v2.sh ${subcmd1} Arabidopsis_thaliana"
    _subcmd1_clean="${subcmd1//-/_}"
    declare -n ref="help_message_${_subcmd1_clean}"
    echo "$ref"
    exit 0
  fi
  ${subcmd1}_genus_species "${cmd_args[@]}"
  ;;
"menu")
  _run_polap_menu
  ;;
download-test-data)
  if [[ "${opt_y_flag}" == false ]]; then
    read -p "Do you want to download test data for polap? (y/N): " confirm
  else
    confirm="yes"
  fi

  if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
    _s=polap-disassemble-test-full.tar.gz
    # local _s=polap-disassemble-test.tar.gz
    if [[ ! -s "${_s}" ]]; then
      # full test data
      curl -L -o "${_s}" "https://figshare.com/ndownloader/files/53457569?private_link=ec1cb394870c7727a2d4"
      #
      # test data
      # curl -L -o "${_s}" "https://figshare.com/ndownloader/files/53457566?private_link=ec1cb394870c7727a2d4"
    fi
    if [[ ! -s "l.fastq" ]]; then
      tar -zxf "${_s}"
      echo "downloaded: l.fastq, s_1.fastq, s_2.fastq"
    else
      echo "You already have: l.fastq"
    fi
  else
    echo "polap test download is canceled."
  fi
  ;;
download-man)
  if [[ "${opt_y_flag}" == false ]]; then
    read -p "Do you want to create the manuscript tables and figures? (y/N): " confirm
  else
    confirm="yes"
  fi

  if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
    echo "Creating the tables and figures ..."
    if [[ -d "polap" ]]; then
      rm -rf polap
    fi
    git clone --quiet https://github.com/goshng/polap.git
    cp -qpr polap/man .
    cd
  elif [[ "${confirm,,}" == "update" || "${confirm,,}" == "u" ]]; then
    cp $HOME/all/polap/github/man/v0.4/Makefile ../man/v0.4/
    cp $HOME/all/polap/github/man/v0.4/manuscript.md ../man/v0.4/
    cp $HOME/all/polap/github/man/v0.4/manuscript-supp.md ../man/v0.4/
  else
    echo "polap man is canceled."
  fi
  ;;
make-man)
  if [[ "${opt_y_flag}" == false ]]; then
    read -p "Do you want to make the manuscript tables and figures? (y/N): " confirm
  else
    confirm="yes"
  fi

  if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
    # Initialize Conda
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate polap-man
    cd ../man/v0.4
    make
    echo ../man/v0.4/manuscript.pdf
    conda deactivate
  else
    echo "polap man is canceled."
  fi
  ;;
*)
  # only print usage if common_case didn't handle it
  if [[ $common_handled -ne 0 ]]; then
    echo "Usage: $0 <subcommand> [species_folder]"
    echo "${help_message}"
    echo "subcommand '$subcmd1' is not recognized."
    exit 1
  fi
  ;;
esac
