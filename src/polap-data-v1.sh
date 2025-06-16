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

if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
  echo "[ERROR] This script must be executed, not sourced: use 'bash $BASH_SOURCE'" >&2
  return 1 2>/dev/null || exit 1
fi
: "${_POLAP_DEBUG:=0}"
export _POLAP_DEBUG
: "${_POLAP_RELEASE:=0}"
export _POLAP_RELEASE

# Data directories
_media1_dir="/media/h1/sra"
_media2_dir="/media/h2/sra"
_media_dir="/media/h2/sra"

# _polap_version="$(${_polap_cmd} --version | awk '{print $2}')"
if [ -z "${_polap_version+x}" ]; then
  _polap_version="0.3.7.3"
fi

help_message=$(
  cat <<HEREDOC
Polap version: ${_polap_version}
Usage: $0 [Options] <subcommand> [species_folder] [[option] ...]

Polap data analysis of subsampling-based plastid genome assembly

Options:
  -c <arg>       Set value for -c option (default: off)
  -t <arg>       Set value for -t option (default: off)
  -y             Enable -y flag
  -v                  Enable verbose mode.
  -f                  Enable -f flag to say YES to profiling.
  -e <name>      Call <name>_genus_species function and exit
  -h, --help     Show this help message

List of subcommands:
  install-conda, setup-conda, install-polap, install-fmlrc, install-cflye
  install-getorganelle, install-pmat, install-man
  download-polap-github, patch-polap, bleeding-edge-polap
  delete-polap-github, uninstall, test-polap, download-test-data, help, system
  mkdir-all, rm-empty, rm, example-data, convert-data, sample-csv
  sra, refs, getorganelle, ptgaul, msbwt,
  local-batch, remote-batch, batch, archive, clean, report, get
  maintable1, supptable1, suppfigure1, suppfigure3, suppfigure4
  download-man, make-man

0. test <SPECIES>
1. send-data-to <SPECIES>
2. mkdir <SPECIES>
3. taxon-reference <SPECIES> : ptgaul-reference.fa with mtDNA related with the species
3. taxon-known <SPECIES> : ptgaul-mtdna.fa with mtDNA of the species
4. ptgaul-reference <SPECIES>
4. ptgaul-known <SPECIES>
5. bandage-ptgaul-known <SPECIES>
5. bandage-ptgaul-reference <SPECIES>
. subsample <SPECIES>
. assemble1-subsample <SPECIES> <COVERAGE>
. bandage1-subsample <SPECIES> <COVERAGE>
. seeds-subsample <SPECIES> <COVERAGE>
. assemble2-subsample <SPECIES> <COVERAGE> [TARGET_ASSEMBLY_ID]
. bandage2-subsample <SPECIES> <COVERAGE> [TARGET_ASSEMBLY_ID]
. run <SPECIES>
write-config polap-data-v1.config

install-conda: install Miniconda3
setup-conda: setup Miniconda3 for Bioconda
install-polap or install: install Polap
install-fmlrc: install ptGAUL's FMLRC short-read polishing tools
patch-polap: update the miniconda3/envs/polap/bin/polap with the github version
Carex_pseudochinensis: test with mtDNA of Carex pseudochinensis
test-polap: test Polap run with a test dataset
download-polap: download Polap
clean: delete analyses
uninstall: uninstall Polap
mkdir: create the 11 species folders.
rm: delete the 11 species folders.
link-fastq: create links to the input data at ${_media_dir} for the 11 folders.
delete-links: delete all links in the current folder and its subfolders.
copy-fastq: copy the input data at ${_media_dir} to the 11 folders.
scopy-fastq <ssh-hostname>: transfer the input data at ${_media_dir} to the remote computers.
zip: compress the 11 species folders.
plot: create the supplementary figures for the -w option tests.
table: create a table for the 11 test data.
sync 0: copy the whole-genome assembly results from other computers to the local.
sync 1: copy the organelle-genome assembly results from other computers to the local.
---
BUG: FIX QT5 problem:
export QT_QPA_PLATFORM=offscreen

How to set a custom CSV:
csv_file=a.csv $0

More help for subcommands:
$(basename $0) help <subcommand>
HEREDOC
)

# Spirodela_polyrhiza
# Taraxacum_mongolicum
# Trifolium_pratense
# Salix_dunnii
# Anthoceros_agrestis
# Anthoceros_angustus
# Brassica_rapa
# Vigna_radiata
# Macadamia_tetraphylla
# Punica_granatum
# Lolium_perenne

# polap-analysis-without-wga
#
# Salix_dunnii: 20 takes too long, 10?
# Taraxacum_mongolicum: 20 has too many seeds, so we do not organelle-genome.
#
S2=(
  'Anthoceros_agrestis'
  'Vigna_radiata'
  'Salix_dunnii'
  'Taraxacum_mongolicum'
  'Brassica_rapa'
  'Spirodela_polyrhiza'
  'Trifolium_pratense'
)

S9=(
  'Spirodela_polyrhiza'
  'Taraxacum_mongolicum'
  'Trifolium_pratense'
  'Anthoceros_agrestis'
  'Anthoceros_angustus'
  'Brassica_rapa'
  'Vigna_radiata'
  'Macadamia_tetraphylla'
)

S=(
  'Spirodela_polyrhiza'
  'Taraxacum_mongolicum'
  'Trifolium_pratense'
  'Salix_dunnii'
  'Anthoceros_agrestis'
  'Anthoceros_angustus'
  'Brassica_rapa'
  'Vigna_radiata'
  'Macadamia_tetraphylla'
  'Punica_granatum'
  'Lolium_perenne'
)

declare -A pmat_fc

# Add key-value pairs
pmat_fc[Spirodela_polyrhiza]=0.1
pmat_fc[Taraxacum_mongolicum]=0.2
pmat_fc[Trifolium_pratense]=0.1
pmat_fc[Salix_dunnii]=0.02
pmat_fc[Anthoceros_agrestis]=0.1
pmat_fc[Anthoceros_angustus]=1.0
pmat_fc[Brassica_rapa]=0.1
pmat_fc[Vigna_radiata]=0.1
pmat_fc[Macadamia_tetraphylla]=0.1
pmat_fc[Punica_granatum]=0.1
pmat_fc[Lolium_perenne]=0.1

_local_host="thorne"

# OTHER STUFF GENERATED BY Argbash
_polap_script_bin_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" || {
  echo "Couldn't determine the script's running directory, which probably matters, bailing out" >&2
  exit 2
}
_POLAPLIB_DIR="${_polap_script_bin_dir}/polaplib"

source "${_POLAPLIB_DIR}/polap-lib-timing.sh"
source "${_POLAPLIB_DIR}/polap-lib-number.sh"
source "${_POLAPLIB_DIR}/polap-lib-data.sh"
source "${_POLAPLIB_DIR}/polap-lib-process.sh"
source <(echo 'export PATH="$PWD/bin:$PATH"')
# source <(echo 'export QT_QPA_PLATFORM=offscreen')
source <(echo 'export QT_QPA_PLATFORM=minimal')

# echo "export QT_QPA_PLATFORM=offscreen"

_polap_data_csv="$(basename "$0" .sh).csv"
_polap_data_data="$(basename "$0" .sh).data"
_polap_data_txt="$(basename "$0" .sh).txt"

_log_echo() {
  if [[ -s "${_brg_outdir}/${_polap_data_txt}" ]]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') [$subcmd1] - $1" >>"${_brg_outdir}/${_polap_data_txt}"
  fi
  # echo "$1"
}

data-example_genus_species() {
  local _brg_data="${1:-${_polap_data_data}}"

  local _data=$(
    cat <<HEREDOC
species,long,short
Carex_pseudochinensis,SRR30757341,SRR30757340
Lolium_perenne,SRR13386519,SRR13386518
Anthoceros_angustus,SRR9696346,SRR9662965
Macadamia_tetraphylla,SRR10424548,SRR10424549
Spirodela_polyrhiza,SRR11472010,SRR11472009
Vigna_radiata,SRR12549541,SRR12549533
Brassica_rapa,ERR6210792,ERR6210790
Trifolium_pratense,SRR15433794,SRR15433795
Anthoceros_agrestis,SRR10190639,SRR10250248
Punica_granatum,SRR24893686,SRR24893685
Taraxacum_mongolicum,SRR19182970,SRR19182971
Salix_dunnii,SRR12893432,SRR12893433
HEREDOC
  )

  echo "${_data}" >"${_brg_data}"
  echo "create ${_brg_data}"
}

############################################################
# main command arguments used before a subcommand
#
# Default values for options
opt_c_arg="off"
opt_t_arg="off"
opt_v_flag=false
opt_y_flag=false
opt_f_flag=false
opt_e_arg=""

print_help() {
  echo "$help_message"
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
  -y)
    opt_y_flag=true
    ;;
  -v)
    opt_v_flag=true
    ;;
  -f)
    opt_f_flag=true
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

_brg_default_target_dir="$HOME/all/manuscript/polap-v0.2/figures/"
if [[ -d "src" ]]; then
  _brg_default_target_dir="$HOME/all/manuscript/polap-v0.2/figures/"
else
  if [[ -d "../man" ]]; then
    _brg_default_target_dir="../man/v0.2/figures"
  else
    mkdir -p ../man/v0.2/figures
  fi
fi

if [[ "${opt_t_arg}" != "off" ]]; then
  _brg_default_target_dir="${opt_t_arg}"
fi

if [[ "${opt_c_arg}" != "off" ]]; then
  csv_file="${opt_c_arg}"
fi

# Input parameter
subcmd1="${1:-help}"

_polap_cmd="${_polap_script_bin_dir}/polap.sh"

##### INSERT_HELP_HERE #####
help_message_polap_batch_assemble=$(
  cat <<HEREDOC

  get ${_local_host}:$PWD/<outdir>-a.tar.gz
HEREDOC
)

help_message_polap_batch_polish=$(
  cat <<HEREDOC

  Batch steps
HEREDOC
)

help_message_polap_msbwt=$(
  cat <<HEREDOC

  get ${_local_host}:$PWD/<outdir>-a.tar.gz
HEREDOC
)

help_message_polap_batch=$(
  cat <<HEREDOC

  Batch steps
  1. polap-batch-assemble <outdir> 30
  2. polap-msbwt <outdir> 30 -> only once
  3. polap-archive <outdir> 30
  4. polap-batch-polish <outdir> 30
  5. polap-clean <outdir> 30

  31 because mt.contig.name-1 only one seed file
  x: means 20 and 30. 20 -> 1x, 30 -> 0.1x
HEREDOC
)

help_message_polap_archive=$(
  cat <<HEREDOC

  get ${_local_host}:$PWD/<outdir>-a.tar.gz
HEREDOC
)

help_message_polap_polish=$(
  cat <<HEREDOC

  get ${_local_host}:$PWD/<outdir>-a.tar.gz
HEREDOC
)

help_message_polap_assemble2=$(
  cat <<HEREDOC

  get ${_local_host}:$PWD/<outdir>-a.tar.gz
HEREDOC
)

help_message_polap_annotate=$(
  cat <<HEREDOC

  get ${_local_host}:$PWD/<outdir>-a.tar.gz
HEREDOC
)

help_message_polap_assemble1=$(
  cat <<HEREDOC

  outdir: species folder delimited by an underscore character
  inum: key appended to the out directory or the first column of a CSV config
  coverage: polap --coverage option for a whole-genome assembly
HEREDOC
)

help_message_get_data=$(
  cat <<HEREDOC

  get ${_local_host}:$PWD/<outdir>-a.tar.gz
HEREDOC
)

help_message_pmat_suppfigure3=$(
  cat <<HEREDOC

  get ${_local_host}:$PWD/<outdir>-a.tar.gz
HEREDOC
)

help_message_data_example=$(
  cat <<HEREDOC

  Create a data file with species name and SRA IDs.
HEREDOC
)

help_message_maintable1=$(
  cat <<HEREDOC

  Create the tables in the main text.
HEREDOC
)

help_message_supptable1=$(
  cat <<HEREDOC

  Create the tables in the main text.
  1. PMAT & nextDenovo: memory and time usage
HEREDOC
)

############################################################
# CSV setting for each analysis
#
declare -A _folder
declare -A _taxon
declare -A _dummy
declare -A _status
declare -A _host
declare -A _long
declare -A _short
declare -A _inref
declare -A _random_seed
declare -A _downsample

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
  while IFS=$',' read -r species long short random_seed downsample inref dummy status; do
    # Skip header line
    [[ "$species" == "species" ]] && continue
    [[ "$species" == \#* ]] && continue

    # Store in associative arrays
    if [[ -z "${species:-}" ]]; then
      continue
    fi
    _long["$species"]="$long"
    _short["$species"]="$short"
    _random_seed["$species"]="$random_seed"
    _downsample["$species"]="$downsample"
    _inref["$species"]="$inref"
    _dummy["$species"]="$dummy"
    _status["$species"]="$status"
  done <"$csv_file"

  # Create an array with unique values from the associative array
  # Sall=($(printf "%s\n" "${_folder[@]}" | sort -u))

  # Extract, clean, sort, and deduplicate keys
  mapfile -t Sall < <(
    for key in "${!_long[@]}"; do
      echo "${key%%-*}"
    done | sort -u
  )
}

read-a-tsv-file-into-associative-arrays

keys_array=($(for key in "${!_long[@]}"; do echo "$key"; done | sort))

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

# Common operations function
common_operations() {
  local long_sra="$1"
  local short_sra="$2"

  source ~/miniconda3/bin/activate polap
  if [[ -s "o1/0/mt.contig.name-1" ]]; then
    mkdir -p o/0
    cp -p o2/*.txt o
    cp -p o1/0/mt.contig.name-1 o/0/
    cp -pr o1/0/30-contigger o/0/
  fi
  if [[ ! -s "${long_sra}.fastq" ]]; then
    ${_polap_cmd} x-ncbi-fetch-sra --sra "$long_sra"
    ${_polap_cmd} x-ncbi-fetch-sra --sra "$short_sra"
    ${_polap_cmd} x-ncbi-fetch-sra-runinfo --sra "$long_sra"
    ${_polap_cmd} x-ncbi-fetch-sra-runinfo --sra "$short_sra"
  fi
  cp -s "${long_sra}.fastq" l.fq
  cp -s "${short_sra}_1.fastq" s1.fq
  cp -s "${short_sra}_2.fastq" s2.fq
  if [[ -s "o/0/mt.contig.name-1" ]]; then
    ${_polap_cmd} reduce-data
  fi
}

##### INSERT_FUNCTION_HERE #####
polap-batch_genus_species() {
  local _brg_outdir="${1}"
  local _brg_inum="${2:-x}"

  if [[ "${_brg_inum}" == "x" ]]; then
    for _brg_inum in 20 10; do
      data-long_genus_species "${_brg_outdir}"
      data-short_genus_species "${_brg_outdir}"

      polap-batch-assemble_genus_species "${_brg_outdir}" "${_brg_inum}"
      if [[ ! -d "${_brg_outdir}/msbwt" ]]; then
        polap-msbwt_genus_species "${_brg_outdir}"
      fi
      polap-batch-polish_genus_species "${_brg_outdir}" "${_brg_inum}"
      polap-archive_genus_species "${_brg_outdir}" "${_brg_inum}"
      polap-clean_genus_species "${_brg_outdir}" "${_brg_inum}"
    done
  else
    data-long_genus_species "${_brg_outdir}"
    data-short_genus_species "${_brg_outdir}"

    polap-batch-assemble_genus_species "${_brg_outdir}" "${_brg_inum}"
    if [[ ! -d "${_brg_outdir}/msbwt" ]]; then
      polap-msbwt_genus_species "${_brg_outdir}"
    fi
    polap-batch-polish_genus_species "${_brg_outdir}" "${_brg_inum}"
    polap-archive_genus_species "${_brg_outdir}" "${_brg_inum}"
    polap-clean_genus_species "${_brg_outdir}" "${_brg_inum}"
  fi
}

polap-batch-polish_genus_species() {
  local _brg_outdir="${1}"
  local _brg_inum="${2:-0}"
  local _brg_msbwt="${3:-x}"
  local target_index="${_brg_outdir}-${_brg_inum}"
  local downsample="${_downsample["$target_index"]}"

  local _outdir_main="${_brg_outdir}"
  rm -rf "${_brg_outdir}/msbwt"
  cp -pr "${_brg_outdir}/msbwts/x" "${_brg_outdir}/msbwt"

  # local _outdir_seeds="${_brg_outdir}"-"${_brg_inum}"-polap-seeds
  # local _i_seeds="${_outdir_seeds}/0"

  for j in {1..9}; do
    _brg_jnum=$((_brg_inum + j))
    local _i_main="${_outdir_main}/${_brg_jnum}"
    if [[ -s "${_i_main}/assembly.fasta" ]]; then
      polap-polish_genus_species "${_brg_outdir}" "${_brg_inum}" "${_brg_jnum}"
    fi
  done
}

# Batch steps
# 1. polap-batch <outdir> 30
# 2. polap-msbwt <outdir> 30
# 3. polap-archive <outdir> 30 31
# 4. polap-batch-polish <outdir> 30 31
polap-msbwt_genus_species() {
  local _brg_outdir="${1}"
  local _brg_inum="${2:-0}"
  _brg_inum=0
  # local _brg_coverage="${3}"

  local target_index="${_brg_outdir}-${_brg_inum}"
  local _brg_outdir_i="${_brg_outdir}/${_brg_inum}"
  local long_sra="${_long["$target_index"]}"
  local short_sra="${_short["$target_index"]}"
  local _outdir_msbwt="${_brg_outdir}"-"${_brg_inum}"-polap-msbwt
  local _outdir_assemble1="${_brg_outdir}"-"${_brg_inum}"-polap-assemble1
  mkdir -p "${_brg_outdir}/msbwts"
  mkdir -p "${_brg_outdir}/tmp"

  source "$(conda info --base)/etc/profile.d/conda.sh"
  if [[ "$CONDA_DEFAULT_ENV" != "polap-fmlrc" ]]; then
    echo "You're not in the polap-fmlrc environment. Chaniging 'polap-fmlrc'..."
    conda activate polap-fmlrc
  fi

  if [[ "$CONDA_DEFAULT_ENV" == "polap-fmlrc" ]]; then

    for k in x; do
      local _timing_polap_msbwt="${_brg_outdir_i}"/timing-polap-msbwt-k-${k}.txt
      local _stdout_polap_msbwt="${_brg_outdir_i}"/stdout-polap-msbwt-k-${k}.txt
      mkdir -p "${_brg_outdir_i}"
      mkdir -p "${_outdir_msbwt}"

      # subsample the short-read data
      if [[ "${k}" == "x" ]]; then
        # rm "${_brg_outdir}/tmp/${short_sra}_1.fastq"
        # rm "${_brg_outdir}/tmp/${short_sra}_2.fastq"
        ln -sf "$(pwd)/${short_sra}_1.fastq" \
          "${_brg_outdir}/tmp/"
        ln -sf "$(pwd)/${short_sra}_2.fastq" \
          "${_brg_outdir}/tmp/"
      else
        local genomesize=$(<"${_outdir_assemble1}/short_expected_genome_size.txt")
        echo data-downsample-short_genus_species "${_brg_outdir}" \
          "${_brg_inum}" "${k}" "${genomesize}"
        data-downsample-short_genus_species "${_brg_outdir}" \
          "${_brg_inum}" "${k}" "${genomesize}"
      fi

      # msbwt
      command time -v ${_polap_cmd} prepare-polishing \
        -a "${_brg_outdir}/tmp/${short_sra}_1.fastq" \
        -b "${_brg_outdir}/tmp/${short_sra}_2.fastq" \
        -o ${_outdir_msbwt} \
        >"${_stdout_polap_msbwt}" \
        2>"${_timing_polap_msbwt}"

      # Record the computer system info
      echo "hostname: $(hostname)" >>"${_timing_polap_msbwt}"
      free -h >>"${_timing_polap_msbwt}"
      lscpu >>"${_timing_polap_msbwt}"

      mv "${_outdir_msbwt}/msbwt" "${_brg_outdir}/msbwts/${k}"
    done

    # local _pmat_dir="${_brg_outdir_i}"/pmat-${_brg_fc}
    # rm -rf "${_pmat_dir}"
    # mkdir "${_pmat_dir}"
    # cp -pr "${_brg_outdir}"-pmat-${_brg_fc}/gfa_result "${_pmat_dir}"
    # rm -rf "${_brg_outdir}"-pmat-${_brg_fc}

    # results
    # echo "results"
    # echo "output: ${_brg_outdir_i}/pmat.out"
    # echo "timing output: ${_brg_outdir_i}/timing-pmat.txt"
    # echo "${_brg_outdir_i}"/pmat/gfa_result/PMAT_mt_master.gfa
    # echo "${_brg_outdir_i}"/pmat/gfa_result/PMAT_pt_master.gfa

    conda deactivate

    # archive it
    # send it back if it is remote
    # if [[ "${_local_host}" == "$(hostname)" ]]; then
    # 	echo "no operation for the local"
    # else
    # 	archive_genus_species ${_brg_outdir}
    # 	scp -p ${_brg_outdir}-a.tar.gz ${_local_host}:$PWD/${_brg_outdir}-a-${_brg_inum}.tar.gz
    # fi

  else
    echo "ERROR: no polap conda environment"
  fi

}

polap-batch-assemble_genus_species() {
  local _brg_outdir="${1}"
  local _brg_inum="${2:-0}"
  local target_index="${_brg_outdir}-${_brg_inum}"
  local downsample="${_downsample["$target_index"]}"

  # prepare the input data
  # data-long_genus_species "${_brg_outdir}"
  # data-short_genus_species "${_brg_outdir}"

  # Download data
  polap-assemble1_genus_species "${_brg_outdir}" "${_brg_inum}" "${downsample}"
  polap-seeds_genus_species "${_brg_outdir}" "${_brg_inum}"

  # polap-msbwt_genus_species "${_brg_outdir}" "${_brg_inum}"

  local _outdir_seeds="${_brg_outdir}"-"${_brg_inum}"-polap-seeds
  local _i_seeds="${_outdir_seeds}/0"
  for j in {1..9}; do
    _brg_jnum=$((_brg_inum + j))
    if [[ -s "${_i_seeds}/mt.contig.name-${j}" ]]; then
      mv "${_i_seeds}/mt.contig.name-${j}" \
        "${_i_seeds}/mt.contig.name-${_brg_jnum}"
      polap-assemble2_genus_species "${_brg_outdir}" "${_brg_inum}" "${_brg_jnum}"
      # polap-polish_genus_species "${_brg_outdir}" "${_brg_inum}" "${_brg_jnum}"
    fi
  done
}

polap-archive_genus_species() {
  local _brg_outdir="${1}"
  local _brg_inum="${2:-0}"
  local _brg_jnum="${3:-1}"

  local target_index="${_brg_outdir}-${_brg_inum}"
  local _brg_outdir_i="${_brg_outdir}/${_brg_inum}"
  local _outdir_main="${_brg_outdir}"
  local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

  # Copy necessary files
  local _outdir_assemble1="${_brg_outdir}"-"${_brg_inum}"-polap-assemble1
  local _outdir_seeds="${_brg_outdir}"-"${_brg_inum}"-polap-seeds
  local _outdir_assemble2="${_brg_outdir}"-"${_brg_inum}"-polap-assemble2
  local _outdir_msbwt="${_brg_outdir}"-0-polap-msbwt
  local _i_main="${_outdir_main}/${_brg_inum}"
  local _i_assemble1="${_outdir_assemble1}/0"
  local _i_seeds="${_outdir_seeds}/0"
  local _i_assemble2="${_outdir_assemble2}/0"
  local _30_contigger_assemble1="${_i_assemble1}/30-contigger"
  local _30_contigger_assemble2="${_i_assemble1}/30-contigger"
  local _30_contigger_main="${_i_main}/30-contigger"

  mkdir -p "${_30_contigger_main}"
  cp -p "${_i_assemble1}"/assembly_info_organelle_annotation_count-all.txt \
    "${_i_main}"
  cp -p "${_i_assemble2}"/mt.contig.name-* "${_i_main}"
  # for i in graph_final.fasta graph_final.gfa edges_stats.txt; do
  for i in graph_final.gfa edges_stats.txt; do
    cp -p "${_30_contigger_assemble2}"/$i "${_30_contigger_main}/"
  done
  cp -p "${_outdir_assemble1}"/*.stats "${_i_main}"
  cp -p "${_outdir_assemble1}"/*.txt "${_i_main}"
  cp -p "${_outdir_seeds}"/polap.log "${_i_main}"/polap-seeds.log
  cp -p "${_outdir_assemble1}"/polap.log "${_i_main}"/polap-assemble1.log
  cp -p "${_outdir_assemble2}"/polap.log "${_i_main}"/polap-assemble2.log
  cp -p "${_outdir_msbwt}"/polap.log "${_i_main}"/polap-msbwt.log
  cp -p "${_i_assemble1}"/*.txt "${_i_main}"

  for j in {1..9}; do
    _brg_jnum=$((_brg_inum + j))
    if [[ -s "${_i_seeds}/mt.contig.name-${_brg_jnum}" ]]; then
      local _j_assemble2="${_outdir_assemble2}/${_brg_jnum}"
      local _j_main="${_outdir_main}/${_brg_jnum}"
      mkdir -p "${_j_main}"
      for i in assembly_graph.gfa assembly_info.txt assembly.fasta; do
        cp -p "${_j_assemble2}"/$i "${_j_main}/"
      done
    fi
  done
}

polap-polish_genus_species() {
  local _brg_outdir="${1}"
  local _brg_inum="${2:-0}"
  local _brg_jnum="${3:-1}"

  local target_index="${_brg_outdir}-${_brg_inum}"
  local _brg_outdir_i="${_brg_outdir}/${_brg_inum}"
  local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

  # Copy necessary files
  local _outdir_main="${_brg_outdir}"
  # local _outdir_assemble1="${_brg_outdir}"-"${_brg_inum}"-polap-assemble1
  # local _outdir_assemble2="${_brg_outdir}"-"${_brg_inum}"-polap-assemble2

  # if [[ -d "${_outdir_assemble1}"/msbwt ]]; then
  # 	mv "${_outdir_assemble1}"/msbwt "${_outdir_assemble2}"
  # fi

  source "$(conda info --base)/etc/profile.d/conda.sh"
  if [[ "$CONDA_DEFAULT_ENV" != "polap" ]]; then
    echo "You're not in the polap environment. Chaniging 'polap'..."
    conda activate polap
  fi

  if [[ "$CONDA_DEFAULT_ENV" == "polap" ]]; then

    local _timing_polap_polish="${_brg_outdir_i}"/timing-polap-polish.txt
    local _stdout_polap_polish="${_brg_outdir_i}"/stdout-polap-polish.txt
    mkdir -p "${_brg_outdir_i}"

    command time -v "${_polap_cmd}" polish \
      -o "${_outdir_main}" \
      -p "${_outdir_main}/${_brg_jnum}/assembly.fasta" \
      -f "${_outdir_main}/${_brg_jnum}/assembly-polished.fa" \
      >"${_stdout_polap_polish}" \
      2>"${_timing_polap_polish}"

    # Record the computer system info
    echo "hostname: $(hostname)" >>"${_timing_polap_polish}"
    free -h >>"${_timing_polap_polish}"
    lscpu >>"${_timing_polap_polish}"

    # local _pmat_dir="${_brg_outdir_i}"/pmat-${_brg_fc}
    # rm -rf "${_pmat_dir}"
    # mkdir "${_pmat_dir}"
    # cp -pr "${_brg_outdir}"-pmat-${_brg_fc}/gfa_result "${_pmat_dir}"
    # rm -rf "${_brg_outdir}"-pmat-${_brg_fc}

    # results
    # echo "results"
    # echo "output: ${_brg_outdir_i}/pmat.out"
    # echo "timing output: ${_brg_outdir_i}/timing-pmat.txt"
    # echo "${_brg_outdir_i}"/pmat/gfa_result/PMAT_mt_master.gfa
    # echo "${_brg_outdir_i}"/pmat/gfa_result/PMAT_pt_master.gfa

    conda deactivate

    # archive it
    # send it back if it is remote
    # if [[ "${_local_host}" == "$(hostname)" ]]; then
    # 	echo "no operation for the local"
    # else
    # 	archive_genus_species ${_brg_outdir}
    # 	scp -p ${_brg_outdir}-a.tar.gz ${_local_host}:$PWD/${_brg_outdir}-a-${_brg_inum}.tar.gz
    # fi

  else
    echo "ERROR: no polap conda environment"
  fi

}

polap-assemble2_genus_species() {
  local _brg_outdir="${1}"
  local _brg_inum="${2:-0}"
  local _brg_jnum="${3:-1}"

  local target_index="${_brg_outdir}-${_brg_inum}"
  local _brg_outdir_i="${_brg_outdir}/${_brg_inum}"
  local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

  # Copy necessary files
  local _outdir_assemble1="${_brg_outdir}"-"${_brg_inum}"-polap-assemble1
  local _outdir_seeds="${_brg_outdir}"-"${_brg_inum}"-polap-seeds
  local _outdir_assemble2="${_brg_outdir}"-"${_brg_inum}"-polap-assemble2
  local _i_seeds="${_outdir_seeds}/0"
  local _i_assemble2="${_outdir_assemble2}/0"
  local _30_contigger_assemble1="${_outdir_assemble1}/0/30-contigger"
  local _30_contigger_assemble2="${_outdir_assemble2}/0/30-contigger"

  mkdir -p "${_30_contigger_assemble2}"
  cp "${_outdir_assemble1}"/lk.fq.gz "${_outdir_assemble2}"
  cp "${_i_seeds}"/mt.contig.name-* "${_i_assemble2}"
  # for i in graph_final.fasta graph_final.gfa edges_stats.txt; do
  for i in graph_final.gfa edges_stats.txt; do
    cp "${_30_contigger_assemble1}"/$i "${_30_contigger_assemble2}/"
  done

  source "$(conda info --base)/etc/profile.d/conda.sh"
  if [[ "$CONDA_DEFAULT_ENV" != "polap" ]]; then
    echo "You're not in the polap environment. Chaniging 'polap'..."
    conda activate polap
  fi

  if [[ "$CONDA_DEFAULT_ENV" == "polap" ]]; then

    local _timing_polap_assemble2="${_brg_outdir_i}"/timing-polap-assemble2.txt
    local _stdout_polap_assemble2="${_brg_outdir_i}"/stdout-polap-assemble2.txt
    mkdir -p "${_brg_outdir_i}"

    command time -v "${_polap_cmd}" assemble2 \
      -o "${_outdir_assemble2}" \
      -j "${_brg_jnum}" \
      >>"${_stdout_polap_assemble2}" \
      2>>"${_timing_polap_assemble2}"

    # Record the computer system info
    echo "hostname: $(hostname)" >>"${_timing_polap_assemble2}"
    free -h >>"${_timing_polap_assemble2}"
    lscpu >>"${_timing_polap_assemble2}"

    # local _pmat_dir="${_brg_outdir_i}"/pmat-${_brg_fc}
    # rm -rf "${_pmat_dir}"
    # mkdir "${_pmat_dir}"
    # cp -pr "${_brg_outdir}"-pmat-${_brg_fc}/gfa_result "${_pmat_dir}"
    # rm -rf "${_brg_outdir}"-pmat-${_brg_fc}

    # results
    # echo "results"
    # echo "output: ${_brg_outdir_i}/pmat.out"
    # echo "timing output: ${_brg_outdir_i}/timing-pmat.txt"
    # echo "${_brg_outdir_i}"/pmat/gfa_result/PMAT_mt_master.gfa
    # echo "${_brg_outdir_i}"/pmat/gfa_result/PMAT_pt_master.gfa

    conda deactivate

    # archive it
    # send it back if it is remote
    # if [[ "${_local_host}" == "$(hostname)" ]]; then
    # 	echo "no operation for the local"
    # else
    # 	archive_genus_species ${_brg_outdir}
    # 	scp -p ${_brg_outdir}-a.tar.gz ${_local_host}:$PWD/${_brg_outdir}-a-${_brg_inum}.tar.gz
    # fi

  else
    echo "ERROR: no polap conda environment"
  fi

}

polap-seeds_genus_species() {
  local _brg_outdir="${1}"
  local _brg_inum="${2:-0}"

  local target_index="${_brg_outdir}-${_brg_inum}"
  local _brg_outdir_i="${_brg_outdir}/${_brg_inum}"
  local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

  # Copy necessary files
  local _outdir_assemble1="${_brg_outdir}"-"${_brg_inum}"-polap-assemble1
  local _outdir_seeds="${_brg_outdir}"-"${_brg_inum}"-polap-seeds
  local _i_assemble1="${_outdir_assemble1}/0"
  local _i_seeds="${_outdir_seeds}/0"
  local _30_contigger_assemble1="${_outdir_assemble1}/0/30-contigger"
  local _30_contigger_seeds="${_outdir_seeds}/0/30-contigger"
  mkdir -p "${_30_contigger_assemble1}"
  mkdir -p "${_30_contigger_seeds}"
  cp "${_i_assemble1}"/assembly_info_organelle_annotation_count-all.txt \
    "${_i_seeds}"
  # for i in graph_final.fasta graph_final.gfa edges_stats.txt; do
  for i in graph_final.gfa edges_stats.txt; do
    cp "${_30_contigger_assemble1}"/$i "${_30_contigger_seeds}/"
  done

  source "$(conda info --base)/etc/profile.d/conda.sh"
  if [[ "$CONDA_DEFAULT_ENV" != "polap" ]]; then
    echo "You're not in the polap environment. Chaniging 'polap'..."
    conda activate polap
  fi

  if [[ "$CONDA_DEFAULT_ENV" == "polap" ]]; then

    local _timing_polap_seeds="${_brg_outdir_i}"/timing-polap-seeds.txt
    local _stdout_polap_seeds="${_brg_outdir_i}"/stdout-polap-seeds.txt
    mkdir -p "${_brg_outdir_i}"

    command time -v "${_polap_cmd}" seeds \
      -o "${_outdir_seeds}" \
      >>"${_stdout_polap_seeds}" \
      2>>"${_timing_polap_seeds}"

    # Record the computer system info
    echo "hostname: $(hostname)" >>"${_timing_polap_seeds}"
    free -h >>"${_timing_polap_seeds}"
    lscpu >>"${_timing_polap_seeds}"

    # local _pmat_dir="${_brg_outdir_i}"/pmat-${_brg_fc}
    # rm -rf "${_pmat_dir}"
    # mkdir "${_pmat_dir}"
    # cp -pr "${_brg_outdir}"-pmat-${_brg_fc}/gfa_result "${_pmat_dir}"
    # rm -rf "${_brg_outdir}"-pmat-${_brg_fc}

    # results
    # echo "results"
    # echo "output: ${_brg_outdir_i}/pmat.out"
    # echo "timing output: ${_brg_outdir_i}/timing-pmat.txt"
    # echo "${_brg_outdir_i}"/pmat/gfa_result/PMAT_mt_master.gfa
    # echo "${_brg_outdir_i}"/pmat/gfa_result/PMAT_pt_master.gfa

    conda deactivate

    # archive it
    # send it back if it is remote
    # if [[ "${_local_host}" == "$(hostname)" ]]; then
    # 	echo "no operation for the local"
    # else
    # 	archive_genus_species ${_brg_outdir}
    # 	scp -p ${_brg_outdir}-a.tar.gz ${_local_host}:$PWD/${_brg_outdir}-a-${_brg_inum}.tar.gz
    # fi

  else
    echo "ERROR: no polap conda environment"
  fi

}

polap-clean_genus_species() {
  local _brg_outdir="${1}"
  local _brg_inum="${2:-0}"

  local target_index="${_brg_outdir}-${_brg_inum}"
  local _brg_outdir_i="${_brg_outdir}/${_brg_inum}"
  local long_sra="${_long["$target_index"]}"
  if [[ -z "$long_sra" ]]; then
    echo "Error: skipping ${_brg_outdir}-${_brg_inum} because it is not in the CSV."
    return
  fi
  local short_sra="${_short["$target_index"]}"

  if [[ -n "${_brg_outdir}" ]]; then
    rm -rf "${_brg_outdir}"-"${_brg_inum}"*
    rm -f "${long_sra}".fastq
    rm -f "${short_sra}"_?.fastq
  fi
}

polap-assemble1_genus_species() {
  local _brg_outdir="${1}"
  local _brg_inum="${2:-0}"
  local _brg_coverage="${3}"

  local target_index="${_brg_outdir}-${_brg_inum}"
  local long_sra="${_long["$target_index"]}"
  local short_sra="${_short["$target_index"]}"
  if [[ -z "${_brg_coverage}" ]]; then
    _brg_coverage="${_downsample["$target_index"]}"
  fi
  local _brg_outdir_i="${_brg_outdir}/t1/${_brg_inum}"
  local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

  # Step 0. Get the archive data file if this runs at a remote
  # where we do not have the file.
  # data-short_genus_species "${_brg_outdir}" "${_brg_inum}"
  # data-long_genus_species "${_brg_outdir}" "${_brg_inum}"

  source "$(conda info --base)/etc/profile.d/conda.sh"
  if [[ "$CONDA_DEFAULT_ENV" != "polap" ]]; then
    echo "You're not in the polap environment. Chaniging 'polap'..."
    conda activate polap
  fi

  if [[ "$CONDA_DEFAULT_ENV" == "polap" ]]; then

    local _timing_polap_assemble1="${_brg_outdir_i}"/timing-polap-assemble1.txt
    local _stdout_polap_assemble1="${_brg_outdir_i}"/stdout-polap-assemble1.txt
    local _outdir="${_brg_outdir}"-"${_brg_inum}"-polap-assemble1
    mkdir -p "${_brg_outdir_i}"

    command time -v "${_polap_cmd}" assemble1 \
      -l "${long_sra}".fastq \
      -a "${short_sra}"_1.fastq \
      -b "${short_sra}"_2.fastq \
      -c "${_brg_coverage}" \
      -o "${_outdir}" \
      >>"${_stdout_polap_assemble1}" \
      2>>"${_timing_polap_assemble1}"

    # Record the computer system info
    echo "hostname: $(hostname)" >>"${_timing_polap_assemble1}"
    free -h >>"${_timing_polap_assemble1}"
    lscpu >>"${_timing_polap_assemble1}"

    # local _pmat_dir="${_brg_outdir_i}"/pmat-${_brg_fc}
    # rm -rf "${_pmat_dir}"
    # mkdir "${_pmat_dir}"
    # cp -pr "${_brg_outdir}"-pmat-${_brg_fc}/gfa_result "${_pmat_dir}"
    # rm -rf "${_brg_outdir}"-pmat-${_brg_fc}

    # results
    # echo "results"
    # echo "output: ${_brg_outdir_i}/pmat.out"
    # echo "timing output: ${_brg_outdir_i}/timing-pmat.txt"
    # echo "${_brg_outdir_i}"/pmat/gfa_result/PMAT_mt_master.gfa
    # echo "${_brg_outdir_i}"/pmat/gfa_result/PMAT_pt_master.gfa

    conda deactivate

    # archive it
    # send it back if it is remote
    # if [[ "${_local_host}" == "$(hostname)" ]]; then
    # 	echo "no operation for the local"
    # else
    # 	archive_genus_species ${_brg_outdir}
    # 	scp -p ${_brg_outdir}-a.tar.gz ${_local_host}:$PWD/${_brg_outdir}-a-${_brg_inum}.tar.gz
    # fi

  else
    echo "ERROR: no polap conda environment"
  fi

}

get-data_genus_species() {
  local _brg_outdir="${1}"
  echo "Preparing archive for ${_brg_outdir} ..."
  # tar zcf "${_brg_outdir}-a.tar.gz" "${_brg_outdir}"
}

data2csv_genus_species() {
  local _brg_data="${1:-${_polap_data_data}}"
  local _brg_csv="${2:-${_polap_data_csv}}"
  local _brg_setting_I="${3:-0}"
  local _brg_setting_S="${4:-0}"
  local _brg_setting_D="${5:-0}"

  # Print header
  echo "species,long,short,random_seed,downsample,dummy,status" >"${_brg_csv}"

  # Read input file line by line
  while IFS=',' read -r species long short; do
    [[ "$species" == "species" ]] && continue

    local base="${species}"
    local species="${species}-${_brg_setting_I}"
    echo "${species},${long},${short},${_brg_setting_S},dummy,done" >>"${_brg_csv}"
  done <"${_brg_data}"

  echo "create ${_brg_csv}"
}
# Compare polap, ptGAUL, GetOrganelle, and PMAT
pmat-suppfigure3_genus_species() {
  local _brg_inum="${1:-0}"
  local _brg_bandage="${2:-no}"
  local _brg_t_dir="${3:-"${_brg_default_target_dir}"}"

  # local _table_name="suppfigure1"
  # if ((_brg_table == 2)); then
  # 	_table_name="suppfigure2"
  # fi
  # local _supptable_md="${_table_name}-${_brg_inum}-${_brg_stage}.md"

  local _tf_mainfigure="F"
  local _supfigure_file="pmat-suppfigure3-${_brg_inum}.md"

  # printf "# Supplementary Figures: Polap, ptGAUL, and GetOrganelle\n\n" \
  # 	>"${_supfigure_file}"
  rm -f "${_supfigure_file}"

  # Start writing to the markdown file
  cat <<EOF >>"${_supfigure_file}"

  | Species (Order) | Polap | PMAT |
|-----------------|-----------------|-----------------|
EOF

  # Extract and sort keys
  sorted_keys=($(for key in "${!_long[@]}"; do echo "$key"; done | sort))

  # Iterate over sorted keys and check if value is "T"
  for _brg_outdir in "${sorted_keys[@]}"; do
    if [[ "${_brg_outdir}" == "Carex_pseudochinensis" ]]; then
      continue
    fi
    key="${_brg_outdir}"
    echo Processing $key ...

    # local _species=${_taxon[$key]}
    local _species="${key//_/ }"
    local _genus=${_species%% *}
    local _order=$(grep "${_genus}" ${_POLAPLIB_DIR}/taxonomy_output.tsv | cut -f 6 | head -1)
    _brg_outdir="$key"

    local _gfa_polap="${_brg_outdir}/o2/1/assembly_graph.gfa"
    local _png_polap="${_brg_outdir}/o2/1/assembly_graph.png"

    local images=()
    local captions=()
    images+=("${_species}")
    captions+=("${_order}")
    images+=(${_png_polap})
    captions+=("Polap")

    local v=${pmat_fc[${_brg_outdir}]}
    local _gfa_pmat="${_brg_outdir}/0/pmat-${v}/gfa_result/PMAT_mt_master.gfa"
    local _png_pmat="${_brg_outdir}/0/pmat-${v}/gfa_result/PMAT_mt_master.png"
    if [[ "${_brg_bandage}" == "yes" ]]; then
      ${_polap_cmd} bandage png \
        ${_gfa_pmat} \
        ${_png_pmat}
    fi
    images+=(${_png_pmat})
    captions+=("PMAT")

    output="$_supfigure_file"

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
      if ((count % 3 == 0)); then
        echo "$row" >>"$output"
        echo "$caption_row" >>"$output"
        echo "|-----------------|-----------------|-----------------|" >>"$output"
        # echo "" >>"$output"
        row="| "
        caption_row="| "
      fi
    done

  done

  # debug local variables
  # for debugging: Inline local printing local var
  # while IFS= read -r line; do
  #   if [[ $line =~ ^declare\ --\ ([^=]+)= ]]; then
  #     var="${BASH_REMATCH[1]}"
  #     printf "%s=%q\n" "$var" "${!var}"
  #   fi
  # done < <(local -p 2>/dev/null)
  # return

  echo ${_supfigure_file}
  cp -p ${_supfigure_file} "${_brg_t_dir}"
  echo cp -p ${_supfigure_file} "${_brg_t_dir}"
}

copy-figures_genus_species() {
  local _brg_t_dir_figures="${1:-"${_brg_default_target_dir}"}"
  rsync -qav --include='*/' --include='*.png' --exclude='*' ./ "${_brg_t_dir_figures}/"
  rsync -qav --include='*/' --include='*.pdf' --exclude='*' ./ "${_brg_t_dir_figures}/"
}

maintable1_genus_species_header() {
  local _brg_table="${1:-1}"

  local _items

  if ((_brg_table == 1)); then
    local _items=(
      "Species"
      "Order"
      "Family"
      "L_size"
      "L_cov"
      "S_size"
      "S_cov"
      "G"
      "D"
      "Mpa0"
      "Tpa0"
      "Mpa1"
      "Tpa1"
      "Cpa2"
      "Mpa2"
      "Tpa2"
      "Mpp1"
      "Tpp1"
      "Mpp2"
      "Tpp2"
      "Mpolap1"
      "Tpolap1"
      "Mpolap2"
      "Tpolap2"
      "Mpolap"
      "Tpolap"
      "Mndp"
      "Tndp"
      "Mpmata"
      "Tpmata"
      "Mpmat"
      "Tpmat"
    )
  elif ((_brg_table == 2)); then
    local _items=(
      "Species"
      "Order"
      "Family"
      "L_size"
      "L_cov"
      "S_size"
      "S_cov"
      "G"
      "D"
      "Mpa0"
      "Tpa0"
      "Mpa1"
      "Tpa1"
      "Cpa2"
      "Mpa2"
      "Tpa2"
      "Mpp1"
      "Tpp1"
      "Mpp2"
      "Tpp2"
      "Mpolap1"
      "Tpolap1"
      "Mpolap2"
      "Tpolap2"
      "Mpolap"
      "Tpolap"
      "Mndp"
      "Tndp"
      "Mpmata"
      "Tpmata"
      "Mpmat"
      "Tpmat"
    )
  elif ((_brg_table == 3)); then
    _items=(
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
      "D"
      "A0"
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
      "T_g"
      "M_p"
      "T_p"
      "M_pmat1"
      "T_pmat1"
      "M_pmat2"
      "T_pmat2"
      "M_t1"
      "M_t2"
      "M_s"
      "M_f"
      "T"
    )
  fi

  printf "%s\t" "${_items[@]::${#_items[@]}-1}"
  printf "%s\n" "${_items[-1]}"
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

# Extract long-read SRA from polap.log
# extract_long_sra() {
# 	local _polap_log="$1"
# 	local _l_sra
# 	_l_sra=$(awk '/long-read/ {match($0, /([A-Z]RR[0-9]+)/, arr); print arr[0]}' "$_polap_log" | sort -u | grep -v '^$')
# 	echo "$_l_sra"
# }
#
extract_long_sra_seqkit_stats() {
  local _fq_stats="$1"

  # Check if input file is readable
  if [[ ! -r "$_fq_stats" ]]; then
    echo "Error: Log file '$_fq_stats' does not exist or is not readable" >&2
    return 1
  fi

  local _l_sra=$(awk 'NR==2 { sub(/\.fastq$/, "", $1); print $1 }' "${_fq_stats}")
  _l_sra=${_l_sra%_1}

  echo "$_l_sra"
}

# extract_long_sra_size_gb() {
# 	local _v1="$1"
# 	local _l_sra_size
#
# 	if [[ -s "${_v1}/long_total_length.txt" ]]; then
# 		_l_sra_size=$(<"${_v1}/long_total_length.txt")
# 	else
# 		_l_sra_size=$(<"${_v1}/l.fq.txt")
# 	fi
#
# 	local _l_sra_size_gb
# 	_l_sra_size_gb=$(_polap_utility_convert_bp "${_l_sra_size}")
# 	echo "$_l_sra_size_gb"
# }
extract_long_sra_size_gb() {
  local _v1="$1"

  # Check that input directory is provided
  if [[ -z "$_v1" ]]; then
    echo "Error: No directory provided to extract_long_sra_size_gb" >&2
    return 1
  fi

  local _long_file="${_v1}/long_total_length.txt"
  local _fallback_file="${_v1}/l.fq.txt"
  local _l_sra_size

  # Check for long_total_length.txt or fallback to l.fq.txt
  if [[ -s "$_long_file" ]]; then
    _l_sra_size=$(<"$_long_file")
  elif [[ -s "$_fallback_file" ]]; then
    _l_sra_size=$(<"$_fallback_file")
  else
    echo "Error: Neither '$_long_file' nor '$_fallback_file' is available or non-empty" >&2
    return 1
  fi

  # Check that value is numeric
  if ! [[ "$_l_sra_size" =~ ^[0-9]+$ ]]; then
    echo "Error: SRA size value is not a valid integer: $_l_sra_size" >&2
    return 1
  fi

  # Ensure the conversion function exists
  if ! command -v _polap_utility_convert_bp &>/dev/null; then
    echo "Error: Conversion function '_polap_utility_convert_bp' not found" >&2
    return 1
  fi

  # Attempt the conversion
  local _l_sra_size_gb
  if ! _l_sra_size_gb=$(_polap_utility_convert_bp "${_l_sra_size}"); then
    echo "Error: _polap_utility_convert_bp failed for input $_l_sra_size" >&2
    return 1
  fi

  echo "$_l_sra_size_gb"
}

extract_seqkit_stats_sum_len() {
  awk 'NR==2 { print $5 }' $1
}

extract_short_sra() {
  local _polap_log="$1"
  local _s_sra
  _s_sra=$(awk '/short-read1/ {match($0, /([A-Z]RR[0-9]+)/, arr); print arr[0]}' "$_polap_log" | sort -u | grep -v '^$')
  echo "$_s_sra"
}

extract_short_sra_size_gb() {
  local _v1="$1"
  local _s_sra_size

  if [[ -s "${_v1}/short_total_length.txt" ]]; then
    _s_sra_size=$(<"${_v1}/short_total_length.txt")
  else
    local _s_sra_size1
    local _s_sra_size2
    _s_sra_size1=$(awk 'NR==2 { print $5 }' "${_v1}"/s1.fq.stats)
    _s_sra_size2=$(awk 'NR==2 { print $5 }' "${_v1}"/s2.fq.stats)
    # _s_sra_size1=$(<"${_v1}/s1.fq.stats")
    # _s_sra_size2=$(<"${_v1}/s2.fq.stats")
    _s_sra_size=$((_s_sra_size1 + _s_sra_size2))
  fi

  local _s_sra_size_gb
  _s_sra_size_gb=$(_polap_utility_convert_bp "${_s_sra_size}")
  echo "$_s_sra_size_gb"
}

extract_short_sra_size() {
  local _v1="$1"
  local _s_sra_size

  if [[ -s "${_v1}/short_total_length.txt" ]]; then
    _s_sra_size=$(<"${_v1}/short_total_length.txt")
  else
    local _s_sra_size1
    local _s_sra_size2
    _s_sra_size1=$(awk 'NR==2 { print $5 }' "${_v1}"/s1.fq.stats)
    _s_sra_size2=$(awk 'NR==2 { print $5 }' "${_v1}"/s2.fq.stats)
    # _s_sra_size1=$(<"${_v1}/s1.fq.stats")
    # _s_sra_size2=$(<"${_v1}/s2.fq.stats")
    _s_sra_size=$((_s_sra_size1 + _s_sra_size2))
  fi

  echo "$_s_sra_size"
}

extract_known_mtdna_accession() {
  local _polap_log="$1"

  # Check if the log file exists and is readable
  if [[ ! -r "$_polap_log" ]]; then
    echo "Error: Cannot read log file '$_polap_log'" >&2
    return 1
  fi

  local _known_mtdna
  _known_mtdna=$(grep 'NCBI accession: ' "$_polap_log" | cut -d: -f4 | tail -n 1 | xargs || true)
  _known_mtdna=${_known_mtdna:-'NA'}

  echo "$_known_mtdna"
}

extract_genome_size() {
  local _brg_outdir="$1"
  local _genome_size="NA"

  if [[ ! -d "$_brg_outdir" ]]; then
    echo "$_genome_size"
    return
  fi

  if [[ -s "${_brg_outdir}/o/short_expected_genome_size.txt" ]]; then
    _genome_size=$(<"${_brg_outdir}/o/short_expected_genome_size.txt")
  elif [[ -s "${_brg_outdir}/short_expected_genome_size.txt" ]]; then
    _genome_size=$(<"${_brg_outdir}/short_expected_genome_size.txt")
  else
    echo "$_genome_size"
    return
  fi

  _genome_size=${_genome_size%%.*}
  echo "$_genome_size"
}

extract_fasta_seqlen_plain() {
  local fasta="$1"

  if [[ ! -f "$fasta" ]]; then
    echo -e "NA\t0"
    return 0
  fi

  awk '
	/^>/ {
		if (seq != "") {
			print header "\t" length(seq)
		}
		header = substr($0, 2)
		seq = ""
		next
	}
	{
		seq = seq $0
	}
	END {
		if (seq != "") {
			print header "\t" length(seq)
		}
	}
	' "$fasta" 2>/dev/null || echo -e "NA\t0"

  return 0
}

extract_fasta_seqlen() {
  local fasta="$1"

  # If file does not exist, print length 0 and return
  if [[ ! -f "$fasta" ]]; then
    echo "NA"
    return 0
  fi

  # If bioawk is available, use it
  if command -v bioawk >/dev/null 2>&1; then
    bioawk -c fastx 'NR==1 {print length($seq); exit}' "$fasta" 2>/dev/null || echo "0"
    return 0
  fi

  # If seqkit is available and succeeds
  if command -v seqkit >/dev/null 2>&1; then
    if seqkit fx2tab -l -n "$fasta" | cut -f2 | xargs 2>/dev/null; then
      return 0
    fi
  fi

  # Fallback: extract headers and print 0 length
  extract_fasta_seqlen_plain "$fasta"
  return 0
}

extract_fasta_numseq() {
  local fasta="$1"

  # Check if file exists
  if [[ ! -f "$fasta" ]]; then
    echo 0
    return 0
  fi

  # Try seqkit if available
  if command -v seqkit >/dev/null 2>&1; then
    seqkit stats "$fasta" 2>/dev/null | awk 'NR > 1 {print $4}' || echo 0
    return 0
  fi

  # Fallback: count '>' lines with grep
  grep -c '^>' "$fasta" 2>/dev/null || echo 0
  return 0
}

extract_time_memory_timing_summary_file() {
  local _v1="$1"
  local _timing_file="${_v1}"

  _memory_gb="NA"
  _total_hours="NA"

  if [[ ! -s "$_timing_file" ]]; then
    echo "Error: timing file not found or empty: $_timing_file" >&2
    return 0
  fi

  _memory_gb=$(grep 'Net increase' "${_summary_file}" | awk -F'[()]' '{print $2}' | sed 's/ GB//')
  _total_hours=$(grep 'Elapsed time' "${_summary_file}" | awk '{print $3}')
  _total_hours=$(_polap_lib_timing-convert_to_hours_or_minutes "${_total_hours}")

  return 0
}

extract_time_memory_timing() {
  local _v1="$1"
  local _timing_file="${_v1}"

  _memory_gb="NA"
  _total_hours="NA"

  if [[ ! -s "$_timing_file" ]]; then
    echo "Error: timing file not found or empty: $_timing_file" >&2
    return 0
  fi

  local _output
  _output=$(_polap_lib_timing-parse-timing "$_timing_file" | tail -1 2>/dev/null)

  # validate output: must contain two fields
  if [[ "$_output" =~ ^[^[:space:]]+[[:space:]]+[^[:space:]]+ ]]; then
    read -r _memory_gb _total_hours <<<"$_output"
  else
    echo "Warning: malformed output from _polap_lib_timing-parse-timing" >&2
    echo "  Got: $_output" >&2
  fi

  return 0
}

extract_cumulative_time_memory_timing() {
  local _v1="$1"
  local _timing_file="${_v1}"

  _memory_gb="NA"
  _total_hours="NA"
  _count_run="NA"

  if [[ ! -s "$_timing_file" ]]; then
    echo "Error: timing file not found or empty: $_timing_file" >&2
    return 0
  fi

  local _output

  _output=$(_polap_lib_timing-parse-multiple-time-only "$_timing_file")
  local _cumulative_time=$(_polap_lib_timing-parse-cumulative-timing "$_timing_file")

  # validate output: must contain two fields
  if [[ "$_output" =~ ^[^[:space:]]+[[:space:]]+[^[:space:]]+ ]]; then
    read -r _memory_gb _count_run <<<"$_output"
  else
    echo "Warning: malformed output from _polap_lib_timing-parse-timing" >&2
    echo "  Got: $_output" >&2
  fi

  _total_hours="${_cumulative_time}"

  return 0
}

extract_time_memory_pmat1() {
  local _v1="$1"
  local _timing_file="${_v1}/timing-nextdenovo.txt"

  _memory_gb_pmat1="NA"
  _total_hours_pmat1="NA"

  if [[ ! -s "$_timing_file" ]]; then
    echo "Error: timing file not found or empty: $_timing_file" >&2
    return 0
  fi

  local _output
  _output=$(_polap_lib_timing-parse-timing "$_timing_file" 2>/dev/null)

  # Accept any two non-whitespace strings
  if [[ "$_output" =~ ^[^[:space:]]+[[:space:]]+[^[:space:]]+ ]]; then
    read -r _memory_gb_pmat1 _total_hours_pmat1 <<<"$_output"
  else
    echo "Warning: malformed output from _polap_lib_timing-parse-timing" >&2
    echo "  Got: $_output" >&2
  fi

  return 0
}

extract_time_memory_pmat2() {
  local _v1="$1"
  local _timing_file="${_v1}/timing-pmat-nextdenovo-fc-0.1.txt"

  _memory_gb_pmat2="NA"
  _total_hours_pmat2="NA"

  if [[ ! -s "$_timing_file" ]]; then
    echo "Error: timing file not found or empty: $_timing_file" >&2
    return 0
  fi

  local _output
  _output=$(_polap_lib_timing-parse-timing "$_timing_file" 2>/dev/null)

  # Accept any two non-whitespace strings
  if [[ "$_output" =~ ^[^[:space:]]+[[:space:]]+[^[:space:]]+ ]]; then
    read -r _memory_gb_pmat2 _total_hours_pmat2 <<<"$_output"
  else
    echo "Warning: malformed output from _polap_lib_timing-parse-timing" >&2
    echo "  Got: $_output" >&2
  fi

  return 0
}

extract_mafft_pident() {
  local _v1_inum="$1"
  local _pident_file="${_v1_inum}/mafft/1/pident.txt"

  local _mafft_pident="NA"

  if [[ -s "$_pident_file" ]]; then
    local _raw
    _raw="$(<"$_pident_file")"

    # Check if it is a valid number (integer or decimal)
    if [[ "$_raw" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
      _mafft_pident="$_raw"
    else
      echo "Warning: invalid number in $_pident_file: '$_raw'" >&2
    fi
  fi

  echo "${_mafft_pident}"

  return 0
}

# A copy of Table 1's row for a given outdir and inum.
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
# Done in pmat-0.1
##################
# Eucalyptus_pauciflora
# Dunaliella_tertiolecta
# Leiosporoceros_dussii
# Anthoceros_agrestis
# Spirodela_polyrhiza
#
maintable1_genus_species_for() {
  local _brg_outdir="${1}"
  local _brg_inum="${2:-10}"
  local _brg_adir="${3:-.}"
  local _brg_table="${4:-1}"
  local _brg_t_dir="${5:-"${_brg_default_target_dir}"}"

  local _key="${_brg_outdir}-${_brg_inum}"
  local downsample="${_downsample["${_key}"]}"
  local _polap_log=${_brg_outdir}/polap.log

  # Taxon names: species, genus, family, and order
  local _species="${_brg_outdir//_/ }"
  local _genus=${_species%% *}
  local _family=$(grep "${_genus}" ${_POLAPLIB_DIR}/taxonomy_output.tsv | cut -f 7 | head -1)
  local _order=$(grep "${_genus}" ${_POLAPLIB_DIR}/taxonomy_output.tsv | cut -f 6 | head -1)

  # Extract long-read SRA ID
  local _l_fq_stats="${_brg_outdir}/${_brg_adir}/${_brg_inum}/l.fq.stats"
  local _l_sra=$(extract_long_sra_seqkit_stats "${_l_fq_stats}")
  local _l_sra_size=$(extract_seqkit_stats_sum_len "${_l_fq_stats}")
  local _l_sra_size_gb=$(_polap_utility_convert_bp "${_l_sra_size}")

  # Extract short-read SRA ID
  local _s1_fq_stats="${_brg_outdir}/${_brg_adir}/${_brg_inum}/s1.fq.stats"
  local _s1_sra=$(extract_long_sra_seqkit_stats "${_s1_fq_stats}")
  local _s1_sra_size=$(extract_seqkit_stats_sum_len "${_s1_fq_stats}")
  local _s1_sra_size_gb=$(_polap_utility_convert_bp "${_s1_sra_size}")
  local _s2_fq_stats="${_brg_outdir}/${_brg_adir}/${_brg_inum}/s2.fq.stats"
  local _s2_sra=$(extract_long_sra_seqkit_stats "${_s2_fq_stats}")
  local _s2_sra_size=$(extract_seqkit_stats_sum_len "${_s2_fq_stats}")
  local _s2_sra_size_gb=$(_polap_utility_convert_bp "${_s2_sra_size}")
  local _s_sra_size=$((_s1_sra_size + _s2_sra_size))
  local _s_sra_size_gb=$(_polap_utility_convert_bp "${_s_sra_size}")

  # Extract the genome size estimate
  local _genome_size_dir="${_brg_outdir}/${_brg_adir}/${_brg_inum}"
  local _genome_size=$(extract_genome_size "${_genome_size_dir}")
  local _genome_size_gb=$(_polap_utility_convert_bp "${_genome_size}")

  # sequencing data coverage
  local _long_coverage=$(echo "scale=2; ${_l_sra_size} / ${_genome_size}" | bc)
  local _short_coverage=$(echo "scale=2; ${_s_sra_size} / ${_genome_size}" | bc)

  # Extract the timing for the corresponing coverage x case
  local _timing_dir="${_brg_outdir}/${_brg_adir}/${_brg_inum}"
  local _memory_gb='NA'
  local _total_hours='NA'
  extract_time_memory_timing "${_timing_dir}/timing-polap-assemble1.txt"
  local _memory_gb_polap_assemble1=$_memory_gb
  local _total_hours_polap_assemble1=$_total_hours

  local _summary_file="${_timing_dir}/summary-polap-analysis-reduce.txt"
  local _memory_gb='NA'
  local _total_hours='NA'
  extract_time_memory_timing_summary_file "${_summary_file}"
  local _memory_gb_polap_reduce=$_memory_gb
  local _total_hours_polap_reduce=$_total_hours

  local _summary_file="${_timing_dir}/summary-polap-analysis-assemble2.txt"
  local _memory_gb='NA'
  local _total_hours='NA'
  extract_time_memory_timing_summary_file "${_summary_file}"
  local _memory_gb_polap_assemble2=$_memory_gb
  local _total_hours_polap_assemble2=$_total_hours

  local _count_run='NA'
  extract_cumulative_time_memory_timing "${_timing_dir}/timing-polap-analysis-assemble2.txt"
  local _memory_gb_polap_assemble2_timing=$_memory_gb
  local _total_hours_polap_assemble2_timing=$_total_hours

  local _timing0_dir="${_brg_outdir}/${_brg_adir}/0"
  local _summary_file="${_timing0_dir}/summary-polap-analysis-msbwt.txt"
  local _memory_gb='NA'
  local _total_hours='NA'
  extract_time_memory_timing_summary_file "${_summary_file}"
  local _memory_gb_polap_msbwt=$_memory_gb
  local _total_hours_polap_msbwt=$_total_hours

  local _summary_file="${_timing_dir}/summary-polap-analysis-polish.txt"
  local _memory_gb='NA'
  local _total_hours='NA'
  extract_time_memory_timing_summary_file "${_summary_file}"
  local _memory_gb_polap_polish=$_memory_gb
  local _total_hours_polap_polish=$_total_hours

  # Total Polap processing time: assemble1
  local _total_hours_polap1=$(_polap_lib_timing-sum_time_values \
    "${_total_hours_polap_reduce}" \
    "${_total_hours_polap_assemble1}")

  # Peak memory for Polap: assemble1
  local _memory_gb_polap1=$(_polap_lib_number-max_float \
    "${_memory_gb_polap_reduce}" \
    "${_memory_gb_polap_assemble1}")

  # Total Polap processing time: assemble2
  local _total_hours_polap2=$(_polap_lib_timing-sum_time_values \
    "${_total_hours_polap_assemble2}" \
    "${_total_hours_polap_msbwt}" \
    "${_total_hours_polap_polish}")

  # Peak memory for Polap: assemble2
  local _memory_gb_polap2=$(_polap_lib_number-max_float \
    "${_memory_gb_polap_assemble2}" \
    "${_memory_gb_polap_msbwt}" \
    "${_memory_gb_polap_polish}")

  # Total Polap processing time
  local _total_hours_polap=$(_polap_lib_timing-sum_time_values \
    "${_total_hours_polap_reduce}" \
    "${_total_hours_polap_assemble1}" \
    "${_total_hours_polap_assemble2}" \
    "${_total_hours_polap_msbwt}" \
    "${_total_hours_polap_polish}")
  # Peak memory for Polap

  local _memory_gb_polap=$(_polap_lib_number-max_float \
    "${_memory_gb_polap_reduce}" \
    "${_memory_gb_polap_assemble1}" \
    "${_memory_gb_polap_assemble2}" \
    "${_memory_gb_polap_msbwt}" \
    "${_memory_gb_polap_polish}")

  # Extract nextDenovo polishing timing
  local _timing_dir="${_brg_outdir}/b1"
  extract_time_memory_timing "${_timing_dir}/timing-nextdenovo-polish.txt"
  local _underestimated_memory_gb_nextdenovo_polish=$_memory_gb
  local _total_hours_nextdenovo_polish=$_total_hours

  local _memory_gb_nextdenovo_polish=$(bash "${_polap_script_bin_dir}/polap-data-v2.sh" nextdenovo-polish-memory "${_brg_outdir}" b1)
  _memory_gb_nextdenovo_polish="${_memory_gb_nextdenovo_polish%GB}"

  # Extract PMAT timing
  local _fc_value=${pmat_fc[${_brg_outdir}]}
  local _timing_dir="${_brg_outdir}/b1"
  extract_time_memory_timing "${_timing_dir}/timing-pmat-nextdenovo-fc-${_fc_value}.txt"
  local _memory_gb_pmat_assembly=$_memory_gb
  local _total_hours_pmat_assembly=$_total_hours

  # Total PMAT processingg time
  local _total_hours_pmat=$(_polap_lib_timing-sum_time_values \
    "${_total_hours_nextdenovo_polish}" \
    "${_total_hours_pmat_assembly}")

  # Peak memory for PMAT
  local _memory_gb_pmat=$(_polap_lib_number-max_float \
    "${_memory_gb_nextdenovo_polish}" \
    "${_memory_gb_pmat_assembly}")

  # debug local variables
  # for debugging: Inline local printing local var
  # while IFS= read -r line; do
  # 	if [[ $line =~ ^declare\ --\ ([^=]+)= ]]; then
  # 		var="${BASH_REMATCH[1]}"
  # 		printf "%s=%q\n" "$var" "${!var}"
  # 	fi
  # done < <(local -p 2>/dev/null)
  # return

  # Known ptDNA NCBI accession
  # local _known_mtdna=$(extract_known_mtdna_accession "${_polap_log}")

  local _items
  if ((_brg_table == 1)); then
    _items=(
      "_${_species}_"
      "${_order}"
      "${_family}"
      "${_l_sra_size_gb}"
      "${_long_coverage}"
      "${_s_sra_size_gb}"
      "${_short_coverage}"
      "${_genome_size_gb}"
      "${downsample}"
      "${_memory_gb_polap_reduce}"
      "${_total_hours_polap_reduce}"
      "${_memory_gb_polap_assemble1}"
      "${_total_hours_polap_assemble1}"
      "${_count_run}"
      "${_memory_gb_polap_assemble2}"
      "${_total_hours_polap_assemble2}"
      "${_memory_gb_polap_msbwt}"
      "${_total_hours_polap_msbwt}"
      "${_memory_gb_polap_polish}"
      "${_total_hours_polap_polish}"
      "${_memory_gb_polap1}"
      "${_total_hours_polap1}"
      "${_memory_gb_polap2}"
      "${_total_hours_polap2}"
      "${_memory_gb_polap}"
      "${_total_hours_polap}"
      "${_memory_gb_nextdenovo_polish}"
      "${_total_hours_nextdenovo_polish}"
      "${_memory_gb_pmat_assembly}"
      "${_total_hours_pmat_assembly}"
      "${_memory_gb_pmat}"
      "${_total_hours_pmat}"
    )
  elif ((_brg_table == 2)); then
    _items=(
      "_${_species}_"
      "${_order}"
      "${_family}"
      "${_l_sra_size_gb}"
      "${_long_coverage}"
      "${_s_sra_size_gb}"
      "${_short_coverage}"
      "${_genome_size_gb}"
      "${downsample}"
      "${_memory_gb_polap_reduce}"
      "${_total_hours_polap_reduce}"
      "${_memory_gb_polap_assemble1}"
      "${_total_hours_polap_assemble1}"
      "${_count_run}"
      "${_memory_gb_polap_assemble2}"
      "${_total_hours_polap_assemble2}"
      "${_memory_gb_polap_msbwt}"
      "${_total_hours_polap_msbwt}"
      "${_memory_gb_polap_polish}"
      "${_total_hours_polap_polish}"
      "${_memory_gb_polap1}"
      "${_total_hours_polap1}"
      "${_memory_gb_polap2}"
      "${_total_hours_polap2}"
      "${_memory_gb_polap}"
      "${_total_hours_polap}"
      "${_memory_gb_nextdenovo_polish}"
      "${_total_hours_nextdenovo_polish}"
      "${_memory_gb_pmat_assembly}"
      "${_total_hours_pmat_assembly}"
      "${_memory_gb_pmat}"
      "${_total_hours_pmat}"
    )
  elif ((_brg_table == 3)); then
    _items=(
      "_${_species}_"
      "${_order}"
      "${_family}"
      "${_brg_inum}"
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
      "${_D}"
      "${_Alpha}"
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
      "${_total_hours_getorganelle}"
      "${_memory_gb_ptgaul}"
      "${_total_hours_ptgaul}"
      "${_memory_gb_pmat1}"
      "${_total_hours_pmat1}"
      "${_memory_gb_pmat2}"
      "${_total_hours_pmat2}"
      "${_memory_gb_prepare_polishing}"
      "${_memory_gb_ptgaul_polishing}"
      "${_memory_gb_subsampling_polishing}"
      "${_memory_gb}"
      "${_total_hours}"
    )
  fi

  printf "%s\t" "${_items[@]::${#_items[@]}-1}" >>"${_table_tsv}"
  printf "%s\n" "${_items[-1]}" >>"${_table_tsv}"
}

# A better version of maintable1
# comparing memory and computing time used by tools
# focusing on table3, which compares the computing time and memory
#
# nextDenovo correction time and memory needs more than command time -v.
# It is because nextDenovo creates jobs in the background, which cannot be
# tracked by command time -v.
#
# inum for polap: 2
# inum for pmat: 7
maintable1_genus_species() {
  local _brg_outdir="${1}"
  local _brg_adir="${2:-.}"
  local _brg_table="${3:-1}"
  local _brg_t_dir="${4:-"${_brg_default_target_dir}"}"

  # local _table_name=$(echo "${FUNCNAME[0]}" | cut -d'_' -f1)

  local _table_name="maintable"${_brg_table}
  local _table_tsv="${_table_name}.tsv"

  maintable1_genus_species_header "${_brg_table}" >"${_table_tsv}"

  if [[ "${_brg_outdir}" == "all" ]]; then
    for outdir in "${S2[@]}"; do
      for _brg_inum in 10 20; do
        if [[ "${outdir}" == "Taraxacum_mongolicum" && "${_brg_inum}" == "20" ]]; then
          continue
        fi
        if [[ "${outdir}" == "Salix_dunnii" && "${_brg_inum}" == "20" ]]; then
          continue
        fi
        maintable1_genus_species_for "${outdir}" "${_brg_inum}" "${@:2}" >>"${_table_tsv}"
        # echo maintable1_genus_species_for "${outdir}" "${_brg_inum}" "${@:2}"
      done
    done
  else
    # maintable1_genus_species_for "$@" >>"${_table_tsv}"
    maintable1_genus_species_for "${_brg_outdir}" 10 "${_brg_adir}" "${_brg_table}" "${_brg_t_dir}" >>"${_table_tsv}"
    # maintable1_genus_species_for "${_brg_outdir}" 10 "${_brg_adir}" "${_brg_table}" "${_brg_t_dir}"
  fi

  echo "Table: ${_table_tsv}"

  source "$(conda info --base)/etc/profile.d/conda.sh"
  if [[ "${CONDA_DEFAULT_ENV:-}" != "polap" ]]; then
    echo "[INFO] Activating conda environment 'polap'..."
    conda activate polap
  fi

  if [[ "${CONDA_DEFAULT_ENV:-}" != "polap" ]]; then
    echo "[ERROR] Failed to enter conda environment 'polap'"
    return
  fi

  if ((_brg_table == 1)); then
    csvtk -t cut -f Species,L_size,L_cov,S_size,S_cov,G,D,Mpa0,Tpa0,Mpa1,Tpa1,Cpa2,Mpa2,Tpa2,Mpp1,Tpp1,Mpp2,Tpp2,Mpolap1,Tpolap1,Mpolap2,Tpolap2,Mpolap,Tpolap,Mndp,Tndp,Mpmata,Tpmata,Mpmat,Tpmat \
      ${_table_tsv} |
      csvtk -t rename -f 1-30 -n Species,L_size,L_cov,S_size,S_cov,G,D,Mpa0,Tpa0,Mpa1,Tpa1,Cpa2,Mpa2,Tpa2,Mpp1,Tpp1,Mpp2,Tpp2,Mpolap1,Tpolap1,Mpolap2,Tpolap2,Mpolap,Tpolap,Mndp,Tndp,Mpmata,Tpmata,Mpmat,Tpmat |
      csvtk -t csv2md -a right -o ${_table_name}-analysis.md
    echo ${_table_name}-analysis.md
  elif ((_brg_table == 2)); then
    csvtk -t cut -f Species,Order,Family,L_cov,S_cov,D,Mpolap1,Tpolap1,Mpolap2,Tpolap2,Mpolap,Tpolap,Mndp,Tndp,Mpmata,Tpmata,Mpmat,Tpmat \
      ${_table_tsv} |
      csvtk -t rename -f 1-18 -n Species,Order,Family,L_cov,S_cov,D,Mpolap1,Tpolap1,Mpolap2,Tpolap2,Mpolap,Tpolap,Mndp,Tndp,Mpmata,Tpmata,Mpmat,Tpmat |
      csvtk -t csv2md -a right -o ${_table_name}-analysis.md
    echo ${_table_name}-analysis.md
  elif ((_brg_table == 3)); then
    csvtk -t cut -f Species,M,T_g,M_g,T_p,M_p,T_pmat1,T_pmat2,M_pmat2,M_t2,M_s,M_f,T \
      ${_table_tsv} |
      csvtk -t rename -f 1-13 -n Species,M,Tg,Mg,Tp,Mp,Tn,Ta,Ma,Mt,Ms,Mf,T |
      csvtk -t csv2md -a right -o ${_table_name}-${_brg_inum}-analysis.md
  fi

  conda deactivate

  # if ((_brg_format == 3)); then
  # 	csvtk -t cut -f Species,A0,D,P,R,Rate,Alpha,Length1,Length2,Pident,N1,Mode,SD,M,M_g,M_p,M_t2,M_s,M_f,T \
  # 		${_table_tsv} |
  # 		csvtk -t rename -f 1-20 -n Species,A0,D,P,R,Rate,Alpha,L1,L2,Pident,N1,Mode,SD,M,Mg,Mp,Mt,Ms,Mf,T |
  # 		csvtk -t csv2md -a right -o ${_table_name}-${_brg_inum}-analysis.md
  # fi

  # csvtk -t cut -f Species,Order,Family,L_SRA,L_size,L_cov,S_SRA,S_size,S_cov ${_table_tsv} |
  # 	csvtk -t rename -f 1-9 -n Species,Order,Family,L_SRA,L_size,L_cov,S_SRA,S_size,S_cov |
  # 	csvtk -t csv2md -a right -o ${_table_name}-${_brg_inum}-data.md

  # echo "  copying to ${_brg_t_dir}"
  # cp -p ${_table_name}-${_brg_inum}-analysis.md "${_brg_t_dir}"
  # cp -p ${_table_name}-${_brg_inum}-data.md "${_brg_t_dir}"

  # csvtk -t csv2md -a right ${_table_tsv} -o ${_table_name}-${_brg_inum}.md

  # Rscript "${_POLAPLIB_DIR}"/polap-data-v2.R --inum ${_brg_inum} \
  # 	>${_table_name}-${_brg_inum}.txt

  # cat "${_table_name}-${_brg_inum}.md"
  # echo "ouput: ${_table_name}-${_brg_inum}-analysis.md"
  # echo "ouput: ${_table_name}-${_brg_inum}-data.md"
  # echo "ouput: ${_table_name}-${_brg_inum}.txt"
  # echo "ouput: ${_table_name}-${_brg_inum}.md"
  # echo "ouput: ${_table_tsv}"

  # echo cp -p ${_table_name}-${_brg_inum}-analysis.md "${_brg_t_dir}"
}

supptable1_genus_species_header() {
  local _brg_table="${1:-1}"

  local _items

  if ((_brg_table == 1)); then
    local _items=(
      "Species"
      "Order"
      "Family"
      "Mndp"
      "Tndp"
      "Mpmata"
      "Tpmata"
    )
  elif ((_brg_table == 2)); then
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
  elif ((_brg_table == 3)); then
    _items=(
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
      "D"
      "A0"
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
      "T_g"
      "M_p"
      "T_p"
      "M_pmat1"
      "T_pmat1"
      "M_pmat2"
      "T_pmat2"
      "M_t1"
      "M_t2"
      "M_s"
      "M_f"
      "T"
    )
  fi

  printf "%s\t" "${_items[@]::${#_items[@]}-1}"
  printf "%s\n" "${_items[-1]}"
}

supptable1_genus_species_for() {
  local _brg_outdir="${1}"
  # local _brg_inum="${2:-10}"
  local _brg_table="${2:-1}"
  local _brg_t_dir="${3:-"${_brg_default_target_dir}"}"
  local _brg_inum="10"

  local _table_name="supptable"${_brg_table}
  local _table_tsv="${_table_name}.tsv"

  local _key="${_brg_outdir}-${_brg_inum}"
  local downsample="${_downsample["${_key}"]}"
  local _polap_log=${_brg_outdir}/polap.log

  # Taxon names: species, genus, family, and order
  local _species="${_brg_outdir//_/ }"
  local _genus=${_species%% *}
  local _family=$(grep "${_genus}" ${_POLAPLIB_DIR}/taxonomy_output.tsv | cut -f 7 | head -1)
  local _order=$(grep "${_genus}" ${_POLAPLIB_DIR}/taxonomy_output.tsv | cut -f 6 | head -1)

  # Extract nextDenovo polishing timing
  local _timing_dir="${_brg_outdir}/b1"
  extract_time_memory_timing "${_timing_dir}/timing-nextdenovo-polish.txt"
  local _underestimated_memory_gb_nextdenovo_polish=$_memory_gb
  local _total_hours_nextdenovo_polish=$_total_hours

  local _memory_gb_nextdenovo_polish=$(bash "${_polap_script_bin_dir}/polap-data-v2.sh" nextdenovo-polish-memory "${_brg_outdir}" b1)
  _memory_gb_nextdenovo_polish="${_memory_gb_nextdenovo_polish%GB}"

  # Extract PMAT timing
  local _fc_value=${pmat_fc[${_brg_outdir}]}
  local _timing_dir="${_brg_outdir}/b1"
  extract_time_memory_timing "${_timing_dir}/timing-pmat-nextdenovo-fc-${_fc_value}.txt"
  local _memory_gb_pmat_assembly=$_memory_gb
  local _total_hours_pmat_assembly=$_total_hours

  # Total PMAT processingg time
  local _total_hours_pmat=$(_polap_lib_timing-sum_time_values \
    "${_total_hours_nextdenovo_polish}" \
    "${_total_hours_pmat_assembly}")

  # Peak memory for PMAT
  local _memory_gb_pmat=$(_polap_lib_number-max_float \
    "${_memory_gb_nextdenovo_polish}" \
    "${_memory_gb_pmat_assembly}")

  # debug local variables
  # for debugging: Inline local printing local var
  # while IFS= read -r line; do
  # 	if [[ $line =~ ^declare\ --\ ([^=]+)= ]]; then
  # 		var="${BASH_REMATCH[1]}"
  # 		printf "%s=%q\n" "$var" "${!var}"
  # 	fi
  # done < <(local -p 2>/dev/null)
  # return

  # Known ptDNA NCBI accession
  # local _known_mtdna=$(extract_known_mtdna_accession "${_polap_log}")

  local _items
  if ((_brg_table == 1)); then
    _items=(
      "_${_species}_"
      "${_order}"
      "${_family}"
      "${_memory_gb_nextdenovo_polish}"
      "${_total_hours_nextdenovo_polish}"
      "${_memory_gb_pmat_assembly}"
      "${_total_hours_pmat_assembly}"
    )
  elif ((_brg_table == 2)); then
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
  elif ((_brg_table == 3)); then
    _items=(
      "_${_species}_"
      "${_order}"
      "${_family}"
      "${_brg_inum}"
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
      "${_D}"
      "${_Alpha}"
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
      "${_total_hours_getorganelle}"
      "${_memory_gb_ptgaul}"
      "${_total_hours_ptgaul}"
      "${_memory_gb_pmat1}"
      "${_total_hours_pmat1}"
      "${_memory_gb_pmat2}"
      "${_total_hours_pmat2}"
      "${_memory_gb_prepare_polishing}"
      "${_memory_gb_ptgaul_polishing}"
      "${_memory_gb_subsampling_polishing}"
      "${_memory_gb}"
      "${_total_hours}"
    )
  fi

  printf "%s\t" "${_items[@]::${#_items[@]}-1}" >>"${_table_tsv}"
  printf "%s\n" "${_items[-1]}" >>"${_table_tsv}"
}

# A better version of supptable1
# comparing memory and computing time used by tools
# focusing on table3, which compares the computing time and memory
#
# nextDenovo correction time and memory needs more than command time -v.
# It is because nextDenovo creates jobs in the background, which cannot be
# tracked by command time -v.
#
# inum for polap: 2
# inum for pmat: 7
supptable1_genus_species() {
  local _brg_outdir="${1:-all}"
  local _brg_table="${2:-1}"
  local _brg_t_dir="${3:-"${_brg_default_target_dir}"}"
  # local _table_name=$(echo "${FUNCNAME[0]}" | cut -d'_' -f1)

  local _table_name="supptable"${_brg_table}
  local _table_tsv="${_table_name}.tsv"

  supptable1_genus_species_header "${_brg_table}" >"${_table_tsv}"

  if [[ "${_brg_outdir}" == "all" ]]; then
    for outdir in "${S[@]}"; do
      supptable1_genus_species_for "${outdir}" "${@:2}" >>"${_table_tsv}"
    done
  else
    # supptable1_genus_species_for "${_brg_outdir}" "${_brg_table}" "${_brg_t_dir}" >>"${_table_tsv}"
    supptable1_genus_species_for "${_brg_outdir}" "${_brg_table}" "${_brg_t_dir}" >>"${_table_tsv}"
  fi

  echo "Table: ${_table_tsv}"

  source "$(conda info --base)/etc/profile.d/conda.sh"
  if [[ "${CONDA_DEFAULT_ENV:-}" != "polap" ]]; then
    echo "[INFO] Activating conda environment 'polap'..."
    conda activate polap
  fi

  if [[ "${CONDA_DEFAULT_ENV:-}" != "polap" ]]; then
    echo "[ERROR] Failed to enter conda environment 'polap'"
    return
  fi

  if ((_brg_table == 1)); then
    csvtk -t cut -f Species,Order,Family,Mndp,Tndp,Mpmata,Tpmata \
      ${_table_tsv} |
      csvtk -t rename -f 1-7 -n Species,Order,Family,Mndp,Tndp,Mpmata,Tpmata |
      csvtk -t csv2md -a right -o ${_table_name}-analysis.md
    echo ${_table_name}-analysis.md
  elif ((_brg_table == 2)); then
    csvtk -t cut -f Species,C,N,P,R,Rate,Alpha,Length2,N1,Mode,SD,M,M_g,M_f,T \
      ${_table_tsv} |
      csvtk -t rename -f 1-15 -n Species,C,N,P,R,Rate,Alpha,L2,N1,Mode,SD,M,Mg,Mf,T |
      csvtk -t csv2md -a right -o ${_table_name}-${_brg_inum}-analysis.md
  elif ((_brg_table == 3)); then
    csvtk -t cut -f Species,M,T_g,M_g,T_p,M_p,T_pmat1,T_pmat2,M_pmat2,M_t2,M_s,M_f,T \
      ${_table_tsv} |
      csvtk -t rename -f 1-13 -n Species,M,Tg,Mg,Tp,Mp,Tn,Ta,Ma,Mt,Ms,Mf,T |
      csvtk -t csv2md -a right -o ${_table_name}-${_brg_inum}-analysis.md
  fi

  conda deactivate

  # if ((_brg_format == 3)); then
  # 	csvtk -t cut -f Species,A0,D,P,R,Rate,Alpha,Length1,Length2,Pident,N1,Mode,SD,M,M_g,M_p,M_t2,M_s,M_f,T \
  # 		${_table_tsv} |
  # 		csvtk -t rename -f 1-20 -n Species,A0,D,P,R,Rate,Alpha,L1,L2,Pident,N1,Mode,SD,M,Mg,Mp,Mt,Ms,Mf,T |
  # 		csvtk -t csv2md -a right -o ${_table_name}-${_brg_inum}-analysis.md
  # fi

  # csvtk -t cut -f Species,Order,Family,L_SRA,L_size,L_cov,S_SRA,S_size,S_cov ${_table_tsv} |
  # 	csvtk -t rename -f 1-9 -n Species,Order,Family,L_SRA,L_size,L_cov,S_SRA,S_size,S_cov |
  # 	csvtk -t csv2md -a right -o ${_table_name}-${_brg_inum}-data.md

  # echo "  copying to ${_brg_t_dir}"
  # cp -p ${_table_name}-${_brg_inum}-analysis.md "${_brg_t_dir}"
  # cp -p ${_table_name}-${_brg_inum}-data.md "${_brg_t_dir}"

  # csvtk -t csv2md -a right ${_table_tsv} -o ${_table_name}-${_brg_inum}.md

  # Rscript "${_POLAPLIB_DIR}"/polap-data-v2.R --inum ${_brg_inum} \
  # 	>${_table_name}-${_brg_inum}.txt

  # cat "${_table_name}-${_brg_inum}.md"
  # echo "ouput: ${_table_name}-${_brg_inum}-analysis.md"
  # echo "ouput: ${_table_name}-${_brg_inum}-data.md"
  # echo "ouput: ${_table_name}-${_brg_inum}.txt"
  # echo "ouput: ${_table_name}-${_brg_inum}.md"
  # echo "ouput: ${_table_tsv}"

  # echo cp -p ${_table_name}-${_brg_inum}-analysis.md "${_brg_t_dir}"
}

# Anthoceros_agrestis
# Brassica_rapa
run_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"

  cd "${output_dir}" || exit
  common_operations "${long_sra}" "${short_sra}"
  if [[ -s "o/0/mt.contig.name-1" ]]; then
    ${_polap_cmd} assemble2
  else
    ${_polap_cmd} assemble
  fi
  cd -
}

subsample_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"
  local random_seed="${_random_seed["$1"]}"
  local genomesize="${_genomesize["$1"]}"

  # subsampling
  local _c
  for _c in 5 10 30; do
    rm -f "${output_dir}/l${_c}x.fq.gz"
    rm -f "${output_dir}/s${_c}x_1.fq.gz"
    rm -f "${output_dir}/s${_c}x_2.fq.gz"

    if [[ "${species_name}" == "Anthoceros agrestis" ||
      "${species_name}" == "Taraxacum mongolicum" ]]; then

      ${_polap_cmd} fastq subsample \
        "${output_dir}/${long_sra}.fastq" \
        "${output_dir}/l${_c}x.fq.gz" \
        -c "${_c}" \
        -o "${output_dir}" \
        --random-seed "${random_seed}" \
        --genomesize "${genomesize}" -v \
        >"${output_dir}/l${_c}x.txt"

      ${_polap_cmd} fastq subsample2 \
        "${output_dir}/${short_sra}_1.fastq" \
        "${output_dir}/${short_sra}_2.fastq" \
        "${output_dir}/s${_c}x_1.fq" \
        "${output_dir}/s${_c}x_2.fq" \
        -c "${_c}" \
        -o "${output_dir}" \
        --random-seed "${random_seed}" \
        --genomesize "${genomesize}" -v \
        >"${output_dir}/s${_c}x.txt"
    else
      ${_polap_cmd} fastq subsample \
        "${output_dir}/${long_sra}.fastq" \
        "${output_dir}/l${_c}x.fq.gz" \
        -c "${_c}" \
        -o "${output_dir}" \
        --random-seed "${random_seed}" \
        --species "${species_name}" -v \
        >"${output_dir}/l${_c}x.txt"

      ${_polap_cmd} fastq subsample2 \
        "${output_dir}/${short_sra}_1.fastq" \
        "${output_dir}/${short_sra}_2.fastq" \
        "${output_dir}/s${_c}x_1.fq" \
        "${output_dir}/s${_c}x_2.fq" \
        -c "${_c}" \
        -o "${output_dir}" \
        --random-seed "${random_seed}" \
        --species "${species_name}" -v \
        >"${output_dir}/s${_c}x.txt"
    fi
  done

  # Brassica_rapa
  # genome size: 364511000
  # long_total_length.txt
}

assemble1-subsample_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"

  # subsampling
  local _c
  for _c in 5 10 30; do
    assemble1-subsample_genus_species_for \
      "${output_dir}" \
      "${_c}"
  done
}

bandage1-subsample_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"

  # subsampling
  local _c
  for _c in 5 10 30; do
    bandage1-subsample_genus_species_for \
      "${output_dir}" \
      "${_c}"
  done
}

seeds-subsample_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"

  # subsampling
  local _c
  for _c in 5 10 30; do
    seeds-subsample_genus_species_for \
      "${output_dir}" \
      "${_c}"
  done
}

assemble2-subsample_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"

  # subsampling
  local _c
  local _i
  for _c in 5 10 30; do
    for _i in {1..3}; do
      local _mtcontigname="${output_dir}/subsample/${_c}x/o1/0/mt.contig.name-${_i}"
      if [[ -s "${_mtcontigname}" ]]; then
        assemble2-subsample_genus_species_for \
          "${output_dir}" \
          "${_c}" "${_i}"
      fi
    done
  done
}

bandage2-subsample_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"

  # subsampling
  local _c
  local _i
  for _c in 5 10 30; do
    for _i in {1..3}; do
      local _mtcontigname="${output_dir}/subsample/${_c}x/o1/0/mt.contig.name-${_i}"
      if [[ -s "${_mtcontigname}" ]]; then
        bandage2-subsample_genus_species_for \
          "${output_dir}" \
          "${_c}" "${_i}"
      fi
    done
  done
}

assemble1-subsample_genus_species_for() {
  local output_dir="$1"
  local coverage="$2"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"
  local subsample_dir="subsample/${coverage}x"
  local subsample1_dir="${subsample_dir}/o1"
  local subsample2_dir="${subsample_dir}/o2"
  local long_fastq="l${coverage}x.fq"
  local short1_fastq="s${coverage}x_1.fq"
  local short2_fastq="s${coverage}x_2.fq"

  cd "${output_dir}" || exit
  mkdir -p "${subsample1_dir}"
  # To gunzip a .gz file to a specific folder
  # while keeping the original .gz file, use:
  if [[ -s "${long_fastq}.gz" ]]; then

    gzip -dc "${long_fastq}.gz" >"${subsample_dir}/${long_fastq}"
    gzip -dc "${short1_fastq}.gz" >"${subsample_dir}/${short1_fastq}"
    gzip -dc "${short2_fastq}.gz" >"${subsample_dir}/${short2_fastq}"
    command time -v polap assemble1 \
      -o "${subsample1_dir}" \
      -l "${subsample_dir}/${long_fastq}" \
      -a "${subsample_dir}/${short1_fastq}" \
      -b "${subsample_dir}/${short2_fastq}" \
      2>timing-polap-assemble1-subsample-${coverage}x.txt

    echo "Execute for the annotation table: polap annotate view -o ${output_dir}/${subsample1_dir}"
    echo "Execute to try to automatically create seeds: polap seeds -o ${output_dir}/${subsample1_dir}"
    echo "Make sure that you have seed contigs in file: ${output_dir}/${subsample1_dir}/0/mt.contig.name-1"
    echo "Use Bandage to view ${output_dir}/${subsample1_dir}/0/30-contigger/graph_final.gfa"
    echo "  to locate edge_<number> at the top of the annotation table."
  else
    echo "WARNING: no long-read data: ${long_fastq}.gz"
    echo "  deleting ${subsample_dir}"
    rm -rf "${subsample_dir}"
  fi

  cd -
}

seeds-subsample_genus_species_for() {
  local output_dir="$1"
  local coverage="$2"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"
  local subsample_dir="subsample/${coverage}x"
  local subsample1_dir="${subsample_dir}/o1"
  local subsample2_dir="${subsample_dir}/o2"
  local long_fastq="l${coverage}x.fq"
  local short1_fastq="s${coverage}x_1.fq"
  local short2_fastq="s${coverage}x_2.fq"

  cd "${output_dir}" || exit

  if [[ -d "${subsample_dir}" ]]; then
    polap seeds \
      -o "${subsample1_dir}"
  else
    echo "WARNING: no such folder: ${subsample_dir}"
  fi

  echo "Execute for the annotation table with seeds: polap seeds view -o ${output_dir}/${subsample1_dir} -j <target_assembly>"
  echo "Execute for organelle-genome assembly: polap assemble2 -o ${output_dir}/${subsample1_dir} -j <target_assembly>"

  cd -
}

assemble2-subsample_genus_species_for() {
  local output_dir="$1"
  local coverage="$2"
  local target_assembly="${3:-1}"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"
  local subsample_dir="subsample/${coverage}x"
  local subsample1_dir="${subsample_dir}/o1"
  local subsample2_dir="${subsample_dir}/o2"
  local long_fastq="l${coverage}x.fq"
  local short1_fastq="s${coverage}x_1.fq"
  local short2_fastq="s${coverage}x_2.fq"

  cd "${output_dir}" || exit

  if [[ -d "${subsample_dir}" ]]; then
    command time -v polap assemble2 \
      -j "${target_assembly}" \
      -o "${subsample1_dir}" \
      2>timing-polap-assemble2-subsample-${coverage}x.txt

    echo "Use Bandage to view ${subsample1_dir}/1/assembly_graph.gfa"
    echo "  to locate and extract mtDNA sequences."
  else
    echo "WARNING: no such folder: ${subsample_dir}"
  fi

  cd -
}

bandage1-subsample_genus_species_for() {
  local output_dir="$1"
  local coverage="$2"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"
  local subsample_dir="subsample/${coverage}x"
  local subsample1_dir="${subsample_dir}/o1"
  local subsample2_dir="${subsample_dir}/o2"
  local long_fastq="l${coverage}x.fq"
  local short1_fastq="s${coverage}x_1.fq"
  local short2_fastq="s${coverage}x_2.fq"

  local _gfa="${output_dir}/${subsample1_dir}/0/30-contigger/graph_final.gfa"
  if [[ -s "${_gfa}" ]]; then
    ${_polap_cmd} bandage png \
      ${_gfa} \
      ${output_dir}/subsample/whole-genome-assembly-${coverage}x.png
    echo "file:	${output_dir}/subsample/whole-genome-assembly-${coverage}x.png"
  else
    echo "Warning: no such file: ${_gfa}"
    echo "FIX QT5 problem:"
    echo "export QT_QPA_PLATFORM=offscreen"
  fi
}

bandage2-subsample_genus_species_for() {
  local output_dir="$1"
  local coverage="$2"
  local target_assembly="${3:-1}"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"
  local subsample_dir="subsample/${coverage}x"
  local subsample1_dir="${subsample_dir}/o1"
  local subsample2_dir="${subsample_dir}/o2"
  local long_fastq="l${coverage}x.fq"
  local short1_fastq="s${coverage}x_1.fq"
  local short2_fastq="s${coverage}x_2.fq"

  local _gfa="${output_dir}/${subsample1_dir}/${target_assembly}/assembly_graph.gfa"
  if [[ -s "${_gfa}" ]]; then
    ${_polap_cmd} bandage png \
      ${_gfa} \
      ${output_dir}/subsample/organelle-genome-assembly-${coverage}x.png
    echo "file: ${output_dir}/subsample/organelle-genome-assembly-${coverage}x.png"
  else
    echo "Warning: no such file: ${_gfa}"
  fi
}

polish-subsample_genus_species() {
  local output_dir="$1"
  local coverage="$2"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"
  local subsample_dir="subsample/${coverage}x"
  local subsample1_dir="${subsample_dir}/o1"
  local subsample2_dir="${subsample_dir}/o2"
  local long_fastq="l${coverage}x.fq"
  local short1_fastq="s${coverage}x_1.fq"
  local short2_fastq="s${coverage}x_2.fq"

  cd "${output_dir}" || exit
  mkdir -p "${subsample1_dir}"
  # To gunzip a .gz file to a specific folder
  # while keeping the original .gz file, use:

  polap polish \
    -o "${subsample1_dir}"

  echo "Execute: polap annotate view -o ${subsample1_dir}"

  cd -
}

write-config_genus_species() {
  local csv_file="$1"

  # Write the output to CSV
  {
    # Print the header
    echo "species,_long,_short,_host,_inref"

    # Loop through species and print their values
    for key in "${!_long[@]}"; do
      echo "$key,${_long[$key]},${_short[$key]},${_host[$key]},${_inref[$key]}"
    done
  } >"$csv_file"
}

# Common operations function
common_operations_subsample2() {
  local long_sra="$1"
  local short_sra="$2"

  source ~/miniconda3/bin/activate polap
  if [[ -s "o1/0/mt.contig.name-1" ]]; then
    mkdir -p o/0
    cp -p o2/*.txt o
    cp -p o1/0/mt.contig.name-1 o/0/
    cp -pr o1/0/30-contigger o/0/
  fi
  if [[ ! -s "${long_sra}.fastq" ]]; then
    ${_polap_cmd} x-ncbi-fetch-sra --sra "$long_sra"
    ${_polap_cmd} x-ncbi-fetch-sra --sra "$short_sra"
    ${_polap_cmd} x-ncbi-fetch-sra-runinfo --sra "$long_sra"
    ${_polap_cmd} x-ncbi-fetch-sra-runinfo --sra "$short_sra"
  fi
  cp -s "${long_sra}.fastq" l.fq
  cp -s "${short_sra}_1.fastq" s1.fq
  cp -s "${short_sra}_2.fastq" s2.fq
  if [[ -s "o/0/mt.contig.name-1" ]]; then
    ${_polap_cmd} reduce-data
  fi
}

run-subsample3_genus_species() {
  local output_dir="$1"
  local coverage="$2"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"

  cd "${output_dir}" || exit
  local long_fastq=l${coverage}x.fq
  gunzip "${long_fastq}.gz"
  common_operations_subsample "${long_fastq}" "${short_sra}"
  if [[ -s "o/0/mt.contig.name-1" ]]; then
    ${_polap_cmd} assemble2
  else
    ${_polap_cmd} assemble
  fi
  cd -
}

# Function for each species
run_anthoceros_agrestis() {
  cd Anthoceros_agrestis || exit
  common_operations "SRR10190639" "SRR10250248"
  if [[ -s "o/0/mt.contig.name-1" ]]; then
    ${_polap_cmd} assemble2
  else
    ${_polap_cmd} assemble
  fi
  cd -
}

run_brassica_rapa() {
  cd Brassica_rapa || exit
  common_operations "ERR6210792" "ERR6210790"
  if [[ -s "o/0/mt.contig.name-1" ]]; then
    ${_polap_cmd} assemble2
  else
    ${_polap_cmd} assemble
  fi
  cd -
}

run_vigna_radiata() {
  cd Vigna_radiata || exit
  common_operations "SRR12549541" "SRR12549533"
  if [[ -s "o/0/mt.contig.name-1" ]]; then
    ${_polap_cmd} assemble2 --polap-reads
  else
    ${_polap_cmd} assemble --polap-reads
  fi
  cd -
}

run_trifolium_pratense() {
  cd Trifolium_pratense || exit
  common_operations "SRR15433794" "SRR15433795"
  if [[ -s "o/0/mt.contig.name-1" ]]; then
    ${_polap_cmd} assemble2
  else
    ${_polap_cmd} assemble
  fi
  cd -
}

run_taraxacum_mongolicum() {
  cd Taraxacum_mongolicum || exit
  common_operations "SRR19182970" "SRR19182971"
  if [[ -s "o/0/mt.contig.name-1" ]]; then
    ${_polap_cmd} assemble2
  else
    ${_polap_cmd} assemble
  fi
  cd -
}

run_spirodela_polyrhiza() {
  cd Spirodela_polyrhiza || exit
  common_operations "SRR11472010" "SRR11472009"
  if [[ -s "o/0/mt.contig.name-1" ]]; then
    ${_polap_cmd} assemble2
  else
    ${_polap_cmd} assemble
  fi
  cd -
}

run_salix_dunnii() {
  cd Salix_dunnii || exit
  common_operations "SRR12893432" "SRR12893433"
  if [[ -s "o/0/mt.contig.name-1" ]]; then
    ${_polap_cmd} assemble2 --polap-reads
  else
    ${_polap_cmd} assemble --polap-reads
  fi
  cd -
}

run_punica_granatum() {
  cd Punica_granatum || exit
  common_operations "SRR24893686" "SRR24893685"
  if [[ -s "o/0/mt.contig.name-1" ]]; then
    ${_polap_cmd} assemble2 --polap-reads
  else
    ${_polap_cmd} assemble --polap-reads
  fi
  cd -
}

run_macadamia_tetraphylla() {
  cd Macadamia_tetraphylla || exit
  common_operations "SRR10424548" "SRR10424549"
  if [[ -s "o/0/mt.contig.name-1" ]]; then
    ${_polap_cmd} assemble2 --polap-reads
  else
    ${_polap_cmd} assemble --polap-reads
  fi
  cd -
}

run_lolium_perenne() {
  cd Lolium_perenne || exit
  common_operations "SRR13386519" "SRR13386518"
  if [[ -s "o/0/mt.contig.name-1" ]]; then
    ${_polap_cmd} assemble2 --polap-reads -w 11000 --max-seeds 50
  else
    ${_polap_cmd} assemble --polap-reads -w 11000 --max-seeds 50
  fi
  cd -
}

run_anthoceros_angustus() {
  cd Anthoceros_angustus || exit
  common_operations "SRR9696346" "SRR9662965"
  if [[ -s "o/0/mt.contig.name-1" ]]; then
    ${_polap_cmd} assemble2
  else
    ${_polap_cmd} assemble
  fi
  cd -
}

run_carex_pseudochinensis() {
  mkdir Carex_pseudochinensis
  cd Carex_pseudochinensis || exit
  common_operations "SRR30757341" "SRR30757340"
  if [[ -s "o/0/mt.contig.name-1" ]]; then
    ${_polap_cmd} assemble2 --polap-reads
  else
    ${_polap_cmd} assemble --polap-reads
  fi
  cd -
}

################################################################################
# Part of genus_species
#
test_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"
  echo host: $(hostname)
  echo output: $output_dir
  echo species: $species_name
  if [[ -s "${long_sra}.fastq" ]]; then
    echo long: $long_sra
  else
    echo long: no such fastq file
    echo run $0 send-data-to at ${_local_host}
    echo ncbitools fetch sra ${long_sra}
  fi
  if [[ -s "${short_sra}_1.fastq" ]] &&
    [[ -s "${short_sra}_2.fastq" ]]; then
    echo short: $short_sra
  else
    echo short: no such fastq file
    echo run $0 send-data-to at ${_local_host}
    echo ncbitools fetch sra ${short_sra}
  fi
}

send-data-to_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"

  local _media_dir="${_media_dir}/${output_dir}"

  if [[ "${_local_host}" == "$(hostname)" ]]; then
    local f="${_media_dir}/${long_sra}.fastq"
    if [[ -s "${f}" ]]; then
      scp "${f}" ${_host["${output_dir}"]}:"$PWD"/
    else
      echo "ERROR: no such file: ${f}"
    fi
    local f="${_media_dir}/${short_sra}_1.fastq"
    if [[ -s "${f}" ]]; then
      scp "${f}" ${_host["${output_dir}"]}:"$PWD"/
    else
      echo "ERROR: no such file: ${f}"
    fi
    local f="${_media_dir}/${short_sra}_2.fastq"
    if [[ -s "${f}" ]]; then
      scp "${f}" ${_host["${output_dir}"]}:"$PWD"/
    else
      echo "ERROR: no such file: ${f}"
    fi
  else
    echo "ERROR: run at the local host."
  fi
}

send-data-to-single_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"

  if [[ "${_local_host}" == "$(hostname)" ]]; then
    local f="${_media_dir}/${long_sra}.fastq"
    if [[ -s "${f}" ]]; then
      scp "${f}" ${_host["${output_dir}"]}:"$PWD"/
    else
      echo "ERROR: no such file: ${f}"
    fi
    local f="${_media_dir}/${short_sra}_1.fastq"
    if [[ -s "${f}" ]]; then
      scp "${f}" ${_host["${output_dir}"]}:"$PWD"/
    else
      echo "ERROR: no such file: ${f}"
    fi
    local f="${_media_dir}/${short_sra}_2.fastq"
    if [[ -s "${f}" ]]; then
      scp "${f}" ${_host["${output_dir}"]}:"$PWD"/
    else
      echo "ERROR: no such file: ${f}"
    fi
  else
    echo "ERROR: run at the local host."
  fi
}

mkdir_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"

  echo "create $output_dir ..."
  mkdir -p $output_dir
}

taxon-reference_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local inref="${_inref["$1"]}"
  # copy_data

  # --taxonomy-rank-ingroup "${inref}" \
  ${_polap_cmd} taxonomy reference \
    --redo \
    --steps-include 1-9 \
    -o "${output_dir}" \
    --taxonomy-rank-ingroup "${inref}" \
    --species "${species_name}"
}

taxon-known_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local inref="${_inref["$1"]}"
  # copy_data

  # --taxonomy-rank-ingroup "${inref}" \
  ${_polap_cmd} taxonomy reference \
    --redo \
    --yes \
    --steps-include 1-9 \
    -o "${output_dir}"/known-mtdna \
    --taxonomy-rank-ingroup "${inref}" \
    --species "${species_name}"
  cp -p "${output_dir}/known-mtdna/ptgaul-reference.fa" \
    "${output_dir}/ptgaul-known.fa"
}

get-ptdna-from-ncbi_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"

  if [[ "${output_dir}" == "Juncus_inflexus" ]]; then
    species_name="Juncus effusus"
    echo "No ptDNA for ${output_dir}, so we use ${species_name}"
  fi
  ${_polap_cmd} get-mtdna \
    --plastid \
    --species "${species_name}" \
    -o ${output_dir}
}

copy-ptdna-of-ncbi-as-reference_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"

  echo "copy ${output_dir}/ptdna-reference.fa"
  cp -p "${output_dir}/00-bioproject/2-mtdna.fasta" \
    "${output_dir}/ptdna-reference.fa"
}

ptgaul_genus_species() {
  local output_dir="$1"
  local _type="${2}"

  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"

  local average_length=$(seqkit stats -Ta ${output_dir}/ptgaul-"${_type}".fa | awk 'NR==2 {print $7}')
  local genome_size=${average_length%.*}

  local _ptgaul_dir="${output_dir}/ptgaul-${_type}"
  mkdir -p ${_ptgaul_dir}

  local m
  local f=3000
  for m in 1000 3000 6000 9000 12000; do
    command time -v bash src/polap-ptGAUL1.sh \
      -o ${output_dir}-ptgaul \
      -r ${output_dir}/ptgaul-${_type}.fa \
      -g "${genome_size}" \
      -l "${output_dir}/${long_sra}.fastq" \
      -m "${m}" \
      -t 8 \
      2>${output_dir}/timing-ptgaul-${_type}-${m}.txt

    rm -rf "${_ptgaul_dir}/m${m}"
    mv "${output_dir}-ptgaul/result_${f}" "${_ptgaul_dir}/m${m}"
    rm -rf "${output_dir}-ptgaul"
  done
}

ptgaul-reference_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  ptgaul_genus_species "${output_dir}" "reference"
}

ptgaul-known_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  ptgaul_genus_species "${output_dir}" "known"
}

bandage-ptgaul_genus_species() {
  local output_dir="$1"
  local _type="${2}"

  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"

  local _ptgaul_dir="${output_dir}/ptgaul-${_type}"

  local m
  local f=3000
  for m in 1000 3000 6000 9000 12000; do
    local _gfa="${_ptgaul_dir}/m${m}/flye_cpONT/assembly_graph.gfa"
    if [[ -s "${_gfa}" ]]; then
      ${_polap_cmd} bandage png \
        ${_gfa} \
        ${_ptgaul_dir}/assembly_graph-${m}.png
    else
      echo "Warning: no such file: ${_gfa}"
    fi
  done
  echo "check ${_ptgaul_dir} for png files."
}

bandage-ptgaul-known_genus_species() {
  local output_dir="$1"
  bandage-ptgaul_genus_species "${output_dir}" "known"
}

bandage-ptgaul-reference_genus_species() {
  local output_dir="$1"
  bandage-ptgaul_genus_species "${output_dir}" "reference"
}

msbwt_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"

  command time -v ${_polap_cmd} prepare-polishing \
    -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
    -o ${output_dir} \
    2>${output_dir}/timing-prepare-polishing.txt
}

extract-ptdna-of-ptgaul_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"

  # extract ptGAUL result
  echo "extract ptDNA from the ptGAUL result with fmlrc polishing"
  command time -v ${_polap_cmd} disassemble ptgaul \
    -v -v -v \
    -o ${output_dir} \
    2>${output_dir}/timing-ptgaul-polishing.txt
  echo "use extract-ptdna-of-ptgaul2 <species_folder> if not working"
}

copy-ptdna-of-ptgaul_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"

  # copy ptGAUL result
  echo "copy ${output_dir}/ptdna-ptgaul.fa"
  _outdir="${output_dir}/result_3000/flye_cpONT/ptdna"
  _arg_final_assembly="${_outdir}/pt.1.fa"
  cp -p ${_arg_final_assembly} ${output_dir}/ptdna-ptgaul.fa
}

taxon-assemble_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"
  # copy_data

  ${_polap_cmd} taxonomy assemble \
    -o ${output_dir} \
    -l ${long_sra}.fastq \
    -a ${short_sra}_1.fastq \
    -b ${short_sra}_2.fastq \
    --species "${species_name}"
}

taxon-sample1_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"
  local ingroup="${_ingroup["$1"]}"
  local outgroup="${_outgroup["$1"]}"
  local allgroup="${_allgroup["$1"]}"
  # copy_data

  ${_polap_cmd} taxonomy sample \
    --steps-include 1-4 \
    -v \
    -o ${output_dir} \
    --taxonomy-rank-ingroup ${ingroup} \
    --taxonomy-rank-allgroup ${allgroup} \
    --species "${species_name}"
}

archive_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"
  # copy_data

  rm -rf ${output_dir}-a
  ${_polap_cmd} disassemble archive \
    -o ${output_dir}
  tar zcf ${output_dir}-a.tar.gz ${output_dir}-a
}

get_genus_species_for() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"

  if [[ "${_local_host}" == "$(hostname)" ]]; then
    case "${output_dir}" in
    "Eucalyptus_pauciflora")
      scp lab01:$PWD/${output_dir}-a.tar.gz .
      ;;
    *)
      scp ${_host["${output_dir}"]}:$PWD/${output_dir}-a.tar.gz .
      ;;
    esac

    tar zxf ${output_dir}-a.tar.gz
    mv ${output_dir}-a ${output_dir}
  else
    echo "ERROR: run at the local host."
  fi
}

get_genus_species() {
  local output_dir="${1:-default}"

  if [[ "${output_dir}" == "default" ]]; then
    for _v1 in "${S[@]}"; do
      get_genus_species_for "${_v1}"
    done
  else
    get_genus_species_for "${output_dir}"
  fi
}

report_genus_species() {
  local output_dir="${1:-default}"

  if [[ "${output_dir}" == "default" ]]; then
    for _v1 in "${S[@]}"; do
      report_genus_species_for "${_v1}"
    done
  else
    report_genus_species_for "${output_dir}"
  fi
}

report_genus_species_for() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"
  # copy_data

  local i=0
  local n
  local p
  IFS=',' read -r -a extracted_array_n <<<"${_compare_n["$1"]}"
  IFS=',' read -r -a extracted_array_p <<<"${_compare_p["$1"]}"
  for n in "${extracted_array_n[@]}"; do
    for p in "${extracted_array_p[@]}"; do
      i=$((i + 1))
      j=$((i + 3))
      for k in {1..3}; do
        ${_polap_cmd} disassemble report ${k} \
          -o ${output_dir} \
          --disassemble-i $i
        ${_polap_cmd} disassemble report ${k} \
          -o ${output_dir} \
          --disassemble-i $j
      done
    done
  done
}

# run_spirodela_polyrhiza
# run_taraxacum_mongolicum
# run_trifolium_pratense
# run_salix_dunnii
# run_anthoceros_agrestis
# run_anthoceros_angustus
# run_brassica_rapa
# run_vigna_radiata
# run_macadamia_tetraphylla
# run_punica_granatum
# run_lolium_perenne
# run_carex_pseudochinensis

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

if [[ "${subcmd1}" == "help" ]]; then
  if [[ "${_arg2}" == "arg2" ]]; then
    subcmd1="help"
  else
    subcmd1="${_arg2}"
    _arg2="arg2"
  fi
fi

# Call common case first
common_handled=1
# _polap_lib_data-execute-common-subcommand "$subcmd1" "${_arg2}" "${_arg3}" "$opt_y_flag"
_polap_lib_data-execute-common-subcommand
common_handled=$?

# Main case statement
case "$subcmd1" in
data2csv)
  if [[ "${_arg2}" == arg2 ]]; then
    echo "Help: ${subcmd1} [data:${_polap_data_data}] [csv:${_polap_data_csv}] [I:0|N] [Seed:0|N] [Downsample:50|N]"
    echo "  $0 ${subcmd1}"
    echo "${help_message_convert_data}"
    _arg2="${_polap_data_data}"
  fi
  [[ "${_arg3}" == arg3 ]] && _arg3=""
  [[ "${_arg4}" == arg4 ]] && _arg4=""
  [[ "${_arg5}" == arg5 ]] && _arg5=""
  [[ "${_arg6}" == arg6 ]] && _arg6=""
  ${subcmd1}_genus_species "${_arg2}" "${_arg3}" "${_arg4}" "${_arg5}" "${_arg6}"
  ;;
system)
  ${subcmd1}_genus_species
  ;;
copy-figures)
  ${subcmd1}_genus_species
  ;;
'send-data-to' | \
  'mkdir' | \
  'taxon-known' | \
  'taxon-reference' | \
  'ptgaul-known' | \
  'ptgaul-reference' | \
  'bandage-ptgaul-known' | \
  'bandage-ptgaul-reference' | \
  'subsample' | \
  'write-config' | \
  'run' | \
  'taxon-sample1' | \
  'taxon-sample2' | \
  'taxon-geseq' | \
  'taxon-orthofinder' | \
  'taxon-phylogeny' | \
  'taxon-tree' | \
  'taxon-species' | \
  'geseq' | \
  'get-ptdna-from-ncbi' | \
  'copy-ptdna-of-ncbi-as-reference' | \
  'msbwt' | \  | \
  'extract-ptdna-of-ptgaul' | \
  'extract-ptdna-of-ptgaul2' | \
  'copy-ptdna-of-ptgaul' | \
  'compare' | \
  'archive' | \
  'get' | \
  'report' | \
  'table1' | \
  'table2' | \
  'suptable1' | \
  'supfigure1' | \
  'test')
  ${subcmd1}_genus_species "${_arg2}"
  ;;
'assemble1-subsample-for' | \
  'bandage1-subsample-for' | \
  'seeds-subsample-for' | \
  'assemble2-subsample-for' | \
  'bandage2-subsample-for' | \
  'polish-subsample' | \
  'test2')
  ${subcmd1}_genus_species_for "${_arg2}" "${_arg3}" "${_arg4}"
  ;;
assemble1-subsample)
  ${subcmd1}_genus_species "${_arg2}"
  ;;
'bandage1-subsample' | \
  'seeds-subsample' | \
  'assemble2-subsample' | \
  'bandage2-subsample' | \
  'polish-subsample' | \
  'test2')
  ${subcmd1}_genus_species "${_arg2}"
  ;;
'chain-subsample')
  assemble1-subsample_genus_species "${_arg2}"
  seeds-subsample_genus_species "${_arg2}"
  bandage1-subsample_genus_species "${_arg2}"
  assemble2-subsample_genus_species "${_arg2}"
  bandage2-subsample_genus_species "${_arg2}"
  ;;
"Anthoceros_agrestis")
  run_anthoceros_agrestis
  ;;
"Brassica_rapa")
  run_brassica_rapa
  ;;
"Vigna_radiata")
  run_vigna_radiata
  ;;
"Trifolium_pratense")
  run_trifolium_pratense
  ;;
"Taraxacum_mongolicum")
  run_taraxacum_mongolicum
  ;;
"Spirodela_polyrhiza")
  run_spirodela_polyrhiza
  ;;
"Salix_dunnii")
  run_salix_dunnii
  ;;
"Punica_granatum")
  run_punica_granatum
  ;;
"Macadamia_tetraphylla")
  run_macadamia_tetraphylla
  ;;
"Lolium_perenne")
  run_lolium_perenne
  ;;
"Anthoceros_angustus")
  run_anthoceros_angustus
  ;;
"Carex_pseudochinensis")
  run_carex_pseudochinensis
  ;;
  ##### INSERT_CASE_HERE #####
polap-batch)
  if [[ -z "${_arg2}" || "${_arg2}" == arg2 ]]; then
    echo "Help: ${subcmd1} <outdir> [inum:x|N]"
    echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
    _subcmd1_clean="${subcmd1//-/_}"
    declare -n ref="help_message_${_subcmd1_clean}"
    echo "$ref"
    exit 0
  fi
  [[ "${_arg3}" == arg3 ]] && _arg3=""
  ${subcmd1}_genus_species "${_arg2}" "${_arg3}"
  ;;
polap-batch-polish)
  if [[ -z "${_arg2}" || "${_arg2}" == arg2 ]]; then
    echo "Help: ${subcmd1} <outdir> [inum:0|N] [msbwt:x|R]"
    echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
    _subcmd1_clean="${subcmd1//-/_}"
    declare -n ref="help_message_${_subcmd1_clean}"
    echo "$ref"
    exit 0
  fi
  [[ "${_arg3}" == arg3 ]] && _arg3=""
  [[ "${_arg4}" == arg4 ]] && _arg4=""
  ${subcmd1}_genus_species "${_arg2}" "${_arg3}" "${_arg4}"
  ;;
polap-msbwt)
  if [[ -z "${_arg2}" || "${_arg2}" == arg2 ]]; then
    echo "Help: ${subcmd1} <outdir> [inum:0|N]"
    echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
    _subcmd1_clean="${subcmd1//-/_}"
    declare -n ref="help_message_${_subcmd1_clean}"
    echo "$ref"
    exit 0
  fi
  [[ "${_arg3}" == arg3 ]] && _arg3=""
  ${subcmd1}_genus_species "${_arg2}" "${_arg3}"
  ;;
polap-batch-assemble)
  if [[ -z "${_arg2}" || "${_arg2}" == arg2 ]]; then
    echo "Help: ${subcmd1} <outdir> [inum:0|N]"
    echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
    _subcmd1_clean="${subcmd1//-/_}"
    declare -n ref="help_message_${_subcmd1_clean}"
    echo "$ref"
    exit 0
  fi
  [[ "${_arg3}" == arg3 ]] && _arg3=""
  ${subcmd1}_genus_species "${_arg2}" "${_arg3}"
  ;;
polap-archive)
  if [[ -z "${_arg2}" || "${_arg2}" == arg2 ]]; then
    echo "Help: ${subcmd1} <outdir> [inum:0|N] [jnum:1|N]"
    echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
    _subcmd1_clean="${subcmd1//-/_}"
    declare -n ref="help_message_${_subcmd1_clean}"
    echo "$ref"
    exit 0
  fi
  [[ "${_arg3}" == arg3 ]] && _arg3=""
  [[ "${_arg4}" == arg4 ]] && _arg4=""
  ${subcmd1}_genus_species "${_arg2}" "${_arg3}" "${_arg4}"
  ;;
polap-polish)
  if [[ -z "${_arg2}" || "${_arg2}" == arg2 ]]; then
    echo "Help: ${subcmd1} <outdir> [inum:0|N] [jnum:1|N]"
    echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
    _subcmd1_clean="${subcmd1//-/_}"
    declare -n ref="help_message_${_subcmd1_clean}"
    echo "$ref"
    exit 0
  fi
  [[ "${_arg3}" == arg3 ]] && _arg3=""
  [[ "${_arg4}" == arg4 ]] && _arg4=""
  ${subcmd1}_genus_species "${_arg2}" "${_arg3}" "${_arg4}"
  ;;
polap-assemble2)
  if [[ -z "${_arg2}" || "${_arg2}" == arg2 ]]; then
    echo "Help: ${subcmd1} <outdir> [inum:0|N] [jnum:1|N]"
    echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
    _subcmd1_clean="${subcmd1//-/_}"
    declare -n ref="help_message_${_subcmd1_clean}"
    echo "$ref"
    exit 0
  fi
  [[ "${_arg3}" == arg3 ]] && _arg3=""
  [[ "${_arg4}" == arg4 ]] && _arg4=""
  ${subcmd1}_genus_species "${_arg2}" "${_arg3}" "${_arg4}"
  ;;
polap-seeds)
  if [[ -z "${_arg2}" || "${_arg2}" == arg2 ]]; then
    echo "Help: ${subcmd1} <outdir> [inum:0|N]"
    echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
    _subcmd1_clean="${subcmd1//-/_}"
    declare -n ref="help_message_${_subcmd1_clean}"
    echo "$ref"
    exit 0
  fi
  [[ "${_arg3}" == arg3 ]] && _arg3=""
  ${subcmd1}_genus_species "${_arg2}" "${_arg3}"
  ;;
polap-assemble1)
  if [[ -z "${_arg2}" || "${_arg2}" == arg2 ]]; then
    echo "Help: ${subcmd1} <outdir> [inum:0|N] [coverage:50|N]"
    echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
    _subcmd1_clean="${subcmd1//-/_}"
    declare -n ref="help_message_${_subcmd1_clean}"
    echo "$ref"
    exit 0
  fi
  [[ "${_arg3}" == arg3 ]] && _arg3=""
  [[ "${_arg4}" == arg4 ]] && _arg4=""
  ${subcmd1}_genus_species "${_arg2}" "${_arg3}" "${_arg4}"
  ;;
polap-clean)
  if [[ -z "${_arg2}" || "${_arg2}" == arg2 ]]; then
    echo "Help: ${subcmd1} <outdir> [inum:0|N]"
    echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
    _subcmd1_clean="${subcmd1//-/_}"
    declare -n ref="help_message_${_subcmd1_clean}"
    echo "$ref"
    exit 0
  fi
  [[ "${_arg3}" == arg3 ]] && _arg3=""
  ${subcmd1}_genus_species "${_arg2}" "${_arg3}"
  ;;
get-data)
  if [[ -z "${_arg2}" || "${_arg2}" == arg2 ]]; then
    echo "Help: ${subcmd1} <outdir>"
    echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
    _subcmd1_clean="${subcmd1//-/_}"
    declare -n ref="help_message_${_subcmd1_clean}"
    echo "$ref"
    exit 0
  fi
  ${subcmd1}_genus_species "${_arg2}"
  ;;
pmat-suppfigure3)
  if [[ -z "${_arg2}" || "${_arg2}" == arg2 ]]; then
    echo "Help: ${subcmd1} [inum:0|N] [bandage:no|yes] [target:${_brg_default_target_dir}]"
    echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
    _subcmd1_clean="${subcmd1//-/_}"
    declare -n ref="help_message_${_subcmd1_clean}"
    echo "$ref"
    exit 0
  fi
  [[ "${_arg3}" == arg3 ]] && _arg3=""
  [[ "${_arg4}" == arg4 ]] && _arg4=""
  ${subcmd1}_genus_species "${_arg2}" "${_arg3}" "${_arg4}"
  ;;
maintable1)
  if [[ "${_arg2}" == arg2 ]]; then
    echo "Help: ${subcmd1} <all|outdir> [analysis:.|t1|string] [table:1|N] [targe_dir]"
    echo "  $0 $subcmd1 all"
    echo "  $0 $subcmd1 Eucalyptus_pauciflora 2"
    _subcmd1_clean="${subcmd1//-/_}"
    declare -n ref="help_message_${_subcmd1_clean}"
    echo "$ref"
    exit 0
  fi
  [[ "${_arg3}" == arg3 ]] && _arg3="."
  [[ "${_arg4}" == arg4 ]] && _arg4="1"
  [[ "${_arg5}" == arg5 ]] && _arg5="${_brg_default_target_dir}"
  ${subcmd1}_genus_species "${_arg2}" "${_arg3}" "${_arg4}" "${_arg5}"
  ;;
supptable1)
  if [[ "${_arg2}" == arg2 ]]; then
    echo "Help: ${subcmd1} <all|outdir> [table:1|N] [targe_dir]"
    echo "  $0 $subcmd1 all"
    echo "  $0 $subcmd1 Eucalyptus_pauciflora 2"
    _subcmd1_clean="${subcmd1//-/_}"
    declare -n ref="help_message_${_subcmd1_clean}"
    echo "$ref"
    exit 0
  fi
  [[ "${_arg3}" == arg3 ]] && _arg3="1"
  [[ "${_arg4}" == arg4 ]] && _arg4="${_brg_default_target_dir}"
  ${subcmd1}_genus_species "${_arg2}" "${_arg3}" "${_arg4}"
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
"test-polap" | "test")
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
"delete-links")
  find . -type l -delete
  ;;
"zip")
  for i in "${S[@]}"; do
    echo "zipping $i ..."
    tar zcf ${i}.tar.gz ${i}
  done
  ;;
"mkdir")
  for i in "${S[@]}"; do
    echo "creating folder $i ..."
    mkdir ${i}
  done
  ;;
"rm")
  for i in "${S[@]}"; do
    echo "deleting folder $i ..."
    rm -rf ${i}
  done
  ;;
"scopy-fastq")
  for i in "${S[@]}"; do
    echo "scure copying $i ..."
    scp -p ${_media_dir}/$i/*.fastq $2:$PWD/$i/
  done
  ;;
"copy-fastq")
  for i in "${S[@]}"; do
    echo "copying $i ..."
    cp -p ${_media_dir}/$i/*.fastq $i/
  done
  ;;
"link-fastq")
  for i in "${S[@]}"; do
    cd $i
    for i in ${_media_dir}/$i/*.fastq; do
      ln -s $i
    done
    cd -
  done
  ;;
"plot")
  pandoc plot.md -o supp-figure-s11-s34.pdf
  pandoc plot.md -o supp-figure-s11-s34.docx
  ;;
"table")
  bash table1.sh >table1.tsv
  pandoc table1.tsv -f tsv -t pdf -o table1.pdf -V geometry:landscape
  pandoc table1.tsv -f tsv -t docx -o table1.docx
  ;;
"sync")
  for i in "${S[@]}"; do
    V1="${i}"
    cd $V1
    if [[ -v _host["${V1}"] ]]; then
      mkdir -p o/${value}/01-contig o/${value}/06-summary o/${value}/07-plot
      rsync -aPh ${_host["${V1}"]}:$PWD/o/${value}/01-contig/ o/${value}/01-contig/
      rsync -aPh ${_host["${V1}"]}:$PWD/o/${value}/06-summary/ o/${value}/06-summary/
      rsync -aPh ${_host["${V1}"]}:$PWD/o/${value}/07-plot/ o/${value}/07-plot/
      rsync -aPh ${_host["${V1}"]}:$PWD/o/${value}/30-contigger o/${value}/
      rsync -aPh ${_host["${V1}"]}:$PWD/o/${value}/*annotation* o/${value}/
      rsync -aPh ${_host["${V1}"]}:$PWD/o/${value}/assembly.fasta o/${value}/
      rsync -aPh ${_host["${V1}"]}:$PWD/o/${value}/assembly_graph.gfa o/${value}/
      rsync -aPh ${_host["${V1}"]}:$PWD/o/00-bioproject o/
      # rsync -aPh o/${value}/mt.0.fasta ${_host["${V1}"]}:$PWD/o/${value}/
      rsync -aPh ${_host["${V1}"]}:$PWD/o/polap.log o/
      # cat o/${value}/07-plot/*.md >>../plot.md
    else
      echo "No such host for $V1"
    fi
    cd -
  done
  ;;
"package")
  for i in "${S[@]}"; do
    V1="${i}"
    rm -rf $V1
    mkdir $V1
    ${_polap_cmd} package -o ../revision/$V1/o --archive $V1/o
    if [[ "$V1" == "Salix_dunnii" ]]; then
      ${_polap_cmd} -o $V1/o test-reads report ptgaul
      ${_polap_cmd} -o $V1/o test-reads report intra
    elif [[ "$V1" == "Lolium_perenne" ]]; then
      ${_polap_cmd} -o $V1/o test-reads report ptgaul --report-x 5000,7000,9000,11000,13000,15000,17000
      ${_polap_cmd} -o $V1/o test-reads report polap --report-x 5000,7000,9000,11000,13000,15000,17000
    else
      ${_polap_cmd} -o $V1/o test-reads report ptgaul
      ${_polap_cmd} -o $V1/o test-reads report polap
    fi
  done
  ;;
sample-csv)
  if [[ "${_arg2}" == arg2 ]]; then
    echo "Help: ${subcmd1} <csv:${_polap_data_csv}> [number:1|2|7|test|all|each] [force:off|on] [inum:0|N]"
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
      head -1 ${_POLAPLIB_DIR}/${_polap_data_csv} >"${_arg2}"
      grep Spirodela_polyrhiza ${_POLAPLIB_DIR}/${_polap_data_csv} >>"${_arg2}"
    elif [[ "${_arg3}" == "2" ]]; then
      head -1 ${_POLAPLIB_DIR}/${_polap_data_csv} >"${_arg2}"
      grep Spirodela_polyrhiza ${_POLAPLIB_DIR}/${_polap_data_csv} >>"${_arg2}"
      grep Eucalyptus_pauciflora ${_POLAPLIB_DIR}/${_polap_data_csv} >>"${_arg2}"
    elif [[ "${_arg3}" == "test" ]]; then
      head -1 ${_POLAPLIB_DIR}/${_polap_data_csv} >"${_arg2}"
      grep test ${_POLAPLIB_DIR}/${_polap_data_csv} >>"${_arg2}"
    elif [[ "${_arg3}" == "7" ]]; then
      head -1 ${_POLAPLIB_DIR}/${_polap_data_csv} >"${_arg2}"
      for item in $(printf "%s\n" "${S7[@]}" | sort); do
        grep $item ${_POLAPLIB_DIR}/${_polap_data_csv} |
          grep -v test >>"${_arg2}"
      done
    elif [[ "${_arg3}" == "each" ]]; then
      extract_and_replace_suffix "${_POLAPLIB_DIR}"/${_polap_data_csv} "-${_arg5}" "${_arg2}"
    elif [[ "${_arg3}" == "all" ]]; then
      grep -v test ${_POLAPLIB_DIR}/${_polap_data_csv} >"${_arg2}"
    fi
    echo "create CSV: ${_arg2}"
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
