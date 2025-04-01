#!/usr/bin/bash

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" || {
  echo "Couldn't determine the script's running directory, which probably matters, bailing out" >&2
  exit 2
}
LIB_DIR="${_POLAPLIB_DIR}/lib"

source "${LIB_DIR}/polap-lib-timing.sh"

################################################################################
# local functions

# Read the config files
read-a-tsv-file-into-associative-arrays() {
  # Define input CSV file
  csv_file="src/polap-data-v4.csv"

  # Read the CSV file (skip header)
  while IFS=$',' read -r species long short host inref min_read range random_seed; do
    # Skip header line
    [[ "$species" == "species" ]] && continue

    # Store in associative arrays
    _long["$species"]="$long"
    _short["$species"]="$short"
    _host["$species"]="$host"
    _min_read["$species"]="$min_read"
    _range["$species"]="$range"
    _random_seed["$species"]="$random_seed"
  done <"$csv_file"
}

# Loop through the associative array
default_output_dir() {
  local _output_dir=""

  for key in "${!_host[@]}"; do
    if [[ "${_host[$key]}" == "$(hostname)" ]]; then
      _output_dir="$key"
      break # Exit loop after finding the first match
    fi
  done
  echo "${_output_dir}"
}

S=(
  'Macadamia_tetraphylla'
  'Punica_granatum'
  'Lolium_perenne'
)

_polap_subcmd=(
  'test'
  'seeds'
  'map'
  'reads'
  'dflye'
  'directional'
  'archive'
  'report'
  'table1'
  'figure1'
)

_local_host="thorne"
if [[ -d "src" ]]; then
  _polap_cmd="src/polap.sh"
else
  _polap_cmd="polap"
fi
_polap_version="0.4.5.1"
_media_dir="/media/h2/sra" # for v2 and v3
_media_dir="/media/h1/sra" # for v1 and v4

help_message=$(
  cat <<HEREDOC
# Polap data analysis for the version 0.4 or directional
#
# Edit two variables:
# _polap_cmd=${_polap_cmd}
# _media_dir=${_media_dir}
# _local_host=${_local_host}
# host=$(hostname)
#
# subcommand can be replaced with the leading number
0. <species_folder>: To verify that this script functions correctly, execute the initial subcommand.
1. <species_folder> seeds:
2. <species_folder> map:
3. <species_folder> reads:
4. <species_folder> dflye:
5. <species_folder> directional:
6. <species_folder> archive:
7. <species_folder> report:
8. <species_folder> table1:
9. <species_folder> figure1:
cmd subcmd:

HEREDOC
)

_arg1=${1:-"default_value"}
# Check if the species folder is provided
if [[ "${_arg1}" == "default_value" ]]; then
  echo "Usage: $0 <species_folder> <subcommand>"
  echo "${help_message}"
  exit 1
fi

# Input parameter
subarg1="${1%/}"
subcmd1="${2:-0}"

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

declare -A _status
declare -A _host
declare -A _long
declare -A _short
declare -A _random_seed
declare -A _ssh
declare -A _range

set +u

read-a-tsv-file-into-associative-arrays

# Order the species names in alphabetical order
# keys_array=("${!_long[@]}")
keys_array=($(for key in "${!_long[@]}"; do echo "$key"; done | sort))

################################################################################
# Part of genus_species
#
test_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"
  local random_seed="${_random_seed["$1"]}"

  echo output: $output_dir
  echo host: $(hostname)
  echo species: $species_name
  echo random seed: $random_seed

  if [[ -s "${long_sra}.fastq" ]]; then
    echo long: $long_sra
  else
    echo long: no such fastq file: $long_sra
  fi
  if [[ -s "${short_sra}_1.fastq" ]] &&
    [[ -s "${short_sra}_2.fastq" ]]; then
    echo short: $short_sra
  else
    echo short: no such fastq file: $short_sra
  fi
}

# whole-genome assembly
# assemble1
# mt.contig.name-1

# create seed contigs
#
# /home/goshng/all/polap/dflye1/Lolium_perenne/o/2
# └── 01-contig
#     ├── contig1.fa
#     ├── contig1.name.txt
#     ├── contig2.fa
#     └── contig2.name.txt
#
# 1 directory, 4 files
seeds_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"
  local random_seed="${_random_seed["$1"]}"

  # mt.contig.name-1 with plus/minus signed edge_<number>[+-]
  ${_polap_cmd} directional-prepare-seeds \
    -o ${output_dir}/o \
    -i 1 -j 2
}

# map the long-read data on the seed contigs
#
# input:
# o/lk.fq.gz
# o/nk.fq.gz
# seed contigs
#
# output:
#
# /home/goshng/all/polap/dflye1/Lolium_perenne/o/2
# └── 01-contig
#     ├── contig1.fa
#     ├── contig1.name.txt
#     ├── contig1.paf
#     ├── contig1.tab
#     ├── contig1_total_length.txt
#     ├── contig2.fa
#     ├── contig2.name.txt
#     ├── contig2.paf
#     └── contig2.tab
#
# 1 directory, 9 files
map_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"
  local random_seed="${_random_seed["$1"]}"

  # o/lk.fq.gz
  # o/nk.fq.gz
  ${_polap_cmd} directional-map-reads \
    -o ${output_dir}/o \
    -i 1 -j 2
}

# select reads mapped on the reference
#
# /home/goshng/all/polap/dflye1/Lolium_perenne/o/2
# ├── 02-reads
# │   └── ptgaul
# │       ├── 0
# │       │   ├── ptgaul.contig1.forward.names.txt
# │       │   └── ptgaul.contig.forward.names.txt
# │       ├── 1
# │       │   ├── ptgaul.contig1.forward.names.txt
# │       │   └── ptgaul.contig.forward.names.txt
# │       └── 2
# │           ├── ptgaul.contig1.forward.names.txt
# │           └── ptgaul.contig.forward.names.txt
# ├── 03-seeds
# │   └── ptgaul
# │       ├── 0.fq.gz
# │       ├── 1.fq.gz
# │       └── 2.fq.gz
# ├── 04-subsample
# │   └── ptgaul
# ├── 05-flye
# │   └── ptgaul
# ├── 06-summary
# │   └── ptgaul
# └── 07-plot
#     └── ptgaul
#
reads_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"
  local random_seed="${_random_seed["$1"]}"

  # -s 3000,6000,2
  # -s 3000,27000,9
  # -s 30000,42000,5
  ${_polap_cmd} directional-test-reads \
    -o ${output_dir}/o \
    -i 1 -j 2 \
    -s 33000,45000,5
}

# assemble the directinoal reads using dflye
dflye_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"
  local random_seed="${_random_seed["$1"]}"

  # mt_fastq=/home/goshng/all/polap/dflye1/Lolium_perenne/o/2/03-seeds/ptgaul/2.fq.gz
  # bin/dflye \
  #   --nano-raw ${mt_fastq} \
  #   --debug \
  #   --stop-after contigger \
  #   --asm-coverage 30 \
  #   --directional-reads \
  #   -g 800000 -o $out_dir -t 56 -m 10000

  ${_polap_cmd} dflye \
    -o ${output_dir}/o \
    -i 1 -j 2
}

# the first step is a typical polap-assemble1 and assemble2 with seed contig selection.
# now, we use the organelle genome assembly as a seed to apply this directional method.
#
# this is later step.
#
directional_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"
  local random_seed="${_random_seed["$1"]}"

  IFS=':' read -r -a extracted_array_n <<<"${_range["$1"]}"
  local extracted_range="${_range["$1"]//:/,}"
  local extracted_min_read="${_min_read["$1"]}"

  ${_polap_cmd} directional \
    -v -v \
    --select-read-range "${extracted_range}" \
    -o ${output_dir}/o \
    --directional-i 1 \
    --min-read-length "${extracted_min_read}" \
    -i 1 -j 3
}

mkdir_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"

  local long_data="${_media_dir}/${long_sra}.fastq.tar.gz"
  local short_data="${_media_dir}/${short_sra}.fastq.tar.gz"

  echo "create ${output_dir} ..."
  mkdir -p "${output_dir}"
  if tar -zxf "$(basename ${long_data})"; then
    rm "$(basename ${long_data})"
    echo "Extraction successful: ${long_sra}. Archive deleted."
  fi
  if tar -zxf "$(basename ${short_data})"; then
    rm "$(basename ${short_data})"
    echo "Extraction successful: ${short_data} Archive deleted."
  fi

  get-ptdna-from-ncbi_genus_species "${output_dir}"
}

clean_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"

  local long_data="${_media_dir}/${long_sra}.fastq.tar.gz"
  local short_data="${_media_dir}/${short_sra}.fastq.tar.gz"

  # Check if the directory exists
  if [[ -d "$output_dir" ]]; then
    echo "Directory '$output_dir' exists."

    # Ask for user confirmation before deleting
    read -p "Do you want to delete it? (y/N): " confirm
    case "$confirm" in
    [yY] | [yY][eE][sS])

      echo "cleaning ${output_dir} ..."
      rm -rf "${output_dir}"
      rm -rf "${output_dir}-a"
      rm -f "${output_dir}-a.tar.gz"
      rm -f "${long_sra}.fastq"
      rm -f "${short_sra}"_?.fastq
      ;;
    *)
      echo "Deletion canceled."
      ;;
    esac
  else
    echo "Directory '$output_dir' does not exist."
  fi

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
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"

  command time -v bash src/polap-ptGAUL1.sh \
    -o ${output_dir}-ptgaul \
    -r ${output_dir}/ptdna-reference.fa \
    -g "${_ptgaul_genomesize["$1"]}" \
    -l ${long_sra}.fastq \
    -t 24 \
    2>${output_dir}/timing-ptgaul.txt

  mv ${output_dir}-ptgaul/result_3000 ${output_dir}/
  rm -rf "${output_dir}-ptgaul"

  msbwt_genus_species "${output_dir}"
}

msbwt_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"

  echo "This polishing preparation takes long ..."
  command time -v ${_polap_cmd} prepare-polishing \
    -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
    -o ${output_dir} \
    2>${output_dir}/timing-prepare-polishing.txt
}

# Extraction of ptDNA from the assembly: a first try
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

# Extraction of ptDNA from the assembly: another try
extract-ptdna-of-ptgaul2_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"

  # extract ptGAUL result
  echo "extract ptDNA from the ptGAUL result with fmlrc polishing"
  ${_polap_cmd} disassemble ptgaul 2 \
    -v -v -v \
    -o ${output_dir}
  command time -v ${_polap_cmd} disassemble ptgaul 3 \
    -v -v -v \
    -o ${output_dir} \
    2>${output_dir}/timing-ptgaul-polishing.txt

  _outdir="${output_dir}/result_3000/flye_cpONT/ptdna"
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
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"

  # copy ptGAUL result
  echo "copy ${output_dir}/ptdna-ptgaul.fa"
  _outdir="${output_dir}/result_3000/flye_cpONT/ptdna"
  _arg_final_assembly="${_outdir}/pt.1.fa"
  cp -p ${_arg_final_assembly} ${output_dir}/ptdna-ptgaul.fa
}

# Case of the compare menu
# --disassemble-c
# --no-disassemble-align-reference
compare_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"
  local random_seed="${_random_seed["$1"]}"
  # copy_data

  local i=0
  local n
  local p
  IFS=':' read -r -a extracted_array_n <<<"${_compare_n["$1"]}"
  IFS=':' read -r -a extracted_array_p <<<"${_compare_p["$1"]}"
  local extracted_r="${_compare_r["$1"]}"
  for n in "${extracted_array_n[@]}"; do
    for p in "${extracted_array_p[@]}"; do
      i=$((i + 1))
      echo "($i) n=$n, p=$p"

      local _d_i
      _d_i="compare-$i"

      command time -v ${_polap_cmd} disassemble \
        -o ${output_dir} \
        -l ${long_sra}.fastq \
        -a ${short_sra}_1.fastq \
        -b ${short_sra}_2.fastq \
        --disassemble-c ${output_dir}/ptdna-ptgaul.fa \
        --disassemble-i "${_d_i}" \
        --disassemble-n $n \
        --disassemble-p $p \
        --disassemble-r ${extracted_r} \
        --disassemble-alpha 1.0 \
        --random-seed "${random_seed}" \
        2>${output_dir}/timing-${_d_i}.txt
    done
  done

}

# Case of the check menu
# --disassemble-c
# --disassemble-align-reference
# --disassemble-simple-polishing
check_genus_species() {
  local output_dir="$1"
  local simple_polishing="${2:-default}"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"
  local random_seed="${_random_seed["$1"]}"
  # copy_data

  local i=0
  local n
  local p
  IFS=':' read -r -a extracted_array_n <<<"${_compare_n["$1"]}"
  IFS=':' read -r -a extracted_array_p <<<"${_compare_p["$1"]}"
  local extracted_r="${_compare_r["$1"]}"
  for n in "${extracted_array_n[@]}"; do
    for p in "${extracted_array_p[@]}"; do
      i=$((i + 1))
      echo "($i) n=$n, p=$p"

      local _d_i="infer-$i"
      local _x_i="check-$i"
      local _s_i="subsample-polish"
      local _stages="--stages-include 3"
      if [[ "${simple_polishing}" != "default" ]]; then
        _s_i="simple-polish"
      fi

      # "${_stages}" \
      command time -v ${_polap_cmd} disassemble \
        ${_stages} \
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
        --disassemble-alpha 1.0 \
        --random-seed "${random_seed}" \
        2>"${output_dir}/timing-${_x_i}-${_s_i}.txt"

      if [[ -d "${output_dir}/disassemble/${_d_i}/3" ]]; then
        rm -rf "${output_dir}/disassemble/${_d_i}/3-check"
        mv "${output_dir}/disassemble/${_d_i}/3" \
          "${output_dir}/disassemble/${_d_i}/3-check"
      fi

      # compare the results
      mkdir -p ${output_dir}/mauve
      mkdir -p ${output_dir}/blast
      if [[ -s "${output_dir}/disassemble/infer-1/pt.subsample-polishing.reference.aligned.1.fa" ]]; then
        ${_polap_cmd} mauve-mtdna -a "${output_dir}/ptdna-ptgaul.fa" -b "${output_dir}/disassemble/infer-1/pt.subsample-polishing.reference.aligned.1.fa" -o "${output_dir}/mauve" >"${output_dir}/mauve/log.txt"
        echo "see ${output_dir}/mauve/log.txt"
        cat "${output_dir}/mauve/log.txt"
        ${_polap_cmd} compare2ptdna -a "${output_dir}/ptdna-ptgaul.fa" -b "${output_dir}/disassemble/infer-1/pt.subsample-polishing.reference.aligned.1.fa" -o "${output_dir}/blast"
        echo "see ${output_dir}/blast/pident.txt"
        cat "${output_dir}/blast/pident.txt"
      else
        echo "ERROR: no such file: ${output_dir}/disassemble/infer-1/pt.subsample-polishing.reference.aligned.1.fa"
      fi

    done
  done

}

mauve_genus_species() {
  local output_dir="$1"
  local simple_polishing="${2:-default}"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"
  local random_seed="${_random_seed["$1"]}"
  # copy_data

  local i=0
  local n
  local p
  IFS=':' read -r -a extracted_array_n <<<"${_compare_n["$1"]}"
  IFS=':' read -r -a extracted_array_p <<<"${_compare_p["$1"]}"
  local extracted_r="${_compare_r["$1"]}"
  for n in "${extracted_array_n[@]}"; do
    for p in "${extracted_array_p[@]}"; do
      i=$((i + 1))
      echo "($i) n=$n, p=$p"

      mkdir -p ${output_dir}/mauve
      mkdir -p ${output_dir}/blast
      if [[ -s "${output_dir}/disassemble/infer-1/pt.subsample-polishing.reference.aligned.1.fa" ]]; then
        ${_polap_cmd} mauve-mtdna -a "${output_dir}/ptdna-ptgaul.fa" -b "${output_dir}/disassemble/infer-1/pt.subsample-polishing.reference.aligned.1.fa" -o "${output_dir}/mauve" >"${output_dir}/mauve/log.txt"
        echo "see ${output_dir}/mauve/log.txt"
        cat "${output_dir}/mauve/log.txt"
        ${_polap_cmd} compare2ptdna -a "${output_dir}/ptdna-ptgaul.fa" -b "${output_dir}/disassemble/infer-1/pt.subsample-polishing.reference.aligned.1.fa" -o "${output_dir}/blast"
        echo "see ${output_dir}/blast/pident.txt"
        cat "${output_dir}/blast/pident.txt"
      else
        echo "ERROR: no such file: ${output_dir}/disassemble/infer-1/pt.subsample-polishing.reference.aligned.1.fa"
      fi

    done
  done

}

# Case of the infer menu
# no --disassemble-c
infer_genus_species() {
  local output_dir="$1"
  local simple_polishing="${2:-default}"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"
  local random_seed="${_random_seed["$1"]}"
  # copy_data

  local i=0
  local n
  local p
  IFS=':' read -r -a extracted_array_n <<<"${_compare_n["$1"]}"
  IFS=':' read -r -a extracted_array_p <<<"${_compare_p["$1"]}"
  local extracted_r="${_compare_r["$1"]}"
  for n in "${extracted_array_n[@]}"; do
    for p in "${extracted_array_p[@]}"; do
      i=$((i + 1))
      echo "($i) n=$n, p=$p"

      local _d_i="infer-$i"
      local _x_i="infer-$i"
      local _s_i="subsample-polish"
      local _stages="--stages-include 1-3"
      if [[ "${simple_polishing}" != "default" ]]; then
        _stages="--stages-include 3"
        _s_i="simple-polish"
      fi

      command time -v ${_polap_cmd} disassemble \
        ${_stages} \
        -o ${output_dir} \
        -l ${long_sra}.fastq \
        -a ${short_sra}_1.fastq \
        -b ${short_sra}_2.fastq \
        ${simple_polishing} \
        --disassemble-i "${_d_i}" \
        --disassemble-n $n \
        --disassemble-p $p \
        --disassemble-r ${extracted_r} \
        --disassemble-alpha 1.0 \
        --random-seed "${random_seed}" \
        2>"${output_dir}/timing-${_x_i}-${_s_i}.txt"

      if [[ -d "${output_dir}/disassemble/${_d_i}/3" ]]; then
        rm -rf "${output_dir}/disassemble/${_d_i}/3-infer"
        mv "${output_dir}/disassemble/${_d_i}/3" \
          "${output_dir}/disassemble/${_d_i}/3-infer"
      fi
    done
  done
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

# the rest of them are test cases
infer12_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"
  # copy_data

  local i=3
  local n
  local p
  IFS=':' read -r -a extracted_array_n <<<"${_compare_n["$1"]}"
  IFS=':' read -r -a extracted_array_p <<<"${_compare_p["$1"]}"
  local extracted_r="${_compare_r["$1"]}"
  for n in "${extracted_array_n[@]}"; do
    for p in "${extracted_array_p[@]}"; do
      i=$((i + 1))
      echo "($i) n=$n, p=$p"

      command time -v ${_polap_cmd} disassemble \
        --stages-include 1-4 \
        -o ${output_dir} \
        -l ${long_sra}.fastq \
        -a ${short_sra}_1.fastq \
        -b ${short_sra}_2.fastq \
        --disassemble-i $i \
        --disassemble-n $n \
        --disassemble-p $p \
        --disassemble-r ${extracted_r} \
        --disassemble-alpha 1.0 \
        2>${output_dir}/timing-infer12-${i}.txt
    done
  done
}

infer1_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"
  # copy_data

  local i=3
  local n
  local p
  IFS=':' read -r -a extracted_array_n <<<"${_compare_n["$1"]}"
  IFS=':' read -r -a extracted_array_p <<<"${_compare_p["$1"]}"
  local extracted_r="${_compare_r["$1"]}"
  for n in "${extracted_array_n[@]}"; do
    for p in "${extracted_array_p[@]}"; do
      i=$((i + 1))
      echo "($i) n=$n, p=$p"

      command time -v ${_polap_cmd} disassemble stage2 \
        --stages-include 1-2 \
        -o ${output_dir} \
        -l ${long_sra}.fastq \
        -a ${short_sra}_1.fastq \
        -b ${short_sra}_2.fastq \
        --disassemble-i $i \
        --disassemble-n $n \
        --disassemble-p $p \
        --disassemble-r ${extracted_r} \
        --disassemble-alpha 1.0 \
        2>${output_dir}/timing-infer2only-${i}.txt
    done
  done
}

infer2_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"
  # copy_data

  local i=3
  local n
  local p
  IFS=':' read -r -a extracted_array_n <<<"${_compare_n["$1"]}"
  IFS=':' read -r -a extracted_array_p <<<"${_compare_p["$1"]}"
  local extracted_r="${_compare_r["$1"]}"
  for n in "${extracted_array_n[@]}"; do
    for p in "${extracted_array_p[@]}"; do
      i=$((i + 1))
      echo "($i) n=$n, p=$p"

      command time -v ${_polap_cmd} disassemble stage2 \
        --stages-include 3-4 \
        -o ${output_dir} \
        -l ${long_sra}.fastq \
        -a ${short_sra}_1.fastq \
        -b ${short_sra}_2.fastq \
        --disassemble-i $i \
        --disassemble-n $n \
        --disassemble-p $p \
        --disassemble-r ${extracted_r} \
        --disassemble-alpha 1.0 \
        2>${output_dir}/timing-infer2only-${i}.txt
    done
  done
}

# polish + stage 6 of disassemble
infer3_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"
  # copy_data

  local i=3
  local n
  local p
  IFS=':' read -r -a extracted_array_n <<<"${_polish_n["$1"]}"
  IFS=':' read -r -a extracted_array_p <<<"${_polish_p["$1"]}"
  local extracted_r="${_compare_r["$1"]}"
  for n in "${extracted_array_n[@]}"; do
    for p in "${extracted_array_p[@]}"; do
      i=$((i + 1))
      echo "($i) n=$n, p=$p"

      command time -v ${_polap_cmd} disassemble \
        --stages-include 5-6 \
        -o ${output_dir} \
        -l ${long_sra}.fastq \
        -a ${short_sra}_1.fastq \
        -b ${short_sra}_2.fastq \
        --disassemble-i $i \
        --disassemble-n $n \
        --disassemble-p $p \
        --disassemble-r ${extracted_r} \
        --disassemble-alpha 1.0 \
        2>${output_dir}/timing-infer3only-${i}.txt
    done
  done
}

# a simple polishing
#
# no need of a known ptDNA; use the first or the selected ptdna.0.fa
# 1. copy the selected j's ptdna.0.fa
# disassemble/4/2/j/52-mtdna/ptdna.0.fa -> disassemble/4/pt.0.fa
# 2. polish the pt.0.fa -> pt.1.fa
simple-polish_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"
  # copy_data

  local i=3
  local n
  local p
  IFS=':' read -r -a extracted_array_n <<<"${_compare_n["$1"]}"
  IFS=':' read -r -a extracted_array_p <<<"${_compare_p["$1"]}"
  for n in "${extracted_array_n[@]}"; do
    for p in "${extracted_array_p[@]}"; do
      i=$((i + 1))
      echo "($i) n=$n, p=$p"

      command time -v ${_polap_cmd} disassemble polishing \
        -o ${output_dir} \
        --disassemble-i $i \
        2>${output_dir}/timing-polishing-${i}.txt
    done
  done
}

# polish with incremental subsampling
#
# just the subsampling polish step
# it is step 5 of disassemble
subsample-polish_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"
  # copy_data

  local i=3
  local n
  local p
  IFS=':' read -r -a extracted_array_n <<<"${_polish_n["$1"]}"
  IFS=':' read -r -a extracted_array_p <<<"${_polish_p["$1"]}"
  for n in "${extracted_array_n[@]}"; do
    for p in "${extracted_array_p[@]}"; do
      i=$((i + 1))
      echo "($i) n=$n, p=$p"

      ${_polap_cmd} disassemble polish \
        -v -v \
        --steps-include 4,6,7,8 \
        -o ${output_dir} \
        -l ${long_sra}.fastq \
        -a ${short_sra}_1.fastq \
        -b ${short_sra}_2.fastq \
        --disassemble-i $i \
        --disassemble-n $n \
        --disassemble-p $p
      # 2>${output_dir}/timing-polish-${i}.txt
    done
  done
}

# WARNING: delete it?
disassemble-check_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"
  # copy_data

  local i=3
  local n
  local p
  IFS=':' read -r -a extracted_array_n <<<"${_compare_n["$1"]}"
  IFS=':' read -r -a extracted_array_p <<<"${_compare_p["$1"]}"
  for n in "${extracted_array_n[@]}"; do
    for p in "${extracted_array_p[@]}"; do
      i=$((i + 1))
      echo "($i) n=$n, p=$p"

      command time -v ${_polap_cmd} disassemble check \
        -o ${output_dir} \
        --disassemble-c ${output_dir}/ptdna-ptgaul.fa \
        --disassemble-i $i \
        2>${output_dir}/timing-check-${i}.txt
    done
  done
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

archive_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"
  # copy_data

  rm -rf "${output_dir}-a"
  rm -f "${output_dir}-a.tar.gz"
  ${_polap_cmd} disassemble archive \
    -o ${output_dir}
  # tar zcf "${output_dir}-a.tar.gz" "${output_dir}-a"
}

get_genus_species_for() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"

  if [[ "${_local_host}" == "$(hostname)" ]]; then

    # Define the directory name (can be passed as an argument)
    dir_to_check="${output_dir}"

    # Check if the directory exists
    if [[ -d "$dir_to_check" ]]; then
      echo "Directory '$dir_to_check' exists."

      # Ask for user confirmation before deleting
      read -p "Do you want to delete it? (y/N): " confirm
      case "$confirm" in
      [yY] | [yY][eE][sS])

        rm -f "${output_dir}-a.tar.gz"
        scp ${_ssh["${output_dir}"]}:$PWD/${output_dir}-a.tar.gz .
        # case "${output_dir}" in
        # "Eucalyptus_pauciflora")
        # 	scp lab01:$PWD/${output_dir}-a.tar.gz .
        # 	;;
        # *)
        # 	scp ${_ssh["${output_dir}"]}:$PWD/${output_dir}-a.tar.gz .
        # 	;;
        # esac

        tar zxf ${output_dir}-a.tar.gz
        mv "$dir_to_check/0-bioproject" "${output_dir}-a/"
        mv "$dir_to_check/bioproject.txt" "${output_dir}-a/"
        rm -rf "$dir_to_check"
        echo "Directory '$dir_to_check' deleted."
        mv ${output_dir}-a ${output_dir}
        echo "${output_dir} is downloaded."
        ;;
      *)
        echo "Deletion canceled."
        ;;
      esac
    else
      echo "Directory '$dir_to_check' does not exist."
    fi
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
  IFS=':' read -r -a extracted_array_n <<<"${_compare_n["$1"]}"
  IFS=':' read -r -a extracted_array_p <<<"${_compare_p["$1"]}"
  for n in "${extracted_array_n[@]}"; do
    for p in "${extracted_array_p[@]}"; do
      i=$((i + 1))
      for k in {1..2}; do
        ${_polap_cmd} disassemble report ${k} infer \
          -o ${output_dir} \
          --disassemble-i infer-$i
        ${_polap_cmd} disassemble report ${k} \
          -o ${output_dir} \
          --disassemble-i compare-$i
      done
    done
  done
}

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
# read -r I P N <<< "$result"
parse_params() {
  local file="$1" # Input file
  local I P N     # Declare variables for the parameters

  # Read the file and extract the values
  while IFS=": " read -r key value; do
    case "$key" in
    "I") I="$value" ;;
    "P") P="$value" ;;
    "N") N="$value" ;;
    esac
  done <"$file"

  # Print the variables (return as output)
  echo "$I $P $N"
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

table1_genus_species() {
  local _table1_file="table1.tsv"
  local _v1
  local _logfile
  local _readmefile
  local _species
  local _l_sra
  local _s_sra
  local _l_sra_size
  local _s_sra_size
  local _l_sra_size_gb
  local _s_sra_size_gb
  local _memory_gb
  local _total_hours
  local _memory_gb_ptgaul
  local _total_hours_ptgaul
  local _I
  local _P
  local _N
  local j

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "Species" \
    "P" \
    "L_SRA" \
    "L_size" \
    "S_SRA" \
    "S_size" \
    "NCBI ptDNA" \
    "Mode" \
    "SD" \
    "M_r" \
    "M_f" \
    "T" \
    >"${_table1_file}"

  for _v1 in "${Stable1[@]}"; do
    _logfile=${_v1}/polap.log
    _readmefile=${_v1}/README
    _species="${_v1/_/ }"
    _l_sra=$(awk '/long-read/ {match($0, /([A-Z]RR[0-9]+)/, arr); print arr[0]}' "$_logfile" | sort -u | grep -v '^$')
    _l_sra_size=$(<${_v1}/long_total_length.txt)
    _l_sra_size_gb=$(_polap_utility_convert_bp "${_l_sra_size}")
    _s_sra=$(awk '/short-read1/ {match($0, /([A-Z]RR[0-9]+)/, arr); print arr[0]}' "$_logfile" | sort -u | grep -v '^$')
    _s_sra_size=$(<${_v1}/short_total_length.txt)
    _s_sra_size_gb=$(_polap_utility_convert_bp "${_s_sra_size}")
    _known_mtdna=$(grep 'NCBI accession:' ${_logfile} | cut -d: -f4 | tail -n 1)

    if [[ -v _host["${_v1}"] ]]; then
      echo host:${_host["${_v1}"]}
    else
      echo "No such host for $_v1"
    fi

    local _ptdna_ptgaul=${_v1}/ptdna-ptgaul.fa
    local _ptdna_reference=${_v1}/ptdna-reference.fa

    for j in {1..1}; do
      local _disassemble_index="compare-${j}"
      local _summary1_ordered_txt=${_v1}/disassemble/${_disassemble_index}/1/summary1-ordered.txt

      # Extract mode value
      # Extract SD value
      # Extract the first index value
      local _mode=$(grep "^#mode:" "$_summary1_ordered_txt" | awk '{print $2}')
      local _sd=$(grep "^#sd:" "$_summary1_ordered_txt" | awk '{print $2}')
      local _first_index=$(grep "^#index:" "$_summary1_ordered_txt" | awk 'NR==1 {print $2}')

      local _params_txt=${_v1}/disassemble/${_disassemble_index}/params.txt
      _ipn=$(parse_params "${_params_txt}")
      read -r _I _P _N <<<"$_ipn"

      # read -r _memory_gb _total_hours < <(parse_timing "${_v1}" "${j}")
      read -r _memory_gb _total_hours < <(_polap_lib_timing-parse-timing "${_v1}/timing-compare-${j}.txt")
      j="ptgaul"
      read -r _memory_gb_ptgaul _total_hours_ptgaul < <(_polap_lib_timing-parse-timing "${_v1}/timing-ptgaul.txt")

      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "_${_species}_" \
        "${_P}" \
        "${_l_sra}" \
        "${_l_sra_size_gb}" \
        "${_s_sra}" \
        "${_s_sra_size_gb}" \
        "${_known_mtdna}" \
        "${_mode}" \
        "${_sd}" \
        "${_memory_gb_ptgaul}" \
        "${_memory_gb}" \
        "${_total_hours}" \
        >>"${_table1_file}"
    done

  done

  csvtk -t csv2md -a right ${_table1_file} -o table1.md
  echo "Check table1.md"
  # pandoc ${_table1_file} -f tsv -t markdown_mmd -o table1.md

  # pandoc table1.tsv -f tsv -t docx -o table1.docx
}

table2_genus_species() {
  local _table_file="table2.tsv"
  local _v1
  local _logfile
  local _readmefile
  local _species
  local _l_sra
  local _s_sra
  local _l_sra_size
  local _s_sra_size
  local _l_sra_size_gb
  local _s_sra_size_gb
  local _memory_gb
  local _total_hours
  local _memory_gb_ptgaul
  local _total_hours_ptgaul
  local _memory_gb_flye1
  local _total_hours_flye1
  local _I
  local _P
  local _N
  local j

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "Species" \
    "P" \
    "L_size" \
    "S_size" \
    "NCBI ptDNA" \
    "Length1" \
    "Length2" \
    "Mode" \
    "SD" \
    "M_r" \
    "M_a" \
    "M_f" \
    "T" \
    >"${_table_file}"

  for _v1 in "${Stable2[@]}"; do
    _logfile=${_v1}/polap.log
    _readmefile=${_v1}/README
    _species="${_v1/_/ }"
    _l_sra=$(awk '/long-read/ {match($0, /([A-Z]RR[0-9]+)/, arr); print arr[0]}' "$_logfile" | sort -u | grep -v '^$')
    _l_sra_size=$(<${_v1}/long_total_length.txt)
    _l_sra_size_gb=$(_polap_utility_convert_bp "${_l_sra_size}")
    _s_sra=$(awk '/short-read1/ {match($0, /([A-Z]RR[0-9]+)/, arr); print arr[0]}' "$_logfile" | sort -u | grep -v '^$')
    _s_sra_size=$(<${_v1}/short_total_length.txt)
    _s_sra_size_gb=$(_polap_utility_convert_bp "${_s_sra_size}")
    _known_mtdna=$(grep 'NCBI accession:' ${_logfile} | cut -d: -f4 | tail -n 1)

    if [[ -v _host["${_v1}"] ]]; then
      echo host:${_host["${_v1}"]}
    else
      echo "No such host for $_v1"
    fi

    local _ptdna_ptgaul=${_v1}/ptdna-ptgaul.fa
    local _ptdna_reference=${_v1}/ptdna-reference.fa
    seq_length1=$(grep -v "^>" "$_ptdna_ptgaul" | tr -d '\n' | wc -c)

    read -r _memory_gb_flye1 _total_hours_flye1 < <(_polap_lib_timing-parse-timing "${_v1}/timing-flye1.txt")

    for j in {1..1}; do
      local _disassemble_index="infer-${j}"
      local _summary1_ordered_txt=${_v1}/disassemble/${_disassemble_index}/1/summary1-ordered.txt
      local fasta_file=${_v1}/disassemble/${_disassemble_index}/pt.subsample-polishing.1.fa

      if [[ -s "${fasta_file}" ]]; then
        # Count the number of sequences (lines starting with '>')
        seq_count=$(grep -c "^>" "$fasta_file")

        # Ensure there is exactly one sequence
        if [[ "$seq_count" -ne 1 ]]; then
          echo "Error: FASTA file does not contain exactly one sequence."
          exit 1
        fi

        # Compute the sequence length (excluding header lines)
        seq_length2=$(grep -v "^>" "$fasta_file" | tr -d '\n' | wc -c)
      else
        seq_length2=0
      fi

      # Extract mode value
      # Extract SD value
      # Extract the first index value
      local _mode=$(grep "^#mode:" "$_summary1_ordered_txt" | awk '{print $2}')
      local _sd=$(grep "^#sd:" "$_summary1_ordered_txt" | awk '{print $2}')
      local _first_index=$(grep "^#index:" "$_summary1_ordered_txt" | awk 'NR==1 {print $2}')

      local _params_txt=${_v1}/disassemble/${_disassemble_index}/params.txt
      _ipn=$(parse_params "${_params_txt}")
      read -r _I _P _N <<<"$_ipn"

      read -r _memory_gb _total_hours < <(_polap_lib_timing-parse-timing "${_v1}/timing-infer-${j}-subsample-polish.txt")
      j="ptgaul"
      read -r _memory_gb_ptgaul _total_hours_ptgaul < <(_polap_lib_timing-parse-timing "${_v1}/timing-ptgaul.txt")

      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "_${_species}_" \
        "${_P}" \
        "${_l_sra_size_gb}" \
        "${_s_sra_size_gb}" \
        "${_known_mtdna}" \
        "${seq_length1}" \
        "${seq_length2}" \
        "${_mode}" \
        "${_sd}" \
        "${_memory_gb_ptgaul}" \
        "${_memory_gb_flye1}" \
        "${_memory_gb}" \
        "${_total_hours}" \
        >>"${_table_file}"
    done

  done

  csvtk -t csv2md -a right ${_table_file} -o table2.md
  echo "Check table2.md"
  # pandoc ${_table_file} -f tsv -t markdown_mmd -o table1.md

  # pandoc table1.tsv -f tsv -t docx -o table1.docx
}

table3_genus_species() {
  local _table_file="table3.tsv"
  local _v1
  local _logfile
  local _readmefile
  local _species
  local _l_sra
  local _s_sra
  local _l_sra_size
  local _s_sra_size
  local _l_sra_size_gb
  local _s_sra_size_gb
  local _memory_gb
  local _total_hours
  local _memory_gb_ptgaul
  local _total_hours_ptgaul
  local _memory_gb_flye1
  local _total_hours_flye1
  local _I
  local _P
  local _N
  local j

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "Species" \
    "P" \
    "L_size" \
    "S_size" \
    "NCBI ptDNA" \
    "Length1" \
    "Length2" \
    "Pident" \
    "Mode" \
    "SD" \
    "M_r" \
    "M_p" \
    "M_f" \
    "T" \
    >"${_table_file}"

  for _v1 in "${Stable3[@]}"; do
    _logfile=${_v1}/polap.log
    _readmefile=${_v1}/README
    _species="${_v1/_/ }"
    _l_sra=$(awk '/long-read/ {match($0, /([A-Z]RR[0-9]+)/, arr); print arr[0]}' "$_logfile" | sort -u | grep -v '^$')
    _l_sra_size=$(<${_v1}/long_total_length.txt)
    _l_sra_size_gb=$(_polap_utility_convert_bp "${_l_sra_size}")
    _s_sra=$(awk '/short-read1/ {match($0, /([A-Z]RR[0-9]+)/, arr); print arr[0]}' "$_logfile" | sort -u | grep -v '^$')
    _s_sra_size=$(<${_v1}/short_total_length.txt)
    _s_sra_size_gb=$(_polap_utility_convert_bp "${_s_sra_size}")
    _known_mtdna=$(grep 'NCBI accession:' ${_logfile} | cut -d: -f4 | tail -n 1)

    if [[ -v _host["${_v1}"] ]]; then
      echo host:${_host["${_v1}"]}
    else
      echo "No such host for $_v1"
    fi

    local _ptdna_ptgaul=${_v1}/ptdna-ptgaul.fa
    local _ptdna_reference=${_v1}/ptdna-reference.fa
    seq_length1=$(grep -v "^>" "$_ptdna_ptgaul" | tr -d '\n' | wc -c)

    # read -r _memory_gb_flye1 _total_hours_flye1 < <(parse_timing "${_v1}" "infer12-4")

    for j in {1..1}; do
      local _disassemble_index="infer-${j}"
      local _summary1_ordered_txt=${_v1}/disassemble/${_disassemble_index}/1/summary1-ordered.txt
      local fasta_file=${_v1}/disassemble/${_disassemble_index}/pt.subsample-polishing.1.fa

      if [[ -s "${fasta_file}" ]]; then
        # Count the number of sequences (lines starting with '>')
        seq_count=$(grep -c "^>" "$fasta_file")

        # Ensure there is exactly one sequence
        if [[ "$seq_count" -ne 1 ]]; then
          echo "Error: FASTA file does not contain exactly one sequence."
          exit 1
        fi

        # Compute the sequence length (excluding header lines)
        seq_length2=$(grep -v "^>" "$fasta_file" | tr -d '\n' | wc -c)
      else
        seq_length2=0
      fi

      # Extract mode value
      # Extract SD value
      # Extract the first index value
      local _mode=$(grep "^#mode:" "$_summary1_ordered_txt" | awk '{print $2}')
      local _sd=$(grep "^#sd:" "$_summary1_ordered_txt" | awk '{print $2}')
      local _first_index=$(grep "^#index:" "$_summary1_ordered_txt" | awk 'NR==1 {print $2}')

      local _params_txt=${_v1}/disassemble/${_disassemble_index}/params.txt
      _ipn=$(parse_params "${_params_txt}")
      read -r _I _P _N <<<"$_ipn"

      read -r _memory_gb _total_hours < <(_polap_lib_timing-parse-timing "${_v1}/timing-infer-${j}-subsample-polish.txt")
      read -r _memory_gb_polishing _total_hours_polishing < <(_polap_lib_timing-parse-timing "${_v1}/timing-check-${j}-subsample-polish.txt")
      read -r _memory_gb_ptgaul _total_hours_ptgaul < <(_polap_lib_timing-parse-timing "${_v1}/timing-ptgaul.txt")

      local _pident="NA"
      if [[ -s "${_v1}/blast/pident.txt" ]]; then
        _pident="$(<${_v1}/blast/pident.txt)"
      fi
      # read -r _memory_gb _total_hours < <(parse_timing "${_v1}" "infer12-${j}")
      # read -r _memory_gb_polishing _total_hours_polishing < <(parse_timing "${_v1}" "infer3only-${j}")

      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "_${_species}_" \
        "${_P}" \
        "${_l_sra_size_gb}" \
        "${_s_sra_size_gb}" \
        "${_known_mtdna}" \
        "${seq_length1}" \
        "${seq_length2}" \
        "${_pident}" \
        "${_mode}" \
        "${_sd}" \
        "${_memory_gb_ptgaul}" \
        "${_memory_gb_polishing}" \
        "${_memory_gb}" \
        "${_total_hours}" \
        >>"${_table_file}"
    done

  done

  csvtk -t csv2md -a right ${_table_file} -o table3.md
  echo "Check table3.md"
  cat table3.md
}

suptable1_genus_species() {
  local _suptable1_file="suptable1.md"
  local _species

  local _v1
  local _logfile
  local _readmefile
  local _l_sra
  local _s_sra
  local _memory_gb
  local _total_hours
  local _memory_gb_ptgaul
  local _total_hours_ptgaul
  local j
  local _percent_identity
  local _I
  local _P
  local _N
  local _label_base
  local _label

  printf "# Supplementary Tables\n\n" \
    >"${_suptable1_file}"

  for _v1 in "${S[@]}"; do
    _logfile=${_v1}/polap.log
    _readmefile=${_v1}/README
    _species="${_v1/_/ }"
    _label_base="${_v1/_/-}"
    _label_base=$(echo "$_label_base" | awk '{print tolower($0)}')

    for j in {1..3}; do
      local _params_txt=${_v1}/disassemble/${j}/params.txt
      _ipn=$(parse_params "${_params_txt}")
      read -r _I _P _N <<<"$_ipn"
      _label=${_label_base}-${_I}

      # stage 1
      printf "\\\blandscape\n\n" \
        >>"${_suptable1_file}"

      printf "Table: Variability in assembly replicates for the maximum subsampling rate of ${_P}%% from the stage 1 of Polap's analysis in _${_species}_. {#tbl:table1-${_label}}\n\n" \
        >>"${_suptable1_file}"

      cat "${_v1}/disassemble/${j}/1/summary1.md" \
        >>"${_suptable1_file}"
      # Table: _Juncus effusus_ Polap's disassemble data analysis. i=3, p=10, n=50 \label{table-juncus-effusus} {#tbl:table-juncus-effusus}

      printf "\n" \
        >>"${_suptable1_file}"

      cat src/polap-data-v2-table1_footnote.tex \
        >>"${_suptable1_file}"

      printf "\\\elandscape\n\n\\\newpage\n\n" \
        >>"${_suptable1_file}"

    done
  done

  for _v1 in "${S[@]}"; do
    _logfile=${_v1}/polap.log
    _readmefile=${_v1}/README
    _species="${_v1/_/ }"
    _label_base="${_v1/_/-}"
    _label_base=$(echo "$_label_base" | awk '{print tolower($0)}')

    for j in {1..3}; do
      local _params_txt=${_v1}/disassemble/${j}/params.txt
      _ipn=$(parse_params "${_params_txt}")
      read -r _I _P _N <<<"$_ipn"
      _label=${_label_base}-${_I}

      # stage 2
      printf "\\\blandscape\n\n" \
        >>"${_suptable1_file}"

      printf "Table: Variability in assembly replicates for the maximum subsampling rate of ${_P}%% from the stage 2 of Polap's analysis in _${_species}_. {#tbl:table2-${_label}}\n\n" \
        >>"${_suptable1_file}"

      cat "${_v1}/disassemble/${j}/2/summary1.md" \
        >>"${_suptable1_file}"
      # Table: _Juncus effusus_ Polap's disassemble data analysis. i=3, p=10, n=50 \label{table-juncus-effusus} {#tbl:table-juncus-effusus}

      printf "\n" \
        >>"${_suptable1_file}"

      cat src/polap-data-v2-table2_footnote.tex \
        >>"${_suptable1_file}"

      printf "\\\elandscape\n\n\\\newpage\n\n" \
        >>"${_suptable1_file}"
    done
  done

  echo "See ${_suptable1_file}"
}

supfigure1_genus_species() {
  local _supfigure1_file="supfigure1.md"
  local _species

  local _v1
  local _logfile
  local _readmefile
  local _l_sra
  local _s_sra
  local _memory_gb
  local _total_hours
  local _memory_gb_ptgaul
  local _total_hours_ptgaul
  local j
  local _percent_identity
  local _I
  local _P
  local _N
  local _label
  local _label_base

  printf "# Supplementary Figures\n\n" \
    >"${_supfigure1_file}"

  for _v1 in "${S[@]}"; do
    _logfile=${_v1}/polap.log
    _readmefile=${_v1}/README
    _species="${_v1/_/ }"
    _label_base="${_v1/_/-}"
    _label_base=$(echo "$_label_base" | awk '{print tolower($0)}')

    for j in {1..3}; do
      local _params_txt=${_v1}/disassemble/${j}/params.txt
      _ipn=$(parse_params "${_params_txt}")
      read -r _I _P _N <<<"$_ipn"
      _label=${_label_base}-${_I}

      # stage 1
      printf "\\\newpage\n\n" \
        >>"${_supfigure1_file}"

      printf "![Selection of a plastid genome using the length distribution for _${_species}_. (A) Trace plot of plastid genome length and the subsample-size upto ${_P} %%. (B) Length distribution for candidate plastid genome of _${_species}_ ](figures/${_label}-summary1-ordered.pdf){#fig:figure-${_label}}\n\n" \
        >>"${_supfigure1_file}"

      cp "${_v1}/disassemble/${j}/1/summary1-ordered.pdf" \
        figures/${_label}-summary1-ordered.pdf

      printf "\n" \
        >>"${_supfigure1_file}"

    done
  done
  echo "See ${_supfigure1_file}"
  echo "Copy PDF files in figures to the manuscript folder."
  echo "  cp figures/* ~/all/manuscript/polap-v0.4/figures/"
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
    echo "0. Test"
    echo "1. Send data to"
    echo "2. Make a starting folder (3,4)"
    echo "3. Get ptDNA from NCBI (4)"
    echo "4. Copy the ptDNA as a reference"
    echo "5. ptGAUL assembly with the reference (6)"
    echo "6. Prepare the short-read polishing"
    echo "7. Extract the ptGAUL assembly (8)"
    echo "8. Copy the ptGAUL assembly"
    echo "9. Infer with POLAP"
    echo "10. Check with POLAP (needs step 5,7)"
    echo "11. (only for 5 datasets) Compare with POLAP (needs step 5,7)"
    echo "A. Archive the ${_species_folder}"
    echo "G. Get the result from the remote"
    echo "C. Clean up the ${_species_folder}"
    echo "M. Back to Main Menu"
    echo "--------------------------------------------------"
    echo your current host: "$(hostname)"
    echo host name: "${_host["$_species_name"]}"
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
          echo "Valid number"
          echo "Invalid choice, please select a valid option."
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
            test_genus_species "${keys_array[$species_key]}"
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
            copy-ptdna-of-ncbi-as-reference_genus_species \
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
          g)
            get_genus_species \
              "${keys_array[$species_key]}"
            ;;
          a)
            archive_genus_species \
              "${keys_array[$species_key]}"
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
  navigate_menu "main"
}

# Main case statement
case "$subcmd1" in
'seeds' | \
  'mkdir' | \
  'map' | \
  'map' | \
  'reads' | \
  'msbwt' | \  | \
  'dflye' | \
  'directional' | \
  'archive' | \
  'report' | \
  'table1' | \
  'figure1' | \
  'copy-ptdna-of-ptgaul' | \
  'compare' | \
  'get' | \
  'table2' | \
  'table3' | \
  'suptable1' | \
  'infer12' | \
  'infer2only' | \
  'infer3only' | \
  'polishing' | \
  'wga' | \
  'polish' | \
  'mauve' | \
  'clean' | \
  'write-config' | \
  'test')
  ${subcmd1}_genus_species ${subarg1}
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
"delete-links")
  find . -type l -delete
  ;;
"zip")
  for i in "${S[@]}"; do
    echo "zipping $i ..."
    tar zcf ${i}.tar.gz ${i}
  done
  ;;
"mkdir-all")
  for i in "${S[@]}"; do
    echo "creating folder $i ..."
    mkdir ${i}
  done
  ;;
"mkdir-media")
  for i in "${S[@]}"; do
    echo "creating folder $i ..."
    mkdir ${_media_dir}/${i}
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
  bash src/polap-v0.4-report-table1.sh >table1.tsv
  pandoc table1.tsv -f tsv -t pdf -o table1.pdf -V geometry:landscape
  pandoc table1.tsv -f tsv -t docx -o table1.docx
  ;;
"sync")
  for i in "${S[@]}"; do
    V1="${i}"
    cd $V1
    if [[ -v _ssh["${V1}"] ]]; then
      mkdir -p o/${_arg2}/01-contig o/${_arg2}/06-summary o/${_arg2}/07-plot
      rsync -aPh ${_ssh["${V1}"]}:$PWD/o/${_arg2}/01-contig/ o/${_arg2}/01-contig/
      rsync -aPh ${_ssh["${V1}"]}:$PWD/o/${_arg2}/06-summary/ o/${_arg2}/06-summary/
      rsync -aPh ${_ssh["${V1}"]}:$PWD/o/${_arg2}/07-plot/ o/${_arg2}/07-plot/
      rsync -aPh ${_ssh["${V1}"]}:$PWD/o/${_arg2}/30-contigger o/${_arg2}/
      rsync -aPh ${_ssh["${V1}"]}:$PWD/o/${_arg2}/*annotation* o/${_arg2}/
      rsync -aPh ${_ssh["${V1}"]}:$PWD/o/${_arg2}/assembly.fasta o/${_arg2}/
      rsync -aPh ${_ssh["${V1}"]}:$PWD/o/${_arg2}/assembly_graph.gfa o/${_arg2}/
      rsync -aPh ${_ssh["${V1}"]}:$PWD/o/00-bioproject o/
      # rsync -aPh o/${_arg2}/mt.0.fasta ${_ssh["${V1}"]}:$PWD/o/${_arg2}/
      rsync -aPh ${_ssh["${V1}"]}:$PWD/o/polap.log o/
      # cat o/${_arg2}/07-plot/*.md >>../plot.md
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
*)
  echo "Species '$subcmd1' is not recognized."
  exit 1
  ;;
esac
