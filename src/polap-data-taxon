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

DEBUG=0
_POLAP_RELEASE=1

# Taxonomy
#
# the 12 datasets
S=(
  'Prunus_mandshurica'
  'Eucalyptus_pauciflora'
  'Eutrema_salsugineum'
  'Linum_usitatissimum'
  'Dioscorea_japonica'
  'Gossypium_herbaceum'
  'Canavalia_ensiformis'
  'Hylodesmum_podocarpum'
  'Crataegus_pinnatifida'
  'Pterocarpus_santalinus'
  'Megaceros_flagellaris'
  'Paraphymatoceros_pearsonii'
)
# the five datasets
S=(
  'Eucalyptus_pauciflora'
  'Eutrema_salsugineum'
  'Hylodesmum_podocarpum'
  'Crataegus_pinnatifida'
  'Megaceros_flagellaris'
)
# the two datasets
S=(
  'Megaceros_flagellaris'
  'Hylodesmum_podocarpum'
  'Anthoceros_agrestis'
)

_local_host="thorne"

# Input parameter
subcmd1="${1:-0}"
# subarg1="${2:-0}"

_polap_subcmd=(
  'test'
  'send-data-to'
  'mkdir'
  'taxon-ref1'
  'taxon-sample1'
  'taxon-sample2'
  'taxon-geseq'
  'taxon-orthofinder'
  'taxon-phylogeny'
  'taxon-tree'
  'taxon-species'
  'archive'
  'get'
  'report'
  'table1'
  'table2'
  'suptable1'
  'supfigure1'
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

# echo "subcmd1 is now: $subcmd1"

if [[ -d "src" ]]; then
  _polap_cmd="src/polap.sh"
else
  _polap_cmd="polap"
fi
_polap_version="0.4.3.7"
_media_dir="/media/h1/run/ptgaul20"
_media_dir="/media/h2/goshng/bioprojects"

help_message=$(
  cat <<HEREDOC
# Polap data analysis for the version 0.4.7 taxonomy
#
# Edit two variables:
# _polap_cmd=${_polap_cmd}
# _media_dir=${_media_dir}
# _local_host=${_local_host}
# host=$(hostname)
#
# subcommand can be replaced with the leading number
0. test <species_folder>
1. send-data-to <species_folder>: scp data files to the remote 
2. mkdir <species_folder>
3. taxon-assemble
4. taxon-sample1
5. taxon-sample2
6. taxon-geseq
7. taxon-orthofinder
8. taxon-phylogeny
9. taxon-tree
10. archive <species_folder>: report the results
11. get <species_folder>: fetch the archive from the remote
12. report [species_folder]: report the results
13. table1: table1
14. table2: table2
15. suptable1: Supplementary tables
16. supfigure1: Supplementary figures
HEREDOC
)

declare -A _host
_host['Megaceros_flagellaris']="kishino"
_host['Hylodesmum_podocarpum']="vincent"

declare -A _long
_long['Megaceros_flagellaris']="SRR25397414"
_long['Hylodesmum_podocarpum']="SRR18714550"
declare -A _short
_short['Megaceros_flagellaris']="SRR25397413"
_short['Hylodesmum_podocarpum']="SRR18714546"
declare -A _inref
_inref['Megaceros_flagellaris']="family"
_inref['Hylodesmum_podocarpum']="family"
declare -A _ingroup
_ingroup['Megaceros_flagellaris']="class"
_ingroup['Hylodesmum_podocarpum']="family"
declare -A _ingroup_sample
_ingroup_sample['Megaceros_flagellaris']="family"
_ingroup_sample['Hylodesmum_podocarpum']="genus"
declare -A _allgroup
_allgroup['Megaceros_flagellaris']="phylum"
_allgroup['Hylodesmum_podocarpum']="class"
declare -A _outgroup_sample
_outgroup_sample['Megaceros_flagellaris']="class"
_outgroup_sample['Hylodesmum_podocarpum']="order"

_taxonomy_min_aa=20

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

if [[ -z "$2" ]]; then
  _arg2=$(default_output_dir)
else
  _arg2="${2%/}"
fi

# Check if the species folder is provided
if [[ -z "$1" ]]; then
  echo "Usage: $0 <subcommand> [species_folder]"
  echo "${help_message}"
  exit 1
fi

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

taxon-ref1_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local inref="${_inref["$1"]}"
  # copy_data

  # --taxonomy-rank-ingroup "${inref}" \
  ${_polap_cmd} taxonomy reference \
    --steps-include 1-6 \
    -v \
    -o ${output_dir} \
    --taxonomy-rank-ingroup "genus" \
    --plastid --redo \
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
    --steps-include 1-6 \
    -v \
    -o ${output_dir} \
    --taxonomy-rank-ingroup ${ingroup} \
    --taxonomy-rank-allgroup ${allgroup} \
    --species "${species_name}"
}

taxon-sample2_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"
  local ingroup="${_ingroup["$1"]}"
  local allgroup="${_allgroup["$1"]}"
  local ingroup_sample="${_ingroup_sample["$1"]}"
  local outgroup_sample="${_outgroup_sample["$1"]}"
  # copy_data

  ${_polap_cmd} taxonomy sample \
    --steps-include 7-12 \
    -o ${output_dir} \
    --taxonomy-ingroup-size 10 \
    --taxonomy-outgroup-size 10 \
    --taxonomy-sample-size-per-rank 1 \
    --taxonomy-rank-ingroup ${ingroup} \
    --taxonomy-rank-ingroup-sample ${ingroup_sample} \
    --taxonomy-rank-allgroup ${allgroup} \
    --taxonomy-rank-outgroup-sample ${outgroup_sample}
}

taxon-geseq_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"
  # copy_data

  ${_polap_cmd} taxonomy geseq \
    --steps-include 1-6 \
    -o ${output_dir} \
    --species "${species_name}"
}

taxon-orthofinder_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"
  # copy_data

  ${_polap_cmd} taxonomy orthofinder \
    --steps-include 1-5 --redo \
    -v \
    -o ${output_dir} \
    --taxonomy-min-aa "${_taxonomy_min_aa}" \
    --species "${species_name}"
}

taxon-phylogeny_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"
  # copy_data

  ${_polap_cmd} taxonomy phylogeny \
    --steps-include 1-2 --redo \
    -o ${output_dir} \
    -l ${long_sra}.fastq \
    -a ${short_sra}_1.fastq \
    -b ${short_sra}_2.fastq \
    --species "${species_name}"
}

taxon-tree_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"
  # copy_data

  ${_polap_cmd} taxonomy tree \
    --steps-include 1-2 \
    -v \
    -o ${output_dir} \
    --species "${species_name}"
}

taxon-species_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  # copy_data

  ${_polap_cmd} taxonomy species \
    --steps-include 1 \
    -o ${output_dir}
}

geseq_genus_species() {
  local output_dir="$1"
  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$1"]}"
  local short_sra="${_short["$1"]}"
  # copy_data

  ${_polap_cmd} taxonomy \
    --stages-include 1 \
    -o ${output_dir} \
    -l ${long_sra}.fastq \
    -a ${short_sra}_1.fastq \
    -b ${short_sra}_2.fastq \
    --taxonomy-ingroup class \
    --taxonomy-outgroup phylum \
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
    local _memory_gb=$(echo "scale=2; $_memory_kbytes / 1048576" | bc)
    # Extract the time portion using grep with regex
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
    "Length" \
    "SD" \
    "M_r" \
    "M_f" \
    "T" \
    >"${_table1_file}"

  for _v1 in "${S[@]}"; do
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

    for j in {1..3}; do
      local _summary1_ordered_txt=${_v1}/disassemble/${j}/1/summary1-ordered.txt

      # Extract mode value
      # Extract SD value
      # Extract the first index value
      local _mode=$(grep "^#mode:" "$_summary1_ordered_txt" | awk '{print $2}')
      local _sd=$(grep "^#sd:" "$_summary1_ordered_txt" | awk '{print $2}')
      local _first_index=$(grep "^#index:" "$_summary1_ordered_txt" | awk 'NR==1 {print $2}')

      local _params_txt=${_v1}/disassemble/${j}/params.txt
      _ipn=$(parse_params "${_params_txt}")
      read -r _I _P _N <<<"$_ipn"

      read -r _memory_gb _total_hours < <(parse_timing "${_v1}" "${j}")
      j="ptgaul"
      read -r _memory_gb_ptgaul _total_hours_ptgaul < <(parse_timing "${_v1}" "${j}")

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
  local _I
  local _P
  local _N
  local j

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
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
    "M_f" \
    "T" \
    >"${_table_file}"

  for _v1 in "${S[@]}"; do
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

    for j in {4..6}; do
      local _summary1_ordered_txt=${_v1}/disassemble/${j}/1/summary1-ordered.txt
      local fasta_file=${_v1}/disassemble/${j}/pt.1.fa

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

      local _params_txt=${_v1}/disassemble/${j}/params.txt
      _ipn=$(parse_params "${_params_txt}")
      read -r _I _P _N <<<"$_ipn"

      read -r _memory_gb _total_hours < <(parse_timing "${_v1}" "${j}")
      j="ptgaul"
      read -r _memory_gb_ptgaul _total_hours_ptgaul < <(parse_timing "${_v1}" "${j}")

      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
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

scopy_data_to_remote() {
  mkdir -p "${output_dir}"
  if [[ -s "${long_sra}.fastq" ]]; then
    scp "${_media_dir}/${long_sra}.fastq" ${_host["${output_dir}"]}:$PWD/
  fi
  if [[ -s "${short_sra}_1.fastq" ]]; then
    scp "${_media_dir}/${short_sra}_1.fastq" ${_host["${output_dir}"]}:$PWD/
  fi
  if [[ -s "${short_sra}_2.fastq" ]]; then
    scp "${_media_dir}/${short_sra}_2.fastq" ${_host["${output_dir}"]}:$PWD/
  fi
}

copy_data() {
  mkdir -p "${output_dir}"
  if [[ ! -s "${long_sra}.fastq" ]]; then
    cp "${_media_dir}/${long_sra}.fastq" .
  fi
  if [[ ! -s "${short_sra}_1.fastq" ]]; then
    cp "${_media_dir}/${short_sra}_1.fastq" .
  fi
  if [[ ! -s "${short_sra}_2.fastq" ]]; then
    cp "${_media_dir}/${short_sra}_2.fastq" .
  fi
  if [[ ! -s "ptdna-${output_dir}.fa" ]]; then
    ${_polap_cmd} get-mtdna --plastid --species "${species_name}" -o ${output_dir}
    cp ${output_dir}/00-bioproject/2-mtdna.fasta ptdna-${output_dir}.fa
  fi
}

# Main case statement
case "$subcmd1" in
'send-data-to' | \
  'mkdir' | \
  'taxon-ref1' | \
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
  'ptgaul' | \
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
  ${subcmd1}_genus_species ${_arg2}
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
    if [[ -v _host["${V1}"] ]]; then
      mkdir -p o/${_arg2}/01-contig o/${_arg2}/06-summary o/${_arg2}/07-plot
      rsync -aPh ${_host["${V1}"]}:$PWD/o/${_arg2}/01-contig/ o/${_arg2}/01-contig/
      rsync -aPh ${_host["${V1}"]}:$PWD/o/${_arg2}/06-summary/ o/${_arg2}/06-summary/
      rsync -aPh ${_host["${V1}"]}:$PWD/o/${_arg2}/07-plot/ o/${_arg2}/07-plot/
      rsync -aPh ${_host["${V1}"]}:$PWD/o/${_arg2}/30-contigger o/${_arg2}/
      rsync -aPh ${_host["${V1}"]}:$PWD/o/${_arg2}/*annotation* o/${_arg2}/
      rsync -aPh ${_host["${V1}"]}:$PWD/o/${_arg2}/assembly.fasta o/${_arg2}/
      rsync -aPh ${_host["${V1}"]}:$PWD/o/${_arg2}/assembly_graph.gfa o/${_arg2}/
      rsync -aPh ${_host["${V1}"]}:$PWD/o/00-bioproject o/
      # rsync -aPh o/${_arg2}/mt.0.fasta ${_host["${V1}"]}:$PWD/o/${_arg2}/
      rsync -aPh ${_host["${V1}"]}:$PWD/o/polap.log o/
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
  echo "Species '$subcmd1' is not recognized as the subcommand of $0."
  exit 1
  ;;
esac
