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

system_genus_species() {
  echo "Host: $(hostname)"
  # echo "  CPU: lscpu"
  lscpu | grep Model | head -1 | tr -s ' '
  # echo "  Memory: free -h"
  free -h | grep Mem | awk '{print "Memory:", $2}'
}

function _polap_lib_data-example-data {
  local _brg_data="${1}"

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

  echo "${_data}" >"${_brg_data}"
  echo "create ${_brg_data}"
}

help_message_bleeding_edge_polap=$(
  cat <<HEREDOC

  Patch polap with the latest github version.
HEREDOC
)

help_message_patch_polap=$(
  cat <<HEREDOC

  Patch polap with a specific version.
  _polap_version=<VERSION> polap patch-polap
HEREDOC
)

help_message_install_getorganelle=$(
  cat <<HEREDOC

  Install conda environments: getorganelle
  conda create -y --name getorganelle bioconda::getorganelle
  conda activate getorganelle
  get_organelle_config.py --add embplant_pt,embplant_mt
HEREDOC
)

help_message_install_cflye=$(
  cat <<HEREDOC

  Install cflye to polap conda environment
  conda install -y goshng::cflye
HEREDOC
)

help_message_install_fmlrc=$(
  cat <<HEREDOC

  Install conda environments: polap-fmlrc
HEREDOC
)

help_message_install_polap=$(
  cat <<HEREDOC

  Install conda environments: polap
  conda create -y --name polap bioconda::polap
HEREDOC
)

help_message_install_tippo=$(
  cat <<HEREDOC

  Install conda environments: tippo
  conda create -y --name tippo bioconda::tippo
HEREDOC
)

help_message_install_oatk=$(
  cat <<HEREDOC

  Install conda environments: oatk

  conda create -y --name oatk bioconda::oatk
  conda activate oatk
  conda install -y biopython
  conda install -y hmmer seqtk mafft parallel entrez-direct
  git clone https://github.com/c-zhou/OatkDB.git
HEREDOC
)

help_message_uninstall=$(
  cat <<HEREDOC

  Uninstall conda environments: polap, polap-fmlrc, getorganelle
  conda remove -n polap --all
  conda remove -n polap-fmlrc --all
  conda remove -n getorganelle --all
  conda env remove -n pmat
HEREDOC
)

##### INSERT_HELP_HERE #####
help_message_merge=$(
  cat <<HEREDOC

  menu title
HEREDOC
)

help_message_hello_world=$(
  cat <<HEREDOC

  menu title
HEREDOC
)

help_message_rsync=$(
  cat <<HEREDOC

  menu title
HEREDOC
)

help_message_oatk_polished=$(
  cat <<HEREDOC

  all: process all 
  outdir: process only outdir
  inum: produce outdir/<inum>
  polished: use hifiasm or nextdenovo-polished data
  fc: oatk -c option value
HEREDOC
)

help_message_hifiasm_polish=$(
  cat <<HEREDOC

  Polish an ONT long-read fastq data file.
  Ref: https://github.com/c-zhou/oatk/issues/16

  Code:
  hifiasm -t64 --ont -o ONT.asm ONT.read.fastq.gz --write-ec

  mkdir temp
  cd temp
  ln -s ../ONT.asm.ec.bin
  ln -s ../ONT.asm.ovlp.source.bin
  ln -s ../ONT.asm.ovlp.reverse.bin
  hifiasm -t64 --ont -o ONT.asm --write-ec /dev/null
HEREDOC
)

help_message_mitohifi=$(
  cat <<HEREDOC

  mitohifi using apptainer
HEREDOC
)

help_message_install_apptainer=$(
  cat <<HEREDOC

  Ref: https://github.com/apptainer/apptainer
HEREDOC
)

help_message_install_novoplasty=$(
  cat <<HEREDOC

  Ref: https://anaconda.org/bioconda/novoplasty
  Ref: https://github.com/ndierckx/NOVOPlasty.git
  NOVOPlasty4.3.1.pl -c config.txt
HEREDOC
)

help_message_install_mitohifi=$(
  cat <<HEREDOC

  Ref: https://github.com/marcelauliano/MitoHiFi.git
  docker pull ghcr.io/marcelauliano/mitohifi:master
  singularity exec --bind /path/to/container_directory:/path/to/container_directory docker://ghcr.io/marcelauliano/mitohifi:master mitohifi.py -h
HEREDOC
)

help_message_data_downsample_long=$(
  cat <<HEREDOC

  Subsample the long-read data and save it under o/tmp.
  We want to see the effect on long-read subsampling on the final result.
HEREDOC
)

help_message_data_downsample_short=$(
  cat <<HEREDOC

  Subsample the short-read data and save it at o/tmp.
  We want to reduce the memory usage by polishing using FMLRC.
HEREDOC
)

help_message_estimate_genomesize=$(
  cat <<HEREDOC

  Estimate the genome size.
HEREDOC
)

help_message_data_long=$(
  cat <<HEREDOC

  Download the long-read data file.
HEREDOC
)

help_message_data_short=$(
  cat <<HEREDOC

  Download the short-read data file.
HEREDOC
)

##### INSERT_FUNCTION_HERE #####
merge_genus_species_for() {
  local _brg_outdir="${1}"
  local _brg_source="${2}"
  local _brg_dest="${3}"
  local _brg_as="${4}"

  # folders in data-v1
  # o1: whole-genome assembly
  # o2: organelle-genome assembly
  # b1: benchmarking with PMAT, tippo, oatk with long-read polishing
  # t1: whole-genome assembly with different depth of input data 10=10x, 20=1x
  # p1: ptgaul assembly for comparison (not presented in a paper)
  # s1: polap analysis with subsampled data (different from t1 where subsampled data used for WGA)
  #     all of the long-read is used in t1 while the subsampled data are used in the organelle-genome
  #     assembly. The results can be similar although we would see smaller data used in OGA with s1.
  #     It seems deirable to perform t1 analysis rather than s1 because we have no reason why we
  #     have to use a subset of the long-read data in OGA. In s1, we subsample data to create
  #     a smaller dataset as an input data for Polap. In t1, we feed the original long-read data
  #     to Polap, which uses -c option to create nk.fq.gz, a ruduced subset by the estimated genome
  #     size.

  # step 1
  # echo "merging ${_brg_outdir} from ${_brg_source} to ${_brg_dest} as ${_brg_as} ..."
  # mkdir -p "${_brg_dest}/${_brg_outdir}/${_brg_as}"
  # rsync -azu --info=progress2 "${_brg_source}/${_brg_outdir}/" \
  #   "${_brg_dest}/${_brg_outdir}/${_brg_as}/"

  # step 2
  for i in p1; do
    mkdir -p "${_brg_dest}/${_brg_outdir}/${i}"
  done
  for i in known-mtdna ptgaul-known ptgaul-reference taxonomy result_3000 ptgaul-known.fa ptgaul-reference.fa ncbi-species-genome-size.txt; do
    mv "${_brg_dest}/${_brg_outdir}/s1/${i}" "${_brg_dest}/${_brg_outdir}/p1/"
  done
  mv "${_brg_dest}/${_brg_outdir}/s1/timing"-ptgaul-*.txt "${_brg_dest}/${_brg_outdir}/p1/"
  mv "${_brg_dest}/${_brg_outdir}/s1/0" "${_brg_dest}/${_brg_outdir}/b1"
  mv "${_brg_dest}/${_brg_outdir}/s1/ptiming" "${_brg_dest}/${_brg_outdir}/t1"
  mv "${_brg_dest}/${_brg_outdir}/s1" "${_brg_dest}/${_brg_outdir}/backup1"
  mv "${_brg_dest}/${_brg_outdir}/backup1/subsample" "${_brg_dest}/${_brg_outdir}/s1"
  mv "${_brg_dest}/${_brg_outdir}/backup1"/s*x.txt "${_brg_dest}/${_brg_outdir}/s1"
  mv "${_brg_dest}/${_brg_outdir}/backup1"/l*x.txt "${_brg_dest}/${_brg_outdir}/s1/"
  mv "${_brg_dest}/${_brg_outdir}/backup1"/*.fq.txt "${_brg_dest}/${_brg_outdir}/s1/"
}

merge_genus_species() {
  local _brg_outdir="${1}"
  local _brg_source="${2:-/media/h2/goshng/aflye1}"
  local _brg_dest="${3:-/media/h3/labshare/goshng/aflye1}"
  local _brg_as="${4}"

  if [[ "${_brg_outdir}" == "all" ]]; then
    for _v1 in "${Sall[@]}"; do
      merge_genus_species_for "${_v1}" "${@:2}"
    done
  elif [[ "${_brg_outdir}" == "each" ]]; then
    for _v1 in "${Sall[@]}"; do
      merge_genus_species_for "${_v1}" "${@:2}"
    done
  else
    merge_genus_species_for "$@"
  fi
}

hello-world_genus_species_for() {
  local _brg_outdir="${1}"
  echo "Preparing archive for ${_brg_outdir} ..."
  # tar zcf "${_brg_outdir}-a.tar.gz" "${_brg_outdir}"
}

hello-world_genus_species() {
  local _brg_outdir="${1-all}"
  local _brg_inum="${2:-0}"
  local _brg_polished="${3:-hifiasm}"
  local _brg_fc="${4:-30}"

  if [[ "${_brg_outdir}" == "all" ]]; then
    for _v1 in "${Sall[@]}"; do
      hello-world_genus_species_for "${_v1}" "${@:2}"
    done
  elif [[ "${_brg_outdir}" == "each" ]]; then
    for _v1 in "${Sall[@]}"; do
      hello-world_genus_species_for "${_v1}" "${@:2}"
    done
  else
    hello-world_genus_species_for "$@"
  fi
}

rsync_genus_species_for() {
  local _brg_outdir="${1}"
  local _brg_inum="${2:-0}"
  local _brg_remote="${3:-siepel}"

  local _brg_outdir_i="${_brg_outdir}/${_brg_inum}"
  rsync -azvu "${_brg_remote}:${PWD}/${_brg_outdir_i}/" "${_brg_outdir_i}/"
}

rsync_genus_species() {
  local _brg_outdir="${1-all}"
  local _brg_inum="${2:-0}"
  local _brg_remote="${3:-siepel}"

  if [[ "${_brg_outdir}" == "all" ]]; then
    for _v1 in "${Sall[@]}"; do
      rsync_genus_species_for "${_v1}" "${@:2}"
    done
  elif [[ "${_brg_outdir}" == "each" ]]; then
    for _v1 in "${Sall[@]}"; do
      rsync_genus_species_for "${_v1}" "${@:2}"
    done
  else
    rsync_genus_species_for "$@"
  fi
}

nextdenovo-check-cns_genus_species() {
  if [[ "${_local_host}" != "$(hostname)" ]]; then
    if ssh ${_local_host} "test -f $PWD/${_brg_outdir_i}/cns.fa"; then
      if [[ ! -s "${_brg_outdir_i}/cns.fa" ]]; then
        scp -p ${_local_host}:$PWD/"${_brg_outdir_i}/cns.fa" "${_brg_outdir_i}"
      fi
    else
      echo "no such file at ${_local_host}: $PWD/${_brg_outdir_i}/cns.fa"
      echo "We need the polished long-read data for oatk."
      return 1
    fi
    if ssh ${_local_host} "test -f $PWD/${_brg_outdir}/short_expected_genome_size.txt"; then
      scp -p ${_local_host}:"$PWD/${_brg_outdir}/short_expected_genome_size.txt" "${_brg_outdir}"
    fi
    if ssh ${_local_host} "test -f $PWD/${_brg_outdir}/o/short_expected_genome_size.txt"; then
      mkdir -p "${_brg_outdir}/o"
      scp -p ${_local_host}:"$PWD/${_brg_outdir}/o/short_expected_genome_size.txt" "${_brg_outdir}/o"
    fi
  else
    if is_file_at_least_1MB "${_brg_outdir_i}/cns.fa"; then
      echo "found: ${_brg_outdir_i}/cns.fa"
    else
      echo "no such file (greater than 1MB) ${_brg_outdir_i}/cns.fa"
      echo "We need the polished long-read data for oatk."
      return 1
    fi
  fi
  return 0
}

hifiasm-check-cns_genus_species() {
  if [[ "${_local_host}" != "$(hostname)" ]]; then
    if ssh ${_local_host} "test -f $PWD/${_brg_outdir_i}/ont.asm.ec.fq"; then
      if [[ ! -s "${_brg_outdir_i}/ont.asm.ec.fq" ]]; then
        scp -p ${_local_host}:$PWD/"${_brg_outdir_i}/ont.asm.ec.fq" "${_brg_outdir_i}"
      fi
    else
      echo "no such file at ${_local_host}: $PWD/${_brg_outdir_i}/ont.asm.ec.fq"
      echo "We need the polished long-read data for oatk."
      return 1
    fi

    if ssh ${_local_host} "test -f $PWD/${_brg_outdir}/short_expected_genome_size.txt"; then
      scp -p ${_local_host}:"$PWD/${_brg_outdir}/short_expected_genome_size.txt" "${_brg_outdir}"
    fi
    if ssh ${_local_host} "test -f $PWD/${_brg_outdir}/o/short_expected_genome_size.txt"; then
      mkdir -p "${_brg_outdir}/o"
      scp -p ${_local_host}:"$PWD/${_brg_outdir}/o/short_expected_genome_size.txt" "${_brg_outdir}/o"
    fi

  else
    if is_file_at_least_1MB "${_brg_outdir_i}/ont.asm.ec.fq"; then
      echo "found: ${_brg_outdir_i}/ont.asm.ec.fq"
    else
      echo "no such file (greater than 1MB) ${_brg_outdir_i}/ont.asm.ec.fq"
      echo "We need the polished long-read data for oatk."
      return 1
    fi
  fi
  return 0
}

check-short-expected-genome-size_genus_species() {
  if [[ "${_local_host}" != "$(hostname)" ]]; then
    if ssh ${_local_host} "test -f $PWD/${_brg_outdir}/short_expected_genome_size.txt"; then
      scp -p ${_local_host}:"$PWD/${_brg_outdir}/short_expected_genome_size.txt" "${_brg_outdir}"
    fi
    if ssh ${_local_host} "test -f $PWD/${_brg_outdir}/o/short_expected_genome_size.txt"; then
      mkdir -p "${_brg_outdir}/o"
      scp -p ${_local_host}:"$PWD/${_brg_outdir}/o/short_expected_genome_size.txt" "${_brg_outdir}/o"
    fi
  else
    if test -s "$PWD/${_brg_outdir}/short_expected_genome_size.txt"; then
      echo "found: ${_brg_outdir}/short_expected_genome_size.txt"
    elif test -s "$PWD/${_brg_outdir}/o/short_expected_genome_size.txt"; then
      echo "found: ${_brg_outdir}/o/short_expected_genome_size.txt"
    else
      echo "We need the expected genome size."
      return 1
    fi
  fi
  return 0
}

oatk-polished_genus_species_for() {
  local _brg_outdir="${1}"
  local _brg_inum="${2:-0}"
  local _brg_polished="${3:-hifiasm}"
  local _brg_fc="${4:-30}"

  local target_index="${_brg_outdir}-${_brg_inum}"
  local long_sra="${_long["$target_index"]}"
  if [[ -z "$long_sra" ]]; then
    echo "Info: skipping oatk-polished on ${_brg_outdir}-${_brg_inum} because it is not in the CSV."
    return
  fi
  local short_sra="${_short["$target_index"]}"
  local _brg_outdir_i="${_brg_outdir}/${_brg_inum}"
  local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

  # echo "long_sra: ${long_sra}"
  # echo "short_sra: ${short_sra}"

  # Step 0. Get the archive data file if this runs at a remote
  # where we do not have the file: cns.fa
  mkdir -p "${_brg_outdir_i}"
  local _brg_corrected_fa="ont.asm.ec.fq"
  if [[ "${_brg_polished}" == h* ]]; then
    _brg_polished="hifiasm"
  elif [[ "${_brg_polished}" == n* ]]; then
    _brg_polished="nextdenovo"
  else
    echo "Error: option polished must either hifiasm or nextdenovo."
    return 1
  fi

  if [[ "${_brg_polished}" == h* ]]; then
    if ! hifiasm-check-cns_genus_species; then
      echo "check ont.asm.ec.fq failed"
      return
    fi
  elif [[ "${_brg_polished}" == n* ]]; then
    _brg_corrected_fa="cns.fa"
    if ! nextdenovo-check-cns_genus_species; then
      echo "check cns.fa failed"
      return
    fi
  else
    echo "Error: option polished must either hifiasm or nextdenovo."
    return 1
  fi

  # Step 2. genome size
  genome_size=$(get_genome_size) || exit 1

  source "$(conda info --base)/etc/profile.d/conda.sh"
  if [[ "$CONDA_DEFAULT_ENV" != "oatk" ]]; then
    echo "You're not in the oatk environment. Chaniging 'oatk'..."
    conda activate oatk
  fi

  if [[ "$CONDA_DEFAULT_ENV" == "oatk" ]]; then

    local fc_list=()

    if [[ "${_brg_fc}" == *,* ]]; then
      IFS=',' read -ra fc_list <<<"$_brg_fc"
    else
      fc_list=("$_brg_fc")
    fi

    for _brg_fc in "${fc_list[@]}"; do
      echo "oatk on ${_brg_outdir}/${_brg_inum} using the ${_brg_polished}-polished data: ${_brg_outdir_i}/${_brg_corrected_fa} with -c ${_brg_fc}"
      local _timing_oatk=${_brg_outdir_i}/timing-oatk-polished-${_brg_polished}-${_brg_fc}.txt
      local _stdout_oatk=${_brg_outdir_i}/sdtout-oatk-polished-${_brg_polished}-${_brg_fc}.txt
      mkdir -p "${_brg_outdir_i}"
      local _outdir_oatk="${_brg_outdir}"-oatk-polished-${_brg_polished}-${_brg_fc}
      mkdir -p "${_outdir_oatk}"

      # oatk -k 1001 -c 30 -t 8 --nhmmscan /bin/nhmmscan -m embryophyta_mito.fam -p embryophyta_pltd.fam -o ddAraThal4 ddAraThal4_organelle.hifi.fa.gz
      #
      # tippo options:
      # -p ont \
      # -t ${_brg_threads} \
      #
      # oatk options:
      # --nhmmscan /bin/nhmmscan \

      # oatk -k 1001 -c 30 -t 8 \
      # 	-m ./OatkDB/v20230921/embryophyta_mito.fam \
      # 	-p ./OatkDB/v20230921/embryophyta_pltd.fam \
      # 	-o oatk-1 \
      # 	cns.fa

      local _oatk_dir="${_brg_outdir_i}"/oatk-polished-${_brg_polished}-${_brg_fc}
      rm -rf "${_oatk_dir}"
      command time -v oatk \
        -k 1001 \
        -c ${_brg_fc} \
        -t 8 \
        -m ./OatkDB/v20230921/embryophyta_mito.fam \
        -p ./OatkDB/v20230921/embryophyta_pltd.fam \
        -o "${_outdir_oatk}"/oatk-1 \
        "${_brg_outdir_i}"/"${_brg_corrected_fa}" \
        >"${_stdout_oatk}" \
        2>"${_timing_oatk}"

      # Record the computer system info
      echo "hostname: $(hostname)" >>"${_timing_oatk}"
      free -h >>"${_timing_oatk}"
      lscpu >>"${_timing_oatk}"

      mv "${_outdir_oatk}" "${_oatk_dir}"
      # rm -rf "${_oatk_dir}"
      # mkdir "${_oatk_dir}"

      # what to copy: consider this
      # cp -pr "${_brg_outdir}"-oatk-${_brg_fc}/gfa_result "${_oatk_dir}"
      # rm -rf "${_brg_outdir}"-oatk-${_brg_fc}
    done

    conda deactivate

  else
    echo "ERROR: no oatk conda environment"
  fi

}

oatk-polished_genus_species() {
  local _brg_outdir="${1:-all}"
  local _brg_inum="${2:-0}"
  local _brg_polished="${3:-hifiasm}"
  local _brg_fc="${4:-30}"

  if [[ "${_brg_outdir}" == "all" ]]; then
    for _v1 in "${Sall[@]}"; do
      oatk-polished_genus_species_for "${_v1}" "${@:2}"
    done
  elif [[ "${_brg_outdir}" == "each" ]]; then
    for _v1 in "${Sall[@]}"; do
      oatk-polished_genus_species_for "${_v1}" "${@:2}"
    done
  else
    oatk-polished_genus_species_for "$@"
  fi
}

hifiasm-polish_genus_species() {
  local _brg_outdir="${1:-all}"
  local _brg_inum="${2:-7}"
  local _brg_fc="${3:-1}"

  local target_index="${_brg_outdir}-${_brg_inum}"
  local long_sra="${_long["$target_index"]}"
  local short_sra="${_short["$target_index"]}"
  if [ -z "${long_sra}" ]; then
    echo "ERROR: no long-read SRA ID: ${long_sra}"
    echo "Suggestion: use -c option for a user-povided CSV."
    return
  fi
  local _brg_outdir_i="${_brg_outdir}/${_brg_inum}"
  local _brg_threads="$(cat /proc/cpuinfo | grep -c processor)"
  # local _memlog_csv="${_brg_outdir_i}/memlog-hifiasm-polish.csv"

  # echo "long_sra: ${long_sra}"
  # echo "short_sra: ${short_sra}"

  # Step 0. Get the archive data file if this runs at a remote
  # where we do not have the file.
  if [[ "${_local_host}" != "$(hostname)" ]]; then
    if [[ ! -s "${_brg_outdir}-a.tar.gz" ]]; then
      scp -p ${_local_host}:$PWD/${_brg_outdir}-a.tar.gz .
    fi
  fi

  mkdir -p "${_brg_outdir_i}"

  if [[ ! -s "${_brg_outdir}/polap.log" ]]; then
    _log_echo "no such file: ${_brg_outdir}/polap.log -> recovering the ${_brg_outdir}"
    recover_genus_species "${_brg_outdir}" "${_brg_inum}"
  fi

  # Step 1. prepare input data
  prepare-long-data_genus_species "${_brg_outdir}" "${_brg_inum}"

  source "$(conda info --base)/etc/profile.d/conda.sh"
  if [[ "$CONDA_DEFAULT_ENV" != "oatk" ]]; then
    echo "You're not in the oatk environment. Chaniging 'oatk'..."
    conda activate oatk
  fi

  if [[ "$CONDA_DEFAULT_ENV" == "oatk" ]]; then
    echo "hifiasm on ${_brg_outdir}/${_brg_inum} using the raw ONT data: ${long_sra}.fastq with option -c ${_brg_fc}"
    local _timing_hifiasm=${_brg_outdir_i}/timing-hifiasm-polish-fc-${_brg_fc}.txt
    local _stdout_hifiasm=${_brg_outdir_i}/stdout-hifiasm-polish-fc-${_brg_fc}.txt
    local _outdir_hifiasm="${_brg_outdir}"-hifiasm-polish-${_brg_fc}
    mkdir -p "${_brg_outdir_i}"
    mkdir -p "${_outdir_hifiasm}"

    # hifiasm -t64 --ont -o ONT.asm ONT.read.fastq.gz --write-ec
    command time -v hifiasm \
      -t${_brg_threads} \
      --ont \
      --write-ec \
      -o "${_outdir_hifiasm}"/ont.asm \
      "${long_sra}.fastq" \
      >"${_stdout_hifiasm}" \
      2>"${_timing_hifiasm}"

    # cat "${work_dir}"/02.cns_align/01.seed_cns.sh.work/seed_cns*/cns.fasta >"${_brg_outdir_i}"/cns.fa

    # Copy the hifiasm run configuration
    # cp "${input_file}" "${_brg_outdir_i}"
    # cp "${hifiasm_cfg}" "${_brg_outdir_i}"

    # Record the computer system info
    echo "hostname: $(hostname)" >>"${_timing_hifiasm}"
    free -h >>"${_timing_hifiasm}"
    lscpu >>"${_timing_hifiasm}"

    local _hifiasm_dir="${_brg_outdir_i}/hifiasm-polish-${_brg_fc}"
    cp -p "${_outdir_hifiasm}"/ont.asm.ec.fq "${_hifiasm_dir}"
    # mv "${_outdir_hifiasm}"/ont.asm.ec.fq "${_hifiasm_dir}"
    # rm -rf "${_outdir_hifiasm}"

    conda deactivate

    # Clean up the hifiasm working folder
    # rm -rf "${work_dir}"
    # rm "${long_sra}".fastq

    if [[ "${_local_host}" != "$(hostname)" ]]; then
      rsync -aq "${_brg_outdir_i}/" "${_local_host}:${PWD}/${_brg_outdir}/"
    fi
  else
    echo "ERROR: no oatk conda environment"
  fi

}

mitohifi_genus_species() {
  local version="${1:-master}"
  shift
  local sif_file="mitohifi-${version}.sif"
  local image_url="docker://ghcr.io/marcelauliano/mitohifi:${version}"
  local bind_dir="${MITOHIFI_BIND_DIR:-$PWD}"
  local tmpdir="${MITOHIFI_TMPDIR:-$(mktemp -d -p "$bind_dir" mitohifi_tmp_XXXX)}"
  local timestamp=$(date '+%Y%m%d_%H%M%S')
  local args=("$@")

  if ! command -v apptainer &>/dev/null; then
    echo "[!] Apptainer not installed." >&2
    return 1
  fi

  if [[ ! -f "$sif_file" ]]; then
    echo "[*] Building SIF image: $sif_file"
    apptainer build "$sif_file" "$image_url" || {
      echo "[!] Failed to build image." >&2
      return 1
    }
  fi

  # Auto-detect input file
  if [[ ! "${args[*]}" =~ -r ]]; then
    local fq
    fq=$(find "$bind_dir" -maxdepth 1 -type f \( -iname "*.fastq" -o -iname "*.fastq.gz" \) | head -n 1)
    if [[ -n "$fq" ]]; then
      echo "[*] Auto-detected input: $fq"
      args+=("-r" "$fq")
    else
      echo "[!] No input FASTQ found and none specified with -r." >&2
      return 1
    fi
  fi

  # Auto-create timestamped output dir if -o not given
  local output_dir
  if [[ "${args[*]}" =~ -o[[:space:]] ]]; then
    output_dir=$(echo "${args[@]}" | sed -n 's/.*-o[[:space:]]\([^[:space:]]*\).*/\1/p')
  else
    output_dir="mitohifi_out_${timestamp}"
    args+=("-o" "$output_dir")
  fi

  mkdir -p "$output_dir"
  local log_file="$output_dir/mitohifi.log"
  local timing_file="$output_dir/mitohifi.timing.txt"

  echo "[*] Running MitoHiFi with output: $output_dir"
  echo "$(date '+%F %T') | mitohifi.py ${args[*]}" >>"$log_file"

  echo "[*] Profiling runtime and memory..."
  {
    /usr/bin/time -v \
      APPTAINER_TMPDIR="$tmpdir" apptainer exec \
      --bind "$bind_dir:$bind_dir" \
      "$sif_file" mitohifi.py "${args[@]}"
  } 2>>"$timing_file"

  local exit_code=$?
  rm -rf "$tmpdir"

  local max_mem time_elapsed
  max_mem=$(grep "Maximum resident set size" "$timing_file" | awk '{print $6}')
  time_elapsed=$(grep "Elapsed (wall clock) time" "$timing_file" | awk -F': ' '{print $2}')
  echo "[*] Summary:"
  echo "    Max memory   : ${max_mem} KB"
  echo "    Elapsed time : ${time_elapsed}"

  echo "$(date '+%F %T') | Max memory: ${max_mem} KB | Time: ${time_elapsed}" >>"$log_file"
  echo "[*] Logs saved to: $log_file and $timing_file"

  if [[ $exit_code -eq 0 ]]; then
    echo "[✔] MitoHiFi completed successfully"
  else
    echo "[!] MitoHiFi failed (exit code $exit_code)"
  fi

  return $exit_code
}

install-apptainer_genus_species() {
  local version="${1:-1.4.0}"
  local deb_url="https://github.com/apptainer/apptainer/releases/download/v${version}/apptainer_${version}_amd64.deb"
  # https://github.com/apptainer/apptainer/releases/download/v1.4.0/apptainer_1.4.0_amd64.deb
  local deb_file="apptainer_${version}_amd64.deb"

  echo "[*] Attempting to install Apptainer v${version} via .deb..."

  if wget -q --show-progress "$deb_url" -O "$deb_file"; then
    if sudo apt install -y "./$deb_file"; then
      echo "[✔] Apptainer installed successfully via .deb"
      rm -f "$deb_file"
      apptainer version
      return 0
    else
      echo "[!] .deb install failed, will build from source..."
      rm -f "$deb_file"
    fi
  else
    echo "[!] Failed to download .deb, will build from source..."
  fi

  return

  echo "[*] Installing build dependencies..."
  sudo apt update
  sudo apt install -y build-essential uuid-dev libgpgme11-dev libseccomp-dev \
    pkg-config squashfs-tools cryptsetup runc git wget

  if ! command -v go >/dev/null 2>&1; then
    echo "[*] Installing Go..."
    wget -q https://go.dev/dl/go1.21.6.linux-amd64.tar.gz
    sudo tar -C /usr/local -xzf go1.21.6.linux-amd64.tar.gz
    export PATH="/usr/local/go/bin:$PATH"
    echo 'export PATH="/usr/local/go/bin:$PATH"' >>~/.bashrc
    rm -f go1.21.6.linux-amd64.tar.gz
  else
    echo "[*] Go is already installed: $(go version)"
  fi

  echo "[*] Cloning Apptainer source..."
  git clone https://github.com/apptainer/apptainer.git
  cd apptainer || {
    echo "[!] Failed to enter apptainer dir"
    return 1
  }
  git checkout "v${version}"

  echo "[*] Building Apptainer..."
  ./mconfig && make -C builddir && sudo make -C builddir install

  echo "[✔] Apptainer installed from source:"
  apptainer version
}

install-novoplasty_genus_species() {
  local _brg_outdir="${1}"
  echo "Preparing archive for ${_brg_outdir} ..."
  # tar zcf "${_brg_outdir}-a.tar.gz" "${_brg_outdir}"
}

install-mitohifi_genus_species() {
  local version="${1:-master}"
  local sif_name="mitohifi-${version}.sif"
  local image_url="docker://ghcr.io/marcelauliano/mitohifi:${version}"

  if [[ -f "$sif_name" ]]; then
    echo "[✔] MitoHiFi SIF already exists: $sif_name"
    return 0
  fi

  echo "[*] Pulling MitoHiFi image version: ${version}"
  echo "[*] Saving as: $sif_name"

  if ! command -v apptainer &>/dev/null; then
    echo "[!] Apptainer is not installed. Please install it first."
    return 1
  fi

  apptainer build "$sif_name" "$image_url" || {
    echo "[!] Failed to build MitoHiFi image"
    return 1
  }

  echo "[✔] MitoHiFi image installed: $sif_name"
}

data-peek_genus_species() {
  local _brg_outdir="$1"
  local _brg_inum="${2:-0}"
  local target_index="${_brg_outdir}-${_brg_inum}"

  local species_name="$(echo ${_brg_outdir} | sed 's/_/ /')"
  local long_sra="${_long["$target_index"]}"
  local short_sra="${_short["$target_index"]}"

  echo "Species: ${species_name}"
  echo "Long SRA ID: ${long_sra}"
  echo "Short SRA ID: ${short_sra}"
}

data-long_genus_species() {
  local _brg_outdir="$1"
  local _brg_inum="${2:-0}"
  local target_index="${_brg_outdir}-${_brg_inum}"

  local species_name="$(echo ${_brg_outdir} | sed 's/_/ /')"
  local long_sra="${_long["$target_index"]}"
  local _brg_outdir_i="${_brg_outdir}/${_brg_inum}"

  # check if long_sra is empty
  if [[ -z "$long_sra" ]]; then
    echo "Error: no long SRA for ${_brg_outdir}-${_brg_inum}"
    return
  fi

  mkdir -p "${_brg_outdir}/${_brg_inum}"

  # Step 1. prepare input fastq data
  local long_data="${_media_dir}/${long_sra}.fastq.tar.gz"
  local long_fastq="${long_sra}.fastq"

  if [[ "${_local_host}" == "$(hostname)" ]]; then
    if [[ -s "${long_fastq}" ]]; then
      echo "found: ${long_fastq}"
    else
      if [[ -s "${long_data}" ]]; then
        tar -zxf ${long_data}
      else
        ln -s "${_media1_dir}/${long_sra}.fastq"
      fi
    fi
  else
    if [[ -s "${long_sra}.fastq" ]]; then
      echo "found: ${long_fastq}"
    else
      if ssh ${_local_host} "test -f ${long_data}"; then
        if [[ -s "${long_sra}".fastq.tar.gz ]]; then
          echo "found: ${long_fastq}.fastq.tar.gz"
        else
          scp ${_local_host}:${long_data} .
        fi
        echo "decompressing ..."
        tar -zxf "${long_sra}".fastq.tar.gz
        rm -f "${long_sra}".fastq.tar.gz
      elif ssh ${_local_host} "test -f ${_media1_dir}/${long_sra}.fastq"; then
        scp ${_local_host}:"${_media1_dir}/${long_sra}.fastq" .
      else
        echo "  downloading long-read SRA ID: ${long_sra} ... be patient!"
        "${_polap_script_bin_dir}"/polap-ncbitools fetch sra "$long_sra"
      fi
    fi
  fi
}

data-short_genus_species() {
  local _brg_outdir="$1"
  local _brg_inum="${2:-0}"
  local target_index="${_brg_outdir}-${_brg_inum}"

  local species_name="$(echo ${_brg_outdir} | sed 's/_/ /')"
  local short_sra="${_short["$target_index"]}"
  local _brg_outdir_i="${_brg_outdir}/${_brg_inum}"

  # check if short_sra is empty
  if [[ -z "$short_sra" ]]; then
    echo "Error: no short SRA for ${_brg_outdir}-${_brg_inum}"
    return
  fi

  mkdir -p "${_brg_outdir}/${_brg_inum}"

  # Step 1. prepare input fastq data
  local short_data="${_media_dir}/${short_sra}.fastq.tar.gz"
  local short_fastq="${short_sra}_1.fastq"

  if [[ "${_local_host}" == "$(hostname)" ]]; then
    if [[ -s "${short_fastq}" ]]; then
      echo "found: ${short_fastq}"
    else
      if [[ -s "${short_data}" ]]; then
        tar -zxf ${short_data}
      else
        ln -s "${_media1_dir}/${short_sra}_1.fastq"
        ln -s "${_media1_dir}/${short_sra}_2.fastq"
      fi
    fi
  else
    if [[ -s "${short_sra}_1.fastq" ]]; then
      echo "found: ${short_fastq}"
    else
      if ssh ${_local_host} "test -f ${short_data}"; then
        if [[ -s "${short_sra}".fastq.tar.gz ]]; then
          echo "found: ${short_fastq}.fastq.tar.gz"
        else
          scp ${_local_host}:${short_data} .
        fi
        echo "decompressing ..."
        tar -zxf "${short_sra}".fastq.tar.gz
      elif ssh ${_local_host} "test -f ${_media1_dir}/${short_sra}_1.fastq"; then
        scp ${_local_host}:"${_media1_dir}/${short_sra}_1.fastq" .
        scp ${_local_host}:"${_media1_dir}/${short_sra}_2.fastq" .
      else
        echo "  downloading short-read SRA ID: ${short_sra} ... be patient!"
        "${_polap_script_bin_dir}"/polap-ncbitools fetch sra "$short_sra"
      fi
    fi
  fi
}

estimate-genomesize_genus_species() {
  local _brg_outdir="${1:-all}"
  local _brg_inum="${2:-0}"

  local target_index="${_brg_outdir}-${_brg_inum}"
  local long_sra="${_long["$target_index"]}"
  local short_sra="${_short["$target_index"]}"

  local _genome_size="0"
  if [[ -s "${_brg_outdir}"/o/short_expected_genome_size.txt ]]; then
    _genome_size=$(<"${_brg_outdir}"/o/short_expected_genome_size.txt)
  elif [[ -s "${_brg_outdir}"/short_expected_genome_size.txt ]]; then
    _genome_size=$(<"${_brg_outdir}"/short_expected_genome_size.txt)
  else
    ${_polap_cmd} find-genome-size \
      -a ${short_sra}_1.fastq \
      -b ${short_sra}_2.fastq \
      -o "${_brg_outdir}"
    if [[ -s "${_brg_outdir}"/short_expected_genome_size.txt ]]; then
      _genome_size=$(<"${_brg_outdir}"/short_expected_genome_size.txt)
    fi
  fi
  _genome_size=${_genome_size%%.*}
  echo "${_genome_size}"
}

data-downsample-long_genus_species() {
  local _brg_outdir="$1"
  local _brg_inum="${2:-0}"
  local _brg_coverage="${3:-50}"
  local _brg_dry="${4:-off}"

  local target_index="${_brg_outdir}-${_brg_inum}"
  local long_sra="${_long["$target_index"]}"
  local random_seed="${_random_seed["$target_index"]}"
  local _brg_outdir_i="${_brg_outdir}/${_brg_inum}"

  local _dry=""
  if [[ "${_brg_dry}" == "on" ]]; then
    _dry="--dry"
  fi

  mkdir -p "${_brg_outdir_i}"

  # Step 1. check input fastq data
  mkdir -p "${_brg_outdir}/tmp"
  if [[ -s "${long_sra}.fastq" ]]; then
    echo "found: ${long_sra}.fastq"
    if [[ -s "${_brg_outdir}/tmp/${long_sra}.fastq" ]]; then
      echo "found and use: ${_brg_outdir}/tmp/${long_sra}.fastq"
    else
      echo mv "${long_sra}.fastq" "${_brg_outdir}"/tmp
      mv "${long_sra}.fastq" "${_brg_outdir}"/tmp
    fi
  else
    if [[ -s "${_brg_outdir}/tmp/${long_sra}.fastq" ]]; then
      echo "found and use: ${_brg_outdir}/tmp/${long_sra}.fastq"
    else
      echo "no such file: ${long_sra}.fastq"
      return
    fi
  fi

  source "$(conda info --base)/etc/profile.d/conda.sh"
  if [[ "$CONDA_DEFAULT_ENV" != "polap" ]]; then
    echo "You're not in the polap environment. Chaniging 'polap'..."
    conda activate polap
  fi

  if [[ "$CONDA_DEFAULT_ENV" == "polap" ]]; then

    # Step 2. Estimate the genome size
    local _genome_size=$(
      estimate-genomesize_genus_species \
        "${_brg_outdir}" \
        "${_brg_inum}" | tail -1
    )
    echo "Genome size estimate: ${_genome_size}"

    # Step 3.
    echo "downsampling ... ${long_sra}"
    echo "input: ${_brg_outdir}/tmp/${long_sra}.fastq"
    echo "output: ${long_sra}.fastq"
    ${_polap_cmd} fastq subsample -v ${_dry} \
      "${_brg_outdir}/tmp/${long_sra}.fastq" \
      "${long_sra}.fastq" \
      -c "${_brg_coverage}" \
      -o "${_brg_outdir_i}" \
      --random-seed "${random_seed}" \
      --genomesize "${_genome_size}" -v \
      >"${_brg_outdir_i}/l${_brg_coverage}x.txt"
    echo "log: ${_brg_outdir_i}/l${_brg_coverage}x.txt"
  fi
}

data-downsample-short_genus_species() {
  local _brg_outdir="$1"
  local _brg_inum="${2:-0}"
  local _brg_coverage="${3:-50}"
  local _brg_dry="${4:-off}"
  local _brg_genomesize="${5}"

  local target_index="${_brg_outdir}-${_brg_inum}"
  local short_sra="${_short["$target_index"]}"
  local random_seed="${_random_seed["$target_index"]}"
  local _brg_outdir_i="${_brg_outdir}/${_brg_inum}"

  local _dry=""
  if [[ "${_brg_dry}" == "on" ]]; then
    _dry="--dry"
  fi

  mkdir -p "${_brg_outdir_i}"

  # Step 1. check input fastq data
  mkdir -p "${_brg_outdir}/tmp"
  rm "${_brg_outdir}/tmp/${short_sra}_1.fastq"
  rm "${_brg_outdir}/tmp/${short_sra}_2.fastq"
  if [[ -s "${short_sra}_1.fastq" ]]; then
    echo "found: ${short_sra}_1.fastq"
  else
    echo "no such file: ${short_sra}_1.fastq"
    return
  fi
  if [[ -s "${short_sra}_2.fastq" ]]; then
    echo "found: ${short_sra}_2.fastq"
  else
    echo "no such file: ${short_sra}_2.fastq"
    return
  fi

  source "$(conda info --base)/etc/profile.d/conda.sh"
  if [[ "$CONDA_DEFAULT_ENV" != "polap" ]]; then
    echo "You're not in the polap environment. Chaniging 'polap'..."
    conda activate polap
  fi

  if [[ "$CONDA_DEFAULT_ENV" == "polap" ]]; then

    local _genome_size=0
    if [[ -z "${_brg_genomesize}" ]]; then
      # Step 2. Estimate the genome size
      _genome_size=$(
        estimate-genomesize_genus_species \
          "${_brg_outdir}" \
          "${_brg_inum}" | tail -1
      )
    else
      _genome_size="${_brg_genomesize}"
    fi
    echo "Genome size estimate: ${_genome_size}"

    # Step 3.
    echo "downsampling ${short_sra} coverage: ${_brg_coverage}x"
    echo "input: ${short_sra}_1.fastq"
    echo "input: ${short_sra}_2.fastq"
    echo "output: ${_brg_outdir}/tmp/${short_sra}_1.fastq"
    echo "output: ${_brg_outdir}/tmp/${short_sra}_2.fastq"

    ${_polap_cmd} fastq subsample2 -v ${_dry} \
      "${short_sra}_1.fastq" \
      "${short_sra}_2.fastq" \
      "${_brg_outdir}/tmp/${short_sra}_1.fastq" \
      "${_brg_outdir}/tmp/${short_sra}_2.fastq" \
      -c "${_brg_coverage}" \
      -o "${_brg_outdir_i}" \
      --random-seed "${random_seed}" \
      --genomesize "${_genome_size}" -v \
      >"${_brg_outdir_i}/s${_brg_coverage}x.txt"
    echo "log: ${_brg_outdir_i}/s${_brg_coverage}x.txt"

  fi
}

clean_input_lines() {
  while IFS= read -r line; do
    # Remove leading spaces
    line="${line#"${line%%[![:space:]]*}"}"

    # Remove trailing ')'
    line="${line%)}"

    echo "$line"
  done
}

function _polap_lib_data-execute-common-subcommand {
  # local subcmd1="$1"
  # local arg2="$2"
  # local arg3="$3"
  # local opt_y_flag="$4"
  local handled=0

  case "$subcmd1" in
  install-conda | setup-conda | uninstall | \
    install-polap | download-polap-github | delete-polap-github | \
    install-fmlrc | patch-polap | bleeding-edge-polap | local-edge-polap | \
    test-polap | \
    install-bandage | install-getorganelle | install-cflye | install-dflye | \
    install-oatk | \
    install-tippo | \
    install-pmat | download-pmat | \
    install-apptainer | \
    install-mitohifi | \
    data-long | data-short | data-peek | \
    data-downsample-long | data-downsample-short | \
    mitohifi | \
    list-subcommands)
    handled=1
    ;;
    ##### INSERT_COMMAND_HERE #####
  merge)
    handled=1
    ;;
  hello-world)
    handled=1
    ;;
  rsync)
    handled=1
    ;;
  oatk-polished)
    handled=1
    ;;
  hifiasm-polish)
    handled=1
    ;;
  mkdir-all | rm-empty | rm)
    handled=1
    ;;
  esac

  case "$subcmd1" in
  data-peek)
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
  data-long)
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
  data-short)
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
  list-subcommands)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 ]]; then
      echo "Help: ${subcmd1} <query> [any|start|end]"
      echo "  ${0} ${subcmd1} polap"
      exit 0
    fi
    [[ "${_arg3}" == arg3 ]] && _arg3="any"
    local _brg_query="${_arg2}"
    local _brg_where="${_arg3}"

    if [[ "${_brg_where}" == "any" ]]; then
      grep ")$" "${_POLAPLIB_DIR}/polap-lib-data.sh" | grep -v "(" | grep "${_brg_query}" | clean_input_lines
      grep ")$" $0 | grep -v "(" | grep "${_brg_query}" | clean_input_lines
    elif [[ "${_brg_where}" == "end" ]]; then
      grep ")$" "${_POLAPLIB_DIR}/polap-lib-data.sh" | grep -v "(" | grep "${_brg_query})$" | clean_input_lines
      grep ")$" $0 | grep -v "(" | grep "${_brg_query})$" | clean_input_lines
    else
      grep ")$" "${_POLAPLIB_DIR}/polap-lib-data.sh" | clean_input_lines | grep -v "(" | grep "^${_brg_query}"
      grep ")$" $0 | clean_input_lines | grep -v "(" | grep "^${_brg_query}"
    fi
    ;;
  install-conda)
    read -p "Do you want to install miniconda3? (y/N): " confirm
    if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
      if [[ -d "$HOME/miniconda3" ]]; then
        echo "ERROR: you already have miniconda3: $HOME/miniconda3"
      else
        mkdir -p ~/miniconda3
        wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
          -O ~/miniconda3/miniconda.sh
        bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
        rm ~/miniconda3/miniconda.sh
        echo After installing Miniconda3, close and reopen your terminal application.
      fi
    else
      echo "miniconda3 installation is canceled."
    fi
    ;;
  setup-conda)
    read -p "Do you want to install miniconda3? (y/N): " confirm
    if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
      # Check if conda is available
      if command -v conda &>/dev/null; then
        source "$(conda info --base)/etc/profile.d/conda.sh"
        conda activate base
        conda config --add channels bioconda
        conda config --add channels conda-forge
        conda config --set channel_priority strict
      else
        echo "Error: conda command not found."
        exit 1
      fi
    else
      echo "miniconda3 bioconda setup is canceled."
    fi
    ;;
  uninstall)
    echo "Removing the following conda environments: polap, polap-fmlrc, getorganelle"
    if [[ "${opt_y_flag}" == false ]]; then
      read -p "Do you want to uninstall all conda environments for polap? (y/N): " confirm
    else
      confirm="yes"
    fi

    if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
      source "$(conda info --base)/etc/profile.d/conda.sh"
      if [[ "$CONDA_DEFAULT_ENV" != "base" ]]; then
        echo "You're not in the base environment. Chaniging 'base'..."
        conda activate base
      fi

      if [[ "$CONDA_DEFAULT_ENV" == "base" ]]; then
        if conda info --envs | awk '{print $1}' | grep -Fxq polap; then
          conda env remove -n polap
          # conda remove -y -n polap --all
        fi
        if conda info --envs | awk '{print $1}' | grep -Fxq polap-fmlrc; then
          conda env remove -n polap-fmlrc
        fi
        if conda info --envs | awk '{print $1}' | grep -Fxq getorganelle; then
          conda env remove -n getorganelle
        fi
        if conda info --envs | awk '{print $1}' | grep -Fxq pmat; then
          conda env remove -n pmat
        fi

        echo "conda deactivate if necessary"
      fi
    else
      echo "Uninstallation is canceled."
      echo "${help_message_uninstall}"
    fi
    ;;
  install-polap)
    if [[ "${opt_y_flag}" == false ]]; then
      read -p "Do you want to install polap? (y/N): " confirm
    else
      confirm="yes"
    fi

    if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
      # Check current conda environment
      # Initialize Conda for non-interactive shells
      source "$(conda info --base)/etc/profile.d/conda.sh"
      if [[ "$CONDA_DEFAULT_ENV" != "base" ]]; then
        echo "You're not in the base environment. Chaniging 'base'..."
        conda activate base
      fi

      if [[ "$CONDA_DEFAULT_ENV" == "base" ]]; then
        echo "You're in the base environment. Creating 'polap'..."
        if conda env list | awk '{print $1}' | grep -qx "polap"; then
          echo "ERROR: Conda environment 'polap' already exists."
        else
          conda create -y --name polap bioconda::polap
        fi
      else
        echo "Error: You're in the '$CONDA_DEFAULT_ENV' environment. Please activate base before running this script."
        exit 1
      fi
    else
      echo "polap installation is canceled."
      echo "${help_message_install_polap}"
    fi
    ;;
  download-polap-github)
    if [[ "${opt_y_flag}" == false ]]; then
      read -p "Do you want to download polap from github? (y/N): " confirm
    else
      confirm="yes"
    fi

    if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
      if [[ -d "polap-${_polap_version}" ]]; then
        echo "ERROR: you already have polap-${_polap_version}"
        echo "  delete it if you want to redownload it."
        exit 0
      fi
      wget -q https://github.com/goshng/polap/archive/refs/tags/${_polap_version}.zip
      unzip -o -q ${_polap_version}.zip
      echo "polap github source is created at polap-${_polap_version}"
    else
      echo "polap download from github is canceled."
    fi
    ;;
  delete-polap-github)
    echo "Deleting the following:"
    echo "  ${_polap_version}.zip and duplicates"
    echo "  polap-${_polap_version}"
    echo "  polap"
    if [[ "${opt_y_flag}" == false ]]; then
      read -p "Do you want to delete all the downloaded polap source? (y/N): " confirm
    else
      confirm="yes"
    fi

    if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
      rm -f ${_polap_version}.zip*
      rm -rf polap-${_polap_version}
      rm -rf polap
      echo "polap download polap-${_polap_version} and its zipped file have been deleted."
    else
      echo "Deleting the polap download from github is canceled."
    fi
    ;;
  install-fmlrc)
    if [[ "${opt_y_flag}" == false ]]; then
      read -p "Do you want to install polap-fmlrc? (y/N): " confirm
    else
      confirm="yes"
    fi

    if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
      # Check current conda environment
      # Initialize Conda for non-interactive shells
      source "$(conda info --base)/etc/profile.d/conda.sh"
      if [[ "$CONDA_DEFAULT_ENV" != "base" ]]; then
        echo "You're not in the base environment. Chaniging 'base'..."
        conda activate base
      fi

      if [[ "$CONDA_DEFAULT_ENV" == "base" ]]; then
        echo "You're in the base environment. Creating 'polap-fmlrc'..."
        if conda env list | awk '{print $1}' | grep -qx "polap-fmlrc"; then
          echo "ERROR: Conda environment 'polap-fmlrc' already exists."
        else
          wget -q https://github.com/goshng/polap/archive/refs/tags/${_polap_version}.zip
          if [[ -s "${_polap_version}.zip" ]]; then
            unzip -o -q ${_polap_version}.zip
            cd polap-${_polap_version}
            conda env create -f src/polaplib/polap-conda-environment-fmlrc.yaml
          else
            echo "Error: no such file: ${_polap_version}.zip - no such polap version"
            echo "Suggestion: _polap_version=0.4.3.7.4 $0 $subcmd1"
          fi
        fi
      else
        echo "Error: You're in the '$CONDA_DEFAULT_ENV' environment. Please activate base before running this script."
        exit 1
      fi
    else
      echo "polap installation is canceled."
      echo "${help_message_install_fmlrc}"
    fi
    ;;
  patch-polap)
    if [[ "${opt_y_flag}" == false ]]; then
      read -p "Do you want to replace conda env polap with the version ${_polap_version}? (y/N): " confirm
    else
      confirm="yes"
    fi

    if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
      if conda env list | awk '{print $1}' | grep -qx "polap"; then
        echo "Updating the Polap Conda environment to version ${_polap_version}."
        wget -q https://github.com/goshng/polap/archive/refs/tags/${_polap_version}.zip
        if [[ -s "${_polap_version}.zip" ]]; then
          unzip -o -q ${_polap_version}.zip
          cd polap-${_polap_version}/src
          bash polaplib/polap-build.sh >../build.sh
          cd ..
          PREFIX="$(conda info --base)/envs/polap" bash build.sh
        else
          echo "Error: no such file: ${_polap_version}.zip - no such polap version"
          echo "Suggestion: _polap_version=0.4.3.7.4 $0 $subcmd1"
        fi
      else
        echo "Error: You do not have polap environment. Please activate polap before running this script."
      fi
    else
      echo "polap patch is canceled."
      echo "${help_message_patch_polap}"
    fi
    ;;
  bleeding-edge-polap)
    if [[ "${opt_y_flag}" == false ]]; then
      read -p "Do you want to update conda env polap with the latest polap github? (y/N): " confirm
    else
      confirm="yes"
    fi

    if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
      if conda env list | awk '{print $1}' | grep -qx "polap"; then
        echo "Updating the Polap Conda environment to its most current state on GitHub ensures access to the latest features and updates."
        if [[ -d "polap" ]]; then
          rm -rf polap
        fi
        git clone --quiet https://github.com/goshng/polap.git
        cd polap/src
        bash polaplib/polap-build.sh >../build.sh
        cd ..
        PREFIX="$(conda info --base)/envs/polap" bash build.sh
      else
        echo "Error: You do not have polap environment. Please activate polap before running this script."
      fi
    else
      echo "polap bleeding-edge version update is canceled."
      echo "${help_message_bleeding_edge_polap}"
    fi
    ;;
  local-edge-polap)
    if [[ "${_arg2}" == arg2 ]]; then
      echo "Help: ${subcmd1} <polap github base path>"
      echo "  $0 ${subcmd1} polap/github"
      exit 0
    fi
    if [[ "${opt_y_flag}" == false ]]; then
      read -p "Do you want to update conda env polap with the local polap source? (y/N): " confirm
    else
      confirm="yes"
    fi

    if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
      if conda env list | awk '{print $1}' | grep -qx "polap"; then
        echo "Updating the Polap Conda environment to the local source at ${_arg2}."
        if [[ -d "${_arg2}/src" ]]; then
          cd "${_arg2}/src"
          bash polaplib/polap-build.sh >../build.sh
          cd ..
          PREFIX="$(conda info --base)/envs/polap" bash build.sh
        else
          echo "Error: The polap base path does not have src directory: ${_arg2}"
        fi
      else
        echo "Error: You do not have polap environment. Please activate polap before running this script."
      fi
    else
      echo "polap local version update is canceled."
    fi
    ;;
  test-polap)
    if [[ "${opt_y_flag}" == false ]]; then
      read -p "Do you want to test polap? (y/N): " confirm
    else
      confirm="yes"
    fi

    if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
      if [[ ! -d "polap-${_polap_version}" ]]; then
        echo "ERROR: no such folder: polap-${_polap_version}"
        echo "  fetch the polap github source"
        echo "$0 download-polap-github"
        exit 0
      fi
      cd polap-${_polap_version}/test
      source "$(conda info --base)/etc/profile.d/conda.sh"
      if [[ "$CONDA_DEFAULT_ENV" != "polap" ]]; then
        echo "You're not in the polap environment. Chaniging 'polap'..."
        conda activate polap
      fi

      if [[ "$CONDA_DEFAULT_ENV" == "polap" ]]; then
        polap assemble --test
      else
        echo "ERROR: no such conda environment: polap"
        echo "$0 install-polap"
      fi
    else
      echo "polap test is canceled."
    fi
    ;;
  install-bandage)
    if [[ "${opt_y_flag}" == false ]]; then
      read -p "Do you want to install bandage? (y/N): " confirm
    else
      confirm="yes"
    fi

    if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
      wget -q https://github.com/rrwick/Bandage/releases/download/v0.8.1/Bandage_Ubuntu_static_v0_8_1.zip
      unzip Bandage_Ubuntu_static_v0_8_1.zip
      mkdir -p bin
      cp Bandage "bin"
      source <(echo 'export PATH="$PWD/bin:$PATH"')
    else
      echo "bandage installation is canceled."
      echo "${help_message_install_bandage}"
    fi
    ;;
  install-getorganelle)
    if [[ "${opt_y_flag}" == false ]]; then
      read -p "Do you want to install getorganelle? (y/N): " confirm
    else
      confirm="yes"
    fi

    if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
      # Check current conda environment
      # Initialize Conda for non-interactive shells
      source "$(conda info --base)/etc/profile.d/conda.sh"
      if [[ "$CONDA_DEFAULT_ENV" != "base" ]]; then
        echo "You're not in the base environment. Chaniging 'base'..."
        conda activate base
      fi

      if [[ "$CONDA_DEFAULT_ENV" == "base" ]]; then
        echo "You're in the base environment. Creating 'getorganelle'..."
        if conda env list | awk '{print $1}' | grep -qx "getorganelle"; then
          echo "ERROR: Conda environment 'getorganelle' already exists."
        else
          conda create -y --name getorganelle bioconda::getorganelle
          conda activate getorganelle
          get_organelle_config.py --add embplant_pt,embplant_mt
        fi
      else
        echo "Error: You're in the '$CONDA_DEFAULT_ENV' environment. Please activate base before running this script."
        exit 1
      fi
    else
      echo "getorganelle installation is canceled."
      echo "${help_message_install_getorganelle}"
    fi
    ;;
  install-cflye)
    if [[ "${opt_y_flag}" == false ]]; then
      read -p "Do you want to install cflye? (y/N): " confirm
    else
      confirm="yes"
    fi

    if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
      # Check current conda environment
      # Initialize Conda for non-interactive shells
      source "$(conda info --base)/etc/profile.d/conda.sh"
      if [[ "$CONDA_DEFAULT_ENV" != "polap" ]]; then
        echo "You're not in the polap environment. Chaniging 'polap'..."
        conda activate polap
      fi

      if [[ "$CONDA_DEFAULT_ENV" == "polap" ]]; then
        echo "You're in the polap environment. Installing 'cflye'..."
        conda install -y goshng::cflye
      else
        echo "Error: You're in the '$CONDA_DEFAULT_ENV' environment. Please activate polap before running this script."
        exit 1
      fi
    else
      echo "cflye or read-coverage filtering version Flye installation is canceled."
    fi
    ;;
  install-dflye)
    if [[ "${opt_y_flag}" == false ]]; then
      read -p "Do you want to install dflye? (y/N): " confirm
    else
      confirm="yes"
    fi

    if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
      # Check current conda environment
      # Initialize Conda for non-interactive shells
      source "$(conda info --base)/etc/profile.d/conda.sh"
      if [[ "$CONDA_DEFAULT_ENV" != "polap" ]]; then
        echo "You're not in the polap environment. Chaniging 'polap'..."
        conda activate polap
      fi

      if [[ "$CONDA_DEFAULT_ENV" == "polap" ]]; then
        echo "You're in the polap environment. Installing 'dflye'..."
        conda install -y goshng::dflye
      else
        echo "Error: You're in the '$CONDA_DEFAULT_ENV' environment. Please activate polap before running this script."
        exit 1
      fi
    else
      echo "dflye or read-coverage filtering version Flye installation is canceled."
    fi
    ;;
  install-tippo)
    if [[ "${opt_y_flag}" == false ]]; then
      read -p "Do you want to install tippo? (y/N): " confirm
    else
      confirm="yes"
    fi

    if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
      # Check current conda environment
      # Initialize Conda for non-interactive shells
      source "$(conda info --base)/etc/profile.d/conda.sh"
      if [[ "$CONDA_DEFAULT_ENV" != "base" ]]; then
        echo "You're not in the base environment. Chaniging 'base'..."
        conda activate base
      fi

      if [[ "$CONDA_DEFAULT_ENV" == "base" ]]; then
        echo "You're in the base environment. Creating 'tippo'..."
        if conda env list | awk '{print $1}' | grep -qx "tippo"; then
          echo "ERROR: Conda environment 'tippo' already exists."
        else
          conda create -n TIPP python=3.8 #please specific the python version to 3.8 :)
          conda activate TIPP
          conda install -y bioconda::tipp
        fi
      else
        echo "Error: You're in the '$CONDA_DEFAULT_ENV' environment. Please activate base before running this script."
        exit 1
      fi
    else
      echo "tippo installation is canceled."
      echo "${help_message_install_tippo}"
    fi
    ;;
  install-oatk)
    if [[ "${opt_y_flag}" == false ]]; then
      read -p "Do you want to install oatk? (y/N): " confirm
    else
      confirm="yes"
    fi

    if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
      # Check current conda environment
      # Initialize Conda for non-interactive shells
      source "$(conda info --base)/etc/profile.d/conda.sh"
      if [[ "$CONDA_DEFAULT_ENV" != "base" ]]; then
        echo "You're not in the base environment. Chaniging 'base'..."
        conda activate base
      fi

      if [[ "$CONDA_DEFAULT_ENV" == "base" ]]; then
        echo "You're in the base environment. Creating 'oatk'..."
        if conda env list | awk '{print $1}' | grep -qx "oatk"; then
          echo "ERROR: Conda environment 'oatk' already exists."
        else
          conda create -y --name oatk bioconda::oatk
          conda activate oatk
          conda install -y biopython
          conda install -y hmmer seqtk mafft parallel entrez-direct
          conda install -y hifiasm
          git clone https://github.com/c-zhou/OatkDB.git
        fi
      else
        echo "Error: You're in the '$CONDA_DEFAULT_ENV' environment. Please activate base before running this script."
        exit 1
      fi
    else
      echo "oatk installation is canceled."
      echo "${help_message_install_oatk}"
    fi
    ;;
  install-bolap)
    if [[ "${opt_y_flag}" == false ]]; then
      read -p "Do you want to install bolap? (y/N): " confirm
    else
      confirm="yes"
    fi

    if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
      # Check current conda environment
      # Initialize Conda for non-interactive shells
      source "$(conda info --base)/etc/profile.d/conda.sh"
      if [[ "$CONDA_DEFAULT_ENV" != "base" ]]; then
        echo "You're not in the base environment. Chaniging 'base'..."
        conda activate base
      fi

      if [[ "$CONDA_DEFAULT_ENV" == "base" ]]; then
        echo "You're in the base environment. Creating 'bolap'..."
        if conda env list | awk '{print $1}' | grep -qx "bolap"; then
          echo "ERROR: Conda environment 'bolap' already exists."
        else
          conda create -y --name bolap
          conda activate bolap
          conda install -y fqtools seqkit seqtk entrez-direct
        fi
      else
        echo "Error: You're in the '$CONDA_DEFAULT_ENV' environment. Please activate base before running this script."
        exit 1
      fi
    else
      echo "bolap installation is canceled."
      echo "${help_message_install_bolap}"
    fi
    ;;
  download-pmat)
    if [[ "${opt_y_flag}" == false ]]; then
      read -p "Do you want to download pmat? (y/N): " confirm
    else
      confirm="yes"
    fi

    if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
      wget https://github.com/bichangwei/PMAT/archive/refs/tags/v1.5.3.tar.gz
      tar -zxvf v1.5.3.tar.gz
      source <(echo 'export PATH="$PWD/PMAT-1.5.3/bin:$PATH"')
      cd PMAT-1.5.3/bin
      chmod +x PMAT
    else
      echo "pmat download is canceled."
    fi
    ;;
  install-pmat)
    if [[ "${opt_y_flag}" == false ]]; then
      read -p "Do you want to install pmat? (y/N): " confirm
    else
      confirm="yes"
    fi

    if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
      # Check current conda environment
      # Initialize Conda for non-interactive shells
      source "$(conda info --base)/etc/profile.d/conda.sh"
      if [[ "$CONDA_DEFAULT_ENV" != "base" ]]; then
        echo "You're not in the base environment. Chaniging 'base'..."
        conda activate base
      fi

      if [[ "$CONDA_DEFAULT_ENV" == "base" ]]; then
        echo "You're in the base environment. Creating 'pmat'..."
        if conda env list | awk '{print $1}' | grep -qx "pmat"; then
          echo "ERROR: Conda environment 'pmat' already exists."
        else
          conda create -y --name pmat apptainer nextdenovo canu blast
          # conda activate pmat
          # conda install -y conda-forge::apptainer
          # conda install -y bioconda::canu
          # conda install -y bioconda::nextdenovo
          # conda install -y bioconda::blast
          # conda install -y bioconda::aspera-cli
        fi
      else
        echo "Error: You're in the '$CONDA_DEFAULT_ENV' environment. Please activate base before running this script."
        exit 1
      fi
    else
      echo "pmat installation is canceled."
    fi
    ;;
  ##### INSERT_CASE_HERE #####
  merge)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 ]]; then
      echo "Help: ${subcmd1} <all|outdir> [source] [dest] [dest_dir]"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    [[ "${_arg3}" == arg3 ]] && _arg3="/media/h2/goshng/aflye1"
    [[ "${_arg4}" == arg4 ]] && _arg4="/media/h3/labshare/goshng/aflye2"
    [[ "${_arg5}" == arg5 ]] && _arg5="s1"
    ${subcmd1}_genus_species "${_arg2}" "${_arg3}" "${_arg4}" "${_arg5}"
    ;;
  hello-world)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 ]]; then
      echo "Help: ${subcmd1} <outdir>"
      echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    [[ "${_arg3}" == arg3 ]] && _arg3=""
    ${subcmd1}_genus_species "${_arg2}"
    ;;
  rsync)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 ]]; then
      echo "Help: ${subcmd1} <outdir|all> [inum:0|N]"
      echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    [[ "${_arg3}" == arg3 ]] && _arg3=""
    ${subcmd1}_genus_species "${_arg2}" "${_arg3}"
    ;;
  oatk-polished)
    if [[ "${_arg2}" == arg2 ]]; then
      echo "Help: ${subcmd1} <all|outdir> [inum:0|N] [polished:hifiasm|nextdenovo] [fc:30|N|N1,N2,N3,...]"
      echo "  polap-data-v2.sh ${subcmd1} all"
      echo "  polap-data-v2.sh ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    [[ "${_arg3}" == arg3 ]] && _arg3="0"
    [[ "${_arg4}" == arg4 ]] && _arg4="hifiasm"
    [[ "${_arg5}" == arg5 ]] && _arg5="30"
    ${subcmd1}_genus_species "${_arg2}" "${_arg3}" "${_arg4}" "${_arg5}"
    ;;
  hifiasm-polish)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 ]]; then
      echo "Help: ${subcmd1} <outdir> [inum:7|N]"
      echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    [[ "${_arg3}" == arg3 ]] && _arg3=""
    ${subcmd1}_genus_species "${_arg2}" "${_arg3}"
    ;;
  mitohifi)
    if [[ "${_arg2}" == "-h" || "${_arg2}" == --h* ]]; then
      echo "Help: ${subcmd1} [master]"
      echo "  ${0} ${subcmd1}"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    [[ "${_arg2}" == arg2 ]] && _arg2="master"
    ${subcmd1}_genus_species "${_arg2}" "${@:3}"
    ;;
  install-apptainer)
    if [[ "${_arg2}" == "-h" || "${_arg2}" == --h* ]]; then
      echo "Help: ${subcmd1} [version:1.4.0]"
      echo "  ${0} ${subcmd1}"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    [[ "${_arg2}" == arg2 ]] && _arg2="1.4.0"
    ${subcmd1}_genus_species "${_arg2}"
    ;;
  install-novoplasty)
    if [[ "${opt_y_flag}" == false ]]; then
      read -p "Do you want to install novoplasty? (y/N): " confirm
    else
      confirm="yes"
    fi

    if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
      # Check current conda environment
      # Initialize Conda for non-interactive shells
      source "$(conda info --base)/etc/profile.d/conda.sh"
      if [[ "$CONDA_DEFAULT_ENV" != "base" ]]; then
        echo "You're not in the base environment. Chaniging 'base'..."
        conda activate base
      fi

      if [[ "$CONDA_DEFAULT_ENV" == "base" ]]; then
        echo "You're in the base environment. Creating 'novoplasty'..."
        if conda env list | awk '{print $1}' | grep -qx "novoplasty"; then
          echo "ERROR: Conda environment 'novoplasty' already exists."
        else
          conda create -y --name novoplasty bioconda::novoplasty
        fi
      else
        echo "Error: You're in the '$CONDA_DEFAULT_ENV' environment. Please activate base before running this script."
        exit 1
      fi
    else
      echo "novoplasty installation is canceled."
      echo "${help_message_install_novoplasty}"
    fi
    ;;
  install-mitohifi)
    if [[ "${_arg2}" == "-h" || "${_arg2}" == --h* ]]; then
      echo "Help: ${subcmd1} [version:master]"
      echo "  ${0} ${subcmd1}"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    [[ "${_arg2}" == arg2 ]] && _arg2="master"
    ${subcmd1}_genus_species "${_arg2}"
    ;;
  estimate-genomesize)
    if [[ "${_arg2}" == arg2 ]]; then
      echo "Help: ${subcmd1} <outdir|all> <inum:N>"
      echo "  polap-data-v2.sh ${subcmd1} all"
      echo "  polap-data-v2.sh ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    [[ "${_arg3}" == arg3 ]] && _arg3=""
    ${subcmd1}_genus_species "${_arg2}" "${_arg3}"
    ;;
  data-downsample-long)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 ]]; then
      echo "Help: ${subcmd1} <outdir> [inum:0|N] [coverage:50|N] [dry:off|on]"
      echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    [[ "${_arg3}" == arg3 ]] && _arg3=""
    [[ "${_arg4}" == arg4 ]] && _arg4=""
    [[ "${_arg5}" == arg5 ]] && _arg5=""
    ${subcmd1}_genus_species "${_arg2}" "${_arg3}" "${_arg4}" "${_arg5}"
    ;;
  data-downsample-short)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 ]]; then
      echo "Help: ${subcmd1} <outdir> [inum:0|N] [coverage:50|N] [dry:off|on] [genomesize]"
      echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    [[ "${_arg3}" == arg3 ]] && _arg3=""
    [[ "${_arg4}" == arg4 ]] && _arg4=""
    [[ "${_arg5}" == arg5 ]] && _arg5=""
    [[ "${_arg6}" == arg6 ]] && _arg6=""
    ${subcmd1}_genus_species "${_arg2}" "${_arg3}" "${_arg4}" "${_arg5}" "${_arg6}"
    ;;
  mkdir-all)
    handled=1
    if [[ "${opt_y_flag}" == false ]]; then
      read -p "Do you want to create all species folders? (y/N): " confirm
    else
      confirm="yes"
    fi

    case "$confirm" in
    [yY] | [yY][eE][sS])
      echo "creating all species folder $i ..."
      for i in "${Sall[@]}"; do
        mkdir ${i}
      done
      ;;
    *)
      echo "Creating all species folders is canceled."
      ;;
    esac
    ;;
  rm-empty)
    handled=1
    find . -type d -empty -delete
    echo "Deleting empty folders ... done."
    ;;
  rm)
    handled=1
    read -p "Do you really want to delete all? (YES/NO): " confirm
    case "$confirm" in
    YES)
      read -p "Really, type in YES? (YES): " confirm
      if [[ "$confirm" == "YES" ]]; then
        for i in "${Sall[@]}"; do
          echo "deleting folder $i ..."
          rm -rf ${i}
        done
      fi
      ;;
    *)
      echo "Deleting all is canceled."
      ;;
    esac
    ;;
  esac

  [[ $handled -eq 1 ]] && return 0 || return 1
}
