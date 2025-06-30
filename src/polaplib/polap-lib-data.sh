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

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  echo "[ERROR] This script must be sourced, not executed: use 'source $BASH_SOURCE'" >&2
  return 1 2>/dev/null || exit 1
fi
: "${_POLAP_DEBUG:=0}"
: "${_POLAP_RELEASE:=0}"

# TODO: Check
# _POLAP_RELEASE
if [[ "${_POLAP_RELEASE}" == "1" ]]; then
  _polap_var_memtracker_time_interval=60
else
  _polap_var_memtracker_time_interval=15
fi

source "${_POLAPLIB_DIR}/polap-lib-conda.sh"
source "${_POLAPLIB_DIR}/polap-lib-timing.sh"
source "${_POLAPLIB_DIR}/polap-lib-unit.sh"
source "${_POLAPLIB_DIR}/polap-lib-array.sh"
source "${_POLAPLIB_DIR}/polap-lib-number.sh"
source "${_POLAPLIB_DIR}/polap-lib-file.sh"
source "${_POLAPLIB_DIR}/polap-lib-process.sh"
source "${_POLAPLIB_DIR}/polap-lib-extract.sh"

: "${_brg_outdir:=.}"
_polap_data_cmd="$(basename "$0" .sh)"
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

# oatk options
declare -a _polap_oatk_options
_polap_oatk_options=(
  30
  20
)

# pmat options
declare -a _polap_pmat_options
_polap_pmat_options=(
  0.1
  1.0
)

# tippo options
declare -a _polap_tippo_options
_polap_tippo_options=(
  onthq
  ont
)

system_genus_species() {
  _polap_lib_timing-get_system_info
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
  conda create -y --name polap-getorganelle bioconda::getorganelle
  conda activate polap-getorganelle
  get_organelle_config.py --add embplant_pt,embplant_mt
HEREDOC
)

help_message_install_cflye=$(
  cat <<HEREDOC

  Install cflye to polap-cflye conda environment
  conda install -y goshng::cflye
HEREDOC
)

help_message_install_polap=$(
  cat <<HEREDOC

  Install conda environments: polap
  conda create -y --name polap bioconda::polap
  conda create -y --name polap bioconda::polap=0.3.7.3
  conda create -y --name polap bioconda::polap=0.4.3.7.7
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
  conda activate polap-oatk
  conda install -y biopython
  conda install -y hmmer seqtk mafft parallel entrez-direct
  git clone https://github.com/c-zhou/OatkDB.git
HEREDOC
)

help_message_uninstall=$(
  cat <<HEREDOC

  Installable tools include:
  conda
  polap
  fmlrc
  getorganelle
  pmat
  tippo
  oatk
  man
  cflye
  dflye

  Conda How-To:
  (better) conda env remove -y -n myenv
  (works the same) conda remove -n myenv --all -y
  e.g., conda env remove -y -n pmat
HEREDOC
)

##### INSERT_HELP_HERE #####
help_message_run_direct_dflye=$(
  cat <<HEREDOC

  menu title
HEREDOC
)

help_message_download_bioproject=$(
  cat <<HEREDOC

  menu title
HEREDOC
)

help_message_archive_run=$(
  cat <<HEREDOC

  menu title
HEREDOC
)

help_message_get_dflye=$(
  cat <<HEREDOC

  menu title
HEREDOC
)

help_message_config_add=$(
  cat <<HEREDOC

  field
  value
  out.csv
HEREDOC
)

help_message_config=$(
  cat <<HEREDOC

  view
  add
HEREDOC
)

help_message_config_view=$(
  cat <<HEREDOC

  menu title
HEREDOC
)

help_message_run_polap_assemble_wga=$(
  cat <<HEREDOC

  menu title
HEREDOC
)

help_message_run_polap_reduce_data=$(
  cat <<HEREDOC

  menu title
HEREDOC
)

help_message_run_polap_prepare_data=$(
  cat <<HEREDOC

  menu title
HEREDOC
)

help_message_run_direct_select_oga=$(
  cat <<HEREDOC

  Select one of oga graphs.
  outdir: species folder, outdir/t2/inum is the actual outdir.
  inum: the number folder. outdir-inum is the key in the CSV config.
  jnum: index of selected oga, out/1/05-flye/ptgaul-intra-base-length/<jnum>
  tdir: t2 can be changed.
HEREDOC
)

help_message_run_direct_bandage_oga=$(
  cat <<HEREDOC

  Create bandage graphs of oga result.
  outdir: species folder, outdir/t2/inum is the actual outdir.
  inum: the number folder. outdir-inum is the key in the CSV config.
  tdir: t2 can be changed.
  bandage: draw bandage graphs.
HEREDOC
)

help_message_run_direct_oga=$(
  cat <<HEREDOC

  Perform an organelle-genome assembly.
  outdir: species folder, outdir/t2/inum is the actual outdir.
  index: the number folder. outdir-inum is the key in the CSV config.
  -t t4: t4 can be changed.
HEREDOC
)

help_message_run_direct_wga=$(
  cat <<HEREDOC

  Perform a whole-genome assembly for seed contig selection.
  outdir: species folder, outdir/t2/inum is the actual outdir.
  index: the number folder. outdir-inum is the key in the CSV config.
  -t t4: t4 can be changed.
HEREDOC
)

help_message_run_direct_flye=$(
  cat <<HEREDOC

  Execute Flye with directional reads. (not tested yet)
  outdir: species folder, outdir/t2/inum is the actual <out>.
  index: the number folder. outdir-inum is the key in the CSV config.

  jnum: index of selected oga, <out>/<junm>
  knum: index of selected oga, <out>/<kunm>
  
  -t t4: t4 can be changed.
HEREDOC
)

help_message_run_direct_read=$(
  cat <<HEREDOC

  Select reads for directional Flye run.
  outdir: species folder, outdir/t2/inum is the actual <out>.
  inum: the number folder. outdir-inum is the key in the CSV config.
  jnum: index of selected oga, <out>/<junm>
  knum: index of selected oga, <out>/<kunm>
  tdir: t2 can be changed.
HEREDOC
)

help_message_run_direct_map=$(
  cat <<HEREDOC

  Map reads on directional seed contigs.
  outdir: species folder, outdir/t2/inum is the actual <out>.
  inum: the number folder. outdir-inum is the key in the CSV config.
  jnum: index of selected oga, <out>/<junm>
  knum: index of selected oga, <out>/<kunm>
  tdir: t2 can be changed.
HEREDOC
)

help_message_run_direct_seed=$(
  cat <<HEREDOC

  Create directional seed contigs.
  outdir: species folder, outdir/t2/inum is the actual <out>.
  inum: the number folder. outdir-inum is the key in the CSV config.
  jnum: index of selected oga, <out>/<junm>
  knum: index of selected oga, <out>/<kunm>
  tdir: t2 can be changed.
HEREDOC
)

help_message_run_direct_dga=$(
  cat <<HEREDOC

  We execute wga, oga, and dga for a directional read assembly.
HEREDOC
)

help_message_install_latex=$(
  cat <<HEREDOC

  menu title
HEREDOC
)

help_message_download_species=$(
  cat <<HEREDOC

  menu title
HEREDOC
)

help_message_setup_completions=$(
  cat <<HEREDOC

  Setup bash command completions:
  aflye
  cflye
  dflye
  taxon

  Source it manually:
  source polap/src/polaplib/polap-data-cflye.complete.sh
HEREDOC
)

help_message_download_test_data_cflye=$(
  cat <<HEREDOC

  Download a test dataset for polap-cflye from:
  https://figshare.com/ndownloader/files/53457569?private_link=ec1cb394870c7727a2d4
  Output1: polap-disassemble-test-full.tar.gz
  Output2: l.fastq, s_1.fastq, s_2.fastq
  Delete these files to redownload or decompress the tar.gz file.
HEREDOC
)

help_message_clean=$(
  cat <<HEREDOC

  Clean-up any output results. <which> can be:
  cflye
HEREDOC
)

help_message_clean_cflye=$(
  cat <<HEREDOC

  Clean-up the output from polap-cflye.
HEREDOC
)

help_message_get_timing_pmat=$(
  cat <<HEREDOC

  menu title
HEREDOC
)

help_message_get_timing=$(
  cat <<HEREDOC

  Get timing for:
  pmat
HEREDOC
)

help_message_recover=$(
  cat <<HEREDOC

  Extract <outdir>-a.tar.gz to <outdir>
HEREDOC
)

help_message_update_local=$(
  cat <<HEREDOC

  Update from a local working directory (e.g., a development version of polap).
HEREDOC
)

help_message_update_github=$(
  cat <<HEREDOC

  Pull and update from the polap GitHub repository.
HEREDOC
)

help_message_update=$(
  cat <<HEREDOC

  Update the tools:
  github
  local
HEREDOC
)

help_message_get_cflye=$(
  cat <<HEREDOC

  Get results cflye1 from a remote host.
  outdir: species folder
  host: remote hostname to download results from
HEREDOC
)

help_message_get=$(
  cat <<HEREDOC

  Subcommands that can follow get are:
  cflye
  timing
HEREDOC
)

help_message_print_all_help=$(
  cat <<HEREDOC

  menu title
HEREDOC
)

help_message_help=$(
  cat <<HEREDOC

  See Also:
  print-help-all     List all help messages.
HEREDOC
)

help_message_install_minimal=$(
  cat <<HEREDOC

  menu title
HEREDOC
)

help_message_install_all=$(
  cat <<HEREDOC

  Install all of the tools including:
  conda
  polap
  fmlrc
  getorganelle
  pmat
  tippo
  oatk
  man
  cflye
  dflye
HEREDOC
)

help_message_remove=$(
  cat <<HEREDOC

  menu title
HEREDOC
)

help_message_run_summary_data=$(
  cat <<HEREDOC

  outdir: species folder
  inum: number folder
HEREDOC
)

help_message_run_polap_check=$(
  cat <<HEREDOC

  Build a plastid genome and check the accuracy via subsampling ONT long-read data.

  simple-polish: default for subsampling-based short-read polishing
  random: off to use CSV's random seed, otherwise on
HEREDOC
)

help_message_run_polap_disassemble=$(
  cat <<HEREDOC

  Build a plastid genome via subsampling ONT long-read data.

  outdir: species folder
  inum: numbber folder
  polishing: default for subsampling-based polishing, simple for ptGAUL's full-data polishing
  random: off to use CSV's random seed, otherwise on
HEREDOC
)

help_message_run_polap=$(
  cat <<HEREDOC

  menu title
HEREDOC
)

help_message_run_oatk_nextdenovo=$(
  cat <<HEREDOC

  outdir: species folder
  inum: number folder
HEREDOC
)

help_message_run_oatk_ont=$(
  cat <<HEREDOC

  outdir: species folder
  inum: number folder
HEREDOC
)

help_message_run_estimate_genomesize=$(
  cat <<HEREDOC

  outdir: species folder
  inum: number folder
HEREDOC
)

help_message_run_nextdenovo_polish=$(
  cat <<HEREDOC

  outdir: species folder
  inum: number folder
HEREDOC
)

help_message_run_pmat=$(
  cat <<HEREDOC

  outdir: species folder
  inum: number folder
HEREDOC
)

help_message_run_tippo=$(
  cat <<HEREDOC

  We use TIPPo v2.4.

  fc: 1
  Execute TIPPo v2.4 on the ONT reads.
  TIPPo.v2.4.pl -f <long_sra>.fastq -p ont
  stderr: o/<inum>/timing-tippo-ont-fc-1.txt
  stdout: o/<inum>/stdout-tippo-ont-fc-1.txt

  fc: 2
  Execute TIPPo v2.3 on the ONT reads.
  TIPPo.v2.3.pl -f <long_sra>.fastq -y --nano-raw -a map-ont

  e.g.,
  TIPPo.v2.3.pl -f ERR2173373.fastq -y --nano-raw -a map-ont
  TIPPo.v2.3.pl -f ERR2173373.fastq -y --nano-raw -a map-ont -m 1500
  TIPPo.v2.4.pl -f ERR2173373.fastq -p ont
  TIPPo.v2.4.pl -f ERR2173373.fastq -p ont -m 1500

  Execute TIPPo v2.4 on the nextDenovo polished reads.
  TIPPo.v2.4.pl -f cns.fa
  stderr: o/<inum>/timing-tippo-nextdenovo-fc-1.txt
  stdout: o/<inum>/stdout-tippo-nextdenovo-fc-1.txt
HEREDOC
)

help_message_run_oatk=$(
  cat <<HEREDOC

  outdir: species folder
  inum: number folder
  type: nextdenovo or ont

  Requirement: git clone https://github.com/c-zhou/OatkDB.git
  Example:
  oatk -k 1001 -c 30 -t 8 --nhmmscan /bin/nhmmscan -m embryophyta_mito.fam -p embryophyta_pltd.fam -o ddAraThal4 ddAraThal4_organelle.hifi.fa.gz
  stderr: o/<inum>/timing-oatk-nextdenovo-fc-1.txt
  stdout: o/<inum>/stdout-oatk-nextdenovo-fc-1.txt
HEREDOC
)

help_message_run_extract_ptdna_ptgaul=$(
  cat <<HEREDOC

  Extract ptDNA from ptGAUL result.
  outdir: species folder
  inum: number folder
HEREDOC
)

help_message_run_polish_ptdna_ptgaul=$(
  cat <<HEREDOC

  Polish ptDNA from ptGAUL result.
HEREDOC
)

help_message_run_ptgaul=$(
  cat <<HEREDOC

  Execute ptGAUL with the long- and short-read data.
HEREDOC
)

help_message_download_ptdna=$(
  cat <<HEREDOC

  Download ptDNA from NCBI.
HEREDOC
)

help_message_run_msbwt=$(
  cat <<HEREDOC

  Execute FMLRC polishing preparation with the short-read data.
HEREDOC
)

help_message_run=$(
  cat <<HEREDOC

  Run the following tools:
  estimate-genomesize
  getorganelle
  msbwt
  ptgaul
  extract-ptdna-ptgaul
  polish-ptdna-ptgaul
  nextdenovo-polish
  pmat
  tippo
  oatk
  polap-disassemble
  polap-disassemble-check
  polap-disassemble-compare
HEREDOC
)

help_message_run_getorganelle=$(
  cat <<HEREDOC

  Execute GetOrganelle with the short-read data.
HEREDOC
)

help_message_install_pmat=$(
  cat <<HEREDOC

  Need: setup pmat
  Ref: wget https://github.com/bichangwei/PMAT/archive/refs/tags/v1.5.3.tar.gz
  conda create -y --name pmat apptainer nextdenovo canu blast
HEREDOC
)

help_message_install_bandage=$(
  cat <<HEREDOC

  Ref: https://github.com/rrwick/Bandage
HEREDOC
)

help_message_download_test_data=$(
  cat <<HEREDOC

  Download the test data.
HEREDOC
)

help_message_download=$(
  cat <<HEREDOC
Usage: $_polap_data_cmd download <subcommand>

Description:
  $_polap_data_cmd is a command-line interface for downloading data.

Subcommands:
  species         Download sequencing data for a given species
  ptdna           Download ptDNA for a given species

Examples:
  $_polap_data_cmd download species <species>

  Download data:
  species <outdir>

  test-data-cflye: polap-cflye
  test-data1 for testing polap-aflye
  test-data2 for testing polap-cflye
  test-data4 for testing polap-dflye
HEREDOC
)

help_message_delete_polap_github=$(
  cat <<HEREDOC

  Delete polap source after installation.
HEREDOC
)

help_message_delete=$(
  cat <<HEREDOC

  polap-github
HEREDOC
)

help_message_install_dflye=$(
  cat <<HEREDOC

  Install dflye to polap-dflye conda environment
  conda install -y goshng::dflye
HEREDOC
)

help_message_setup_conda=$(
  cat <<HEREDOC

  conda config --add channels bioconda
  conda config --add channels conda-forge
  conda config --set channel_priority strict
HEREDOC
)

help_message_setup=$(
  cat <<HEREDOC

  Setup a tool after installation before using it.
  Tools that need a setup include:
  conda
  polap
  pmat
  csv
  completions

  Use setup-<tool> to look up the help of each tool.
HEREDOC
)

help_message_list=$(
  cat <<HEREDOC

  List subcommands.
HEREDOC
)

help_message_install=$(
  cat <<HEREDOC

  Installable tools include:
  conda
  polap
  fmlrc
  getorganelle
  pmat
  tippo
  oatk
  man
  cflye
  dflye

  Some other tools:
  latex
  nvim

  Simply use to install most of the above tools:
  all

  Appened a version to the tool if necessary.
  install polap=0.4.3.7.7
HEREDOC
)

help_message_install_fmlrc=$(
  cat <<HEREDOC

  You could do:
  git clone --quiet https://github.com/goshng/polap.git
  conda env create -f polap/polaplib/polap-conda-environment-fmlrc.yaml
HEREDOC
)

help_message_install_conda=$(
  cat <<HEREDOC

  menu title
HEREDOC
)

help_message_polap_analysis_oga=$(
  cat <<HEREDOC

  p1 <= polap-data-v1.sh
  mkdir Test_species
  cp src/polaplib/polap-data-vt.csv polap-data-v1.csv
  p1 polap-analysis-assemble2 Test_species 10 t1
  p1 polap-analysis-polish Test_species 10 t1

  Note: msbwt is done before oga
  p1 polap-analysis-msbwt Test_species t1 [0]
HEREDOC
)

help_message_polap_analysis_wga=$(
  cat <<HEREDOC

  p1 <= polap-data-v1.sh
  mkdir Test_species
  cp src/polaplib/polap-data-vt.csv polap-data-v1.csv
  p1 polap-analysis-data Test_species
  p1 polap-analysis-reduce Test_species 10 t1 [30]
  p1 polap-analysis-assemble1 Test_species 10 t1 [on]
HEREDOC
)

help_message_nextdenovo_polish=$(
  cat <<HEREDOC

  Polish long-read data using nextDenovo.
HEREDOC
)

help_message_polap_analysis_without_wga=$(
  cat <<HEREDOC

  Copy WGA results rather than doing it.
  Precondition: copy t1 to t1_ using ranger and delete t1.
  p1 <= polap-data-v1.sh
  mkdir Test_species
  cp src/polaplib/polap-data-vt.csv polap-data-v1.csv
  p1 polap-analysis-data Test_species
  p1 polap-analysis-reduce Test_species 10 t1 [30]
  # not doing it: p1 polap-analysis-assemble1 Test_species 10 t1 [on]
  p1 polap-analysis-assemble2 Test_species 10 t1
  p1 polap-analysis-msbwt Test_species t1 [0]
  p1 polap-analysis-polish Test_species 10 t1
HEREDOC
)

help_message_polap_analysis_data=$(
  cat <<HEREDOC

  Prepare the long- and short-read fastq data at the current folder.
  It calls data-short and data-long subcommands:
  p1 data-short Test_species
  p1 data-long Test_species
HEREDOC
)

help_message_polap_analysis_polish=$(
  cat <<HEREDOC

  menu title
HEREDOC
)

help_message_polap_analysis_msbwt=$(
  cat <<HEREDOC

  NOTE: coverage > 0 case is not implemented yet.
  Prepare the FMLRC polishing using all of the short-read data.
  outdir: species folder delimited by an underscore character
  coverage: subsampling the short-read data by the genome size estimate multiplied by this
HEREDOC
)

help_message_polap_analysis_assemble2=$(
  cat <<HEREDOC

  Perform a organelle-genome assembly usnig the seed contigs
  Precondition: polap-analysis-assemble1 must have been executed.
  outdir: species folder delimited by an underscore character
  inum: appended to the 'outdir' directory to match the the first column
    of a CSV config as the key
  analysis: default is . -> o/<inum>, otherwise o/string/<inum>
HEREDOC
)

help_message_polap_analysis=$(
  cat <<HEREDOC

  p1 <= polap-data-v1.sh
  mkdir Test_species
  cp src/polaplib/polap-data-vt.csv polap-data-v1.csv
  p1 polap-analysis-data Test_species
  p1 polap-analysis-reduce Test_species 10 t1 [30]
  p1 polap-analysis-assemble1 Test_species 10 t1 [on]
  p1 polap-analysis-assemble2 Test_species 10 t1
  p1 polap-analysis-msbwt Test_species t1 [0]
  p1 polap-analysis-polish Test_species 10 t1
HEREDOC
)

help_message_polap_analysis_assemble1=$(
  cat <<HEREDOC

  Perform a whole-genome assembly, annotation, and seed selection
  Precondition: input data are needed by polap-analysis-data subcommand.
  Precondition: nk.fq.gz
  Precondition: polap-analysis-reduce must have been executed.
  outdir: species folder delimited by an underscore character
  inum: appended to the 'outdir' directory to match the the first column
    of a CSV config as the key
  analysis: default is . -> o/<inum>, otherwise o/string/<inum>
  plastid: on for plastid seeds
HEREDOC
)

help_message_polap_analysis_reduce=$(
  cat <<HEREDOC

  Prepare the rudeced input data for assemble1 menu.
  Precondition: input data are needed by polap-analysis-data subcommand.
  outdir: species folder delimited by an underscore character
  inum: appended to the 'outdir' directory to match the the first column
    of a CSV config as the key
  analysisdir: default is . -> o/<inum>, otherwise o/string/<inum>
  coverage: polap --coverage option for a whole-genome assembly,
    use downsample of CSV if 0
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

has_help() {
  for arg in "$@"; do
    [[ "$arg" == "-h" || "$arg" == "--help" ]] && return 0
  done
  return 1
}

# setup command auto complete for polap-data-cflye
# 2025-06-07
setup-polap_data_completion() {
  local type="cflye"
  local quiet=0
  local force=0
  local remove=0
  local dry_run=0
  local log_file=""

  # Parse options
  while [[ "$1" ]]; do
    case "$1" in
    --quiet) quiet=1 ;;
    --force) force=1 ;;
    --remove) remove=1 ;;
    --dry-run) dry_run=1 ;;
    --log)
      shift
      log_file="$1"
      ;;
    --type | -t)
      shift
      type="$1"
      ;;
    esac
    shift
  done

  # Derived values
  local bashrc="$HOME/.bashrc"
  local completion_dir="$HOME/.bash_completion.d"
  local completion_path="$completion_dir/_polap-data-${type}"
  local start_marker="# >>> polap-data-${type} initialize >>>"
  local end_marker="# <<< polap-data-${type} initialize <<<"

  # Logging helper
  log_msg() {
    [[ -n "$log_file" ]] && echo "$1" >>"$log_file"
    [[ "$quiet" -eq 0 ]] && echo "$1"
  }

  local script_dir
  script_dir="$(cd -- "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"

  local source_candidates=(
    "$script_dir/polaplib/polap-data-${type}.complete.sh"
    "$script_dir/polap-data-${type}.complete.sh"
  )

  local source_file=""
  for f in "${source_candidates[@]}"; do
    [[ -f "$f" ]] && source_file="$f" && break
  done

  if [[ "$remove" -eq 1 ]]; then
    if [[ "$dry_run" -eq 1 ]]; then
      log_msg "[DRY-RUN] Would remove block from $bashrc"
      log_msg "[DRY-RUN] Would delete $completion_path"
      return 0
    fi
    grep -q "$start_marker" "$bashrc" && {
      if [[ "$quiet" -eq 1 || "$force" -eq 1 ]]; then
        sed -i "/$start_marker/,/$end_marker/d" "$bashrc"
      else
        read -rp "[CONFIRM] Remove block for $type from $bashrc? [y/N] " confirm
        [[ "${confirm,,}" =~ ^y ]] && sed -i "/$start_marker/,/$end_marker/d" "$bashrc"
      fi
    }
    [[ -f "$completion_path" ]] && {
      if [[ "$quiet" -eq 1 || "$force" -eq 1 ]]; then
        rm -f "$completion_path"
        complete -r "polap-data-${type}" 2>/dev/null
      else
        read -rp "[CONFIRM] Delete $completion_path? [y/N] " confirm
        [[ "${confirm,,}" =~ ^y ]] && {
          rm -f "$completion_path"
          complete -r "polap-data-${type}" 2>/dev/null
        }
      fi
    }
    log_msg "[INFO] Reloading $bashrc"
    source "$bashrc"
    return 0
  fi

  if [[ ! -f "$source_file" ]]; then
    log_msg "[ERROR] Could not find polap-data-${type}.complete.sh"
    printf '%s\n' "${source_candidates[@]}" | while read -r c; do log_msg "  - $c"; done
    return 1
  fi

  if [[ "$dry_run" -eq 1 ]]; then
    log_msg "[DRY-RUN] Would copy $source_file to $completion_path"
    log_msg "[DRY-RUN] Would insert block in $bashrc"
    return 0
  fi

  mkdir -p "$completion_dir"
  if [[ "$force" -eq 1 || "$quiet" -eq 1 ]]; then
    cp -f "$source_file" "$completion_path"
  else
    cp -i "$source_file" "$completion_path" || return 1
  fi
  log_msg "[INFO] Copied to $completion_path"

  # Replace block
  sed -i "/$start_marker/,/$end_marker/d" "$bashrc"
  cat >>"$bashrc" <<EOF
$start_marker
[[ -f "$completion_path" ]] && source "$completion_path"
$end_marker
EOF
  log_msg "[INFO] Added block to $bashrc"

  # Validate before sourcing
  if [[ ! -f "$completion_path" ]]; then
    log_msg "[ERROR] Completion file missing: $completion_path"
    return 1
  fi

  if ! bash -n "$completion_path"; then
    log_msg "[ERROR] Syntax error in: $completion_path"
    return 1
  fi

  source "$completion_path"
  log_msg "[INFO] Registered completion for polap-data-${type}"

  source "$bashrc"
  log_msg "[INFO] Reloaded $bashrc"
}

##### INSERT_FUNCTION_HERE #####
download-bioproject_genus_species() {
  local _brg_outdir="${1}"

  "${_polap_cmd}" get-bioproject \
    -o "${_brg_outdir}" \
    --bioproject "${_brg_outdir}"
}

archive-run_genus_species() {
  local _brg_outdir="$1"
  local _brg_sindex="${2:-0}"
  local _brg_adir _brg_title _brg_target _brg_rundir _brg_outdir_i
  local _timing_txt _stdout_txt _memlog_file _summary_file

  _brg_outdir _brg_sindex _brg_adir _brg_title \
    brg_common_setup \
    _brg_target _brg_rundir _brg_outdir_i \
    _timing_txt _stdout_txt _memlog_file _summary_file

  echo "archiving ${_brg_rundir} in ${_brg_outdir_i} ..."
  rsync -azuq --max-size=5M \
    "${_brg_rundir}"/ \
    "${_brg_outdir_i}"/
  rsync -aq "${_brg_rundir}"/0/30-contigger/ \
    "${_brg_outdir_i}"/0/30-contigger/
}

get-dflye_genus_species() {
  local _brg_outdir="$1"
  local _brg_sindex="${2:-0}"
  local _brg_host="${3:-0}"
  local _brg_adir _brg_title _brg_target _brg_rundir _brg_outdir_i
  local _timing_txt _stdout_txt _memlog_file _summary_file

  brg_common_setup \
    _brg_outdir _brg_sindex _brg_adir _brg_title \
    _brg_target _brg_rundir _brg_outdir_i \
    _timing_txt _stdout_txt _memlog_file _summary_file

  # local _brg_host="${_brg_sindex}"

  if [[ "${_brg_outdir}" == "-h" || "${_brg_outdir}" == "--help" ]]; then
    echo "$help_message_get_dflye"
    return
  fi

  if [[ -z "${_brg_outdir}" ]]; then
    echo "[ERROR] outdir is required."
    return
  fi

  if [[ ! -d "${_brg_outdir}" ]]; then
    echo "[ERROR] no such folder: ${_brg_outdir}"
    return
  fi

  if [[ -z "${_brg_host}" ]]; then
    echo "[ERROR] hostname is required".
    return
  fi

  # rsync from the remote
  echo Getting dflye1 results of ${_brg_outdir} from ${_brg_host} ...
  rsync -azuq --max-size=5M \
    "${_brg_host}:${PWD}/${_brg_outdir}"/ \
    "${_brg_outdir}"/
  rsync -azuq \
    "${_brg_host}:${PWD}/${_brg_outdir_i}"/0/30-contigger/ \
    "${_brg_outdir_i}"/0/30-contigger/
  # copy and compress cns.fa
  # rsync -azq "${_brg_host}:${PWD}/${_brg_outdir}/t1/0/cns.fa" "${_brg_outdir}/t1/0/"
  # rsync -azq "${_brg_host}:${PWD}/${_brg_outdir}/t1/0/msbwt/comp_msbwt.npy" "${_brg_outdir}/t1/0/msbwt/"
}

config-add_genus_species() {
  local first_arg="$1"
  local second_arg="${2:-0}"
  local third_arg="${3:-out.csv}"
  print_species_field_summary --add-field="${first_arg}=${second_arg}" \
    --values --out="${third_arg}"
  echo cp -p "${third_arg}" "${_POLAPLIB_DIR}/${_polap_data_csv}"
}

config_genus_species() {
  local first_arg="$1"
  local remaining_args=("${@:2}")

  # Remove trailing slash from the first element
  # remaining_args[0]="${remaining_args[0]%/}"

  config-${first_arg}_genus_species "${remaining_args[@]:-}"
}

config-view_genus_species() {

  # print_species_field_summary --add-field=fruit=banana --values
  print_species_field_summary --values
}

run-polap-assemble-wga_genus_species() {
  local _brg_outdir="$1"
  local _brg_sindex="${2:-0}"
  local _brg_adir _brg_title _brg_target _brg_rundir _brg_outdir_i
  local _timing_txt _stdout_txt _memlog_file _summary_file

  brg_common_setup \
    _brg_outdir _brg_sindex _brg_adir _brg_title \
    _brg_target _brg_rundir _brg_outdir_i \
    _timing_txt _stdout_txt _memlog_file _summary_file

  _polap_lib_conda-ensure_conda_env polap || exit 1

  command time -v "${_polap_cmd}" flye1 \
    -o "${_brg_rundir}" \
    >"${_stdout_txt}" \
    2>"${_timing_txt}"

  command time -v "${_polap_cmd}" edges-stats \
    -o "${_brg_rundir}" \
    >>"${_stdout_txt}" \
    2>>"${_timing_txt}"

  command time -v "${_polap_cmd}" annotate \
    -o "${_brg_rundir}" \
    >>"${_stdout_txt}" \
    2>>"${_timing_txt}"

  local plastid="${_plastid["$_brg_target"]}"
  if [[ "${plastid}" == "false" ]]; then
    command time -v "${_polap_cmd}" seeds \
      -o "${_brg_rundir}" \
      >>"${_stdout_txt}" \
      2>>"${_timing_txt}"
  else
    command time -v "${_polap_cmd}" seeds \
      -o "${_brg_rundir}" \
      --plastid \
      >>"${_stdout_txt}" \
      2>>"${_timing_txt}"
  fi

  conda deactivate
}

run-polap-reduce-data_genus_species() {
  local _brg_outdir="$1"
  local _brg_sindex="${2:-0}"
  local _brg_adir _brg_title _brg_target _brg_rundir _brg_outdir_i
  local _timing_txt _stdout_txt _memlog_file _summary_file

  brg_common_setup \
    _brg_outdir _brg_sindex _brg_adir _brg_title \
    _brg_target _brg_rundir _brg_outdir_i \
    _timing_txt _stdout_txt _memlog_file _summary_file

  local long_sra="${_long["$_brg_target"]}"
  local short_sra="${_short["$_brg_target"]}"
  local coverage="${_coverage["$_brg_target"]}"

  _polap_lib_conda-ensure_conda_env polap || exit 1

  "${_polap_cmd}" total-length-long \
    -l "${long_sra}".fastq \
    -o "${_brg_rundir}"

  "${_polap_cmd}" total-length-short \
    -a "${short_sra}"_1.fastq \
    -b "${short_sra}"_2.fastq \
    -o "${_brg_rundir}"

  local genomesize="${_genomesize["$_brg_target"]}"
  if [[ "${genomesize}" -eq 0 ]]; then
    "${_polap_cmd}" find-genome-size \
      -a "${short_sra}"_1.fastq \
      -b "${short_sra}"_2.fastq \
      -o "${_brg_rundir}"
  else
    echo "${genomesize}" >"${_brg_rundir}/short_expected_genome_size.txt"
  fi

  "${_polap_cmd}" reduce-data \
    -l "${long_sra}".fastq \
    -a "${short_sra}"_1.fastq \
    -b "${short_sra}"_2.fastq \
    -c "${coverage}" \
    -o "${_brg_rundir}"

  conda deactivate
}

run-polap-prepare-data_genus_species() {
  local _brg_outdir="$1"
  local _brg_sindex="${2:-0}"
  local _brg_adir _brg_title _brg_target _brg_rundir _brg_outdir_i
  local _timing_txt _stdout_txt _memlog_file _summary_file

  brg_common_setup \
    _brg_outdir _brg_sindex _brg_adir _brg_title \
    _brg_target _brg_rundir _brg_outdir_i \
    _timing_txt _stdout_txt _memlog_file _summary_file

  if [[ -v _long["$_brg_target"] ]]; then
    local long_sra="${_long["$_brg_target"]}"
  else
    echo "Error: ${_brg_target} because it is not in the CSV."
    return
  fi

  mkdir -p "${_brg_rundir}"

  data-long_genus_species "${_brg_outdir}"
  data-short_genus_species "${_brg_outdir}"
}

run-direct-bandage-oga_genus_species() {
  local _brg_outdir="${1}"
  local _brg_inum="${2:-0}"
  local _brg_adir="${3:-t2}"
  local _brg_bandage="${4:-bandage}"
  local _run_title="run-direct-bandage-oga"

  local target_index="${_brg_outdir}-${_brg_inum}"
  # Check the key exists in the CSV config.
  local long_sra="${_long["$target_index"]}"
  if [[ -z "$long_sra" ]]; then
    echo "Error: skipping ${_brg_outdir}-${_brg_inum} because it is not in the CSV."
    return
  fi
  local extracted_range="${_range["$target_index"]//:/,}"
  local extracted_min_read="${_min_read["$target_index"]}"
  local _outdir="${_brg_outdir}"-"${_brg_inum}"
  local _brg_outdir_i="${_brg_outdir}/${_brg_adir}/${_brg_inum}"
  local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

  # debug local variables
  # for debugging: Inline local printing local var
  while IFS= read -r line; do
    if [[ $line =~ ^declare\ --\ ([^=]+)= ]]; then
      var="${BASH_REMATCH[1]}"
      printf "%s=%q\n" "$var" "${!var}"
    fi
  done < <(local -p 2>/dev/null)

  echo "A: ${extracted_range}"

  # We fix the j to be 1.
  # Test_species-0/t2/0/1/05-flye/ptgaul-intra-base-length/X-graph_final.gfa
  # outdir/1/05-flye/ptgaul-intra-base-length/X-graph_final.gfa
  # outdir/1/01-contig/ptgaul-intra-base-length.txt
  local index_txt="${_outdir}"/1/01-contig/ptgaul-intra-base-length.txt
  cat "${index_txt}"
  local w_values=()
  read -a w_values <"${index_txt}"

  for i in "${w_values[@]}"; do
    echo "$i"
  done

  local _base_figure="."
  local _fc=1
  local _species=Species
  local _brg_csv="1.csv"

  for ((i = 0; i < ${#w_values[@]}; i++)); do
    local _gfa_infer="${_outdir}/1/05-flye/ptgaul-intra-base-length/$i-graph_final.gfa"
    local _png_infer="${_outdir}/1/05-flye/ptgaul-intra-base-length/$i-graph_final.png"
    if [[ -s "${_gfa_infer}" ]]; then
      # echo "gfa file: ${_gfa_infer}"
      if [[ "${_brg_bandage}" == "bandage" ]]; then
        ${_polap_cmd} bandage png \
          ${_gfa_infer} \
          ${_png_infer}
      fi
      printf "%s,%s,%s,%s\n" "${_run_title}" "${_species}" "${_fc}" "${_base_figure}/${_png_infer}" >>"${_brg_csv}"

    else
      printf "%s,%s,%s,%s\n" "${_run_title}" "${_species}" "${_fc}" "${_base_figure}/na.png" >>"${_brg_csv}"

    fi
  done

  return
}

brg_common_setup() {
  local -n outdir_ref="$1"
  local -n sindex_ref="$2"
  local -n adir_ref="$3"
  local -n title_ref="$4"
  local -n target_ref="$5"
  local -n rundir_ref="$6"
  local -n outdir_i_ref="$7"
  local -n timing_txt_ref="$8"
  local -n stdout_txt_ref="$9"
  local -n memlog_file_ref="${10}"
  local -n summary_file_ref="${11}"

  adir_ref="${opt_t_arg:-t4}"
  title_ref="${FUNCNAME[1]#run-}"
  title_ref="${title_ref%%_*}"
  target_ref="${outdir_ref}-${sindex_ref}"
  rundir_ref="${target_ref}"
  outdir_i_ref="${outdir_ref}/${adir_ref}/${sindex_ref}"

  mkdir -p "${outdir_i_ref}"
  timing_txt_ref="${outdir_i_ref}/timing-${title_ref}.txt"
  stdout_txt_ref="${outdir_i_ref}/stdout-${title_ref}.txt"
  memlog_file_ref="${outdir_i_ref}/memlog-${title_ref}.csv"
  summary_file_ref="${outdir_i_ref}/summary-${title_ref}.txt"
}

# create rundir.
run-direct-wga_genus_species() {
  local _brg_outdir="$1"
  local _brg_sindex="${2:-0}"
  local _brg_adir _brg_title _brg_target _brg_rundir _brg_outdir_i
  local _timing_txt _stdout_txt _memlog_file _summary_file

  brg_common_setup \
    _brg_outdir _brg_sindex _brg_adir _brg_title \
    _brg_target _brg_rundir _brg_outdir_i \
    _timing_txt _stdout_txt _memlog_file _summary_file

  run-polap-prepare-data_genus_species "$_brg_outdir" "$_brg_sindex"
  run-polap-reduce-data_genus_species "$_brg_outdir" "$_brg_sindex"
  run-polap-assemble-wga_genus_species "$_brg_outdir" "$_brg_sindex"

  echo "Create seed contigs: ${_brg_rundir}/0/mt.contig.name-1"
}

# use rundir.
run-direct-oga_genus_species() {
  local _brg_outdir="$1"
  local _brg_sindex="${2:-0}"
  local _brg_inum="${3:-0}"
  local _brg_adir _brg_title _brg_target _brg_rundir _brg_outdir_i
  local _timing_txt _stdout_txt _memlog_file _summary_file

  brg_common_setup \
    _brg_outdir _brg_sindex _brg_adir _brg_title \
    _brg_target _brg_rundir _brg_outdir_i \
    _timing_txt _stdout_txt _memlog_file _summary_file

  if [[ -v _long["$_brg_target"] ]]; then
    local long_sra="${_long["$_brg_target"]}"
  else
    echo "Error: ${_brg_target} because it is not in the CSV."
    return
  fi

  local range0="${_range["$_brg_target"]//:/,}"
  local range1="${_range1["$_brg_target"]//:/,}"
  local range2="${_range2["$_brg_target"]//:/,}"
  # if [[ "${_brg_inum}" == "0" ]]; then
  #   local range="${ranger0}"
  # elif [[ "${_brg_inum}" == "1" ]]; then
  #   local range="${ranger1}"
  # fi
  # declare -n ref="ranger${_brg_inum}"
  # local range="$ref"
  local -n range="range${_brg_inum}"

  local min_read="${_min_read["$_brg_target"]}"
  local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

  _polap_lib_conda-ensure_conda_env polap || exit 1

  rm -f "${_timing_txt}"
  rm -f "${_stdout_txt}"

  echo "[INFO] Starting an organelle-genome assembly pipeline on ${_brg_outdir}-${_brg_sindex} on a range of omega values: ${range} and minimum read length of ${min_read}"

  # Start memory logger
  _polap_lib_process-start_memtracker "${_memlog_file}" \
    "${_polap_var_memtracker_time_interval}"

  local _inum="${_brg_inum}"
  for j in {1..9}; do
    _brg_jnum=$((_inum + j))
    if [[ -s "${_brg_rundir}/${_inum}/mt.contig.name-${_brg_jnum}" ]]; then

      ############################################################
      # the output dir should not exist.
      #
      if [[ -d "${_brg_rundir}/${_brg_jnum}" ]]; then
        while true; do
          read -r -p "Folder [${_brg_rundir}/${_brg_jnum}] already exists. Do you want to replace it? [y/n] " yn
          case $yn in
          [Yy]*)
            rm -rf ${_brg_rundir}/${_brg_jnum}
            echo "  folder ${_brg_rundir}/${_brg_jnum} is deleted."
            break
            ;;
          [Nn]*)
            echo "  folder ${_brg_rundir}/${_brg_jnum} is not deleted."
            echo "  the subcommand is cancelled."
            return
            exit $EXIT_FAIL
            ;;
          *) echo "Please answer yes or no." ;;
          esac
        done
      fi

      command time -v "${_polap_cmd}" assemble-wrange \
        -o "${_brg_rundir}" \
        -i "${_inum}" \
        -j "${_brg_jnum}" \
        -s "${range}" \
        -m "${min_read}" \
        >>"${_stdout_txt}" \
        2>>"${_timing_txt}"
      echo "Create seed contigs: ${_brg_rundir}/${_brg_jnum}/05-flye/ptgaul-intra-base-length/j/mt.contig.name"
    fi
  done

  # Summarize results after job (with previously defined summary function)
  _polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"

  conda deactivate

  _polap_lib_timing-get_system_info >>"${_timing_txt}"
}

# We do two things:
# 1. Copy an OGA folder to a run folder
# 2. Manually create directional seed contigs
#
# Select one of the oga assemblies
# there are a few methods in OGA depending on how we choose seeds.
# This is too much complicated.
# ptgaul-intra-base-length <- by test-reads
# polap-reads <- select-reads
# ptgaul-reads <- select-reads
# we may need to specify this method in this subcommand.
# Check with the test-reads or select-reads polap subcommand.
# We just copy these folder to run folder.
#
# We use:
# ptgaul-intra-base-length
run-direct-select-oga_genus_species() {
  local _brg_outdir="$1"
  local _brg_sindex="${2:-0}"
  local _brg_inum="${3:-1}"
  local _brg_jnum="${4:-0}"

  local _brg_adir _brg_title _brg_target _brg_rundir _brg_outdir_i
  local _timing_txt _stdout_txt _memlog_file _summary_file

  brg_common_setup \
    _brg_outdir _brg_sindex _brg_adir _brg_title \
    _brg_target _brg_rundir _brg_outdir_i \
    _timing_txt _stdout_txt _memlog_file _summary_file

  local i="${_brg_inum}"
  local j="${_brg_jnum}"
  local k=$((i + 1))
  # We fix the j to be 1.
  # Test_species-0/t2/0/1/05-flye/ptgaul-intra-base-length/X-graph_final.gfa
  # outdir/1/05-flye/ptgaul-intra-base-length/X-graph_final.gfa
  # outdir/1/01-contig/ptgaul-intra-base-length.txt
  local index_txt="${_brg_rundir}"/$i/01-contig/ptgaul-intra-base-length.txt
  cat "${index_txt}"
  local w_values=()
  read -a w_values <"${index_txt}"

  for w in "${w_values[@]}"; do
    echo "w: $w"
  done

  rm -rf "${_brg_rundir}/$i/30-contigger"
  cp -pr "${_brg_rundir}/$i/05-flye/ptgaul-intra-base-length/$j/30-contigger" "${_brg_rundir}/$i"
  cp -p "${_brg_rundir}/$i/05-flye/ptgaul-intra-base-length/$j/30-contigger/graph_final.gfa" \
    "${_brg_rundir}/$i/assembly_graph.gfa"
  cp -p "${_brg_rundir}/$i/05-flye/ptgaul-intra-base-length/$j/mt.contig.name" \
    "${_brg_rundir}/$i/mt.contig.name-$k"
  echo "Create seed contigs: ${_brg_rundir}/$i/mt.contig.name-#"
  echo ptgaul-intra-base-length/$j >"${_brg_rundir}/$i/index.txt"
  echo "select ptgaul-intra-base-length/$j in $i"
}

# On the run folder we execute:
# seed, map, read, and dflye.
# dga.
# seed and map can be done in one step.
# read might take a while.
# we will rename these:
# seed -> prepare-seeds
# map -> map-reads
# read -> select-reads (a range of omega values)
# flye -> execute-dflye
# We could combine these all into direct-dga.
#
# map-reads could take much time.
# prepare-seeds and map-reads can be combined -> direct-prepare-reads
# select-reads and execute-dflye can be one -> direct-assemble-reads
# if we do not mind of the map-reads processing time, we could combine these
# all into a single subcommand, say direct-dga, on the input of the directional
# seed contigs from direct-oga.
#
# We still need to start with oga.
# We need to be able to execute dflye.
run-direct-seed_genus_species() {
  local _brg_outdir="$1"
  local _brg_sindex="${2:-0}"
  local _brg_inum="${3:-1}"
  local _brg_jnum="${4:-2}"
  local _brg_adir _brg_title _brg_target _brg_rundir _brg_outdir_i
  local _timing_txt _stdout_txt _memlog_file _summary_file

  brg_common_setup \
    _brg_outdir _brg_sindex _brg_adir _brg_title \
    _brg_target _brg_rundir _brg_outdir_i \
    _timing_txt _stdout_txt _memlog_file _summary_file

  _polap_lib_conda-ensure_conda_env polap || exit 1

  # mt.contig.name-1 with plus/minus signed edge_<number>[+-]
  ${_polap_cmd} directional-prepare-seeds \
    -o ${_brg_rundir} \
    -i "${_brg_inum}" -j "${_brg_jnum}"

  conda deactivate
}

run-direct-map_genus_species() {
  local _brg_outdir="$1"
  local _brg_sindex="${2:-0}"
  local _brg_inum="${3:-1}"
  local _brg_jnum="${4:-2}"
  local _brg_adir _brg_title _brg_target _brg_rundir _brg_outdir_i
  local _timing_txt _stdout_txt _memlog_file _summary_file

  brg_common_setup \
    _brg_outdir _brg_sindex _brg_adir _brg_title \
    _brg_target _brg_rundir _brg_outdir_i \
    _timing_txt _stdout_txt _memlog_file _summary_file

  _polap_lib_conda-ensure_conda_env polap || exit 1

  # mt.contig.name-1 with plus/minus signed edge_<number>[+-]
  ${_polap_cmd} directional-map-reads \
    -o ${_brg_rundir} \
    -i "${_brg_inum}" -j "${_brg_jnum}"

  conda deactivate
}

run-direct-read_genus_species() {
  local _brg_outdir="$1"
  local _brg_sindex="${2:-0}"
  local _brg_inum="${3:-1}"
  local _brg_jnum="${4:-2}"
  local _brg_adir _brg_title _brg_target _brg_rundir _brg_outdir_i
  local _timing_txt _stdout_txt _memlog_file _summary_file

  brg_common_setup \
    _brg_outdir _brg_sindex _brg_adir _brg_title \
    _brg_target _brg_rundir _brg_outdir_i \
    _timing_txt _stdout_txt _memlog_file _summary_file

  if [[ -v _long["$_brg_target"] ]]; then
    local long_sra="${_long["$_brg_target"]}"
  else
    echo "Error: ${_brg_target} because it is not in the CSV."
    return
  fi

  local range0="${_range["$_brg_target"]//:/,}"
  local range1="${_range1["$_brg_target"]//:/,}"
  local range2="${_range2["$_brg_target"]//:/,}"
  local -n range="range${_brg_inum}"

  local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

  _polap_lib_conda-ensure_conda_env polap || exit 1

  rm -f "${_timing_txt}"
  rm -f "${_stdout_txt}"

  echo "[INFO] Starting an organelle-genome assembly pipeline on ${_brg_outdir}-${_brg_sindex} on a range of omega values: ${range}"

  _polap_lib_conda-ensure_conda_env polap || exit 1

  ${_polap_cmd} directional-select-reads \
    -o ${_brg_rundir} \
    -i "${_brg_inum}" -j "${_brg_jnum}" \
    -s "${range}"

  conda deactivate
}

run-direct-flye_genus_species() {
  local _brg_outdir="$1"
  local _brg_sindex="${2:-0}"
  local _brg_inum="${3:-1}"
  local _brg_jnum="${4:-2}"
  local _brg_knum="${5:--1}"
  local _brg_adir _brg_title _brg_target _brg_rundir _brg_outdir_i
  local _timing_txt _stdout_txt _memlog_file _summary_file

  brg_common_setup \
    _brg_outdir _brg_sindex _brg_adir _brg_title \
    _brg_target _brg_rundir _brg_outdir_i \
    _timing_txt _stdout_txt _memlog_file _summary_file

  if [[ -v _long["$_brg_target"] ]]; then
    local long_sra="${_long["$_brg_target"]}"
  else
    echo "Error: ${_brg_target} because it is not in the CSV."
    return
  fi

  local range0="${_range["$_brg_target"]//:/,}"
  local range1="${_range1["$_brg_target"]//:/,}"
  local range2="${_range2["$_brg_target"]//:/,}"
  local -n range="range${_brg_inum}"
  local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

  echo "[INFO] Starting an organelle-genome assembly pipeline on ${_brg_outdir}-${_brg_sindex} on a range of omega values: ${range}"

  _polap_lib_conda-ensure_conda_env polap || exit 1

  if [[ "${_brg_knum}" == "-1" ]]; then
    ${_polap_cmd} directional-flye-reads \
      -o ${_brg_rundir} \
      -i "${_brg_inum}" -j "${_brg_jnum}" \
      --no-directional \
      -s "${range}"
  else
    local _next_knum=$((_brg_knum + 1))
    rm -rf "${_brg_rundir}/${_brg_jnum}/05-flye/ptgaul/${_brg_knum}"
    ${_polap_cmd} directional-flye-reads \
      -o ${_brg_rundir} \
      -i "${_brg_inum}" -j "${_brg_jnum}" \
      --no-directional \
      -s "${range}" --start-index "${_brg_knum}" --end-index "${_next_knum}"
  fi

  conda deactivate
}

run-direct-dflye_genus_species() {
  local _brg_outdir="$1"
  local _brg_sindex="${2:-0}"
  local _brg_inum="${3:-1}"
  local _brg_jnum="${4:-2}"
  local _brg_knum="${5:--1}"
  local _brg_adir _brg_title _brg_target _brg_rundir _brg_outdir_i
  local _timing_txt _stdout_txt _memlog_file _summary_file

  brg_common_setup \
    _brg_outdir _brg_sindex _brg_adir _brg_title \
    _brg_target _brg_rundir _brg_outdir_i \
    _timing_txt _stdout_txt _memlog_file _summary_file

  if [[ -v _long["$_brg_target"] ]]; then
    local long_sra="${_long["$_brg_target"]}"
  else
    echo "Error: ${_brg_target} because it is not in the CSV."
    return
  fi

  local range0="${_range["$_brg_target"]//:/,}"
  local range1="${_range1["$_brg_target"]//:/,}"
  local range2="${_range2["$_brg_target"]//:/,}"
  local -n range="range${_brg_inum}"
  local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

  echo "[INFO] Starting an organelle-genome assembly pipeline on ${_brg_outdir}-${_brg_sindex} on a range of omega values: ${range}"

  _polap_lib_conda-ensure_conda_env polap-dflye || exit 1

  if [[ "${_brg_knum}" == "-1" ]]; then
    ${_polap_cmd} directional-flye-reads \
      -o ${_brg_rundir} \
      -i "${_brg_inum}" -j "${_brg_jnum}" \
      --directional \
      -s "${range}"
  else
    local _next_knum=$((_brg_knum + 1))
    rm -rf "${_brg_rundir}/${_brg_jnum}/05-dflye/ptgaul/${_brg_knum}"
    ${_polap_cmd} directional-flye-reads \
      -o ${_brg_rundir} \
      -i "${_brg_inum}" -j "${_brg_jnum}" \
      --directional \
      -s "${range}" --start-index "${_brg_knum}" --end-index "${_next_knum}"
  fi

  conda deactivate
}

# Given a OGA and directional seeds to assemble mtDNA.
# dflye has still bugs that is seg. fault from time to time.
# we use just the directinal reads mapped.
run-direct-dga_genus_species() {
  local _brg_outdir="$1"
  local _brg_sindex="${2:-0}"
  local _brg_inum="${3:-1}"
  local _brg_jnum="${4:-2}"

  # run-direct-select-oga_genus_species
  run-direct-seed_genus_species "${_brg_outdir}" "${_brg_sindex}" "${_brg_inum}" "${_brg_jnum}"
  run-direct-map_genus_species "${_brg_outdir}" "${_brg_sindex}" "${_brg_inum}" "${_brg_jnum}"
  run-direct-read_genus_species "${_brg_outdir}" "${_brg_sindex}" "${_brg_inum}" "${_brg_jnum}"
  run-direct-flye_genus_species "${_brg_outdir}" "${_brg_sindex}" "${_brg_inum}" "${_brg_jnum}"
  # run-direct-dflye_genus_species "${_brg_outdir}" "${_brg_sindex}" "${_brg_inum}" "${_brg_jnum}"
}

install-latex_genus_species() {
  local missing=()

  if ! command -v pdflatex >/dev/null 2>&1; then
    missing+=("pdflatex")
  fi

  if ! command -v biber >/dev/null 2>&1; then
    missing+=("biber")
  fi

  if [[ ${#missing[@]} -eq 0 ]]; then
    echo "[INFO] LaTeX and biber are already installed."
    return 0
  fi

  echo "[INFO] Missing: ${missing[*]}. Installing LaTeX packages..."
  sudo apt-get update
  sudo apt-get install -y \
    texlive \
    texlive-latex-recommended \
    texlive-xetex \
    texlive-fonts-recommended \
    texlive-fonts-extra \
    texlive-lang-all \
    biber
}

download-species_genus_species() {
  local _brg_outdir="${1}"

  # implemented by polap-analysis-data
  polap-analysis-data_genus_species "${_brg_outdir}"
}

setup-completions_genus_species() {
  local _brg_outdir="${1:-cflye}"

  if has_help "${args[@]}"; then
    echo "$help_message_setup_completions"
    return
  fi

  if [[ "${_brg_outdir}" == "cflye" ]]; then
    setup-polap_data_completion
  elif [[ "${_brg_outdir}" == "dflye" ]]; then
    setup-polap_data_completion -t dflye
  fi
}

download-test-data-cflye_genus_species() {
  local _brg_outdir="${1}"

  if [[ "${_brg_outdir}" == "-h" || "${_brg_outdir}" == "--help" ]]; then
    echo "$help_message_download_test_data_cflye"
    return
  fi

  if [[ "${opt_y_flag}" == false ]]; then
    read -p "Do you want to download test data for polap-cflye? (y/N): " confirm
  else
    confirm="yes"
  fi

  if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
    local _s=polap-disassemble-test-full.tar.gz
    # local _s=polap-disassemble-test.tar.gz
    if [[ ! -s "${_s}" ]]; then
      # full test data
      # curl -L -o "${_s}" "https://figshare.com/ndownloader/files/53457569?private_link=ec1cb394870c7727a2d4"
      wget -O "${_s}" "https://figshare.com/ndownloader/files/53457569?private_link=ec1cb394870c7727a2d4"
      #
      # test data
      # curl -L -o "${_s}" "https://figshare.com/ndownloader/files/53457566?private_link=ec1cb394870c7727a2d4"
    fi
    if [[ ! -s "l.fastq" ]]; then
      tar -zxf "${_s}"
      echo "downloaded: l.fastq, s_1.fastq, s_2.fastq"
    else
      echo "You already have: l.fastq"
      echo "  delete l.fastq if you decompress ${_s}"
    fi
  else
    echo "polap test download is canceled."
  fi
}

clean_genus_species() {
  local first_arg="$1"
  local remaining_args=("${@:2}")

  clean-${first_arg}_genus_species "${remaining_args[@]}"
}

clean-cflye_genus_species() {
  local _brg_outdir="${1}"

  if [[ "${_brg_outdir}" == "-h" || "${_brg_outdir}" == "--help" ]]; then
    echo "$help_message_clean_cflye"
    return
  fi

  _brg_outdir="${_brg_outdir%/}"

  if [[ ! -d "${_brg_outdir}" ]]; then
    echo "[ERROR] no such folder: ${_brg_outdir}"
    return 1
  fi

  local target_index="${_brg_outdir}-0"
  local long_sra="${_long["$target_index"]}"
  if [[ -z "$long_sra" ]]; then
    echo "Error: skipping ${_brg_outdir}-${_brg_inum} because it is not in the CSV."
    return
  fi
  local short_sra="${_short["$target_index"]}"

  echo "Deleting ouput files and folders from ${_brg_outdir} ..."
  rm -rf "${_brg_outdir}"
  rm -f "${long_sra}.fastq"
  rm -f "${short_sra}_1.fastq"
  rm -f "${short_sra}_2.fastq"
  rm -f "${short_sra}.fastq.tar.gz"
  rm -f "${_brg_outdir}_0_input.fofn"
  rm -f "${_brg_outdir}-0-nextdenovo-polish.cfg"
  rm -rf "${_brg_outdir}-polap-disassemble"
}

get-timing-pmat_genus_species() {
  local _brg_outdir="${1}"
  local _brg_inum="${2:-0}"
  local _run_title="pmat-nextdenovo"

  local _brg_outdir_t="${_brg_outdir}/${opt_t_arg}"
  local _brg_outdir_i="${_brg_outdir_t}/${_brg_inum}"

  local _memlog_file="${_brg_outdir_i}/memlog-${_run_title}.csv"
  local _summary_file="${_brg_outdir_i}/summary-${_run_title}.txt"

  _polap_lib_process-analyze_memtracker_log "${_memlog_file}" "${_summary_file}"
}

get-timing_genus_species() {
  local args=("$@")

  local item="${args[0]}"
  local rest=("${args[@]:1}")

  case "$item" in
  pmat | oatk | tippo | getorganelle | ptgaul | polap)
    get-timing-${item}_genus_species ${rest[@]}
    ;;
  -h | --help)
    echo $help_message_get_timing
    ;;
  *)
    echo "[ERROR] Invalid timing tool: $item" >&2
    exit 1
    ;;
  esac
}

update-local_genus_species() {
  local local_path="$1"

  if [[ "${local_path}" == "-h" || "${local_path}" == "--help" ]]; then
    echo "$help_message_update_local"
    return
  fi

  local repo_url="https://github.com/goshng/polap.git"
  local github_path="$PWD"
  local conda_env="polap"
  local conda_base
  conda_base=$(conda info --base)
  local target_path="$conda_base/envs/$conda_env/bin/polap"

  if [[ "${opt_y_flag}" == false ]]; then
    echo "[INFO] polap local path: ${local_path}"
    read -p "Do you want to update conda env polap with the local polap source? (y/N): " confirm
  else
    confirm="yes"
  fi

  if ! [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
    echo "polap local version update is canceled."
    return
  fi

  echo "Updating the Polap Conda environment to the local source at ${local_path}."

  # Check if conda environment exists
  if [[ ! -d "$conda_base/envs/$conda_env" ]]; then
    echo "[ERROR] Conda environment '$conda_env' not found." >&2
    return 1
  fi

  # Create a temporary directory
  echo "[INFO] Checking the local repository: $local_path"
  if [[ ! -d "${local_path}" ]]; then
    echo "[ERROR] no such local path: ${local_path}"
    return 1
  fi
  if [[ ! -d "${local_path}/src" ]]; then
    echo "[ERROR] no src in the local path: ${local_path}"
    return 1
  fi

  cd $local_path/src
  bash polaplib/polap-build.sh >../build.sh
  cd ..
  PREFIX="$(conda info --base)/envs/polap" bash build.sh
  # echo "Error: You do not have polap environment."

  echo "[INFO] Updated polap in conda env '$conda_env' from the local version at ${local_path}"
}

update-github_genus_species() {
  local ref="${1:-main}" # tag or commit hash (default: main)
  local repo_url="https://github.com/goshng/polap.git"
  local conda_env="polap"
  local conda_base
  conda_base=$(conda info --base)
  local target_path="$conda_base/envs/$conda_env/bin/polap"

  if [[ "${opt_y_flag}" == false ]]; then
    read -p "Do you want to update conda env '$conda_env' from GitHub ($ref)? (y/N): " confirm
  else
    confirm="yes"
  fi

  if ! [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
    echo "polap update canceled."
    echo "${help_message_bleeding_edge_polap}"
    return
  fi

  echo "[INFO] Updating polap to ref '$ref' from GitHub"

  if [[ ! -d "$conda_base/envs/$conda_env" ]]; then
    echo "[ERROR] Conda environment '$conda_env' not found." >&2
    return 1
  fi

  local temp_dir
  temp_dir=$(mktemp -d)
  echo "[INFO] Cloning repository into $temp_dir"

  git clone "$repo_url" "$temp_dir" || {
    echo "[ERROR] Failed to clone repository." >&2
    rm -rf "$temp_dir"
    return 1
  }

  cd "$temp_dir"
  git checkout "$ref" || {
    echo "[ERROR] Failed to checkout ref '$ref'" >&2
    cd /
    rm -rf "$temp_dir"
    return 1
  }

  cd src
  bash polaplib/polap-build.sh >../build.sh
  cd ..
  PREFIX="$conda_base/envs/$conda_env" bash build.sh

  echo "[INFO] Updated polap in conda env '$conda_env' to ref '$ref'"

  cd /
  rm -rf "$temp_dir"
}

update_genus_species() {
  local args=("$@")

  local item="${args[0]}"
  local rest=("${args[@]:1}")
  update-${item}_genus_species ${rest[@]}
}

get-cflye_genus_species() {
  local _brg_outdir="${1}"
  local _brg_host="${2}"

  if [[ "${_brg_outdir}" == "-h" || "${_brg_outdir}" == "--help" ]]; then
    echo "$help_message_get_cflye"
    return
  fi

  if [[ -z "${_brg_outdir}" ]]; then
    echo "[ERROR] outdir is required."
    return
  fi

  if [[ ! -d "${_brg_outdir}" ]]; then
    echo "[ERROR] no such folder: ${_brg_outdir}"
    return
  fi

  if [[ -z "${_brg_host}" ]]; then
    echo "[ERROR] hostname is required".
    return
  fi

  # rsync from the remote
  echo Getting cflye1 results of ${_brg_outdir} from ${_brg_host} ...
  rsync -azuq --max-size=5M \
    "${_brg_host}:${PWD}/${_brg_outdir}"/ \
    "${_brg_outdir}"/
  # copy and compress cns.fa
  rsync -azq "${_brg_host}:${PWD}/${_brg_outdir}/t1/0/cns.fa" "${_brg_outdir}/t1/0/"
  rsync -azq "${_brg_host}:${PWD}/${_brg_outdir}/t1/0/msbwt/comp_msbwt.npy" "${_brg_outdir}/t1/0/msbwt/"
  # gzip "${_brg_outdir}/t1/0/cns.fa"
}

get_genus_species() {
  local args=("$@")

  local item="${args[0]}"
  local rest=("${args[@]:1}")
  get-${item}_genus_species ${rest[@]}
}

print-help-all_genus_species() {
  local files=("polap-data-v2.sh" "polaplib/polap-lib-data.sh")
  local var_name
  local var_suffix
  local file

  for file in "${files[@]}"; do
    file="${_polap_script_bin_dir}/$file"

    if [[ -f "$file" ]]; then
      # Source the file in a subshell to avoid polluting current shell
      (
        source "$file"
        for var_name in ${!help_message_*}; do
          var_suffix="${var_name#help_message_}"
          echo ">>> $var_suffix"
          printf "%s\n\n" "${!var_name}"
        done
      )
    else
      echo "[WARN] File not found: $file" >&2
    fi
  done
}

help_genus_species() {
  local _brg_outdir="${1}"

  local subcmd1="${_brg_outdir}"
  local _subcmd1_clean="${subcmd1//-/_}"

  local target_var="help_message_${_subcmd1_clean}"

  if declare -p "$target_var" &>/dev/null; then
    declare -n ref="$target_var"
    echo "$ref"
  else
    echo "[ERROR] No such polap help or subcommand: $subcmd1" >&2
  fi

}

setup-fmlrc2_genus_species() {
  local dest_dir="$HOME/.cargo"

  # Add to PATH if not already present
  if ! grep -q "$dest_dir/bin" ~/.bashrc; then
    echo "export PATH=\"$dest_dir/bin:\$PATH\"" >>~/.bashrc
    echo "[INFO] Added $HOME/.cargo to PATH in ~/.bashrc" >&2

    # Final message
    echo "[INFO] fmlrc2 installed. Run 'source ~/.bashrc' or restart your terminal." >&2
    echo "[INFO] Test with: msbwt2-build -h" >&2
    echo "[INFO] Test with: fmlrc2 -h" >&2
  fi
}

_install-polap-fmlrc2_conda_env() {
  local env_name="polap-fmlrc2"
  local tools=("msbwt2-build" "msbwt2-convert")
  local cargo_bin="$HOME/.cargo/bin"

  echo "[INFO] Creating or updating Conda environment: $env_name" >&2
  conda create -y -n "$env_name" -c bioconda -c conda-forge fmlrc2 rust || {
    echo "[ERROR] Failed to create environment: $env_name" >&2
    return 1
  }

  # Activate environment
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "$env_name"

  # Install msbwt2
  if ! command -v msbwt2-build &>/dev/null || ! command -v msbwt2-convert &>/dev/null; then
    echo "[INFO] Installing msbwt2 using cargo..." >&2
    cargo install msbwt2 || {
      echo "[ERROR] cargo install msbwt2 failed" >&2
      return 1
    }
  fi

  # Symlink into conda's bin directory
  local conda_bin="${CONDA_PREFIX}/bin"
  mkdir -p "$conda_bin"
  for tool in "${tools[@]}"; do
    if [[ -x "${cargo_bin}/${tool}" ]]; then
      ln -sf "${cargo_bin}/${tool}" "${conda_bin}/${tool}"
      echo "[INFO] Linked ${tool} to ${conda_bin}" >&2
    else
      echo "[WARNING] ${tool} not found in ${cargo_bin}" >&2
    fi
  done

  echo "[SUCCESS] Environment '$env_name' is ready with fmlrc2 and msbwt2" >&2
  return 0
}

install-fmlrc2_genus_species() {
  if [[ "${opt_y_flag}" == false ]]; then
    read -p "Do you want to install polap-fmlrc2? (y/N): " confirm
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
      echo "You're in the base environment. Creating 'polap-fmlrc2'..."
      if conda env list | awk '{print $1}' | grep -qx "polap-fmlrc2"; then
        echo "ERROR: Conda environment 'polap-fmlrc2' already exists."
      else
        _install-polap-fmlrc2_conda_env
        setup-fmlrc2_genus_species
      fi
    else
      echo "Error: You're in the '$CONDA_DEFAULT_ENV' environment. Please activate base before running this script."
      exit 1
    fi
  else
    echo "polap-fmlrc2 installation is canceled."
    echo "${help_message_install_fmlrc2}"
  fi
}

install-minimal_genus_species() {
  local tools_to_install=(
    polap
    fmlrc
    man
    cflye
  )

  for item in "${tools_to_install[@]}"; do
    install-${item}_genus_species
  done
}

install-all_genus_species() {
  local _version="${1}"

  local tools_to_install=(
    fmlrc
    fmlrc2
    getorganelle
    pmat
    tippo
    oatk
    man
    cflye
    dflye
  )

  if [[ "${_brg_outdir}" == "-h" || "${_brg_outdir}" == "--help" ]]; then
    echo "$help_message_install_all"
    return
  fi

  install-polap_genus_species "${_version}"
  for item in "${tools_to_install[@]}"; do
    install-${item}_genus_species
  done
}

setup-nvim_genus_species() {
  # required
  mv ~/.config/nvim{,.bak} >/dev/null 2>&1

  # optional but recommended
  mv ~/.local/share/nvim{,.bak} >/dev/null 2>&1
  mv ~/.local/state/nvim{,.bak} >/dev/null 2>&1
  mv ~/.cache/nvim{,.bak} >/dev/null 2>&1

  # Clone the starter
  git clone https://github.com/LazyVim/starter ~/.config/nvim

  # Remove the .git folder, so you can add it to your own repo later
  rm -rf ~/.config/nvim/.git

  # Start Neovim!
  echo nvim

  if ! grep -q "alias v=" ~/.bashrc; then
    echo "alias v='nvim'" >>~/.bashrc
    echo "export EDITOR='nvim'" >>~/.bashrc
    echo "[INFO] Added alias v for nvim in ~/.bashrc" >&2
  fi
}

finalize-pmat_genus_species() {
  local url_basename="PMAT-1.5.3"
  local url="https://github.com/bichangwei/PMAT/archive/refs/tags/v1.5.3.tar.gz"
  local dest_dir="$HOME/bin/pmat"
  local tmp_dir

  # Create a temporary working directory
  tmp_dir=$(mktemp -d)
  echo "[INFO] Downloading PMAT to $tmp_dir..." >&2

  # Download and extract
  # Download using wget
  wget -q --show-progress -O "$tmp_dir/pmat.tar.gz" "$url" || {
    echo "[ERROR] Failed to download PMAT binary." >&2
    return 1
  }

  # Download using curl
  # curl -L "$url" -o "$tmp_dir/pmat.tar.gz" || {
  # 	echo "[ERROR] Failed to download PMAT binary." >&2
  # 	return 1
  # }

  tar -xzf "$tmp_dir/pmat.tar.gz" -C "$tmp_dir" || {
    echo "[ERROR] Failed to extract PMAT archive." >&2
    return 1
  }

  # Move to /opt (requires sudo)
  echo "[INFO] Installing to $dest_dir..." >&2
  rm -rf "$dest_dir"
  if [[ ! -d "$HOME/bin" ]]; then

    if [[ "${opt_y_flag}" == false ]]; then
      read -p "Proceed to create a folder at $HOME/bin? (y/N): " confirm
    else
      confirm="yes"
    fi

    if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
      mkdir -p "$HOME/bin"
    else
      echo "[ERROR] Failed to create $HOME/bin"
      return 1
    fi
  fi
  mv "$tmp_dir/${url_basename}" "$dest_dir"

  # Add to PATH if not already present
  if ! grep -q "$dest_dir/bin" ~/.bashrc; then
    echo "export PATH=\"$dest_dir/bin:\$PATH\"" >>~/.bashrc
    echo "[INFO] Added PMAT to PATH in ~/.bashrc" >&2
  fi

  # Clean up
  rm -rf "$tmp_dir"

  # Final message
  echo "[INFO] PMAT installed. Run 'source ~/.bashrc' or restart your terminal." >&2
  echo "[INFO] Test with: PMAT --version" >&2
}

setup-pmat_genus_species() {
  echo "Install fuse2fs gocryptfs ..."
  sudo apt update
  sudo apt -y install fuse2fs gocryptfs

  echo "Setup apptainer ..."
  _polap_lib_conda-ensure_conda_env polap-pmat || exit 1
  apptainer exec $HOME/bin/pmat/container/runAssembly.sif echo "Fail: PMAT setup"
  echo 'kernel.unprivileged_userns_clone=1' | sudo tee /etc/sysctl.d/90-userns.conf
  echo 'kernel.apparmor_restrict_unprivileged_userns=0' | sudo tee /etc/sysctl.d/80-apparmor-userns.conf
  sudo sysctl --system
  apptainer exec $HOME/bin/pmat/container/runAssembly.sif echo "Success: PMAT setup"
  conda deactivate
}

setup-polap_genus_species() {
  local _brg_outdir="${1:-na}"
  local bashrc="$HOME/.bashrc"
  # local bashrc="$PWD/bashrc.txt"
  local start_marker="# >>> polap initialize >>>"
  local end_marker="# <<< polap initialize <<<"

  # Customize this path if you move the scripts
  local TEMPLATE_DIR="polap/src"

  if [[ "$_brg_outdir" == "--remove" ]]; then
    if grep -q "$start_marker" "$bashrc"; then
      sed -i "/$start_marker/,/$end_marker/d" "$bashrc"
      echo "[INFO] Removed polap initialize block from $bashrc" >&2
    else
      echo "[INFO] No polap initialize block found in $bashrc" >&2
    fi
    return 0
  fi

  local block_content
  block_content=$(
    cat <<EOF
$start_marker
export PATH="\$HOME/bin/pmat/bin:\$PATH"
export PATH="\$HOME/.cargo/bin:\$PATH"
alias p1='bash $TEMPLATE_DIR/polap-data-v1.sh'
alias p2='bash $TEMPLATE_DIR/polap-data-v2.sh'
alias pl='bash $TEMPLATE_DIR/polap.sh'
alias p='bash $TEMPLATE_DIR/polap-data-v2.sh'
alias p4='bash $TEMPLATE_DIR/polap-data-v4.sh'
$end_marker
EOF
  )

  sed -i "/$start_marker/,/$end_marker/d" "$bashrc"
  echo "$block_content" >>"$bashrc"
  echo "[INFO] Added polap initialize block to $bashrc using template dir: $TEMPLATE_DIR" >&2
}

setup-bandage_genus_species() {
  local bashrc="$HOME/.bashrc"
  # local bashrc="$PWD/bashrc.txt"
  local start_marker="# >>> polap-bandage initialize >>>"
  local end_marker="# <<< polap-bandage initialize <<<"

  # Customize this path if you move the scripts
  local TEMPLATE_DIR="polap/src"

  if [[ "$1" == "--remove" ]]; then
    if grep -q "$start_marker" "$bashrc"; then
      sed -i "/$start_marker/,/$end_marker/d" "$bashrc"
      echo "[INFO] Removed polap-bandage initialize block from $bashrc" >&2
    else
      echo "[INFO] No polap-bandage initialize block found in $bashrc" >&2
    fi
    return 0
  fi

  local block_content
  block_content=$(
    cat <<EOF
$start_marker
export PATH="\$HOME/bin:\$PATH"
$end_marker
EOF
  )

  sed -i "/$start_marker/,/$end_marker/d" "$bashrc"
  echo "$block_content" >>"$bashrc"
  echo "[INFO] Added polap-bandage initialize block to $bashrc using template dir: $TEMPLATE_DIR" >&2
}

install-nvim_genus_species() {
  local url_basename="nvim-linux-x86_64"
  local url="https://github.com/neovim/neovim/releases/latest/download/${url_basename}.tar.gz"
  local dest_dir="/opt/nvim"
  local tmp_dir

  # Create a temporary working directory
  tmp_dir=$(mktemp -d)
  echo "[INFO] Downloading Neovim to $tmp_dir..." >&2

  # Download and extract
  wget -q --show-progress -O "$tmp_dir/nvim.tar.gz" "$url" || {
    echo "[ERROR] Failed to download Neovim binary." >&2
    return 1
  }
  # curl -L "$url" -o "$tmp_dir/nvim.tar.gz" || {
  # 	echo "[ERROR] Failed to download Neovim binary." >&2
  # 	return 1
  # }

  tar -xzf "$tmp_dir/nvim.tar.gz" -C "$tmp_dir" || {
    echo "[ERROR] Failed to extract Neovim archive." >&2
    return 1
  }

  # Move to /opt (requires sudo)
  echo "[INFO] Installing to $dest_dir..." >&2
  sudo rm -rf "$dest_dir"
  sudo mv "$tmp_dir/${url_basename}" "$dest_dir"

  # Add to PATH if not already present
  if ! grep -q "$dest_dir/bin" ~/.bashrc; then
    echo "export PATH=\"$dest_dir/bin:\$PATH\"" >>~/.bashrc
    echo "[INFO] Added Neovim to PATH in ~/.bashrc" >&2
  fi

  # Clean up
  rm -rf "$tmp_dir"

  # Final message
  echo "[INFO] Neovim installed. Run 'source ~/.bashrc' or restart your terminal." >&2
  echo "[INFO] Test with: nvim --version" >&2
}

uninstall_genus_species() {
  local args=("$@")

  if has_help "${args[@]}"; then
    echo "$help_message_uninstall"
    return
  fi

  for item in "$@"; do
    uninstall-${item}_genus_species
  done
}

uninstall-polap_genus_species() {
  echo "Removing the following conda environments: polap"
  if [[ "${opt_y_flag}" == false ]]; then
    read -p "Do you want to uninstall polap conda environments? (y/N): " confirm
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
        conda env remove -y -n polap
        # conda remove -y -n polap --all
      fi
      echo "conda deactivate if necessary"
    fi
  else
    echo "Uninstallation of polap is canceled."
    echo "${help_message_uninstall}"
  fi
}

uninstall-getorganelle_genus_species() {
  echo "Removing the following conda environments: polap-getorganelle"
  if [[ "${opt_y_flag}" == false ]]; then
    read -p "Do you want to uninstall polap-getorganelle conda environments for polap? (y/N): " confirm
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
      if conda info --envs | awk '{print $1}' | grep -Fxq polap-getorganelle; then
        conda env remove -y -n polap-getorganelle
        # conda remove -y -n polap-getorganelle --all
      fi
      echo "conda deactivate if necessary"
    fi
  else
    echo "Uninstallation of getorganelle is canceled."
    echo "${help_message_uninstall_getorganelle}"
  fi
}

uninstall-fmlrc2_genus_species() {
  echo "Removing the following conda environments: polap-fmlrc2"
  if [[ "${opt_y_flag}" == false ]]; then
    read -p "Do you want to uninstall polap-fmlrc2 conda environments for polap? (y/N): " confirm
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
      if conda info --envs | awk '{print $1}' | grep -Fxq polap-fmlrc2; then
        conda env remove -y -n polap-fmlrc2
        # conda remove -y -n polap-fmlrc2 --all
      fi
      echo "conda deactivate if necessary"
    fi
  else
    echo "Uninstallation of fmlrc2 is canceled."
    echo "${help_message_uninstall_fmlrc2}"
  fi
}

uninstall-fmlrc_genus_species() {
  echo "Removing the following conda environments: polap-fmlrc"
  if [[ "${opt_y_flag}" == false ]]; then
    read -p "Do you want to uninstall polap-fmlrc conda environments for polap? (y/N): " confirm
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
      if conda info --envs | awk '{print $1}' | grep -Fxq polap-fmlrc; then
        conda env remove -y -n polap-fmlrc
        # conda remove -y -n polap-fmlrc --all
      fi
      echo "conda deactivate if necessary"
    fi
  else
    echo "Uninstallation of fmlrc is canceled."
    echo "${help_message_uninstall_fmlrc}"
  fi
}

uninstall-dflye_genus_species() {
  echo "Removing the following conda environments: polap-dflye"
  if [[ "${opt_y_flag}" == false ]]; then
    read -p "Do you want to uninstall polap-dflye conda environments for polap? (y/N): " confirm
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
      if conda info --envs | awk '{print $1}' | grep -Fxq polap-dflye; then
        conda env remove -y -n polap-dflye
        # conda remove -y -n polap-dflye --all
      fi
      echo "conda deactivate if necessary"
    fi
  else
    echo "Uninstallation of dflye is canceled."
    echo "${help_message_uninstall_dflye}"
  fi
}

uninstall-cflye_genus_species() {
  echo "Removing the following conda environments: polap-cflye"
  if [[ "${opt_y_flag}" == false ]]; then
    read -p "Do you want to uninstall polap-cflye conda environments for polap? (y/N): " confirm
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
      if conda info --envs | awk '{print $1}' | grep -Fxq polap-cflye; then
        conda env remove -y -n polap-cflye
        # conda remove -y -n polap-cflye --all
      fi
      echo "conda deactivate if necessary"
    fi
  else
    echo "Uninstallation of cflye is canceled."
    echo "${help_message_uninstall_cflye}"
  fi
}

uninstall-oatk_genus_species() {
  echo "Removing the following conda environments: polap-oatk"
  if [[ "${opt_y_flag}" == false ]]; then
    read -p "Do you want to uninstall polap-oatk conda environments for polap? (y/N): " confirm
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
      if conda info --envs | awk '{print $1}' | grep -Fxq polap-oatk; then
        conda env remove -y -n polap-oatk
        # conda remove -y -n polap-oatk --all
      fi
      echo "conda deactivate if necessary"
    fi
  else
    echo "Uninstallation of oatk is canceled."
    echo "${help_message_uninstall_oatk}"
  fi
}

uninstall-tippo_genus_species() {
  echo "Removing the following conda environments: polap-tippo"
  if [[ "${opt_y_flag}" == false ]]; then
    read -p "Do you want to uninstall polap-tippo conda environments for polap? (y/N): " confirm
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
      if conda info --envs | awk '{print $1}' | grep -Fxq polap-tippo; then
        conda env remove -y -n polap-tippo
        # conda remove -y -n polap-tippo --all
      fi
      echo "conda deactivate if necessary"
    fi
  else
    echo "Uninstallation of tippo is canceled."
    echo "${help_message_uninstall_tippo}"
  fi
}

uninstall-pmat_genus_species() {
  echo "Removing the following conda environments: polap-pmat"
  if [[ "${opt_y_flag}" == false ]]; then
    read -p "Do you want to uninstall polap-pmat conda environments for polap? (y/N): " confirm
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
      if conda info --envs | awk '{print $1}' | grep -Fxq polap-pmat; then
        conda env remove -y -n polap-pmat
        rm -rf $HOME/pmat
        # conda remove -y -n polap-pmat --all
      fi
      echo "conda deactivate if necessary"
    fi
  else
    echo "Uninstallation of pmat is canceled."
    echo "${help_message_uninstall_pmat}"
  fi
}

uninstall-efg_genus_species() {
  echo "Removing the following conda environments: polap-efg"
  if [[ "${opt_y_flag}" == false ]]; then
    read -p "Do you want to uninstall polap-efg conda environments for polap? (y/N): " confirm
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
      if conda info --envs | awk '{print $1}' | grep -Fxq polap-efg; then
        conda env remove -y -n polap-efg
        # conda remove -y -n polap-efg --all
      fi
      echo "conda deactivate if necessary"
    fi
  else
    echo "Uninstallation of efg is canceled."
    echo "${help_message_uninstall_efg}"
  fi
}

run-summary-data_genus_species() {
  local _brg_outdir="${1}"
  local _brg_inum="${2:-0}"

  # Set the run title
  local full_name="${FUNCNAME[0]}"
  local middle_part="${full_name#run-}"
  local _run_title="${middle_part%%_*}"

  # Directory without a slash
  _brg_outdir="${_brg_outdir%/}"

  # Folders
  local _brg_outdir_t="${_brg_outdir}/${opt_t_arg}"
  local _brg_outdir_i="${_brg_outdir_t}/${_brg_inum}"
  local _brg_rundir="${_brg_outdir_i}/${_run_title}"
  local _run_dir="${_brg_outdir}-${_run_title}"
  local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

  mkdir -p "${_run_dir}"
  mkdir -p "${_brg_rundir}"

  local target_index="${_brg_outdir}-${_brg_inum}"
  local long_sra="${_long["$target_index"]}"
  if [[ -z "$long_sra" ]]; then
    echo "Error: skipping ${_brg_outdir}-${_brg_inum} because it is not in the CSV."
    return
  fi
  local short_sra="${_short["$target_index"]}"

  # Activate a conda environment
  _polap_lib_conda-ensure_conda_env polap || exit 1

  # Files for memtracker
  local _stdout_txt="${_brg_outdir_i}/stdout-${_run_title}.txt"
  local _timing_txt="${_brg_outdir_i}/timing-${_run_title}.txt"
  local _memlog_file="${_brg_outdir_i}/memlog-${_run_title}.csv"
  local _summary_file="${_brg_outdir_i}/summary-${_run_title}.txt"

  rm -f "${_stdout_txt}"
  rm -f "${_timing_txt}"

  # Start memory logger
  _polap_lib_process-start_memtracker "${_memlog_file}" \
    "${_polap_var_memtracker_time_interval}"

  command time -v seqkit stats -T \
    ${long_sra}.fastq \
    -o "${_run_dir}"/l.fq.stats \
    >>"${_stdout_txt}" \
    2>>"${_timing_txt}"

  command time -v seqkit stats -T \
    ${short_sra}_1.fastq \
    -o "${_run_dir}"/s1.fq.stats \
    >>"${_stdout_txt}" \
    2>>"${_timing_txt}"

  command time -v seqkit stats -T \
    ${short_sra}_2.fastq \
    -o "${_run_dir}"/s2.fq.stats \
    >>"${_stdout_txt}" \
    2>>"${_timing_txt}"

  # Save system info
  _polap_lib_timing-get_system_info >>"${_timing_txt}"

  _polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"

  conda deactivate

  # Save results
  rsync -azuq "${_run_dir}"/ "${_brg_rundir}"/

  # Clean-up
  rm -rf "${_run_dir}"
}

install-man_genus_species() {
  if [[ "${opt_y_flag}" == false ]]; then
    read -p "Do you want to install man in the polap-man conda environment? (y/N): " confirm
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
      echo "You're in the base environment. Creating 'polap-man'..."
      if conda env list | awk '{print $1}' | grep -qx "polap-man"; then
        echo "ERROR: Conda environment 'polap-man' already exists."
      else
        conda create -y --name polap-man pandoc pandoc-crossref yq jq
      fi
    else
      echo "Error: You're in the '$CONDA_DEFAULT_ENV' environment. Please activate base before running this script."
      exit 1
    fi
  else
    echo "man installation is canceled."
    echo "Ref: https://pandoc.org"
    echo "${help_message_install_man}"
  fi
}

# Case of the check menu
# --disassemble-c
# --disassemble-align-reference
# --disassemble-simple-polishing
#
# check <outdir> <inum> [simple|polish]
#
run-polap-disassemble-check_genus_species() {
  local _brg_outdir="$1"
  local _brg_inum="${2:-0}"
  local simple_polishing="${3:-default}"
  local _brg_random="${4:-off}"

  # Set the run title
  local full_name="${FUNCNAME[0]}"
  local middle_part="${full_name#run-}"
  local _run_title="${middle_part%%-check*}"
  # local _run_title="${middle_part%%_*}-${simple_polishing}"

  # Directory without a slash
  _brg_outdir="${_brg_outdir%/}"

  # Folders
  local _brg_outdir_t="${_brg_outdir}/${opt_t_arg}"
  local _brg_outdir_i="${_brg_outdir_t}/${_brg_inum}"
  local _brg_rundir="${_brg_outdir_i}/${_run_title}"
  local _run_dir="${_brg_outdir}-${_run_title}"
  local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

  mkdir -p ${_brg_rundir}
  if [[ ! -d "${_run_dir}" ]]; then
    echo "Error: run polap-disassemble before running this command."
    return
  fi

  local target_index="${_brg_outdir}-${_brg_inum}"

  local species_name="$(echo ${_brg_outdir} | sed 's/_/ /')"
  local long_sra="${_long["$target_index"]}"
  local short_sra="${_short["$target_index"]}"
  local random_seed="${_random_seed["$target_index"]}"

  if [ -z "${long_sra}" ]; then
    echo "ERROR: no long-read SRA ID: ${long_sra}"
    echo "Suggestion: use -c option for a user-povided CSV."
    return
  fi

  if [[ "${_brg_random}" == "on" ]]; then
    random_seed=0
  fi

  local i=0
  local n
  local p
  local extracted_n="${_compare_n["$target_index"]}"
  local extracted_p="${_compare_p["$target_index"]}"
  local extracted_r="${_compare_r["$target_index"]}"
  local extracted_memory="${_memory["$target_index"]}"
  local extracted_downsample="${_downsample["$target_index"]}"
  local extracted_alpha="${_disassemble_alpha["$target_index"]}"
  local extracted_delta="${_disassemble_delta["$target_index"]}"

  mkdir -p "${_brg_outdir_i}"

  # Files for memtracker
  local _stdout_txt="${_brg_outdir_i}/stdout-${_run_title}-check-${simple_polishing}.txt"
  local _timing_txt="${_brg_outdir_i}/timing-${_run_title}-check-${simple_polishing}.txt"
  local _memlog_file="${_brg_outdir_i}/memlog-${_run_title}-check-${simple_polishing}.csv"
  local _summary_file="${_brg_outdir_i}/summary-${_run_title}-check-${simple_polishing}.txt"

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

  # Start memory logger
  _polap_lib_process-start_memtracker "${_memlog_file}" \
    "${_polap_var_memtracker_time_interval}"

  local i=0
  local n="${extracted_n}"
  local p="${extracted_p}"

  i=$((i + 1))

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

  _log_echo "Analysis (check-${_s_i}): ${_brg_outdir} at ${_brg_inum}"
  _log_echo "($i) n=$n, p=$p, memory=${extracted_memory}G, downsample=${extracted_downsample}x"

  # "${_stages}" \
  command time -v ${_polap_cmd} disassemble \
    ${_stages} \
    --downsample ${extracted_downsample} \
    -i ${_brg_inum} \
    -o ${_run_dir} \
    -l ${long_sra}.fastq \
    -a ${short_sra}_1.fastq \
    -b ${short_sra}_2.fastq \
    --disassemble-c ${_brg_outdir_i}/ptdna-ptgaul.fa \
    --disassemble-align-reference \
    ${simple_polishing} \
    --disassemble-i "${_d_i}" \
    --disassemble-n $n \
    --disassemble-p $p \
    --disassemble-r ${extracted_r} \
    --disassemble-memory ${extracted_memory} \
    --random-seed "${random_seed}" \
    >"${_stdout_txt}" \
    2>"${_timing_txt}"

  # Save results separately from 3-check
  if [[ -d "${_run_dir}/${_brg_inum}/disassemble/${_d_i}/3" ]]; then
    rm -rf "${_run_dir}/${_brg_inum}/disassemble/${_d_i}/3-check"
    mv "${_run_dir}/${_brg_inum}/disassemble/${_d_i}/3" \
      "${_run_dir}/${_brg_inum}/disassemble/${_d_i}/3-check"
  fi

  # Compare the results using mauve, blast, and mafft
  local mauve_dir="${_run_dir}/mauve/${i}"
  local blast_dir="${_run_dir}/blast/${i}"
  local mafft_dir="${_run_dir}/mafft/${i}"
  mkdir -p "${mauve_dir}"
  mkdir -p "${blast_dir}"
  mkdir -p "${mafft_dir}"
  if [[ -s "${_run_dir}/${_brg_inum}/disassemble/${_d_i}/pt.subsample-polishing.reference.aligned.1.fa" ]]; then

    ${_polap_cmd} mauve-mtdna -a "${_brg_outdir_i}/ptdna-ptgaul.fa" \
      -b "${_run_dir}/${_brg_inum}/disassemble/${_d_i}/pt.subsample-polishing.reference.aligned.1.fa" \
      -o "${mauve_dir}" \
      >"${mauve_dir}/log.txt"
    # echo "see ${mauve_dir}/log.txt"
    # cat "${mauve_dir}/log.txt"

    ${_polap_cmd} compare2ptdna -a "${_brg_outdir_i}/ptdna-ptgaul.fa" \
      -b "${_run_dir}/${_brg_inum}/disassemble/${_d_i}/pt.subsample-polishing.reference.aligned.1.fa" \
      -o "${blast_dir}"
    # echo "see ${blast_dir}/pident.txt"
    # cat "${blast_dir}/pident.txt"

    # mafft/1. ptGAUL vs. subsample-polishing
    i=1
    ${_polap_cmd} mafft-mtdna \
      -a "${_brg_outdir_i}/ptdna-ptgaul.fa" \
      -b "${_run_dir}/${_brg_inum}/disassemble/${_d_i}/pt.subsample-polishing.reference.aligned.1.fa" \
      -o "${mafft_dir}" \
      >"${mafft_dir}/log.txt"
    # echo "see ${mafft_dir}/pident.txt"
    # cat "${mafft_dir}/pident.txt"

  else
    echo "ERROR: no such file: ${_brg_outdir_i}/disassemble/${_d_i}/pt.subsample-polishing.reference.aligned.1.fa"
  fi

  # End with summary of the system usage
  _polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"

  # Save system info
  _polap_lib_timing-get_system_info >>"${_timing_txt}"

  conda deactivate

  # Save results
  rsync -azuq --max-size=5M "${_run_dir}"/ "${_brg_rundir}"/
  # rm -rf "${_run_dir}"
}

run-polap-disassemble-compare_genus_species() {
  local _brg_outdir="$1"
  local _brg_inum="${2:-0}"
  local simple_polishing="${3:-default}"
  local _brg_random="${4:-off}"

  # Set the run title
  local full_name="${FUNCNAME[0]}"
  local middle_part="${full_name#run-}"
  local _run_title="${middle_part%%-compare*}"
  # local _run_title="${middle_part%%_*}-${simple_polishing}"

  # Directory without a slash
  _brg_outdir="${_brg_outdir%/}"

  # Folders
  local _brg_outdir_t="${_brg_outdir}/${opt_t_arg}"
  local _brg_outdir_i="${_brg_outdir_t}/${_brg_inum}"
  local _brg_rundir="${_brg_outdir_i}/${_run_title}"
  local _run_dir="${_brg_outdir}-${_run_title}"
  local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

  mkdir -p ${_brg_rundir}
  if [[ ! -d "${_run_dir}" ]]; then
    echo "Error: run polap-disassemble before running this command."
    return
  fi

  local target_index="${_brg_outdir}-${_brg_inum}"

  local species_name="$(echo ${_brg_outdir} | sed 's/_/ /')"
  local long_sra="${_long["$target_index"]}"
  local short_sra="${_short["$target_index"]}"
  local random_seed="${_random_seed["$target_index"]}"

  if [ -z "${long_sra}" ]; then
    echo "ERROR: no long-read SRA ID: ${long_sra}"
    echo "Suggestion: use -c option for a user-povided CSV."
    return
  fi

  if [[ "${_brg_random}" == "on" ]]; then
    random_seed=0
  fi

  local i=0
  local n
  local p
  local extracted_n="${_compare_n["$target_index"]}"
  local extracted_p="${_compare_p["$target_index"]}"
  local extracted_r="${_compare_r["$target_index"]}"
  local extracted_memory="${_memory["$target_index"]}"
  local extracted_downsample="${_downsample["$target_index"]}"
  local extracted_alpha="${_disassemble_alpha["$target_index"]}"
  local extracted_delta="${_disassemble_delta["$target_index"]}"

  mkdir -p "${_brg_outdir_i}"

  _polap_lib_conda-ensure_conda_env polap || exit 1

  # Files for memtracker
  local _stdout_txt="${_brg_outdir_i}/stdout-${_run_title}-check-${simple_polishing}.txt"
  local _timing_txt="${_brg_outdir_i}/timing-${_run_title}-check-${simple_polishing}.txt"
  local _memlog_file="${_brg_outdir_i}/memlog-${_run_title}-check-${simple_polishing}.csv"
  local _summary_file="${_brg_outdir_i}/summary-${_run_title}-check-${simple_polishing}.txt"

  local i=0
  local n="${extracted_n}"
  local p="${extracted_p}"

  i=$((i + 1))
  _log_echo "Analysis (compare): ${_brg_outdir} at ${_brg_inum}"
  _log_echo "($i) n=$n, p=$p, memory=${extracted_memory}G, downsample=${extracted_downsample}x"

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

  # Compare the results using mauve, blast, and mafft
  local mafft_dir="${_run_dir}/mafft/${i}"
  mkdir -p "${mafft_dir}"
  if [[ -s "${_run_dir}/${_brg_inum}/disassemble/${_d_i}/pt.subsample-polishing.reference.aligned.1.fa" ]]; then

    # mafft/1. ptGAUL vs. subsample-polishing
    i=1
    ${_polap_cmd} mafft-mtdna \
      -a "${_brg_outdir_i}/ptdna-ptgaul.fa" \
      -b "${_run_dir}/${_brg_inum}/disassemble/${_d_i}/pt.subsample-polishing.reference.aligned.1.fa" \
      -o "${mafft_dir}" \
      >"${mafft_dir}/log.txt"
    # echo "see ${mafft_dir}/pident.txt"
    # cat "${mafft_dir}/pident.txt"

    # mafft/2. simple-polishing vs. subsample-polishing
    i=2
    local mafft_dir="${_run_dir}/mafft/${i}"
    mkdir -p "${mafft_dir}"
    ${_polap_cmd} mafft-mtdna \
      -a "${_run_dir}/${_brg_inum}/disassemble/${_d_i}/pt.simple-polishing.reference.aligned.1.fa" \
      -b "${_run_dir}/${_brg_inum}/disassemble/${_d_i}/pt.subsample-polishing.reference.aligned.1.fa" \
      -o "${mafft_dir}" \
      >"${mafft_dir}/log.txt"

    # mafft/3. ptGAUL vs. simple-polishing
    i=3
    local mafft_dir="${_run_dir}/mafft/${i}"
    mkdir -p "${mafft_dir}"
    ${_polap_cmd} mafft-mtdna \
      -a "${_brg_outdir_i}/ptdna-ptgaul.fa" \
      -b "${_run_dir}/${_brg_inum}/disassemble/${_d_i}/pt.simple-polishing.reference.aligned.1.fa" \
      -o "${mafft_dir}" \
      >"${mafft_dir}/log.txt"

  else
    echo "ERROR: no such file: ${_brg_outdir_i}/disassemble/${_d_i}/pt.subsample-polishing.reference.aligned.1.fa"
  fi

  conda deactivate

  # Save results
  rsync -azuq --max-size=5M "${_run_dir}"/ "${_brg_rundir}"/
  echo rm -rf "${_run_dir}"
}

# Case of the infer menu
# no --disassemble-c
#
# infer <outdir> <inum> [default|polishing|simple]
#
# default: all stages
# polish: stage 3 only
# simple: stage 3 only but simple polishing
run-polap-disassemble_genus_species() {
  local _brg_outdir="$1"
  local _brg_inum="${2:-0}"
  local simple_polishing="${3:-default}"
  local _brg_random="${4:-off}"

  if [[ "${_brg_outdir}" == "-h" || "${_brg_outdir}" == "--help" ]]; then
    echo "$help_message_run_polap_disassemble"
    return
  fi

  # Set the run title
  local full_name="${FUNCNAME[0]}"
  local middle_part="${full_name#run-}"
  local _run_title="${middle_part%%_*}"
  # local _run_title="${middle_part%%_*}-${simple_polishing}"

  # Directory without a slash
  _brg_outdir="${_brg_outdir%/}"

  # Folders
  local _brg_outdir_t="${_brg_outdir}/${opt_t_arg}"
  local _brg_outdir_i="${_brg_outdir_t}/${_brg_inum}"
  local _brg_rundir="${_brg_outdir_i}/${_run_title}"
  local _run_dir="${_brg_outdir}-${_run_title}"
  local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

  if [[ "${simple_polishing}" == "default" ]]; then
    rm -rf ${_run_dir}
  elif [[ ! -d "${_run_dir}" ]]; then
    echo "Error: run polap-disassemble default before running this command."
    return
  fi
  mkdir -p ${_brg_rundir}
  mkdir -p ${_run_dir}

  local target_index="${_brg_outdir}-${_brg_inum}"

  local species_name="$(echo ${_brg_outdir} | sed 's/_/ /')"
  local long_sra="${_long["$target_index"]}"
  local short_sra="${_short["$target_index"]}"
  local random_seed="${_random_seed["$target_index"]}"

  if [ -z "${long_sra}" ]; then
    echo "ERROR: no long-read SRA ID: ${long_sra}"
    echo "Suggestion: use -c option for a user-povided CSV."
    return
  fi

  if [[ "${_brg_random}" == "on" ]]; then
    random_seed=0
  fi

  local i=0
  local n
  local p
  local extracted_n="${_compare_n["$target_index"]}"
  local extracted_p="${_compare_p["$target_index"]}"
  local extracted_r="${_compare_r["$target_index"]}"
  local extracted_memory="${_memory["$target_index"]}"
  local extracted_downsample="${_downsample["$target_index"]}"
  local extracted_alpha="${_disassemble_alpha["$target_index"]}"
  local extracted_delta="${_disassemble_delta["$target_index"]}"

  mkdir -p "${_brg_outdir_i}"

  # Files for memtracker
  local _stdout_txt="${_brg_outdir_i}/stdout-${_run_title}-${simple_polishing}.txt"
  local _timing_txt="${_brg_outdir_i}/timing-${_run_title}-${simple_polishing}.txt"
  local _memlog_file="${_brg_outdir_i}/memlog-${_run_title}-${simple_polishing}.csv"
  local _summary_file="${_brg_outdir_i}/summary-${_run_title}-${simple_polishing}.txt"

  _polap_lib_conda-ensure_conda_env polap || exit 1

  # Start memory logger
  _polap_lib_process-start_memtracker "${_memlog_file}" \
    "${_polap_var_memtracker_time_interval}"

  local i=0
  local n="${extracted_n}"
  local p="${extracted_p}"

  i=$((i + 1))

  #
  # default: all stages
  # polish: stage 3 only
  # simple: stage 3 only but simple polishing
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
  else
    simple_polishing=""
  fi

  _log_echo "Analysis (inference-${_s_i}): ${_brg_outdir} at ${_brg_inum}"
  _log_echo "($i) n=$n, p=$p, r=${extracted_r} memory=${extracted_memory}G, downsample=${extracted_downsample}x"

  # infer <outdir> <inum> [default|polishing|simple]

  # NOTE: "${_stages}" is a bug.
  # use it without quotations.
  command time -v ${_polap_cmd} disassemble \
    ${_stages} \
    --downsample ${extracted_downsample} \
    -i ${_brg_inum} \
    -o ${_run_dir} \
    -l ${long_sra}.fastq \
    -a ${short_sra}_1.fastq \
    -b ${short_sra}_2.fastq \
    ${simple_polishing} \
    --disassemble-i "${_d_i}" \
    --disassemble-n $n \
    --disassemble-p $p \
    --disassemble-r ${extracted_r} \
    --disassemble-memory ${extracted_memory} \
    --disassemble-alpha ${extracted_alpha} \
    --disassemble-delta ${extracted_delta} \
    --random-seed "${random_seed}" \
    >"${_stdout_txt}" \
    2>"${_timing_txt}"

  # Save results separately from 3-check
  if [[ -d "${_run_dir}/${_brg_inum}/disassemble/${_d_i}/3" ]]; then
    rm -rf "${_run_dir}/${_brg_inum}/disassemble/${_d_i}/3-infer"
    mv "${_run_dir}/${_brg_inum}/disassemble/${_d_i}/3" \
      "${_run_dir}/${_brg_inum}/disassemble/${_d_i}/3-infer"
  fi

  # End with summary of the system usage
  _polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"

  # Save system info
  _polap_lib_timing-get_system_info >>"${_timing_txt}"

  conda deactivate

  # Save results
  rsync -azuq --max-size=5M "${_run_dir}"/ "${_brg_rundir}"/
  rm -rf "${_run_dir}"
}

run-polap_genus_species() {
  local _brg_outdir="${1}"
  local _brg_inum="${2:-0}"

  # Set the run title
  local full_name="${FUNCNAME[0]}"
  local middle_part="${full_name#run-}"
  local _run_title="${middle_part%%_*}"
}

run-oatk-nextdenovo_genus_species() {
  local _brg_outdir="${1}"
  local _brg_inum="${2:-0}"

  local full_name="${FUNCNAME[0]}"
  local middle_part="${full_name#run-}"
  local _run_title="${middle_part%%_*}"

  run-oatk_genus_species "${_brg_outdir}" "${_brg_inum}" "nextdenovo"
}

run-oatk-ont_genus_species() {
  local _brg_outdir="${1}"
  local _brg_inum="${2:-0}"

  local full_name="${FUNCNAME[0]}"
  local middle_part="${full_name#run-}"
  local _run_title="${middle_part%%_*}"

  run-oatk_genus_species "${_brg_outdir}" "${_brg_inum}" "ont"
}

run-estimate-genomesize_genus_species() {
  local _brg_outdir="${1}"
  local _brg_inum="${2:-0}"

  # Set the run title
  local full_name="${FUNCNAME[0]}"
  local middle_part="${full_name#run-}"
  local _run_title="${middle_part%%_*}"
  local _run_title="estimate-genomesize"

  # Directory without a slash
  _brg_outdir="${_brg_outdir%/}"

  # Folders
  local _brg_outdir_t="${_brg_outdir}/${opt_t_arg}"
  local _brg_outdir_i="${_brg_outdir_t}/${_brg_inum}"
  local _brg_rundir="${_brg_outdir_i}/${_run_title}"
  local _run_dir="${_brg_outdir}-${_run_title}"
  local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

  # Initialize the run folders
  rm -rf "${_run_dir}"
  mkdir -p "${_run_dir}"
  mkdir -p "${_brg_rundir}"

  # Input data
  local target_index="${_brg_outdir}-${_brg_inum}"
  local long_sra="${_long["$target_index"]}"
  local short_sra="${_short["$target_index"]}"
  if [ -z "${long_sra}" ]; then
    echo "ERROR: no long-read SRA ID: ${long_sra}"
    echo "Suggestion: use -c option for a user-povided CSV."
    return
  fi

  # Files for memtracker
  local _stdout_txt="${_brg_outdir_i}/stdout-${_run_title}.txt"
  local _timing_txt="${_brg_outdir_i}/timing-${_run_title}.txt"
  local _memlog_file="${_brg_outdir_i}/memlog-${_run_title}.csv"
  local _summary_file="${_brg_outdir_i}/summary-${_run_title}.txt"

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

  # Start memory logger
  _polap_lib_process-start_memtracker "${_memlog_file}" \
    "${_polap_var_memtracker_time_interval}"

  command time -v ${_polap_cmd} find-genome-size \
    -a ${short_sra}_1.fastq \
    -b ${short_sra}_2.fastq \
    -o "${_run_dir}" \
    >"${_stdout_txt}" \
    2>"${_timing_txt}"

  # End with summary of the system usage
  _polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"

  # Save system info
  _polap_lib_timing-get_system_info >>"${_timing_txt}"

  conda deactivate

  # Save results
  rsync -azuq --max-size=5M "${_run_dir}"/ "${_brg_rundir}"/
  rm -rf "${_run_dir}"
}

run-nextdenovo-polish_genus_species() {
  local _brg_outdir="${1}"
  local _brg_inum="${2:-0}"
  # local _brg_adir="${3:-.}"
  local _run_title="nextdenovo-polish"

  _brg_outdir="${_brg_outdir%/}"

  local target_index="${_brg_outdir}-${_brg_inum}"
  local long_sra="${_long["$target_index"]}"
  local short_sra="${_short["$target_index"]}"
  if [ -z "${long_sra}" ]; then
    echo "ERROR: no long-read SRA ID: ${long_sra}"
    echo "Suggestion: use -c option for a user-povided CSV."
    return
  fi

  # Folders
  # brg: outdir_t -> outdir_i -> rundir
  # run_dir
  local _brg_outdir_t="${_brg_outdir}/${opt_t_arg}"
  local _brg_outdir_i="${_brg_outdir_t}/${_brg_inum}"
  local _brg_rundir="${_brg_outdir_i}/${_run_title}"
  local _run_dir="${_brg_outdir}-${_run_title}"
  local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

  # Step 0. Get the archive data file if this runs at a remote
  # where we do not have the file.
  # if [[ "${_local_host}" != "$(hostname)" ]]; then
  # 	if [[ ! -s "${_brg_outdir}-a.tar.gz" ]]; then
  # 		scp -p ${_local_host}:$PWD/${_brg_outdir}-a.tar.gz .
  # 	fi
  # fi

  rm -rf "${_run_dir}"
  mkdir -p "${_run_dir}"
  mkdir -p "${_brg_rundir}"

  # if [[ ! -s "${_brg_outdir}/polap.log" ]]; then
  # 	_log_echo "no such file: ${_brg_outdir}/polap.log -> recovering the ${_brg_outdir}"
  # 	recover_genus_species "${_brg_outdir}" "${_brg_inum}"
  # fi

  # NOTE: fofn file needs to be at the present working directory
  # cfg file needs to be at the same directory as well.
  local input_file="${_brg_outdir}_${_brg_inum}_input.fofn"
  echo "${long_sra}.fastq" >"${input_file}"
  cp "${input_file}" "${_brg_outdir_i}"
  local nextdenovo_cfg="${_brg_outdir}-${_brg_inum}-${_run_title}.cfg"

  # NOTE: we may need to precondition a long-read data and a genome size estimate.
  # Step 1. prepare input data
  # prepare-long-data_genus_species "${_brg_outdir}" "${_brg_inum}"
  data-long_genus_species "${_brg_outdir}"

  # Craeat nextdenovo config.
  # NOTE: use one that already exists if there is one.
  # We could manually edit the cfg instead of using the automatically created one
  # because it is often difficult to be optimal conf for NextDenovo.
  # In short, NextDenovo run is not easy to optimize.
  if [[ -s "${nextdenovo_cfg}" ]]; then
    echo "found: ${nextdenovo_cfg}"
  else

    # Step 2. genome size
    # genome_size=$(get_genome_size) || exit 1
    local _estimate_genomesize_rundir="${_brg_outdir_i}/estimate-genomesize"
    if [[ -s "${_estimate_genomesize_rundir}" ]]; then
      local genome_size=$(<"${_estimate_genomesize_rundir}"/short_expected_genome_size.txt)
    else
      echo "Error: no such file: ${_estimate_genomesize_rundir}"
      echo "Suggestion: run estimate-genomesize before run-nextdenovo-polish"
      return
    fi

    generate_nextdenovo_cfg "${genome_size}" "${input_file}" "${_run_dir}" "${nextdenovo_cfg}"
    cp "${nextdenovo_cfg}" "${_brg_outdir_i}"
  fi

  _polap_lib_conda-ensure_conda_env polap-pmat || exit 1

  # Files for memtracker
  local _stdout_txt="${_brg_outdir_i}/stdout-${_run_title}.txt"
  local _timing_txt="${_brg_outdir_i}/timing-${_run_title}.txt"
  local _memlog_file="${_brg_outdir_i}/memlog-${_run_title}.csv"
  local _summary_file="${_brg_outdir_i}/summary-${_run_title}.txt"

  _log_echo "nextDenovo error-correction on ${_brg_outdir} at dir:${_brg_outdir}/${_brg_adir}"

  # Start memory logger
  _polap_lib_process-start_memtracker "${_memlog_file}" \
    "${_polap_var_memtracker_time_interval}"

  # Execute NextDenovo
  command time -v nextDenovo "${nextdenovo_cfg}" \
    >"${_stdout_txt}" \
    2>"${_timing_txt}"

  _polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"

  cat "${_run_dir}"/02.cns_align/01.seed_cns.sh.work/seed_cns*/cns.fasta >"${_brg_outdir_i}"/cns.fa

  # Save system info
  _polap_lib_timing-get_system_info >>"${_timing_txt}"

  conda deactivate

  # Clean up the nextDenovo working folder
  # if [[ "${_local_host}" != "$(hostname)" ]]; then
  #   echo rsync -aq "${_brg_outdir_i}/" "${_local_host}:${PWD}/${_brg_outdir_i}/"
  # fi
  rm -rf "${_run_dir}"
}

run-pmat_genus_species() {
  local _brg_outdir="${1}"
  local _brg_inum="${2:-0}"
  local _brg_type="${3:-nextdenovo}"
  local _run_title="pmat-${_brg_type}"

  source <(echo 'export PATH="$PWD/PMAT-1.5.3/bin:$PATH"')

  # Directory without a slash
  _brg_outdir="${_brg_outdir%/}"

  # Folders
  # brg: outdir_t -> outdir_i -> rundir
  # run_dir
  local _brg_outdir_t="${_brg_outdir}/${opt_t_arg}"
  local _brg_outdir_i="${_brg_outdir_t}/${_brg_inum}"
  local _brg_rundir="${_brg_outdir_i}/${_run_title}"
  local _run_dir="${_brg_outdir}-${_run_title}"
  local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

  # Create folders
  mkdir -p "${_brg_rundir}"
  rm -rf "${_run_dir}"
  mkdir -p "${_run_dir}"

  local target_index="${_brg_outdir}-${_brg_inum}"
  local long_sra="${_long["$target_index"]}"
  if [[ -z "$long_sra" ]]; then
    echo "Error: no long SRA for ${_brg_outdir}-${_brg_inum}"
    echo "Info: skipping pmat-nextdenovo on ${_brg_outdir}-${_brg_inum}"
    return
  fi

  # Step 1. input data
  local _brg_input_data="${_brg_outdir_i}"/cns.fa
  if [[ "${_brg_type}" == "ont" ]]; then
    _brg_input_data="${long_sra}.fastq"
  fi

  # Step 2. genome size
  local _estimate_genomesize_rundir="${_brg_outdir_i}/estimate-genomesize"
  if [[ -s "${_estimate_genomesize_rundir}" ]]; then
    local genome_size=$(<"${_estimate_genomesize_rundir}"/short_expected_genome_size.txt)
  else
    echo "Error: no such file: ${_estimate_genomesize_rundir}"
    echo "Suggestion: run estimate-genomesize before run-nextdenovo-polish"
    return
  fi

  _polap_lib_conda-ensure_conda_env polap-pmat || exit 1

  # local fc_list=()
  # local _brg_fc="0.1,1.0"
  # if [[ "${_brg_fc}" == *,* ]]; then
  #   IFS=',' read -ra fc_list <<<"$_brg_fc"
  # else
  #   fc_list=("$_brg_fc")
  # fi

  local _fc=0.1
  for _fc in "${_polap_pmat_options[@]}"; do
    _log_echo "pmat on ${_brg_rundir} using the nextDenovo-polished data: ${_brg_outdir_i}/cns.fa with -fc ${_fc}"
    # local formatted_fc=$(printf "%02d" "${_fc}")

    # https://github.com/c-zhou/pmat

    # NOTE: pmat threads 56 -> seg. fault.
    # local _pmat_threads=8

    # Files for memtracker
    local _stdout_txt="${_brg_outdir_i}/stdout-${_run_title}-${_fc}.txt"
    local _timing_txt="${_brg_outdir_i}/timing-${_run_title}-${_fc}.txt"
    local _memlog_file="${_brg_outdir_i}/memlog-${_run_title}-${_fc}.csv"
    local _summary_file="${_brg_outdir_i}/summary-${_run_title}-${_fc}.txt"

    rm -f "${_stdout_txt}"
    rm -f "${_timing_txt}"

    # Start memory logger
    _polap_lib_process-start_memtracker "${_memlog_file}" \
      "${_polap_var_memtracker_time_interval}"

    rm -rf "${_run_dir}"

    # command time -v PMAT autoMito \
    command time -v timeout 6h PMAT autoMito \
      -i "${_brg_input_data}" \
      -o "${_run_dir}" \
      -st ont \
      -g ${genome_size} \
      --task p1 \
      --type all \
      -cpu ${_brg_threads} \
      -fc ${_fc} \
      -m \
      >"${_stdout_txt}" \
      2>"${_timing_txt}"

    # Save results
    mkdir -p "${_brg_rundir}/${_fc}"
    # rsync -azuq "${_run_dir}"/ "${_brg_rundir}/${_fc}"/
    rsync -azuq --max-size=5M \
      "${_run_dir}"/ "${_brg_rundir}/${_fc}"/

    # Save system info
    _polap_lib_timing-get_system_info >>"${_timing_txt}"

    _polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"

  done

  conda deactivate

  # Save results
  # rsync -azuq "${_run_dir}"/ "${_brg_rundir}"/

  # Clean-up
  rm -rf "${_run_dir}"
}

run-tippo_genus_species() {
  local _brg_outdir="${1}"
  local _brg_inum="${2:-0}"
  local _brg_type="${3:-nextdenovo}"

  local _run_title="tippo-${_brg_type}"

  # Directory without a slash
  _brg_outdir="${_brg_outdir%/}"

  # Folders
  # brg: outdir_t -> outdir_i -> rundir
  # run_dir
  local _brg_outdir_t="${_brg_outdir}/${opt_t_arg}"
  local _brg_outdir_i="${_brg_outdir_t}/${_brg_inum}"
  local _brg_rundir="${_brg_outdir_i}/${_run_title}"
  local _run_dir="${_brg_outdir}-${_run_title}"
  local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

  # Create folders
  mkdir -p "${_brg_rundir}"
  rm -rf "${_run_dir}"
  mkdir -p "${_run_dir}"

  local target_index="${_brg_outdir}-${_brg_inum}"
  local long_sra="${_long["$target_index"]}"
  if [[ -z "$long_sra" ]]; then
    echo "Error: no long SRA for ${_brg_outdir}-${_brg_inum}"
    echo "Info: skipping tippo-nextdenovo on ${_brg_outdir}-${_brg_inum}"
    return
  fi

  local _brg_input_data="${_brg_outdir_i}"/cns.fa
  if [[ "${_brg_type}" == "ont" ]]; then
    _brg_input_data="${long_sra}.fastq"
  fi

  _polap_lib_conda-ensure_conda_env polap-tippo || exit 1

  # local fc_list=()
  # # local _brg_fc="hifi,clr,ont,onthq"
  # local _brg_fc="onthq,ont"
  # if [[ "${_brg_fc}" == *,* ]]; then
  #   IFS=',' read -ra fc_list <<<"$_brg_fc"
  # else
  #   fc_list=("$_brg_fc")
  # fi

  # TIPPo has no option for output directory
  local _tippo_actual_outdir="cns.fa.organelle"

  for _fc in "${_polap_tippo_options[@]}"; do
    _log_echo "tippo on ${_brg_rundir} using the nextDenovo-polished data: ${_brg_outdir_i}/cns.fa with -p ${_fc}"

    # https://github.com/Wenfei-Xian/TIPP

    # Files for memtracker
    local _stdout_txt="${_brg_outdir_i}/stdout-${_run_title}-${_fc}.txt"
    local _timing_txt="${_brg_outdir_i}/timing-${_run_title}-${_fc}.txt"
    local _memlog_file="${_brg_outdir_i}/memlog-${_run_title}-${_fc}.csv"
    local _summary_file="${_brg_outdir_i}/summary-${_run_title}-${_fc}.txt"

    rm -f "${_stdout_txt}"
    rm -f "${_timing_txt}"

    # Start memory logger
    _polap_lib_process-start_memtracker "${_memlog_file}" \
      "${_polap_var_memtracker_time_interval}"

    rm -rf "${_tippo_actual_outdir}"

    # if [[ "${_brg_fc}" == "1" ]]; then
    command time -v timeout 6h TIPPo.v2.4.pl \
      -f "${_brg_input_data}" \
      -t ${_brg_threads} \
      -p ${_fc} \
      >"${_stdout_txt}" \
      2>"${_timing_txt}"

    # Save results
    mkdir -p "${_brg_rundir}/${_fc}"
    rsync -azuq "${_tippo_actual_outdir}"/ "${_brg_rundir}/${_fc}"/

    # Save system info
    _polap_lib_timing-get_system_info >>"${_timing_txt}"

    _polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"

  done

  conda deactivate

  # Output files
  # cns.fa.chloroplast.fasta.filter.800.round1.edge_*.edge_*.edge_*.organelle.chloroplast.fasta
  # cns.fa.mitochondrial.fasta.filter.fasta.flye/assembly_graph.gfa

  # Clean-up
  rm -rf "${_run_dir}"
  rm -rf "${_tippo_actual_outdir}"
}

run-oatk_genus_species() {
  local _brg_outdir="${1}"
  local _brg_inum="${2:-0}"
  local _brg_type="${3:-nextdenovo}"
  local _run_title="oatk-${_brg_type}"

  if [[ "${_brg_outdir}" == "-h" || "${_brg_outdir}" == "--help" ]]; then
    echo "$help_message_run_oatk"
    return
  fi

  # Directory without a slash
  _brg_outdir="${_brg_outdir%/}"

  # Folders
  # brg: outdir_t -> outdir_i -> rundir
  # run_dir
  local _brg_outdir_t="${_brg_outdir}/${opt_t_arg}"
  local _brg_outdir_i="${_brg_outdir_t}/${_brg_inum}"
  local _brg_rundir="${_brg_outdir_i}/${_run_title}"
  local _run_dir="${_brg_outdir}-${_run_title}"
  local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

  # Create folders
  mkdir -p "${_brg_rundir}"
  rm -rf "${_run_dir}"
  mkdir -p "${_run_dir}"

  local target_index="${_brg_outdir}-${_brg_inum}"
  local long_sra="${_long["$target_index"]}"
  if [[ -z "$long_sra" ]]; then
    echo "Error: no long SRA for ${_brg_outdir}-${_brg_inum}"
    echo "Info: skipping oatk-nextdenovo on ${_brg_outdir}-${_brg_inum}"
    return
  fi

  local _brg_input_data="${_brg_outdir_i}"/cns.fa
  if [[ "${_brg_type}" == "ont" ]]; then
    _brg_input_data="${long_sra}.fastq"
  fi

  _polap_lib_conda-ensure_conda_env polap-oatk || exit 1

  # local fc_list=()
  # local _brg_fc="30,20"
  # if [[ "${_brg_fc}" == *,* ]]; then
  #   IFS=',' read -ra fc_list <<<"$_brg_fc"
  # else
  #   fc_list=("$_brg_fc")
  # fi

  local _fc=1
  for _fc in "${_polap_oatk_options[@]}"; do
    _log_echo "oatk on ${_brg_rundir} using the nextDenovo-polished data: ${_brg_outdir_i}/cns.fa with -c ${_fc}"
    local formatted_fc=$(printf "%02d" "${_fc}")

    # Files for memtracker
    local _stdout_txt="${_brg_outdir_i}/stdout-${_run_title}-${formatted_fc}.txt"
    local _timing_txt="${_brg_outdir_i}/timing-${_run_title}-${formatted_fc}.txt"
    local _memlog_file="${_brg_outdir_i}/memlog-${_run_title}-${formatted_fc}.csv"
    local _summary_file="${_brg_outdir_i}/summary-${_run_title}-${formatted_fc}.txt"

    rm -f "${_stdout_txt}"
    rm -f "${_timing_txt}"

    # Start memory logger
    _polap_lib_process-start_memtracker "${_memlog_file}" \
      "${_polap_var_memtracker_time_interval}"

    # https://github.com/c-zhou/oatk

    # NOTE: oatk threads 56 -> seg. fault.
    local _oatk_threads=8
    command time -v timeout 6h oatk \
      -k 1001 \
      -c ${_fc} \
      -t ${_oatk_threads} \
      -m ./OatkDB/v20230921/embryophyta_mito.fam \
      -p ./OatkDB/v20230921/embryophyta_pltd.fam \
      -o "${_run_dir}/oatk-${_brg_type}-${formatted_fc}" \
      "${_brg_input_data}" \
      >"${_stdout_txt}" \
      2>"${_timing_txt}"

    # Save system info
    _polap_lib_timing-get_system_info >>"${_timing_txt}"

    _polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"

  done

  conda deactivate

  # Save results
  rsync -azuq "${_run_dir}"/ "${_brg_rundir}"/

  # Clean-up
  rm -rf "${_run_dir}"
}

run-extract-ptdna-ptgaul_genus_species() {
  local _brg_outdir="$1"
  local _brg_inum="${2:-0}"

  if [[ "${_brg_outdir}" == "-h" || "${_brg_outdir}" == "--help" ]]; then
    echo "$help_message_run_extract_ptdna_ptgaul"
    return
  fi

  # Directory without a slash
  _brg_outdir="${_brg_outdir%/}"

  # Config CSV
  local _run_title="extract-ptdna-ptgaul"

  # Folders
  local _brg_outdir_t="${_brg_outdir}/${opt_t_arg}"
  local _brg_outdir_i="${_brg_outdir_t}/${_brg_inum}"
  local _brg_rundir="${_brg_outdir_i}/${_run_title}"
  local _run_dir="${_brg_outdir}-${_run_title}"
  local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

  # Files for memtracker
  local _stdout_txt="${_brg_outdir_i}/stdout-${_run_title}.txt"
  local _timing_txt="${_brg_outdir_i}/timing-${_run_title}.txt"
  local _memlog_file="${_brg_outdir_i}/memlog-${_run_title}.csv"
  local _summary_file="${_brg_outdir_i}/summary-${_run_title}.txt"

  # Initialize Conda
  _polap_lib_conda-ensure_conda_env polap || exit 1

  # extract ptGAUL result
  _log_echo "extract ptDNA from the ptGAUL result"
  # echo command time -v ${_polap_cmd} disassemble ptgaul \
  # 	-o ${_brg_outdir_i}

  # Too fast to catch
  _polap_lib_process-start_memtracker "${_memlog_file}" 3
  # "${_polap_var_memtracker_time_interval}"

  # command time -v ${_polap_cmd} disassemble ptgaul \
  #   -o ${_brg_outdir_i} \
  #   >"${_stdout_txt}" \
  #   2>"${_timing_txt}"

  command time -v ${_polap_cmd} ptgaul extract-ptdna \
    -o ${_brg_outdir_i} \
    >"${_stdout_txt}" \
    2>"${_timing_txt}"

  conda deactivate

  # End with summary of the system usage
  _polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"

  # Save system info
  _polap_lib_timing-get_system_info >>"${_timing_txt}"

  # copy ptGAUL result
  local _outdir="${_brg_outdir_i}/ptgaul/flye_cpONT/ptdna"
  local _brg_unpolished_fasta="${_outdir}/circular_path_1_concatenated.fa"
  if [[ ! -s "${_brg_unpolished_fasta}" ]]; then
    _log_echo "[Fail] extraction of a draft ptDNA"
    echo "Check (tail): ${_brg_outdir_i}/polap.log"
  fi
}

run-polish-ptdna-ptgaul_genus_species() {
  local _brg_outdir="$1"
  local _brg_inum="${2:-0}"

  if [[ "${_brg_outdir}" == "-h" || "${_brg_outdir}" == "--help" ]]; then
    echo "$help_message_run_polish_ptdna_ptgaul"
    return
  fi

  # Directory without a slash
  _brg_outdir="${_brg_outdir%/}"

  # Config CSV
  local _run_title="polish-ptdna-ptgaul"

  # Folders
  local _brg_outdir_t="${_brg_outdir}/${opt_t_arg}"
  local _brg_outdir_i="${_brg_outdir_t}/${_brg_inum}"
  local _brg_rundir="${_brg_outdir_i}/${_run_title}"
  local _run_dir="${_brg_outdir}-${_run_title}"
  local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

  # Files for memtracker
  local _stdout_txt="${_brg_outdir_i}/stdout-${_run_title}.txt"
  local _timing_txt="${_brg_outdir_i}/timing-${_run_title}.txt"
  local _memlog_file="${_brg_outdir_i}/memlog-${_run_title}.csv"
  local _summary_file="${_brg_outdir_i}/summary-${_run_title}.txt"

  # Initialize Conda
  _polap_lib_conda-ensure_conda_env polap || exit 1

  # polish ptGAUL result
  _log_echo "polish ptDNA from the ptGAUL result with fmlrc polishing"
  # echo command time -v ${_polap_cmd} disassemble ptgaul \
  # 	-o ${_brg_outdir_i}

  # It could be too fast to catch the system info.
  _polap_lib_process-start_memtracker "${_memlog_file}" 10
  # "${_polap_var_memtracker_time_interval}"

  # command time -v ${_polap_cmd} disassemble ptgaul \
  #   -o ${_brg_outdir_i} \
  #   >"${_stdout_txt}" \
  #   2>"${_timing_txt}"

  command time -v ${_polap_cmd} ptgaul polish-ptdna \
    -o ${_brg_outdir_i} \
    >"${_stdout_txt}" \
    2>"${_timing_txt}"

  conda deactivate

  # End with summary of the system usage
  _polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"

  # Save system info
  _polap_lib_timing-get_system_info >>"${_timing_txt}"

  # copy ptGAUL result
  local _outdir="${_brg_outdir_i}/ptgaul/flye_cpONT/ptdna"
  local _brg_final_assembly="${_outdir}/pt.1.fa"
  if [[ -s "${_brg_final_assembly}" ]]; then
    cp -p ${_brg_final_assembly} ${_brg_outdir_i}/ptdna-ptgaul.fa
    # _log_echo "[Success] extraction of a polished ptDNA"
  else
    _log_echo "[Fail] extraction of a polished ptDNA"
  fi
}

run-ptgaul_genus_species() {
  local _brg_outdir="$1"
  local _brg_inum="${2:-0}"

  # Directory without a slash
  _brg_outdir="${_brg_outdir%/}"

  # Config CSV
  local _run_title="ptgaul"
  local target_index="${_brg_outdir}-${_brg_inum}"
  local long_sra="${_long["$target_index"]}"
  local extracted_ptgaul_genomesize="${_ptgaul_genomesize["$target_index"]}"

  # Check the config
  if [[ -z "$long_sra" ]]; then
    echo "Error: skipping ${_brg_outdir}-${_brg_inum} because it is not in the CSV."
    return
  fi

  # Folders
  local _brg_outdir_t="${_brg_outdir}/${opt_t_arg}"
  local _brg_outdir_i="${_brg_outdir_t}/${_brg_inum}"
  local _brg_rundir="${_brg_outdir_i}/${_run_title}"
  local _run_dir="${_brg_outdir}-${_run_title}"
  local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

  # Create the output folder
  mkdir -p "${_brg_rundir}"

  # Initialize Conda
  _polap_lib_conda-ensure_conda_env polap || exit 1
  # source "$(conda info --base)/etc/profile.d/conda.sh"
  # if [[ "${CONDA_DEFAULT_ENV:-}" != "polap" ]]; then
  # 	if [[ "$_POLAP_DEBUG" == "1" ]]; then
  # 		echo "[INFO] Activating conda environment 'polap'..."
  # 	fi
  # 	conda activate polap
  # fi
  #
  # if [[ "${CONDA_DEFAULT_ENV:-}" != "polap" ]]; then
  # 	echo "[ERROR] Failed to enter conda environment 'polap'"
  # 	return
  # fi

  # Files for memtracker
  local _stdout_txt="${_brg_outdir_i}/stdout-${_run_title}.txt"
  local _timing_txt="${_brg_outdir_i}/timing-${_run_title}.txt"
  local _memlog_file="${_brg_outdir_i}/memlog-${_run_title}.csv"
  local _summary_file="${_brg_outdir_i}/summary-${_run_title}.txt"

  # Start memory logger
  # GetOrganelle takes not much time to finish, and we use 10 seconds of interval for memlog
  _polap_lib_process-start_memtracker "${_memlog_file}" \
    "${_polap_var_memtracker_time_interval}"

  mkdir -p "${_brg_outdir}/timing"
  command time -v bash ${_POLAPLIB_DIR}/polap-ptGAUL1.sh \
    -o ${_run_dir} \
    -r ${_brg_outdir_i}/ncbi-ptdna/ptdna-reference.fa \
    -g "${extracted_ptgaul_genomesize}" \
    -l ${long_sra}.fastq \
    -t ${_brg_threads} \
    >"${_stdout_txt}" \
    2>"${_timing_txt}"

  # End with summary of the system usage
  _polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"

  # Save system info
  _polap_lib_timing-get_system_info >>"${_timing_txt}"

  conda deactivate

  # Save some results
  rsync -azuq --max-size=5M \
    "${_run_dir}"/result_3000/ \
    "${_brg_rundir}"/

  rm -rf "${_run_dir}"
}

# NOTE: we use _brg_rundir not _run_dir to download the data.
download-ptdna_genus_species() {
  local _brg_outdir="$1"
  local _brg_inum="${2:-0}"

  # Directory without a slash
  _brg_outdir="${_brg_outdir%/}"

  local _run_title="ncbi-ptdna"
  local target_index="${_brg_outdir}-${_brg_inum}"
  local _brg_outdir_t="${_brg_outdir}/${opt_t_arg}"
  local _brg_outdir_i="${_brg_outdir_t}/${_brg_inum}"
  local _brg_rundir="${_brg_outdir_i}/${_run_title}"
  local _run_dir="${_brg_outdir}-${_run_title}"
  local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"
  local _species=${_taxon[$target_index]}
  local _ref_species=${_ptgaul_ref[$target_index]}
  local species_name="${_species//_/ }"

  mkdir -p "${_brg_rundir}"

  # Initialize Conda
  _polap_lib_conda-ensure_conda_env polap || exit 1

  if [[ "${_ref_species}" != "NA" ]]; then
    species_name="${_ref_species}"
    echo "No ptDNA for ${_brg_outdir}, so we use ${species_name}"
  fi

  if [[ -s "${_brg_outdir_i}/ptdna-reference.fa" ]]; then
    echo "found: ptDNA reference: ${_brg_outdir}/ptdna-reference.fa"
  else

    ${_polap_cmd} get-mtdna \
      --plastid \
      --species "${species_name}" \
      -o ${_brg_rundir}

    if [[ -s "${_brg_rundir}/00-bioproject/2-mtdna.fasta" ]]; then
      echo "copy ${_brg_rundir}/ptdna-reference.fa"
      cp -p "${_brg_rundir}/00-bioproject/2-mtdna.fasta" \
        "${_brg_rundir}/ptdna-reference.fa"
    else
      echo "No such file: ${_brg_rundir}/ptdna-reference.fa"
      local _genus_name=$(echo ${species_name} | awk '{print $1}')
      echo "  trying to search NCBI plastid genomes for genus name only: ${_genus_name}"
      ${_polap_cmd} get-mtdna \
        --plastid \
        --species "${_genus_name}" \
        -o ${_brg_rundir}
      if [[ -s "${_brg_rundir}/00-bioproject/2-mtdna.fasta" ]]; then
        echo "copy ${_brg_rundir}/ptdna-reference.fa"
        cp -p "${_brg_rundir}/00-bioproject/2-mtdna.fasta" \
          "${_brg_rundir}/ptdna-reference.fa"
      else
        echo "  we could not find one even in the genus level."
        echo "No such file: ${_brg_rundir}/ptdna-reference.fa"
      fi
    fi
  fi

  conda deactivate
}

run-msbwt_genus_species() {
  # echo "${FUNCNAME[0]}: $@"
  local _brg_outdir="$1"
  local _brg_inum="${2:-0}"

  # Directory without a slash
  _brg_outdir="${_brg_outdir%/}"

  local _run_title="msbwt"
  local target_index="${_brg_outdir}-${_brg_inum}"
  local short_sra="${_short["$target_index"]}"
  local _brg_outdir_t="${_brg_outdir}/${opt_t_arg}"
  local _brg_outdir_i="${_brg_outdir_t}/${_brg_inum}"
  local _brg_rundir="${_brg_outdir_i}/${_run_title}"
  local _run_dir="${_brg_outdir}-${_run_title}"
  local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

  # Check the input short-read data
  mkdir -p "${_brg_outdir_i}"/msbwt

  # Initialize Conda
  _polap_lib_conda-ensure_conda_env polap-fmlrc || exit 1
  # _polap_lib_conda-ensure_conda_env polap-fmlrc2 || exit 1

  # Files for memtracker
  local _stdout_txt="${_brg_outdir_i}/stdout-${_run_title}.txt"
  local _timing_txt="${_brg_outdir_i}/timing-${_run_title}.txt"
  local _memlog_file="${_brg_outdir_i}/memlog-${_run_title}.csv"
  local _summary_file="${_brg_outdir_i}/summary-${_run_title}.txt"

  # Start memory logger
  # GetOrganelle takes not much time to finish, and we use 10 seconds of interval for memlog
  _polap_lib_process-start_memtracker "${_memlog_file}" \
    "${_polap_var_memtracker_time_interval}"

  command time -v ${_polap_cmd} prepare-polishing \
    -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
    -o ${_run_dir} \
    >"${_stdout_txt}" \
    2>"${_timing_txt}"

  conda deactivate

  # End with summary of the system usage
  _polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"

  # Save system info
  _polap_lib_timing-get_system_info >>"${_timing_txt}"

  # Save some results
  rsync -azuq "${_run_dir}"/msbwt/ "${_brg_outdir_i}"/msbwt/
  rm -rf "${_run_dir}"
}

run_genus_species() {
  local first_arg="$1"
  local remaining_args=("${@:2}")

  # Remove trailing slash from the first element
  remaining_args[0]="${remaining_args[0]%/}"

  run-${first_arg}_genus_species "${remaining_args[@]}"
}

run-getorganelle_genus_species() {
  # echo "${FUNCNAME[0]}: $@"
  local _brg_outdir="$1"
  local _brg_inum="${2:-0}"

  # Directory without a slash
  _brg_outdir="${_brg_outdir%/}"

  local _run_title="getorganelle"
  local target_index="${_brg_outdir}-${_brg_inum}"
  local short_sra="${_short["$target_index"]}"
  local _brg_outdir_t="${_brg_outdir}/${opt_t_arg}"
  local _brg_outdir_i="${_brg_outdir_t}/${_brg_inum}"
  local _brg_rundir="${_brg_outdir_i}/${_run_title}"
  local _run_dir="${_brg_outdir}-${_run_title}"
  local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

  # Check the input short-read data
  mkdir -p "${_brg_outdir_i}"

  # Initialize Conda
  _polap_lib_conda-ensure_conda_env polap-getorganelle || exit 1

  # Files for memtracker
  local _stdout_txt="${_brg_outdir_i}/stdout-${_run_title}.txt"
  local _timing_txt="${_brg_outdir_i}/timing-${_run_title}.txt"
  local _memlog_file="${_brg_outdir_i}/memlog-${_run_title}.csv"
  local _summary_file="${_brg_outdir_i}/summary-${_run_title}.txt"

  # Start memory logger
  # GetOrganelle takes not much time to finish, and we use 10 seconds of interval for memlog
  _polap_lib_process-start_memtracker "${_memlog_file}" 10

  command time -v get_organelle_from_reads.py \
    -1 ${short_sra}_1.fastq \
    -2 ${short_sra}_2.fastq \
    -o ${_run_dir} \
    -t ${_brg_threads} \
    -F embplant_pt \
    >"${_stdout_txt}" \
    2>"${_timing_txt}"

  conda deactivate

  # End with summary of the system usage
  _polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"

  # Save system info
  _polap_lib_timing-get_system_info >>"${_timing_txt}"

  # Save some results
  mkdir -p "${_brg_rundir}"
  rsync -azuq \
    --max-size=5M \
    "${_run_dir}"/ \
    "${_brg_rundir}"/
  rm -rf "${_run_dir}"

}

install-pmat_genus_species() {
  if [[ "${opt_y_flag}" == false ]]; then
    read -p "Do you want to install pmat in the polap-pmat conda environment? (y/N): " confirm
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
      echo "You're in the base environment. Creating 'polap-pmat'..."
      if conda env list | awk '{print $1}' | grep -qx "polap-pmat"; then
        echo "ERROR: Conda environment 'polap-pmat' already exists."
      else
        conda create -y --name polap-pmat apptainer nextdenovo canu blast
        finalize-pmat_genus_species
      fi
    else
      echo "Error: You're in the '$CONDA_DEFAULT_ENV' environment. Please activate base before running this script."
      exit 1
    fi
  else
    echo "pmat installation is canceled."
    echo "Ref: https://github.com/bichangwei/PMAT"
    echo "Ref: wget https://github.com/bichangwei/PMAT/archive/refs/tags/v1.5.3.tar.gz"
    echo "${help_message_install_pmat}"
  fi
}

install-bandage_genus_species() {
  if [[ "${opt_y_flag}" == false ]]; then
    read -p "Do you want to install bandage? (y/N): " confirm
  else
    confirm="yes"
  fi

  if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
    wget -q https://github.com/rrwick/Bandage/releases/download/v0.8.1/Bandage_Ubuntu_static_v0_8_1.zip
    unzip Bandage_Ubuntu_static_v0_8_1.zip
    mkdir -p $HOME/bin
    cp Bandage "$HOME/bin"
    # source <(echo 'export PATH="$PWD/bin:$PATH"')
  else
    echo "bandage installation is canceled."
    echo "${help_message_install_bandage}"
  fi
}

download-test-data2_genus_species() {
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
}

download_genus_species() {
  local first_arg="$1"
  local remaining_args=("${@:2}")

  if ! declare -F download-${first_arg}_genus_species >/dev/null; then
    echo "No such download: ${first_arg}"
    return 1
  fi

  download-${first_arg}_genus_species "${remaining_args[@]}"
}

delete-polap-github_genus_species() {
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
}

delete_genus_species() {
  local args=("$@")

  for item in "$@"; do
    delete-${item}_genus_species
  done
}

install-dflye_genus_species() {
  if [[ "${opt_y_flag}" == false ]]; then
    read -p "Do you want to install dflye in the polap-dflye conda environment? (y/N): " confirm
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
      echo "You're in the base environment. Creating 'polap-dflye'..."
      if conda env list | awk '{print $1}' | grep -qx "polap-dflye"; then
        echo "ERROR: Conda environment 'polap-dflye' already exists."
      else
        conda create -y --name polap-dflye goshng::dflye
        conda activate polap-dflye
        conda install -y flye=2.9.5
      fi
    else
      echo "Error: You're in the '$CONDA_DEFAULT_ENV' environment. Please activate base before running this script."
      exit 1
    fi
  else
    echo "dflye or read-coverage filtering version Flye installation is canceled."
    echo "Ref: https://anaconda.org/goshng/dflye"
    echo "Execute: dflye"
  fi
}

install-cflye_genus_species() {
  if [[ "${opt_y_flag}" == false ]]; then
    read -p "Do you want to install cflye in the polap-cflye conda environment? (y/N): " confirm
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
      echo "You're in the base environment. Creating 'polap-cflye'..."
      if conda env list | awk '{print $1}' | grep -qx "polap-cflye"; then
        echo "ERROR: Conda environment 'polap-cflye' already exists."
      else
        conda create -y --name polap-cflye goshng::cflye
      fi
    else
      echo "Error: You're in the '$CONDA_DEFAULT_ENV' environment. Please activate base before running this script."
      exit 1
    fi
  else
    echo "cflye or read-coverage filtering version Flye installation is canceled."
    echo "Ref: https://anaconda.org/goshng/cflye"
    echo "Execute: cflye"
  fi
}

install-polap_genus_species() {
  local _version="${1}"

  if [[ "${opt_y_flag}" == false ]]; then
    if [[ -n "${_version}" ]]; then
      read -p "Do you want to install polap ${_version}? (y/N): " confirm
    else
      read -p "Do you want to install polap? (y/N): " confirm
    fi

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
        if [[ -n "${_version}" ]]; then
          conda create -y --name polap bioconda::polap="${_version}"
        else
          conda create -y --name polap bioconda::polap
        fi
      fi
    else
      echo "Error: You're in the '$CONDA_DEFAULT_ENV' environment. Please activate base before running this script."
      exit 1
    fi
  else
    echo "polap installation is canceled."
    echo "${help_message_install_polap}"
  fi
}

install-fmlrc_genus_species() {
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
        local v3=0.3.7.3

        wget -q https://github.com/goshng/polap/archive/refs/tags/${v3}.zip
        if [[ -s "${v3}.zip" ]]; then
          unzip -o -q ${v3}.zip
          conda env create -f polap-${v3}/src/polap-conda-environment-fmlrc.yaml
        else
          echo "Error: no such file: ${v3}.zip"
          # echo "Suggestion: _polap_version=0.4.3.7.4 $0 $subcmd1"
        fi

      fi
    else
      echo "Error: You're in the '$CONDA_DEFAULT_ENV' environment. Please activate base before running this script."
      exit 1
    fi
  else
    echo "polap-fmlrc installation is canceled."
    echo "${help_message_install_fmlrc}"
  fi
}

setup-conda_genus_species() {
  if [[ "${opt_y_flag}" == false ]]; then
    read -p "Do you want to setup miniconda3 for conda-forge and bioconda? (y/N): " confirm
  else
    confirm="yes"
  fi

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
}

setup_genus_species() {
  local args=("$@")
  local first_arg="$1"
  local remaining_args=("${@:2}")

  if has_help "${args[@]}"; then
    echo "$help_message_setup"
    return
  fi

  setup-${first_arg}_genus_species "${remaining_args[@]:-}"
}

list_genus_species() {
  local _brg_query="${1}"
  local _brg_where="${2:-any}"

  if [[ "${_brg_where}" == "any" ]]; then
    grep ")$" "${_POLAPLIB_DIR}/polap-lib-data.sh" |
      grep -v "(" | grep "${_brg_query}" | clean_input_lines | sort | uniq
    grep ")$" $0 | grep -v "(" | grep "${_brg_query}" | clean_input_lines | sort | uniq
  elif [[ "${_brg_where}" == "end" ]]; then
    grep ")$" "${_POLAPLIB_DIR}/polap-lib-data.sh" |
      grep -v "(" | grep "${_brg_query})$" | clean_input_lines | sort | uniq
    grep ")$" $0 | grep -v "(" | grep "${_brg_query})$" | clean_input_lines | sort | uniq
  else
    grep ")$" "${_POLAPLIB_DIR}/polap-lib-data.sh" |
      clean_input_lines | grep -v "(" | grep "^${_brg_query}" | sort | uniq
    grep ")$" $0 | clean_input_lines | grep -v "(" | grep "^${_brg_query}" | sort | uniq
  fi
}

install_genus_species() {
  local args=("$@")

  if has_help "${args[@]}"; then
    echo "$help_message_install"
    return
  fi

  for item in "${args[@]}"; do
    local tool="${item%%=*}"
    local version=""
    [[ "$item" == *"="* ]] && version="${item#*=}"

    install-${tool}_genus_species "$version"
  done
}

install-efg_genus_species() {
  if [[ "${opt_y_flag}" == false ]]; then
    read -p "Do you want to install efg in the polap-efg conda environment? (y/N): " confirm
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
      echo "You're in the base environment. Creating 'polap-efg'..."
      if conda env list | awk '{print $1}' | grep -qx "polap-efg"; then
        echo "ERROR: Conda environment 'polap-efg' already exists."
      else
        echo conda create -y --name polap-efg efg
      fi
    else
      echo "Error: You're in the '$CONDA_DEFAULT_ENV' environment. Please activate base before running this script."
      exit 1
    fi
  else
    echo "efg installation is canceled."
    echo "Check: https://github.com/maickrau/efg"
    echo "conda install -c bioconda efg"
    echo "Execute: efg"
  fi
}

install-abc_genus_species() {
  if [[ "${opt_y_flag}" == false ]]; then
    read -p "Do you want to install abc in the polap-abc conda environment? (y/N): " confirm
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
      echo "You're in the base environment. Creating 'polap-abc'..."
      if conda env list | awk '{print $1}' | grep -qx "polap-abc"; then
        echo "ERROR: Conda environment 'polap-abc' already exists."
      else
        echo conda create -y --name polap-abc abc
      fi
    else
      echo "Error: You're in the '$CONDA_DEFAULT_ENV' environment. Please activate base before running this script."
      exit 1
    fi
  else
    echo "abc installation is canceled."
    echo "Check: https://github.com/maickrau/abc"
    echo "conda install -c bioconda abc"
    echo "Execute: abc"
  fi
}

install-getorganelle_genus_species() {
  if [[ "${opt_y_flag}" == false ]]; then
    read -p "Do you want to install getorganelle in the polap-getorganelle conda environment? (y/N): " confirm
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
      echo "You're in the base environment. Creating 'polap-getorganelle'..."
      if conda env list | awk '{print $1}' | grep -qx "polap-getorganelle"; then
        echo "ERROR: Conda environment 'polap-getorganelle' already exists."
      else
        conda create -y --name polap-getorganelle getorganelle
        conda activate polap-getorganelle
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
}

install-tippo_genus_species() {
  if [[ "${opt_y_flag}" == false ]]; then
    read -p "Do you want to install tippo in the polap-tippo conda environment? (y/N): " confirm
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
      echo "You're in the base environment. Creating 'polap-tippo'..."
      if conda env list | awk '{print $1}' | grep -qx "polap-tippo"; then
        echo "ERROR: Conda environment 'polap-tippo' already exists."
      else
        #please specific the python version to 3.8 :)
        conda create -y --name polap-tippo python=3.8
        conda activate polap-tippo
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
}

install-oatk_genus_species() {
  if [[ "${opt_y_flag}" == false ]]; then
    read -p "Do you want to install oatk in the polap-oatk conda environment? (y/N): " confirm
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
      echo "You're in the base environment. Creating 'polap-oatk'..."
      if conda env list | awk '{print $1}' | grep -qx "polap-oatk"; then
        echo "ERROR: Conda environment 'polap-oatk' already exists."
      else
        conda create -y --name polap-oatk bioconda::oatk
        conda activate polap-oatk
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
}

install-conda_genus_species() {
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
}

install-mbg_genus_species() {
  if [[ "${opt_y_flag}" == false ]]; then
    read -p "Do you want to install mbg in the polap-mbg conda environment? (y/N): " confirm
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
      echo "You're in the base environment. Creating 'polap-mbg'..."
      if conda env list | awk '{print $1}' | grep -qx "polap-mbg"; then
        echo "ERROR: Conda environment 'polap-mbg' already exists."
      else
        conda create -y --name polap-mbg mbg
      fi
    else
      echo "Error: You're in the '$CONDA_DEFAULT_ENV' environment. Please activate base before running this script."
      exit 1
    fi
  else
    echo "mbg installation is canceled."
    echo "Check: https://github.com/maickrau/MBG"
    echo "conda install -c bioconda mbg"
    echo "Execute: MBG"
  fi
}

polap-analysis-oga_genus_species_for() {
  local _brg_outdir="${1}"
  local _brg_inum="${2:-0}"
  local _brg_adir="${3:-.}"
  local _brg_plastid="${4:-off}"
  local _brg_title="polalp-analysis-oga"

  local target_index="${_brg_outdir}-${_brg_inum}"
  local _brg_outdir_i="${_brg_outdir}/${_brg_adir}/${_brg_inum}"
  local long_sra="${_long["$target_index"]}"
  if [[ -z "$long_sra" ]]; then
    echo "Error: skipping ${_brg_outdir}-${_brg_inum} because it is not in the CSV."
    return
  fi

  mkdir -p "${_brg_outdir_i}"
  local _memlog_file="${_brg_outdir_i}/memlog-${_brg_title}.csv"
  local _summary_file="${_brg_outdir_i}/summary-${_brg_title}.txt"

  echo "[INFO] Starting polap analysis pipeline on ${_brg_outdir}-${_brg_inum} at dir:${_brg_outdir}/${_brg_adir}"

  # Start memory logger
  _polap_lib_process-start_memtracker "${_memlog_file}" \
    "${_polap_var_memtracker_time_interval}"

  polap-analysis-assemble2_genus_species "${_brg_outdir}" "${_brg_inum}" "${_brg_adir}"
  local _brg_inum_backup="${_brg_inum}"
  # polap-analysis-msbwt_genus_species "${_brg_outdir}" "${_brg_adir}" 0
  # echo before:polap-analysis-polish_genus_species "${_brg_outdir}" "${_brg_inum}" "${_brg_adir}"
  # _brg_inum="${_brg_inum_backup}"
  # echo after:polap-analysis-polish_genus_species "${_brg_outdir}" "${_brg_inum}" "${_brg_adir}"
  polap-analysis-polish_genus_species "${_brg_outdir}" "${_brg_inum}" "${_brg_adir}"
  # polap-clean_genus_species "${_brg_outdir}" 0
  # polap-clean_genus_species "${_brg_outdir}" "${_brg_inum}"

  # Summarize results after job (with previously defined summary function)
  _polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"

  echo "[INFO] End polap analysis pipeline on ${_brg_outdir}-${_brg_inum} at dir:${_brg_outdir}/${_brg_adir}"
}

polap-analysis-oga_genus_species() {
  local _brg_outdir="${1:-all}"
  local _brg_inum="${2:-0}"
  local _brg_adir="${3:-.}"
  local _brg_plastid="${4:-off}"

  if [[ "${_brg_outdir}" == "all" ]]; then
    for _v1 in "${Sall[@]}"; do
      polap-analysis-oga_genus_species_for "${_v1}" "${@:2}"
    done
  else
    polap-analysis-oga_genus_species_for "$@"
  fi
}

polap-analysis-wga_genus_species_for() {
  local _brg_outdir="${1}"
  local _brg_inum="${2:-0}"
  local _brg_adir="${3:-.}"
  local _brg_plastid="${4:-off}"
  local _brg_test="${5:-off}"
  local _brg_title="polalp-analysis-wga"

  local _brg_outdir="${_brg_outdir%/}"
  local target_index="${_brg_outdir}-${_brg_inum}"
  local _brg_outdir_i="${_brg_outdir}/${_brg_adir}/${_brg_inum}"
  local long_sra="${_long["$target_index"]}"
  if [[ -z "$long_sra" ]]; then
    echo "Error: skipping ${_brg_outdir}-${_brg_inum} because it is not in the CSV."
    return
  fi

  mkdir -p "${_brg_outdir_i}"
  local _memlog_file="${_brg_outdir_i}/memlog-${_brg_title}.csv"
  local _summary_file="${_brg_outdir_i}/summary-${_brg_title}.txt"

  echo "[INFO] Starting polap analysis WGA pipeline on ${_brg_outdir}-${_brg_inum} at dir:${_brg_outdir}/${_brg_adir}"

  # Start memory logger
  _polap_lib_process-start_memtracker "${_memlog_file}" \
    "${_polap_var_memtracker_time_interval}"

  polap-analysis-data_genus_species "${_brg_outdir}"
  polap-analysis-reduce_genus_species "${_brg_outdir}" "${_brg_inum}" "${_brg_adir}" 0
  polap-analysis-assemble1_genus_species "${_brg_outdir}" "${_brg_inum}" "${_brg_adir}" "${_brg_plastid}"

  # Summarize results after job (with previously defined summary function)
  _polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"

  echo "[INFO] End polap analysis WGA pipeline on ${_brg_outdir}-${_brg_inum} at dir:${_brg_outdir}/${_brg_adir}"
}

polap-analysis-wga_genus_species() {
  local _brg_outdir="${1:-all}"
  local _brg_inum="${2:-0}"
  local _brg_adir="${3:-.}"
  local _brg_plastid="${4:-off}"
  local _brg_test="${5:-off}"

  local _brg_outdir="${_brg_outdir%/}"

  if [[ "${_brg_outdir}" == "all" ]]; then
    for _v1 in "${Sall[@]}"; do
      polap-analysis-wga_genus_species_for "${_v1}" "${@:2}"
    done
  else
    polap-analysis-wga_genus_species_for "$@"
  fi
}

# duplicates
# polap-analysis-wga_genus_species() {
#
wga_genus_species() {
  local _brg_outdir="$1"
  local _brg_inum="${2:-0}"

  local target_index="${_brg_outdir}-0"

  local species_name="$(echo $1 | sed 's/_/ /')"
  local long_sra="${_long["$target_index"]}"
  local short_sra="${_short["$target_index"]}"

  rm -rf "${_brg_outdir}/0"

  ${_polap_cmd} assemble1 \
    -o ${_brg_outdir} \
    -l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
    --stopafter data

  rm -rf "${_brg_outdir}/0"

  command time -v ${_polap_cmd} flye1 polishing \
    -o ${_brg_outdir} \
    -l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
    2>${_brg_outdir}/timing-flye1.txt
}

# extract the archive at the remote
recover_genus_species() {
  local _brg_outdir="$1"

  # Directory without a slash
  _brg_outdir="${_brg_outdir%/}"

  if [[ -s "${_brg_outdir}-a.tar.gz" ]]; then
    _log_echo "Deleting ${_brg_outdir} ..."
    rm -rf "${_brg_outdir}"
    tar -zxf "${_brg_outdir}-a.tar.gz"
    mv "${_brg_outdir}-a" "${_brg_outdir}"
    _log_echo "we have recovered ${_brg_outdir}"
  else
    mkdir -p "${_brg_outdir}"
  fi
}

generate_nextdenovo_cfg() {
  local genome_size="$1"
  local input_file="${2:-input.fofn}"
  local work_dir="${3:-01_rundir}"
  local output_file="${4:-nextdenovo.cfg}"

  local _ncpu="$(cat /proc/cpuinfo | grep -c processor)"
  local mem_gb=$(free -g | awk '/^Mem:/{print $2}')
  local parallel_jobs=2
  local pa_correction=2
  local _nthreads=4
  if ((mem_gb > 260)); then
    parallel_jobs=8
    _nthreads=$((_ncpu / parallel_jobs - 1))
  elif ((mem_gb > 130)); then
    parallel_jobs=6
    _nthreads=$((_ncpu / parallel_jobs - 1))
  elif ((mem_gb > 70)); then
    parallel_jobs=3
    _nthreads=$((_ncpu / parallel_jobs - 1))
  fi
  local pa_correction=$parallel_jobs

  cat >"$output_file" <<EOF
[General]
job_type = local
job_prefix = nextDenovo
task = correct
rewrite = yes
deltmp = yes
parallel_jobs = ${parallel_jobs}
input_type = raw
read_type = ont
input_fofn = ${input_file}
workdir = ${work_dir}

[correct_option]
read_cutoff = 1k
genome_size = ${genome_size}
pa_correction = ${pa_correction}
sort_options = -m 20g -t ${_nthreads}
correction_options = -p ${_nthreads} -b
minimap2_options_raw = -t ${_nthreads}

[assemble_option]
minimap2_options_cns = -t ${_nthreads}
nextgraph_options = -a 1
EOF
}

polap-analysis-without-wga_genus_species_for() {
  local _brg_outdir="${1}"
  local _brg_inum="${2:-0}"
  local _brg_adir="${3:-.}"
  local _brg_plastid="${4:-off}"
  local _brg_title="polap-analysis-without-wga"

  local target_index="${_brg_outdir}-${_brg_inum}"
  local _brg_outdir_i="${_brg_outdir}/${_brg_adir}/${_brg_inum}"
  local _brg_outdir_source="${_brg_outdir}/${_brg_adir}_"
  local _brg_outdir_i_source="${_brg_outdir}/${_brg_adir}_/${_brg_inum}"
  local _outdir="${_brg_outdir}"-"${_brg_inum}"

  mkdir -p "${_brg_outdir_i}"
  local _memlog_file="${_brg_outdir_i}/memlog-${_brg_title}.csv"
  local _summary_file="${_brg_outdir_i}/summary-${_brg_title}.txt"

  if [[ "${_local_host}" != "$(hostname)" ]]; then
    mkdir -p "${_brg_outdir_source}"
    rsync -azuq ${_local_host}:"${PWD}"/"${_brg_outdir_source}"/ "${_brg_outdir_source}"/
  fi

  echo "[INFO] Starting polap analysis pipeline without WGA on ${_brg_outdir}-${_brg_inum} at dir:${_brg_outdir}/${_brg_adir}"

  # Start memory logger
  _polap_lib_process-start_memtracker "${_memlog_file}" \
    "${_polap_var_memtracker_time_interval}"

  polap-analysis-data_genus_species "${_brg_outdir}"
  polap-analysis-reduce_genus_species "${_brg_outdir}" "${_brg_inum}" "${_brg_adir}"

  # alternative to assemble1
  mkdir -p "${_outdir}"/0/30-contigger
  rsync -azuq "${_brg_outdir_i_source}"/30-contigger/ "${_outdir}"/0/30-contigger/
  cp -p "${_brg_outdir_i_source}"/mt.contig.name-? "${_outdir}"/0/
  cp -p "${_brg_outdir_i_source}"/*-polap-assemble1.txt "${_outdir}"/

  polap-analysis-assemble2_genus_species "${_brg_outdir}" "${_brg_inum}" "${_brg_adir}"
  local _brg_inum_backup="${_brg_inum}"
  polap-analysis-msbwt_genus_species "${_brg_outdir}" "${_brg_adir}" 0
  echo before:polap-analysis-polish_genus_species "${_brg_outdir}" "${_brg_inum}" "${_brg_adir}"
  _brg_inum="${_brg_inum_backup}"
  echo after:polap-analysis-polish_genus_species "${_brg_outdir}" "${_brg_inum}" "${_brg_adir}"
  polap-analysis-polish_genus_species "${_brg_outdir}" "${_brg_inum}" "${_brg_adir}"

  # polap-clean_genus_species "${_brg_outdir}" 0
  # polap-clean_genus_species "${_brg_outdir}" "${_brg_inum}"

  # Summarize results after job (with previously defined summary function)
  _polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"

  echo "[INFO] End polap analysis pipeline without WGA on ${_brg_outdir}-${_brg_inum} at dir:${_brg_outdir}/${_brg_adir}"

  if [[ "${_local_host}" != "$(hostname)" ]]; then
    rsync -azuq "${_brg_outdir_source}"/ ${_local_host}:"${PWD}"/"${_brg_outdir_source}"/
  fi
}

polap-analysis-without-wga_genus_species() {
  local _brg_outdir="${1:-all}"
  local _brg_inum="${2:-0}"
  local _brg_adir="${3:-.}"
  local _brg_plastid="${4:-off}"

  if [[ "${_brg_outdir}" == "all" ]]; then
    for _v1 in "${Sall[@]}"; do
      polap-analysis-without-wga_genus_species_for "${_v1}" "${@:2}"
    done
  else
    polap-analysis-without-wga_genus_species_for "$@"
  fi
}

polap-analysis_genus_species_for() {
  local _brg_outdir="${1}"
  local _brg_inum="${2:-0}"
  local _brg_adir="${3:-.}"
  local _brg_plastid="${4:-off}"
  local _brg_title="polalp-analysis"

  local target_index="${_brg_outdir}-${_brg_inum}"
  local _brg_outdir_i="${_brg_outdir}/${_brg_adir}/${_brg_inum}"
  local long_sra="${_long["$target_index"]}"
  if [[ -z "$long_sra" ]]; then
    echo "Error: skipping ${_brg_outdir}-${_brg_inum} because it is not in the CSV."
    return
  fi

  mkdir -p "${_brg_outdir_i}"
  local _memlog_file="${_brg_outdir_i}/memlog-${_brg_title}.csv"
  local _summary_file="${_brg_outdir_i}/summary-${_brg_title}.txt"

  echo "[INFO] Starting polap analysis pipeline on ${_brg_outdir}-${_brg_inum} at dir:${_brg_outdir}/${_brg_adir}"

  # Start memory logger
  _polap_lib_process-start_memtracker "${_memlog_file}" \
    "${_polap_var_memtracker_time_interval}"

  polap-analysis-data_genus_species "${_brg_outdir}"
  polap-analysis-reduce_genus_species "${_brg_outdir}" "${_brg_inum}" "${_brg_adir}"
  polap-analysis-assemble1_genus_species "${_brg_outdir}" "${_brg_inum}" "${_brg_adir}" "${_brg_plastid}"
  polap-analysis-assemble2_genus_species "${_brg_outdir}" "${_brg_inum}" "${_brg_adir}"
  local _brg_inum_backup="${_brg_inum}"
  polap-analysis-msbwt_genus_species "${_brg_outdir}" "${_brg_adir}" 0
  echo before:polap-analysis-polish_genus_species "${_brg_outdir}" "${_brg_inum}" "${_brg_adir}"
  _brg_inum="${_brg_inum_backup}"
  echo after:polap-analysis-polish_genus_species "${_brg_outdir}" "${_brg_inum}" "${_brg_adir}"
  polap-analysis-polish_genus_species "${_brg_outdir}" "${_brg_inum}" "${_brg_adir}"
  polap-clean_genus_species "${_brg_outdir}" 0
  polap-clean_genus_species "${_brg_outdir}" "${_brg_inum}"

  # Summarize results after job (with previously defined summary function)
  _polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"

  echo "[INFO] End polap analysis pipeline on ${_brg_outdir}-${_brg_inum} at dir:${_brg_outdir}/${_brg_adir}"
}

polap-analysis_genus_species() {
  local _brg_outdir="${1:-all}"
  local _brg_inum="${2:-0}"
  local _brg_adir="${3:-.}"
  local _brg_plastid="${4:-off}"

  if [[ "${_brg_outdir}" == "all" ]]; then
    for _v1 in "${Sall[@]}"; do
      polap-analysis_genus_species_for "${_v1}" "${@:2}"
    done
  else
    polap-analysis_genus_species_for "$@"
  fi
}

polap-analysis-data_genus_species_for() {
  local _brg_outdir="${1}"

  # Directory without a slash
  _brg_outdir="${_brg_outdir%/}"

  data-long_genus_species "${_brg_outdir}"
  data-short_genus_species "${_brg_outdir}"
}

polap-analysis-data_genus_species() {
  local _brg_outdir="${1:-all}"

  if [[ "${_brg_outdir}" == "all" ]]; then
    for _v1 in "${Sall[@]}"; do
      polap-analysis-data_genus_species_for "${_v1}" "${@:2}"
    done
  elif [[ "${_brg_outdir}" == "each" ]]; then
    for _v1 in "${Sall[@]}"; do
      polap-analysis-data_genus_species_for "${_v1}" "${@:2}"
    done
  else
    polap-analysis-data_genus_species_for "$@"
  fi
}

polap-analysis-polish_genus_species() {
  local _brg_outdir="${1}"
  local _brg_inum="${2:-0}"
  local _brg_adir="${3:-.}"
  local _brg_title="polap-analysis-polish"

  local target_index="${_brg_outdir}-${_brg_inum}"

  # local _outdir_data1="${_brg_outdir}"-"${_brg_inum}"-"polap-analysis-reduce"
  # local _outdir_assemble1="${_brg_outdir}"-"${_brg_inum}"-"polap-analysis-assemble1"
  # local _outdir_assemble2="${_brg_outdir}"-"${_brg_inum}"-"polap-analysis-assemble2"

  local _outdir="${_brg_outdir}"-"${_brg_inum}"
  local _brg_outdir_i="${_brg_outdir}/${_brg_adir}/${_brg_inum}"
  local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

  _polap_lib_conda-ensure_conda_env polap-fmlrc || exit 1
  # source "$(conda info --base)/etc/profile.d/conda.sh"
  # if [[ "${CONDA_DEFAULT_ENV:-}" != "polap-fmlrc" ]]; then
  # 	echo "[INFO] Activating conda environment 'polap-fmlrc'..."
  # 	conda activate polap-fmlrc
  # fi
  #
  # if [[ "${CONDA_DEFAULT_ENV:-}" != "polap-fmlrc" ]]; then
  # 	echo "[ERROR] Failed to enter conda environment 'polap-fmlrc'"
  # 	return
  # fi

  local _timing_txt="${_brg_outdir_i}/timing-${_brg_title}.txt"
  local _stdout_txt="${_brg_outdir_i}/stdout-${_brg_title}.txt"
  mkdir -p "${_brg_outdir_i}"
  local _memlog_file="${_brg_outdir_i}/memlog-${_brg_title}.csv"
  local _summary_file="${_brg_outdir_i}/summary-${_brg_title}.txt"
  # local _msbwt_basedir="${_brg_outdir}/${_brg_adir}"
  rm -f "${_timing_txt}"
  rm -f "${_stdout_txt}"

  echo "[INFO] Starting polishing pipeline on ${_brg_outdir}-${_brg_inum}"

  # Start memory logger
  _polap_lib_process-start_memtracker "${_memlog_file}" \
    "${_polap_var_memtracker_time_interval}"

  local _inum=0
  for j in {1..9}; do
    _brg_jnum=$((_inum + j))
    local _brg_outdir_j="${_brg_outdir}/${_brg_adir}/${_brg_inum}/${_brg_jnum}"
    if [[ -s "${_brg_outdir_j}/assembly.fasta" ]]; then
      command time -v "${_polap_cmd}" polish \
        -o "${_brg_outdir}/${_brg_adir}" \
        -p "${_brg_outdir_j}/assembly.fasta" \
        -f "${_brg_outdir_j}/assembly-polished.fa" \
        >>"${_stdout_txt}" \
        2>>"${_timing_txt}"
    else
      echo "Error: no such file: ${_brg_outdir_j}/assembly.fasta"
    fi
  done

  # Record the computer system info
  echo "hostname: $(hostname)" >>"${_timing_txt}"
  free -h >>"${_timing_txt}"
  lscpu >>"${_timing_txt}"

  # Summarize results after job (with previously defined summary function)
  _polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"

  conda deactivate
}

polap-analysis-msbwt_genus_species() {
  local _brg_outdir="${1}"
  local _brg_adir="${2:-.}"
  local _brg_coverage="${3:-0}"
  local _brg_title="polap-analysis-msbwt"
  local _brg_inum=0
  # local _brg_coverage="${3}"

  local target_index="${_brg_outdir}-${_brg_inum}"
  local _brg_outdir_i="${_brg_outdir}/${_brg_inum}"
  local long_sra="${_long["$target_index"]}"
  if [[ -z "$long_sra" ]]; then
    echo "Error: skipping ${_brg_outdir}-${_brg_inum} because it is not in the CSV."
    return
  fi
  local short_sra="${_short["$target_index"]}"
  local _outdir="${_brg_outdir}"-"${_brg_coverage}"
  local _brg_outdir_i="${_brg_outdir}/${_brg_adir}/${_brg_inum}"
  local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

  # Check msbwt
  local _msbwt="${_brg_outdir}/${_brg_adir}"/msbwt/comp_msbwt.npy
  local filesize
  if [[ -s "${_msbwt}" ]]; then
    # Get the file size in bytes
    filesize=$(stat --format=%s "${_msbwt}")
  else
    filesize=0
  fi

  if ((filesize > 1024)); then
    # Check if the file size is greater than 100 KB (100 * 1024 bytes)
    echo "[INFO] skipping the short-read polishing preparation."
    return
  fi

  _polap_lib_conda-ensure_conda_env polap-fmlrc || exit 1
  # source "$(conda info --base)/etc/profile.d/conda.sh"
  # if [[ "${CONDA_DEFAULT_ENV:-}" != "polap-fmlrc" ]]; then
  # 	echo "[INFO] Activating conda environment 'polap-fmlrc'..."
  # 	conda activate polap-fmlrc
  # fi
  #
  # if [[ "${CONDA_DEFAULT_ENV:-}" != "polap-fmlrc" ]]; then
  # 	echo "[ERROR] Failed to enter conda environment 'polap-fmlrc'"
  # 	return
  # fi

  local _timing_txt="${_brg_outdir_i}/timing-${_brg_title}.txt"
  local _stdout_txt="${_brg_outdir_i}/stdout-${_brg_title}.txt"
  mkdir -p "${_brg_outdir_i}"
  local _memlog_file="${_brg_outdir_i}/memlog-${_brg_title}.csv"
  local _summary_file="${_brg_outdir_i}/summary-${_brg_title}.txt"
  # mkdir -p "${_brg_outdir}/msbwts"
  mkdir -p "${_brg_outdir}/tmp"

  echo "[INFO] Starting FMLRC msbwt preparation pipeline on ${_brg_outdir}-${_brg_inum}"

  # subsample the short-read data
  if [[ "${_brg_coverage}" == "0" ]]; then
    # rm "${_brg_outdir}/tmp/${short_sra}_1.fastq"
    # rm "${_brg_outdir}/tmp/${short_sra}_2.fastq"
    ln -sf "$(pwd)/${short_sra}_1.fastq" \
      "${_brg_outdir}/tmp/"
    ln -sf "$(pwd)/${short_sra}_2.fastq" \
      "${_brg_outdir}/tmp/"
  else
    echo "Error: I do not know which genome size estimate is used."
    return
    local genomesize=$(<"${_outdir_assemble1}/short_expected_genome_size.txt")
    echo data-downsample-short_genus_species "${_brg_outdir}" \
      "${_brg_inum}" "${k}" "${genomesize}"
    data-downsample-short_genus_species "${_brg_outdir}" \
      "${_brg_inum}" "${k}" "${genomesize}"
  fi

  # Start memory logger
  _polap_lib_process-start_memtracker "${_memlog_file}" \
    "${_polap_var_memtracker_time_interval}"

  # msbwt
  command time -v ${_polap_cmd} prepare-polishing \
    -a "${_brg_outdir}/tmp/${short_sra}_1.fastq" \
    -b "${_brg_outdir}/tmp/${short_sra}_2.fastq" \
    -o ${_outdir} \
    >"${_stdout_txt}" \
    2>"${_timing_txt}"

  # Summarize results after job (with previously defined summary function)
  _polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"

  echo "hostname: $(hostname)" >>"${_timing_txt}"
  free -h >>"${_timing_txt}"
  lscpu >>"${_timing_txt}"

  # Record the computer system info
  conda deactivate

  rsync -azuq "${_outdir}"/msbwt/ "${_brg_outdir}/${_brg_adir}"/msbwt/
}

polap-analysis-assemble2_genus_species() {
  local _brg_outdir="${1}"
  local _brg_inum="${2:-0}"
  local _brg_adir="${3:-.}"
  local _brg_title="polap-analysis-assemble2"

  local target_index="${_brg_outdir}-${_brg_inum}"

  local _outdir_data1="${_brg_outdir}"-"${_brg_inum}"-"polap-analysis-reduce"
  # Check the key exists in the CSV config.
  local long_sra="${_long["$target_index"]}"
  if [[ -z "$long_sra" ]]; then
    echo "Error: skipping ${_brg_outdir}-${_brg_inum} because it is not in the CSV."
    return
  fi
  local _outdir="${_brg_outdir}"-"${_brg_inum}"
  local _brg_outdir_i="${_brg_outdir}/${_brg_adir}/${_brg_inum}"
  local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

  # Prepare the input data files
  rsync -azuq "${_brg_outdir_i}"/ "${_outdir}"/

  _polap_lib_conda-ensure_conda_env polap || exit 1

  local _timing_txt="${_brg_outdir_i}/timing-${_brg_title}.txt"
  local _stdout_txt="${_brg_outdir_i}/stdout-${_brg_title}.txt"
  mkdir -p "${_brg_outdir_i}"
  local _memlog_file="${_brg_outdir_i}/memlog-${_brg_title}.csv"
  local _summary_file="${_brg_outdir_i}/summary-${_brg_title}.txt"
  rm -f "${_timing_txt}"
  rm -f "${_stdout_txt}"

  echo "[INFO] Starting an organelle-genome assembly pipeline on ${_brg_outdir}-${_brg_inum}"

  # Start memory logger
  _polap_lib_process-start_memtracker "${_memlog_file}" \
    "${_polap_var_memtracker_time_interval}"

  local _inum=0
  for j in {1..9}; do
    _brg_jnum=$((_inum + j))
    if [[ -s "${_outdir}/${_inum}/mt.contig.name-${_brg_jnum}" ]]; then
      command time -v "${_polap_cmd}" assemble2 \
        -o "${_outdir}" \
        -i "${_inum}" \
        -j "${_brg_jnum}" \
        >>"${_stdout_txt}" \
        2>>"${_timing_txt}"
    fi
  done

  # Summarize results after job (with previously defined summary function)
  _polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"

  conda deactivate

  _polap_lib_timing-get_system_info >>"${_timing_txt}"

  rsync -azuq --max-size=100M "${_outdir}"/ "${_brg_outdir_i}"/
  # rm -rf "${_outdir_assemble2}"
}

polap-analysis-assemble1_genus_species() {
  local _brg_outdir="${1}"
  local _brg_inum="${2:-0}"
  local _brg_adir="${3:-.}"
  local _brg_plastid="${4:-off}"
  local _brg_test="${5}"
  local _brg_title="polap-analysis-assemble1"

  local target_index="${_brg_outdir}-${_brg_inum}"
  # Check the key exists in the CSV config.
  local long_sra="${_long["$target_index"]}"
  if [[ -z "$long_sra" ]]; then
    echo "Error: skipping ${_brg_outdir}-${_brg_inum} because it is not in the CSV."
    return
  fi
  local _outdir="${_brg_outdir}"-"${_brg_inum}"
  local _brg_outdir_i="${_brg_outdir}/${_brg_adir}/${_brg_inum}"
  local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

  # Prepare the input data files
  rsync -azuq "${_brg_outdir_i}"/ "${_outdir}"/

  _polap_lib_conda-ensure_conda_env polap || exit 1
  # source "$(conda info --base)/etc/profile.d/conda.sh"
  # if [[ "${CONDA_DEFAULT_ENV:-}" != "polap" ]]; then
  # 	echo "[INFO] Activating conda environment 'polap'..."
  # 	conda activate polap
  # fi
  #
  # if [[ "${CONDA_DEFAULT_ENV:-}" != "polap" ]]; then
  # 	echo "[ERROR] Failed to enter conda environment 'polap'"
  # 	return
  # fi

  local _timing_txt="${_brg_outdir_i}/timing-${_brg_title}.txt"
  local _stdout_txt="${_brg_outdir_i}/stdout-${_brg_title}.txt"
  mkdir -p "${_brg_outdir_i}"
  local _memlog_file="${_brg_outdir_i}/memlog-${_brg_title}.csv"
  local _summary_file="${_brg_outdir_i}/summary-${_brg_title}.txt"

  echo "[INFO] Starting a whole-genome assembly pipeline on ${_brg_outdir}-${_brg_inum}"

  # Start memory logger
  _polap_lib_process-start_memtracker "${_memlog_file}" \
    "${_polap_var_memtracker_time_interval}"

  command time -v "${_polap_cmd}" flye1 \
    -o "${_outdir}" \
    >"${_stdout_txt}" \
    2>"${_timing_txt}"

  command time -v "${_polap_cmd}" edges-stats \
    -o "${_outdir}" \
    >>"${_stdout_txt}" \
    2>>"${_timing_txt}"

  command time -v "${_polap_cmd}" annotate \
    -o "${_outdir}" \
    >>"${_stdout_txt}" \
    2>>"${_timing_txt}"

  if [[ "${_brg_plastid}" == "off" ]]; then
    command time -v "${_polap_cmd}" seeds \
      -o "${_outdir}" \
      >>"${_stdout_txt}" \
      2>>"${_timing_txt}"
  else
    command time -v "${_polap_cmd}" seeds \
      -o "${_outdir}" \
      --plastid \
      >>"${_stdout_txt}" \
      2>>"${_timing_txt}"
  fi

  # Record the computer system info
  echo "hostname: $(hostname)" >>"${_timing_txt}"
  free -h >>"${_timing_txt}"
  lscpu >>"${_timing_txt}"

  # Summarize results after job (with previously defined summary function)
  _polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"

  conda deactivate

  rsync -azuq "${_outdir}"/ "${_brg_outdir_i}"/
  # rm -rf "${_outdir_assemble1}"
}

# workdir is not copied back to the outdir_i.
polap-analysis-reduce_genus_species() {
  local _brg_outdir="${1}"
  local _brg_inum="${2:-0}"
  local _brg_adir="${3:-.}"

  local _brg_coverage="${4:-0}"
  local _brg_test="${5}"
  local _brg_title="polap-analysis-reduce"

  local target_index="${_brg_outdir}-${_brg_inum}"
  # Check the key exists in the CSV config.
  local long_sra="${_long["$target_index"]}"
  if [[ -z "$long_sra" ]]; then
    echo "Error: skipping ${_brg_outdir}-${_brg_inum} because it is not in the CSV."
    return
  fi
  local short_sra="${_short["$target_index"]}"
  if [[ "${_brg_coverage}" == "0" ]]; then
    _brg_coverage="${_downsample["$target_index"]}"
  fi
  local _brg_outdir_i="${_brg_outdir}/${_brg_adir}/${_brg_inum}"
  local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

  _polap_lib_conda-ensure_conda_env polap || exit 1

  local _timing_txt="${_brg_outdir_i}/timing-${_brg_title}.txt"
  local _stdout_txt="${_brg_outdir_i}/stdout-${_brg_title}.txt"
  # local _outdir_data1="${_brg_outdir}"-"${_brg_inum}"-"${_brg_title}"
  local _outdir="${_brg_outdir}"-"${_brg_inum}"
  mkdir -p "${_brg_outdir_i}"

  local _memlog_file="${_brg_outdir_i}/memlog-${_brg_title}.csv"
  local _summary_file="${_brg_outdir_i}/summary-${_brg_title}.txt"
  local _log_file="${_brg_outdir_i}/log-${_brg_title}.txt"

  echo "[INFO] Starting data reduction pipeline on ${_brg_outdir}-${_brg_inum} with ${_brg_coverage}x depth"

  # Start memory logger
  _polap_lib_process-start_memtracker "${_memlog_file}" \
    "${_polap_var_memtracker_time_interval}"

  # Run your command
  command time -v "${_polap_cmd}" assemble1 \
    -l "${long_sra}".fastq \
    -a "${short_sra}"_1.fastq \
    -b "${short_sra}"_2.fastq \
    --stopafter data \
    --redo \
    -c "${_brg_coverage}" \
    -o "${_outdir}" \
    >"${_stdout_txt}" \
    2>"${_timing_txt}"

  # Summarize results after job (with previously defined summary function)
  _polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"

  conda deactivate

  # Record the computer system info
  # echo "hostname: $(hostname)" >>"${_timing_txt}"
  # free -h >>"${_timing_txt}"
  # lscpu >>"${_timing_txt}"
  _polap_lib_timing-get_system_info >>"${_timing_txt}"

  # rsync -azuq --max-size=100M "${_outdir}"/ "${_brg_outdir_i}"/
  # rm -rf "${_outdir_data1}"
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
    if _polap_lib_file-is_at_least_1MB "${_brg_outdir_i}/cns.fa"; then
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
    if _polap_lib_file-is_at_least_1MB "${_brg_outdir_i}/ont.asm.ec.fq"; then
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
  echo Not tested: prepare-long-data_genus_species "${_brg_outdir}" "${_brg_inum}"
  exit 1
  prepare-long-data_genus_species "${_brg_outdir}" "${_brg_inum}"

  source "$(conda info --base)/etc/profile.d/conda.sh"
  if [[ "$CONDA_DEFAULT_ENV" != "oatk" ]]; then
    echo "You're not in the oatk environment. Chaniging 'oatk'..."
    conda activate polap-oatk
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

  if [[ "$_POLAP_DEBUG" == "1" ]]; then
    for key in "${!_long[@]}"; do
      echo "$key => ${_long[$key]}"
    done
  fi

  # check if long_sra is empty
  if [[ -z "$long_sra" ]]; then
    _log_echo "Error: no long SRA for ${_brg_outdir}-${_brg_inum}"
    return
  fi

  mkdir -p "${_brg_outdir}/${_brg_inum}"

  # Step 1. prepare input fastq data
  local long_data="${_media_dir}/${long_sra}.fastq.tar.gz"
  local long_fastq="${long_sra}.fastq"

  if [[ "${_local_host}" == "$(hostname)" ]]; then
    if [[ -s "${long_fastq}" ]]; then
      _log_echo "found: ${long_fastq}"
    else
      if [[ -s "${long_data}" ]]; then
        tar -zxf ${long_data}
      else
        ln -s "${_media1_dir}/${long_sra}.fastq"
      fi
    fi
  else
    if [[ -s "${long_sra}.fastq" ]]; then
      _log_echo "Found: ${long_fastq}"
    else
      if ssh ${_local_host} "test -f ${long_data}" 2>/dev/null; then
        if [[ -s "${long_sra}".fastq.tar.gz ]]; then
          _log_echo "Found: ${long_fastq}.fastq.tar.gz"
        else
          _log_echo "copying ${long_data} from ${_local_host}"
          scp ${_local_host}:${long_data} .
        fi
        _log_echo "decompressing ${long_sra}.fastq.tar.gz"
        tar -zxf "${long_sra}".fastq.tar.gz
        rm -f "${long_sra}".fastq.tar.gz
      elif ssh ${_local_host} "test -f ${_media1_dir}/${long_sra}.fastq" 2>/dev/null; then
        scp ${_local_host}:"${_media1_dir}/${long_sra}.fastq" .
      else
        _polap_lib_conda-ensure_conda_env polap || exit 1
        _log_echo "  downloading long-read SRA ID: ${long_sra} ... be patient!"
        bash "${_polap_script_bin_dir}"/polap-ncbitools fetch sra "$long_sra"
        conda deactivate
      fi
    fi
  fi
  if [[ "${_POLAP_DEBUG}" == "1" ]]; then
    _log_echo "[INFO] long-read: ${long_sra} is ready."
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
    _log_echo "Error: no short SRA for ${_brg_outdir}-${_brg_inum}"
    return
  fi

  mkdir -p "${_brg_outdir}/${_brg_inum}"

  # Step 1. prepare input fastq data
  local short_data="${_media_dir}/${short_sra}.fastq.tar.gz"
  local short_fastq="${short_sra}_1.fastq"

  if [[ "${_local_host}" == "$(hostname)" ]]; then
    if [[ -s "${short_fastq}" ]]; then
      _log_echo "found: ${short_fastq}"
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
      _log_echo "Found: ${short_fastq}"
    else
      if ssh ${_local_host} "test -f ${short_data}" 2>/dev/null; then
        if [[ -s "${short_sra}".fastq.tar.gz ]]; then
          _log_echo "Found: ${short_fastq}.fastq.tar.gz"
        else
          _log_echo "copying ${short_data} from ${_local_host}"
          scp ${_local_host}:${short_data} .
        fi
        _log_echo "decompressing ${short_sra}.fastq.tar.gz"
        tar -zxf "${short_sra}".fastq.tar.gz
        rm -f "${short_sra}".fastq.tar.gz
      elif ssh ${_local_host} "test -f ${_media1_dir}/${short_sra}_1.fastq" 2>/dev/null; then
        scp ${_local_host}:"${_media1_dir}/${short_sra}_1.fastq" .
        scp ${_local_host}:"${_media1_dir}/${short_sra}_2.fastq" .
      else
        _polap_lib_conda-ensure_conda_env polap || exit 1
        _log_echo "  downloading short-read SRA ID: ${short_sra} ... be patient!"
        bash "${_polap_script_bin_dir}"/polap-ncbitools fetch sra "$short_sra"
        conda deactivate
      fi
    fi
  fi
  if [[ "${_POLAP_DEBUG}" == "1" ]]; then
    _log_echo "[INFO] short-read: ${short_sra} is ready."
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
  local subcmd1="$1"
  local opt_y_flag="$2"
  local -n cmd_args_ref="$3" # name reference to array

  # set the global variables
  # _brg_outdir
  # _brg_sindex
  if [[ -v cmd_args_ref[0] ]]; then
    local _arg2=${cmd_args_ref[0]}
  else
    local _arg2=""
  fi
  if [[ -v cmd_args_ref[1] ]]; then
    local _arg3=${cmd_args_ref[1]}
  else
    local _arg3=""
  fi
  _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

  case "$subcmd1" in
  download-polap-github | \
    patch-polap | bleeding-edge-polap | local-edge-polap | \
    test-polap | \
    install-apptainer | \
    install-mitohifi | \
    data-long | data-short | data-peek | \
    data-downsample-long | data-downsample-short | \
    mitohifi | \
    list-subcommands)
    handled=1
    ;;
    ##### INSERT_COMMAND_HERE #####
  run-direct-dflye)
    handled=1
    ;;
  download-bioproject)
    handled=1
    ;;
  archive-run)
    handled=1
    ;;
  get-dflye)
    handled=1
    ;;
  config-add)
    handled=1
    ;;
  config)
    handled=1
    ;;
  config-view)
    handled=1
    ;;
  run-polap-assemble-wga)
    handled=1
    ;;
  run-polap-reduce-data)
    handled=1
    ;;
  run-polap-prepare-data)
    handled=1
    ;;
  run-direct-select-oga)
    handled=1
    ;;
  run-direct-bandage-oga)
    handled=1
    ;;
  run-direct-oga)
    handled=1
    ;;
  run-direct-wga)
    handled=1
    ;;
  run-direct-flye)
    handled=1
    ;;
  run-direct-read)
    handled=1
    ;;
  run-direct-map)
    handled=1
    ;;
  run-direct-seed)
    handled=1
    ;;
  run-direct)
    handled=1
    ;;
  install-latex)
    handled=1
    ;;
  download-species)
    handled=1
    ;;
  setup-completions)
    handled=1
    ;;
  download-test-data-cflye)
    handled=1
    ;;
  clean)
    handled=1
    ;;
  clean-cflye)
    handled=1
    ;;
  get-timing-pmat)
    handled=1
    ;;
  get-timing)
    handled=1
    ;;
  update-local)
    handled=1
    ;;
  update-github)
    handled=1
    ;;
  update)
    handled=1
    ;;
  get-cflye)
    handled=1
    ;;
  get)
    handled=1
    ;;
  print-help-all)
    handled=1
    ;;
  help)
    handled=1
    ;;
  setup-fmlrc2)
    handled=1
    ;;
  install-fmlrc2)
    handled=1
    ;;
  install-minimal)
    handled=1
    ;;
  install-all)
    handled=1
    ;;
  setup-nvim)
    handled=1
    ;;
  setup-pmat)
    handled=1
    ;;
  setup-polap)
    handled=1
    ;;
  setup-bandage)
    handled=1
    ;;
  install-nvim)
    handled=1
    ;;
  remove)
    handled=1
    ;;
  uninstall-fmlrc)
    handled=1
    ;;
  uninstall-dflye)
    handled=1
    ;;
  uninstall-cflye)
    handled=1
    ;;
  uninstall-oatk)
    handled=1
    ;;
  uninstall-tippo)
    handled=1
    ;;
  uninstall-pmat)
    handled=1
    ;;
  uninstall-getorganelle)
    handled=1
    ;;
  uninstall-polap)
    handled=1
    ;;
  uninstall)
    handled=1
    ;;
  run-summary-data)
    handled=1
    ;;
  install-man)
    handled=1
    ;;
  run-polap-disassemble-compare)
    handled=1
    ;;
  run-polap-disassemble-check)
    handled=1
    ;;
  run-polap-disassemble)
    handled=1
    ;;
  run-polap)
    handled=1
    ;;
  run-oatk-nextdenovo)
    handled=1
    ;;
  run-oatk-ont)
    handled=1
    ;;
  run-estimate-genomesize)
    handled=1
    ;;
  run-nextdenovo-polish)
    handled=1
    ;;
  run-pmat)
    handled=1
    ;;
  run-tippo)
    handled=1
    ;;
  run-oatk)
    handled=1
    ;;
  run-extract-ptdna-ptgaul)
    handled=1
    ;;
  run-polish-ptdna-ptgaul)
    handled=1
    ;;
  run-ptgaul)
    handled=1
    ;;
  download-ptdna)
    handled=1
    ;;
  run-msbwt)
    handled=1
    ;;
  run)
    handled=1
    ;;
  run-getorganelle)
    handled=1
    ;;
  install-pmat)
    handled=1
    ;;
  install-bandage)
    handled=1
    ;;
  download)
    handled=1
    ;;
  delete-polap-github)
    handled=1
    ;;
  delete)
    handled=1
    ;;
  install-getorganelle)
    handled=1
    ;;
  install-oatk)
    handled=1
    ;;
  install-tippo)
    handled=1
    ;;
  install-dflye)
    handled=1
    ;;
  install-cflye)
    handled=1
    ;;
  install-polap)
    handled=1
    ;;
  install-fmlrc)
    handled=1
    ;;
  setup-conda)
    handled=1
    ;;
  setup)
    handled=1
    ;;
  list)
    handled=1
    ;;
  recover)
    handled=1
    ;;
  install)
    handled=1
    ;;
  install-conda)
    handled=1
    ;;
  install-mbg)
    handled=1
    ;;
  polap-analysis-oga)
    handled=1
    ;;
  polap-analysis-wga)
    handled=1
    ;;
  nextdenovo-polish)
    handled=1
    ;;
  polap-analysis-without-wga)
    handled=1
    ;;
  polap-analysis)
    handled=1
    ;;
  polap-analysis-data)
    handled=1
    ;;
  polap-analysis-polish)
    handled=1
    ;;
  polap-analysis-msbwt)
    handled=1
    ;;
  polap-analysis-assemble2)
    handled=1
    ;;
  polap-analysis-assemble1)
    handled=1
    ;;
  polap-analysis-reduce)
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
  mkdir-all)
    handled=1
    ;;
  rm-empty)
    handled=1
    ;;
  rm)
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
  install-getorganelle)
    ${subcmd1}_genus_species
    ;;
  install-tippo)
    ${subcmd1}_genus_species
    ;;
  install-oatk)
    ${subcmd1}_genus_species
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
    ##### INSERT_CASE_HERE #####
  run-direct-dflye)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <outdir> [index:0|N] [inum:0|N] [jnum:1|N] [knum:0|N]"
      echo "  $(basename ${0}) ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  download-bioproject)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <bioproject accession ID>"
      echo "  $(basename ${0}) ${subcmd1} PRJNA990649"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  archive-run)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <outdir> [index:0|N]"
      echo "  $(basename ${0}) ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  get-dflye)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <outdir> <index:0|N> <hostname>"
      echo "  $(basename ${0}) ${subcmd1} Arabidopsis_thaliana hostname1"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  config-add)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <outdir> [inum:0|N]"
      echo "  $(basename ${0}) ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  config)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <view|add>"
      echo "  $(basename ${0}) ${subcmd1} view"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  config-view)
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  run-polap-assemble-wga)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <outdir> [index:0|N]"
      echo "  $(basename ${0}) ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  run-polap-reduce-data)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <outdir> [index:0|N]"
      echo "  $(basename ${0}) ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  run-polap-prepare-data)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <outdir> [index:0|N]"
      echo "  $(basename ${0}) ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  run-direct-select-oga)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <outdir> [inum:0|N] [jnum:0] [adir:t2]"
      echo "  $(basename ${0}) ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  run-direct-bandage-oga)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <outdir> [inum:0|N] [adir:t2] [bandage|off]"
      echo "  $(basename ${0}) ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  run-direct-oga)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <outdir> [index:0|N]"
      echo "  $(basename ${0}) ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  run-direct-wga)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <outdir> [index:0|N]"
      echo "  $(basename ${0}) ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  run-direct-flye)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <outdir> [index:0|N] [inum:0|N] [jnum:1|N]"
      echo "  $(basename ${0}) ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  run-direct-read)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <outdir> [inum:0|N] [jnum:1|N] [knum:2|N] [adir:t2]"
      echo "  $(basename ${0}) ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  run-direct-map)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <outdir> [inum:0|N] [jnum:1|N] [knum:2|N] [adir:t2]"
      echo "  $(basename ${0}) ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  run-direct-seed)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <outdir> [inum:0|N] [jnum:1|N] [knum:2|N] [adir:t2]"
      echo "  $(basename ${0}) ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  run-direct-dga)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <outdir> [inum:0|N]"
      echo "  $(basename ${0}) ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  download-species)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <outdir>"
      echo "  $(basename ${0}) ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    polap-analysis-data_genus_species "${cmd_args_ref[@]}"
    ;;
  setup-completions)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1}"
      echo "  $(basename ${0}) ${subcmd1} cflye"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  download-test-data-cflye)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1}"
      echo "  $(basename ${0}) ${subcmd1}"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  clean)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <which> <outdir>"
      echo "  $(basename ${0}) ${subcmd1} cflye Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  clean-cflye)
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  get-timing-pmat)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <outdir> [inum:0|N]"
      echo "  $(basename ${0}) ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  get-timing)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <tool> <outdir> [inum:0|N]"
      echo "  $(basename ${0}) ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  update-local)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <local_github_path>"
      echo "  git clone --quiet https://github.com/goshng/polap.git"
      echo "  $(basename ${0}) ${subcmd1} polap"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  update-github)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} [main|tag|hash]"
      echo "  $(basename ${0}) ${subcmd1} 0.4.3.7.8"
      echo "  $(basename ${0}) ${subcmd1} befbb51"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  update)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} TOOL"
      echo "  $(basename ${0}) ${subcmd1} github"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  get-cflye)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <outdir> <hostname>"
      echo "  $(basename ${0}) ${subcmd1} Arabidopsis_thaliana hostname1"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  get)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} TOOL..."
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  print-help-all)
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  help)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} [command]"
      echo "  $(basename ${0}) ${subcmd1} example"
      echo "  $(basename ${0}) ${subcmd1} install"
      _subcmd1_clean="help_message_${subcmd1//-/_}"
      declare -n ref="${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  setup-fmlrc2)
    ${subcmd1}_genus_species
    ;;
  install-latex)
    ${subcmd1}_genus_species
    ;;
  install-fmlrc2)
    ${subcmd1}_genus_species
    ;;
  install-minimal)
    ${subcmd1}_genus_species
    ;;
  install-all)
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  setup-nvim)
    ${subcmd1}_genus_species
    ;;
  setup-pmat)
    ${subcmd1}_genus_species
    ;;
  setup-polap)
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  setup-bandage)
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  remove)
    uninstall_genus_species "${cmd_args_ref[@]}"
    ;;
  uninstall-fmlrc2)
    ${subcmd1}_genus_species
    ;;
  uninstall-fmlrc)
    ${subcmd1}_genus_species
    ;;
  uninstall-dflye)
    ${subcmd1}_genus_species
    ;;
  uninstall-cflye)
    ${subcmd1}_genus_species
    ;;
  uninstall-oatk)
    ${subcmd1}_genus_species
    ;;
  uninstall-tippo)
    ${subcmd1}_genus_species
    ;;
  uninstall-pmat)
    ${subcmd1}_genus_species
    ;;
  uninstall-getorganelle)
    ${subcmd1}_genus_species
    ;;
  uninstall-polap)
    ${subcmd1}_genus_species
    ;;
  run-summary-data)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <outdir> [inum:0|N]"
      echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  install-man)
    ${subcmd1}_genus_species
    ;;
  run-polap-disassemble-compare)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <outdir> [inum:0|N]"
      echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  run-polap-disassemble-check)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <outdir> [inum:0|N]"
      echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  run-polap-disassemble)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <outdir> [inum:0|N] [polishing:default|simple] [random:off|on|random]"
      echo "  $(basename ${0}) ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  run-polap)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <outdir> [inum:0|N]"
      echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  run-oatk-nextdenovo)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <outdir> [inum:0|N]"
      echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  run-oatk-ont)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <outdir> [inum:0|N]"
      echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  run-estimate-genomesize)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <outdir> [inum:0|N]"
      echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  run-nextdenovo-polish)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <outdir> [inum:0|N]"
      echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  run-pmat)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <outdir> [inum:0|N]"
      echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  run-tippo)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <outdir> [inum:0|N]"
      echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  run-oatk)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <outdir> [inum:0|N]"
      echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  run-extract-ptdna-ptgaul)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <outdir>"
      echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  run-polish-ptdna-ptgaul)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <outdir>"
      echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  download-ptdna)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <outdir>"
      echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${_arg2}"
    ;;
  run-msbwt)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <outdir> [inum:0|N]"
      echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  run)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <command> ..."
      echo "  $(basename ${0}) ${subcmd1} getorganelle Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  run-getorganelle)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <outdir> [inum:0|N]"
      echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    # [[ "${_arg3}" == arg3 ]] && _arg3=""
    # ${subcmd1}_genus_species "${_arg2}" "${_arg3}"
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  run-ptgaul)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <outdir>"
      echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  install-pmat)
    if [[ "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1}"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species
    ;;
  install-bandage)
    ${subcmd1}_genus_species
    ;;
  download-test-data2)
    ${subcmd1}_genus_species
    ;;
  download)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} <data>"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  delete-polap-github)
    ${subcmd1}_genus_species
    ;;
  delete)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} TOOL..."
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  install-dflye)
    ${subcmd1}_genus_species
    ;;
  install-cflye)
    ${subcmd1}_genus_species
    ;;
  install-polap)
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  install-fmlrc)
    ${subcmd1}_genus_species
    ;;
  install-nvim)
    ${subcmd1}_genus_species
    ;;
  setup-conda)
    if [[ "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1}"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species
    ;;
  setup)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} TOOL..."
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  list)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} QUERY [any|start|end]"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  uninstall)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} TOOL..."
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    # [[ -z "${_arg3}" ]] && _arg3=""
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  install)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 || "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1} TOOL..."
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    # [[ -z "${_arg3}" ]] && _arg3=""
    ${subcmd1}_genus_species "${cmd_args_ref[@]}"
    ;;
  install-conda)
    ${subcmd1}_genus_species
    ;;
  install-mbg)
    ${subcmd1}_genus_species
    ;;
  polap-analysis-oga)
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
  polap-analysis-wga)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 ]]; then
      echo "Help: ${subcmd1} <outdir> [inum:0|N] [analysis:.|string] [plastid:off|on] [test:off|on|test]"
      echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    [[ "${_arg3}" == arg3 ]] && _arg3="0"
    [[ "${_arg4}" == arg4 ]] && _arg4="."
    [[ "${_arg5}" == arg5 ]] && _arg5="off"
    [[ "${_arg6}" == arg6 ]] && _arg6="off"
    ${subcmd1}_genus_species "${_arg2}" "${_arg3}" "${_arg4}" "${_arg5}" "${_arg6}"
    ;;
  nextdenovo-polish)
    if [[ "${_arg2}" == arg2 ]]; then
      echo "Help: ${subcmd1} <outdir> [inum:0|N] [analysis:.|string]"
      echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    [[ "${_arg3}" == arg3 ]] && _arg3="0"
    [[ "${_arg4}" == arg4 ]] && _arg4="."
    ${subcmd1}_genus_species "${_arg2}" "${_arg3}" "${_arg4}"
    ;;
  polap-analysis-without-wga)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 ]]; then
      echo "Help: ${subcmd1} <outdir> [inum:0|N] [analysis:.|string] [plastid:off|on]"
      echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    [[ "${_arg3}" == arg3 ]] && _arg3="0"
    [[ "${_arg4}" == arg4 ]] && _arg4="."
    [[ "${_arg5}" == arg5 ]] && _arg5="off"
    ${subcmd1}_genus_species "${_arg2}" "${_arg3}" "${_arg4}" "${_arg5}"
    ;;
  polap-analysis)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 ]]; then
      echo "Help: ${subcmd1} <outdir> [inum:0|N] [analysis:.|string] [plastid:off|on]"
      echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    [[ "${_arg3}" == arg3 ]] && _arg3="0"
    [[ "${_arg4}" == arg4 ]] && _arg4="."
    [[ "${_arg5}" == arg5 ]] && _arg5="off"
    ${subcmd1}_genus_species "${_arg2}" "${_arg3}" "${_arg4}" "${_arg5}"
    ;;
  polap-analysis-data)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 ]]; then
      echo "Help: ${subcmd1} <all|outdir>"
      echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
      echo "  ${0} ${subcmd1} all"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    ${subcmd1}_genus_species "${_arg2}"
    ;;
  polap-analysis-polish)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 ]]; then
      echo "Help: ${subcmd1} <outdir> [inum:0|N] [analysis:.|string]"
      echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    [[ "${_arg3}" == arg3 ]] && _arg3="0"
    [[ "${_arg4}" == arg4 ]] && _arg4="."
    ${subcmd1}_genus_species "${_arg2}" "${_arg3}" "${_arg4}"
    ;;
  polap-analysis-msbwt)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 ]]; then
      echo "Help: ${subcmd1} <outdir> [analysis:.|string] [coverage:0|N]"
      echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    [[ "${_arg3}" == arg3 ]] && _arg3="."
    [[ "${_arg4}" == arg4 ]] && _arg4="0"
    ${subcmd1}_genus_species "${_arg2}" "${_arg3}" "${_arg4}"
    ;;
  polap-analysis-assemble2)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 ]]; then
      echo "Help: ${subcmd1} <outdir> [inum:0|N] [analysis:.|string]"
      echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    [[ "${_arg3}" == arg3 ]] && _arg3="0"
    [[ "${_arg4}" == arg4 ]] && _arg4="."
    ${subcmd1}_genus_species "${_arg2}" "${_arg3}" "${_arg4}"
    ;;
  polap-analysis-assemble1)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 ]]; then
      echo "Help: ${subcmd1} <outdir> [inum:0|N] [analysis:.|string] [plastid:off|on] [test:off|on|test]"
      echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    [[ "${_arg3}" == arg3 ]] && _arg3="0"
    [[ "${_arg4}" == arg4 ]] && _arg4="."
    [[ "${_arg5}" == arg5 ]] && _arg5="off"
    ${subcmd1}_genus_species "${_arg2}" "${_arg3}" "${_arg4}" "${_arg5}"
    ;;
  polap-analysis-reduce)
    if [[ -z "${_arg2}" || "${_arg2}" == arg2 ]]; then
      echo "Help: ${subcmd1} <outdir> [inum:0|N] [analysis:.|t1|string] [coverage:0|N] [test:off|on|test]"
      echo "  ${0} ${subcmd1} Arabidopsis_thaliana"
      _subcmd1_clean="${subcmd1//-/_}"
      declare -n ref="help_message_${_subcmd1_clean}"
      echo "$ref"
      exit 0
    fi
    [[ "${_arg3}" == arg3 ]] && _arg3="0"
    [[ "${_arg4}" == arg4 ]] && _arg4="."
    [[ "${_arg5}" == arg5 ]] && _arg5="0"
    [[ "${_arg6}" == arg6 ]] && _arg6="off"
    ${subcmd1}_genus_species "${_arg2}" "${_arg3}" "${_arg4}" "${_arg5}" "${_arg6}"
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
  recover)
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
    if [[ "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1}"
      echo "  Delete all empty folders at the current folder."
      exit 0
    fi
    handled=1
    find . -type d -empty -delete
    echo "Deleting empty folders ... done."
    ;;
  rm)
    if [[ "${_arg2}" == "-h" || "${_arg2}" == "--help" ]]; then
      echo "Help: ${subcmd1}"
      echo "  Delete all species folders at the current folder."
      exit 0
    fi
    handled=1
    read -p "Do you really want to delete all species folders? (YES/NO): " confirm
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
      echo "Deleting all species folders is canceled."
      ;;
    esac
    ;;
  esac

  [[ $handled -eq 1 ]] && return 0 || return 1
}
