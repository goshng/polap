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
# Polap has its main script called polap.sh and its testing scripts including
# polap-data-v1.sh and polap-data-v2.sh. There are more such scripts with
# similar names. But, they are extensively tested yet although they are used
# for testing purpose.
# These testing scripts, e.g., polap-data-v2.sh, have their own subcommands.
# Common subcommands used by these testing scripts are placed in this script.
# Such common subcommands include installation and execution of other tools.
# See Also:
# A helper script named polap-data-menu.sh can be used to add a new subcommand
# for the testing scripts.
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
	_polap_var_memtracker_time_interval=120
	_polap_var_oatk_memtracker_time_interval=12
else
	_polap_var_memtracker_time_interval=60
	_polap_var_oatk_memtracker_time_interval=6
fi

source "${_POLAPLIB_DIR}/polap-lib-conda.sh"
source "${_POLAPLIB_DIR}/polap-lib-timing.sh"
source "${_POLAPLIB_DIR}/polap-lib-unit.sh"
source "${_POLAPLIB_DIR}/polap-lib-array.sh"
source "${_POLAPLIB_DIR}/polap-lib-number.sh"
source "${_POLAPLIB_DIR}/polap-lib-file.sh"
source "${_POLAPLIB_DIR}/polap-lib-process.sh"
source "${_POLAPLIB_DIR}/polap-lib-extract.sh"
source "${_POLAPLIB_DIR}/run-polap-function-utilities.sh"

: "${_brg_outdir:=.}"
_polap_data_cmd="$(basename "$0" .sh)"
_polap_data_csv="$(basename "$0" .sh).csv"
_polap_data_data="$(basename "$0" .sh).data"
_polap_data_txt="$(basename "$0" .sh).txt"
_polap_data_md="man"

# _log_echo() {
# 	if [[ -d "${_brg_outdir}" ]]; then
# 		echo "$(date '+%Y-%m-%d %H:%M:%S') [$subcmd1] - $1" >>"${_brg_outdir}/${_polap_data_txt}"
# 	else
# 		echo "$(date '+%Y-%m-%d %H:%M:%S') [$subcmd1] - $1" >>"./${_polap_data_txt}"
# 	fi
# 	echo "$(date '+%Y-%m-%d %H:%M:%S') [$subcmd1] - $1"
# }

# ------------------------------
# Core logger
# ------------------------------
# usage: __log <level> <message...>
__log() {
	local level="$1"
	shift
	local msg="$*"

	# If silent mode, suppress all logs
	# (( _silent == 1 )) && return 0

	# Print only if verbosity is high enough
	local func="${FUNCNAME[2]}"
	local file="$(basename "${BASH_SOURCE[2]}")"
	local line="${BASH_LINENO[1]}"
	local tag="[$func@$file:$line]"
	# local log_line="$(date '+%Y-%m-%d %H:%M:%S') $tag [$subcmd1] - $msg"
	local log_line="$(date '+%H:%M:%S') $tag [$subcmd1] - $msg"

	if [[ -d "${_brg_outdir}" ]]; then
		echo "$log_line" >>"${_brg_outdir}/${_polap_data_txt}"
	else
		echo "$log_line" >>"./${_polap_data_txt}"
	fi

	if ((_brg_verbose > level)); then
		echo "$log_line"
	fi
}

# ------------------------------
# Convenience wrappers
# ------------------------------
_log_echo0() { __log 0 "$@"; } # Always (unless -q)
_log_echo1() { __log 1 "$@"; } # Needs -v
_log_echo2() { __log 2 "$@"; } # Needs -v -v

# Default _log_echo is the most verbose
_log_echo() { _log_echo2 "$@"; }

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
declare -a _polap_tippo_ont_options
_polap_tippo_ont_options=(
	onthq
)
# ont

declare -a _polap_tippo_hifi_options
_polap_tippo_hifi_options=(
	hifi
)

declare -a _polap_tippo_clr_options
_polap_tippo_clr_options=(
	clr
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
help_message_bolap=$(
	cat <<HEREDOC
bolap help
HEREDOC
)

help_message_do_mv=$(
	cat <<HEREDOC

  menu title
HEREDOC
)

help_message_run_polap_disassemble_hifi=$(
	cat <<HEREDOC

  menu title
HEREDOC
)

help_message_get_conf_long_sra_id=$(
	cat <<HEREDOC

  menu title
HEREDOC
)

help_message_run_polap_annotate_check=$(
	cat <<HEREDOC

  Compare ptDNA and mtDNA from annotating and selecting long-reads
  with ptGAUL's results.
HEREDOC
)

x_help_message_run_polap_syncassemble=$(
	cat <<HEREDOC

  Assemble organelle genomes by annotating and selecting long-reads.

  This is the main run function for assembling plant organelle genomes
  using ONT reads.
HEREDOC
)

help_message_demo=$(
	cat <<EOF
Name:
  bolap demo - execute bolap's demonstration

Synopsis:
  bolap demo

Description:
  polap's demo shows how one can use different bolap/polap commands.

Examples:
  read:
    bolap -c a.csv demo read

  syncassemble:
    bolap demo syncasm

  syncassemble:
    bolap polap-syncassemble Brassica_napus 1

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
)

help_message_simulate=$(
	cat <<EOF
Name:
  bolap simulate - execute polap's syncassemble

Synopsis:
  bolap simulate <SPECIES>

Description:
  This simulates ONT or HiFi sequencing data.

Examples:
  Execute polap :
    bolap download-ptdna Actinidia_arguta
    bolap download-mtdna Actinidia_arguta
    bolap simulate Actinidia_arguta

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
)

help_message_run_polap_syncassemble=$(
	cat <<EOF
Name:
  bolap run-polap-syncassemble - execute polap's syncassemble

Synopsis:
  bolap run polap-syncassemble <SPECIES>

Description:
  polap's syncassemble filters out nuclear-derived reads.

Examples:
  Execute polap readassemble:
    bolap run polap-syncassemble Carex_pseudochinensis

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
)

help_message_run_readassemble=$(
	cat <<HEREDOC

  Assemble organelle genomes by annotating and selecting long-reads.

  This is the main run function for assembling plant organelle genomes
  using ONT reads.
HEREDOC
)

help_message_run_readassemble_pt=$(
	cat <<HEREDOC

  Assemble ptDNA by annotating and selecting long-reads.
HEREDOC
)

help_message_run_readassemble_nt=$(
	cat <<HEREDOC

  Assemble mtDNA by annotating and selecting long-reads.

  This is the main run function for assembling plant mitochondrial genome
  using ONT reads.
HEREDOC
)

help_message_run_readassemble_mt=$(
	cat <<HEREDOC

  Assemble mtDNA by annotating and selecting long-reads.

  This is the main run function for assembling plant mitochondrial genome
  using ONT reads.
HEREDOC
)

help_message_run_polap_readassemble=$(
	cat <<HEREDOC

  use polap-readassemble command.

  Assemble organelle genomes by annotating and selecting long-reads.

  This is the main run function for assembling plant organelle genomes
  using ONT reads.
HEREDOC
)

help_message_run_polap_readassemble_pt=$(
	cat <<HEREDOC

  use polap-readassemble command.

  Assemble ptDNA by annotating and selecting long-reads.
HEREDOC
)

help_message_run_polap_readassemble_mt=$(
	cat <<HEREDOC

  use polap-readassemble command.

  Assemble mtDNA by annotating and selecting long-reads.

  This is the main run function for assembling plant mitochondrial genome
  using ONT reads.
HEREDOC
)

help_message_sync=$(
	cat <<HEREDOC

  menu title
HEREDOC
)

help_message_run_readassemble_ont_pt=$(
	cat <<HEREDOC

  Assemble reads selected by plastid gene annotation.

  output: o/annotate-read/mt.1.fa
  output: o/annotate-read/mt.1.gfa
  output: o/annotate-read/pt.1.fa
  output: o/annotate-read/pt.1.gfa
HEREDOC
)

help_message_run_readassemble_ont_mt=$(
	cat <<HEREDOC

  Assemble reads selected by mitochondrial gene annotation.

  output: o/annotate-read/mt.1.fa
  output: o/annotate-read/mt.1.gfa
  output: o/annotate-read/pt.1.fa
  output: o/annotate-read/pt.1.gfa
HEREDOC
)

help_message_readassemble_annotate_pt=$(
	cat <<HEREDOC

  Annotate ONT reads with organelle genes.

  output: o/annotate-read/pt
HEREDOC
)

help_message_readassemble_annotate_mt=$(
	cat <<HEREDOC

  Annotate ONT reads with organelle genes.

  output: o/annotate-read/mt
HEREDOC
)

help_message_readassemble_annotate_nt=$(
	cat <<HEREDOC

  Annotate ONT reads with organelle genes.

  output: o/annotate-read/mt
HEREDOC
)

help_message_run_pca_oga=$(
	cat <<HEREDOC

  menu title
HEREDOC
)

help_message_run_pca_wga=$(
	cat <<HEREDOC

  menu title
HEREDOC
)

help_message_run_pca_pt=$(
	cat <<HEREDOC

  menu title
HEREDOC
)

help_message_run_pca_ont_kmer=$(
	cat <<HEREDOC

  We apply PCA and UMAP to ONT data.
  Downsample the data before the application.
HEREDOC
)

help_message_run_pca_ont_pt=$(
	cat <<HEREDOC

  We use flye v2.9.6 without genome size or assembly coverage option.
HEREDOC
)

help_message_run_pca_kmer=$(
	cat <<HEREDOC

  menu title
HEREDOC
)

help_message_run_pca_kmer_short=$(
	cat <<HEREDOC

  menu title
HEREDOC
)

help_message_run_pca_count_filter=$(
	cat <<HEREDOC

  menu title
HEREDOC
)

help_message_run_pca_count=$(
	cat <<HEREDOC

  menu title
HEREDOC
)

help_message_run_hifi_oga=$(
	cat <<HEREDOC

  menu title
HEREDOC
)

help_message_sample_fastq=$(
	cat <<HEREDOC

  menu title
HEREDOC
)

help_message_get_aflye=$(
	cat <<HEREDOC

  menu title
HEREDOC
)

help_message_polap_cmd_find_genome_size=$(
	cat <<HEREDOC

  polap command title
HEREDOC
)

help_message_polap_cmd_xxx=$(
	cat <<HEREDOC

  polap command title
HEREDOC
)

help_message_polap_cmd_total_length_long=$(
	cat <<HEREDOC

  polap command title
HEREDOC
)

help_message_polap_cmd=$(
	cat <<HEREDOC

  menu title
HEREDOC
)

help_message_polap_cmd_clean_menus=$(
	cat <<HEREDOC

  menu title
HEREDOC
)

help_message_update_dev=$(
	cat <<HEREDOC

  cd .. && rsync -aq thorne:$PWD/github/ github/ && cd -
HEREDOC
)

help_message_polap_cmd_init=$(
	cat <<HEREDOC

  menu title
HEREDOC
)

help_message_download_sra=$(
	cat <<HEREDOC

  Download SRA data.
HEREDOC
)

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

  view long,short chilo
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

help_message_man=$(
	cat <<HEREDOC

  These are implemented in each polap-data script.

  man init
  man update
  man table
  man figure
  man pdf
  man final
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
Name:
  bolap run tipo - execute TIPPo v2.4

Synopsis:
  bolap run tippo <SPECIES>

Description:
  polap's syncassemble filters out nuclear-derived reads.

  We use TIPPo v2.4.

Examples:
  Execute polap readassemble:
    bolap run tippo Carex_pseudochinensis

  Execute TIPPo v2.4 on the ONT reads.
    TIPPo.v2.4.pl -f <long_sra>.fastq -p ont -m 1500

  Execute TIPPo v2.4 on the nextDenovo polished reads.
    TIPPo.v2.4.pl -f cns.fa

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
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

  Subsample the long-read data and save the original data file under o/tmp.
  We want to see the effect on long-read subsampling on the final result.
  Because we save the original data file at o/tmp, we can keep calling
  this command to regenerate a down-sampled data file.
  If the coverage is greater than the input data, the down-sampled data
  file is the same as the original data.
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

#
# END: HELP
################################################################################

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
do-mv_genus_species_for() {
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"
	local _brg_type="${3:-nt}"

	local _brg_inum=0
	local _brg_adir _brg_title _brg_target _brg_rundir _brg_outdir_i
	local _timing_txt _stdout_txt _memlog_file _summary_file

	brg_common_setup \
		_brg_outdir _brg_sindex _brg_adir _brg_title \
		_brg_target _brg_rundir _brg_outdir_i \
		_timing_txt _stdout_txt _memlog_file _summary_file

	local platform="${_platform["$_brg_target"]}"
	if [[ "$platform" != "ONT" ]]; then
		mkdir -p "${_brg_outdir_i}/estimate-genomesize"
		if [[ -s "${_brg_outdir}/short_expected_genome_size.txt" ]]; then
			cp -p "${_brg_outdir}/short_expected_genome_size.txt" \
				"${_brg_outdir_i}/estimate-genomesize"
		else
			echo "no genome size: ${_brg_outdir}"
		fi
	fi

	# debug local variables
	# for debugging: Inline local printing local var
	# while IFS= read -r line; do
	# 	if [[ $line =~ ^declare\ --\ ([^=]+)= ]]; then
	# 		var="${BASH_REMATCH[1]}"
	# 		printf "%s=%q\n" "$var" "${!var}"
	# 	fi
	# done < <(local -p 2>/dev/null)

}

do-mv_genus_species() {
	local _brg_outdir="${1:-some}"

	if [[ "${_brg_outdir}" == "all" ]]; then
		for _v1 in "${Sall[@]}"; do
			do-mv_genus_species_for "${_v1}" "${@:2}"
		done
	elif [[ "${_brg_outdir}" == "some" ]]; then
		for _v1 in "${Ssome[@]}"; do
			do-mv_genus_species_for "${_v1}" "${@:2}"
		done
	elif [[ "${_brg_outdir}" == "each" ]]; then
		for _v1 in "${Skeys[@]}"; do
			do-mv_genus_species_for "${_v1}" "${@:2}"
		done
	else
		do-mv_genus_species_for "$@"
	fi
}

run-polap-disassemble-hifi_genus_species() {
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"
	local _brg_type="${3:-nt}"

	local _brg_inum=0
	local _brg_adir _brg_title _brg_target _brg_rundir _brg_outdir_i
	local _timing_txt _stdout_txt _memlog_file _summary_file

	brg_common_setup \
		_brg_outdir _brg_sindex _brg_adir _brg_title \
		_brg_target _brg_rundir _brg_outdir_i \
		_timing_txt _stdout_txt _memlog_file _summary_file

	# Extra folders
	local _brg_outdir_t="${_brg_outdir}/${opt_t_arg}"
	local _brg_titledir="${_brg_outdir_i}/${_brg_title}-${_brg_type}"
	local _brg_runtitledir="${_brg_rundir}-${_brg_title}-${_brg_type}"

	if [[ -v _long["$_brg_target"] ]]; then
		local long_sra="${_long["$_brg_target"]}"
	else
		echo "Error: ${_brg_target} because it is not in the CSV."
		return
	fi

	local _brg_coverage="50m"

	local platform="${_platform["$_brg_target"]}"
	local option_data_type="--nano-raw"

	if [[ "${platform}" == "PACBIO_SMRT" ]]; then
		option_data_type="--pacbio-hifi"
		if [[ "${_brg_type}" == "pt" ]]; then
			lib-polap-readassemble-hifi-pt-brg-coverage
		else
			lib-polap-readassemble-hifi-brg-coverage
		fi
	elif [[ "${platform}" == "ONT" ]]; then
		lib-polap-readassemble-brg-coverage
	fi

	_log_echo "Do: disassemble ${_brg_type}DNA using seed reads of coverage (${_brg_coverage}) with ${option_data_type}: ${_brg_rundir}"

	if [[ ! -s "${long_sra}".fastq ]]; then
		data-long_genus_species "${_brg_outdir}"
	fi

	rm -rf "${_brg_rundir}"
	if [[ "${platform}" == "PACBIO_SMRT" ]]; then
		if [[ ! -s "${_brg_outdir}/short_expected_genome_size.txt" ]]; then
			write-genomesize-10x "${_brg_outdir}"
		fi
		mkdir -p "${_brg_rundir}"
		cp "${_brg_outdir}/short_expected_genome_size.txt" "${_brg_rundir}"
	fi

	data-downsample-long_genus_species "${_brg_outdir}" "${_brg_sindex}" "${_brg_coverage}"

	# local resolved_fastq=$(<"${_brg_outdir}/${_brg_inum}/l.fastq.path.txt")
	local resolved_fastq="${long_sra}.fastq"
	# echo "  - FASTQ: $resolved_fastq"

	# Start memory logger
	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"

	_polap_lib_conda-ensure_conda_env polap || exit 1

	if [[ "${_brg_type}" == "pt" ]]; then
		_log_echo "Not Tested Yet!"
		return
		${_polap_cmd} disassemble \
			"${option_data_type}" \
			--plastid \
			-l "${resolved_fastq}" \
			-o "${_brg_rundir}"
	elif [[ "${_brg_type}" == "animal" ]]; then
		_log_echo "Not Tested Yet!"
		return
		${_polap_cmd} disassemble \
			"${option_data_type}" \
			--animal \
			-l "${resolved_fastq}" \
			-o "${_brg_rundir}"
	else
		# elif [[ "${_brg_type}" == "mt" ]]; then
		if [[ "${_local_host}" != "$(hostname)" ]]; then
			mkdir -p "${_brg_rundir}"
			scp "${_local_host}:$PWD/${_brg_rundir}/pt.0.gfa" "${_brg_rundir}/"
			# sync_genus_species "${_brg_outdir}" "${_brg_sindex}" --pull
		else
			cp "$HOME/all/polap/hifi1/${_brg_rundir}/pt.0.gfa" "${_brg_rundir}/"
		fi
		# else

		# assemble ptDNA first
		if [[ ! -s "${_brg_rundir}/pt.0.gfa" ]]; then
			${_polap_cmd} readassemble \
				"${option_data_type}" \
				--plastid \
				-l "${resolved_fastq}" \
				-o "${_brg_rundir}"
		fi

		if [[ -s "${_brg_rundir}/pt.0.gfa" ]]; then

			# filter by ptDNA reference
			if [[ ! -s "${_brg_rundir}/kmer/ref-filtered.fastq" ]]; then
				_log_echo "filter by ptDNA reference"
				${_polap_cmd} filter reference hifi \
					-l "${resolved_fastq}" \
					--reference "${_brg_rundir}/pt.0.gfa" \
					-o "${_brg_rundir}"
			fi
			local p=$(<"${_brg_rundir}/kmer/chloroplast-ratio.txt")
			_log_echo "chloroplast-raito: $p"

			# plant mt
			# --disassemble-a 100m \
			# --disassemble-b 1500m \
			# --disassemble-m 1.5m \
			local n=10
			${_polap_cmd} disassemble-hifi \
				"${option_data_type}" \
				-l "${_brg_rundir}/kmer/ref-filtered.fastq" \
				--disassemble-p $p \
				--disassemble-m 1m \
				--disassemble-memory 100 \
				--disassemble-n $n \
				--disassemble-stop-after stage1 \
				-o "${_brg_rundir}"

			# check the last index or select one index for seed contigs
			n=$(find "${_brg_rundir}/0/disassemble/1/1" -maxdepth 1 -type d -printf "%f\n" |
				grep -E '^[0-9]+$' |
				sort -n |
				tail -n 1)

			n=$((n + 1))

			# create seed contigs
			${_polap_cmd} seed-mito \
				-o "${_brg_rundir}/0/disassemble/1/1" \
				-i $((n - 1)) -j $n

			for ((i = 1; i < 5; i++)); do
				cp "${_brg_rundir}/0/disassemble/1/1/$((n - 1))/mt.contig.name-$n" \
					"${_brg_rundir}/0/disassemble/1/1/$((n - 1))/mt.contig.name-$((n + i))"
			done

			${_polap_cmd} assemble-rate \
				"${option_data_type}" \
				-o "${_brg_rundir}/0/disassemble/1/1" \
				-l "${resolved_fastq}" \
				-i $((n - 1)) -j $n -w 3000
			ln -s "0/disassemble/1/1/$n/assembly_graph.gfa" \
				"${_brg_rundir}/nt.2.gfa"

			${_polap_cmd} assemble-omega \
				"${option_data_type}" \
				-o "${_brg_rundir}/0/disassemble/1/1" \
				-l "${resolved_fastq}" \
				-i $((n - 1)) -j $((n + 1))
			ln -s "0/disassemble/1/1/$((n + 1))/assembly_graph.gfa" \
				"${_brg_rundir}/nt.3.gfa"

			# use the filtered data only
			resolved_fastq="${_brg_rundir}/kmer/ref-filtered.fastq"

			${_polap_cmd} assemble-rate \
				"${option_data_type}" \
				-o "${_brg_rundir}/0/disassemble/1/1" \
				-l "${resolved_fastq}" \
				-i $((n - 1)) -j $((n + 2)) -w 3000
			ln -s "0/disassemble/1/1/$((n + 2))/assembly_graph.gfa" \
				"${_brg_rundir}/nt.4.gfa"

			${_polap_cmd} assemble-omega \
				"${option_data_type}" \
				-o "${_brg_rundir}/0/disassemble/1/1" \
				-l "${resolved_fastq}" \
				-i $((n - 1)) -j $((n + 3))
			ln -s "0/disassemble/1/1/$((n + 3))/assembly_graph.gfa" \
				"${_brg_rundir}/nt.5.gfa"

		else
			_log_echo "No ptDNA to filter!"
		fi
	fi

	conda deactivate

	# Summarize results after job (with previously defined summary function)
	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"
	_polap_lib_timing-get_system_info >>"${_timing_txt}"

	# Save some results
	rsync -azuq --max-size=5M \
		"${_brg_rundir}/annotate-read-${_brg_type}"/ \
		"${_brg_titledir}"/

	if [[ "${_local_host}" != "$(hostname)" ]]; then
		sync_genus_species "${_brg_outdir}" "${_brg_sindex}" --main-push
		if [[ "${_brg_type}" == "nt" ]]; then
			rm -rf "${_brg_rundir}" "${_brg_outdir}"
			rm -f "${long_sra}.fastq" "${long_sra}-10x.fastq.tar.gz"
		fi
	fi
}

run-polap-downsample_genus_species() {
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"
	source "${_POLAPLIB_DIR}/polap-variables-data.sh"

	if [[ "$#" -gt 1 ]]; then
		shift 2
	else
		shift
	fi

	# defaults
	local _brg_type="miniasm"
	local _brg_omega="1500"
	local _brg_downsample="no-downsample"
	local _brg_coverage="10g"

	# parse options
	while [[ $# -gt 0 ]]; do
		case "$1" in
		--type)
			_brg_type="$2"
			shift 2
			;;
		-w)
			_brg_omega="$2"
			shift 2
			;;
		-f)
			shift
			;;
		-d)
			_brg_downsample="downsample"
			shift
			;;
		--coverage)
			_brg_coverage="$2"
			shift 2
			;;
		--dry-run)
			_brg_dry="on"
			shift
			;;
		--)
			shift
			break
			;;
		*)
			_log_echo "Unknown option: $1"
			return 1
			;;
		esac
	done

	# redefine
	local _brg_rundir="${_brg_outdir_i}/${_brg_title}-1-${_brg_type}"

	mkdir -p "${_brg_outdir_i}"

	local platform="${_platform["$_brg_target"]}"

	_log_echo1 "Do: downsample ${_brg_type}DNA using seed reads of coverage (${_brg_coverage}): ${_brg_rundir}"

	# rm -rf "${_brg_target}"
	mkdir -p "${_brg_target}/tmp"

	# Start memory logger
	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"

	_polap_lib_conda-ensure_conda_env polap || exit 1

	# check out input data
	# ln -sf "${_brg_tmpdir}/l.fq" "${long_sra}.fastq"
	# _brg_input_data="${long_sra}.fastq"
	if [[ ! -s "${long_sra}".fastq ]]; then
		data-long_genus_species "${_brg_outdir}"
	fi

	# plant mt
	_log_echo1 ${_polap_cmd} fastq-sample-to \
		-g "${_brg_coverage}" \
		-l "${long_sra}.fastq" \
		--outfile "${_brg_target}/tmp/l.fq"
	${_polap_cmd} fastq-sample-to \
		-g "${_brg_coverage}" \
		-l "${long_sra}.fastq" \
		--outfile "${_brg_target}/tmp/l.fq"

	_log_echo1 mv "${_brg_target}/tmp/l.fq" "${long_sra}.fastq"
	if ! mv "${_brg_target}/tmp/l.fq" "${long_sra}.fastq" 2>/dev/null; then
		_log_echo1 "skip: mv failed because the two files are the same"
	fi

	_log_echo1 seqkit stats -Ta "${long_sra}.fastq" -o "${long_sra}.fastq.seqkit.stats.Ta.txt"
	seqkit stats -Ta "${long_sra}.fastq" -o "${long_sra}.fastq.seqkit.stats.Ta.txt"
	cut -f1-5 "${long_sra}.fastq.seqkit.stats.Ta.txt"
	conda deactivate

	# Summarize results after job (with previously defined summary function)
	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"
	_polap_lib_timing-get_system_info >>"${_timing_txt}"

	if [[ -d "${_brg_target}" ]]; then
		_log_echo1 rm -rf "${_brg_target}"
		rm -rf "${_brg_target}"
	fi

}

run-polap-aflye_genus_species() {
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"
	source "${_POLAPLIB_DIR}/polap-variables-data.sh"

	if [[ "$#" -gt 1 ]]; then
		shift 2
	else
		shift
	fi

	# defaults
	local _brg_type="aflye"
	local _brg_omega="1500"
	local _brg_downsample="no-downsample"
	local _brg_redo="off"
	local _brg_cleanup="on"

	# parse options
	while [[ $# -gt 0 ]]; do
		case "$1" in
		--type)
			_brg_type="$2"
			shift 2
			;;
		-w)
			_brg_omega="$2"
			shift 2
			;;
		-f)
			shift
			;;
		--redo)
			_brg_redo="on"
			shift
			;;
		--no-cleanup)
			_brg_cleanup="off"
			shift
			;;
		-d)
			_brg_downsample="downsample"
			shift
			;;
		--dry-run)
			_brg_dry="on"
			shift
			;;
		--)
			shift
			break
			;;
		*)
			_log_echo "Unknown option: $1"
			return 1
			;;
		esac
	done

	# redefine
	local _brg_rundir="${_brg_outdir_i}/${_brg_title}-1-${_brg_type}"

	local mtn="${_mtn["$_brg_target"]}"
	local ptn="${_ptn["$_brg_target"]}"
	local short_sra="${_short["$_brg_target"]}"

	# _log_echo0 "mtn: ${mtn}"
	# _log_echo0 "ptn: ${ptn}"
	# return 0

	_log_echo0 "_brg_redo: ${_brg_redo}"
	_log_echo0 "-s ${_brg_rundir}/mt.1.gfa"
	if [[ "${_brg_redo}" == "off" && -s "${_brg_rundir}/mt.1.gfa" ]]; then
		_log_echo0 "We have already done the assembly. Use --redo if you redo it."
		return 0
	fi

	mkdir -p "${_brg_outdir_i}"

	local platform="${_platform["$_brg_target"]}"

	local _brg_coverage="50m"

	local option_data_type="--nano-raw"
	if [[ "${platform}" == "PACBIO_SMRT" ]]; then
		option_data_type="--pacbio-hifi"
		if [[ "${_brg_type}" == "pt" ]]; then
			lib-polap-readassemble-hifi-pt-brg-coverage
		else
			lib-polap-readassemble-hifi-brg-coverage
		fi
	elif [[ "${platform}" == "PACBIO_CLR" ]]; then
		option_data_type="--pacbio-raw"
		if [[ "${_brg_type}" == "pt" ]]; then
			lib-polap-readassemble-clr-pt-brg-coverage
		else
			lib-polap-readassemble-clr-brg-coverage
		fi
	elif [[ "${platform}" == "ONT" ]]; then
		if [[ "${_brg_type}" == "pt" ]]; then
			lib-polap-readassemble-pt-brg-coverage
		else
			lib-polap-readassemble-brg-coverage
		fi
	fi

	_log_echo1 "Do: assemble pt/mt DNA using seed reads of coverage (${_brg_coverage}) with ${option_data_type}: ${_brg_rundir}"

	rm -rf "${_brg_target}"
	mkdir -p "${_brg_target}"

	# Start memory logger
	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"

	_polap_lib_conda-ensure_conda_env polap || exit 1

	# check out input data
	# ln -sf "${_brg_tmpdir}/l.fq" "${long_sra}.fastq"
	# _brg_input_data="${long_sra}.fastq"
	if [[ ! -s "${long_sra}".fastq ]]; then
		data-long_genus_species "${_brg_outdir}"
	fi

	if [[ ! -s "${short_sra}_1.fastq" ]]; then
		data-short_genus_species "${_brg_outdir}"
	fi

	local resolved_fastq="${long_sra}.fastq"

	_log_echo2 "plant mt $_brg_verbose_str"

	${_polap_cmd} assemble1 \
		"${option_data_type}" \
		-l "${resolved_fastq}" \
		-a "${short_sra}_1.fastq" \
		-b "${short_sra}_2.fastq" \
		${_brg_verbose_str} \
		-o "${_brg_target}"

	_log_echo0 ${_polap_cmd} assemble2 \
		"${option_data_type}" \
		-l "${resolved_fastq}" \
		-a "${short_sra}_1.fastq" \
		-b "${short_sra}_2.fastq" \
		${_brg_verbose_str} \
		-o "${_brg_target}"

	conda deactivate

	# Summarize results after job (with previously defined summary function)
	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"
	_polap_lib_timing-get_system_info >>"${_timing_txt}"

	# Save some results
	rsync -azuq --max-size=5M \
		"${_brg_target}/" "${_brg_rundir}/"

	# if [[ -s "${long_sra}".fastq ]]; then
	# 	_log_echo1 rm -f "${long_sra}".fastq
	# 	rm -f "${long_sra}".fastq
	# fi
	if [[ "${_brg_cleanup}" == "on" ]]; then
		if [[ -d "${_brg_target}" ]]; then
			_log_echo1 rm -rf "${_brg_target}"
			rm -rf "${_brg_target}"
		fi
	fi
}

# two positional arguments: outdir and index
# optiotns
# type: -t pt|mt|nt
# omega: -w 1500
# downsample: -d
run-polap-readassemble_genus_species() {
	local _brg_outdir="${1:-$_brg_outdir}"
	local _brg_sindex="${2:-$_brg_sindex}"

	local bolap_cmd="${FUNCNAME%%_*}"

	help_message=$(
		cat <<EOF
Name:
  bolap - $bolap_cmd

Synopsis:
  bolap $bolap_cmd

Description:
  bolap

Examples:
  Execute polap readassemble:
    bolap readassemble -s Vigna_radiata

  Execute polap disassemble:
    bolap disassemble -s Vigna_radiata

  Execute polap syncassemble:
    bolap syncassemble -s Vigna_radiata

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	# defaults
	local _brg_type="miniasm"
	local _brg_omega="1500"
	local _brg_downsample="no-downsample"
	local _brg_redo="off"
	local _brg_cleanup="on"

	parse_commandline() {
		set -- "${_brg_unknown_opts[@]}"

		# source "${_POLAPLIB_DIR}/polap-cmd-version.sh" # '.' means 'source'
		while test $# -gt 0; do
			_key="$1"
			case "$_key" in
			--type)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_type="$2"
					shift || true
				fi
				;;
			-w)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_omega="$2"
					shift || true
				fi
				;;
			--redo)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_redo="on"
					shift || true
				fi
				;;
			--no-cleanup)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_cleanup="off"
					shift || true
				fi
				;;
			-d)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_downsample="downsample"
					shift || true
				fi
				;;
			--dry-run)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_dry="on"
					shift || true
				fi
				;;
			-*)
				_log_echo0 "[INFO] no such options: $1"
				;;
			*)
				break
				;;
			esac
			shift || true
		done
	}

	declare -n ref="help_message"
	if [[ "${_brg_help}" == "on" ]]; then
		local manfile=$(_bolap_lib_man-convert_help_message "$ref" "${bolap_cmd}")
		man "$manfile"
		rm -f "$manfile"
		return
	fi

	source "${_POLAPLIB_DIR}/polap-variables-data.sh"
	local long_sra="${_long["$_brg_target"]}"

	parse_commandline

	# redefine
	# local _brg_rundir="${_brg_outdir_i}/${_brg_title}-1-${_brg_type}"

	_log_echo0 "target: $_brg_target"
	# local mtn="${_mtn["$_brg_target"]}"
	# local ptn="${_ptn["$_brg_target"]}"
	local mtn=3
	local ptn=2

	# _log_echo0 "mtn: ${mtn}"
	# _log_echo0 "ptn: ${ptn}"
	# return 0

	_log_echo0 "_brg_redo: ${_brg_redo}"
	_log_echo0 "-s ${_brg_rundir}/mt.1.gfa"
	if [[ "${_brg_redo}" == "off" && -s "${_brg_rundir}/mt.1.gfa" ]]; then
		_log_echo0 "We have already done the assembly. Use --redo if you redo it."
		return 0
	fi

	mkdir -p "${_brg_outdir_i}"

	local platform="${_platform["$_brg_target"]}"

	local _brg_coverage="50m"

	local option_data_type="--nano-raw"
	# [[ "${platform}" == "ONT" ]]
	# pt: 10g
	# mt: 10g or more
	# if [[ "${_brg_type}" == "pt" ]]; then
	# 	lib-polap-readassemble-pt-brg-coverage
	# else
	# 	lib-polap-readassemble-brg-coverage
	# fi

	_log_echo1 "Asemble organelle genome sequences with ${option_data_type}"

	# Always redo
	rm -rf "${_brg_rundir}"
	mkdir -p "${_brg_rundir}"

	# Always redo
	rm -rf "${_brg_target}"
	mkdir -p "${_brg_target}"

	# Start memory logger
	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"

	_polap_lib_conda-ensure_conda_env polap || exit 1

	# check out input data
	# ln -sf "${_brg_tmpdir}/l.fq" "${long_sra}.fastq"
	# _brg_input_data="${long_sra}.fastq"
	if [[ ! -s "${long_sra}".fastq ]]; then
		data-long_genus_species "${_brg_outdir}"
	fi

	local resolved_fastq="${long_sra}.fastq"

	# if [[ "${_brg_type}" == "pt" ]]; then
	# 	${_polap_cmd} readassemble \
	# 		"${option_data_type}" \
	# 		--plastid \
	# 		--no-reduction-reads \
	# 		-w "${_brg_omega}" \
	# 		--downsample "${_brg_coverage}" \
	# 		-l "${resolved_fastq}" \
	# 		-o "${_brg_target}"
	# elif [[ "${_brg_type}" == "animal" ]]; then
	# 	${_polap_cmd} readassemble \
	# 		"${option_data_type}" \
	# 		--animal \
	# 		--downsample "${_brg_coverage}" \
	# 		-l "${resolved_fastq}" \
	# 		-o "${_brg_target}"

	# elif [[ "${_brg_type}" == "miniasm" ]]; then
	_log_echo2 "plant mt $_brg_verbose_str"

	_log_echo0 ${_polap_cmd} readassemble \
		"${option_data_type}" \
		-l "${resolved_fastq}" \
		--readassemble-mtn "${mtn}" \
		${_brg_verbose_str} \
		-o "${_brg_target}"

	${_polap_cmd} readassemble \
		"${option_data_type}" \
		-l "${resolved_fastq}" \
		--readassemble-mtn "${mtn}" \
		${_brg_verbose_str} \
		-o "${_brg_target}"

	# elif [[ "${_brg_type}" == "nt" ]]; then
	# 	# plant mt
	# 	${_polap_cmd} readassemble \
	# 		"${option_data_type}" \
	# 		--downsample "${_brg_coverage}" \
	# 		-l "${resolved_fastq}" \
	# 		-o "${_brg_target}"
	# elif [[ "${_brg_type}" == "nt-no-noncoding" ]]; then
	# 	# plant mt
	# 	${_polap_cmd} readassemble \
	# 		"${option_data_type}" \
	# 		--no-noncoding \
	# 		--downsample "${_brg_coverage}" \
	# 		-l "${resolved_fastq}" \
	# 		-o "${_brg_target}"
	# else
	# 	# plant mt
	# 	# use hifi100k.sh and ont100k.sh
	# 	# assemble ptDNA first
	# 	${_polap_cmd} readassemble \
	# 		"${option_data_type}" \
	# 		--no-noncoding \
	# 		--readassemble-t 300000 \
	# 		--readassemble-n 100 \
	# 		-w "${_brg_omega}" \
	# 		--downsample "${_brg_coverage}" \
	# 		-l "${resolved_fastq}" \
	# 		-o "${_brg_target}"
	# fi

	conda deactivate

	# Summarize results after job (with previously defined summary function)
	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"
	_polap_lib_timing-get_system_info >>"${_timing_txt}"

	# Save some results
	rsync -azuq --max-size=5M \
		"${_brg_target}/" "${_brg_rundir}/"

	# if [[ -s "${long_sra}".fastq ]]; then
	# 	_log_echo1 rm -f "${long_sra}".fastq
	# 	rm -f "${long_sra}".fastq
	# fi
	if [[ "${_brg_cleanup}" == "on" ]]; then
		if [[ -d "${_brg_target}" ]]; then
			_log_echo1 rm -rf "${_brg_target}"
			rm -rf "${_brg_target}"
		fi
	fi
}

# 2025-10-20
run-polap-miniassemble_genus_species() {
	# if [[ -n "${_brg_outdir}" ]]; then
	# 	_log_echo0 "brg_outdir 1: $_brg_outdir"
	# fi
	# local _brg_outdir="${1:-$_brg_outdir}"
	# local _brg_sindex="${2:-$_brg_sindex}"
	# if [[ -n "${_brg_outdir}" ]]; then
	# 	_log_echo0 "brg_outdir 2: $_brg_outdir"
	# fi
	#
	# local _brg_outdir
	# local _brg_sindex

	local bolap_cmd="${FUNCNAME%%_*}"

	help_message=$(
		cat <<EOF
Name:
  bolap - $bolap_cmd

Synopsis:
  bolap $bolap_cmd

Description:
  bolap

Examples:
  Execute polap miniassemble:
    bolap miniassemble -s Vigna_radiata

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	# defaults
	local _brg_type="miniasm"
	local _brg_omega="1500"
	local _brg_downsample="no-downsample"
	local _brg_redo="off"
	local _brg_cleanup="on"

	parse_commandline() {
		# set -- "${_brg_unknown_opts[@]}"
		set -- "${_brg_args[@]}"
		# _log_echo0 "parse: $@"

		# source "${_POLAPLIB_DIR}/polap-cmd-version.sh" # '.' means 'source'
		while test $# -gt 0; do
			_key="$1"
			case "$_key" in
			-s)
				test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
				_brg_outdir="${2%/}"
				_brg_outdir_list+=("$_brg_outdir") # append to array
				shift
				;;
			-i)
				test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
				_brg_sindex="$2"
				shift
				;;
			--type)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_type="$2"
					shift || true
				fi
				;;
			-w)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_omega="$2"
					shift || true
				fi
				;;
			--redo)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_redo="on"
					shift || true
				fi
				;;
			--no-cleanup)
				# _log_echo0 "--no-cleanup"
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_cleanup="off"
					shift || true
				fi
				;;
			-d)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_downsample="downsample"
					shift || true
				fi
				;;
			--dry-run)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_dry="on"
					shift || true
				fi
				;;
			-*)
				_log_echo0 "[INFO] no such options: $1"
				;;
			*) ;;
			esac
			shift || true
		done
	}

	# declare -n ref="help_message"
	# if [[ "${_brg_help}" == "on" ]]; then
	# 	local manfile=$(_bolap_lib_man-convert_help_message "$ref" "${bolap_cmd}")
	# 	man "$manfile"
	# 	rm -f "$manfile"
	# 	return
	# fi

	_polap_lib_help-maybe-show "$bolap_cmd" help_message || return 0

	source "${_POLAPLIB_DIR}/polap-variables-data.sh"
	local long_sra="${_long["$_brg_target"]}"

	parse_commandline
	_log_echo1 "target: $_brg_target"

	# local mtn="${_mtn["$_brg_target"]}"
	# local ptn="${_ptn["$_brg_target"]}"
	local mtn=3
	local ptn=2

	mkdir -p "${_brg_outdir_i}"

	local platform="${_platform["$_brg_target"]}"

	local _brg_coverage="50m"

	local option_data_type="--nano-raw"

	_log_echo1 "Asemble organelle genomes using miniasm as a reference generator with ${option_data_type}"

	# Always redo
	if [[ "${_brg_cleanup}" == "on" ]]; then
		rm -rf "${_brg_rundir}"
	fi
	mkdir -p "${_brg_rundir}"

	# Always redo
	if [[ "${_brg_cleanup}" == "on" ]]; then
		rm -rf "${_brg_target}"
	fi
	mkdir -p "${_brg_target}"

	# Start memory logger
	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"

	_polap_lib_conda-ensure_conda_env polap || exit 1

	# check out input data
	# ln -sf "${_brg_tmpdir}/l.fq" "${long_sra}.fastq"
	# _brg_input_data="${long_sra}.fastq"
	if [[ ! -s "${long_sra}".fastq ]]; then
		# data-long_genus_species "${_brg_outdir}"
		_log_echo0 "No input data: ${long_sra}.fastq"
		_log_echo0 "run data long -s ${_brg_outdir}"
		_log_echo0 "data downsample long -s ${_brg_outdir}"
		return 0
	fi

	local resolved_fastq="${long_sra}.fastq"

	# elif [[ "${_brg_type}" == "miniasm" ]]; then
	_log_echo2 "plant mt $_brg_verbose_str"

	_log_echo1 ${_polap_cmd} miniassemble \
		"${option_data_type}" \
		-l "${resolved_fastq}" \
		--readassemble-mtn "${mtn}" \
		${_brg_verbose_str} \
		-o "${_brg_target}"

	${_polap_cmd} miniassemble \
		"${option_data_type}" \
		-l "${resolved_fastq}" \
		--readassemble-mtn "${mtn}" \
		${_brg_verbose_str} \
		-o "${_brg_target}"

	conda deactivate

	# Summarize results after job (with previously defined summary function)
	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"
	_polap_lib_timing-get_system_info >>"${_timing_txt}"

	# Save some results
	rsync -azuq --max-size=5M \
		"${_brg_target}/" "${_brg_rundir}/"

	if [[ "${_brg_cleanup}" == "on" ]]; then
		if [[ -d "${_brg_target}" ]]; then
			_log_echo1 rm -rf "${_brg_target}"
			rm -rf "${_brg_target}"
		fi
	fi
}

# 2025-10-30
run-polap-assemble_genus_species() {

	local bolap_cmd="${FUNCNAME%%_*}"

	help_message=$(
		cat <<EOF
Name:
  bolap - $bolap_cmd

Synopsis:
  bolap $bolap_cmd

Description:
  bolap

Examples:
  Execute polap assemble:
    bolap run-polap-assemble -s Vigna_radiata

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	# defaults
	local _brg_type="miniasm"
	local _brg_omega="1500"
	local _brg_downsample="no-downsample"
	local _brg_redo="off"
	local _brg_cleanup="on"

	parse_commandline() {
		# set -- "${_brg_unknown_opts[@]}"
		set -- "${_brg_args[@]}"
		# _log_echo0 "parse: $@"

		# source "${_POLAPLIB_DIR}/polap-cmd-version.sh" # '.' means 'source'
		while test $# -gt 0; do
			_key="$1"
			case "$_key" in
			-s)
				test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
				_brg_outdir="${2%/}"
				_brg_outdir_list+=("$_brg_outdir") # append to array
				shift
				;;
			-i)
				test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
				_brg_sindex="$2"
				shift
				;;
			--type)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_type="$2"
					shift || true
				fi
				;;
			-w)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_omega="$2"
					shift || true
				fi
				;;
			--redo)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_redo="on"
					shift || true
				fi
				;;
			--no-cleanup)
				# _log_echo0 "--no-cleanup"
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_cleanup="off"
					shift || true
				fi
				;;
			-d)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_downsample="downsample"
					shift || true
				fi
				;;
			--dry-run)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_dry="on"
					shift || true
				fi
				;;
			esac
			shift || true
		done
	}

	# declare -n ref="help_message"
	# if [[ "${_brg_help}" == "on" ]]; then
	# 	local manfile=$(_bolap_lib_man-convert_help_message "$ref" "${bolap_cmd}")
	# 	man "$manfile"
	# 	rm -f "$manfile"
	# 	return
	# fi

	_polap_lib_help-maybe-show "$bolap_cmd" help_message || return 0

	source "${_POLAPLIB_DIR}/polap-variables-data.sh"
	local long_sra="${_long["$_brg_target"]}"
	local short_sra="${_short["$_brg_target"]}"

	local option_data_type="--nano-raw"

	parse_commandline
	_log_echo1 "target: $_brg_target"

	mkdir -p "${_brg_outdir_i}"

	local platform="${_platform["$_brg_target"]}"

	_log_echo1 "Finalize organelle genomes using miniasm, oatk-pathfinder, fmlrc2"

	# Always redo
	if [[ "${_brg_cleanup}" == "on" ]]; then
		rm -rf "${_brg_rundir}"
	fi
	mkdir -p "${_brg_rundir}"

	# Always redo
	if [[ "${_brg_cleanup}" == "on" ]]; then
		rm -rf "${_brg_target}"
	fi
	mkdir -p "${_brg_target}"

	# Start memory logger
	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"

	# check out input data
	# ln -sf "${_brg_tmpdir}/l.fq" "${long_sra}.fastq"
	# _brg_input_data="${long_sra}.fastq"
	if [[ ! -s "${long_sra}".fastq ]]; then
		# data-long_genus_species "${_brg_outdir}"
		_log_echo0 "No input data: ${long_sra}.fastq"
		_log_echo0 "run data long -s ${_brg_outdir}"
		_log_echo0 "data downsample long -s ${_brg_outdir}"
		return 0
	fi

	${_polap_cmd} assemble \
		"${option_data_type}" \
		-l "${long_sra}.fastq" \
		-a "${short_sra}_1.fastq" \
		-b "${short_sra}_2.fastq" \
		${_brg_verbose_str} \
		-o "${_brg_target}"

	# Summarize results after job (with previously defined summary function)
	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"
	_polap_lib_timing-get_system_info >>"${_timing_txt}"

	# Save some results
	rsync -azuq --max-size=5M \
		"${_brg_target}/" "${_brg_rundir}/"

	if [[ "${_brg_cleanup}" == "on" ]]; then
		if [[ -d "${_brg_target}" ]]; then
			_log_echo1 rm -rf "${_brg_target}"
			rm -rf "${_brg_target}"
		fi
	fi
}

run-polap-coverage_genus_species() {

	local bolap_cmd="${FUNCNAME%%_*}"

	help_message=$(
		cat <<EOF
Name:
  bolap - $bolap_cmd

Synopsis:
  bolap $bolap_cmd

Description:
  bolap

Examples:
  Execute polap assemble:
    bolap run-polap-assemble -s Vigna_radiata

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	# defaults
	local _brg_type="miniasm"
	local _brg_omega="1500"
	local _brg_downsample="no-downsample"
	local _brg_redo="off"
	local _brg_cleanup="on"

	parse_commandline() {
		# set -- "${_brg_unknown_opts[@]}"
		set -- "${_brg_args[@]}"
		# _log_echo0 "parse: $@"

		# source "${_POLAPLIB_DIR}/polap-cmd-version.sh" # '.' means 'source'
		while test $# -gt 0; do
			_key="$1"
			case "$_key" in
			-s)
				test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
				_brg_outdir="${2%/}"
				_brg_outdir_list+=("$_brg_outdir") # append to array
				shift
				;;
			-i)
				test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
				_brg_sindex="$2"
				shift
				;;
			--type)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_type="$2"
					shift || true
				fi
				;;
			-w)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_omega="$2"
					shift || true
				fi
				;;
			--redo)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_redo="on"
					shift || true
				fi
				;;
			--no-cleanup)
				# _log_echo0 "--no-cleanup"
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_cleanup="off"
					shift || true
				fi
				;;
			-d)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_downsample="downsample"
					shift || true
				fi
				;;
			--dry-run)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_dry="on"
					shift || true
				fi
				;;
			esac
			shift || true
		done
	}

	_polap_lib_help-maybe-show "$bolap_cmd" help_message || return 0

	source "${_POLAPLIB_DIR}/polap-variables-data.sh"
	local long_sra="${_long["$_brg_target"]}"
	local short_sra="${_short["$_brg_target"]}"

	local option_data_type="--nano-raw"

	parse_commandline
	_log_echo1 "target: $_brg_target"

	mkdir -p "${_brg_outdir_i}"

	local platform="${_platform["$_brg_target"]}"

	_log_echo1 "Finalize organelle genomes using miniasm, oatk-pathfinder, fmlrc2"

	# Always redo
	if [[ "${_brg_cleanup}" == "on" ]]; then
		rm -rf "${_brg_rundir}"
	fi
	mkdir -p "${_brg_rundir}"

	# Always redo
	if [[ "${_brg_cleanup}" == "on" ]]; then
		rm -rf "${_brg_target}"
	fi
	mkdir -p "${_brg_target}"

	# Start memory logger
	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"

	# check out input data
	# ln -sf "${_brg_tmpdir}/l.fq" "${long_sra}.fastq"
	# _brg_input_data="${long_sra}.fastq"
	if [[ ! -s "${long_sra}".fastq ]]; then
		# data-long_genus_species "${_brg_outdir}"
		_log_echo0 "No input data: ${long_sra}.fastq"
		_log_echo0 "run data long -s ${_brg_outdir}"
		_log_echo0 "data downsample long -s ${_brg_outdir}"
		return 0
	fi

	if [[ -s "${_brg_outdir_i}/polap-assemble/mt.1.fasta" ]]; then
		${_polap_cmd} coverage \
			-l "${long_sra}.fastq" \
			--infile "${_brg_outdir_i}/polap-assemble/mt.1.fasta" \
			-o "${_brg_target}"
	fi

	if [[ -s "${_brg_outdir_i}/polap-assemble/pt.1.fasta" ]]; then
		${_polap_cmd} coverage \
			--plastid \
			-l "${long_sra}.fastq" \
			--infile "${_brg_outdir_i}/polap-assemble/pt.1.fasta" \
			-o "${_brg_target}"
	fi

	# Summarize results after job (with previously defined summary function)
	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"
	_polap_lib_timing-get_system_info >>"${_timing_txt}"

	# Save some results
	rsync -azuq --max-size=5M \
		"${_brg_target}/" "${_brg_rundir}/"

	if [[ "${_brg_cleanup}" == "on" ]]; then
		if [[ -d "${_brg_target}" ]]; then
			_log_echo1 rm -rf "${_brg_target}"
			rm -rf "${_brg_target}"
		fi
	fi
}

run-polap-extract-using-oatk_genus_species() {

	local bolap_cmd="${FUNCNAME%%_*}"

	help_message=$(
		cat <<EOF
Name:
  bolap - $bolap_cmd

Synopsis:
  bolap $bolap_cmd

Description:
  bolap

Examples:
  Execute what?:
    bolap syncassemble -s Vigna_radiata

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	# defaults
	local _brg_type="mtpt"
	local _brg_omega="1500"
	local _brg_downsample="no-downsample"
	local _brg_redo="off"
	local _brg_cleanup="on"

	parse_commandline() {
		set -- "${_brg_args[@]}"

		# source "${_POLAPLIB_DIR}/polap-cmd-version.sh" # '.' means 'source'
		while test $# -gt 0; do
			_key="$1"
			case "$_key" in
			-s)
				test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
				_brg_outdir="${2%/}"
				_brg_outdir_list+=("$_brg_outdir") # append to array
				shift
				;;
			-i)
				test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
				_brg_sindex="$2"
				shift
				;;
			--type)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_type="$2"
					shift || true
				fi
				;;
			-w)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_omega="$2"
					shift || true
				fi
				;;
			--redo)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_redo="on"
					shift || true
				fi
				;;
			--no-cleanup)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_cleanup="off"
					shift || true
				fi
				;;
			-d)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_downsample="downsample"
					shift || true
				fi
				;;
			--dry-run)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_dry="on"
					shift || true
				fi
				;;
			-*)
				_log_echo0 "[INFO] no such options: $1"
				;;
			*)
				break
				;;
			esac
			shift || true
		done
	}

	# --help
	_polap_lib_help-maybe-show "$bolap_cmd" help_message || return 0

	parse_commandline

	# _brg_outdir: -s
	# _brg_sindex: -i
	source "${_POLAPLIB_DIR}/polap-variables-data.sh"
	local long_sra="${_long["$_brg_target"]}"

	# main
	#
	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"

	# redefine
	local _brg_miniasm="${_brg_outdir_i}/polap-miniassemble"
	local _brg_rundir="${_brg_outdir_i}/${_brg_title}"

	mkdir -p "${_brg_outdir_i}"

	local platform="${_platform["$_brg_target"]}"

	mkdir -p "${_brg_target}"

	# copy mt.1.gfa
	if [[ ! -s "${_brg_miniasm}/mt.1.gfa" ]]; then
		_log_echo0 "[ERROR] no MT assembly to extract: mt.1.gfa"
		return 0
	fi

	if [[ ! -s "${_brg_miniasm}/pt.1.gfa" ]]; then
		_log_echo0 "[ERROR] no PT assembly to extract: pt.1.gfa"
		return 0
	fi

	${_polap_cmd} extract \
		--infile1 "${_brg_miniasm}/pt.1.gfa" \
		--infile2 "${_brg_miniasm}/mt.1.gfa" \
		-o "${_brg_target}" -v

	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"

	# Save some results
	rsync -azuq --max-size=5M \
		"${_brg_target}/" "${_brg_rundir}/"

	if [[ "${_brg_cleanup}" == "on" ]]; then
		if [[ -d "${_brg_target}" ]]; then
			_log_echo1 rm -rf "${_brg_target}"
			rm -rf "${_brg_target}"
		fi
	fi
}

# 2025-10-29
run-polap-polish-longshort_genus_species() {

	local bolap_cmd="${FUNCNAME%%_*}"

	help_message=$(
		cat <<EOF
Name:
  bolap - $bolap_cmd

Synopsis:
  bolap $bolap_cmd

Description:
  bolap

Examples:
  Execute what?:
    bolap syncassemble -s Vigna_radiata

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	# defaults
	local _brg_type="mtpt"
	local _brg_omega="1500"
	local _brg_downsample="no-downsample"
	local _brg_redo="off"
	local _brg_cleanup="on"

	parse_commandline() {
		set -- "${_brg_args[@]}"

		# source "${_POLAPLIB_DIR}/polap-cmd-version.sh" # '.' means 'source'
		while test $# -gt 0; do
			_key="$1"
			case "$_key" in
			-s)
				test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
				_brg_outdir="${2%/}"
				_brg_outdir_list+=("$_brg_outdir") # append to array
				shift
				;;
			-i)
				test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
				_brg_sindex="$2"
				shift
				;;
			--type)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_type="$2"
					shift || true
				fi
				;;
			-w)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_omega="$2"
					shift || true
				fi
				;;
			--redo)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_redo="on"
					shift || true
				fi
				;;
			--no-cleanup)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_cleanup="off"
					shift || true
				fi
				;;
			-d)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_downsample="downsample"
					shift || true
				fi
				;;
			--dry-run)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_dry="on"
					shift || true
				fi
				;;
			-*)
				_log_echo0 "[INFO] no such options: $1"
				;;
			*)
				break
				;;
			esac
			shift || true
		done
	}

	# --help
	_polap_lib_help-maybe-show "$bolap_cmd" help_message || return 0

	parse_commandline

	# _brg_outdir: -s
	# _brg_sindex: -i
	source "${_POLAPLIB_DIR}/polap-variables-data.sh"
	local long_sra="${_long["$_brg_target"]}"
	local short_sra="${_short["$_brg_target"]}"

	# main
	#
	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"

	# redefine
	local _brg_extract="${_brg_outdir_i}/polap-assemble/extract"
	local _brg_rundir="${_brg_outdir_i}/${_brg_title}"

	mkdir -p "${_brg_outdir_i}"

	local platform="${_platform["$_brg_target"]}"

	rm -rf "${_brg_target}"
	mkdir -p "${_brg_target}"

	# copy mt.1.gfa

	if [[ ! -s "${_brg_extract}/oatk.pltd.ctg.fasta" ]]; then
		_log_echo0 "[ERROR] no draft plastid sequence to polish"
		return 0
	fi

	if [[ ! -s "${_brg_extract}/oatk.mito.ctg.fasta" ]]; then
		_log_echo0 "[ERROR] no draft mitochondrial sequence to polish"
		return 0
	fi

	if [[ ! -s "${short_sra}_1.fastq" ]]; then
		_log_echo0 "[ERROR] no short-read1 data to use for polishing"
		return 0
	fi

	if [[ ! -s "${short_sra}_2.fastq" ]]; then
		_log_echo0 "[ERROR] no short-read2 data to use for polishing"
		return 0
	fi

	if [[ ! -s "${long_sra}.fastq" ]]; then
		_log_echo0 "[ERROR] no long-read data to use for polishing"
		return 0
	fi

	# _log_echo0 "START: ptDNA polish"
	# ${_polap_cmd} polish-longshort \
	# 	--plastid \
	# 	-l "${long_sra}.fastq" \
	# 	-a "${short_sra}_1.fastq" \
	# 	-b "${short_sra}_2.fastq" \
	# 	--infile1 "${_brg_extract}/oatk.mito.ctg.fasta" \
	# 	--infile2 "${_brg_extract}/oatk.pltd.ctg.fasta" \
	# 	--outfile "${_brg_target}/pt.1.fasta" \
	# 	-v \
	# 	-o "${_brg_target}"
	# _log_echo0 "END: ptDNA polish"

	${_polap_cmd} polish-longshort \
		-l "${long_sra}.fastq" \
		-a "${short_sra}_1.fastq" \
		-b "${short_sra}_2.fastq" \
		--infile1 "${_brg_extract}/oatk.mito.ctg.fasta" \
		--infile2 "${_brg_extract}/oatk.pltd.ctg.fasta" \
		--outfile "${_brg_target}/mt.1.fasta" \
		-v \
		-o "${_brg_target}"

	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"

	# Save some results
	rsync -azuq --max-size=5M \
		"${_brg_target}/" "${_brg_rundir}/"

	_log_echo0 "mtDNA polished: ${_brg_rundir}/mt.1.fasta"
	_log_echo0 "ptDNA polished: ${_brg_rundir}/pt.1.fasta"

	if [[ "${_brg_cleanup}" == "on" ]]; then
		if [[ -d "${_brg_target}" ]]; then
			_log_echo1 rm -rf "${_brg_target}"
			rm -rf "${_brg_target}"
		fi
	fi
}

run-polap-polish-long_genus_species() {
	local _brg_outdir="${1:-$_brg_outdir}"
	local _brg_sindex="${2:-$_brg_sindex}"

	local bolap_cmd="${FUNCNAME%%_*}"

	help_message=$(
		cat <<EOF
Name:
  bolap - $bolap_cmd

Synopsis:
  bolap $bolap_cmd

Description:
  bolap

Examples:
  Execute what?:
    bolap syncassemble -s Vigna_radiata

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	# defaults
	local _brg_type="mtpt"
	local _brg_omega="1500"
	local _brg_downsample="no-downsample"
	local _brg_redo="off"
	local _brg_cleanup="on"

	parse_commandline() {
		set -- "${_brg_unknown_opts[@]}"

		# source "${_POLAPLIB_DIR}/polap-cmd-version.sh" # '.' means 'source'
		while test $# -gt 0; do
			_key="$1"
			case "$_key" in
			--type)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_type="$2"
					shift || true
				fi
				;;
			-w)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_omega="$2"
					shift || true
				fi
				;;
			--redo)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_redo="on"
					shift || true
				fi
				;;
			--no-cleanup)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_cleanup="off"
					shift || true
				fi
				;;
			-d)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_downsample="downsample"
					shift || true
				fi
				;;
			--dry-run)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_dry="on"
					shift || true
				fi
				;;
			-*)
				_log_echo0 "[INFO] no such options: $1"
				;;
			*)
				break
				;;
			esac
			shift || true
		done
	}

	declare -n ref="help_message"
	if [[ "${_brg_help}" == "on" ]]; then
		local manfile=$(_bolap_lib_man-convert_help_message "$ref" "${bolap_cmd}")
		man "$manfile"
		rm -f "$manfile"
		return
	fi

	source "${_POLAPLIB_DIR}/polap-variables-data.sh"
	local long_sra="${_long["$_brg_target"]}"

	parse_commandline

	# main
	#
	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"

	# redefine
	local _brg_miniasm="${_brg_outdir_i}/polap-readassemble"
	local _brg_rundir="${_brg_outdir_i}/${_brg_title}"

	mkdir -p "${_brg_outdir_i}"

	local platform="${_platform["$_brg_target"]}"

	mkdir -p "${_brg_target}"

	# copy mt.1.gfa
	cp -p "${_brg_miniasm}/mt.1.gfa" \
		"${_brg_target}"

	# Note: we need to edit mt.1.txt manusally.
	if [[ -s "${_brg_miniasm}/mt.1.txt" ]]; then
		cp -p "${_brg_miniasm}/mt.1.txt" \
			"${_brg_target}"
	fi

	cp -p "${_brg_miniasm}/pt.1.gfa" \
		"${_brg_target}/pt.1.gfa"

	if [[ -s "${_brg_miniasm}/mt.1.txt" ]]; then
		${_polap_cmd} convert gfa2fasta \
			"${_brg_target}/mt.1.gfa" \
			"${_brg_target}/mt.1.fa" \
			--ids "${_brg_target}/mt.1.txt"
	else
		${_polap_cmd} convert gfa2fasta \
			"${_brg_target}/mt.1.gfa" \
			"${_brg_target}/mt.1.fa"
	fi

	${_polap_cmd} convert gfa2fasta \
		"${_brg_target}/pt.1.gfa" \
		"${_brg_target}/pt.1.fa"

	${_polap_cmd} polish2 \
		-l "${long_sra}.fastq" \
		-p "${_brg_target}/pt.1.fa" \
		-f "${_brg_target}/pt.1.fasta"

	${_polap_cmd} polish2 \
		-l "${long_sra}.fastq" \
		-p "${_brg_target}/mt.1.fa" \
		-f "${_brg_target}/mt.1.fasta"

	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"

	# ${_polap_cmd} isomer \
	# 	-l "${long_sra}.fastq" \
	# 	-o "${_brg_target}" \
	# 	--pt-ref "${_brg_target}/pt.1.fasta" \
	# 	--mt-ref "${_brg_target}/mt.1.fasta"

	# Save some results
	rsync -azuq --max-size=5M \
		"${_brg_target}/" "${_brg_rundir}/"

	if [[ "${_brg_cleanup}" == "on" ]]; then
		if [[ -d "${_brg_target}" ]]; then
			_log_echo1 rm -rf "${_brg_target}"
			rm -rf "${_brg_target}"
		fi
	fi
}

run-polap-mtpt_genus_species() {
	local _brg_outdir="${1:-$_brg_outdir}"
	local _brg_sindex="${2:-$_brg_sindex}"

	# Debug trace with timestamps, file:line, and function
	if [[ "${POLAP_DEBUG:-0}" -eq 1 ]]; then
		export PS4='+$(date "+%F %T") ${BASH_SOURCE##*/}:${LINENO}:${FUNCNAME[0]}: '
		set -o errtrace  # inherit ERR in functions
		set -o functrace # inherit DEBUG/RETURN traps
		set -x           # trace execution
		trap 'st=$?; echo "ERR ${BASH_SOURCE[0]}:${LINENO}: ${BASH_COMMAND}" >&2; exit $st' ERR
	fi

	local bolap_cmd="${FUNCNAME%%_*}"

	help_message=$(
		cat <<EOF
Name:
  bolap - $bolap_cmd

Synopsis:
  bolap $bolap_cmd

Description:
  bolap

Examples:
  Execute what?:
    bolap syncassemble -s Vigna_radiata

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	# defaults
	local _brg_type="mtpt"
	local _brg_omega="1500"
	local _brg_downsample="no-downsample"
	local _brg_redo="off"
	local _brg_cleanup="on"

	parse_commandline() {
		set -- "${_brg_unknown_opts[@]}"

		# source "${_POLAPLIB_DIR}/polap-cmd-version.sh" # '.' means 'source'
		while test $# -gt 0; do
			_key="$1"
			case "$_key" in
			--type)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_type="$2"
					shift || true
				fi
				;;
			-w)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_omega="$2"
					shift || true
				fi
				;;
			--redo)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_redo="on"
					shift || true
				fi
				;;
			--no-cleanup)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_cleanup="off"
					shift || true
				fi
				;;
			-d)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_downsample="downsample"
					shift || true
				fi
				;;
			--dry-run)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_dry="on"
					shift || true
				fi
				;;
			-*)
				_log_echo0 "[INFO] no such options: $1"
				;;
			*)
				break
				;;
			esac
			shift || true
		done
	}

	declare -n ref="help_message"
	if [[ "${_brg_help:-off}" == "on" ]]; then
		local manfile=$(_bolap_lib_man-convert_help_message "$ref" "${bolap_cmd}")
		man "$manfile"
		rm -f "$manfile"
		return
	fi

	source "${_POLAPLIB_DIR}/polap-variables-data.sh"
	local long_sra="${_long["$_brg_target"]}"

	parse_commandline

	# main

	# redefine
	local _brg_miniasm="${_brg_outdir_i}/polap-readassemble"
	local _brg_polap_assemble="${_brg_outdir_i}/polap-assemble"
	local _brg_rundir="${_brg_outdir_i}/${_brg_title}"

	mkdir -p "${_brg_outdir_i}"

	local platform="${_platform["$_brg_target"]}"

	mkdir -p "${_brg_target}"

	# copy mt.1.gfa
	cp -p "${_brg_polap_assemble}/mt.1.fasta" \
		"${_brg_target}"
	cp -p "${_brg_polap_assemble}/pt.1.fasta" \
		"${_brg_target}"

	${_polap_cmd} mtpt \
		-o "${_brg_target}" \
		-l "${long_sra}.fastq" \
		--pt-ref "${_brg_target}/pt.1.fasta" \
		--mt-ref "${_brg_target}/mt.1.fasta"

	# Save some results
	rsync -azuq --max-size=5M \
		"${_brg_target}/" "${_brg_rundir}/"

	if [[ "${_brg_cleanup}" == "on" ]]; then
		if [[ -d "${_brg_target}" ]]; then
			_log_echo1 rm -rf "${_brg_target}"
			rm -rf "${_brg_target}"
		fi
	fi
}

run-polap-readassemble-animal-mt_genus_species() {
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"

	run-polap-readassemble_genus_species "${_brg_outdir}" "${_brg_sindex}" \
		-t "animal"
}

run-polap-readassemble-nt_genus_species() {
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"

	run-polap-readassemble_genus_species "${_brg_outdir}" "${_brg_sindex}" \
		-t "nt"
}

run-polap-readassemble-mt_genus_species() {
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"

	run-polap-readassemble_genus_species "${_brg_outdir}" "${_brg_sindex}" \
		-t "mt"
}

run-polap-readassemble-miniasm_genus_species() {
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"

	run-polap-readassemble_genus_species "${_brg_outdir}" "${_brg_sindex}" \
		-t "miniasm"
}

run-polap-readassemble-pt_genus_species() {
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"

	run-polap-readassemble_genus_species "${_brg_outdir}" "${_brg_sindex}" \
		-t "pt" -w 1500
}

run-polap-syncasm_genus_species() {
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"
	source "${_POLAPLIB_DIR}/polap-variables-data.sh"

	mkdir -p "${_brg_rundir}"

	local platform="${_platform["$_brg_target"]}"

	local option_data_type="--nano-raw"
	if [[ "${platform}" == "PACBIO_SMRT" ]]; then
		option_data_type="--pacbio-hifi"
	elif [[ "${platform}" == "PACBIO_CLR" ]]; then
		option_data_type="--pacbio-raw"
	elif [[ "${platform}" == "ONT" ]]; then
		local option_data_type="--nano-raw"
	fi

	ln -sf "${_brg_tmpdir}/l.fq" "${long_sra}.fastq"
	local _brg_input_data="${long_sra}.fastq"
	local resolved_fastq="${long_sra}.fastq"

	_log_echo "Do: assemble mtDNA using syncasm as a seed generator with ${option_data_type}"

	# Start memory logger
	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"

	_polap_lib_conda-ensure_conda_env polap || exit 1

	echo ${_polap2_cmd} syncasm \
		"${option_data_type}" \
		-l "${resolved_fastq}" \
		-o "${_brg_target}"

	conda deactivate

	# Summarize results after job (with previously defined summary function)
	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"
	_polap_lib_timing-get_system_info >>"${_timing_txt}"

	# Save some results
	rsync -azuq --max-size=5M "${_brg_target}/" "${_brg_rundir}/"
}

run-polap-syncassemble_genus_species() {
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"
	source "${_POLAPLIB_DIR}/polap-variables-data.sh"

	mkdir -p "${_brg_rundir}"

	local platform="${_platform["$_brg_target"]}"

	local option_data_type="--nano-raw"
	if [[ "${platform}" == "PACBIO_SMRT" ]]; then
		option_data_type="--pacbio-hifi"
	elif [[ "${platform}" == "PACBIO_CLR" ]]; then
		option_data_type="--pacbio-raw"
	elif [[ "${platform}" == "ONT" ]]; then
		local option_data_type="--nano-raw"
	fi

	ln -sf "${_brg_tmpdir}/l.fq" "${long_sra}.fastq"
	local _brg_input_data="${long_sra}.fastq"
	local resolved_fastq="${long_sra}.fastq"

	_log_echo "Do: assemble pt/mtDNA using syncfilter reads with ${option_data_type}"

	# Start memory logger
	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"

	_polap_lib_conda-ensure_conda_env polap || exit 1

	${_polap_cmd} syncassemble \
		"${option_data_type}" \
		-l "${resolved_fastq}" \
		-o "${_brg_target}"

	conda deactivate

	# Summarize results after job (with previously defined summary function)
	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"
	_polap_lib_timing-get_system_info >>"${_timing_txt}"

	# Save some results
	rsync -azuq --max-size=5M "${_brg_target}/" "${_brg_rundir}/"
}

run-polap-oatk_genus_species() {
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"
	source "${_POLAPLIB_DIR}/polap-variables-data.sh"

	mkdir -p "${_brg_rundir}"
	# rm -rf "${_brg_target}/oatk"
	# rm -rf "${_brg_target}/annotate-read-oatk"

	local platform="${_platform["$_brg_target"]}"

	local option_data_type="--nano-raw"
	if [[ "${platform}" == "PACBIO_SMRT" ]]; then
		option_data_type="--pacbio-hifi"
	elif [[ "${platform}" == "PACBIO_CLR" ]]; then
		option_data_type="--pacbio-raw"
	elif [[ "${platform}" == "ONT" ]]; then
		local option_data_type="--nano-raw"
	fi

	ln -sf "${_brg_tmpdir}/l.fq" "${long_sra}.fastq"
	local _brg_input_data="${long_sra}.fastq"

	_log_echo "Do: assemble pt/mtDNA using oatk's components with ${option_data_type}"

	# Start memory logger
	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"

	_polap_lib_conda-ensure_conda_env polap || exit 1

	${_polap_cmd} oatk \
		"${option_data_type}" \
		-l "${_brg_input_data}" \
		-o "${_brg_target}"

	conda deactivate

	# Summarize results after job (with previously defined summary function)
	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"
	_polap_lib_timing-get_system_info >>"${_timing_txt}"

	# Save some results
	rsync -azuq --max-size=5M "${_brg_target}/" "${_brg_rundir}/"
}

get-conf-long-sra-id_genus_species() {
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"

	local _brg_inum=0
	local _brg_adir _brg_title _brg_target _brg_rundir _brg_outdir_i
	local _timing_txt _stdout_txt _memlog_file _summary_file

	brg_common_setup \
		_brg_outdir _brg_sindex _brg_adir _brg_title \
		_brg_target _brg_rundir _brg_outdir_i \
		_timing_txt _stdout_txt _memlog_file _summary_file

	if [[ -v _long["$_brg_target"] ]]; then
		local long_sra="${_long["$_brg_target"]}"
		echo "${long_sra}"
	fi
}

# input: _brg_sindex
# output: _brg_coverage
#
lib-polap-readassemble-brg-coverage() {
	if [[ "${_brg_sindex}" == "0" ]]; then
		_brg_coverage="10g"
	elif [[ "${_brg_sindex}" == "1" ]]; then
		_brg_coverage="1g"
	elif [[ "${_brg_sindex}" == "2" ]]; then
		_brg_coverage="200m"
	else
		_brg_coverage="50m"
	fi
}

lib-polap-readassemble-pt-brg-coverage() {
	if [[ "${_brg_sindex}" == "0" ]]; then
		_brg_coverage="3g"
	elif [[ "${_brg_sindex}" == "1" ]]; then
		_brg_coverage="1g"
	elif [[ "${_brg_sindex}" == "2" ]]; then
		_brg_coverage="200m"
	else
		_brg_coverage="50m"
	fi
}

lib-polap-readassemble-hifi-brg-coverage() {
	if [[ "${_brg_sindex}" == "0" ]]; then
		_brg_coverage="10"
	elif [[ "${_brg_sindex}" == "1" ]]; then
		_brg_coverage="2"
	elif [[ "${_brg_sindex}" == "2" ]]; then
		_brg_coverage="1"
	else
		_brg_coverage="0.5"
	fi
}

lib-polap-readassemble-hifi-pt-brg-coverage() {
	if [[ "${_brg_sindex}" == "0" ]]; then
		_brg_coverage="3"
	elif [[ "${_brg_sindex}" == "1" ]]; then
		_brg_coverage="1"
	elif [[ "${_brg_sindex}" == "2" ]]; then
		_brg_coverage="0.5"
	else
		_brg_coverage="0.1"
	fi
}

lib-polap-readassemble-clr-brg-coverage() {
	if [[ "${_brg_sindex}" == "0" ]]; then
		_brg_coverage="10"
	elif [[ "${_brg_sindex}" == "1" ]]; then
		_brg_coverage="2"
	elif [[ "${_brg_sindex}" == "2" ]]; then
		_brg_coverage="1"
	else
		_brg_coverage="0.5"
	fi
}

lib-polap-readassemble-clr-pt-brg-coverage() {
	if [[ "${_brg_sindex}" == "0" ]]; then
		_brg_coverage="30"
	elif [[ "${_brg_sindex}" == "1" ]]; then
		_brg_coverage="1"
	elif [[ "${_brg_sindex}" == "2" ]]; then
		_brg_coverage="0.5"
	else
		_brg_coverage="0.1"
	fi
}

run-readassemble_genus_species() {
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"

	_brg_outdir="${_brg_outdir%/}"

	run-readassemble-mt_genus_species "${_brg_outdir}" "${_brg_sindex}" iterate
	run-readassemble-nt_genus_species "${_brg_outdir}" "${_brg_sindex}" iterate
	run-readassemble-pt_genus_species "${_brg_outdir}" "${_brg_sindex}" iterate

	local _brg_rundir="${_brg_outdir}-${_brg_sindex}"
	rm -rf "${_brg_rundir}"
}

# variables
#
# _brg_adir=t5
# _brg_inum=0
# _brg_outdir=Ophrys_lutea
# _brg_outdir_i=Ophrys_lutea/t5/0
# _brg_outdir_t=Ophrys_lutea/t5
# _brg_rundir=Ophrys_lutea-0
# _brg_runtitledir=Ophrys_lutea-0-polap-readassemble-mt
# _brg_sindex=0
# _brg_target=Ophrys_lutea-0
# _brg_title=polap-readassemble-mt
# _brg_titledir=Ophrys_lutea/t5/0/polap-readassemble-mt
# _memlog_file=Ophrys_lutea/t5/0/memlog-polap-readassemble-mt.csv
# _stdout_txt=Ophrys_lutea/t5/0/stdout-polap-readassemble-mt.txt
# _summary_file=Ophrys_lutea/t5/0/summary-polap-readassemble-mt.txt
# _timing_txt=Ophrys_lutea/t5/0/timing-polap-readassemble-mt.txt
# long_sra=ERR5167480
#
# 2025-08-01
run-readassemble-mt_genus_species() {
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"
	local _brg_iterate="${3:-iterate}"

	local _brg_inum=0
	local _brg_adir _brg_title _brg_target _brg_rundir _brg_outdir_i
	local _timing_txt _stdout_txt _memlog_file _summary_file

	brg_common_setup \
		_brg_outdir _brg_sindex _brg_adir _brg_title \
		_brg_target _brg_rundir _brg_outdir_i \
		_timing_txt _stdout_txt _memlog_file _summary_file

	# Extra folders
	local _brg_outdir_t="${_brg_outdir}/${opt_t_arg}"
	local _brg_titledir="${_brg_outdir_i}/${_brg_title}"
	local _brg_runtitledir="${_brg_rundir}-${_brg_title}"

	if [[ -v _long["$_brg_target"] ]]; then
		local long_sra="${_long["$_brg_target"]}"
	else
		echo "Error: ${_brg_target} because it is not in the CSV."
		return
	fi

	local _brg_coverage="50m"
	lib-polap-readassemble-brg-coverage

	if [[ ! -s "${long_sra}".fastq ]]; then
		data-long_genus_species "${_brg_outdir}"
	fi

	data-downsample-long_genus_species "${_brg_outdir}" "${_brg_sindex}" "${_brg_coverage}"
	readassemble-annotate-mt_genus_species "${_brg_outdir}" "${_brg_sindex}"
	if [[ "${_brg_iterate}" == "iterate" ]]; then
		readassemble-ont-mt-iterate_genus_species "${_brg_outdir}" "${_brg_sindex}"
	else
		readassemble-ont-mt_genus_species "${_brg_outdir}" "${_brg_sindex}"
	fi

	# Save some results
	rsync -azuq --max-size=5M \
		"${_brg_rundir}"/annotate-read-mt/ \
		"${_brg_titledir}"/

	if [[ "${_local_host}" != "$(hostname)" ]]; then
		sync_genus_species "${_brg_outdir}" "${_brg_sindex}" --main-push
	fi

	# rm -rf "${_brg_rundir}"
}

run-readassemble-nt_genus_species() {
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"
	local _brg_iterate="${3:-iterate}"

	local _brg_inum=0
	local _brg_adir _brg_title _brg_target _brg_rundir _brg_outdir_i
	local _timing_txt _stdout_txt _memlog_file _summary_file

	brg_common_setup \
		_brg_outdir _brg_sindex _brg_adir _brg_title \
		_brg_target _brg_rundir _brg_outdir_i \
		_timing_txt _stdout_txt _memlog_file _summary_file

	# Extra folders
	local _brg_outdir_t="${_brg_outdir}/${opt_t_arg}"
	local _brg_titledir="${_brg_outdir_i}/${_brg_title}"
	local _brg_runtitledir="${_brg_rundir}-${_brg_title}"

	if [[ -v _long["$_brg_target"] ]]; then
		local long_sra="${_long["$_brg_target"]}"
	else
		echo "Error: ${_brg_target} because it is not in the CSV."
		return
	fi

	local _brg_coverage="50m"
	lib-polap-readassemble-brg-coverage

	if [[ -s "${long_sra}-10x".fq ]]; then
		echo "${long_sra}-10x".fq >"${_brg_outdir}/${_brg_inum}/l.fastq.path.txt"
	else
		if [[ ! -s "${long_sra}".fastq ]]; then
			data-long_genus_species "${_brg_outdir}"
		fi

	fi
	data-downsample-long_genus_species "${_brg_outdir}" "${_brg_sindex}" "${_brg_coverage}"
	readassemble-annotate-nt_genus_species "${_brg_outdir}" "${_brg_sindex}"
	if [[ "${_brg_iterate}" == "iterate" ]]; then
		readassemble-ont-nt-iterate_genus_species "${_brg_outdir}" "${_brg_sindex}"
	else
		readassemble-ont-nt_genus_species "${_brg_outdir}" "${_brg_sindex}"
	fi

	# Save some results
	rsync -azuq --max-size=5M \
		"${_brg_rundir}"/annotate-read-nt/ \
		"${_brg_titledir}"/

	if [[ "${_local_host}" != "$(hostname)" ]]; then
		sync_genus_species "${_brg_outdir}" "${_brg_sindex}" --main-push
	fi

	# rm -rf "${_brg_rundir}"
}

run-readassemble-pt_genus_species() {
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"
	local _brg_iterate="${3:-iterate}"

	local _brg_inum=0
	local _brg_adir _brg_title _brg_target _brg_rundir _brg_outdir_i
	local _timing_txt _stdout_txt _memlog_file _summary_file

	brg_common_setup \
		_brg_outdir _brg_sindex _brg_adir _brg_title \
		_brg_target _brg_rundir _brg_outdir_i \
		_timing_txt _stdout_txt _memlog_file _summary_file

	# Extra folders
	local _brg_outdir_t="${_brg_outdir}/${opt_t_arg}"
	local _brg_titledir="${_brg_outdir_i}/${_brg_title}"
	local _brg_runtitledir="${_brg_rundir}-${_brg_title}"

	if [[ -v _long["$_brg_target"] ]]; then
		local long_sra="${_long["$_brg_target"]}"
	else
		echo "Error: ${_brg_target} because it is not in the CSV."
		return
	fi

	local _brg_coverage="50m"
	lib-polap-readassemble-brg-coverage

	if [[ ! -s "${long_sra}".fastq ]]; then
		data-long_genus_species "${_brg_outdir}"
	fi

	data-downsample-long_genus_species "${_brg_outdir}" "${_brg_sindex}" "${_brg_coverage}"
	readassemble-annotate-pt_genus_species "${_brg_outdir}" "${_brg_sindex}"
	if [[ "${_brg_iterate}" == "iterate" ]]; then
		readassemble-ont-pt-iterate_genus_species "${_brg_outdir}" "${_brg_sindex}"
	else
		readassemble-ont-pt_genus_species "${_brg_outdir}" "${_brg_sindex}"
	fi

	# Save some results
	rsync -azuq --max-size=5M \
		"${_brg_rundir}"/annotate-read-pt/ \
		"${_brg_titledir}"/

	if [[ "${_local_host}" != "$(hostname)" ]]; then
		sync_genus_species "${_brg_outdir}" "${_brg_sindex}" --main-push
	fi

	# rm -rf "${_brg_rundir}"
}

run-readassemble-aninmal-mt_genus_species() {
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"

	local _brg_inum=0
	local _brg_adir _brg_title _brg_target _brg_rundir _brg_outdir_i
	local _timing_txt _stdout_txt _memlog_file _summary_file

	brg_common_setup \
		_brg_outdir _brg_sindex _brg_adir _brg_title \
		_brg_target _brg_rundir _brg_outdir_i \
		_timing_txt _stdout_txt _memlog_file _summary_file

	# Extra folders
	local _brg_outdir_t="${_brg_outdir}/${opt_t_arg}"
	local _brg_titledir="${_brg_outdir_i}/${_brg_title}"
	local _brg_runtitledir="${_brg_rundir}-${_brg_title}"

	if [[ -v _long["$_brg_target"] ]]; then
		local long_sra="${_long["$_brg_target"]}"
	else
		echo "Error: ${_brg_target} because it is not in the CSV."
		return
	fi

	if [[ ! -s "${long_sra}".fastq ]]; then
		data-long_genus_species "${_brg_outdir}"
	fi

	if [[ "${_brg_sindex}" == "0" ]]; then
		_brg_coverage="10g"
	elif [[ "${_brg_sindex}" == "1" ]]; then
		_brg_coverage="1g"
	elif [[ "${_brg_sindex}" == "2" ]]; then
		_brg_coverage="200m"
	fi

	data-downsample-long_genus_species "${_brg_outdir}" "${_brg_sindex}" "${_brg_coverage}"
	readassemble-annotate-mt_genus_species "${_brg_outdir}" "${_brg_sindex}"
	readassemble-ont-animal-mt_genus_species "${_brg_outdir}" "${_brg_sindex}"

	# Save some results
	rsync -azuq --max-size=5M \
		"${_brg_rundir}"/annotate-read-mt/ \
		"${_brg_titledir}"/

	rm -rf "${_brg_rundir}"
}

sync_genus_species() {
	local _brg_outdir="${1:-$_brg_outdir}"
	local _brg_sindex="${2:-$_brg_sindex}"

	local direction="${3:---push}" # Default: --both
	local maxsize="${4:-5M}"       # Default: --both
	local dryrun="${5:-}"          # Optional: --dry-run

	local _brg_inum=0
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

	# local dimred="${_dimred["$_brg_target"]}"
	local platform="${_platform["$_brg_target"]}"
	local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

	local host="${_host["$_brg_target"]}"
	# local remote_host="${host}"
	local remote_host="thorne"
	local folder="${_brg_rundir}"

	# Determine full local and remote paths
	local local_dir="$(cd "$(dirname "$folder")" && pwd)"
	local folder_name="$(basename "$folder")"
	local local_path="${local_dir}/${folder_name}"
	local remote_path="${remote_host}:${local_path}"
	local main_path="${_local_host}:${local_path}"

	# Prepare rsync options
	local rsync_opts=(-az --info=progress2 --no-inc-recursive)
	[[ "$dryrun" == "--dry-run" ]] && rsync_opts+=(--dry-run)
	rsync_opts+=(--max-size="${maxsize}")

	# Ensure remote directory exists
	echo "📁 Ensuring '$local_path' exists on $remote_host..."
	ssh "$remote_host" "mkdir -p '$local_path'"
	mkdir -p "$_brg_rundir"

	# to transfer small datasets
	local local_outdir="${local_dir}/${_brg_outdir}"
	local remote_outdir="${remote_host}:${local_outdir}"
	local main_outdir="${_local_host}:${local_outdir}"
	ssh "$remote_host" "mkdir -p '$local_outdir'"
	local rsync_opts2=(-az --info=progress2 --no-inc-recursive)
	rsync_opts2+=(--max-size="${maxsize}")

	# Perform sync
	case "$direction" in
	--pull)
		echo "Pulling from $remote_host:$local_path to $local_path..."
		rsync "${rsync_opts[@]}" "$remote_path/" "$local_path/"
		rsync "${rsync_opts2[@]}" "$remote_outdir/" "$_brg_outdir/"
		;;
	--push)
		echo "Pushing from $local_path to $remote_host:$local_path..."
		rsync "${rsync_opts[@]}" "$local_path/" "$remote_path/"
		rsync "${rsync_opts2[@]}" "$_brg_outdir/" "$remote_outdir/"
		;;
	--both)
		echo "Syncing both directions..."
		echo "Pulling from $remote_host..."
		rsync "${rsync_opts[@]}" "$remote_path/" "$local_path/"
		echo "Pushing to $remote_host..."
		rsync "${rsync_opts[@]}" "$local_path/" "$remote_path/"
		;;
	--main-push)
		echo "Pushing from $local_path to $_local_host:$local_path..."
		rsync "${rsync_opts[@]}" "$local_path/" "$main_path/"
		rsync "${rsync_opts2[@]}" "$_brg_outdir/" "$main_outdir/"
		;;
	*)
		echo "❌ Unknown direction: $direction"
		echo "Usage: $0 [--pull <folder >[index:0] | --push | --both | --main-push] [5M] [--dry-run]"
		exit 1
		;;
	esac

	echo "Done."
}

simulate_genus_species() {
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"

	source "${_POLAPLIB_DIR}/polap-variables-data.sh"

	mkdir -p "${_brg_target}"

	local _brg_rundir_ncbi_mtdna="${_brg_outdir_i}/ncbi-mtdna"
	local _brg_rundir_ncbi_ptdna="${_brg_outdir_i}/ncbi-ptdna"

	local platform="${_platform["$_brg_target"]}"
	if [[ "${platform}" == "PACBIO_SMRT" ]]; then
		${_polap_cmd} simulate hifi-reference \
			-p "${_brg_rundir_ncbi_ptdna}/ptdna-reference.fa" \
			-m "${_brg_rundir_ncbi_mtdna}/mtdna-reference.fa" \
			-o "${_brg_target}"
	elif [[ "${platform}" == "ONT" ]]; then
		${_polap_cmd} simulate ont-reference \
			-p "${_brg_rundir_ncbi_ptdna}/ptdna-reference.fa" \
			-m "${_brg_rundir_ncbi_mtdna}/mtdna-reference.fa" \
			-o "${_brg_target}"
	fi
}

readassemble-ont-pt_genus_species() {
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"

	local _brg_inum=0
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

	# local dimred="${_dimred["$_brg_target"]}"
	local platform="${_platform["$_brg_target"]}"
	local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

	rm -f "${_timing_txt}"
	rm -f "${_stdout_txt}"

	echo "------------------------------------------------------------------------------"
	echo "[INFO] Seed a DNA annotate organelle-genome assembly pipeline on ${_brg_outdir}-${_brg_sindex}"
	echo " - platform: ${platform}"

	mkdir -p "${_brg_rundir}"

	# prepare data:
	# version 1 and 2: from local
	# version 3: from -1.00x data
	data-long_genus_species "${_brg_outdir}" "${_brg_sindex}"

	# Start memory logger
	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"

	_polap_lib_conda-ensure_conda_env polap || exit 1

	local resolved_fastq=$(<"${_brg_outdir}/${_brg_inum}/l.fastq.path.txt")
	echo " - FASTQ: $resolved_fastq"

	local annotatedir="${_brg_rundir}"/annotate-read-pt

	local option_data_type="--nano-raw"
	if [[ "${platform}" == "PACBIO_SMRT" ]]; then
		_log_echo "Do: assemble ptDNA using seed reads with pacbio hifi: ${_brg_rundir}"
		option_data_type="--pacbio-hifi"
	elif [[ "${platform}" == "ONT" ]]; then
		_log_echo "Do: assemble ptDNA using seed reads with ONT raw: ${_brg_rundir}"
	fi

	${_polap_cmd} assemble-annotated-read \
		"${option_data_type}" \
		--plastid \
		-l "${resolved_fastq}" \
		-o "${_brg_rundir}"

	if [[ -s "${annotatedir}/pt/assembly_graph.gfa" ]]; then
		# extract ptDNA from the pt gfa
		${_polap_cmd} pt \
			--infile "${annotatedir}/pt/assembly_graph.gfa" \
			-o "${annotatedir}/pt/ptdna"

		${_polap_cmd} bandage png \
			"${annotatedir}/pt/assembly_graph.gfa" \
			"${annotatedir}/pt/assembly_graph.png"

		ln -sf "pt/ptdna/pt.0.fa" \
			"${annotatedir}/pt.0.fa"

		ln -sf "pt/assembly_graph.gfa" \
			"${annotatedir}/pt.0.gfa"

		ln -sf "pt/assembly_graph.png" \
			"${annotatedir}/pt.0.png"

		# NOTE: annotate for seeding
		# select connected components of the first contigs only
		# a new seeds menu
		${_polap_cmd} annotate \
			--quiet \
			-o "${annotatedir}" \
			-i pt

		# NOTE: plastid seed
		${_polap_cmd} seed-plastid \
			-o "${annotatedir}" \
			-i pt -j pt1
	else
		_log_echo "No PT assembly 0"
	fi

	# use all connected reads

	# NOTE: assemble-rate
	# repeat so that the sampling rate is between 0.1 and 0.5
	# map and subsample to compute the subsampling rate
	# we need the original input long-read data
	# we may need a new assemble command for this.
	# assemble-rate
	${_polap_cmd} assemble-rate \
		"${option_data_type}" \
		--plastid \
		-o "${annotatedir}" \
		-l "${resolved_fastq}" \
		-i pt -j pt1 -w 1500

	if [[ -s "${annotatedir}/pt1/assembly_graph.gfa" ]]; then
		${_polap_cmd} pt \
			--infile "${annotatedir}/pt1/assembly_graph.gfa" \
			-o "${annotatedir}/pt1/ptdna"

		${_polap_cmd} bandage png \
			"${annotatedir}/pt1/assembly_graph.gfa" \
			"${annotatedir}/pt1/assembly_graph.png"

		# using polap-py-find-plastid-gfa2fasta.py
		ln -sf "pt1/ptdna/pt.0.fa" \
			"${annotatedir}/pt.1.fa"

		# using polap-py-find-plastid-gfa.py
		ln -sf "pt1/ptdna/pt2.0.fa" \
			"${annotatedir}/pt2.1.fa"

		ln -sf "pt1/assembly_graph.gfa" \
			"${annotatedir}/pt.1.gfa"

		ln -sf "pt1/assembly_graph.png" \
			"${annotatedir}/pt.1.png"

		echo ptDNA assembly: "${annotatedir}/pt.1.gfa"
		echo ptDNA sequence: "${annotatedir}/pt.1.fa"
	else
		_log_echo "No PT assembly 1"
	fi

	conda deactivate

	# Summarize results after job (with previously defined summary function)
	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"
	_polap_lib_timing-get_system_info >>"${_timing_txt}"
}

# iterate the seed and assemble
readassemble-ont-pt-iterate_genus_species() {
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"

	local _brg_inum=0
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

	# local dimred="${_dimred["$_brg_target"]}"
	local platform="${_platform["$_brg_target"]}"
	local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

	rm -f "${_timing_txt}"
	rm -f "${_stdout_txt}"

	echo "------------------------------------------------------------------------------"
	echo "[INFO] Seed a DNA annotate organelle-genome assembly pipeline on ${_brg_outdir}-${_brg_sindex}"
	echo " - platform: ${platform}"

	mkdir -p "${_brg_rundir}"

	# prepare data:
	# version 1 and 2: from local
	# version 3: from -1.00x data
	data-long_genus_species "${_brg_outdir}" "${_brg_sindex}"

	# Start memory logger
	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"

	_polap_lib_conda-ensure_conda_env polap || exit 1

	local resolved_fastq=$(<"${_brg_outdir}/${_brg_inum}/l.fastq.path.txt")
	echo " we use FASTQ file ready: $resolved_fastq"

	local annotatedir="${_brg_rundir}"/annotate-read-pt

	local option_data_type="--nano-raw"
	if [[ "${platform}" == "PACBIO_SMRT" ]]; then
		_log_echo "Do: assemble ptDNA iteratively using seed reads with pacbio hifi: ${_brg_rundir}"
		option_data_type="--pacbio-hifi"
	elif [[ "${platform}" == "ONT" ]]; then
		_log_echo "Do: assemble ptDNA iteratively using seed reads with ONT raw: ${_brg_rundir}"
	fi

	if [[ -s "${annotatedir}/pt/assembly_graph.gfa" ]]; then
		_log_echo "found and use: ${annotatedir}/pt/assembly_graph.gfa"
		_log_echo "skip: assemble-annotated-read"
	else
		${_polap_cmd} assemble-annotated-read \
			"${option_data_type}" \
			--plastid \
			-l "${resolved_fastq}" \
			-o "${_brg_rundir}"
	fi

	if [[ -s "${annotatedir}/pt/assembly_graph.gfa" ]]; then
		${_polap_cmd} pt \
			--infile "${annotatedir}/pt/assembly_graph.gfa" \
			-o "${annotatedir}/pt/ptdna"

		${_polap_cmd} bandage png \
			"${annotatedir}/pt/assembly_graph.gfa" \
			"${annotatedir}/pt/assembly_graph.png"

		ln -sf "pt/ptdna/pt.0.fa" \
			"${annotatedir}/pt.0.fa"

		ln -sf "pt/assembly_graph.gfa" \
			"${annotatedir}/pt.0.gfa"

		ln -sf "pt/assembly_graph.png" \
			"${annotatedir}/pt.0.png"

	else
		_log_echo "No pt assembly 0"
	fi

	ln -s pt "${annotatedir}"/pt0

	local i
	for i in {0..1}; do
		local j=$((i + 1))

		# NOTE: annotate for seeding
		# select connected components of the pt contigs only
		${_polap_cmd} annotate \
			--quiet \
			-o "${annotatedir}" \
			-i pt$i

		# NOTE: mito seed
		${_polap_cmd} seed-plastid \
			-o "${annotatedir}" \
			-i pt$i -j pt$j

		if [[ ! -s "${annotatedir}/pt$i/mt.contig.name-pt$j" ]]; then
			_log_echo "No pt seed for pt$j"
			break
		fi

		# -l "${resolved_fastq}" \
		# -l "${annotatedir}"/pt.fq \
		${_polap_cmd} assemble-rate \
			"${option_data_type}" \
			--plastid \
			-o "${annotatedir}" \
			-l "${resolved_fastq}" \
			-i pt$i -j pt$j -w 1500

		if [[ -s "${annotatedir}/pt$j/assembly_graph.gfa" ]]; then
			${_polap_cmd} pt \
				--infile "${annotatedir}/pt$j/assembly_graph.gfa" \
				-o "${annotatedir}/pt$j/ptdna"

			${_polap_cmd} bandage png \
				"${annotatedir}/pt$j/assembly_graph.gfa" \
				"${annotatedir}/pt$j/assembly_graph.png"

			ln -sf "pt$j/ptdna/pt.0.fa" \
				"${annotatedir}/pt.$j.fa"

			ln -sf "pt$j/assembly_graph.gfa" \
				"${annotatedir}/pt.$j.gfa"

			ln -sf "pt$j/assembly_graph.png" \
				"${annotatedir}/pt.$j.png"

			echo ptDNA assembly: "${annotatedir}/pt.$j.gfa"
		else
			_log_echo "No pt assembly $j"
		fi

	done

	conda deactivate

	# Summarize results after job (with previously defined summary function)
	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"
	_polap_lib_timing-get_system_info >>"${_timing_txt}"
}

readassemble-ont-mt_genus_species() {
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"

	local _brg_inum=0
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

	# local dimred="${_dimred["$_brg_target"]}"
	local platform="${_platform["$_brg_target"]}"
	local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

	rm -f "${_timing_txt}"
	rm -f "${_stdout_txt}"

	echo "------------------------------------------------------------------------------"
	echo "[INFO] Seed a DNA annotate organelle-genome assembly pipeline on ${_brg_outdir}-${_brg_sindex}"
	echo " - platform: ${platform}"

	mkdir -p "${_brg_rundir}"

	# prepare data:
	# version 1 and 2: from local
	# version 3: from -1.00x data
	data-long_genus_species "${_brg_outdir}" "${_brg_sindex}"

	# Start memory logger
	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"

	_polap_lib_conda-ensure_conda_env polap || exit 1

	local resolved_fastq=$(<"${_brg_outdir}/${_brg_inum}/l.fastq.path.txt")
	echo " we use FASTQ file ready: $resolved_fastq"

	local annotatedir="${_brg_rundir}"/annotate-read-mt

	local option_data_type="--nano-raw"
	if [[ "${platform}" == "PACBIO_SMRT" ]]; then
		_log_echo "Do: assemble mtDNA using seed reads with pacbio hifi: ${_brg_rundir}"
		option_data_type="--pacbio-hifi"
	elif [[ "${platform}" == "ONT" ]]; then
		_log_echo "Do: assemble mtDNA using seed reads with ONT raw: ${_brg_rundir}"
	fi

	if [[ -s "${annotatedir}/mt/assembly_graph.gfa" ]]; then
		_log_echo "found and use: ${annotatedir}/mt/assembly_graph.gfa"
		_log_echo "skip: assemble-annotated-read"
	else
		${_polap_cmd} assemble-annotated-read \
			"${option_data_type}" \
			-l "${resolved_fastq}" \
			-o "${_brg_rundir}"
	fi

	if [[ -s "${annotatedir}/mt/assembly_graph.gfa" ]]; then
		${_polap_cmd} mt \
			--infile "${annotatedir}/mt/assembly_graph.gfa" \
			-o "${annotatedir}/mt/mtdna"

		${_polap_cmd} bandage png \
			"${annotatedir}/mt/assembly_graph.gfa" \
			"${annotatedir}/mt/assembly_graph.png"

		ln -sf "mt/assembly_graph.gfa" \
			"${annotatedir}/mt.0.gfa"

		ln -sf "mt/assembly_graph.png" \
			"${annotatedir}/mt.0.png"
	else
		_log_echo "No MT assembly 0"
	fi

	if [[ "${option_data_type}" == "--pacbio-hifi" ]]; then
		${_polap_cmd} filter reference hifi \
			"${option_data_type}" \
			-o "${annotatedir}" \
			-l "${resolved_fastq}" \
			--reference "${annotatedir_pt}"/pt.3.gfa
		local resolved_fastq="${annotatedir}/kmer/ref-filtered.fastq"
	fi

	# NOTE: annotate for seeding
	# select connected components of the MT contigs only
	${_polap_cmd} annotate \
		--quiet \
		-o "${annotatedir}" \
		-i mt

	# NOTE: mito seed
	${_polap_cmd} seed-mito \
		-o "${annotatedir}" \
		-i mt -j mt1

	# -l "${resolved_fastq}" \
	${_polap_cmd} assemble-rate \
		"${option_data_type}" \
		-o "${annotatedir}" \
		-l "${annotatedir}"/mt.fq \
		-i mt -j mt1 -w 3000

	if [[ -s "${annotatedir}/mt1/assembly_graph.gfa" ]]; then
		${_polap_cmd} bandage png \
			"${annotatedir}/mt1/assembly_graph.gfa" \
			"${annotatedir}/mt1/assembly_graph.png"

		ln -sf "mt1/assembly_graph.gfa" \
			"${annotatedir}/mt.1.gfa"

		ln -sf "mt1/assembly_graph.png" \
			"${annotatedir}/mt.1.png"

		echo mtDNA assembly: "${annotatedir}/mt.1.gfa"
	else
		_log_echo "No MT assembly 1"
	fi

	conda deactivate

	# Summarize results after job (with previously defined summary function)
	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"
	_polap_lib_timing-get_system_info >>"${_timing_txt}"
}

# iterate the seed and assemble
readassemble-ont-mt-iterate_genus_species() {
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"

	local _brg_inum=0
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

	# local dimred="${_dimred["$_brg_target"]}"
	local platform="${_platform["$_brg_target"]}"
	local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

	rm -f "${_timing_txt}"
	rm -f "${_stdout_txt}"

	echo "------------------------------------------------------------------------------"
	echo "[INFO] Seed a DNA annotate organelle-genome assembly pipeline on ${_brg_outdir}-${_brg_sindex}"
	echo " - platform: ${platform}"

	mkdir -p "${_brg_rundir}"

	# prepare data:
	# version 1 and 2: from local
	# version 3: from -1.00x data
	data-long_genus_species "${_brg_outdir}" "${_brg_sindex}"

	# Start memory logger
	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"

	_polap_lib_conda-ensure_conda_env polap || exit 1

	local resolved_fastq=$(<"${_brg_outdir}/${_brg_inum}/l.fastq.path.txt")
	echo " we use FASTQ file ready: $resolved_fastq"

	local annotatedir="${_brg_rundir}"/annotate-read-mt

	local option_data_type="--nano-raw"
	if [[ "${platform}" == "PACBIO_SMRT" ]]; then
		_log_echo "Do: assemble mtDNA iteratively using seed reads with pacbio hifi: ${_brg_rundir}"
		option_data_type="--pacbio-hifi"
	elif [[ "${platform}" == "ONT" ]]; then
		_log_echo "Do: assemble mtDNA iteratively using seed reads with ONT raw: ${_brg_rundir}"
	fi

	if [[ -s "${annotatedir}/mt/assembly_graph.gfa" ]]; then
		_log_echo "found and use: ${annotatedir}/mt/assembly_graph.gfa"
		_log_echo "skip: assemble-annotated-read"
	else
		${_polap_cmd} assemble-annotated-read \
			"${option_data_type}" \
			-l "${resolved_fastq}" \
			-o "${_brg_rundir}"
	fi

	if [[ -s "${annotatedir}/mt/assembly_graph.gfa" ]]; then
		${_polap_cmd} mt \
			--infile "${annotatedir}/mt/assembly_graph.gfa" \
			-o "${annotatedir}/mt/mtdna"

		${_polap_cmd} bandage png \
			"${annotatedir}/mt/assembly_graph.gfa" \
			"${annotatedir}/mt/assembly_graph.png"

		ln -sf "mt/assembly_graph.gfa" \
			"${annotatedir}/mt.0.gfa"

		ln -sf "mt/assembly_graph.png" \
			"${annotatedir}/mt.0.png"

	else
		_log_echo "No MT assembly 0"
	fi

	if [[ "${option_data_type}" == "--pacbio-hifi" ]]; then
		${_polap_cmd} filter reference hifi \
			"${option_data_type}" \
			-o "${annotatedir}" \
			-l "${resolved_fastq}" \
			--reference "${annotatedir_pt}"/pt.3.gfa
		local resolved_fastq="${annotatedir}/kmer/ref-filtered.fastq"
	fi

	ln -s mt "${annotatedir}"/mt0

	local i
	for i in {0..5}; do
		local j=$((i + 1))

		# NOTE: annotate for seeding
		# select connected components of the MT contigs only
		${_polap_cmd} annotate \
			--quiet \
			-o "${annotatedir}" \
			-i mt$i

		# NOTE: mito seed
		${_polap_cmd} seed-mito \
			-o "${annotatedir}" \
			-i mt$i -j mt$j

		if [[ ! -s "${annotatedir}/mt$i/mt.contig.name-mt$j" ]]; then
			_log_echo "No MT seed for mt$j"
			break
		fi

		# -l "${resolved_fastq}" \
		${_polap_cmd} assemble-rate \
			"${option_data_type}" \
			-o "${annotatedir}" \
			-l "${annotatedir}"/mt.fq \
			-i mt$i -j mt$j -w 3000

		if [[ -s "${annotatedir}/mt$j/assembly_graph.gfa" ]]; then
			${_polap_cmd} bandage png \
				"${annotatedir}/mt$j/assembly_graph.gfa" \
				"${annotatedir}/mt$j/assembly_graph.png"

			ln -sf "mt$j/assembly_graph.gfa" \
				"${annotatedir}/mt.$j.gfa"

			ln -sf "mt$j/assembly_graph.png" \
				"${annotatedir}/mt.$j.png"

			echo mtDNA assembly: "${annotatedir}/mt.$j.gfa"
		else
			_log_echo "No MT assembly $j"
		fi

	done

	conda deactivate

	# Summarize results after job (with previously defined summary function)
	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"
	_polap_lib_timing-get_system_info >>"${_timing_txt}"
}

readassemble-ont-nt_genus_species() {
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"

	local _brg_inum=0
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

	# local dimred="${_dimred["$_brg_target"]}"
	local platform="${_platform["$_brg_target"]}"
	local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

	rm -f "${_timing_txt}"
	rm -f "${_stdout_txt}"

	echo "------------------------------------------------------------------------------"
	echo "[INFO] Seed a DNA annotate organelle-genome assembly pipeline on ${_brg_outdir}-${_brg_sindex}"
	echo " - platform: ${platform}"

	mkdir -p "${_brg_rundir}"

	# prepare data:
	# version 1 and 2: from local
	# version 3: from -1.00x data
	data-long_genus_species "${_brg_outdir}" "${_brg_sindex}"
	if [[ -s "${long_sra}-10x".fq ]]; then
		echo "${long_sra}-10x".fq >"${_brg_outdir}/${_brg_inum}/l.fastq.path.txt"
	fi

	# Start memory logger
	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"

	_polap_lib_conda-ensure_conda_env polap || exit 1

	local resolved_fastq=$(<"${_brg_outdir}/${_brg_inum}/l.fastq.path.txt")
	echo " we use FASTQ file ready: $resolved_fastq"

	local annotatedir="${_brg_rundir}"/annotate-read-nt
	local annotatedir_pt="${_brg_rundir}"/annotate-read-pt

	local option_data_type="--nano-raw"
	if [[ "${platform}" == "PACBIO_SMRT" ]]; then
		_log_echo "Do: assemble mtDNA including noncoding seed reads with pacbio hifi: ${_brg_rundir}"
		option_data_type="--pacbio-hifi"
	elif [[ "${platform}" == "ONT" ]]; then
		_log_echo "Do: assemble mtDNA including noncoding seed reads with ONT raw: ${_brg_rundir}"
	fi

	if [[ -s "${annotatedir}/mt/assembly_graph.gfa" ]]; then
		_log_echo "found and use: ${annotatedir}/mt/assembly_graph.gfa"
		_log_echo "skip: assemble-annotated-read"
	else
		# ${_polap_cmd} assemble-annotated-read \
		${_polap_cmd} assemble-annotated-read-nt \
			"${option_data_type}" \
			-l "${resolved_fastq}" \
			-o "${_brg_rundir}"
	fi

	if [[ -s "${annotatedir}/mt/assembly_graph.gfa" ]]; then
		${_polap_cmd} mt \
			--infile "${annotatedir}/mt/assembly_graph.gfa" \
			-o "${annotatedir}/mt/mtdna"

		${_polap_cmd} bandage png \
			"${annotatedir}/mt/assembly_graph.gfa" \
			"${annotatedir}/mt/assembly_graph.png"

		ln -sf "mt/assembly_graph.gfa" \
			"${annotatedir}/mt.0.gfa"

		ln -sf "mt/assembly_graph.png" \
			"${annotatedir}/mt.0.png"
	else
		_log_echo "No MT assembly 0"
	fi

	if [[ "${option_data_type}" == "--pacbio-hifi" ]]; then
		${_polap_cmd} filter reference hifi \
			"${option_data_type}" \
			-o "${annotatedir}" \
			-l "${resolved_fastq}" \
			--reference "${annotatedir_pt}"/pt.3.gfa
		local resolved_fastq="${annotatedir}/kmer/ref-filtered.fastq"
	fi

	# NOTE: annotate for seeding
	# select connected components of the MT contigs only
	${_polap_cmd} annotate \
		--quiet \
		-o "${annotatedir}" \
		-i mt

	# NOTE: mito seed
	${_polap_cmd} seed-mito \
		-o "${annotatedir}" \
		-i mt -j mt1

	# -l "${resolved_fastq}" \
	# -l "${annotatedir}"/mt.fq \
	if [[ "${platform}" == "PACBIO_SMRT" ]]; then
		${_polap_cmd} assemble-rate \
			"${option_data_type}" \
			-o "${annotatedir}" \
			-l "${resolved_fastq}" \
			--reference "${annotatedir_pt}"/pt.3.gfa \
			-i mt -j mt1 -w 3000
	elif [[ "${platform}" == "ONT" ]]; then
		${_polap_cmd} assemble-rate \
			"${option_data_type}" \
			-o "${annotatedir}" \
			-l "${annotatedir}"/mt.fq \
			-i mt -j mt1 -w 3000
	fi

	if [[ -s "${annotatedir}/mt1/assembly_graph.gfa" ]]; then
		${_polap_cmd} bandage png \
			"${annotatedir}/mt1/assembly_graph.gfa" \
			"${annotatedir}/mt1/assembly_graph.png"

		ln -sf "mt1/assembly_graph.gfa" \
			"${annotatedir}/mt.1.gfa"

		ln -sf "mt1/assembly_graph.png" \
			"${annotatedir}/mt.1.png"

		echo mtDNA assembly: "${annotatedir}/mt.1.gfa"
	else
		_log_echo "No MT assembly 1"
	fi

	conda deactivate

	# Summarize results after job (with previously defined summary function)
	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"
	_polap_lib_timing-get_system_info >>"${_timing_txt}"
}

# iterate the seed and assemble
readassemble-ont-nt-iterate_genus_species() {
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"

	local _brg_inum=0
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

	# local dimred="${_dimred["$_brg_target"]}"
	local platform="${_platform["$_brg_target"]}"
	local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

	rm -f "${_timing_txt}"
	rm -f "${_stdout_txt}"

	echo "------------------------------------------------------------------------------"
	echo "[INFO] Seed a DNA annotate organelle-genome assembly pipeline on ${_brg_outdir}-${_brg_sindex}"
	echo " - platform: ${platform}"

	mkdir -p "${_brg_rundir}"

	# prepare data:
	# version 1 and 2: from local
	# version 3: from -1.00x data
	data-long_genus_species "${_brg_outdir}" "${_brg_sindex}"
	if [[ -s "${long_sra}-10x".fq ]]; then
		echo "${long_sra}-10x".fq >"${_brg_outdir}/${_brg_inum}/l.fastq.path.txt"
	fi

	# Start memory logger
	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"

	_polap_lib_conda-ensure_conda_env polap || exit 1

	local resolved_fastq=$(<"${_brg_outdir}/${_brg_inum}/l.fastq.path.txt")
	echo " we use FASTQ file ready: $resolved_fastq"

	local annotatedir="${_brg_rundir}"/annotate-read-nt
	local annotatedir_pt="${_brg_rundir}"/annotate-read-pt

	local option_data_type="--nano-raw"
	if [[ "${platform}" == "PACBIO_SMRT" ]]; then
		_log_echo "Do: assemble mtDNA iteratively including noncoding seed reads with pacbio hifi: ${_brg_rundir}"
		option_data_type="--pacbio-hifi"
	elif [[ "${platform}" == "ONT" ]]; then
		_log_echo "Do: assemble mtDNA iteratively including noncoding seed reads with ONT raw: ${_brg_rundir}"
	fi

	if [[ -s "${annotatedir}/mt/assembly_graph.gfa" ]]; then
		_log_echo "found and use: ${annotatedir}/mt/assembly_graph.gfa"
		_log_echo "skip: assemble-annotated-read"
	else
		${_polap_cmd} assemble-annotated-read-nt \
			"${option_data_type}" \
			-l "${resolved_fastq}" \
			-o "${_brg_rundir}"
	fi

	if [[ -s "${annotatedir}/mt/assembly_graph.gfa" ]]; then
		${_polap_cmd} mt \
			--infile "${annotatedir}/mt/assembly_graph.gfa" \
			-o "${annotatedir}/mt/mtdna"

		${_polap_cmd} bandage png \
			"${annotatedir}/mt/assembly_graph.gfa" \
			"${annotatedir}/mt/assembly_graph.png"

		ln -sf "mt/assembly_graph.gfa" \
			"${annotatedir}/mt.0.gfa"

		ln -sf "mt/assembly_graph.png" \
			"${annotatedir}/mt.0.png"

	else
		_log_echo "No MT assembly 0"
	fi

	ln -s mt "${annotatedir}"/mt0

	# TODO: remove pt reads from the input if not --plastid
	if [[ "${option_data_type}" == "--pacbio-hifi" ]]; then
		${_polap_cmd} filter reference hifi \
			"${option_data_type}" \
			-o "${annotatedir}" \
			-l "${resolved_fastq}" \
			--reference "${annotatedir_pt}"/pt.3.gfa
		local resolved_fastq="${annotatedir}/kmer/ref-filtered.fastq"
	fi

	local i
	for i in {0..5}; do
		local j=$((i + 1))

		# NOTE: annotate for seeding
		# select connected components of the MT contigs only
		${_polap_cmd} annotate \
			--quiet \
			-o "${annotatedir}" \
			-i mt$i

		# NOTE: mito seed
		${_polap_cmd} seed-mito \
			-o "${annotatedir}" \
			-i mt$i -j mt$j

		if [[ ! -s "${annotatedir}/mt$i/mt.contig.name-mt$j" ]]; then
			_log_echo "No MT seed for mt$j"
			break
		fi

		# -l "${resolved_fastq}" \
		# -l "${annotatedir}"/mt.fq \
		if [[ "${platform}" == "PACBIO_SMRT" ]]; then
			${_polap_cmd} assemble-rate \
				"${option_data_type}" \
				-o "${annotatedir}" \
				-l "${resolved_fastq}" \
				--reference "${annotatedir_pt}"/pt.3.gfa \
				-i mt$i -j mt$j -w 3000
		elif [[ "${platform}" == "ONT" ]]; then
			${_polap_cmd} assemble-rate \
				"${option_data_type}" \
				-o "${annotatedir}" \
				-l "${annotatedir}"/mt.fq \
				-i mt$i -j mt$j -w 3000
		fi

		if [[ -s "${annotatedir}/mt$j/assembly_graph.gfa" ]]; then
			${_polap_cmd} bandage png \
				"${annotatedir}/mt$j/assembly_graph.gfa" \
				"${annotatedir}/mt$j/assembly_graph.png"

			ln -sf "mt$j/assembly_graph.gfa" \
				"${annotatedir}/mt.$j.gfa"

			ln -sf "mt$j/assembly_graph.png" \
				"${annotatedir}/mt.$j.png"

			echo mtDNA assembly: "${annotatedir}/mt.$j.gfa"
		else
			_log_echo "No MT assembly $j"
		fi

	done

	conda deactivate

	# Summarize results after job (with previously defined summary function)
	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"
	_polap_lib_timing-get_system_info >>"${_timing_txt}"
}

# animal mtDNA assembly
# copy of run-assemble-ont-pt
readassemble-ont-animal-mt_genus_species() {
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"

	local _brg_inum=0
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

	# local dimred="${_dimred["$_brg_target"]}"
	local platform="${_platform["$_brg_target"]}"
	local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

	rm -f "${_timing_txt}"
	rm -f "${_stdout_txt}"

	echo "="
	echo "------------------------------------------------------------------------------"
	echo "[INFO] Seed a DNA annotate organelle-genome assembly pipeline on ${_brg_outdir}-${_brg_sindex}"
	echo " - platform: ${platform}"

	mkdir -p "${_brg_rundir}"

	# prepare data:
	# version 1 and 2: from local
	# version 3: from -1.00x data
	data-long_genus_species "${_brg_outdir}" "${_brg_sindex}"

	# Start memory logger
	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"

	_polap_lib_conda-ensure_conda_env polap || exit 1

	local resolved_fastq=$(<"${_brg_outdir}/${_brg_inum}/l.fastq.path.txt")
	echo " we use FASTQ file ready: $resolved_fastq"

	local annotatedir="${_brg_rundir}"/annotate-read-mt

	if [[ "${platform}" == "PACBIO_SMRT" ]]; then
		${_polap_cmd} flye1 --pacbio-hifi \
			-l "${_brg_rundir}"/lk.fq.gz \
			-g $genomesize \
			-o "${_brg_rundir}"
	elif [[ "${platform}" == "ONT" ]]; then
		echo "assemble using seed reads for mtDNA: ${_brg_rundir}"

		${_polap_cmd} assemble-annotated-read \
			--animal \
			-l "${resolved_fastq}" \
			-o "${_brg_rundir}"

		# extract mtDNA from the mt gfa
		${_polap_cmd} mt \
			--infile "${annotatedir}/mt/assembly_graph.gfa" \
			-o "${annotatedir}/mt/mtdna"

		${_polap_cmd} bandage png \
			"${annotatedir}/mt/assembly_graph.gfa" \
			"${annotatedir}/mt/assembly_graph.png"

		ln -sf "mt/mtdna/mt.0.fa" \
			"${annotatedir}/mt.0.fa"

		ln -sf "mt/assembly_graph.gfa" \
			"${annotatedir}/mt.0.gfa"

		ln -sf "mt/assembly_graph.png" \
			"${annotatedir}/mt.0.png"

		# NOTE: annotate for seeding
		# select connected components of the first contigs only
		# a new seeds menu
		${_polap_cmd} annotate \
			--quiet \
			-o "${annotatedir}" \
			-i mt

		# NOTE: plastid seed
		${_polap_cmd} seed-plastid \
			-o "${annotatedir}" \
			-i mt -j mt1

		# NOTE: assemble-rate
		# repeat so that the sampling rate is between 0.1 and 0.5
		# map and subsample to compute the subsampling rate
		# we need the original input long-read data
		# we may need a new assemble command for this.
		# assemble-rate
		${_polap_cmd} assemble-rate \
			-o "${annotatedir}" \
			-l "${resolved_fastq}" \
			-i mt -j mt1 -w 3000

		${_polap_cmd} mt \
			--infile "${annotatedir}/mt1/assembly_graph.gfa" \
			-o "${annotatedir}/mt1/mtdna"

		${_polap_cmd} bandage png \
			"${annotatedir}/mt1/assembly_graph.gfa" \
			"${annotatedir}/mt1/assembly_graph.png"

		ln -sf "mt1/mtdna/mt.0.fa" \
			"${annotatedir}/mt.1.fa"

		ln -sf "mt1/assembly_graph.gfa" \
			"${annotatedir}/mt.1.gfa"

		ln -sf "mt1/assembly_graph.png" \
			"${annotatedir}/mt.1.png"

		echo mtDNA assembly: "${annotatedir}/mt.1.gfa"
		echo mtDNA sequence: "${annotatedir}/mt.1.fa"
	fi

	conda deactivate

	# Summarize results after job (with previously defined summary function)
	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"
	_polap_lib_timing-get_system_info >>"${_timing_txt}"
}

readassemble-annotate-pt_genus_species() {
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"

	local _brg_inum=0
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

	# local dimred="${_dimred["$_brg_target"]}"
	local platform="${_platform["$_brg_target"]}"
	local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

	rm -f "${_timing_txt}"
	rm -f "${_stdout_txt}"

	echo "------------------------------------------------------------------------------"
	echo "[INFO] Start a DNA annotate organelle-genome assembly pipeline on ${_brg_outdir}-${_brg_sindex}"
	echo " - platform: ${platform}"

	mkdir -p "${_brg_rundir}"

	# prepare data:
	# version 1 and 2: from local
	# version 3: from -1.00x data
	data-long_genus_species "${_brg_outdir}" "${_brg_sindex}"

	# Start memory logger
	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"

	_polap_lib_conda-ensure_conda_env polap || exit 1

	local resolved_fastq=$(<"${_brg_outdir}/${_brg_inum}/l.fastq.path.txt")
	echo " - FASTQ: $resolved_fastq"

	rm -rf "${_brg_rundir}/annotate-read-pt"

	if [[ "${platform}" == "PACBIO_SMRT" ]]; then
		${_polap_cmd} annotate-read-pt \
			--pacbio-hifi \
			-l "${resolved_fastq}" \
			--annotate-read-min-mapq 1 \
			--annotate-read-min-identity 0.15 \
			-o "${_brg_rundir}"
	elif [[ "${platform}" == "ONT" ]]; then
		${_polap_cmd} annotate-read-pt \
			-l "${resolved_fastq}" \
			--annotate-read-min-mapq 1 \
			--annotate-read-min-identity 0.15 \
			-o "${_brg_rundir}"
	fi

	conda deactivate

	# Summarize results after job (with previously defined summary function)
	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"
	_polap_lib_timing-get_system_info >>"${_timing_txt}"
}

readassemble-annotate-mt_genus_species() {
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"

	local _brg_inum=0
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

	# local dimred="${_dimred["$_brg_target"]}"
	local platform="${_platform["$_brg_target"]}"
	local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

	rm -f "${_timing_txt}"
	rm -f "${_stdout_txt}"

	echo "------------------------------------------------------------------------------"
	echo "[INFO] Annotate ONT reads on ${_brg_outdir}-${_brg_sindex}"
	echo " - platform: ${platform}"

	mkdir -p "${_brg_rundir}"

	# prepare data:
	# version 1 and 2: from local
	# version 3: from -1.00x data
	data-long_genus_species "${_brg_outdir}" "${_brg_sindex}"

	# Start memory logger
	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"

	_polap_lib_conda-ensure_conda_env polap || exit 1

	local resolved_fastq=$(<"${_brg_outdir}/${_brg_inum}/l.fastq.path.txt")
	echo " we use FASTQ file: $resolved_fastq"

	echo " delete and start with an empty folder: ${_brg_rundir}/annotate-read-mt"
	rm -rf "${_brg_rundir}/annotate-read-mt"

	if [[ "${platform}" == "PACBIO_SMRT" ]]; then
		${_polap_cmd} annotate-read-mt \
			--pacbio-hifi \
			-l "${resolved_fastq}" \
			--annotate-read-min-mapq 1 \
			--annotate-read-min-identity 0.15 \
			-o "${_brg_rundir}"
	elif [[ "${platform}" == "ONT" ]]; then
		${_polap_cmd} annotate-read-mt \
			-l "${resolved_fastq}" \
			--annotate-read-min-mapq 1 \
			--annotate-read-min-identity 0.15 \
			-o "${_brg_rundir}"
	fi

	conda deactivate

	# Summarize results after job (with previously defined summary function)
	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"
	_polap_lib_timing-get_system_info >>"${_timing_txt}"
}

readassemble-annotate-nt_genus_species() {
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"

	local _brg_inum=0
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

	# local dimred="${_dimred["$_brg_target"]}"
	local platform="${_platform["$_brg_target"]}"
	local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

	rm -f "${_timing_txt}"
	rm -f "${_stdout_txt}"

	echo "------------------------------------------------------------------------------"
	echo "[INFO] Annotate ONT reads on ${_brg_outdir}-${_brg_sindex}"
	echo " - platform: ${platform}"

	mkdir -p "${_brg_rundir}"

	# prepare data:
	# version 1 and 2: from local
	# version 3: from -1.00x data
	data-long_genus_species "${_brg_outdir}" "${_brg_sindex}"
	if [[ -s "${long_sra}-10x".fq ]]; then
		echo "${long_sra}-10x".fq >"${_brg_outdir}/${_brg_inum}/l.fastq.path.txt"
	fi

	# Start memory logger
	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"

	_polap_lib_conda-ensure_conda_env polap || exit 1

	local resolved_fastq=$(<"${_brg_outdir}/${_brg_inum}/l.fastq.path.txt")
	echo " we use FASTQ file: $resolved_fastq"

	echo " delete and start with an empty folder: ${_brg_rundir}/annotate-read-nt"
	rm -rf "${_brg_rundir}/annotate-read-nt"

	if [[ "${platform}" == "PACBIO_SMRT" ]]; then
		${_polap_cmd} annotate-read-nt \
			--pacbio-hifi \
			-l "${resolved_fastq}" \
			--annotate-read-min-mapq 1 \
			--annotate-read-min-identity 0.15 \
			-o "${_brg_rundir}"
	elif [[ "${platform}" == "ONT" ]]; then
		${_polap_cmd} annotate-read-nt \
			-l "${resolved_fastq}" \
			--annotate-read-min-mapq 1 \
			--annotate-read-min-identity 0.15 \
			-o "${_brg_rundir}"
	fi

	conda deactivate

	# Summarize results after job (with previously defined summary function)
	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"
	_polap_lib_timing-get_system_info >>"${_timing_txt}"
}

run-pca-oga_genus_species() {
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

	local min_read="${_min_read["$_brg_target"]}"
	local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

	rm -f "${_timing_txt}"
	rm -f "${_stdout_txt}"

	echo "="
	echo "------------------------------------------------------------------------------"
	echo "[INFO] Start an whole-genome HiFi assembly pipeline on ${_brg_outdir}-${_brg_sindex}"

	mkdir -p "${_brg_rundir}"

	# prepare data:
	# version 1 and 2: from local
	# version 3: from -1.00x data
	data-long_genus_species "${_brg_outdir}" "${_brg_sindex}"

	# Start memory logger
	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"

	_polap_lib_conda-ensure_conda_env polap || exit 1

	local resolved_fastq=$(<"${_brg_outdir}/${_brg_inum}/l.fastq.path.txt")
	echo "here FASTQ file ready: $resolved_fastq"

	# extract ptDNA from the long-read polish
	${_polap_cmd} search gfa plastid \
		-i "${_brg_rundir}"/0/assembly_graph.gfa \
		--outfile "${_brg_rundir}"

	# identify ptDNA from the 30-contigger
	${_polap_cmd} search gfa plastid \
		-i "${_brg_rundir}"/0/30-contiggre/graphy_final.gfa \
		--outfile "${_brg_rundir}"

	# select mtDNA connected component
	# mtDNA gfa -> mt.0.fa -> mt.1.fa
	${_polap_cmd} assemble2 \
		-o "${_brg_rundir}"

	conda deactivate

	# Summarize results after job (with previously defined summary function)
	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"
	_polap_lib_timing-get_system_info >>"${_timing_txt}"
}

run-pca-wga_genus_species() {
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"
	local _brg_inum="${3:-0}"
	local _brg_knum="${4:-6}"

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

	local dimred="${_dimred["$_brg_target"]}"
	local readsel="${_readsel["$_brg_target"]}"
	local min_read="${_min_read["$_brg_target"]}"
	local platform="${_platform["$_brg_target"]}"
	local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

	rm -f "${_timing_txt}"
	rm -f "${_stdout_txt}"

	echo "="
	echo "------------------------------------------------------------------------------"
	echo "[INFO] Start an whole-genome HiFi assembly pipeline on ${_brg_outdir}-${_brg_sindex}"
	echo " - platform: ${platform}"
	echo " - dimension reduction: ${dimred}"
	echo " - read selection: ${readsel}"

	mkdir -p "${_brg_rundir}"

	# prepare data:
	# version 1 and 2: from local
	# version 3: from -1.00x data
	data-long_genus_species "${_brg_outdir}" "${_brg_sindex}"

	# Start memory logger
	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"

	_polap_lib_conda-ensure_conda_env polap || exit 1

	local resolved_fastq=$(<"${_brg_outdir}/${_brg_inum}/l.fastq.path.txt")
	echo "here FASTQ file ready: $resolved_fastq"

	# ${_polap_cmd} init \
	# 	-o "${_brg_rundir}"
	#
	# ${_polap_cmd} filter pca \
	# 	-l "${resolved_fastq}" \
	# 	-o "${_brg_rundir}"

	# o/kmer/k2.pca.tsv
	# o/kmer/k2.umap.tsv
	local K="${_brg_knum}"
	if [[ "$K" == "0" ]]; then
		K=2
	fi

	local organelle_fastq="${resolved_fastq}".k${K}.deselected.fastq
	local organelle_read_names="${_brg_rundir}"/kmer/k${K}.${dimred}.deselected.tsv
	if [[ "${readsel}" == "ellipse" ]]; then
		organelle_read_names="${_brg_rundir}"/kmer/k${K}.${dimred}.ellipse.deselected.tsv
	fi

	cut -f5 "${organelle_read_names}" |
		seqtk subseq "${resolved_fastq}" - \
			>"${organelle_fastq}"

	# genome size estimate
	# pl total-length-long -l DRR503528.fastq
	# pl find-genome-size-for-pacbio -l DRR503528.fastq
	# pl reduce-data -l DRR503528.fastq
	# pl flye1 --pacbio-hifi --blast --blast-mt 0.0001 --blast-pt 0.0001

	${_polap_cmd} total-length-long \
		-l "${organelle_fastq}" \
		-o "${_brg_rundir}"

	if [[ "${platform}" == "PACBIO_SMRT" ]]; then
		${_polap_cmd} find-genome-size-for-pacbio \
			-l "${organelle_fastq}" \
			-o "${_brg_rundir}"
	elif [[ "${platform}" == "ONT" ]]; then
		echo "ONT: fixed genome size: 3Mb"
		echo 3000000 >"${_brg_rundir}"/short_expected_genome_size.txt
		# ${_polap_cmd} find-genome-size-for-pacbio \
		# 	-l "${organelle_fastq}" \
		# 	-o "${_brg_rundir}"
	fi

	${_polap_cmd} reduce-data \
		-l "${organelle_fastq}" \
		-o "${_brg_rundir}"

	mkdir -p "${_brg_rundir}/pt"
	if [[ "${platform}" == "PACBIO_SMRT" ]]; then
		${_polap_cmd} assemble1 --pacbio-hifi \
			-l "${organelle_fastq}" \
			-o "${_brg_rundir}/pt"
		# flye --pacbio-hifi "${resolved_fastq}".k${K}.pca.deselected.fastq \
		# 	-t "${_brg_threads}" --asm-coverage 30 -g 3m \
		# 	--out-dir "${_brg_rundir}"/0
	elif [[ "${platform}" == "ONT" ]]; then
		# ${_polap_cmd} assemble1 --nano-raw \
		# 	-l "${organelle_fastq}" \
		# 	-o "${_brg_rundir}/pt"
		flye --nano-raw "${organelle_fastq}" \
			-t "${_brg_threads}" --asm-coverage 30 -g 3m \
			--out-dir "${_brg_rundir}"/pt/0
	fi

	local genomesize=$(<"${_brg_rundir}"/short_expected_genome_size.txt)
	if ((genomesize < 1000000)); then
		genomesize=1000000
	fi

	if [[ "${platform}" == "PACBIO_SMRT" ]]; then
		${_polap_cmd} flye1 --pacbio-hifi \
			-l "${_brg_rundir}"/lk.fq.gz \
			-g $genomesize \
			-o "${_brg_rundir}"
		# flye --pacbio-hifi "${resolved_fastq}".k${K}.pca.deselected.fastq \
		# 	-t "${_brg_threads}" --asm-coverage 30 -g 3m \
		# 	--out-dir "${_brg_rundir}"/0
	elif [[ "${platform}" == "ONT" ]]; then
		${_polap_cmd} flye1 --nano-raw \
			-l "${_brg_rundir}"/lk.fq.gz \
			-g $genomesize \
			-o "${_brg_rundir}"
		# flye --nano-raw "${resolved_fastq}".k${K}.pca.deselected.fastq \
		# 	-t "${_brg_threads}" --asm-coverage 30 -g 3m \
		# 	--out-dir "${_brg_rundir}"/0
	fi

	# identify ptDNA gfa -> pt.0.fa -> pt.1.fa
	# _disassemble-step11 ~ step14
	# 1. annotate
	# 2. seeds
	# 3. polap-py-find-plastid-gfa2fasta.py
	# 4. select one of the 4: polap-py-compare2ptdna.py
	#
	# remove ptDNA gfa from the wga gfa
	#
	# create mt.contig.seed-1 for mtDNA
	# oga for mtDNA
	#
	# select mtDNA connected component
	# mtDNA gfa -> mt.0.fa -> mt.1.fa

	conda deactivate

	# Summarize results after job (with previously defined summary function)
	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"
	_polap_lib_timing-get_system_info >>"${_timing_txt}"
}

run-pca-pt_genus_species() {
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"
	local _brg_inum="${3:-0}"
	local _brg_knum="${4:-5}"
	local _brg_genomesize="${5:-0}"

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

	local dimred="${_dimred["$_brg_target"]}"
	local readsel="${_readsel["$_brg_target"]}"
	local min_read="${_min_read["$_brg_target"]}"
	local platform="${_platform["$_brg_target"]}"
	local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

	rm -f "${_timing_txt}"
	rm -f "${_stdout_txt}"

	echo "="
	echo "------------------------------------------------------------------------------"
	echo "[INFO] Start an whole-genome HiFi assembly pipeline on ${_brg_outdir}-${_brg_sindex}"
	echo " - platform: ${platform}"
	echo " - dimension reduction: ${dimred}"
	echo " - read selection: ${readsel}"

	mkdir -p "${_brg_rundir}"

	# prepare data:
	# version 1 and 2: from local
	# version 3: from -1.00x data
	data-long_genus_species "${_brg_outdir}" "${_brg_sindex}"

	# Start memory logger
	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"

	_polap_lib_conda-ensure_conda_env polap || exit 1

	local resolved_fastq=$(<"${_brg_outdir}/${_brg_inum}/l.fastq.path.txt")
	echo "here FASTQ file ready: $resolved_fastq"

	# ${_polap_cmd} init \
	# 	-o "${_brg_rundir}"
	#
	# ${_polap_cmd} filter pca \
	# 	-l "${resolved_fastq}" \
	# 	-o "${_brg_rundir}"

	# o/kmer/k2.pca.tsv
	# o/kmer/k2.umap.tsv
	local K="${_brg_knum}"
	if [[ "$K" == "0" ]]; then
		K=2
	fi

	local organelle_fastq="${resolved_fastq}".k${K}.deselected.fastq
	local organelle_read_names="${_brg_rundir}"/kmer/k${K}.${dimred}.deselected.tsv
	if [[ "${readsel}" == "ellipse" ]]; then
		organelle_read_names="${_brg_rundir}"/kmer/k${K}.${dimred}.ellipse.deselected.tsv
	fi

	# if [[ ! -s "${organelle_fastq}" ]]; then
	cut -f5 "${organelle_read_names}" |
		seqtk subseq "${resolved_fastq}" - \
			>"${organelle_fastq}"
	# fi

	# genome size estimate
	# pl total-length-long -l DRR503528.fastq
	# pl find-genome-size-for-pacbio -l DRR503528.fastq
	# pl reduce-data -l DRR503528.fastq
	# pl flye1 --pacbio-hifi --blast --blast-mt 0.0001 --blast-pt 0.0001

	${_polap_cmd} total-length-long \
		-l "${organelle_fastq}" \
		-o "${_brg_rundir}"

	if [[ "${platform}" == "PACBIO_SMRT" ]]; then
		${_polap_cmd} find-genome-size-for-pacbio \
			-l "${organelle_fastq}" \
			-o "${_brg_rundir}"
	elif [[ "${platform}" == "ONT" ]]; then

		if [[ "${_brg_genomesize}" == "0" ]]; then
			# local long_len=$(<"${_brg_rundir}"/long_total_length.txt)
			# local reduced_len=$((long_len / 50))
			# echo $reduced_len >"${_brg_rundir}"/short_expected_genome_size.txt
			echo "ONT: genome size: auto flye 2.9.6 is used"
		else
			echo "ONT: genome size: auto flye 2.9.6 is used"
			# echo "${_brg_genomesize}" >"${_brg_rundir}"/short_expected_genome_size.txt

		fi
	fi

	# ${_polap_cmd} reduce-data \
	# 	-l "${organelle_fastq}" \
	# 	-o "${_brg_rundir}"

	# local genomesize=$(<"${_brg_rundir}"/short_expected_genome_size.txt)

	if [[ "${platform}" == "PACBIO_SMRT" ]]; then
		${_polap_cmd} flye1 --pacbio-hifi \
			-l "${_brg_rundir}"/lk.fq.gz \
			-g $genomesize \
			-o "${_brg_rundir}"
	elif [[ "${platform}" == "ONT" ]]; then
		# ${_polap_cmd} flye1 --nano-raw \
		# 	-g $genomesize \
		# 	-o "${_brg_rundir}"
		_polap_lib_conda-ensure_conda_env polap-flye || exit 1
		flye --nano-raw "${organelle_fastq}" \
			-t "${_brg_threads}" \
			--out-dir "${_brg_rundir}"/0
		conda deactivate
	fi

	# identify ptDNA gfa -> pt.0.fa -> pt.1.fa
	# _disassemble-step11 ~ step14
	# 1. annotate
	# 2. seeds
	# 3. polap-py-find-plastid-gfa2fasta.py
	# 4. select one of the 4: polap-py-compare2ptdna.py
	#
	# remove ptDNA gfa from the wga gfa
	#
	# create mt.contig.seed-1 for mtDNA
	# oga for mtDNA
	#
	# select mtDNA connected component
	# mtDNA gfa -> mt.0.fa -> mt.1.fa

	conda deactivate

	# Summarize results after job (with previously defined summary function)
	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"
	_polap_lib_timing-get_system_info >>"${_timing_txt}"
}

run-pca-ont-pt_genus_species() {
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"
	local _brg_knum="${3:-5}"
	local _brg_dimred="${4:-umap}"
	local _brg_genomesize="${5:-0}"

	local _brg_inum=0
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

	# local dimred="${_dimred["$_brg_target"]}"
	dimred="${_brg_dimred}"
	local readsel="${_readsel["$_brg_target"]}"
	local min_read="${_min_read["$_brg_target"]}"
	local platform="${_platform["$_brg_target"]}"
	local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

	rm -f "${_timing_txt}"
	rm -f "${_stdout_txt}"

	echo "="
	echo "------------------------------------------------------------------------------"
	echo "[INFO] Start an whole-genome HiFi assembly pipeline on ${_brg_outdir}-${_brg_sindex}"
	echo " - platform: ${platform}"
	echo " - dimension reduction: ${dimred}"
	echo " - read selection: ${readsel}"

	mkdir -p "${_brg_rundir}"

	# prepare data:
	# version 1 and 2: from local
	# version 3: from -1.00x data
	data-long_genus_species "${_brg_outdir}" "${_brg_sindex}"

	# Start memory logger
	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"

	_polap_lib_conda-ensure_conda_env polap || exit 1

	local resolved_fastq=$(<"${_brg_outdir}/${_brg_inum}/l.fastq.path.txt")
	echo "here FASTQ file ready: $resolved_fastq"

	# ${_polap_cmd} init \
	# 	-o "${_brg_rundir}"
	#
	# ${_polap_cmd} filter pca \
	# 	-l "${resolved_fastq}" \
	# 	-o "${_brg_rundir}"

	# o/kmer/k2.pca.tsv
	# o/kmer/k2.umap.tsv
	local K="${_brg_knum}"
	if [[ "$K" == "0" ]]; then
		K=2
	fi

	local organelle_fastq="${resolved_fastq}".k${K}.deselected.fastq
	local organelle_read_names="${_brg_rundir}"/kmer/k${K}.${dimred}.deselected.tsv
	local nuclear_fastq="${resolved_fastq}".k${K}.selected.fastq
	local nuclear_read_names="${_brg_rundir}"/kmer/k${K}.${dimred}.selected.tsv
	if [[ "${readsel}" == "ellipse" ]]; then
		organelle_read_names="${_brg_rundir}"/kmer/k${K}.${dimred}.ellipse.deselected.tsv
	fi

	# if [[ ! -s "${organelle_fastq}" ]]; then
	cut -f5 "${organelle_read_names}" |
		seqtk subseq "${resolved_fastq}" - \
			>"${organelle_fastq}"

	cut -f5 "${nuclear_read_names}" |
		seqtk subseq "${resolved_fastq}" - \
			>"${nuclear_fastq}"

	# conda deactivate
	#
	# # Summarize results after job (with previously defined summary function)
	# _polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"
	# _polap_lib_timing-get_system_info >>"${_timing_txt}"
	# return
	# # fi

	# genome size estimate
	# pl total-length-long -l DRR503528.fastq
	# pl find-genome-size-for-pacbio -l DRR503528.fastq
	# pl reduce-data -l DRR503528.fastq
	# pl flye1 --pacbio-hifi --blast --blast-mt 0.0001 --blast-pt 0.0001

	${_polap_cmd} total-length-long \
		-l "${organelle_fastq}" \
		-o "${_brg_rundir}"

	if [[ "${platform}" == "PACBIO_SMRT" ]]; then
		${_polap_cmd} find-genome-size-for-pacbio \
			-l "${organelle_fastq}" \
			-o "${_brg_rundir}"
	elif [[ "${platform}" == "ONT" ]]; then

		if [[ "${_brg_genomesize}" == "0" ]]; then
			# local long_len=$(<"${_brg_rundir}"/long_total_length.txt)
			# local reduced_len=$((long_len / 50))
			# echo $reduced_len >"${_brg_rundir}"/short_expected_genome_size.txt
			echo "ONT: genome size: auto flye 2.9.6 is used"
		else
			echo "ONT: genome size: auto flye 2.9.6 is used"
			# echo "${_brg_genomesize}" >"${_brg_rundir}"/short_expected_genome_size.txt

		fi
	fi

	# ${_polap_cmd} reduce-data \
	# 	-l "${organelle_fastq}" \
	# 	-o "${_brg_rundir}"

	# local genomesize=$(<"${_brg_rundir}"/short_expected_genome_size.txt)

	if [[ "${platform}" == "PACBIO_SMRT" ]]; then
		${_polap_cmd} flye1 --pacbio-hifi \
			-l "${_brg_rundir}"/lk.fq.gz \
			-g $genomesize \
			-o "${_brg_rundir}"
	elif [[ "${platform}" == "ONT" ]]; then
		# ${_polap_cmd} flye1 --nano-raw \
		# 	-g $genomesize \
		# 	-o "${_brg_rundir}"
		# ${_polap_cmd} blast-pt \
		# 	--blast-evalue 0 \
		# 	--nucleotide \
		# 	-l "${organelle_fastq}" \
		# 	-o "${_brg_rundir}"/b

		_polap_lib_conda-ensure_conda_env polap-flye || exit 1
		flye --nano-raw "${organelle_fastq}" \
			-t "${_brg_threads}" \
			--out-dir "${_brg_rundir}"/0
		flye --nano-raw "${resolved_fastq}" \
			-t "${_brg_threads}" \
			--out-dir "${_brg_rundir}"/a
		conda deactivate
	fi

	conda deactivate

	# Summarize results after job (with previously defined summary function)
	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"
	_polap_lib_timing-get_system_info >>"${_timing_txt}"
}

run-pca-ont-kmer_genus_species() {
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"

	# local _brg_knum="${3:-6}"
	local _brg_inum="0"
	# local _brg_knum="${5:-6}"

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

	local platform="${_platform["$_brg_target"]}"
	local short="${_short["$_brg_target"]}"
	local min_read="${_min_read["$_brg_target"]}"
	local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

	rm -f "${_timing_txt}"
	rm -f "${_stdout_txt}"

	# if [[ "${platform}" == "ONT" ]]; then
	# 	if [[ "${_brg_knum}" == "6" ]]; then
	# 		_brg_knum=5
	# 	fi
	# fi

	echo "------------------------------------------------------------------------------"
	echo "[INFO] Start an whole-genome HiFi assembly pipeline on ${_brg_outdir}-${_brg_sindex}"
	echo " - platform: ${platform}"
	echo " - K: 2 ~ 8 upto k as much memory as allowed"

	mkdir -p "${_brg_rundir}"

	# prepare data:
	# version 1 and 2: from local
	# version 3: from -1.00x data
	# data-downsample-long_genus_species "${_brg_outdir}" "${_brg_sindex}" 1
	data-long_genus_species "${_brg_outdir}" "${_brg_sindex}"

	# genome size estimate from a short-read dataset
	# downsample

	# Start memory logger
	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"

	_polap_lib_conda-ensure_conda_env polap || exit 1

	local resolved_fastq=$(<"${_brg_outdir}/${_brg_inum}/l.fastq.path.txt")
	echo "FASTQ file ready: $resolved_fastq"

	# local resolved_short_fastq=$(<"${_brg_outdir}/${_brg_inum}/s.fastq.path.txt")
	# echo "FASTQ file ready: $resolved_short_fastq"

	${_polap_cmd} init \
		-o "${_brg_rundir}"

	if [[ "${platform}" == "PACBIO_SMRT" ]]; then
		${_polap_cmd} filter pca \
			-k "${_brg_knum}" \
			-l "${resolved_fastq}" \
			-o "${_brg_rundir}"
	elif [[ "${platform}" == "ONT" ]]; then
		mkdir -p "${_brg_rundir}"/kmer/figures
		for k in {2..8}; do
			# for k in {7..9}; do
			${_polap_cmd} filter pca \
				-k "${k}" \
				-l "${resolved_fastq}" \
				-o "${_brg_rundir}"
			cp "${_brg_rundir}"/kmer/k*.density_heatmap.pdf \
				"${_brg_rundir}"/kmer/figures
		done
	fi

	# o/kmer/k2.pca.tsv
	# o/kmer/k2.umap.tsv

	conda deactivate

	# Summarize results after job (with previously defined summary function)
	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"
	_polap_lib_timing-get_system_info >>"${_timing_txt}"

	echo "see: ${_brg_rundir}/kmer"
}

run-pca-kmer_genus_species() {
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"
	local _brg_inum="${3:-0}"
	local _brg_knum="${4:-6}"

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

	local platform="${_platform["$_brg_target"]}"
	local short="${_short["$_brg_target"]}"
	local min_read="${_min_read["$_brg_target"]}"
	local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

	rm -f "${_timing_txt}"
	rm -f "${_stdout_txt}"

	if [[ "${platform}" == "ONT" ]]; then
		if [[ "${_brg_knum}" == "6" ]]; then
			_brg_knum=5
		fi
	fi

	echo "="
	echo "------------------------------------------------------------------------------"
	echo "[INFO] Start an whole-genome HiFi assembly pipeline on ${_brg_outdir}-${_brg_sindex}"
	echo " - platform: ${platform}"
	echo " - K: ${_brg_knum}"

	mkdir -p "${_brg_rundir}"

	# prepare data:
	# version 1 and 2: from local
	# version 3: from -1.00x data
	data-long_genus_species "${_brg_outdir}" "${_brg_sindex}"
	# if [[ "${platform}" == "ONT" ]]; then
	# 	data-short_genus_species "${_brg_outdir}" "${_brg_sindex}"
	# 	if [[ ! -s "${short}.fastq" ]]; then
	# 		data-downsample-short_genus_species "${_brg_outdir}" "${_brg_sindex}" 1 "off" ""
	# 		echo "${short}_1.fastq" + "${short}_2.fastq" = "${short}.fastq"
	# 		cat "${_brg_outdir}/tmp/${short}_1.fastq" "${_brg_outdir}/tmp/${short}_2.fastq" >"${short}.fastq"
	# 		# rm "${short}_1.fastq" "${short}_2.fastq"
	# 	fi
	# fi

	# Start memory logger
	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"

	_polap_lib_conda-ensure_conda_env polap || exit 1

	local resolved_fastq=$(<"${_brg_outdir}/${_brg_inum}/l.fastq.path.txt")
	echo "FASTQ file ready: $resolved_fastq"

	# local resolved_short_fastq=$(<"${_brg_outdir}/${_brg_inum}/s.fastq.path.txt")
	# echo "FASTQ file ready: $resolved_short_fastq"

	${_polap_cmd} init \
		-o "${_brg_rundir}"

	if [[ "${platform}" == "PACBIO_SMRT" ]]; then
		${_polap_cmd} filter pca \
			-k "${_brg_knum}" \
			-l "${resolved_fastq}" \
			-o "${_brg_rundir}"
	elif [[ "${platform}" == "ONT" ]]; then
		${_polap_cmd} filter pca \
			-k "${_brg_knum}" \
			-l "${resolved_fastq}" \
			-o "${_brg_rundir}"

		# ${_polap_cmd} filter short \
		# 	-k "${_brg_knum}" \
		# 	-l "${short}.fastq" \
		# 	-o "${_brg_rundir}"
	fi

	# o/kmer/k2.pca.tsv
	# o/kmer/k2.umap.tsv

	conda deactivate

	# Summarize results after job (with previously defined summary function)
	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"
	_polap_lib_timing-get_system_info >>"${_timing_txt}"

	echo "see: ${_brg_rundir}/kmer"
}

run-pca-kmer-short_genus_species() {
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"
	local _brg_inum="${3:-0}"
	local _brg_knum="${4:-6}"

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

	local platform="${_platform["$_brg_target"]}"
	local short="${_short["$_brg_target"]}"
	local min_read="${_min_read["$_brg_target"]}"
	local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

	rm -f "${_timing_txt}"
	rm -f "${_stdout_txt}"

	if [[ "${platform}" == "ONT" ]]; then
		if [[ "${_brg_knum}" == "6" ]]; then
			_brg_knum=5
		fi
	fi

	echo "="
	echo "------------------------------------------------------------------------------"
	echo "[INFO] Start an whole-genome HiFi assembly pipeline on ${_brg_outdir}-${_brg_sindex}"
	echo " - platform: ${platform}"
	echo " - K: ${_brg_knum}"

	mkdir -p "${_brg_rundir}"

	# prepare data:
	# version 1 and 2: from local
	# version 3: from -1.00x data
	# data-long_genus_species "${_brg_outdir}" "${_brg_sindex}"
	if [[ "${platform}" == "ONT" ]]; then
		data-short_genus_species "${_brg_outdir}" "${_brg_sindex}"
		if [[ ! -s "${short}.fastq" ]]; then
			data-downsample-short_genus_species "${_brg_outdir}" "${_brg_sindex}" 1 "off" ""
			echo "${short}_1.fastq" + "${short}_2.fastq" = "${short}.fastq"
			cat "${_brg_outdir}/tmp/${short}_1.fastq" "${_brg_outdir}/tmp/${short}_2.fastq" >"${short}.fastq"
			rm "${short}_1.fastq" "${short}_2.fastq"
		fi
	fi

	# Start memory logger
	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"

	_polap_lib_conda-ensure_conda_env polap || exit 1

	local resolved_fastq=$(<"${_brg_outdir}/${_brg_inum}/l.fastq.path.txt")
	echo "FASTQ file ready: $resolved_fastq"

	# local resolved_short_fastq=$(<"${_brg_outdir}/${_brg_inum}/s.fastq.path.txt")
	# echo "FASTQ file ready: $resolved_short_fastq"

	${_polap_cmd} init \
		-o "${_brg_rundir}"

	if [[ "${platform}" == "PACBIO_SMRT" ]]; then
		${_polap_cmd} filter pca \
			-k "${_brg_knum}" \
			-l "${resolved_fastq}" \
			-o "${_brg_rundir}"
	elif [[ "${platform}" == "ONT" ]]; then
		# ${_polap_cmd} filter pca \
		# 	-k "${_brg_knum}" \
		# 	-l "${resolved_fastq}" \
		# 	-o "${_brg_rundir}"

		${_polap_cmd} filter short \
			-k "${_brg_knum}" \
			-l "${short}.fastq" \
			-o "${_brg_rundir}"
	fi

	# o/kmer/k2.pca.tsv
	# o/kmer/k2.umap.tsv

	conda deactivate

	# Summarize results after job (with previously defined summary function)
	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"
	_polap_lib_timing-get_system_info >>"${_timing_txt}"

	echo "see: ${_brg_rundir}/kmer"
}

run-pca-count-filter_genus_species() {
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"
	local _brg_inum="${3:-0}"
	local _brg_ranger="${4:-0}"

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

	# local range1="${_range1["$_brg_target"]//:/,}"
	# local range2="${_range2["$_brg_target"]//:/,}"
	local first_item=0
	local second_item=0
	local third_item=0
	if [[ "${_brg_ranger}" != "0" ]]; then
		local delimiter=","
		local range1="${_range1["$_brg_target"]//:/,}"
		local -n range="range${_brg_ranger}"
		local first_item=$(echo "$range" | cut -d "$delimiter" -f 1)
		local second_item=$(echo "$range" | cut -d "$delimiter" -f 2)
		local third_item=$(echo "$range" | cut -d "$delimiter" -f 3)
	fi

	local min_read="${_min_read["$_brg_target"]}"
	local platform="${_platform["$_brg_target"]}"
	local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

	rm -f "${_timing_txt}"
	rm -f "${_stdout_txt}"

	echo "="
	echo "------------------------------------------------------------------------------"
	echo "[INFO] Start an whole-genome HiFi assembly pipeline on ${_brg_outdir}-${_brg_sindex}"

	mkdir -p "${_brg_rundir}"

	# prepare data:
	# version 1 and 2: from local
	# version 3: from -1.00x data
	data-long_genus_species "${_brg_outdir}" "${_brg_sindex}"

	# Start memory logger
	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"

	_polap_lib_conda-ensure_conda_env polap || exit 1

	local resolved_fastq=$(<"${_brg_outdir}/${_brg_inum}/l.fastq.path.txt")
	echo "here FASTQ file ready: $resolved_fastq"

	echo "platform: ${platform}"

	${_polap_cmd} init \
		-o "${_brg_rundir}"

	if [[ "${platform}" == "PACBIO_SMRT" ]]; then
		if [[ "${first_item}" == "${second_item}" ]]; then
			${_polap_cmd} count filter \
				--pacbio-hifi \
				--x-min 0 --x-max 1500 \
				-l "${resolved_fastq}" \
				-o "${_brg_rundir}"
		else
			${_polap_cmd} count filter \
				--pacbio-hifi \
				--x-min "${first_item}" --x-max "${second_item}" -c "${third_item}" \
				-l "${resolved_fastq}" \
				-o "${_brg_rundir}"
		fi
	elif [[ "${platform}" == "ONT" ]]; then
		if [[ "${first_item}" == "${second_item}" ]]; then
			${_polap_cmd} count filter \
				--nano-raw \
				-l "${resolved_fastq}" \
				-o "${_brg_rundir}"
		else
			${_polap_cmd} count filter \
				--nano-raw \
				--x-min "${first_item}" --x-max "${second_item}" -c "${third_item}" \
				-l "${resolved_fastq}" \
				-o "${_brg_rundir}"
		fi

	fi

	# o/kmer/k2.pca.tsv
	# o/kmer/k2.umap.tsv

	conda deactivate

	# Summarize results after job (with previously defined summary function)
	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"
	_polap_lib_timing-get_system_info >>"${_timing_txt}"

}

run-pca-count_genus_species() {
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"
	local _brg_inum="${3:-0}"
	local _brg_ranger="${4:-0}"

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

	# local range1="${_range1["$_brg_target"]//:/,}"
	# local range2="${_range2["$_brg_target"]//:/,}"
	local first_item=0
	local second_item=0
	if [[ "${_brg_ranger}" != "0" ]]; then
		local delimiter=","
		local range1="${_range1["$_brg_target"]//:/,}"
		local -n range="range${_brg_ranger}"
		local first_item=$(echo "$range" | cut -d "$delimiter" -f 1)
		local second_item=$(echo "$range" | cut -d "$delimiter" -f 2)
	fi

	local min_read="${_min_read["$_brg_target"]}"
	local platform="${_platform["$_brg_target"]}"
	local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

	rm -f "${_timing_txt}"
	rm -f "${_stdout_txt}"

	echo "="
	echo "------------------------------------------------------------------------------"
	echo "[INFO] Start an whole-genome HiFi assembly pipeline on ${_brg_outdir}-${_brg_sindex}"

	mkdir -p "${_brg_rundir}"

	# prepare data:
	# version 1 and 2: from local
	# version 3: from -1.00x data
	data-long_genus_species "${_brg_outdir}" "${_brg_sindex}"

	# Start memory logger
	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"

	_polap_lib_conda-ensure_conda_env polap || exit 1

	local resolved_fastq=$(<"${_brg_outdir}/${_brg_inum}/l.fastq.path.txt")
	echo "here FASTQ file ready: $resolved_fastq"

	echo "platform: ${platform}"

	${_polap_cmd} init \
		-o "${_brg_rundir}"

	if [[ "${platform}" == "PACBIO_SMRT" ]]; then
		if [[ "${first_item}" == "${second_item}" ]]; then
			${_polap_cmd} count view \
				--pacbio-hifi \
				--x-min 0 --x-max 1500 \
				-l "${resolved_fastq}" \
				-o "${_brg_rundir}"
		else
			${_polap_cmd} count view \
				--pacbio-hifi \
				--x-min "${first_item}" --x-max "${second_item}" \
				-l "${resolved_fastq}" \
				-o "${_brg_rundir}"
		fi
	elif [[ "${platform}" == "ONT" ]]; then
		if [[ "${first_item}" == "${second_item}" ]]; then
			${_polap_cmd} count view \
				--nano-raw \
				--x-min 0 --x-max 1500 \
				-l "${resolved_fastq}" \
				-o "${_brg_rundir}"
		else
			${_polap_cmd} count view \
				--nano-raw \
				--x-min "${first_item}" --x-max "${second_item}" \
				-l "${resolved_fastq}" \
				-o "${_brg_rundir}"
		fi

	fi

	# o/kmer/k2.pca.tsv
	# o/kmer/k2.umap.tsv

	conda deactivate

	# Summarize results after job (with previously defined summary function)
	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"
	_polap_lib_timing-get_system_info >>"${_timing_txt}"

}

run-download-sra-clr_genus_species() {
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"

	local _brg_adir _brg_title _brg_target _brg_rundir _brg_outdir_i
	local _timing_txt _stdout_txt _memlog_file _summary_file
	local long_sra
	source "${_POLAPLIB_DIR}/polap-variables-data.sh"

	echo "------------------------------------------------------------------------------"
	echo "[INFO] Download organelle-genome CLR assembly pipeline on ${_brg_outdir}-${_brg_sindex}"

	mkdir -p "${_brg_rundir}"
	mkdir -p "${_brg_outdir}"

	download-sra_genus_species "${long_sra}"

	rm -rf "${long_sra}"

	_polap_lib_conda-ensure_conda_env polap || exit 1

	# ${_polap_cmd} init \
	# 	-o "${_brg_rundir}"

	if [[ -s "${_brg_rundir}/long_total_length.txt" ]]; then
		echo "found: ${_brg_rundir}/long_total_length.txt"
	else
		${_polap_cmd} total-length-long -l "${long_sra}".fastq \
			-o "${_brg_rundir}"
	fi

	if [[ -s "${_brg_rundir}/short_expected_genome_size.txt" ]]; then
		echo "found: ${_brg_rundir}/short_expected_genome_size.txt"
	else
		${_polap_cmd} find-genome-size-for-pbclr -l "${long_sra}".fastq \
			-o "${_brg_rundir}"
	fi

	if [[ -s "${_brg_rundir}"/short_expected_genome_size.txt ]]; then
		local _genome_size=$(<"${_brg_rundir}/short_expected_genome_size.txt")
	else
		echo "ERROR: no such file: ${_brg_outdir}/short_expected_genome_size.txt"
		return
	fi
	cp -p "${_brg_rundir}/short_expected_genome_size.txt" "${_brg_outdir}"
	cp -p "${_brg_rundir}/long_total_length.txt" "${_brg_outdir}"

	local _sample_coverage=10
	seqkit seq -m 2000 "${long_sra}".fastq -o "${long_sra}".2k.fastq
	mv "${long_sra}".2k.fastq "${long_sra}".fastq

	# downsample to the target coverage
	# if [[ ! -s "${long_sra}-${_sample_coverage}x.fq" ]]; then
	# 	sample-fastq_genus_species "${long_sra}".fastq "${_brg_rundir}" "${_sample_coverage}"
	# fi

	conda deactivate

	# Transfer data back to the storage for later use
	# mv "${long_sra}-${_sample_coverage}"x.fq "${long_sra}.fastq"
	# tar zcf "${long_sra}-${_sample_coverage}x.fastq.tar.gz" "${long_sra}.fastq"
	tar zcf "${long_sra}.2k.fastq.tar.gz" "${long_sra}.fastq"
	scp -qp "${long_sra}.2k.fastq.tar.gz" "${_local_host}:/media/h1/sra/"

	rsync -azuq --max-size=5M \
		"${_brg_outdir}"/ \
		"${_local_host}:${PWD}/${_brg_outdir}"/

	# Clean-up
	rm -f "${long_sra}.fastq"
	rm -f "${long_sra}.2k.fastq.tar.gz"

	rm -rf "${_brg_rundir}"
	rm -rf "${_brg_outdir}"
}

run-download-sra-hifi_genus_species() {
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

	echo "------------------------------------------------------------------------------"
	echo "[INFO] Start an organelle-genome HiFi assembly pipeline on ${_brg_outdir}-${_brg_sindex}"

	mkdir -p "${_brg_rundir}"

	for i in long_total_length.txt short_expected_genome_size.txt; do
		scp "${_local_host}":$PWD/${_brg_outdir}/t5/0/$i "${_brg_rundir}"
	done

	download-sra_genus_species "${long_sra}"
	rm -rf "${long_sra}"

	_polap_lib_conda-ensure_conda_env polap || exit 1

	# ${_polap_cmd} init \
	# 	-o "${_brg_rundir}"

	if [[ -s "${_brg_rundir}/long_total_length.txt" ]]; then
		echo "found: ${_brg_rundir}/long_total_length.txt"
	else
		${_polap_cmd} total-length-long -l "${long_sra}".fastq \
			-o "${_brg_rundir}"
	fi

	if [[ -s "${_brg_rundir}/short_expected_genome_size.txt" ]]; then
		echo "found: ${_brg_rundir}/short_expected_genome_size.txt"
	else
		${_polap_cmd} find-genome-size-for-pacbio -l "${long_sra}".fastq \
			-o "${_brg_rundir}"
	fi

	if [[ -s "${_brg_rundir}"/short_expected_genome_size.txt ]]; then
		local _genome_size=$(<"${_brg_rundir}/short_expected_genome_size.txt")
	else
		echo "ERROR: no such file: ${_brg_outdir}/short_expected_genome_size.txt"
		return
	fi

	local _sample_coverage=10

	# downsample to the target coverage
	if [[ ! -s "${long_sra}-${_sample_coverage}x.fq" ]]; then
		sample-fastq_genus_species "${long_sra}".fastq "${_brg_rundir}" "${_sample_coverage}"
	fi

	conda deactivate

	scp -qp "${long_sra}-${_sample_coverage}"x.fq "${_local_host}:$PWD"/

	rm -f "${long_sra}.fastq"
	rm -rf "${_brg_rundir}"

	rm -rf "${_brg_outdir}"
	rm -f "${long_sra}-${_sample_coverage}"x.fq
}

run-hifi-oga_genus_species() {
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"
	local _brg_inum="${3:-0}"
	local _brg_ranger="${4:-0}"
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

	local range0="${_range0["$_brg_target"]//:/,}"
	local range1="${_range1["$_brg_target"]//:/,}"
	local range2="${_range2["$_brg_target"]//:/,}"
	local -n range="range${_brg_ranger}"

	local min_read="${_min_read["$_brg_target"]}"
	local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

	rm -f "${_timing_txt}"
	rm -f "${_stdout_txt}"

	echo "="
	echo "------------------------------------------------------------------------------"
	echo "[INFO] Start an organelle-genome HiFi assembly pipeline on ${_brg_outdir}-${_brg_sindex}"

	local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

	mkdir -p "${_brg_rundir}"

	_create_range_float "${range}" \
		"${_brg_rundir}/range.txt"

	download-sra_genus_species "${long_sra}"
	rm -rf "${long_sra}"

	# Start memory logger
	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"

	_polap_lib_conda-ensure_conda_env polap || exit 1

	${_polap_cmd} init \
		-o "${_brg_rundir}"

	if [[ -s "${_brg_rundir}/${long_sra}".fastq.txt ]]; then
		echo "found: ${_brg_rundir}/${long_sra}.fastq.txt"
	else
		${_polap_cmd} fastq-inspect "${long_sra}".fastq \
			"${_brg_rundir}/${long_sra}".fastq.txt \
			-o "${_brg_rundir}"
	fi

	local seq_type=$(<"${_brg_rundir}/${long_sra}".fastq.txt)
	if [[ "${seq_type}" == "HiFi" ]]; then
		echo "[INFO] SRA ${long_sra}: HiFi data"
		echo "${_brg_outdir}:SRA:${long_sra}" >>report.hifi.txt
	else
		echo "ERROR: SRA ${long_sra}: no HiFi data: ${seq_type}"
		echo "${_brg_outdir}:SRA:${long_sra}" >>report.no-hifi.txt
		ssh ${_local_host} "touch $PWD/no-hifi-${_brg_outdir}"
		rm -f "${long_sra}.fastq"
		_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"
		conda deactivate
		return
	fi

	if [[ -s "${_brg_rundir}/long_total_length.txt" ]]; then
		echo "found: ${_brg_rundir}/long_total_length.txt"
	else
		${_polap_cmd} total-length-long -l "${long_sra}".fastq \
			-o "${_brg_rundir}"
	fi

	if [[ -s "${_brg_rundir}/short_expected_genome_size.txt" ]]; then
		echo "found: ${_brg_rundir}/short_expected_genome_size.txt"
	else
		${_polap_cmd} find-genome-size-for-pacbio -l "${long_sra}".fastq \
			-o "${_brg_rundir}"
	fi

	# if [[ ! -s "${_brg_rundir}/lk.fq.gz" ]]; then
	# 	${_polap_cmd} reduce-data -l "${long_sra}".fastq \
	# 		-o "${_brg_rundir}"
	# fi

	local is_complete=0

	read -a restored_array <"${_brg_rundir}/range.txt"
	local array_length=${#restored_array[@]}

	# 10 Mb at least
	# 100 x 10 Mb / Genome size => sampling percent minimum
	if [[ -s "${_brg_rundir}"/short_expected_genome_size.txt ]]; then
		local _genome_size=$(<"${_brg_rundir}/short_expected_genome_size.txt")
	else
		echo "ERROR: no such file: ${_brg_outdir}/short_expected_genome_size.txt"
		return
	fi
	if ((_genome_size < 500000000)); then
		local _coverage=$(echo "scale=2; 500000000 / $_genome_size" | bc)
		restored_array[0]="$_coverage"
	else
		restored_array[0]="1.00"
	fi

	# Iterate over the array using an index
	for ((i = 0; i < array_length; i++)); do
		local _sample_coverage="${restored_array[i]}"

		# downsample to the target coverage
		if [[ ! -s "${long_sra}-${_sample_coverage}x.fq" ]]; then
			sample-fastq_genus_species "${long_sra}".fastq "${_brg_rundir}" "${_sample_coverage}"
		fi

		${_polap_cmd} blast-pt \
			-l "${long_sra}-${_sample_coverage}"x.fq \
			-o "${_brg_rundir}/${_sample_coverage}"

		if [[ ! -s "${_brg_rundir}/${_sample_coverage}/blast/pt/30-contigger/graph_final.gfa" ]]; then
			echo "No plastid assembly from ${_sample_coverage}x data"
			ssh ${_local_host} "touch $PWD/no-plastid-assembly-${_brg_outdir}-${_sample_coverage}x"
			continue
		fi

		${_polap_cmd} annotate \
			-o "${_brg_rundir}/${_sample_coverage}/blast" \
			-i pt

		${_polap_cmd} seeds \
			--plastid \
			-o "${_brg_rundir}/${_sample_coverage}/blast" \
			-i pt -j 1

		# dummy file to avoid errors in assemble2 polap command
		# ln -sf "../../lk.fq.gz" \
		# 	"${_brg_rundir}/${_sample_coverage}/blast/lk.fq.gz"
		if [[ ! -s "${_brg_rundir}/${_sample_coverage}/blast/lk.fq.gz" ]]; then
			echo 1 >"${_brg_rundir}/${_sample_coverage}/blast/lk.fq.gz"
		fi

		${_polap_cmd} assemble2 \
			--pacbio-hifi \
			--coverage-oga 300 \
			-l "${long_sra}".fastq \
			-o "${_brg_rundir}/${_sample_coverage}/blast" \
			-i pt -j 1

		local gfa="${_brg_rundir}/${_sample_coverage}/blast/1/assembly_graph.gfa"
		if [[ -s "${gfa}" ]]; then
			is_complete=1
			ln -sf "${_sample_coverage}/blast/1/assembly_graph.gfa" \
				"${_brg_rundir}/assembly_graph.gfa"
		fi

		# now, remove reads mapped on the ptDNA from the input data
		# create ptdna4.fa with 4 possible ptDNA sequences from the ptDNA gfa
		# polap command subtract reads mapped on ptdna4.fa from the input long-read data
		# sample 1x of the subtracted data, wga, and oga or blast-mt
		#
		# local assembly_fasta="${_brg_rundir}/${_sample_coverage}/blast/1/assembly.fasta"
		# ${_polap_cmd} subtract \
		# 	-l "${long_sra}".fastq \
		# 	-p "${assembly_fasta}" \
		# 	-o "${_brg_rundir}/${_sample_coverage}"
		#
		# data without ptDNA
		# local subtracted_fq="${_brg_rundir}/${_sample_coverage}"/subtract/ptgaul-reads/0.fq.gz
		#
		# if [[ ! -s "${long_sra}-${_sample_coverage}x.fq" ]]; then
		# 	sample-fastq_genus_species "${subtracted_fq}" "${_brg_rundir}" "${_sample_coverage}"
		# fi
		#
		# echo "use the sampled subtracted_fq to assemble mtDNA"

		echo "${_brg_rundir}"/blast/pt
		# pacbio hifi assembly - clr and raw are not working yet
	done

	# Summarize results after job (with previously defined summary function)
	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"

	conda deactivate

	_polap_lib_timing-get_system_info >>"${_timing_txt}"

	if [[ "${is_complete}" == "1" ]]; then
		mkdir -p "${_brg_outdir_i}"
		rsync -azuq --max-size=5M \
			"${_brg_rundir}"/ \
			"${_brg_outdir_i}"/

		ssh ${_local_host} "mkdir -p $PWD/${_brg_outdir}"

		rsync -azuq --max-size=5M \
			"${_brg_outdir}"/ \
			"${_local_host}:${PWD}/${_brg_outdir}"/

		scp -qp "${long_sra}-${_sample_coverage}"x.fq "${_local_host}:$PWD"/
		scp -qp "${long_sra}-${_sample_coverage}"x_pt.fq.gz "${_local_host}:$PWD"/
		scp -qp "${long_sra}.fastq.tsv" "${_local_host}:$PWD"/

		rm -f "${long_sra}.fastq"
		rm -rf "${_brg_rundir}"
	else
		echo "[INFO] no assembly for SRA ${long_sra}"
		if [[ ! -s "${_brg_rundir}/${_sample_coverage}/blast/pt/30-contigger/graph_final.gfa" ]]; then
			echo "No plastid assembly from ${_sample_coverage}x data"
			ssh ${_local_host} "touch $PWD/no-plastid-assembly-${_brg_outdir}-${_sample_coverage}x"
		fi
	fi

	rm -rf "${_brg_outdir}"
	rm -f "${long_sra}-${_sample_coverage}"x_pt.fq.gz
	rm -f "${long_sra}-${_sample_coverage}"x.fq
}

sample-fastq_genus_species() {
	local _brg_infile="$1"
	local _brg_outdir="${2:-o}"
	local _brg_coverage="${3:-10}"

	local _brg_infilename="${_brg_infile%%.*}"

	# source "$(conda info --base)/etc/profile.d/conda.sh"
	# if [[ "$CONDA_DEFAULT_ENV" != "polap" ]]; then
	# 	echo "You're not in the polap environment. Chaniging 'polap'..."
	# 	conda activate polap
	# fi
	_polap_lib_conda-ensure_conda_env polap || exit 1

	if [[ "$CONDA_DEFAULT_ENV" == "polap" ]]; then

		if [[ -s "${_brg_outdir}"/short_expected_genome_size.txt ]]; then
			echo "found: ${_brg_outdir}/short_expected_genome_size.txt"
			local _genome_size=$(<"${_brg_outdir}/short_expected_genome_size.txt")
		else
			echo "ERROR: no such file: ${_brg_outdir}/short_expected_genome_size.txt"
			return
		fi

		echo "Genome size estimate: ${_genome_size}"

		# Step 3.
		echo "downsampling ..."
		echo "input: ${_brg_infile}"
		echo "output: ${_brg_infilename}-${_brg_coverage}x.fq"
		${_polap_cmd} fastq subsample \
			"${_brg_infile}" \
			"${_brg_infilename}-${_brg_coverage}x.fq" \
			-c "${_brg_coverage}" \
			-o "${_brg_outdir}" \
			--genomesize "${_genome_size}" \
			>"${_brg_outdir}/${_brg_infilename}-${_brg_coverage}x.txt"
		echo "log: ${_brg_outdir}/${_brg_infilename}-${_brg_coverage}x.txt"
		echo "log: ${_brg_infilename}-${_brg_coverage}x.fq"
	else
		echo "ERROR: no conda environment of polap"
	fi
}

get-aflye_genus_species() {
	local _brg_outdir="${1}"
	local _brg_host="${2}"

	if [[ "${_brg_outdir}" == "-h" || "${_brg_outdir}" == "--help" ]]; then
		echo "$help_message_get_aflye"
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
	echo Getting aflye1 results of ${_brg_outdir} from ${_brg_host} ...
	rsync -azuq --max-size=5M \
		"${_brg_host}:${PWD}/${_brg_outdir}"/ \
		"${_brg_outdir}"/
	rsync -azq \
		"${_brg_host}:${PWD}/${_brg_outdir}/0/30-contigger/graph_final.gfa" \
		"${_brg_outdir}/0/30-contigger"
}

polap-cmd-xxx_genus_species() {
	_polap_lib_conda-ensure_conda_env polap || exit 1

	"${_polap_cmd}" xxx

	conda deactivate
}

polap-cmd-find-genome-size_genus_species() {
	local _brg_outdir="${1}"
	local _brg_long="${2}"

	_polap_lib_conda-ensure_conda_env polap || exit 1

	"${_polap_cmd}" find-genome-size-for-pacbio -o "${_brg_outdir}" -l "${_brg_long}"

	conda deactivate
}

polap-cmd-total-length-long_genus_species() {
	local _brg_outdir="${1}"
	local _brg_long="${2}"

	_polap_lib_conda-ensure_conda_env polap || exit 1

	"${_polap_cmd}" total-length-long -o "${_brg_outdir}" -l "${_brg_long}"

	conda deactivate
}

polap-cmd_genus_species() {
	local args=("$@")

	local item="${args[0]}"
	local rest=("${args[@]:1}")
	polap-cmd-${item}_genus_species ${rest[@]}
}

polap-cmd-clean-menus_genus_species() {
	_polap_lib_conda-ensure_conda_env polap || exit 1

	"${_polap_cmd}" clean-menus

	conda deactivate
}

update-dev_genus_species() {
	local _brg_outdir="${1:-$HOME/all/polap}"

	cd "${_brg_outdir}"
	rsync -aq thorne:$PWD/github/ github/
	cd -
}

polap-cmd-init_genus_species() {
	local _brg_outdir="${1}"

	_polap_lib_conda-ensure_conda_env polap || exit 1

	"${_polap_cmd}" init -o "${_brg_outdir}"

	conda deactivate
}

download-sra_genus_species() {
	local args=("$@")

	if has_help "${args[@]}"; then
		echo "$help_message_download_sra"
		return
	fi

	_polap_lib_conda-ensure_conda_env polap || exit 1
	for item in "$@"; do
		if [[ -s "${item}.fastq" ]]; then
			_log_echo "found: ${item}.fastq"
		elif [[ -s "${item}_1.fastq" ]]; then
			_log_echo "found: ${item}_1.fastq"
		else
			_log_echo "  downloading long-read SRA ID: ${item} ... be patient!"
			bash "${_polap_script_bin_dir}"/polap-ncbitools fetch sra "$item"
		fi
	done
	conda deactivate
}

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

config-info_genus_species() {
	local _brg_csv="${1:-polap-data-read.csv}"

	echo "CSV: ${_brg_csv}"
	echo "columns: species,long,short,version,platform,type,folder,bioproject,species2,base,nseq,size,libray,type,model,taxon,coverage,down,dummy,genomesize,host,inref,min_read,plastid,range1,range2,range3,dimred,readsel,bench_ptgaul,bench_pmat,bench_tippo,bench_oatk,ptgaul_ref,mt_ref,random,random_seed,status"
}

config-add_genus_species() {
	local _brg_species="$1"
	local _brg_sindex="${2:-0}"
	local _brg_long="${3}"
	local _brg_short="${4}"
	local _brg_platform="${5}"
	local _brg_csv="${6:-polap-data-read.csv}"

	local _brg_target="${_brg_species}-${_brg_sindex}"
	# source "${_POLAPLIB_DIR}/polap-variables-data.sh"

	local confirm="no"
	if [[ -v _long["${_brg_target}"] ]]; then
		echo "Warning: ${_brg_target} is in the CSV."
		if [[ "${opt_y_flag}" == false ]]; then
			read -p "Do you want to add one more of ${_brg_target}? (y/N): " confirm
		else
			confirm="yes"
		fi
	else
		confirm="yes"
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

	if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
		# echo cp -p "${third_arg}" "${_POLAPLIB_DIR}/${_polap_data_csv}"
		local species_name="${_brg_species//_/ }"

		echo "species,long,short,version,platform,type,folder,bioproject,species2,base,nseq,size,libray,type,model,taxon,coverage,down,dummy,genomesize,host,inref,min_read,plastid,range1,range2,range3,dimred,readsel,bench_ptgaul,bench_pmat,bench_tippo,bench_oatk,ptgaul_ref,mt_ref,random,random_seed,status" >"${_brg_csv}"
		printf "%s-%d,%s,%s,3,%s,NA,%s,PRJEB89667,%s,46922076059,1297,11436,WGS,GENOMIC,Sequel II,%s,1,10,dummy,0,lab01,family,3000,TRUE,1:1:1,0.1:1:10,3000:6000:2,pca,manual,160000,0.1,onthq,30,NA,%s,367,367,done\n" \
			"${_brg_species}" \
			"${_brg_sindex}" \
			"${_brg_long}" \
			"${_brg_short}" \
			"${_brg_platform}" \
			"${_brg_species}" \
			"${_brg_species}" \
			"${species_name}" \
			"${species_name}" \
			>>"${_brg_csv}"
	fi
}

config-add-field_genus_species() {
	local first_arg="$1"
	local second_arg="${2:-0}"
	local third_arg="${3:-out.csv}"
	print_species_field_summary --add-field="${first_arg}=${second_arg}" \
		--values --out="${third_arg}"
	echo cp -p "${third_arg}" "${_POLAPLIB_DIR}/${_polap_data_csv}"
}

config-update-field_genus_species() {
	local first_arg="$1"
	local second_arg="${2:-0}"
	local third_arg="${3:-0}"
	local csv_arg="${4:-out.csv}"
	print_species_field_summary --set-field="${first_arg}:${second_arg}:${third_arg}" \
		--values --out="${csv_arg}"
	echo cp -p "${csv_arg}" "${_POLAPLIB_DIR}/${_polap_data_csv}"
}

config_genus_species() {
	local first_arg="$1"
	local remaining_args=("${@:2}")

	# Remove trailing slash from the first element
	# remaining_args[0]="${remaining_args[0]%/}"

	config-${first_arg}_genus_species "${remaining_args[@]:-}"
}

config-view_genus_species() {
	# _log_echo0 "outdir: $_brg_outdir"
	# _log_echo0 "sindex: $_brg_sindex"

	local bolap_cmd="${FUNCNAME##*_}"

	help_message=$(
		cat <<EOF
Name:
  bolap - config view

Synopsis:
  bolap config view --fields long,short --match Anthoceros 

Description:
  bolap

Examples:
  View config:
    bolap config view --fields long,short --match Anthoceros 

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	local _fields="all"
	local _match="any"

	parse_commandline() {
		set -- "${_brg_unknown_opts[@]}"

		# source "${_POLAPLIB_DIR}/polap-cmd-version.sh" # '.' means 'source'
		while test $# -gt 0; do
			_key="$1"
			case "$_key" in
			--fields)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_fields="$2"
					shift || true
				fi
				;;
			--match)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_match="$2"
					shift || true
				fi
				;;
			-*)
				_log_echo1 "[INFO] no such options: $1"
				;;
			*)
				break
				;;
			esac
			shift || true
		done
	}

	declare -n ref="help_message"
	if [[ "${_brg_help}" == "on" ]]; then
		local manfile=$(_bolap_lib_man-convert_help_message "$ref" "${bolap_cmd}")
		man "$manfile"
		rm -f "$manfile"
		return
	fi

	parse_commandline

	echo "reading data config ..." >&2

	# print_species_field_summary --add-field=fruit=banana --values
	if [[ "${_fields}" == "all" ]]; then
		if [[ "${_match}" == "any" ]]; then
			print_species_field_summary --values
		else
			print_species_field_summary --values --match="${_match}"
		fi
	else
		if [[ "${_match}" == "any" ]]; then
			print_species_field_summary --values --fields="${_fields}"
		else
			print_species_field_summary --values --match="${_match}" --fields="${_fields}"
		fi
	fi

}

#-------------------------------------------------------------------------------
# dataset
dataset-import_genus_species() {
	local bolap_cmd="${FUNCNAME%%_*}"

	help_message=$(
		cat <<EOF
Name:
  bolap - dataset import

Synopsis:
  bolap $bolap_cmd --infile CSV

Description:
  Prepare a CSV file with columns including species, long, and short. Additional fields are optional.

Examples:
  Import in.csv as a json dataset:
    bolap $bolap_cmd --infile in.csv

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	local infile=""
	parse_commandline() {
		set -- "${_brg_args[@]}"

		# source "${_POLAPLIB_DIR}/polap-cmd-version.sh" # '.' means 'source'
		while test $# -gt 0; do
			_key="$1"
			case "$_key" in
			--infile)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					infile="$2"
					shift || true
				fi
				;;
			-*)
				_log_echo0 "[ERROR] no such options: $1"
				;;
			esac
			shift || true
		done
	}

	declare -n ref="help_message"
	if [[ "${_brg_help}" == "on" ]]; then
		local manfile=$(_bolap_lib_man-convert_help_message "$ref" "${bolap_cmd}")
		man "$manfile"
		rm -f "$manfile"
		return
	fi

	parse_commandline

	_polap_lib_conda-ensure_conda_env polap || exit 1
	if [[ -s "$infile" ]]; then
		# Basic import (key defaults to 'species', stores everything as strings)
		dataset-import-csv "$infile"
	else
		_log_echo0 "No such file: $infile"
	fi
	conda deactivate

}

dataset-view_genus_species() {
	local bolap_cmd="${FUNCNAME%%_*}"

	help_message=$(
		cat <<EOF
Name:
  bolap - dataset view

Synopsis:
  bolap $bolap_cmd

Description:
  bolap

Examples:
  View the json dataset:
    bolap $bolap_cmd --fields long,short

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	local fields=""
	parse_commandline() {
		set -- "${_brg_args[@]}"

		# source "${_POLAPLIB_DIR}/polap-cmd-version.sh" # '.' means 'source'
		while test $# -gt 0; do
			_key="$1"
			case "$_key" in
			--fields)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					fields="$2"
					shift || true
				fi
				;;
			-*)
				_log_echo0 "[ERROR] no such options: $1"
				;;
			esac
			shift || true
		done
	}

	declare -n ref="help_message"
	if [[ "${_brg_help}" == "on" ]]; then
		local manfile=$(_bolap_lib_man-convert_help_message "$ref" "${bolap_cmd}")
		man "$manfile"
		rm -f "$manfile"
		return
	fi

	parse_commandline

	_polap_lib_conda-ensure_conda_env polap || exit 1
	if [[ -z "${fields:-}" ]]; then
		dataset-view-table
	else
		dataset-view-table --fields="$fields"
	fi
	# echo dataset-view-table --fields="$fields" >&2
	conda deactivate

}

demo_genus_species() {
	local first_arg="$1"
	local remaining_args=("${@:2}")

	local title="${FUNCNAME[0]%%_*}"

	# Remove trailing slash from the first element
	if [[ -n "${remaining_args[0]}" ]]; then
		remaining_args[0]="${remaining_args[0]%/}"
	fi

	if declare -F ${title}-${first_arg}_genus_species >/dev/null; then
		${title}-${first_arg}_genus_species "${remaining_args[@]}"
	else
		echo "No such method: ${title}-${first_arg}"
		local manfile=$(_bolap_lib_man-convert_help_message "$help_message_demo" "demo")
		man "$manfile"
		rm -f "$manfile"
		return
	fi
}

demo-read_genus_species() {
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"
	source "${_POLAPLIB_DIR}/polap-variables-data.sh"
	local _brg_inum_0="0"

	# local _brg_outdir="Salix_purpurea"
	# local _brg_sindex="0"
	# local _brg_inum_0="0"

	# create a config
	if [[ ! -s "a.csv" ]]; then
		_log_echo0 "bolap config add ${_brg_outdir} ${_brg_sindex} SRR21824870 NA a.csv -> a.csv"
		config-add_genus_species "${_brg_outdir}" "${_brg_sindex}" SRR21824870 NA a.csv
		echo "rerun it!"
		return
	fi
	if [[ ! -s "${_brg_outdir}/tmp/l.fq" ]]; then
		_log_echo0 "bolap -c a.csv run data-long ${_brg_outdir} ${_brg_sindex}"
		run-data-long_genus_species "${_brg_outdir}" "${_brg_sindex}"
	fi
	_log_echo0 "bolap -c a.csv run summary-data ${_brg_outdir} ${_brg_inum_0}"
	# run-summary-data_genus_species "${_brg_outdir}" "${_brg_inum_0}"
	_log_echo0 "bolap -c a.csv run estimate-genomesize ${_brg_outdir} ${_brg_inum_0}"
	# run-estimate-genomesize_genus_species "${_brg_outdir}" "${_brg_inum_0}"
	_log_echo0 "bolap -c a.csv download ptdna ${_brg_outdir}"
	# download-ptdna_genus_species "${_brg_outdir}"
	_log_echo0 "bolap -c a.csv download mtdna ${_brg_outdir}"
	# download-mtdna_genus_species "${_brg_outdir}"
	_log_echo0 "bolap -c a.csv run ptgaul ${_brg_outdir}"
	# run-ptgaul_genus_species "${_brg_outdir}" "${_brg_inum_0}"
	_log_echo0 "bolap -c a.csv run mtgaul ${_brg_outdir}"
	# run-mtgaul_genus_species "${_brg_outdir}" "${_brg_inum_0}"
	_log_echo0 "bolap -c a.csv run tippo ${_brg_outdir}"
	run-tippo_genus_species "${_brg_outdir}" "${_brg_inum_0}"
	_log_echo0 "bolap -c a.csv run oatk ${_brg_outdir}"
	run-oatk_genus_species "${_brg_outdir}" "${_brg_inum_0}"
	_log_echo0 "bolap -c a.csv run polap-readassemble-pt" "${_brg_outdir}" "${_brg_sindex}"
	run-polap-readassemble-pt_genus_species "${_brg_outdir}" "${_brg_sindex}"
	_log_echo0 "bolap -c a.csv run polap-syncassemble ${_brg_outdir}"
	run-polap-syncassemble_genus_species "${_brg_outdir}" "${_brg_sindex}"
}

demo-syncfilter_genus_species_for() {
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"
	source "${_POLAPLIB_DIR}/polap-variables-data.sh"
	local _brg_inum_0="0"

	# local _brg_outdir="Salix_purpurea"
	# local _brg_sindex="0"

	rm -rf "${_brg_outdir}/tmp"
	rm -rf "${_brg_target}"

	# create a config
	if [[ ! -s "${_brg_outdir}/tmp/l.fq" ]]; then
		_log_echo0 "bolap -c a.csv run data-long ${_brg_outdir} ${_brg_sindex}"
		run-data-long_genus_species "${_brg_outdir}" "${_brg_sindex}"
	fi
	_log_echo0 "bolap -c a.csv run summary-data ${_brg_outdir} ${_brg_inum_0}"
	run-summary-data_genus_species "${_brg_outdir}" "${_brg_inum_0}"
	download-mtdna_genus_species "${_brg_outdir}" "${_brg_inum_0}"
	_log_echo0 "bolap -c a.csv run estimate-genomesize ${_brg_outdir} ${_brg_inum_0}"
	# run-estimate-genomesize_genus_species "${_brg_outdir}" "${_brg_inum_0}"
	_log_echo0 "bolap -c a.csv download ptdna ${_brg_outdir}"
	# download-ptdna_genus_species "${_brg_outdir}"
	_log_echo0 "bolap -c a.csv download mtdna ${_brg_outdir}"
	# download-mtdna_genus_species "${_brg_outdir}"
	_log_echo0 "bolap -c a.csv run ptgaul ${_brg_outdir}"
	# run-ptgaul_genus_species "${_brg_outdir}" "${_brg_inum_0}"
	_log_echo0 "bolap -c a.csv run mtgaul ${_brg_outdir}"
	# run-mtgaul_genus_species "${_brg_outdir}" "${_brg_inum_0}"
	_log_echo0 "bolap -c a.csv run tippo ${_brg_outdir}"
	# run-tippo_genus_species "${_brg_outdir}" "${_brg_inum_0}"
	_log_echo0 "bolap -c a.csv run oatk ${_brg_outdir}"
	# run-oatk_genus_species "${_brg_outdir}" "${_brg_inum_0}"
	_log_echo0 "bolap -c a.csv run polap-readassemble-pt" "${_brg_outdir}" "${_brg_sindex}"
	run-polap-readassemble-pt_genus_species "${_brg_outdir}" "${_brg_sindex}"
	_log_echo0 "bolap -c a.csv run polap-syncassemble ${_brg_outdir}"
	run-polap-syncassemble_genus_species "${_brg_outdir}" "${_brg_sindex}"

	rm -rf "${_brg_outdir}/tmp"
	rm -rf "${_brg_target}"
}

demo-syncfilter_genus_species() {
	local _brg_outdir="${1:-some}"

	if [[ "${_brg_outdir}" == "all" ]]; then
		for _v1 in "${Sall[@]}"; do
			demo-syncfilter_genus_species_for "${_v1}" "${@:2}"
		done
	elif [[ "${_brg_outdir}" == "test" ]]; then
		for _v1 in "${Stest[@]}"; do
			demo-syncfilter_genus_species_for "${_v1}" "${@:2}"
		done
	elif [[ "${_brg_outdir}" == "some" ]]; then
		for _v1 in "${Ssome[@]}"; do
			demo-syncfilter_genus_species_for "${_v1}" "${@:2}"
		done
	elif [[ "${_brg_outdir}" == "each" ]]; then
		for _v1 in "${Skeys[@]}"; do
			demo-syncfilter_genus_species_for "${_v1}" "${@:2}"
		done
	else
		demo-syncfilter_genus_species_for "$@"
	fi
}

demo-sim-syncassemble_genus_species() {
	# local _brg_outdir="$1"
	# local _brg_sindex="${2:-0}"
	# source "${_POLAPLIB_DIR}/polap-variables-data.sh"

	local _brg_outdir="Salix_purpurea"
	local _brg_sindex="0"
	local _brg_inum_0="0"

	# create a config
	if [[ ! -s "a.csv" ]]; then
		_log_echo0 "bolap config add ${_brg_outdir} ${_brg_sindex} SRR21824870 NA a.csv -> a.csv"
		config-add_genus_species "${_brg_outdir}" "${_brg_sindex}" SRR21824870 NA ONT a.csv
		echo "rerun it!"
		return
	fi
	_log_echo0 "bolap -c a.csv download ptdna ${_brg_outdir}"
	# download-ptdna_genus_species "${_brg_outdir}"
	_log_echo0 "bolap -c a.csv download mtdna ${_brg_outdir}"
	# download-mtdna_genus_species "${_brg_outdir}"
	_log_echo0 "bolap -c a.csv simulate ${_brg_outdir} ${_brg_sindex}"
	# simulate_genus_species "${_brg_outdir}" "${_brg_sindex}"
	_log_echo0 "bolap -c a.csv run tippo ${_brg_outdir}"
	# run-tippo_genus_species "${_brg_outdir}" "${_brg_inum_0}"
	_log_echo0 "bolap -c a.csv run oatk ${_brg_outdir}"
	# run-oatk_genus_species "${_brg_outdir}" "${_brg_inum_0}"
	_log_echo0 "bolap -c a.csv run polap-readassemble-pt" "${_brg_outdir}" "${_brg_sindex}"
	# run-polap-readassemble-pt_genus_species "${_brg_outdir}" "${_brg_sindex}"
	_log_echo0 "bolap -c a.csv run polap-syncassemble ${_brg_outdir}"
	run-polap-syncassemble_genus_species "${_brg_outdir}" "${_brg_sindex}"
}

config-test_genus_species() {

	# print_species_field_summary --add-field=fruit=banana --values
	print_species_field_summary --values --fields=long,short
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

	rsync -azuq --max-size=5M "${_brg_rundir}"/ "${_brg_outdir_i}/${_brg_title}"/
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
	local _brg_ranger="${4:-0}"
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
	local -n range="range${_brg_ranger}"

	local min_read="${_min_read["$_brg_target"]}"
	local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

	_polap_lib_conda-ensure_conda_env polap || exit 1

	rm -f "${_timing_txt}"
	rm -f "${_stdout_txt}"

	echo "[INFO] Start an organelle-genome assembly pipeline on ${_brg_outdir}-${_brg_sindex} on a range of omega values: ${range} and minimum read length of ${min_read}"

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
	local _brg_range="${5:-0}"
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
	local -n range="range${_brg_range}"

	local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

	_polap_lib_conda-ensure_conda_env polap || exit 1

	rm -f "${_timing_txt}"
	rm -f "${_stdout_txt}"

	echo "[INFO] Start an organelle-genome assembly pipeline on ${_brg_outdir}-${_brg_sindex} on a range of omega values: ${range}"

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
	local _brg_range="${5:-0}"
	local _brg_knum="${6:--1}"
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
	local -n range="range${_brg_range}"
	local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

	echo "[INFO] Start an organelle-genome assembly pipeline on ${_brg_outdir}-${_brg_sindex} on a range of omega values: ${range}"

	_polap_lib_conda-ensure_conda_env polap || exit 1

	if [[ "${_brg_knum}" == "-1" ]]; then
		rm -rf "${_brg_rundir}/${_brg_jnum}/05-flye/ptgaul"
		mkdir -p "${_brg_rundir}/${_brg_jnum}/05-flye/ptgaul"
		${_polap_cmd} directional-flye-reads \
			-o ${_brg_rundir} \
			-i "${_brg_inum}" -j "${_brg_jnum}" \
			--no-directional \
			--no-coverage-check \
			-v -v \
			-s "${range}"
	else
		local _next_knum=$((_brg_knum + 1))
		rm -rf "${_brg_rundir}/${_brg_jnum}/05-flye/ptgaul/${_brg_knum}"
		${_polap_cmd} directional-flye-reads \
			-o ${_brg_rundir} \
			-i "${_brg_inum}" -j "${_brg_jnum}" \
			--no-directional \
			--no-coverage-check \
			-s "${range}" --start-index "${_brg_knum}" --end-index "${_next_knum}"
	fi

	conda deactivate
}

run-direct-dflye_genus_species() {
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"
	local _brg_inum="${3:-1}"
	local _brg_jnum="${4:-2}"
	local _brg_range="${5:-0}"
	local _brg_knum="${6:--1}"
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
	local -n range="range${_brg_range}"
	local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

	echo "[INFO] Start an organelle-genome assembly pipeline on ${_brg_outdir}-${_brg_sindex} on a range of omega values: ${range}"

	_polap_lib_conda-ensure_conda_env polap-dflye || exit 1

	if [[ "${_brg_knum}" == "-1" ]]; then
		rm -rf "${_brg_rundir}/${_brg_jnum}/05-dflye/ptgaul"
		mkdir -p "${_brg_rundir}/${_brg_jnum}/05-dflye/ptgaul"
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
	local _brg_range="${5:-0}"

	# run-direct-select-oga_genus_species
	run-direct-seed_genus_species "${_brg_outdir}" "${_brg_sindex}" "${_brg_inum}" "${_brg_jnum}"
	run-direct-map_genus_species "${_brg_outdir}" "${_brg_sindex}" "${_brg_inum}" "${_brg_jnum}"
	run-direct-read_genus_species "${_brg_outdir}" "${_brg_sindex}" "${_brg_inum}" "${_brg_jnum}" "${_brg_range}"
	run-direct-flye_genus_species "${_brg_outdir}" "${_brg_sindex}" "${_brg_inum}" "${_brg_jnum}" "${_brg_range}"
	# run-direct-dflye_genus_species "${_brg_outdir}" "${_brg_sindex}" "${_brg_inum}" "${_brg_jnum}"
}

conda_genus_species() {
	source "${_POLAPLIB_DIR}/polap-variables-data.sh"

	local bolap_cmd="${FUNCNAME%%_*}"

	help_message=$(
		cat <<EOF
Name:
  bolap - $bolap_cmd

Synopsis:
  bolap $bolap_cmd

Description:
  Clean-up conda environments. Recreate them.
    bolap $bolap_cmd --cleanup
    bolap $bolap_cmd --export polap
    bolap $bolap_cmd --export-all
    bolap $bolap_cmd --recreate
    bolap $bolap_cmd --delete polap
    bolap $bolap_cmd --create polap
    bolap $bolap_cmd --delete polap --create polap
    bolap $bolap_cmd --list
    bolap $bolap_cmd --table

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	local _export=""
	local _export_all=0
	local _delete=""
	local _cleanup=0
	local _recreate=0
	local _create=""
	local _list=0
	local _table=0

	parse_commandline() {
		set -- "${_brg_args[@]}"

		# source "${_POLAPLIB_DIR}/polap-cmd-version.sh" # '.' means 'source'
		while test $# -gt 0; do
			_key="$1"
			# _log_echo0 "$#: $_key"
			case "$_key" in
			--export-all)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_export_all=1
				fi
				;;
			--export)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_export="$2"
				fi
				shift || true
				;;
			--delete)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_delete="$2"
				fi
				shift || true
				;;
			--cleanup)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_cleanup="1"
				fi
				;;
			--list)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_list="1"
				fi
				;;
			--table)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_table="1"
				fi
				;;
			--recreate)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_recreate="1"
				fi
				;;
			--create)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_create="$2"
				fi
				shift || true
				;;
			# -*)
			# 	_log_echo0 "[ERROR] no such options: $1"
			# 	;;
			# *)
			# 	_log_echo0 "[POS] arg: $1"
			# 	;;
			esac
			shift || true
		done
		# _log_echo0 "p:cleanup: $_cleanup"
		# _log_echo0 "p:export: $_export"
		# _log_echo0 "p:recreate: $_recreate"
	}

	declare -n ref="help_message"
	if [[ "${_brg_help}" == "on" ]]; then
		local manfile=$(_bolap_lib_man-convert_help_message "$ref" "${bolap_cmd}")
		man "$manfile"
		rm -f "$manfile"
		return
	fi

	parse_commandline

	# _log_echo0 "cleanup: $_cleanup"
	# _log_echo0 "export: $_export"
	# _log_echo0 "table: $_table"
	# _log_echo0 "delete: $_delete"
	# _log_echo0 "export-all: $_export_all"
	# _log_echo0 "recreate: $_recreate"
	# _log_echo0 "create: $_create"

	_polap_lib_conda-ensure_conda_env base || exit 1

	if [[ "$_cleanup" == "1" ]]; then
		if _polap_lib_dialog-yes-no "Delete all polap conda environments?" "y"; then
			make -f "${_POLAPLIB_DIR}"/Makefile.env clean-envs
		else
			echo "[SKIP] User declined."
		fi
	fi

	if [[ -n "$_delete" ]]; then
		if _polap_lib_dialog-yes-no "Delete conda environments: $_delete?" "y"; then
			make -f "${_POLAPLIB_DIR}"/Makefile.env clean-envs ENVS="$_delete"
		else
			echo "[SKIP] User declined."
		fi
	fi

	if [[ "$_export_all" == "1" ]]; then
		if _polap_lib_dialog-yes-no "Export all conda environments?" "y"; then
			make -f "${_POLAPLIB_DIR}"/Makefile.env export-full
		else
			echo "[SKIP] User declined."
		fi
	fi

	if [[ -n "$_export" ]]; then
		if _polap_lib_dialog-yes-no "Export a conda environments: $_export?" "y"; then
			make -f "${_POLAPLIB_DIR}"/Makefile.env export ENVS="$_export"
		else
			echo "[SKIP] User declined."
		fi
	fi

	if [[ -n "$_create" ]]; then
		if _polap_lib_dialog-yes-no "Recreate polap conda environments: $_create?" "y"; then
			make -f "${_POLAPLIB_DIR}"/Makefile.env recreate ENVS="$_create"
			# make -f "${_POLAPLIB_DIR}"/Makefile.env recreate-"$_create"
		else
			echo "[SKIP] User declined."
		fi
	fi

	if [[ "$_recreate" == "1" ]]; then
		if _polap_lib_dialog-yes-no "Recreate all polap conda environments?" "y"; then
			make -f "${_POLAPLIB_DIR}"/Makefile.env recreate-skip
		else
			echo "[SKIP] User declined."
		fi
	fi

	if [[ "$_list" == "1" ]]; then
		conda env list | grep polap
	fi

	if [[ "$_table" == "1" ]]; then
		if _polap_lib_dialog-yes-no "List all conda packages (it takes time; be patient!)?" "y"; then

			# Concatenate conda env yaml files.
			local output_file="man/md/conda-packages.yml"
			>"$output_file"

			# Find all environment.yml files (depth = 2)
			find "${_POLAPLIB_DIR}/envs" -mindepth 2 -maxdepth 2 -type f -name "environment.yml" | sort | while read -r file; do
				dir_name=$(basename "$(dirname "$file")")
				echo "# ───────────────────────────────────────────────" >>"$output_file"
				echo "# From: $dir_name/environment.yml" >>"$output_file"
				echo "# ───────────────────────────────────────────────" >>"$output_file"
				cat "$file" >>"$output_file"
				echo -e "\n" >>"$output_file"
			done

			python "${_POLAPLIB_DIR}/scripts/conda-packages.py" --simple-channel \
				"${output_file}" >man/md/conda-packages.csv

			pandoc -f csv -t markdown_mmd man/md/conda-packages.csv -o man/md/conda-packages.md

			echo "[INFO] man/md/conda-packages.csv"
			echo "[INFO] man/md/conda-packages.md"

			# make -f "${_POLAPLIB_DIR}"/Makefile.env recreate-skip
		else
			echo "[SKIP] User declined."
		fi
	fi

	conda deactivate
}

conda-clean_genus_species() {
	source "${_POLAPLIB_DIR}/polap-variables-data.sh"

	local bolap_cmd="${FUNCNAME%%_*}"

	help_message=$(
		cat <<EOF
Name:
  bolap - $bolap_cmd

Synopsis:
  bolap $bolap_cmd

Description:
  It deletes all polap-related environments.

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	declare -n ref="help_message"
	if [[ "${_brg_help}" == "on" ]]; then
		local manfile=$(_bolap_lib_man-convert_help_message "$ref" "${bolap_cmd}")
		man "$manfile"
		rm -f "$manfile"
		return
	fi

	make -f "${_POLAPLIB_DIR}"/Makefile.env clean-envs
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

download-test-data-hifi_genus_species() {

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

download-test-data-read_genus_species() {

	if [[ "${_brg_outdir}" == "-h" || "${_brg_outdir}" == "--help" ]]; then
		echo "$help_message_download_test_data_cflye"
		return
	fi

	if [[ "${opt_y_flag}" == false ]]; then
		read -p "Do you want to download test data for polap-read? (y/N): " confirm
	else
		confirm="yes"
	fi

	if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
		local _s=polap-readassemble-test-full.tar.gz
		# local _s=polap-disassemble-test.tar.gz
		if [[ ! -s "${_s}" ]]; then
			# full test data
			# curl -L -o "${_s}" "https://figshare.com/ndownloader/files/53457569?private_link=ec1cb394870c7727a2d4"
			wget -O "${_s}" "https://figshare.com/ndownloader/files/57379570?private_link=2ee3e102cadcf0d28bbb"
			#
			# test data
			# curl -L -o "${_s}" "https://figshare.com/ndownloader/files/53457566?private_link=ec1cb394870c7727a2d4"
		fi
		if [[ ! -s "anthoceros_angustus_long.fq" ]]; then
			tar -zxf "${_s}"
			echo "downloaded: anthoceros_angustus_long.fq, anthoceros_angustus_short_1.fq, anthoceros_angustus_short_2.fq"
		else
			echo "You already have: anthoceros_angustus_long.fq"
		fi
	else
		echo "polap read-assemble test download is canceled."
	fi
}

cleanup-all_genus_species() {
	local first_arg="${1:-read}"
	local _brg_sindex="${2:-$_brg_sindex}"

	local odir
	for odir in "${_brg_outdir_list[@]}"; do
		echo "Deleting all folders from ${odir} ..."
		local _brg_outdir="$odir"
		source "${_POLAPLIB_DIR}/polap-variables-data.sh"
		local long_sra="${_long["$_brg_target"]}"
		local short_sra="${_short["$_brg_target"]}"
		set -x
		rm -rf "${odir}"
		rm -rf "${odir}-${_brg_sindex}"
		rm -f "$long_sra"*
		rm -f "$short_sra"*
		set +x
	done
}

clean_genus_species() {
	local first_arg="${1:-read}"

	local outdir
	for outdir in "${_brg_outdir_list[@]}"; do
		clean-"${first_arg}"_genus_species_for "${outdir}"
	done
}

clean-read_genus_species_for() {
	local odir="${1}"

	source "${_POLAPLIB_DIR}/polap-variables-data.sh"

	if [[ "${odir}" == "-h" || "${odir}" == "--help" ]]; then
		echo "$help_message_clean_read"
		return 0
	fi

	odir="${odir%/}"

	if [[ ! -d "${odir}" ]]; then
		echo "[ERROR] no such folder: ${odir}"
		return 1
	fi

	local target_index="${odir}-0"

	local long_sra="${_long["$target_index"]}"
	if [[ -z "$long_sra" ]]; then
		echo "Error: skipping ${odir}-${_brg_inum} because it is not in the CSV."
		return 1
	fi

	local long_sra1=$(dataset-get-field "$_brg_target" long)
	local short_sra1=$(dataset-get-field "$_brg_target" short)

	echo "Deleting some folders from ${_brg_target} ..."
	set -x
	rm -f "${odir}/tmp"/*.gz
	rm -rf "${_brg_target}"
	rm -f "$long_sra1".{fastq,fq,mt.gfa}
	rm -f "$long_sra1".*
	rm -f "$short_sra1".{fastq,fq,mt.gfa}
	rm -f "$short_sra1".*
	rm -f 2-oga.gfa
	set +x
}

clean-hifi_genus_species() {
	local _brg_outdir="${1}"

	if [[ "${_brg_outdir}" == "-h" || "${_brg_outdir}" == "--help" ]]; then
		echo "$help_message_clean_hifi"
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

	echo "Deleting input/output files and folders from ${_brg_outdir} ..."
	rm -rf "${_brg_outdir}"
	rm -rf "${_brg_outdir}-0"
	rm -rf "${_brg_outdir}-1"
	rm -rf "${_brg_outdir}-2"
	rm -f "${long_sra}.fastq"
	rm -f "${short_sra}_1.fastq"
	rm -f "${short_sra}_2.fastq"
	rm -f "${short_sra}.fastq.tar.gz"
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
	local local_path="${1:-polap}"

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

	if declare -F update-${item}_genus_species >/dev/null; then
		update-${item}_genus_species ${rest[@]}
	else
		echo "$help_message_update"
	fi
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
	local _brg_outdir="${1:-example}"

	local subcmd1="${_brg_outdir}"
	local _subcmd1_clean="${subcmd1//-/_}"

	local target_var="help_message_${_subcmd1_clean}"

	if declare -p "$target_var" &>/dev/null; then
		declare -n ref="$target_var"
		echo "$ref"
	else
		echo "[ERROR] No such bolap help or subcommand: $subcmd1" >&2
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
				echo "[INFO] Conda environment 'polap-fmlrc2' already exists."
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

# 2025-12-03
# 2025-12-04
install-hifi_genus_species() {
	local bolap_cmd="${FUNCNAME%%_*}"

	help_message=$(
		cat <<EOF
Name:
  bolap - $bolap_cmd

Synopsis:
  bolap $bolap_cmd

Description:
  Install all conda environments for bolap hifi analysis.

Examples:
  Install all for bolap hifi analysis:
    bolap install hifi

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	# defaults
	local _brg_ver=""

	_polap_lib_help-maybe-show "$bolap_cmd" help_message || return 0

	polap_setup_parse_commandline

	local _version="${_brg_ver}"

	local tools_to_install=(
		ncbitools
		fmlrc
		fmlrc2
		oatk
		himt
		man
		tippo
		pmat2
		bandage
	)

	set +u
	for item in "${tools_to_install[@]}"; do
		install-${item}_genus_species
	done
	set -u
}

setup-nvim_genus_species() {
	set +e

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
	set -e
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

setup-pmat2-apptainer_genus_species() {
	echo "Install fuse2fs gocryptfs ..."
	sudo apt update
	sudo apt -y install fuse2fs gocryptfs

	echo "Setup apptainer ..."
	_polap_lib_conda-ensure_conda_env polap-pmat2 || exit 1
	apptainer exec $CONDA_PREFIX/bin/container/runAssembly.sif echo "Fail: PMAT setup"
	echo 'kernel.unprivileged_userns_clone=1' | sudo tee /etc/sysctl.d/90-userns.conf
	echo 'kernel.apparmor_restrict_unprivileged_userns=0' | sudo tee /etc/sysctl.d/80-apparmor-userns.conf
	sudo sysctl --system
	apptainer exec $CONDA_PREFIX/bin/container/runAssembly.sif echo "Success: PMAT setup"
	set +u
	conda deactivate
	set -u
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
alias pl='bash $TEMPLATE_DIR/polap.sh'
alias bl='bash $TEMPLATE_DIR/bolap.sh'
alias pl2='bash $TEMPLATE_DIR/polap2.sh'
alias bl2='bash $TEMPLATE_DIR/bolap2.sh'
alias p1='bash $TEMPLATE_DIR/polap-data-v1.sh'
alias p2='bash $TEMPLATE_DIR/polap-data-v2.sh'
alias p3='bash $TEMPLATE_DIR/polap-data-v3.sh'
alias p4='bash $TEMPLATE_DIR/polap-data-v4.sh'
alias p5='bash $TEMPLATE_DIR/polap-data-v5.sh'
alias p6='bash $TEMPLATE_DIR/polap-data-v6.sh'
alias p7='bash $TEMPLATE_DIR/polap-data-v7.sh'
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
export QT_QPA_PLATFORM=minimal
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
	echo "Removing the conda environments: polap"
	if [[ "${opt_y_flag-}" == false ]]; then
		read -p "Do you want to uninstall polap conda environments? (y/N): " confirm
	else
		confirm="yes"
	fi

	if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
		_polap_lib_conda-ensure_conda_env base || exit 1
		if [[ "$CONDA_DEFAULT_ENV" == "base" ]]; then
			conda env remove -y -n polap
		else
			echo "Error: You're in the '$CONDA_DEFAULT_ENV' environment. Please activate base before running this script."
			exit 1
		fi
	else
		echo "Uninstallation of polap is canceled."
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
	local _brg_outdir="${1:-$_brg_outdir}"
	local _brg_sindex="${2:-$_brg_sindex}"

	# _log_echo0 "_brg_outdir: $_brg_outdir"
	# _log_echo0 "_brg_sindex: $_brg_sindex"
	local bolap_cmd="${FUNCNAME##*_}"

	help_message=$(
		cat <<EOF
Name:
  bolap - summary-data

Synopsis:
  bolap run-summary-data

Description:
  bolap

Examples:
  Execute:
    bolap run-summary-data -s Vigna_radiata

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	local _data_long_cleanup=false
	local _data_long_redo=false

	parse_commandline() {
		set -- "${_brg_unknown_opts[@]}"

		# source "${_POLAPLIB_DIR}/polap-cmd-version.sh" # '.' means 'source'
		while test $# -gt 0; do
			_key="$1"
			case "$_key" in
			--data-long-cleanup)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'."
				else
					_data_long_cleanup="$2"
					shift || true
				fi
				;;
			--data-long-redo)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'."
				else
					_data_long_redo="$2"
					shift || true
				fi
				;;
			-*)
				_log_echo0 "[ERROR] no such options: $1"
				;;
			*)
				break
				;;
			esac
			shift || true
		done
	}

	declare -n ref="help_message"
	if [[ "${_brg_help}" == "on" ]]; then
		local manfile=$(_bolap_lib_man-convert_help_message "$ref" "${bolap_cmd}")
		man "$manfile"
		rm -f "$manfile"
		return
	fi

	parse_commandline

	# source "${_POLAPLIB_DIR}/polap-variables-data.sh"
	# local long_sra="${_long["$_brg_target"]}"

	# redownload
	source "${_POLAPLIB_DIR}/polap-variables-data.sh"

	mkdir -p "${_brg_rundir}"

	local short_sra="${_short["$_brg_target"]}"

	_polap_lib_conda-ensure_conda_env polap || exit 1

	_log_echo0 _polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"
	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"

	rm -f "${_stdout_txt}"
	rm -f "${_timing_txt}"

	if [[ -s "${_brg_tmpdir}/l.fq.gz" ]]; then
		_polap_lib_seqkit-stats \
			--in "${_brg_tmpdir}/l.fq.gz" \
			--out "${_brg_rundir}/l.fq.seqkit.stats.ta.txt"

		ln -sf \
			$(realpath --relative-to "${_brg_tmpdir}" "${_brg_rundir}/l.fq.seqkit.stats.ta.txt") \
			"${_brg_tmpdir}/l.fq.seqkit.stats.ta.txt"

	fi

	local i
	for i in 1 2; do
		if [[ -s "${_brg_tmpdir}/s_$i.fq.gz" ]]; then
			_polap_lib_seqkit-stats \
				--in "${_brg_tmpdir}/s_$i.fq.gz" \
				--out "${_brg_rundir}/s_$i.fq.seqkit.stats.ta.txt"

			ln -sf \
				$(realpath --relative-to "${_brg_tmpdir}" "${_brg_rundir}/s_$i.fq.seqkit.stats.ta.txt") \
				"${_brg_tmpdir}/s_$i.fq.seqkit.stats.ta.txt"
		fi
	done

	# Save system info
	# _polap_lib_timing-get_system_info >>"${_timing_txt}"

	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"

	conda deactivate
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
				echo "[INFO] Conda environment 'polap-man' already exists."
			else
				conda create -y --name polap-man pandoc pandoc-crossref yq jq \
					r-ggplot2 r-patchwork r-dplyr r-readr r-optparse
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
	run-data_genus_species
	run-downsample_genus_species
	run-polap-assemble_genus_species
	run-polap-coverage_genus_species

	if [[ "${_local_host}" != "$(hostname)" ]]; then
		sync_genus_species
	fi
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
	local _brg_outdir="$1"
	local _brg_sindex="${2:-0}"
	source "${_POLAPLIB_DIR}/polap-variables-data.sh"

	# Initialize the run folders
	rm -rf "${_brg_target}"
	mkdir -p "${_brg_target}"

	# Input data
	local short_sra="${_short["$_brg_target"]}"

	_polap_lib_conda-ensure_conda_env polap || exit 1

	# Start memory logger
	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"

	if [[ -s "${_brg_tmpdir}/s_1.fq" ]]; then
		command time -v ${_polap_cmd} find-genome-size \
			-a ${_brg_tmpdir}/s_1.fq \
			-b ${_brg_tmpdir}/s_2.fq \
			-o "${_brg_target}" \
			>"${_stdout_txt}" \
			2>"${_timing_txt}"
	elif [[ -s "${_brg_tmpdir}/l.fq" ]]; then
		command time -v ${_polap_cmd} find-genome-size-for-pacbio \
			-l ${_brg_tmpdir}/l.fq \
			-o "${_brg_target}" \
			>"${_stdout_txt}" \
			2>"${_timing_txt}"
	else
		echo "No data to estimate the genome size" >"${_stdout_txt}"
	fi

	# End with summary of the system usage
	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"

	# Save system info
	_polap_lib_timing-get_system_info >>"${_timing_txt}"

	conda deactivate

	# Save results
	rsync -azuq --max-size=5M "${_brg_target}"/ "${_brg_rundir}"/
	rm -rf "${_brg_target}"
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
		command time -v timeout ${_brg_timeout} PMAT autoMito \
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

# 2025-09-04
#
# Usage: /home/goshng/miniconda3/envs/polap-tippo/bin/TIPPo.v2.4.pl [options]
#   -h: Show this help message.
#   -f: Long reads (required).
#   -d: reference sequence (default: None).
#   -g: chloroplast or organelle (default: organelle).
#   -t: Threads for tiara, flye, KMC3 and readskmercount.
#   -n: Number of reads in each downsample for chloroplast.
#   -r: Number of random downsamplings (default: 5).
#   -p: Sequence technology - 'hifi','clr','ont', 'onthq' (default: hifi)
#   -i: Assume the presence of the inverted repeats in the chloroplast genome (default: 1).
#   -l: lower kmer count - lkc (default: 0.3).
#   -c: high kmer count - hkc (default: 5).
#   -m: minimum overlap in repeat graph construction (default:800)
#   --trf: remove the reads are tandem repeats, only avaliable for reference-free and hifi/onthq reads
#   -v: version.
run-tippo_genus_species() {
	local _brg_outdir="${1:-$_brg_outdir}"
	local _brg_sindex="${2:-$_brg_sindex}"

	source "${_POLAPLIB_DIR}/polap-variables-data.sh"

	# Always redo
	rm -rf "${_brg_target}"
	mkdir -p "${_brg_target}"

	# Always redo
	rm -rf "${_brg_rundir}"
	mkdir -p "${_brg_rundir}"

	# ln -sf "${_brg_tmpdir}/l.fq" "${long_sra}.fastq"
	_brg_input_data="${long_sra}.fastq"

	if [[ ! -s "${_brg_input_data}" ]]; then
		_log_echo0 "[ERROR] no inut data: ${_brg_input_data}"
		return 1
	fi

	# TIPPo has no option for output directory
	local _tippo_actual_outdir="${_brg_input_data}.organelle"

	# https://github.com/Wenfei-Xian/TIPP
	rm -f "${_stdout_txt}"
	rm -f "${_timing_txt}"

	# Start memory logger
	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"

	local platform="${_platform["$_brg_target"]}"

	if [[ "$platform" == "ONT" ]]; then
		local _fc="ont"
	elif [[ "$platform" == "ONTHQ" ]]; then
		local _fc="onthq"
	elif [[ "${platform}" == "PACBIO_SMRT" ]]; then
		local _fc="hifi"
	elif [[ "${platform}" == "PACBIO_CLR" ]]; then
		local _fc="clr"
	else
		_log_echo "No such platform: $platform"
		return 2
	fi

	rm -rf "${_tippo_actual_outdir}"

	_polap_lib_conda-ensure_conda_env polap-tippo || exit 1

	# -g organelle
	# -g chloroplast \
	command time -v timeout ${_brg_timeout} TIPPo.v2.4.pl \
		-f "${_brg_input_data}" \
		-t ${_brg_threads} \
		-g organelle \
		-p ${_fc} \
		>"${_stdout_txt}" \
		2>"${_timing_txt}" || rc=$?

	if [[ "${rc:-0}" -ne 0 ]]; then
		printf '[ERR] oatk failed (rc=%d)\n' "$rc" >&2
		# decide: exit "$rc" / continue
	fi

	conda deactivate

	# Save system info
	_polap_lib_timing-get_system_info >>"${_timing_txt}"

	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"

	# Save results
	rsync -azuq --max-size=5M \
		"${_tippo_actual_outdir}/" "${_brg_rundir}/"

	# Output files
	# cns.fa.chloroplast.fasta.filter.800.round1.edge_*.edge_*.edge_*.organelle.chloroplast.fasta
	# cns.fa.mitochondrial.fasta.filter.fasta.flye/assembly_graph.gfa

	# Clean-up
	rm -rf "${_brg_target}"
	rm -rf "${_tippo_actual_outdir}"
	rm -f "${long_sra}.fastq.fasta"
}

# 2025-09-04
run-oatk_genus_species() {
	local _brg_outdir="${1:-$_brg_outdir}"
	local _brg_sindex="${2:-$_brg_sindex}"

	source "${_POLAPLIB_DIR}/polap-variables-data.sh"

	# Create folders
	rm -rf "${_brg_target}"
	mkdir -p "${_brg_target}"
	mkdir -p "${_brg_outdir_i}"

	# ln -sf "${_brg_tmpdir}/l.fq" "${long_sra}.fastq"
	_brg_input_data="${long_sra}.fastq"

	# check the input
	# prepare if there is no such one

	_polap_lib_conda-ensure_conda_env polap-oatk || exit 1

	local _fc
	# for _fc in 30 25 20 15 10; do
	for _fc in 30 20 10; do

		local formatted_fc=$(printf "%02d" "${_fc}")

		_log_echo "oatk on ${_brg_rundir} with -c ${_fc}"
		local formatted_fc=$(printf "%02d" "${_fc}")
		local stdout_txt_oatk="${_stdout_txt%.txt}-${formatted_fc}.txt"
		local timing_txt_oatk="${_timing_txt%.txt}-${formatted_fc}.txt"
		local memlog_file_oatk="${_memlog_file%.csv}-${formatted_fc}.csv"
		local summary_file_oatk="${_summary_file%.txt}-${formatted_fc}.txt"

		rm -f "${stdout_txt_oatk}"
		rm -f "${timing_txt_oatk}"

		# Start memory logger
		_polap_lib_process-start_memtracker "${memlog_file_oatk}" \
			"${_polap_var_oatk_memtracker_time_interval}"

		# https://github.com/c-zhou/oatk
		# NOTE: oatk threads 56 -> seg. fault.
		# Either way, it is still very fast.
		local rc
		local _oatk_threads=8
		command time -v timeout ${_brg_timeout} oatk \
			-k 1001 \
			-c ${_fc} \
			-t ${_oatk_threads} \
			-m ./OatkDB/v20230921/embryophyta_mito.fam \
			-p ./OatkDB/v20230921/embryophyta_pltd.fam \
			-o "${_brg_target}/oatk-${formatted_fc}" \
			"${_brg_input_data}" \
			>"${stdout_txt_oatk}" \
			2>"${timing_txt_oatk}" || rc=$?

		if [[ "${rc:-0}" -ne 0 ]]; then
			printf '[ERR] oatk failed (rc=%d)\n' "$rc" >&2
			# decide: exit "$rc" / continue
		fi

		# Save system info
		_polap_lib_timing-get_system_info >>"${timing_txt_oatk}"

		_polap_lib_process-end_memtracker "${memlog_file_oatk}" "${summary_file_oatk}" "no_verbose"

	done

	conda deactivate

	# Save results
	rsync -azuq --max-size=5M \
		"${_brg_target}/" "${_brg_rundir}/"

	# Clean-up
	rm -rf "${_brg_target}"
}

run-oatk-hifi_genus_species() {
	local _brg_outdir="${1}"
	local _brg_sindex="${2:-0}"
	local _brg_type="${3:-hifi}"
	local _run_title="oatk-${_brg_type}"

	# Directory without a slash
	_brg_outdir="${_brg_outdir%/}"

	if [[ "${_brg_outdir}" == "-h" || "${_brg_outdir}" == "--help" ]]; then
		echo "$help_message_run_oatk"
		return
	fi

	local _brg_inum=0
	local _brg_adir _brg_title _brg_target _brg_rundir _brg_outdir_i
	local _timing_txt _stdout_txt _memlog_file _summary_file
	brg_common_setup \
		_brg_outdir _brg_sindex _brg_adir _brg_title \
		_brg_target _brg_rundir _brg_outdir_i \
		_timing_txt _stdout_txt _memlog_file _summary_file
	# Extra folders
	local _brg_outdir_t="${_brg_outdir}/${opt_t_arg}"
	local _brg_titledir="${_brg_outdir_i}/${_brg_title}"
	local _brg_runtitledir="${_brg_rundir}-${_brg_title}"

	if [[ -v _long["$_brg_target"] ]]; then
		local long_sra="${_long["$_brg_target"]}"
	else
		echo "Error: ${_brg_target} because it is not in the CSV."
		return
	fi

	# Folders
	# brg: outdir_t -> outdir_i -> rundir
	# run_dir
	local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

	# debug local variables
	# for debugging: Inline local printing local var
	# while IFS= read -r line; do
	# 	if [[ $line =~ ^declare\ --\ ([^=]+)= ]]; then
	# 		var="${BASH_REMATCH[1]}"
	# 		printf "%s=%q\n" "$var" "${!var}"
	# 	fi
	# done < <(local -p 2>/dev/null)
	# return

	# Create folders
	mkdir -p "${_brg_titledir}"
	rm -rf "${_brg_runtitledir}"
	mkdir -p "${_brg_runtitledir}"

	local _brg_input_data="${long_sra}.fastq"

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
		_log_echo "oatk on ${_brg_rundir} with -c ${_fc}"
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
		command time -v timeout ${_brg_timeout} oatk \
			-k 1001 \
			-c ${_fc} \
			-t ${_oatk_threads} \
			-m ./OatkDB/v20230921/embryophyta_mito.fam \
			-p ./OatkDB/v20230921/embryophyta_pltd.fam \
			-o "${_brg_runtitledir}/oatk-${_brg_type}-${formatted_fc}" \
			"${_brg_input_data}" \
			>"${_stdout_txt}" \
			2>"${_timing_txt}"

		# Save system info
		_polap_lib_timing-get_system_info >>"${_timing_txt}"

		_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"

	done

	conda deactivate

	# Save results
	rsync -azuq "${_brg_runtitledir}"/ "${_brg_titledir}"/

	# Clean-up
	rm -rf "${_brg_runtitledir}"
}

# Run Tiara and return mito FASTQ path
# Args:
#   $1 = input reads (FASTQ)
#   $2 = work dir (e.g. "${_brg_target}")
#   $3 = name of shell variable to set with mito FASTQ path
#
# Side products:
#   ${work}/tiara/tiara.out.tsv
#   ${work}/tiara/tiara.mito.fq.gz
#   ${work}/tiara/tiara.plastid.fq.gz
_polap_oatk_step_tiara() {
	local in_reads="$1"
	local work="$2"
	local __outvar="$3"

	local _stdout_txt="${_stdout_txt_base}-tiara.txt"
	local _timing_txt="${_timing_txt_base}-tiara.txt"

	_polap_lib_conda-ensure_conda_env polap-tiara || return 1

	local tiara_dir="${work}/tiara"
	mkdir -p "$tiara_dir"

	local tiara_out="${tiara_dir}/tiara.out.tsv"
	local tiara_mito_fq="${tiara_dir}/tiara.mito.fq.gz"
	local tiara_plastid_fq="${tiara_dir}/tiara.plastid.fq.gz"
	local tiara_fa="${tiara_dir}/tiara.input.fa"

	echo "[INFO] Tiara: converting FASTQ to FASTA for $in_reads"
	seqkit fq2fa "$in_reads" |
		seqkit replace -p "\s.+" -r "" >"$tiara_fa"

	echo "[INFO] Tiara classification on $in_reads"
	command time -v tiara \
		-i "$tiara_fa" \
		-o "$tiara_out" \
		--tf mit pla \
		-t "${_brg_threads:-8}" \
		>"${_stdout_txt}" \
		2>"${_timing_txt}" || rc=$?

	if [[ "${rc:-0}" -ne 0 ]]; then
		printf '[ERR] tiara failed (rc=%d)\n' "$rc" >&2
	fi

	# mitogenome reads: class_snd_stage == "mitochondrion"
	awk 'BEGIN{FS="\t"} NR>1 && $3=="mitochondrion" {print $1}' "$tiara_out" |
		seqkit grep -f - "$in_reads" \
			-o "$tiara_mito_fq"

	# plastome reads: class_snd_stage == "plastid"
	awk 'BEGIN{FS="\t"} NR>1 && $3=="plastid" {print $1}' "$tiara_out" |
		seqkit grep -f - "$in_reads" \
			-o "$tiara_plastid_fq"

	rm -f "$tiara_fa"
	_polap_lib_conda-ensure_conda_env_deactivate || true

	# return mito FASTQ path in the requested variable
	printf -v "$__outvar" '%s' "$tiara_mito_fq"
}

# Run HiMT filter and return filtered FASTQ path
# Args:
#   $1 = input reads (FASTQ)
#   $2 = work dir (e.g. "${_brg_target}")
#   $3 = name of shell variable to set with filtered FASTQ path
_polap_oatk_step_himt() {
	local in_reads="$1"
	local work="$2"
	local __outvar="$3"

	local _stdout_txt="${_stdout_txt_base}-himt-filter.txt"
	local _timing_txt="${_timing_txt_base}-himt-filter.txt"

	_polap_lib_conda-ensure_conda_env polap-himt || return 1

	local himt_dir="${work}/himt"
	mkdir -p "$himt_dir"

	local himt_out_fa="${himt_dir}/extract.fa"
	local himt_out="${himt_dir}/himt.fq.gz"

	echo "[INFO] HiMT filter on $in_reads"
	command time -v himt filter \
		--input_file "$in_reads" \
		--output_dir "$himt_dir" \
		-s plant \
		--thread "${_brg_threads:-8}" \
		>"${_stdout_txt}" \
		2>"${_timing_txt}" || rc=$?

	if [[ "${rc:-0}" -ne 0 ]]; then
		printf '[ERR] himt filter failed (rc=%d)\n' "$rc" >&2
	fi

	_polap_lib_conda-ensure_conda_env_deactivate || conda deactivate || true

	seqkit seq -ni "${himt_out_fa}" |
		seqkit grep -f - "${in_reads}" -o "${himt_out}"

	printf -v "$__outvar" '%s' "$himt_out"
}

# 2025-12-07
################################################################################
# Run Oatk on HiFi reads with Tiara + HiMT filtering in configurable order
#
# Patterns:
#   t-h : Tiara -> HiMT -> Oatk        -> prefix oatk-t-h-XY
#   h-t : HiMT -> Tiara -> Oatk        -> prefix oatk-h-t-XY
#   t-x : Tiara -> Oatk                -> prefix oatk-t-x-XY
#   x-h : HiMT  -> Oatk                -> prefix oatk-x-h-XY
#   x-x : Oatk only                    -> prefix oatk-x-x-XY
#
# Usage:
#   run-oatk-tiara-himt_genus_species OUTDIR [INUM] [PATTERN] [C_LIST]
#     OUTDIR  : species folder (e.g. Erythranthe_laciniata)
#     INUM    : integer index (default 0)
#     PATTERN : t-h | h-t | t-x | x-h | x-x  (default t-h)
#     C_LIST  : comma-separated list of c-values, 2-digit padded (01..48)
#               default: 01,02,...,48
################################################################################
run-oatk-tiara-himt_genus_species() {
	local bolap_cmd="${FUNCNAME%%_*}"

	help_message=$(
		cat <<EOF
Name:
  bolap ${bolap_cmd}

Synopsis:
  bolap ${bolap_cmd}

Description:
  Run Oatk on HiFi reads with Tiara + HiMT filtering in configurable order

Options:
  -s STR
    species folder (e.g. Erythranthe_laciniata)

  -i INT
    integer index (default 0)

  --pattern
    PATTERN : t-h | h-t | t-x | x-h | x-x  (default t-h)

  --c-list
    C_LIST  : comma-separated list of c-values, 2-digit padded (01..48)
              default: 01,02,...,48

Examples:
  Subtitle:
    bolap ${bolap_cmd} -s species_name -i 0 --pattern t-h --c-list 01,02,03

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	parse_commandline() {
		# echo "${_brg_args[@]}" >&2
		set -- "${_brg_args[@]}"

		while test $# -gt 0; do
			_key="$1"
			case "$_key" in
			--pattern)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_pattern="$2"
					shift || true
				fi
				;;
			--c-list)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_c_list="$2"
					shift || true
				fi
				;;
			esac
			shift || true
		done
	}

	_polap_lib_help-maybe-show "$bolap_cmd" help_message || return 0

	# local _brg_outdir="${1:-${_brg_outdir}}"
	# local _brg_sindex="${2:-${_brg_sindex:-0}}"
	local _pattern="t-h"
	local _c_list=""
	# local _c_list="01..48"

	bolap_s_parse_commandline
	parse_commandline

	get_padded_sequence_pure() {
		# Extract start and end numbers using parameter expansion
		local start=${1%%..*}
		local end=${1##*..}
		local output=""

		# Loop from start to end
		for ((i = start; i <= end; i++)); do
			# Format the number with 0 padding
			printf -v num "%02d" $i

			# Append to output with comma
			output="${output}${num},"
		done

		# Print output with the trailing comma removed
		echo "${output%,}"
	}

	# Default c_list: 01..48
	if [[ -z "$_c_list" ]]; then
		local i
		for i in $(seq 2 49); do
			printf -v c "%02d" "$i"
			if [[ -z "$_c_list" ]]; then
				_c_list="$c"
			else
				_c_list="${_c_list},$c"
			fi
		done
	elif [[ "$_c_list" == *..* ]]; then
		# Usage Example
		_c_list=$(get_padded_sequence_pure "$_c_list")
	fi

	# Normalize pattern
	case "$_pattern" in
	t-h | h-t | t-x | x-h | x-x) ;;
	*)
		echo "[ERR] run-oatk-tiara-himt_genus_species: unknown pattern: $_pattern" >&2
		return 2
		;;
	esac

	# Boiler-plate: defines _brg_outdir_i etc
	local _brg_title="oatk-tiara-himt-${_pattern}"
	source "${_POLAPLIB_DIR}/polap-variables-data.sh"
	local _brg_rundir="${_brg_outdir_i}/${_brg_title}"

	mkdir -p "${_brg_outdir_i}"

	# Input reads (HiFi)
	local reads="${long_sra}.fastq"
	if [[ ! -s "$reads" ]]; then
		echo "[ERR] run-oatk-tiara-himt_genus_species: no reads found for ${_brg_outdir}/${_brg_adir}/${_brg_sindex}" >&2
		return 2
	fi

	# Working folder for this pipeline
	# local work="${_brg_target}/oatk-tiara-himt"
	local work="${_brg_target}"
	mkdir -p "$work"
	local tiara_reads=""
	local himt_reads=""
	local oatk_reads=""

	# ------------------------------------------------------------------#
	# Decide the sequence:
	#   t-h : Tiara -> HiMT
	#   h-t : HiMT  -> Tiara
	#   t-x : Tiara only
	#   x-h : HiMT  only
	#   x-x : no filtering (raw reads)
	# ------------------------------------------------------------------#
	case "$_pattern" in
	t-h)
		_polap_oatk_step_tiara "$reads" "$work" tiara_reads || return 1
		# tiara_reads="$work/tiara/tiara.mito.fq.gz"
		_polap_oatk_step_himt "$tiara_reads" "$work" himt_reads || return 1
		oatk_reads="$himt_reads"
		;;

	h-t)
		_polap_oatk_step_himt "$reads" "$work" himt_reads || return 1
		_polap_oatk_step_tiara "$himt_reads" "$work" tiara_reads || return 1
		oatk_reads="$tiara_reads"
		;;

	t-x)
		_polap_oatk_step_tiara "$reads" "$work" tiara_reads || return 1
		oatk_reads="$tiara_reads"
		;;

	x-h)
		_polap_oatk_step_himt "$reads" "$work" himt_reads || return 1
		oatk_reads="$himt_reads"
		;;

	x-x)
		oatk_reads="$reads"
		;;

	*)
		echo "[ERR] unknown pattern '$_pattern' (expected one of t-h h-t t-x x-h x-x)" >&2
		return 2
		;;
	esac

	# ------------------------------------------------------------------#
	# From here on, use $oatk_reads as the input to Oatk
	# ------------------------------------------------------------------#
	echo "[INFO] Oatk will use reads: $oatk_reads"

	# -------------------------------------------------------------------------- #
	# Step 3: Oatk for each c in _c_list
	# Output prefix: oatk-<pattern>-<c>
	# utg GFA/FASTA etc should be written into oatk/ so that downstream tools see them
	# -------------------------------------------------------------------------- #
	_polap_lib_conda-ensure_conda_env polap-oatk || return 1

	local oatk_dir="${work}/oatk"
	mkdir -p "$oatk_dir"

	local c
	IFS=',' read -r -a _c_arr <<<"$_c_list"
	for c in "${_c_arr[@]}"; do
		local prefix="oatk-${_pattern}-${c}"
		local outprefix="${oatk_dir}/${prefix}"

		c_num=$((10#$c))
		local _stdout_txt="${_stdout_txt_base}-c${c}.txt"
		local _timing_txt="${_timing_txt_base}-c${c}.txt"

		echo "[INFO] Running Oatk ($prefix) on $oatk_reads (c=$c_num)"

		local rc=0
		local _oatk_threads=8
		command time -v oatk \
			-k 1001 \
			-c ${c_num} \
			-t ${_oatk_threads} \
			-m ./OatkDB/v20230921/embryophyta_mito.fam \
			-p ./OatkDB/v20230921/embryophyta_pltd.fam \
			-o "${outprefix}" \
			"${oatk_reads}" \
			>"${_stdout_txt}" \
			2>"${_timing_txt}" || rc=$?

		if [[ "${rc:-0}" -ne 0 ]]; then
			printf '[ERR] oatk failed (rc=%d)\n' "${rc:-0}" >&2
		fi

	done

	_polap_lib_conda-ensure_conda_env_deactivate || true

	# Save results
	# rsync -azuq --max-size=5M \
	# 	"${_brg_target}/" "${_brg_rundir}/"
	rsync -azuq "${_brg_target}/" "${_brg_rundir}/"

	# Clean-up
	rm -rf "${_brg_target}"
}

################################################################################
# Run HiMT "assess" on all oatk-* patterns and c values
# Produces 5 * 48 assess folders:
#   <outdir>/v5/0/assess/oatk-t-h-01, ..., oatk-x-x-48
################################################################################
run-assess-oatk-tiara-himt_genus_species() {
	local bolap_cmd="${FUNCNAME%%_*}"

	help_message=$(
		cat <<EOF
Name:
  bolap ${bolap_cmd}

Synopsis:
  bolap ${bolap_cmd}

Description:
  text.

Examples:
  Subtitle:
    bolap ${bolap_cmd} -s Genus_species

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	_polap_lib_help-maybe-show "$bolap_cmd" help_message || return 0

	bolap_s_parse_commandline

	source "${_POLAPLIB_DIR}/polap-variables-data.sh"

	local assess_root="${_brg_outdir_i}/assess"
	mkdir -p "$assess_root"

	_polap_lib_conda-ensure_conda_env polap-himt || return 1

	rm -f "${_stdout_txt}"
	rm -f "${_timing_txt}"

	local pattern c
	for pattern in t-h h-t t-x x-h x-x; do
		for i in $(seq 2 49); do
			printf -v c "%02d" "$i"
			local oatk_dir="${_brg_outdir_i}/oatk-tiara-himt-${pattern}"
			local prefix="oatk-${pattern}-${c}"
			local gfa="${oatk_dir}/oatk/${prefix}.mito.gfa"
			local png="${oatk_dir}/oatk/${prefix}.mito.png"
			local mito_fa="${oatk_dir}/oatk/${prefix}.mito.ctg.fasta"
			local outdir="${assess_root}/${prefix}"

			if [[ ! -s "$mito_fa" ]]; then
				# echo "[WARN] Missing mito assembly for $prefix ($mito_fa); skipping assess" >&2
				continue
			fi

			rm -rf "$outdir"
			mkdir -p "$outdir"
			echo "[INFO] HiMT assess for $prefix"
			if [[ -s "${mito_fa}" ]]; then

				himt assess \
					--input_file "${mito_fa}" \
					--output_dir "$outdir" \
					>>"${_stdout_txt}" \
					2>>"${_timing_txt}" || rc=$?
			else
				echo "[ERROR] no fasta: ${mito_fa}" >&2
			fi

			if [[ "${rc:-0}" -ne 0 ]]; then
				printf "[ERR] himt assess on ${prefix} failed (rc=%d)\n" "${rc:-0}" >&2
			fi

			# convert HiMT assess html to report
			local assess_himt_html="${outdir}/himt_mitochondrial.html"
			local assess_report_prefix="${outdir}/report/prefix"
			if [[ -s "${assess_himt_html}" ]]; then
				${_polap_cmd} convert himt2prefix \
					"${assess_himt_html}" "${assess_report_prefix}"
			else
				echo "[ERROR] no himt html: ${assess_himt_html}" >&2
			fi

			if [[ -s "${gfa}" ]]; then
				# echo ${_polap_cmd} bandage png "${gfa}" "${png}" >&2
				${_polap_cmd} bandage png "${gfa}" "${png}"
			else
				echo "[ERROR] no gfa: ${gfa}" >&2
			fi
		done
	done

	_polap_lib_conda-ensure_conda_env_deactivate || true
}

# 2025-11-30
run-himt_genus_species() {
	local _brg_outdir="${1:-$_brg_outdir}"
	local _brg_sindex="${2:-$_brg_sindex}"

	source "${_POLAPLIB_DIR}/polap-variables-data.sh"

	# Create folders
	rm -rf "${_brg_target}"
	mkdir -p "${_brg_target}"
	mkdir -p "${_brg_outdir_i}"

	# ln -sf "${_brg_tmpdir}/l.fq" "${long_sra}.fastq"
	_brg_input_data="${long_sra}.fastq"

	# check the input
	# prepare if there is no such one

	_polap_lib_conda-ensure_conda_env polap-himt || exit 1

	# local _fc=30
	# local formatted_fc=$(printf "%02d" "${_fc}")

	_log_echo "himt on ${_brg_rundir}"
	# local formatted_fc=$(printf "%02d" "${_fc}")

	rm -f "${_stdout_txt}"
	rm -f "${_timing_txt}"

	# Start memory logger
	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_oatk_memtracker_time_interval}"

	local rc
	local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"
	local _himt_threads="${_brg_threads}"

	command time -v timeout ${_brg_timeout} himt \
		assemble \
		-t ${_himt_threads} \
		-n 4 \
		-b 3 \
		-s plant \
		--input_file "${_brg_input_data}" \
		--output_dir "${_brg_target}" \
		>"${_stdout_txt}" \
		2>"${_timing_txt}" || rc=$?

	if [[ "${rc:-0}" -ne 0 ]]; then
		printf '[ERR] oatk failed (rc=%d)\n' "$rc" >&2
		# decide: exit "$rc" / continue
	fi

	# Save system info
	_polap_lib_timing-get_system_info >>"${_timing_txt}"

	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"

	conda deactivate

	# Save results
	rsync -azuq --max-size=5M \
		"${_brg_target}/" "${_brg_rundir}/"

	# Clean-up
	rm -rf "${_brg_target}"
}

# 2025-12-01
# HiMT assess -> report for all 4 pipelines
# gfa -> png using Bandage
run-assess-himt_genus_species() {
	local _brg_outdir="${1:-$_brg_outdir}"
	local _brg_sindex="${2:-$_brg_sindex}"

	source "${_POLAPLIB_DIR}/polap-variables-data.sh"

	# Create folders
	rm -rf "${_brg_target}"
	mkdir -p "${_brg_target}"
	mkdir -p "${_brg_outdir_i}"

	# ln -sf "${_brg_tmpdir}/l.fq" "${long_sra}.fastq"
	_brg_input_data="${long_sra}.fastq"

	# check the input
	# prepare if there is no such one

	_polap_lib_conda-ensure_conda_env polap-himt || exit 1

	# local _fc=30
	# local formatted_fc=$(printf "%02d" "${_fc}")

	_log_echo "assess himt on ${_brg_rundir}"
	# local formatted_fc=$(printf "%02d" "${_fc}")

	rm -f "${_stdout_txt}"
	rm -f "${_timing_txt}"

	# Start memory logger
	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_oatk_memtracker_time_interval}"

	local rc
	local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"
	local _himt_threads="${_brg_threads}"

	# assess oatk
	local _fc
	# for _fc in 30 25 20 15 10 5; do
	for _fc in 30 20 10; do
		local formatted_fc=$(printf "%02d" "${_fc}")
		# local gfa="${_brg_outdir_i}/oatk/oatk-${formatted_fc}.mito.gfa"
		local gfa="${_brg_outdir_i}/oatk/oatk-${formatted_fc}.mito.gfa"
		local png="${_brg_outdir_i}/oatk/oatk-${formatted_fc}.mito.png"
		local fasta="${_brg_outdir_i}/oatk/oatk-${formatted_fc}.mito.ctg.fasta"
		local assess_gfa="${_brg_outdir_i}/assess/oatk-${formatted_fc}"
		himt \
			assess \
			--input_file "${fasta}" \
			--output_dir "${assess_gfa}" \
			>>"${_stdout_txt}" \
			2>>"${_timing_txt}" || rc=$?

		if [[ "${rc:-0}" -ne 0 ]]; then
			printf "[ERR] himt assess on oatk-${formatted_fc} failed (rc=%d)\n" "$rc" >&2
		fi

		# convert HiMT assess html to report
		local assess_himt_html="${_brg_outdir_i}/assess/oatk-${formatted_fc}/himt_mitochondrial.html"
		local assess_report_prefix="${_brg_outdir_i}/assess/oatk-${formatted_fc}/report/prefix"
		${_polap_cmd} convert himt2prefix \
			"${assess_himt_html}" "${assess_report_prefix}"

		if [[ -s "${gfa}" ]]; then
			${_polap_cmd} bandage png "${gfa}" "${png}"
		fi

	done

	# assess himt
	local gfa="${_brg_outdir_i}/himt/himt_mitochondrial.gfa"
	local png="${_brg_outdir_i}/himt/himt_mitochondrial.png"
	local assess_gfa="${_brg_outdir_i}/assess/himt"
	himt \
		assess \
		--input_file "${gfa}" \
		--output_dir "${assess_gfa}" \
		>"${_stdout_txt}" \
		2>"${_timing_txt}" || rc=$?

	if [[ "${rc:-0}" -ne 0 ]]; then
		printf '[ERR] himt assess on himt failed (rc=%d)\n' "$rc" >&2
	fi

	# convert HiMT assess html to report
	local assess_himt_html="${_brg_outdir_i}/assess/himt/himt_mitochondrial.html"
	local assess_report_prefix="${_brg_outdir_i}/assess/himt/report/prefix"
	${_polap_cmd} convert himt2prefix \
		"${assess_himt_html}" "${assess_report_prefix}"

	if [[ -s "${gfa}" ]]; then
		${_polap_cmd} bandage png "${gfa}" "${png}"
	fi

	# assess tippo
	local gfa="${_brg_outdir_i}/tippo/${long_sra}.fastq.mitochondrial.fasta.filter.fasta.flye/assembly_graph.gfa"
	local png="${_brg_outdir_i}/tippo/${long_sra}.fastq.mitochondrial.fasta.filter.fasta.flye/assembly_graph.png"
	local assess_gfa="${_brg_outdir_i}/assess/tippo"
	himt \
		assess \
		--input_file "${gfa}" \
		--output_dir "${assess_gfa}" \
		>>"${_stdout_txt}" \
		2>>"${_timing_txt}" || rc=$?

	if [[ "${rc:-0}" -ne 0 ]]; then
		printf '[ERR] himt assess on tippo failed (rc=%d)\n' "$rc" >&2
	fi

	# convert HiMT assess html to report
	local assess_himt_html="${_brg_outdir_i}/assess/tippo/himt_mitochondrial.html"
	local assess_report_prefix="${_brg_outdir_i}/assess/tippo/report/prefix"
	${_polap_cmd} convert himt2prefix \
		"${assess_himt_html}" "${assess_report_prefix}"

	if [[ -s "${gfa}" ]]; then
		${_polap_cmd} bandage png "${gfa}" "${png}"
	fi

	local gfa="${_brg_outdir_i}/pmat2/gfa_result/PMAT_mt_main.gfa"
	local png="${_brg_outdir_i}/pmat2/gfa_result/PMAT_mt_main.png"
	local assess_gfa="${_brg_outdir_i}/assess/pmat2"
	himt \
		assess \
		--input_file "${gfa}" \
		--output_dir "${assess_gfa}" \
		>>"${_stdout_txt}" \
		2>>"${_timing_txt}" || rc=$?

	if [[ "${rc:-0}" -ne 0 ]]; then
		printf '[ERR] himt assess on pmat2 failed (rc=%d)\n' "$rc" >&2
	fi

	# convert HiMT assess html to report
	local assess_himt_html="${_brg_outdir_i}/assess/pmat2/himt_mitochondrial.html"
	local assess_report_prefix="${_brg_outdir_i}/assess/pmat2/report/prefix"
	${_polap_cmd} convert himt2prefix \
		"${assess_himt_html}" "${assess_report_prefix}"

	if [[ -s "${gfa}" ]]; then
		${_polap_cmd} bandage png "${gfa}" "${png}"
	fi

	# Save system info
	_polap_lib_timing-get_system_info >>"${_timing_txt}"

	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"

	conda deactivate

	# Clean-up
	rm -rf "${_brg_target}"
}

# 2025-11-30
run-pmat2_genus_species() {
	local _brg_outdir="${1:-$_brg_outdir}"
	local _brg_sindex="${2:-$_brg_sindex}"

	source "${_POLAPLIB_DIR}/polap-variables-data.sh"

	# Create folders
	rm -rf "${_brg_target}"
	mkdir -p "${_brg_target}"
	mkdir -p "/tmp/${_brg_target}"
	mkdir -p "${_brg_outdir_i}"

	# ln -sf "${_brg_tmpdir}/l.fq" "${long_sra}.fastq"
	_brg_input_data="${long_sra}.fastq"

	# check the input
	# prepare if there is no such one

	_polap_lib_conda-ensure_conda_env polap-pmat2 || exit 1

	# local _fc=30
	# local formatted_fc=$(printf "%02d" "${_fc}")

	_log_echo "pmat2 on ${_brg_rundir}"
	# local formatted_fc=$(printf "%02d" "${_fc}")

	rm -f "${_stdout_txt}"
	rm -f "${_timing_txt}"

	# Start memory logger
	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_oatk_memtracker_time_interval}"

	local rc
	local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"
	local _pmat2_threads="${_brg_threads}"

	# pmat2 options
	# -k 31 \
	# -p 1 \
	# -G mt \
	# -x 0 \
	# -F 1 \
	local pmat_cmd=$(which PMAT)
	local pmat_out=$(realpath "${_brg_target}")
	command time -v timeout ${_brg_timeout} $pmat_cmd autoMito \
		-i "${_brg_input_data}" \
		-o "/tmp/${_brg_target}" \
		-t hifi \
		-m \
		-T ${_brg_threads} \
		>"${_stdout_txt}" \
		2>"${_timing_txt}" || rc=$?

	if [[ "${rc:-0}" -ne 0 ]]; then
		printf '[ERR] pmat2 failed (rc=%d)\n' "$rc" >&2
		# decide: exit "$rc" / continue
	fi

	# Save system info
	_polap_lib_timing-get_system_info >>"${_timing_txt}"

	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"

	set +u
	conda deactivate
	set -u

	# Save results
	rsync -azuq --max-size=5M \
		"/tmp/${_brg_target}/" "${_brg_rundir}/"

	# Clean-up
	rm -rf "${_brg_target}"
	rm -rf "/tmp/${_brg_target}"
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
	local _brg_outdir="${1:-$_brg_outdir}"
	local _brg_sindex="${2:-$_brg_sindex}"
	local _brg_type="${3:-pt}"

	local _brg_title="${_brg_type}gaul"
	source "${_POLAPLIB_DIR}/polap-variables-data.sh"

	# Always redo
	rm -f "${_timing_txt}"
	rm -f "${_stdout_txt}"

	# Always redo
	rm -rf "${_brg_target}"
	mkdir -p "${_brg_target}"

	# Initialize Conda
	_polap_lib_conda-ensure_conda_env polap || exit 1

	# Start memory logger
	# GetOrganelle takes not much time to finish, and we use 10 seconds of interval for memlog
	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"

	local genome_size=$(seqkit stats -T "${_brg_outdir_i}/ncbi-${_brg_type}dna/${_brg_type}dna-reference.fa" | awk 'NR==2 {print $7/1}')

	command time -v bash ${_POLAPLIB_DIR}/polap-ptGAUL1.sh \
		-o ${_brg_target} \
		-r "${_brg_outdir_i}/ncbi-${_brg_type}dna/${_brg_type}dna-reference.fa" \
		-g "${genome_size}" \
		-l ${_brg_tmpdir}/l.fq.gz \
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
		"${_brg_target}"/result_3000/ \
		"${_brg_rundir}"/

	rm -rf "${_brg_target}"
}

run-mtgaul_genus_species() {
	local _brg_outdir="${1:-$_brg_outdir}"
	local _brg_sindex="${2:-$_brg_sindex}"
	run-ptgaul_genus_species "${_brg_outdir}" "${_brg_sindex}" "mt"
}

run-download-mtdna_genus_species() {
	download-mtdna_genus_species "$@"
}

download-mtdna_genus_species() {
	local _brg_outdir="${1:-$_brg_outdir}"
	local _brg_sindex="${2:-$_brg_sindex}"
	local _brg_type="${3:-mtdna}"

	local _brg_title="ncbi-${_brg_type}"
	source "${_POLAPLIB_DIR}/polap-variables-data.sh"

	local _species=${_taxon[$_brg_target]}
	local species_name="${_species//_/ }"

	mkdir -p "${_brg_rundir}"

	# Initialize Conda
	_polap_lib_conda-ensure_conda_env polap || exit 1

	if [[ -s "${_brg_rundir}/${_brg_type}-reference.fa" ]]; then
		echo "found: ${_brg_type} reference: ${_brg_outdir_i}/${_brg_type}-reference.fa"
	else

		${_polap_cmd} get-mtdna \
			$([[ "${_brg_type}" == "ptdna" ]] && echo "--plastid") \
			--species "${species_name}" \
			-o "${_brg_rundir}"

		if [[ -s "${_brg_rundir}/00-bioproject/2-mtdna.fasta" ]]; then
			echo "${_brg_type} of ${species_name}: ${_brg_rundir}/${_brg_type}-reference.fa"
			cp -p "${_brg_rundir}/00-bioproject/2-mtdna.fasta" \
				"${_brg_rundir}/${_brg_type}-reference.fa"
		else
			local _genus_name=$(echo ${species_name} | awk '{print $1}')
			echo "No ${_brg_type} of the species name: ${species_name}"
			echo "  trying to search NCBI plastid genomes for genus name only: ${_genus_name}"

			${_polap_cmd} get-mtdna \
				$([[ "${_brg_type}" == "ptdna" ]] && echo "--plastid") \
				--species "${species_name}" \
				-o "${_brg_rundir}"

			if [[ -s "${_brg_rundir}/00-bioproject/2-mtdna.fasta" ]]; then
				echo "${_brg_type} of ${_genus_name}: ${_brg_rundir}/${_brg_type}-reference.fa"
				cp -p "${_brg_rundir}/00-bioproject/2-mtdna.fasta" \
					"${_brg_rundir}/${_brg_type}-reference.fa"
			else
				echo "  we could not find one even in the genus level."
				echo "No such file: ${_brg_rundir}/${_brg_type}-reference.fa"
			fi
		fi
	fi

	conda deactivate
}

download-ptdna_genus_species() {
	local _brg_outdir="${1:-$_brg_outdir}"
	local _brg_sindex="${2:-$_brg_sindex}"
	download-mtdna_genus_species "${_brg_outdir}" "${_brg_sindex}" "ptdna"
}

run-download-ptdna_genus_species() {
	download-ptdna_genus_species "$@"
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
	benchmark_genus_species
	# local first_arg="$1"
	# local remaining_args=("${@:2}")

	# # Remove trailing slash from the first element
	# if [[ -n "${remaining_args[0]}" ]]; then
	# 	remaining_args[0]="${remaining_args[0]%/}"
	# fi

	# if declare -F run-${first_arg}_genus_species >/dev/null; then
	# 	run-${first_arg}_genus_species "${remaining_args[@]}"
	# else
	# 	echo "$help_message_run"
	# fi
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

install-pmat2_genus_species() {
	local want_env="polap-pmat2"
	local -a pkgs=(
		zlib gcc_linux-64 gxx_linux-64
		apptainer nextdenovo canu blast
	)

	# Confirm (honors opt_y_flag if you set it elsewhere)
	local confirm
	if [[ "${opt_y_flag-}" == "true" ]]; then
		confirm="yes"
	else
		read -r -p "Do you want to install pmat2 in the ${want_env} conda environment? (y/N): " confirm
	fi

	if [[ "${confirm,,}" == "y" || "${confirm,,}" == "yes" ]]; then
		# Create/upgrade env (channels first for R/bioconda harmony)
		set +u
		_polap_lib_conda-create-env \
			"$want_env" "${pkgs[@]}" \
			--channel conda-forge --channel bioconda -y || return 1

		# Quick sanity check that key tools are on PATH
		local t
		for t in blastn; do
			if ! command -v "$t" >/dev/null 2>&1; then
				echo "WARNING: '$t' not found in '$want_env' PATH." >&2
			fi
		done
		set -u

		# wget https://github.com/aiPGAB/PMAT2/archive/refs/tags/v2.1.5.tar.gz
		# tar -zxvf v2.1.5.tar.gz
		# cd PMAT2-2.1.5
		# make
		# ./PMAT --help

	else
		echo "pmat2 installation is canceled."
		echo "conda install -c conda-forge zlib gcc_linux-64 gxx_linux-64"
		echo "Check: https://github.com/aiPGAB/PMAT2"
		echo "wget https://github.com/aiPGAB/PMAT2/archive/refs/tags/v2.1.5.tar.gz"
		echo "tar -zxvf PMAT2-2.1.5.tar.gz"
		echo "cd PMAT2-2.1.5"
		echo "make"
		echo "./PMAT --help"
	fi
}

setup-pmat2_genus_species() {

	set +u
	_polap_lib_conda-ensure_conda_env polap-pmat2 || exit 1
	set -u

	if [[ "$CONDA_DEFAULT_ENV" == "polap-pmat2" ]]; then
		rm -f v2.1.5.tar.gz
		wget https://github.com/aiPGAB/PMAT2/archive/refs/tags/v2.1.5.tar.gz
		rm -rf PMAT2-2.1.5
		tar -zxvf v2.1.5.tar.gz
		cd PMAT2-2.1.5

		make -f "${_POLAPLIB_DIR}/Makefile.pmat2" \
			CFLAGS="-g -O3 -Wall -Wno-unused-function -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-result -I${CONDA_PREFIX}/include" \
			LDFLAGS="-L${CONDA_PREFIX}/lib"

		install -Dm755 PMAT "$CONDA_PREFIX/bin/PMAT"
		cp -r Conserved_PCGs_db "$CONDA_PREFIX/bin/"
		cp -r container "$CONDA_PREFIX/bin/"
		cp -r lib "$CONDA_PREFIX/bin/"
	else
		echo "Error: You're in the '$CONDA_DEFAULT_ENV' environment. Please activate base before running this script."
		exit 1
	fi
}

test-pmat2_genus_species() {
	_polap_lib_conda-ensure_conda_env polap-pmat2 || exit 1
	if [[ "$CONDA_DEFAULT_ENV" == "polap-pmat2" ]]; then
		$CONDA_PREFIX/bin/PMAT --help
	else
		echo "Error: You're in the '$CONDA_DEFAULT_ENV' environment. Please activate base before running this script."
		exit 1
	fi
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
				echo "[INFO] Conda environment 'polap-pmat' already exists."
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

install-himt_genus_species() {
	local want_env="polap-himt"
	local -a pkgs=(
		shuyuan_tang::himt
	)

	# Confirm (honors opt_y_flag if you set it elsewhere)
	local confirm
	if [[ "${opt_y_flag-}" == "true" ]]; then
		confirm="yes"
	else
		read -r -p "Do you want to install himt in the ${want_env} conda environment? (y/N): " confirm
	fi

	if [[ "${confirm,,}" == "y" || "${confirm,,}" == "yes" ]]; then
		# Create/upgrade env (channels first for R/bioconda harmony)
		_polap_lib_conda-create-env \
			"$want_env" "${pkgs[@]}" \
			--channel conda-forge --channel bioconda -y || return 1

		# Quick sanity check that key tools are on PATH
		local t
		for t in python; do
			if ! command -v "$t" >/dev/null 2>&1; then
				echo "WARNING: '$t' not found in '$want_env' PATH." >&2
			fi
		done
	else
		echo "himt installation is canceled."
		echo "https://anaconda.org/channels/shuyuan_tang/packages/himt/overview"
		echo "Check: https://github.com/tang-shuyuan/HiMT"
		echo "conda install shuyuan_tang::himt -c bioconda"
		echo "Execute: himt"
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
				echo "[INFO] Conda environment 'polap-dflye' already exists."
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
				echo "[INFO] Conda environment 'polap-cflye' already exists."
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

x-install-polap_genus_species() {
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
				echo "[INFO] Conda environment 'polap' already exists."
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

polap_setup_parse_commandline() {
	set -- "${_brg_args[@]}"

	# source "${_POLAPLIB_DIR}/polap-cmd-version.sh" # '.' means 'source'
	while test $# -gt 0; do
		_key="$1"
		case "$_key" in
		--ver)
			test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
			_brg_ver="$2"
			shift || true
			;;
		esac
		shift || true
	done
}

# 2025-12-03
install-polap_genus_species() {
	local bolap_cmd="${FUNCNAME%%_*}"

	help_message=$(
		cat <<EOF
Name:
  bolap - $bolap_cmd

Synopsis:
  bolap $bolap_cmd

Description:
  Install polap in polap conda environment.

Examples:
  Install polap:
    bolap install polap

  Install polap of version 0.4:
    bolap install polap --ver 0.4

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	# defaults
	local _brg_ver=""

	_polap_lib_help-maybe-show "$bolap_cmd" help_message || return 0

	polap_setup_parse_commandline

	local _version="${_brg_ver}"

	local want_env="polap"

	if [[ -n "${_version}" ]]; then
		local -a pkgs=(
			polap="${_version}"
		)
	else
		local -a pkgs=(
			polap
		)
	fi

	# Confirm (honors opt_y_flag if you set it elsewhere)
	local confirm
	if [[ "${opt_y_flag-}" == "true" ]]; then
		confirm="yes"
	else
		if [[ -n "${_version}" ]]; then
			read -r -p "Do you want to install polap v${_version} in ${want_env} conda environment? (y/N): " confirm
		else
			read -r -p "Do you want to install polap in ${want_env} conda environment? (y/N): " confirm
		fi
	fi

	if [[ "${confirm,,}" == "y" || "${confirm,,}" == "yes" ]]; then
		# Create/upgrade env (channels first for R/bioconda harmony)
		if [[ "${opt_y_flag-}" == "true" ]]; then
			_polap_lib_conda-create-env \
				"$want_env" "${pkgs[@]}" \
				--channel conda-forge --channel bioconda -y || return 1
		else
			_polap_lib_conda-create-env \
				"$want_env" "${pkgs[@]}" \
				--channel conda-forge --channel bioconda || return 1
		fi

		# Quick sanity check that key tools are on PATH
		local t
		for t in polap; do
			if ! command -v "$t" >/dev/null 2>&1; then
				echo "WARNING: '$t' not found in '$want_env' PATH." >&2
			fi
		done
	else
		echo "polap installation is canceled."
		echo "Check: https://github.com/goshng/polap"
		echo "conda install -c bioconda polap"
		echo "Execute: polap"
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
				echo "[INFO] Conda environment 'polap-fmlrc' already exists."
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
	if [[ "${opt_y_flag-}" == false ]]; then
		read -p "Do you want to setup miniconda3 for conda-forge and bioconda? (y/N): " confirm
	else
		confirm="yes"
	fi

	if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
		# Check if conda is available
		_polap_lib_conda-ensure_conda_env base || exit 1

		if [[ "$CONDA_DEFAULT_ENV" == "base" ]]; then
			if command -v conda &>/dev/null; then
				source "$(conda info --base)/etc/profile.d/conda.sh"
				conda activate base
				conda config --add channels bioconda
				conda config --add channels conda-forge
				conda config --set channel_priority strict

				# Terms of Service have not been accepted for the following channels.
				conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
				conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
			else
				echo "Error: conda command not found."
				exit 1
			fi
		else
			echo "Error: You're in the '$CONDA_DEFAULT_ENV' environment. Please activate base before running this script."
			exit 1
		fi
	else
		echo "miniconda3 bioconda setup is canceled."
	fi
}

show-conda_genus_species() {
	# Check if conda is available
	_polap_lib_conda-ensure_conda_env base || exit 1

	if [[ "$CONDA_DEFAULT_ENV" == "base" ]]; then
		if command -v conda &>/dev/null; then
			conda config --show channels
		else
			echo "Error: conda command not found."
			exit 1
		fi
	else
		echo "Error: You're in the '$CONDA_DEFAULT_ENV' environment. Please activate base before running this script."
		exit 1
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

	if declare -F setup-${first_arg}_genus_species >/dev/null; then
		setup-${first_arg}_genus_species "${remaining_args[@]:-}"
	else
		echo "$help_message_setup"
	fi
}

# list_genus_species() {
# 	local _brg_query="${1}"
# 	local _brg_where="${2:-any}"
#
# 	if [[ "${_brg_where}" == "any" ]]; then
# 		grep ")$" "${_POLAPLIB_DIR}/polap-lib-data.sh" |
# 			grep -v "(" | grep "${_brg_query}" | clean_input_lines | sort | uniq
# 		grep ")$" "${_POLAPLIB_DIR}/polap-lib-data.sh" |
# 			grep -v "(" | grep "${_brg_query}" | clean_input_lines | sort | uniq
# 		grep ")$" $0 | grep -v "(" | grep "${_brg_query}" | clean_input_lines | sort | uniq
# 	elif [[ "${_brg_where}" == "end" ]]; then
# 		grep ")$" "${_POLAPLIB_DIR}/polap-lib-data.sh" |
# 			grep -v "(" | grep "${_brg_query})$" | clean_input_lines | sort | uniq
# 		grep ")$" $0 | grep -v "(" | grep "${_brg_query})$" | clean_input_lines | sort | uniq
# 	else
# 		grep ")$" "${_POLAPLIB_DIR}/polap-lib-data.sh" |
# 			clean_input_lines | grep -v "(" | grep "^${_brg_query}" | sort | uniq
# 		grep ")$" $0 | clean_input_lines | grep -v "(" | grep "^${_brg_query}" | sort | uniq
# 	fi
# }

# ────────────────────────────────────────────────────────────────
# clean_input_lines
#   - trims whitespace
#   - collapses multiple spaces/tabs into one TAB
#   - removes empty lines
# Usage:  some_command | clean_input_lines
# ────────────────────────────────────────────────────────────────
clean_input_lines() {
	awk '
	{
		# trim leading/trailing whitespace
		gsub(/^[ \t]+/, "", $0)
		gsub(/[ \t]+$/, "", $0)
		# collapse runs of spaces/tabs to a single TAB
		gsub(/[ \t]+/, "\t", $0)
		# print only if line is not empty
		if (length($0) > 0) print
	}'
}

# ------------------------------------------------------------------------------
# list_genus_species
# - Finds genus/species tokens from two places:
#   (1) data lines that end with ')' (your original behavior)
#   (2) function headers of the form  XXX_genus_species() {   and extracts 'XXX'
# Usage:
#   list_genus_species "<query>" [any|end|start]
# ------------------------------------------------------------------------------
list_genus_species() {

	# Where to look
	local _files=(
		"${_POLAPLIB_DIR}/polap-lib-data.sh"
		"${_POLAPLIB_DIR}/polap-data-read.sh"
		"$0"
	)

	local bolap_cmd="${FUNCNAME%%_*}"

	help_message=$(
		cat <<EOF
Name:
  bolap - search commands

Synopsis:
  bolap ${bolap_cmd} --query STR --where STR

Description:
  --query STR

  --where any|end|start

Examples:
  Search command:
    bolap ${bolap_cmd} --query data --where any

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	local _brg_query=".*"
	local _brg_where="any"

	parse_commandline() {
		set -- "${_brg_unknown_opts[@]}"

		# source "${_POLAPLIB_DIR}/polap-cmd-version.sh" # '.' means 'source'
		while test $# -gt 0; do
			_key="$1"
			case "$_key" in
			--query)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_query="$2"
					shift || true
				fi
				;;
			--where)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_where="$2"
					shift || true
				fi
				;;
			-*)
				_log_echo0 "[ERROR] no such options: $1"
				;;
			*)
				break
				;;
			esac
			shift || true
		done
	}

	declare -n ref="help_message"
	if [[ "${_brg_help}" == "on" ]]; then
		local manfile=$(_bolap_lib_man-convert_help_message "$ref" "${bolap_cmd}")
		man "$manfile"
		rm -f "$manfile"
		return
	fi

	parse_commandline

	# Build a regex according to where-to-match semantics
	local _re
	case "${_brg_where}" in
	any) _re="${_brg_query}" ;;
	end) _re="${_brg_query}\$" ;;
	start) _re="^${_brg_query}" ;;
	*) _re="${_brg_query}" ;;
	esac

	# If you have a cleaner, use it; otherwise pass-through.
	local _cleaner="cat"
	if declare -F clean_input_lines >/dev/null 2>&1; then
		_cleaner="clean_input_lines"
	fi

	# Collect candidates to a single stream, then uniq/sort once.
	{
		# ---------- Source A: original “data lines” ----------
		# lines that end with ')', but do not contain '('
		# (Your pipelines run twice on polap-lib-data.sh; removed the duplicate.)
		for f in "${_files[@]}"; do
			# guard missing file(s)
			[[ -f "$f" ]] || continue

			# Mode: any/end/start done via _re
			# We first normalize lines via cleaner (if any), then apply the regex.
			if [[ "${_brg_where}" == "start" ]]; then
				# Start-anchored match on the content (after cleaning)
				grep ')$' "$f" | grep -v '(' | ${_cleaner} | grep -E "^${_brg_query}" || true
			elif [[ "${_brg_where}" == "end" ]]; then
				# End-anchored match
				grep ')$' "$f" | grep -v '(' | ${_cleaner} | grep -E "${_brg_query})$" || true
			else
				# Anywhere
				grep ')$' "$f" | grep -v '(' | ${_cleaner} | grep -E "${_brg_query}" || true
			fi
		done

		# ---------- Source B: function headers ----------
		# Extract 'XXX' from lines like:   XXX_genus_species() {
		# Accept leading spaces; allow [A-Za-z0-9_]+ in the XXX part.
		# Then apply the same where-to-match rule to the extracted XXX.
		# for f in "${_files[@]}"; do
		# 	[[ -f "$f" ]] || continue
		# 	awk -v re="${_re}" '
		# 		BEGIN { OFS="\n" }
		# 		# Match function headers:
		# 		#   ^\s* (NAME) _genus_species () {   with flexible spaces
		# 		# Capture group 1 = NAME (prefix before _genus_species)
		# 		match($0, /^[[:space:]]*([[:alnum:]_]+)_genus_species[[:space:]]*\(\)[[:space:]]*\{/, m) {
		# 			name = m[1]
		# 			if (name ~ re) print name
		# 		}
		# 	' "$f"
		# done

		# 		for f in "${_files[@]}"; do
		# 			[[ -f "$f" ]] || continue
		# 			# 	awk -v re="${_re}" '
		# 			#   BEGIN { OFS="\n" }
		# 			#   # Match function headers like:
		# 			#   #   name_genus_species() {            # classic
		# 			#   #   function name_genus_species() {   # with "function" keyword
		# 			#   # NAME may include letters, digits, underscore, and hyphen.
		# 			#   match($0, /^[[:space:]]*(function[[:space:]]+)?([[:alnum:]_-]+)_genus_species[[:space:]]*\(\)[[:space:]]*\{/, m) {
		# 			#     name = m[2]
		# 			#     if (name ~ re) print name
		# 			#   }
		# 			# ' "$f"
		#
		# 			awk -v re="${_re}" '
		#   function norm(s){ gsub(/-/, "_", s); return s }
		#   match($0, /^[[:space:]]*(function[[:space:]]+)?([[:alnum:]_-]+)_genus_species[[:space:]]*\(\)[[:space:]]*\{/, m) {
		#     name = m[2]
		#     if (name ~ re || norm(name) ~ re) print name
		#   }
		# ' "$f"
		# 		done

		for f in "${_files[@]}"; do
			[[ -f "$f" ]] || continue
			awk -v re="${_re}" -f "${_POLAPLIB_DIR}/scripts/extract_genus_species_prefixes.awk" "$f"
		done
	} | sort -u
}

v0_list_genus_species() {
	# Files to search
	local _files=("${_POLAPLIB_DIR}/polap-lib-data.sh" "$0")
	local bolap_cmd="${FUNCNAME##*_}"

	help_message=$(
		cat <<EOF
Name:
  bolap ${bolap_cmd} - search commands

Synopsis:
  bolap ${bolap_cmd} [--query REGEX] [--where any|start|end]

Description:
  Enumerates available command prefixes (the part before "_genus_species")
  from polap-lib-data.sh and the current script. Matches are filtered by
  a regular expression and anchoring mode.

Examples:
  bolap ${bolap_cmd} --query data --where any
  bolap ${bolap_cmd} data              # legacy shorthand (same as above)
  bolap ${bolap_cmd} --where start --query '^read'
EOF
	)

	local _brg_query=".*"
	local _brg_where="any"

	# Local, isolated unknowns (don’t rely on global collector)
	local -a _unknown=()
	# Back-compat: lone positional ->  --query
	if [[ $# -eq 1 && "$1" != -* ]]; then
		set -- --query "$1" --where any
	fi

	while (($#)); do
		case "$1" in
		--query)
			shift
			_brg_query="${1:-.*}"
			[[ $# -gt 0 ]] && shift || true
			;;
		--where)
			shift
			_brg_where="${1:-any}"
			[[ $# -gt 0 ]] && shift || true
			;;
		--help | -h)
			declare -n ref="help_message"
			local manfile=$(_bolap_lib_man-convert_help_message "$ref" "${bolap_cmd}")
			man "$manfile"
			rm -f "$manfile"
			return 0
			;;
		-*)
			_unknown+=("$1")
			shift
			;;
		*) # extra positionals ignored (future use)
			shift ;;
		esac
	done

	# Build anchor once
	local _re
	case "${_brg_where}" in
	start) _re="^${_brg_query}" ;;
	end) _re="${_brg_query}\$" ;;
	*) _re="${_brg_query}" ;;
	esac

	# AWK to extract function headers once; grep for data-lines once.
	{
		# for f in "${_files[@]}"; do
		# 	[[ -f "$f" ]] || continue
		#
		# 	# Source A: data lines (whatever your project encodes as "commands" in data).
		# 	# Keep your original heuristic if those lines end with ')', without '('.
		# 	# Match anywhere vs anchored via _re.
		# 	if grep -qE . "$f"; then
		# 		# Clean (optional) then filter
		# 		if declare -F clean_input_lines >/dev/null 2>&1; then
		# 			clean_input_lines <"$f" | grep -E ')$' | grep -v '(' | grep -E "${_re}" || true
		# 		else
		# 			grep -E ')$' "$f" | grep -v '(' | grep -E "${_re}" || true
		# 		fi
		# 	fi
		#
		# 	# Source B: function headers NAME_genus_species() {
		# 	awk -v re="${_re}" '
		#       match($0, /^[[:space:]]*([[:alnum:]_]+)_genus_species[[:space:]]*\(\)[[:space:]]*\{/, m) {
		#         name = m[1]
		#         if (name ~ re) print name
		#       }' "$f"
		# done

		for f in "${_files[@]}"; do
			[[ -f "$f" ]] || continue
			awk -v re="${_re}" '
				BEGIN { OFS="\n" }
				# Match function headers:
				#   ^\s* (NAME) _genus_species () {   with flexible spaces
				# Capture group 1 = NAME (prefix before _genus_species)
				match($0, /^[[:space:]]*([[:alnum:]_]+)_genus_species[[:space:]]*\(\)[[:space:]]*\{/, m) {
					name = m[1]
					if (name ~ re) print name
				}
			' "$f"
		done
	} | sort -u
}

install_genus_species() {
	local args=("$@")

	# _log_echo "1"
	# echo "2"

	if has_help "${args[@]}"; then
		echo "$help_message_install"
		return
	fi
	# echo "3"

	for item in "${args[@]}"; do
		local tool="${item%%=*}"
		local version=""
		[[ "$item" == *"="* ]] && version="${item#*=}"

		install-${tool}_genus_species "$version"
	done
}

install-ont_genus_species() {
	if [[ "${opt_y_flag}" == false ]]; then
		read -p "Do you want to install ont in the polap-ont conda environment? (y/N): " confirm
	else
		confirm="yes"
	fi

	if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
		# Check current conda environment
		# Initialize Conda for non-interactive shells
		# source "$(conda info --base)/etc/profile.d/conda.sh"
		# if [[ "$CONDA_DEFAULT_ENV" != "base" ]]; then
		# 	echo "You're not in the base environment. Chaniging 'base'..."
		# 	conda activate base
		# fi
		_polap_lib_conda-ensure_conda_env base || exit 1

		if [[ "$CONDA_DEFAULT_ENV" == "base" ]]; then
			echo "You're in the base environment. Creating 'polap-ont'..."
			if conda env list | awk '{print $1}' | grep -qx "polap-ont"; then
				echo "[INFO] Conda environment 'polap-ont' already exists."
			else
				conda create -y --name polap-ont porechop_abi
				conda activate polap-ont
				conda install -y bbmap meryl seqtk seqkit yacrd
			fi
		else
			echo "Error: You're in the '$CONDA_DEFAULT_ENV' environment. Please activate base before running this script."
			exit 1
		fi
	else
		echo "ont installation is canceled."
		echo "Check: https://github.com/bonsai-team/Porechop_ABI"
		echo "conda install -c bioconda porechop_abi"
		echo "Execute: ont"
	fi
}

install-sim_genus_species() {
	if [[ "${opt_y_flag}" == false ]]; then
		read -p "Do you want to install sim in the polap-sim conda environment? (y/N): " confirm
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
			echo "You're in the base environment. Creating 'polap-sim'..."
			if conda env list | awk '{print $1}' | grep -qx "polap-sim"; then
				echo "[INFO] Conda environment 'polap-sim' already exists."
			else
				conda create -y --name polap-sim
				conda activate polap-sim
				conda install -y seqkit seqtk
				conda install -y badread nanosim pbsim3
				conda install -y busco
			fi
		else
			echo "Error: You're in the '$CONDA_DEFAULT_ENV' environment. Please activate base before running this script."
			exit 1
		fi
	else
		echo "sim installation is canceled."
		echo "Check: https://github.com/maickrau/sim"
		echo "conda install -c bioconda sim"
		echo "Execute: sim"
	fi
}

# Wrapper that uses your _polap_lib_conda-* helpers
# Installs “abc” into env 'polap-abc' with a sensible default stack.
install-template_genus_species() {
	local want_env="polap-template"
	local -a pkgs=(
		python=3.11
		samtools mummer4 blast
		biopython pysam
		r-base r-data.table r-ggplot2 r-gridextra
		template
	)

	# Confirm (honors opt_y_flag if you set it elsewhere)
	local confirm
	if [[ "${opt_y_flag-}" == "true" ]]; then
		confirm="yes"
	else
		read -r -p "Do you want to install template in the ${want_env} conda environment? (y/N): " confirm
	fi

	if [[ "${confirm,,}" == "y" || "${confirm,,}" == "yes" ]]; then
		# Create/upgrade env (channels first for R/bioconda harmony)
		_polap_lib_conda-create-env \
			"$want_env" "${pkgs[@]}" \
			--channel conda-forge --channel bioconda -y || return 1

		# Quick sanity check that key tools are on PATH
		local t
		for t in python samtools nucmer blastn R template; do
			if ! command -v "$t" >/dev/null 2>&1; then
				echo "WARNING: '$t' not found in '$want_env' PATH." >&2
			fi
		done
	else
		echo "template installation is canceled."
		echo "Check: https://github.com/maickrau/template"
		echo "conda install -c bioconda template"
		echo "Execute: template"
	fi
}

install-tiara_genus_species() {
	local bolap_cmd="${FUNCNAME%%_*}"

	help_message=$(
		cat <<EOF
Name:
  bolap ${bolap_cmd}

Synopsis:
  bolap ${bolap_cmd}

Description:
  Install tiara.

Examples:
  Basic usage:
    tiara -i sample_input.fasta -o out.txt

  Advanced:
    tiara -i sample_input.fasta -o out.txt --tf mit pla pro -t 4 -p 0.65 0.60 --probabilities

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	local want_env="polap-tiara"
	local -a pkgs=(
		tiara
		seqtk seqkit
	)

	# Confirm (honors opt_y_flag if you set it elsewhere)
	local confirm
	if [[ "${opt_y_flag-}" == "true" ]]; then
		confirm="yes"
	else
		read -r -p "Do you want to install tiara in the ${want_env} conda environment? (y/N): " confirm
	fi

	if [[ "${confirm,,}" == "y" || "${confirm,,}" == "yes" ]]; then
		# Create/upgrade env (channels first for R/bioconda harmony)
		_polap_lib_conda-create-env \
			"$want_env" "${pkgs[@]}" \
			--channel conda-forge --channel bioconda -y || return 1

		# Quick sanity check that key tools are on PATH
		local t
		for t in tiara; do
			if ! command -v "$t" >/dev/null 2>&1; then
				echo "WARNING: '$t' not found in '$want_env' PATH." >&2
			fi
		done
	else
		echo "tiara installation is canceled."
		echo "Check: https://github.com/ibe-uw/tiara"
		echo "conda install -c conda-forge tiara"
		echo "Unfortunately currently it does work only for python 3.7 and 3.8."
		echo "License: MIT license"
		echo "Execute: tiara"
	fi
}

install-path_genus_species() {
	local want_env="polap-path"
	local -a pkgs=(
		python=3.11
	)

	# Confirm (honors opt_y_flag if you set it elsewhere)
	local confirm
	if [[ "${opt_y_flag-}" == "true" ]]; then
		confirm="yes"
	else
		read -r -p "Do you want to install path in the ${want_env} conda environment? (y/N): " confirm
	fi

	if [[ "${confirm,,}" == "y" || "${confirm,,}" == "yes" ]]; then
		# Create/upgrade env (channels first for R/bioconda harmony)
		_polap_lib_conda-create-env \
			"$want_env" "${pkgs[@]}" \
			--channel conda-forge --channel bioconda -y || return 1

		# Quick sanity check that key tools are on PATH
		local t
		for t in python samtools nucmer blastn R path; do
			if ! command -v "$t" >/dev/null 2>&1; then
				echo "WARNING: '$t' not found in '$want_env' PATH." >&2
			fi
		done
	else
		echo "path installation is canceled."
		echo "Check: https://github.com/maickrau/path"
		echo "conda install -c bioconda path"
		echo "Execute: path"
	fi
}

install-polish_genus_species() {
	local want_env="polap-polish"
	local -a pkgs=(
		seqtk seqkit
		minimap2 gfatools
		pysam
		fmlrc2 ropebwt2 bwa-mem2 samtools polypolish
		meryl merqury r-base python
		conda-forge::r-data.table
		bowtie2
		cmake=3.26.* ninja gcc=12 gxx=12 zlib git make
	)

	# Confirm (honors opt_y_flag if you set it elsewhere)
	local confirm
	if [[ "${opt_y_flag-}" == "true" ]]; then
		confirm="yes"
	else
		read -r -p "Do you want to install polish in the ${want_env} conda environment? (y/N): " confirm
	fi

	if [[ "${confirm,,}" == "y" || "${confirm,,}" == "yes" ]]; then
		# Create/upgrade env (channels first for R/bioconda harmony)
		_polap_lib_conda-create-env \
			"$want_env" "${pkgs[@]}" \
			--channel conda-forge --channel bioconda -y || return 1

		# Quick sanity check that key tools are on PATH
		local t
		for t in merqury.sh; do
			if ! command -v "$t" >/dev/null 2>&1; then
				echo "WARNING: '$t' not found in '$want_env' PATH." >&2
			fi
		done

		setup-racon_genus_species

		# # Install racon
		# rm -rf racon
		# git clone --recursive https://github.com/lbcb-sci/racon.git
		# cd racon
		# git fetch --tags --all
		#
		# # Pick a tag automatically (prefer 1.5.x if present; else latest semver tag)
		# TAG="$(git tag -l | grep -E '^v?1\.5(\.[0-9]+)?$' | sort -V | tail -1)"
		# if [ -z "$TAG" ]; then
		# 	TAG="$(git tag -l | grep -E '^v?[0-9]+\.[0-9]+(\.[0-9]+)?$' | sort -V | tail -1)"
		# fi
		#
		# if [ -z "$TAG" ]; then
		# 	echo "[ERROR] No semver-like tags found; building from current default branch."
		# else
		# 	echo "[INFO] Using racon tag: $TAG"
		# 	git checkout -B "build-$TAG" "$TAG"
		# fi
		#
		# # Resync submodules to match the tag
		# git submodule sync --recursive || true
		# git submodule update --init --recursive
		#
		# # Clean build dir
		# rm -rf build && mkdir build && cd build
		#
		# # Configure with pinned compilers, zlib, and policies (edlib SIMD off)
		# cmake -G Ninja \
		# 	-DCMAKE_BUILD_TYPE=Release \
		# 	-DCMAKE_INSTALL_PREFIX="$CONDA_PREFIX" \
		# 	-DCMAKE_C_COMPILER="$CONDA_PREFIX/bin/gcc" \
		# 	-DCMAKE_CXX_COMPILER="$CONDA_PREFIX/bin/g++" \
		# 	-DCMAKE_CXX_STANDARD=17 \
		# 	-Dracon_build_tests=OFF \
		# 	-DEDLIB_FORCE_DISABLE_SIMD=ON \
		# 	-DZLIB_ROOT="$CONDA_PREFIX" \
		# 	-DZLIB_INCLUDE_DIR="$CONDA_PREFIX/include" \
		# 	-DZLIB_LIBRARY="$CONDA_PREFIX/lib/libz.so" \
		# 	-DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
		# 	-DCMAKE_POLICY_DEFAULT_CMP0169=OLD \
		# 	-DCMAKE_POLICY_DEFAULT_CMP0074=NEW \
		# 	-Wno-dev \
		# 	..
		#
		# ninja -j"$(nproc)"
		# ninja install
		#
		# echo "Installed: $CONDA_PREFIX/bin/racon"
		# "$CONDA_PREFIX/bin/racon" --version
		echo "[INFO] we have created conda env $CONDA_PREFIX" >&2
	else
		echo "polap-polish env is not created"
	fi
}

setup-busco_genus_species() {
	${_polap_cmd} data busco
}

setup-racon_genus_species() {

	_polap_lib_conda-ensure_conda_env polap-polish || exit 1

	# Install racon
	rm -rf racon
	git clone --recursive https://github.com/lbcb-sci/racon.git
	cd racon
	git fetch --tags --all

	# Pick a tag automatically (prefer 1.5.x if present; else latest semver tag)
	TAG="$(git tag -l | grep -E '^v?1\.5(\.[0-9]+)?$' | sort -V | tail -1)"
	if [ -z "$TAG" ]; then
		TAG="$(git tag -l | grep -E '^v?[0-9]+\.[0-9]+(\.[0-9]+)?$' | sort -V | tail -1)"
	fi

	if [ -z "$TAG" ]; then
		echo "[ERROR] No semver-like tags found; building from current default branch."
	else
		echo "[INFO] Using racon tag: $TAG"
		git checkout -B "build-$TAG" "$TAG"
	fi

	# Resync submodules to match the tag
	git submodule sync --recursive || true
	git submodule update --init --recursive

	# Clean build dir
	rm -rf build && mkdir build && cd build

	# Configure with pinned compilers, zlib, and policies (edlib SIMD off)
	cmake -G Ninja \
		-DCMAKE_BUILD_TYPE=Release \
		-DCMAKE_INSTALL_PREFIX="$CONDA_PREFIX" \
		-DCMAKE_C_COMPILER="$CONDA_PREFIX/bin/gcc" \
		-DCMAKE_CXX_COMPILER="$CONDA_PREFIX/bin/g++" \
		-DCMAKE_CXX_STANDARD=17 \
		-Dracon_build_tests=OFF \
		-DEDLIB_FORCE_DISABLE_SIMD=ON \
		-DZLIB_ROOT="$CONDA_PREFIX" \
		-DZLIB_INCLUDE_DIR="$CONDA_PREFIX/include" \
		-DZLIB_LIBRARY="$CONDA_PREFIX/lib/libz.so" \
		-DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
		-DCMAKE_POLICY_DEFAULT_CMP0169=OLD \
		-DCMAKE_POLICY_DEFAULT_CMP0074=NEW \
		-Wno-dev \
		..

	ninja -j"$(nproc)"
	ninja install

	echo "Installed: $CONDA_PREFIX/bin/racon"
	"$CONDA_PREFIX/bin/racon" --version
	echo "[INFO] we have installed racon in the conda env $CONDA_PREFIX." >&2

	conda deactivate
}

v1_install-polish_genus_species() {
	local want_env="polap-polish-not-working-in-other-linux-box"
	local -a pkgs=(
		minimap2 racon gfatools fmlrc2 ropebwt2 bwa-mem2 samtools polypolish
		seqtk
		meryl merqury r-base python
	)

	# Confirm (honors opt_y_flag if you set it elsewhere)
	local confirm
	if [[ "${opt_y_flag-}" == "true" ]]; then
		confirm="yes"
	else
		read -r -p "Do you want to install polish in the ${want_env} conda environment? (y/N): " confirm
	fi

	if [[ "${confirm,,}" == "y" || "${confirm,,}" == "yes" ]]; then
		# Create/upgrade env (channels first for R/bioconda harmony)
		_polap_lib_conda-create-env \
			"$want_env" "${pkgs[@]}" \
			--channel conda-forge --channel bioconda -y || return 1

		# Quick sanity check that key tools are on PATH
		local t
		for t in merqury.sh; do
			if ! command -v "$t" >/dev/null 2>&1; then
				echo "WARNING: '$t' not found in '$want_env' PATH." >&2
			fi
		done
	else
		echo "polap-polish env is not created"
	fi
}

install-abc-old-version_genus_species() {
	if [[ "${opt_y_flag}" == false ]]; then
		read -p "Do you want to install abc in the polap-abc conda environment? (y/N): " confirm
	else
		confirm="yes"
	fi

	if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then

		# --- Conda bootstrap (works under `set -u`) ---
		# Try to make `conda activate` work in a non-interactive shell
		if ! command -v conda >/dev/null 2>&1; then
			for d in "$HOME/miniconda3" "$HOME/mambaforge" "$HOME/anaconda3" "/opt/conda"; do
				if [[ -f "$d/etc/profile.d/conda.sh" ]]; then
					# shellcheck disable=SC1090
					source "$d/etc/profile.d/conda.sh"
					break
				fi
			done
		fi
		# Enable `conda activate` regardless of how Conda is installed
		command -v conda >/dev/null 2>&1 && eval "$(conda shell.bash hook)" || {
			echo "ERROR: Could not initialize Conda in this shell." >&2
			exit 1
		}

		# Helper: robust current env name (safe with `set -u`)
		_current_env() {
			if [[ -n "${CONDA_PREFIX-}" ]]; then
				basename -- "$CONDA_PREFIX"
			elif [[ -n "${CONDA_DEFAULT_ENV-}" ]]; then
				printf '%s' "$CONDA_DEFAULT_ENV"
			else
				printf ''
			fi
		}

		want_env="polap-abc"
		cur_env="$(_current_env)"

		# We don't need to hop to 'base' first; Conda can switch directly.
		if [[ "$cur_env" != "$want_env" && -n "$cur_env" ]]; then
			echo "You're in '$cur_env'. Switching to '$want_env' (no stacking)."
		fi

		# Create env if missing
		if ! conda env list | awk '{print $1}' | grep -qx "$want_env"; then
			echo "Creating Conda environment '$want_env'..."
			# Prefer creating with all packages at once; add -c conda-forge for R stack
			conda create -y -n "$want_env" \
				python=3.11 \
				samtools mummer4 blast \
				biopython pysam \
				r-base r-data.table r-ggplot2 r-gridextra
		else
			echo "Conda environment '$want_env' already exists."
		fi

		# Activate and verify
		conda activate "$want_env" || {
			echo "ERROR: Failed to activate '$want_env'." >&2
			exit 1
		}

		echo "Activated env: $(_current_env)"

		# (optional) quick sanity check that core tools are on PATH
		for tool in python samtools nucmer blastn R; do
			command -v "$tool" >/dev/null 2>&1 || {
				echo "WARNING: '$tool' not found in '$want_env' PATH." >&2
			}
		done

	else
		echo "abc installation is canceled."
		echo "Check: https://github.com/maickrau/abc"
		echo "conda install -c bioconda abc"
		echo "Execute: abc"
	fi
}

install-evo_genus_species() {
	local want_env="polap-evo"
	local -a pkgs=(
		python=3.11
		# core mappers / consensus
		minimap2=2.*
		racon=1.*
		polypolish=0.*
		bowtie2=2.*
		bwa-mem2=2.*
		samtools=1.*
		# MTPT detection / annotation / repeats
		blast=2.*
		mummer4=4.*
		hmmer=3.*
		bedtools=2.*
		seqkit=2.*
		# clustering / alignment
		mmseqs2=15.*
		cd-hit=4.*
		mafft=7.*
		# cp phylogenomics
		mauve=2.*
		iqtree=2.*
		# coverage
		mosdepth=0.*
		# R stack
		r-base=4.*
		r-data.table
		r-ggplot2
		r-ape
		r-phangorn
		r-phytools
		r-castor
		r-quantreg
		bioconductor-ggtree
		# Python stack
		pysam
		pandas
		numpy
		biopython
		matplotlib
	)

	# Confirm (honors opt_y_flag if you set it elsewhere)
	local confirm
	if [[ "${opt_y_flag-}" == "true" ]]; then
		confirm="yes"
	else
		read -r -p "Do you want to install evo in the ${want_env} conda environment? (y/N): " confirm
	fi

	if [[ "${confirm,,}" == "y" || "${confirm,,}" == "yes" ]]; then
		# Create/upgrade env (channels first for R/bioconda harmony)
		_polap_lib_conda-create-env \
			"$want_env" "${pkgs[@]}" \
			--channel conda-forge --channel bioconda -y || return 1

		# Quick sanity check that key tools are on PATH
		local t
		for t in iqtree; do
			if ! command -v "$t" >/dev/null 2>&1; then
				echo "WARNING: '$t' not found in '$want_env' PATH." >&2
			fi
		done
	else
		echo "evo installation is canceled."
		echo "Check: https://github.com/maickrau/evo"
		echo "conda install -c bioconda evo"
		echo "Execute: evo"
	fi
}

install-dev_genus_species() {
	if [[ "${opt_y_flag}" == false ]]; then
		read -p "Do you want to install dev in the polap-dev conda environment? (y/N): " confirm
	else
		confirm="yes"
	fi

	if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
		_polap_lib_conda-ensure_conda_env base || exit 1

		if [[ "$CONDA_DEFAULT_ENV" == "base" ]]; then
			echo "You're in the base environment. Creating 'polap-dev'..."
			if conda env list | awk '{print $1}' | grep -qx "polap-dev"; then
				echo "[INFO] Conda environment 'polap-dev' already exists."
			else
				echo conda create -y --name polap-dev
				echo conda activate polap-dev
				echo conda install -y cd-hit
			fi
		else
			echo "Error: You're in the '$CONDA_DEFAULT_ENV' environment. Please activate base before running this script."
			exit 1
		fi
	else
		echo "dev installation is canceled."
		echo https://github.com/weizhongli/cdhit
	fi
}

install-ncbitools_genus_species() {
	if [[ "${opt_y_flag}" == false ]]; then
		read -p "Do you want to install ncbitools in the polap-ncbitools conda environment? (y/N): " confirm
	else
		confirm="yes"
	fi

	if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
		# Check current conda environment
		# Initialize Conda for non-interactive shells
		_polap_lib_conda-ensure_conda_env base || exit 1

		# source "$(conda info --base)/etc/profile.d/conda.sh"
		# if [[ "$CONDA_DEFAULT_ENV" != "base" ]]; then
		# 	echo "You're not in the base environment. Chaniging 'base'..."
		# 	conda activate base
		# fi
		#
		if [[ "$CONDA_DEFAULT_ENV" == "base" ]]; then
			echo "You're in the base environment. Creating 'polap-ncbitools'..."
			if conda env list | awk '{print $1}' | grep -qx "polap-ncbitools"; then
				echo "[INFO] Conda environment 'polap-ncbitools' already exists."
			else
				# conda create -y --name polap-ncbitools blast entrez-direct sra-tools
				# conda activate polap-ncbitools
				# conda install -y seqtk seqkit
				# conda install -y busco
				# conda install -y conda-forge::pv
				# _polap_lib_conda-ensure_conda_env_deactivate

				conda create -y --name polap-ncbitools blast entrez-direct sra-tools \
					seqtk seqkit busco conda-forge::pv
			fi
		else
			echo "Error: You're in the '$CONDA_DEFAULT_ENV' environment. Please activate base before running this script."
			exit 1
		fi
		_polap_lib_conda-ensure_conda_env_deactivate
		# set +u
		# conda deactivate
		# set -u
	else
		echo "ncbitools installation is canceled."
		echo "conda install -c bioconda blast entrez-direct sra-tools"
		echo "https://www.ivarch.com/programs/pv.shtml"
	fi
}

install-graphaligner_genus_species() {
	if [[ "${opt_y_flag}" == false ]]; then
		read -p "Do you want to install graphaligner in the polap-graphaligner conda environment? (y/N): " confirm
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
			echo "You're in the base environment. Creating 'polap-graphaligner'..."
			if conda env list | awk '{print $1}' | grep -qx "polap-graphaligner"; then
				echo "[INFO] Conda environment 'polap-graphaligner' already exists."
			else
				conda create -y --name polap-graphaligner graphaligner
			fi
		else
			echo "Error: You're in the '$CONDA_DEFAULT_ENV' environment. Please activate base before running this script."
			exit 1
		fi
	else
		echo "graphaligner installation is canceled."
		echo "Check: https://github.com/maickrau/GraphAligner"
		echo "conda install -c bioconda graphaligner"
		echo "Execute: GraphAligner"
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
				echo "[INFO] Conda environment 'polap-getorganelle' already exists."
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
				echo "[INFO] Conda environment 'polap-tippo' already exists."
			else
				#please specific the python version to 3.8 :)
				conda create -y --name polap-tippo python=3.8
				conda activate polap-tippo
				set +u
				conda install -y bioconda::tipp
				set -u
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
				echo "[INFO] Conda environment 'polap-oatk' already exists."
			else
				conda create -y --name polap-oatk bioconda::oatk
				conda activate polap-oatk
				conda install -y biopython
				conda install -y hmmer seqtk mafft parallel entrez-direct
				conda install -y hifiasm
				conda install -y minimap2 seqtk seqkit samtools gfatools
				# conda install -y rust-mdbg ntcard racon medaka
				conda install -y rust-mdbg ntcard racon
				conda install -y filtlong
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

setup-oatk_genus_species() {

	git clone https://github.com/c-zhou/OatkDB.git

	echo "[INFO] create: OatkDB" >&2
}

install-conda_genus_species() {
	if [[ "${opt_y_flag}" == false ]]; then
		read -p "Do you want to install miniconda3? (y/N): " confirm
	else
		confirm="yes"
	fi

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
				echo "[INFO] Conda environment 'polap-mbg' already exists."
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

	echo "[INFO] Start an organelle-genome assembly pipeline on ${_brg_outdir}-${_brg_inum}"

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

# 2025-10-31
# not used much
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

check_and_prepare_fastq() {
	local long_sra="$1"
	local media_dir="$2"
	local media1_dir="$3"
	local remote_host="$4"
	local dry_run="$5"
	local __result_var="$6" # output: resolved file path

	local result=""
	local dirs=("$media_dir" "$media1_dir" "$PWD")

	# echo "${remote_host}" >&2
	# return 0

	# 1. Try local working directory
	# for ext in ".fastq" "-*.fq" ".fq.gz" ".fastq.gz"; do
	for ext in ".fastq" ".fq" ".fq.gz" ".fastq.gz"; do
		match=$(ls ${long_sra}${ext} 2>/dev/null | head -n 1)
		if [[ -s "$match" ]]; then
			result="$match"
			eval $__result_var="'$result'"
			return 0
		fi
	done

	for dir in "${dirs[@]}"; do
		local archive="${long_sra}.fastq.tar.gz"
		local full_archive_path="${dir}/${archive}"

		# 2a. Local archive
		if [[ -f "$full_archive_path" ]]; then
			echo "Found local archive: $full_archive_path"
			if [[ "$dry_run" != true ]]; then
				echo "  decompressing ..."
				if command -v pv >/dev/null 2>&1; then
					pv -p -t -e -r -b "$full_archive_path" | tar xzf - -C .
				else
					tar -zxf "$full_archive_path"
				fi
			fi
			if [[ -s "${long_sra}.fastq" ]]; then
				result="${long_sra}.fastq"
			else
				result="no-such-file.fastq"
			fi
			# result=$(tar -tzf "$full_archive_path" | head -n 1)
			eval $__result_var="'$result'"
			return 0
		fi

		# 2b. Remote archive
		if [[ -n "$remote_host" ]]; then
			if ssh "$remote_host" "test -f '$full_archive_path'" 2>/dev/null; then
				echo "Found remote archive: $remote_host:$full_archive_path"
				if [[ "$dry_run" != true ]]; then
					# scp "$remote_host:$full_archive_path" . && tar -zxf "$archive"
					scp "$remote_host:$full_archive_path" .
					if [[ -f "$archive" ]]; then
						echo "  decompressing ..."
						tar -zxf "$archive"
					fi
				fi
				result=$(tar -tzf "$archive" | head -n 1)
				eval $__result_var="'$result'"
				return 0
			fi
		fi

		# 2c. Local FASTQ
		match_plain=$(find "$dir" -maxdepth 1 -type f \
			\( -name "${long_sra}*.fq" -o -name "${long_sra}*.fastq" -o -name "${long_sra}*.fq.gz" -o -name "${long_sra}*.fastq.gz" -o -name "${long_sra}*.fastq.tar.gz" \) |
			head -n 1)
		if [[ -n "$match_plain" && -s "$match_plain" ]]; then
			echo "Found local: $match_plain"
			if [[ "$dry_run" != true ]]; then
				ln -s "$match_plain" .
			fi
			if [[ "$match_plain" == *.tar.gz ]]; then
				tar -zxf "$match_plain"
			fi
			result=$(basename "$match_plain")
			eval $__result_var="'$result'"
			return 0
		fi

		# 2d. Remote FASTQ
		if [[ -n "$remote_host" ]]; then
			remote_match=$(ssh "$remote_host" "find '$dir' -maxdepth 1 -type f \
        \\( -name '${long_sra}*.fq' -o -name '${long_sra}*.fastq' -o -name '${long_sra}*.fastq.tar.gz' \\) | head -n 1" 2>/dev/null)
			if [[ -n "$remote_match" ]]; then
				echo "Found remote FASTQ: $remote_match"
				if [[ "$dry_run" != true ]]; then
					scp "$remote_host:$remote_match" .
				fi
				if [[ "$remote_match" == *.tar.gz ]]; then
					echo "  decompressing ..."
					tar -zxf $(basename "$remote_match")
				fi
				result=$(basename "$remote_match")
				eval $__result_var="'$result'"
				return 0
			fi
		fi
	done

	_polap_lib_conda-ensure_conda_env polap || exit 1
	_log_echo "  downloading long-read SRA ID: ${long_sra} ... be patient!"
	bash "${_polap_script_bin_dir}"/polap-ncbitools fetch sra "$long_sra"
	rm -rf "${long_sra}"
	conda deactivate

	# 1. Try local working directory
	# for ext in ".fastq" "-*.fq" ".fq.gz" ".fastq.gz"; do
	for ext in ".fastq" ".fq" ".fq.gz" ".fastq.gz"; do
		match=$(ls ${long_sra}${ext} 2>/dev/null | head -n 1)
		if [[ -s "$match" ]]; then
			result="$match"
			eval $__result_var="'$result'"
			return 0
		fi
	done

	echo "❌ Could not locate FASTQ for $long_sra"
	return 1
}

data-long_genus_species() {
	# local _brg_outdir="$1"
	# local _brg_inum="${2:-0}"
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

	resolved_fastq=""
	if check_and_prepare_fastq "$long_sra" "$_media_dir" "$_media1_dir" "$_local_host" false resolved_fastq; then
		echo "$resolved_fastq" >"${_brg_outdir}/${_brg_inum}/l.fastq.path.txt"
		echo "$resolved_fastq"
	else
		echo "FASTQ file not found or extracted"
	fi

	if [[ "${_POLAP_DEBUG}" == "1" ]]; then
		_log_echo "[INFO] long-read: ${long_sra} is ready."
	fi
}

x_data-long_genus_species() {
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

# test
# template
run-a-b-c_genus_species() {
	_log_echo0 "outdir: $_brg_outdir"
	_log_echo0 "sindex: $_brg_sindex"

	source "${_POLAPLIB_DIR}/polap-variables-data.sh"

	_log_echo0 dataset-get-field "$_brg_target" long
	local long_sra1=$(dataset-get-field "$_brg_target" long)
	_log_echo0 "long SRA1: $long_sra1"

	# dataset-cache-fields long
	# _log_echo0 "${_dataset["$_brg_target|long"]}"
	# local long_sra="${_dataset["$_brg_target|long"]}"
	# _log_echo0 "long SRA: $long_sra"

	local bolap_cmd="${FUNCNAME%%_*}"

	help_message=$(
		cat <<EOF
Name:
  bolap - $bolap_cmd

Synopsis:
  bolap $bolap_cmd

Description:
  bolap

Examples:
  Execute polap readassemble:
    bolap readassemble -s Vigna_radiata

  Execute polap disassemble:
    bolap disassemble -s Vigna_radiata

  Execute polap syncassemble:
    bolap syncassemble -s Vigna_radiata

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	local _abc_arg1=1
	local _abc_arg2=1

	parse_commandline() {
		set -- "${_brg_unknown_opts[@]}"

		# source "${_POLAPLIB_DIR}/polap-cmd-version.sh" # '.' means 'source'
		while test $# -gt 0; do
			_key="$1"
			case "$_key" in
			--abc-arg1)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_abc_arg1="$2"
					shift || true
				fi
				;;
			--abc-arg2)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_abc_arg2="$2"
					shift || true
				fi
				;;
			-*)
				_log_echo0 "[ERROR] no such options: $1"
				;;
			*)
				break
				;;
			esac
			shift || true
		done
	}

	declare -n ref="help_message"
	if [[ "${_brg_help}" == "on" ]]; then
		local manfile=$(_bolap_lib_man-convert_help_message "$ref" "${bolap_cmd}")
		man "$manfile"
		rm -f "$manfile"
		return
	fi

	_log_echo0 "A: abc-arg1: $_abc_arg1"
	parse_commandline

	_log_echo0 "abc-arg1: $_abc_arg1"
	_log_echo0 "abc-arg2: $_abc_arg2"
}

# 2025-10-25
run-polap-assemble1_genus_species() {
	local bolap_cmd="${FUNCNAME%%_*}"

	help_message=$(
		cat <<EOF
Name:
  bolap - $bolap_cmd

Synopsis:
  bolap $bolap_cmd

Description:
  bolap

Examples:
  Execute polap assemble1:
    bolap run data-long -s Vigna_radiata
    bolap run data-short -s Vigna_radiata
    bolap run summary-data -s Vigna_radiata
    bolap data downsample long -s Vigna_radiata --coverage 100g
    bolap data downsample short -s Vigna_radiata --coverage 100g
    bolap run polap assemble1 -s Vigna_radiata

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	# defaults
	local _brg_type="miniasm"
	local _brg_omega="1500"
	local _brg_downsample="no-downsample"
	local _brg_redo="off"
	local _brg_cleanup="on"

	parse_commandline() {
		# set -- "${_brg_unknown_opts[@]}"
		set -- "${_brg_args[@]}"
		# _log_echo0 "parse: $@"

		# source "${_POLAPLIB_DIR}/polap-cmd-version.sh" # '.' means 'source'
		while test $# -gt 0; do
			_key="$1"
			case "$_key" in
			-s)
				test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
				_brg_outdir="${2%/}"
				_brg_outdir_list+=("$_brg_outdir") # append to array
				shift
				;;
			-i)
				test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
				_brg_sindex="$2"
				shift
				;;
			--type)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_type="$2"
					shift || true
				fi
				;;
			-w)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_omega="$2"
					shift || true
				fi
				;;
			--redo)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_redo="on"
					shift || true
				fi
				;;
			--no-cleanup)
				# _log_echo0 "--no-cleanup"
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_cleanup="off"
					shift || true
				fi
				;;
			-d)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_downsample="downsample"
					shift || true
				fi
				;;
			--dry-run)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_dry="on"
					shift || true
				fi
				;;
			esac
			shift || true
		done
	}

	_polap_lib_help-maybe-show "$bolap_cmd" help_message || return 0

	# required option: -s
	#
	if [[ "${_brg_outdir}" == "." ]]; then
		_log_echo0 "-s Genus_species"
		return 0
	fi

	source "${_POLAPLIB_DIR}/polap-variables-data.sh"
	local long_sra="${_long["$_brg_target"]}"
	local short_sra="${_short["$_brg_target"]}"

	parse_commandline
	_log_echo1 "target: $_brg_target"

	mkdir -p "${_brg_outdir_i}"

	local platform="${_platform["$_brg_target"]}"

	local _brg_coverage="50m"

	local option_data_type="--nano-raw"

	_log_echo1 "Asemble organelle genomes using miniasm as a reference generator with ${option_data_type}"

	# Always redo
	rm -rf "${_brg_rundir}"
	mkdir -p "${_brg_rundir}"

	# Always redo
	rm -rf "${_brg_target}"
	mkdir -p "${_brg_target}"

	# Start memory logger
	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"

	# check out input data
	# ln -sf "${_brg_tmpdir}/l.fq" "${long_sra}.fastq"
	# _brg_input_data="${long_sra}.fastq"
	if [[ ! -s "${long_sra}".fastq ]]; then
		# data-long_genus_species "${_brg_outdir}"
		_log_echo0 "No input data: ${long_sra}.fastq"
		_log_echo0 "run data long -s ${_brg_outdir}"
		_log_echo0 "data downsample long -s ${_brg_outdir}"
		return 0
	fi

	local lfq="${long_sra}.fastq"
	local sfq1="${short_sra}_1.fastq"
	local sfq2="${short_sra}_2.fastq"

	# elif [[ "${_brg_type}" == "miniasm" ]]; then
	_log_echo2 "plant mt $_brg_verbose_str"

	${_polap_cmd} assemble1 \
		"${option_data_type}" \
		-l "${lfq}" \
		-a "${sfq1}" \
		-b "${sfq2}" \
		${_brg_verbose_str} \
		-o "${_brg_target}"

	# Summarize results after job (with previously defined summary function)
	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"
	_polap_lib_timing-get_system_info >>"${_timing_txt}"

	# Save some results
	# rsync -azuq --max-size=5M \
	# 	"${_brg_target}/" "${_brg_rundir}/"
	#
	# if [[ "${_brg_cleanup}" == "on" ]]; then
	# 	if [[ -d "${_brg_target}" ]]; then
	# 		_log_echo1 rm -rf "${_brg_target}"
	# 		rm -rf "${_brg_target}"
	# 	fi
	# fi
}

# select seed contigs
#

# 2025-10-26
run-polap-assemble2_genus_species() {
	local bolap_cmd="${FUNCNAME%%_*}"

	help_message=$(
		cat <<EOF
Name:
  bolap - $bolap_cmd

Synopsis:
  bolap $bolap_cmd

Description:
  bolap

Examples:
  Execute polap assemble2:
    bolap run data-long -s Vigna_radiata
    bolap run data-short -s Vigna_radiata
    bolap run summary-data -s Vigna_radiata
    bolap data downsample long -s Vigna_radiata --coverage 100g
    bolap data downsample short -s Vigna_radiata --coverage 100g
    bolap run polap assemble1 -s Vigna_radiata

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	# defaults
	local _brg_type="miniasm"
	local _brg_omega="1500"
	local _brg_downsample="no-downsample"
	local _brg_redo="off"
	local _brg_cleanup="on"

	parse_commandline() {
		# set -- "${_brg_unknown_opts[@]}"
		set -- "${_brg_args[@]}"
		# _log_echo0 "parse: $@"

		# source "${_POLAPLIB_DIR}/polap-cmd-version.sh" # '.' means 'source'
		while test $# -gt 0; do
			_key="$1"
			case "$_key" in
			-s)
				test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
				_brg_outdir="${2%/}"
				_brg_outdir_list+=("$_brg_outdir") # append to array
				shift
				;;
			-i)
				test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
				_brg_sindex="$2"
				shift
				;;
			--type)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_type="$2"
					shift || true
				fi
				;;
			-w)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_omega="$2"
					shift || true
				fi
				;;
			--redo)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_redo="on"
					shift || true
				fi
				;;
			--no-cleanup)
				# _log_echo0 "--no-cleanup"
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_cleanup="off"
					shift || true
				fi
				;;
			-d)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_downsample="downsample"
					shift || true
				fi
				;;
			--dry-run)
				if test $# -lt 1; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_dry="on"
					shift || true
				fi
				;;
			-*)
				_log_echo0 "[INFO] no such options: $1"
				;;
			*) ;;
			esac
			shift || true
		done
	}

	_polap_lib_help-maybe-show "$bolap_cmd" help_message || return 0

	# required option: -s
	#
	if [[ "${_brg_outdir}" == "." ]]; then
		_log_echo0 "-s Genus_species"
		return 0
	fi

	source "${_POLAPLIB_DIR}/polap-variables-data.sh"
	local long_sra="${_long["$_brg_target"]}"
	local short_sra="${_short["$_brg_target"]}"

	parse_commandline
	_log_echo1 "target: $_brg_target"

	mkdir -p "${_brg_outdir_i}"

	local platform="${_platform["$_brg_target"]}"

	local _brg_coverage="50m"

	local option_data_type="--nano-raw"

	_log_echo1 "Asemble organelle genomes using miniasm as a reference generator with ${option_data_type}"

	# Always redo
	rm -rf "${_brg_rundir}"
	mkdir -p "${_brg_rundir}"

	# Always redo
	rm -rf "${_brg_target}"
	mkdir -p "${_brg_target}"

	# Start memory logger
	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"

	# check out input data
	# ln -sf "${_brg_tmpdir}/l.fq" "${long_sra}.fastq"
	# _brg_input_data="${long_sra}.fastq"
	if [[ ! -s "${long_sra}".fastq ]]; then
		# data-long_genus_species "${_brg_outdir}"
		_log_echo0 "No input data: ${long_sra}.fastq"
		_log_echo0 "run data long -s ${_brg_outdir}"
		_log_echo0 "data downsample long -s ${_brg_outdir}"
		return 0
	fi

	local lfq="${long_sra}.fastq"
	local sfq1="${short_sra}_1.fastq"
	local sfq2="${short_sra}_2.fastq"

	# elif [[ "${_brg_type}" == "miniasm" ]]; then
	_log_echo2 "plant mt $_brg_verbose_str"

	_log_echo0 ${_polap_cmd} assemble2 \
		"${option_data_type}" \
		-l "${lfq}" \
		-a "${sfq1}" \
		-b "${sfq2}" \
		${_brg_verbose_str} \
		-o "${_brg_target}"

	# ${_polap_cmd} assemble1 \
	# 	"${option_data_type}" \
	# 	-l "${lfq}" \
	# 	-a "${sfq1}" \
	# 	-b "${sfq2}" \
	# 	${_brg_verbose_str} \
	# 	-o "${_brg_target}"

	# Summarize results after job (with previously defined summary function)
	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"
	_polap_lib_timing-get_system_info >>"${_timing_txt}"

	# Save some results
	# rsync -azuq --max-size=5M \
	# 	"${_brg_target}/" "${_brg_rundir}/"
	#
	# if [[ "${_brg_cleanup}" == "on" ]]; then
	# 	if [[ -d "${_brg_target}" ]]; then
	# 		_log_echo1 rm -rf "${_brg_target}"
	# 		rm -rf "${_brg_target}"
	# 	fi
	# fi
}

run-data_genus_species() {
	run-data-long_genus_species
	run-data-short_genus_species
}

run-data-long_genus_species() {
	local bolap_cmd="${FUNCNAME##*_}"

	help_message=$(
		cat <<EOF
Name:
  bolap - prepare long-read data at a species directory

Synopsis:
  bolap run-data-long

Description:
  bolap

Examples:
  Execute polap readassemble:
    bolap run-data-long -s Vigna_radiata

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	local _data_long_cleanup=false
	local _data_long_redo=false
	local _remote=""

	parse_commandline() {
		# set -- "$@"
		set -- "${_brg_args[@]}"

		# source "${_POLAPLIB_DIR}/polap-cmd-version.sh" # '.' means 'source'
		while test $# -gt 0; do
			_key="$1"
			case "$_key" in
			--cleanup)
				_data_long_cleanup="true"
				;;
			--redo)
				_data_long_redo="true"
				;;
			--remote)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'."
				else
					_remote="$2"
					shift || true
				fi
				;;
			-s)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'."
				else
					_brg_outdir="$2"
					shift || true
				fi
				;;
			-i)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'."
				else
					_brg_sindex="$2"
					shift || true
				fi
				;;
			esac
			shift || true
		done
	}

	declare -n ref="help_message"
	if [[ "${_brg_help}" == "on" ]]; then
		local manfile=$(_bolap_lib_man-convert_help_message "$ref" "${bolap_cmd}")
		man "$manfile"
		rm -f "$manfile"
		return
	fi

	parse_commandline "$@"

	# echo "s:${_brg_outdir}" >&2
	# echo "i:${_brg_sindex}" >&2

	source "${_POLAPLIB_DIR}/polap-variables-data.sh"

	# check long-read in the following order
	# where long_sra=l
	# l.fastq
	# l.fq
	# l.fastq.gz
	# l.fq.gz
	# l.fastq.tar.gz
	# l.fq.tar.gz
	# storage: check each dir for the filenames as above
	# $media_dir
	# $media1_dir
	# $media2_dir
	# remote: or the dispatcher or _local_host's storage
	# if not check NCBI SRA database
	#
	# save it outdir/tmp/l.fq

	_polap_lib_conda-ensure_conda_env polap-ncbitools || exit 1

	local args=(
		-l "${long_sra}"
		-o "${_brg_outdir}"
	)

	# Safe defaults under `set -u`
	_data_long_redo="${_data_long_redo:-false}"
	_data_long_cleanup="${_data_long_cleanup:-false}"

	if [[ "${_brg_redo,,}" != "false" ]]; then
		args+=(--redo)
	fi

	if [[ "${_brg_cleanup,,}" != "false" ]]; then
		args+=(--cleanup)
	fi

	if [[ -n "${_remote}" ]]; then
		args+=(--remote ${_remote})
	fi

	bash "${_POLAPLIB_DIR}/polap-bash-run-data-long.sh" "${args[@]}"

	_polap_lib_conda-ensure_conda_env_deactivate
	# set +u
	# conda deactivate
	# set -u

	# redownload

}

# Version: v0.1.0
run-data-short_genus_species() {
	# _log_echo0 "outdir: $_brg_outdir"
	# _log_echo0 "sindex: $_brg_sindex"
	# return 0

	local bolap_cmd="${FUNCNAME##*_}"

	help_message=$(
		cat <<EOF
Name:
  bolap ${bolap_cmd} - prepare short-read data at a species directory

Synopsis:
  bolap ${bolap_cmd} -s Genus_species [-i 0] [--data-short-redo BOOL] [--data-short-cleanup BOOL]

Description:
  Prepare paired short-read FASTQ as:
    sfolder/tmp/s1.fq
    sfolder/tmp/s2.fq   (if present)
  Record SRA IDs to:
    sfolder/tmp/s.sra.txt
    sfolder/v< tier >/< inum >/summary-data/s.sra.txt
  Emit seqkit stats (-Ta) to:
    sfolder/v< tier >/< inum >/summary-data/s1.fq.seqkit.stats.ta.txt
    sfolder/v< tier >/< inum >/summary-data/s2.fq.seqkit.stats.ta.txt (if present)

Options:
  -s, --species   Species folder (e.g., Vigna_radiata)
  -i, --inum      Instance index (default: 0)
  --tier          Tier string (default: v6)
  --threads N     Threads for seqkit (default: 4)
  --data-short-redo BOOL
  --data-short-cleanup BOOL

Notes:
  Short-read SRA ID is taken from mapping: _short["<sfolder>-<inum>"].
EOF
	)

	local _data_long_cleanup=false
	local _data_long_redo=false
	local _remote=""

	parse_commandline() {
		set -- "${_brg_args[@]}"

		# source "${_POLAPLIB_DIR}/polap-cmd-version.sh" # '.' means 'source'
		while test $# -gt 0; do
			_key="$1"
			case "$_key" in
			--cleanup)
				_data_long_cleanup="true"
				;;
			--redo)
				_data_long_redo="true"
				;;
			--remote)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'."
				else
					_remote="$2"
					shift || true
				fi
				;;
			-*)
				_log_echo1 "[INFO] no such options: $1"
				;;
			*)
				:
				;;
			esac
			shift || true
		done
	}

	if [[ "${_brg_help:-off}" == "on" ]]; then
		declare -n ref="help_message"
		local manfile=$(_bolap_lib_man-convert_help_message "$ref" "${bolap_cmd}")
		man "$manfile"
		rm -f "$manfile"
		return 0
	fi

	# _log_echo0 "outdir: $_brg_outdir"
	# _log_echo0 "sindex: $_brg_sindex"
	parse_commandline

	# _log_echo0 "outdir: $_brg_outdir"
	# _log_echo0 "sindex: $_brg_sindex"

	# Load variables (expects _short map, _media_dir*, _local_host, etc.)
	source "${_POLAPLIB_DIR}/polap-variables-data.sh"

	if [[ -z "${_brg_outdir}" ]]; then
		_log_echo0 "[ERROR] --species is required"
		return 2
	fi

	local short_sra="${_short["$_brg_target"]:-}"

	if [[ -z "$short_sra" ]]; then
		_log_echo0 "[ERROR] No short-read SRA for ${_brg_target}"
		return 3
	fi

	# Build args for the worker script
	# local -a args=(-s "${short_sra}" -o "${_brg_outdir}" -i "${_brg_sindex}" --tier "${_tier}" --threads "${_threads}")
	local -a args=(-s "${short_sra}" -o "${_brg_outdir}")

	if [[ -n "${_remote}" ]]; then
		args+=(--remote ${_remote})
	fi

	_polap_lib_conda-ensure_conda_env polap-ncbitools || return 1
	bash "${_POLAPLIB_DIR}/polap-bash-run-data-short.sh" "${args[@]}"
	conda deactivate
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
				_log_echo "extract ${short_sra}"
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
				rm -rf "${short_sra}"
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
	elif [[ -s "${_brg_outdir}/${opt_t_arg}/0/short_expected_genome_size.txt" ]]; then
		_genome_size=$(<"${_brg_outdir}/${opt_t_arg}/0/short_expected_genome_size.txt")
	else

		"${_polap_cmd}" find-genome-size-for-pacbio -o "${_brg_outdir}" -l "${_brg_long}"

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

estimate-genomesize-platform_genus_species() {
	local _brg_outdir="${1:-all}"
	local _brg_inum="${2:-0}"

	local target_index="${_brg_outdir}-${_brg_inum}"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local platform="${_platform["$_brg_target"]}"

	local _genome_size="0"
	if [[ -s "${_brg_outdir}"/o/short_expected_genome_size.txt ]]; then
		_genome_size=$(<"${_brg_outdir}"/o/short_expected_genome_size.txt)
	elif [[ -s "${_brg_outdir}"/short_expected_genome_size.txt ]]; then
		_genome_size=$(<"${_brg_outdir}"/short_expected_genome_size.txt)
	elif [[ -s "${_brg_outdir}/${opt_t_arg}/0/short_expected_genome_size.txt" ]]; then
		_genome_size=$(<"${_brg_outdir}/${opt_t_arg}/0/short_expected_genome_size.txt")
	else
		if [[ "${platform}" == "PACBIO_SMRT" ]]; then
			"${_polap_cmd}" find-genome-size-for-pacbio \
				-o "${_brg_outdir}" \
				-l "${long_sra}.fastq"
		elif [[ "${platform}" == "ONT" ]]; then
			${_polap_cmd} find-genome-size \
				-a ${short_sra}_1.fastq \
				-b ${short_sra}_2.fastq \
				-o "${_brg_outdir}"
		fi
		if [[ -s "${_brg_outdir}"/short_expected_genome_size.txt ]]; then
			_genome_size=$(<"${_brg_outdir}"/short_expected_genome_size.txt)
		fi
	fi
	_genome_size=${_genome_size%%.*}
	echo "${_genome_size}"
}

write-genomesize-10x() {
	local _brg_outdir="${1}"
	local _brg_inum="${2:-0}"

	local target_index="${_brg_outdir}-${_brg_inum}"
	local long_sra="${_long["$target_index"]}"

	_polap_lib_conda-ensure_conda_env polap || exit 1

	local l=$(seqkit stats -Ta "${long_sra}.fastq" | csvtk cut -t -f sum_len | csvtk del-header)
	local _genome_size=$((l / 10))
	echo "${_genome_size}" >"${_brg_outdir}/short_expected_genome_size.txt"

	conda deactivate
}

run-downsample-long-data_genus_species() {
	local _brg_outdir="$1"
	local _brg_sindex="$2"
	source "${_POLAPLIB_DIR}/polap-variables-data.sh"
	shift 2

	# defaults
	local _brg_coverage="0.1"
	local _brg_random_seed="-1"
	local _brg_dry="off"

	# parse options
	while [[ $# -gt 0 ]]; do
		case "$1" in
		-c)
			_brg_coverage="$2"
			shift 2
			;;
		-s)
			_brg_random_seed="$2"
			shift 2
			;;
		--dry-run)
			_brg_dry="on"
			shift
			;;
		--)
			shift
			break
			;;
		*)
			echo "Unknown option: $1" >&2
			return 1
			;;
		esac
	done

	# load variables
	source "${_POLAPLIB_DIR}/polap-variables-data.sh"
	mkdir -p "${_brg_rundir}"

	# Step 1. check input fastq data
	if [[ ! -s "${_brg_tmpdir}/l.fq" ]]; then
		return 1
	fi

	_polap_lib_conda-ensure_conda_env polap || exit 1

	if [[ "$CONDA_DEFAULT_ENV" == "polap" ]]; then

		# echo "Genome size estimate: $(_polap_lib_unit-convert_bp ${_genome_size})"

		# coverage: 0.1x -> genome size 0.1
		# coverage: 0.1g -> 0.1 Gb
		if [[ "${_brg_coverage}" =~ [gGkKmM]$ ]]; then
			${_polap_cmd} fastq-sample-to -v \
				-l "${_brg_tmpdir}/l.fq" \
				--outfile "${_brg_tmpdir}/ld.fq" \
				-g "${_brg_coverage}" \
				--random-seed "${random_seed}" \
				>"${_brg_rundir}/l-${_brg_coverage}.txt"
		else

			# Step 2. Estimate the genome size
			local _genome_size="0"
			if [[ -v _platform["$_brg_target"] ]]; then
				_genome_size=$(
					estimate-genomesize-platform_genus_species \
						"${_brg_outdir}" \
						"${_brg_sindex}" | tail -1
				)
			else
				_genome_size=$(
					estimate-genomesize_genus_species \
						"${_brg_outdir}" \
						"${_brg_sindex}" | tail -1
				)
			fi
			if [[ ! "$_genome_size" =~ ^[0-9]+$ ]]; then
				echo "$_genome_size is NOT an integer"
				return 2
			fi

			${_polap_cmd} fastq subsample --redo -v \
				"${_brg_tmpdir}/l.fq" \
				"${_brg_tmpdir}/ld.fq" \
				-c "${_brg_coverage}" \
				-o "${_brg_rundir}" \
				--random-seed "${random_seed}" \
				--genomesize "${_genome_size}" \
				>"${_brg_rundir}/l-${_brg_coverage}x.txt"
		fi
	fi

	conda deactivate
}

run-downsample_genus_species() {
	run-downsample-long_genus_species
	run-downsample-short_genus_species
}

run-downsample-long_genus_species() {
	data-downsample-long_genus_species
}

data-downsample-long_genus_species() {
	# Debug trace with timestamps, file:line, and function
	if [[ "${POLAP_DEBUG:-0}" -eq 1 ]]; then
		export PS4='+$(date "+%F %T") ${BASH_SOURCE##*/}:${LINENO}:${FUNCNAME[0]}: '
		set -o errtrace  # inherit ERR in functions
		set -o functrace # inherit DEBUG/RETURN traps
		set -x           # trace execution
		trap 'st=$?; echo "ERR ${BASH_SOURCE[0]}:${LINENO}: ${BASH_COMMAND}" >&2; exit $st' ERR
	fi

	local _brg_outdir="${1:-$_brg_outdir}"
	local _brg_sindex="${2:-$_brg_sindex}"

	source "${_POLAPLIB_DIR}/polap-variables-data.sh"
	local long_sra="${_long["$_brg_target"]}"

	local bolap_cmd="${FUNCNAME%%_*}"

	help_message=$(
		cat <<EOF
Name:
  bolap - downsample long-read data 

Synopsis:
  bolap $bolap_cmd

Description:
  bolap

  --coverage STR [default: 30g]

Examples:
  Downsample the fastq to 10g:
    bolap data downsample-long -s Vigna_radiata --coverage 10g

  No downsample the fastq:
    bolap data downsample-long -s Vigna_radiata --coverage 0g

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	local _brg_random_seed="-1"
	local _brg_coverage="30g"
	local _brg_dry="off"

	local long_sra="${_long["$_brg_target"]}"
	local _brg_random_seed="${_random_seed["$_brg_target"]}"

	if [[ "${_brg_random_seed}" != "-1" ]]; then
		random_seed="${_brg_random_seed}"
	fi

	local _dry=""
	if [[ "${_brg_dry}" == "on" ]]; then
		_dry="--dry"
	fi

	parse_commandline() {
		set -- "${_brg_args[@]}"
		# set -- "${_brg_unknown_opts[@]}"

		# source "${_POLAPLIB_DIR}/polap-cmd-version.sh" # '.' means 'source'
		while test $# -gt 0; do
			_key="$1"
			case "$_key" in
			--coverage)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_coverage="$2"
					shift || true
				fi
				;;
			--random-seed)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_random_seed="$2"
					shift || true
				fi
				;;
			esac
			shift || true
		done
	}

	_polap_lib_help-maybe-show "$bolap_cmd" help_message || return 0

	parse_commandline

	_log_echo0 "Downsample long-read: ${_brg_coverage}"

	mkdir -p "${_brg_outdir_i}"
	mkdir -p "${_brg_target}"

	local l_fq_gz="${_brg_outdir}/tmp/l.fq.gz"
	local s_1_fq_gz="${_brg_outdir}/tmp/s_1.fq.gz"
	local s_2_fq_gz="${_brg_outdir}/tmp/s_2.fq.gz"

	# check if the files exist.
	if ! [[ -e "$l_fq_gz" && -s "$l_fq_gz" ]]; then
		echo "$l_fq_gz missing, empty, or broken link"
		echo "run data-long"
		return 0
	fi

	# _polap_lib_conda-ensure_conda_env polap || exit 1

	_log_echo0 _polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"

	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"

	if [[ "${_brg_coverage}" =~ [gGkKmM]$ ]]; then
		_log_echo1 ${_polap_cmd} fastq-sample-to -v \
			-l "${l_fq_gz}" \
			--outfile "${long_sra}.fastq" \
			-g "${_brg_coverage}" \
			--random-seed "${random_seed}"

		${_polap_cmd} fastq-sample-to -v \
			-l "${l_fq_gz}" \
			--outfile "${long_sra}.fastq" \
			-g "${_brg_coverage}" \
			--random-seed "${random_seed}" \
			>"${_brg_outdir_i}/l-${_brg_coverage}.txt"

		echo ${_polap_cmd} fastq-sample-to >"${_brg_target}/${long_sra}.fq.downsample.txt"
	else
		${_polap_cmd} fastq subsample --redo -v \
			"${l_fq_gz}" \
			"${long_sra}.fastq" \
			-c "${_brg_coverage}" \
			-o "${_brg_outdir_i}" \
			--random-seed "${random_seed}" \
			--genomesize "${_genome_size}" \
			>"${_brg_outdir_i}/l-${_brg_coverage}x.txt"
		echo ${_polap_cmd} fastq subsample >"${_brg_target}/${long_sra}.fq.downsample.txt"
	fi

	# seqkit stats -Ta
	echo "input: ${l_fq_gz}" >>"${_brg_target}/${long_sra}.fq.downsample.txt"
	echo "coverage: ${_brg_coverage}" >>"${_brg_target}/${long_sra}.fq.downsample.txt"

	_polap_lib_conda-ensure_conda_env polap || exit 1
	seqkit stats -Ta "${long_sra}.fastq" \
		-o "${_brg_target}/${long_sra}.fastq.seqkit.stats.ta.tsv"
	conda deactivate

	_log_echo0 "SRA: ${long_sra}"

	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"

	rsync -azuq --max-size=5M \
		"${_brg_target}/" "${_brg_rundir}/"

	if [[ -d "${_brg_target}" ]]; then
		_log_echo1 rm -rf "${_brg_target}"
		rm -rf "${_brg_target}"
	fi
}

run-downsample-short_genus_species() {
	data-downsample-short_genus_species
}

# Use the coverage 10g or 10 Gb.
#
# outdir/tmp/s_1.fq.gz
# outdir/tmp/s_2.fq.gz
#
# -> SRAID_1.fastq
# -> SRAID_2.fastq
data-downsample-short_genus_species() {
	local _brg_outdir="${1:-$_brg_outdir}"
	local _brg_sindex="${2:-$_brg_sindex}"

	source "${_POLAPLIB_DIR}/polap-variables-data.sh"
	local short_sra="${_short["$_brg_target"]}"

	local bolap_cmd="${FUNCNAME%%_*}"

	help_message=$(
		cat <<EOF
Name:
  bolap - downsample short-read data 

Synopsis:
  bolap $bolap_cmd

Description:
  bolap

Examples:
  Genome size:
    bolap data downsample-long -s Vigna_radiata --coverage 10g

  Downsample the fastq to 10g:
    bolap data downsample-long -s Vigna_radiata --coverage 10g

  No downsample the fastq to 10g:
    bolap data downsample-long -s Vigna_radiata --coverage 0

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	local _brg_random_seed="-1"
	local _brg_coverage="30g"
	local _brg_dry="off"

	local long_sra="${_long["$_brg_target"]}"
	local _brg_random_seed="${_random_seed["$_brg_target"]}"

	if [[ "${_brg_random_seed}" != "-1" ]]; then
		random_seed="${_brg_random_seed}"
	fi

	local _dry=""
	if [[ "${_brg_dry}" == "on" ]]; then
		_dry="--dry"
	fi

	parse_commandline() {
		set -- "${_brg_args[@]}"

		# source "${_POLAPLIB_DIR}/polap-cmd-version.sh" # '.' means 'source'
		while test $# -gt 0; do
			_key="$1"
			case "$_key" in
			--coverage)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_coverage="$2"
					shift || true
				fi
				;;
			esac
			shift || true
		done
	}

	_polap_lib_help-maybe-show "$bolap_cmd" help_message || return 0

	parse_commandline

	_log_echo0 "Downsample short-read: ${_brg_coverage}"

	mkdir -p "${_brg_outdir_i}"
	mkdir -p "${_brg_target}"

	local l_fq_gz="${_brg_outdir}/tmp/l.fq.gz"
	local s_1_fq_gz="${_brg_outdir}/tmp/s_1.fq.gz"
	local s_2_fq_gz="${_brg_outdir}/tmp/s_2.fq.gz"

	# check if the files exist.
	if ! [[ -e "$s_1_fq_gz" && -s "$s_1_fq_gz" ]]; then
		echo "$s_1_fq_gz missing, empty, or broken link"
		echo "run data-short"
		return 0
	fi

	_polap_lib_process-start_memtracker "${_memlog_file}" \
		"${_polap_var_memtracker_time_interval}"

	_log_echo1 "Sampling paired-end data"
	_log_echo1 ${_polap_cmd} fastq sample2-to \
		"${_brg_outdir}/tmp/s_1.fq.gz" \
		"${_brg_outdir}/tmp/s_2.fq.gz" \
		"${short_sra}_1.fastq" \
		"${short_sra}_2.fastq" \
		-g "${_brg_coverage}" \
		-o "${_brg_target}" \
		--random-seed "${random_seed}" \
		">${_brg_outdir_i}/s${_brg_coverage}x.txt"
	${_polap_cmd} fastq sample2-to \
		"${_brg_outdir}/tmp/s_1.fq.gz" \
		"${_brg_outdir}/tmp/s_2.fq.gz" \
		"${short_sra}_1.fastq" \
		"${short_sra}_2.fastq" \
		-g "${_brg_coverage}" \
		-o "${_brg_target}" \
		--random-seed "${random_seed}" \
		>"${_brg_outdir_i}/s${_brg_coverage}x.txt"

	echo "log: ${_brg_outdir_i}/s${_brg_coverage}x.txt"

	# seqkit stats -Ta
	echo "input: ${s_1_fq_gz}" >>"${_brg_target}/${short_sra}.fq.downsample.txt"
	echo "coverage: ${_brg_coverage}" >>"${_brg_target}/${short_sra}.fq.downsample.txt"

	_polap_lib_conda-ensure_conda_env polap || exit 1
	seqkit stats -Ta "${short_sra}_1.fastq" \
		-o "${_brg_target}/${short_sra}_1.fastq.seqkit.stats.ta.tsv"
	seqkit stats -Ta "${short_sra}_2.fastq" \
		-o "${_brg_target}/${short_sra}_2.fastq.seqkit.stats.ta.tsv"
	conda deactivate

	_polap_lib_process-end_memtracker "${_memlog_file}" "${_summary_file}" "no_verbose"

	rsync -azuq --max-size=5M \
		"${_brg_target}/" "${_brg_rundir}/"

	if [[ -d "${_brg_target}" ]]; then
		_log_echo1 rm -rf "${_brg_target}"
		rm -rf "${_brg_target}"
	fi
}

download-polap-github_genus_species() {
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
}

patch-polap_genus_species() {
	local bolap_cmd="${FUNCNAME%%_*}"

	help_message=$(
		cat <<EOF
Name:
  bolap - $bolap_cmd

Synopsis:
  bolap $bolap_cmd

Description:
  Patch the installed polap in conda environment with a github source version.

Examples:
  Patch polap:
    bolap patch-polap

  Patch polap of version 0.4:
    bolap install polap --ver 0.4

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	# defaults
	local _brg_ver=""

	_polap_lib_help-maybe-show "$bolap_cmd" help_message || return 0

	polap_setup_parse_commandline

	local _version="${_brg_ver}"

	if [[ "${opt_y_flag-}" == false ]]; then
		read -p "Do you want to replace conda env polap with the version ${_version}? (y/N): " confirm
	else
		confirm="yes"
	fi

	if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then

		_polap_lib_conda-ensure_conda_env polap || exit 1

		if [[ "$CONDA_DEFAULT_ENV" == "polap" ]]; then

			echo "[INFO] downloading polap v${_version} ..."
			rm -f "${_version}.zip"
			wget -q https://github.com/goshng/polap/archive/refs/tags/${_version}.zip
			if [[ -s "${_version}.zip" ]]; then
				rm -rf "polap-${_version}"
				unzip -o -q "${_version}.zip"
				cd "polap-${_version}/src"
				bash polaplib/polap-build.sh >../build.sh
				cd ..
				echo "[INFO] updating polap with v${_version} ..."
				PREFIX="$(conda info --base)/envs/polap" bash build.sh
				echo "[INFO] polap is patched with v${_version}"
			else
				echo "Error: no such polap version ${_version}"
				echo "Suggestion: $0 $bolap_cmd --ver 0.3.7.3"
				echo "Suggestion: $0 $bolap_cmd --ver 0.4.3.7.9"
				echo "Suggestion: $0 $bolap_cmd --ver 0.5.1.1"
			fi

		else
			echo "Error: You do not have polap environment."
			echo "Error: You're in the '$CONDA_DEFAULT_ENV' environment."
			exit 1
		fi

	else
		echo "polap patch is canceled."
	fi
}

bleeding-edge-polap_genus_species() {
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

}

local-edge-polap_genus_species() {

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
}

test-polap_genus_species() {
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

}

install-bolap_genus_species() {

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
				echo "[INFO] Conda environment 'bolap' already exists."
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
}

install-novoplasty_genus_species() {
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
				echo "[INFO] Conda environment 'novoplasty' already exists."
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
}

mkdir-all_genus_species() {
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
}

rm-empty_genus_species() {
	find . -type d -empty -delete
	echo "Deleting empty folders ... done."
}

rm_genus_species() {
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
}

xxx_genus_species() {
	echo "bolap subcommand xxx"
}

test_genus_species_for() {
	local _brg_outdir="${1:-all}"
	local _brg_sindex="${2:-0}"
	source "${_POLAPLIB_DIR}/polap-variables-data.sh"

	# local short_sra="${_short["$_brg_target"]}"
	${_polap_cmd} test log

	echo "Key: $_brg_target"
	if [[ "${_POLAP_RELEASE}" == "0" ]]; then
		echo "  long_sra: ${long_sra}"
	fi
}

test_genus_species() {
	local _brg_outdir="${1:-some}"
	local _brg_sindex="${2:-0}"
	source "${_POLAPLIB_DIR}/polap-variables-data.sh"

	if [[ "${_brg_outdir}" == "all" ]]; then
		for _v1 in "${Sall[@]}"; do
			test_genus_species_for "${_v1}" "${@:2}"
		done
	elif [[ "${_brg_outdir}" == "some" ]]; then
		for _v1 in "${Ssome[@]}"; do
			test_genus_species_for "${_v1}" "${_brg_sindex}"
		done
	elif [[ "${_brg_outdir}" == "test" ]]; then
		for _v1 in "${Stest[@]}"; do
			test_genus_species_for "${_v1}" "${_brg_sindex}"
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

################################################################################
# BEGIN: common manuscript functions
#
man_genus_species() {
	local first_arg="${1:-help}"
	local remaining_args=("${@:2}")

	man-${first_arg}_genus_species "${remaining_args[@]}"
}

man-figure-assembly_genus_species() {
	local args=("$@")

	_polap_lib_conda-ensure_conda_env polap || exit 1
	make -f "${_POLAPLIB_DIR}/../Makefile" manifests
	make -f "${_POLAPLIB_DIR}/../Makefile" figure-assemblies-two-pages
	conda deactivate
}

man-minimap2_genus_species() {
	local args=("$@")

	if [[ ! -s minimap2.1 ]]; then
		curl -L https://raw.githubusercontent.com/lh3/minimap2/master/minimap2.1 \
			-o minimap2.1
	fi
	man ./minimap2.1
}

man-update_genus_species() {
	local args=("$@")
	local first_arg="${1:-na}"
	local remaining_args=("${@:2}")

	if declare -F man-update-${first_arg}_genus_species >/dev/null; then
		man-update-${first_arg}_genus_species "${remaining_args[@]}"
	else
		echo "$help_message_man_update"
	fi
	return
}

man-table_genus_species() {
	local args=("$@")
	local first_arg="${1:-na}"
	local remaining_args=("${@:2}")

	case "${first_arg}" in
	benchmark)
		man-table-${first_arg}_genus_species "${remaining_args[@]}"
		;;
	dataset)
		man-table-${first_arg}_genus_species "${remaining_args[@]}"
		;;
	*)
		echo "$help_message_man_table"
		;;
	esac
}

#
# END: common manuscript functions
################################################################################

clean_input_lines() {
	while IFS= read -r line; do
		# Remove leading spaces
		line="${line#"${line%%[![:space:]]*}"}"

		# Remove trailing ')'
		line="${line%)}"

		echo "$line"
	done
}

dev_genus_species() {
	local bolap_cmd="${FUNCNAME%%_*}"

	help_message_main=$(
		cat <<EOF
Name:
  bolap - development

Synopsis:
  bolap $bolap_cmd

Description:
  Edit bolap-parsing.sh for bolap -h

  Edit polap-data-read.sh for bolap help

  Edit polap-lib-data.sh to add a new command

  Makefile.read to add man commands

  rsync -a -v thorne:/home/goshng/all/polap/github/src/ /home/goshng/all/polap/github/src/


Examples:
  Setup:
    mkdir -p $HOME/all/polap
    cd $HOME/all/polap
    git clone https://github.com/goshng/polap.git github
    alias pl='bash polap/src/polap.sh'
    alias bl='bash polap/src/bolap.sh'
    cd $HOME/all/polap/hifi1
    ln -s $HOME/all/polap/github polap
    cd .. && rsync -aq thorne:\$PWD/github/ github/ && rsync -aq thorne:.polap/ \$HOME/.polap/ && cd -

See also:
  Makefile.read

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	local _brg_topic="main"

	parse_commandline() {
		set -- "${_brg_args[@]}"

		while test $# -gt 0; do
			_key="$1"
			case "$_key" in
			--topic)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_brg_topic="$2"
					shift || true
				fi
				;;
			-*)
				_log_echo0 "[ERROR] no such options: $1"
				;;
			esac
			shift || true
		done
	}

	_polap_lib_help-maybe-show "$bolap_cmd" help_message || return

	parse_commandline

	# echo "$ref"
	local varname="help_message_${_brg_topic}"

	if [[ -v $varname ]]; then
		declare -n ref="$varname"
		manfile=$(_bolap_lib_man-convert_help_message "$ref" "$_brg_topic")
		man "$manfile"
		rm -f "$manfile" >&2
	else
		echo "[WARN] No help message defined for topic $_brg_topic" >&2
	fi
}

use_genus_species() {
	local bolap_cmd="${FUNCNAME%%_*}"

	help_message=$(
		cat <<EOF
Name:
  bolap ${bolap_cmd} - choose polap's data analysis

Synopsis:
  bolap ${bolap_cmd} DATA

Description:
  bolap ${bolap_cmd} allows to choose polap's data analysis; DATA can be aflye, cflye, read.
Each DATA has its own benchmark and man bolap subcommand.

Examples:
  Choose read polap's data analysis.
    bolap ${bolap_cmd} --data read --folder \$HOME

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	local _brg_topic=""
	local data="read"
	local folder="${HOME}"

	parse_commandline() {
		# echo "${_brg_args[@]}" >&2
		set -- "${_brg_args[@]}"

		while test $# -gt 0; do
			_key="$1"
			case "$_key" in
			--folder)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					folder="$2"
					shift || true
				fi
				;;
			--data)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					data="$2"
					shift || true
				fi
				;;
			-*)
				_log_echo0 "[ERROR] no such options: $1"
				;;
			esac
			shift || true
		done
	}

	# echo "H: $_brg_help" >&2

	declare -n ref="help_message"
	if [[ "${_brg_help}" == "on" ]]; then
		echo "$_brg_help" >&2
		local manfile=$(_bolap_lib_man-convert_help_message "$ref" "${bolap_cmd}")
		man "$manfile"
		rm -f "$manfile"
		return
	fi

	# echo "H1: $data" >&2
	parse_commandline
	# echo "H2: $data" >&2

	_log_echo "${bolap_cmd}_genus_species" "${data}" "${folder}"
	echo "${data}" >"${folder}/.bolaprc"
	echo "bolap use data: ${data}"

	# echo "$ref: $_brg_topic" >&2
}

run-save-sequence-id_genus_species() {
	local _brg_outdir="${1:-$_brg_outdir}"
	local _brg_sindex="${2:-$_brg_sindex}"

	source "${_POLAPLIB_DIR}/polap-variables-data.sh"

	local bolap_cmd="${FUNCNAME%%_*}"

	help_message=$(
		cat <<EOF
Name:
  bolap ${bolap_cmd} - template

Synopsis:
  bolap ${bolap_cmd}

Description:
  Save sequence IDs in Genus_species/tmp/l.fq.gz 

Examples:
  Subtitle:
    bolap ${bolap_cmd} -s Genus_species

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	local _brg_topic=""
	local data="read"
	local folder="${HOME}"

	parse_commandline() {
		# echo "${_brg_args[@]}" >&2
		set -- "${_brg_args[@]}"

		while test $# -gt 0; do
			_key="$1"
			case "$_key" in
			--folder)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					folder="$2"
					shift || true
				fi
				;;
			--data)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					data="$2"
					shift || true
				fi
				;;
			-s)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'."
				else
					_brg_outdir="$2"
					shift || true
				fi
				;;
			-i)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'."
				else
					_brg_sindex="$2"
					shift || true
				fi
				;;
			esac
			shift || true
		done
	}

	_polap_lib_help-maybe-show "$bolap_cmd" help_message || return 0

	parse_commandline

	source "${_POLAPLIB_DIR}/polap-variables-data.sh"

	_polap_lib_conda-ensure_conda_env polap || exit 1

	if [[ -s "${_brg_tmpdir}/l.fq.gz" ]]; then
		local sra_id=$(cat "${_brg_tmpdir}/l.sra.txt")
		if [[ -s "${_brg_outdir_i}/${sra_id}.txt.gz" ]]; then
			_log_echo0 "[INFO] file exist: ${_brg_outdir_i}/${sra_id}.txt.gz"
			# cp "${_brg_outdir_i}/${sra_id}.txt.gz" reads-id/
			# echo "${_brg_outdir} ${sra_id}" >>check-reads-species.txt
		else
			seqkit seq -ni "${_brg_tmpdir}/l.fq.gz" \
				-o "${_brg_outdir_i}/${sra_id}.txt.gz"
		fi
	else
		_log_echo0 "[ERROR] no such file: ${_brg_tmpdir}/l.fq.gz"
	fi

	conda deactivate
}

# species list (HiFi/root)
bolap_get_species_list() {
	local SET="${1:-auto}"

	case "$SET" in
	auto)
		find . -maxdepth 1 -type d -printf '%f\n' |
			grep -E '^[A-Z][a-z]+_[a-z]+$' |
			sort -u
		;;
	some)
		cat <<'__SOME__'
Biancaea_sappan
__SOME__
		;;
	file:*)
		local f="${SET#file:}"
		[[ -s "$f" ]] || {
			echo "[ERR] list not found: $f" >&2
			exit 2
		}
		sed -e 's/^[[:space:]]\+//' -e 's/[[:space:]]\+$//' -e '/^$/d' "$f"
		;;
	*) echo "$SET" ;;
	esac
}

batch_genus_species() {

	local bolap_cmd="${FUNCNAME%%_*}"

	help_message=$(
		cat <<EOF
Name:
  bolap ${bolap_cmd} - template

Synopsis:
  bolap ${bolap_cmd}

Description:
  Batch process of a bolap subcommand.

Examples:
  Get all data:
    bolap ${bolap_cmd} --cmd run-data-long --type auto

  Get sequence IDs for all:
    bolap ${bolap_cmd} --cmd run-save-sequence-id --type some

  Export bolap analysis:
    bolap ${bolap_cmd} --cmd export --type some

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	local _batch_command=""
	local _batch_data_type="auto"

	parse_commandline() {
		set -- "${_brg_args[@]}"

		while test $# -gt 0; do
			_key="$1"
			case "$_key" in
			--cmd)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_batch_command="$2"
					shift || true
				fi
				;;
			--type)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					_batch_data_type="$2"
					shift || true
				fi
				;;
			esac
			shift || true
		done
	}

	_polap_lib_help-maybe-show "$bolap_cmd" help_message || return 0

	parse_commandline

	if [[ -z "$_batch_command" ]]; then
		_log_echo0 "[ERROR] missing --cmd"
		return 0
	fi

	# get all species folder names.
	# 1. all in the current folder
	# 2. species name text file
	while read -r sp; do
		[[ -z "${sp// /}" ]] && continue
		_brg_args=(-s $sp -i 0)
		if [[ "${_batch_command}" == "benchmark" ]]; then
			_brg_args=(-f 6h -y -s $sp -i 0)
		fi
		echo "${_batch_command}_genus_species" "${_brg_args[@]}" >&2
		"${_batch_command}_genus_species"
	done < <(bolap_get_species_list "${_batch_data_type}")
}

bolap_s_parse_commandline() {
	# echo "${_brg_args[@]}" >&2
	set -- "${_brg_args[@]}"

	while test $# -gt 0; do
		_key="$1"
		case "$_key" in
		-s)
			if test $# -lt 2; then
				_log_echo0 "[ERROR] Missing value for the optional argument '$_key'."
			else
				_brg_outdir="$2"
				shift || true
			fi
			;;
		-i)
			if test $# -lt 2; then
				_log_echo0 "[ERROR] Missing value for the optional argument '$_key'."
			else
				_brg_sindex="$2"
				shift || true
			fi
			;;
		esac
		shift || true
	done
}

tar_gz_with_progress() {
	local src="$1"
	local out="$2"

	if [[ -z "$src" || -z "$out" ]]; then
		echo "Usage: tar_gz_with_progress <source_folder> <output.tar.gz>"
		return 1
	fi

	if ! command -v tar >/dev/null; then
		echo "Error: tar not found."
		return 1
	fi

	# Size in bytes for progress
	local size
	size=$(du -sb "$src" | awk '{print $1}')

	if command -v pv >/dev/null; then
		tar cf - "$src" |
			pv -s "$size" |
			gzip >"$out"
	else
		tar czf "$out" "$src"
	fi

	echo "Done: $out"
}

export_genus_species() {
	local bolap_cmd="${FUNCNAME%%_*}"

	help_message=$(
		cat <<EOF
Name:
  bolap ${bolap_cmd} - compress a spceies folder

Synopsis:
  bolap ${bolap_cmd}

Description:
  text.

Examples:
  Subtitle:
    bolap ${bolap_cmd} -s Genus_species

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	_polap_lib_help-maybe-show "$bolap_cmd" help_message || return 0

	bolap_s_parse_commandline

	source "${_POLAPLIB_DIR}/polap-variables-data.sh"

	# echo "s: ${_brg_outdir}" >&2
	# echo "i: ${_brg_sindex}" >&2

	tar_gz_with_progress "${_brg_outdir}" "${_brg_outdir}.tar.gz"
}

version_genus_species() {
	local bolap_cmd="${FUNCNAME%%_*}"

	help_message=$(
		cat <<EOF
Name:
  bolap ${bolap_cmd} - display versions of different tools

Synopsis:
  bolap ${bolap_cmd}

Description:
  Display the versions of some tools used.

Examples:
  Display versions:
    bolap ${bolap_cmd}

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	_polap_lib_help-maybe-show "$bolap_cmd" help_message || return 0

	_polap_lib_conda-ensure_conda_env polap-oatk || exit 1
	echo "Oatk $(oatk --version)"
	conda deactivate

	_polap_lib_conda-ensure_conda_env polap-tippo || exit 1
	echo "TIPPo $(TIPPo.v2.4.pl -v)"
	conda deactivate

	_polap_lib_conda-ensure_conda_env polap-himt || exit 1
	himt --version
	conda deactivate

	_polap_lib_conda-ensure_conda_env polap-pmat2 || exit 1
	local pmat_cmd=$(which PMAT)
	$pmat_cmd --version
	set +u
	conda deactivate
	set -u
}

run-data-long-10x_genus_species() {
	local bolap_cmd="${FUNCNAME%%_*}"

	help_message=$(
		cat <<EOF
Name:
  bolap ${bolap_cmd} - template

Synopsis:
  bolap ${bolap_cmd}

Description:
  Use read names in polap/data/hifi to prepare 10x input data.

Examples:
  Subtitle:
    bolap ${bolap_cmd} -s Genus_species

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	_polap_lib_help-maybe-show "$bolap_cmd" help_message || return 0

	bolap_s_parse_commandline

	source "${_POLAPLIB_DIR}/polap-variables-data.sh"

	_polap_lib_conda-ensure_conda_env polap || exit 1

	if [[ -s "polap/data/hifi/${long_sra}.txt.gz" ]]; then

		echo "[INFO] creating 10x ${long_sra} in ${_brg_tmpdir}/l.fq.gz" >&2
		gzip -dc "polap/data/hifi/${long_sra}.txt.gz" >"${long_sra}.txt"

		echo "[WARN] backing up ${_brg_tmpdir}/lb.fq.gz" >&2
		if [[ ! -e "${_brg_tmpdir}/lb.fq.gz" ]]; then
			mv "${_brg_tmpdir}/l.fq.gz" "${_brg_tmpdir}/lb.fq.gz"
		fi

		seqkit grep -f "${long_sra}.txt" \
			"${_brg_tmpdir}/lb.fq.gz" \
			-o "${long_sra}.fastq"

		gzip "${long_sra}.fastq"
		mv "${long_sra}.fastq.gz" "${_brg_tmpdir}/l.fq.gz"
		rm "${long_sra}.txt"
	else
		echo "[ERROR] no read name file: polap/data/hifi/${long_sra}.txt.gz" >&2
	fi

	_polap_lib_conda-ensure_conda_env_deactivate
}

bolap_cmd_parse_commandline() {
	set -- "${_brg_args[@]}"

	while test $# -gt 0; do
		_key="$1"
		case "$_key" in
		--cmd)
			if test $# -lt 2; then
				_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
			else
				cmd="$2"
				shift || true
			fi
			;;
		esac
		shift || true
	done
}

help-himt_genus_species() {
	local bolap_cmd="${FUNCNAME%%_*}"

	help_message=$(
		cat <<EOF
Name:
  bolap ${bolap_cmd}

Synopsis:
  bolap ${bolap_cmd} [--cmd assemble|assess|filter|compare|refassemble]

Description:
  Display the help message of HiMT.

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	_polap_lib_help-maybe-show "$bolap_cmd" help_message || return 0

	local cmd=""
	bolap_cmd_parse_commandline

	_polap_lib_conda-ensure_conda_env polap-himt || exit 1

	if ! himt $cmd -h 2>/dev/null; then
		echo "[ERROR] use --cmd with proper option value; use --help"
	fi
	# if [[ -n "$cmd" ]]; then
	# 	himt $cmd -h
	# else
	# 	himt -h
	# fi

	_polap_lib_conda-ensure_conda_env_deactivate
}

help-pmat2_genus_species() {
	local bolap_cmd="${FUNCNAME%%_*}"

	help_message=$(
		cat <<EOF
Name:
  bolap ${bolap_cmd}

Synopsis:
  bolap ${bolap_cmd} [--cmd autoMito|graphBuild]

Description:
  Display the help message of PMAT2.

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	_polap_lib_help-maybe-show "$bolap_cmd" help_message || return 0

	local cmd=""
	bolap_cmd_parse_commandline

	_polap_lib_conda-ensure_conda_env polap-pmat2 || exit 1

	local pmat_cmd=$(which PMAT)
	if ! "${pmat_cmd}" $cmd -h 2>/dev/null; then
		echo "[ERROR] use --cmd with proper option value; use --help"
	fi

	_polap_lib_conda-ensure_conda_env_deactivate
}

help-tippo_genus_species() {
	local bolap_cmd="${FUNCNAME%%_*}"

	help_message=$(
		cat <<EOF
Name:
  bolap ${bolap_cmd}

Synopsis:
  bolap ${bolap_cmd}

Description:
  Display the help message of TIPPo.

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	_polap_lib_help-maybe-show "$bolap_cmd" help_message || return 0

	_polap_lib_conda-ensure_conda_env polap-tippo || exit 1

	TIPPo.v2.4.pl -h

	_polap_lib_conda-ensure_conda_env_deactivate
}

help-oatk_genus_species() {
	local bolap_cmd="${FUNCNAME%%_*}"

	help_message=$(
		cat <<EOF
Name:
  bolap ${bolap_cmd}

Synopsis:
  bolap ${bolap_cmd}

Description:
  Display the help message of Oatk.

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	_polap_lib_help-maybe-show "$bolap_cmd" help_message || return 0

	_polap_lib_conda-ensure_conda_env polap-oatk || exit 1

	oatk -h

	_polap_lib_conda-ensure_conda_env_deactivate
}

# 2025-12-05
# improve this template bolap command.
# start with it for creating a new bolap subcommand.
template_genus_species() {
	local bolap_cmd="${FUNCNAME%%_*}"

	help_message=$(
		cat <<EOF
Name:
  bolap ${bolap_cmd} - template

Synopsis:
  bolap ${bolap_cmd}

Description:
  text.

Examples:
  Subtitle:
    bolap ${bolap_cmd} template

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	local _brg_topic=""
	local data="read"
	local folder="${HOME}"

	parse_commandline() {
		# echo "${_brg_args[@]}" >&2
		set -- "${_brg_args[@]}"

		while test $# -gt 0; do
			_key="$1"
			case "$_key" in
			--folder)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					folder="$2"
					shift || true
				fi
				;;
			--data)
				if test $# -lt 2; then
					_log_echo0 "[ERROR] Missing value for the optional argument '$_key'." 1
				else
					data="$2"
					shift || true
				fi
				;;
			esac
			shift || true
		done
	}

	_polap_lib_help-maybe-show "$bolap_cmd" help_message || return 0

	parse_commandline
}
