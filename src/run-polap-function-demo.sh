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
source "$script_dir/run-polap-function-include.sh"
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

source "$script_dir/run-polap-function-miscellaneous.sh"

# Common operations function
common_operations() {
	local long_sra="$1"
	local short_sra="$2"

	if [[ -s "o1/0/mt.contig.name-1" ]]; then
		mkdir -p o/0
		cp -p o2/*.txt o
		cp -p o1/0/mt.contig.name-1 o/0/
		cp -pr o1/0/30-contigger o/0/
	fi

	if [[ ! -s "${long_sra}.fastq" ]]; then
		_arg_sra="${long_sra}"
		_run_polap_x-ncbi-fetch-sra
		_arg_sra="${short_sra}"
		_run_polap_x-ncbi-fetch-sra
	fi
	cp -s "${long_sra}.fastq" l.fq
	cp -s "${short_sra}_1.fastq" s1.fq
	cp -s "${short_sra}_2.fastq" s2.fq
	if [[ -s "o/0/mt.contig.name-1" ]]; then
		_run_polap_reduce-data
	fi
}

# Function for each species
run_spirodela_polyrhiza() {
	cd Spirodela_polyrhiza || exit
	common_operations "SRR11472010" "SRR11472009"
	if [[ -s "o/0/mt.contig.name-1" ]]; then
		_run_polap_assemble2
	else
		_run_polap_assemble
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

function _run_polap_demo {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "$script_dir/polap-variables-common.sh" # '.' means 'source'

	local S=('Spirodela_polyrhiza' 'Taraxacum_mongolicum' 'Trifolium_pratense' 'Salix_dunnii' 'Anthoceros_agrestis' 'Anthoceros_angustus' 'Brassica_rapa' 'Vigna_radiata' 'Macadamia_tetraphylla' 'Punica_granatum' 'Lolium_perenne')

	local _m2=${_arg_menu[1]%/}

	if [[ -d "src" ]]; then
		_polap_cmd="src/polap.sh"
	else
		_polap_cmd="polap"
	fi
	local _polap_version="0.3.7.3"
	local _media_dir="/media/h1/run/mtdna"

	local -A _host
	_host['Spirodela_polyrhiza']="kishino"
	_host['Taraxacum_mongolicum']="lab01"
	_host['Trifolium_pratense']="vincent"
	_host['Salix_dunnii']="marybeth"
	_host['Anthoceros_agrestis']="kishino"
	_host['Anthoceros_angustus']="kishino"
	_host['Brassica_rapa']="lab01"
	_host['Vigna_radiata']="marybeth"
	_host['Macadamia_tetraphylla']="vincent"
	_host['Punica_granatum']="marybeth"
	_host['Lolium_perenne']="siepel"

	# Print help message if requested
	help_message=$(
		cat <<HEREDOC
# Polap demo
#
# <species_folder>: run polap on the species_folder
# install-conda: install Miniconda3
# setup-conda: setup Miniconda3 for Bioconda
# install-polap or install: install Polap
# install-fmlrc: install ptGAUL's FMLRC short-read polishing tools
# patch-polap: update the miniconda3/envs/polap/bin/polap with the github version
# Carex_pseudochinensis: test with mtDNA of Carex pseudochinensis
# test-polap: test Polap run with a test dataset
# download-polap: download Polap
# clean: delete analyses
# uninstall: uninstall Polap
# mkdir: create the 11 species folders.
# rm: delete the 11 species folders.
# link-fastq: create links to the input data at ${_media_dir} for the 11 folders.
# delete-links: delete all links in the current folder and its subfolders.
# copy-fastq: copy the input data at ${_media_dir} to the 11 folders.
# scopy-fastq <ssh-hostname>: transfer the input data at ${_media_dir} to the remote computers.
# zip: compress the 11 species folders.
# plot: create the supplementary figures for the -w option tests.
# table: create a table for the 11 test data.
# sync 0: copy the whole-genome assembly results from other computers to the local.
# sync 1: copy the organelle-genome assembly results from other computers to the local.
#
# Arguments:
#   -i ${_arg_inum}: source Flye (usually whole-genome) assembly number
#
# Inputs:
#   ${_polap_var_ga_annotation_all}
#
# Outputs:
#   ${_polap_var_mtcontigname}
#
# Menu2 (subcommands):
#   install-conda
#   setup-conda
#   install-polap
#   install-fmlrc
#   patch-polap or patch
#   test-polap or test
#   mkdir
#   rm
#   zip
#   table1
#
# Menu2:
#   Taraxacum_mongolicum
#   Trifolium_pratense
#   Salix_dunnii
#   Anthoceros_agrestis
#   Anthoceros_angustus
#   Brassica_rapa
#   Vigna_radiata
#   Macadamia_tetraphylla
#   Punica_granatum
#   Lolium_perenne)
Example: $(basename $0) ${_arg_menu[0]} Spirodela_polyrhiza
Example: $(basename $0) ${_arg_menu[0]} install-fmlrc
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	case "${_m2}" in
	"Spirodela_polyrhiza" | \
		"Taraxacum_mongolicum" | \
		"Trifolium_pratense" | \
		"Salix_dunnii" | \
		"Anthoceros_agrestis" | \
		"Anthoceros_angustus" | \
		"Brassica_rapa" | \
		"Vigna_radiata" | \
		"Punica_granatum" | \
		"Macadamia_tetraphylla" | \
		"Lolium_perenne" | \
		"Carex_pseudochinensis")
		bash ${script_dir}/polap-data-v1.sh "${_m2}"
		;;
	"install-conda")
		mkdir -p $HOME/miniconda3
		wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O $HOME/miniconda3/miniconda.sh
		bash $HOME/miniconda3/miniconda.sh -b -u -p $HOME/miniconda3
		rm $HOME/miniconda3/miniconda.sh
		_polap_log0 After installing Miniconda3, close and reopen your terminal application.
		;;
	"setup-conda")
		source $HOME/miniconda3/bin/activate
		conda config --add channels bioconda
		conda config --add channels conda-forge
		conda config --set channel_priority strict
		_polap_log0 Execute: source $HOME/miniconda3/bin/activate
		;;
	"install-polap")
		conda create --name polap polap
		;;
	"install-fmlrc" | "install")
		conda env create -f ${script_dir}/polap-conda-environment-fmlrc.yaml
		;;
	"git-clone")
		git clone https://github.com/goshng/polap.git
		;;
	"clean")
		rm -f ${_polap_version}.zip
		rm -rf polap-${_polap_version}
		;;
	"patch-polap" | "patch")
		rm -rf polap
		git clone -q https://github.com/goshng/polap.git
		cd polap/src
		bash polap-build.sh >../build.sh
		cd ..
		PREFIX="$HOME/miniconda3/envs/polap" bash build.sh
		;;
	"test-polap" | "test")
		rm -rf polap
		git clone -q https://github.com/goshng/polap.git
		cd polap/test
		source $HOME/miniconda3/bin/activate polap
		polap assemble --test
		;;
	"uninstall")
		_polap_log0 source $HOME/miniconda3/bin/deactivate
		_polap_log0 conda remove -n polap --all
		_polap_log0 conda remove -n polap-fmlrc --all
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
		pandoc ${script_dir}/doc/plot.md -o supp-figure-s11-s32.pdf
		pandoc ${script_dir}/doc/plot.md -o supp-figure-s11-s32.docx
		;;
	"table1")
		bash ${script_dir}/polap-report-table1.sh >table1.tsv
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
	*)
		echo "Menu2 $species_folder is not recognized."
		return 0
		;;
	esac

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}
