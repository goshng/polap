#!/bin/bash

S=('Spirodela_polyrhiza' 'Taraxacum_mongolicum' 'Trifolium_pratense' 'Salix_dunnii' 'Anthoceros_agrestis' 'Anthoceros_angustus' 'Brassica_rapa' 'Vigna_radiata' 'Macadamia_tetraphylla' 'Punica_granatum' 'Lolium_perenne')

# Input parameter
species_folder="$1"

if [[ -d "src" ]]; then
	_polap_cmd="src/polap.sh"
else
	_polap_cmd="polap"
fi
_polap_version="0.3.7.3"
_media_dir="/media/h1/run/mtdna"

help_message=$(
	cat <<HEREDOC
# Polap data analysis for the revision 1
#
# Edit two variables:
# _polap_cmd=${_polap_cmd}
# _media_dir=${_media_dir}
#
# Argument:
<species_folder>: run polap on the species_folder
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
HEREDOC
)

# Check if the species folder is provided
if [[ -z "$species_folder" ]]; then
	echo "Usage: $0 <species_folder>"
	echo "       $0 <arg1> <arg2>"
	echo "${help_message}"
	exit 1
fi

value="${2:-1}"

declare -A _host
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
		${_polap_cmd} assemble2 --polap-reads -w 11000
	else
		${_polap_cmd} assemble --polap-reads -w 11000
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

# Main case statement
case "$species_folder" in
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
	unzip ${_polap_version}.zip
	cd polap-${_polap_version}
	;;
"clean")
	rm -f ${_polap_version}.zip
	rm -rf polap-${_polap_version}
	;;
"install-fmlrc")
	wget https://github.com/goshng/polap/archive/refs/tags/${_polap_version}.zip
	unzip ${_polap_version}.zip
	cd polap-${_polap_version}
	conda env create -f src/polap-conda-environment-fmlrc.yaml
	;;
"patch-polap")
	wget https://github.com/goshng/polap/archive/refs/tags/${_polap_version}.zip
	unzip ${_polap_version}.zip
	cd polap-${_polap_version}
	bash src/polap-build.sh >build.sh
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
		src/polap.sh package -o ../revision/$V1/o --archive $V1/o
		if [[ "$V1" == "Salix_dunnii" ]]; then
			src/polap.sh -o $V1/o test-reads report ptgaul
			src/polap.sh -o $V1/o test-reads report intra
		elif [[ "$V1" == "Lolium_perenne" ]]; then
			src/polap.sh -o $V1/o test-reads report ptgaul --report-x 5000,7000,9000,11000,13000,15000,17000
			src/polap.sh -o $V1/o test-reads report polap --report-x 5000,7000,9000,11000,13000,15000,17000
		else
			src/polap.sh -o $V1/o test-reads report ptgaul
			src/polap.sh -o $V1/o test-reads report polap
		fi
	done
	;;
*)
	echo "Species '$species_folder' is not recognized."
	exit 1
	;;
esac
