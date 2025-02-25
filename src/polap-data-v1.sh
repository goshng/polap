#!/bin/bash

S=(
	'Spirodela_polyrhiza'
	'Taraxacum_mongolicum'
	'Trifolium_pratense'
	'Salix_dunnii'
	'Anthoceros_agrestis'
	'Anthoceros_angustus' 'Brassica_rapa'
	'Vigna_radiata'
	'Macadamia_tetraphylla'
	'Punica_granatum'
	'Lolium_perenne'
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
	'ptgaul'
	'run'
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
# _polap_version="0.3.7.3"
_polap_version="0.4.3.7"
_media_dir="/media/h1/run/ptgaul20"
_media_dir="/media/h2/goshng/bioprojects"
_media_dir="/media/h1/run/mtdna"

# Input parameter
species_folder="$1"

help_message=$(
	cat <<HEREDOC
# Polap data analysis for the revision 1
#
# Edit two variables:
# _polap_cmd=${_polap_cmd}
# _media_dir=${_media_dir}
#
# Argument:
0. test
1. send-data-to
2. mkdir
3. taxon-ref1
4. ptgaul
5. run
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
if [[ -z "$1" ]]; then
	echo "Usage: $0 <subcommand> [species_folder]"
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

declare -A _long
declare -A _short
_long['Anthoceros_agrestis']="SRR10190639"
_short['Anthoceros_agrestis']="SRR10250248"
_long['Brassica_rapa']="ERR6210792"
_short['Brassica_rapa']="ERR6210790"
_long['Vigna_radiata']="SRR12549541"
_short['Vigna_radiata']="SRR12549533"
_long['Trifolium_pratense']="SRR15433794"
_short['Trifolium_pratense']="SRR15433795"
_long['Taraxacum_mongolicum']="SRR19182970"
_short['Taraxacum_mongolicum']="SRR19182971"
_long['Spirodela_polyrhiza']="SRR11472010"
_short['Spirodela_polyrhiza']="SRR11472009"
_long['Salix_dunnii']="SRR12893432"
_short['Salix_dunnii']="SRR12893433"
_long['Punica_granatum']="SRR24893686"
_short['Punica_granatum']="SRR24893685"
_long['Macadamia_tetraphylla']="SRR10424548"
_short['Macadamia_tetraphylla']="SRR10424549"
_long['Lolium_perenne']="SRR13386519"
_short['Lolium_perenne']="SRR13386518"
_long['Anthoceros_angustus']="SRR9696346"
_short['Anthoceros_angustus']="SRR9662965"
_long['Carex_pseudochinensis']="SRR30757341"
_short['Carex_pseudochinensis']="SRR30757340"

declare -A _inref
_inref['Anthoceros_agrestis']="genus"
_inref['Brassica_rapa']="genus"
_inref['Vigna_radiata']="genus"
_inref['Trifolium_pratense']="genus"
_inref['Taraxacum_mongolicum']="genus"
_inref['Spirodela_polyrhiza']="family"
_inref['Salix_dunnii']="genus"
_inref['Punica_granatum']="family"
_inref['Macadamia_tetraphylla']="genus"
_inref['Lolium_perenne']="family"
_inref['Anthoceros_angustus']="genus"
_inref['Carex_pseudochinensis']="genus"
declare -A _ingroup
_ingroup['Spirodela_polyrhiza']="class"
_ingroup['Taraxacum_mongolicum']="genus"
declare -A _outgroup
_outgroup['Spirodela_polyrhiza']="class"
_outgroup['Taraxacum_mongolicum']="family"
declare -A _allgroup
_allgroup['Spirodela_polyrhiza']="phylum"
_allgroup['Taraxacum_mongolicum']="family"

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

taxon-ref1_genus_species() {
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

	local average_length=$(seqkit stats -Ta ${output_dir}/ptgaul-reference.fa | awk 'NR==2 {print $7}')
	local genome_size=${average_length%.*}

	command time -v bash src/ptGAUL1.sh \
		-o ${output_dir}-ptgaul \
		-r ${output_dir}/ptgaul-reference.fa \
		-g "${genome_size}" \
		-l "${output_dir}/${long_sra}.fastq" \
		-t 24 \
		2>${output_dir}/timing-ptgaul.txt

	rm -rf ${output_dir}/result_3000
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

# Main case statement
case "$subcmd1" in
'send-data-to' | \
	'mkdir' | \
	'taxon-ref1' | \
	'ptgaul' | \
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
*)
	echo "Species '$species_folder' is not recognized."
	exit 1
	;;
esac
