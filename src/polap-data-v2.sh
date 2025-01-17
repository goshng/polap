#!/usr/bin/bash

S=('Juncus_effusus' 'Juncus_inflexus' 'Juncus_roemerianus' 'Juncus_validus' 'Eucalyptus_pauciflora' 'Arctostaphylos_glauca' 'Lepidium_sativum' 'Chaetoceros_muellerii' 'Potentilla_micrantha' 'Durio_zibethinus' 'Beta_vulgaris' 'Eleocharis_dulcis' 'Leucanthemum_vulgare' 'Oryza_glaberrima' 'Cenchrus_americanus' 'Digitaria_exilis' 'Podococcus_acaulis' 'Raphia_textilis' 'Phytelephas_aequatorialis' 'Picea_glauca')

# Input parameter
species_folder="${1%/}"

if [[ -d "src" ]]; then
	_polap_cmd="src/polap.sh"
else
	_polap_cmd="polap"
fi
_polap_version="0.4.1.9"
_media_dir="/media/h1/run/ptgaul20"

help_message=$(
	cat <<HEREDOC
# Polap data analysis for the version 0.4 or disassemble (subsampling method)
#
# Edit two variables:
# _polap_cmd=${_polap_cmd}
# _media_dir=${_media_dir}
#
# Argument:
mkdir: create the 16 species folders.

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

# Function for each species
# kishino
run_Juncus_effusus() {
	local output_dir="$(echo $FUNCNAME | sed s/run_//)"
	local species_name="$(echo $FUNCNAME | sed 's/run_//' | sed 's/_/ /')"
	local long_sra="SRR14298760"
	local short_sra="SRR14298746"
	# copy_data

	local i=0
	local n
	local p
	for n in 10 30 100; do
		for p in 1 5 10; do
			i=$((i + 1))
			if [[ -d "${output_dir}/disassemble/${i}" ]]; then
				echo "exists: $i, $p, $n"
			else
				command time -v ${_polap_cmd} disassemble -o ${output_dir} \
					-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
					--disassemble-i $i \
					--disassemble-p $p \
					--disassemble-n $n 2>${output_dir}/timing-${i}.txt
			fi
		done
	done

	# command time -v ${_polap_cmd} disassemble -o ${output_dir}-a \
	# 	-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
	# 	--disassemble-a 100m \
	# 	--disassemble-b 1g \
	# 	--disassemble-n 100 \
	# 	--disassemble-best \
	# 	--disassemble-compare-to-fasta ptdna-${output_dir}-best.fa -v

	# Elapsed (wall clock) time (h:mm:ss or m:ss): 5:27:43
	# Maximum resident set size (kbytes): 32578112
	# rm -rf ${output_dir}/disassemble/x
	# command time -v ${_polap_cmd} disassemble -o ${output_dir} \
	# 	-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
	# 	--disassemble-compare-to-fasta ptdna-${output_dir}.fa \
	# 	--disassemble-s 181m --disassemble-alpha 3.25 \
	# 	--disassemble-n 30 \
	# 	--disassemble-stop-after assemble -v
}

run_Juncus_effusus-a() {
	local output_dir="$(echo $FUNCNAME | sed s/run_//)"
	local species_name="$(echo $FUNCNAME | sed 's/run_//' | sed 's/_/ /')"
	local long_sra="SRR14298760"
	local short_sra="SRR14298746"

	local i=0
	local n
	local p
	# for n in 10 30 100; do
	# 	for p in 1 5 10; do
	# --stages-include 1-6 \
	# --disassemble-compare-to-fasta ptdna-${output_dir}-ptgaul.fa \
	for n in 10; do
		for p in 1; do
			i=$((i + 1))
			if [[ -d "${output_dir}/disassemble/${i}" ]]; then
				${_polap_cmd} disassemble -o ${output_dir} \
					--disassemble-compare-to-fasta ptdna-Juncus_effusus-ptgaul.fa \
					--stages-include 1 \
					--disassemble-best \
					--disassemble-polish \
					-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
					--disassemble-i $i \
					--disassemble-p $p \
					--disassemble-n $n
			else
				echo "exists: $i, $p, $n"
			fi
		done
	done

}

# vincent
run_Juncus_inflexus() {
	local output_dir="$(echo $FUNCNAME | sed s/run_//)"
	local species_name="$(echo $FUNCNAME | sed 's/run_//' | sed 's/_/ /')"
	species_name="Juncus effusus"
	local long_sra="SRR14298751"
	local short_sra="SRR14298745"

	# copy_data
	# cp -s ptdna-Juncus_effusus.fa ptdna-${output_dir}.fa

	# Elapsed (wall clock) time (h:mm:ss or m:ss): 20:45:16
	# Maximum resident set size (kbytes): 41721844
	# rm -rf ${output_dir}/disassemble/0
	# command time -v ${_polap_cmd} disassemble -o ${output_dir} \
	# 	-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
	# 	--disassemble-p 5 \
	# 	--disassemble-n 100

	# step 2
	bash ptgaul/ptGAUL.sh -o ${output_dir}-ptgaul -r ptdna-${output_dir}.fa -g 180000 -l ${long_sra}.fastq -t 24
	cp -pr "${output_dir}-ptgaul/result_3000" "${output_dir}"

	# local i=0
	# local n
	# local p
	# for n in 10 30 100; do
	# 	for p in 1 5 10; do
	# 		i=$((i + 1))
	# 		if [[ -d "${output_dir}/disassemble/${i}" ]]; then
	# 			echo "exists: $i, $p, $n"
	# 		else
	# 			command time -v ${_polap_cmd} disassemble -o ${output_dir} \
	# 				-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
	# 				--disassemble-i $i \
	# 				--disassemble-p $p \
	# 				--disassemble-n $n 2>${output_dir}/timing-${i}.txt
	# 		fi
	# 	done
	# done

	# rm -rf ${output_dir}/disassemble/x
	# command time -v ${_polap_cmd} disassemble -o ${output_dir} \
	# 	-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
	# 	--disassemble-compare-to-fasta ptdna-${output_dir}.fa \
	# 	--disassemble-s 181m --disassemble-alpha 3.25 \
	# 	--disassemble-n 30 \
	# 	--disassemble-stop-after assemble -v
}

run_Juncus_inflexus-a() {
	local output_dir="$(echo $FUNCNAME | sed s/run_//)"
	local species_name="$(echo $FUNCNAME | sed 's/run_//' | sed 's/_/ /')"
	species_name="Juncus effusus"
	local long_sra="SRR14298751"
	local short_sra="SRR14298745"

	# copy_data
	# cp -s ptdna-Juncus_effusus.fa ptdna-${output_dir}.fa

	# Elapsed (wall clock) time (h:mm:ss or m:ss): 20:45:16
	# Maximum resident set size (kbytes): 41721844
	# rm -rf ${output_dir}/disassemble/0
	# command time -v ${_polap_cmd} disassemble -o ${output_dir} \
	# 	-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
	# 	--disassemble-p 5 \
	# 	--disassemble-n 100

	# rm -rf ${output_dir}/disassemble/x
	# command time -v ${_polap_cmd} disassemble -o ${output_dir} \
	# 	-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
	# 	--disassemble-compare-to-fasta ptdna-${output_dir}.fa \
	# 	--disassemble-s 181m --disassemble-alpha 3.25 \
	# 	--disassemble-n 30 \
	# 	--disassemble-stop-after assemble -v
}

# siepel
run_Juncus_roemerianus() {
	local output_dir="$(echo $FUNCNAME | sed s/run_//)"
	local species_name="$(echo $FUNCNAME | sed 's/run_//' | sed 's/_/ /')"
	local long_sra="SRR21976090"
	local short_sra="SRR21976092"
	copy_data

	# Elapsed (wall clock) time (h:mm:ss or m:ss): 9:20:00
	# Maximum resident set size (kbytes): 27453508
	# rm -rf ${output_dir}/disassemble/0
	# command time -v ${_polap_cmd} disassemble -o ${output_dir} \
	# 	-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
	# 	--disassemble-p 5 \
	# 	--disassemble-n 100

	# ${_polap_cmd} disassemble -o ${output_dir} \
	# 	-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
	# 	--disassemble-i 1 \
	# 	--stages-include 6

	local i=0
	local n
	local p
	# step 1
	# for n in 10 30 100; do
	# 	for p in 1 5 10; do
	# 		i=$((i + 1))
	# 		if [[ -d "${output_dir}/disassemble/${i}" ]]; then
	# 			echo "exists: $i, $p, $n"
	# 		else
	# 			command time -v ${_polap_cmd} disassemble -o ${output_dir} \
	# 				-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
	# 				--disassemble-i $i \
	# 				--disassemble-p $p \
	# 				--disassemble-n $n 2>${output_dir}/timing-${i}.txt
	# 		fi
	# 	done
	# done

	# step 2
	bash ptgaul/ptGAUL.sh -o ${output_dir}-ptgaul -r ptdna-${output_dir}.fa -g 200000 -l ${long_sra}.fastq -t 24
	cp -pr "${output_dir}-ptgaul/result_3000" "${output_dir}"

	# step 3
	# for n in 10 30 100; do
	# 	for p in 1 5 10; do
	# 		i=$((i + 1))
	# 		if [[ -d "${output_dir}/disassemble/${i}" ]]; then
	# 			${_polap_cmd} disassemble -o ${output_dir} \
	# 				--stages-include 1-6 \
	# 				--disassemble-best \
	# 				--disassemble-compare-to-fasta ptdna-${output_dir}-ptgaul.fa \
	# 				-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
	# 				--disassemble-i $i \
	# 				--disassemble-p $p \
	# 				--disassemble-n $n
	# 		else
	# 			echo "exists: $i, $p, $n"
	# 		fi
	# 	done
	# done

	# rm -rf ${output_dir}/disassemble/x
	# command time -v ${_polap_cmd} disassemble -o ${output_dir} \
	# 	-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
	# 	--disassemble-compare-to-fasta ptdna-${output_dir}.fa \
	# 	--disassemble-s 600m --disassemble-alpha 3.25 \
	# 	--disassemble-n 30 \
	# 	--disassemble-stop-after assemble \
	# 	-v
}

run_Juncus_roemerianus-a() {
	local output_dir="$(echo $FUNCNAME | sed s/run_//)"
	local species_name="$(echo $FUNCNAME | sed 's/run_//' | sed 's/_/ /')"
	local long_sra="SRR21976090"
	local short_sra="SRR21976092"
	copy_data

	# Elapsed (wall clock) time (h:mm:ss or m:ss): 9:20:00
	# Maximum resident set size (kbytes): 27453508
	# rm -rf ${output_dir}/disassemble/0
	# command time -v ${_polap_cmd} disassemble -o ${output_dir} \
	# 	-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
	# 	--disassemble-p 5 \
	# 	--disassemble-n 100

	# ${_polap_cmd} disassemble -o ${output_dir} \
	# 	-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
	# 	--disassemble-i 1 \
	# 	--stages-include 6

	local i=0
	local n
	local p
	# step 1
	# for n in 10 30 100; do
	# 	for p in 1 5 10; do
	# 		i=$((i + 1))
	# 		if [[ -d "${output_dir}/disassemble/${i}" ]]; then
	# 			echo "exists: $i, $p, $n"
	# 		else
	# 			command time -v ${_polap_cmd} disassemble -o ${output_dir} \
	# 				-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
	# 				--disassemble-i $i \
	# 				--disassemble-p $p \
	# 				--disassemble-n $n 2>${output_dir}/timing-${i}.txt
	# 		fi
	# 	done
	# done

	# step 2
	# bash ptgaul/ptGAUL.sh -o ${output_dir}-ptgaul -r ptdna-${output_dir}.fa -g 200000 -l ${long_sra}.fastq -t 24
	# cp -pr "${output_dir}-ptgaul/result_3000" "${output_dir}"

	# step 3
	for n in 10 30 100; do
		for p in 1 5 10; do
			i=$((i + 1))
			if [[ -d "${output_dir}/disassemble/${i}" ]]; then
				${_polap_cmd} disassemble -o ${output_dir} \
					--stages-include 1-6 \
					--disassemble-best \
					--disassemble-compare-to-fasta ptdna-${output_dir}-ptgaul.fa \
					-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
					--disassemble-i $i \
					--disassemble-p $p \
					--disassemble-n $n
			else
				echo "exists: $i, $p, $n"
			fi
		done
	done

	# rm -rf ${output_dir}/disassemble/x
	# command time -v ${_polap_cmd} disassemble -o ${output_dir} \
	# 	-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
	# 	--disassemble-compare-to-fasta ptdna-${output_dir}.fa \
	# 	--disassemble-s 600m --disassemble-alpha 3.25 \
	# 	--disassemble-n 30 \
	# 	--disassemble-stop-after assemble \
	# 	-v
}

# marybeth
run_Juncus_validus() {
	local output_dir="$(echo $FUNCNAME | sed s/run_//)"
	local species_name="$(echo $FUNCNAME | sed 's/run_//' | sed 's/_/ /')"
	local long_sra="SRR21976089"
	local short_sra="SRR21976091"
	# copy_data

	# Elapsed (wall clock) time (h:mm:ss or m:ss): 20:37:27
	# Maximum resident set size (kbytes): 21494240
	# rm -rf ${output_dir}/disassemble/0
	# command time -v ${_polap_cmd} disassemble -o ${output_dir} \
	# 	-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
	# 	--disassemble-p 5 \
	# 	--disassemble-n 100

	# rm -rf ${output_dir}/disassemble/1/3
	# ${_polap_cmd} disassemble -o ${output_dir} \
	# 	-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
	# 	--stages-include 6 \
	# 	--disassemble-i 1

	local i=0
	local n
	local p
	# step 1
	# for n in 10 30 100; do
	# 	for p in 1 5 10; do
	# 		i=$((i + 1))
	# 		if [[ -d "${output_dir}/disassemble/${i}" ]]; then
	# 			echo "exists: $i, $p, $n"
	# 		else
	# 			command time -v ${_polap_cmd} disassemble -o ${output_dir} \
	# 				-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
	# 				--disassemble-i $i \
	# 				--disassemble-p $p \
	# 				--disassemble-n $n 2>${output_dir}/timing-${i}.txt
	# 		fi
	# 	done
	# done

	# step 2
	# bash ptgaul/ptGAUL.sh -o ${output_dir}-ptgaul -r ptdna-${output_dir}.fa -l ${long_sra}.fastq -t 24
	# cp -pr "${output_dir}-ptgaul/result_3000" "${output_dir}"

	# step 3
	for n in 10 30 100; do
		for p in 1 5 10; do
			i=$((i + 1))
			if [[ -d "${output_dir}/disassemble/${i}" ]]; then
				${_polap_cmd} disassemble -o ${output_dir} \
					--stages-include 1-6 \
					--disassemble-best \
					--disassemble-compare-to-fasta ptdna-${output_dir}-ptgaul.fa \
					-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
					--disassemble-i $i \
					--disassemble-p $p \
					--disassemble-n $n
			else
				echo "exists: $i, $p, $n"
			fi
		done
	done

	# rm -rf ${output_dir}/disassemble/x
	# command time -v ${_polap_cmd} disassemble -o ${output_dir} \
	# 	-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
	# 	--disassemble-compare-to-fasta ptdna-${output_dir}.fa \
	# 	--disassemble-s 200m --disassemble-alpha 2 \
	# 	--disassemble-n 2 \
	# 	--disassemble-stop-after assemble \
	# 	-v
}

run_Juncus_validus-a() {
	local output_dir="$(echo $FUNCNAME | sed s/run_//)"
	local species_name="$(echo $FUNCNAME | sed 's/run_//' | sed 's/_/ /')"
	local long_sra="SRR21976089"
	local short_sra="SRR21976091"
	# copy_data

	# Elapsed (wall clock) time (h:mm:ss or m:ss): 20:37:27
	# Maximum resident set size (kbytes): 21494240
	# rm -rf ${output_dir}/disassemble/0
	# command time -v ${_polap_cmd} disassemble -o ${output_dir} \
	# 	-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
	# 	--disassemble-p 5 \
	# 	--disassemble-n 100

	# rm -rf ${output_dir}/disassemble/1/3
	# ${_polap_cmd} disassemble -o ${output_dir} \
	# 	-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
	# 	--stages-include 6 \
	# 	--disassemble-i 1

	local i=0
	local n
	local p
	# step 1
	# for n in 10 30 100; do
	# 	for p in 1 5 10; do
	# 		i=$((i + 1))
	# 		if [[ -d "${output_dir}/disassemble/${i}" ]]; then
	# 			echo "exists: $i, $p, $n"
	# 		else
	# 			command time -v ${_polap_cmd} disassemble -o ${output_dir} \
	# 				-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
	# 				--disassemble-i $i \
	# 				--disassemble-p $p \
	# 				--disassemble-n $n 2>${output_dir}/timing-${i}.txt
	# 		fi
	# 	done
	# done

	# step 2
	# bash ptgaul/ptGAUL.sh -o ${output_dir}-ptgaul -r ptdna-${output_dir}.fa -l ${long_sra}.fastq -t 24
	# cp -pr "${output_dir}-ptgaul/result_3000" "${output_dir}"

	# step 3
	for n in 10 30 100; do
		for p in 1 5 10; do
			i=$((i + 1))
			if [[ -d "${output_dir}/disassemble/${i}" ]]; then
				${_polap_cmd} disassemble -o ${output_dir} \
					--stages-include 1-6 \
					--disassemble-best \
					--disassemble-compare-to-fasta ptdna-${output_dir}-ptgaul.fa \
					-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
					--disassemble-i $i \
					--disassemble-p $p \
					--disassemble-n $n
			else
				echo "exists: $i, $p, $n"
			fi
		done
	done

	# rm -rf ${output_dir}/disassemble/x
	# command time -v ${_polap_cmd} disassemble -o ${output_dir} \
	# 	-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
	# 	--disassemble-compare-to-fasta ptdna-${output_dir}.fa \
	# 	--disassemble-s 200m --disassemble-alpha 2 \
	# 	--disassemble-n 2 \
	# 	--disassemble-stop-after assemble \
	# 	-v
}

# lab01
run_Eucalyptus_pauciflora() {
	local output_dir="$(echo $FUNCNAME | sed s/run_//)"
	local species_name="$(echo $FUNCNAME | sed 's/run_//' | sed 's/_/ /')"
	local long_sra="SRR7153095"
	local short_sra="SRR7161123"
	copy_data

	# command time -v ${_polap_cmd} disassemble -o ${output_dir} \
	# rm -rf ${output_dir}/disassemble/0
	# Elapsed (wall clock) time (h:mm:ss or m:ss): 10:33:34
	# Maximum resident set size (kbytes): 21504672
	# command time -v ${_polap_cmd} disassemble -o ${output_dir} \
	# 	-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
	# 	--disassemble-p 5 \
	# 	--disassemble-n 100

	# ${_polap_cmd} disassemble -o ${output_dir} \
	# 	-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
	# 	--stages-include 4-6 \
	# 	--disassemble-i 1

	# 2
	# Elapsed (wall clock) time (h:mm:ss or m:ss): 1:08:49
	# Maximum resident set size (kbytes): 21504364

	local i=0
	local n
	local p
	# step 1
	# for n in 10 30 100; do
	# 	for p in 1 5 10; do
	# 		i=$((i + 1))
	# 		if [[ -d "${output_dir}/disassemble/${i}" ]]; then
	# 			echo "exists: $i, $p, $n"
	# 		else
	# 			command time -v ${_polap_cmd} disassemble -o ${output_dir} \
	# 				-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
	# 				--disassemble-i $i \
	# 				--disassemble-p $p \
	# 				--disassemble-n $n 2>${output_dir}/timing-${i}.txt
	# 		fi
	# 	done
	# done

	# step 2
	# bash ptgaul/ptGAUL.sh -o ${output_dir}-ptgaul -r ptdna-${output_dir}.fa -l ${long_sra}.fastq -t 24
	# cp -pr "${output_dir}-ptgaul/result_3000" "${output_dir}"

	# command time -v ${_polap_cmd} disassemble -o ${output_dir} \
	# 	-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
	# 	--disassemble-i 2 \
	# 	--disassemble-p 5 \
	# 	--disassemble-n 10

	# Stage 1: v0.4.1.7
	# ${_polap_cmd} disassemble -o ${output_dir} \
	# 	-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
	# 	--disassemble-a 50m \
	# 	--disassemble-b 1g \
	# 	--disassemble-n 30 \
	# 	--disassemble-compare-to-fasta ptdna-${output_dir}.fa -v

	# rm -rf ${output_dir}/disassemble/x
	# command time -v ${_polap_cmd} disassemble -o ${output_dir} \
	# 	-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
	# 	--disassemble-compare-to-fasta ptdna-${output_dir}.fa \
	# 	--disassemble-s 82m --disassemble-alpha 0.25 \
	# 	--disassemble-n 10 \
	# 	--disassemble-stop-after assemble \
	# 	-v
}

run_Eucalyptus_pauciflora-a() {
	local output_dir="$(echo $FUNCNAME | sed s/run_//)"
	local species_name="$(echo $FUNCNAME | sed 's/run_//' | sed 's/_/ /')"
	local long_sra="SRR7153095"
	local short_sra="SRR7161123"
	copy_data

	# command time -v ${_polap_cmd} disassemble -o ${output_dir} \
	# rm -rf ${output_dir}/disassemble/0
	# Elapsed (wall clock) time (h:mm:ss or m:ss): 10:33:34
	# Maximum resident set size (kbytes): 21504672
	# command time -v ${_polap_cmd} disassemble -o ${output_dir} \
	# 	-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
	# 	--disassemble-p 5 \
	# 	--disassemble-n 100

	# ${_polap_cmd} disassemble -o ${output_dir} \
	# 	-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
	# 	--stages-include 4-6 \
	# 	--disassemble-i 1

	# 2
	# Elapsed (wall clock) time (h:mm:ss or m:ss): 1:08:49
	# Maximum resident set size (kbytes): 21504364

	local i=0
	local n
	local p
	# step 1
	# for n in 10 30 100; do
	# 	for p in 1 5 10; do
	# 		i=$((i + 1))
	# 		if [[ -d "${output_dir}/disassemble/${i}" ]]; then
	# 			echo "exists: $i, $p, $n"
	# 		else
	# 			command time -v ${_polap_cmd} disassemble -o ${output_dir} \
	# 				-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
	# 				--disassemble-i $i \
	# 				--disassemble-p $p \
	# 				--disassemble-n $n 2>${output_dir}/timing-${i}.txt
	# 		fi
	# 	done
	# done

	# step 2
	# bash ptgaul/ptGAUL.sh -o ${output_dir}-ptgaul -r ptdna-${output_dir}.fa -l ${long_sra}.fastq -t 24
	# cp -pr "${output_dir}-ptgaul/result_3000" "${output_dir}"

	for n in 10 30 100; do
		for p in 1 5 10; do
			i=$((i + 1))
			if [[ -d "${output_dir}/disassemble/${i}" ]]; then
				${_polap_cmd} disassemble -o ${output_dir} \
					--stages-include 1-6 \
					--disassemble-best \
					--disassemble-compare-to-fasta ptdna-${output_dir}-ptgaul.fa \
					-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
					--disassemble-i $i \
					--disassemble-p $p \
					--disassemble-n $n
			else
				echo "exists: $i, $p, $n"
			fi
		done
	done

	# command time -v ${_polap_cmd} disassemble -o ${output_dir} \
	# 	-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
	# 	--disassemble-i 2 \
	# 	--disassemble-p 5 \
	# 	--disassemble-n 10

	# Stage 1: v0.4.1.7
	# ${_polap_cmd} disassemble -o ${output_dir} \
	# 	-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
	# 	--disassemble-a 50m \
	# 	--disassemble-b 1g \
	# 	--disassemble-n 30 \
	# 	--disassemble-compare-to-fasta ptdna-${output_dir}.fa -v

	# rm -rf ${output_dir}/disassemble/x
	# command time -v ${_polap_cmd} disassemble -o ${output_dir} \
	# 	-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
	# 	--disassemble-compare-to-fasta ptdna-${output_dir}.fa \
	# 	--disassemble-s 82m --disassemble-alpha 0.25 \
	# 	--disassemble-n 10 \
	# 	--disassemble-stop-after assemble \
	# 	-v
}

run_test() {
	local output_dir="$(echo $FUNCNAME | sed s/run_//)"
	local species_name="$(echo $FUNCNAME | sed 's/run_//' | sed 's/_/ /')"
	local long_sra="SRR7153095"
	local short_sra="SRR7161123"
	local long_sra="test"
	local short_sra="test"

	# command time -v ${_polap_cmd} disassemble -o ${output_dir} \
	# rm -rf ${output_dir}/disassemble/0
	# i:0
	# command time -v ${_polap_cmd} disassemble -o ${output_dir} \
	# 	-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
	# 	--disassemble-r 3
	# i:1
	# command time -v ${_polap_cmd} disassemble -o ${output_dir} \

	# rm -rf ${output_dir}/disassemble/1
	# rm -rf ${output_dir}/disassemble/1/2
	${_polap_cmd} disassemble -o ${output_dir} \
		-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
		--disassemble-i 3 \
		--disassemble-p 90 \
		--disassemble-n 5 \
		--disassemble-r 3

	# rm -rf ${output_dir}/disassemble/1/1
	# ${_polap_cmd} disassemble -o ${output_dir} \
	# 	-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
	# 	--stages-include 1 \
	# 	--disassemble-i 1 \
	# 	--disassemble-p 50 \
	# 	--disassemble-n 5 \
	# 	--disassemble-r 2

	# rm -rf ${output_dir}/disassemble/1/2
	# ${_polap_cmd} disassemble -o ${output_dir} \
	# 	-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
	# 	--stages-include 2 \
	# 	--disassemble-i 1 \
	# 	--disassemble-s 60m \
	# 	--disassemble-alpha 0.25 \
	# 	--disassemble-r 2

	# ${_polap_cmd} disassemble -o ${output_dir} \
	# 	-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
	# 	--stages-include 3 \
	# 	--disassemble-i 1 \
	# 	--disassemble-s 60m \
	# 	--genomesize 520115 \
	# 	--random-seed 27046 -v -v -v

	# Stage 1: v0.4.1.7
	# ${_polap_cmd} disassemble -o ${output_dir} \
	# 	-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
	# 	--disassemble-a 50m \
	# 	--disassemble-b 1g \
	# 	--disassemble-n 30 \
	# 	--disassemble-compare-to-fasta ptdna-${output_dir}.fa -v

	# rm -rf ${output_dir}/disassemble/x
	# command time -v ${_polap_cmd} disassemble -o ${output_dir} \
	# 	-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
	# 	--disassemble-compare-to-fasta ptdna-${output_dir}.fa \
	# 	--disassemble-s 82m --disassemble-alpha 0.25 \
	# 	--disassemble-n 10 \
	# 	--disassemble-stop-after assemble \
	# 	-v
}

run_Arctostaphylos_glauca() {
	local output_dir="$(echo $FUNCNAME | sed s/run_//)"
	local long_sra="SRR14883332"
	# if [[ ! -s "ptdna-${output_dir}.fa" ]]; then
	# 	${_polap_cmd} get-mtdna --plastid --species "Arctostaphylos glauca" -o ${output_dir}
	# 	cp ${output_dir}/00-bioproject/2-mtdna.fasta ptdna-${output_dir}.fa
	# fi
	${_polap_cmd} disassemble --flye-pacbio-hifi -o ${output_dir} -l ${long_sra}.fastq --genomesize 200k -v -v -v
}

run_Lepidium_sativum() {
	local output_dir="$(echo $FUNCNAME | sed s/run_//)"
	local species_name="$(echo $FUNCNAME | sed 's/run_//' | sed 's/_/ /')"
	local long_sra="SRR18079816"
	if [[ ! -s "ptdna-${output_dir}.fa" ]]; then
		${_polap_cmd} get-mtdna --plastid --species "${species_name}" -o ${output_dir}
		cp ${output_dir}/00-bioproject/2-mtdna.fasta ptdna-${output_dir}.fa
	fi
	# long-read data size:
	${_polap_cmd} disassemble -o ${output_dir} -l ${long_sra}.fastq --disassemble-compare-to-fasta ptdna-${output_dir}.fa --disassemble-alpha 0.01 --genomesize 200k --disassemble-s 555m -v
}

run_Chaetoceros_muellerii() {
	local output_dir="$(echo $FUNCNAME | sed s/run_//)"
	local species_name="$(echo $FUNCNAME | sed 's/run_//' | sed 's/_/ /')"
	local long_sra="SRR13238610"
	if [[ ! -s "ptdna-${output_dir}.fa" ]]; then
		${_polap_cmd} get-mtdna --plastid --species "${species_name}" -o ${output_dir}
		cp ${output_dir}/00-bioproject/2-mtdna.fasta ptdna-${output_dir}.fa
	fi
	# long-read data size:
	${_polap_cmd} disassemble -o ${output_dir} -l ${long_sra}.fastq --disassemble-compare-to-fasta ptdna-${output_dir}.fa --disassemble-alpha 0.01 --genomesize 200k --disassemble-s 555m -v
}

run_Potentilla_micrantha() {
	local output_dir="$(echo $FUNCNAME | sed s/run_//)"
	local species_name="$(echo $FUNCNAME | sed 's/run_//' | sed 's/_/ /')"
	local long_sra="ERR338629"
	if [[ ! -s "ptdna-${output_dir}.fa" ]]; then
		${_polap_cmd} get-mtdna --plastid --species "${species_name}" -o ${output_dir}
		cp ${output_dir}/00-bioproject/2-mtdna.fasta ptdna-${output_dir}.fa
	fi
	# long-read data size:
	${_polap_cmd} disassemble -o ${output_dir} -l ${long_sra}.fastq --disassemble-compare-to-fasta ptdna-${output_dir}.fa --disassemble-alpha 0.01 --genomesize 200k --disassemble-s 555m -v
}

run_Durio_zibethinus() {
	local output_dir="$(echo $FUNCNAME | sed s/run_//)"
	local species_name="$(echo $FUNCNAME | sed 's/run_//' | sed 's/_/ /')"
	local long_sra="SRR11547303"
	if [[ ! -s "ptdna-${output_dir}.fa" ]]; then
		${_polap_cmd} get-mtdna --plastid --species "${species_name}" -o ${output_dir}
		cp ${output_dir}/00-bioproject/2-mtdna.fasta ptdna-${output_dir}.fa
	fi
	# long-read data size:
	${_polap_cmd} disassemble -o ${output_dir} -l ${long_sra}.fastq --disassemble-compare-to-fasta ptdna-${output_dir}.fa --disassemble-alpha 0.01 --genomesize 200k --disassemble-s 555m -v
}

run_Beta_vulgaris() {
	local output_dir="$(echo $FUNCNAME | sed s/run_//)"
	local species_name="$(echo $FUNCNAME | sed 's/run_//' | sed 's/_/ /')"
	local long_sra="SRR1980665"
	if [[ ! -s "ptdna-${output_dir}.fa" ]]; then
		${_polap_cmd} get-mtdna --plastid --species "${species_name}" -o ${output_dir}
		cp ${output_dir}/00-bioproject/2-mtdna.fasta ptdna-${output_dir}.fa
	fi
	# long-read data size:
	${_polap_cmd} disassemble -o ${output_dir} -l ${long_sra}.fastq --disassemble-compare-to-fasta ptdna-${output_dir}.fa --disassemble-alpha 0.01 --genomesize 200k --disassemble-s 555m -v
}

run_Leucanthemum_vulgare() {
	local output_dir="$(echo $FUNCNAME | sed s/run_//)"
	local long_sra="SRR10948618"
	if [[ ! -s "ptdna-${output_dir}.fa" ]]; then
		${_polap_cmd} get-mtdna --plastid --species "Leucanthemum vulgare" -o ${output_dir}
		cp ${output_dir}/00-bioproject/2-mtdna.fasta ptdna-${output_dir}.fa
	fi
	${_polap_cmd} disassemble -o ${output_dir} -l ${long_sra}.fastq --disassemble-compare-to-fasta ptdna-${output_dir}.fa --disassemble-alpha 0.01 --genomesize 200k --disassemble-s 66m -v
}

run_Oryza_glaberrima() {
	local output_dir="$(echo $FUNCNAME | sed s/run_//)"
	local species_name="$(echo $FUNCNAME | sed 's/run_//' | sed 's/_/ /')"
	local long_sra="SRR8989349"
	if [[ ! -s "ptdna-${output_dir}.fa" ]]; then
		${_polap_cmd} get-mtdna --plastid --species "${species_name}" -o ${output_dir}
		cp ${output_dir}/00-bioproject/2-mtdna.fasta ptdna-${output_dir}.fa
	fi
	${_polap_cmd} disassemble -o ${output_dir} -l ${long_sra}.fastq --disassemble-compare-to-fasta ptdna-${output_dir}.fa --disassemble-alpha 0.01 --genomesize 200k --disassemble-s 317m -v
}

run_Cenchrus_americanus() {
	local output_dir="$(echo $FUNCNAME | sed s/run_//)"
	local species_name="$(echo $FUNCNAME | sed 's/run_//' | sed 's/_/ /')"
	local long_sra="SRR8989348"
	if [[ ! -s "ptdna-${output_dir}.fa" ]]; then
		${_polap_cmd} get-mtdna --plastid --species "${species_name}" -o ${output_dir}
		cp ${output_dir}/00-bioproject/2-mtdna.fasta ptdna-${output_dir}.fa
	fi
	${_polap_cmd} disassemble -o ${output_dir} -l ${long_sra}.fastq --disassemble-compare-to-fasta ptdna-${output_dir}.fa --disassemble-alpha 0.01 --genomesize 200k -v
}

run_Digitaria_exilis() {
	local output_dir="$(echo $FUNCNAME | sed s/run_//)"
	local species_name="$(echo $FUNCNAME | sed 's/run_//' | sed 's/_/ /')"
	local long_sra="SRR8989347"
	if [[ ! -s "ptdna-${output_dir}.fa" ]]; then
		${_polap_cmd} get-mtdna --plastid --species "${species_name}" -o ${output_dir}
		cp ${output_dir}/00-bioproject/2-mtdna.fasta ptdna-${output_dir}.fa
	fi
	# long-read data size: 555m
	# rm -rf ${output_dir}/disassemble
	${_polap_cmd} disassemble -o ${output_dir} -l ${long_sra}.fastq --disassemble-compare-to-fasta ptdna-${output_dir}.fa --disassemble-alpha 0.01 --genomesize 200k -v
}

run_Podococcus_acaulis() {
	local output_dir="$(echo $FUNCNAME | sed s/run_//)"
	local species_name="$(echo $FUNCNAME | sed 's/run_//' | sed 's/_/ /')"
	local long_sra="SRR8989346"
	if [[ ! -s "ptdna-${output_dir}.fa" ]]; then
		${_polap_cmd} get-mtdna --plastid --species "${species_name}" -o ${output_dir}
		cp ${output_dir}/00-bioproject/2-mtdna.fasta ptdna-${output_dir}.fa
	fi
	# long-read data size:
	# ${_polap_cmd} disassemble -o ${output_dir} -l ${long_sra}.fastq --disassemble-compare-to-fasta ptdna-${output_dir}.fa --disassemble-alpha 0.01 --genomesize 200k -v
	# ${_polap_cmd}
	seqkit seq -m 1000 ${long_sra}.fastq -o ${output_dir}/lk.fq.gz
	src/polap.sh flye1 -g 1m -o Podococcus_acaulis -l Podococcus_acaulis/lk.fq.gz
	src/polap.sh annotate -o Podococcus_acaulis
	src/polap.sh annotate -o Podococcus_acaulis view pt-table
	src/polap.sh -o Podococcus_acaulis seeds bandage
	src/polap.sh -o Podococcus_acaulis assemble2
	src/polap.sh -o Podococcus_acaulis annotate -i 1 view pt-table
}

run_Raphia_textilis() {
	local output_dir="$(echo $FUNCNAME | sed s/run_//)"
	local species_name="$(echo $FUNCNAME | sed 's/run_//' | sed 's/_/ /')"
	local long_sra="SRR8989345"
	if [[ ! -s "ptdna-${output_dir}.fa" ]]; then
		${_polap_cmd} get-mtdna --plastid --species "${species_name}" -o ${output_dir}
		cp ${output_dir}/00-bioproject/2-mtdna.fasta ptdna-${output_dir}.fa
	fi
	# long-read data size:
	# ${_polap_cmd} disassemble -o ${output_dir} -l ${long_sra}.fastq --disassemble-compare-to-fasta ptdna-${output_dir}.fa --disassemble-alpha 0.01 --genomesize 200k --disassemble-s 555m -v
	# ${_polap_cmd} get-dna-by-accession NC_020365 ptdna-${output_dir}.fa
	# seqkit seq -m 1000 ${long_sra}.fastq -o ${output_dir}/lk.fq.gz
	# ${_polap_cmd} flye1 -g 1m -o ${output_dir} -l ${output_dir}/lk.fq.gz
	# ${_polap_cmd} annotate -o ${output_dir}
	${_polap_cmd} annotate -o ${output_dir} view pt-table
}

run_Phytelephas_aequatorialis() {
	local output_dir="$(echo $FUNCNAME | sed s/run_//)"
	local species_name="$(echo $FUNCNAME | sed 's/run_//' | sed 's/_/ /')"
	local long_sra="SRR8989344"
	# no reference ptDNA
	# other ptDNA: ON248677.1 newer ptDNA not NC_029957
	if [[ ! -s "ptdna-${output_dir}.fa" ]]; then
		${_polap_cmd} get-mtdna --plastid --species "${species_name}" -o ${output_dir}
		cp ${output_dir}/00-bioproject/2-mtdna.fasta ptdna-${output_dir}.fa
	fi
	# long-read data size: 505m
	${_polap_cmd} disassemble -o ${output_dir} -l ${long_sra}.fastq --disassemble-compare-to-fasta ptdna-${output_dir}.fa --disassemble-alpha 0.01 --genomesize 200k -v
	# seqkit seq -m 1000 ${long_sra}.fastq -o ${output_dir}/lk.fq.gz
	# ${_polap_cmd} flye1 -g 1m -o ${output_dir} -l ${output_dir}/lk.fq.gz
	# ${_polap_cmd} annotate -o ${output_dir}
	# ${_polap_cmd} annotate -o ${output_dir} view pt-table
	# ${_polap_cmd} seeds bandage -o ${output_dir}
	# ${_polap_cmd} assemble2 -o ${output_dir} -w 1000 --plastid
	# ${_polap_cmd} annotate -o ${output_dir} -i 1
	# ${_polap_cmd} annotate -o ${output_dir} -i 1 view pt-table
}

run_Picea_glauca() {
	local output_dir="$(echo $FUNCNAME | sed s/run_//)"
	local species_name="$(echo $FUNCNAME | sed 's/run_//' | sed 's/_/ /')"
	local long_sra="SRR21038753"
	if [[ ! -s "ptdna-${output_dir}.fa" ]]; then
		${_polap_cmd} get-mtdna --plastid --species "${species_name}" -o ${output_dir}
		cp ${output_dir}/00-bioproject/2-mtdna.fasta ptdna-${output_dir}.fa
	fi
	# long-read data size:
	${_polap_cmd} disassemble --flye-pacbio-corr -o ${output_dir} -l ${long_sra}.fastq --disassemble-compare-to-fasta ptdna-${output_dir}.fa --genomesize 200k -v -v -v
}

# Main case statement
case "$species_folder" in
"Juncus_effusus" | \
	"Juncus_effusus-a" | \
	"Juncus_inflexus" | \
	"Juncus_inflexus-a" | \
	"Juncus_roemerianus" | \
	"Juncus_roemerianus-a" | \
	"Juncus_validus" | \
	"Juncus_validus-a" | \
	"Eucalyptus_pauciflora" | \
	"Eucalyptus_pauciflora-a" | \
	"Arctostaphylos_glauca" | \
	"Lepidium_sativum" | \
	"Chaetoceros_muellerii" | \
	"Potentilla_micrantha" | \
	"Durio_zibethinus" | \
	"Beta_vulgaris" | \
	"Eleocharis_dulcis" | \
	"Leucanthemum_vulgare" | \
	"Oryza_glaberrima" | \
	"Cenchrus_americanus" | \
	"Digitaria_exilis" | \
	"Podococcus_acaulis" | \
	"Raphia_textilis" | \
	"Phytelephas_aequatorialis" | \
	"test" | \
	"Picea_glauca")
	run_${species_folder}
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

exit

# Juncus effusus
# SRR14298760 35261850229 WGS GENOMIC SINGLE OXFORD_NANOPORE Juncus effusus
# SRR14298746 28996069500 WGS GENOMIC PAIRED ILLUMINA Juncus effusus

src/polap.sh disassemble -o j.effusus -l SRR14298760.fastq -a SRR14298746_1.fastq -a SRR14298746_2.fastq

# Eucalyptus pauciflora
rm -rf epauciflora
src/polap.sh disassemble -o epauciflora -l SRR7153095.fastq --disassemble-s-max 5m --disassemble-compare-to-fasta ptdna-epauciflora.fa --disassemble-min-memory 9 -v -v >1.txt

exit

# src/polap.sh x-get-sra-info --sra <SRA>

S=(
	SRR14883332
	SRR18079816
	SRR13238610
	ERR338629
	SRR11547303
	SRR1980665
	SRR7153095
	SRR10948618
	SRR8989349
	SRR8989348
	SRR8989347
	SRR8989346
	SRR8989345
	SRR8989344
	SRR21038753
)

for i in "${S[@]}"; do
	# src/polap.sh x-ncbi-fetch-sra --sra "$i"
	# sleep 5
	rm -rf "$i"
done

exit
#!/usr/bin/bash

# src/polap.sh x-get-sra-info --sra <SRA>

S=(SRR7153095
	SRR7161123
	SRR21976089
	SRR21976090
	SRR21976091
	SRR21976092
	SRR14298760
	SRR14298751
	SRR14298746
	SRR14298745

	SRR14883332
	SRR18079816
	SRR13238610
	ERR338629
	SRR11547303
	SRR1980665
	SRR7153095
	SRR10948618
	SRR8989349
	SRR8989348
	SRR8989347
	SRR8989346
	SRR8989345
	SRR8989344
	SRR21038753)

for i in "${S[@]}"; do
	src/polap.sh x-get-sra-info --sra "$i"
	sleep 5
done

exit

# Table 2 of Zhou2023
SRR14298760 35261850229 WGS GENOMIC SINGLE OXFORD_NANOPORE Juncus effusus
SRR14298746 28996069500 WGS GENOMIC PAIRED ILLUMINA Juncus effusus

SRR14298751 37591299287 WGS GENOMIC SINGLE OXFORD_NANOPORE Juncus inflexus
SRR14298745 25023621900 WGS GENOMIC PAIRED ILLUMINA Juncus inflexus

SRR21976090 2302422897 WGS GENOMIC SINGLE OXFORD_NANOPORE Juncus roemerianus
SRR21976092 79461161000 WGS GENOMIC PAIRED ILLUMINA Juncus roemerianus

SRR21976089 561799410 WGS GENOMIC SINGLE OXFORD_NANOPORE Juncus validus
SRR21976091 78356215000 WGS GENOMIC PAIRED ILLUMINA Juncus validus

# Table 1 of Zhou2023: with short-read data as well
SRR7153095 6040851677 WGS GENOMIC SINGLE OXFORD_NANOPORE Eucalyptus pauciflora
SRR7161123 6376665372 WGS GENOMIC PAIRED ILLUMINA Eucalyptus pauciflora

# Table 1 of Zhou2023
SRR14883332 27304352670 WGS GENOMIC SINGLE PACBIO_SMRT Arctostaphylos glauca
SRR18079816 2615947808 WGA GENOMIC SINGLE PACBIO_SMRT Lepidium sativum
SRR13238610 792414508 Targeted-Capture OTHER SINGLE PACBIO_SMRT Chaetoceros muellerii
ERR338629 54492250 WGS GENOMIC SINGLE PACBIO_SMRT Potentilla micrantha
SRR11547303 5987567965 OTHER GENOMIC SINGLE PACBIO_SMRT Durio zibethinus
SRR1980665 296752589 WGS GENOMIC SINGLE PACBIO_SMRT Beta vulgaris subsp. vulgaris
# Eleocharis dulcis - available by asking authors (Ref: NC_047447.1)
SRR7153095 6040851677 WGS GENOMIC SINGLE OXFORD_NANOPORE Eucalyptus pauciflora
SRR10948618 66965150 AMPLICON GENOMIC SINGLE OXFORD_NANOPORE Leucanthemum vulgare
SRR8989349 317851242 Targeted-Capture GENOMIC SINGLE OXFORD_NANOPORE Oryza glaberrima
SRR8989348 525852568 Targeted-Capture GENOMIC SINGLE OXFORD_NANOPORE Cenchrus americanus
SRR8989347 555524185 Targeted-Capture GENOMIC SINGLE OXFORD_NANOPORE Digitaria exilis
SRR8989346 645897952 Targeted-Capture GENOMIC SINGLE OXFORD_NANOPORE Podococcus acaulis
SRR8989345 202292169 Targeted-Capture GENOMIC SINGLE OXFORD_NANOPORE Raphia textilis
SRR8989344 505797880 Targeted-Capture GENOMIC SINGLE OXFORD_NANOPORE Phytelephas aequatorialis

SRR21038753 5928239580 WGS GENOMIC SINGLE OXFORD_NANOPORE Picea glauca

# Eucalyptus pauciflora
src/polap.sh x-ncbi-fetch-sra --sra SRR7153095
src/polap.sh x-ncbi-fetch-sra --sra SRR7161123

#
src/polap.sh x-ncbi-fetch-sra --sra SRR21976089
src/polap.sh x-ncbi-fetch-sra --sra SRR21976090

src/polap.sh x-ncbi-fetch-sra --sra SRR21976091
# J. validus (SRR21976091)
# J. roamerianus (SRR21976092)
src/polap.sh x-ncbi-fetch-sra --sra SRR21976092

#
src/polap.sh x-ncbi-fetch-sra --sra SRR14298760
src/polap.sh x-ncbi-fetch-sra --sra SRR14298751

src/polap.sh x-ncbi-fetch-sra --sra SRR14298746
src/polap.sh x-ncbi-fetch-sra --sra SRR14298745

src/polap.sh x-ncbi-fetch-sra-runinfo --sra "$long_sra"

exit

#Download Long read data or Illumina data
module add sratoolkit/3.0.0

fasterq-dump --split-files SRR14883332
#Download reference genome
#Long read data
esearch -db nucleotide -query "NC_035584.1" | efetch -format fasta >NC_035584.1.fasta

# Replace SRA accession numbers (SRR14883332) and reference number (NC_035584.1) with the following numbers.
# Long read data for ptGAUL validation:

# Arctostaphylos glauca (SRR14883332) (Ref: NC_035584.1/NC_042713.1/NC_047438.1)
SRR14883332

# Lepidium sativum (SRR18079816) (Ref: NC_047178.1)
SRR18079816

# Chaetoceros muellerii (SRR13238610) (Ref: MW004650.1)
SRR13238610

# Potentilla micrantha (ERR338629) (Ref: NC_015206.1)
ERR338629

# Durio zibethinus (SRR11547303) (Ref: MT321069)
SRR11547303
# Beta vulgaris (SRR1980665) (Ref: KR230391.1)
SRR1980665
# Eleocharis dulcis - available by asking authors (Ref: NC_047447.1)

# Eucalyptus pauciflora (SRR7153095) (Ref: MZ670598.1/HM347959.1/NC_014570.1/AY780259.1/ NC_039597.1)
SRR7153095

# Leucanthemum vulgare (SRR10948618) (Ref: NC_047460.1)
SRR10948618

# Oryza glaberrima (SRR8989349) (Ref: NC_024175.1)
SRR8989349

# Cenchrus americanus (SRR8989348) (Ref: NC_024171.1)
SRR8989348

# Digitaria exilis (SRR8989347) (Ref: NC_024176.1)
SRR8989347

# Podococcus acaulis (SRR8989346) (Ref: NC_027276.1)
SRR8989346

# Raphia textilis (SRR8989345) (Ref: NC_020365.1)
SRR8989345

# Phytelephas aequatorialis (SRR8989344) (Ref: NC_029957.1)
SRR8989344

# Picea glauca (SRR21038753) (Ref: NC_021456.1).
SRR21038753

# Illumina data of Juncus for GetOrganelle:

# J. inflexus (SRR14298745)
# J. effusus (SRR13309655 (Lu et al. 2021) and SRR14298746)
# J. roamerianus (SRR21976092)
# J. validus (SRR21976091)

# example of using GetOrganelle for J. effusus
fasterq-dump --split-files SRR13309655

# load GetOrganelle
get_organelle_from_reads.py -1 SRR13309655.1_1.fastq -2 SRR13309655.1_2.fastq -o SSR55 -R 15 -k 21,39,45,65,85,105,115,125 -F embplant_pt -t 20
get_organelle_from_reads.py -1 SRR13309655.1_1.fastq -2 SRR13309655.1_2.fastq -o SSR55_w075 -R 15 -k 21,39,45,65,85,105,115,125 -F embplant_pt -t 20 -w 0.75

# To run the GetOrganelle command, when the input data is too large, choose a fraction of the large dataset using seqtk module.
# example: J. validus SRR21976091_1.fastq and SRR21976091_2.fastq have about 87G, respectively.
module add seqtk
seqtk sample -s100 SRR21976091_1.fastq 0.05 >J_validus_s1.fastq
seqtk sample -s100 SRR21976091_2.fastq 0.05 >J_validus_s2.fastq
# long-read only
# rm -rf epauciflora/disassemble
# rm -rf epauciflora1/disassemble
rm -rf epauciflora1/disassemble/x
src/polap.sh disassemble -o epauciflora1 -l SRR7153095.fastq -g 250k --disassemble-s 130m --disassemble-compare-to-fasta ptdna-epauciflora.fa --disassemble-min-memory 9 -v -v -v >1.txt
exit

rm -rf epauciflora
src/polap.sh disassemble -o epauciflora -l SRR7153095.fastq --disassemble-s-max 5m --disassemble-compare-to-fasta ptdna-epauciflora.fa --disassemble-min-memory 9 -v -v >1.txt
exit

# Eucalyptus pauciflora
src/polap.sh disassemble -o epauciflora -l SRR7153095.fastq -a SRR7161123_1.fastq -b SRR7161123_2.fastq --disassemble-compare-to-fasta ptdna-epauciflora.fa --disassemble-min-memory 9 -v -v -v

exit

polap.sh disassemble view
3:
polap.sh disassemble view 2
polap.sh disassemble -o o2 --anotherdir o --disassemble-s 314588465 --disassemble-alpha 3.5 -l SRR7153095.fastq -a SRR7161123_1.fastq -b SRR7161123_2.fastq

For ptDNA of Juncus validus
polap.sh get-mtdna -o jvalidus --plastid --species "Juncus validus"
polap.sh disassemble -o jvalidus -l SRR21976089.fastq -a SRR21976091_1.fastq -b SRR21976091_2.fastq --disassemble-min-memory 9 -v --disassemble-a 50mb --disassemble-n 20
src/polap.sh disassemble -l SRR21976089.fastq -a SRR21976091_1.fastq -b SRR21976091_2.fastq -o jvalidus --disassemble-min-memory 9 -v --disassemble-a 50mb --disassemble-n 20
