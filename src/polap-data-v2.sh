#!/usr/bin/bash

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" || {
	echo "Couldn't determine the script's running directory, which probably matters, bailing out" >&2
	exit 2
}
LIB_DIR="${script_dir}/lib"

source "${LIB_DIR}/polap-lib-timing.sh"

_log_echo() {
	echo "$(date '+%Y-%m-%d %H:%M:%S') [$subarg1] - $1" >>"${output_dir}/polap-data-v2.txt"
	echo "$1"
}

Sall=(
	Anthoceros_agrestis
	Arabidopsis_thaliana
	Brassica_carinata
	Canavalia_ensiformis
	Cinchona_pubescens
	Codonopsis_lanceolata
	Cucumis_sativus_var_hardwickii
	Dioscorea_japonica
	Dunaliella_tertiolecta
	Eucalyptus_pauciflora
	Euonymus_alatus
	Gossypium_herbaceum
	Juncus_effusus
	Juncus_inflexus
	Juncus_roemerianus
	Juncus_validus
	Leiosporoceros_dussii
	Macadamia_jansenii
	Musa_acuminata_subsp_malaccensis
	Notothylas_orbicularis
	Ophrys_lutea
	Oryza_rufipogon
	Phaeomegaceros_chiloensis
	Pisum_sativum
	Populus_x_sibirica
	Prunus_mandshurica
	Pterocarpus_santalinus
	Solanum_lycopersicum
	Spirodela_polyrhiza
	Vaccinium_vitis-idaea
	Vitis_vinifera
)

Snot=(
	Delphinium_montanum
)

Stable5=(
	'Juncus_effusus'
	'Juncus_inflexus'
	'Juncus_roemerianus'
	'Juncus_validus'
	'Eucalyptus_pauciflora'
)

_local_host="thorne"

# Input parameter
subcmd1="${1:-help}"

_polap_subcmd=(
	'batch'
	'recover'
	'mkdir'
	'refs'
	'getorganelle'
	'ptgaul'
	'msbwt'
	'extract-ptdna-of-ptgaul'
	'coverage'
	'infer'
	'check'
	'compare'
	'best'
	'bandage1'
	'bandage2'
	'copy-figure'
	'downsample'
	'use-downsample'
	'mauve'
	'restart'
	'write-config'
	'wga'
	'archive'
	'get'
	'report'
	'table1'
	'table2'
	'suptable1'
	'supfigure1'
	'supfigure2'
	'clean'
	'menu'
	'help'
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

if [[ -d "src" ]]; then
	_polap_cmd="src/polap.sh"
else
	_polap_cmd="polap"
fi
_polap_version="0.4.3.7"
_media1_dir="/media/h1/sra"
_media2_dir="/media/h2/sra"
_media_dir="/media/h2/sra"

help_message=$(
	cat <<HEREDOC
# Polap data analysis for the version 0.4 or disassemble (subsampling method)
#
# Edit two variables:
# _polap_cmd=${_polap_cmd}
# _media_dir=${_media_dir}
# _local_host=${_local_host}
# host=$(hostname)
#
# subcommand can be replaced with the leading number
0. batch <species_folder> [number]: To execute all subcommands.
1. sra <species_folder>: scp data files to the remote 
1. send <species_folder>: scp data files to the remote 
2. recover <species_folder>
2. mkdir <species_folder>
3. refs <species_folder>: get-ptdna-from-ncbi <species_folder>
4. getorganelle <species_folder>
5. ptgaul <species_folder>: ptGAUL analysis on the data
6. msbwt <species_folder>: prepare short-read polishing
7. extract-ptdna-of-ptgaul <species_folder>: extract ptDNA from ptGAUL's result
8. coverage <species_folder> [number]: ptGAUL's result
9. infer <species_folder> [number] [--disassemble-simple-polishing]: assemble the ptDNA without a reference ptDNA
9. infer2 <species_folder> [number] [--disassemble-simple-polishing]: assemble the ptDNA without a reference ptDNA
10. check <species_folder> [number] [index:0] [--disassemble-simple-polishing]: assemble ptDNA with subsampling by comparing it with the ptGAUL assembly
11. compare <species_folder>: compare ptDNAs with the ptGAUL ptDNA
12. archive <species_folder>: archive the results
13. get <species_folder>: fetch the archive from the remote
14. report [species_folder]: report the results
15. table1: table1
16. table2: table2
17. table3: table3
17. table4: table4
18. suptable1: Supplementary tables
19. supfigure1: Supplementary figures
20. wga <species_folder>: whole-genome assembly
21. simple-polish <species_folder>: simple polish the ptDNA using the short-read data
22. subsample-polish <species_folder>: subsampling polish the ptDNA using the short-read data
23. mauve <species_folder>: subsampling polish the ptDNA using the short-read data
24. clean <species_folder>: [DANGER] delete the species folder
write-config src/polap-data-v2.csv: write config
#
# call the following menu to create tables and figures
#
bandage1 <species_folder> 2
bandage2 <species_folder> 2
copy-figure: copy all figures to the target directory
table1 2
table1 0
table2 2
table2 0
suptable1 2 off 1
suptable1 2 off 2
suptable1 2 off 3-infer
suptable1 2 on 1
suptable1 2 on 2
suptable1 2 on 3-infer
supfigure1 2 no off
supfigure1 2 no on
supfigure2 2 off <- for sup
supfigure2 2 on
copy-figure: copy all figures to the target directory

# files
table1-1.md
table1-2.md
table2.md
suptable1.md
supfigure2.md
supfigure1.md

HEREDOC
)

declare -A _folder
declare -A _memory
declare -A _inum
declare -A _table1
declare -A _table2
declare -A _dummy
declare -A _downsample
declare -A _status
declare -A _host
declare -A _long
declare -A _short
declare -A _ptgaul_genomesize
declare -A _compare_n
declare -A _compare_p
declare -A _compare_r
declare -A _random_seed
declare -A _ssh
declare -A _polish_n
declare -A _polish_p

set +u

# Read the config files
read-a-tsv-file-into-associative-arrays() {
	# Define input TSV file
	csv_file="${script_dir}/polap-data-v2.csv"

	# Read the TSV file (skip header)
	while IFS=$',' read -r species folder long short host ptgaul_genomesize compare_n compare_p compare_r polish_n polish_p random_seed ssh memory downsample inum table1 table2 dummy status; do
		# Skip header line
		[[ "$species" == "species" ]] && continue

		# Store in associative arrays
		_folder["$species"]="$folder"
		_long["$species"]="$long"
		_short["$species"]="$short"
		_host["$species"]="$host"
		_ptgaul_genomesize["$species"]="$ptgaul_genomesize"
		_compare_n["$species"]="$compare_n"
		_compare_p["$species"]="$compare_p"
		_compare_r["$species"]="$compare_r"
		_polish_n["$species"]="$polish_n"
		_polish_p["$species"]="$polish_p"
		_random_seed["$species"]="$random_seed"
		_ssh["$species"]="$ssh"
		_memory["$species"]="$memory"
		_downsample["$species"]="$downsample"
		_inum["$species"]="$inum"
		_table1["$species"]="$table1"
		_table2["$species"]="$table2"
		_dummy["$species"]="$dummy"
		_status["$species"]="$status"
	done <"$csv_file"
}

read-a-tsv-file-into-associative-arrays

keys_array=($(for key in "${!_long[@]}"; do echo "$key"; done | sort))

_arg2=${2:-"default_value"}
if [[ "${_arg2}" != "default_value" ]]; then
	_arg2="${2%/}"
fi

_arg3=${3:-0}
_arg4=${4:-0}

_arg1=${1:-"default_value"}
# Check if the species folder is provided
if [[ "${_arg1}" == "default_value" ]]; then
	echo "Usage: $0 <subcommand> [species_folder]"
	echo "${help_message}"
	exit 1
fi

################################################################################
# Part of genus_species
#
batch_genus_species() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local do_simple_polishing="${3:-off}"
	local target_index="${output_dir}-${isuffix}"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local ssh_remote="${_ssh["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	# rm -rf "${output_dir}"
	# tar -zxf "${output_dir}-a.tar.gz"
	# mv "${output_dir}-a" "${output_dir}"

	echo host: $(hostname)
	echo ssh-remote: $ssh_remote
	echo output: $output_dir
	echo inum: $extracted_inum
	echo species: $species_name
	echo random seed: $random_seed
	local long_data="${_media_dir}/${long_sra}.fastq.tar.gz"
	local short_data="${_media_dir}/${short_sra}.fastq.tar.gz"

	if [[ "${_local_host}" == "$(hostname)" ]]; then
		if [[ -s "${long_data}" ]]; then
			echo long: $(du -h ${long_data})
			echo short: $(du -h ${short_data})
		else
			echo long: $(du -h "${_media1_dir}/${long_sra}.fastq")
			echo short1: $(du -h "${_media1_dir}/${short_sra}_1.fastq")
			echo short2: $(du -h "${_media1_dir}/${short_sra}_2.fastq")
		fi
		return 0
	fi

	if [[ -s "${long_sra}.fastq.tar.gz" ]] ||
		[[ -s "${long_sra}.fastq" ]]; then

		read -p "Do you want to execute the batch procedure? (y/N): " confirm
		case "$confirm" in
		[yY] | [yY][eE][sS])

			if [[ -s "${output_dir}/polap.log" ]]; then
				echo "Folder ${output_dir} is not empty"
			else
				echo "Folder ${output_dir} is empty"
				recover_genus_species "${output_dir}" "${isuffix}"
			fi

			mkdir -p "${output_dir_i}"

			if [[ -s "${long_sra}.fastq" ]] &&
				[[ -s "${short_sra}_1.fastq" ]] &&
				[[ -s "${short_sra}_2.fastq" ]]; then
				_log_echo "Found: sequencing data"
			else
				mkdir_genus_species "${output_dir}" "${isuffix}"
				if [[ -s "${long_sra}.fastq" ]] &&
					[[ -s "${short_sra}_1.fastq" ]] &&
					[[ -s "${short_sra}_2.fastq" ]]; then
					_log_echo "Success: sequencing data"
				else
					_log_echo "Fail: sequencing data"
					return 1
				fi
			fi

			if [[ -d "${output_dir}/getorganelle" ]] &&
				find "${output_dir}/getorganelle" -maxdepth 1 -type f -name 'embplant_pt.*.gfa' -size +0c | grep -q .; then
				_log_echo "Found: GetOrganelle assembled ptDNA"
			else
				getorganelle_genus_species "${output_dir}"
				if find "${output_dir}/getorganelle" -maxdepth 1 -type f -name 'embplant_pt.*.gfa' -size +0c | grep -q .; then
					_log_echo "Success: GetOrganelle assembled ptDNA"
				else
					_log_echo "Fail: GetOrganelle assembled ptDNA"
					return 1
				fi
			fi

			if [[ -s "${output_dir}/ptdna-reference.fa" ]]; then
				_log_echo "Found: reference ptDNA"
			else
				get-ptdna-from-ncbi_genus_species "${output_dir}"
				if [[ -s "${output_dir}/ptdna-reference.fa" ]]; then
					_log_echo "Success: reference ptDNA"
				else
					_log_echo "Fail: reference ptDNA"
					return 1
				fi
			fi

			if [[ -s "${output_dir}/ptgaul/flye_cpONT/assembly_graph.gfa" ]]; then
				_log_echo "Found: ptGAUL assembly"
			else
				ptgaul_genus_species "${output_dir}"
				if [[ -s "${output_dir}/ptgaul/flye_cpONT/assembly_graph.gfa" ]]; then
					_log_echo "Success: ptGAUL assembly"
				else
					_log_echo "Fail: ptGAUL assembly"
					return 1
				fi
			fi

			if [[ -s "${output_dir}/msbwt/comp_msbwt.npy" ]]; then
				_log_echo "Found: FMLRC msbwt"
			else
				msbwt_genus_species "${output_dir}"
				if [[ -s "${output_dir}/msbwt/comp_msbwt.npy" ]]; then
					_log_echo "Success: FMLRC msbwt"
				else
					_log_echo "Fail: FMLRC msbwt"
					return 1
				fi
			fi

			if [[ -s "${output_dir}/ptdna-ptgaul.fa" ]]; then
				_log_echo "Found: ptGAUL polished genome"
			else
				extract-ptdna-of-ptgaul_genus_species "${output_dir}"
				copy-ptdna-of-ptgaul_genus_species "${output_dir}"
				if [[ -s "${output_dir}/ptdna-ptgaul.fa" ]]; then
					_log_echo "Success: ptGAUL polished genome"
				else
					_log_echo "Fail: ptGAUL polished genome"
					return 1
				fi
			fi

			if [[ -s "${output_dir_i}/disassemble/infer-1/pt.subsample-polishing.1.fa" ]]; then
				_log_echo "Found: infer case"
			else
				infer_genus_species "${output_dir}" "${isuffix}"
				if [[ -s "${output_dir_i}/disassemble/infer-1/pt.subsample-polishing.1.fa" ]]; then
					_log_echo "Success: infer case"
				else
					_log_echo "Fail: infer case"
					return 1
				fi
			fi

			if [[ "${do_simple_polishing}" == "on" ]]; then

				if [[ -s "${output_dir_i}/disassemble/infer-1/pt.simple-polishing.1.fa" ]]; then
					_log_echo "Found: infer case - simple polishing"
				else
					infer_genus_species "${output_dir}" "${isuffix}" --disassemble-simple-polishing
					if [[ -s "${output_dir_i}/disassemble/infer-1/pt.simple-polishing.1.fa" ]]; then
						_log_echo "Success: infer case - simple polishing"
					else
						_log_echo "Fail: infer case - simple polishing"
						_log_echo "  potenial error: fmlrc not enough memory"
						return 1
					fi
				fi
			fi

			if [[ -s "${output_dir}/ptgaul/flye_cpONT/assembly_graph.gfa" ]]; then

				if [[ -s "${output_dir_i}/disassemble/infer-1/pt.subsample-polishing.reference.aligned.1.fa" ]]; then
					_log_echo "Found: check case"
				else
					check_genus_species "${output_dir}" "${isuffix}"
					if [[ -s "${output_dir_i}/disassemble/infer-1/pt.subsample-polishing.reference.aligned.1.fa" ]]; then
						_log_echo "Success: check case"
					else
						_log_echo "Fail: check case"
						return 1
					fi
				fi

				if [[ "${do_simple_polishing}" == "on" ]]; then
					if [[ -s "${output_dir_i}/disassemble/infer-1/pt.simple-polishing.reference.aligned.1.fa" ]]; then
						_log_echo "Found: check case - simple polishing"
					else
						check_genus_species "${output_dir}" "${isuffix}" --disassemble-simple-polishing
						if [[ -s "${output_dir_i}/disassemble/infer-1/pt.simple-polishing.reference.aligned.1.fa" ]]; then
							_log_echo "Success: check case - simple polishing"
						else
							_log_echo "Fail: check case - simple polishing"
							_log_echo "  potenial error: fmlrc not enough memory"
							return 1
						fi
					fi
				fi

				if printf '%s\n' "${Stable5[@]}" | grep -qx "${output_dir}"; then
					if [[ -s "${output_dir_i}/disassemble/compare-1/pt.simple-polishing.1.fa" ]]; then
						_log_echo "Found: compare case"
					else
						compare_genus_species "${output_dir}" "${isuffix}"
						if [[ -s "${output_dir_i}/disassemble/compare-1/pt.simple-polishing.1.fa" ]]; then
							_log_echo "Success: compare case"
						else
							_log_echo "Fail: compare case"
							return 1
						fi
					fi
				fi
			fi

			;;
		*)
			echo "Batch procedure is canceled."
			;;
		esac
	fi

}

send-data_genus_species() {
	local output_dir="$1"
	local inum="${2:-0}"
	local target_index="${output_dir}-${inum}"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local ssh_name="${_ssh["$target_index"]}"

	echo host: $(hostname)
	echo remote dir: ${ssh_name}:$PWD
	echo output: ${output_dir}
	echo target_index: ${target_index}
	echo species: $species_name
	echo long: $long_sra
	echo short: $short_sra

	local long_data="${_media_dir}/${long_sra}.fastq.tar.gz"
	local short_data="${_media_dir}/${short_sra}.fastq.tar.gz"

	if [[ "${_local_host}" == "$(hostname)" ]]; then
		if [[ -s "${long_data}" ]]; then
			echo long: $(du -h ${long_data})
			echo short: $(du -h ${short_data})
		else
			echo long: $(du -h "${_media1_dir}/${long_sra}.fastq")
			echo short1: $(du -h "${_media1_dir}/${short_sra}_1.fastq")
			echo short2: $(du -h "${_media1_dir}/${short_sra}_2.fastq")
		fi
	else
		echo "You can send the data from the central host: ${_local_host}"
		return 1
	fi

	if [[ "${_local_host}" == "$(hostname)" ]]; then

		ssh ${ssh_name} "mkdir -p $PWD/${output_dir}"

		if [[ -s "${long_data}" ]]; then
			scp "${long_data}" ${ssh_name}:$PWD/
		else
			scp \
				"${_media1_dir}/${long_sra}.fastq" \
				${ssh_name}:$PWD/
		fi

		if [[ -s "${short_data}" ]]; then
			scp "${short_data}" ${ssh_name}:$PWD/
		else
			scp \
				"${_media1_dir}/${short_sra}_1.fastq" \
				"${_media1_dir}/${short_sra}_2.fastq" \
				${ssh_name}:$PWD/
		fi
	else
		echo "ERROR: run at the local host."
	fi
}

# send the input data to the remote
send_genus_species() {
	local output_dir="$1"
	local inum="${2:-0}"
	local target_index="${output_dir}-${inum}"

	if [[ "${_local_host}" == "$(hostname)" ]]; then
		archive_genus_species "${output_dir}"
		scp "${output_dir}-a.tar.gz" ${_ssh["${target_index}"]}:$PWD/
	else
		echo "ERROR: run at the local host."
	fi

	send-data_genus_species "${output_dir}" "${inum}"
}

# extract the archive at the remote
recover_genus_species() {
	local output_dir="$1"

	if [[ -s "${output_dir}-a.tar.gz" ]]; then
		echo "Deleting ${output_dir} ..."
		rm -rf "${output_dir}"
		tar -zxf "${output_dir}-a.tar.gz"
		mv "${output_dir}-a" "${output_dir}"
		echo "  we have recreated ${output_dir}!"
		echo "Next: $0 mkdir ${output_dir}"
	else
		mkdir -p "${output_dir}"
	fi
}

# create input files
mkdir_genus_species() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local target_index="${output_dir}-${isuffix}"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local ssh_remote="${_ssh["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	local long_data="${_media_dir}/${long_sra}.fastq.tar.gz"
	local short_data="${_media_dir}/${short_sra}.fastq.tar.gz"

	echo "create ${output_dir} ..."
	mkdir -p "${output_dir}/timing"
	if [[ -s "${long_sra}.fastq" ]]; then
		echo "  found: long SRA: ${long_sra}.fastq"
	else
		if tar -zxf "$(basename ${long_data})"; then
			rm "$(basename ${long_data})"
			echo "Extraction successful: ${long_sra}. Archive deleted."
		fi
	fi
	if [[ -s "${short_sra}_1.fastq" ]] && [[ -s "${short_sra}_2.fastq" ]]; then
		echo "  found: short SRA1: ${short_sra}_1.fastq"
		echo "  found: short SRA1: ${short_sra}_2.fastq"
	else
		if tar -zxf "$(basename ${short_data})"; then
			rm "$(basename ${short_data})"
			echo "Extraction successful: ${short_sra} Archive deleted."
		fi
	fi

	echo "Next: $0 refs ${output_dir} [number] to download reference ptDNAs from NCBI"
	echo "Next: $0 coverage ${output_dir} [number] to overview your data"
}

sra_genus_species() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local target_index="${output_dir}-${isuffix}"

	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"

	echo "create ${output_dir} ..."
	mkdir -p "${output_dir}/timing"
	run-polap-ncbitools fetch sra "$long_sra"
	run-polap-ncbitools fetch sra "$short_sra"

	echo "Next: $0 refs ${output_dir} [number] to download reference ptDNAs from NCBI"
	echo "Next: $0 coverage ${output_dir} [number] to overview your data"
}

get-ptdna-from-ncbi_genus_species() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local target_index="${output_dir}-${isuffix}"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local ssh_remote="${_ssh["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	if [[ "${output_dir}" == "Juncus_inflexus" ]]; then
		species_name="Juncus effusus"
		echo "No ptDNA for ${output_dir}, so we use ${species_name}"
	fi

	if [[ -s "${output_dir}/ptdna-reference.fa" ]]; then
		echo "found: ptDNA reference: ${output_dir}/ptdna-reference.fa"
	else

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
	fi
}

copy-ptdna-of-ncbi-as-reference_genus_species() {
	local output_dir="$1"

	echo "copy ${output_dir}/ptdna-reference.fa"
	cp -p "${output_dir}/00-bioproject/2-mtdna.fasta" \
		"${output_dir}/ptdna-reference.fa"
}

clean_genus_species() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local target_index="${output_dir}-${isuffix}"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local ssh_remote="${_ssh["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

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

getorganelle_genus_species_for() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local target_index="${output_dir}-${isuffix}"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local ssh_remote="${_ssh["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	if [[ "${_local_host}" == "$(hostname)" ]]; then
		decompress-short_genus_species "${output_dir}"
	else
		echo "  run at the local host for the decompressing the data"
	fi

	if [[ -s "${short_sra}_1.fastq" ]] &&
		[[ -s "${short_sra}_2.fastq" ]]; then
		# Initialize Conda
		source $HOME/miniconda3/etc/profile.d/conda.sh
		conda activate getorganelle

		# Example
		# -1 Arabidopsis_simulated.1.fq.gz -2 Arabidopsis_simulated.2.fq.gz -t 1 -o Arabidopsis_simulated.plastome -F embplant_pt -R 10
		# default
		# -R 15 \
		mkdir -p "${output_dir}/timing"
		command time -v get_organelle_from_reads.py \
			-1 ${short_sra}_1.fastq \
			-2 ${short_sra}_2.fastq \
			-o ${output_dir}-getorganelle \
			-t 24 \
			-F embplant_pt \
			2>${output_dir}/timing/timing-getorganelle.txt

		rsync -aq \
			--max-size=5M \
			"${output_dir}-getorganelle/" \
			"${output_dir}/getorganelle/"
		rm -rf "${output_dir}-getorganelle"

		conda deactivate

		if [[ "${_local_host}" == "$(hostname)" ]]; then
			delete-short_genus_species "${output_dir}"
		else
			echo "  run at the local host for deleting the data"
		fi
	fi
}

decompress-short_genus_species() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local target_index="${output_dir}-${isuffix}"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local ssh_remote="${_ssh["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	local short_data="${_media_dir}/${short_sra}.fastq.tar.gz"

	if [[ -s "${short_data}" ]]; then
		echo "decompressing the short-read data ..."
		tar zxf "${short_data}"
	else
		echo "copying the short-read data ..."
		cp -p "/media/h1/sra/${short_sra}_1.fastq" .
		cp -p "/media/h1/sra/${short_sra}_2.fastq" .
	fi

	if [[ -s "${short_sra}_1.fastq" ]] &&
		[[ -s "${short_sra}_2.fastq" ]]; then
		echo short: $short_sra
	else
		echo short: no such fastq file: $short_sra
		echo run $0 send-data-to at ${_local_host}
	fi
}

delete-short_genus_species() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local target_index="${output_dir}-${isuffix}"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local ssh_remote="${_ssh["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	if [[ -s "${short_sra}_1.fastq" ]] &&
		[[ -s "${short_sra}_2.fastq" ]]; then
		echo deleteing short-read data: $short_sra
		rm -f "${short_sra}_1.fastq"
		rm -f "${short_sra}_2.fastq"
	else
		echo short: no such fastq file: $short_sra
		echo run $0 send-data-to at ${_local_host}
	fi
}

getorganelle_genus_species() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local target_index="${output_dir}-${isuffix}"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local ssh_remote="${_ssh["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	if [[ "${output_dir}" == "default" ]]; then
		for _v1 in "${Sall[@]}"; do
			if [[ -d "${_v1}" ]]; then
				getorganelle_genus_species_for "${_v1}" "${isuffix}"
			else
				echo "No such folder: ${_v1}"
			fi
		done
	else
		getorganelle_genus_species_for "${output_dir}" "${isuffix}"
	fi
}

ptgaul_genus_species() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local target_index="${output_dir}-${isuffix}"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local ssh_remote="${_ssh["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"
	local extracted_ptgaul_genomesize="${_ptgaul_genomesize["$target_index"]}"

	mkdir -p "${output_dir}/timing"
	command time -v bash ${script_dir}/ptGAUL1.sh \
		-o ${output_dir}-ptgaul \
		-r ${output_dir}/ptdna-reference.fa \
		-g "${extracted_ptgaul_genomesize}" \
		-l ${long_sra}.fastq \
		-t 24 \
		2>${output_dir}/timing/timing-ptgaul.txt

	mv ${output_dir}-ptgaul/result_3000 ${output_dir}/ptgaul
	rm -rf "${output_dir}-ptgaul"

	msbwt_genus_species "${output_dir}"
}

msbwt_genus_species() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local target_index="${output_dir}-${isuffix}"

	local short_sra="${_short["$target_index"]}"

	mkdir -p "${output_dir}/timing"
	echo "This polishing preparation takes long ..."
	command time -v ${_polap_cmd} prepare-polishing \
		-a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
		-o ${output_dir} \
		2>${output_dir}/timing/timing-prepare-polishing.txt
}

# Extraction of ptDNA from the assembly: a first try
extract-ptdna-of-ptgaul_genus_species() {
	local output_dir="$1"

	mkdir -p "${output_dir}/timing"
	# extract ptGAUL result
	echo "extract ptDNA from the ptGAUL result with fmlrc polishing"
	command time -v ${_polap_cmd} disassemble ptgaul \
		-o ${output_dir} \
		2>${output_dir}/timing/timing-ptgaul-polishing.txt
	echo "use extract-ptdna-of-ptgaul2 <species_folder> if not working"
}

# Extraction of ptDNA from the assembly: another try
# not using it much
extract-ptdna-of-ptgaul2_genus_species() {
	local output_dir="$1"

	mkdir -p "${output_dir}/timing"
	# extract ptGAUL result
	echo "extract ptDNA from the ptGAUL result with fmlrc polishing"
	${_polap_cmd} disassemble ptgaul 2 \
		-o ${output_dir}
	command time -v ${_polap_cmd} disassemble ptgaul 3 \
		-v -v -v \
		-o ${output_dir} \
		2>${output_dir}/timing/timing-ptgaul-polishing.txt

	_outdir="${output_dir}/ptgaul/flye_cpONT/ptdna"
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

	# copy ptGAUL result
	echo "copy ${output_dir}/ptdna-ptgaul.fa"
	_outdir="${output_dir}/ptgaul/flye_cpONT/ptdna"
	_arg_final_assembly="${_outdir}/pt.1.fa"
	cp -p ${_arg_final_assembly} ${output_dir}/ptdna-ptgaul.fa
}

# Case of the infer menu
# no --disassemble-c
coverage_genus_species() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local target_index="${output_dir}-${isuffix}"

	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	local extracted_downsample="${_downsample["$target_index"]}"

	local i=1
	mkdir -p "${output_dir_i}"
	echo "Analysis (coverage): ${output_dir} at ${isuffix}"
	echo "  downsample=${extracted_downsample}x"

	local _d_i="infer-$i"
	local _x_i="coverage-$i"
	local _stages="--stages-include 0"

	# NOTE: "${_stages}" is a bug.
	# use it without quotations.
	command time -v ${_polap_cmd} disassemble \
		${_stages} \
		--downsample ${extracted_downsample} \
		-i ${extracted_inum} \
		-o ${output_dir} \
		-l ${long_sra}.fastq \
		-a ${short_sra}_1.fastq \
		-b ${short_sra}_2.fastq \
		--disassemble-i "${_d_i}" \
		--random-seed "${random_seed}" \
		2>"${output_dir_i}/timing-${_x_i}.txt"

	echo "See ${output_dir_i}/{lx.txt,sx.txt} for the basic statistics of your data"
	echo "Next: $0 infer ${output_dir} [number] to assemble the ptDNA genome"
}

delete-disassemble_genus_species() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local target_index="${output_dir}-${isuffix}"

	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	local i=1
	echo "Deleting disassemble: ${output_dir_i}/disassemble"
	read -p "Do you want to delete it? (y/N): " confirm

	case "$confirm" in
	[yY] | [yY][eE][sS])
		rm -rf "${output_dir_i}/disassemble"
		;;
	*)
		echo "Deletion canceled."
		;;
	esac

	echo "See ${output_dir_i}/{lx.txt,sx.txt} for the basic statistics of your data"
	echo "Next: $0 infer ${output_dir} [number] to assemble the ptDNA genome"
}

# Case of the infer menu
# no --disassemble-c
infer_genus_species() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local simple_polishing="${3:-default}"
	local target_index="${output_dir}-${isuffix}"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local ssh_remote="${_ssh["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	local i=0
	local n
	local p
	IFS=':' read -r -a extracted_array_n <<<"${_compare_n["$target_index"]}"
	IFS=':' read -r -a extracted_array_p <<<"${_compare_p["$target_index"]}"
	local extracted_r="${_compare_r["$target_index"]}"
	local extracted_memory="${_memory["$target_index"]}"
	local extracted_downsample="${_downsample["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	mkdir -p "${output_dir_i}"

	for n in "${extracted_array_n[@]}"; do
		for p in "${extracted_array_p[@]}"; do
			i=$((i + 1))
			echo "Analysis (inference): ${output_dir} at ${isuffix}"
			echo "($i) n=$n, p=$p, r=${extracted_r} memory=${extracted_memory}G, downsample=${extracted_downsample}x"

			local _d_i="infer-$i"
			local _x_i="infer-$i"
			local _s_i="subsample-polish"
			local _stages="--stages-include 0-3"
			if [[ "${simple_polishing}" != "default" ]]; then
				_stages="--stages-include 3"
				_s_i="simple-polish"
			fi

			# NOTE: "${_stages}" is a bug.
			# use it without quotations.
			command time -v ${_polap_cmd} disassemble \
				${_stages} \
				--downsample ${extracted_downsample} \
				-i ${extracted_inum} \
				-o ${output_dir} \
				-l ${long_sra}.fastq \
				-a ${short_sra}_1.fastq \
				-b ${short_sra}_2.fastq \
				${simple_polishing} \
				--disassemble-i "${_d_i}" \
				--disassemble-n $n \
				--disassemble-p $p \
				--disassemble-r ${extracted_r} \
				--disassemble-memory ${extracted_memory} \
				--disassemble-alpha 1.0 \
				--random-seed "${random_seed}" \
				2>"${output_dir_i}/timing-${_x_i}-${_s_i}.txt"

			if [[ -d "${output_dir_i}/disassemble/${_d_i}/3" ]]; then
				rm -rf "${output_dir_i}/disassemble/${_d_i}/3-infer"
				mv "${output_dir_i}/disassemble/${_d_i}/3" \
					"${output_dir_i}/disassemble/${_d_i}/3-infer"
			fi
		done
	done
	echo "Next: $0 check ${output_dir} [number] to compare the assembly of the ptDNA genome with ptGAUL's"
}

infer2_genus_species() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local simple_polishing="${3:-default}"
	local target_index="${output_dir}-${isuffix}"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local ssh_remote="${_ssh["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	local i=0
	local n
	local p
	IFS=':' read -r -a extracted_array_n <<<"${_compare_n["$target_index"]}"
	IFS=':' read -r -a extracted_array_p <<<"${_compare_p["$target_index"]}"
	local extracted_r="${_compare_r["$target_index"]}"
	local extracted_memory="${_memory["$target_index"]}"
	local extracted_downsample="${_downsample["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	mkdir -p "${output_dir_i}"

	for n in "${extracted_array_n[@]}"; do
		for p in "${extracted_array_p[@]}"; do
			i=$((i + 2))
			echo "Analysis (inference): ${output_dir} at ${isuffix}"
			echo "($i) n=$n, p=$p, r=${extracted_r} memory=${extracted_memory}G, downsample=${extracted_downsample}x"

			local _d_i="infer-$i"
			local _x_i="infer-$i"
			local _s_i="subsample-polish"
			local _stages="--stages-include 0-3"
			if [[ "${simple_polishing}" != "default" ]]; then
				_stages="--stages-include 3"
				_s_i="simple-polish"
			fi

			# NOTE: "${_stages}" is a bug.
			# use it without quotations.
			command time -v ${_polap_cmd} disassemble \
				${_stages} \
				--downsample ${extracted_downsample} \
				-i ${extracted_inum} \
				-o ${output_dir} \
				-l ${long_sra}.fastq \
				-a ${short_sra}_1.fastq \
				-b ${short_sra}_2.fastq \
				${simple_polishing} \
				--disassemble-i "${_d_i}" \
				--disassemble-n $n \
				--disassemble-p $p \
				--disassemble-r ${extracted_r} \
				--disassemble-memory ${extracted_memory} \
				--disassemble-alpha 1.5 \
				--disassemble-beta 0.05 \
				--random-seed "${random_seed}" \
				2>"${output_dir_i}/timing-${_x_i}-${_s_i}.txt"

			if [[ -d "${output_dir_i}/disassemble/${_d_i}/3" ]]; then
				rm -rf "${output_dir_i}/disassemble/${_d_i}/3-infer"
				mv "${output_dir_i}/disassemble/${_d_i}/3" \
					"${output_dir_i}/disassemble/${_d_i}/3-infer"
			fi
		done
	done
	echo "Next: $0 check ${output_dir} [number] to compare the assembly of the ptDNA genome with ptGAUL's"
}

# Case of the check menu
# --disassemble-c
# --disassemble-align-reference
# --disassemble-simple-polishing
check_genus_species() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local icount="${3:-0}"
	local simple_polishing="${3:-default}"
	local target_index="${output_dir}-${isuffix}"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local ssh_remote="${_ssh["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	local i="${icount}"
	local n
	local p
	IFS=':' read -r -a extracted_array_n <<<"${_compare_n["$target_index"]}"
	IFS=':' read -r -a extracted_array_p <<<"${_compare_p["$target_index"]}"
	local extracted_r="${_compare_r["$target_index"]}"
	local extracted_memory="${_memory["$target_index"]}"
	local extracted_downsample="${_downsample["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	mkdir -p "${output_dir_i}"

	for n in "${extracted_array_n[@]}"; do
		for p in "${extracted_array_p[@]}"; do
			i=$((i + 1))
			echo "Analysis (check): ${output_dir} at ${isuffix}"
			echo "($i) n=$n, p=$p, memory=${extracted_memory}G, downsample=${extracted_downsample}x"

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
				--downsample ${extracted_downsample} \
				-i ${extracted_inum} \
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
				--disassemble-memory ${extracted_memory} \
				--disassemble-alpha 1.0 \
				--random-seed "${random_seed}" \
				2>"${output_dir_i}/timing-${_x_i}-${_s_i}.txt"

			# compare the results

			if [[ -d "${output_dir_i}/disassemble/${_d_i}/3" ]]; then
				rm -rf "${output_dir_i}/disassemble/${_d_i}/3-check"
				mv "${output_dir_i}/disassemble/${_d_i}/3" \
					"${output_dir_i}/disassemble/${_d_i}/3-check"
			fi

			local mauve_dir="${output_dir_i}/mauve/${i}"
			local blast_dir="${output_dir_i}/blast/${i}"
			local mafft_dir="${output_dir_i}/mafft/${i}"
			mkdir -p "${mauve_dir}"
			mkdir -p "${blast_dir}"
			mkdir -p "${mafft_dir}"
			if [[ -s "${output_dir_i}/disassemble/${_d_i}/pt.subsample-polishing.reference.aligned.1.fa" ]]; then

				${_polap_cmd} mauve-mtdna -a "${output_dir}/ptdna-ptgaul.fa" \
					-b "${output_dir_i}/disassemble/${_d_i}/pt.subsample-polishing.reference.aligned.1.fa" \
					-o "${mauve_dir}" \
					>"${mauve_dir}/log.txt"
				echo "see ${mauve_dir}/log.txt"
				cat "${mauve_dir}/log.txt"

				${_polap_cmd} compare2ptdna -a "${output_dir}/ptdna-ptgaul.fa" \
					-b "${output_dir_i}/disassemble/${_d_i}/pt.subsample-polishing.reference.aligned.1.fa" \
					-o "${blast_dir}"
				echo "see ${blast_dir}/pident.txt"
				cat "${blast_dir}/pident.txt"

				${_polap_cmd} mafft-mtdna -a "${output_dir}/ptdna-ptgaul.fa" \
					-b "${output_dir_i}/disassemble/${_d_i}/pt.subsample-polishing.reference.aligned.1.fa" \
					-o "${mafft_dir}" \
					>"${mafft_dir}/log.txt"
				echo "see ${mafft_dir}/pident.txt"
				cat "${mafft_dir}/pident.txt"
			else
				echo "ERROR: no such file: ${output_dir_i}/disassemble/${_d_i}/pt.subsample-polishing.reference.aligned.1.fa"
			fi

		done
	done
}

# Case of the compare menu
# --disassemble-c
# --no-disassemble-align-reference
compare_genus_species() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local simple_polishing="${3:-default}"
	local target_index="${output_dir}-${isuffix}"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local ssh_remote="${_ssh["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	local i=0
	local n
	local p
	IFS=':' read -r -a extracted_array_n <<<"${_compare_n["$target_index"]}"
	IFS=':' read -r -a extracted_array_p <<<"${_compare_p["$target_index"]}"
	local extracted_r="${_compare_r["$target_index"]}"
	local extracted_memory="${_memory["$target_index"]}"
	local extracted_downsample="${_downsample["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	mkdir -p "${output_dir_i}"

	for n in "${extracted_array_n[@]}"; do
		for p in "${extracted_array_p[@]}"; do
			i=$((i + 1))
			echo "Analysis (check stage 1 or compare): ${output_dir} at ${isuffix}"
			echo "($i) n=$n, p=$p, memory=${extracted_memory}G, downsample=${extracted_downsample}x"

			local _d_i
			_d_i="compare-$i"
			local _stages="--stages-include 1-3" # to avoid another down-sampling

			command time -v ${_polap_cmd} disassemble \
				${_stages} \
				--downsample ${extracted_downsample} \
				-i ${extracted_inum} \
				-o ${output_dir} \
				-l ${long_sra}.fastq \
				-a ${short_sra}_1.fastq \
				-b ${short_sra}_2.fastq \
				--disassemble-c ${output_dir}/ptdna-ptgaul.fa \
				--disassemble-i "${_d_i}" \
				--disassemble-n $n \
				--disassemble-p $p \
				--disassemble-r ${extracted_r} \
				--disassemble-memory ${extracted_memory} \
				--disassemble-alpha 1.0 \
				--random-seed "${random_seed}" \
				2>${output_dir_i}/timing-${_d_i}.txt
		done
	done
}

best-genus_species() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local simple_polishing="${3:-default}"
	local target_index="${output_dir}-${isuffix}"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local ssh_remote="${_ssh["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	${_polap_cmd} disassemble best \
		-i ${extracted_inum} \
		-o ${output_dir} \
		--disassemble-i infer-1

}

bandage1_genus_species() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local simple_polishing="${3:-default}"
	local target_index="${output_dir}-${isuffix}"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local ssh_remote="${_ssh["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"
	local images=()
	local captions=()

	${_polap_cmd} disassemble bandage \
		-i ${extracted_inum} \
		-o ${output_dir} \
		--disassemble-i infer-1
	local _gfa_polap="${output_dir}/${isuffix}/disassemble/infer-1/pt.1.gfa"
	local _png_polap="${output_dir}/${isuffix}/disassemble/infer-1/pt.1.png"

	images+=(${_png_polap})
	captions+=("Polap")

	local _gfa_ptgaul="${output_dir}/ptgaul/flye_cpONT/assembly_graph.gfa"
	local _png_ptgaul="${output_dir}/ptgaul.png"
	${_polap_cmd} bandage png \
		${_gfa_ptgaul} \
		${_png_ptgaul}
	images+=(${_png_ptgaul})
	captions+=("ptGAUL")

	echo "ptGAUL ptDNA gfa: ${_gfa_ptgaul}"
	echo "ptGAUL ptDNA png: ${_png_ptgaul}"

	local _gfa_getorganelle=$(find "${output_dir}/getorganelle" -type f -name 'embplant_pt*.gfa' | head -n 1)
	local _png_getorganelle="${output_dir}/getorganelle.png"
	${_polap_cmd} bandage png \
		${_gfa_getorganelle} \
		${_png_getorganelle}
	images+=(${_png_getorganelle})
	captions+=("GetOrganelle")

	echo "GetOrganelle ptDNA gfa: ${_gfa_getorganelle}"
	echo "GetOrganelle ptDNA png: ${_png_getorganelle}"

	local _figure_md="figure1.md"

	# echo "Warning: no such file: ${_gfa}"
	# echo "FIX QT5 problem:"
	# problem solved
	# echo "export QT_QPA_PLATFORM=offscreen"
	# Output Markdown file
	#
	output="${_figure_md}"

	# Start writing to the markdown file
	cat <<EOF >"$output"
---
title: "${species_name}"
geometry: margin=1in
---

# Image Grid

| Three Organelle genome assemblies |  |  |
|-----------------|-----------------|-----------------|
EOF

	# Generate the image table with subcaptions
	count=0
	row="| "
	caption_row="| "

	for ((i = 0; i < ${#images[@]}; i++)); do
		row+="![${captions[i]}](figures/${images[i]}){width=100px} | "
		caption_row+="**${captions[i]}** | "
		((count++))

		# End row if 3 images are added
		if ((count % 3 == 0)); then
			echo "$row" >>"$output"
			echo "|-----------------|-----------------|-----------------|" >>"$output"
			echo "$caption_row" >>"$output"
			# echo "" >>"$output"
			row="| "
			caption_row="| "
		fi
	done
}

# figure1's all bandage assembly_graph
# figure2.md
bandage2_genus_species() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local simple_polishing="${3:-default}"
	local target_index="${output_dir}-${isuffix}"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local ssh_remote="${_ssh["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"
	local extracted_n="${_compare_n["$target_index"]}"
	local i
	local width=13
	local images=()
	local captions=()

	local _figure2_md="figure2.md"
	local k=0
	for ((i = 0; i < extracted_n; i++)); do

		local _gfa_infer="${output_dir}/${extracted_inum}/disassemble/infer-1/1/${i}/30-contigger/graph_final.gfa"
		local _png_infer="${output_dir}/${extracted_inum}/disassemble/infer-1/1/${i}/30-contigger/graph_final.png"
		if [[ -s "${_gfa_infer}" ]]; then
			echo "gfa file: ${_gfa_infer}"
			${_polap_cmd} bandage png \
				${_gfa_infer} \
				${_png_infer}
			printf "| ![polap %s](figures/%s){ width=%s%% } " "${i}" "${_png_infer}" "${width}" >>"${_figure2_md}"
			images+=("figures/${_png_infer}")
			captions+=(${i})
			((k++))
		else
			echo "no such file: ${_gfa_infer}"
		fi
	done

	printf "Figure. Polap inference - %s\n" "${species_name}" >>"${_figure2_md}"

	#!/bin/bash

	# Define an array of PNG files and their subcaptions
	# images=("fig1.png" "fig2.png" "fig3.png" "fig4.png" "fig5.png" "fig6.png" "fig7.png" "fig8.png")
	# captions=("Subcaption 1" "Subcaption 2" "Subcaption 3" "Subcaption 4" "Subcaption 5" "Subcaption 6" "Subcaption 7" "Subcaption 8")

	# Output Markdown file
	output="${_figure2_md}"

	# Start writing to the markdown file
	cat <<EOF >"$output"
---
title: "${species_name}"
geometry: margin=1in
---

# Image Grid

| Organelle genome assemblies |  |  |
|-----------------|-----------------|-----------------|
EOF

	# Generate the image table with subcaptions
	count=0
	row="| "
	caption_row="| "

	for ((i = 0; i < ${#images[@]}; i++)); do
		row+="![${captions[i]}](${images[i]}){width=100px} | "
		caption_row+="**${captions[i]}** | "
		((count++))

		# End row if 3 images are added
		if ((count % 3 == 0)); then
			echo "$row" >>"$output"
			echo "|-----------------|-----------------|-----------------|" >>"$output"
			echo "$caption_row" >>"$output"
			# echo "" >>"$output"
			row="| "
			caption_row="| "
		fi
	done

	# Handle the last incomplete row (if any)
	remaining=$((3 - count % 3))
	if ((remaining < 3)); then
		# Fill empty image columns
		for ((i = 0; i < remaining; i++)); do
			row+=" | "
			caption_row+=" | "
		done
		echo "$row" >>"$output"
		echo "|-----------------|-----------------|-----------------|" >>"$output"
		echo "$caption_row" >>"$output"
	fi

	echo "" >>"${output}"
	printf "Figure. Polap inference - %s\n" "${species_name}" >>"${output}"

}

copy-figure_genus_species() {
	rsync -av --include='*/' --include='*.png' --exclude='*' ./ ../../manuscript/polap-v0.4/figures/
	rsync -av --include='*/' --include='*.pdf' --exclude='*' ./ ../../manuscript/polap-v0.4/figures/
}

# downsample <folder> <isuffix> [--dry]
# e.g., downsample Pisum_sativum 0 --dry
downsample_genus_species() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local _dry="${3:-default_value}"
	local target_index="${output_dir}-${isuffix}"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local ssh_remote="${_ssh["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"
	local extracted_downsample="${_downsample["$target_index"]}"

	_c=${extracted_downsample}
	if [[ "${_dry}" == "default_value" ]]; then
		_dry=""
	fi

	local _d_i
	_d_i="downsample"

	mkdir -p "${output_dir}"

	local long_data="${_media_dir}/${long_sra}.fastq.tar.gz"
	local short_data="${_media_dir}/${short_sra}.fastq.tar.gz"

	if [[ "${_local_host}" == "$(hostname)" ]]; then
		if [[ -s "${long_data}" ]]; then
			tar -zxf ${long_data}
			tar -zxf ${short_data}
		else
			ln -s "${_media1_dir}/${long_sra}.fastq"
			ln -s "${_media1_dir}/${short_sra}_1.fastq"
			ln -s "${_media1_dir}/${short_sra}_2.fastq"
		fi
	# else
	# 	echo "You can send the data from the central host: ${_local_host}"
	# 	return 1
	fi

	if [[ -s "$output_dir/short_expected_genome_size.txt" ]]; then
		echo "Found: $output_dir/short_expected_genome_size.txt"
	else
		${_polap_cmd} find-genome-size \
			-a ${short_sra}_1.fastq \
			-b ${short_sra}_2.fastq \
			-o "${output_dir}"
	fi
	local _genome_size=$(<"$output_dir/short_expected_genome_size.txt")
	echo "Genome size estimate: ${_genome_size}"

	${_polap_cmd} fastq subsample -v ${_dry} \
		"${long_sra}.fastq" \
		"${output_dir}/l${_c}x.fq" \
		-c "${_c}" \
		-o "${output_dir}" \
		--random-seed "${random_seed}" \
		--genomesize "${_genome_size}" -v \
		>"${output_dir}/l${_c}x.txt"

	${_polap_cmd} fastq subsample2 -v ${_dry} \
		"${short_sra}_1.fastq" \
		"${short_sra}_2.fastq" \
		"${output_dir}/s${_c}x_1.fq" \
		"${output_dir}/s${_c}x_2.fq" \
		-c "${_c}" \
		-o "${output_dir}" \
		--random-seed "${random_seed}" \
		--genomesize "${_genome_size}" -v \
		>"${output_dir}/s${_c}x.txt"

}

# use-downsample <folder> [inum]
# replace the input data with the downsampled data
use-downsample_genus_species() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local target_index="${output_dir}-${isuffix}"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local ssh_remote="${_ssh["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"
	local extracted_downsample="${_downsample["$target_index"]}"

	_c=${extracted_downsample}

	# mv "${long_sra}.fastq" "${output_dir}/l.fq"
	rm -f "${long_sra}.fastq"
	# gunzip -f "${output_dir}/l${_c}x.fq.gz"
	mv "${output_dir}/l${_c}x.fq" "${long_sra}.fastq"

	# mv "${short_sra}_1.fastq" "${output_dir}/s_1.fq"
	rm -f "${short_sra}_1.fastq"
	# gunzip -f "${output_dir}/s${_c}x_1.fq.gz"
	mv "${output_dir}/s${_c}x_1.fq" "${short_sra}_1.fastq"

	# mv "${short_sra}_2.fastq" "${output_dir}/s_2.fq"
	rm -f "${short_sra}_2.fastq"
	# gunzip -f "${output_dir}/s${_c}x_2.fq.gz"
	mv "${output_dir}/s${_c}x_2.fq" "${short_sra}_2.fastq"
}

mauve_genus_species() {
	local output_dir="${1:-default}"
	local isuffix="${2:-0}"
	local target_index="${output_dir}-${isuffix}"

	if [[ "${output_dir}" == "default" ]]; then

		# Extract and sort keys
		sorted_keys=($(for key in "${!_table1[@]}"; do echo "$key"; done | sort))

		# Iterate over sorted keys and check if value is "T"
		for key in "${sorted_keys[@]}"; do
			if [[ "${_table1[$key]}" == "F" ]]; then
				continue
			fi
			local _v1=${_folder[$key]}

			mauve_genus_species_for "${_v1}" "${isuffix}"
		done
	else
		mauve_genus_species_for "${output_dir}" "${isuffix}"
	fi
}

mauve_genus_species_for() {
	local output_dir="$1"
	local isuffix="${2:-0}"
	local target_index="${output_dir}-${isuffix}"

	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	local i=0
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	mkdir -p "${output_dir_i}"

	i=$((i + 1))
	echo "($i)"

	local mauve_dir="${output_dir_i}/mauve/${i}"
	local mafft_dir="${output_dir_i}/mafft/${i}"
	local blast_dir="${output_dir_i}/blast/${i}"
	mkdir -p "${mauve_dir}"
	mkdir -p "${mafft_dir}"
	mkdir -p "${blast_dir}"
	if [[ -s "${output_dir_i}/disassemble/infer-${i}/pt.subsample-polishing.reference.aligned.1.fa" ]]; then
		${_polap_cmd} mafft-mtdna -a "${output_dir}/ptdna-ptgaul.fa" \
			-b "${output_dir_i}/disassemble/infer-${i}/pt.subsample-polishing.reference.aligned.1.fa" \
			-o "${mafft_dir}" \
			>"${mafft_dir}/log.txt"
		echo "see ${mafft_dir}/pident.txt"
		cat "${mafft_dir}/pident.txt"

		# ${_polap_cmd} mauve-mtdna -a "${output_dir}/ptdna-ptgaul.fa" \
		# 	-b "${output_dir_i}/disassemble/infer-${i}/pt.subsample-polishing.reference.aligned.1.fa" \
		# 	-o "${mauve_dir}" \
		# 	>"${mauve_dir}/log.txt"
		# echo "see ${mauve_dir}/log.txt"
		# cat "${mauve_dir}/log.txt"
		#
		# ${_polap_cmd} compare2ptdna -a "${output_dir}/ptdna-ptgaul.fa" \
		# 	-b "${output_dir_i}/disassemble/infer-${i}/pt.subsample-polishing.reference.aligned.1.fa" \
		# 	-o "${blast_dir}"
		# echo "see ${blast_dir}/pident.txt"
		# cat "${blast_dir}/pident.txt"
	else
		echo "ERROR: no such file: ${output_dir_i}/disassemble/infer-${i}/pt.subsample-polishing.reference.aligned.1.fa"
	fi

}

restart_genus_species() {
	local output_dir="$1"

	local target_dir="${output_dir}-t"
	mkdir "${target_dir}"

	for i in \
		00-bioproject \
		0-bioproject \
		getorganelle \
		ptgaul \
		msbwt \
		timing \
		polap.log; do
		mv "${output_dir}/$i" "${target_dir}"
	done
	mv "${output_dir}"/*.fa "${target_dir}"

	mv "${output_dir}" "${output_dir}-backup"
	mv "${target_dir}" "${output_dir}"
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
	# copy_data

	rm -rf "${output_dir}-a"
	rm -f "${output_dir}-a.tar.gz"
	${_polap_cmd} disassemble archive \
		--max-filesize 5M \
		-o ${output_dir}
	# tar zcf "${output_dir}-a.tar.gz" "${output_dir}-a"
}

get_genus_species() {
	local output_dir="${1:-default}"
	local isuffix="${2:-0}"
	local target_index="${output_dir}-${isuffix}"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local ssh_remote="${_ssh["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	if [[ "${output_dir}" == "default" ]]; then
		for _v1 in "${S[@]}"; do
			get_genus_species_for "${_v1}" "${isuffix}"
		done
	else
		get_genus_species_for "${output_dir}" "${isuffix}"
	fi
}

get_genus_species_for() {
	local output_dir="${1:-default}"
	local isuffix="${2:-0}"
	local target_index="${output_dir}-${isuffix}"

	local species_name="$(echo ${output_dir} | sed 's/_/ /')"
	local long_sra="${_long["$target_index"]}"
	local short_sra="${_short["$target_index"]}"
	local random_seed="${_random_seed["$target_index"]}"
	local ssh_remote="${_ssh["$target_index"]}"
	local extracted_inum="${_inum["$target_index"]}"
	local output_dir_i="${output_dir}/${extracted_inum}"

	if [[ "${_local_host}" == "$(hostname)" ]]; then

		# Define the directory name (can be passed as an argument)
		dir_to_check="${output_dir}"

		# Check if the directory exists
		if [[ -d "$dir_to_check" ]]; then
			echo "Directory '$dir_to_check' exists."
			# Ask for user confirmation before deleting
			read -p "Do you want to delete it? (y/N): " confirm
		else
			echo "Directory '$dir_to_check' does not exist."
			confirm="yes"
		fi

		case "$confirm" in
		[yY] | [yY][eE][sS])

			rm -f "${output_dir}-a.tar.gz"
			scp -p ${ssh_remote}:$PWD/${output_dir}-a.tar.gz .

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
		echo "ERROR: run at the local host."
	fi
}

# TODO
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

# TODO
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

################################################################################
# create tables and figures
# some helper functions first and then the main table and figure functions
#

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
# read -r I P N R <<< "$result"
parse_params() {
	local file="$1" # Input file
	local I P N R   # Declare variables for the parameters

	# Read the file and extract the values
	while IFS=": " read -r key value; do
		case "$key" in
		"I") I="$value" ;;
		"P") P="$value" ;;
		"N") N="$value" ;;
		"R") R="$value" ;;
		esac
	done <"$file"

	# Print the variables (return as output)
	echo "$I $P $N $R"
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

# extract short-read coverage
get_short_read_coverage() {
	local file="$1"
	if [[ ! -f "$file" ]]; then
		echo "Error: File not found!"
		return 1
	fi

	local coverage
	# coverage=$(grep -oP 'short-read coverage:\s*\K\d+' "$file")
	coverage=$(grep -oP 'short-read coverage:\s*\K[\d.]+(?=x)' "$file" | tail -1)

	if [[ -n "$coverage" ]]; then
		coverage=$(echo "$coverage" | tr -d '[:space:]') # Remove spaces
		printf "%.1f\n" "$coverage"
	else
		echo "Error: Coverage information not found!"
		return 1
	fi
}

get_long_read_coverage() {
	local file="$1"
	if [[ ! -f "$file" ]]; then
		echo "Error: File not found!"
		return 1
	fi

	local coverage
	# coverage=$(grep -oP 'long-read coverage:\s*\K\d+' "$file")
	coverage=$(grep -oP 'long-read coverage:\s*\K[\d.]+(?=x)' "$file" | tail -1)

	if [[ -n "$coverage" ]]; then
		coverage=$(echo "$coverage" | tr -d '[:space:]') # Remove spaces
		printf "%.1f\n" "$coverage"
	else
		echo "Error: Coverage information not found!"
		return 1
	fi
}

get_target_read_coverage() {
	local file="$1"
	if [[ ! -f "$file" ]]; then
		echo "Error: File not found!"
		return 1
	fi

	local coverage
	# coverage=$(grep -oP 'target coverage:\s*\K\d+' "$file")
	coverage=$(grep -oP 'target coverage:\s*\K[\d.]+(?=x)' "$file" | tail -1)

	if [[ -n "$coverage" ]]; then
		coverage=$(echo "$coverage" | tr -d '[:space:]') # Remove spaces
		printf "%.1f\n" "$coverage"
	else
		echo "Error: Coverage information not found!"
		return 1
	fi
}

#
# table1 2
#
# output:
# table1.md
# table1.tsv
# table1-1.md
# table1-2.md
table1_genus_species() {
	local _arg_inum="${1:-2}"
	local _table="table1"
	local _table_file="${_table}.tsv"
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
	local _memory_gb_getorganelle
	local _total_hours_getorganelle
	local _memory_gb_flye1
	local _total_hours_flye1
	local _I
	local _P
	local _N
	local j

	printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
		"Species" \
		"Order" \
		"Family" \
		"I" \
		"C" \
		"L_SRA" \
		"L_size" \
		"L_cov" \
		"S_SRA" \
		"S_size" \
		"S_cov" \
		"G" \
		"N" \
		"P" \
		"R" \
		"Rate" \
		"Size" \
		"Alpha" \
		"ptDNA" \
		"Length1" \
		"Length2" \
		"Pident" \
		"N1" \
		"Mode" \
		"SD" \
		"N2" \
		"M" \
		"M_g" \
		"M_p" \
		"M_t1" \
		"M_t2" \
		"M_s" \
		"M_f" \
		"T" \
		>"${_table_file}"

	# Extract and sort keys
	sorted_keys=($(for key in "${!_table1[@]}"; do echo "$key"; done | sort))

	# Iterate over sorted keys and check if value is "T"
	for key in "${sorted_keys[@]}"; do
		if [[ "${_table1[$key]}" == "F" ]]; then
			continue
		fi
		local _v1=${_folder[$key]}
		local _extracted_inum=${_inum[$key]}
		if [[ "${_extracted_inum}" != "${_arg_inum}" ]]; then
			continue
		fi

		# done

		_logfile=${_v1}/polap.log
		_readmefile=${_v1}/README
		_species="${_v1//_/ }"
		_genus=${_species%% *}
		_order=$(grep "${_genus}" taxonomy_output.tsv | cut -f 6 | head -1)
		_family=$(grep "${_genus}" taxonomy_output.tsv | cut -f 7 | head -1)
		echo "${_species}: ${_order}: ${_family}"

		_l_sra=$(awk '/long-read/ {match($0, /([A-Z]RR[0-9]+)/, arr); print arr[0]}' "$_logfile" | sort -u | grep -v '^$')
		if [[ -s "${_v1}/long_total_length.txt" ]]; then
			_l_sra_size=$(<"${_v1}/long_total_length.txt")
		else
			_l_sra_size=$(<"${_v1}/l.fq.txt")
		fi
		_l_sra_size_gb=$(_polap_utility_convert_bp "${_l_sra_size}")

		_s_sra=$(awk '/short-read1/ {match($0, /([A-Z]RR[0-9]+)/, arr); print arr[0]}' "$_logfile" | sort -u | grep -v '^$')
		if [[ -s "${_v1}/short_total_length.txt" ]]; then
			_s_sra_size=$(<"${_v1}/short_total_length.txt")
		else
			_s_sra_size1=$(<"${_v1}/s1.fq.txt")
			_s_sra_size2=$(<"${_v1}/s2.fq.txt")
			_s_sra_size=$((_s_sra_size1 + _s_sra_size2))
		fi
		_s_sra_size_gb=$(_polap_utility_convert_bp "${_s_sra_size}")
		_known_mtdna=$(grep 'NCBI accession:' ${_logfile} | cut -d: -f4 | tail -n 1)

		if [[ -s "${_v1}/short_expected_genome_size.txt" ]]; then
			_genome_size=$(<"${_v1}/short_expected_genome_size.txt")
		else
			echo "ERROR: no such file: ${_v1}/short_expected_genome_size.txt"
			exit 1
		fi

		# if [[ -v _host["${_v1}"] ]]; then
		# 	echo species:"${_v1}"
		# else
		# 	echo "No such host for $_v1"
		# fi

		local _ptdna_ptgaul=${_v1}/ptdna-ptgaul.fa
		local _ptdna_reference=${_v1}/ptdna-reference.fa
		seq_length1=$(grep -v "^>" "$_ptdna_ptgaul" | tr -d '\n' | wc -c)

		# read -r _memory_gb_flye1 _total_hours_flye1 < <(parse_timing "${_v1}" "infer12-4")

		local j=1
		local _i="${_extracted_inum}"
		local _v1_inum="${_v1}/${_i}"
		if [[ ! -d "${_v1_inum}" ]]; then
			continue
		fi
		local _disassemble_index="infer-${j}"
		local _summary1_ordered_txt=${_v1_inum}/disassemble/${_disassemble_index}/1/summary1-ordered.txt
		local fasta_file=${_v1_inum}/disassemble/${_disassemble_index}/pt.subsample-polishing.1.fa

		local target_index="${_v1}-${_i}"
		local extracted_memory="${_memory["$target_index"]}"

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
		local _n1=$(grep "^#n:" "$_summary1_ordered_txt" | awk 'NR==1 {print $2}')

		local _params_txt=${_v1_inum}/disassemble/${_disassemble_index}/params.txt
		_ipn=$(parse_params "${_params_txt}")
		read -r _I _P _N _R <<<"$_ipn"

		local _summary2_ordered_txt=${_v1_inum}/disassemble/${_disassemble_index}/2/summary1-ordered.txt
		local _n2=$(grep "^#n:" "$_summary2_ordered_txt" | awk 'NR==1 {print $2}')
		local output=$(awk -F'\t' 'NR==2 {print $1, $2, $4, $11}' "${_summary2_ordered_txt}")
		local _summary2_index
		local _summary2_size
		local _summary2_size_gb
		local _summary2_rate_rounded
		local _summary2_alpha
		local _summary2_alpha_formatted
		read -r _summary2_index _summary2_size _summary2_rate _summary2_alpha <<<"$output"
		_summary2_size_gb=$(_polap_utility_convert_bp "${_summary2_size}")
		local _summary2_rate_decimal=$(printf "%.10f" "$_summary2_rate")
		_summary2_rate_rounded=$(echo "scale=4; $_summary2_rate_decimal / 1" | bc)
		_summary2_alpha_formatted=$(echo "scale=2; $_summary2_alpha / 1" | bc | awk '{printf "%.2f\n", $1}')
		# _summary2_alpha_formatted=$(echo "$_summary2_alpha_formatted" | awk '{printf "%.10g\n", $1}')

		read -r _memory_gb_getorganelle _total_hours_getorganelle < <(_polap_lib_timing-parse-timing "${_v1}/timing/timing-getorganelle.txt")
		read -r _memory_gb_ptgaul _total_hours_ptgaul < <(_polap_lib_timing-parse-timing "${_v1}/timing/timing-ptgaul.txt")
		read -r _memory_gb_prepare_polishing _total_hours_prepare_polishing < <(_polap_lib_timing-parse-timing "${_v1}/timing/timing-prepare-polishing.txt")
		read -r _memory_gb_ptgaul_polishing _total_hours_ptgaul_polishing < <(_polap_lib_timing-parse-timing "${_v1}/timing/timing-ptgaul-polishing.txt")
		read -r _memory_gb_subsampling_polishing _total_hours_subsampling_polishing < <(_polap_lib_timing-parse-timing "${_v1_inum}/timing-check-${j}-subsample-polish.txt")
		read -r _memory_gb _total_hours < <(_polap_lib_timing-parse-timing "${_v1_inum}/timing-infer-${j}-subsample-polish.txt")

		local _mafft_pident="NA"
		if [[ -s "${_v1_inum}/mafft/${j}/pident.txt" ]]; then
			_mafft_pident="$(<${_v1_inum}/mafft/${j}/pident.txt)"
		else
			_mafft_pident=0
		fi

		local _blast_pident="NA"
		if [[ -s "${_v1_inum}/blast/${j}/pident.txt" ]]; then
			_blast_pident="$(<${_v1_inum}/blast/${j}/pident.txt)"
		else
			_blast_pident=0
		fi

		local _mauve_lcb_coverage="NA"
		if [[ -s "${_v1_inum}/mauve/${j}/log.txt" ]]; then
			_mauve_lcb_coverage=$(awk '{print $2}' "${_v1_inum}/mauve/${j}/log.txt")
		else
			_mauve_lcb_coverage=0
		fi

		local _pident="${_mafft_pident}"
		# if (($(echo "$_blast_pident < $_mauve_lcb_coverage" | bc -l))); then
		# 	_pident="${_mauve_lcb_coverage}"
		# else
		# 	_pident="${_blast_pident}"
		# fi

		local _short_coverage=$(get_short_read_coverage "${_v1_inum}/sx.txt")
		local _long_coverage=$(get_long_read_coverage "${_v1_inum}/lx.txt")
		local _target_coverage=$(get_target_read_coverage "${_v1_inum}/lx.txt")

		# read -r _memory_gb _total_hours < <(parse_timing "${_v1}" "infer12-${j}")
		# read -r _memory_gb_polishing _total_hours_polishing < <(parse_timing "${_v1}" "infer3only-${j}")

		printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
			"_${_species}_" \
			"${_order}" \
			"${_family}" \
			"${_i}" \
			"${_target_coverage}" \
			"${_l_sra}" \
			"${_l_sra_size_gb}" \
			"${_long_coverage}" \
			"${_s_sra}" \
			"${_s_sra_size_gb}" \
			"${_short_coverage}" \
			"${_genome_size}" \
			"${_N}" \
			"${_P}" \
			"${_R}" \
			"${_summary2_rate_rounded}" \
			"${_summary2_size_gb}" \
			"${_summary2_alpha_formatted}" \
			"${_known_mtdna}" \
			"${seq_length1}" \
			"${seq_length2}" \
			"${_pident}" \
			"${_n1}" \
			"${_mode}" \
			"${_sd}" \
			"${_n2}" \
			"${extracted_memory}" \
			"${_memory_gb_getorganelle}" \
			"${_memory_gb_ptgaul}" \
			"${_memory_gb_prepare_polishing}" \
			"${_memory_gb_ptgaul_polishing}" \
			"${_memory_gb_subsampling_polishing}" \
			"${_memory_gb}" \
			"${_total_hours}" \
			>>"${_table_file}"
	done

	# numfmt --field=7,8 --header=1 --delimiter=\t --grouping >${_table}-1.md
	csvtk -t cut -f Species,C,N,P,R,Rate,Alpha,Length1,Length2,Pident,N1,Mode,SD,M,M_g,M_p,M_t1,M_t2,M_s,M_f,T \
		${_table_file} |
		csvtk -t rename -f 1-21 -n Species,C,N,P,R,Rate,Alpha,L1,L2,Pident,N1,Mode,SD,M,Mg,Mp,Mt1,Mt2,Ms,Mf,T |
		csvtk -t csv2md -a right -o ${_table}-analysis-${_arg_inum}.md

	csvtk -t cut -f Species,Order,Family,L_SRA,L_size,L_cov,S_SRA,S_size,S_cov ${_table_file} |
		csvtk -t rename -f 1-9 -n Species,Order,Family,L_SRA,L_size,L_cov,S_SRA,S_size,S_cov |
		csvtk -t csv2md -a right -o ${_table}-data.md

	cp -p ${_table}-analysis-${_arg_inum}.md ~/all/manuscript/polap-v0.4/
	cp -p ${_table}-data.md ~/all/manuscript/polap-v0.4/

	csvtk -t csv2md -a right ${_table_file} -o ${_table}.md

	echo "Check ${_table}.md"
	cat "${_table}.md"

	Rscript src/polap-data-v2.R
}

# previously table3

# table2.md
# table2.tsv
table2_genus_species() {
	local _arg_inum="${1:-2}"
	local _table="table2"
	local _table_file="${_table}.tsv"
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
	local _memory_gb_getorganelle
	local _total_hours_getorganelle
	local _memory_gb_flye1
	local _total_hours_flye1
	local _I
	local _P
	local _N
	local j

	printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
		"Species" \
		"N" \
		"C" \
		"L_size" \
		"L_cov" \
		"S_size" \
		"S_cov" \
		"P" \
		"R" \
		"Rate" \
		"Size" \
		"Alpha" \
		"Length2" \
		"N1" \
		"Mode" \
		"SD" \
		"N2" \
		"M" \
		"M_g" \
		"M_f" \
		"T" \
		>"${_table_file}"

	# Extract and sort keys
	sorted_keys=($(for key in "${!_table2[@]}"; do echo "$key"; done | sort))

	# Iterate over sorted keys and check if value is "T"
	for key in "${sorted_keys[@]}"; do
		if [[ "${_table2[$key]}" == "F" ]]; then
			continue
		fi
		local _v1=${_folder[$key]}
		local _extracted_inum=${_inum[$key]}
		if [[ "${_extracted_inum}" != "${_arg_inum}" ]]; then
			continue
		fi

		# echo "${_v1}"

		# done

		_logfile=${_v1}/polap.log
		_readmefile=${_v1}/README
		_species="${_v1/_/ }"

		_l_sra=$(awk '/long-read/ {match($0, /([A-Z]RR[0-9]+)/, arr); print arr[0]}' "$_logfile" | sort -u | grep -v '^$')
		if [[ -s "${_v1}/long_total_length.txt" ]]; then
			_l_sra_size=$(<"${_v1}/long_total_length.txt")
		else
			_l_sra_size=$(<"${_v1}/l.fq.txt")
		fi
		_l_sra_size_gb=$(_polap_utility_convert_bp "${_l_sra_size}")

		_s_sra=$(awk '/short-read1/ {match($0, /([A-Z]RR[0-9]+)/, arr); print arr[0]}' "$_logfile" | sort -u | grep -v '^$')
		if [[ -s "${_v1}/short_total_length.txt" ]]; then
			_s_sra_size=$(<"${_v1}/short_total_length.txt")
		else
			_s_sra_size1=$(<"${_v1}/s1.fq.txt")
			_s_sra_size2=$(<"${_v1}/s2.fq.txt")
			_s_sra_size=$((_s_sra_size1 + _s_sra_size2))
		fi
		_s_sra_size_gb=$(_polap_utility_convert_bp "${_s_sra_size}")
		_known_mtdna=$(grep 'NCBI accession:' ${_logfile} | cut -d: -f4 | tail -n 1)

		# if [[ -v _host["${_v1}"] ]]; then
		# 	echo species:"${_v1}"
		# else
		# 	echo "No such host for $_v1"
		# fi

		local _ptdna_ptgaul=${_v1}/ptdna-ptgaul.fa
		local _ptdna_reference=${_v1}/ptdna-reference.fa
		# seq_length1=$(grep -v "^>" "$_ptdna_ptgaul" | tr -d '\n' | wc -c)

		# read -r _memory_gb_flye1 _total_hours_flye1 < <(parse_timing "${_v1}" "infer12-4")

		local j=1
		local _i="${_extracted_inum}"
		# if [[ "${_i}" != "0" ]]; then
		# 	continue
		# fi
		local _v1_inum="${_v1}/${_i}"
		if [[ ! -d "${_v1_inum}" ]]; then
			continue
		fi
		local _disassemble_index="infer-${j}"
		local _summary1_ordered_txt=${_v1_inum}/disassemble/${_disassemble_index}/1/summary1-ordered.txt
		local fasta_file=${_v1_inum}/disassemble/${_disassemble_index}/pt.subsample-polishing.1.fa

		local target_index="${_v1}-${_i}"
		local extracted_memory="${_memory["$target_index"]}"

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
		local _n1=$(grep "^#n:" "$_summary1_ordered_txt" | awk 'NR==1 {print $2}')

		local _params_txt=${_v1_inum}/disassemble/${_disassemble_index}/params.txt
		_ipn=$(parse_params "${_params_txt}")
		read -r _I _P _N _R <<<"$_ipn"

		local _summary2_ordered_txt=${_v1_inum}/disassemble/${_disassemble_index}/2/summary1-ordered.txt
		local _n2=$(grep "^#n:" "$_summary2_ordered_txt" | awk 'NR==1 {print $2}')
		local output=$(awk -F'\t' 'NR==2 {print $1, $2, $4, $11}' "${_summary2_ordered_txt}")
		local _summary2_index
		local _summary2_size
		local _summary2_size_gb
		local _summary2_rate_rounded
		local _summary2_alpha
		local _summary2_alpha_formatted
		read -r _summary2_index _summary2_size _summary2_rate _summary2_alpha <<<"$output"
		_summary2_size_gb=$(_polap_utility_convert_bp "${_summary2_size}")
		local _summary2_rate_decimal=$(printf "%.10f" "$_summary2_rate")
		_summary2_rate_rounded=$(echo "scale=4; $_summary2_rate_decimal / 1" | bc)
		_summary2_alpha_formatted=$(echo "scale=2; $_summary2_alpha / 1" | bc | awk '{printf "%.2f\n", $1}')
		# _summary2_alpha_formatted=$(echo "$_summary2_alpha_formatted" | awk '{printf "%.10g\n", $1}')

		read -r _memory_gb _total_hours < <(_polap_lib_timing-parse-timing "${_v1_inum}/timing-infer-${j}-subsample-polish.txt")
		read -r _memory_gb_polishing _total_hours_polishing < <(_polap_lib_timing-parse-timing "${_v1}/timing/timing-prepare-polishing.txt")
		read -r _memory_gb_subsampling_polishing _total_hours_subsampling_polishing < <(_polap_lib_timing-parse-timing "${_v1_inum}/timing-infer-${j}-subsample-polish.txt")
		read -r _memory_gb_getorganelle _total_hours_getorganelle < <(_polap_lib_timing-parse-timing "${_v1}/timing/timing-getorganelle.txt")

		local _blast_pident="NA"
		if [[ -s "${_v1_inum}/blast/${j}/pident.txt" ]]; then
			_blast_pident="$(<${_v1_inum}/blast/${j}/pident.txt)"
		else
			_blast_pident=0
		fi

		local _mauve_lcb_coverage="NA"
		if [[ -s "${_v1_inum}/mauve/${j}/log.txt" ]]; then
			_mauve_lcb_coverage=$(awk '{print $2}' "${_v1_inum}/mauve/${j}/log.txt")
		else
			_mauve_lcb_coverage=0
		fi

		local _pident
		if (($(echo "$_blast_pident < $_mauve_lcb_coverage" | bc -l))); then
			_pident="${_mauve_lcb_coverage}"
		else
			_pident="${_blast_pident}"
		fi

		local _short_coverage=$(get_short_read_coverage "${_v1_inum}/sx.txt")
		local _long_coverage=$(get_long_read_coverage "${_v1_inum}/lx.txt")
		local _target_coverage=$(get_target_read_coverage "${_v1_inum}/lx.txt")

		# read -r _memory_gb _total_hours < <(parse_timing "${_v1}" "infer12-${j}")
		# read -r _memory_gb_polishing _total_hours_polishing < <(parse_timing "${_v1}" "infer3only-${j}")

		printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
			"_${_species}_" \
			"${_i}" \
			"${_target_coverage}" \
			"${_l_sra_size_gb}" \
			"${_long_coverage}" \
			"${_s_sra_size_gb}" \
			"${_short_coverage}" \
			"${_P}" \
			"${_R}" \
			"${_summary2_rate_rounded}" \
			"${_summary2_size_gb}" \
			"${_summary2_alpha_formatted}" \
			"${seq_length2}" \
			"${_n1}" \
			"${_mode}" \
			"${_sd}" \
			"${_n2}" \
			"${extracted_memory}" \
			"${_memory_gb_getorganelle}" \
			"${_memory_gb}" \
			"${_total_hours}" \
			>>"${_table_file}"
	done

	# table2 - 1
	# csvtk -t cut -f Species,L_size,L_cov,S_size,S_cov,ptDNA,Length1 \
	# 	${_table_file} | csvtk -t csv2md -a right -o ${_table}-1.md
	# cp -p ${_table}-1.md ~/all/manuscript/polap-v0.4/

	# csvtk -t cut -f Species,C,P,Rate,Size,Alpha,Pident,N1,Mode,SD,N2,M,M_g,M_p,M_c,M_f,T \
	# 	${_table_file} | csvtk -t csv2md -a right -o ${_table}-2.md
	# cp -p ${_table}-2.md ~/all/manuscript/polap-v0.4/
	#
	#

	csvtk -t cut -f Species,C,N,P,R,Rate,Alpha,Length2,N1,Mode,SD,M,M_g,M_f,T \
		${_table_file} |
		csvtk -t rename -f 1-15 -n Species,C,N,P,R,Rate,Alpha,L2,N1,Mode,SD,M,Mg,Mf,T |
		csvtk -t csv2md -a right -o ${_table}-analysis-${_arg_inum}.md

	# csvtk -t csv2md -a right ${_table_file} -o ${_table}-${_arg_inum}.md
	cp -p ${_table}-analysis-${_arg_inum}.md ~/all/manuscript/polap-v0.4/
	echo "Check ${_table}-analysis-${_arg_inum}.md"
	cat "${_table}-analysis-${_arg_inum}.md"
}

# stage 1 table
#
# suptable1.md
# 2 off 1
# 2 off 2
# 2 off 3-infer
# 2 on 1
# 2 on 2
# 2 on 3-infer
# 0 on
suptable1_genus_species() {
	local _arg_inum="${1:-2}"
	local _main="${2:-off}"
	local _arg_stage="${3:-1}"
	if [[ "${_main}" == "default_value" ]]; then
		_main="off"
	fi
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

	# printf "# Supplementary Tables\n\n" \
	# 	>"${_suptable_file}"

	local _i=${_arg_inum}
	# for j in 1 2 3-infer; do
	j="${_arg_stage}"

	local _suptable_file="suptable1-${_arg_inum}-${_arg_stage}.md"

	rm -f "${_suptable_file}"

	# Extract and sort keys
	sorted_keys=($(for key in "${!_table1[@]}"; do echo "$key"; done | sort))

	# Iterate over sorted keys and check if value is "T"
	for key in "${sorted_keys[@]}"; do
		if [[ "${_table1[$key]}" == "F" ]]; then
			continue
		fi
		local _v1=${_folder[$key]}
		if [[ "${_main}" == "on" ]]; then
			if [[ "${_v1}" != "Eucalyptus_pauciflora" ]]; then
				continue
			fi
		else
			if [[ "${_v1}" == "Eucalyptus_pauciflora" ]]; then
				continue
			fi
		fi

		extracted_inum=${_inum[$key]}
		if [[ "${extracted_inum}" != "${_arg_inum}" ]]; then
			continue
		fi

		_logfile=${_v1}/polap.log
		_readmefile=${_v1}/README
		_species="${_v1/_/ }"
		_label_base="${_v1/_/-}"
		_label_base=$(echo "$_label_base" | awk '{print tolower($0)}')

		# local _i=0
		# for j in 1 2 3-infer; do
		local _v1_inum="${_v1}/${_i}"
		local _params_txt=${_v1_inum}/disassemble/infer-1/params.txt
		_ipn=$(parse_params "${_params_txt}")
		read -r _I _P _N _R <<<"$_ipn"
		_label=${_label_base}

		# stage 1
		# printf "\\\blandscape\n\n" \
		# 	>>"${_suptable_file}"

		case "$j" in
		1)
			printf "Table: Plastid genome assemblies with the increasing subsample size upto the maximum subsampling rate of ${_P}%% from Stage 1 of Polap's subsampling-based analysis for the dataset of _${_species}_. {#tbl:suptable1-${_label}}\n\n" \
				>>"${_suptable_file}"
			;;
		2)
			printf "Table: Plastid genome assemblies with a fixed subsample size and subsampling-rate from Stage 2 of Polap's subsampling-based analysis for the dataset of _${_species}_. {#tbl:suptable2-${_label}}\n\n" \
				>>"${_suptable_file}"
			;;
		3-infer)
			printf "Table: Plastid genome assemblies with the subsampling-based short-read polishing from Stage 3 of Polap's subsampling-based analysis for the dataset of _${_species}_. {#tbl:suptable3-${_label}}\n\n" \
				>>"${_suptable_file}"
			;;
		*)
			echo default
			;;
		esac

		cat "${_v1_inum}/disassemble/infer-1/$j/summary1.md" \
			>>"${_suptable_file}"
		# Table: _Juncus effusus_ Polap's disassemble data analysis. i=3, p=10, n=50 \label{table-juncus-effusus} {#tbl:table-juncus-effusus}

		printf "\n" \
			>>"${_suptable_file}"

		case "$j" in
		1)
			cat src/polap-data-v2-suptable1_footnote.tex \
				>>"${_suptable_file}"
			;;
		2)
			cat src/polap-data-v2-suptable2_footnote.tex \
				>>"${_suptable_file}"
			;;
		3-infer)
			cat src/polap-data-v2-suptable3_footnote.tex \
				>>"${_suptable_file}"
			;;
		*)
			echo default
			;;
		esac

		# printf "\\\elandscape\n\n\\\newpage\n\n" \
		# 	>>"${_suptable_file}"

		printf "\n\n\\\newpage\n\n" \
			>>"${_suptable_file}"
	done
	# done

	if [[ "${_main}" == "on" ]]; then
		echo "See main-${_suptable_file}"
		cp -p ${_suptable_file} main-${_suptable_file}
		cp -p ${_suptable_file} ~/all/manuscript/polap-v0.4/main-${_suptable_file}
	else
		echo "See ${_suptable_file}"
		cp -p ${_suptable_file} ~/all/manuscript/polap-v0.4/
	fi
}

# figures of Polap, ptGAUL, and GetOrganelle's assemblies
#
# supfigure1.md
# 2 no off <- inum:2, no png extraction, supp figures
# 2 no on
# 2 yes off
# 2 yes on
supfigure1_genus_species() {
	local _arg_inum="${1:-2}"
	local _draw_bandage="${2:-default_value}"
	local _main="${3:-off}"
	local _supfigure_file="supfigure1.md"

	# printf "# Supplementary Figures: traceplot and bandage graphs\n\n" \
	# 	>"${_supfigure_file}"
	rm -f "${_supfigure_file}"

	# Extract and sort keys
	sorted_keys=($(for key in "${!_table1[@]}"; do echo "$key"; done | sort))

	# Iterate over sorted keys and check if value is "T"
	for key in "${sorted_keys[@]}"; do
		if [[ "${_table1[$key]}" == "F" ]]; then
			continue
		fi
		local _v1=${_folder[$key]}
		if [[ "${_main}" == "on" ]]; then
			if [[ "${_v1}" != "Eucalyptus_pauciflora" ]]; then
				continue
			fi
		else
			if [[ "${_v1}" == "Eucalyptus_pauciflora" ]]; then
				continue
			fi
		fi
		_species="${_v1/_/ }"
		_label_base="${_v1/_/-}"
		_label_base=$(echo "$_label_base" | awk '{print tolower($0)}')

		# bandage graph
		extracted_inum=${_inum[$key]}
		if [[ "${extracted_inum}" != "${_arg_inum}" ]]; then
			continue
		fi

		output_dir=${_v1}
		echo inum: ${extracted_inum}
		echo output_dir: ${output_dir}

		local j=1
		local _v1_inum="${_v1}/${extracted_inum}"
		local _params_txt=${_v1_inum}/disassemble/infer-1/params.txt
		_ipn=$(parse_params "${_params_txt}")
		read -r _I _P _N _R <<<"$_ipn"
		_label=${_label_base}

		# stage 1
		# printf "![Selection of a plastid genome using the length distribution for _${_species}_. (A) Trace plot of plastid genome length and the subsample-size upto ${_P} %%. (B) Length distribution for candidate plastid genome of _${_species}_ ](figures/${_v1_inum}/disassemble/infer-1/1/summary1-ordered.pdf){#fig:supfigure1-${_label}}\n\n" \
		# 	>>"${_supfigure_file}"

		printf "\\\newpage\n\n" \
			>>"${_supfigure_file}"

		printf "\n\n" \
			>>"${_supfigure_file}"

		#########################################################
		# all figures
		local extracted_n="${_compare_n["$key"]}"
		local i
		local width=13
		local images=()
		local captions=()

		local k=0
		for ((i = 0; i < extracted_n; i++)); do

			local _gfa_infer="${output_dir}/${extracted_inum}/disassemble/infer-1/1/${i}/30-contigger/graph_final.gfa"
			local _png_infer="${output_dir}/${extracted_inum}/disassemble/infer-1/1/${i}/30-contigger/graph_final.png"
			if [[ -s "${_gfa_infer}" ]]; then
				echo "gfa file: ${_gfa_infer}"
				if [[ "${_draw_bandage}" == "yes" ]]; then
					${_polap_cmd} bandage png \
						${_gfa_infer} \
						${_png_infer}
				fi
				# printf "| ![polap %s](figures/%s){ width=%s%% } " "${i}" "${_png_infer}" "${width}" >>"${_supfigure_file}"
				images+=("figures/${_png_infer}")
				captions+=(${i})
				((k++))
			else
				echo "no such file: ${_gfa_infer}"
			fi
		done

		#!/bin/bash

		# Define an array of PNG files and their subcaptions
		# images=("fig1.png" "fig2.png" "fig3.png" "fig4.png" "fig5.png" "fig6.png" "fig7.png" "fig8.png")
		# captions=("Subcaption 1" "Subcaption 2" "Subcaption 3" "Subcaption 4" "Subcaption 5" "Subcaption 6" "Subcaption 7" "Subcaption 8")

		# Output Markdown file
		output="${_supfigure_file}"

		cat <<EOF >>"$output"

Table: Plastid genome assembly graphs generated from Stage 1 of the subsampling-based method for the dataset of ${_species}. The graphs were generated with Bandage software. Each number corresponds to the iteration index of Stage 1 ([@tbl:suptable1-${_label_base}]). {#tbl:supfigure3-${_label_base}} 

EOF

		# Start writing to the markdown file
		cat <<EOF >>"$output"

| Organelle genome assemblies |  |  |
|-----------------|-----------------|-----------------|
EOF

		# Generate the image table with subcaptions
		count=0
		row="| "
		caption_row="| "

		for ((i = 0; i < ${#images[@]}; i++)); do
			row+="![${captions[i]}](${images[i]}){width=100px} | "
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

		# Handle the last incomplete row (if any)
		remaining=$((3 - count % 3))
		if ((remaining < 3)); then
			# Fill empty image columns
			for ((i = 0; i < remaining; i++)); do
				row+=" | "
				caption_row+=" | "
			done
			echo "$row" >>"$output"
			echo "$caption_row" >>"$output"
			echo "|-----------------|-----------------|-----------------|" >>"$output"
		fi

		echo "" >>"${output}"

	done

	if [[ "${_main}" == "on" ]]; then
		echo "main-${_supfigure_file}"
		cp -p ${_supfigure_file} main-${_supfigure_file}
		cp -p ${_supfigure_file} ~/all/manuscript/polap-v0.4/main-${_supfigure_file}
	else
		echo "${_draw_bandage}"
		echo "${_supfigure_file}"
		cp -p ${_supfigure_file} ~/all/manuscript/polap-v0.4/
	fi
}

# figures of Polap, ptGAUL, and GetOrganelle's assemblies

# supfigure2.md
# 2 no
# 2 yes to extract bandage graph png
supfigure2_genus_species() {
	local _arg_inum="${1:-2}"
	local _draw_bandage="${2:-default_value}"
	local _supfigure_file="supfigure2.md"

	# printf "# Supplementary Figures: Polap, ptGAUL, and GetOrganelle\n\n" \
	# 	>"${_supfigure_file}"
	rm -f "${_supfigure_file}"

	# Start writing to the markdown file
	cat <<EOF >>"${_supfigure_file}"

| Species | Polap | ptGAUL | GetOrganelle |
|-----------------|-----------------|-----------------|-----------------|
EOF

	# Extract and sort keys
	sorted_keys=($(for key in "${!_table1[@]}"; do echo "$key"; done | sort))

	# Iterate over sorted keys and check if value is "T"
	for key in "${sorted_keys[@]}"; do
		if [[ "${_table1[$key]}" == "F" ]]; then
			continue
		fi
		local _v1=${_folder[$key]}
		_species="${_v1/_/ }"
		_label_base="${_v1/_/-}"
		_label_base=$(echo "$_label_base" | awk '{print tolower($0)}')

		# bandage graph
		extracted_inum=${_inum[$key]}
		if [[ "${extracted_inum}" != "${_arg_inum}" ]]; then
			continue
		fi

		output_dir=${_v1}
		echo inum: ${extracted_inum}
		echo output_dir: ${output_dir}
		echo species: ${_species}

		if [[ "${_draw_bandage}" == "yes" ]]; then
			${_polap_cmd} disassemble bandage \
				-i ${extracted_inum} \
				-o ${output_dir} \
				--disassemble-i infer-1
		fi
		local _gfa_polap="${output_dir}/${extracted_inum}/disassemble/infer-1/pt.1.gfa"
		local _png_polap="${output_dir}/${extracted_inum}/disassemble/infer-1/pt.1.png"

		local images=()
		local captions=()
		images+=("${_species}")
		captions+=(" ")
		images+=(${_png_polap})
		captions+=("Polap")

		local _gfa_ptgaul="${output_dir}/ptgaul/flye_cpONT/assembly_graph.gfa"
		local _png_ptgaul="${output_dir}/ptgaul.png"
		if [[ "${_draw_bandage}" == "yes" ]]; then
			${_polap_cmd} bandage png \
				${_gfa_ptgaul} \
				${_png_ptgaul}
		fi
		images+=(${_png_ptgaul})
		captions+=("ptGAUL")

		local _gfa_getorganelle=$(find "${output_dir}/getorganelle" -type f -name 'embplant_pt*.gfa' | head -n 1)
		local _png_getorganelle="${output_dir}/getorganelle.png"
		if [[ "${_draw_bandage}" == "yes" ]]; then
			${_polap_cmd} bandage png \
				${_gfa_getorganelle} \
				${_png_getorganelle}
		fi
		images+=(${_png_getorganelle})
		captions+=("GetOrganelle")

		# echo "Warning: no such file: ${_gfa}"
		# echo "FIX QT5 problem:"
		# problem solved
		# echo "export QT_QPA_PLATFORM=offscreen"
		# Output Markdown file
		#
		output="${_supfigure_file}"

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
			if ((count % 4 == 0)); then
				echo "$row" >>"$output"
				echo "$caption_row" >>"$output"
				echo "|-----------------|-----------------|-----------------|-----------------|" >>"$output"
				# echo "" >>"$output"
				row="| "
				caption_row="| "
			fi
		done

	done

	cp -p ${_supfigure_file} ~/all/manuscript/polap-v0.4/
	echo ${_supfigure_file}
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
		echo "0. Batch"
		echo "1. Send data to the remote: ${_remote_name}"
		echo "2. Make a starting folder (3)"
		echo "3. Get ptDNA from NCBI"
		echo "4. GetOrganelle assembly with the short-read data"
		echo "5. ptGAUL assembly with the reference (6)"
		echo "6. Prepare the short-read polishing"
		echo "7. Extract the ptGAUL assembly (8)"
		echo "8. Copy the ptGAUL assembly"
		echo "9. Infer with POLAP"
		echo "10. Check with POLAP (needs step 5,7)"
		echo "11. (only for 5 datasets) Compare with POLAP (needs step 5,7)"
		echo "L. List the folder: ${_species_folder}"
		echo "S. Estimate the coverage ${_species_folder}"
		echo "A. Archive the ${_species_folder}"
		echo "G. Get the result from the remote"
		echo "C. Clean up the ${_species_folder}"
		echo "M. Back to Main Menu"
		echo "--------------------------------------------------"
		echo current host: "$(hostname)"
		echo rement host name: "${_host["$_species_name"]}"
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
						batch_genus_species "${keys_array[$species_key]}"
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
						getorganelle_genus_species \
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
					l)
						ls -l \
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
					s)
						downsample_genus_species \
							"${keys_array[$species_key]}" \
							5 --dry
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
	# navigate_menu "main"
	navigate_menu "bioproject"
}

# Main case statement
case "$subcmd1" in
'copy-figure')
	${subcmd1}_genus_species
	;;
'refs')
	get-ptdna-from-ncbi_genus_species "${_arg2}"
	;;
'send-data' | \
	'send' | \
	'recover' | \
	'mkdir' | \
	'sra' | \
	'get-ptdna-from-ncbi' | \
	'copy-ptdna-of-ncbi-as-reference' | \
	'getorganelle' | \
	'ptgaul' | \
	'msbwt' | \  | \
	'extract-ptdna-of-ptgaul' | \
	'extract-ptdna-of-ptgaul2' | \
	'copy-ptdna-of-ptgaul' | \
	'archive' | \
	'get' | \
	'report' | \
	'polishing' | \
	'wga' | \
	'polish' | \
	'clean' | \
	'restart' | \
	'clean-infer' | \
	'write-config' | \
	'table1' | \
	'table2' | \
	'single-argument')
	${subcmd1}_genus_species "${_arg2}"
	;;
'batch' | \
	'supfigure2' | \
	'compare' | \
	'use-downsample' | \
	'best' | \
	'mauve' | \
	bandage* | \
	'coverage' | \
	'delete-disassemble' | \
	'two-arguments')
	${subcmd1}_genus_species "${_arg2}" "${_arg3}"
	;;
'infer' | \
	'infer2' | \
	'suptable1' | \
	'supfigure1' | \
	'check')
	${subcmd1}_genus_species "${_arg2}" "${_arg3}" "${_arg4}"
	# ${subcmd1}_genus_species "${_arg2}" --disassemble-simple-polishing
	;;
'downsample')
	${subcmd1}_genus_species "${_arg2}" "${_arg3}" "${_arg4}"
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
"install-getorganelle")
	conda create --name getorganelle bioconda::getorganelle
	conda activate getorganelle
	get_organelle_config.py --add embplant_pt,embplant_mt
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
	for i in "${Sall[@]}"; do
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
