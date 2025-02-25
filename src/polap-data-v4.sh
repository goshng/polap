#!/usr/bin/bash

# S=('Juncus_effusus' 'Juncus_inflexus' 'Juncus_roemerianus' 'Juncus_validus' 'Eucalyptus_pauciflora' 'Arctostaphylos_glauca' 'Lepidium_sativum' 'Chaetoceros_muellerii' 'Potentilla_micrantha' 'Durio_zibethinus' 'Beta_vulgaris' 'Eleocharis_dulcis' 'Leucanthemum_vulgare' 'Oryza_glaberrima' 'Cenchrus_americanus' 'Digitaria_exilis' 'Podococcus_acaulis' 'Raphia_textilis' 'Phytelephas_aequatorialis' 'Picea_glauca')
S=('Juncus_effusus' 'Juncus_inflexus' 'Juncus_roemerianus' 'Juncus_validus' 'Eucalyptus_pauciflora')
_local_host="thorne"

# Input parameter
subcmd1="${1:-0}"
subarg1="${2:-0}"

_polap_subcmd=(
	'test'
	'send-data-to'
	'mkdir'
	'get-ptdna-from-ncbi'
	'copy-ptdna-of-ncbi-as-reference'
	'ptgaul'
	'msbwt'
	'extract-ptdna-of-ptgaul'
	'copy-ptdna-of-ptgaul'
	'compare'
	'archive'
	'get'
	'report'
	'table1'
	'table2'
	'suptable1'
	'supfigure1'
	'infer'
	'infer2only'
	'polish'
	'check'
	'wga'
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
_polap_version="0.5.1.1"
_media_dir="/media/h1/run/ptgaul20"

help_message=$(
	cat <<HEREDOC
# Polap data analysis for the version 0.5 or directional (directional method)
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
3. get-ptdna-from-ncbi <species_folder>
4. copy-ptdna-of-ncbi-as-reference <species_folder>: copy the NCBI's ptDNA to a reference
5. ptgaul <species_folder>: ptGAUL analysis on the data
6. msbwt <species_folder>: prepare short-read polishing
7. extract-ptdna-of-ptgaul <species_folder>: extract ptDNA from ptGAUL's result
8. copy-ptdna-of-ptgaul <species_folder>: ptGAUL's result
9. compare <species_folder>: compare ptDNAs with the ptGAUL ptDNA
10. archive <species_folder>: report the results
11. get <species_folder>: fetch the archive from the remote
12. report [species_folder]: report the results
13. table1: table1
14. table2: table2
15. suptable1: Supplementary tables
16. supfigure1: Supplementary figures
17. infer <species_folder>: assemble the ptDNA without a reference ptDNA
18. infer2only <species_folder>: assemble stage2 only the ptDNA without a reference ptDNA
19. polish <species_folder>: polish the ptDNA using the short-read data
20. check <species_folder>: compare the polished infered ptDNA and the known one
21. wga <species_folder>: whole-genome assembly
HEREDOC
)

declare -A _host
_host['Juncus_effusus']="kishino"
_host['Eucalyptus_pauciflora']="lab01"
_host['Juncus_inflexus']="vincent"
_host['Juncus_roemerianus']="siepel"
_host['Juncus_validus']="marybeth"
_host['Eucalyptus_pauciflora']="user1-X99-PR8-H"

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

declare -A _long
_long['Juncus_effusus']="SRR14298760"
_long['Juncus_inflexus']="SRR14298751"
_long['Juncus_roemerianus']="SRR21976090"
_long['Juncus_validus']="SRR21976089"
_long['Eucalyptus_pauciflora']="SRR7153095"
declare -A _short
_short['Juncus_effusus']="SRR14298746"
_short['Juncus_inflexus']="SRR14298745"
_short['Juncus_roemerianus']="SRR21976092"
_short['Juncus_validus']="SRR21976091"
_short['Eucalyptus_pauciflora']="SRR7161123"
declare -A _ptgaul_genomesize
_ptgaul_genomesize['Juncus_effusus']="180000"
_ptgaul_genomesize['Juncus_inflexus']="180000"
_ptgaul_genomesize['Juncus_roemerianus']="200000"
_ptgaul_genomesize['Juncus_validus']="160000"
_ptgaul_genomesize['Eucalyptus_pauciflora']="160000"
declare -A _compare_n
_compare_n['Juncus_effusus']="50"
_compare_n['Juncus_inflexus']="50"
_compare_n['Juncus_roemerianus']="50"
_compare_n['Juncus_validus']="50"
_compare_n['Eucalyptus_pauciflora']="50"
declare -A _compare_p
_compare_p['Juncus_effusus']="1,5,10"
_compare_p['Juncus_inflexus']="1,5,10"
_compare_p['Juncus_roemerianus']="5,10,20"
_compare_p['Juncus_validus']="5,10,20"
_compare_p['Eucalyptus_pauciflora']="5,10,20"
declare -A _compare_r
_compare_n['Juncus_effusus']="5"
_compare_n['Juncus_inflexus']="5"
_compare_n['Juncus_roemerianus']="10"
_compare_n['Juncus_validus']="5"
_compare_n['Eucalyptus_pauciflora']="5"

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
	fi
	if [[ -s "${short_sra}_1.fastq" ]] &&
		[[ -s "${short_sra}_2.fastq" ]]; then
		echo short: $short_sra
	else
		echo short: no such fastq file
		echo run $0 send-data-to at ${_local_host}
	fi
}

send-data-to_genus_species() {
	local output_dir="$1"
	local species_name="$(echo $1 | sed 's/_/ /')"
	local long_sra="${_long["$1"]}"
	local short_sra="${_short["$1"]}"

	if [[ "${_local_host}" == "$(hostname)" ]]; then
		if [[ -s "${long_sra}.fastq" ]]; then
			scp "${_media_dir}/${long_sra}.fastq" ${_host["${output_dir}"]}:$PWD/
		fi
		if [[ -s "${short_sra}_1.fastq" ]]; then
			scp "${_media_dir}/${short_sra}_1.fastq" ${_host["${output_dir}"]}:$PWD/
		fi
		if [[ -s "${short_sra}_2.fastq" ]]; then
			scp "${_media_dir}/${short_sra}_2.fastq" ${_host["${output_dir}"]}:$PWD/
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

	command time -v bash src/ptGAUL1.sh \
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

compare_genus_species() {
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
			echo "($i) n=$n, p=$p"

			command time -v ${_polap_cmd} disassemble \
				-o ${output_dir} \
				-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
				--disassemble-c ${output_dir}/ptdna-ptgaul.fa \
				--disassemble-i $i \
				--disassemble-n $n \
				--disassemble-p $p \
				--disassemble-alpha 1.0 \
				2>${output_dir}/timing-${i}.txt
		done
	done

}

infer_genus_species() {
	local output_dir="$1"
	local species_name="$(echo $1 | sed 's/_/ /')"
	local long_sra="${_long["$1"]}"
	local short_sra="${_short["$1"]}"
	# copy_data

	local i=3
	local n
	local p
	IFS=',' read -r -a extracted_array_n <<<"${_compare_n["$1"]}"
	IFS=',' read -r -a extracted_array_p <<<"${_compare_p["$1"]}"
	local extracted_r="${_compare_r["$1"]}"
	for n in "${extracted_array_n[@]}"; do
		for p in "${extracted_array_p[@]}"; do
			i=$((i + 1))
			echo "($i) n=$n, p=$p"

			command time -v ${_polap_cmd} disassemble \
				-o ${output_dir} \
				-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
				--disassemble-i $i \
				--disassemble-n $n \
				--disassemble-p $p \
				--disassemble-r ${extracted_r} \
				--disassemble-alpha 1.0 \
				2>${output_dir}/timing-${i}.txt
		done
	done
}

infer2only_genus_species() {
	local output_dir="$1"
	local species_name="$(echo $1 | sed 's/_/ /')"
	local long_sra="${_long["$1"]}"
	local short_sra="${_short["$1"]}"
	# copy_data

	local i=3
	local n
	local p
	IFS=',' read -r -a extracted_array_n <<<"${_compare_n["$1"]}"
	IFS=',' read -r -a extracted_array_p <<<"${_compare_p["$1"]}"
	local extracted_r="${_compare_r["$1"]}"
	for n in "${extracted_array_n[@]}"; do
		for p in "${extracted_array_p[@]}"; do
			i=$((i + 1))
			echo "($i) n=$n, p=$p"

			command time -v ${_polap_cmd} disassemble stage2 \
				-o ${output_dir} \
				-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
				--disassemble-i $i \
				--disassemble-n $n \
				--disassemble-p $p \
				--disassemble-r ${extracted_r} \
				--disassemble-alpha 1.0 \
				2>${output_dir}/timing-${i}.txt
		done
	done
}

polish_genus_species() {
	local output_dir="$1"
	local species_name="$(echo $1 | sed 's/_/ /')"
	local long_sra="${_long["$1"]}"
	local short_sra="${_short["$1"]}"
	# copy_data

	local i=3
	local n
	local p
	IFS=',' read -r -a extracted_array_n <<<"${_compare_n["$1"]}"
	IFS=',' read -r -a extracted_array_p <<<"${_compare_p["$1"]}"
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

check_genus_species() {
	local output_dir="$1"
	local species_name="$(echo $1 | sed 's/_/ /')"
	local long_sra="${_long["$1"]}"
	local short_sra="${_short["$1"]}"
	# copy_data

	local i=3
	local n
	local p
	IFS=',' read -r -a extracted_array_n <<<"${_compare_n["$1"]}"
	IFS=',' read -r -a extracted_array_p <<<"${_compare_p["$1"]}"
	for n in "${extracted_array_n[@]}"; do
		for p in "${extracted_array_p[@]}"; do
			i=$((i + 1))
			echo "($i) n=$n, p=$p"

			command time -v ${_polap_cmd} disassemble check \
				-o ${output_dir} \
				--disassemble-c ${output_dir}/ptdna-ptgaul.fa \
				--disassemble-i $i \
				--polish \
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

	command time -v ${_polap_cmd} assemble1 \
		-o ${output_dir} \
		-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
		2>${output_dir}/timing-wga.txt
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

backup_genus_species() {
	local output_dir="$1"
	local species_name="$(echo $1 | sed 's/_/ /')"
	local long_sra="${_long["$1"]}"
	local short_sra="${_short["$1"]}"

	tar zcf ${output_dir}.tar.gz ${output_dir}
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

################################################################################
# Part of Juncus_effusus
#
# Function for each species
# kishino
scopy_Juncus_effusus() {
	local output_dir="$(echo $FUNCNAME | sed s/scopy_//)"
	local species_name="$(echo $FUNCNAME | sed 's/scopy_//' | sed 's/_/ /')"
	local long_sra="SRR14298760"
	local short_sra="SRR14298746"
	scopy_data_to_remote
}

ptgaul_Juncus_effusus() {
	local output_dir="$(echo $FUNCNAME | sed s/ptgaul_//)"
	local species_name="$(echo $FUNCNAME | sed 's/ptgaul_//' | sed 's/_/ /')"
	local long_sra="SRR14298760"
	local short_sra="SRR14298746"
	bash src/ptGAUL1.sh \
		-o ${output_dir}-ptgaul \
		-r ptdna-${output_dir}.fa \
		-g 180000 \
		-l ${long_sra}.fastq \
		-t 24
}

compare_Juncus_effusus() {
	local output_dir="$(echo $FUNCNAME | sed s/run_//)"
	local species_name="$(echo $FUNCNAME | sed 's/run_//' | sed 's/_/ /')"
	local long_sra="SRR14298760"
	local short_sra="SRR14298746"
	# copy_data

	local i=0
	local n
	local p
	for n in 30 50; do
		for p in 1 5; do
			i=$((i + 1))
			if [[ -d "${output_dir}/disassemble/${i}" ]]; then
				echo "exists: $i, $p, $n"
			else
				command time -v ${_polap_cmd} disassemble \
					-o ${output_dir} \
					-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
					--disassemble-compare-to-fasta ptdna-Juncus_effusus-ptgaul.fa \
					--stages-include 1 \
					--disassemble-best \
					--disassemble-polish \
					--no-contigger \
					--disassemble-i $i \
					--disassemble-p $p \
					--disassemble-n $n 2>{output_dir}/timing-${i}.txt
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

run_Juncus_effusus() {
	local output_dir="$(echo $FUNCNAME | sed s/run_//)"
	local species_name="$(echo $FUNCNAME | sed 's/run_//' | sed 's/_/ /')"
	local long_sra="SRR14298760"
	local short_sra="SRR14298746"
	# copy_data

	local i=0
	local n
	local p
	for n in 30 50; do
		for p in 1 5; do
			i=$((i + 1))
			if [[ -d "${output_dir}/disassemble/${i}" ]]; then
				echo "exists: $i, $p, $n"
			else
				command time -v ${_polap_cmd} disassemble \
					-o ${output_dir} \
					-l ${long_sra}.fastq -a ${short_sra}_1.fastq -b ${short_sra}_2.fastq \
					--disassemble-compare-to-fasta ptdna-Juncus_effusus-ptgaul.fa \
					--stages-include 1 \
					--disassemble-best \
					--disassemble-polish \
					--no-contigger \
					--disassemble-i $i \
					--disassemble-p $p \
					--disassemble-n $n 2>{output_dir}/timing-${i}.txt
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
				${_polap_cmd} disassemble \
					-o ${output_dir} \
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

################################################################################
# Part of Juncus_inflexus
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

################################################################################
# Part of Juncus_roemerianus
#
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

################################################################################
# Part of Juncus_validus
#
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

################################################################################
# Part of Eucalyptus_pauciflora
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
case "$subcmd1" in
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
	"Picea_glauca")
	run_${subcmd1}
	;;
'send-data-to' | \
	'mkdir' | \
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
	'infer' | \
	'infer2only' | \
	'polish' | \
	'check' | \
	'wga' | \
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
	echo "Species '$subcmd1' is not recognized."
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
