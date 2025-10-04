#!/bin/bash
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
# This is called by polap-data-v3.sh not polap-data-v2 or v4.
# FIXME: 
# v3 is polap-data-taxon, which we have not finished with that yet.
# TEST-SCC: not yet. I do not know what we should do with this.
################################################################################
# This is called by polap-data-v3.sh not polap-data-v2 or v4.

S=('Juncus_effusus' 'Juncus_inflexus' 'Juncus_roemerianus' 'Juncus_validus' 'Eucalyptus_pauciflora')
S=('Juncus_effusus' 'Juncus_roemerianus' 'Juncus_validus' 'Eucalyptus_pauciflora')
S=('Juncus_effusus' 'Eucalyptus_pauciflora')

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

b="ptgaul20"

# Define the header
printf "Species\tP\tN\tLong_SRA\tShort_SRA\tKnown_mtDNA\tAlignment Coverage\tTime(hrs)\tMemory(GB)\n" >table1.tsv

for i in "${S[@]}"; do
	V1="${i// /_}-a"
	LOGFILE=${b}/${V1}/polap.log
	READMEFILE=${b}/${V1}/README

	# _taxonomy=$(grep taxonomy: ${LOGFILE} | tail -1 | cut -d: -f4 | head -n 1)
	# _class=$(grep taxonomy: ${LOGFILE} | tail -1 | cut -d: -f4 | cut -f10 | head -n 1)
	# _order=$(grep taxonomy: ${LOGFILE} | tail -1 | cut -d: -f4 | cut -f15 | head -n 1)
	# _species=$(grep taxonomy: ${LOGFILE} | tail -1 | cut -d: -f4 | cut -f24 | head -n 1)
	# _l_sra=$(grep 'Long SRA:' ${LOGFILE} | cut -d: -f4 | head -n 1)
	# _s_sra=$(grep 'Short SRA:' ${LOGFILE} | cut -d: -f4 | head -n 1)
	# _length=$(grep 'mtDNA total length:' ${LOGFILE} | cut -d: -f4 | tail -n 1)
	# _identity=$(grep 'mauve_lcb_length_coverage:' ${LOGFILE} | cut -d: -f4 | tail -n 1)
	# _known_mtdna=$(grep 'NCBI accession:' ${LOGFILE} | cut -d: -f4 | tail -n 1)
	# _species=$(grep 'species:' ${LOGFILE} | cut -d: -f4 | head -n 1)

	#   Viridiplanta
	# Elapsed
	#         Elapsed (wall clock) time (h:mm:ss or m:ss): 43:33.67
	#         Maximum resident set size (kbytes): 38957276
	_species="${i/_/ }"
	# _l_sra=$(grep 'Long SRA:' ${LOGFILE} | cut -d: -f4 | head -n 1)
	# _s_sra=$(grep 'Short SRA:' ${LOGFILE} | cut -d: -f4 | head -n 1)
	_l_sra=$(awk '/long-read/ {match($0, /([A-Z]RR[0-9]+)/, arr); print arr[0]}' "$LOGFILE" | sort -u | grep -v '^$')
	_s_sra=$(awk '/short-read1/ {match($0, /([A-Z]RR[0-9]+)/, arr); print arr[0]}' "$LOGFILE" | sort -u | grep -v '^$')
	_known_mtdna=$(grep 'NCBI accession:' ${LOGFILE} | cut -d: -f4 | tail -n 1)

	# for j in {1..9}; do
	for j in 6; do

		_timing_file=${b}/${V1}/timing-${j}.txt
		params_txt=${b}/${V1}/disassemble/${j}/params.txt
		percent_identity=$(<"${b}/${V1}/disassemble/${j}/c/coverage.txt")
		# percent_identity2=$(printf "%.5f" "$percent_identity")
		# echo $percent_identity2

		ipn=$(parse_params "${params_txt}")
		read -r I P N <<<"$ipn"

		if [[ -s "${_timing_file}" ]]; then
			_time_wga=$(grep 'Elapsed' "${_timing_file}" | head -1)
			_memory_wga=$(grep 'Maximum resident set size' "${_timing_file}" | head -1)
			# Extract the number in kilobytes
			memory_kbytes=$(echo "$_memory_wga" | grep -oE "[0-9]+")
			# Convert kilobytes to gigabytes
			memory_gb=$(echo "scale=2; $memory_kbytes / 1048576" | bc)
			# Extract the time portion using grep with regex
			# time_only=$(echo "$_time_wga" | grep -oE "[0-9]+(:[0-9]{2}){1,2}")
			time_only=$(grep "Elapsed (wall clock) time (h:mm:ss or m:ss):" "$_timing_file" | awk -F': ' '{print $2}')

			total_hours=$(convert_to_hours "${time_only}")
		else
			memory_gb=0
			total_hours=0
		fi

		# Assuming all variables are already defined
		echo "----- Summary Report -----"
		printf "%-20s: %s\n" "Species" "${_species}"
		# printf "%-20s: %s\n" "Taxonomy" "${_taxonomy}"
		# printf "%-20s: %s\n" "Class" "${_class}"
		# printf "%-20s: %s\n" "Order" "${_order}"
		# printf "%-20s: %s\n" "Known mtDNA" "${_known_mtdna}"
		# printf "%-20s: %s\n" "Long SRA" "${_l_sra}"
		# printf "%-20s: %s\n" "Short SRA" "${_s_sra}"
		# printf "%-20s: %s\n" "Length" "${_length}"
		# if [[ -n "${_identity}" ]]; then
		# 	printf "%-20s: %s\n" "Identity" "${_identity}"
		# else
		# 	printf "%-20s: %s\n" "Identity" N/A
		# fi
		printf "%-20s: %s\n" "Total Time" "${total_hours} hrs"
		printf "%-20s: %s\n" "Memory Usage" "${memory_gb} GB"
		echo "-------------------------"

		# Print each species data in TSV format
		# printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
		# 	"$_species" "$_class" "$_known_mtdna" "$_l_sra" "$_s_sra" \
		# 	"$_length" "$_identity" "$total_hours" "$memory_gb" >>table1.tsv
		printf "%s\t%s\t%s\t%s\t%s\t%s\t%.7f\t%s\t%s\n" \
			"$_species" \
			"$P" \
			"$N" \
			"$_l_sra" "$_s_sra" \
			"$_known_mtdna" \
			"$percent_identity" \
			"$total_hours" "$memory_gb" >>table1.tsv
	done

done
