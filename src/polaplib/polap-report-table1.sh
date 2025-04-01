#!/bin/bash
S=('Spirodela polyrhiza' 'Taraxacum mongolicum' 'Trifolium pratense' 'Salix dunnii' 'Anthoceros agrestis' 'Anthoceros angustus' 'Brassica rapa' 'Vigna radiata' 'Macadamia tetraphylla' 'Punica granatum' 'Lolium perenne')
S=('Spirodela polyrhiza'
	'Taraxacum mongolicum'
	'Trifolium pratense'
	'Salix dunnii'
	'Anthoceros agrestis'
	'Anthoceros angustus'
	'Brassica rapa'
	'Vigna radiata'
	'Macadamia tetraphylla'
	'Punica granatum'
	'Lolium perenne')
# S=('Spirodela polyrhiza' 'Taraxacum mongolicum' 'Trifolium pratense' 'Anthoceros agrestis' 'Anthoceros angustus' 'Brassica rapa' 'Vigna radiata' 'Macadamia tetraphylla' 'Punica granatum' 'Lolium perenne')

# Define the header
printf "Species\tClass\tKnown_mtDNA\tLong_SRA\tShort_SRA\tLength\tIdentity\tTime(hrs)\tMemory(GB)\n" >table1.tsv

for i in "${S[@]}"; do
	V1="${i// /_}"
	LOGFILE=./$V1/o2/polap.log
	READMEFILE=./$V1/README
	_taxonomy=$(grep taxonomy: ${LOGFILE} | tail -1 | cut -d: -f4 | head -n 1)
	_class=$(grep taxonomy: ${LOGFILE} | tail -1 | cut -d: -f4 | cut -f10 | head -n 1)
	_order=$(grep taxonomy: ${LOGFILE} | tail -1 | cut -d: -f4 | cut -f15 | head -n 1)
	_species=$(grep taxonomy: ${LOGFILE} | tail -1 | cut -d: -f4 | cut -f24 | head -n 1)
	_l_sra=$(grep 'Long SRA:' ${LOGFILE} | cut -d: -f4 | head -n 1)
	_s_sra=$(grep 'Short SRA:' ${LOGFILE} | cut -d: -f4 | head -n 1)
	_length=$(grep 'mtDNA total length:' ${LOGFILE} | cut -d: -f4 | tail -n 1)
	_identity=$(grep 'mauve_lcb_length_coverage:' ${LOGFILE} | cut -d: -f4 | tail -n 1)
	_known_mtdna=$(grep 'NCBI accession:' ${LOGFILE} | cut -d: -f4 | tail -n 1)
	_species=$(grep 'species:' ${LOGFILE} | cut -d: -f4 | head -n 1)

	#   Viridiplanta
	# Elapsed
	#         Elapsed (wall clock) time (h:mm:ss or m:ss): 43:33.67
	#         Maximum resident set size (kbytes): 38957276
	_timing_file=timing/run-${V1}.sh
	_time_wga=$(grep 'Elapsed' "${_timing_file}" | head -1)
	_memory_wga=$(grep 'Maximum resident set size' "${_timing_file}" | head -1)

	# Extract the number in kilobytes
	memory_kbytes=$(echo "$_memory_wga" | grep -oE "[0-9]+")

	# Convert kilobytes to gigabytes
	memory_gb=$(echo "scale=2; $memory_kbytes / 1048576" | bc)

	# Extract the time portion using grep with regex
	time_only=$(echo "$_time_wga" | grep -oE "[0-9]+(:[0-9]{2}){1,2}")

	# Split the time into components based on the presence of hours
	IFS=":" read -r part1 part2 part3 <<<"$time_only"

	# Determine whether part1 is hours or minutes based on the presence of a third component
	if [ -z "$part3" ]; then
		hours=0
		minutes=$part1
		seconds=$part2
	else
		hours=$part1
		minutes=$part2
		seconds=$part3
	fi

	# Convert the time to total hours as an integer
	total_hours=$(echo "($hours + $minutes / 60 + $seconds / 3600)/1" | bc)

	# Display the result

	# Assuming all variables are already defined

	echo "----- Summary Report -----"
	printf "%-20s: %s\n" "Species" "${_species}"
	printf "%-20s: %s\n" "Taxonomy" "${_taxonomy}"
	printf "%-20s: %s\n" "Class" "${_class}"
	printf "%-20s: %s\n" "Order" "${_order}"
	printf "%-20s: %s\n" "Known mtDNA" "${_known_mtdna}"
	printf "%-20s: %s\n" "Long SRA" "${_l_sra}"
	printf "%-20s: %s\n" "Short SRA" "${_s_sra}"
	printf "%-20s: %s\n" "Length" "${_length}"
	if [[ -n "${_identity}" ]]; then
		printf "%-20s: %s\n" "Identity" "${_identity}"
	else
		printf "%-20s: %s\n" "Identity" N/A
	fi
	printf "%-20s: %s\n" "Total Time" "${total_hours} hrs"
	printf "%-20s: %s\n" "Memory Usage" "${memory_gb} GB"
	echo "-------------------------"

	# Print each species data in TSV format
	printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
		"$_species" "$_class" "$_known_mtdna" "$_l_sra" "$_s_sra" \
		"$_length" "$_identity" "$total_hours" "$memory_gb" >>table1.tsv

done
