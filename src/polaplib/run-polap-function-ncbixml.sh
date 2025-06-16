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
# TODO: not related with polap but related with taxon module.
# TODO: rename: ncbixml -> tree or phylogeny
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

# Define the get_family function
get_family() {
	local taxid="$1" # Input TaxID
	local taxon="$2" # Input order or family
	local family=""  # Variable to store the family name

	if [[ -z "${taxon}" ]]; then
		taxon="order"
	fi

	# Traverse the taxonomy hierarchy
	while [ "$taxid" != "1" ]; do
		# Get the current TaxID details from the taxonomy file
		entry=$(grep -P "^$taxid\t" taxonomy_full.tsv)
		if [ -z "$entry" ]; then
			echo "TaxID $taxid not found in the taxonomy database." >&2
			return 1
		fi

		rank=$(echo "$entry" | cut -f3) # Extract rank
		name=$(echo "$entry" | cut -f4) # Extract scientific name

		# Check if the rank is "family"
		if [ "$rank" == "${taxon}" ]; then
			family="$name"
			break
		fi

		# Update TaxID to the Parent_TaxID
		taxid=$(echo "$entry" | cut -f2)
	done

	# Output the family name or indicate that it's not found
	if [ -n "$family" ]; then
		echo "$family"
		return 0
	else
		echo "Family not found for the given TaxID." >&2
		return 1
	fi
}

get_taxonomy_info() {
	local species_name="$1"
	local seq_id="$2"
	local taxid=""
	local order=""
	local family=""

	echo $seq_id
	# Find the Taxonomy ID for the species name
	taxid=$(grep -P "\t${species_name}$" taxonomy_full.tsv | cut -f1)

	if [ -z "$taxid" ]; then
		echo "Error: Species '$species_name' not found in the taxonomy database." >&2
		return 1
	fi

	# Traverse the taxonomy hierarchy
	while [ "$taxid" != "1" ]; do
		# Get the current TaxID details
		entry=$(grep -P "^$taxid\t" taxonomy_full.tsv)
		if [ -z "$entry" ]; then
			echo "Error: TaxID $taxid not found in the taxonomy database." >&2
			return 1
		fi

		rank=$(echo "$entry" | cut -f3) # Extract rank
		name=$(echo "$entry" | cut -f4) # Extract scientific name

		# Save the order and family if found
		if [ "$rank" == "order" ]; then
			order="$name"
		elif [ "$rank" == "family" ]; then
			family="$name"
		fi

		# Stop if both order and family are found
		if [ -n "$order" ] && [ -n "$family" ]; then
			break
		fi

		# Update TaxID to the Parent_TaxID
		taxid=$(echo "$entry" | cut -f2)
	done

	# Output the results
	# echo -e "Species Name: $species_name"
	# echo -e "Taxonomy ID: $taxid"
	# echo -e "Order: $order"
	# echo -e "Family: $family"
	echo -e "\t${taxid}\t${species_name}\t${order}\t${family}"
}

function _run_polap_ncbixml {
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

	# Print help message if requested
	help_message=$(
		cat <<HEREDOC

${_polap_var_ncbi}
${_polap_var_ncbi_sequence_tsv}
Example: $(basename $0) ${_arg_menu[0]}
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
	fi

	if [[ "${_arg_menu[1]}" == "1" ]]; then
		_polap_log0 Go to https://www.ncbi.nlm.nih.gov/sites/batchentrez to fetch sequences in batch mode.
		_polap_log0 You could copy the content of the ${_polap_var_ncbi_accessions}
		_polap_log0 and paste it to the query at https://www.ncbi.nlm.nih.gov/nuccore
		_polap_log0 Then, use Send to: button to select TinySeq XML to save sequences and others.
		_polap_log0 Create a flie: sequence.fasta.xml
		mkdir "${_polap_var_ncbi}"
	fi

	if [[ "${_arg_menu[1]}" == "2" ]]; then

		# Input and output file names
		if [[ "${_arg_menu[2]}" == "outfile" ]]; then
			local input_file="sequence.fasta.xml"
		else
			local input_file="${_arg_menu[2]}"
		fi

		output_file="${_polap_var_ncbi_sequence_tsv}"

		# Write the TSV header
		echo -e "accver\ttaxid\torgname\tdefline\tlength" >"$output_file"

		# Extract the fields and write to the TSV file
		xmllint --format "$input_file" |
			awk '
/<TSeq>/,/<\/TSeq>/ {
    if ($1 ~ /<TSeq_accver>/) accver=gensub(/.*<TSeq_accver>(.*)<\/TSeq_accver>.*/, "\\1", "g")
    if ($1 ~ /<TSeq_taxid>/) taxid=gensub(/.*<TSeq_taxid>(.*)<\/TSeq_taxid>.*/, "\\1", "g")
    if ($1 ~ /<TSeq_orgname>/) orgname=gensub(/.*<TSeq_orgname>(.*)<\/TSeq_orgname>.*/, "\\1", "g")
    if ($1 ~ /<TSeq_defline>/) defline=gensub(/.*<TSeq_defline>(.*)<\/TSeq_defline>.*/, "\\1", "g")
    if ($1 ~ /<TSeq_length>/) len=gensub(/.*<TSeq_length>(.*)<\/TSeq_length>.*/, "\\1", "g")
    if ($1 ~ /<\/TSeq>/) {
        printf("%s\t%s\t%s\t%s\t%s\n", accver, taxid, orgname, defline, len) >> "'$output_file'"
        accver=""; taxid=""; orgname=""; defline=""; len=""
    }
}'
		_polap_log0_head "${output_file}"

	fi

	if [[ "${_arg_menu[1]}" == "3-download-db" ]]; then
		# wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
		tar -xvzf taxdump.tar.gz
		# Extract relevant data from nodes.dmp and names.dmp
		awk -F "\t|\t" 'BEGIN {OFS="\t"} {print $1, $3, $5}' nodes.dmp >nodes.tsv
		awk -F "\t|\t" 'BEGIN {OFS="\t"} $7 == "scientific name" {print $1, $3}' names.dmp >names.tsv

		# Merge the two files on TaxID
		join -t $'\t' -1 1 -2 1 <(sort -k1,1 nodes.tsv) <(sort -k1,1 names.tsv) >taxonomy_full.tsv

		_polap_log0 "File: taxonomy_full.tsv"
	fi

	if [[ "${_arg_menu[1]}" == "4-get-taxon" ]]; then

		# Add Class-level and Order-level TaxIDs
		while IFS=$'\t' read -r acc taxid organism desc len; do
			# Fetch taxonomy details for the current TaxID
			local _family=$(get_family "$taxid" family)
			local _order=$(get_family "$taxid" order)
			echo -e "${acc}\t${taxid}\t${organism}\t${_order}\t${_family}"
			_order=""
			_family=""
		done < <(tail -n +2 "${_polap_var_ncbi_sequence_tsv}") >"${_polap_var_ncbi_sequence_taxon_tsv}"

		_polap_log0 "output: ${_polap_var_ncbi_sequence_taxon_tsv}"
	fi

	if [[ "${_arg_menu[1]}" == "5-sample" ]]; then
		if [[ "${_arg_menu[2]}" == "outfile" ]]; then
			# sample 2 or 3 species sequences for each family or order
			Rscript "${_POLAPLIB_DIR}/run-polap-r-taxonomy.R" \
				sample \
				-t ${_polap_var_ncbi_sequence_taxon_tsv} \
				-o ${_polap_var_ncbi_sampled_accession}
		else
			Rscript "${_POLAPLIB_DIR}/run-polap-r-taxonomy.R" \
				sample \
				-f \
				-t ${_polap_var_ncbi_sequence_taxon_tsv} \
				-o ${_polap_var_ncbi_sampled_accession}
		fi
		_polap_log0 file: ${_polap_var_ncbi_sampled_accession}
	fi

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_polap_geseq {
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

	# Print help message if requested
	help_message=$(
		cat <<HEREDOC
# Prepare files for IQ-TREE inference.
#
# Arguments:
#
# Inputs:
#
# Outputs:
#
# See:
Example: $(basename $0) ${_arg_menu[0]} -a <other-out-folder> --start-index 1 --species "Eucalyptus pauciflora"
Example: $(basename $0) ${_arg_menu[0]} step1
Example: $(basename $0) ${_arg_menu[0]} step2
Example: $(basename $0) ${_arg_menu[0]} step3
Example: $(basename $0) ${_arg_menu[0]} step4 --start-index 1
Example: $(basename $0) ${_arg_menu[0]} step5 --species "Eucalyptus pauciflora"
Example: $(basename $0) ${_arg_menu[0]} step6 -a <other-out-folder> --species "Eucalyptus pauciflora"
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
	fi

	_polap_log0 "start: ------------------------------ ${_arg_menu[1]}"
	if [[ "${_arg_menu[1]}" == "step0" ]] || [[ "${_arg_menu[1]}" == "infile" ]]; then
		_polap_log0 "step0: ------------------------------"

		# Define the folder containing the .zip file
		folder="${_polap_var_outdir}"

		# Get the .zip file name in the folder
		zip_file=$(find "$folder" -maxdepth 1 -type f -name "job-results*.zip" | head -n 1)

		if [[ -n "$zip_file" ]]; then
			# Extract the filename without the path
			zip_filename=$(basename "$zip_file")

			# Print the filename
			_polap_log0 "Zip file found: $zip_filename"

			# Extract the .zip file in the same folder
			unzip "$zip_file" -d "$folder"
			_polap_log0 "Contents extracted in $folder"
		else
			_polap_log0 "No .zip file found in $folder"
		fi
	fi

	if [[ "${_arg_menu[1]}" == "step1" ]] || [[ "${_arg_menu[1]}" == "infile" ]]; then
		_polap_log0 "step1: ------------------------------"

		mkdir -p "${_polap_var_geseq}"
		_polap_log0 "copy the gff3 file to folder: ${_polap_var_geseq}"

		# Define the target folder to copy the files
		local target_folder="${_polap_var_geseq}"

		# Iterate over each job-results folder
		for folder in "${_polap_var_outdir}"/job-results-*/; do
			# Check for multi-GFF3.gff3 file
			multi_file=$(find "$folder" -type f -name "*multi-GFF3.gff3" | head -n 1)

			if [[ -n "$multi_file" ]]; then
				# If a multi-GFF3.gff3 file exists, copy it
				_polap_log0 "Found multi-GFF3 file: $multi_file."
				_polap_log0 "  Copying it to $target_folder."
				_polap_log0_cmd cp "$multi_file" "${_polap_var_geseq_gff}"
			else
				# Otherwise, look for the general GFF3.gff3 file
				gff3_file=$(find "$folder" -type f -name "*GFF3.gff3" | head -n 1)
				if [[ -n "$gff3_file" ]]; then
					_polap_log0 "No multi-GFF3 file in $folder."
					_polap_log0 "  Copying GFF3 file: $gff3_file to $target_folder."
					_polap_log0_cmd cp "$gff3_file" "${_polap_var_geseq_gff}"
				else
					_polap_log0 "No GFF3 files found in $folder."
				fi
			fi
		done
		if [[ -s "${_polap_var_geseq_gff}" ]]; then
			sed -i 's/1\tI/0\tI/' "${_polap_var_geseq_gff}"
			_polap_log0 "correct the conda start position: ${_polap_var_geseq_gff}"
		fi
	fi

	if [[ "${_arg_menu[1]}" == "step2" ]] || [[ "${_arg_menu[1]}" == "infile" ]]; then
		_polap_log0 "step2: ------------------------------"
		# extract
		gffread -g ${_polap_var_mt1fa} -y ${_polap_var_geseq_faa} ${_polap_var_geseq_gff}
		gffread -g ${_polap_var_mt1fa} -x ${_polap_var_geseq_fna} ${_polap_var_geseq_gff}
		_polap_log0 "${_polap_var_geseq_faa}"
		_polap_log0 "${_polap_var_geseq_fna}"
	fi

	if [[ "${_arg_menu[1]}" == "step3" ]] || [[ "${_arg_menu[1]}" == "infile" ]]; then
		_polap_log0 "step3: ------------------------------"
		_polap_log0 "remove sequences with non-stanadrad amino acids: ${_polap_var_geseq_faa}"
		_polap_log0_pipe "seqkit grep -s -r -n -v -p '[^ACDEFGHIKLMNPQRSTVWY*]' \
      ${_polap_var_geseq_faa} >${_polap_var_geseq_filtered_faa}"
	fi

	if [[ "${_arg_menu[1]}" == "step4" ]] || [[ "${_arg_menu[1]}" == "infile" ]]; then
		_polap_log0 "step4: ------------------------------"
		# Usage: ./modify_fasta.sh input.fasta output.fasta MG37

		_polap_log0 "modify the sequence IDs to conform to NCBI's: ${_polap_var_geseq_formated_faa}"
		# Input arguments
		local INPUT_FASTA="${_polap_var_geseq_filtered_faa}"
		local OUTPUT_FASTA="${_polap_var_geseq_formated_faa}"
		PREFIX="AA"
		local NUMBER="${_arg_start_index}"
		if [[ "${_arg_start_index}" -eq 0 ]]; then
			local NUMBER=1
		fi
		# Generate AA000001 (6 digits padded with zeros)
		local ID1=$(printf "%s%06d" "$PREFIX" "$NUMBER")
		# Generate AA01 (2 digits padded with zeros)
		local ID2=$(printf "%s%02d" "$PREFIX" "$NUMBER")
		local PREFIX1=${ID1}
		local PREFIX2=${ID2}

		# Temporary file for processing
		TEMP_FILE=$(mktemp)

		# Process the FASTA file
		awk -v prefix1="$PREFIX1" -v prefix2="$PREFIX2" '
  BEGIN {
    counter = 1
  }
  /^>/ {
    # Extract the gene name from the header
    split($0, arr, "cds-blatx_")
    gene_id = arr[2]

    # Generate the new ID
    padded_counter = sprintf("%04d", counter)
    new_id = "lcl|" prefix1 ".1_prot_" prefix2 padded_counter ".1_" counter " [gene=" gene_id "]"
    
    # Print the new header
    print ">" new_id

    # Increment the counter
    counter++
  }
  /^[^>]/ {
    # Print the sequence line as is
    print
  }
' "$INPUT_FASTA" >"$TEMP_FILE"

		# Save the processed output
		mv "$TEMP_FILE" "$OUTPUT_FASTA"

		rm -f "${_polap_var_geseq}/${PREFIX}"*.faa
		cp "${OUTPUT_FASTA}" "${_polap_var_geseq}/${ID1}.faa"
		_polap_log0 "Modified FASTA file written to $OUTPUT_FASTA"
		_polap_log0 "Modified FASTA file written to ${_polap_var_geseq}/${ID1}.faa"
	fi

	if [[ "${_arg_menu[1]}" == "step5" ]] || [[ "${_arg_menu[1]}" == "infile" ]]; then
		_polap_log0 "step5: ------------------------------"
		local species_name="${_arg_species}"

		first_file=$(basename "$(ls ${_polap_var_geseq}/AA*.faa 2>/dev/null | head -n 1)" .faa)
		# if [ -n "$first_file" ]; then
		# 	local file_base="${first_file%.faa}"
		# else
		# 	echo "No file matches the pattern."
		# fi

		_polap_log0 $(get_taxonomy_info "$species_name" "$first_file")
	fi

	# copy the AA file to taxonomy analysis folder: orthofinder
	if [[ "${_arg_menu[1]}" == "step6" ]] || [[ "${_arg_menu[1]}" == "infile" ]]; then
		_polap_log0 "step6: ------------------------------"
		local species_name="${_arg_species}"
		local taxonomy_folder="${_arg_short_read1}"
		local taxonomy_orthofinder_folder="${taxonomy_folder}/orthofinder"

		_polap_log0_cmd mkdir -p "${taxonomy_folder}/orthofinder"

		first_file=$(basename "$(ls ${_polap_var_geseq}/AA*.faa 2>/dev/null | head -n 1)" .faa)
		# if [ -n "$first_file" ]; then
		# 	local file_base="${first_file%.faa}"
		# else
		# 	echo "No file matches the pattern."
		# fi

		_polap_log0_cmd cp "${_polap_var_geseq}/$first_file.faa" "${taxonomy_orthofinder_folder}"

		echo $(get_taxonomy_info "$species_name" "$first_file") >>"${taxonomy_folder}/custom.taxa.txt"
	fi

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_polap_orthofinder {
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

	# Print help message if requested
	help_message=$(
		cat <<HEREDOC
# Prepare files for IQ-TREE inference.
#
# Arguments:
#
# Inputs:
#
# Outputs:
#
# See:
Example: $(basename $0) ${_arg_menu[0]} -f sequence.txt -r 20 -x 20
Example: $(basename $0) ${_arg_menu[0]} step1 -f sequence.txt
Example: $(basename $0) ${_arg_menu[0]} step2 -r 20
Example: $(basename $0) ${_arg_menu[0]} step4
Example: $(basename $0) ${_arg_menu[0]} step5 -x 20
Example: $(basename $0) ${_arg_menu[0]} step6 -x 20
Example: $(basename $0) ${_arg_menu[0]} step3 <- copy your AA fasta files done in geseq
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
	fi

	# sed -i 's/1\tI/0\tI/' ccapricornis-20230702-142430_chr1_GFF3.gff3

	if [[ -d "${_polap_var_orthofinder}" ]]; then
		_polap_log0 "ERROR: no such folder: ${_polap_var_orthofinder}"
	fi

	if [[ "${_arg_menu[1]}" == "step1" ]] || [[ "${_arg_menu[1]}" == "infile" ]]; then
		# local SPECIES=species.of
		local SPECIES="${_polap_var_ncbi_sampled_accession}"
		local FILTERFASTA="sequence1514.txt"
		local FILTERFASTA="${_arg_final_assembly}" # -f sequence1514.txt

		mkdir -p "${_polap_var_orthofinder}"
		while read -r g; do
			seqkit grep -r -p "l/|$g" $FILTERFASTA -o "${_polap_var_orthofinder}/$g.fa"
		done < <(cut -d',' -f1 $SPECIES | tail -n +2)
	fi

	# filter by the number of amino acid sequences
	if [[ "${_arg_menu[1]}" == "step2" ]] || [[ "${_arg_menu[1]}" == "infile" ]]; then
		# Get folder path and threshold length from arguments
		FASTA_FOLDER="${_polap_var_orthofinder}"
		LENGTH_THRESHOLD="${_arg_pair_min}" # -r 20

		# Validate that the threshold is a positive integer
		if ! [[ "$LENGTH_THRESHOLD" =~ ^[0-9]+$ ]]; then
			echo "Error: Length threshold must be a positive integer."
			exit 1
		fi

		# Loop through each FASTA file in the folder
		for FILE in "$FASTA_FOLDER"/*.fa; do
			# Get the maximum sequence length in the file using seqkit stats
			MAX_LENGTH=$(seqkit stats -T "$FILE" | awk 'NR>1 {print $4}') # 5th column is the max length

			# Check if the max length is below the threshold
			if [ "$MAX_LENGTH" -lt "$LENGTH_THRESHOLD" ]; then
				rm "$FILE"
			fi
		done
	fi

	if [[ "${_arg_menu[1]}" == "step3" ]]; then
		_polap_log0 "Copy your AA fasta files to ${_polap_var_orthofinder}"
	fi

	if [[ "${_arg_menu[1]}" == "step4" ]] || [[ "${_arg_menu[1]}" == "infile" ]]; then
		orthofinder -f "${_polap_var_orthofinder}"
	fi

	if [[ "${_arg_menu[1]}" == "step5" ]] || [[ "${_arg_menu[1]}" == "infile" ]]; then
		local N="{_arg_bridge_min}" # -x 20
		local of="${_polap_var_orthofinder}"

		rm -rf "$of"/d
		rm -rf "$of"/i
		mkdir -p "$of"/d
		mkdir -p "$of"/i
		for ((i = 0; i <= N; i++)); do
			printf -v j "%07d" $i
			cp "$(find $of -name Orthogroup_Sequences)"/OG"${j}".fa "$of"/d
		done

		for i in $of/d/*.fa; do muscle -quiet -in $i -out $of/i/$(basename ${i%.fa}.afa); done

		_polap_log0 "Check: $of/i for alignment"
		_polap_log0 "find . -name Orthogroups.GeneCount.tsv | xargs open"
	fi

	if [[ "${_arg_menu[1]}" == "step6" ]] || [[ "${_arg_menu[1]}" == "infile" ]]; then
		local N="{_arg_bridge_min}" # -x 20
		local OG="${_polap_var_orthofinder}"
		local phylo="${_polap_var_phylogenies}"

		rm -rf "$phylo"
		mkdir -p "$phylo"/d
		mkdir -p "$phylo"/i
		for ((i = 0; i <= N; i++)); do
			printf -v j "%07d" $i
			cp "$(find $OG -name Orthogroup_Sequences)"/OG"${j}".fa "$phylo"/d
		done

		for i in $phylo/d/*.fa; do
			sed -i -e "s/\(>[a-zA-Z0-9]\+\)_[a-zA-Z0-9]\+/\1/" $i
		done

		for i in $phylo/d/*.fa; do
			# cat $i | perl -pe's/l\|([A-Za-z0-9]+)\.\d+_(\S+)/>$1/' | perl -pe's/>(target)_\S+/>$1/' > ${i%.fa}.fa1;
			cat $i | perl -pe's/^>[a-z]{3}\|([A-Za-z0-9_]+)\.{1}.+/>$1/' >${i%.fa}.fa1
		done

		for i in $phylo/d/*.fa1; do seqkit rmdup -n $i -o ${i%.fa1}.faa; done
		for i in $phylo/d/*.faa; do muscle -quiet -in $i -out ${i%.faa}.afa; done
		mkdir -p $phylo/s
		cp $phylo/d/*.afa $phylo/s
		iqtree -p $phylo/s -T AUTO -B 1000
	fi

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_polap_tree {
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

	# Print help message if requested
	help_message=$(
		cat <<HEREDOC
# Prepare files for IQ-TREE inference.
#
# Arguments:
#
# Inputs:
#
# Outputs:
#
# See:
Example: $(basename $0) ${_arg_menu[0]}
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
	fi

	while IFS=$'\t' read -r col1 col2 col3 col4 col5; do
		clean_col1="${col1%.*}"
		# Replace spaces in columns with underscores
		concatenated_column=$(echo "$col3"_"$col4"_"$col5" | tr ' ' '_')

		# Print the transformed line
		echo -e "$clean_col1\t$concatenated_column"
	done <"${_polap_var_ncbi_sequence_taxon_tsv}" >"${_polap_var_ncbi_sequence_taxon_map}"

	# nw_rename short_IDs.nw id2longname.map
	nw_rename "${_polap_var_phylogenies_tree}" \
		"${_polap_var_ncbi_sequence_taxon_map}" \
		>"${_polap_var_phylogenies_tree_taxa}"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_polap_ncbi {
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

	# Print help message if requested
	help_message=$(
		cat <<HEREDOC

Example: $(basename $0) ${_arg_menu[0]}
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	local _v=$(_polap_lib_ncbi-query-genome-size "${_arg_species}")
	_polap_log0 "Genome size of ${_arg_species}: ${_v}"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}
