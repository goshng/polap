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

# Define the get_family function
_taxonomy_get_ranks() {
	local taxid="$1" # Input TaxID
	local _tsv="$2"  #

	# Variable to store the family name
	local genus=""
	local family=""
	local order=""
	local class=""
	local phylum=""
	local kindom=""

	# Traverse the taxonomy hierarchy
	while [ "$taxid" != "1" ]; do
		# Get the current TaxID details from the taxonomy file
		entry=$(grep -P "^$taxid\t" "${_tsv}")
		if [ -z "$entry" ]; then
			echo "TaxID $taxid not found in the taxonomy database." >&2
			return 1
		fi

		rank=$(echo "$entry" | cut -f3) # Extract rank
		name=$(echo "$entry" | cut -f4) # Extract scientific name

		# Check if the rank is "family"
		if [ "$rank" == "genus" ]; then
			genus="$name"
		elif [ "$rank" == "family" ]; then
			family="$name"
		elif [ "$rank" == "order" ]; then
			order="$name"
		elif [ "$rank" == "class" ]; then
			class="$name"
		elif [ "$rank" == "phylum" ]; then
			phylum="$name"
		elif [ "$rank" == "kindom" ]; then
			kindom="$name"
		fi

		# Update TaxID to the Parent_TaxID
		taxid=$(echo "$entry" | cut -f2)
	done

	# Output the family name or indicate that it's not found
	if [ -n "$family" ]; then
		echo "$phylum $class $order $family $genus"
		return 0
	else
		echo "Family not found for the given TaxID." >&2
		return 1
	fi
}

# Define the get_family function
_taxonomy_get_family() {
	local taxid="$1" # Input TaxID
	local taxon="$2" # Input order or family
	local _tsv="$3"
	local family="" # Variable to store the family name

	# Traverse the taxonomy hierarchy
	while [ "$taxid" != "1" ]; do
		# Get the current TaxID details from the taxonomy file
		entry=$(grep -P "^$taxid\t" "${_tsv}")
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

_taxonomy_get_taxonomy_info() {
	local species_name="$1"
	local _tsv="$2" # taxonomy_full.tsv
	local taxid=""
	local phylum=""
	local class=""
	local order=""
	local family=""

	# Find the Taxonomy ID for the species name
	taxid=$(grep -P "\t${species_name}$" "${_tsv}" | cut -f1)

	if [ -z "$taxid" ]; then
		echo "Error: Species '$species_name' not found in the taxonomy database." >&2
		return 1
	fi

	# Traverse the taxonomy hierarchy
	while [ "$taxid" != "1" ]; do
		# Get the current TaxID details
		entry=$(grep -P "^$taxid\t" "${_tsv}")
		if [ -z "$entry" ]; then
			echo "Error: TaxID $taxid not found in the taxonomy database." >&2
			return 1
		fi

		rank=$(echo "$entry" | cut -f3) # Extract rank
		name=$(echo "$entry" | cut -f4) # Extract scientific name

		# Save the order and family if found
		if [ "$rank" == "phylum" ]; then
			phylum="$name"
		elif [ "$rank" == "class" ]; then
			class="$name"
		elif [ "$rank" == "order" ]; then
			order="$name"
		elif [ "$rank" == "family" ]; then
			family="$name"
		fi

		# Stop if both order and family are found
		if [[ -n "$phylum" ]] &&
			[[ -n "$class" ]] &&
			[[ -n "$order" ]] &&
			[[ -n "$family" ]]; then
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
	echo -e "${taxid}\t${species_name}\t${phylum}\t${class}\t${order}\t${family}"
}

_taxonomy_get_taxon_group() {
	local _classification_rank="$1"
	local _taxon_group="${_taxonomy_genus}"

	case "$_classification_rank" in
	'phylum')
		_taxon_group="${_taxonomy_phylum}"
		;;
	'class')
		_taxon_group="${_taxonomy_class}"
		;;
	'order')
		_taxon_group="${_taxonomy_order}"
		;;
	'family')
		_taxon_group="${_taxonomy_family}"
		;;
	*)
		_taxon_group="${_taxonomy_genus}"
		;;
	esac

	echo "$_taxon_group"
}

function _run_polap_taxonomy-assemble {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_/Menu_/)"

	source "$script_dir/polap-variables-common.sh"
	local _taxonomy_dir="${_arg_outdir}/taxonomy"
	local _assemble="${_taxonomy_dir}/01-assemble"
	local _sample="${_taxonomy_dir}/02-sample"
	local _geseq="${_taxonomy_dir}/03-geseq"
	local _orthofinder="${_taxonomy_dir}/04-orthofinder"
	local _phylogenies="${_taxonomy_dir}/05-phylogenies"
	local _tree="${_taxonomy_dir}/06-tree"

	if [[ -z "${_arg_steps_include}" ]]; then
		_arg_steps_include="0-4"
		_arg_steps_exclude=""
	fi

	_include="${_arg_steps_include}"
	_exclude="${_arg_steps_exclude}" # Optional range or list of steps to exclude
	_step_array=()
	_step_array=($(_polap_parse_steps "${_include}" "${_exclude}"))

	# main
	if _polap_contains_step 1 "${_step_array[@]}"; then
		_polap_log0 "  step 1: copy the mtDNA or ptDNA assembly in FASTA format"
	fi

	# Disable debugging if previously enabled
	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_polap_taxonomy-species {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_/Menu_/)"

	source "$script_dir/polap-variables-common.sh"
	local _taxonomy_dir="${_arg_outdir}/taxonomy"
	local _assemble="${_taxonomy_dir}/01-assemble"
	local _sample="${_taxonomy_dir}/02-sample"
	local _geseq="${_taxonomy_dir}/03-geseq"
	local _orthofinder="${_taxonomy_dir}/04-orthofinder"
	local _phylogenies="${_taxonomy_dir}/05-phylogenies"
	local _tree="${_taxonomy_dir}/06-tree"

	# input1: species_list.txt
	#
	if [[ ! -f "taxdump.tar.gz" ]]; then
		wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
		tar -xvzf taxdump.tar.gz
	fi

	# Extract TaxID for each species from names.dmp
	grep -P '\tscientific name\t' names.dmp | cut -f1,3 >taxid_name.tsv
	# awk -F'\t\\|\t' '$4 == "scientific name"' names.dmp | awk -F'\t\\|\t' '{print $1 "\t" $2}' >taxid_name.tsv
	# Extract Taxonomy Hierarchy from nodes.dmp
	awk -F'\t\\|\t' '{print $1 "\t" $3 "\t" $5}' nodes.dmp >taxid_parent_rank.tsv

	# Create a temporary lookup from names.dmp
	# awk -F'\t\|\t' '$4 == "scientific name"' names.dmp | awk -F'\t\|\t' '{print $1 "\t" $2}' >taxid_name.tsv
	# awk -F'\t\|\t' '{print $1 "\t" $3 "\t" $5}' nodes.dmp >taxid_parent_rank.tsv

	# Output header
	echo -e "Species\tTaxID\tKingdom\tPhylum\tClass\tOrder\tFamily" >taxonomy_output.tsv

	# Function to trace taxonomy lineage
	get_lineage() {
		local taxid=$1
		local lineage="NA\tNA\tNA\tNA\tNA"
		local kingdom phylum class order family
		while [ "$taxid" != "1" ] && [ -n "$taxid" ]; do
			rank=$(awk -v id="$taxid" '$1 == id {print $3}' taxid_parent_rank.tsv)
			parent=$(awk -v id="$taxid" '$1 == id {print $2}' taxid_parent_rank.tsv)
			name=$(awk -v id="$taxid" '$1 == id {print $2}' taxid_name.tsv)

			case $rank in
			kingdom) kingdom=$name ;;
			phylum) phylum=$name ;;
			class) class=$name ;;
			order) order=$name ;;
			family) family=$name ;;
			esac
			taxid=$parent
		done
		echo -e "${kingdom:-NA}\t${phylum:-NA}\t${class:-NA}\t${order:-NA}\t${family:-NA}"
	}

	# Process species list
	while IFS= read -r species; do
		taxid=$(awk -F'\t' -v sp="$species" '$2 == sp {print $1}' taxid_name.tsv)

		if [ -z "$taxid" ]; then
			echo -e "$species\tNot_Found\tNA\tNA\tNA\tNA\tNA" >>taxonomy_output.tsv
			continue
		fi

		lineage=$(get_lineage "$taxid")
		echo -e "$species\t$taxid\t$lineage" >>taxonomy_output.tsv
	done <species_list.txt

	# Disable debugging if previously enabled
	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

#
function _taxonomy_timer_set {
	# Progress
	# Get terminal width
	_terminal_width=$(tput cols)
	_polap_set_start_time
	_total_iterations=$(cat "${_sequence_outgroup_tsv}" | wc -l)
	return 0
}

function _taxonomy_timer_count {
	time_per_iteration=$((($(date +%s) - _polap_var_start_time) / i))
	remaining_iterations=$((_total_iterations - i - 1))
	remaining_time=$((remaining_iterations * time_per_iteration))
	status="    iteration $((i + 1))/$_total_iterations, remaining time: $(_polap_get_time_format ${remaining_time})"
	printf "\r%-${_terminal_width}s" " " >&3
	_polap_log0_ne "\r$status"
}

function _taxonomy-extract-tsv-from-tinyxml {
	local _sequence_tsv="${1}"
	local _sequence_xml="${2}"

	# Write the TSV header
	echo -e "accver\ttaxid\torgname\tdefline\tlength" >"${_sequence_tsv}"

	# Extract the fields and write to the TSV file
	xmllint --format "${_sequence_xml}" |
		awk '
/<TSeq>/,/<\/TSeq>/ {
    if ($1 ~ /<TSeq_accver>/) accver=gensub(/.*<TSeq_accver>(.*)<\/TSeq_accver>.*/, "\\1", "g")
    if ($1 ~ /<TSeq_taxid>/) taxid=gensub(/.*<TSeq_taxid>(.*)<\/TSeq_taxid>.*/, "\\1", "g")
    if ($1 ~ /<TSeq_orgname>/) orgname=gensub(/.*<TSeq_orgname>(.*)<\/TSeq_orgname>.*/, "\\1", "g")
    if ($1 ~ /<TSeq_defline>/) defline=gensub(/.*<TSeq_defline>(.*)<\/TSeq_defline>.*/, "\\1", "g")
    if ($1 ~ /<TSeq_length>/) len=gensub(/.*<TSeq_length>(.*)<\/TSeq_length>.*/, "\\1", "g")
    if ($1 ~ /<\/TSeq>/) {
        printf("%s\t%s\t%s\t%s\t%s\n", accver, taxid, orgname, defline, len) >> "'${_sequence_tsv}'"
        accver=""; taxid=""; orgname=""; defline=""; len=""
    }
}'
}

function _taxonomy-extract-taxon-from-tsv-using-taxonomy {
	local _terminal_width
	local _sequence_tsv="${1}"
	local _taxonomy_full_tsv="${2}"
	local _sequence_taxon_tsv="${3}"

	# _taxonomy_timer_set
	_terminal_width=$(tput cols)
	_polap_set_start_time
	_total_iterations=$(cat "${_sequence_tsv}" | wc -l)

	_polap_log1 "    It takes time ... please, wait ... total ${_total_iterations}"
	i=0
	rm -f "${_sequence_taxon_tsv}"
	# Add Class-level and Order-level TaxIDs
	tail -n +2 "${_sequence_tsv}" |
		while IFS=$'\t' read -r acc taxid organism desc len; do
			i=$((i + 1))
			# Fetch taxonomy details for the current TaxID
			read -r _phylum _class _order _family _genus \
				<<<$(_taxonomy_get_ranks "$taxid" "${_taxonomy_full_tsv}")

			# local _genus=$(_taxonomy_get_family "$taxid" genus "${_taxonomy_full_tsv}")
			# local _family=$(_taxonomy_get_family "$taxid" family "${_taxonomy_full_tsv}")
			# local _order=$(_taxonomy_get_family "$taxid" order "${_taxonomy_full_tsv}")
			# local _class=$(_taxonomy_get_family "$taxid" class "${_taxonomy_full_tsv}")
			# local _phylum=$(_taxonomy_get_family "$taxid" phylum "${_taxonomy_full_tsv}")
			echo -e "${acc}\t${taxid}\t${organism}\t${_phylum}\t${_class}\t${_order}\t${_family}\t${_genus}" \
				>>"${_sequence_taxon_tsv}"
			_phylum=""
			_class=""
			_order=""
			_family=""
			_genus=""

			# Progress
			# alculate elapsed time and remaining time
			# _taxonomy_timer_count
			time_per_iteration=$((($(date +%s) - _polap_var_start_time) / i))
			remaining_iterations=$((_total_iterations - i - 1))
			remaining_time=$((remaining_iterations * time_per_iteration))
			status="    iteration $((i + 1))/$_total_iterations, remaining time: $(_polap_get_time_format ${remaining_time})"
			# if [ "${_arg_verbose}" -eq "1" ]; then
			printf "\r%-${_terminal_width}s" " " >&3
			# fi
			_polap_log1_ne "\r$status"
			# _polap_log0 "$status"
		done
	# done < <(tail -n +2 "${_sequence_tsv}")

	_polap_log0 ""
}

# select mt or ptDNA reference
function _run_polap_taxonomy-reference {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_/Menu_/)"

	source "$script_dir/polap-variables-common.sh"
	local _taxonomy_dir="${_arg_outdir}/taxonomy"
	local _reference="${_taxonomy_dir}/01-reference"
	local _assemble="${_taxonomy_dir}/01-assemble"
	local _sample="${_taxonomy_dir}/02-sample"
	local _geseq="${_taxonomy_dir}/03-geseq"
	local _orthofinder="${_taxonomy_dir}/04-orthofinder"
	local _phylogenies="${_taxonomy_dir}/05-phylogenies"
	local _tree="${_taxonomy_dir}/06-tree"

	_polap_log3_cmd mkdir -p "${_reference}"

	if [[ -z "${_arg_steps_include}" ]]; then
		_arg_steps_include="1-9"
		_arg_steps_exclude=""
	fi

	_include="${_arg_steps_include}"
	_exclude="${_arg_steps_exclude}" # Optional range or list of steps to exclude
	_step_array=()
	_step_array=($(_polap_parse_steps "${_include}" "${_exclude}"))

	# main
	_polap_log0 "Taxonomy sample of ${_arg_species}"

	local _taxonomy_kingdom
	local _taxonomy_phylum
	local _taxonomy_class
	local _taxonomy_order
	local _taxonomy_family
	local _taxonomy_genus
	local _taxonomy_xml
	local _taxonomy_ranks_species
	_taxonomy_xml="${_reference}/taxonomy.xml"
	_taxonomy_ranks_species="${_reference}/taxonomy_ranks_species.txt"
	if _polap_contains_step 1 "${_step_array[@]}"; then
		_polap_log0 "  step 1: determine taxonomy ranks -> ${_taxonomy_ranks_species}"
		_polap_log1 "    input1: species: ${_arg_species}"
		_polap_log1 "    output1: ${_taxonomy_xml}"
		_polap_log1 "    output2: ${_taxonomy_ranks_species}"

		if [[ -s "${_taxonomy_xml}" ]] && [[ "${_arg_redo}" == "off" ]]; then
			_polap_log1 "    found1: ${_taxonomy_xml}"
		else
			esearch -db taxonomy -query "${_arg_species}" |
				efetch -format xml >"${_taxonomy_xml}"
		fi
		_taxonomy_kingdom=$(xtract -input "${_taxonomy_xml}" -pattern Taxon -block "LineageEx/Taxon" -if Rank -equals "kingdom" -element ScientificName)
		_taxonomy_phylum=$(xtract -input "${_taxonomy_xml}" -pattern Taxon -block "LineageEx/Taxon" -if Rank -equals "phylum" -element ScientificName)
		_taxonomy_class=$(xtract -input "${_taxonomy_xml}" -pattern Taxon -block "LineageEx/Taxon" -if Rank -equals "class" -element ScientificName)
		_taxonomy_order=$(xtract -input "${_taxonomy_xml}" -pattern Taxon -block "LineageEx/Taxon" -if Rank -equals "order" -element ScientificName)
		_taxonomy_family=$(xtract -input "${_taxonomy_xml}" -pattern Taxon -block "LineageEx/Taxon" -if Rank -equals "family" -element ScientificName)
		_taxonomy_genus=$(xtract -input "${_taxonomy_xml}" -pattern Taxon -block "LineageEx/Taxon" -if Rank -equals "genus" -element ScientificName)
		printf "%s\n%s\n%s\n%s\n%s\n%s\n" \
			"kingdom: ${_taxonomy_kingdom}" \
			"phylum: ${_taxonomy_phylum}" \
			"class: ${_taxonomy_class}" \
			"order: ${_taxonomy_order}" \
			"family: ${_taxonomy_family}" \
			"genus: ${_taxonomy_genus}" \
			>"${_taxonomy_ranks_species}"
		_polap_log2_cat "${_taxonomy_ranks_species}"
	fi

	local _taxdump="${_reference}/taxdump.tar.gz"
	local _taxdump_dir="${_reference}/taxdump"
	mkdir -p "${_taxdump_dir}"
	if _polap_contains_step 2 "${_step_array[@]}"; then
		_polap_log0 "  step 2: download taxonomy dump file"
		_polap_log1 "    output1: ${_taxdump}"

		if [[ -s "${_taxdump}" ]] && [[ "${_arg_redo}" == "off" ]]; then
			_polap_log1 "    found1: ${_taxdump}"
		else
			rm -f "${_taxdump}"
			wget -q -O "${_taxdump}" \
				ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
			tar zxf "${_taxdump}" -C "${_taxdump_dir}"
		fi
	fi

	local _taxonomy_full_tsv="${_reference}/taxonomy_full.tsv"
	if _polap_contains_step 3 "${_step_array[@]}"; then
		_polap_log0 "  step 3: extract TSV from the taxonomy dump file"
		_polap_log1 "    input1: ${_taxdump}"
		_polap_log1 "    output1: ${_taxonomy_full_tsv}"

		if [[ -s "${_taxonomy_full_tsv}" ]] && [[ "${_arg_redo}" == "off" ]]; then
			_polap_log1 "    found1: ${_taxonomy_full_tsv}"
		else

			# Extract relevant data from nodes.dmp and names.dmp
			awk -F "\t|\t" 'BEGIN {OFS="\t"} {print $1, $3, $5}' \
				"${_taxdump_dir}"/nodes.dmp >"${_taxdump_dir}"/nodes.tsv
			awk -F "\t|\t" 'BEGIN {OFS="\t"} $7 == "scientific name" {print $1, $3}' \
				"${_taxdump_dir}"/names.dmp >"${_taxdump_dir}"/names.tsv

			# Merge the two files on TaxID
			join -t $'\t' -1 1 -2 1 \
				<(sort -k1,1 "${_taxdump_dir}"/nodes.tsv) \
				<(sort -k1,1 "${_taxdump_dir}"/names.tsv) \
				>"${_taxonomy_full_tsv}"
		fi

		_polap_log3_head "${_taxonomy_full_tsv}"
	fi

	local _taxon_ingroup
	local _sequence_ingroup_xml="${_reference}/sequence.ingroup.fasta.xml"
	local _accessions_ingroup="${_reference}/accessions_ingroup.txt"
	if _polap_contains_step 5 "${_step_array[@]}"; then
		_polap_log0 "  step 5: download ingroup accessions -> ${_accessions_ingroup}"
		_polap_log1 "    input1: classification rank: ${_arg_taxonomy_rank_ingroup}"
		_polap_log1 "    output1: ${_accessions_ingroup}"

		if [[ -z "${_arg_taxonomy_rank_ingroup}" ]]; then
			die "ERROR: no value for option --taxonomy-rank-ingroup"
		fi

		_taxon_ingroup=$(_taxonomy_get_taxon_group "${_arg_taxonomy_rank_ingroup}")

		local _ncbi_query_ingroup
		if [[ "${_arg_plastid}" == "off" ]]; then
			_ncbi_query_ingroup="(mitochondrion[Title] AND complete[Title] AND genome[Title]) AND ${_taxon_ingroup}[Organism]"
		else
			_ncbi_query_ingroup="(chloroplast[Title] AND complete[Title] AND genome[Title]) AND ${_taxon_ingroup}[Organism]"
		fi
		if [[ -s "${_accessions_ingroup}" ]] && [[ "${_arg_redo}" == "off" ]]; then
			_polap_log1 "    found1: ${_accessions_ingroup}"
		else
			esearch \
				-db nuccore \
				-query "${_ncbi_query_ingroup}" |
				efetch -format acc >"${_accessions_ingroup}"
		fi

		_polap_log2_head "${_accessions_ingroup}"

		_polap_log0 "Download ${_sequence_ingroup_xml}"
		_polap_log0 "  1. visit https://www.ncbi.nlm.nih.gov/nuccore"
		_polap_log0 "  2. search: ${_ncbi_query_ingroup}"
		_polap_log0 "  3. use Send to: button to select TinySeq XML to save."
		_polap_log0 "  4. download flie: sequence.fasta.xml"
		_polap_log0 "  5. move it to ${_sequence_ingroup_xml}"
		_polap_log0 "  $ mv sequence.fasta.xml ${_sequence_ingroup_xml}"
	fi

	# pause

	# list unique sequences by length
	#
	local _ingroup_fasta="${_reference}/ingroup.fa"
	local _inref_fasta="${_reference}/inref.fa"
	local _inref2_fasta="${_reference}/inref2.fa"
	if _polap_contains_step 8 "${_step_array[@]}"; then
		_polap_log0 "  step 8: extract sequences"
		_polap_log1 "    input1: ${_accessions_ingroup}"
		_polap_log1 "    output1: ${_inref2_fasta}"

		esearch -db nucleotide -query \
			"$(awk '{printf "%s OR ", $1}' ${_accessions_ingroup} | sed 's/ OR $//')" |
			efetch -format fasta >"${_ingroup_fasta}"

		read first second <<<"${_arg_species}"

		seqkit grep -vnr -p "${second}" \
			"${_ingroup_fasta}" \
			-o "${_inref_fasta}"
		seqkit rmdup -s "${_inref_fasta}" -o "${_inref2_fasta}"

		_polap_log0 "-----------------------------------------------"
		_polap_log0 "List of reference sequences for ${_arg_species}"
		_polap_log0 "-----------------------------------------------"
		seqkit seq -n "${_inref2_fasta}" >&3
	fi

	local _selected_accessions="${_reference}/selected_accessions.txt"
	if _polap_contains_step 9 "${_step_array[@]}"; then
		_polap_log0 "  step 9: choose one per genus"
		_polap_log1 "    input1: ${_inref2_fasta}"
		_polap_log1 "    output1: ${_arg_outdir}/ptgaul-reference.fa"

		# Input FASTA file
		local input_fasta="${_inref2_fasta}"
		local output_fasta="${_arg_outdir}/ptgaul-reference.fa"

		# Count the number of sequences using seqkit
		local sequence_count=$(seqkit stats $input_fasta | awk 'NR==2 {print $4}')

		# Decision based on the number of sequences
		if [ "$sequence_count" -le 10 ]; then
			_polap_log1 "    The FASTA file ${_inref2_fasta} contains 10 or fewer sequences."
			cp -p "${_inref2_fasta}" "${output_fasta}"
		else
			if [[ "${_arg_taxonomy_rank_ingroup}" == "family" ]]; then

				# Extract accession and genus, then randomly select one per genus
				seqkit seq -n $input_fasta |
					awk '{print $1, $2}' |
					sort |
					awk '!seen[$2]++' |
					shuf |
					awk '{print $1}' >"${_selected_accessions}"

			else

				# Select one species only per species
				# Extract accession and genus, then randomly select one per genus
				seqkit seq -n $input_fasta |
					awk '{print $1, $2, $3}' |
					awk '{species=$2"_"$3; print $1, species}' |
					sort |
					awk '!seen[$2]++' |
					shuf |
					awk '{print $1}' >"${_selected_accessions}"

			fi
			# Filter the original FASTA file to include only selected accessions
			seqkit grep -f "${_selected_accessions}" $input_fasta -o $output_fasta
		fi

		_polap_log0 "-----------------------------------------------"
		_polap_log0 "List of reference sequences for ${_arg_species}"
		_polap_log0 "-----------------------------------------------"
		seqkit seq -n "${output_fasta}" >&3
	fi

	# Disable debugging if previously enabled
	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_polap_taxonomy-sample {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_/Menu_/)"

	source "$script_dir/polap-variables-common.sh"
	local _taxonomy_dir="${_arg_outdir}/taxonomy"
	local _assemble="${_taxonomy_dir}/01-assemble"
	local _sample="${_taxonomy_dir}/02-sample"
	local _geseq="${_taxonomy_dir}/03-geseq"
	local _orthofinder="${_taxonomy_dir}/04-orthofinder"
	local _phylogenies="${_taxonomy_dir}/05-phylogenies"
	local _tree="${_taxonomy_dir}/06-tree"

	_polap_log3_cmd mkdir -p "${_sample}"

	if [[ -z "${_arg_steps_include}" ]]; then
		_arg_steps_include="1-4"
		_arg_steps_exclude=""
	fi

	_include="${_arg_steps_include}"
	_exclude="${_arg_steps_exclude}" # Optional range or list of steps to exclude
	_step_array=()
	_step_array=($(_polap_parse_steps "${_include}" "${_exclude}"))

	# main
	_polap_log0 "Taxonomy sample of ${_arg_species}"

	local _taxonomy_kingdom
	local _taxonomy_phylum
	local _taxonomy_class
	local _taxonomy_order
	local _taxonomy_family
	local _taxonomy_genus
	local _taxonomy_xml
	local _taxonomy_ranks_species
	_taxonomy_xml="${_sample}/taxonomy.xml"
	_taxonomy_ranks_species="${_sample}/taxonomy_ranks_species.txt"
	if _polap_contains_step 1 "${_step_array[@]}"; then
		_polap_log0 "  step 1: determine taxonomy ranks -> ${_taxonomy_ranks_species}"
		_polap_log1 "    input1: species: ${_arg_species}"
		_polap_log1 "    output1: ${_taxonomy_xml}"
		_polap_log1 "    output2: ${_taxonomy_ranks_species}"

		if [[ -s "${_taxonomy_xml}" ]] && [[ "${_arg_redo}" == "off" ]]; then
			_polap_log1 "    found1: ${_taxonomy_xml}"
		else
			esearch -db taxonomy -query "${_arg_species}" |
				efetch -format xml >"${_taxonomy_xml}"
		fi
		_taxonomy_kingdom=$(xtract -input "${_taxonomy_xml}" -pattern Taxon -block "LineageEx/Taxon" -if Rank -equals "kingdom" -element ScientificName)
		_taxonomy_phylum=$(xtract -input "${_taxonomy_xml}" -pattern Taxon -block "LineageEx/Taxon" -if Rank -equals "phylum" -element ScientificName)
		_taxonomy_class=$(xtract -input "${_taxonomy_xml}" -pattern Taxon -block "LineageEx/Taxon" -if Rank -equals "class" -element ScientificName)
		_taxonomy_order=$(xtract -input "${_taxonomy_xml}" -pattern Taxon -block "LineageEx/Taxon" -if Rank -equals "order" -element ScientificName)
		_taxonomy_family=$(xtract -input "${_taxonomy_xml}" -pattern Taxon -block "LineageEx/Taxon" -if Rank -equals "family" -element ScientificName)
		_taxonomy_genus=$(xtract -input "${_taxonomy_xml}" -pattern Taxon -block "LineageEx/Taxon" -if Rank -equals "genus" -element ScientificName)
		printf "%s\n%s\n%s\n%s\n%s\n%s\n" \
			"kingdom: ${_taxonomy_kingdom}" \
			"phylum: ${_taxonomy_phylum}" \
			"class: ${_taxonomy_class}" \
			"order: ${_taxonomy_order}" \
			"family: ${_taxonomy_family}" \
			"genus: ${_taxonomy_genus}" \
			>"${_taxonomy_ranks_species}"
		_polap_log2_cat "${_taxonomy_ranks_species}"
	fi

	local _taxdump="${_sample}/taxdump.tar.gz"
	local _taxdump_dir="${_sample}/taxdump"
	mkdir -p "${_taxdump_dir}"
	if _polap_contains_step 2 "${_step_array[@]}"; then
		_polap_log0 "  step 2: download taxonomy dump file"
		_polap_log1 "    output1: ${_taxdump}"

		if [[ -s "${_taxdump}" ]] && [[ "${_arg_redo}" == "off" ]]; then
			_polap_log1 "    found1: ${_taxdump}"
		else
			rm -f "${_taxdump}"
			wget -O "${_taxdump}" \
				ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
			tar zxf "${_taxdump}" -C "${_taxdump_dir}"
		fi
	fi

	local _taxonomy_full_tsv="${_sample}/taxonomy_full.tsv"
	if _polap_contains_step 3 "${_step_array[@]}"; then
		_polap_log0 "  step 3: extract TSV from the taxonomy dump file"
		_polap_log1 "    input1: ${_taxdump}"
		_polap_log1 "    output1: ${_taxonomy_full_tsv}"

		if [[ -s "${_taxonomy_full_tsv}" ]] && [[ "${_arg_redo}" == "off" ]]; then
			_polap_log1 "    found1: ${_taxonomy_full_tsv}"
		else

			# Extract relevant data from nodes.dmp and names.dmp
			awk -F "\t|\t" 'BEGIN {OFS="\t"} {print $1, $3, $5}' \
				"${_taxdump_dir}"/nodes.dmp >"${_taxdump_dir}"/nodes.tsv
			awk -F "\t|\t" 'BEGIN {OFS="\t"} $7 == "scientific name" {print $1, $3}' \
				"${_taxdump_dir}"/names.dmp >"${_taxdump_dir}"/names.tsv

			# Merge the two files on TaxID
			join -t $'\t' -1 1 -2 1 \
				<(sort -k1,1 "${_taxdump_dir}"/nodes.tsv) \
				<(sort -k1,1 "${_taxdump_dir}"/names.tsv) \
				>"${_taxonomy_full_tsv}"
		fi

		_polap_log3_head "${_taxonomy_full_tsv}"
	fi

	local _taxon_allgroup
	local _accessions_allgroup="${_sample}/accessions_allgroup.txt"
	if _polap_contains_step 4 "${_step_array[@]}"; then
		_polap_log0 "  step 4: download all or ingroup/outgroup accessions -> ${_accessions_allgroup}"
		_polap_log1 "    input1: classification rank: ${_arg_taxonomy_rank_allgroup}"
		_polap_log1 "    output1: ${_accessions_allgroup}"

		_taxon_allgroup=$(_taxonomy_get_taxon_group "${_arg_taxonomy_rank_allgroup}")

		local _ncbi_query_allgroup
		_ncbi_query_allgroup="(mitochondrion[Title] AND complete[Title] AND genome[Title]) AND ${_taxon_allgroup}[Organism]"
		if [[ -s "${_accessions_allgroup}" ]] && [[ "${_arg_redo}" == "off" ]]; then
			_polap_log1 "    found1: ${_accessions_allgroup}"
		else
			esearch \
				-db nuccore \
				-query "${_ncbi_query_allgroup}" |
				efetch -format acc >"${_accessions_allgroup}"
		fi
		_polap_log3_head "${_accessions_allgroup}"
	fi

	local _taxon_ingroup
	local _sequence_ingroup_xml="${_sample}/sequence.ingroup.fasta.xml"
	local _accessions_ingroup="${_sample}/accessions_ingroup.txt"
	if _polap_contains_step 5 "${_step_array[@]}"; then
		_polap_log0 "  step 5: download ingroup accessions -> ${_accessions_ingroup}"
		_polap_log1 "    input1: classification rank: ${_arg_taxonomy_rank_ingroup}"
		_polap_log1 "    output1: ${_accessions_ingroup}"

		_taxon_ingroup=$(_taxonomy_get_taxon_group "${_arg_taxonomy_rank_ingroup}")

		local _ncbi_query_ingroup
		_ncbi_query_ingroup="(mitochondrion[Title] AND complete[Title] AND genome[Title]) AND ${_taxon_ingroup}[Organism]"
		if [[ -s "${_accessions_ingroup}" ]] && [[ "${_arg_redo}" == "off" ]]; then
			_polap_log1 "    found1: ${_accessions_ingroup}"
		else
			esearch \
				-db nuccore \
				-query "${_ncbi_query_ingroup}" |
				efetch -format acc >"${_accessions_ingroup}"
		fi

		_polap_log2_head "${_accessions_ingroup}"

		_polap_log0 "Download ${_sequence_ingroup_xml}"
		_polap_log0 "  1. visit https://www.ncbi.nlm.nih.gov/nuccore"
		_polap_log0 "  2. search: ${_ncbi_query_ingroup}"
		_polap_log0 "  3. use Send to: button to select TinySeq XML to save."
		_polap_log0 "  4. download flie: sequence.fasta.xml"
		_polap_log0 "  5. move it to ${_sequence_ingroup_xml}"
		_polap_log0 "  $ mv sequence.fasta.xml ${_sequence_ingroup_xml}"
	fi

	local _taxon_outgroup
	local _sequence_outgroup_xml="${_sample}/sequence.outgroup.fasta.xml"
	local _accessions_outgroup="${_sample}/accessions_outgroup.txt"
	if _polap_contains_step 6 "${_step_array[@]}"; then
		_polap_log0 "  step 6: download outgroup accessions without ingroup -> ${_accessions_outgroup}"
		_polap_log1 "    input1: classification rank of all: ${_arg_taxonomy_rank_allgroup}"
		_polap_log1 "    input2: classification rank of ingroup: ${_arg_taxonomy_rank_ingroup}"
		_polap_log1 "    output1: ${_accessions_outgroup}"

		_taxon_ingroup=$(_taxonomy_get_taxon_group "${_arg_taxonomy_rank_ingroup}")
		_taxon_allgroup=$(_taxonomy_get_taxon_group "${_arg_taxonomy_rank_allgroup}")

		local _ncbi_query_outgroup
		_ncbi_query_outgroup="(mitochondrion[Title] AND complete[Title] AND genome[Title]) AND ${_taxon_allgroup}[Organism] NOT ${_taxon_ingroup}[Organism]"
		if [[ -s "${_accessions_outgroup}" ]] && [[ "${_arg_redo}" == "off" ]]; then
			_polap_log1 "    found1: ${_accessions_outgroup}"
		else
			esearch \
				-db nuccore \
				-query "${_ncbi_query_outgroup}" |
				efetch -format acc >"${_accessions_outgroup}"
		fi
		_polap_log3_head "${_accessions_outgroup}"

		_polap_log0 "Download ${_sequence_outgroup_xml}"
		_polap_log0 "  1. visit https://www.ncbi.nlm.nih.gov/nuccore"
		_polap_log0 "  2. search: ${_ncbi_query_outgroup}"
		_polap_log0 "  3. use Send to: button to select TinySeq XML to save."
		_polap_log0 "  4. download flie: sequence.fasta.xml"
		_polap_log0 "  5. move it to ${_sequence_outgroup_xml}"
		_polap_log0 "  $ mv sequence.fasta.xml ${_sequence_outgroup_xml}"
	fi

	# pause

	local _sequence_ingroup_xml="${_sample}/sequence.ingroup.fasta.xml"
	local _sequence_ingroup_tsv="${_sample}/sequence.ingroup.tsv"
	if _polap_contains_step 7 "${_step_array[@]}"; then
		_polap_log0 "  step 7: extract taxonomy info from the tinyseq XML -> ${_sequence_ingroup_tsv}"
		_polap_log1 "    input1: ingroup mtDNA xml: ${_sequence_ingroup_xml}"
		_polap_log1 "    output1: ${_sequence_ingroup_tsv}"

		if [[ ! -s "${_sequence_ingroup_tsv}" ]] || [[ "${_arg_redo}" == "on" ]]; then
			_taxonomy-extract-tsv-from-tinyxml \
				"${_sequence_ingroup_tsv}" \
				"${_sequence_ingroup_xml}"
		fi
		_polap_log3_head "${_sequence_ingroup_tsv}"
	fi

	local _sequence_outgroup_xml="${_sample}/sequence.outgroup.fasta.xml"
	local _sequence_outgroup_tsv="${_sample}/sequence.outgroup.tsv"
	if _polap_contains_step 8 "${_step_array[@]}"; then
		_polap_log0 "  step 8: extract taxonomy info from the tinyseq XML -> ${_sequence_outgroup_tsv}"
		_polap_log1 "    input1: outgroup mtDNA xml: ${_sequence_outgroup_xml}"
		_polap_log1 "    output1: ${_sequence_outgroup_tsv}"

		if [[ ! -s "${_sequence_outgroup_tsv}" ]] || [[ "${_arg_redo}" == "on" ]]; then
			_taxonomy-extract-tsv-from-tinyxml \
				"${_sequence_outgroup_tsv}" \
				"${_sequence_outgroup_xml}"
		fi
		_polap_log3_head "${_sequence_outgroup_tsv}"
	fi

	local _sequence_ingroup_taxon_tsv="${_sample}/sequence.ingroup.taxon.tsv"
	if _polap_contains_step 9 "${_step_array[@]}"; then
		_polap_log0 "  step 9: select ingroup taxonomy info -> ${_sequence_ingroup_taxon_tsv}"
		_polap_log1 "    input1: ${_sequence_ingroup_tsv}"
		_polap_log1 "    input2: ${_taxonomy_full_tsv}"
		_polap_log1 "    output1: ${_sequence_ingroup_taxon_tsv}"

		if [[ -s "${_sequence_ingroup_taxon_tsv}" ]] && [[ "${_arg_redo}" == "off" ]]; then
			_polap_log1 "    found1: ${_sequence_ingroup_taxon_tsv}"
		else
			_taxonomy-extract-taxon-from-tsv-using-taxonomy \
				"${_sequence_ingroup_tsv}" \
				"${_taxonomy_full_tsv}" \
				"${_sequence_ingroup_taxon_tsv}"
		fi

		_polap_log3_head "${_sequence_ingroup_taxon_tsv}"
	fi

	local _sequence_outgroup_taxon_tsv="${_sample}/sequence.outgroup.taxon.tsv"
	if _polap_contains_step 10 "${_step_array[@]}"; then
		_polap_log0 "  step 10: select outgroup taxonomy info -> ${_sequence_outgroup_taxon_tsv}"
		_polap_log1 "    input1: ${_sequence_outgroup_tsv}"
		_polap_log1 "    input2: ${_taxonomy_full_tsv}"
		_polap_log1 "    output1: ${_sequence_outgroup_taxon_tsv}"

		if [[ -s "${_sequence_outgroup_taxon_tsv}" ]] && [[ "${_arg_redo}" == "off" ]]; then
			_polap_log1 "    found1: ${_sequence_outgroup_taxon_tsv}"
		else
			_taxonomy-extract-taxon-from-tsv-using-taxonomy \
				"${_sequence_outgroup_tsv}" \
				"${_taxonomy_full_tsv}" \
				"${_sequence_outgroup_taxon_tsv}"
		fi

		_polap_log3_head "${_sequence_outgroup_taxon_tsv}"
	fi

	local _accessions_ingroup_sampled="${_sample}/accesssions_ingroup_sampled.txt"
	if _polap_contains_step 11 "${_step_array[@]}"; then
		_polap_log0 "  step 11: sample taxa in the ingroup"
		_polap_log1 "    input1: ingroup size: ${_arg_taxonomy_ingroup_size}"
		_polap_log1 "    input2: rank sample size: ${_arg_taxonomy_sample_size_per_rank}"
		_polap_log1 "    input3: ingroup classification rank: ${_arg_taxonomy_rank_ingroup_sample}"
		_polap_log1 "    input4: ${_sequence_ingroup_taxon_tsv}"
		_polap_log1 "    output1: ${_accessions_ingroup_sampled}"

		_polap_log3_pipe "Rscript ${script_dir}/run-polap-r-taxonomy.R \
			sample \
			--size ${_arg_taxonomy_ingroup_size} \
			-n ${_arg_taxonomy_sample_size_per_rank} \
			-r ${_arg_taxonomy_rank_ingroup_sample} \
			-t ${_sequence_ingroup_taxon_tsv} \
			-o ${_accessions_ingroup_sampled} \
		  >${_polap_output_dest} 2>&1"

		_polap_log3_head "${_accessions_ingroup_sampled}"
	fi

	local _accessions_outgroup_sampled="${_sample}/accesssions_outgroup_sampled.txt"
	if _polap_contains_step 12 "${_step_array[@]}"; then
		_polap_log0 "  step 12: sample taxa in the outgroup"
		_polap_log1 "    input1: outgroup size: ${_arg_taxonomy_outgroup_size}"
		_polap_log1 "    input2: rank sample size: ${_arg_taxonomy_sample_size_per_rank}"
		_polap_log1 "    input3: sample rank: ${_arg_taxonomy_rank_outgroup_sample}"
		_polap_log1 "    input4: ${_sequence_outgroup_taxon_tsv}"
		_polap_log1 "    output1: ${_accessions_outgroup_sampled}"

		_polap_log3_pipe "Rscript ${script_dir}/run-polap-r-taxonomy.R \
			sample \
			--size ${_arg_taxonomy_outgroup_size} \
			-n ${_arg_taxonomy_sample_size_per_rank} \
			-r ${_arg_taxonomy_rank_outgroup_sample} \
			-t ${_sequence_outgroup_taxon_tsv} \
			-o ${_accessions_outgroup_sampled} \
		  >${_polap_output_dest} 2>&1"

		_polap_log3_head "${_accessions_outgroup_sampled}"
	fi

	# Disable debugging if previously enabled
	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_polap_taxonomy-geseq {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_/Menu_/)"

	local _geseq

	source "$script_dir/polap-variables-common.sh"
	local _taxonomy_dir="${_arg_outdir}/taxonomy"
	local _assemble="${_taxonomy_dir}/01-assemble"
	local _sample="${_taxonomy_dir}/02-sample"
	local _geseq="${_taxonomy_dir}/03-geseq"
	local _orthofinder="${_taxonomy_dir}/04-orthofinder"
	local _phylogenies="${_taxonomy_dir}/05-phylogenies"
	local _tree="${_taxonomy_dir}/06-tree"

	local _input_short_reads="${_arg_outdir}/s1.fq"

	# Display help message
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	# main
	#
	# select default steps
	if [[ -z "${_arg_steps_include}" ]]; then
		_arg_steps_include="1-6"
		_arg_steps_exclude=""
	fi
	_include="${_arg_steps_include}"
	_exclude="${_arg_steps_exclude}" # Optional range or list of steps to exclude
	_step_array=()
	_step_array=($(_polap_parse_steps "${_include}" "${_exclude}"))

	# main
	_polap_log0 "Taxonomy geseq of ${_arg_species}"

	# Make sure that mt.1.fa exists at the current folder.
	local _mt1fa="${_arg_outdir}/mt.1.fa"
	if [[ ! -s "${_mt1fa}" ]]; then
		die "ERROR: no such file: ${_mt1fa}"
	fi

	# Search for zip files that match the pattern
	local zip_files=(job-results*.zip)

	# Count the number of matching files
	local file_count=${#zip_files[@]}

	# Check the count and handle accordingly
	if [ $file_count -eq 0 ]; then
		echo "Error: No zip files starting with 'job-results' found."
		exit 1
	elif [ $file_count -gt 1 ]; then
		echo "Error: Multiple zip files found:"
		for file in "${zip_files[@]}"; do
			echo " - $file"
		done
		exit 1
	else
		_polap_log1 "    Only one zip file found - ${zip_files[0]}"
	fi

	mkdir -p "${_geseq}"

	local _geseq_results="${_geseq}/geseq-results"
	rm -rf "${_geseq_results}"
	if _polap_contains_step 1 "${_step_array[@]}"; then
		_polap_log0 "  step 1: extract geseq results"
		_polap_log1 "    output1: ${_geseq_results}/"

		# Get the .zip file name in the {_arg_outdir}
		zip_file=$(find "${_arg_outdir}" -maxdepth 1 -type f -name "job-results*.zip" | head -n 1)

		if [[ -n "$zip_file" ]]; then
			# Extract the filename without the path
			zip_filename=$(basename "$zip_file" .zip)

			_polap_log2 "    zip file found: $zip_filename.zip"
			_polap_log2 "    contents extracted in ${_arg_outdir}"
			_polap_log3_cmd unzip "$zip_file" -d "${_arg_outdir}"
			_polap_log3_cmd mv "${_arg_outdir}/${zip_filename}" "${_geseq_results}"
		else
			die "No .zip file found in ${_arg_outdir}"
		fi
	fi

	local _geseq_gff="${_geseq}/mt.1.gff"
	local _geseq_faa="${_geseq}/mt.1.faa"
	local _geseq_fna="${_geseq}/mt.1.fna"
	local _geseq_filtered_faa="${_geseq}/mt.1.filtered.faa"
	local _geseq_formated_faa="${_geseq}/mt.1.formated.faa"
	if _polap_contains_step 2 "${_step_array[@]}"; then
		_polap_log0 "  step 2: correct the gff3 and copy it to folder: ${_geseq}"
		_polap_log1 "    input1: ${_geseq_results}/"
		_polap_log1 "    output1: ${_geseq_gff}"

		# Define the target folder to copy the files
		local target_folder="${_geseq}"

		# Iterate over each job-results folder
		folder="${_geseq_results}"
		# Check for multi-GFF3.gff3 file
		multi_file=$(find "$folder" -type f -name "*multi-GFF3.gff3" | head -n 1)

		if [[ -n "$multi_file" ]]; then
			# If a multi-GFF3.gff3 file exists, copy it
			_polap_log2 "  found multi-GFF3 file: $multi_file."
			_polap_log2 "  copying it to $target_folder."
			_polap_log3_cmd cp -p "$multi_file" "${_geseq_gff}"
		else
			# Otherwise, look for the general GFF3.gff3 file
			gff3_file=$(find "$folder" -type f -name "*GFF3.gff3" | head -n 1)
			if [[ -n "$gff3_file" ]]; then
				_polap_log2 "    no multi-GFF3 file in $folder."
				_polap_log2 "    copying GFF3 file: $gff3_file to $target_folder."
				_polap_log3_cmd cp -p "$gff3_file" "${_geseq_gff}"
			else
				_polap_log0 "    no GFF3 files found in $folder."
			fi
		fi

		if [[ -s "${_geseq_gff}" ]]; then
			sed -i 's/1\tI/0\tI/' "${_geseq_gff}"
			_polap_log2 "    correct the conda start position: ${_geseq_gff}"
		fi
	fi

	if _polap_contains_step 3 "${_step_array[@]}"; then
		_polap_log0 "  step 3: extract CDS and AA sequences using ${_geseq_gff}"
		_polap_log1 "    input1: ${_mt1fa}"
		_polap_log1 "    input2: ${_geseq_gff}"
		_polap_log1 "    output1: ${_geseq_faa}"
		_polap_log1 "    output2: ${_geseq_fna}"

		_polap_log3_cmd rm -f "${_geseq_faa}" "${_geseq_fna}"
		_polap_log3_cmd gffread -g "${_mt1fa}" -y "${_geseq_faa}" "${_geseq_gff}"
		_polap_log3_cmd gffread -g "${_mt1fa}" -x "${_geseq_fna}" "${_geseq_gff}"
	fi

	if _polap_contains_step 4 "${_step_array[@]}"; then
		_polap_log0 "  step 4: remove sequences with non-stanadrad amino acids: ${_geseq_faa}"
		_polap_log1 "    input1: ${_geseq_faa}"
		_polap_log1 "    output1: ${_geseq_filtered_faa}"

		_polap_log3_cmd rm -f "${_geseq_filtered_faa}"
		_polap_log3_pipe "seqkit grep -s -r -n -v -p '[^ACDEFGHIKLMNPQRSTVWY*]' \
      ${_geseq_faa} >${_geseq_filtered_faa}"
	fi

	# Usage: ./modify_fasta.sh input.fasta output.fasta MG37
	if _polap_contains_step 5 "${_step_array[@]}"; then
		_polap_log0 "  step 5: modify the sequence IDs to conform to NCBI's"
		_polap_log1 "    input1: ${_geseq_filtered_faa}"
		_polap_log1 "    output1: ${_geseq_formated_faa}"

		# Input arguments
		local INPUT_FASTA="${_geseq_filtered_faa}"
		local OUTPUT_FASTA="${_geseq_formated_faa}"
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

		_polap_log1 "    output2: output1 file written to ${_geseq}/${ID1}.faa"

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

		rm -f "${_geseq}/${PREFIX}"*.faa
		cp -p "${OUTPUT_FASTA}" "${_geseq}/${ID1}.faa"
	fi

	if _polap_contains_step 6 "${_step_array[@]}"; then
		_polap_log0 "  step 6: NCBI search and save AA sequences in FASTA format"
		if [[ "${_arg_plastid}" == "off" ]]; then
			_ncbi_query="(mitochondrion[Title] AND complete[Title] AND genome[Title]) AND Viridiplantae[Organism]"
		else
			_ncbi_query="(chloroplast[Title] AND complete[Title] AND genome[Title]) AND Viridiplantae[Organism]"
		fi
		_polap_log0 1. visit https://www.ncbi.nlm.nih.gov/nuccore
		_polap_log0 2. search: ${_ncbi_query}
		_polap_log0 3. use Send to: button to select Coding Sequences in FASTA Protein
		_polap_log0 4. create a flie: sequence.txt
		_polap_log0 5. put sequence.txt in the base directory: e.g., ~/all/polap/taxon1
	fi

	# Disable debugging if previously enabled
	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_polap_taxonomy-orthofinder {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_/Menu_/)"

	source "$script_dir/polap-variables-common.sh"
	local _taxonomy_dir="${_arg_outdir}/taxonomy"
	local _assemble="${_taxonomy_dir}/01-assemble"
	local _sample="${_taxonomy_dir}/02-sample"
	local _geseq="${_taxonomy_dir}/03-geseq"
	local _orthofinder="${_taxonomy_dir}/04-orthofinder"
	local _phylogenies="${_taxonomy_dir}/05-phylogenies"
	local _tree="${_taxonomy_dir}/06-tree"

	local _input_short_reads="${_arg_outdir}/s1.fq"

	# Display help message
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	# main
	#
	# select default steps
	if [[ -z "${_arg_steps_include}" ]]; then
		_arg_steps_include="1-5"
		_arg_steps_exclude=""
	fi
	_include="${_arg_steps_include}"
	_exclude="${_arg_steps_exclude}" # Optional range or list of steps to exclude
	_step_array=()
	_step_array=($(_polap_parse_steps "${_include}" "${_exclude}"))

	# Make sure that mt.1.fa exists at the current folder.
	local _aa_all_fasta="sequence.txt"
	if [[ ! -s "${_aa_all_fasta}" ]]; then
		die "ERROR: no such file: ${_aa_all_fasta}"
	fi

	if [[ "${_arg_redo}" == "on" ]]; then
		_polap_log3_cmd rm -rf "${_orthofinder}"
	fi
	_polap_log3_cmd mkdir -p "${_orthofinder}"

	local _taxonomy_full_tsv="${_sample}/taxonomy_full.tsv"
	if _polap_contains_step 1 "${_step_array[@]}"; then
		_polap_log0 "  step 1: copy the AA sequences to ${_orthofinder}"
		_polap_log1 "    output1: ${_orthofinder}/*.faa"
		_polap_log1 "    output2: ${_orthofinder}/custom.taxa.txt"

		_polap_log3_cmd cp -p "${_geseq}"/AA*.faa "${_orthofinder}"
		echo $(_taxonomy_get_taxonomy_info "${_arg_species}" "${_taxonomy_full_tsv}") >"${_orthofinder}/custom.taxa.txt"
	fi

	local _accessions_ingroup_sampled="${_sample}/accesssions_ingroup_sampled.txt"
	local _accessions_outgroup_sampled="${_sample}/accesssions_outgroup_sampled.txt"
	if _polap_contains_step 2 "${_step_array[@]}"; then
		_polap_log0 "  step 2: extract AA fasta files per species"
		_polap_log1 "    input1: sequence.txt"
		_polap_log1 "    output1: ${_orthofinder}/*.fa"

		while read -r g; do
			seqkit grep -r -p "l/|$g" sequence.txt -o "${_orthofinder}/$g.fa"
		done < <(cut -d',' -f1 ${_accessions_ingroup_sampled} | tail -n +2)
		while read -r g; do
			seqkit grep -r -p "l/|$g" sequence.txt -o "${_orthofinder}/$g.fa"
		done < <(cut -d',' -f1 ${_accessions_outgroup_sampled} | tail -n +2)
	fi

	# filter by the number of amino acid sequences
	if _polap_contains_step 3 "${_step_array[@]}"; then
		_polap_log0 "  step 3: filter out AA fasta files with too few genes."
		_polap_log1 "    input1: ${_orthofinder}/*.fa"
		_polap_log1 "    output1: delete some of ${_orthofinder}/*.fa"

		# Get folder path and threshold length from arguments
		FASTA_FOLDER="${_orthofinder}"
		LENGTH_THRESHOLD="${_arg_taxonomy_min_aa}"

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
				_polap_log2 "    deleting $FILE"
				_polap_log3_cmd rm -f "$FILE"
			fi
		done
	fi

	if _polap_contains_step 4 "${_step_array[@]}"; then
		_polap_log0 "  step 4: execute orthofinder on ${_orthofinder}"
		_polap_log1 "    input1: ${_orthofinder}/*.fa"
		_polap_log1 "    output1: ${_orthofinder}/OrthoFinder"
		orthofinder -f "${_orthofinder}" -n 1
	fi

	local _of_genecount="${_orthofinder}/OrthoFinder/Results_1/Orthogroups/Orthogroups.GeneCount.tsv"
	local _of_genecount_max="${_orthofinder}/genecount_max.txt"
	if _polap_contains_step 5 "${_step_array[@]}"; then
		_polap_log0 "  step 5: determine the top OG numbers for phylogeny on ${_orthofinder}"
		_polap_log1 "    input1: ${_of_genecount}"
		_polap_log1 "    output1: ${_of_genecount_max}"
		Rscript --vanilla "${script_dir}/run-polap-r-orthofinder.R" \
			-t ${_of_genecount} >"${_of_genecount_max}"
	fi

	# Disable debugging if previously enabled
	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_polap_taxonomy-phylogeny {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_/Menu_/)"

	source "$script_dir/polap-variables-common.sh"
	local _taxonomy_dir="${_arg_outdir}/taxonomy"
	local _assemble="${_taxonomy_dir}/01-assemble"
	local _sample="${_taxonomy_dir}/02-sample"
	local _geseq="${_taxonomy_dir}/03-geseq"
	local _orthofinder="${_taxonomy_dir}/04-orthofinder"
	local _phylogenies="${_taxonomy_dir}/05-phylogenies"
	local _tree="${_taxonomy_dir}/06-tree"

	local _input_short_reads="${_arg_outdir}/s1.fq"

	# Display help message
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	# main
	#
	# select default steps
	if [[ -z "${_arg_steps_include}" ]]; then
		_arg_steps_include="1-2"
		_arg_steps_exclude=""
	fi
	_include="${_arg_steps_include}"
	_exclude="${_arg_steps_exclude}" # Optional range or list of steps to exclude
	_step_array=()
	_step_array=($(_polap_parse_steps "${_include}" "${_exclude}"))

	if [[ "${_arg_redo}" == "on" ]]; then
		_polap_log3_cmd rm -rf "${_phylogenies}"
	fi
	mkdir -p "${_phylogenies}"

	local _of_genecount_max="${_orthofinder}/genecount_max.txt"
	local _phylogenies_tree="${_phylogenies}/s.treefile"
	local _phylogenies_tree_taxa="${_phylogenies}/taxa.treefile"
	if _polap_contains_step 1 "${_step_array[@]}"; then
		_polap_log0 "  step 1: align the AA sequences per gene"
		_polap_log1 "    input1: ${_orthofinder}/OrthoFinder"
		_polap_log1 "    output1: ${_orthofinder}/d"
		_polap_log1 "    output2: ${_orthofinder}/i"

		# _polap_log1 "    output1: ${_orthofinder}/*.fa"
		local N="${_arg_bridge_min}" # -x 20
		local N=$(<"${_of_genecount_max}")
		local _of="${_orthofinder}"

		rm -rf "${_of}"/d
		rm -rf "${_of}"/i
		mkdir -p "${_of}"/d
		mkdir -p "${_of}"/i
		for ((i = 0; i <= N; i++)); do
			printf -v j "%07d" $i
			cp "$(find ${_of} -name Orthogroup_Sequences)"/OG"${j}".fa "${_of}"/d
		done

		for i in ${_of}/d/*.fa; do muscle -quiet -in $i -out ${_of}/i/$(basename ${i%.fa}.afa); done

		_polap_log0 "Check: ${_of}/i for alignment"
		_polap_log0 "find . -name Orthogroups.GeneCount.tsv | xargs open"
		_polap_log0 "find . -name Orthogroups.GeneCount.tsv"
	fi

	if _polap_contains_step 2 "${_step_array[@]}"; then
		_polap_log0 "  step 2: tree"
		_polap_log1 "    input1: ${_orthofinder}/i"
		local N="${_arg_bridge_min}" # -x 20
		local N=$(<"${_of_genecount_max}")
		local OG="${_orthofinder}"
		local _of="${_orthofinder}"
		local _phy="${_phylogenies}"

		rm -rf "${_phy}"
		mkdir -p "${_phy}"/d
		mkdir -p "${_phy}"/i
		mkdir -p "${_phy}"/s
		for ((i = 0; i <= N; i++)); do
			printf -v j "%07d" $i
			cp "$(find ${_of} -name Orthogroup_Sequences)"/OG"${j}".fa "${_phy}"/d
		done

		for i in ${_phy}/d/*.fa; do
			sed -i -e "s/\(>[a-zA-Z0-9]\+\)_[a-zA-Z0-9]\+/\1/" $i
		done

		for i in ${_phy}/d/*.fa; do
			# cat $i | perl -pe's/l\|([A-Za-z0-9]+)\.\d+_(\S+)/>$1/' | perl -pe's/>(target)_\S+/>$1/' > ${i%.fa}.fa1;
			cat $i | perl -pe's/^>[a-z]{3}\|([A-Za-z0-9_]+)\.{1}.+/>$1/' >${i%.fa}.fa1
		done

		for i in ${_phy}/d/*.fa1; do seqkit rmdup -n $i -o ${i%.fa1}.faa; done
		for i in ${_phy}/d/*.faa; do muscle -quiet -in $i -out ${i%.faa}.afa; done
		cp ${_phy}/d/*.afa ${_phy}/s
		iqtree -p ${_phy}/s -T AUTO -B 1000
	fi

	# Disable debugging if previously enabled
	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_polap_taxonomy-tree {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_/Menu_/)"

	source "$script_dir/polap-variables-common.sh"
	local _assemble="${_taxonomy_dir}/01-assemble"
	local _sample="${_taxonomy_dir}/02-sample"
	local _geseq="${_taxonomy_dir}/03-geseq"
	local _orthofinder="${_taxonomy_dir}/04-orthofinder"
	local _phylogenies="${_taxonomy_dir}/05-phylogenies"
	local _tree="${_taxonomy_dir}/06-tree"

	# Display help message
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	# main
	#
	# select default steps
	if [[ -z "${_arg_steps_include}" ]]; then
		_arg_steps_include="1-4"
		_arg_steps_exclude=""
	fi
	_include="${_arg_steps_include}"
	_exclude="${_arg_steps_exclude}" # Optional range or list of steps to exclude
	_step_array=()
	_step_array=($(_polap_parse_steps "${_include}" "${_exclude}"))

	mkdir -p "${_tree}"

	local _sequence_ingroup_taxon_tsv="${_sample}/sequence.ingroup.taxon.tsv"
	local _sequence_outgroup_taxon_tsv="${_sample}/sequence.outgroup.taxon.tsv"
	local _sequence_taxon_tsv="${_tree}/sequence.taxon.tsv"
	local _sequence_taxon_map="${_tree}/sequence.taxon.map"

	if _polap_contains_step 1 "${_step_array[@]}"; then
		_polap_log0 "  step 1: "
		_polap_log1 "    input1: ${_sequence_ingroup_taxon_tsv}"
		_polap_log1 "    input2: ${_sequence_outgroup_taxon_tsv}"
		_polap_log1 "    output1: ${_sequence_taxon_tsv}"
		cat "${_sequence_ingroup_taxon_tsv}" "${_sequence_outgroup_taxon_tsv}" \
			>"${_sequence_taxon_tsv}"
	fi

	local _phylogenies_tree="${_phylogenies}/s.treefile"
	local _phylogenies_tree_taxa="${_tree}/s.1.treefile"
	if _polap_contains_step 2 "${_step_array[@]}"; then
		_polap_log0 "  step 2: tree rename tip"
		_polap_log1 "    input1: ${_sequence_taxon_tsv}"
		_polap_log1 "    input2: ${_phylogenies_tree}"
		_polap_log1 "    output1: ${_sequence_taxon_map}"
		_polap_log1 "    output2: ${_phylogenies_tree_taxa}"

		while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7; do
			clean_col1="${col1%.*}"
			# Replace spaces in columns with underscores
			concatenated_column=$(echo "$col3"_"$col4"_"$col5"_"$col6"_"$col7" | tr ' ' '_')

			# Print the transformed line
			echo -e "$clean_col1\t$concatenated_column"
		done <"${_sequence_taxon_tsv}" >"${_sequence_taxon_map}"

		# nw_rename short_IDs.nw id2longname.map
		nw_rename "${_phylogenies_tree}" \
			"${_sequence_taxon_map}" \
			>"${_phylogenies_tree_taxa}"
	fi

	if _polap_contains_step 3 "${_step_array[@]}"; then
		_polap_log0 "  step 3: tree draw gene presence and absence table"

	fi

	if _polap_contains_step 4 "${_step_array[@]}"; then
		_polap_log0 "  step 4: tree draw tree with classification ranks"

	fi

	# phylogeny draw using ggtree kind
	# proteins presence-absence table
	# sequence alignment per OG or protein group
	#

	# Disable debugging if previously enabled
	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_polap_taxonomy {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_/Menu_/)"

	local _taxonomy_dir

	source "$script_dir/polap-variables-common.sh"
	_taxonomy_dir="${_arg_outdir}/taxonomy"

	help_message=$(
		cat <<HEREDOC
Organel genome analysis - assembly, taxon sampling, alignment, and phylogeny

Inputs
------

- long-read data: ${_arg_long_reads} (default: l.fq)
- short-read data 1 : ${_arg_short_read1} (no default)
- short-read data 2 (optional): ${_arg_short_read2} (no default)

--taxonomy-s: the threshold to change alpha
--taxonomy-stop-after: stage1, stage2, stage3

Outputs
-------

- plastid genome assembly: ${_arg_final_assembly}
- trace plots for the features of plastid genome assemblies
  ${_arg_outdir}/taxonomy/0/summary1-ordered.pdf
  ${_arg_outdir}/taxonomy/0/summary1-ordered.txt
  ${_arg_outdir}/taxonomy/x/summary1-ordered.pdf
  ${_arg_outdir}/taxonomy/x/summary1-ordered.txt

Arguments
---------

Stage 1:
-o ${_arg_outdir}
-l ${_arg_long_reads}: a long-read fastq data file
-a ${_arg_short_read1}: a short-read fastq data file 1
-b ${_arg_short_read2} (optional): a short-read fastq data file 2
--taxonomy-compare-to-fasta <FASTA>
--start-index <INDEX>: from <INDEX>
--end-index <INDEX>: to <INDEX> - 1
--steps-include <STPES>: STEPS can be 1,2,3 or 1-15
--steps-exclude <STPES>: STEPS can be 1,2,3 or 1-15
-t ${_arg_threads}: the number of CPU cores
--random-seed <arg>: 5-digit number
--no-contigger for a single final step: use the long-read polished gfa


--stages-include <STAGES> for internal use
--taxonomy-stop-after stage1 for external interface
--taxonomy-stop-after stage2 for external interface

--genomesize <NUMBER>: the number of steps is equal to 1.



Stage 2:
--taxonomy-s: the threshold to change alpha

Test stage: with neither of these two options, one sequence is selected as 
a final genome.
--species <SPECIES>
--taxonomy-compare-to-fasta <FASTA>

Menus
-----

- help: display this help message
- redo: overwrite previous results
- view: show some results
- reset or clean: delete output files 

Usages
------
$(basename "$0") ${_arg_menu[0]} -l ${_arg_long_reads}
$(basename "$0") ${_arg_menu[0]} -l ${_arg_long_reads} -a ${_arg_short_read1}
$(basename "$0") ${_arg_menu[0]} -l ${_arg_long_reads} -a ${_arg_short_read1} -b ${_arg_short_read2}

Examples
--------
1:
$(basename "$0") x-ncbi-fetch-sra --sra SRR7153095
$(basename "$0") x-ncbi-fetch-sra --sra SRR7161123
$(basename "$0") get-mtdna --plastid --species "Eucalyptus pauciflora"
cp o/00-bioproject/2-mtdna.fasta ptdna-Eucalyptus_pauciflora-known.fa

bash ptgaul/ptGAUL.sh -o o-ptgaul -r ptdna-Eucalyptus_pauciflora-known.fa -g 180000 -l long.fastq
cp -pr o-ptgaul/result_3000 o

2:
Stage 1 - long-read and short-read data
$(basename "$0") taxonomy -l SRR7153095.fastq -a SRR7161123_1.fastq -b SRR7161123_2.fastq --taxonomy-i 1 --taxonomy-p 1 --taxonomy-n 10

$(basename "$0") taxonomy archive

$(basename "$0") taxonomy -o o-a --taxonomy-best --taxonomy-compare-to-fasta ptdna-epauciflora.fa --taxonomy-i 1 -l SRR7153095.fastq -a SRR7161123_1.fastq -b SRR7161123_2.fastq 

$(basename "$0") taxonomy view
3:
$(basename "$0") taxonomy view 2
$(basename "$0") taxonomy -o o2 --anotherdir o --taxonomy-s 314588465 --taxonomy-alpha 3.5 -l SRR7153095.fastq -a SRR7161123_1.fastq -b SRR7161123_2.fastq

For ptDNA of Juncus validus
$(basename "$0") get-mtdna -o jvalidus --plastid --species "Juncus validus"
$(basename "$0") taxonomy -o jvalidus -l SRR21976089.fastq -a SRR21976091_1.fastq -b SRR21976091_2.fastq --taxonomy-min-memory 9 -v --taxonomy-a 50mb --taxonomy-n 20
src/polap.sh taxonomy -l SRR21976089.fastq -a SRR21976091_1.fastq -b SRR21976091_2.fastq -o jvalidus --taxonomy-min-memory 9 -v --taxonomy-a 50mb --taxonomy-n 20

Stage 1 - long-read only
$(basename "$0") taxonomy -l SRR7153095.fastq --taxonomy-compare-to-fasta ptdna-epauciflora.fa --taxonomy-min-memory 9

Menu examples:
$(basename "$0") taxonomy view
$(basename "$0") taxonomy view 0
$(basename "$0") taxonomy view x
$(basename "$0") taxonomy view best 0
$(basename "$0") taxonomy view best x
$(basename "$0") taxonomy report
$(basename "$0") taxonomy report coverage <- with known ptDNA
$(basename "$0") taxonomy best 0 46
$(basename "$0") taxonomy best x 2
$(basename "$0") taxonomy archive
$(basename "$0") taxonomy archive polishing <- backup msbwt as well
$(basename "$0") taxonomy polishing

$(basename "$0") taxonomy assemble
$(basename "$0") taxonomy sample
$(basename "$0") taxonomy geseq
$(basename "$0") taxonomy orthofinder
$(basename "$0") taxonomy iqtree
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return 0
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		# Disable debugging if previously enabled
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	# menu: report
	if [[ "${_arg_menu[1]}" == "report" ]]; then

		# Disable debugging if previously enabled
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	# menu: archive
	# archive the taxonomy analysis
	if [[ "${_arg_menu[1]}" == "archive" ]]; then
		_polap_log0 "archive ${_arg_outdir} to ${_arg_archive}"

		# Disable debugging if previously enabled
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	# menu: ptgaul
	# ptgaul
	if [[ "${_arg_menu[1]}" == "ptgaul" ]]; then

		# Disable debugging if previously enabled
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	# menu: example
	if [[ "${_arg_menu[1]}" == "example" ]]; then
		if [[ "${_arg_menu[2]}" == "1" ]]; then
			_polap_log0 "example 1"
		fi

		if [[ "${_arg_menu[2]}" == "2" ]]; then
			_polap_log0 "example 1"
		fi

		# Disable debugging if previously enabled
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	# main
	#
	# select stages
	if [[ -z "${_arg_stages_include}" ]]; then
		_arg_stages_include="0-4"
		_arg_stages_exclude=""
	fi
	if [[ "${_arg_menu[1]}" == "stage1" ]]; then
		_arg_stages_include="0-1"
		_arg_stages_exclude=""
	elif [[ "${_arg_menu[1]}" == "stage2" ]]; then
		_arg_stages_include="2-4"
		_arg_stages_exclude=""
	fi

	_include="${_arg_stages_include}"
	_exclude="${_arg_stages_exclude}" # Optional range or list of steps to exclude
	_stage_array=()
	_stage_array=($(_polap_parse_steps "${_include}" "${_exclude}"))

	_polap_log0 "Start taxonomy"

	if [[ "${_arg_menu[1]}" == "assemble" ]]; then
		_polap_log1 "  input1: long-read: ${_arg_long_reads}"
		[[ -s "${_arg_long_reads}" ]] || return ${_POLAP_ERR_CMD_OPTION_LONGREAD}
		if [[ "${_arg_short_read1_is}" == "off" ]]; then
			_polap_log1 "  input2: no short-read1"
		else
			_polap_log1 "  input2: short-read1: ${_arg_short_read1}"
			# [[ -s "${_arg_short_read1}" ]] || return ${_POLAP_ERR_CMD_OPTION_SHORTREAD}
		fi
		if [[ "${_arg_short_read2_is}" == "off" ]]; then
			_polap_log1 "  input3: no short-read2"
		else
			_polap_log1 "  input3: short-read2: ${_arg_short_read2}"
			# [[ -s "${_arg_short_read2}" ]] || return ${_POLAP_ERR_CMD_OPTION_SHORTREAD}
		fi

		_run_polap_taxonomy-assemble

		# Disable debugging if previously enabled
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	if [[ "${_arg_menu[1]}" == "species" ]]; then

		_run_polap_taxonomy-species

		# Disable debugging if previously enabled
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	if [[ "${_arg_menu[1]}" == "reference" ]]; then

		_run_polap_taxonomy-reference

		# Disable debugging if previously enabled
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	if [[ "${_arg_menu[1]}" == "sample" ]]; then

		_run_polap_taxonomy-sample

		# Disable debugging if previously enabled
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	if [[ "${_arg_menu[1]}" == "geseq" ]]; then

		_run_polap_taxonomy-geseq

		# Disable debugging if previously enabled
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	if [[ "${_arg_menu[1]}" == "orthofinder" ]]; then

		_run_polap_taxonomy-orthofinder

		# Disable debugging if previously enabled
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	if [[ "${_arg_menu[1]}" == "phylogeny" ]]; then

		_run_polap_taxonomy-phylogeny

		# Disable debugging if previously enabled
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	if [[ "${_arg_menu[1]}" == "tree" ]]; then

		_run_polap_taxonomy-tree

		# Disable debugging if previously enabled
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	# Disable debugging if previously enabled
	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_polap_taxon {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "$script_dir/polap-variables-common.sh" # '.' means 'source'

	# Print help message if requested
	help_message=$(
		cat <<HEREDOC
# Taxonomy
#
# Arguments:
#
# Inputs:
#   sequence.fasta.xml : TinySeq XML of NCBI nucleotide database
# Outputs:
#
# See:
# Menu:
# 1-accession Brassicales
# 1-accession [Magnoliopsida]
# 2-extract [sequence.fasta.xml]
# 3-download-db
# 4-get-taxon
# 5-sample [family]
Example: $(basename $0) ${_arg_menu[0]}
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

	# not a good way because of a single file download
	if [[ "${_arg_menu[1]}" == "download" ]]; then
		# Define the file containing the list of accession numbers (one per line)
		ACCESSION_FILE="s70.txt"

		mkdir aa cds

		# Loop through each accession number in the file
		while IFS= read -r ACCESSION; do
			_polap_log0 "Processing accession: $ACCESSION"

			esearch -db nuccore -query "$ACCESSION[ACCN]" </dev/null |
				elink -target protein |
				efetch -format fasta_cds_aa >"aa/${ACCESSION}.faa"

			sleep 5

			esearch -db nuccore -query "$ACCESSION[ACCN]" </dev/null |
				elink -target protein |
				efetch -format fasta_cds_na >"cds/${ACCESSION}.fna"

			sleep 5

		done <"$ACCESSION_FILE"
	fi

	# csvtk cut -f1 s70.csv | csvtk del-header > s70.txt
	# s70.txt
	# https://www.ncbi.nlm.nih.gov/sites/batchentrez
	# progressivemauve --output=s70.xmfa s70.fa
	# extract LCB with a certain value
	# stripSubsetLCBs s70.xmfa s70.xmfa.bbcols s70.core.xmfa 500 70
	# perl ../../../proj/src/xmfa2fasta.pl -f core.xmfa >core.fa
	# raxmlHPC-PTHREADS -T 56 -s core.fa -w "$PWD" -n core_raxml -m GTRCAT -p 789
	# -m ?
	# -p ?
	#
	# add a new mtDNA to the s70 set to draw the tree.
	#
	# open tree database
	# The Open Tree of Life (OTOL) provides a comprehensive and
	# synthesized phylogeny for all life forms, including plants.
	# https://tree.opentreeoflife.org/
	# https://tree.opentreeoflife.org/opentree/argus/ottol@99252/Magnoliopsida
	# Magnoliopsida or angiosperm
	# https://opentreeoflife.github.io/
	#

	if [[ "${_arg_menu[1]}" == "1-accession" ]]; then
		if [[ "${_arg_menu[2]}" == "outfile" ]]; then
			local SPECIES="Magnoliopsida"
		else
			local SPECIES="${_arg_menu[2]}"
		fi
		# NCBI:
		mkdir -p "${_polap_var_ncbi}"
		esearch \
			-db nuccore \
			-query "(mitochondrion[Title] AND complete[Title] AND genome[Title]) AND ${SPECIES}[Organism]" |
			efetch -format acc >"${_polap_var_ncbi_accessions}"
		_polap_log0_head "${_polap_var_ncbi_accessions}"
		_polap_log0 Go to https://www.ncbi.nlm.nih.gov/sites/batchentrez to fetch sequences in batch mode.
		_polap_log0 You could copy the content of the ${_polap_var_ncbi_accessions}
		_polap_log0 and paste it to the query at https://www.ncbi.nlm.nih.gov/nuccore
		_polap_log0 Then, use Send to: button to select TinySeq XML to save sequences and others.
		_polap_log0 Create a flie: sequence.fasta.xml
	fi

	if [[ "${_arg_menu[1]}" == "2-extract" ]]; then

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
			Rscript "${script_dir}/run-polap-r-taxonomy.R" \
				sample \
				-t ${_polap_var_ncbi_sequence_taxon_tsv} \
				-o ${_polap_var_ncbi_sampled_accession}
		else
			Rscript "${script_dir}/run-polap-r-taxonomy.R" \
				sample \
				-f \
				-t ${_polap_var_ncbi_sequence_taxon_tsv} \
				-o ${_polap_var_ncbi_sampled_accession}
		fi
		_polap_log0 file: ${_polap_var_ncbi_sampled_accession}
	fi

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_polap_geseq {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "$script_dir/polap-variables-common.sh" # '.' means 'source'

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
		[ "$DEBUG" -eq 1 ] && set +x
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
			_polap_log0 "correct the codon start position: ${_polap_var_geseq_gff}"
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
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_polap_orthofinder {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "$script_dir/polap-variables-common.sh" # '.' means 'source'

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
		[ "$DEBUG" -eq 1 ] && set +x
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
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_polap_tree {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "$script_dir/polap-variables-common.sh" # '.' means 'source'

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
		[ "$DEBUG" -eq 1 ] && set +x
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
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}
