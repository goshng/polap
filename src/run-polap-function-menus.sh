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
[[ -n "${!_POLAP_INCLUDE_}" ]] && return 0
set -u
declare "$_POLAP_INCLUDE_=1"
#
################################################################################

################################################################################
# Makes menu commands as empty files.
################################################################################
function _run_polap_make-menus() { # makes menu commands as empty files.
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	help_message=$(
		cat <<HEREDOC
# Makes menu files for easier command completion.
#
# Arguments: none
# Inputs: none
# Outputs: empty files with filenames of menus.
Example: $(basename $0) ${_arg_menu[0]}
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	_polap_log0 "creating empty menu files ..."

	cat "$script_dir"/*.sh |
		grep "^function _run_polap" |
		grep run_polap |
		grep -v run_polap_x |
		sed 's/function _run_polap_//' |
		sed 's/() {//' |
		awk '{print $1}' |
		parallel touch {}

	local _other_menus=(
		"ptgaul-intra-base-length"
		"single-intra-base-length"
		"polap-rw-base-length"
		"ptgaul-reads"
		"intra-reads"
		"polap-reads"
	)
	for _menu in "${_other_menus[@]}"; do
		touch "${_menu}"
	done

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
# Makes menu commands as empty files.
# Creates menus prefixed with x.
################################################################################
function _run_polap_make-menus-all() { # makes all menu commands as empty files including development version
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	help_message=$(
		cat <<HEREDOC
# Makes menu files including development version for easier command completion.
#
# Arguments: none
# Inputs: none
# Outputs: empty files with filenames of menus.
Example: $(basename $0) ${_arg_menu[0]}
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	_polap_log0 "creating empty menu files including development versions ..."

	cat "$script_dir"/*.sh |
		grep "^function _run_polap" |
		grep run_polap |
		sed 's/function _run_polap_//' |
		sed 's/() {//' |
		awk '{print $1}' |
		parallel touch {}

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
# Deletes menu commands.
# Leaves make-menus command.
################################################################################
function _run_polap_clean-menus() { # deletes menu command empty files
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	help_message=$(
		cat <<HEREDOC
# Delete empty menu files that are no longer required.
#
# Arguments: none
# Inputs: none
# Outputs: no empty menu files
Example: $0 ${_arg_menu[0]}
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	_polap_log0 "cleaning up the empty menu files ..."

	cat "$script_dir"/*.sh |
		grep "^function _run_polap" |
		grep run_polap |
		sed 's/function _run_polap_//' |
		sed 's/() {//' |
		awk '{print $1}' |
		parallel rm -f {}

	local _other_menus=(
		"ptgaul-intra-base-length"
		"single-intra-base-length"
		"polap-rw-base-length"
		"ptgaul-reads"
		"intra-reads"
		"polap-reads"
	)
	for _menu in "${_other_menus[@]}"; do
		rm -f "${_menu}"
	done

	# Leave the make-menus empty file.
	touch make-menus

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
# Lists menu of POLAP.
################################################################################
function _run_polap_list() { # List POLAP menus.
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	help_message=$(
		cat <<HEREDOC
# List menus. A POLAP menu is a first argument of the command.
#
# You need to execute make-menus menu if nothing is displayed.
Example: $0 ${_arg_menu[0]} [all|main|assemble|assemble1|annotate|seeds|assemble2|mtdna|polish]
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	case "${_arg_menu[1]}" in
	all)
		cat "$script_dir"/*.sh |
			grep "^function _run_polap" |
			grep run_polap |
			grep -v run_polap_x |
			sed 's/function _run_polap_//' |
			sed 's/() {//' |
			sort >&3
		;;
	dev)
		if [[ "${_arg_menu[2]}" = "outfile" ]]; then
			local _n=0
		else
			local _n="${_arg_menu[2]}"
		fi
		grep -n "^function _run_polap" "$script_dir"/*polap*.sh |
			sed 's/function _run_polap_//' |
			sed 's/() {//' |
			awk -F'/' '{print $NF}' |
			sort >&3
		;;
	x)
		cat "$script_dir"/*.sh |
			grep "^function _run_polap" |
			grep run_polap |
			sed 's/function _run_polap_//' |
			sed 's/() {//' |
			sort >&3
		;;
	main)
		_polap_log0 assemble
		_polap_log0 assemble1
		_polap_log0 annotate
		_polap_log0 seeds
		_polap_log0 assemble2
		_polap_log0 flye-polishing
		_polap_log0 prepare-polishing
		_polap_log0 polish
		;;
	assemble)
		_polap_log0 assemble1
		_polap_log0 annotate
		_polap_log0 seeds
		_polap_log0 assemble2
		;;
	assemble1)
		_polap_log0 "reset"
		_polap_log0 "summary-reads"
		_polap_log0 "total-length-long"
		_polap_log0 "find-genome-size"
		_polap_log0 "reduce-data"
		_polap_log0 "flye1"
		;;
	annotate)
		_polap_log0 "edges-stats"
		_polap_log0 "blast-genome"
		_polap_log0 "count-gene"
		;;
	seeds)
		_polap_log0 "seeds-1-annotation"
		_polap_log0 "seeds-2-depth-range"
		_polap_log0 "seeds-4-annotation-depth"
		_polap_log0 "seeds"
		;;
	assemble2)
		_polap_log0 "select-reads"
		_polap_log0 "flye2"
		;;
	mtdna)
		_polap_log0 "select-mtdna"
		;;
	polish)
		_polap_log0 "flye-polishing"
		_polap_log0 "prepare-polishing"
		_polap_log0 "polish"
		;;
	*)
		_polap_log0 "assemble"
		_polap_log0 "  assemble1"
		_polap_log0 "    reset"
		_polap_log0 "    summary-reads"
		_polap_log0 "    total-length-long"
		_polap_log0 "    find-genome-size"
		_polap_log0 "    reduce-data"
		_polap_log0 "    flye1"
		_polap_log0 "    edges-stats"
		_polap_log0 "  annotate"
		_polap_log0 "    blast-genome"
		_polap_log0 "    count-gene"
		_polap_log0 "  seeds"
		_polap_log0 "    seeds-gene"
		_polap_log0 "  assemble2"
		_polap_log0 "    map-reads"
		_polap_log0 "    test-reads"
		_polap_log0 "    select-reads"
		_polap_log0 "    flye2"
		_polap_log0 "polish"
		_polap_log0 "  flye-polishing"
		_polap_log0 "  prepare-polishing"
		_polap_log0 "mtdna"
		_polap_log0 "  mauve-mtdna"
		;;
	esac

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
#
################################################################################
function _run_polap_menu() { # Interactive menu interface
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "$script_dir/polap-variables-common.sh"
	source "$script_dir/run-polap-function-utilities.sh"

	help_message=$(
		cat <<HEREDOC
# Menu
Example: $(basename "$0") ${_arg_menu[0]}
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap__polap_log00 "${help_message}" && return

	if [[ "${_arg_menu[1]}" != "infile" ]]; then
		_polap_log0 "ODIR: ${_arg_menu[1]}"
		ODIR="${_arg_menu[1]}"
	else
		_polap_log0 "ODIR: ${ODIR}"
	fi

	source "$script_dir/polap-variables-main.sh"
	source "$script_dir/polap-variables-base.sh"
	source "$script_dir/polap-variables-bioproject.sh"
	source "$script_dir/polap-variables-ga.sh"
	source "$script_dir/polap-variables-wga.sh"
	source "$script_dir/polap-variables-oga.sh"
	source "$script_dir/polap-variables-mtcontig.sh"

	local _index_bioproject=1

	# Function to display the main menu
	show_main_menu() {
		_polap_log0 "===================="
		_polap_log0 "     Main Menu      "
		_polap_log0 "===================="
		_polap_log0 "1. BioProject Folder"
		_polap_log0 "2. BioProject-species"
		_polap_log0 "3. Seed Contigs"
		_polap_log0 "R. Run Polap"
		_polap_log0 "T. Test"
		_polap_log0 "E. Exit"
		_polap_log0 "--------------------------------------------------"
		_polap_log0 "POLAP Output Folder (${_index_bioproject}): ${ODIR}"
		_polap_log0 "===================="
		_polap_log0_n "Enter your choice [2] (or just press ENTER): "
	}

	# Function to display the species operations submenu
	show_bioproject_menu() {
		_polap_log0 "===================="
		_polap_log0 "  BioProject"
		_polap_log0 "===================="
		_polap_log0 "1. Create BioProject"
		_polap_log0 "2. Clean-up just BioProject folder"
		_polap_log0 "M. Back to Main Menu"
		_polap_log0 "--------------------------------------------------"
		_polap_log0 "POLAP Output Folder (${_index_bioproject}): ${ODIR}"
		_polap_log0 "===================="
		_polap_log0_n "Enter your choice: "
	}

	# Function to display the system operations submenu
	#
	# BioProject folder has species name appended
	#
	#
	# Download BioProject runinfo (creating folders for all species)
	# Download
	show_species_menu() {
		_polap_log0 "===================="
		_polap_log0 "  BioProject-species"
		_polap_log0 "===================="
		_polap_log0 "1. Check Species"
		_polap_log0 "2. Check known mtDNA"
		_polap_log0 "3. Custom species name for checking known mtDNA"
		_polap_log0 "4. Show the known mtDNA"
		_polap_log0 "5. Select seed contigs"
		_polap_log0 "6. View seed contigs"
		_polap_log0 "7. Detail of seed contigs"
		_polap_log0 "list. List the content of genome assembly"
		_polap_log0 "status. Check the status of genome assembly"
		_polap_log0 "O. Check the status of organelle genome assembly"
		_polap_log0 "C. Clean OGA"
		_polap_log0 "D. Display SRA"
		_polap_log0 "L. Download long-read SRA"
		_polap_log0 "S. Download short-read SRA"
		_polap_log0 "W. Whole-genome assembly"
		_polap_log0 "M. Back to Main Menu"
		_polap_log0 "--------------------------------------------------"
		_polap_log0 "POLAP Output Folder (${_index_bioproject}): ${ODIR}"
		_polap_log0 "===================="
		_polap_log0_n "Enter your choice: "
	}

	# Function to display the file operations submenu
	show_seeds_menu() {
		_polap_log0 "===================="
		_polap_log0 "  Seed Selection"
		_polap_log0 "===================="
		_polap_log0 "1. Display seeds selection"
		_polap_log0 "2. Display seeds for a particular approach"
		_polap_log0 "M. Back to Main Menu"
		_polap_log0 "--------------------------------------------------"
		_polap_log0 "POLAP Output Folder (${_index_bioproject}): ${ODIR}"
		_polap_log0 "===================="
		_polap_log0_n "Enter your choice: "
	}

	# Function to display the file operations submenu
	show_contigs_menu() {
		_polap_log0 "===================="
		_polap_log0 "  Seed Contigs"
		_polap_log0 "===================="
		_polap_log0 "0. Auto"
		_polap_log0 "1. Annotation"
		_polap_log0 "2. Annotation plus copy CDF depth and graph"
		_polap_log0 "3. Annotation plus mixture depth and graph"
		_polap_log0 "4. Annotation plus custom depth range"
		_polap_log0 "5. Depth & Annotation plus copy CDF depth and graph"
		_polap_log0 "6. Depth & Annotation plus mixture depth and graph"
		_polap_log0 "7. Depth & Annotation plus custom depth range"
		_polap_log0 "8. View each"
		_polap_log0 "9. View all mt.contig.name"
		_polap_log0 "A. View all mt contig table"
		_polap_log0 "C. Create mt.contig.name from Bandage"
		_polap_log0 "X. Assemble OGA"
		_polap_log0 "D. Clean OGA"
		_polap_log0 "M. Back to Main Menu"
		_polap_log0 "--------------------------------------------------"
		_polap_log0 "POLAP Output Folder (${_index_bioproject}): ${ODIR}"
		_polap_log0 "===================="
		_polap_log0_n "Enter your choice: "
	}

	# Function to display the file operations submenu
	show_polap_menu() {
		_polap_log0 "===================="
		_polap_log0 "  Run POLAP"
		_polap_log0 "===================="
		_polap_log0 "1. Annotation"
		_polap_log0 "M. Back to Main Menu"
		_polap_log0 "--------------------------------------------------"
		_polap_log0 "POLAP Output Folder (${_index_bioproject}): ${ODIR}"
		_polap_log0 "===================="
		_polap_log0_n "Enter your choice: "
	}

	# Function to display the file operations submenu
	show_test_menu() {
		_polap_log0 "===================="
		_polap_log0 "  POLAP Test"
		_polap_log0 "===================="
		_polap_log0 "1. Annotation"
		_polap_log0 "M. Back to Main Menu"
		_polap_log0 "--------------------------------------------------"
		_polap_log0 "POLAP Output Folder (${_index_bioproject}): ${ODIR}"
		_polap_log0 "===================="
		_polap_log0_n "Enter your choice: "
	}

	# Function to navigate between menus
	navigate_menu() {
		local current_menu="$1"
		while true; do
			date +"%Y-%m-%d %H:%M:%S" >&3
			case "$current_menu" in
			main)
				show_main_menu
				read -r choice
				if [[ -z "${choice}" ]]; then
					choice=2
				fi
				_polap_log0_only "${choice}"
				_polap_log0_only "${choice}"
				# Convert the input to lowercase for case-insensitive comparison
				local choice="${choice,,}"
				case $choice in
				1) current_menu="bioproject" ;;
				2)
					# Get the list of folders (directories only)
					folders=(PRJ*-*/)
					folders=(o "${folders[@]}")

					# Check if there are folders available
					if [ ${#folders[@]} -gt 0 ]; then
						# Display the list of folders with numbers
						_polap_log0 "Please select a folder:"
						for i in "${!folders[@]}"; do
							_polap_log0 "$((i + 1)). ${folders[i]}"
						done

						# Prompt user for input
						_polap_log0_n "Enter the folder number you want to select [0] (or just press ENTER): "
						read -r choice
						if [[ -z "${choice}" ]]; then
							choice=0
						fi
						_polap_log0_only "${choice}"
						local _index_bioproject="${choice}"

						# Validate user input
						if [[ $choice -lt 0 || $choice -gt ${#folders[@]} ]]; then
							_polap_log0 "Invalid selection. Please select a valid number."
							current_menu="main" # Go back to Main Menu
						else
							if [[ "${choice}" -gt 0 ]]; then
								# Get the selected folder
								selected_folder="${folders[$((choice - 1))]}"

								# Do something with the selected folder (e.g., change directory)
								_polap_log0 "You selected: $selected_folder"
								ls -l "$selected_folder" >&2
								_menu_selected_folder="${selected_folder%/}"
								ODIR="${selected_folder%/}"
								source "$script_dir/polap-variables-main.sh"
								source "$script_dir/polap-variables-base.sh"
								source "$script_dir/polap-variables-bioproject.sh"
								source "$script_dir/polap-variables-ga.sh"
								source "$script_dir/polap-variables-wga.sh"
								source "$script_dir/polap-variables-oga.sh"
								source "$script_dir/polap-variables-mtcontig.sh"
								_polap_log0 "Now in BioProject: ${_menu_selected_folder}"
							else
								_polap_log0 "Stay at the current folder: ${ODIR}"
							fi
							current_menu="species"
						fi
					else
						_polap_log0 "No folders found!"
						current_menu="main" # Go back to Main Menu
					fi
					;;

				3) current_menu="contigs" ;;
				4) current_menu="seeds" ;;
				r) current_menu="assemble" ;;
				t) current_menu="test" ;;
				e | q)
					_polap_log0 "Exiting..."
					exit 0
					;;
				*) _polap_log0 "Invalid choice, please select a valid option." ;;
				esac
				;;
			bioproject)
				show_bioproject_menu
				read -r choice
				_polap_log0_only "${choice}"
				# Convert the input to lowercase for case-insensitive comparison
				local choice="${choice,,}"
				case $choice in
				1)
					# Prompt user for input
					_polap_log0_n "Type in BioProject Accession ID (or just press ENTER): "
					read -r choice
					if [[ -z "${choice}" ]]; then
						choice="PRJNA996135"
					fi
					_polap_log0_only "${choice}"
					local _accession_bioproject="${choice}"
					$0 get-bioproject -o ${_accession_bioproject} \
						--bioproject ${_accession_bioproject}
					_polap_log0 "  deleting the bioproject main folder ... leaving other species"
					if [[ "${_accession_bioproject}" != *"*"* && "${_accession_bioproject}" != *"/"* ]]; then
						_polap_log3_cmd rm -rf ${_accession_bioproject}
					else
						_polap_log0 "Dangerous input, contains * or /"
					fi
					;;
				2)
					# Prompt user for input
					_polap_log0_n "Type in BioProject Accession ID (or just press ENTER): "
					read -r choice
					if [[ -z "${choice}" ]]; then
						choice="PRJNA996135"
					fi
					_polap_log0_only "${choice}"
					local _accession_bioproject="${choice}"
					_polap_log0 "  deleting the bioproject main folder ... leaving other species"
					if [[ "${_accession_bioproject}" != *"*"* && "${_accession_bioproject}" != *"/"* ]]; then
						_polap_log3_cmd rm -rf ${_accession_bioproject}
					else
						_polap_log0 "Dangerous input, contains * or /"
					fi
					;;
				m | e) current_menu="main" ;; # Go back to Main Menu
				*) _polap_log0 "Invalid choice, please select a valid option." ;;
				esac
				;;
			species)
				show_species_menu
				read -r choice
				_polap_log0_only "${choice}"
				case $choice in
				1)
					# src/polap.sh -o "${ODIR}" get-mtdna file
					_polap_log0_cat "${_polap_var_bioproject_species}"
					;;
				2)
					_polap_log0 "..."
					src/polap.sh -o "${ODIR}" get-mtdna --log-stderr
					if [ -s "${_polap_var_bioproject_mtdna_fasta2}" ]; then
						seqkit stats "${_polap_var_bioproject_mtdna_fasta2}" >&3
					else
						_polap_log0 "No mtDNA sequence found for the species: $(<${_polap_var_bioproject_species})"
					fi
					;;
				3)
					_polap_log0_n "Enter spceice name you want to search for in NCBI: "
					read -r choice
					_polap_log0_only "${choice}"
					src/polap.sh -o "${ODIR}" get-mtdna --species "${choice}" --log-stderr
					if [ -s "${_polap_var_bioproject_mtdna_fasta2}" ]; then
						seqkit stats "${_polap_var_bioproject_mtdna_fasta2}"
					else
						_polap_log0 "No mtDNA sequence found for the species: ${choice}"
					fi
					;;
				4)
					if [ -s "${_polap_var_bioproject_mtdna_fasta2}" ]; then
						seqkit stats "${_polap_var_bioproject_mtdna_fasta2}" >&2
					else
						_polap_log0 "No mtDNA sequence found for the species: $(<${_polap_var_bioproject_species})"
					fi
					;;
				5)
					src/polap.sh -o "${ODIR}" seeds-1-annotation --log-stderr
					current_menu="contigs"
					;;
				6)
					src/polap.sh -o "${ODIR}" seeds view --log-stderr
					current_menu="contigs"
					;;
				ls)
					ls -l "${ODIR}" >&3
					;;
				status)
					_polap_log0_file "${_polap_var_wga_contigger_gfa}"
					if [[ -s "${_polap_var_wga_contigger_gfa}" ]]; then
						_polap_log0 "WGA Assembled: ${ODIR}"
					else
						_polap_log0 "WGA not assembled yet: ${ODIR}"
					fi
					;;
				o)
					for i in {1..8}; do
						_polap_log0 "--------------------------------------------------------------------------------"
						_polap_log0 "organelle-genome assembly ${i}:"
						INUM="${i}"
						source "$script_dir/polap-variables-ga.sh"
						_polap_log0_file "${_polap_var_contigger_gfa}"
						if [[ -s "${_polap_var_contigger_gfa}" ]]; then
							_polap_log0 "OGA Assembled: ${_polap_var_contigger_gfa}"
						else
							_polap_log0 "OGA not assembled yet!"
						fi
					done
					;;
				c)
					_polap_log0 "..."
					_polap_log0_cmd rm -rf ${_polap_var_wga}/?
					_polap_log0_cmd rm -rf ${_polap_var_wga}/1-custom.depth.range.txt
					_polap_log0_cmd rm -rf ${_polap_var_wga}/2-custom.depth.range.txt
					_polap_log0_cmd rm -rf ${ODIR}/{1..9}
					_polap_log0_cmd rm -f ${_polap_var_wga}/mt.contig.name-?
					;;
				d)
					_polap_log0 "..."
					local SRA_L=$(cut -f1 ${ODIR}/0-bioproject/1-sra-long-read.tsv)
					local SRA_S=$(cut -f1 ${ODIR}/0-bioproject/1-sra-short-read.tsv)
					_polap_log0 "Long-read SRA ID: ${SRA_L}"
					_polap_log0 "Short-read SRA ID: ${SRA_S}"
					;;
				l)
					_polap_log0_n "getting the file size ... "
					local SRA_L=$(cut -f1 ${ODIR}/0-bioproject/1-sra-long-read.tsv)
					_polap_log0_n "long-read SRA ID: ${SRA_L}: "
					_polap_log0 $(get_sra_size ${SRA_L})
					if confirm "Do you want to download SRA: ${SRA_L}"; then
						_polap_log0 "downloading the long-read ..."
						directory_name=$(dirname "$0")
						${directory_name}/run-polap-ncbitools fetch sra ${SRA_L}
						_polap_log3_cmd mv ${SRA_L}.fastq ${ODIR}/
						_polap_log0 "moving ${SRA_L}.fastq to ${ODIR}"
					else
						_polap_log0 "You cancelled the operation."
					fi
					;;
				s)
					_polap_log0_n "getting the file size ... "
					local SRA_S=$(cut -f1 ${ODIR}/0-bioproject/1-sra-short-read.tsv)
					_polap_log0_n "short-read SRA ID: ${SRA_S}: "
					_polap_log0 $(get_sra_size ${SRA_S})
					if confirm "Do you want to download SRA: ${SRA_S}"; then
						_polap_log0 "downloading the short-read ..."
						directory_name=$(dirname "$0")
						${directory_name}/run-polap-ncbitools fetch sra ${SRA_S}
						_polap_log3_cmd mv ${SRA_L}_?.fastq ${ODIR}/
						_polap_log0 "moving " ${SRA_S}_?.fastq "to ${ODIR}"
					else
						_polap_log0 "You cancelled the operation."
					fi
					_polap_log0 "moving" ${SRA_S}_?.fastq "to ${ODIR}"
					;;
				w)
					$0 assemble-bioproject -o ${ODIR} --log-stderr
					;;
				o)
					$0 assemble2 -o ${ODIR} --log-stderr
					;;
				m | e) current_menu="main" ;; # Go back to Main Menu
				*) _polap_log0 "Invalid choice, please select a valid option." ;;
				esac
				;;
			contigs)
				show_contigs_menu
				read -r choice
				_polap_log0_only "${choice}"
				case $choice in
				0)
					_polap_log0_cmd src/polap.sh -o "${ODIR}" seeds auto -j 8 --log-stderr
					;;
				1)
					src/polap.sh -o "${ODIR}" seeds1 --log-stderr
					if [[ -s "${_polap_var_mtcontig_table}" ]]; then
						_polap_log0_cat "${_polap_var_mtcontig_table}"
					else
						_polap_log0 "No such file: ${_polap_var_mtcontig_table}"
					fi
					;;
				2)
					src/polap.sh -o "${ODIR}" seeds2 --log-stderr
					;;
				3)
					src/polap.sh -o "${ODIR}" seeds3 --log-stderr
					;;
				4)
					src/polap.sh -o "${ODIR}" seeds4 create --log-stderr
					;;
				5)
					src/polap.sh -o "${ODIR}" seeds5 --log-stderr
					;;
				6)
					src/polap.sh -o "${ODIR}" seeds6 --log-stderr
					;;
				7)
					src/polap.sh -o "${ODIR}" seeds7 create --log-stderr
					;;
				8)
					# Prompt user for input
					_polap_log0_n "Enter the seed-contig you want to select [1] (or just press ENTER): "
					read -r choice
					if [[ -z "${choice}" ]]; then
						choice=1
					fi
					_polap_log0_only "${choice}"
					local _index_bioproject="${choice}"
					src/polap.sh -o "${ODIR}" seeds view "${choice}" --log-stderr
					;;
				a)
					for i in {1..8}; do
						_polap_log0 "--------------------------------------------------------------------------------"
						_polap_log0 "Seed contig selection ${i}:"
						src/polap.sh -o "${ODIR}" seeds view "${i}" --log-stderr
					done
					;;
				9)
					src/polap.sh -o "${ODIR}" seeds view --log-stderr
					current_menu="contigs"
					;;
				c)
					# Prompt user for input
					_polap_log0_n "Enter the source number [0]: "
					read -r choice
					if [[ -z "${choice}" ]]; then
						choice=0
					fi
					local _source_assembly=$choice
					_polap_log0_only "${choice}"

					_polap_log0_n "Enter the target number [1]: "
					read -r choice
					if [[ -z "${choice}" ]]; then
						choice=1
					fi
					local _target_assembly=$choice
					_polap_log0_only "${choice}"

					_polap_log0_n "Enter the Bandage edge names: "
					read -r choice
					_polap_log0 "  creating ${_source_assembly}/mt.contig.name-${_target_assembly} ..."
					echo "${choice}" | tr -d ' ' | tr ',' '\n' >"${ODIR}/${_source_assembly}/mt.contig.name-${_target_assembly}"
					;;
				x)
					# Prompt user for input
					_polap_log0_n "Enter the seed contigs you use for OGA assembly [8] (or just press ENTER): "
					read -r choice
					if [[ -z "${choice}" ]]; then
						choice=8
					fi
					_polap_log0_only "${choice}"
					local _j_number="${choice}"

					_polap_log0_n "Enter the minimum rwx [3000] (or just press ENTER): "
					read -r choice
					if [[ -z "${choice}" ]]; then
						choice=3000
					fi
					_polap_log0_only "${choice}"
					local _rwx_number="${choice}"
					src/polap.sh assemble2 -j "${_j_number}" --rwx "${_rwx_number}" -o "${ODIR}" --log-stderr
					;;
				d)
					_polap_log0 "..."
					_polap_log0_cmd rm -rf ${_polap_var_wga}/?
					_polap_log0_cmd rm -rf ${_polap_var_wga}/1-custom.depth.range.txt
					_polap_log0_cmd rm -rf ${_polap_var_wga}/2-custom.depth.range.txt
					_polap_log0_cmd rm -rf ${ODIR}/{1..9}
					_polap_log0_cmd rm -f ${_polap_var_wga}/mt.contig.name-?
					;;
				m | e) current_menu="main" ;; # Go back to Main Menu
				*) _polap_log0 "Invalid choice, please select a valid option." ;;
				esac
				;;
			polap)
				show_polap_menu
				read -r choice
				_polap_log0_only "${choice}"
				# Convert the input to lowercase for case-insensitive comparison
				local choice="${choice,,}"
				case $choice in
				1)
					_polap_log0 "Not implemented yet!"
					# src/polap.sh -o "${ODIR}" get-mtdna file
					;;
				2)
					_polap_log0 "File permissions..."
					ls -l
					;;
				m | e) current_menu="main" ;; # Go back to Main Menu
				*) _polap_log0 "Invalid choice, please select a valid option." ;;
				esac
				;;
			test)
				show_test_menu
				read -r choice
				_polap_log0_only "${choice}"
				# Convert the input to lowercase for case-insensitive comparison
				local choice="${choice,,}"
				case $choice in
				1)
					_polap_log0 "Not implemented yet!"
					# src/polap.sh -o "${ODIR}" get-mtdna file
					;;
				2)
					_polap_log0 "File permissions..."
					ls -l
					;;
				m | e) current_menu="main" ;; # Go back to Main Menu
				*) _polap_log0 "Invalid choice, please select a valid option." ;;
				esac
				;;
			*)
				_polap_log0 "Unknown menu. Returning to main menu."
				current_menu="main"
				;;
			esac
			_polap_log0
		done
	}

	# Start at the main menu
	navigate_menu "main"

	_polap_log3 "Function end: $(_polap_log0 $FUNCNAME | sed s/_run_polap_//)"\.
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}
