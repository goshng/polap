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
[[ -n "${!_POLAP_INCLUDE_}" ]] && return 0
declare "$_POLAP_INCLUDE_=1"
#
################################################################################

################################################################################
# Makes menu commands as empty files.
################################################################################
function _run_polap_make-menus() {
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
	[[ ${_arg_menu[1]} == "help" ]] && _polap_echo0 "${help_message}" && exit $EXIT_SUCCESS

	_polap_log0 "creating empty menu files ..."

	cat "$script_dir"/*.sh |
		grep "^function _run_polap" |
		grep run_polap |
		grep -v run_polap_x |
		sed 's/function _run_polap_//' | sed 's/() {//' |
		parallel touch {}

	_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}

################################################################################
# Makes menu commands as empty files.
# Creates menus prefixed with x.
################################################################################
function _run_polap_make-menus-all() {
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
	[[ ${_arg_menu[1]} == "help" ]] && _polap_echo0 "${help_message}" && exit $EXIT_SUCCESS

	_polap_log0 "creating empty menu files including development versions ..."

	cat "$script_dir"/*.sh |
		grep "^function _run_polap" |
		grep run_polap |
		sed 's/function _run_polap_//' | sed 's/() {//' |
		parallel touch {}

	_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}

################################################################################
# Deletes menu commands.
# Leaves make-menus command.
################################################################################
function _run_polap_clean-menus() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	help_message=$(
		cat <<HEREDOC
# Clean up the menu files.
#
# Arguments: none
# Inputs: none
# Outputs: empty files with filenames of menus.
Example: $(basename $0) ${_arg_menu[0]}
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && _polap_echo0 "${help_message}" && exit $EXIT_SUCCESS

	_polap_log0 "cleaning up the empty menu files ..."

	cat "$script_dir"/*.sh |
		grep "^function _run_polap" |
		grep run_polap |
		sed 's/function _run_polap_//' | sed 's/() {//' |
		parallel rm -f {}

	# Leave the make-menus empty file.
	touch make-menus

	_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}

################################################################################
# Lists menu of POLAP.
# You need to execute make-menus menu if nothing is displayed.
################################################################################
function _run_polap_list() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	help_message=$(
		cat <<HEREDOC
# Lists menu of POLAP.
#
# You need to execute make-menus menu if nothing is displayed.
Example: $(basename "$0") ${_arg_menu[0]} [all|main|assemble1|annotate|assemble2|polish]
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && _polap_echo0 "${help_message}" && exit $EXIT_SUCCESS

	if [[ ${_arg_menu[1]} == "all" ]]; then
		find . -maxdepth 1 -type f -empty -exec basename {} \; |
			sort >&2
	elif [[ ${_arg_menu[1]} == "main" ]]; then
		_polap_log0 assemble1
		_polap_log0 annotate
		_polap_log0 assemble2
		_polap_log0 flye-polishing
		_polap_log0 prepare-polishing
		_polap_log0 polish
	elif [[ ${_arg_menu[1]} == "assemble1" ]]; then
		_polap_log0 reset
		_polap_log0 summary-reads
		_polap_log0 total-length-long
		_polap_log0 find-genome-size
		_polap_log0 reduce-data
		_polap_log0 flye1
	elif [[ ${_arg_menu[1]} == "annotate" ]]; then
		_polap_log0 blast-genome
		_polap_log0 count-gene
	elif [[ ${_arg_menu[1]} == "assemble2" ]]; then
		_polap_log0 select-reads
		_polap_log0 flye2
	elif [[ ${_arg_menu[1]} == "polish" ]]; then
		_polap_log0 flye-polishing
		_polap_log0 prepare-polishing
		_polap_log0 polish
	else
		_polap_log0 "${help_message}"
	fi

	_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}

################################################################################
#
################################################################################
function _run_polap_menu() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	help_message=$(
		cat <<HEREDOC
# Menu
#
Example: $(basename "$0") ${_arg_menu[0]}
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && _polap__polap_log00 "${help_message}" && exit $EXIT_SUCCESS

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

	local _index_bioproject=0

	# Function to display the main menu
	show_main_menu() {
		_polap_log0 "===================="
		_polap_log0 "     Main Menu      "
		_polap_log0 "===================="
		_polap_log0 "1. BioProject"
		_polap_log0 "2. Species"
		_polap_log0 "3. Seeds"
		_polap_log0 "4. Seed Contigs"
		_polap_log0 "5. Annotation"
		_polap_log0 "E. Exit"
		_polap_log0 "--------------------------------------------------"
		_polap_log0 "POLAP Output Folder (${_index_bioproject}): ${ODIR}"
		_polap_log0 "===================="
		_polap_log0_n "Enter your choice [1] (or just press ENTER): "
	}

	# Function to display the system operations submenu
	show_bioproject_menu() {
		_polap_log0 "===================="
		_polap_log0 "  BioProject "
		_polap_log0 "===================="
		_polap_log0 "1. Check Species"
		_polap_log0 "2. Check known mtDNA"
		_polap_log0 "3. Custom species name for checking known mtDNA"
		_polap_log0 "4. Show the known mtDNA"
		_polap_log0 "5. Select seed contigs"
		_polap_log0 "6. View seed contigs"
		_polap_log0 "7. Detail of seed contigs"
		_polap_log0 "ls. List the content of genome assembly"
		_polap_log0 "S. Check the status of genome assembly"
		_polap_log0 "D. Clean OGA"
		_polap_log0 "M. Back to Main Menu"
		_polap_log0 "--------------------------------------------------"
		_polap_log0 "POLAP Output Folder (${_index_bioproject}): ${ODIR}"
		_polap_log0 "===================="
		_polap_log0_n "Enter your choice: "
	}

	# Function to display the species operations submenu
	show_species_menu() {
		_polap_log0 "===================="
		_polap_log0 "  Species"
		_polap_log0 "===================="
		_polap_log0 "1. View species file"
		_polap_log0 "2. Show file permissions (ls -l)"
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
		_polap_log0 "1. Annotation"
		_polap_log0 "2. Annotation plus copy CDF depth and graph"
		_polap_log0 "3. Annotation plus mixture depth and graph"
		_polap_log0 "4. Annotation plus custom depth range"
		_polap_log0 "5. Depth & Annotation plus copy CDF depth and graph"
		_polap_log0 "6. Depth & Annotation plus mixture depth and graph"
		_polap_log0 "7. Depth & Annotation plus custom depth range"
		_polap_log0 "8. View each"
		_polap_log0 "9. View all mt.contig.name"
		_polap_log0 "C. Create mt.contig.name from Bandage"
		_polap_log0 "D. Clean OGA"
		_polap_log0 "M. Back to Main Menu"
		_polap_log0 "--------------------------------------------------"
		_polap_log0 "POLAP Output Folder (${_index_bioproject}): ${ODIR}"
		_polap_log0 "===================="
		_polap_log0_n "Enter your choice: "
	}

	# Function to display the file operations submenu
	show_assemble_menu() {
		_polap_log0 "===================="
		_polap_log0 "  Assembly"
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
			case "$current_menu" in
			main)
				show_main_menu
				read -r choice
				if [[ -z "${choice}" ]]; then
					choice=1
				fi
				_polap_log0_only "${choice}"
				_polap_log0_only "${choice}"
				# Convert the input to lowercase for case-insensitive comparison
				local choice="${choice,,}"
				case $choice in
				1)
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
							current_menu="bioproject"
						fi
					else
						_polap_log0 "No folders found!"
						current_menu="main" # Go back to Main Menu
					fi
					;;

				2) current_menu="species" ;;
				3) current_menu="seeds" ;;
				4) current_menu="contigs" ;;
				5) current_menu="assemble" ;;
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
					src/polap.sh -o "${ODIR}" select-contigs-1-annotation --log-stderr
					current_menu="contigs"
					;;
				6)
					src/polap.sh -o "${ODIR}" select-contigs view --log-stderr
					current_menu="contigs"
					;;
				ls)
					ls -l "${ODIR}" >&3
					;;
				s)
					_polap_log0_file "${_polap_var_wga_contigger_gfa}"
					if [[ -s "${_polap_var_wga_contigger_gfa}" ]]; then
						_polap_log0 "WGA Assembled: ${ODIR}"
					else
						_polap_log0 "WGA not assembled yet: ${ODIR}"
					fi
					;;
				d)
					_polap_log0 "..."
					_polap_log0_cmd rm -rf ${ODIR}/0/?
					_polap_log0_cmd rm -rf ${ODIR}/0/1-custom.depth.range.txt
					_polap_log0_cmd rm -rf ${ODIR}/0/2-custom.depth.range.txt
					_polap_log0_cmd rm -rf ${ODIR}/{1..9}
					_polap_log0_cmd rm -f ${ODIR}/0/mt.contig.name-?
					;;
				m | e) current_menu="main" ;; # Go back to Main Menu
				*) _polap_log0 "Invalid choice, please select a valid option." ;;
				esac
				;;
			species)
				show_species_menu
				read -r choice
				_polap_log0_only "${choice}"
				# Convert the input to lowercase for case-insensitive comparison
				local choice="${choice,,}"
				case $choice in
				1)
					src/polap.sh -o "${ODIR}" get-mtdna file
					;;
				2)
					_polap_log0 "File permissions..."
					ls -l
					;;
				m | e) current_menu="main" ;; # Go back to Main Menu
				*) _polap_log0 "Invalid choice, please select a valid option." ;;
				esac
				;;
			seeds)
				show_seeds_menu
				read -r choice
				_polap_log0_only "${choice}"
				case $choice in
				1)
					_polap_log0 "..."
					src/polap.sh -o "${ODIR}" select-contigs view
					;;
				2) current_menu="contigs" ;;
				m | e) current_menu="main" ;; # Go back to Main Menu
				*) _polap_log0 "Invalid choice, please select a valid option." ;;
				esac
				;;
			contigs)
				show_contigs_menu
				read -r choice
				_polap_log0_only "${choice}"
				case $choice in
				1)
					src/polap.sh -o "${ODIR}" select-contigs-1-annotation --log-stderr
					if [[ -s "${_polap_var_mtcontig_table}" ]]; then
						_polap_log0_cat "${_polap_var_mtcontig_table}"
					else
						_polap_log0 "No such file: ${_polap_var_mtcontig_table}"
					fi
					;;
				2)
					src/polap.sh -o "${ODIR}" select-contigs-2-depth-range --log-stderr
					;;
				3)
					src/polap.sh -o "${ODIR}" select-contigs-3-depth-range --log-stderr
					;;
				4)
					src/polap.sh -o "${ODIR}" select-contigs-4-annotation-depth create --log-stderr
					;;
				5)
					src/polap.sh -o "${ODIR}" sc5 --log-stderr
					;;
				6)
					src/polap.sh -o "${ODIR}" sc6 --log-stderr
					;;
				7)
					src/polap.sh -o "${ODIR}" sc7 create --log-stderr
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
					src/polap.sh -o "${ODIR}" select-contigs view "${choice}" --log-stderr
					;;
				9)
					src/polap.sh -o "${ODIR}" select-contigs view --log-stderr
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
				d)
					_polap_log0 "..."
					_polap_log0_cmd rm -rf ${ODIR}/0/?
					_polap_log0_cmd rm -rf ${ODIR}/0/1-custom.depth.range.txt
					_polap_log0_cmd rm -rf ${ODIR}/0/2-custom.depth.range.txt
					_polap_log0_cmd rm -rf ${ODIR}/{1..9}
					_polap_log0_cmd rm -f ${ODIR}/0/mt.contig.name-?
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

	_polap_log2 "Function end: $(_polap_log0 $FUNCNAME | sed s/_run_polap_//)"\.
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}
