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
if [[ -n "${!_POLAP_INCLUDE_}" ]]; then
	set -u
	return 0
fi
set -u
declare "$_POLAP_INCLUDE_=1"
#
################################################################################

# Function to parse and process steps
_polap_parse_steps() {
	local include="$1"
	local exclude="$2"
	local step_array=()
	local exclude_array=()

	# Parse the include input
	if [[ "$include" =~ ^[0-9]+-[0-9]+$ ]]; then
		# Range method (e.g., 2-4)
		local start=${include%-*}
		local end=${include#*-}
		for ((step = start; step <= end; step++)); do
			step_array+=("$step")
		done
	elif [[ "$include" =~ ^[0-9]+(,[0-9]+)*$ ]]; then
		# Comma-separated list (e.g., 1,3,4)
		IFS=',' read -r -a step_array <<<"$include"
	else
		return 1
	fi

	# Parse the exclude input if provided
	if [[ -n "$exclude" ]]; then
		if [[ "$exclude" =~ ^[0-9]+-[0-9]+$ ]]; then
			# Range method (e.g., 2-3)
			local start=${exclude%-*}
			local end=${exclude#*-}
			for ((step = start; step <= end; step++)); do
				exclude_array+=("$step")
			done
		elif [[ "$exclude" =~ ^[0-9]+(,[0-9]+)*$ ]]; then
			# Comma-separated list (e.g., 3,4)
			IFS=',' read -r -a exclude_array <<<"$exclude"
		else
			return 1
		fi
	fi

	# # Remove excluded steps from the step array
	# for exclude_step in "${exclude_array[@]}"; do
	# 	step_array=("${step_array[@]/$exclude_step/}")
	# done
	#
	# # Filter out empty strings explicitly

	# Remove excluded steps from the step array explicitly
	local filtered_array=()
	for step in "${step_array[@]}"; do
		if ! [[ " ${exclude_array[@]} " =~ " $step " ]]; then
			# if [[ -n "$step" ]]; then
			filtered_array+=("$step")
		fi
	done
	step_array=("${filtered_array[@]}")

	# Return the step_array
	echo "${step_array[@]}"
	return 0
}

# Check if a value exists in an array
_polap_contains_step() {
	local step="$1"
	shift
	local step_array=("$@")
	# local s
	for s in "${step_array[@]}"; do
		if [[ "$s" == "$step" ]]; then
			return 0 # Found the step
		fi
	done
	return 1 # Step not found
}

# Function to control execution steps
function _run_polap_test-steps {
	local include="${_arg_steps_include}"
	local exclude="${_arg_steps_exclude}" # Optional range or list of steps to exclude
	local step_array=()
	local exclude_array=()

	help_message=$(
		cat <<HEREDOC
Test execute steps.

Arguments
---------

--steps-include <STEPS>: e.g., STEPS 1,2 or 3-5
--steps-exclude <STEPS>: e.g., STEPS 2 or 2-3
STEPS format: use a range (e.g., 2-4) or list (e.g., 1,3,4).
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return 0
	if [[ "${_arg_menu[1]}" == "infile" ]]; then
		if [[ -z "${_arg_steps_include}" ]]; then
			_polap_echo0 "${help_message}" && return 0
		fi
	fi

	# Capture steps and check for errors
	# NOTE: local variable first then use it to get the return value from
	# a bash function. If local A=($(func1)), $? or the return of func1 function
	# is always 0.
	local step_array
	include=""
	step_array=($(_polap_parse_steps "$include" "$exclude"))
	local rvalue=$?
	_polap_log0 "rvalue: $rvalue"
	if [[ "$rvalue" -ne 0 ]]; then
		_polap_log0 "ERROR: Error parsing steps."
		return "${_POLAP_ERR_CMD_OPTION_STEPS}"
	fi

	_polap_log0 "Executing steps: ${step_array[@]}"

	# Check for specific steps and execute selectively
	if _polap_contains_step 1 "${step_array[@]}"; then
		_polap_log0 "  step 1: Initializing..."
		# Add commands for step 1
	fi

	if _polap_contains_step 2 "${step_array[@]}"; then
		_polap_log0 "  step 2: Processing data..."
		# Add commands for step 2
	fi

	if _polap_contains_step 3 "${step_array[@]}"; then
		_polap_log0 "  step 3: Validating results..."
		# Add commands for step 3
	fi

	if _polap_contains_step 4 "${step_array[@]}"; then
		_polap_log0 "  step 4: Finalizing..."
		# Add commands for step 4
	fi

	if _polap_contains_step 5 "${step_array[@]}"; then
		_polap_log0 "  step 5: Finalizing..."
		# Add commands for step 4
	fi

	if _polap_contains_step 6 "${step_array[@]}"; then
		_polap_log0 "  step 6: Finalizing..."
		# Add commands for step 4
	fi

	if _polap_contains_step 7 "${step_array[@]}"; then
		_polap_log0 "  step 7: Finalizing..."
		# Add commands for step 4
	fi

	if _polap_contains_step 8 "${step_array[@]}"; then
		_polap_log0 "  step 8: Finalizing..."
		# Add commands for step 4
	fi

	if _polap_contains_step 9 "${step_array[@]}"; then
		_polap_log0 "  step 9: Finalizing..."
		# Add commands for step 4
	fi

	if _polap_contains_step 10 "${step_array[@]}"; then
		_polap_log0 "  step 10: Finalizing..."
		# Add commands for step 4
	fi

	if _polap_contains_step 11 "${step_array[@]}"; then
		_polap_log0 "  step 11: Finalizing..."
		# Add commands for step 4
	fi

	if _polap_contains_step 12 "${step_array[@]}"; then
		_polap_log0 "  step 12: Finalizing..."
		# Add commands for step 4
	fi

	if _polap_contains_step 13 "${step_array[@]}"; then
		_polap_log0 "  step 13: Finalizing..."
		# Add commands for step 4
	fi

	if _polap_contains_step 14 "${step_array[@]}"; then
		_polap_log0 "  step 14: Finalizing..."
		# Add commands for step 4
	fi

	if _polap_contains_step 15 "${step_array[@]}"; then
		_polap_log0 "  step 15: Finalizing..."
		# Add commands for step 4
	fi

	_polap_log0 "Execution complete."
}
