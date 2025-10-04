# Template: safe helper function
# my_helper <required_arg1> <required_arg2> [optional_arg3]
my_helper() {
	# ---- Inputs ----
	local arg1="$1"
	local arg2="$2"
	local opt3="${3:-}" # optional, default empty (or set a default value)

	# ---- Pre-flight checks ----
	if [[ -z "$arg1" ]]; then
		note0 "ERR[my_helper]: arg1 is empty"
		return 1
	fi
	if [[ -z "$arg2" ]]; then
		note0 "ERR[my_helper]: arg2 is empty"
		return 1
	fi
	# If arg1 is a file, ensure it exists and is non-empty
	if [[ ! -s "$arg1" ]]; then
		note0 "ERR[my_helper]: input file not found or empty: '$arg1'"
		return 1
	fi

	# ---- Guard clause (already done) ----
	# Example: if output exists, skip work and return success
	# local out="path/to/output"
	# if [[ -s "$out" ]]; then
	#   note1 "my_helper: reusing existing output: $out"
	#   return 0
	# fi

	# ---- Do the work (external tool) ----
	# Each critical command is inside an if ! ...; then block
	# so we can log and return non-zero explicitly.
	# Example:
	# if ! external_tool --in "$arg1" --out "$out" --flag "$opt3"; then
	#   note0 "ERR[my_helper]: external_tool failed for '$arg1'"
	#   return 1
	# fi

	# ---- Post-check results ----
	# if [[ ! -s "$out" ]]; then
	#   note0 "ERR[my_helper]: expected output missing: $out"
	#   return 1
	# fi

	# ---- Success ----
	return 0
}
