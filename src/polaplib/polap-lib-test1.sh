#!/usr/bin/env bash
# Version: v0.1.0
# Library file sourced by polap-bash-test1.sh; defines nested functions.

# Expect the logger to be sourced already by the caller.

polap_lib_test1_entry() {
	local _mode="${1:-bash}"
	polap_log_info "lib: entry -> level2 (mode=${_mode})"
	polap_lib_level2 "${_mode}"
	polap_log_info "lib: level2 returned"
}

polap_lib_level2() {
	local _mode="${1:-bash}"
	polap_log_info "lib: in level2 -> level3 (mode=${_mode})"
	polap_lib_level3 "${_mode}"
	polap_log_info "lib: level3 returned"
}

polap_lib_level3() {
	local _mode="${1:-bash}"
	polap_log_info "lib: level3 dispatch -> ${_mode}"

	if [[ "${_mode}" == "bash" ]]; then
		polap_log_err "Simulated Bash failure via 'false'"
		false # triggers ERR; stack should be printed by trap
	elif [[ "${_mode}" == "python" ]]; then
		polap_log_info "RUN: python3 ${_POLAPLIB_DIR}/polap-py-test1.py"
		python3 "${_POLAPLIB_DIR}/polap-py-test1.py"
		polap_log_info "Python returned (unexpected on error path)"
	elif [[ "${_mode}" == "r" ]]; then
		polap_log_info "RUN: Rscript --vanilla ${_POLAPLIB_DIR}/polap-r-test1.R"
		Rscript --vanilla "${_POLAPLIB_DIR}/polap-r-test1.R"
		polap_log_info "R returned (unexpected on error path)"
	else
		polap_log_err "Unknown mode: ${_mode}"
		return 2
	fi
}
