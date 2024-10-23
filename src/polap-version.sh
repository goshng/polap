local _polap_version=0.3.5

print_version_history() {
	local _message=$(
		cat <<HEREDOC
POLAP - Plant organelle DNA long-read assembly pipeline.
version ${_polap_version}

----------
v0.3.6
----------
- Option --flye-asm-coverage is added so that --coverage is used only for POLAP
- Option --seed is added so that random sampling can be done. Previous versions
- Semi-automatic seed contig selection

----------
v0.2.6
----------
- Bioconda package is available.
HEREDOC
	)

	echo "${_message}"
}
