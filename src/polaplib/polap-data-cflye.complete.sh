# A bash command completion source file for polap-data-cflye.
# It is not working yet.

# polap-data-cflye.complete.sh
# Bash completion for polap-data-cflye
_polap_data_cflye_completions() {
	local cur prev first
	cur="${COMP_WORDS[COMP_CWORD]}"
	prev="${COMP_WORDS[COMP_CWORD - 1]}"
	first="${COMP_WORDS[1]}"

	local top_opts="-h --help -y -c -t -m -e --version"
	local commands="install setup update list search run remove uninstall build assemble download mkdir benchmark clean delete rm get archive man help print-help-all recover"

	local install_opts="conda polap fmlrc getorganelle pmat tippo oatk man cflye dflye all"
	local run_opts="estimate-genomesize getorganelle msbwt ptgaul extract-ptdna-ptgaul polish-ptdna-ptgaul nextdenovo-polish pmat tippo oatk polap-disassemble polap-disassemble-check polap-disassemble-compare"
	local setup_opts="conda polap pmat csv completions"
	local update_opts="github local"
	local man_opts="table-benchmark figure-sheet figure-alpha figure-delta pdf"
	local download_opts="species ptdna man test-data test-data-cflye"
	local get_opts="aflye cflye dflye taxon"

	# Sub-subcommand argument: suggest *.fastq for getorganelle
	if [[ "$first" == "run" && "$prev" == "getorganelle" ]]; then
		COMPREPLY=($(compgen -f -X '!*.fastq' -- "$cur"))
		return 0
	fi

	# download species -> show folders
	if [[ "$first" == "download" && "$prev" == "species" ]]; then
		COMPREPLY=($(compgen -d -- "$cur"))
		return 0
	fi

	if [[ "$first" == "get" && "$prev" == "cflye" ]]; then
		COMPREPLY=($(compgen -d -- "$cur"))
		return 0
	fi

	if [[ "$first" == "run" && "$prev" == "polap-disassemble" ]]; then
		local species_file="./species.tsv"
		if [[ -f "$species_file" ]]; then
			local species_list
			mapfile -t species_list < <(cut -f1 "$species_file")
			COMPREPLY=($(compgen -W "${species_list[*]}" -- "$cur"))
		else
			COMPREPLY=($(compgen -d -- "$cur"))
		fi
		return 0
	fi

	if [[ "$first" == "benchmark" || "$first" == "archive" || "$first" == "recover" ]]; then
		local species_file="./species.tsv"
		if [[ -f "$species_file" ]]; then
			local species_list
			mapfile -t species_list < <(cut -f1 "$species_file")
			COMPREPLY=($(compgen -W "${species_list[*]}" -- "$cur"))
		else
			COMPREPLY=($(compgen -d -- "$cur"))
		fi
		return 0
	fi

	# Command-specific completions
	case "$first" in
	install)
		COMPREPLY=($(compgen -W "$install_opts" -- "$cur"))
		return 0
		;;
	run)
		COMPREPLY=($(compgen -W "$run_opts" -- "$cur"))
		return 0
		;;
	setup)
		COMPREPLY=($(compgen -W "$setup_opts" -- "$cur"))
		return 0
		;;
	update)
		COMPREPLY=($(compgen -W "$update_opts" -- "$cur"))
		return 0
		;;
	man)
		COMPREPLY=($(compgen -W "$man_opts" -- "$cur"))
		return 0
		;;
	download | mkdir)
		COMPREPLY=($(compgen -W "$download_opts" -- "$cur"))
		return 0
		;;
	get)
		COMPREPLY=($(compgen -W "$get_opts" -- "$cur"))
		return 0
		;;
	*)
		if [[ $COMP_CWORD -eq 1 ]]; then
			COMPREPLY=($(compgen -W "$top_opts $commands" -- "$cur"))
			compopt +o default 2>/dev/null
		elif [[ "$cur" == -* ]]; then
			COMPREPLY=($(compgen -W "$top_opts" -- "$cur"))
		fi
		return 0
		;;
	esac
}

alias p2='bash polap/src/polap-data-v2.sh'
complete -F _polap_data_cflye_completions p2
complete -F _polap_data_cflye_completions polap-data-cflye
