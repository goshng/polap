#!/bin/bash

# A bash command completion source file for polap3 command.
# It is not working yet.

_polap3_completions() {
	local cur prev subcommand
	cur="${COMP_WORDS[COMP_CWORD]}"
	prev="${COMP_WORDS[COMP_CWORD - 1]}"
	subcommand="${COMP_WORDS[1]}"

	# List of subcommands
	local subcommands="annotate archive assemble assemble1 assemble2 assemble-bioproject bioproject-postprocess bioproject-prepare blast-genome blast-mtdna clean-menus cleanup compare-mtdna copy-sra-bioproject count-gene depth-distribution edges-stats find-genome-size flye1 flye2 flye-polishing get-bioproject get-bioproject-sra get-mtdna get-taxonomy-species gm init intra-reads list log make-menus make-menus-all map-reads mauve-mtdna menu plot-mtdna polap-reads polap-rw-base-length polish prepare-polishing prepare-seeds ptgaul-intra-base-length ptgaul-reads reduce-data report-assembly report-table seeds seeds-gene select-mtdna select-mtdna-org select-reads single-intra-base-length source-menus summary-reads template test-reads total-length-long version"

	# Top-level completion for subcommands
	if [[ $COMP_CWORD -eq 1 ]]; then
		COMPREPLY=($(compgen -W "${subcommands}" -- "$cur"))
		return 0
	fi

	# Define options for each subcommand
	case "$subcommand" in
	annotate)
		local opts="-o -i view help redo all mt table pt-table no-depth"
		;;
	init)
		local opts="-o --help view help redo"
		;;
	summary-reads)
		local opts="-o -l -a -b view help redo"
		;;
	total-length-long)
		local opts="-o -l view help redo"
		;;
	find-genome-size)
		local opts="-o -a -b view help redo"
		;;
	reduce-data)
		local opts="-o -l -m --no-reduction-reads --random-seed view help redo split"
		;;
	flye1)
		local opts="-o -t -g --flye-asm-coverage --test view help redo"
		;;
	flye2)
		local opts="-o -i -j view help redo"
		;;
	edges-stats)
		local opts="-o -i view help redo"
		;;
	seeds | map-reads)
		local opts="-o -i -j view help redo"
		;;
	test-reads | select-reads)
		local opts="-o -i -j -w -c view help redo"
		;;
	prepare-polishing)
		local opts="-o -a -b view help redo"
		;;
	polish)
		local opts="-o -p -f view help redo"
		;;
	blast-genome | count-gene)
		local opts="-o -i view help redo"
		;;
	flye-polishing)
		local opts="-o -j view help redo"
		;;
	*)
		local opts="-o view help redo"
		;;
	esac

	# Check if the current word matches an option; otherwise, show folder contents
	if [[ "$cur" == -* ]]; then
		COMPREPLY=($(compgen -W "${opts}" -- "$cur"))
	else
		COMPREPLY=($(compgen -f -- "$cur"))
	fi

	return 0
}

# Register the autocompletion function with the command
complete -F _polap3_completions polap3
