planemo tool_init --force \
	--macros \
	--id "polap_test-reads" \
	--name "Polap test-reads" \
	--description ": test reads on seed contigs from a Flye genome assembly" \
	--requirement polap@0.3.7.3 \
	--example_command "polap test-reads" \
	--example_input o/lk.fq.gz \
	--example_input o/0/mt.contig.name-1 \
	--example_input o/0/30-contigger/graph_final.fasta \
	--example_output o/1/06-summary/ptgaul-reads/0.gfa \
	--test_case \
	--doi 10.1371/journal.pbio.1001241 \
	--help_from_command "polap test-reads --help"
