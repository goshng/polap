planemo tool_init --force \
	--macros \
	--id "polap_reduce_data" \
	--name "Polap reduce_data" \
	--description ": reduce the long-read data for a Flye genome assembly" \
	--requirement polap@0.3.7.3 \
	--example_command "polap reduce_data" \
	--example_input o/short_expected_genome_size.txt \
	--example_input o/long_total_length.txt \
	--example_input l.fq \
	--example_output o/lk.fq.gz \
	--example_output o/nk.fq.gz \
	--test_case \
	--doi 10.1371/journal.pbio.1001241 \
	--help_from_command "polap reduce_data --help"
