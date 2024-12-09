planemo tool_init --force \
	--macros \
	--id "polap_total_length_long" \
	--name "Polap total_length_long" \
	--description ": compute the total number of nucleotides of long-read data." \
	--requirement polap@0.3.7.3 \
	--example_command "polap total-length-long -l l.fq" \
	--example_input l.fq \
	--example_output o/long_total_length.txt \
	--test_case \
	--doi 10.1371/journal.pbio.1001241 \
	--help_from_command "polap total_length_long --help"

echo planemo t --galaxy_root ../../../../galaxy --update_test_data
echo planemo s --host 0.0.0.0 --port 8081 --galaxy_root ../../../../galaxy
