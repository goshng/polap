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

################################################################################
# The idea is as follows. One intriguing yet somewhat challenging aspect is the use of directional reads in genome assembly. This approach could be useful when the results of plant mitochondrial genome assemblies are difficult to interpret.
#
# Let us consider repetitive sequences, specifically direct repeats and inverted repeats. In an assembly graph, multiple sequences are linked by edges. Following a path along these edges corresponds to interpreting that path as the assembled genome. If the assembly graph is simple, identifying this path is straightforward. A plastid genome (chloroplast genome), for instance, may form a simplest circular structure if a single sequence is connected end to end. More commonly, however, there are three sequences connected by four edges, composed of the large subunit, small subunit, and an inverted repeat. In the case of an inverted repeat, if we regard it as two strands (similar to a DNA double helix), it can essentially be treated as a single sequence. In reality, chloroplast genomes often consist of four sequences arranged in a circular form.
#
# However, if an assembled mitochondrial genome is made up of multiple sequences linked in a complex arrangement, tracing a path through the graph—like one would do with a simpler chloroplast genome—may not be easy. One possible approach is to search for the simplest path in the assembly graph, but confirming whether such a path actually exists is not trivial. Instead, we aim to investigate whether we can simplify the genome assembly graph by using that plausible path as seed contigs.
#
# The aforementioned direct and inverted repeats are included in the assembled genome if the assembly path traverses them more than once. However, direct repeats lie in the same orientation along the path, while inverted repeats are in the opposite orientation. If we only follow the completed path in one direction, direct repeats remain repetitive sequences distinguishable by their flanking regions, whereas inverted repeats, despite being repetitive, differ entirely in orientation. This implies that in a unidirectional assembly process, direct repeats may still be indistinguishable, while inverted repeats might be resolved.
#
# The issue is that raw nucleotide sequences are typically generated without preserving directional information, making it rare to find genome assembly tools that handle directional reads. One exception is directional RNA sequencing, which does preserve orientation and thus factors it into transcript assembly. However, because standard genome assembly does not generate directionally specific sequences, such tools have generally not been necessary. For nuclear genomes, incorporating directional information into assembly likely offers little benefit. But for smaller genomes, such as plant mitochondria or possibly bacterial genomes—where assembly can be complicated by repetitive sequences—taking directionality into account could be advantageous in certain cases. In essence, direct repeats would still remain in the graph, but repeats like inverted repeats might be resolved, simplifying the assembly graph overall.
#
# There are two specific challenges: (1) generating directional read data, and (2) developing or adapting a tool capable of performing genome assembly using these directional reads. We generate directional read data by selecting sequences that map in the same orientation to seed contigs. Among available genome assembly tools, we modified Flye to utilize this directional data in the assembly process.
#
# As a proof of concept, we tested our method on a plant mitochondrial genome that was otherwise difficult to interpret, and found that it produced a somewhat simpler assembly graph. Going forward, it seems promising to apply this method to bacterial genome assembly as well.
################################################################################
function _run_polap_dissemble {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "$script_dir/polap-variables-common.sh" # '.' means 'source'

	help_message=$(
		cat <<HEREDOC
Plastid genome assembly using long-read data without reference

Inputs
------

- long-read data: ${_arg_long_reads}
- short-read data: ${_arg_short_read1}

Outputs
-------

- plastid genome assembly
- trace plots

Arguments
---------
-o ${_arg_outdir}
-l ${_arg_long_reads}: a long-read fastq data file
-a ${_arg_short_read1}: a short-read fastq data file

-m ${_arg_min_read_length}: the long-read sequence length threshold
-t ${_arg_threads}: the number of CPU cores
-c ${_arg_coverage}: the coverage option
--random-seed <arg>: 5-digit number

view
----

Usages
------
$(basename "$0") ${_arg_menu[0]} -l ${_arg_long_reads} -a ${_arg_short_read1}
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		if [[ -s "${_polap_var_wga_contigger_edges_gfa}" ]]; then
			_polap_log0_file "${_polap_var_wga_contigger_edges_gfa}"
		else
			_polap_log0 "No such file: ${_polap_var_wga_contigger_edges_gfa}"
		fi
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
		exit $EXIT_SUCCESS
	fi

	_polap_log0 "starting the whole-genome assembly on ${_arg_outdir} ..."
	_polap_log1 "  output1: ${_polap_var_outdir_s1_fq_stats}"
	_polap_log1 "  output1: ${_polap_var_outdir_s2_fq_stats}"
	_polap_log1 "  output2: ${_polap_var_outdir_long_total_length}"
	_polap_log1 "  output3: ${_polap_var_outdir_genome_size}"
	_polap_log1 "  output4: ${_polap_var_outdir_lk_fq_gz}"
	_polap_log1 "  output5: ${_polap_var_outdir_nk_fq_gz}"
	_polap_log1 "  output6: ${_polap_var_wga_contigger_edges_gfa}"
	_polap_log1 "  output7: ${_polap_var_wga_contigger_edges_stats}"

	# Skip flye1 if you want
	if [[ "${_arg_flye}" == "on" ]]; then
		if [[ -d "${_polap_var_wga}" ]]; then
			if [[ "${_arg_redo}" = "on" ]]; then
				_polap_log3_cmd rm -rf "${_polap_var_wga}"
				_polap_log3_cmd mkdir -p "${_polap_var_wga}"
			else
				if confirm "Do you want to do the whole-genome assembly, which will delete ${_polap_var_wga}?"; then
					_polap_log0 "  deleting and creating ${_polap_var_wga} ..."
					_polap_log3_cmd rm -rf "${_polap_var_wga}"
					_polap_log3_cmd mkdir -p "${_polap_var_wga}"
				else
					_polap_log0 "You have cancelled the whole-genome assembly."
					return
				fi
			fi
		else
			_polap_log0 "  creating ${_polap_var_wga} ..."
			_polap_log3_cmd mkdir -p "${_polap_var_wga}"
		fi
	fi

	if [ -s "${_polap_var_outdir_long_total_length}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log2 "  skipping total-length-long ..."
	else
		_run_polap_total-length-long
	fi

	if [ -s "${_polap_var_outdir_genome_size}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log2 "  skipping find-genome-size ..."
	else
		_run_polap_find-genome-size
	fi

	# prepare-polishing early on to delete short-read data files
	if [ -s "${_polap_var_outdir_msbwt}" ]; then
		_polap_log1 "  skipping the preparation of short-read polishing ..."
	else
		if [ -s "${_polap_var_outdir_msbwt_tar_gz}" ]; then
			_polap_log1 "  decompressing ${_polap_var_outdir_msbwt_tar_gz} ... later when we polish it with the short-read data."
			# tar zxf "${_polap_var_outdir_msbwt_tar_gz}"
		else
			_polap_log0 "  Do the preparation of short-read polishing ... early"
			check_file_existence "${_arg_short_read1}"
			_run_polap_prepare-polishing
		fi
	fi

	if [ -s "${_polap_var_outdir_nk_fq_gz}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log2 "  skipping reduce-data ..."
	else
		_run_polap_reduce-data
	fi

	if [ -s "${_polap_var_outdir_lk_fq_stats}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log2 "  skipping summary-reads ..."
	else
		_run_polap_summary-reads
	fi

	if [ -s "${_polap_var_wga_contigger_edges_gfa}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log2 "  skipping flye1 ..."
	else
		if [[ "${_arg_flye}" == "on" ]]; then
			_run_polap_flye1
		else
			_polap_log2 "  skipping flye1 ..."
		fi
	fi

	if [ -s "${_polap_var_ga_contigger_edges_stats}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log2 "  skipping edges-stats ..."
	else
		_run_polap_edges-stats
		if [ $? -eq $RETURN_FAIL ]; then
			_polap_log0 "ERROR: edges-stats step failed."
			return $RETURN_FAIL
		fi
	fi

	if [ -s "${_polap_var_ga_annotation_all}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log2 "  skipping annotation ..."
	else
		_run_polap_annotate
		if [ $? -eq $RETURN_FAIL ]; then
			_polap_log0 "ERROR: annotation step failed."
			return $RETURN_FAIL
		fi
	fi

	_polap_log1 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}
