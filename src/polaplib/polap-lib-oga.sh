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
# Convert numbers between different units.
################################################################################

################################################################################
# Ensure that the current script is sourced only once
source "${_POLAPLIB_DIR}/run-polap-function-include.sh"
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

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
	echo "[ERROR] This script must be sourced, not executed: use 'source $BASH_SOURCE'" >&2
	return 1 2>/dev/null || exit 1
fi
: "${_POLAP_DEBUG:=0}"
: "${_POLAP_RELEASE:=0}"

################################################################################
# Function to convert base pairs to the highest appropriate unit
# Example usage
# bp=31846726397
# convert_bp $bp
################################################################################
function _polap_lib_oga-estimate-read-sampling-rate {
	# polap-cmd-oga: test-reads
	local _pread_sel="ptgaul-reads"
	local _read_names="ptgaul"

	_polap_utility_get_contig_length \
		"${_polap_var_oga_contig}/contig.fa" \
		"${_polap_var_oga_contig}/contig_total_length.txt"
	local CONTIG_LENGTH=$(<"${_polap_var_oga_contig}/contig_total_length.txt")
	_polap_log2_cat "${_polap_var_oga_contig}/contig_total_length.txt"

	_polap_log1 "  determines which long-read data to use ..."
	local _source_long_reads_fq=""
	_polap_oga_determine-long-read-file _source_long_reads_fq

	local size_source_lfq="${_polap_var_oga_contig}"/l.txt
	if [[ ! -s "${size_source_lfq}" ]]; then
		_polap_lib_fastq-total-length-of "${_source_long_reads_fq}" "${size_source_lfq}"
	fi
	local len_source_lfq=$(<"${size_source_lfq}")

	local subsample_lfq="${_polap_var_oga_contig}"/l.subsample.fq
	if [[ ! -s "${subsample_lfq}" ]]; then
		_polap_lib_fastq-sample-to "${_source_long_reads_fq}" "${subsample_lfq}" 1g
	fi

	local size_subsample_lfq="${_polap_var_oga_contig}"/l.subsample.txt
	if [[ ! -s "${size_subsample_lfq}" ]]; then
		_polap_lib_fastq-total-length-of "${subsample_lfq}" "${size_subsample_lfq}"
	fi
	local len_subsample_lfq=$(<"${size_subsample_lfq}")

	local rate_lfq=$(echo "scale=9; ${len_subsample_lfq}/$len_source_lfq" | bc)

	_polap_log2 "  creating folders for read selection type: ${_pread_sel}"
	_polap_log3_cmd mkdir -p "${_polap_var_oga_reads}/${_pread_sel}"
	_polap_log3_cmd mkdir -p "${_polap_var_oga_seeds}/${_pread_sel}"
	_polap_log3_cmd mkdir -p "${_polap_var_oga_subsample}/${_pread_sel}"

	local i=0
	for ((i = 0; i < 5; i++)); do
		_polap_log3_cmd mkdir -p "${_polap_var_oga_reads}/${_pread_sel}/${i}"

		_polap_lib_fastq-sample "${_source_long_reads_fq}" "${subsample_lfq}" "${rate_lfq}"

		_polap_log1 "  mapping long-read data on the seed contigs using minimap2 ..."
		_polap_log2 "    input1: ${_polap_var_oga_contig}/contig.fa"
		_polap_log2 "    input2: ${_source_long_reads_fq}"
		_polap_log2 "    input2: ${subsample_lfq}"
		_polap_log2 "    output: ${_polap_var_oga_contig}/contig.paf"
		# -p 0.5 -N 1 --secondary=no
		# -p 0.5 -N 1 \
		_polap_log3_pipe "minimap2 -cx \
      ${_arg_minimap2_data_type} \
      ${_polap_var_oga_contig}/contig.fa \
      ${subsample_lfq} \
      -t ${_arg_threads} \
      -o ${_polap_var_oga_contig}/contig.paf \
      >${_polap_output_dest} 2>&1"

		_polap_log1 "  converting PAF to TAB ..."
		_polap_log2 "    input1: ${_polap_var_oga_contig}/contig.paf"
		_polap_log2 "    output: ${_polap_var_oga_contig}/contig.tab"
		_polap_log3_cmd bash "${_POLAPLIB_DIR}/polap-bash-minimap2-paf2tab.sh" "${_arg_min_read_length}" "${_polap_var_oga_contig}/contig.paf" "${_polap_var_oga_contig}/contig.tab"

		_polap_log3_pipe "Rscript --vanilla ${_POLAPLIB_DIR}/polap-r-pairs.R \
	    -m ${_polap_var_mtcontigname} \
		  -t ${_polap_var_oga_contig}/contig.tab \
		  --out ${_polap_var_oga_reads}/${_pread_sel}/${i} \
      -w ${_arg_single_min} \
		  -r ${_arg_pair_min} \
		  -x ${_arg_bridge_min} \
      --all \
		  >${_polap_output_dest} 2>&1"

		_polap_log1 "  selecting long reads for ${_pread_sel}"
		_polap_log2 "    input1: ${_source_long_reads_fq}"
		_polap_log2 "    input2: ${_polap_var_oga_reads}/${_pread_sel}/${i}/${_read_names}.names"
		_polap_log2 "    output: ${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz"

		_polap_log3_pipe "seqtk subseq \
		    ${_source_long_reads_fq} \
		    ${_polap_var_oga_reads}/${_pread_sel}/${i}/${_read_names}.names |\
		    gzip >${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz"

		_polap_log1 "  sampling reads ..."
		_polap_log2 "    input1: ${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz"
		local _contig_length_bp=$(_polap_utility_convert_bp ${CONTIG_LENGTH})
		_polap_log2 "    input2 (seed contig size): ${_contig_length_bp}"
		_polap_utility_get_contig_length \
			"${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz" \
			"${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.len"
		local _seeds_length=$(<"${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.len")
		_polap_log2_cat "${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.len"
		local _seeds_length_bp=$(_polap_utility_convert_bp ${_seeds_length})
		_polap_log2 "    result1 (total size of reads mapped contigs): ${_seeds_length_bp}"
		local _expected_organelle_coverage=$((_seeds_length / CONTIG_LENGTH))
		_polap_log2 "    result2 (expected organelle coverage): ${_expected_organelle_coverage}x"

		# sampling rate
		local _rate=$(echo "scale=9; ${_arg_coverage_oga}/$_expected_organelle_coverage" | bc)

		_polap_log0 "rate_lfq (<1.0): ${rate_lfq}"
		_polap_log0 "rate (0.1 ~ 0.5): ${_rate}"

		if (($(echo "${rate_lfq} > 1.0" | bc -l))); then
			echo "1.0" >"${_polap_var_oga_contig}/rate_lfq.txt"
			echo ${i} >"${_polap_var_oga_contig}/index.txt"
			return
		fi

		if (($(echo "$_rate < 0.1" | bc -l))); then
			rate_lfq=$(echo "scale=9; ${rate_lfq}/2" | bc)
		elif (($(echo "$_rate > 0.5" | bc -l))); then
			rate_lfq=$(echo "scale=9; ${rate_lfq}*2" | bc)
		else
			echo "${rate_lfq}" >"${_polap_var_oga_contig}/rate_lfq.txt"
			echo ${i} >"${_polap_var_oga_contig}/index.txt"
			return
		fi

	done

	echa 4 >"${_polap_var_oga_contig}/index.txt"
}

# This function is used to run the assembly rate estimation and read selection
# for organelle genomes, specifically for the removal of plastid DNA.
# It estimates the sampling rate of long reads mapped to the seed contig and
# selects reads based on that rate.
#
function _polap_lib_oga-estimate-read-sampling-rate-remove-ptdna {
	# polap-cmd-oga: test-reads
	local _pread_sel="ptgaul-reads"
	local _read_names="ptgaul"

	mkdir -p "${_polap_var_oga_contig}"
	_polap_log1 "  extracts contig sequeces from the assembly: ${_polap_var_ga_contigger_edges_fasta}"
	_polap_log2 "    input1: ${_polap_var_ga_contigger_edges_fasta}"
	_polap_log2 "    input2: ${_polap_var_mtcontigname}"
	_polap_log2 "    output: ${_polap_var_oga_contig}/contig.fa"
	_polap_log3_pipe "seqkit grep \
    -f ${_polap_var_mtcontigname} \
		${_polap_var_ga_contigger_edges_fasta} \
		-o ${_polap_var_oga_contig}/contig.fa \
		2>${_polap_output_dest}"

	_polap_utility_get_contig_length \
		"${_polap_var_oga_contig}/contig.fa" \
		"${_polap_var_oga_contig}/contig_total_length.txt"
	local CONTIG_LENGTH=$(<"${_polap_var_oga_contig}/contig_total_length.txt")
	_polap_log2_cat "${_polap_var_oga_contig}/contig_total_length.txt"

	_polap_log1 "  determines which long-read data to use ..."
	local _source_long_reads_fq=""
	_polap_oga_determine-long-read-file _source_long_reads_fq
	# _source_long_reads_fq="${annotatedir}"/mt.fq

	#############################################################################
	# remove ptDNA-origin long reads
	#
	# local ptdna="${_arg_unpolished_fasta}"
	# _polap_log3_cmdout minimap2 -x \
	# 	${_arg_minimap2_data_type} \
	# 	${ptdna} \
	# 	${_source_long_reads_fq} \
	# 	-t ${_arg_threads} \
	# 	-o ${_polap_var_oga_contig}/ptdna.paf
	#
	# _polap_log3_cmd bash "${_POLAPLIB_DIR}/polap-bash-minimap2-paf2tab.sh" \
	# 	500 \
	# 	"${_polap_var_oga_contig}/ptdna.paf" \
	# 	"${_polap_var_oga_contig}/ptdna.tab"
	#
	# _polap_log3_cmdout Rscript --vanilla \
	# 	"${_POLAPLIB_DIR}/polap-r-select-read-using-paf.R" \
	# 	--min-length 500 \
	# 	--min-identity "${_arg_remove_plastid_min_identity}" \
	# 	${_polap_var_oga_contig}/ptdna.tab \
	# 	${_polap_var_oga_contig}/ptdna.txt
	#
	# seqkit grep -v \
	# 	-f ${_polap_var_oga_contig}/ptdna.txt \
	# 	${_source_long_reads_fq} \
	# 	>"${_source_long_reads_fq}".tmp
	#
	# cp -p "${_source_long_reads_fq}" "${_source_long_reads_fq}".old
	# mv "${_source_long_reads_fq}".tmp "${_source_long_reads_fq}"

	local size_source_lfq="${_polap_var_oga_contig}"/l.txt
	if [[ ! -s "${size_source_lfq}" ]]; then
		_polap_lib_fastq-total-length-of "${_source_long_reads_fq}" "${size_source_lfq}"
	fi
	local len_source_lfq=$(<"${size_source_lfq}")

	local subsample_lfq="${_polap_var_oga_contig}"/l.subsample.fq
	if [[ ! -s "${subsample_lfq}" ]]; then
		_polap_lib_fastq-sample-to "${_source_long_reads_fq}" "${subsample_lfq}" 1g
	fi

	local size_subsample_lfq="${_polap_var_oga_contig}"/l.subsample.txt
	if [[ ! -s "${size_subsample_lfq}" ]]; then
		_polap_lib_fastq-total-length-of "${subsample_lfq}" "${size_subsample_lfq}"
	fi
	local len_subsample_lfq=$(<"${size_subsample_lfq}")

	local rate_lfq=$(echo "scale=9; ${len_subsample_lfq}/$len_source_lfq" | bc)

	_polap_log2 "  creating folders for read selection type: ${_pread_sel}"
	_polap_log3_cmd mkdir -p "${_polap_var_oga_reads}/${_pread_sel}"
	_polap_log3_cmd mkdir -p "${_polap_var_oga_seeds}/${_pread_sel}"
	_polap_log3_cmd mkdir -p "${_polap_var_oga_subsample}/${_pread_sel}"

	local i=0
	for ((i = 0; i < 5; i++)); do
		_polap_log3_cmd mkdir -p "${_polap_var_oga_reads}/${_pread_sel}/${i}"

		_polap_lib_fastq-sample "${_source_long_reads_fq}" "${subsample_lfq}" "${rate_lfq}"

		_polap_log1 "  mapping long-read data on the seed contigs using minimap2 ..."
		_polap_log2 "    input1: ${_polap_var_oga_contig}/contig.fa"
		_polap_log2 "    input2: ${_source_long_reads_fq}"
		_polap_log2 "    input2: ${subsample_lfq}"
		_polap_log2 "    output: ${_polap_var_oga_contig}/contig.paf"
		# -p 0.5 -N 1 --secondary=no
		# -p 0.5 -N 1 \
		_polap_log3_pipe "minimap2 -cx \
      ${_arg_minimap2_data_type} \
      ${_polap_var_oga_contig}/contig.fa \
      ${subsample_lfq} \
      -t ${_arg_threads} \
      -o ${_polap_var_oga_contig}/contig.paf \
      >${_polap_output_dest} 2>&1"

		_polap_log1 "  converting PAF to TAB ..."
		_polap_log2 "    input1: ${_polap_var_oga_contig}/contig.paf"
		_polap_log2 "    output: ${_polap_var_oga_contig}/contig.tab"
		_polap_log3_cmd bash "${_POLAPLIB_DIR}/polap-bash-minimap2-paf2tab.sh" "${_arg_min_read_length}" "${_polap_var_oga_contig}/contig.paf" "${_polap_var_oga_contig}/contig.tab"

		_polap_log3_pipe "Rscript --vanilla ${_POLAPLIB_DIR}/polap-r-pairs.R \
	    -m ${_polap_var_mtcontigname} \
		  -t ${_polap_var_oga_contig}/contig.tab \
		  --out ${_polap_var_oga_reads}/${_pread_sel}/${i} \
      -w ${_arg_single_min} \
		  -r ${_arg_pair_min} \
		  -x ${_arg_bridge_min} \
      --all \
		  >${_polap_output_dest} 2>&1"

		_polap_log1 "  selecting long reads for ${_pread_sel}"
		_polap_log2 "    input1: ${_source_long_reads_fq}"
		_polap_log2 "    input2: ${_polap_var_oga_reads}/${_pread_sel}/${i}/${_read_names}.names"
		_polap_log2 "    output: ${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz"

		_polap_log3_pipe "seqtk subseq \
		    ${_source_long_reads_fq} \
		    ${_polap_var_oga_reads}/${_pread_sel}/${i}/${_read_names}.names |\
		    gzip >${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz"

		_polap_log1 "  sampling reads ..."
		_polap_log2 "    input1: ${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz"
		local _contig_length_bp=$(_polap_utility_convert_bp ${CONTIG_LENGTH})
		_polap_log2 "    input2 (seed contig size): ${_contig_length_bp}"
		_polap_utility_get_contig_length \
			"${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz" \
			"${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.len"
		local _seeds_length=$(<"${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.len")
		_polap_log2_cat "${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.len"
		local _seeds_length_bp=$(_polap_utility_convert_bp ${_seeds_length})
		_polap_log2 "    result1 (total size of reads mapped contigs): ${_seeds_length_bp}"
		local _expected_organelle_coverage=$((_seeds_length / CONTIG_LENGTH))
		_polap_log2 "    result2 (expected organelle coverage): ${_expected_organelle_coverage}x"

		# sampling rate
		local _rate=$(echo "scale=9; ${_arg_coverage_oga}/$_expected_organelle_coverage" | bc)

		_polap_log0 "rate_lfq (<1.0): ${rate_lfq}"
		_polap_log0 "rate (0.1 ~ 0.5): ${_rate}"

		if (($(echo "${rate_lfq} > 1.0" | bc -l))); then
			echo "1.0" >"${_polap_var_oga_contig}/rate_lfq.txt"
			echo ${i} >"${_polap_var_oga_contig}/index.txt"
			return
		fi

		if (($(echo "$_rate < 0.1" | bc -l))); then
			rate_lfq=$(echo "scale=9; ${rate_lfq}/2" | bc)
		elif (($(echo "$_rate > 0.5" | bc -l))); then
			rate_lfq=$(echo "scale=9; ${rate_lfq}*2" | bc)
		else
			echo "${rate_lfq}" >"${_polap_var_oga_contig}/rate_lfq.txt"
			echo ${i} >"${_polap_var_oga_contig}/index.txt"
			return
		fi

	done

	echa 4 >"${_polap_var_oga_contig}/index.txt"
}

function _polap_lib_oga-select-reads-with-sampling-rate {
	# polap-cmd-oga: test-reads
	local _pread_sel="ptgaul-reads"
	local _read_names="ptgaul"

	_polap_utility_get_contig_length \
		"${_polap_var_oga_contig}/contig.fa" \
		"${_polap_var_oga_contig}/contig_total_length.txt"
	local CONTIG_LENGTH=$(<"${_polap_var_oga_contig}/contig_total_length.txt")
	_polap_log2_cat "${_polap_var_oga_contig}/contig_total_length.txt"

	_polap_log1 "  determines which long-read data to use ..."
	local subsample_lfq="${_polap_var_oga_contig}"/l.subsample.fq
	local _source_long_reads_fq="${subsample_lfq}"

	# local rate_lfq=$(echo "scale=9; ${len_subsample_lfq}/$len_source_lfq" | bc)

	_polap_log2 "  creating folders for read selection type: ${_pread_sel}"
	_polap_log3_cmd mkdir -p "${_polap_var_oga_reads}/${_pread_sel}"
	_polap_log3_cmd mkdir -p "${_polap_var_oga_seeds}/${_pread_sel}"
	_polap_log3_cmd mkdir -p "${_polap_var_oga_subsample}/${_pread_sel}"

	local i=0
	{
		_polap_log3_cmd mkdir -p "${_polap_var_oga_reads}/${_pread_sel}/${i}"

		_polap_log3_pipe "Rscript --vanilla ${_POLAPLIB_DIR}/polap-r-pairs.R \
	    -m ${_polap_var_mtcontigname} \
		  -t ${_polap_var_oga_contig}/contig.tab \
		  --out ${_polap_var_oga_reads}/${_pread_sel}/${i} \
      -w ${_arg_single_min} \
		  -r ${_arg_pair_min} \
		  -x ${_arg_bridge_min} \
      --all \
		  >${_polap_output_dest} 2>&1"

		_polap_log1 "  selecting long reads for ${_pread_sel}"
		_polap_log2 "    input1: ${_source_long_reads_fq}"
		_polap_log2 "    input2: ${_polap_var_oga_reads}/${_pread_sel}/${i}/${_read_names}.names"
		_polap_log2 "    output: ${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz"

		_polap_log3_pipe "seqtk subseq \
		    ${_source_long_reads_fq} \
		    ${_polap_var_oga_reads}/${_pread_sel}/${i}/${_read_names}.names |\
		    gzip >${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz"

		_polap_log1 "  sampling reads ..."
		_polap_log2 "    input1: ${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz"
		local _contig_length_bp=$(_polap_utility_convert_bp ${CONTIG_LENGTH})
		_polap_log2 "    input2 (seed contig size): ${_contig_length_bp}"
		_polap_utility_get_contig_length \
			"${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz" \
			"${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.len"
		local _seeds_length=$(<"${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.len")
		_polap_log2_cat "${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.len"
		local _seeds_length_bp=$(_polap_utility_convert_bp ${_seeds_length})
		_polap_log2 "    result1 (total size of reads mapped contigs): ${_seeds_length_bp}"
		local _expected_organelle_coverage=$((_seeds_length / CONTIG_LENGTH))
		_polap_log2 "    result2 (expected organelle coverage): ${_expected_organelle_coverage}x"

		# sampling rate
		local _rate=$(echo "scale=9; ${_arg_coverage_oga}/$_expected_organelle_coverage" | bc)

		_polap_log0 "rate: ${_rate}"

		if [[ "$_expected_organelle_coverage" -gt "${_arg_coverage_oga}" ]]; then

			# sampling rate
			local _rate=$(echo "scale=9; ${_arg_coverage_oga}/$_expected_organelle_coverage" | bc)
			echo "${_rate}" "${_polap_var_oga_subsample}/${_pread_sel}/${i}.txt"

			if [[ "${_arg_coverage_check}" == "on" ]]; then
				_polap_log0 "  long-read data reduction by rate of ${_rate} <= COV[${_arg_coverage_oga}] / long-read organelle coverage[$_expected_organelle_coverage]"

				_polap_log1 "    sampling long-read data by ${_rate} ... wait ..."
				_polap_lib_random-get
				local _random_seed=${_polap_var_random_number}
				# local _random_seed=11
				_polap_log1 "    random seed for reducing long reads mapped on potential seed contigs: ${_random_seed}"
				rm -f "${_polap_var_oga_subsample}/${_pread_sel}/${i}.fq.gz"
				_polap_log3_pipe "seqkit sample \
          -p ${_rate} \
          -s ${_random_seed} \
			    ${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz \
          -o ${_polap_var_oga_subsample}/${_pread_sel}/${i}.fq.gz \
          2>${_polap_output_dest}"
				_polap_log3_pipe "echo ${_random_seed} >${_polap_var_oga_subsample}/${_pread_sel}/${i}.random.seed.${_random_seed}"
			else
				_polap_log0 "    no reduction of the long-read data because of the option --no-coverage-check: expected coverage: ${_expected_organelle_coverage}"
				_polap_log3_cmd ln -s $(realpath "${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz") "${_polap_var_oga_subsample}/${_pread_sel}/${i}.fq.gz"
			fi
		else
			_polap_log0 "    no reduction of the long-read data because $_expected_organelle_coverage < ${_arg_coverage_oga}"
			_polap_log3_cmd ln -s $(realpath "${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz") "${_polap_var_oga_subsample}/${_pread_sel}/${i}.fq.gz"

		fi

		# _polap_log1 "  flye assembly for ${_pread_sel}"
		# _polap_log2 "    input1: ${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz"
		# _polap_log2 "    output: ${_polap_var_oga_flye}/${_pread_sel}/${i}/30-contigger/graph_final.gfa"
		# # if [[ "${_arg_plastid}" == "on" ]]; then
		# # 	CONTIG_LENGTH=$((CONTIG_LENGTH * 3))
		# # fi
		# local _command1="flye \
		#     ${_arg_flye_data_type} \
		#     ${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz \
		#   --out-dir ${_polap_var_oga_flye}/${_pread_sel}/${i} \
		#   --threads ${_arg_threads}"
		# if [[ "${_arg_flye_asm_coverage}" -gt 0 ]]; then
		# 	_command1+=" \
		#   --asm-coverage ${_arg_flye_asm_coverage} \
		#   --genome-size $CONTIG_LENGTH"
		# fi
		# if [[ "${_arg_menu[2]}" == "polishing" ]]; then
		# 	_command1+=" \
		#   --resume"
		# else
		# 	_command1+=" \
		#   --stop-after contigger"
		# fi
		# _command1+=" \
		#   2>${_polap_output_dest}"
		#
		# if [[ "${_arg_flye}" == "on" ]]; then
		# 	_polap_log3_pipe "${_command1}"
		# else
		# 	_polap_log0 "No flye run in test-reads"
		# fi

	}
}

function _polap_lib_oga-flye-select-reads {
	# polap-cmd-oga: test-reads
	local _pread_sel="ptgaul-reads"
	local _read_names="ptgaul"

	_polap_utility_get_contig_length \
		"${_polap_var_oga_contig}/contig.fa" \
		"${_polap_var_oga_contig}/contig_total_length.txt"
	local CONTIG_LENGTH=$(<"${_polap_var_oga_contig}/contig_total_length.txt")
	_polap_log2_cat "${_polap_var_oga_contig}/contig_total_length.txt"

	# local rate_lfq=$(echo "scale=9; ${len_subsample_lfq}/$len_source_lfq" | bc)

	_polap_log2 "  creating folders for read selection type: ${_pread_sel}"
	_polap_log3_cmd mkdir -p "${_polap_var_oga_flye}/${_pread_sel}"

	local i=$(<"${_polap_var_oga_contig}/index.txt")
	{

		_polap_log1 "  flye assembly for ${_pread_sel}"
		_polap_log2 "    input1: ${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz"
		_polap_log2 "    output: ${_polap_var_oga}"
		# if [[ "${_arg_plastid}" == "on" ]]; then
		# 	CONTIG_LENGTH=$((CONTIG_LENGTH * 3))
		# fi
		local _command1="flye \
      ${_arg_flye_data_type} \
      ${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz \
		  --out-dir ${_polap_var_oga} \
		  --threads ${_arg_threads}"
		if [[ "${_arg_flye_asm_coverage}" -gt 0 ]]; then
			_command1+=" \
		  --asm-coverage ${_arg_flye_asm_coverage} \
		  --genome-size $CONTIG_LENGTH"
		fi
		_command1+=" \
		  2>${_polap_output_dest}"

		if [[ "${_arg_flye}" == "on" ]]; then
			_polap_log3_pipe "${_command1}"
		else
			_polap_log0 "No flye run in test-reads"
		fi

	}
}

# not tested yet
function _polap_lib_oga-select-reads {
	# use the ptgaul mapping
	#
	# polap-cmd-oga: test-reads
	local _pread_sel="ptgaul-reads"
	local _read_names="ptgaul"

	_polap_utility_get_contig_length \
		"${_polap_var_oga_contig}/contig.fa" \
		"${_polap_var_oga_contig}/contig_total_length.txt"
	local CONTIG_LENGTH=$(<"${_polap_var_oga_contig}/contig_total_length.txt")
	_polap_log2_cat "${_polap_var_oga_contig}/contig_total_length.txt"

	_polap_log1 "  determines which long-read data to use ..."
	local _source_long_reads_fq=""
	_polap_oga_determine-long-read-file _source_long_reads_fq

	_polap_log2 "  creating folders for read selection type: ${_pread_sel}"
	_polap_log3_cmd mkdir -p "${_polap_var_oga_reads}/${_pread_sel}"
	_polap_log3_cmd mkdir -p "${_polap_var_oga_seeds}/${_pread_sel}"
	_polap_log3_cmd mkdir -p "${_polap_var_oga_subsample}/${_pread_sel}"

	{
		local i=0
		_polap_log3_cmd mkdir -p "${_polap_var_oga_reads}/${_pread_sel}/${i}"

		_polap_log3_pipe "Rscript --vanilla ${_POLAPLIB_DIR}/polap-r-pairs.R \
	    -m ${_polap_var_mtcontigname} \
		  -t ${_polap_var_oga_contig}/contig.tab \
		  --out ${_polap_var_oga_reads}/${_pread_sel}/${i} \
      -w ${_arg_single_min} \
		  -r ${_arg_pair_min} \
		  -x ${_arg_bridge_min} \
      --all \
		  >${_polap_output_dest} 2>&1"

		_polap_log1 "  selecting long reads for ${_pread_sel}"
		_polap_log2 "    input1: ${_source_long_reads_fq}"
		_polap_log2 "    input2: ${_polap_var_oga_reads}/${_pread_sel}/${i}/${_read_names}.names"
		_polap_log2 "    output: ${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz"

		_polap_log3_pipe "seqtk subseq \
		    ${_source_long_reads_fq} \
		    ${_polap_var_oga_reads}/${_pread_sel}/${i}/${_read_names}.names |\
		    gzip >${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz"

		_polap_log1 "  sampling reads ..."
		_polap_log2 "    input1: ${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz"
		local _contig_length_bp=$(_polap_utility_convert_bp ${CONTIG_LENGTH})
		_polap_log2 "    input2 (seed contig size): ${_contig_length_bp}"
		_polap_utility_get_contig_length \
			"${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz" \
			"${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.len"
		local _seeds_length=$(<"${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.len")
		_polap_log2_cat "${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.len"
		local _seeds_length_bp=$(_polap_utility_convert_bp ${_seeds_length})
		_polap_log2 "    result1 (total size of reads mapped contigs): ${_seeds_length_bp}"
		local _expected_organelle_coverage=$((_seeds_length / CONTIG_LENGTH))
		_polap_log2 "    result2 (expected organelle coverage): ${_expected_organelle_coverage}x"

		if [[ "$_expected_organelle_coverage" -gt "${_arg_coverage_oga}" ]]; then

			# sampling rate
			local _rate=$(echo "scale=9; ${_arg_coverage_oga}/$_expected_organelle_coverage" | bc)
			echo "${_rate}" "${_polap_var_oga_subsample}/${_pread_sel}/${i}.txt"

			if [[ "${_arg_coverage_check}" == "on" ]]; then
				_polap_log0 "  long-read data reduction by rate of ${_rate} <= COV[${_arg_coverage_oga}] / long-read organelle coverage[$_expected_organelle_coverage]"

				_polap_log1 "    sampling long-read data by ${_rate} ... wait ..."
				_polap_lib_random-get
				local _random_seed=${_polap_var_random_number}
				# local _random_seed=11
				_polap_log1 "    random seed for reducing long reads mapped on potential seed contigs: ${_random_seed}"
				rm -f "${_polap_var_oga_subsample}/${_pread_sel}/${i}.fq.gz"
				_polap_log3_pipe "seqkit sample \
          -p ${_rate} \
          -s ${_random_seed} \
			    ${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz \
          -o ${_polap_var_oga_subsample}/${_pread_sel}/${i}.fq.gz \
          2>${_polap_output_dest}"
				_polap_log3_pipe "echo ${_random_seed} >${_polap_var_oga_subsample}/${_pread_sel}/${i}.random.seed.${_random_seed}"
			else
				_polap_log0 "    no reduction of the long-read data because of the option --no-coverage-check: expected coverage: ${_expected_organelle_coverage}"
				_polap_log3_cmd ln -s $(realpath "${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz") "${_polap_var_oga_subsample}/${_pread_sel}/${i}.fq.gz"
			fi
		else
			_polap_log0 "    no reduction of the long-read data because $_expected_organelle_coverage < ${_arg_coverage_oga}"
			_polap_log3_cmd ln -s $(realpath "${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz") "${_polap_var_oga_subsample}/${_pread_sel}/${i}.fq.gz"

		fi

		#
	}
}
