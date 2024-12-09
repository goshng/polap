#!/usr/bin/bash

_arg_short_read1=$1
_arg_short_read2=$2
_arg_outdir=o
_polap_output_dest=/dev/null

if [[ ${_arg_short_read1} = *.fastq || ${_arg_short_read1} = *.fq ]]; then
	cat "${_arg_short_read1}" "${_arg_short_read2:-/dev/null}" |
		awk 'NR % 4 == 2' | sort | tr NT TN |
		ropebwt2 -LR 2>"${_polap_output_dest}" |
		tr NT TN |
		msbwt convert "${_arg_outdir}"/msbwt \
			>/dev/null 2>&1
elif [[ ${_arg_short_read1} = *.fq.gz ]] || [[ ${_arg_short_read1} = *.fastq.gz ]]; then
	zcat "${_arg_short_read1}" "${_arg_short_read2}" |
		awk 'NR % 4 == 2' | sort | tr NT TN |
		ropebwt2 -LR 2>"${_polap_output_dest}" |
		tr NT TN |
		msbwt convert "${_arg_outdir}"/msbwt \
			>/dev/null 2>&1
else
	echo "ERROR: short-read1: no fastq or fq file: ${_arg_short_read1}"
	echo "ERROR: short-read2: no fastq or fq file: ${_arg_short_read2}"
fi
