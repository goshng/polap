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

function _run_polap_simulate {
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

	help_message=$(
		cat <<'EOF'
Name:
  polap simulate - simulate sequencing data

Synopsis:
  polap simulate [options]

Description:
  polap simulate uses sequencing data simulation tools such as pbsim3
  to generate sequencing data.
  It creates o/reads.fq.gz.

Options:
  -g INT
    Genome size

  -p FASTA
    Plastid reference genome

  -m FASTA
    Mitochondrial reference genome

Examples:
  Simulate PacBio HiFi sequencing data from Actinidia arguta genome (10 Mb):
    polap get-mtdna --species "Actinidia arguta"
    cp o/00-bioproject/2-mtdna.fasta mito_ref.fasta
    polap get-mtdna --species "Actinidia arguta" --plastid
    cp o/00-bioproject/2-mtdna.fasta plastid_ref.fasta
    polap simulate hifi -g 10000000 -p plastid_ref.fasta -m mito_ref.fasta

  Simulate PacBio HiFi data for syncassemble test:
    polap simulate syncam wgs.fq.gz

  Simulate data with pt/mt reference:
    polap simulate hifi-reference -p ptdna-reference.fa -m mtdna-reference.fa
    polap simulate ont-reference -p ptdna-reference.fa -m mtdna-reference.fa
    polap simulate onthq-reference -p ptdna-reference.fa -m mtdna-reference.fa
    polap simulate onthq30-reference -p ptdna-reference.fa -m mtdna-reference.fa
    polap simulate illumina-reference -p ptdna-reference.fa -m mtdna-reference.fa

  Simulate ONT sequencing data from Anthoceros agrestis
    polap simulate ont-reference --species "Anthoceros agrestis"

  Save ~/.polap/profiles/sim1.yaml
    preset: sim1
    species: Anthoceros agrestis
    outdir: sim1
  Then,
    polap simulate ont-reference --preset sim1
    polap simulate ont-reference-nanosim --preset sim1 -l l.fastq

  Simulate nuclear genome:
    polap simulate nuclear 

TODO:
  Check the required options values.
  Add ONT and other data type simulation.
  Add NUMT/NUPT contamination in the nuclear reads.

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (1998–2018)

Author:
  Sang Chul Choi
EOF
	)

	# Display help message
	if [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]]; then
		local manfile=$(_polap_lib_man-convert_help_message "$help_message" "${_arg_menu[0]}")
		man "$manfile" >&3
		rm -f "$manfile"
		return
	fi

	# Display help message
	# [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
	fi

	# FIXME: check the three options
	# -g
	# -p
	# -m
	if [[ "${_arg_menu[1]}" == "hifi" ]]; then
		_polap_simulate_hifi
	elif [[ "${_arg_menu[1]}" == "ont" ]]; then
		_polap_simulate_ont
	elif [[ "${_arg_menu[1]}" == "hifi-numt" ]]; then
		_polap_simulate_hifi_numt
	elif [[ "${_arg_menu[1]}" == "syncasm" ]]; then
		_polap_simulate_syncassemble_hifi "${_arg_menu[2]}"
	elif [[ "${_arg_menu[1]}" == "hifi-reference" ]]; then
		_polap_simulate_hifi_reference
	elif [[ "${_arg_menu[1]}" == "ont-reference" ]]; then
		_polap_simulate_ont_reference
	elif [[ "${_arg_menu[1]}" == "ont-reference-nanosim" ]]; then
		_polap_simulate_ont_reference_nanosim
	elif [[ "${_arg_menu[1]}" == "onthq-reference" ]]; then
		_polap_simulate_onthq_reference
	elif [[ "${_arg_menu[1]}" == "onthq30-reference" ]]; then
		_polap_simulate_onthq30_reference
	elif [[ "${_arg_menu[1]}" == "illumina-reference" ]]; then
		_polap_simulate_illumina_reference
	elif [[ "${_arg_menu[1]}" == "nuclear" ]]; then
		_polap_lib_conda-ensure_conda_env polap-sim || exit 1
		python "${_POLAPLIB_DIR}/polap-py-make-nuclear-from-busco-aa.py" \
			-o "${_arg_outdir}/nuclear" \
			--genome-size 10000000
		conda deactivate
	else
		_polap_log0 "ERROR: no such simulation: ${_arg_menu[1]}"
	fi

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

_polap_simulate_hifi() {

	# Simulate PacBio HiFi using PBSIM2
	# wget https://github.com/PacificBiosciences/ccs/releases/download/v6.4.0/ccs.tar.gz
	# tar zxf ccs.tar.gz
	# git clone https://github.com/yukiteruono/pbsim3.git
	local simdir="${_arg_outdir}"/sim
	mkdir -p "${simdir}"

	if [[ ! -s "${simdir}"/nuclear_ref.fasta ]]; then
		python3 -c "import random; print('>nuclear_sim'); print(''.join(random.choices('ACGT', k=${_arg_genomesize})))" >"${simdir}"/nuclear_ref.fasta
	fi

	if [[ ! -s "${simdir}"/pbsim_plastid_0001.fq.gz ]]; then
		pbsim --strategy wgs \
			--method qshmm \
			--qshmm "${_POLAPLIB_DIR}"/pbsim3/data/QSHMM-RSII.model \
			--depth 500 \
			--genome "${_arg_unpolished_fasta}" \
			--pass-num 10 \
			--id-prefix plastid \
			--prefix "${simdir}"/pbsim_plastid

		"${_POLAPLIB_DIR}"/ccs \
			"${simdir}"/pbsim_plastid_0001.bam \
			"${simdir}"/pbsim_plastid_0001.fq.gz
	fi

	if [[ ! -s "${simdir}"/pbsim_mito_0001.fq.gz ]]; then
		pbsim --strategy wgs \
			--method qshmm \
			--qshmm "${_POLAPLIB_DIR}"/pbsim3/data/QSHMM-RSII.model \
			--depth 50 \
			--genome "${_arg_min_read_length}" \
			--pass-num 10 \
			--id-prefix mito \
			--prefix "${simdir}"/pbsim_mito

		"${_POLAPLIB_DIR}"/ccs \
			"${simdir}"/pbsim_mito_0001.bam \
			"${simdir}"/pbsim_mito_0001.fq.gz
	fi

	if [[ ! -s "${simdir}"/pbsim_nuclear_0001.fq.gz ]]; then
		pbsim --strategy wgs \
			--method qshmm \
			--qshmm "${_POLAPLIB_DIR}"/pbsim3/data/QSHMM-RSII.model \
			--depth 10 \
			--genome "${simdir}"/nuclear_ref.fasta \
			--pass-num 10 \
			--id-prefix nuclear \
			--prefix "${simdir}"/pbsim_nuclear

		"${_POLAPLIB_DIR}"/ccs \
			"${simdir}"/pbsim_nuclear_0001.bam \
			"${simdir}"/pbsim_nuclear_0001.fq.gz
	fi

	zcat "${simdir}"/*.fq.gz >"${_arg_outdir}"/reads.fq

	_polap_log0 "output: ${_arg_outdir}"/reads.fq
}

_polap_simulate_ont() {

	# Simulate PacBio HiFi using PBSIM2
	# wget https://github.com/PacificBiosciences/ccs/releases/download/v6.4.0/ccs.tar.gz
	# tar zxf ccs.tar.gz
	# git clone https://github.com/yukiteruono/pbsim3.git
	local simdir="${_arg_outdir}"/sim
	mkdir -p "${simdir}"

	if [[ ! -s "${simdir}"/nuclear_ref.fasta ]]; then
		python3 -c "import random; print('>nuclear_sim'); print(''.join(random.choices('ACGT', k=${_arg_genomesize})))" >"${simdir}"/nuclear_ref.fasta
	fi

	if [[ ! -s "${simdir}"/pbsim_plastid_0001.fq.gz ]]; then
		pbsim --strategy wgs \
			--method qshmm \
			--qshmm "${_POLAPLIB_DIR}"/pbsim3/data/QSHMM-ONT.model \
			--depth 500 \
			--genome "${_arg_unpolished_fasta}" \
			--id-prefix plastid \
			--prefix "${simdir}"/pbsim_plastid
	fi

	if [[ ! -s "${simdir}"/pbsim_mito_0001.fq.gz ]]; then
		pbsim --strategy wgs \
			--method qshmm \
			--qshmm "${_POLAPLIB_DIR}"/pbsim3/data/QSHMM-ONT.model \
			--depth 50 \
			--genome "${_arg_min_read_length}" \
			--id-prefix mito \
			--prefix "${simdir}"/pbsim_mito
	fi

	if [[ ! -s "${simdir}"/pbsim_nuclear_0001.fq.gz ]]; then
		pbsim --strategy wgs \
			--method qshmm \
			--qshmm "${_POLAPLIB_DIR}"/pbsim3/data/QSHMM-ONT.model \
			--depth 10 \
			--genome "${simdir}"/nuclear_ref.fasta \
			--id-prefix nuclear \
			--prefix "${simdir}"/pbsim_nuclear
	fi

	zcat "${simdir}"/*.fq.gz >"${_arg_outdir}"/reads.fq

	_polap_log0 "output: ${_arg_outdir}"/reads.fq
}

# ref:
# plastid ref:
# mito ref:
_polap_simulate_hifi_numt() {
	# === Inputs ===
	NUCLEAR_REF="reference/nuclear_ref.fasta"
	INSERT_REF_PLASTID="reference/plastid_ref.fasta"
	INSERT_REF_MITO="reference/mito_ref.fasta"

	NUM_INSERTS=5
	INSERT_LEN=1000
	DEPTH=20
	OUTDIR="sim_output"
	PBSIM_BIN="pbsim"
	QSHMM="pbsim3/data/QSHMM-RSII.model"

	# === Output files ===
	BED_FILE="${OUTDIR}/numt_nupt.bed"
	FASTA_FILE="${OUTDIR}/numt_nupt_sequences.fasta"
	MASKED="${OUTDIR}/nuclear_masked.fasta"
	INJECTED="${OUTDIR}/nuclear_with_inserts.fasta"
	READS="${OUTDIR}/pbsim_nuclear_with_numt_nupt_0001.fq.gz"
	READ_NAMES="${OUTDIR}/reads_from_inserts.txt"

	mkdir -p "$OUTDIR"

	# === Step 1: Create BED and FASTA for NUMT/NUPT
	echo "🧬 Extracting and generating NUMT/NUPT fragments..."
	NUCLEN=$(awk '/^>/ {if (seqlen) print seqlen; seqlen=0; next} {seqlen+=length($0)} END{print seqlen}' "$NUCLEAR_REF")

	START_POS=$(shuf -i 1000-$((NUCLEN - INSERT_LEN - 1000)) -n $NUM_INSERTS)
	INSERT_IDS=()
	>"$BED_FILE"
	>"$FASTA_FILE"

	for i in $(seq 1 $NUM_INSERTS); do
		pos=$(echo "$START_POS" | sed -n "${i}p")
		end=$((pos + INSERT_LEN))
		chr=$(grep '^>' "$NUCLEAR_REF" | head -n1 | sed 's/^>//')

		label=$(if ((i <= NUM_INSERTS / 2)); then echo "nupt"; else echo "numt"; fi)
		insert_id="${label}_${i}_pos${pos}"
		INSERT_IDS+=("$insert_id")

		ref_source=$(if [[ "$label" == "nupt" ]]; then echo "$INSERT_REF_PLASTID"; else echo "$INSERT_REF_MITO"; fi)
		frag=$(seqkit subseq -r 1:$INSERT_LEN "$ref_source" | seqkit seq -n -i | sed -n 2p)

		echo -e "${chr}\t${pos}\t${end}" >>"$BED_FILE"
		echo -e ">${insert_id}\n${frag}" >>"$FASTA_FILE"
	done

	# === Step 2: Mask nuclear genome
	bedtools maskfasta -fi "$NUCLEAR_REF" -bed "$BED_FILE" -fo "$MASKED" -mc N

	# === Step 3: Inject fragments
	python "${_POLAPLIB_DIR}"/polap-py-inject-sequences.py "$MASKED" "$INSERTS" "$INJECTED"
	# 	python3 - "$MASKED" "$FASTA_FILE" "$INJECTED" <<'EOF'
	# from Bio import SeqIO
	# import sys
	# masked_fasta, inserts_fasta, output_fasta = sys.argv[1:]
	# inserts = list(SeqIO.parse(inserts_fasta, "fasta"))
	# insert_index = 0
	# def inject(seq):
	#     global insert_index
	#     new_seq, i = [], 0
	#     while i < len(seq):
	#         if seq[i] == "N":
	#             frag = str(inserts[insert_index].seq)
	#             new_seq.append(frag)
	#             i += len(frag)
	#             insert_index += 1
	#         else:
	#             new_seq.append(seq[i])
	#             i += 1
	#     return ''.join(new_seq)
	# with open(output_fasta, "w") as out:
	#     for record in SeqIO.parse(masked_fasta, "fasta"):
	#         record.seq = inject(str(record.seq))
	#         SeqIO.write(record, out, "fasta")
	# EOF

	# === Step 4: Simulate reads
	"$PBSIM_BIN" --strategy wgs \
		--method qshmm \
		--qshmm "$QSHMM" \
		--depth "$DEPTH" \
		--length-mean 9000 \
		--length-sd 10 \
		--genome "$INJECTED" \
		--prefix "${OUTDIR}/pbsim_nuclear_with_numt_nupt"

	# === Step 5: Align and track NUMT/NUPT-derived reads
	minimap2 -ax map-pb "$INJECTED" "$READS" |
		samtools view -b -o "${OUTDIR}/aln.bam" -
	samtools sort -o "${OUTDIR}/aln.sorted.bam" "${OUTDIR}/aln.bam"
	samtools index "${OUTDIR}/aln.sorted.bam"

	bedtools intersect -a <(samtools view -F 4 "${OUTDIR}/aln.sorted.bam" |
		awk '{print $1"\t"$3"\t"$4"\t"($4 + length($10))}') \
		-b "$BED_FILE" -wa | cut -f1 | sort | uniq >"$READ_NAMES"

	echo "✅ Pipeline complete"
	echo "📁 Injected genome: $INJECTED"
	echo "📁 Simulated reads: $READS"
	echo "📁 Reads from injected regions: $READ_NAMES"

}

_polap_simulate_syncassemble_hifi() {
	local READS_OUT="${1:-out.fq.gz}"

	local HIFI_RL_MEAN=10000
	local HIFI_RL_SD=1500
	local HIFI_SUB=0.002
	local HIFI_INS=0.0003
	local HIFI_DEL=0.0003
	local ONT_RL_MEAN=10000
	local ONT_RL_SD=3000
	local ONT_SUB=0.01
	local ONT_INS=0.01
	local ONT_DEL=0.01
	local OUTDIR="${_arg_outdir}"

	# simulate or real
	python3 "${_POLAPLIB_DIR}/polap-py-simulate-wgs-hifi.py" \
		--nuclear-size "10m" --nuclear-cov "5" \
		--mito-size "500k" --mito-cov "50" \
		--plastid-size "150k" --plastid-cov "250" \
		--rl-mean "$HIFI_RL_MEAN" \
		--rl-sd "$HIFI_RL_SD" \
		--sub "$HIFI_SUB" \
		--ins "$HIFI_INS" \
		--dele "$HIFI_DEL" \
		--emit-refs -o "$READS_OUT" \
		--summary "$OUTDIR/wgs_mix.summary.tsv"

}

_polap_simulate_hifi_reference() {
	# local READS_OUT="${1:-out.fq.gz}"

	_polap_lib_conda-ensure_conda_env polap-sim || exit 1
	python "${_POLAPLIB_DIR}/polap-py-simulate-wgs-reads.py" \
		--platform hifi \
		-o "${_arg_outdir}"/sim-hifi \
		--mt-fasta "${_arg_min_read_length}" \
		--pt-fasta "${_arg_unpolished_fasta}" \
		--nuc-depth 5 --mt-depth 50 --pt-depth 250
	conda deactivate
}

_polap_simulate_ont_reference() {
	local simdir="${_arg_outdir}/sim-ont"

	_polap_lib_conda-ensure_conda_env polap-sim || exit 1
	_polap_lib_simulate-ont-with-reference \
		--mt "${_arg_min_read_length}" \
		--pt "${_arg_unpolished_fasta}" \
		--nuc-size "${_arg_sim_nuc_size}" \
		${_arg_sim_nuc:+--nuc "$_arg_sim_nuc"} \
		--depth-nuc "${_arg_sim_depth_nuc}" \
		--depth-mt "${_arg_sim_depth_mt}" \
		--depth-pt "${_arg_sim_depth_pt}" \
		--mode badread \
		--tech duplex \
		--len-mean "${_arg_sim_len_mean}" \
		--len-sd 8000 \
		--threads 48 \
		--verbose \
		--out "${simdir}/simulate"

	conda deactivate
}

_polap_simulate_ont_reference_nanosim() {
	local simdir="${_arg_outdir}/sim"

	# mtDNA: fetch all, list, pick the longest as representative
	# _polap_lib_conda-ensure_conda_env polap || exit 1
	# local organelle_type
	# for organelle_type in pt mt; do
	# 	polap_fetch_select_organelle \
	# 		--species "${_arg_species}" \
	# 		--type ${organelle_type} \
	# 		--out-all "${simdir}"/${organelle_type}.all.fa \
	# 		--out-tsv "${simdir}"/${organelle_type}.list.tsv \
	# 		--out-rep-fa "${simdir}"/${organelle_type}.rep.fa \
	# 		--out-rep-id "${simdir}"/${organelle_type}.rep.id \
	# 		--strategy longest
	# done
	# conda deactivate

	# See the following options:
	# badread options:
	# --err-model
	# -q-model
	_polap_lib_conda-ensure_conda_env polap-sim || exit 1
	_polap_lib_simulate-ont-with-reference \
		--mt "${simdir}"/mt.rep.fa \
		--pt "${simdir}"/pt.rep.fa \
		--nuc-size 100000 \
		--depth-nuc 1x \
		--depth-mt 5x \
		--depth-pt 5x \
		--mode nanosim-nu \
		--in-fastq "${_arg_long_reads}" \
		--verbose \
		--out "${simdir}/simulate"
	conda deactivate

}

_polap_simulate_onthq_reference() {

	python "${_POLAPLIB_DIR}/polap-py-simulate-wgs-reads.py" \
		--platform onthq -o "${_arg_outdir}"/onthq \
		--mt-fasta "${_arg_min_read_length}" \
		--pt-fasta "${_arg_unpolished_fasta}" \
		--nuc-depth 5 --mt-depth 50 --pt-depth 250
}

_polap_simulate_onthq30_reference() {

	python "${_POLAPLIB_DIR}/polap-py-simulate-wgs-reads.py" \
		--platform onthq30 -o "${_arg_outdir}"/onthq30 \
		--mt-fasta "${_arg_min_read_length}" \
		--pt-fasta "${_arg_unpolished_fasta}" \
		--nuc-depth 5 --mt-depth 50 --pt-depth 250 \
		--sub 0.0007 --ins 0.0012 --dele 0.0012 --hp-rate 0.03 --qconst 32
}

_polap_simulate_illumina_reference() {

	python "${_POLAPLIB_DIR}/polap-py-simulate-wgs-reads.py" \
		--platform illumina -o "${_arg_outdir}"/illumina \
		--mt-fasta "${_arg_min_read_length}" \
		--pt-fasta "${_arg_unpolished_fasta}" \
		--nuc-depth 5 --mt-depth 50 --pt-depth 250
}
