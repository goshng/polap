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

################################################################################
# mafft alignment
################################################################################
function _run_polap_mafft-mtdna {
  # Enable debugging if _POLAP_DEBUG is set
  [ "$_POLAP_DEBUG" -eq 1 ] && set -x
  _polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

  # Set verbosity level: stderr if verbose >= 2, otherwise discard output
  local _polap_output_dest="/dev/null"
  [ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

  # Set paths for bioproject data
  source "${_POLAPLIB_DIR}/polap-variables-common.sh"

  # Help message
  local help_message=$(
    cat <<HEREDOC
# Compare the known mitochondrial DNA (mtDNA) sequence with the newly assembled 
# one for verification purposes using mafft.
#
# We utilize progressiveMauve to align two sequences and 
# calculate the total length of the LCB, which is then divided by the length 
# of a known mtDNA sequence. This approach provides a basic method for 
# comparing two sequences, albeit not highly sophisticated.
#
# Arguments:
#   -i ${_arg_inum}
# Inputs:
#   a.fasta: known mtDNA in fasta format
#   b.fasta: another DNA sequence in fasta format
# Outputs:
#   the ratio of LCB total length divided by the known mtDNA sequence
Example: $(basename "$0") ${_arg_menu[0]} -a o/00-bioproject/2-mtdna.fasta -b o/1/assembly.fasta
HEREDOC
  )

  # Display help message
  [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

  _mafft_dir="${_arg_outdir}"
  mkdir -p "${_mafft_dir}"
  cat "${_arg_short_read1}" >"${_mafft_dir}/in.fa"
  cat "${_arg_short_read2}" >>"${_mafft_dir}/in.fa"

  _polap_log3_pipe "mafft \
    --auto \
    ${_mafft_dir}/in.fa \
    >${_mafft_dir}/out.mafft \
    2>${_polap_output_dest}"

  _polap_log3_pipe "Rscript --vanilla ${_POLAPLIB_DIR}/run-polap-r-mafft.R \
    --input ${_mafft_dir}/out.mafft \
    --out ${_mafft_dir}/out.txt 2>${_polap_output_dest}"

  pident_mafft=$(grep -oP '(?<=Percent Identity: )\S+' "${_mafft_dir}/out.txt")
  pident_mafft=${pident_mafft%\%}

  echo "$pident_mafft" >"${_mafft_dir}/pident.txt"

  _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
  # Disable debugging if previously enabled
  [ "$_POLAP_DEBUG" -eq 1 ] && set +x
  return 0
}

################################################################################
# mauve alignment
################################################################################
function _run_polap_mauve-mtdna {
  # Enable debugging if _POLAP_DEBUG is set
  [ "$_POLAP_DEBUG" -eq 1 ] && set -x
  _polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

  # Set verbosity level: stderr if verbose >= 2, otherwise discard output
  local _polap_output_dest="/dev/null"
  [ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

  # Set paths for bioproject data
  source "${_POLAPLIB_DIR}/polap-variables-common.sh"

  # Help message
  local help_message=$(
    cat <<HEREDOC
# Compare the known mitochondrial DNA (mtDNA) sequence with the newly assembled 
# one for verification purposes.
#
# We utilize progressiveMauve to align two sequences and 
# calculate the total length of the LCB, which is then divided by the length 
# of a known mtDNA sequence. This approach provides a basic method for 
# comparing two sequences, albeit not highly sophisticated.
#
# Arguments:
#   -i ${_arg_inum}
# Inputs:
#   a.fasta: known mtDNA in fasta format
#   b.fasta: another DNA sequence in fasta format
# Outputs:
#   the ratio of LCB total length divided by the known mtDNA sequence
Example: $(basename "$0") ${_arg_menu[0]} -a o/00-bioproject/2-mtdna.fasta -b o/1/assembly.fasta
HEREDOC
  )

  # Display help message
  [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

  _polap_log3_pipe "progressiveMauve \
    --output=${_arg_outdir}/mt.xmfa \
    ${_arg_short_read1} \
    ${_arg_short_read2} \
    >${_polap_output_dest} 2>&1"

  # awk 'NR > 1 {sum += $2 - $1} END {print sum}' "${_arg_outdir}/mt.xmfa.backbone" >"${_arg_outdir}/mt.identity.length.txt"
  awk 'NR > 1 {sum += ($1 > $2 ? $1 - $2 : $2 - $1)} END {print sum}' "${_arg_outdir}/mt.xmfa.backbone" >"${_arg_outdir}/mt.identity.length.txt"

  local _alen=$(<"${_arg_outdir}/mt.identity.length.txt")
  _polap_utility_get_contig_length \
    "${_arg_short_read1}" \
    "${_arg_outdir}/mt.reference.length.txt"
  local _blen=$(<"${_arg_outdir}/mt.reference.length.txt")
  local _percent_identity=$(echo "scale=3; ${_alen}/${_blen}" | bc)
  _polap_log0 "mauve_lcb_length_coverage: ${_percent_identity}"

  _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
  # Disable debugging if previously enabled
  [ "$_POLAP_DEBUG" -eq 1 ] && set +x
  return 0
}

function _run_polap_compare2ptdna {
  # Enable debugging if _POLAP_DEBUG is set
  [ "$_POLAP_DEBUG" -eq 1 ] && set -x
  _polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

  # Set verbosity level: stderr if verbose >= 2, otherwise discard output
  local _polap_output_dest="/dev/null"
  [ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

  # Set paths for bioproject data
  source "${_POLAPLIB_DIR}/polap-variables-common.sh"

  # Help message
  local help_message=$(
    cat <<HEREDOC
# Compare the plastid DNA (ptDNA) sequence with the newly assembled 
# one for verification purposes.
#
# Inputs:
#   a.fasta: known mtDNA in fasta format
#   b.fasta: another DNA sequence in fasta format
# Outputs:
#   the ratio of LCB total length divided by the known mtDNA sequence
Example: $(basename "$0") ${_arg_menu[0]} -a o/00-bioproject/2-mtdna.fasta -b o/1/assembly.fasta
HEREDOC
  )

  # Display help message
  [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

  local _ptdir="${_arg_outdir}"
  _polap_log3_pipe "python ${_POLAPLIB_DIR}/run-polap-py-compare2ptdna.py \
    --seq1 ${_arg_short_read1} \
    --seq2 ${_arg_short_read2} \
		--out ${_ptdir} \
		2>${_polap_output_dest}"

  # Disable debugging if previously enabled
  [ "$_POLAP_DEBUG" -eq 1 ] && set +x
  return 0
}

################################################################################
# Blast the final mtDNA sequence against mitochondrial and plastid genes.
################################################################################
function _run_polap_blast-mtdna {
  # Enable debugging if _POLAP_DEBUG is set
  [ "$_POLAP_DEBUG" -eq 1 ] && set -x

  # Set verbosity level: stderr if verbose >= 2, otherwise discard output
  local _polap_output_dest="/dev/null"
  [ "${_arg_verbose}" -ge 2 ] && _polap_output_dest="/dev/stderr"

  # Grouped file path declarations
  source "${_POLAPLIB_DIR}/polap-variables-common.sh"

  # Set variables for file paths
  local MTAA="${_POLAPLIB_DIR}/polap-mt.1.c70.3.faa"
  local PTAA="${_POLAPLIB_DIR}/polap-pt.2.c70.3.faa"
  local _polap_var_assembly="${_polap_var_ga}"
  local _polap_var_chloroplot="${_polap_var_assembly}/60-chloroplot"

  # Help message
  local help_message=$(
    cat <<HEREDOC
# Blast the final mtDNA sequence against the POLAP organelle gene set.
# Arguments:
#   -f ${_arg_final_assembly}
# Inputs:
#   ${_arg_final_assembly}: mtDNA sequence
# Outputs:
#
Example: $(basename "$0") ${_arg_menu[0]} [-i|--inum <arg>] -f mt.1.fa -o [${_arg_outdir}]
HEREDOC
  )

  # Display help message if needed
  if [[ ${_arg_menu[1]} == "help" ]]; then
    echoerr "${help_message}"
    exit $EXIT_SUCCESS
  fi

  # Check if the input FASTA file exists
  check_file_existence "${FA}"

  # Clean up and create the chloroplot directory
  if [ -d "${_polap_var_chloroplot}" ]; then
    rm -rf "${_polap_var_chloroplot}"
    _polap_log2 "LOG: ${_polap_var_chloroplot} is deleted."
  fi
  mkdir -p "${_polap_var_chloroplot}"
  _polap_log2 "LOG: ${_polap_var_chloroplot} is created."

  # Create the BLAST database for the mtDNA sequences
  makeblastdb -dbtype nucl \
    -in "${FA}" \
    -out "${_polap_var_chloroplot}/dna" \
    >/dev/null 2>&1
  _polap_log2 "LOG: BLASTDB of the contig sequences: ${_polap_var_chloroplot}/dna"

  # BLAST mtDNA gene annotation
  _polap_log1 "LOG: executing tblastn ... be patient!"
  _polap_log2 "LOG: BLAST of the mitochondrial proteins against ${_polap_var_chloroplot}/dna ..."

  tblastn -query "$MTAA" \
    -db "${_polap_var_chloroplot}/dna" \
    -out "${_polap_var_chloroplot}/mtaa.blast" \
    -evalue 1e-30 \
    -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle salltitles' \
    -num_threads "${_arg_threads}" \
    >/dev/null 2>&1

  _polap_log2 "INFO: tblastn complete"

  # Process the tblastn results
  "${_POLAPLIB_DIR}"/run-polap-r-genes.R \
    "${_polap_var_chloroplot}/mtaa.blast" \
    "${_polap_var_chloroplot}/mtaa.blast.bed" \
    >/dev/null 2>&1

  sort -k1,1 -k2,2n \
    "${_polap_var_chloroplot}/mtaa.blast.bed" \
    >"${_polap_var_chloroplot}/mtaa.blast.sorted.bed"

  # Create directory for gene bed files
  mkdir "${_polap_var_chloroplot}/mtaa.bed"

  "${_POLAPLIB_DIR}"/run-polap-r-genes-bed4.R \
    "${_polap_var_chloroplot}/mtaa.blast" \
    "${_polap_var_chloroplot}/mtaa.blast.bed4"

  _polap_log2 "LOG: counting mitochondrial genes in the contigs ..."
  bedtools merge -i "${_polap_var_chloroplot}/mtaa.blast.sorted.bed" >"${_polap_var_chloroplot}/mtaa.blast.sorted.bed.txt"

  _polap_log2_file "${_polap_var_chloroplot}/mtaa.blast.sorted.bed.txt"

  # Process individual genes and their annotations
  local _polap_var_i=1
  while IFS= read -r gene; do
    printf "%s\n" "$gene" >"${_polap_var_chloroplot}/${_polap_var_i}.gene.bed"

    bedtools intersect -wa \
      -a "${_polap_var_chloroplot}/mtaa.blast.bed4" \
      -b "${_polap_var_chloroplot}/${_polap_var_i}.gene.bed" \
      >"${_polap_var_chloroplot}/${_polap_var_i}.bed4"

    "${_POLAPLIB_DIR}"/run-polap-r-blast-mtdna-1-determine-gene.R \
      "${_polap_var_chloroplot}/${_polap_var_i}.bed4" \
      "${_POLAPLIB_DIR}/polap-mt.1.c70.3.faa.name" \
      "${_POLAPLIB_DIR}/polap-mtgenes.txt" \
      "${_polap_var_chloroplot}/${_polap_var_i}.gene" \
      "${_polap_var_chloroplot}/${_polap_var_i}.bed4.description" \
      "${_polap_var_chloroplot}/${_polap_var_i}.bed4.count" \
      >/dev/null 2>&1

    _polap_var_i=$((_polap_var_i + 1))

  done <"${_polap_var_chloroplot}/mtaa.blast.sorted.bed.txt"

  # Combine the gene annotations into one file
  paste \
    <(cat "${_polap_var_chloroplot}"/*.gene.bed) \
    <(cat "${_polap_var_chloroplot}"/*.gene) \
    >"${_polap_var_chloroplot}/annotation.bed"

  _polap_log1_file "Annotation: ${_polap_var_chloroplot}/annotation.bed"
  _polap_log1_file "Ambiguous gene annotations: ${_polap_var_chloroplot}/*.check"
  ls "${_polap_var_chloroplot}"/*.check \
    >"$_polap_output_dest"

  echoerr "NEXT: $(basename $0) gene-table-mtdna -o ${_arg_outdir} [-i ${_arg_inum}]"

  # Disable debugging if previously enabled
  [ "$_POLAP_DEBUG" -eq 1 ] && set +x
  return 0
}

################################################################################
# Compares the known mtDNA sequence and the assembled one.
################################################################################
function _run_polap_compare-mtdna {
  # Enable debugging if _POLAP_DEBUG is set
  [ "$_POLAP_DEBUG" -eq 1 ] && set -x
  _polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

  # Set verbosity level: stderr if verbose >= 2, otherwise discard output
  local _polap_output_dest="/dev/null"
  [ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

  # Set paths for bioproject data
  source "${_POLAPLIB_DIR}/polap-variables-common.sh"
  source "${_POLAPLIB_DIR}/polap-variables-common.sh"

  # Help message
  local help_message=$(
    cat <<HEREDOC
# Compares the known mtDNA sequence and the assembled one.
#
# Arguments:
#   -i ${_arg_inum}
# Inputs:
#   ${_polap_var_compare_mtdna3}
# Outputs:
Example: $(basename "$0") ${_arg_menu[0]} [-i|--inum <arg>] -f mt.1.fa -o [${_arg_outdir}]
HEREDOC
  )

  # Display help message
  [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

  check_folder_existence "${_arg_outdir}"

  local n1=$(cut -f1 "${_polap_var_project_mtdna_fasta2_accession}")
  local l1=$(seqkit stats -T "${_polap_var_project_mtdna_fasta2}" |
    csvtk cut -t -f sum_len |
    csvtk del-header)

  # Run blastn between known mtDNA and assembled mtDNA (first round)
  blastn -query "${_polap_var_compare_mtdna1}" \
    -subject "${_polap_var_project_mtdna_fasta2}" \
    -outfmt "6 qseqid qstart sseqid sstart sstrand" \
    >"${_polap_var_compare_oga_blastn1}"

  blastn -query "${_polap_var_compare_mtdna1}" \
    -subject "${_polap_var_project_mtdna_fasta2}" \
    >"${_polap_var_compare_oga_blastn1}.full"

  _polap_log2_file "${_polap_var_compare_oga_blastn1}"

  # Determine the strand orientation
  if [ -s "${_polap_var_compare_oga_blastn1}" ]; then
    local _polap_var_project_strand=$(cut -f5 "${_polap_var_compare_oga_blastn1}" | head -n 1)
  else
    _polap_log1 "No hit in ${_polap_var_compare_oga_blastn1}"
    local n1=$(cut -f1 "${_polap_var_project_mtdna_fasta2_accession}")
    printf "%s\t%d\t0\t0\n" ${n1} ${l1} >"${_polap_var_compare_mtdna_compare}"
    return
  fi

  # Reverse sequence if the strand is negative
  if [ "${_polap_var_project_strand}" = "plus" ]; then
    cp "${_polap_var_compare_mtdna1}" "${_polap_var_compare_mtdna2}"
  else
    seqkit seq -t dna -v -p -r "${_polap_var_compare_mtdna1}" -o "${_polap_var_compare_mtdna2}"
  fi

  # Run blastn (second round)
  blastn -query "${_polap_var_compare_mtdna2}" -subject "${_polap_var_project_mtdna_fasta2}" \
    -outfmt "6 qseqid qstart sseqid sstart sstrand" \
    >"${_polap_var_compare_oga_blastn2}"

  _polap_log2_file "${_polap_var_compare_oga_blastn2}"

  # Restart sequence alignment at the lowest start position
  local _polap_var_project_restart_position=$(sort -n -k4 "${_polap_var_compare_oga_blastn2}" | head -n 1 | cut -f2)
  _polap_log2 "LOG: restart position: ${_polap_var_project_restart_position}"
  seqkit restart -i "${_polap_var_project_restart_position}" "${_polap_var_compare_mtdna2}" -o "${_polap_var_compare_mtdna3}"

  # Run blastn (third round)
  blastn -query "${_polap_var_compare_mtdna3}" -subject "${_polap_var_project_mtdna_fasta2}" \
    -outfmt "6 qseqid qstart sseqid sstart sstrand length pident" \
    >"${_polap_var_compare_oga_blastn3}"

  _polap_log2_file "${_polap_var_compare_oga_blastn3}"

  # Analyze the length of the match using an R script
  "${_POLAPLIB_DIR}"/run-polap-r-assemble-bioproject-3-length-match.R \
    "${_polap_var_compare_oga_blastn3}" \
    "${_polap_var_compare_oga_blastn3_length}" \
    2>"$_polap_output_dest"

  _polap_log2_file "${_polap_var_compare_oga_blastn3_length}"

  local l2=$(<"${_polap_var_compare_oga_blastn3_length}")
  local c1=$(echo "scale=3; ${l2}/${l1}" | bc)
  _polap_log2 "Length of ${_polap_var_project_mtdna_fasta2}: ${l1}"
  _polap_log2 "Length of match alignment: ${l2}"
  _polap_log1 "length coverage: ${c1}"
  printf "%s\t%d\t%d\t%f\n" ${n1} ${l1} ${l2} ${c1} >"${_polap_var_compare_mtdna_compare}"

  _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
  # Disable debugging if previously enabled
  [ "$_POLAP_DEBUG" -eq 1 ] && set +x
  return 0
}

source "${_POLAPLIB_DIR}/run-polap-function-utilities.sh"

################################################################################
# Downloads the mtDNA sequence for a given species name.
################################################################################
function _run_polap_get-mtdna {
  # Enable debugging if _POLAP_DEBUG is set
  [ "$_POLAP_DEBUG" -eq 1 ] && set -x
  _polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

  # Set verbosity level: stderr if verbose >= 2, otherwise discard output
  local _polap_output_dest="/dev/null"
  [ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

  source "${_POLAPLIB_DIR}/polap-variables-common.sh"

  if ! run_check_ncbitools; then
    error_polap_conda
    exit $EXIT_ERROR
  fi

  # Help message
  local help_message=$(
    cat <<HEREDOC
# Download the organelle-genome sequence in FASTA format from NCBI.
#
# Note: A downloaded fasta file may have multiple sequences.
#
# 1. NCBI search for a mitochondrial genome:
#   "(mitochondrion[Title] AND complete[Title] AND genome[Title]) AND <species>[Organism]"
# 2. NCBI search for a plastid genome:
#   "(chloroplast[Title] AND complete[Title] AND genome[Title]) AND <species>[Organism]"
#
# Arguments:
#   --species "scientific name" (highest priority)
#   or
#   -o ${_arg_outdir} (next priority)
# Outputs:
#   ${_polap_var_project_mtdna_fasta1}
#   ${_polap_var_project_mtdna_fasta2}
#   ${_polap_var_project_mtdna_fasta2_accession}
# Preconditions:
#   get-bioproject
Example: $(basename "$0") ${_arg_menu[0]} --species "Anthoceros agrestis"
Example: $(basename "$0") ${_arg_menu[0]} -o o
Example: $(basename "$0") ${_arg_menu[0]} view
HEREDOC
  )

  # Display help message
  [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

  if [[ ${_arg_menu[1]} == "file" ]]; then
    if [[ -s "${_polap_var_project_species}" ]]; then
      _polap_log0_file "${_polap_var_project_species}"
      _polap_log0_cat "${_polap_var_project_species}"
    else
      _polap_log0 "No species file."
    fi
    exit $EXIT_SUCCESS
  fi

  if [[ ${_arg_menu[1]} == "view" ]]; then
    if [[ -s "${_polap_var_project_mtdna_fasta2_accession}" ]]; then
      _polap_log0_cat "${_polap_var_project_mtdna_fasta2_accession}"
      if [ -s "${_polap_var_project_mtdna_fasta1}" ]; then
        seqkit fx2tab -n -l "${_polap_var_project_mtdna_fasta1}" >&3
        seqkit stats "${_polap_var_project_mtdna_fasta1}" >&3
      else
        _polap_log1 "No such file: ${_polap_var_project_mtdna_fasta1}"
        _polap_log0 "No organelle genome sequence"
      fi
      if [ -s "${_polap_var_project_mtdna_fasta2}" ]; then
        seqkit stats "${_polap_var_project_mtdna_fasta2}" >&3
      else
        _polap_log1 "No such file: ${_polap_var_project_mtdna_fasta2}"
      fi

    else
      _polap_log0 "No result yet."
    fi
    exit $EXIT_SUCCESS
  fi

  if [[ "${_arg_plastid}" == "off" ]]; then
    _polap_log0 "getting the mitochondrial organelle-genome sequence of a plant species ..."
  else
    _polap_log0 "getting the plastid organelle-genome sequence of a plant species ..."
  fi

  # Determine species name
  _polap_log1 "  step 1: determine the species name"
  local SPECIES=""
  if [ -n "${_arg_species}" ]; then
    _polap_log1 "  option --species: ${_arg_species}"
    SPECIES="${_arg_species}"
    _polap_log3_cmd mkdir -p "${_polap_var_project}"
  elif [ -s "${_polap_var_project_species}" ]; then
    if [[ -s "${_polap_var_project_txt}" ]]; then
      local bioproject_id=$(<"${_polap_var_project_txt}")
      local n=$(wc -l <"${_polap_var_project_species}")
      if [[ "${n}" -gt 1 ]]; then
        _polap_log0_file "${_polap_var_project_species}"
        _polap_log0_cat "${_polap_var_project_species}"
        die "ERROR: you have multiple species names in BioProject: ${bioproject_id}"
      fi
      SPECIES=$(<"${_polap_var_project_species}")
      _polap_log0 "  bioproject's species: $SPECIES"
    else
      die "ERROR: you have not run the menu: get-bioproject on the output folder [${_arg_outdir}], yet!"
    fi
  elif [ -n "${_arg_bioproject}" ]; then
    _run_polap_get-bioproject
    if [ -s "${_polap_var_project_species}" ]; then
      local bioproject_id=$(<"${_polap_var_project_txt}")
      local n=$(wc -l <"${_polap_var_project_species}")
      if [[ "${n}" -gt 1 ]]; then
        _polap_log0_file "${_polap_var_project_species}"
        _polap_log0_cat "${_polap_var_project_species}"
        die "ERROR: you have multiple species names in BioProject: ${bioproject_id}"
      fi
      SPECIES=$(<"${_polap_var_project_species}")
      _polap_log1 "  bioproject's species: $SPECIES"
    else
      die "  no species name is provided in the BioProject: ${_arg_bioproject}"
    fi
  else
    die "  use --species or --bioproject option."
  fi

  # Download the mitochondrial genome sequence for the given species
  _polap_log1 "  step 2: download the mitochondrial genome sequence for the given species"
  if [ "${_arg_plastid}" = "off" ]; then
    _polap_log1 "  downloading mitochondrial complete genomes of ${SPECIES} ..."
    esearch \
      -db nuccore \
      -query "(mitochondrion[Title] AND complete[Title] AND genome[Title]) AND ${SPECIES}[Organism]" |
      efetch -format fasta >"${_polap_var_project_mtdna_fasta1}"
  else
    _polap_log1 "  downloading chloroplast complete genomes of ${SPECIES} ..."
    esearch \
      -db nuccore \
      -query "(chloroplast[Title] AND complete[Title] AND genome[Title]) AND ${SPECIES}[Organism]" |
      efetch -format fasta >"${_polap_var_project_mtdna_fasta1}"
    if [[ ! -s "${_polap_var_project_mtdna_fasta1}" ]]; then
      esearch \
        -db nuccore \
        -query "(plastid[Title] AND complete[Title] AND genome[Title]) AND ${SPECIES}[Organism]" |
        efetch -format fasta >"${_polap_var_project_mtdna_fasta1}"
    fi
  fi
  _polap_log2_file "${_polap_var_project_mtdna_fasta1}"

  # Check if the fasta file was successfully downloaded
  _polap_log1 "  step 3: check if the fasta file was successfully downloaded"
  if [ -s "${_polap_var_project_mtdna_fasta1}" ]; then
    seqkit fx2tab --length --name --header-line \
      "${_polap_var_project_mtdna_fasta1}" \
      >"${_polap_var_project_mtdna_fasta1_stats}"

    local n=$(seqkit stats -T "${_polap_var_project_mtdna_fasta1}" |
      csvtk cut -t -f num_seqs |
      csvtk del-header)

    if [[ "${n}" -gt 1 ]]; then
      _polap_log2 "  multiple (${n}) accession numbers for the organelle genome sequences"
    fi

    seqkit head -n 1 "${_polap_var_project_mtdna_fasta1}" \
      -o "${_polap_var_project_mtdna_fasta2}"

    seqkit seq -ni "${_polap_var_project_mtdna_fasta2}" \
      -o "${_polap_var_project_mtdna_fasta2_accession}"

    local accession=$(<"${_polap_var_project_mtdna_fasta2_accession}")
    _polap_log2_cat "${_polap_var_project_mtdna_fasta1_stats}"
    _polap_log0 "Species: ${SPECIES}"
    _polap_log0 "NCBI accession: ${accession}"
  else
    echo "no mtDNA" >"${_polap_var_project_mtdna_fasta2_accession}"
    _polap_log0 "No mtDNA sequence found for the species: ${SPECIES}"
    rm -f \
      "${_polap_var_project_mtdna_fasta1_stats}" \
      "${_polap_var_project_mtdna_fasta1}" \
      "${_polap_var_project_mtdna_fasta2}"
  fi

  _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
  # Disable debugging if previously enabled
  [ "$_POLAP_DEBUG" -eq 1 ] && set +x
  return 0
}

function _run_polap_get-dna-by-accession {
  # Enable debugging if _POLAP_DEBUG is set
  [ "$_POLAP_DEBUG" -eq 1 ] && set -x
  _polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

  # Set verbosity level: stderr if verbose >= 2, otherwise discard output
  local _polap_output_dest="/dev/null"
  [ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

  source "${_POLAPLIB_DIR}/polap-variables-common.sh"

  if ! run_check_ncbitools; then
    error_polap_conda
    exit $EXIT_ERROR
  fi

  # Help message
  local help_message=$(
    cat <<HEREDOC
Obtain the DNA sequence by downloading it from NCBI in FASTA format 
using its corresponding accession ID.
Example: $(basename "$0") ${_arg_menu[0]} <accession> <output_file>
Example: $(basename "$0") ${_arg_menu[0]} NC_027276 ptdna-Podococcus_acaulis.fa
HEREDOC
  )

  # Display help message
  [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

  _polap_log0 "getting the DNA sequence using accesion ${_arg_menu[1]}"

  # Download the mitochondrial genome sequence for the given species
  esearch \
    -db nuccore \
    -query "${_arg_menu[1]}[ACCN]" |
    efetch -format fasta >"${_arg_menu[2]}"

  _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
  # Disable debugging if previously enabled
  [ "$_POLAP_DEBUG" -eq 1 ] && set +x
  return 0
}

################################################################################
# Gene table for importing to the Chloroplot R package.
################################################################################
function _run_polap_plot-mtdna {
  # Enable debugging if _POLAP_DEBUG is set
  [ "$_POLAP_DEBUG" -eq 1 ] && set -x

  # Set verbosity level: stderr if verbose >= 2, otherwise discard output
  local _polap_output_dest="/dev/null"
  [ "${_arg_verbose}" -ge 2 ] && _polap_output_dest="/dev/stderr"

  # Grouped file path declarations
  source "${_POLAPLIB_DIR}/polap-variables-common.sh"

  # Set working directories and file paths
  local _polap_var_assembly="${_polap_var_ga}"
  local _polap_var_chloroplot="${_polap_var_assembly}/60-chloroplot"

  # Help message
  local help_message=$(
    cat <<HEREDOC
# Counts genes annotated on a genome assembly and plots the mtDNA genome.
# Arguments:
#   -f ${_arg_final_assembly}
#   -i ${_arg_inum}: a Flye genome assembly number
# Inputs:
#   ${_arg_final_assembly}: mtDNA sequence
# Outputs:
Example: $(basename "$0") ${_arg_menu[0]} [-i|--inum <arg>] [-f ${_arg_final_assembly}]
HEREDOC
  )

  # Display help message
  [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

  _polap_log1 "LOG: Plotting mitochondrial DNA genome for Flye assembly ${_arg_inum}..."

  # Run the R script to generate the mtDNA plot
  "${_POLAPLIB_DIR}/run-polap-r-plot-mtdna.R " \
    "${_polap_var_chloroplot}/annotation.bed" \
    "${FA}" \
    2>"$_polap_output_dest"

  # Output file information
  _polap_log0 "FILE: mt.3.pdf has been created."

  _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
  # Disable debugging if previously enabled
  [ "$_POLAP_DEBUG" -eq 1 ] && set +x
  return 0
}

function _run_polap_mtdna-edges {
  # Enable debugging if _POLAP_DEBUG is set
  [ "$_POLAP_DEBUG" -eq 1 ] && set -x
  _polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

  # Set verbosity level: stderr if verbose >= 2, otherwise discard output
  local _polap_output_dest="/dev/null"
  [ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

  # Set variables
  source "${_POLAPLIB_DIR}/polap-variables-common.sh"

  # Help message
  local help_message=$(
    cat <<HEREDOC
# Selects mtDNA sequences from a GFA.
#
# Arguments:
#   -j ${_arg_jnum}: organelle assembly number
#
# Inputs:
#   ${_polap_var_oga_assembly_graph_gfa}
#
# Outputs:
#   ${_polap_var_oga_mt_fasta}
#
Example: $(basename $0) ${_arg_menu[0]} -j <arg>
HEREDOC
  )

  # Display help message
  [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

  # Create necessary directories
  rm -rf "${_polap_var_mtdna}"
  mkdir -p "${_polap_var_mtdna}"

  # Prompt the user for input
  read -p "Enter edges using Bandage (e.g., edge_265, edge_520, edge_425): " edges

  # Convert the input string into new lines (replace ", " with newlines)
  local formatted_edges=$(echo "$edges" | tr ', ' '\n' | sed '/^ *$/d')
  echo "${formatted_edges}" >"${_polap_var_mtdna_names}"

  _polap_log1_file "Input GFA: ${_polap_var_oga_assembly_graph_gfa}"
  # _polap_log1_file "Input Annotation Table: ${_polap_var_ga_annotation_all}"
  #
  _polap_log3_pipe "gfatools view \
		-l @${_polap_var_mtdna_names} \
		${_polap_var_oga_assembly_graph_gfa} \
		2>$_polap_output_dest \
		>${_polap_var_mtdna_sub_gfa}"

  _polap_log0 "subgfa: ${_polap_var_mtdna_sub_gfa}"

  # Convert GFA to FASTA
  _polap_log3_pipe "gfatools gfa2fa ${_polap_var_mtdna_sub_gfa} \
		>${_polap_var_mtdna_fasta} \
		2>$_polap_output_dest"

  cp ${_polap_var_mtdna_fasta} ${_polap_var_oga_mt_fasta}

  _polap_log0 "FASTA of the input GFA: ${_polap_var_oga_mt_fasta}"

  _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
  # Disable debugging if previously enabled
  [ "$_POLAP_DEBUG" -eq 1 ] && set +x
  return 0
}

################################################################################
# Selects mtDNA sequences from a GFA.
#
# annotation: uses 30-contigger/graph_final.gfa
#
# extraction: uses assembly_graph.gfa
################################################################################
function _run_polap_select-mtdna {
  # Enable debugging if _POLAP_DEBUG is set
  [ "$_POLAP_DEBUG" -eq 1 ] && set -x
  _polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

  # Set verbosity level: stderr if verbose >= 2, otherwise discard output
  local _polap_output_dest="/dev/null"
  [ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

  # Set variables
  # CHECK: local function
  source "${_POLAPLIB_DIR}/polap-variables-common.sh"

  # local FDIR="${_arg_outdir}/${_arg_inum}"
  # local _polap_var_mtdna="$FDIR/mtdna"
  #
  # # File paths
  # local _polap_var_assembly_graph_gfa="$FDIR/assembly_graph.gfa"
  # local _polap_var_ga_annotation_all="$FDIR/assembly_info_organelle_annotation_count-all.txt"
  # local _polap_var_mt_fasta="$FDIR/mt.0.fasta"
  # local _polap_var_mt_edges="$FDIR/mt.0.edges"
  # local _polap_var_mtdna_1_gfa_all="${_polap_var_mtdna}/1-gfa.all.gfa"
  # local _polap_var_mtdna_gfa_links="${_polap_var_mtdna}/1-gfa.links.tsv"
  # local _polap_var_mtdna_gfa_links_edges="${_polap_var_mtdna}/1-gfa.links.edges.txt"
  # local _polap_var_mtdna_gfa_links_circular_path="${_polap_var_mtdna}/2-gfa.links.circular.path.txt"
  # local _polap_var_mtdna_circular_path="${_polap_var_mtdna}/3-circular.path.txt"
  # local _polap_var_mtdna_edge_fasta="${_polap_var_mtdna}/5-edge.fasta"
  # local _polap_var_mtdna_fasta="${_polap_var_mtdna}/4-gfa.fasta"

  # Help message
  local help_message=$(
    cat <<HEREDOC
# Selects mtDNA sequences from a GFA.
#
# Arguments:
#   -j ${_arg_jnum}: organelle assembly number
#
# Inputs:
#   ${_polap_var_oga_assembly_graph_gfa}
#
# Outputs:
#   ${_polap_var_oga_mt_fasta}
#
Example: $(basename $0) ${_arg_menu[0]} -j <arg>
HEREDOC
  )

  # Display help message
  [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

  if [[ "${_arg_menu[1]}" = "bandage" ]]; then
    # Prompt the user for input
    read -p "Enter edges using Bandage (e.g., edge_265, edge_520, edge_425): " edges

    # Convert the input string into new lines (replace ", " with newlines)
    local formatted_edges=$(echo "$edges" | tr ', ' '\n' | sed '/^ *$/d')
    echo "${formatted_edges}" >"${_polap_var_mtdna_names}"

    return $RETURN_SUCCESS
  fi

  # Check for required files
  # check_file_existence "${_polap_var_assembly_graph_gfa}"
  # check_file_existence "${_polap_var_ga_annotation_all}"

  # Create necessary directories
  # _polap_log3_cmd rm -rf "${_polap_var_mtdna}"
  # _polap_log3_cmd mkdir -p "${_polap_var_mtdna}"

  _polap_log1_file "Input GFA: ${_polap_var_oga_assembly_graph_gfa}"
  # _polap_log1_file "Input Annotation Table: ${_polap_var_ga_annotation_all}"
  #
  _polap_log3_pipe "gfatools view \
		-l @${_polap_var_mtdna_names} \
		${_polap_var_oga_assembly_graph_gfa} \
		2>$_polap_output_dest \
		>${_polap_var_mtdna_sub_gfa}"

  _polap_log0 "subgfa: ${_polap_var_mtdna_sub_gfa}"

  # Convert GFA to FASTA
  _polap_log3_pipe "gfatools gfa2fa ${_polap_var_mtdna_sub_gfa} \
		>${_polap_var_mtdna_fasta} \
		2>$_polap_output_dest"

  _polap_log2_file "FASTA of the input GFA: ${_polap_var_mtdna_fasta}"

  return

  # Step 1: List connected components (GFA processing)
  _polap_log3_pipe "gfatools view -S ${_polap_var_assembly_graph_gfa} \
		>${_polap_var_mtdna_1_gfa_all} \
		2>$_polap_output_dest"

  _polap_log2_file "GFA content: ${_polap_var_mtdna_1_gfa_all}"

  # Extract links from GFA and store them
  grep "^L" "${_polap_var_mtdna_1_gfa_all}" >"${_polap_var_mtdna_gfa_links}"

  local _polap_var_number_links_gfa=$(wc -l <"${_polap_var_mtdna_gfa_links}")
  if [ "${_polap_var_number_links_gfa}" -eq 0 ]; then
    echoerr "LOG: No circular sequences are found."
    return
  fi

  _polap_log2_file "GFA Links: ${_polap_var_mtdna_gfa_links}"

  # Process the GFA links
  "${_POLAPLIB_DIR}"/run-polap-r-select-mtdna-1-nx-gfa-links.R \
    "${_polap_var_mtdna_gfa_links}" \
    "${_polap_var_mtdna_gfa_links_edges}" \
    2>"$_polap_output_dest"

  _polap_log2_file "GFA Links Edges: ${_polap_var_mtdna_gfa_links_edges}"

  # Step 2: Find circular paths using Python script
  python "${_POLAPLIB_DIR}"/run-polap-py-select-mtdna-2-nx-simple-cycles.py \
    "${_polap_var_mtdna_gfa_links_edges}" \
    "${_polap_var_mtdna_gfa_links_circular_path}" \
    2>"$_polap_output_dest"
  # python "${_POLAPLIB_DIR}"/run-polap-py-select-mtdna-2-nx-find-circular-path.py \
  # 	"${_polap_var_mtdna_gfa_links_edges}" \
  # 	"${_polap_var_mtdna_gfa_links_circular_path}" \
  # 	2>"$_polap_output_dest"

  _polap_log2_file "Circular Path: ${_polap_var_mtdna_gfa_links_circular_path}"

  # return
  #
  # Process circular path
  # process_circular_path "${_polap_var_mtdna_gfa_links_circular_path}" "${_polap_var_mtdna_circular_path}"
  # _polap_log2_file "Processed Circular Path: ${_polap_var_mtdna_circular_path}"

  # Step 3: Concatenate sequences from circular path
  for _polap_file in "${_polap_var_mtdna_gfa_links_circular_path}"_*.tsv; do
    # Extract the number from the filename using parameter expansion
    local index=$(echo "${_polap_file}" | sed -E 's/.*_(.*)\.tsv/\1/')
    # _polap_log0 "${index}"
    concatenate_sequences "${_polap_file}" "${_polap_var_mtdna_fasta}" "${_polap_var_mtdna_edge_fasta}-${index}"
    finalize_mtdna_fasta "${_polap_var_mtdna_edge_fasta}-${index}" "${_polap_var_mt_fasta}-${index}.fasta"
    _polap_log2_file "Final mtDNA FASTA: ${_polap_var_mt_fasta}-${index}.fasta"
    cp "${_polap_file}" "${_polap_var_mt_edges}-${index}.edges"
    seqkit stats "${_polap_var_mt_fasta}-${index}.fasta" 1>&2
  done

  largest_file=$(ls -S "${_polap_var_mt_fasta}"-*.fasta | head -n 1)
  local index=$(echo "${largest_file}" | sed -E 's/.*-(.*)\.fasta/\1/')
  cp "${largest_file}" "${_polap_var_mt_fasta}"
  cp "${_polap_var_mt_edges}-${index}.edges" "${_polap_var_mt_edges}"

  # concatenate_sequences "${_polap_var_mtdna_circular_path}" "${_polap_var_mtdna_fasta}" "${_polap_var_mtdna_edge_fasta}"

  # Finalize the mtDNA FASTA
  # finalize_mtdna_fasta "${_polap_var_mtdna_edge_fasta}" "${_polap_var_mt_fasta}"
  # _polap_log2_file "Final mtDNA FASTA: ${_polap_var_mt_fasta}"

  # cp "${_polap_var_mtdna_circular_path}" "${_polap_var_mt_edges}"

  _polap_log1_file "Output mtDNA FASTA: ${_polap_var_mt_fasta}"

  _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
  # Disable debugging if previously enabled
  [ "$_POLAP_DEBUG" -eq 1 ] && set +x
  return 0
}

################################################################################
# Before 2024-09-27
################################################################################
function _run_polap_select-mtdna-org {
  # Enable debugging if _POLAP_DEBUG is set
  [ "$_POLAP_DEBUG" -eq 1 ] && set -x
  _polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

  # Set verbosity level: stderr if verbose >= 2, otherwise discard output
  local _polap_output_dest="/dev/null"
  [ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

  # Set variables
  # CHECK: local function
  source "${_POLAPLIB_DIR}/polap-variables-common.sh"

  # local FDIR="${_arg_outdir}/${_arg_inum}"
  # local _polap_var_mtdna="$FDIR/mtdna"
  #
  # # File paths
  # local _polap_var_assembly_graph_gfa="$FDIR/assembly_graph.gfa"
  # local _polap_var_ga_annotation_all="$FDIR/assembly_info_organelle_annotation_count-all.txt"
  # local _polap_var_mt_fasta="$FDIR/mt.0.fasta"
  # local _polap_var_mt_edges="$FDIR/mt.0.edges"
  # local _polap_var_mtdna_1_gfa_all="${_polap_var_mtdna}/1-gfa.all.gfa"
  # local _polap_var_mtdna_gfa_links="${_polap_var_mtdna}/1-gfa.links.tsv"
  # local _polap_var_mtdna_gfa_links_edges="${_polap_var_mtdna}/1-gfa.links.edges.txt"
  # local _polap_var_mtdna_gfa_links_circular_path="${_polap_var_mtdna}/2-gfa.links.circular.path.txt"
  # local _polap_var_mtdna_circular_path="${_polap_var_mtdna}/3-circular.path.txt"
  # local _polap_var_mtdna_edge_fasta="${_polap_var_mtdna}/5-edge.fasta"
  # local _polap_var_mtdna_fasta="${_polap_var_mtdna}/4-gfa.fasta"

  # Help message
  local help_message=$(
    cat <<HEREDOC
# Selects mtDNA sequences from a GFA.
#
# Arguments:
#   -i ${_arg_inum}: organelle assembly number
#
# Inputs:
#   ${_polap_var_assembly_graph_gfa}
#   ${_polap_var_ga_annotation_all}
#
# Outputs:
#   ${_polap_var_mt_fasta}
#
Example: $(basename $0) ${_arg_menu[0]} [-i|--inum <arg>]
HEREDOC
  )

  # Display help message
  [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && echo "${help_message}" >&2 && exit $EXIT_SUCCESS

  # Check for required files
  check_file_existence "${_polap_var_assembly_graph_gfa}"
  check_file_existence "${_polap_var_ga_annotation_all}"

  # Create necessary directories
  rm -rf "${_polap_var_mtdna}"
  mkdir -p "${_polap_var_mtdna}"

  _polap_log1_file "Input GFA: ${_polap_var_assembly_graph_gfa}"
  _polap_log1_file "Input Annotation Table: ${_polap_var_ga_annotation_all}"

  # Convert GFA to FASTA
  gfatools gfa2fa "${_polap_var_assembly_graph_gfa}" \
    >"${_polap_var_mtdna_fasta}" \
    2>"$_polap_output_dest"

  _polap_log2_file "FASTA of the input GFA: ${_polap_var_mtdna_fasta}"

  # Step 1: List connected components (GFA processing)
  gfatools view -S "${_polap_var_assembly_graph_gfa}" \
    >"${_polap_var_mtdna_1_gfa_all}" \
    2>"$_polap_output_dest"

  _polap_log2_file "GFA content: ${_polap_var_mtdna_1_gfa_all}"

  # Extract links from GFA and store them
  grep "^L" "${_polap_var_mtdna_1_gfa_all}" >"${_polap_var_mtdna_gfa_links}"

  local _polap_var_number_links_gfa=$(wc -l <"${_polap_var_mtdna_gfa_links}")
  if [ "${_polap_var_number_links_gfa}" -eq 0 ]; then
    echoerr "LOG: No circular sequences are found."
    return
  fi

  _polap_log2_file "GFA Links: ${_polap_var_mtdna_gfa_links}"

  # Process the GFA links
  "${_POLAPLIB_DIR}"/run-polap-r-select-mtdna-1-nx-gfa-links.R \
    "${_polap_var_mtdna_gfa_links}" \
    "${_polap_var_mtdna_gfa_links_edges}" \
    2>"$_polap_output_dest"

  _polap_log2_file "GFA Links Edges: ${_polap_var_mtdna_gfa_links_edges}"

  # Step 2: Find circular paths using Python script
  # python "${_POLAPLIB_DIR}"/run-polap-py-select-mtdna-2-nx-simple-cycles.py \
  # 	"${_polap_var_mtdna_gfa_links_edges}" \
  # 	"${_polap_var_mtdna_gfa_links_circular_path}" \
  # 	2>"$_polap_output_dest"
  python "${_POLAPLIB_DIR}"/run-polap-py-select-mtdna-2-nx-find-circular-path.py \
    "${_polap_var_mtdna_gfa_links_edges}" \
    "${_polap_var_mtdna_gfa_links_circular_path}" \
    2>"$_polap_output_dest"

  _polap_log2_file "Circular Path: ${_polap_var_mtdna_gfa_links_circular_path}"

  # Process circular path
  process_circular_path "${_polap_var_mtdna_gfa_links_circular_path}" "${_polap_var_mtdna_circular_path}"
  _polap_log2_file "Processed Circular Path: ${_polap_var_mtdna_circular_path}"

  # Step 3: Concatenate sequences from circular path
  concatenate_sequences "${_polap_var_mtdna_circular_path}" "${_polap_var_mtdna_fasta}" "${_polap_var_mtdna_edge_fasta}"

  # Finalize the mtDNA FASTA
  finalize_mtdna_fasta "${_polap_var_mtdna_edge_fasta}" "${_polap_var_mt_fasta}"
  _polap_log2_file "Final mtDNA FASTA: ${_polap_var_mt_fasta}"

  cp "${_polap_var_mtdna_circular_path}" "${_polap_var_mt_edges}"

  _polap_log1_file "Output mtDNA FASTA: ${_polap_var_mt_fasta}"
  seqkit stats "${_polap_var_mt_fasta}" 1>&2

  _polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
  # Disable debugging if previously enabled
  [ "$_POLAP_DEBUG" -eq 1 ] && set +x
  return 0
}

# Helper function to process the circular path and format the edge data
process_circular_path() {
  local input_file=$1
  local output_file=$2

  tail -n +2 "${input_file}" |
    while IFS=$'\t' read -r col1 col2; do
      local number1=$(echo "$col1" | grep -o '^[0-9]*')
      local sign1=$(echo "$col1" | grep -o '[+-]$')
      local new_col1="edge_${number1}\t${sign1}"
      echo -e "$new_col1"
    done >"${output_file}"
}

# Helper function to concatenate sequences based on circular path
concatenate_sequences() {
  local circular_path_file=$1
  local fasta_file=$2
  local output_fasta=$3

  # Create an empty file for the concatenated sequence
  >"${output_fasta}"

  # Loop through the circular path IDs and concatenate sequences
  while read -r id strand; do
    echoerr "Processing ID: $id with strand: $strand"
    if [[ "$strand" == "+" ]]; then
      seqkit grep -p "$id" "${fasta_file}" | seqkit seq -t dna -v >>"${output_fasta}"
    elif [[ "$strand" == "-" ]]; then
      seqkit grep -p "$id" "${fasta_file}" | seqkit seq -t dna -v -r -p >>"${output_fasta}"
    fi
  done <"${circular_path_file}"
}

# Helper function to finalize mtDNA FASTA by concatenating the sequences into one entry
finalize_mtdna_fasta() {
  local input_fasta=$1
  local output_fasta=$2

  echo ">concatenated_sequence" >"${output_fasta}"
  seqkit fx2tab "${input_fasta}" | cut -f2 | tr -d '\n' >>"${output_fasta}"
}
