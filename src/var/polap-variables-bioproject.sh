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

source "$script_dir/polap-variables-common.sh"

local _polap_var_output="${ODIR}"
local _polap_var_bioproject_txt="${_polap_var_output}/bioproject.txt"
local _polap_var_bioproject_runinfo_all="${_polap_var_output}/bioproject.runinfo"
local _polap_var_bioproject="${_polap_var_output}/0-bioproject"
local _polap_var_bioproject_runinfo="${_polap_var_bioproject}/1-runinfo.tsv"
local _polap_var_bioproject_sra_per_species="${_polap_var_bioproject}/1-runinfo.per.species.tsv"
local _polap_var_bioproject_sra_long_read="${_polap_var_bioproject}/1-sra-long-read.tsv"
local _polap_var_bioproject_sra_short_read="${_polap_var_bioproject}/1-sra-short-read.tsv"
local _polap_var_bioproject_species="${_polap_var_bioproject}/1-species.txt"
local _polap_var_bioproject_taxon_id="${_polap_var_bioproject}/1-taxon-id.txt"
local _polap_var_bioproject_taxonomy="${_polap_var_bioproject}/1-taxonomy.txt"
local _polap_var_bioproject_passed="${_polap_var_bioproject}/1-passed.txt"

# local SRA=$(cut -f1 "${_polap_var_bioproject_sra_long_read}")
# local _polap_var_base_sra_long_fastq="${_polap_var_output}/${SRA}.fastq"

local _polap_var_bioproject_mtdna_fasta1="${_polap_var_bioproject}/1-mtdna.fasta"
local _polap_var_bioproject_mtdna_fasta1_stats="${_polap_var_bioproject}/1-mtdna.fasta.stats"
local _polap_var_bioproject_mtdna_fasta2="${_polap_var_bioproject}/2-mtdna.fasta"
local _polap_var_bioproject_mtdna_fasta2_accession="${_polap_var_bioproject}/2-mtdna.accession"

local _polap_var_bioproject_blastn1="${_polap_var_bioproject}/3-blastn1.txt"
local _polap_var_bioproject_blastn2="${_polap_var_bioproject}/3-blastn2.txt"
local _polap_var_bioproject_blastn3="${_polap_var_bioproject}/3-blastn3.txt"
local _polap_var_bioproject_blastn3_length="${_polap_var_bioproject}/3-blastn3.length.txt"
