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

local _polap_var_base_jellyfish_out="${_polap_var_output}/jellyfish_out"
local _polap_var_base_jellyfish_out_histo="${_polap_var_output}/jellyfish_out.histo"
local _polap_var_base_genome_size="${_polap_var_output}/short_expected_genome_size.txt"
local _polap_var_base_long_total_length="${_polap_var_output}/long_total_length.txt"
local _polap_var_base_nk_fq_gz="${_polap_var_output}/nk.fq.gz"
local _polap_var_base_l_fq_gz="${_polap_var_output}/l.fq.gz"
local _polap_var_base_nk_fq_stats="${_polap_var_output}/nk.fq.stats"
local _polap_var_base_fq_stats="${_polap_var_output}/fq.stats"
local _polap_var_base_msbwt="${_polap_var_output}/msbwt/comp_msbwt.npy"
local _polap_var_base_msbwt_tar_gz="${_polap_var_output}/msbwt.tar.gz"
local _polap_var_bioproject_txt="${_polap_var_output}/bioproject.txt"
local _polap_var_bioproject_runinfo_all="${_polap_var_output}/bioproject.runinfo"
