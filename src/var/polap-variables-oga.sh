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

local _polap_var_oga="${ODIR}/${INUM}"
local _polap_var_mtdna1="${_polap_var_oga}/mt.1.fa"
local _polap_var_mtdna2="${_polap_var_oga}/mt.2.fa"
local _polap_var_mtdna3="${_polap_var_oga}/mt.3.fa"
local _polap_var_mtdna_compare="${_polap_var_oga}/mt.compare.txt"
local _polap_var_oga_blastn1="${_polap_var_oga}/3-blastn1.txt"
local _polap_var_oga_blastn2="${_polap_var_oga}/3-blastn2.txt"
local _polap_var_oga_blastn3="${_polap_var_oga}/3-blastn3.txt"
local _polap_var_oga_blastn3_length="${_polap_var_oga}/3-blastn3.length.txt"
