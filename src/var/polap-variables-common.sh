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

# wga -> annotation, wga
# oga ->

# Base variables used across multiple scripts
local _polap_var_output="${ODIR}"
local _polap_var_base="${ODIR}"
local _polap_var_wga="${ODIR}/0"
local _polap_var_oga="${ODIR}/${INUM}"
# local _polap_var_ga="${ODIR}/${INUM}"

local _polap_var_wga_contigger="${_polap_var_wga}/30-contigger"
local _polap_var_wga_contigger_gfa="${_polap_var_wga_contigger}/graph_final.gfa"
local _polap_var_wga_annotation="${_polap_var_wga}/assembly_info_organelle_annotation_count-all.txt"

# Add other common variables here as needed
