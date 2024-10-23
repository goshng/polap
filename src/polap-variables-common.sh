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

# Global option variables
#
# ODIR
# INUM
# JNUM

# Base variables used across multiple scripts
local _polap_var_output="${ODIR}"
local _polap_var_base="${_polap_var_output}"
local _polap_var_bioproject="${_polap_var_output}/0-bioproject"
local _polap_var_wga="${_polap_var_output}/0"
local _polap_var_ga="${_polap_var_output}/${INUM}"
local _polap_var_contigger="${_polap_var_ga}/30-contigger"
local _polap_var_ann="${_polap_var_ga}/50-annotation"
local _polap_var_mtcontigs="${_polap_var_ga}/${JNUM}/mtcontigs"
# local _polap_var_mtcontigs="${_polap_var_ga}/60-mtcontigs"
local _polap_var_mtdna="${_polap_var_ga}/70-mtdna"
local _polap_var_compare="${_polap_var_ga}/80-compare"
local MTCONTIGNAME="${_polap_var_ga}/mt.contig.name-${JNUM}"
local _polap_var_mtcontigname="${_polap_var_ga}/mt.contig.name-${JNUM}"
local _polap_var_links="${_polap_var_mtcontigs}/4-gfa.links"
local _polap_var_oga="${_polap_var_output}/${JNUM}"
local _polap_var_oga_contigger="${_polap_var_oga}/30-contigger"
# local _polap_var_oga="${_polap_var_ga}"
local _polap_var_seeds="${_polap_var_oga}/seeds"