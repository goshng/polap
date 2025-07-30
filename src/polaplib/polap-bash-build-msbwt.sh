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
# A bash script used in polap-cmd-disassemble.sh
# It is used to prepare a MSBWT file before polishing.
# command time -v bash "${_POLAPLIB_DIR}"/polap-bash-build-msbwt.sh
# command time -v fmlrc -p "${_arg_threads}"

_short_read1="${1}"
_msbwt_dir="${2}"

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate polap-fmlrc

rm -rf "${_msbwt_dir}"
cat "${_short_read1}" |
	awk 'NR % 4 == 2' | sort | tr NT TN |
	ropebwt2 -LR |
	tr NT TN |
	msbwt convert "${_msbwt_dir}" \
		>/dev/null

conda deactivate
