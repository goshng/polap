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

# variables
INUM=${_arg_inum}
JNUM=${_arg_jnum}
if [ "${_arg_archive_is}" = "off" ]; then
	_arg_archive="${ODIR}-a"
fi

# tuning variables for optimal performance
# COV=${_arg_coverage} -> --asm-coverage option meaning
if test -z "$DEBUG"; then
	DEBUG=0
fi

SECONDS=0

if [ "${_arg_log_is}" = "off" ]; then
	LOG_FILE="${ODIR}/${_arg_log}"
	[[ ! -d "${ODIR}" ]] && mkdir -p "${ODIR}"
else
	LOG_FILE="${_arg_log}"
fi

################################################################################
# for the magic logit function
################################################################################
function logit() {
	while read; do
		# echo "$(date) $REPLY" >> ${LOG_FILE}
		# -1 means "current time"
		# printf "[%(%Y-%m-%d %T)T] %s\n" -1 "$REPLY" >> ${LOG_FILE}
		printf "[%s] %s\n" "$(date +"%Y-%m-%d %T")" "$REPLY" >>${LOG_FILE}
	done
}
