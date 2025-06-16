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
# This script is included in the early part of the polap.sh script.
# It initializes the output folder.
# It also defines the log function.
# Put something here that must be initialized early on.
################################################################################

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  echo "[ERROR] This script must be sourced, not executed: use 'source $BASH_SOURCE'" >&2
  return 1 2>/dev/null || exit 1
fi
: "${_POLAP_DEBUG:=0}"
: "${_POLAP_RELEASE:=0}"

# variables
if [ "${_arg_archive_is}" = "off" ]; then
  _arg_archive="${_arg_outdir}-a"
fi

# tuning variables for optimal performance
# COV=${_arg_coverage} -> --asm-coverage option meaning
# if test -z "$_POLAP_DEBUG"; then
#   _POLAP_DEBUG=0
# fi

SECONDS=0

[[ ! -d "${_arg_outdir}" ]] && mkdir -p "${_arg_outdir}"
[[ ! -d "${_arg_outdir}/tmp" ]] && mkdir -p "${_arg_outdir}/tmp"
[[ ! -d "${_arg_outdir}/log" ]] && mkdir -p "${_arg_outdir}/log"
if [ "${_arg_log_is}" = "off" ]; then
  LOG_FILE="${_arg_outdir}/${_arg_log}"
else
  LOG_FILE="${_arg_log}"
fi

################################################################################
# for the magic logit function
################################################################################
function logit {
  while read; do
    # echo "$(date) $REPLY" >> ${LOG_FILE}
    # -1 means "current time"
    # printf "[%(%Y-%m-%d %T)T] %s\n" -1 "$REPLY" >> ${LOG_FILE}
    printf "[%s] %s\n" "$(date +"%Y-%m-%d %T")" "$REPLY" >>${LOG_FILE}
  done
}
