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

# include and execute other BASH and R scripts
WDIR="$(dirname "$0")"
WDIR="${BASH_SOURCE%/*}"
if [[ ! -d "$WDIR" ]]; then
	WDIR="$PWD"
fi
WDIR=${script_dir}

# variables for input data file names for flexible data processing.
LR=${_arg_long_reads}       # long-read data file
SR1=${_arg_short_read1}     # paired short-read data file 1
SR2=${_arg_short_read2}     # paired short-read data file 2
PA=${_arg_unpolished_fasta} # assembled draft sequence extracted from bandage
FA=${_arg_final_assembly}   # polished sequence

# variables for output
INUM=${_arg_inum}
JNUM=${_arg_jnum}
FDIR="$ODIR"/0 # flye 1st output
if [ "${_arg_archive_is}" = "off" ]; then
	_arg_archive="${ODIR}-a"
fi

# tuning variables for optimal performance
LRNK="$ODIR/nk.fq.gz"
MR=${_arg_min_read_length}
MPAIR=${_arg_pair_min}     # 3000 for MT, 1000 for PT
MBRIDGE=${_arg_bridge_min} # used to be 3000,
MSINGLE=${_arg_single_min} # not used deprecated
COV=${_arg_coverage}
# NT=$(cat /proc/cpuinfo | grep -c processor)
NT=${_arg_threads}
if test -z "$DEBUG"; then
	DEBUG=0
fi
CIRCULARIZE=${_arg_circularize} # "--circularize"
SPECIES=${_arg_species}

################################################################
# Variables
SRA=${_arg_sra}
SRALONG=""
SRASHORT=""
RESUME=${_arg_resume}
ALL_ANNOTATE="--selective-annotate"
FLYE_CONTIGGER="--contigger"
USE_EDGES="--no-use-edges"
NO_REDUCTION_READS=${_arg_reduction_reads}
NO_COVERAGE_CHECK=${_arg_coverage_check}

# Constants
EXIT_SUCCESS=0
EXIT_FAIL=1
EXIT_ERROR=2
RETURN_SUCCESS=0
RETURN_FAIL=1

SECONDS=0

if [ "${_arg_log_is}" = "off" ]; then
	LOG_FILE="${ODIR}/polap.log"
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
