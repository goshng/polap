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
# Convert numbers between different units.
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
# Function to convert base pairs to the highest appropriate unit
# Example usage
# bp=31846726397
# convert_bp $bp
################################################################################
function _polap_lib_genomesize {
	#

	local g=$(<"${_arg_outdir}/short_expected_genome_size.txt")
	echo $g
}

# find_genome_size: estimate genome size from short- or long-reads.
# - Short reads (-a [-b]) use Jellyfish (default) or ntCard when --mem-gb ≤ 2.
# - Long reads (-l) use ntCard. Genome size is inferred with GenomeScope2.
#
# Requirements: jellyfish, ntcard, genomescope2, seqkit, awk, sort, cut, tail.
#
# Usage examples:
#   find_genome_size -a R1.fq.gz -b R2.fq.gz -o out -t 8
#   find_genome_size -l hifi.fastq.gz --tech hifi -o out -t 8
#   find_genome_size -a R1.fq.gz -o out --mem-gb 2   # low-RAM path (ntCard)
#
# Output:
#   <out>/genome_size.txt                  : numeric genome size (bp, plain number)
#   intermediates in <out>/
#
# Options:
#   -a <file>        Short-read R1 (required for SR mode)
#   -b <file>        Short-read R2 (optional)
#   -l <file>        Long-read FASTQ (required for LR mode if no -a)
#   -o <dir>         Output directory (default: ./genomesize_out)
#   -t <int>         Threads (default: 4)
#   --tech <str>     Long-read technology: hifi|pbclr|ont (default: hifi)
#   --mem-gb <num>   Memory cap in GB; ≤2 forces ntCard on SR (default: 2)
#   --k-sr <int>     k-mer for SR (default: 21)
#   --k-hifi <int>   k-mer for HiFi (default: 21)
#   --k-clr <int>    k-mer for CLR/ONT (default: 41)
#   --lr-minlen <n>  Min read length prefilter for LR (default: 1000)
#   --redo           Recompute even if outputs exist
#   -h|--help        Show help
#
_polap_lib_genomesize-estimate() {
	# defaults
	local R1="" R2="" LR=""
	local OUTDIR="./genomesize_out"
	local THREADS=4
	local TECH="hifi" # hifi|pbclr|ont
	local MEM_GB=2
	local K_SR=21
	local K_HIFI=21
	local K_CLR=41
	local LR_MINLEN=1000
	local REDO=1

	# parse options (short with getopts; long manually)
	local argv=()
	while (("$#")); do
		case "$1" in
		-a)
			R1="$2"
			shift 2
			;;
		-b)
			R2="$2"
			shift 2
			;;
		-l)
			LR="$2"
			shift 2
			;;
		-o)
			OUTDIR="$2"
			shift 2
			;;
		-t)
			THREADS="$2"
			shift 2
			;;
		--tech)
			TECH="$2"
			shift 2
			;;
		--mem-gb)
			MEM_GB="$2"
			shift 2
			;;
		--k-sr)
			K_SR="$2"
			shift 2
			;;
		--k-hifi)
			K_HIFI="$2"
			shift 2
			;;
		--k-clr)
			K_CLR="$2"
			shift 2
			;;
		--lr-minlen)
			LR_MINLEN="$2"
			shift 2
			;;
		--no-redo)
			REDO=0
			shift
			;;
		-h | --help)
			sed -n '1,70p' "$0" | sed -n '/^# find_genome_size:/,$p' | sed '/^find_genome_size()/q'
			return 0
			;;
		--)
			shift
			break
			;;
		-*)
			echo "[ERR] Unknown option: $1" >&2
			return 2
			;;
		*)
			argv+=("$1")
			shift
			;;
		esac
	done

	# helper: die
	_die() {
		echo "[ERR] $*" >&2
		return 1
	}

	# helper: run GenomeScope2 and extract a single number
	_hist_to_genomesize() {
		local HIST="$1" K="$2" OUTTXT="$3"
		local OUTDIR_GS="${HIST}.gscope2"
		mkdir -p "${OUTDIR_GS}" || return 1
		genomescope2 -i "${HIST}" -o "${OUTDIR_GS}" -k "${K}" -p 1 >"${OUTTXT}.genomescope2.txt" || return 1
		cut -d: -f6 "${OUTTXT}.genomescope2.txt" | tail -n 1 >"${OUTTXT}" || return 1
	}

	# sanity
	mkdir -p "${OUTDIR}" || _die "Cannot create outdir: ${OUTDIR}"
	local GS_TXT="${OUTDIR}/genome_size.txt"

	# MODE: SR if -a present; otherwise LR if -l present
	if [[ -n "${R1}" ]]; then
		[[ -s "${R1}" ]] || _die "missing file: ${R1}"
		local K="${K_SR}"
		if ((REDO == 0)) && [[ -s "${GS_TXT}" ]]; then
			echo "[INFO] Found ${GS_TXT}; use --redo to recompute."
			cat "${GS_TXT}"
			return 0
		fi

		if (($(printf "%.0f" "${MEM_GB}") <= 2)); then
			# low-RAM: ntCard -> GenomeScope2
			echo "[INFO] SR low-RAM mode: ntCard -> GenomeScope2 (k=${K}, mem=${MEM_GB}G)"
			local NTP="${OUTDIR}/ntc${K}"
			local NTHIST="${NTP}_k${K}.hist"
			if ((REDO == 1)) || [[ ! -s "${NTHIST}" ]]; then
				if [[ -n "${R2}" && -s "${R2}" ]]; then
					ntcard -k${K} -t "${THREADS}" -p "${NTP}" "${R1}" "${R2}" || _die "ntcard failed"
				else
					ntcard -k${K} -t "${THREADS}" -p "${NTP}" "${R1}" || _die "ntcard failed"
				fi
			fi
			[[ -s "${NTHIST}" ]] || _die "ntCard histogram not found: ${NTHIST}"
			local CLEAN="${OUTDIR}/k${K}.clean.hist"
			awk '($1~/^[0-9]+$/)&&($2~/^[0-9]+$/){print $1"\t"$2}' "${NTHIST}" | sort -k1,1n >"${CLEAN}" || _die "clean hist failed"
			_hist_to_genomesize "${CLEAN}" "${K}" "${GS_TXT}" || _die "GenomeScope2 failed"
		else
			# normal-RAM: Jellyfish -> GenomeScope2 (cap -s by mem)
			local JF="${OUTDIR}/reads.jf"
			local HIST="${OUTDIR}/reads.histo"
			# heuristic hash size: ~0.6×mem (GiB), min 0.5G
			local JF_S=$(awk -v m="${MEM_GB}" 'BEGIN{s=m*0.6;if(s<0.5)s=0.5;printf "%.1fG",s}')
			echo "[INFO] SR Jellyfish mode: -m ${K} -s ${JF_S} -t ${THREADS}"
			if ((REDO == 1)) || [[ ! -s "${JF}" ]]; then
				if [[ -n "${R2}" && -s "${R2}" ]]; then
					jellyfish count -t "${THREADS}" -C -m "${K}" -s "${JF_S}" -o "${JF}" --min-qual-char=? "${R1}" "${R2}" || _die "jellyfish count failed"
				else
					jellyfish count -t "${THREADS}" -C -m "${K}" -s "${JF_S}" -o "${JF}" --min-qual-char=? "${R1}" || _die "jellyfish count failed"
				fi
			fi
			if ((REDO == 1)) || [[ ! -s "${HIST}" ]]; then
				jellyfish histo -t "${THREADS}" -o "${HIST}" "${JF}" || _die "jellyfish histo failed"
			fi
			_hist_to_genomesize "${HIST}" "${K}" "${GS_TXT}" || _die "GenomeScope2 failed"
		fi

		echo "[OK] genome size -> ${GS_TXT}"
		cat "${GS_TXT}"
		return 0

	elif [[ -n "${LR}" ]]; then
		# Long-read path: ntCard (with prefilter) -> GenomeScope2
		[[ -s "${LR}" ]] || _die "missing file: ${LR}"
		local K="${K_HIFI}"
		[[ "${TECH}" == "pbclr" || "${TECH}" == "ont" ]] && K="${K_CLR}"
		if ((REDO == 0)) && [[ -s "${GS_TXT}" ]]; then
			echo "[INFO] Found ${GS_TXT}; use --redo to recompute."
			cat "${GS_TXT}"
			return 0
		fi
		echo "[INFO] LR mode (tech=${TECH}) ntCard -> GenomeScope2 (k=${K})"

		local LR_MIN="${OUTDIR}/lr.min${LR_MINLEN}.fq"
		if ((REDO == 1)) || [[ ! -s "${LR_MIN}" ]]; then
			seqkit seq -m "${LR_MINLEN}" "${LR}" -o "${LR_MIN}" || _die "seqkit prefilter failed"
		fi

		local NTP="${OUTDIR}/ntc${K}"
		local NTHIST="${NTP}_k${K}.hist"
		if ((REDO == 1)) || [[ ! -s "${NTHIST}" ]]; then
			ntcard -k${K} -t "${THREADS}" -p "${NTP}" "${LR_MIN}" || _die "ntcard failed"
		fi
		[[ -s "${NTHIST}" ]] || _die "ntCard histogram not found: ${NTHIST}"

		local CLEAN="${OUTDIR}/k${K}.clean.hist"
		awk '($1~/^[0-9]+$/)&&($2~/^[0-9]+$/){print $1"\t"$2}' "${NTHIST}" | sort -k1,1n >"${CLEAN}" || _die "clean hist failed"
		_hist_to_genomesize "${CLEAN}" "${K}" "${GS_TXT}" || _die "GenomeScope2 failed"

		echo "[OK] genome size -> ${GS_TXT}"
		cat "${GS_TXT}"
		return 0

	else
		_die "Provide either short reads (-a [ -b ]) OR long reads (-l). Try -h."
	fi
}
