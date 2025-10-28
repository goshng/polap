#!/usr/bin/env bash
# Version: v0.1.0
# polap-bash-fq2gfa.sh
# Function: polap_fq_to_gfa  — input FASTQ(.gz), output a single merged GFA via seqkit+minimap2+miniasm
#
# Requirements:
#   - seqkit (for stats & splitting)     https://bioinf.shenwei.me/seqkit/             # docs
#   - minimap2 (ava-ont/ava-pb presets)  https://lh3.github.io/minimap2/minimap2.html  # manual
#   - miniasm                             https://nanoporetech.com/resource-centre/minimap-and-miniasm-fast-mapping-and-de-novo-assembly-noisy-long-sequences-0
#   - (optional) gfak (GFAKluge)          https://edawson.github.io/gfakluge/          # for robust GFA merge
#   - (optional) gfatools                 https://github.com/lh3/gfatools
#
# Usage (examples):
#   source polap-bash-fq2gfa.sh
#   polap_fq_to_gfa reads.fq.gz -o work -t 16 --mode ont --gb 1
#   polap_fq_to_gfa reads.fq      -o work -t 32 --mode hifi --gb 0.5
#
# Outputs:
#   work/gfa/*.gfa                      (per-shard miniasm outputs)
#   work/all_shards.prefixed.gfa        (final merged GFA; printed on stdout as the path)

set -euo pipefail

polap_fq_to_gfa() {
	local in_fq="${1:?FASTQ(.gz) required}"
	shift
	local outdir="work"
	local threads=8
	local mode="ont" # ont|hifi
	local gb="10"    # gigabases per shard (float ok)
	local keep_shards=0

	local assembler="miniasm" # miniasm|raven

	# parse args
	while [[ $# -gt 0 ]]; do
		case "$1" in
		-o | --outdir)
			outdir="${2:?}"
			shift 2
			;;
		-t | --threads)
			threads="${2:?}"
			shift 2
			;;
		--mode)
			mode="${2:?}"
			shift 2
			;;
		--gb)
			gb="${2:?}"
			shift 2
			;;
		--keep-shards)
			keep_shards=1
			shift
			;;
		--assembler)
			assembler="${2:?}"
			shift 2
			;;
		*)
			echo "Unknown option: $1" >&2
			return 2
			;;
		esac
	done

	mkdir -p "${outdir}"/{shards,paf,gfa,log}
	local shards="${outdir}/shards"

	# 2) Split FASTQ into shards of ~reads_per_shard sequences
	seqkit split2 -l "$gb"g -O "${shards}" "$in_fq" 2>/dev/null

	# Normalize shard names to part_*.fq.gz
	find "${shards}" -type f -name "*.part_*.fq.gz" | while read -r f; do
		b=$(basename "$f")
		# turn X.part_001.fq.gz -> part_001.fq.gz
		bn="part_${b##*.part_}"
		ln="${shards}/${bn}"
		[[ -e "$ln" ]] || ln -s "$b" "$ln"
	done

	# 3) Per-shard overlaps + assembly
	local mm_preset="-x ava-ont" # all-vs-all overlaps preset for long reads.  [oai_citation:2‡GitHub Pages](https://lh3.github.io/minimap2/minimap2.html?utm_source=chatgpt.com)
	case "$mode" in
	ont) mm_preset="-x ava-ont" ;;
	hifi) mm_preset="-x ava-pb" ;;
	*)
		echo "ERROR: --mode must be ont|hifi" >&2
		return 2
		;;
	esac

	#
	# shopt -s nullglob
	# for fqgz in "${shards}"/part_*.fq.gz; do
	for fqgz in $(compgen -G "${shards}/part_*.fq.gz"); do

		base=$(basename "$fqgz" .fq.gz)
		local paf="${outdir}/paf/${base}.paf"
		local paf_err="${outdir}/paf/${base}.err"
		local gfa="${outdir}/gfa/${base}.gfa"
		local gfa_err="${outdir}/gfa/${base}.err"
		echo "[INFO] Assembling shard ${base} ..." | tee -a "${outdir}/log/steps.log"

		if [[ "$assembler" == "raven" ]]; then
			# Raven can read gz FASTQ directly; write graph to GFA.
			# Use --disable-polishing for speed if you only need seed contigs.
			raven --threads "${threads}" --disable-polishing \
				--graph-output "${gfa}" \
				--output /dev/null \
				"$fqgz"
		else
			# Default: minimap2 (all-vs-all) + miniasm
			# Decompress to a temp plain fq for miniasm -f pairing
			local tmpfq="${outdir}/${base}.tmp.fq"
			gzip -cd "$fqgz" >"$tmpfq"

			minimap2 $mm_preset -t "${threads}" "$tmpfq" "$tmpfq" \
				>"$paf" \
				2>"paf_err"

			miniasm -f "$tmpfq" "$paf" >"$gfa" 2>"$gfa_err"

			rm -f "$tmpfq"
		fi
	done

	# 4) Merge GFAs safely:
	#    Prefer gfak (ids+merge). Otherwise awk-prefix IDs and concatenate.
	local merged="${outdir}/miniasm.gfa"
	if command -v gfak >/dev/null 2>&1; then
		# GFAKluge: coordinate ID spaces & merge.  [oai_citation:4‡Edawson](https://edawson.github.io/gfakluge/?utm_source=chatgpt.com)
		tmpd="$(mktemp -d)"
		for g in "${outdir}"/gfa/*.gfa; do
			b=$(basename "$g" .gfa)
			gfak ids -p "${b}_" "$g" >"${tmpd}/${b}.pref.gfa"
		done
		gfak merge "${tmpd}"/*.pref.gfa >"$merged"
		rm -rf "$tmpd"
	else
		# portable fallback: prefix IDs via awk, add one header
		awk 'BEGIN{print "H\tVN:Z:1.0"}' >"$merged"
		for g in "${outdir}"/gfa/*.gfa; do
			b=$(basename "$g" .gfa)
			awk -v P="${b}_" 'BEGIN{FS=OFS="\t"}
        /^S\t/{$2=P $2}
        /^L\t/{$2=P $2; $4=P $4}
        /^P\t/{$2=P $2}
        /^H\t/ {next}
        {print}
      ' "$g" >>"$merged"
		done
	fi

	echo "$merged"
	# Optional: keep shard FASTQs if requested; else remove the symlinked normalized names
	if [[ $keep_shards -eq 0 ]]; then
		find "${shards}" -type l -name 'part_*.fq.gz' -delete || true
	fi
}

# If invoked as a script, allow direct run:
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
	polap_fq_to_gfa "$@"
fi
