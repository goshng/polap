#!/usr/bin/env bash
# Version: v0.2.0
# polap-bash-gfa-extract-organelle.sh
#
# Goal: Extract plastid/mitochondrial sequences from a Flye assembly_graph.gfa
# using gfatk (preferred) or OATK path_to_fasta when an explicit path is given.
#
# Usage:
#   polap-bash-gfa-extract-organelle.sh \
#     --gfa assembly_graph.gfa \
#     --outdir outdir \
#     [--target mito|pltd|both] \
#     [--linearize] \
#     [--oatk-path "u1+,u2-,u3+"] [--oatk-path-file path.txt] \
#     [--prefix SAMPLE] \
#     [--threads N]
#
# Outputs (in --outdir):
#   <prefix>.pltd.gfa, <prefix>.pltd.segs.fa, <prefix>.pltd.linear.fa (if --linearize)
#   <prefix>.mito.gfa, <prefix>.mito.segs.fa, <prefix>.mito.linear.fa (if --linearize)
#   <prefix>.path_to_fasta.fa (when using OATK with an explicit path)
#
# Requirements:
#   - gfatk in PATH (install helper provided separately)
#   - (optional) OATK’s path_to_fasta in PATH for explicit path export
#   - GNU awk, coreutils
#
# Notes:
#   - gfatk has organelle heuristics:
#       extract-chloro (plastid), extract-mito (mitochondrial)  [reads Flye GFA v1]
#   - If heuristics are not decisive for your dataset, supply an explicit node path
#     (from Bandage or your notes) via --oatk-path or --oatk-path-file for OATK.
#
set -euo pipefail
IFS=$'\n\t'

# ---- defaults
gfa=""
outdir="gfa-organelle"
target="both" # mito, pltd, both
linearize=false
oatk_path=""
oatk_path_file=""
prefix="sample"
threads="${THREADS:-4}"

show_help() {
	sed -n '1,120p' "$0" | sed 's/^# \{0,1\}//'
}

# ---- parse args
while (($#)); do
	case "$1" in
	--gfa)
		gfa="$2"
		shift 2
		;;
	--outdir)
		outdir="$2"
		shift 2
		;;
	--target)
		target="$2"
		shift 2
		;;
	--linearize)
		linearize=true
		shift
		;;
	--oatk-path)
		oatk_path="$2"
		shift 2
		;;
	--oatk-path-file)
		oatk_path_file="$2"
		shift 2
		;;
	--prefix)
		prefix="$2"
		shift 2
		;;
	--threads)
		threads="$2"
		shift 2
		;;
	-h | --help)
		show_help
		exit 0
		;;
	*)
		echo "ERR: unknown option: $1" >&2
		exit 2
		;;
	esac
done

# ---- sanity
if [[ -z "$gfa" || ! -s "$gfa" ]]; then
	echo "ERR: --gfa is required and must exist: $gfa" >&2
	exit 2
fi
mkdir -p "$outdir"

need() { command -v "$1" >/dev/null 2>&1 || {
	echo "ERR: missing tool: $1" >&2
	exit 127
}; }
need awk
if [[ -n "$oatk_path" || -n "$oatk_path_file" ]]; then
	need path_to_fasta
fi
# gfatk is preferred, but if absent and no OATK path given, we fail early.
if ! command -v gfatk >/dev/null 2>&1; then
	if [[ -z "$oatk_path" && -z "$oatk_path_file" ]]; then
		echo "ERR: gfatk not found. Either install gfatk or supply --oatk-path{,-file}." >&2
		exit 127
	fi
fi

# ---- helpers
scripts_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/scripts"
validator="$scripts_dir/validate_oatk_path.awk"

_ts() { date +"[%Y-%m-%d %H:%M:%S]"; }
log() { echo "$(_ts) $*"; }

# ---- gfatk path
run_gfatk_branch() {
	local which="$1" # pltd or mito
	local subcmd out_gfa out_segs out_linear
	case "$which" in
	pltd) subcmd="extract-chloro" ;;
	mito) subcmd="extract-mito" ;;
	*)
		echo "ERR: internal: unknown organelle $which" >&2
		exit 3
		;;
	esac

	out_gfa="${outdir}/${prefix}.${which}.gfa"
	out_segs="${outdir}/${prefix}.${which}.segs.fa"
	out_linear="${outdir}/${prefix}.${which}.linear.fa"

	log "gfatk $subcmd …"
	gfatk "$subcmd" "$gfa" >"$out_gfa"

	log "gfatk fasta (segments) …"
	gfatk fasta "$out_gfa" >"$out_segs"

	if $linearize; then
		# try to coerce a single linear path (longest legal path in subgraph)
		log "gfatk linear …"
		gfatk linear "$out_gfa" >"$out_linear" || {
			echo "WARN: gfatk linear failed; keeping only segment FASTA for $which" >&2
		}
	fi

	log "done $which → $out_gfa ; $out_segs ${linearize:+; $out_linear}"
}

# ---- OATK explicit path mode
run_oatk_path2fa() {
	local pathstr="$1"
	local outf="${outdir}/${prefix}.path_to_fasta.fa"
	log "OATK path_to_fasta on explicit path"
	# Validate path format (e.g., u1+,u2-,u3+)
	if ! awk -f "$validator" <(echo "$pathstr") >/dev/null; then
		echo "ERR: invalid path format for path_to_fasta (expect comma-delimited node+orient)" >&2
		echo "    example: edge_1+,edge_3-,edge_9+" >&2
		exit 2
	fi
	path_to_fasta "$gfa" "$pathstr" >"$outf"
	log "done OATK → $outf"
}

# ---- main
log "Start: target=$target  linearize=$linearize  prefix=$prefix"
case "$target" in
pltd) command -v gfatk >/dev/null 2>&1 && run_gfatk_branch pltd ;;
mito) command -v gfatk >/dev/null 2>&1 && run_gfatk_branch mito ;;
both) command -v gfatk >/dev/null 2>&1 && {
	run_gfatk_branch pltd
	run_gfatk_branch mito
} ;;
*)
	echo "ERR: --target must be mito|pltd|both" >&2
	exit 2
	;;
esac

# Optional explicit path export (useful when Bandage gave you a path)
if [[ -n "$oatk_path_file" ]]; then
	oatk_path="$(tr -d '\n\r\t ' <"$oatk_path_file")"
fi
if [[ -n "$oatk_path" ]]; then
	run_oatk_path2fa "$oatk_path"
fi

log "All done."
