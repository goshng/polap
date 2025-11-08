#!/usr/bin/env bash
# polaplib/tooltest/seqkit/run.php
# Version: v0.1.0 — test seqkit 'stats' systematically with tiny fixtures.
#
# Tests:
#   1) smoke         : binary present, usable
#   2) fasta_basic   : stats on small FASTA, check num_seqs & sum_len vs our own calc
#   3) fastq_basic   : stats on small FASTQ, check format + num_seqs + sum_len
#   4) stdin_support : pipe FASTA to stdin and ensure it parses (file==-, row produced)
#   5) bad_input     : nonexistent file should fail (non-zero exit)
#
# Output is quiet: each test prints one PASS/FAIL line. On shell errors, your
# polap-lib-failsafe (if enabled in the parent) will show the failing shell line.

set -Eeuo pipefail

tdir="$(cd "$(dirname "${Bash_SOURCE[0]:-${0}}")" && pwd)"
# Harness: quiet asserts & helpers
# shellcheck source=polaplib/tooltest/tooltest-lib.sh
# source "$(dirname "$tdir")/tooltest-lib.sh"
source "$(dirname "$tdir")/common.sh"

seqkit_path="$(command -v seqkit || true)"

smoke() {
	[[ -n "$seqkit_path" ]] || {
		echo "[smoke] FAIL: seqkit not found in PATH" >&2
		return 3
	}
	# Print version for log context (no assertions on version number by default)
	local ver
	ver="$("$seqkit_path" version 2>/dev/null || "$seqkit_path" --version || true)"
	echo "[smoke] seqkit: ${ver:-<no version reported>}"
}

# Parse a TSV header to get column positions by name.
# usage: col_index <colname> <file>  → echoes 1-based index or empty
col_index() {
	local name="$1" file="$2"
	awk -F'\t' -v tgt="$name" 'NR==1{for(i=1;i<=NF;i++) if($i==tgt) {print i; exit 0}}' "$file"
}

fasta_basic() {
	local in="${tdir}/fixtures/toy.fa"
	polap_assert_file "$in"

	local out="${TOUT:-$tdir/_out}/fasta.stats.tsv"
	mkdir -p "$(dirname "$out")"

	# Compute expected sums ourselves
	local exp_n exp_sum
	exp_n="$(grep -c '^>' "$in")"
	exp_sum="$(awk '/^>/{next} {s+=length($0)} END{print s+0}' "$in")"

	# Run seqkit
	seqkit stats -T "$in" >"$out"
	polap_assert_file "$out"

	# column indices
	local i_file i_format i_num i_sum
	i_file=$(col_index "file" "$out")
	i_format=$(col_index "format" "$out")
	i_num=$(col_index "num_seqs" "$out")
	i_sum=$(col_index "sum_len" "$out")

	[[ -n "$i_file" && -n "$i_num" && -n "$i_sum" ]] || {
		echo "[fasta_basic] FAIL: unexpected header in $out" >&2
		return 1
	}

	# Extract row 2 values and compare
	local got_file got_format got_num got_sum
	read -r got_file got_format got_num got_sum < <(awk -F'\t' -v a="$i_file" -v b="$i_format" -v c="$i_num" -v d="$i_sum" 'NR==2{print $a, $b, $c, $d}' "$out")

	[[ "$got_num" == "$exp_n" ]] || {
		echo "[fasta_basic] FAIL: num_seqs=$got_num != $exp_n" >&2
		return 1
	}
	[[ "$got_sum" == "$exp_sum" ]] || {
		echo "[fasta_basic] FAIL: sum_len=$got_sum != $exp_sum" >&2
		return 1
	}
	# Format should be FASTA (case-insensitive)
	echo "$got_format" | grep -qi 'fasta' || {
		echo "[fasta_basic] FAIL: format=$got_format not FASTA" >&2
		return 1
	}

	echo "[fasta_basic] OK ($got_file num_seqs=$got_num sum_len=$got_sum)"
}

fastq_basic() {
	local in="${tdir}/fixtures/toy.fq"
	polap_assert_file "$in"

	local out="${TOUT:-$tdir/_out}/fastq.stats.tsv"
	mkdir -p "$(dirname "$out")"

	# Expected sums
	local exp_n exp_sum
	exp_n="$(grep -c '^@' "$in")"
	# Sum lengths across sequence lines (2nd of each 4-line block)
	exp_sum="$(awk 'NR%4==2{ s+=length($0) } END{ print s+0 }' "$in")"

	seqkit stats -T "$in" >"$out"
	polap_assert_file "$out"

	local i_file i_format i_num i_sum
	# header → indices
	i_file=$(col_index "file" "$out")
	i_format=$(col_index "format" "$out")
	i_num=$(col_index "num_seqs" "$out")
	i_sum=$(col_index "sum_len" "$out")

	# ensure we at least have the core columns; "format" can be absent on older seqkit
	[[ -n "$i_file" && -n "$i_num" && -n "$i_sum" ]] || {
		echo "[fastq_basic] FAIL: unexpected header in $out" >&2
		return 1
	}

	# extract first data row (NR==2) by column indices
	read -r got_file got_format got_num got_sum < <(
		awk -F'\t' \
			-v a="${i_file}" -v b="${i_format:-0}" -v c="${i_num}" -v d="${i_sum}" \
			'NR==2{
         f  = (a>0)?$a:"";
         fm = (b>0)?$b:"";
         n  = (c>0)?$c:"";
         s  = (d>0)?$d:"";
         print f, fm, n, s
       }' "$out"
	)
	echo "$got_format" | grep -qi 'fastq' || {
		echo "[fastq_basic] FAIL: format=$got_format not FASTQ" >&2
		return 1
	}
	[[ "$got_num" == "$exp_n" ]] || {
		echo "[fastq_basic] FAIL: num_seqs=$got_num != $exp_n" >&2
		return 1
	}
	[[ "$got_sum" == "$exp_sum" ]] || {
		echo "[fastq_basic] FAIL: sum_len=$got_sum != $exp_sum" >&2
		return 1
	}

	echo "[fastq_basic] OK ($got_file num_seqs=$got_num sum_len=$got_sum)"
}

stdin_support() {
	local in="${tdir}/fixtures/toy.fa"
	polap_assert_file "$in"
	local out="${TOUT:-$tdir/_out}/stdin.tsv"
	mkdir -p "$(dirname "$out")"

	# Pipe via stdin and request '-' as filename
	cat "$in" | seqkit stats -T - >"$out"
	polap_assert_file "$out"
	# first data row's file column should be '-' (seqkit convention)
	local i_file
	i_file=$(col_index "file" "$out")
	local got
	got="$(awk -F'\t' -v a="$i_file" 'NR==2{print $a}' "$out")"
	[[ "$got" == "-" ]] || {
		echo "[stdin] FAIL: expected file '-' got '$got'" >&2
		return 1
	}
	echo "[stdin] OK (reads from stdin)"
}

bad_input() {
	if seqkit stats /no/such/file 2>/dev/null; then
		echo "[bad_input] FAIL: expected non-zero exit" >&2
		return 1
	fi
	echo "[bad_input] OK (nonexistent file handled)"
}

# ---- Run selected tests ----
main() {
	# Optional output dir override from runner: export TOUT=... or use default
	export TOUT="${TOUT:-$tdir/_out}"
	mkdir -p "$TOUT"

	# 1. smoke
	smoke || exit $?
	# 2. fasta basic
	fasta_basic || exit $?
	# 3. fastq basic
	fastq_basic || exit $?
	# 4. stdin support
	stdin_support || exit $?
	# 5. bad input case
	bad_input || exit $?

	echo "[seqkit] ALL TESTS PASSED"
}

main "$@"
