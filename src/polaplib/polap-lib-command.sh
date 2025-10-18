#!/usr/bin/env bash
# polap-lib-command.sh
# Version: v0.1.1
# SPDX-License-Identifier: GPL-3.0-or-later
set -euo pipefail

# Safely join tokens until a token that starts with '-' or equals '0' or '--'
concat_until_zero_or_anydash() {
	local sep="${1:- }"
	shift || true
	local out="" tok
	for tok in "$@"; do
		[[ "$tok" == "0" || "$tok" == "--" || "$tok" == -* ]] && break
		[[ -n "$out" ]] && out+="$sep"
		out+="$tok"
	done
	printf '%s\n' "$out"
}

# Print tokens after first dashed token (or after '--'), stop at '0'
get_after_dash_until_zero() {
	local seen_dash=0 tok
	for tok in "$@"; do
		[[ "$tok" == "0" ]] && break
		if [[ "$tok" == "--" || "$tok" == -* ]]; then
			seen_dash=1
			[[ "$tok" == "--" ]] && continue
			printf '%s\n' "$tok"
			continue
		fi
		((seen_dash)) && printf '%s\n' "$tok"
	done
}

# Very small suggester (Levenshtein distance <= 2) in pure awk
_polap_suggest_closest() {
	local needle="$1"
	shift || true
	awk -v n="$needle" '
    function min(a,b){return a<b?a:b}
    function lev(a,b,  i,j,c,la,lb,d,prev,cur){
      la=length(a); lb=length(b)
      for(i=0;i<=lb;i++) cur[i]=i
      for(i=1;i<=la;i++){
        split("",prev); for(j=0;j<=lb;j++) prev[j]=cur[j]
        cur[0]=i
        for(j=1;j<=lb;j++){
          c = (substr(a,i,1)==substr(b,j,1)) ? 0 : 1
          cur[j]=min(min(prev[j]+1, cur[j-1]+1), prev[j-1]+c)
        }
      }
      return cur[lb]
    }
    { if (lev(n,$0) <= 2) print $0 }
  ' <<<"$(printf '%s\n' "$@")"
}
