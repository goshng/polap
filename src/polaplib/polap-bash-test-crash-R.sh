#!/usr/bin/env bash
# Version: v0.1.0
# Demo: Bash → R crash with R-FAILED + traceback.
set -Eeuo pipefail
_POLAPLIB_DIR="${_POLAPLIB_DIR:-"$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"}"

R_PROFILE="${_POLAPLIB_DIR}/scripts/polap_r_profile.R" \
	R_PROFILE_USER="${_POLAPLIB_DIR}/scripts/polap_r_profile.R" \
	Rscript "${_POLAPLIB_DIR}/scripts/test_crash.R"
