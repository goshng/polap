# Require bash
[[ -n "${BASH_VERSION:-}" ]] || {
	echo "[bolap-complete] bash required" >&2
	return 1 2>/dev/null || exit 1
}

# Release: bolap in PATH
bolap_completion() {
	# Enable completion in current shell
	eval "$(bolap completion)"
}

# Dev: use your local bolap.sh via $BOLAP, keep your alias 'bl' untouched
bl_completion() {
	if [[ -z "${BOLAP:-}" ]]; then
		echo "[bolap-complete] Set BOLAP=/abs/path/to/bolap.sh first." >&2
		return 1
	fi
	eval "$(bash "$BOLAP" completion)"
}

# Optional: auto-enable on source if requested
# if [[ "${BOLAP_AUTOCOMP:-0}" == 1 ]]; then
if command -v bolap >/dev/null 2>&1; then
	bolap_completion
elif [[ -n "${BOLAP:-}" ]]; then
	bl_completion
fi
# fi
