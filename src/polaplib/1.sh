grep -oP '^\s*function\s+_run_polap_[a-zA-Z0-9_]+\s*\(\)?|^\s*_run_polap_[a-zA-Z0-9_]+\s*\(\)' *.sh |
	sed -E 's/^\s*(function\s+)?_run_polap_([a-zA-Z0-9_]+)\s*\(\)?.*/\2/' |
	sort -u
