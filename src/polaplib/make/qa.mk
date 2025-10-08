# polaplib/make/qa.mk
BASH_SCRIPTS ?= $(shell git ls-files 'polaplib/*.sh' 'polaplib/scripts/*.sh' 2>/dev/null)
PY_SCRIPTS   ?= $(shell git ls-files 'polaplib/scripts/*.py' 2>/dev/null)
R_SCRIPTS    ?= $(shell git ls-files 'polaplib/scripts/*.R' 2>/dev/null)

.PHONY: qa fmt lint test

qa: lint test

fmt:
	@if command -v black >/dev/null 2>&1; then black $(PY_SCRIPTS); fi
	@if command -v ruff >/dev/null 2>&1; then ruff --fix $(PY_SCRIPTS); fi
	# add shfmt if desired
	@if command -v shfmt >/dev/null 2>&1; then shfmt -w $(BASH_SCRIPTS); fi

lint:
	@if command -v shellcheck >/dev/null 2>&1; then shellcheck $(BASH_SCRIPTS); fi
	@if command -v ruff >/dev/null 2>&1; then ruff $(PY_SCRIPTS); fi
	@if command -v flake8 >/dev/null 2>&1; then flake8 $(PY_SCRIPTS); fi
	@if command -v lintr >/dev/null 2>&1; then \
	  for r in $(R_SCRIPTS); do Rscript -e "lintr::lint('$$r')" || true; done; \
	fi

test:
	@if command -v pytest >/dev/null 2>&1; then pytest -q; else echo "[WARN] pytest not found, skipping"; fi
