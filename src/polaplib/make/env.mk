# polaplib/make/env.mk
POLAPLIB_DIR ?= $(abspath $(CURDIR))
CONDA_ENV    ?= polap
REQ_CMDS     ?= seqkit csvtk Rscript python3 awk sed grep

.PHONY: env.check env.create env.versions

env.check:
	@ok=1; \
	for c in $(REQ_CMDS); do \
	  if ! command -v $$c >/dev/null 2>&1; then \
	    echo "[ERR] Missing command: $$c" >&2; ok=0; fi; \
	done; \
	[ $$ok -eq 1 ] || { echo "[ERR] Missing required tools"; exit 2; }

env.create:
	@[ -f environment.yml ] || { echo "[ERR] environment.yml not found"; exit 2; }
	@mamba env create -f environment.yml -n $(CONDA_ENV) || conda env create -f environment.yml -n $(CONDA_ENV)

env.versions:
	@echo "== Tool versions =="; \
	for c in $(REQ_CMDS); do \
	  if command -v $$c >/dev/null 2>&1; then \
	    echo -n "$$c: "; $$c --version 2>&1 | head -n1; \
	  fi; \
	done
