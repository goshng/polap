# polaplib/tooltest/Makefile
# Robust test harness for tool-specific run.sh scripts.
# - Discovers tests as subdirectories containing a run.sh.
# - Runs all or a filtered subset, forwards extra args/environment.
# - Captures logs per test, prints PASS/FAIL summary.
# - Integrates with polap-lib-failsafe.sh (enable via POLAP_FAILSAFE=1).

SHELL := /usr/bin/env bash
.DEFAULT_GOAL := run

# --- Paths & Env --------------------------------------------------------------

# Absolute path to this directory (polaplib/tooltest)
THIS_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

# Where to find polap libraries (set by caller or fallback to parent)
POLAPLIB ?= $(abspath $(THIS_DIR)/..)
export POLAPLIB

# Where to put test outputs and logs
OUT_DIR ?= $(abspath $(THIS_DIR)/_out)
LOG_DIR ?= $(abspath $(THIS_DIR)/_logs)
export TOUT ?= $(OUT_DIR)
export LOG_FILE ?=

# Failsafe & verbosity (see polap-lib-failsafe.sh)
POLAP_FAILSAFE ?= 1
POLAP_VERBOSE  ?= 1
export POLAP_FAILSAFE POLAP_VERBOSE

# Extra arguments forwarded to each run.sh (e.g., ARGS="--stack" or " --code 9")
ARGS ?=

# --- Discovery ----------------------------------------------------------------
# Discover all test dirs that contain run.sh (exclude hidden dirs)
ALL_TESTS := $(sort $(patsubst %/,%,$(filter-out ./_%, $(dir $(wildcard $(THIS_DIR)*/run.sh)))))

# Allow override / filtering, e.g.:
#   make run TESTS="minimap2 seqkit"
#   make run PATTERN=seq
TESTS ?=
PATTERN ?=
ifeq ($(strip $(TESTS)),)
  ifneq ($(strip $(PATTERN)),)
    FILTERED := $(filter $(PATTERN)%,$(ALL_TESTS))
  else
    FILTERED := $(ALL_TESTS)
  endif
else
  FILTERED := $(TESTS)
endif

# --- Helper functions ---------------------------------------------------------

define run_one
  @echo "────────────────────────────────────────────────────────"
  @echo "▶  Running: $(1)"
  @mkdir -p "$(OUT_DIR)/$(1)" "$(LOG_DIR)"
  @echo "[info] TOUT=$(OUT_DIR)/$(1)"
  @echo "[info] log → $(LOG_DIR)/$(1).log"
  @ ( \
      export TOUT="$(OUT_DIR)/$(1)"; \
      export POLAPLIB="$(POLAPLIB)"; \
      export POLAP_FAILSAFE=$(POLAP_FAILSAFE); \
      export POLAP_VERBOSE=$(POLAP_VERBOSE); \
      export LOG_FILE="$(LOG_DIR)/$(1).log"; \
      cd "$(THIS_DIR)/$(1)" && \
      bash ./run.sh $(ARGS); \
    ); \
    status="$$?"; \
    if [ "$$status" -eq 0 ]; then \
      echo "✔ PASS $(1)"; \
    else \
      echo "✖ FAIL($(status)) $(1) — see $(LOG_DIR)/$(1).log"; \
    fi; \
    echo "$$status" > "$(LOG_DIR)/$(1).rc"; \
    exit "$$status"
endef

# --- Phony targets ------------------------------------------------------------

.PHONY: help list run all clean clean-logs clean-out \
        $(ALL_TESTS) run-%

help:
	@echo "polap tooltest harness"
	@echo
	@echo "Usage:"
	@echo "  make run                 # run all discovered tests"
	@echo "  make run TESTS='t1 t2'   # run only listed tests"
	@echo "  make run PATTERN=seq     # run tests whose dir starts with 'seq'"
	@echo "  make list                # list discovered tests"
	@echo "  make clean               # remove _out and _logs"
	@echo
	@echo "Env vars:"
	@echo "  POLAPLIB=$(POLAPLIB)"
	@echo "  OUT_DIR=$(OUT_DIR)"
	@echo "  LOG_DIR=$(LOG_DIR)"
	@echo "  POLAP_FAILSAFE=$(POLAP_FAILSAFE)  (set to 0 to disable bash failsafe)"
	@echo "  POLAP_VERBOSE=$(POLAP_VERBOSE)"
	@echo "  ARGS=$(ARGS)            (extra args forwarded to each run.sh)"
	@echo

list:
	@echo "Discovered tests under $(THIS_DIR):"
	@$(foreach t,$(ALL_TESTS),echo "  - $(t)";)

# Run all filtered tests with summary; continue on failure, summarize rc.
run:
	@echo "Found tests: $(ALL_TESTS)"
	@echo "Selected:    $(FILTERED)"
	@mkdir -p "$(OUT_DIR)" "$(LOG_DIR)"
	@passes=0; fails=0; total=0; \
	for t in $(FILTERED); do \
	  total=$$((total+1)); \
	  echo "────────────────────────────────────────────────────────"; \
	  echo "▶ Running $$t"; \
	  ( export TOUT="$(OUT_DIR)/$$t"; \
	    export POLAPLIB="$(POLAPLIB)"; \
	    export POLAP_FAILSAFE=$(POLAP_FAILSAFE); \
	    export POLAP_VERBOSE=$(POLAP_VERBOSE); \
	    export LOG_FILE="$(LOG_DIR)/$$t.log"; \
	    cd "$(THIS_DIR)/$$t" && bash ./run.sh $(ARGS) ) ; \
	  rc="$$?"; \
	  if [ "$$rc" -eq 0 ]; then \
	    echo "✔ PASS $$t"; passes=$$((passes+1)); \
	  else \
	    echo "✖ FAIL($$rc) $$t — see $(LOG_DIR)/$$t.log"; fails=$$((fails+1)); \
	  fi; \
	  echo "$$rc" > "$(LOG_DIR)/$$t.rc"; \
	done; \
	echo "────────────────────────────────────────────────────────"; \
	echo "Summary: $$passes passed / $$((passes+fails)) total (fails=$$fails)"; \
	if [ "$$fails" -gt 0 ]; then exit 1; fi

# Allow: make minimap2  (runs polaplib/tooltest/minimap2/run.sh)
$(ALL_TESTS):
	$(call run_one,$@)

# Allow: make run-seqkit (alias for: run only seqkit)
run-%:
	$(call run_one,$*)

clean:
	@echo "Removing $(OUT_DIR) and $(LOG_DIR)…"
	@rm -rf "$(OUT_DIR)" "$(LOG_DIR)"

clean-logs:
	@echo "Removing logs: $(LOG_DIR)…"
	@rm -rf "$(LOG_DIR)"

all: run

