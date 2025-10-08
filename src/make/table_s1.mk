# make/table_s1.mk (real)

TABLE_S1_TSV ?= $(TOP)/md/tableS1-dataset-summary.tsv
TABLE_S1_MD  := $(TABLE_S1_TSV:.tsv=.md)

MANIFESTS    ?= $(foreach T,$(TIERS),$(TOP)/md/manifest-$(T)-$(SET).json)
ABS_MANIFESTS:= $(foreach m,$(MANIFESTS),$(abspath $(m)))

.PHONY: table-s1
table-s1: $(TABLE_S1_TSV) $(TABLE_S1_MD)

$(TABLE_S1_TSV): $(MANIFESTS)
	@mkdir -p $(dir $@)
	@if [ -s "$(PLATFORM_TSV)" ]; then \
	  echo "[INFO] Using platform map: $(PLATFORM_TSV)"; \
	  bash $(POLAPLIB_DIR)/polap-bash-make-table-s1.sh \
	    --out "$@" --markdown \
	    $(foreach m,$(ABS_MANIFESTS),--manifest $(m)) \
	    --platform-map "$(PLATFORM_TSV)" \
	    --dedup "$(DEDUP)" --select "$(SELECT)" \
	    --summary-scopes "$(SUMMARY_SCOPES)" --summary-stat "$(SUMMARY_STAT)"; \
	else \
	  echo "[WARN] No platform map at $(PLATFORM_TSV); platform=NA in summaries"; \
	  bash $(POLAPLIB_DIR)/polap-bash-make-table-s1.sh \
	    --out "$@" --markdown \
	    $(foreach m,$(ABS_MANIFESTS),--manifest $(m)) \
	    --dedup "$(DEDUP)" --select "$(SELECT)" \
	    --summary-scopes "$(SUMMARY_SCOPES)" --summary-stat "$(SUMMARY_STAT)"; \
	fi

$(TABLE_S1_MD): $(TABLE_S1_TSV)
	@true
