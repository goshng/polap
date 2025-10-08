# make/plots.mk — complexity plots (segments vs total length)

GF_TSV         ?= $(TOP)/md/graph-features.tsv
COMPLEX_ALL    ?= $(OUTDIR)/complexity-$(SET).pdf
COMPLEX_ALL_LOG?= $(OUTDIR)/complexity-$(SET)-log.pdf

.PHONY: complexity complexity-log
complexity: $(COMPLEX_ALL)
complexity-log: $(COMPLEX_ALL_LOG)

# If you need to build GF_TSV from your pipeline, add a rule here

$(COMPLEX_ALL): $(GF_TSV)
	@mkdir -p $(dir $@)
	@Rscript --vanilla $(POLAPLIB_DIR)/scripts/make_complexity_plot.R \
	  --graph "$(GF_TSV)" \
	  --out "$@" \
	  --title "Assembly complexity (set=$(SET))" \
	  --label-outliers

$(COMPLEX_ALL_LOG): $(GF_TSV)
	@mkdir -p $(dir $@)
	@Rscript --vanilla $(POLAPLIB_DIR)/scripts/make_complexity_plot.R \
	  --graph "$(GF_TSV)" \
	  --out "$@" \
	  --title "Assembly complexity (log x, set=$(SET))" \
	  --logx --label-outliers
