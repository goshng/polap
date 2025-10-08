# polaplib/make/gfa_stats.mk
GFA_LIST      ?= $(shell git ls-files '*/*.gfa' '*/*.gfa.gz' 2>/dev/null)
GFA_FEATURES  ?= md/graph-features.tsv

.PHONY: gfa.features
gfa.features: $(GFA_FEATURES)

$(GFA_FEATURES): $(GFA_LIST)
	@mkdir -p $(dir $@)
	@echo -e "species\ttag\tn_segments\tn_links\ttotal_len\tN50\tmax_seg\tis_circular" > "$@"
	@for g in $(GFA_LIST); do \
	  sp=$$(echo "$$g" | awk -F/ '{print $$1}'); \
	  tag=$$(basename "$$g" | sed 's/\..*//'); \
	  python3 $(POLAPLIB_DIR)/scripts/gfa_features_from_paths.py --assemblies <( \
	    echo -e "species\torganelle\tgfa\tpng\tcircular_count\n$$sp\t$${tag%.*}\t$$g\tNA\tNA" \
	  ) --out md/_tmp.tsv >/dev/null 2>&1 || true; \
	  if [ -s md/_tmp.tsv ]; then \
	    tail -n +2 md/_tmp.tsv | \
	      awk -v s="$$sp" -v t="$$tag" 'BEGIN{FS=OFS="\t"}{print s,t,$$3,$$4,$$5,$$6,$$7,$$8}' >> "$@"; \
	  fi; \
	done
	@rm -f md/_tmp.tsv || true
