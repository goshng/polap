# make/manifests.mk

MANIFESTS := $(foreach T,$(TIERS),$(MANIFEST_DIR)/manifest-$(T)-$(SET).json)

.PHONY: manifests
manifests: $(MANIFESTS)

$(MANIFEST_DIR)/manifest-%-$(SET).json:
	@mkdir -p $(MANIFEST_DIR)
	@echo "[INFO] Building manifest tier=$* set=$(SET)"
	@bash -o pipefail -e -x "$(POLAPLIB_DIR)/polap-bash-make-manifest.sh" \
	  --set "$(SET)" \
	  --tier "$*" \
	  --inum "$(INUM)" \
	  --out "$@" \
	  --pretty \
	  2> "$(MANIFEST_DIR)/manifest-$*-$(SET).log" \
	|| { \
	  echo ""; \
	  echo "[ERR] Manifest step failed for tier=$* set=$(SET). See:"; \
	  echo "      $(MANIFEST_DIR)/manifest-$*-$(SET).log"; \
	  echo ""; \
	  tail -n 50 "$(MANIFEST_DIR)/manifest-$*-$(SET).log" || true; \
	  exit 1; \
	}

