# polaplib/Makefile
# Minimal, self-contained. Tests that `make -C polaplib version.bump` works.

VERSION_FILE ?= $(CURDIR)/VERSION
DIST_DIR     ?= $(CURDIR)/dist

.PHONY: version.bump pack echo-vars

version.bump:
	@# Create VERSION if missing, else bump patch
	@if [ ! -f "$(VERSION_FILE)" ]; then \
		echo "0.1.0" > "$(VERSION_FILE)"; \
		echo "[INFO] created $(VERSION_FILE) as 0.1.0"; \
	else \
		V=$$(cat "$(VERSION_FILE)"); \
		M=$${V%%.*}; R=$${V#*.}; m=$${R%%.*}; p=$${R#*.}; \
		p=$$((p+1)); \
		echo "$${M}.$${m}.$${p}" > "$(VERSION_FILE)"; \
		echo "[INFO] bumped version $$V -> $$(cat "$(VERSION_FILE)")"; \
	fi
	@echo "[OK] $(VERSION_FILE) = $$(cat "$(VERSION_FILE)")"

pack:
	@mkdir -p "$(DIST_DIR)"
	@tar -czf "$(DIST_DIR)/polaplib-$$(cat "$(VERSION_FILE)").tar.gz" ./
	@echo "[OK] wrote $(DIST_DIR)/polaplib-$$(cat "$(VERSION_FILE)").tar.gz"

echo-vars:
	@echo "CURDIR=$(CURDIR)"
	@echo "VERSION_FILE=$(VERSION_FILE)"
	@echo "DIST_DIR=$(DIST_DIR)"
