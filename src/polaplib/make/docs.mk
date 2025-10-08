# polaplib/make/docs.mk
DOCS_DIR ?= $(POLAPLIB_DIR)/docs
QMD_MAIN ?= $(DOCS_DIR)/index.qmd
DOC_PDF  ?= $(DOCS_DIR)/polaplib-docs.pdf

.PHONY: docs docs.clean
docs:
	@if command -v quarto >/dev/null 2>&1; then \
	  quarto render "$(QMD_MAIN)" --to pdf; \
	  echo "[INFO] Wrote docs under $(DOCS_DIR)"; \
	else \
	  echo "[WARN] Quarto not found; skipping docs"; \
	fi

docs.clean:
	@rm -rf $(DOCS_DIR)/_site $(DOCS_DIR)/*.pdf
