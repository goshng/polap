# polaplib/make/bandage.mk
# Converts any *.gfa (or .gfa.gz) into a PNG next to it via polap.sh bandage
# Usage: make path/to/graph.png

%.png: %.gfa
	@mkdir -p $(dir $@)
	@if command -v polap.sh >/dev/null 2>&1; then \
	  polap.sh bandage png "$<" "$@" || cp -f $(POLAPLIB_DIR)/man/figures/na.png "$@"; \
	else \
	  echo "[WARN] polap.sh not found; copying NA image"; \
	  cp -f $(POLAPLIB_DIR)/man/figures/na.png "$@"; \
	fi

%.png: %.gfa.gz
	@mkdir -p $(dir $@)
	@if command -v polap.sh >/dev/null 2>&1; then \
	  gunzip -c "$<" > "$<.tmp.gfa" && polap.sh bandage png "$<.tmp.gfa" "$@" || cp -f $(POLAPLIB_DIR)/man/figures/na.png "$@"; \
	  rm -f "$<.tmp.gfa"; \
	else \
	  echo "[WARN] polap.sh not found; copying NA image"; \
	  cp -f $(POLAPLIB_DIR)/man/figures/na.png "$@"; \
	fi
