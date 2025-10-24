# polaplib/make/pack.mk
# Version bump + tarball creation for polaplib
# NOTE: every recipe line MUST start with a TAB character.

# Run each recipe in a single shell (so heredoc works reliably)
.ONESHELL:

# Where to write/read version and output tarballs
VERSION_FILE ?= $(POLAPLIB_DIR)/VERSION
DIST_DIR     ?= $(POLAPLIB_DIR)/dist

.PHONY: version.bump pack tarball

version.bump:
	@# Create VERSION if missing, else bump patch x.y.z -> x.y.(z+1)
	if [ ! -f "$(VERSION_FILE)" ]; then
		echo "0.1.0" > "$(VERSION_FILE)"
		echo "[INFO] created $(VERSION_FILE) as 0.1.0"
	else
		V="$$(cat "$(VERSION_FILE)")"
		M="$${V%%.*}"
		R="$${V#*.}"
		m="$${R%%.*}"
		p="$${R#*.}"
		p="$$(($$p+1))"
		echo "$${M}.$${m}.$${p}" > "$(VERSION_FILE)"
		echo "[INFO] bumped version $$V -> $$(cat "$(VERSION_FILE)")"
	fi

pack:
	@mkdir -p "$(DIST_DIR)"
	tar -czf "$(DIST_DIR)/polaplib-$$(cat $(VERSION_FILE)).tar.gz" polaplib
	echo "[INFO] wrote $(DIST_DIR)/polaplib-$$(cat $(VERSION_FILE)).tar.gz"

tarball: pack
