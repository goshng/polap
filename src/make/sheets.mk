# make/sheets.mk — atlas sheets (pt & mt) from Bandage PNGs

SHEET_PT ?= $(OUTDIR)/figureS1-$(SET)-pt.pdf
SHEET_MT ?= $(OUTDIR)/figureS1-$(SET)-mt.pdf

.PHONY: sheets pngs
sheets: $(SHEET_PT) $(SHEET_MT)

pngs:
	@bash $(TOP)/polaplib/scripts/generate_pngs.sh --set "$(SET)" --type pt
	@bash $(TOP)/polaplib/scripts/generate_pngs.sh --set "$(SET)" --type mt

$(SHEET_PT): pngs
	@mkdir -p $(dir $@)
	@bash $(TOP)/polap-bash-figure-sheet.sh --set "$(SET)" --type pt --out "$@" --rows 5 --cols 4

$(SHEET_MT): pngs
	@mkdir -p $(dir $@)
	@bash $(TOP)/polap-bash-figure-sheet.sh --set "$(SET)" --type mt --out "$@" --rows 5 --cols 4
