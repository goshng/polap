CSL = molecular-biology-and-evolution-2.csl
BIBLIOGRAPHY = j-polap.bib
FIGURES= \
mainfigure3-2.md \
maintable1-0-analysis.md \
maintable1-1-analysis.md \
maintable1-2-analysis.md \
maintable3-2-analysis.md \
maintable1-2-data.md \
maintable1-3-analysis.md \
maintable1-4-analysis.md \
suppfigure1-2-1-Eucalyptus_pauciflora.md \
suppfigure3-2.md \
supptable1-0-x-Eucalyptus_pauciflora.md \
supptable1-2-x-Anthoceros_agrestis.md \
supptable1-2-x-Codonopsis_lanceolata.md \
supptable1-2-x-Eucalyptus_pauciflora.md \
supptable1-2-x-Juncus_roemerianus.md \
supptable1-3-x-Eucalyptus_pauciflora.md \
alpha0.pdf \
delta.pdf

FIGURE_PATHS = $(addprefix figures/,$(FIGURES))

all: pdf

bibsync:
	bash sync_bib.sh

update:
	unison -auto -batch manuscript

backup:
	cp manuscript.md manuscript-`date +%m%d%H`.md

docx:
	bash replace_include.sh manuscript-supp.md >expanded.md
	cat manuscript.md expanded.md | \
	pandoc -s -f markdown -t docx \
		--filter pandoc-crossref \
		--lua-filter toggle-section.lua \
		--citeproc \
		--bibliography=$(BIBLIOGRAPHY) \
		--csl=$(CSL) \
		-o manuscript.docx

ref:
	pandoc --from=bibtex --to=csljson j-polap.bib -o references.yaml
	yq -p=json -o=yaml '{"references": .}' references.yaml > j-polap.yaml

pdf:
	bash replace_include.sh manuscript-supp.md >expanded.md
	cat manuscript.md expanded.md > 1.md
	pandoc 1.md -s -f markdown -t pdf \
		--filter pandoc-crossref \
		--citeproc \
		--pdf-engine=xelatex \
		--metadata-file=j-polap.yaml \
		--csl=$(CSL) \
		-o manuscript1.pdf
	# pdfunite manuscript1.pdf figures/sheet_benchmark-2.csv.pdf manuscript.pdf

tex:
	bash replace_include.sh manuscript-supp.md >expanded.md
	cat manuscript.md expanded.md | \
	pandoc -s -f markdown -t latex \
		--filter pandoc-crossref \
		--citeproc \
		--pdf-engine=xelatex \
		--bibliography=$(BIBLIOGRAPHY) \
		--csl=$(CSL) \
		-o manuscript.tex

cover:
	pandoc cover-letter.md -s -f markdown -t pdf \
		--filter pandoc-crossref \
		--lua-filter toggle-section.lua \
		--citeproc \
		--pdf-engine=xelatex \
		--metadata-file=j-polap.yaml \
		--csl=$(CSL) \
		-o cover-letter.pdf

print:
	cp manuscript.pdf backup/manuscript-`date +%m%d%H%M%S`.pdf
	cp manuscript.pdf /Users/goshng/Dropbox/Apps/onyx-knote/print/manuscript-`date +%m%d%H%M%S`.pdf
	# cp manuscript.docx /Users/goshng/Dropbox/Apps/onyx-knote/print/manuscript-`date +%m%d%H`.docx

manuscript:
	unison /Users/goshng/all/polap ssh://thorne/all/polap -ignore 'Name mtest' -ignore 'Name o' -ignore 'Name files' -ignore 'Name .DS_Store' -ignore 'Name .git' -ignore 'Name test' -ignore 'Name .*' -ignore 'Name run' -auto -silent
	cp /Users/goshng/Dropbox/Documents/notes/manuscript/00-Bibliography/j-polap.bib .
	cp /Users/goshng/all/notes/note/noexport-subtrees.lua .
	pandoc manuscript.md -s -f markdown -t plain \
		--lua-filter=noexport-subtrees.lua \
		--citeproc \
		--pdf-engine=xelatex \
		--bibliography=j-polap.bib \
		--csl=american-journal-of-botany.csl \
		--wrap=preserve \
		-o manuscript.txt
	pandoc manuscript.md -s -f markdown -t pdf \
		--lua-filter=noexport-subtrees.lua \
		--citeproc \
		--pdf-engine=xelatex \
		--bibliography=j-polap.bib \
		--csl=american-journal-of-botany.csl \
		-o manuscript.pdf
	cp manuscript.pdf /Users/goshng/Dropbox/Apps/onyx-knote/print/manuscript-`date +%m%d%H`.pdf
	unison /Users/goshng/all/polap ssh://thorne/all/polap -ignore 'Name mtest' -ignore 'Name o' -ignore 'Name files' -ignore 'Name .DS_Store' -ignore 'Name .git' -ignore 'Name test' -ignore 'Name .*' -ignore 'Name run' -auto -silent
	# open manuscript.pdf


figures:
	cp $(FIGURE_PATHS) .

mv:
	mv main*.md sup*.md *.pdf table*.md manuscript-????????.* t

clean:
	rm -f manuscript.aux manuscript.bbl manuscript.bcf manuscript.blg manuscript.docx
	rm -f manuscript.log manuscript.pdf manuscript.run.xml manuscript.tex
	rm -f supplementary.aux supplementary.bbl supplementary.bcf supplementary.blg 
	rm -f supplementary.log supplementary.pdf supplementary.run.xml supplementary.tex
	rm -f clean.yaml expanded.md j-polap.json
