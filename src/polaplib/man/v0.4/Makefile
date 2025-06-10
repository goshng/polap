CSL = molecular-biology-and-evolution-2.csl
BIBLIOGRAPHY = j-polap.bib

all: pdf

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
	pdfunite manuscript1.pdf figures/sheet_benchmark-some-2.csv.pdf manuscript.pdf

test:
	cat manuscript-supp.md \
		| sed 's/some-0.md/test-0.md/g' \
		| sed 's/some-1.md/test-1.md/g' \
		| sed 's/some-2.md/test-2.md/g' \
		| sed 's/some-3.md/test-3.md/g' \
		| sed 's/some-4.md/test-4.md/g' \
		| sed 's/Eucalyptus_pauciflora-0.md/Taxon_genus-0.md/g' \
		| sed 's/Eucalyptus_pauciflora-1.md/Taxon_genus-1.md/g' \
		| sed 's/Eucalyptus_pauciflora-2.md/Taxon_genus-2.md/g' \
		| sed 's/Eucalyptus_pauciflora-3.md/Taxon_genus-3.md/g' \
		| sed 's/Eucalyptus_pauciflora-4.md/Taxon_genus-4.md/g' \
		| sed 's/some-1.pdf/test-1.pdf/g' \
		| sed 's/some-2.pdf/test-2.pdf/g' \
		| sed 's/some-3.pdf/test-3.pdf/g' \
		| sed 's/some-4.pdf/test-4.pdf/g' \
		| sed 's/some-0.pdf/test-0.pdf/g'>expanded1.md
	bash replace_include.sh expanded1.md >expanded.md
	cat manuscript.md expanded.md > 1.md
	pandoc 1.md -s -f markdown -t pdf \
		--filter pandoc-crossref \
		--citeproc \
		--pdf-engine=xelatex \
		--metadata-file=j-polap.yaml \
		--csl=$(CSL) \
		-o manuscript1.pdf
	pdfunite manuscript1.pdf figures/sheet_benchmark-test-2.csv.pdf manuscript-test.pdf

clean:
	rm -f manuscript.aux manuscript.bbl manuscript.bcf manuscript.blg manuscript.docx
	rm -f manuscript.log manuscript.pdf manuscript.run.xml manuscript.tex
	rm -f supplementary.aux supplementary.bbl supplementary.bcf supplementary.blg 
	rm -f supplementary.log supplementary.pdf supplementary.run.xml supplementary.tex
	rm -f clean.yaml expanded.md j-polap.json
	rm -f 1.md
