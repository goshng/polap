all: bioconda

update2proj:
	cp -u -p src/polap.sh src/polap
	cp -u -p src/* ../proj/src/

update4proj:
	cp -u -p ../proj/src/polap src/
	cp -u -p ../proj/src/polap-conda-environment-fmlrc.yaml src/
	cp -u -p ../proj/src/polap-conda-environment.yaml src/
	cp -u -p ../proj/src/polap-mt.1.c70.3.faa src/
	cp -u -p ../proj/src/polap-parsing.sh src/
	cp -u -p ../proj/src/polap-pt.2.c70.3.faa src/
	cp -u -p ../proj/src/polap.sh src/
	cp -u -p ../proj/src/run-polap-genes.R src/
	cp -u -p ../proj/src/run-polap-jellyfish.R src/
	cp -u -p ../proj/src/run-polap-mtcontig.R src/
	cp -u -p ../proj/src/run-polap-pairs.R src/

docs:
	pandoc --standalone -t html man.md -o man.html
	pandoc --standalone -t html README.md -o README.html

bioconda:
	cp src/bioconda-recipes-polap/* ~/all/tools/conda/bioconda-recipes/recipes/polap/
	cp src/bioconda-recipes-polap-fmlrc/* ~/all/tools/conda/bioconda-recipes/recipes/polap-fmlrc/

back:
	cp ~/all/tools/conda/bioconda-recipes/recipes/polap/* src/bioconda-recipes-polap/
	cp ~/all/tools/conda/bioconda-recipes/recipes/polap-fmlrc/* src/bioconda-recipes-polap-fmlrc/

clean:
	rm -f man.html README.pdf README.txt README.html
