# How-To Add Polap commands

Copy the polap-cmd-template.sh to add a new polap command.

# Bolap man

2025-12-01

A few reporting tables and figures are generated.
Bolap analysis types each require a specific set of source files, including `Makefile.type` and `polap-data-type.sh`, where _type_ can be hifi, read, aflye, or cflye.
None of these analysis types have been thoroughly tested because Bolap has not yet been released.
Each type is associated with a particular release version: aflye corresponds to v0.3, cflye to v0.4, hifi to v0.5, and read is planned for v0.6.
I previously worked on `polap-data-read.sh` for the Bolap command, but I realized that completing this component would require more development time.
Therefore, I have returned to the hifi analysis path and am focusing on reviewing the recently released plant mitogenome pipelines, including PMAT2, TIPPo, Oatk, and HiMT.

# manifest for reporting

We use polap-bash-make-manifest.sh to collect results in a structured json file.

First, we list all files prepared in a species folder using:

bash polap/src/polaplib/polap-dev-list-files.sh Breynia_androgyna >polap/src/polaplib/report-files-hifi.txt

Then, edit the report-files-...txt file.

Edit polap-bash-make-manifest.sh according to the report-files-...txt document.

We have polap-bash-hifi-make-manifest.sh and polap-bash-read-make-manifest.sh for each bolap analysis type because results would be different.

# About this polaplib

The `src` folder has the main bash shell scripts including `polap.sh` and
`polap-data-v2.sh`. The `polaplib` folder contains other scripts including those
in bash, R, or python. It could have other scripts as well as long as they can
be executed by bash shell.

# Testing code

## python and VScode

https://stevenmortimer.com/setting-up-vs-code-for-python-development-like-rstudio/
Backup these

- user-settings.json
- keybindings.json

And replace them with ones in the site above. This allows us to use VSCode more
like a RStudio-like environment. We can execute a section of code to test, which
is handy to use.

Use Remote-SSH to connect a Linux computer using VSCode and then install
extentions.

See Also:

- hello-world.py
- https://www.youtube.com/watch?v=PwGKhvqJCQM

## bash

Use `bats` to test bash scripts.

# Modules

Polap has grown to include different parts including the original plant
mitochondrial genome assembly pipeline.

- plant mitochondrial genome assembly using ONT and Illumina data
- plant plastid genome assembly using ONT and Illumina data

We are considering other modules:

- taxonomy: a pipeline to infer trees using organelle genome
- directional: a pipeline to use seeds with direction

# Files

## polap data files

Change `_POLAP_LIB` to use polap-test.sh to execute some of the scripts.
polap-test.sh
input
output

## man

- man/v0.2: manual for polap assemble (see v0.4)
- man/v0.3: manual for polap taxonomy (not yet)
- man/v0.4: manual for polap disassemble

## git commit

We used to append the commit hash to the version so that we could tract the code
version. We now use post-commit in the git directiory: .git/hooks.
These are just backups of those because it is stored per git folder.

- git_hooks_post-commit
- git_hooks_pre-commit

## Figures and Tables

We use the TSV file format to create tables using `pandoc`. If this is hard to
use, we use LaTeX to create tables. For figures, we mostly use R scripts with
`tidyverse`-family R packages including `dplyr`, `readr`, and the like.

- polap-bash-figure-latex.sh -> v0.4 LaTeX for plastid genome assembly graph figures
- polap-bash-report-table1.sh -> v0.2 main Table 1

## Polap build and bioconda

Bioconda is the important part of polap. We have two files for bioconda
installation.

- meta.yaml
- build.sh

We used to use `polap-build.sh` to create `build.sh` for including new scripts
in the installation. Now, we have `polaplib` folder for script library. We may
not need `polap-build.sh`.

We also have conda yaml files that we used before. Keep these for the record.

- polap-conda-environment-fmlrc.yaml
- polap-conda-environment.yaml

## Name changes

We used to prepend `run` to script file names. We may want to remove such `run`
from the script file names.

- run-polap-function-xxx.sh -> polap-cmd-xxx.sh
- run-polap-py-xxx.sh -> polap-py-xxx.sh
- run-polap-bash-xxx.sh -> polap-bash-xxx.sh
- run-polap-r-xxx.sh -> polap-r-xxx.sh

We have other bash shell scripts with the following names:

- polap-function-xxx.sh
- polap-lib-xxx.sh

Both of them have shell script functions used in polap-cmd-xxx.sh or
run-polap-function-xxx.sh.
More general library files are those with `polap-lib-xxx.sh` and some of
functions extracted from `polap-cmd-xxx.sh` are those with
`polap-function-xxx.sh`.

## Cammond-line completion

This could be useful so that users are forced to use specific options and files.
We have not implemented this yet because it is experimental.

- polap-command-completion.sh
- polap-data-cflye.complete.sh -> polap-command-completion-polap-data-cflye.sh

## polap bash shell script global variables

There are variables and functions often used in polap's script.

- polap-constant.sh

## CSV data files

Data analysis using polap needs configurations including file paths and options
used. These are set in CSV-format files:

- polap-data-aflye.csv
- polap-data-cflye.csv
- polap-data-dflye.csv
- polap-data-taxon.csv
- polap-data-v1.csv
- polap-data-v2.csv
- polap-data-v3.csv
- polap-data-v4.csv
- polap-data-vt.csv

Keep this for a simple data format.

- polap-data-v2.data

## polap data analysis

We currently use `polap.sh` and `polap-data-v2.sh`. The others are not actively
used for polap disassemble subcommand. However, we would use `polap-data-v1.sh`
for the v0.3.7.3 polap. We also use `polap-data-menu.sh` to add subcommands to
data analysis scripts such as `polap-data-v2.sh` and `polaplib/polap-lib-data.sh`.

- bolap.sh -> Benchmarking script (not developed yet)
- polap-batch-v0.sh -> template batch script
- polap-batch-v2.sh -> a batch script (not used due to polap-data-v2.sh)
- polap-data-menu.sh -> dev script to add a new subcommand to data analysis
- polap-data-v0.sh -> data analysis template script
- polap-data-v1.sh -> v0.3.7.3 data analysis
- polap-data-v2.sh -> v0.4.3 data analysis
- polap-data-v3.sh -> taxonomy module data analysis
- polap-data-v4.sh -> directional module data analysis
- polap.sh -> the main polap script

## polap test scripts

A shell script that tests R, python, bash scripts in polaplib folder. It needs
`input` and `output` folders.

- polap-bash-test.sh

## polap bash variables

A polap is a giant bash shell script with global variables. So, we need a naming
scheme so that we can access those variables from other functions or scripts.
Files are named as systematically as possible. They are placed in either
`polap-variables-common.sh` or `polap-variables-mtcontigs.sh`.

- polap-variables-common.sh
- polap-variables-mtcontigs.sh

Before it enters the main part of the polap script, it sources
`polap-variables-main.sh` so that some preconditions are met. Edit this
`polap-variables-main.sh` script to change that behavior such as creating
default output folders.

- polap-variables-main.sh

## polap-lib-xxx

These more general bash shell script libraries can be used.

Function names can be prepended by the library script name: e.g.,

```bash
_polap_lib_array-csv_to_array() {
```

- polap-lib-array.sh
- polap-lib-conda.sh
- polap-lib-config.sh
- polap-lib-data.sh
- polap-lib-debug.sh
- polap-lib-errors.sh
- polap-lib-extract.sh
- polap-lib-fastq.sh
- polap-lib-filepath.sh
- polap-lib-file.sh
- polap-lib-log.sh
- polap-lib-ncbi.sh
- polap-lib-number.sh
- polap-lib-process.sh
- polap-lib-random.sh
- polap-lib-steps.sh
- polap-lib-timing.sh
- polap-lib-unit.sh

## polap data files for plant organelle gene annotation

polap-mt.1.c70.3.faa
polap-mt.1.c70.3.faa.name
polap-pt.2.c70.3.faa
polap-pt.2.c70.3.faa.name
polap-mtgenes.txt

## polap command-line parsing

- polap-parsing.sh

## External modified scripts

We had to keep it because of the genome size option that is not fed into Flye's
execution. The following script has the genome size option that is used by Flye.

- polap-ptGAUL1.sh

## Python scripts

Python scripts are used for finding connected components, extracting plastid
genome sequences, unique MT/PT contig names, compare two ptDNA sequences from
GFA.

Note: we will rename run-polap-py-xxx.py as polap-py-xxx.py.

- hello-world.py
- polap-py-compare2ptdna.py
- polap-py-find-cc-with-seeds.py
- polap-py-find-plastid-gfa2fasta.py
- polap-py-unique-mtcontigs.py
- polap-py-compare2ptdna.py
- polap-py-find-cc.py
- polap-py-find-cc-with-seeds.py
- polap-py-find-plastid-gfa2fasta.py
- polap-py-unique-mtcontigs.py

## v0.4 figures

We convert a TSV file named table-some-2.tsv to boxplots for displaying time,
memory, and time with nextdenovo-polishing.

- polaplib/polap-r-disassemble-man-benchmark-boxplots.R

## Templates

We have template scripts that we copy to start coding scripts.

- polap-cmd-template.sh
- polap-r-template.R
- run-polap-function-template.sh
- run-polap-r-template.R

We may need a python template script.

- polap-py-template.py

## Playground scripts

We may have scripts for testing any codes.

- polap-playground.sh
- polap-playground.py
- polap-playground.R
- polap-playground.c
- polap-playground.cpp
- polap-playground.java
- polap-playground.js

## Text and LaTeX files

We could use the polap data analysis command to create figures in PDF format:

- polap-data-cflye man figure-sheet some 2 bandage # Requires latex
- polap-data-cflye man figure-sheet-pmat some # Requires latex
- polap-data-cflye man figure-sheet-oatk some # Requires latex
- polap-data-cflye man figure-sheet-tippo some # Requires latex
- polap-data-cflye man pdf # Requires latex

It converts data files to a LaTeX file, which turns into a PDF file. Some of the
following text file contents are added to the LaTeX file.

- polap_cflye_sheet_benchmark.txt
- polap_cflye_sheet_oatk.txt
- polap_cflye_sheet_pmat.txt
- polap_cflye_sheet_polap-0.txt
- polap_cflye_sheet_polap-1.txt
- polap_cflye_sheet_polap-2.txt
- polap_cflye_sheet_polap-4.txt
- polap_cflye_sheet_polap.txt
- polap_cflye_sheet_tippo.txt

## NCBI SRA and BioProject

We start with the NCBI's SRA database to collect data.

```bash
polap download-runinfo-bioproject <bioproject-accessions.txt>
```

where `bioproject-accessions.txt` has BioProject accession ID one per line.
It downloads runinfo files from the NCBI. Search the NCBI for a list of
BioProject accessions.

e.g.,

- NCBI:BioProject:"Viridiplantae"[Organism] AND oxford[All Fields]
- NCBI:BioProject:"Viridiplantae"[Organism] AND pacbio[All Fields]
- NCBI:BioProject:"animalia"[Organism] AND oxford[All Fields]

polap's get-bioproject would download the runinfo of a BioProject ID.
So, we need another subcommand to parse a runinfo to extract what we want.
A BioProject could have multiple species data. So, parse a runinfo to
list sequencing data per species. See, the get-bioproject subcommand and
`polap-r-get-bioproject.R` script.

## Pairwise alignments

- polaplib/polap-cmd-mtdna.sh

## Demo

We rename this:

- polaplib/run-polap-function-demo.sh -> polap-cmd-demo.sh

However, `polap-data-v1.sh` or `polap-data-aflye` replaces it.

## Menus or subcommands

We use `list` subcommands to list them.

## module mtdna

The mtdna module has somewhat mixed features implemented in a single bash
script.

- polap-cmd-mtdna.sh

## Other subcommands

Subcommands or functions that are not hosted in any other bash scripts.

- polaplib/run-polap-function-miscellaneous.sh
- polaplib/run-polap-function-utilities.sh

## main modules

- assemble
- disassemble
- directional
- taxonomy

## files that could be deleted in the future but not yet

polaplib/polap-package-common.sh
polaplib/polap-package-mtcontigs.sh
polaplib/run-polap-function-annotate-contig.sh

## user subcommands

We need to document the help messages for these subcommands:

- init
- list
- make-menus
- clean-menus
- annotate
- archive
- assemble
- assemble1
- assemble2
- bandage
- blast-genome
- cleanup
- copy-sra-bioproject
- count-gene
- depth-distribution
- disassemble
- download-runinfo-bioproject
- edges-stats
- fastq
- find-genome-size
- flye1
- flye2
- flye-polishing
- get-bioproject
- get-bioproject-sra
- get-dna-by-accession
- get-mtdna
- get-taxonomy-species
- install
- log
- mafft-mtdna
- map-reads
- menu
- polish
- polish-disassemble
- prepare-polishing
- prepare-seeds
- reduce-data
- report-assembly
- report-table
- seeds
- select-reads
- summary-reads
- test-reads
- total-length-long
- total-length-short
- compare2ptdna
- compare-mtdna
- version

## subcommands for development and testing

- get-revision1
- mauve-mtdna
- make-menus-all
- demo
- assemble-bioproject
- assemble-wrange
- rng
- test
- choose-seed
- test-steps
- prepare-polishing_sort_efficient
- polish_v2
- prepare-polishing_v2
- source-menus
- step-disassemble
- step-disassemble-archive-cfile
- step-disassemble-seeds-graph
- bioproject-postprocess
- bioproject-prepare
- archive-rsync-template
- seeds-graph
- ptgaul
- template

## subcommands that will be deprecated

- annotate-contig
- count-gene-contig
- blast-genome-contig
