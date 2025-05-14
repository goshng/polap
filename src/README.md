# 2025-05-13

Write a function for adding a subcommand

# 2025-04-30

## pmat run errors

fuse2fs not found
gocryptfs not found

## solution

intsall apptainer to the system-wide using apt install

# 2025-04-29

We need to limit polap bioconda package to a certain python version.

requirements:
host: - python >=3.8,<3.11
run: - python >=3.8,<3.11

# 2025-04-23

Polap has proven to be an invaluable tool for my research in developing the organelle pipeline, consistently delivering development process that have significantly impacted my project's progress.
I require renaming the existing file names to ensure consistency and clarity throughout the system.

- run-polap-function-directional.sh -> polap-subcmd-directional.sh
- run-polap-r-data-v2-alpha0.R -> polap-r-data-v2-alpha0.R
- run-polap-py-xyz.py -> polap-py-xyz.py
- polap-lib-timing.sh: keep it that way.
- polap-data-v2.sh: keep it that way.

# 2025-04-22

directional feature
run-polap-function-directional -> the subcommand
run-polap-function-dga -> functions defined and they are called in the subcommand
run-polap-function-nextdenovo : a new machine

# 2025-04-17

sequencing data types

PacBio HiFi data
PacBio CLR data

https://github.com/bichangwei/PMAT/issues/5
https://github.com/bichangwei/PMAT/issues/5#issuecomment-1722112839
Thank you very much for using PMAT. It is primarily based on highly accurate long-read data for subsequent assembly. When the input is ONT data, it is necessary to call Nextdenovo or Canu for error correction, and this step will take a lot of time, especially for Canu. If you have assembled your nuclear genome using ONT data, you can use the corrected data as the input for PMAT, and assemble its mitochondrial genome by the following command in the latest PMAT:
PMAT autoMito -i CORRECTED_ONT_DATA -o OUTPUT -st ONT -g 2.3G -tk p1 -cpu 32

https://nextdenovo.readthedocs.io/en/latest/QSTART.html#quick-start

# 2024-12-28

## polap-constants.sh

error code

# 2024-11-16

Need to clean up this README

# 2024-10-20

Learn something about Flye:
https://github.com/mikolmogorov/Flye/blob/flye/docs/USAGE.md

not-used-snippets/ : moved scripts not used or combined.

about log

log level 0: user interaction or main steps, e.g., skipping ...
log level 1: I/O for the main steps
log level 2: details of I/O
log level 3: function-level and command executions

# 2024-10-13

## Sequence length from Fasta

https://www.biostars.org/p/118954/

# 2024-10-02

## Compare POLAP vith PMAT

conda install conda-forge::apptainer
conda install bioconda::canu
conda install bioconda::nextdenovo

https://nextdenovo.readthedocs.io/en/latest/OPTION.html

# 2024-05-01 (roughly)

## Semi-auto seed contig selection

conda install networkx
conda install gftools
conda install pandas
conda install r::r-mixr
