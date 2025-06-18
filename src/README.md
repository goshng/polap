# polap_disassemble-seeds

function \_disassemble-step12 {

input files:
local \_contigger_edges_gfa="${_outdir}/30-contigger/graph_final.gfa"
	local _contigger_edges_fasta="${\_outdir}/30-contigger/graph_final.fasta"
local \_ga_annotation_all="${\_outdir}/assembly_info_organelle_annotation_count-all.txt"

output file:
local \_mtcontigname="${\_outdir}/mt.contig.name"

    		polap_disassemble-seeds "${_contigger_edges_gfa}" \
    			"${_ga_annotation_all}" \
    			"${_mtcontigname}"

polap_disassemble-seeds() {

    	_run_polap_step-disassemble-seeds-graph \
    _polap_log3_pipe "python ${_POLAPLIB_DIR}/run-polap-py-unique-mtcontigs.py \

function \_run_polap_step-disassemble-seeds-graph {

    	_polap_log1 "  step 1: clean-up the mtcontigs folder: ${_mtcontigs}"
    	_polap_log3_cmd rm -rf "${_mtcontigs}"
    	_polap_log3_cmd mkdir -p "${_mtcontigs}"

    	_polap_log1 "  step 2: determine the depth range using length CDF"
    	_polap_disassemble_seeds_determine-depth-range \
    		"${_ga_annotation_all}" \
    		"${_knum}" \
    		"${_mtcontigs_1_custom_depth_range}"

    	_polap_log1 "  step 3: pre-select contigs based on organelle gene annotation"
    	_polap_log2 "    input1: ${_ga_annotation_all}"
    	_polap_log2 "    input2: ${_mtcontigs_depth_range_preselection}"
    	_polap_log2 "    output: ${_mtcontigs_preselection}"
    	_polap_disassemble_seeds_preselect-contigs \
    		"${_ga_annotation_all}" \
    		"${_mtcontigs_depth_range_preselection}" \
    		"${_mtcontigs_preselection}"

    	_polap_log1 "  step 4: filter GFA by the depth range"
    	_polap_log2 "    input1: ${_ga_contigger_edges_gfa}"
    	_polap_log2 "    input2: ${_mtcontigs_depth_range_preselection}"
    	_polap_log2 "    output: ${_mtcontigs_gfa_depthfiltered_gfa}"
    	_polap_disassemble_seeds_depthfilter-gfa \
    		"${_ga_contigger_edges_gfa}" \
    		"${_mtcontigs_depth_range_preselection}" \
    		"${_mtcontigs_gfa_depthfiltered_gfa}"

    	_polap_log1 "  step 5: find connected components with the preselected contigs"
    	_polap_log2 "    input1: ${_mtcontigs_gfa_depthfiltered_gfa}"
    	_polap_log3_pipe "python ${_POLAPLIB_DIR}/run-polap-py-find-cc-with-seeds.py \
      --gfa ${_mtcontigs_gfa_depthfiltered_gfa} \
    		--nodes ${_mtcontigs_preselection} \
      --output ${_mtcontigs_gfa_depthfiltered_cc_seed}"

\_polap_disassemble_seeds_depthfilter-gfa() {

    _polap_log1 "    step 4-1: create GFA without sequence data using gfatools view"

    _polap_log1 "    step 4-2: extracte sequence part of GFA: ${_mtcontigs_gfa_seq_part}"

    _polap_log1 "    step 4-3: filter GFA sequence part using depth range"

    _polap_log1 "    step 4-4: subsetting GFA using the depth-filtered GFA sequence part with gfatools view"

# Files need attention

1. polaplib/run-polap-r-data-v2-alpha0.R - command-line processing: use argparser

2. simplify run-polap-function-archive.sh: this would remove the followings:

- polaplib/run-polap-r-contig2edge.R

look for this to search for potentially deleted ones.
TAG: SCC delete

# Add a new subcommand to polap

How to create a polap menu or subcommand:
polap works with a subcommand. If you want to add a hello-world subcommand.

1. Copy polaplib/polap-cmd-template.sh to create polaplib/polap-cmd-hello-world.sh
2. Replace the function named \_run_polap_template with \_run_polap_hello-world.
3. Find the line of 'ADD A NEW SUBCOMMAND SOURCE HERE' and add this line.
   source "${\_POLAPLIB_DIR}/polap-cmd-hello-world.sh"
4. A subcommand has 'help' menu; 'polap hello-world help' should display a message.
5. A subcommand can be create using \_arg_menu[1] variable.

# Script file structure

polap.sh
polaplib/polap-parsing.sh - commmand options
polaplib/polap-variables-main.sh - variables defined before entering the subcommand

polap-lib-xxx.sh - polap bash library
polap-cmd-xxx.sh - subcommand scripts (New)
run-polap-function-xxx.sh - subcommand scripts (Old)

## Template files

polap-cmd-template.sh
polap-r-template.sh

# 2025-06-15

conda build?
(bioconda) $ bioconda-utils build --packages dflye
anaconda login (if needed)
anaconda upload /home/goshng/miniconda3/envs/bioconda/conda-bld/linux-64

# 2025-06-13

To add one to polaplib/polap-lib-data.sh as a common subcommand to all polap-data-vX.sh:
bash polap-data-menu.sh subcommand-name
To add one to polap-data-v4.sh as a subcommand:
bash polap-data-menu.sh polap-data-v4.sh subcommand-name

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
