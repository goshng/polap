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

# FIXME

╰─ grep \_polap_var_annotation_table polap-variables-\*.sh ─╯
polap-variables-mtcontig.sh:local \_polap_var_annotation_table="${_polap_var_ga}/assembly_info_organelle_annotation_count-all.txt"
polap-variables-mtdna.sh:local _polap_var_annotation_table="${\_polap_var_ga}/assembly_info_organelle_annotation_count-all.txt"

╰─ grep \_polap_var_assembly_graph_final_gfa polap-variables-\*.sh ─╯
polap-variables-ga.sh:local \_polap_var_assembly_graph_final_gfa="${_polap_var_contigger}/graph_final.gfa"
polap-variables-mtcontig.sh:local _polap_var_assembly_graph_final_gfa="${\_polap_var_contigger}/graph_final.gfa"

╰─ grep \_polap_var_wga_annotation polap-variables-\*.sh ─╯
polap-variables-ga.sh:local \_polap_var_wga_annotation="${_polap_var_wga}/assembly_info_organelle_annotation_count-all.txt"
polap-variables-wga.sh:local _polap_var_wga_annotation="${\_polap_var_wga}/assembly_info_organelle_annotation_count-all.txt"

╰─ grep \_polap_var_wga_contigger polap-variables-\*.sh ─╯
polap-variables-ga.sh:local \_polap_var_wga_contigger="${_polap_var_wga}/30-contigger"
polap-variables-ga.sh:local _polap_var_wga_contigger_gfa="${\_polap_var_wga_contigger}/graph_final.gfa"
polap-variables-wga.sh:local \_polap_var_wga_contigger="${_polap_var_wga}/30-contigger"
polap-variables-wga.sh:local _polap_var_wga_contigger_gfa="${\_polap_var_wga_contigger}/graph_final.gfa"
polap-variables-wga.sh:local \_polap_var_wga_contigger_contigs_stats="${_polap_var_wga_contigger}/contigs_stats.txt"
polap-variables-wga.sh:local _polap_var_wga_contigger_contigs_fasta="${\_polap_var_wga_contigger}/contigs.fasta"
