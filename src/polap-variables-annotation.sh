local _polap_var_wga="${ODIR}/0"
local _polap_var_wga_contigger="${_polap_var_wga}/30-contigger"
local _polap_var_wga_contigger_gfa="${_polap_var_wga_contigger}/graph_final.gfa"
local _polap_var_wga_annotation="${_polap_var_wga}/assembly_info_organelle_annotation_count-all.txt"

ANUM=$INUM

local _polap_var_ga="$ODIR"/$INUM

assembly_contigs_stats="$FDIR"/30-contigger/contigs_stats.txt
assembly_contigs_fasta="$FDIR"/30-contigger/contigs.fasta

local _polap_var_ann="${_polap_var_ga}"/50-annotation
local _polap_var_ann_CONTIGFILE="${_polap_var_ann}"/contig.fasta
local _polap_var_ann_CONTIGDB="${_polap_var_ann}"/contig
local _polap_var_ann_MTAABLAST="${_polap_var_ann}"/mtaa.blast
local _polap_var_ann_MTAABED="${_polap_var_ann}"/mtaa.bed
local _polap_var_ann_MTGENECOUNT="${_polap_var_ann}"/mt.gene.count
local _polap_var_ann_PTAABLAST="${_polap_var_ann}"/ptaa.blast
local _polap_var_ann_PTAABED="${_polap_var_ann}"/ptaa.bed
local _polap_var_ann_PTGENECOUNT="${_polap_var_ann}"/pt.gene.count
local _polap_var_ann_CONTIGNAME="${_polap_var_ann}"/contig.name
