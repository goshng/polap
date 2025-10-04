---
title: "Polap Documentation"
version: "0.5.3.1"
header-includes: |
  \usepackage{gensymb}
  \usepackage[left]{lineno}
  <!-- \linenumbers -->
  \raggedright
  \usepackage[backend=biber,style=authoryear]{biblatex}
  \setlength{\parskip}{1.2\baselineskip}
  \usepackage{pdflscape}
  \newcommand{\blandscape}{\begin{landscape}}
  \newcommand{\elandscape}{\end{landscape}}
  \usepackage{graphicx}
  \usepackage{array}
  \usepackage{booktabs}
  \usepackage{caption}
  \renewcommand{\arraystretch}{0.5}
  \setlength{\tabcolsep}{2pt}
  \AtBeginEnvironment{tabular}{\tiny}
  \AtBeginEnvironment{longtable}{\scriptsize}
  <!-- \widowpenalties 1 10000 -->
  <!-- \raggedbottom -->
classoption:
  - 12pt
geometry:
  - top=25mm
  - left=25mm
linestretch: 2
indent: true
colorlinks: urlcolor, citecolor, toccolor, linkcolor
lof: false
lot: false
thanks: false
toc: false
chapters: True
chaptersDepth: 1
chapDelim: ""
figureTitle: |
  Figure
figPrefix:
  - "Figure"
  - "Figures"
tableTitle: |
  Table
tblPrefix:
  - "Table"
  - "Tables"
showDetails: false
---

<!-- https://github.com/jgm/pandoc/issues/3148 -->

# Title

Plastid genome assembly by read selection using homology search

Address

\newpage

# Abstract {-}

Content of abstract.

\newpage

# Introduction {-}

Plastid genoms can be assembled using either short- or long-read sequencing data.

`GetOrganelle` [@Jin2020] pioneered plastid genome assembly pipelines using the assembler
`SPAdes` [@Bankevich2012].
They use a dataset of plastid genes to select short reads for plastid genome assembly.

`ptGAUL` [@Zhou2023] uses known but closely-related reference plastid genome to recruit long reads for plastid genome assembly using `Flye` genome assembler [@Kolmogorov2019].
Similarly, `CLAW` [@Phillips2024] devoleped a pipeline using a reference sequences.

Reference sequences could be generated to be used as seeds for recruiting reads potentially originated in a target genome.
Although `PMAT` [@Bi2024] targets on plant mitochondrial genome assemblies, one could generate such plastid genome reference from raw sequencing data to use as seeds.

`TIPPo` [@Xian2025] employed `Tiara` [@Karlicki2022] to select plant organelle derived reads for plastid and mitochondrial genome assemblies.
Instead of aligning sequencing reads on a reference sequence dataset, they directly discriminate plant plastid and mitochondrial from nuclear derived reads.
`Tiara` compiled a set of reference sequences from various species across wide taxa to train a neural netwrk model to predict the origin of sequences.

`Oatk` [@Zhou2024] uses an efficient assembler `syncasm` specialized in organelle genome assembly.

`HiMT` [@Tang2025] provide a user-friendly tool for plant organelle genome assemblies.

Some of these tools require highly accurate sequence quality like PacBio HiFi data.
`HiMT` uses k-mer to estimate depths of reads to select plant organelle reads.
`Oatk` needs very long k-mer to build sparse deBruijn graph for assembling organelle genomes quickly.
`TIPPo` depends on `Tiara`, which use a very short k-mer to make input data its neural network model.
`PMAT` is more efficient for PacBio HiFi data and can be used with low-quality ONT data because the low-quality data are better-off to be error-corrected before using as input.

Although PacBio HiFi (HiFi) sequencing data take the main stream in the sequencing field, Oxford Nanopore (ONT) sequencing is also advancing.
Even though ONT data quality is improving, it was still lower than PacBio HiFi.

For low-quality data, most of the k-mer based method are not efficient or not working properly in selecting long reads for plant organelle genome assembly.
Reference sequences could be generated from long-read data and be chosen for seeds for organelle genome just like `PMAT` approach.
Furthermore, plastid genomes occur in very high copy number and genome assembly of long-read data can produce plastid genome.
However, whole-genome sequencing would take longer than organelle genome assembly because input data consist of mostly nuclear DNA molecules.
This is the reason why `TIPPo` annd `HiMT` were developed using `Flye` as their backend assembler.
These two tools are very efficient to select plant organelle derived long reads especially PacBio HiFi reads.

These two tools could be used for selecting ONT long-reads if ONT long-read quality is improving.
However, if the long-read quality is not high enough, the selecting approach is not applicable.
For HiFi sequencing data, `Oatk` with the assembler `syncasm` could be more efficient.
So, `TIPPo` and `HiMT` might easily adopt `syncasm` for their main assembler.

The k-mer based method decrease the nucleotide complexity with the extent of the k size.
For selecting reads using reference sequences, we could still used alignment based method.
The alignment-based method could be slower but it could deliver organelle genome assemblies.

`Polap` and `HiMT` use amino acid sequences from plant organelle genes in their pipeline.
`GetOrganelle` and `ptGAUL` use nucleotide sequences in their pipelines.
`TIPPo` directly applies seed selecting scheme without reference-generating.

We aim to take these approaches together to assemble plastid genomes to show that one could select long-reads for plastid genome assembly.
We use CDS sequences of the amino acid sequences that we used previously in `Polap` to select long reads.
Instead of BLAST as pairwise alignment tool used in `Polap` and `HiMT`, we use `minimap2` to select long reads.
There are many different ways of constructing a set of DNA sequences to select reads for plant organelle genome.
Because we have used the dataset in `Polap` for selecting seed contigs, it is one case that is applied for selecting long-read data for plant organelle genome.
For mitochondrial genome assemblies, this method may not be applicable or it may require a specific set of reference DNA sequences.
Here, we showcase a read-selection approach for selecting low-quality ONT reads for plastid genome assembly.

- `BWA` [@Li2009]
- `Minimap2` [@Li2018]
- `FMLRC` [@Wang2018]
- `GeSeq` [@Tillich2017]
- Subsampling [@Efron1987]
- `NextDenovo` [@Hu2024]

# Results

## Comparison with other plastid assembly pipelines

<!-- polap man table-benchmark-summary some 0 polap -->

We used 38 ONT and 98 PacBio HiFi long-read datasets.

For PacBio HiFi datasets, all of them have similar success rates of plastid genome assembly.
`Oatk` assembled plastid genomes is 78 PacBio HiFi datasets.
It misassembled in 11 PacBio HiFi datasets.
`TIPPo` assembled plastid genomes in 86 PacBio HiFi datasets.
`ptGAUL` assembled plastid genomes in 40 PacBio HiFi datasets.
It misassembled in 37 PacBio HiFi datasets.
We suspect that this large number of mis-assembly should be due to the misspecification of `Flye` data type option because `ptGAUL` does not have an option for `Flye` data type.
`Polap` assembled plastid genomes in 91 PacBio HiFi datasets.

`Oatk` did not assemble plastid genomes using the ONT datasets.
`TIPPo` assembled plastid genomes in 6 ONT datasets.
`ptGAUL` assembled plastid genomes in 33 ONT datasets.
It misassembled in one ONT dataset.
`Polap` assembled plastid genomes in 34 ONT datasets.

(Supporting Materials – Plastid genome assemblies using the six pipelines)

(@fig:time-some-0)

(@fig:memory-some-0)

(@tbl:benchmark-data-some-0)

(@tbl:benchmark-computer-some-0)

Code

Assembly graph figure

# Discussion

`Polap` (Plant Organelle Long-read Assembly Pipeline v0.4.3.7), which includes the subsampling-based plastid genome assembly feature, is available under the GNU General Public License version 3.0 at [http://github.com/goshng/polap](http://github.com/goshng/polap).

# Materials and Methods

Download all reference mito genome genbank data for plant mito genome.
Remove coding genes and collect the remaining parts for a set of noncoding sequences.

We collected a set of CDS sequences for the plant organelle genes, for which we used amino acid sequences to select seed contigs using BLAST.

Because plastid-derived reads are abundant, we downsample the ONT reads to 3 gigabases.
For HiFi data, we estimate the genome size and downsample to 3x coverage.

We map the reads on the CDS sequences using `minimap2` and count alignment segments for the reads to rank by the number of plastid genes mapped on the reads.
We collect reads upto 30 Mb of cumulative length distribution and sample 10 Mb from the collection to form an input data for Flye v2.9.6 assembly.
From the assembly, we choose plastid-like contigs as seed for mapping reads to select plastid-derived reads. We then execute the next round of assembly.
We repeat the 2nd round of assembly several times until we find plastid-like genome assembly graph.

<!-- (@tbl:benchmark-data-some-2) -->

<!-- (@fig:mainfigure1) -->

- `Flye` [@Kolmogorov2019]
- `JellyFish` [@Marcais2011]
- `BLAST` [@Altschul1997]
- `SeqKit` [@Shen2016]
- `MAFFT` [@Katoh2013]

<!-- (@fig:mainfigure1) -->

- `Bandage` [@Wick2015]

- `Canu` [@Koren2017]
- `NextDenovo` better than `Canu` [@Wick2021]

<!-- (@tbl:benchmark-polap-some-0) -->

<!-- (@tbl:benchmark-computer) -->

# Supplementary Material

Supplementary material, including 10 tables and three figures, is appended to the main text of this manuscript. A `BASH` script for executing the pipeline used to generate the results presented in the manuscript is also included.

# Acknowledgements {-}

We thank Jeffrey L. Thorne for improving the presentation of this work.

# Author Contributions

S.C.C. developed the Polap pipeline and prepared the manuscript.

# Conflict of Interest

The author declare no conflicts.

# Data availability

`Polap` (Plant Organelle Long-read Assembly Pipeline v0.4.3.7) is available under the GNU General Public License version 3.0 at [http://github.com/goshng/polap](http://github.com/goshng/polap).
The results presented in this manuscript are available at Figshare: [https://figshare.com/s/ec1cb394870c7727a2d4](https://figshare.com/s/ec1cb394870c7727a2d4).

# References {-}

::: {#refs}
:::

<!-- # Species -->
<!---->
<!-- - _Anthoceros agrestis_ -->
<!-- - _Arabidopsis thaliana_ -->
<!-- - _Canavalia ensiformis_ -->
<!-- - _Cinchona pubescens_ -->
<!-- - _Codonopsis lanceolata_ -->
<!-- - _Cucumis sativus_var_hardwickii_ -->
<!-- - _Dioscorea japonica_ -->
<!-- - _Dunaliella tertiolecta_ -->
<!-- - _Eucalyptus pauciflora_ -->
<!-- - _Euonymus alatus_ -->
<!-- - _Gossypium herbaceum_ -->
<!-- - _Juncus effusus_ -->
<!-- - _Juncus inflexus_ -->
<!-- - _Juncus roemerianus_ -->
<!-- - _Juncus validus_ -->
<!-- - _Leiosporoceros dussii_ -->
<!-- - _Macadamia jansenii_ -->
<!-- - _Musa acuminata_subsp_malaccensis_ -->
<!-- - _Ophrys lutea_ -->
<!-- - _Oryza rufipogon_ -->
<!-- - _Pisum sativum_ -->
<!-- - _Populus x_sibirica_ -->
<!-- - _Prunus mandshurica_ -->
<!-- - _Pterocarpus santalinus_ -->
<!-- - _Solanum lycopersicum_ -->
<!-- - _Spirodela polyrhiza_ -->
<!-- - _Vaccinium vitis-idaea_ -->
<!-- - _Vitis vinifera_ -->
<!-- - _Notothylas orbicularis_ -->
<!-- - _Phaeomegaceros chiloensis_ -->
