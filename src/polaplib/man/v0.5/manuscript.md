---
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

Author Last and first names

Address

\newpage

# Abstract {-}

Content of abstract.

\newpage

# Introduction {-}

- `BWA` [@Li2009]
- `Minimap2` [@Li2018]
- `SPAdes` [@Bankevich2012]
- `Flye` [@Kolmogorov2019]
- `FMLRC` [@Wang2018]
- `GeSeq` [@Tillich2017]
- `GetOrganelle` [@Jin2020]
- `ptGAUL` [@Zhou2023]
- `CLAW` [@Phillips2024]
- `PMAT` [@Bi2024]
- `Oatk` [@Zhou2024]
- `TIPPo` [@Xian2025]
- Subsampling [@Efron1987]
- `NextDenovo` [@Hu2024]

# Materials and Methods

(@tbl:benchmark-data-some-2)

(@fig:mainfigure1)

- `Flye` [@Kolmogorov2019]
- `JellyFish` [@Marcais2011]
- `BLAST` [@Altschul1997]
- `SeqKit` [@Shen2016]
- `MAFFT` [@Katoh2013]

(@fig:mainfigure1)

- `Bandage` [@Wick2015]

- `Canu` [@Koren2017]
- `NextDenovo` better than `Canu` [@Wick2021]

(@tbl:benchmark-polap-some-0)

(@tbl:benchmark-computer)

# Results

## Comparison with other plastid assembly pipelines

(@tbl:benchmark-polap-some-1)

(Supporting Materials – Plastid genome assemblies using the six pipelines)

(@tbl:benchmark-time-2)

(@tbl:benchmark-memory-2)

## Subsampling-based plastome assemblies

(@tbl:benchmark-polap-some-1)

(@tbl:benchmark-polap-some-0)

(@tbl:benchmark-polap-some-4)

## Three-stage of subsampling-based assembly

- (@tbl:polap-disassemble-Eucalyptus_pauciflora-2)

- (@tbl:polap-disassemble-Eucalyptus_pauciflora-0)

- (@tbl:polap-disassemble-Eucalyptus_pauciflora-3)

# Discussion

`Polap` (Plant Organelle Long-read Assembly Pipeline v0.4.3.7), which includes the subsampling-based plastid genome assembly feature, is available under the GNU General Public License version 3.0 at [http://github.com/goshng/polap](http://github.com/goshng/polap).

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
