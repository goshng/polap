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

# Title: On the plastid genome assembly by subsampling low-quality Oxford Nanopore long-read data

Authors

Affiliations

_Corresponding author_:
Author and mailing address

Email: author@email.edu

Tel.

Fax.

_Running head_:
Plastid genome assembly via long-read subsampling

_Keywords_:
ptGAUL, Flye, Minimap2, TIPPo, PMAT, Oatk, long-read assembly

Word count of the abstract: xxx

Word count of the text: xxxx

\newpage

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
- `Bandage` [@Wick2015]
- `Canu` [@Koren2017]
- `NextDenovo` better than `Canu` [@Wick2021]

(@tbl:benchmark-polap-some-0)

(@tbl:benchmark-computer)

# Results

(@tbl:benchmark-polap-some-1)

(Supporting Materials – Plastid genome assemblies using the six pipelines)

<!--(@tbl:benchmark-time-2)-->

<!--(@tbl:benchmark-memory-2)-->

(@tbl:benchmark-polap-some-1)

(@tbl:benchmark-polap-some-0)

(@tbl:benchmark-polap-some-4)

- (@tbl:polap-disassemble-Eucalyptus-pauciflora-2)

- (@tbl:polap-disassemble-Eucalyptus-pauciflora-0)

- (@tbl:polap-disassemble-Eucalyptus-pauciflora-3)

# Discussion

content

# Supplementary Materials

content

# Acknowledgements {-}

content

# Author Contributions

content

# Data availability

content

# References {-}

::: {#refs}
:::
