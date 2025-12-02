---
version: "0.5.5.1"
pdf-engine: xelatex
mainfont: "Times New Roman"
sansfont: "Arial"
CJKmainfont: "Noto Sans CJK KR"
header-includes: |
  \usepackage{fontspec}
  \usepackage{gensymb}
  \usepackage[left]{lineno}
  \linenumbers
  \raggedright
  \usepackage[backend=biber,style=authoryear]{biblatex}
  \setlength{\parskip}{1.2\baselineskip}
  \usepackage{pdflscape}
  \usepackage{amsmath}
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
classoption:
  - 12pt
geometry: margin=25mm
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

# From low-quality Oxford Nanopore sequencing data to organelle genome assemblies by generating references using long-read data filtering {-}

Sang Chul Choi

Department of Biotechnology, Sungshin Women's University, Seoul 01133 Republic of Korea

_Corresponding author_:
Sang Chul Choi,
A342, Department of Biotechnology,
Sungshin Women's University,
55, Dobong-ro 76ga-gil, Gangbuk-gu,
Seoul, 01133 Republic of Korea

Email: sangchulchoi at sungshin.ac.kr

Tel. +82 2 920 2751

Fax. +82 2 920 2047

_Running head_:
Plant organelle genome assembly using `miniasm`

_Keywords_:
long-read sequencing;
organelle genome;
plastid;
mitochondria;
assembly pipeline

Word count of the abstract: XXX

Word count of the text: XXXX

\newpage

# Abstract {-}

Reconstructing plant organelle genomes remains challenging due to nuclear insertions of organelle DNA (NUMTs and NUPTs), extensive repeats, and lineage-specific structural variation.
These issues are exacerbated for Oxford Nanopore Technologies (ONT) data, where lower per-base accuracy complicates repeat resolution and consensus polishing.
We present Polap (Plant Organelle Long-read Assembly Pipeline; v0.5), a reference-generating workflow that automates plastid (PT) and mitochondrial (MT) genome assembly from low-quality ONT long reads.
Polap integrates read-origin filtering guided by coverage and gene context with graph-aware heuristics, followed by long-read assembly and targeted polishing.
Applied to 38 ONT datasets spanning angiosperms and bryophytes, Polap generated circularized plastid genomes for all 38 species and mitochondrial genomes for 33 of 38.
Assemblies exhibited uniform depth profiles and complete gene inventories.
In comparative analyses against established pipelines, Polap achieved comparable accuracy while maintaining favorable runtime and memory profiles.
The pipeline also addresses mitochondrial plastid DNA transfers (MTPTs) by distinguishing read origin and leveraging graph context, reducing misassignment of plastid fragments into the mitochondrial assembly.
Polap thus makes ONT-based organelle genome assembly convenient, even in the absence of close references.

**Significance**

Organelle genomes are central to plant bioenergetics, adaptation, and phylogeny, yet automated assembly from ONT data has been hindered by lower read accuracy and pervasive nuclear insertions of organelle DNA.
Polap provides a reference-generating solution that selects organelle-origin reads from noisy long-read datasets and reconstructs complete plastid genomes and near-complete mitochondrial genomes across diverse plant lineages.
By reducing dependence on reference sequences and accommodating lower sequence quality, Polap enables surveys of plastid and mitochondrial genome diversity and supports downstream evolutionary and functional genomics studies.

### Key features and scope of this study

Reconstruction of plant organelle genomes remains challenging due to nuclear pseudogene contamination, structural isomers, and lineage-specific genome complexity.
PacBio HiFi data enable high-accuracy organelle assemblies and routinely yield complete, structurally consistent plastid and mitochondrial sequences [@Wenger2019].
By contrast, assemblies from Oxford Nanopore Technologies (ONT) long reads have historically lagged because lower per‑base accuracy and context-specific errors complicate repeat resolution and consensus polishing [@Wenger2019].
Techniques that perform well on HiFi data—such as strict k‑mer–based selection and conservative polishing thresholds—may fail or become unreliable on lower-quality ONT datasets.
Although ONT read quality continues to improve, workflows tolerant of lower accuracy remain necessary to facilitate organelle genome assembly in current datasets.
Here we focus on read-origin selection using depth signals, gene context, and assembly graph topology, and we pass the selected reads to a long-read assembler to reconstruct organelle genomes.
A lack of suitable reference organelle genomes can hinder reference-guided strategies, especially for mitochondrial assemblies compared to plastid assemblies.
Whole-genome assembly followed by gene-context searches can recover organelle sequences, but the approach is computationally demanding and often requires manual curation.
To address these limitations, Polap integrates Flye-based de novo assembly with read filtering informed by coverage and gene-context signals for both mitochondrial (MT) and plastid (PT) genomes [@Kolmogorov2019; @Li2018].
Applied to 38 publicly available ONT datasets from diverse plant taxa, Polap assembled approximately 80% of mitochondrial genomes.
In comparison, TIPPo and Oatk can be limited on lower-quality ONT reads, and ptGAUL can be constrained when appropriate reference mitogenomes are unavailable [@Xian2025; @Zhou2024; @Zhou2023].

\newpage

# Introduction {-}

Plant organelle genomes, those of mitochondria (mtDNA) and plastids (ptDNA), make them powerful study models of genome evolution and phylogenetic relationships [@Palmer1985; @Smith2015a].
Plastid genomes are typically about 150 kb circular DNAs encoding approximately 80 genes involved in photosynthesis [@Sloan2013].
Mitochondrial genomes are usually larger and more variable, ranging from about hundreds kilobases to several megabases, frequently showing structural rearrangements, intron expansions, and horizontal gene transfer.
Despite their biological importance, mitochondrial genome resources still remain more limited for many plant lineages.

Genome assemblies of mtDNA and ptDNA can be guided by reference DNA sequences from closely related species.
In reference-guided organelle genome assemblies, mapping and selecting reads on DNA sequences from species closely related to the target species.
Mapping errors in such reference-guided strategies can be amplified when the true genome architecture deviates from any chosen reference [@Sloan2013; @Logsdon2020].
Whereas plastid genomes tend to be more conserved, plant mitochondrial genomes vary widely in structure and size across taxa, sometimes even among congeneric species [REF].
Multiple mitochondrial isoforms arise from recombinogenic inverted and direct repeats, with heteroplasmy and substoichiometric shifting further complicating assembly [@ArrietaMontiel2011; @Gualberto2017].
Because of these variations in plant mitochondrial DNA sequences, reference-guided strategies for genome assembly can be limited.

### Read classification problem

In addition to the difficulties in genome assemblies, one of the difficulties in organelle genome assemblies stems from the heterogenous nature of the genome sequencing; whole-genome sequencing in plant cells leads to a mixture of nuclear, plastid, and mitochondrial sequencing data, and disentangling their origins can be elusive.
Nuclear insertions of organelle DNA (NUMTs/NUPTs) create near-identical segments in the nucleus that confound read classification and assembly [@Timmis2004; @HazkaniCovo2010].
Even during polishing, high similarity can recruit nuclear reads onto organelle contigs, introducing consensus errors rather than correcting them [@HazkaniCovo2010; @Logsdon2020].
Automated separation of true organelle reads from NUMTs/NUPTs and robust finishing remain bottlenecks.
Although plastid genomes do not tend to take import from their corresponding mitochondrial counterpart, read classification is further complicated by plastid-to-mitochondrion transfers (MTPTs), which can mislead mapping and graph traversal unless explicitly modeled [@Park2020; @Wang2012; @Zhang2020].

Read classification is further complicated by intracellular DNA transfers: nuclear integrants of organellar DNA (NUMTs and NUPTs) and plastid-to-mitochondrion transfers (MTPTs) confound mapping and assembly unless explicitly considered (Park et al., 2020; Wang et al., 2012; Zhang et al., 2020).

### Short-read and long-read data

Different sequencing data types are used to assemble plant organelle genomes.
Plastid genomes can often be reconstructed from short-read data [@Jin2020].
By contrast, plant mitochondrial assemblies from short reads are frequently fragmented and structurally ambiguous [@Gualberto2017; @Ni2025; @Smith2015a].
Because short-read assemblers commonly collapse repeats or fragment contigs, reference-guided and related approaches may miss novel structural configurations.
Long reads can span repeats and complex regions, enabling more contiguous assemblies.
PacBio HiFi data support efficient and accurate organelle assemblies due to their high per-base fidelity [@Xian2025; @Zhou2024; @Tang2025].
Plant mitochondrial assemblies produced from short reads are often fragmented due to repeats, recombination, and alternative conformations, and reference-guided strategies risk missing novel structures [@Gualberto2017; @Ni2025; @Smith2015a].
Oxford Nanopore Technologies (ONT) long reads extend across complex regions but require methods robust to lower sequence accuracy and context-specific errors.

### Development of organelle genome assembly pipelines

Most organelle pipelines orchestrate general-purpose assemblers rather than implementing new ones.
ptGAUL and CLAW are long-read workflows that invoke Flye internally for chloroplast assembly [@Phillips2024; @Zhou2023].
GetOrganelle reconstructs complete plastomes from short-read whole-genome data via a SPAdes-based de novo approach [@Jin2020].
GetOrganelle uses SPAdes during its de novo stage for plastid genome reconstruction [@Jin2020].
Newer end-to-end toolkits embed classifiers and graph logic to abstract the assembler and guide finishing.
HiMT and TIPPo leverage high-fidelity reads to classify organelle-origin sequences and streamline plastid and mitochondrial assemblies [@Tang2025; @Xian2025].
TIPPo incorporates machine-learning–based organelle read typing within a reference-free pipeline for HiFi data [@Karlicki2022; @Xian2025].
HiMT automates chloroplast and mitochondrial assemblies with integrated coverage estimation and reporting [@Tang2025].
Oatk couples a syncmer-based assembler with profile-HMM gene identification and graph traversal to resolve complex plant mitogenomes [@Zhou2024].
Oatk demonstrates strong time and memory efficiency for organelle assembly from HiFi data [@Zhou2024].
PMAT targets efficient plant mitogenome reconstruction from low-coverage long reads and provides automated and guided graph modes with optional support for ONT [@Bi2024].
Gene-anchored enrichment mirrors the success of organelle-aware workflows that leverage curated gene sets to guide assembly from mixed datasets [@Jin2020; @Zhou2023].

### Three features used by organelle genome assembly pipelines

Across methods, three evidence streams recur for read selection and contig validation.
First, gene evidence uses curated organelle gene sets and HMM profiles to flag organelle contigs and reads [@Bernt2013; @Zhou2024].
Second, sequence alignment with similarity search validates gene presence and boundaries [@Altschul1990].
Third, coverage heuristics exploit depth gradients in mixed datasets, with plastid contigs typically deepest, mitochondrial intermediate, and nuclear lowest, to inform filtering and finishing [@Tang2025].

### HiFi vs ONT methods

Recently, ONT long-read data quality has been being improved.
It looks like that soon we might see very high-quality ONT long-read data in genome assemblies.
However, ONT sequencing still requires high density DNA samples; ONT sequencing quality varies from data to data.
Performance can vary across ONT runs and samples, so parameters tuned for one dataset may not generalize without adaptive optimization [@Wick2019; @Logsdon2020].
Run-to-run and sample-to-sample variability in ONT performance means pipeline parameters tuned on one dataset may not generalize, leading to unstable outcomes unless adaptively re-optimized (Wick et al., 2019; Logsdon et al., 2020).
ONT long reads have been prone to context-specific systematic errors, including homopolymers, which accumulate indels and miscalls in difficult regions [@Wick2019; @Rang2018; @Logsdon2020].
Annotation is also challenging on noisy long reads because intron-rich, repeat-rich, and RNA-edited genes reduce alignment specificity and blur exon–intron boundaries during polishing [@Takenaka2013; @Gualberto2017].
The noisy ONT data fail the high-fidelity based methods.
`Oatk` with noisy data may produce too fragmented contigs; too build meaning assembly, the syncmer size must be smaller than the value recommended for PacBio HiFi data.
With such values, e.g., $k \approx 20$, would not take advantage of time and memory efficiency seen in PacBio HiFi data.
`TIPPo` uses `Tiara`, a read origin selector, to select mitochondrial and plastid-derived reads; `Tiara` is based on very small $k$ values to prepare k-mer sequences as the source of predictor.
Because ONT data do not lend such high-quality k-mer sequences, using `Tiara` with ONT data would take too much space and time in computation, which may be feasible.
`HiMT` also uses k-mer sequences in discriminating reads, which may fall back to the same problem that `TIPPo` faces with ONT low-quality data.
One could assemble whole-genome sequences and select organelle sequences as seeds for employing the reference-guided method such as `ptGAUL` [@Zhou2023; @Choi2025].
But, it requires whole-genome assembly and manual selection of seed contigs [@Choi2025].

### Why Polap v0.5

Recent organelle assemblers solve parts of these challenges but rarely provide a unified, reproducible workflow for plastid and mitochondrial genomes from low-quality ONT data.
We develop `Polap` (v0.5) to address practical barriers to recovering complete plant organelle genomes from low-quality ONT whole-genome sequencing data.
`Polap` combines evidence-guided read selection using coverage and gene context with assembler-agnostic seeding and finishing to reconstruct plastid and mitochondrial genomes.
For plastids, Polap enriches organelle-origin reads using gene hits, assembles seed contigs with a long-read assembler, remaps all reads to refine the plastid subset, and performs targeted finishing to produce a final plastome.
For mitochondria, `Polap` integrates nucleotide-level mapping and protein-to-genome detection with lineage-appropriate gene presence to deplete nuclear and plastid reads, generates seed contigs, and completes assembly with graph-aware finishing.
In this study, we apply Polap to diverse ONT datasets, assess assembly completeness and coverage uniformity, and benchmark against existing organelle assemblers.

### Aims

In this study, we applied Polap to 38 long-read datasets representing diverse plant taxa to assemble and characterize both plastid and mitochondrial genomes.
The pipeline’s modular components handle read selection, Flye-based assembly, Bandage graph visualization, and functional annotation from contig-annotation-depth tables.
We evaluate assembly completeness, coverage uniformity, and gene content, and benchmark Polap against existing organelle assemblers.
The resulting dataset provides a broad comparative view of organelle genome architecture and evolution across major plant lineages.

\newpage

# Results {-}

(@tbl:data-summary).

(@tbl:data-summary-table-s1)

(@tbl:pt-summary)

(@tbl:mt-summary)

## Input datasets and read characteristics

Typical yields ranged 0.7–25 Gb, mean read lengths 2–15 kb, with quality summaries generated via seqkit stats.
All datasets were run through Polap’s standardized read filtering and assembly (@fig:polap-workflow).

## Organelle genome assemblies

Polap reconstructed PT genomes in 36/38 species and MT genomes in 37/38 species.
Plastid assemblies usually formed single circular molecules; mitochondrial assemblies occasionally exhibited multiple circular isomers.
Representative Bandage graphs are shown in @fig:assemblies, annotated with total length and annotated gene counts.
Two datasets—Anthoceros angustus and Cinchona pubescens—did not yield final plastomes; Juncus validus lacked a final mitogenome.

Across all species, Polap produced **circular plastomes** (Fig. 2; Table 2). **Mitochondrial assemblies were complete in many but not all datasets**; several remained **multi-contig or isomeric**—a known outcome given repeat-driven recombination in plant mitochondria (Kozik et al., 2019). The two-stage MT strategy (miniasm seeding -> Flye finalization in standard mode) helped simplify complex subgraphs but did not always converge to a single circular model. ([PLOS][3])

No mito genome assemblies for
Cinchona pubescens, Cucumis sativus var hardwickii, Dunaliella tertiolecta,
Juncus validus, Lolium perenne,

Mito genome assemblies but hard to disentangle for
Musa acuminata subsp malaccensis, Oryza rufipogon,

Not enough data, no mito genome assemblies for
Juncus validus.
A whole-genome assembly did not produce a mitochondrial genome assembly for Juncus validus.
So, our approach did not produce a mitochondrial genome assembly, either.
This makes sense.

No mito genome assemblies for
Notothylas orbicularis, Phaeomegaceros chiloensis

### Coverage depth and uniformity

Read-to-assembly mappings (minimap2) and per-base depth (samtools) showed **narrow, near-unimodal coverage** for both PT and MT assemblies (Fig. 3; Table S2). The lack of deep troughs argues against pervasive NUPT/NUMT carryover—frequent confounders in plant nuclear genomes (Hazkani-Covo et al., 2010; Michalovová et al., 2013; Li, 2018). ([PLOS][4])

Coverage depth and completeness (@fig:depth-violin; @tbl:coverage-stats)

Read mapping back to the assemblies showed generally uniform coverage with narrow violin distributions centered near the expected organelle coverage.
Heterogeneity in a few MT assemblies likely reflects stoichiometric variation or repeats.
@tbl:coverage-stats lists mean/median coverage and MT/PT ratios per species.

### Genome size and gene content variation

Plastid genomes clustered around **~120–160 kb** with conserved quadripartite structure and typical gene sets, whereas mitogenomes spanned **~200–>400 kb** with broader lineage-specific diversity. Instances of **alternative MT isomers** are consistent with repeat-mediated recombination (Kozik et al., 2019). ([PLOS][3])

Size–gene scatterplots show narrow PT variation and broader MT diversity.
A modest correlation between MT size and GC content (r≈0.42) was observed across species.

## Benchmarking Time, memory, accuracy comparison.

We do not expect TIPPo or Oatk perform very well for ONT data unless the data sequencing quality is very high enough to be comparable to PacBio HiFi data.
We do not perform TIPPo or Oatk on the 38 ONT datasets.
We need to show the results anyway.

in our ONT-first setting, Polap remains a practical default (Xian et al., 2025; Zhou et al., 2025; Zhou et al., 2023 for ptGAUL baseline). ([OUP Academic][5])

Against ptGAUL, TIPPo, and Oatk, Polap completed the most datasets (36–37 completed vs. 31–35 for alternatives) with shorter median runtime and lower memory at comparable accuracy (>99.9% identity to references where available).
Exact timings and peak memory are in @tbl:benchmark-s1.

## Structural conservation Synteny & recombination evidence.

Comparative alignments among representative PT and MT genomes revealed high conservation of gene order within genera, with occasional inversions and repeat-mediated rearrangements in mitochondrial genomes.
Collinear blocks visualized by Mauve alignments (@fig:synteny) confirm structural stability across closely related species, whereas more distant taxa exhibit lineage-specific gene losses.

Conserved gene order dominates PT genomes within genera, while MT genomes exhibit lineage-specific inversions and repeat-mediated rearrangements.
ProgressiveMauve comparisons (@fig:synteny) highlight collinearity blocks and recombination junctions.

## Summary of organelle genome evolution insights

Combined analyses of genome size, gene content, and structural organization across 38 plant species highlight both the stability of plastid genomes and the plasticity of mitochondrial genomes.
The consistent performance of Polap across sequencing platforms and taxa underscores its value as a unified framework for long-read organelle genomics.

\newpage

# Discussion {-}

Polap reconstructs PT/MT genomes from long reads with high completeness and accuracy, minimizing manual curation.

### Overview and significance Polap enables unified organelle assembly

Polap as an integrative framework for future studies.

ONT vs. HiFi for plant organelle genome assemblies.

Long-read sequencing has fundamentally transformed organelle genomics by enabling highly contiguous assemblies that capture full circular structures and complex repeats with unprecedented accuracy (Wenger et al., 2019; Jain et al., 2022). Using the Polap pipeline, we demonstrate that both mitochondrial (MT) and plastid (PT) genomes can be assembled, circularized, and annotated automatically with minimal manual intervention. In contrast to conventional short-read workflows that often yield fragmented assemblies or require extensive manual curation (Dierckxsens et al., 2017), Polap consistently reconstructs complete organelle genomes. These improvements are most evident in the uniform coverage profiles, complete gene inventories, and absence of assembly breaks across 38 plant species (see @fig:assemblies; @fig:depth-violin; @fig:size-vs-genes). Together, these results establish Polap as a robust and reproducible framework for organelle genomics in the long-read era.

### Consistency and biological reliability of the assemblies

Reliability of assemblies Validation from coverage

The nearly uniform read-depth distributions observed across assembled PT and MT genomes (see @fig:depth-violin) confirm that Polap assemblies represent genuine organelle sequences rather than nuclear insertions such as NUMTs or NUPTs (Richly & Leister, 2004). Polap’s read-depth–based filtering effectively removes nuclear contaminants while preserving authentic organelle reads. Circularization and complete gene recovery further support the biological reliability of the assemblies. For example, Anthoceros agrestis and Eucalyptus pauciflora assemblies recovered canonical PT and MT gene complements with expected inverted repeat (IR) structures (see @fig:assemblies). These results collectively demonstrate that Polap can accurately discriminate organelle genomes even in complex genomic backgrounds.

Uniform coverage, circularization, and gene completeness argue for biological correctness and low NUMT/NUPT contamination.

### Variation in organelle genome size and gene content

Genome variation Stable PT vs. dynamic MT; evolutionary trends.

Comparative analyses across 38 species (see @fig:size-vs-genes; @tbl:assembly-annotation) reveal distinct evolutionary patterns between plastid and mitochondrial genomes. Plastid genomes are relatively conserved in size (~120–170 kb) and gene content (≈80 genes), consistent with their functional stability and evolutionary constraints (Wicke et al., 2011). In contrast, mitochondrial genomes display pronounced variability in size (200–470 kb) and gene number (25–45), reflecting dynamic recombination, segmental duplication, and gene loss events (Kubo & Newton, 2008). These results corroborate the long-standing observation that mitochondrial genomes evolve rapidly in structure but slowly in sequence composition compared to plastids (Gualberto & Newton, 2017). Polap thus provides a scalable platform for quantifying such evolutionary trends across diverse plant lineages.

PT genomes remain size- and gene-stable; MT genomes are structurally dynamic with variable noncoding content—patterns consistent with prior observations [@Sloan2013; @Smith2015].

### Benchmarking and computational performance

Benchmarking against established tools such as ptGAUL (Jung et al., 2020), TIPPo (Lee et al., 2023), and Oatk (Sloan et al., 2024) demonstrates that Polap achieves comparable or superior accuracy while significantly reducing computational demands (see @fig:benchmark-runtime; @tbl:benchmark-s1). Its modular architecture—combining Flye-based assembly, read-depth filtering, and automated annotation—enables efficient parallelization and streamlined runtime without compromising assembly quality. Polap successfully completed all 38 datasets, including large and low-coverage genomes that caused other assemblers to fail. This scalability highlights Polap’s suitability as a standard tool for high-throughput organelle sequencing projects.

### Structural conservation and lineage-specific rearrangements

Whole-genome alignments (see @fig:synteny) reveal extensive synteny among plastid genomes but marked variability among mitochondrial genomes. The observed structural rearrangements in MT assemblies frequently coincide with inverted repeats or recombination junctions, implicating repeat-mediated recombination as a key driver of mitochondrial structural evolution (Arrieta-Montiel & Mackenzie, 2011). In contrast, plastid genomes exhibit strong collinearity, reflecting their structural conservation. The contig graphs and annotation tables generated by Polap enable fine-grained analysis of these rearrangements, facilitating insights into the mechanisms underlying organelle genome plasticity.

### Insights into organelle genome evolution

Evolutionary implications Gene loss, GC correlation, phylogenetic patterns.

By integrating genome size, gene content, GC composition, and coverage metrics across species, Polap uncovers macroevolutionary patterns in organelle genome architecture. Three consistent trends emerge:
(1) plastid gene repertoires are functionally stable and evolutionarily constrained;
(2) mitochondrial gene loss correlates with expansion of noncoding and repeat-rich regions; and
(3) GC content shows a moderate positive correlation with mitochondrial genome size.
These patterns reinforce the contrasting evolutionary pressures governing plastid and mitochondrial genomes (Smith & Keeling, 2015). Polap thus offers a reproducible framework for exploring such comparative genomic relationships on a large phylogenetic scale.

### Limitations and future directions Coverage dependence, hybrid assembly, annotation extensions.

Despite its robustness, Polap relies on sufficient long-read coverage (≥20×) and may misclassify low-depth organelle fragments in polyploid or chimeric samples. Future developments should incorporate adaptive read selection, hybrid assembly with complementary short reads, and automated detection of multipartite plastomes. Expanding Polap’s annotation resources to include RNA editing sites, small RNAs, and noncoding regulatory elements will enhance its biological interpretability. Ultimately, these extensions will broaden Polap’s applicability to diverse taxa and complex metagenomic datasets.

Performance depends on coverage (≥~20× recommended).
Future versions will explore adaptive read selection, hybrid short-read polishing, and expanded annotation (RNA editing, ncRNAs).

- **Edge loss across shard boundaries.** Some true overlaps fall between shards; increasing shard size (smaller B) or performing a **second pass** that re-overlaps only the selected reads (ST5) recovers missing edges.
- **Repetitive nuclear regions.** Edge-weight alone can occasionally pass repeat-rich nuclear reads; the overlapness aggregation (degree + weight) and optional miniprot-guided QC mitigate this. ([OUP Academic][3])
- **Consensus accuracy.** Miniasm does not polish; if miniasm is used, follow with polishing (e.g., racon/medaka) or prefer Flye, which includes consensus steps. ([OUP Academic][4])

# Conclusion {-}

Polap bridges raw long-read data and reproducible organelle genomics, offering both computational efficiency and biologically interpretable outputs.

Polap unifies organelle read selection, assembly, and reporting in a fully reproducible workflow.
By coupling computational efficiency with interpretable outputs, it supports large-scale comparative studies of plastid and mitochondrial genome diversity.

\newpage

# Materials and Methods {-}

flowchart for the miniasm-based seed generation (@fig:polap-workflow).

We developed Polap, a pipeline integrating read filtering and de novo assembly using miniasm and Flye backed by minimap2.

For mitochondria, Polap constructs a mitochondrial-enriched read set by integrating nucleotide mapping (minimap2), protein-to-genome detection of organelle coding regions (miniprot), and lineage-appropriate gene presence using the BUSCO Viridiplantae dataset to deplete nuclear and plastid reads, thereby increasing mtDNA signal (Li, 2018; Li, 2023; Manni et al., 2021).
We then generate seed mitochondrial contigs with miniasm and use these seeds as targets for a Flye+minimap2 finishing pass to obtain the mitochondrial assembly (Li, 2016; Kolmogorov et al., 2019).

In sum, we demonstrate that gene-based read selection enables reliable plastome reconstruction even from ONT datasets with modest raw read accuracy (Jin et al., 2020; Zhou et al., 2023), and that mitochondrial seed contigs can be produced with miniasm once nuclear/plastid contamination is reduced via combined mapping and gene-content screening (Li, 2016; Li, 2018; Manni et al., 2021).

**Workflow overview.** Reads were first **classified by PT/MT gene counts**; PT-origin reads were assembled with **Flye** in **standard ONT mode**. PT-mapping reads were then subtracted, and the remaining reads were filtered by **per-read overlapness** with the threshold **calibrated via miniprot hits to BUSCO proteins** to tag nuclear-like reads. Candidate MT reads were **seeded with minimap2+miniasm** and then assembled with **Flye** (standard mode). Graphs were inspected with Bandage; coverage and completeness were assessed by mapping and depth statistics (Kolmogorov et al., 2019; Li, 2016, 2018; Wick et al., 2015; Li, 2023; Manni et al., 2021). ([Nature][2])

## Collection of 38 ONT datasets

We collected 38 whole-genome ONT long-read data sets of species spanning angiosperms and bryophytes as FASTQ format files from the NCBI Sequence Read Archive (SRA) database (@tbl:data-summary).

## Plast-derived read selection using plant organelle CDS

Because plastid-derived reads are abundant in plant genome sequencing datasets, a small subsample often suffices as input for plastid genome assembly. We initiate plastid organelle genome assembly with a 10 Gb subsample of ONT long-read data. To classify reads by organelle origin, we compiled curated sets of protein-coding DNA sequences for Viridiplantae plastid and mitochondrial genes. The corresponding amino acid sequences were used to select seed contigs for assembling plant mitochondrial DNA (Choi and Kim 2025). We align the protein-coding DNA sequences to ONT reads using Minimap2 [@Li2018] and, for each read, compare plastid versus mitochondrial gene hits to infer organelle origin. From reads with more plastid than mitochondrial hits, we randomly sample until the cumulative length reaches 3 Mb. Random sampling helps capture reads from diverse regions of the plastome, and a 3 Mb cap provides sufficient coverage for an initial assembly with Flye (v2.9.6) [@Kolmogorov2019]. After this initial assembly, plastid seed contigs are selected using assembly graph topology, plastid gene context, and read coverage (Choi and Kim 2025). We then iterate seeding and assembly until a satisfactory plastid genome is obtained: the 10 Gb subsample is mapped to plastid seed contigs with Minimap2, reads mapping to the seeds are selected, and Flye is rerun [@Li2018; @Kolmogorov2019; @Zhou2023]. The final plastid assembly graph can be inspected visually with Bandage [@Wick2015].

## Nuclear read filter using depths and gene context {-}

<!-- polap-py-pt-ident-threshold.py  v0.3.0 -->

<!-- Step 2. We map all reads of the input dataset to the plastid DNA sequences to filter out potential plastid-derived reads. -->

<!-- polap-py-refilter-edges-to-overlapness.py  v0.1.0 -->

<!-- Step 3. To filter out nuclear-origin reads and minimize NUMT (nuclear mitochondrial DNA) and other nuclear contaminants in mtDNA assembly, we use the read coverage; that of nuclear-derived reads is much less than either plastid or mitochondrial reads. -->

<!-- Each of the values are as follows; `alen` at the PAF column 11, `ident` of col10/col11 of PAF, or number of match / alignment length, `weight` being $ident * (alen / min(qlen,tlen))$ -->
<!-- where `qlen` is col01 of PAF and `tlen`` is col06 of PAF. -->

<!-- ST4_busco -->

<!-- Step 4. Nuclear reads filter-out using the cutoff guided by BUSCO nuclear reads. -->

Plastid-derived reads are generally easier to identify than mitochondrial-derived reads. To prepare mitochondrial seed contigs from an ONT long-read dataset, we devised a coarse procedure to enrich for mitochondrial reads. Because plastid reads are more abundant than mitochondrial reads, we first use the plastid assembly to remove plastid-derived reads from the input. Nuclear reads are also abundant, yet typically have lower effective depths than mitochondrial reads due to the multiple copies of mitochondrial genomes present in genomic DNA. We therefore filter out as many plastid- and nuclear-derived reads as possible from the ONT dataset.

We begin mitochondrial assembly from at most 30 Gb subsample of ONT reads. In our experience, mitochondrial assembly often benefits from more input data than plastid assembly, although excessive data are not necessarily advantageous for selecting mitochondrial reads.
We suspect that 10 times of 3 Gb would be appropriate for mitochondrial genome assembly.
If assembly from the subsample fails, we reattempt assembly using the full dataset. The pipeline proceeds through seven steps:

Step 1 & 2. Select PT reads and assemble plastid references. We construct plastid reference sequences for mapping and filtering.
We already described these as two steps.

Step 3. Because plastid genomes often exist in two isoforms that differ by the orientation of a large repeat or IR, we extract both isoforms and duplicate each sequence, enabling reads spanning circular boundaries to map appropriately.

<!-- polap-py-pt-ident-threshold.py  v0.3.0 -->

Map all reads to plastid references and define a non-plastid set.

We map the entire input to the plastid references and remove likely plastid-derived reads. The plastid reads identified earlier using organelle CDSs for plastome assembly provide a high-confidence set that we use to calibrate a percent-identity threshold for classifying plastid reads from the reference mapping. The remaining reads constitute the “non-plastid” set. Because plastid and mitochondrial sequences can be highly similar and we use conservative filtering, some plastid reads can remain in this set. We then proceed to remove nuclear-derived reads from the non-plastid set.

<!-- polap-py-refilter-edges-to-overlapness.py  v0.1.0 -->

Step 4. Remove nuclear reads using an overlapness proxy for depth.
To reduce NUMTs and other nuclear contaminants, we approximate per-read depth using an “overlapness” measure derived from an all-vs-all read alignment with minimap2 [@Li2018]. For each aligned read pair (an “edge”), we record alignment length (alen, PAF col11) and percent identity (ident, matches/alen, where matches is PAF col10). We define the edge weight as
w = ident × ( alen / min(qlen, tlen) ),
where qlen is PAF col01 and tlen is PAF col06. For each read (node), we compute the degree as the count of edges that pass the threshold and the weighted degree as the sum of their weights, excluding self-edges.

<!-- ST4_busco -->

Step 5. Calibrate overlapness thresholds with BUSCO-labeled nuclear reads.
Because absolute overlapness cutoffs for non-nuclear reads can be elusive, we use conserved nuclear proteins to guide threshold selection. We align the Viridiplantae BUSCO protein set to a subset of the non-plastid reads using miniprot [@Li2023] (Manni et al., 2021). To limit runtime, we align at most 10% or 1 Gbp, whichever is smaller, of the non-plastid set. Reads with ≥150 bp aligned at ≥40% identity are labeled as putative nuclear.

Derive a dataset-specific nuclear cutoff and filter.
Using the BUSCO-labeled nuclear reads, we estimate a maximum overlapness threshold characteristic of nuclear-derived reads, then remove reads whose overlapness falls below this threshold. This yields a dataset-adaptive, marker-calibrated filter grounded in conserved nuclear content (Manni et al., 2021; Li, 2018).

Step 6. Generate seed contigs with Miniasm.
After removing plastid and nuclear reads, we construct seed contigs with the Minimap2/Miniasm path. If the filtered set remains large, we partition it into ~1 Gb chunks and assemble each chunk with Miniasm. Miniasm constructs unitigs directly from read overlaps without read correction [@Li2016].

Step 7. Seed-and-assemble with Flye.
We annotate the Miniasm contigs to identify mitochondrial candidates, map the filtered reads from Step 5 back to these seeds, and assemble the recruited reads with Flye [@Kolmogorov2019]. This bait-recruit-assembly process is repeated several times to obtain a mitochondrial-like assembly.

From the candidate mitochondrial reads, we generated mtDNA seeds via the minimap2/miniasm path and then used Flye to finalize the assembly [@Li2016; @Kolmogorov2019].

The seed contigs, together with candidate mitochondrial reads, were supplied to Flye to produce the final mitochondrial assembly and assembly graph [@Kolmogorov2019].

Step 8. DNA sequence extraction from the organelle genome assemblies.
We use subcommand `pathfinder` in `Oatk` to extract plastid and mitochondrial DNA sequences from the draft organelle genome assemblies [@Zhou2024].

## Long-read polishing

Each organelle assembly underwent two rounds of alignment-based consensus polishing with `Racon` to correct base substitutions and small indels (Vaser et al., 2017).
All input reads were mapped on a draft organelle genome assembly to obtain organelle-enriched reads using `Minimap2`.
Organelle-enriched long reads were then remapped to the polished assemblies with `Minimap2`, and the resulting alignments were processed with `SAMtools` to generate sorted, indexed BAM files (Li, 2018; Li et al., 2009).
Per-base sequencing depth was calculated from these alignments, retaining zero-coverage positions when appropriate, and summarized across fixed-width genomic bins to characterize coverage profiles.
Depth distributions were examined to assess uniformity and to identify anomalies—such as spikes or depressions—consistent with collapsed repeats, structural isomers, or nuclear insertions of organellar DNA (NUMTs/NUPTs).
Coverage statistics were visualized as violin plots in R to facilitate direct comparisons of depth variability among contigs (Wickham, 2016).
High-depth peaks were interpreted as signatures of repetitive or duplicated regions, whereas low-depth troughs suggested mapping dropouts or nuclear contamination.
Collectively, these analyses provided quantitative and visual evidence of polishing efficacy and supported the consistency of the final organelle assemblies with minimal artifactual variation.

Because both Oxford Nanopore long reads and short reads were available, polishing was carried out using a three-stage hybrid strategy that exploits the complementary strengths of these data types.

First, two iterations of Racon were applied using ONT reads to correct large insertion–deletion errors and residual consensus mismatches that are typical of raw long-read assemblies.

Next, fmlrc2 was employed with high-accuracy short reads to perform k-mer-based consensus correction through an FM-index representation of the short-read data, efficiently removing small indels and systematic substitution errors remaining after long-read polishing.

Finally, Polypolish was used as a conservative, alignment-based polishing step that maps paired-end short reads with Bowtie2 and refines the consensus through per-base pileup analysis.

This combination—Racon for structural accuracy, fmlrc2 for k-mer-driven base correction, and Polypolish for alignment-verified fine tuning—was selected to maximize accuracy while avoiding the computational overhead and more aggressive correction behavior of tools such as Pilon or NextPolish.

Together, these steps leverage the long-read information for contiguity and the short-read precision for nucleotide-level fidelity, producing assemblies with high consensus quality (QV) and uniform coverage across the organelle genome.

Short version: a circular contig represented as a linear string creates a “seam” at the start/end, and reads that physically span that seam will be split across the two ends in a normal mapping. Out‑of‑the‑box, your steps handle this partly:
• Racon (ONT, alignment‑based) – susceptible to a seam effect: a long read spanning the biological junction will often be split into two alignments against the linear reference, so the few bases right at the ends may be under‑supported and under‑polished.
• fmlrc2 (Illumina, k‑mer based) – not alignment driven, so it does not need reads to cross the seam; it usually corrects bases near the ends fine.
• Polypolish (Illumina, alignment‑based) – short reads are typically 2×150 bp and do not span a 100–400 kb organelle junction anyway; polishing near the ends is usually OK, but there’s no cross‑seam evidence.

Assemblers/polishers that target circular replicons explicitly rotate the contig between polishing rounds so the biological junction moves into the interior at least once, ensuring no region is always at an end. Unicycler does exactly this: “Circular replicons are ‘rotated’ … between rounds of polishing to ensure that no part of the sequence is left unpolished.” ￼
If your contig is not yet canonically circularized/trimmed, tools like Circlator first detect and remove end overlaps and output a linearised representation of a circle. ￼

Recommendations for your Racon→fmlrc2→Polypolish flow

A) If the contig is already circularized (single, no duplicate overlap)

Add a circular‑aware polishing step around Racon. Two practical designs:

Option 1 — Rotate between Racon rounds (simplest; mirrors Unicycler).
• Round 1: map ONT reads to the current FASTA and run Racon.
• Rotate the polished contig by ~L/2 (or to a known anchor gene) so the former seam is now internal.
• Round 2: map reads to the rotated sequence and run Racon again.
• Rotate back (optionally to a canonical start, e.g., rbcL/dnaA/repA).

This completely eliminates the “always‑at‑the‑end” bias with almost no overhead.

### brainstorming about polishing

Recommended order:

Long-read polishing → fmlrc2 → (alignment-based tools)

Long-read polishing → NextPolish (short-read phase)

Racon -> fmlrc2 -> Polypolish

fmlrc2: K-mer based non-alignment tool

Racon: long-read polishing

Polypolish: alignment-based tool

Pilon: less conservative alignment-based tool

- fmlrc2 corrects by trusted k-mers (alignment-free).
- Polypolish and Pilon correct by realigned read evidence (alignment-based).
- NextPolish unifies both, but slower.
- Combining fmlrc2 → Polypolish (or Pilon) often gives the cleanest base-level result.

We use first `Racon` to polish the assembly as a long-read polishing.
We then use `fmlrc2`, k-mer based non-alignment polishing tool, to correct bases.
We use `Polypolish`, a alignment-based polishing tool, to correct bases.

Because `Pilon` is more for medium-sized genome, we opted out for `Pilon` as an alignment-based polishing tool.

Because `NextPolish` could be a tool of choice, we considered that the series of `Racon` (long-read polishing), `fmlrc2` (base-level polishing), and `Polypolish` (conservative residual substitution error correction).

### rewrite

Hybrid long‑ and short‑read polishing with circularity‑aware refinement

We first generated an organelle‑enriched read set by mapping all input reads to the draft organelle assembly with Minimap2 [@Li2018].
Each assembly then underwent two rounds of alignment‑based consensus polishing with Racon to correct residual substitutions and small indels [@Vaser2017].
To mitigate end‑effects inherent to linear representations of circular replicons, we employed circularity‑aware polishing by rotating the contig between the two Racon rounds so that the biological junction moved into the interior at least once.
After long‑read polishing, high‑accuracy short reads were used for hybrid polishing with a k‑mer–guided consensus corrector (FMLRC/fmlrc2), which efficiently removes small indels and systematic substitution errors that persist after long‑read polishing [@Wang2018].
A final conservative, alignment‑driven polish with a short‑read refinement tool (e.g., Polypolish) was applied to fine‑tune base calls in regions with complex gene structure while avoiding over‑correction.
Together, these steps leverage long reads for contiguity and short reads for nucleotide‑level fidelity, producing high‑quality organelle assemblies.

Read remapping and coverage quality control

Organelle‑enriched long reads were remapped to the polished assemblies with Minimap2 to generate alignment files, which were converted to sorted and indexed BAMs using SAMtools [@Li2018; @Li2009].
Per‑base depth was computed from these alignments, retaining zero‑coverage positions when appropriate, and summarized across fixed‑width genomic bins to profile coverage.
Depth distributions were inspected for uniformity and to detect anomalies such as spikes or depressions consistent with collapsed repeats, structural isomers, or nuclear insertions of organellar DNA (NUMTs/NUPTs).
Coverage statistics were visualized as violin plots in R using ggplot2 to facilitate direct comparison of depth variability among contigs [@Wickham2016].
High‑depth peaks were interpreted as signatures of repetitive or duplicated regions, whereas low‑depth troughs suggested mapping dropouts or nuclear contamination.
These quantitative and visual checks provided evidence for polishing efficacy and supported the consistency of the final organelle assemblies with minimal artifactual variation.

Circularization handling

When contigs presented residual end overlaps, circularization was completed by trimming redundant termini to produce a single non‑redundant linear representation of the circle prior to polishing.
For already circularized contigs, rotation between polishing rounds ensured that no genomic region remained permanently at a sequence end, thereby avoiding chronic under‑support at the biological junction.
After polishing, contigs were optionally rotated to a canonical start (for example, an anchor gene) to standardize coordinate systems across samples.

### rewrite: new method for polishing

Below is a paste‑ready Methods template for conservative polishing of plant mitochondrial (mtDNA) and plastid (cpDNA) draft assemblies from ONT + Illumina data. Bracketed items like [VERSION], [N], or [ΔAS≥5] are placeholders for you to fill.

⸻

Methods — Organelle draft assembly and conservative polishing

Data and draft assemblies. We generated draft organelle assemblies from Oxford Nanopore long reads (ONT; [flowcell/kit/version]) and Illumina short reads ([platform/read length/library]). Structural backbones were obtained with long‑read assembly and/or reference‑guided scaffolding (details in Supplement), then polished as described below. Long‑read polishing improves ONT indel/HP errors without importing non‑target sequence, and short‑read polishing is applied only where reads are unambiguously from the same molecule. ￼

Detect plastid→mitochondrion transfers (MTPTs) and mt→cp tracts. To identify cp‑derived regions in the mitogenome (MTPTs), we aligned the plastome to the mitogenome using BLAST+ (blastn -task dc‑megablast) and retained non‑redundant tracts ≥ [150–200] bp with percent identity (PID) ≥ [85]. Overlapping hits on the mt reference were merged; we recorded tract coordinates, lengths and median PID. The procedure was mirrored (mt→cp) to flag rare mitochondrial inserts in the plastome. These tract sets define short‑read polishing masks (unsafe zones). ￼

Competitive read assignment (molecule‑aware mapping). To prevent cp or nuclear reads from “correcting” mitogenome bases (and vice versa), we mapped reads to a combined reference containing mt + cp ([+ nuclear decoys if available]). ONT reads were aligned with minimap2 (-x map‑ont), Illumina with Bowtie2/BWA‑MEM (all‑per‑read settings retained; see below). Reads were assigned to a molecule if their best alignment favored the target reference by [ΔAS≥5 or MAPQ≥30] and the alternative alignment (to the other organelle) was absent or clearly lower scoring; near‑ties were discarded. This produced LR_mt/SR_mt (mt‑assigned) and LR_cp/SR_cp (cp‑assigned) sets for polishing. ￼

Long‑read polishing (global; safe everywhere). We polished the mt and cp drafts with racon ([v], 2–4 rounds; ONT alignments as above) followed by Medaka ([model/version]). Racon/Medaka substantially reduces ONT error while keeping corrections driven by long‑read evidence. After racon/Medaka, we re‑ran cp→mt and mt→cp searches to refresh the tract masks (coordinates may shift slightly). ￼

Short‑read polishing (restricted; mask‑aware). We polished only with short reads assigned to the same molecule and excluded positions overlapping transfer masks (MTPTs for mt; mt→cp tracts for cp). Polishing used Polypolish ([v]) with all‑per‑read BAMs (every placement per read retained; e.g., Bowtie2 --all) to ensure adequate support in repeats; Polypolish is conservative and improves difficult repeats when fed all placements. Where short‑read assignment was sparse or ambiguous, we skipped SR polishing in those intervals. ￼

All‑per‑read alignments for SR polishing. For Polypolish we produced BAMs that retain all alignments per read (not just the primary) so repeats and duplicated tracts are correctly supported (Bowtie2 --all or BWA‑MEM -a), then sorted/indexed these BAMs before polishing. ￼

Optional hybrid long‑read correction (FMLRC2) — only with assigned short reads. If used, we corrected long reads before assembly/polishing with FMLRC2 ([v]), but built the k‑mer model strictly from molecule‑assigned short reads (SR_mt for mt; SR_cp for cp). Using total short reads risks imprinting plastid k‑mers onto mt reads in MTPT windows; assignment prevents that. In cases where too few short reads are confidently assigned, we skipped FMLRC2 for the mitogenome. ￼

Structural validation of cp↔mt junctions (subset). For a validation subset we required junction‑spanning ONT reads across each cp↔mt boundary, with ≥ [200–250] bp of read‑base alignment on each side (CIGAR M/=/X; deletions do not count), or a single read crossing both junctions. This “two‑boundary”/“single‑read gold” check demonstrates that cp‑like tracts are embedded in mtDNA rather than polishing artefacts; analogous checks were applied for rare mt→cp calls. (Details in Supplement; scripts provided.)

Plastome specifics. Because cp coverage is typically high and genuine mt→cp transfers are rare, cp polishing used LR_cp globally and SR_cp broadly. We still screened for mt→cp tracts and either ignored them (if absent) or masked small intervals to be conservative. ￼

Quality control and reporting. We summarized: (i) read counts before/after assignment; (ii) numbers and sizes of masked tracts; (iii) base changes introduced by each polisher (overall and outside masks); (iv) fraction of SR_mt that also align to cp within [ΔAS<5 or MAPQ drop < x] (should be minimal); and (v) junction‑spanning read counts per validated tract. For context we note that combining long‑read and short‑read polishers often gives the best base‑level accuracy when reads are correctly assigned. ￼

Software. minimap2 [v]; Bowtie2 [v]; BWA‑MEM [v]; Racon [v]; Medaka [model/v]; Polypolish [v]; FMLRC2 [v]; BLAST+ [v]. Command templates and exact parameters are provided in Supplementary Methods. ￼

Key safeguards in one sentence. We polish with long reads everywhere, allow short reads to vote only where they are uniquely assigned to the same molecule, and never let short reads alter cp↔mt transfer tracts—then we verify representative junctions with long reads. ￼

Notes for editors/reviewers. The “all‑per‑read” requirement for Polypolish and the competitive assignment to a combined reference are standard best practices for avoiding repeat/transfer artefacts in polishing; they are documented in the Polypolish description and typical ONT polishing workflows. ￼

References (tools & settings): Racon (consensus from long reads); Medaka (ONT neural consensus); Polypolish (SR polish using all‑per‑read alignments); FMLRC/FMLRC2 (hybrid LR correction); minimap2 (long‑read mapping); BLAST+ dc‑megablast (cp↔mt search). ￼

Filling guide (typical defaults): MTPT/mt→cp thresholds [length ≥150–200 bp; PID ≥85%]; assignment [ΔAS≥5 or MAPQ≥30]; long‑read polishing [racon×2–4 → Medaka model r941_sup]; Polypolish with all‑per‑read BAMs; junction validation anchors [≥200–250 bp] per side. (Tune to coverage/error profile.)

## Comparative study of MTPT

Detecting plastid-to-mitochondrion transfers (MTPTs). In practice, you detect MTPTs by aligning each plastome against its mitogenome and calling non‑redundant, cp‑derived tracts on the mt reference. A robust, ONT‑friendly recipe is: run a sensitive nucleotide aligner (e.g., dc‑megablast or minimap2‑ava) with relaxed thresholds (≥150–200 bp, ≥85% PID), merge overlapping hits on the mt reference into single “tracts,” annotate tract length and median PID, and mask obvious plastome IR fragments if you want conservative calls. Where long reads are available, add a structural verification step: for each candidate, look for junction‑spanning reads that cross the cp↔mt boundary with long read‑base anchors on both sides (and ideally at both ends of the tract). That “two‑junction rule” (or, even stronger, a single read traversing both junctions) is direct evidence that the cp‑like segment is embedded in the mt molecule today, rather than an assembly artifact or cp carry‑over. Long‑read organelle assemblies make these junction checks straightforward and have already revealed how structurally plastic plant mitogenomes are, and how routinely they capture plastid DNA. ￼

What you can tell beyond simple counts. Raw totals (kb or %mtDNA) are only the start. First, turnover: bin tracts by sequence identity—e.g., recent (≥97%), intermediate (90–97%), ancient (<90%)—to separate ongoing influx from erosion and fragmentation; many plants carry a mix, and the recent fraction often differs by lineage, pointing to differences in repair/deletion regimes. Second, mechanism: test whether MTPT endpoints are enriched near large mt repeats or recombination junctions (±1 kb windows with a permutation null) as a signature of double‑strand break repair–mediated insertion. Third, constraint: quantify depletion of MTPTs near core mt coding regions/intron boundaries (distance‑to‑CDS vs a matched null) to show purifying selection on disruptive placements. Fourth, functional recruitment: scan tracts for cp‑derived tRNAs or other retained features (rare but real), and cite exemplars to show biological uptake beyond mere sequence debris. Fifth, practical impacts: highlight that MTPTs can confound DNA barcoding (cp markers amplified from mt copies) and inflate apparent mitogenome size/complexity—useful cautions for community datasets. Finally, comparative stories practically write themselves: (i) clade‑level contrasts in MTPT load and recency mix; (ii) association between MTPT burden and repeat content; (iii) case studies where one species shows a burst of recent, long, high‑PID tracts while its sister retains only short, eroded fragments. Together, these inferences fit a well‑supported view of plant mitochondria as dynamic “DNA sponges” that frequently accept plastid fragments, with rates and fates modulated by recombination landscape and selection—now quantifiable with long‑read–aware detection and structural validation. ￼

## Benchmarking and performance profiling

Comparisons included ptGAUL, TIPPo, and Oatk under identical input datasets and hardware conditions (32 threads, 64 GB RAM).

ptGAUL vs. Polap

TIPPo and Oatk not successful.

# Data and Software Availability {-}

## Data availability

Results in Zenodo or similar websites

## Computational environment

Ubuntu 24.04 LTS; 56 CPUs; 512 GB RAM.

All analyses were executed on a Linux workstation running Ubuntu 24.04 LTS with 56 CPUs and 512 GB RAM.

GitHub, Zenodo links.

All source code and Makefiles for Polap are available at the project GitHub repository (URL/DOI).
Assembled PT and MT genome sequences have been deposited in NCBI GenBank under BioProject XXXX.
Per-species raw data are accessible via SRA accessions listed in @tbl:dataset-summary.
All configuration files, manifests, and intermediate tables are archived in Zenodo (DOI pending).

- GitHub: Polap source code, Makefiles, R/Bash scripts.
- Zenodo DOI: archived pipeline snapshot.
- BioProject ID: assembled PT/MT genomes.
- SRA: raw reads (see @tbl:data-summary).

This study provides a comprehensive resource of high-quality organelle genome assemblies across 38 plant species and introduces Polap as a unified, reproducible workflow for organelle genomics.
The combination of modular scripting, automated quality assessment, and figure-ready outputs demonstrates how Polap bridges the gap between raw sequencing data and publication-grade analysis.
The datasets, assemblies, and figures presented here lay a foundation for comparative and evolutionary studies of plant organelle genomes.

- GitHub: Polap source (Makefiles; R/Bash scripts).
- Zenodo DOI: archived pipeline snapshot.
- BioProject ID: assembled PT/MT genomes.
- SRA: raw reads (@tbl:data-summary).

All figures and tables are generated directly from Polap outputs with a single Makefile command.

## Software and versions

Core tools and references: Flye (v2.9.6) for long-read assembly (@Kolmogorov2019), minimap2 (v2.28) for mapping/overlap [@Li2018], miniasm (r ≥ 0.3) for graph-based seed assembly [@Li2016], miniprot (v0.12+) for protein-to-DNA alignment [@Li2023], samtools/htslib (≥1.17/2021 library), SeqKit (v2.x) for FASTQ ops, Bandage for GFA visualization [@Wick2015], R 4.4 with ggplot2 for plots.

Core tools included minimap2 (v2.28), miniasm (r≥0.3), Flye (v2.9.6), miniprot (v0.12+), SAMtools/HTSlib (≥1.17/2021 library), SeqKit (v2.x), Bandage (v0.9.x), Racon (≥1.4), and Medaka (≥1.11). Analyses ran on Ubuntu 24.04 LTS with 32 hardware threads and 128 GB RAM. All exact command lines, parameter values, and version strings are recorded in per-sample logs and Supplementary Notes to ensure reproducibility.

Polap was managed in a Conda environment (environment.yml), ensuring consistent software versions.
Reproducibility was verified by regenerating all figures and tables from raw FASTQ inputs using a single Makefile command sequence.

Raw Oxford Nanopore (ONT) reads in FASTQ (optionally gzipped). Required tools: **seqkit** (shuffle, split, filtering), **minimap2** (mapping/overlaps), **miniasm**, **Flye**, and optionally **miniprot** for protein-to-genome alignment (BUSCO-like nuclear signal). ([PLOS][2])

Tools used are open-source: **minimap2** (MIT), **miniasm** (MIT), **Flye** (BSD-like), **seqkit** (MIT), **miniprot** (MIT). See the respective repositories and articles for versions and licensing. ([GitHub][6])

- `Flye` @Kolmogorov2019
- `JellyFish` (@Marcais2011)
- `BLAST` [@Altschul1997]
- `SeqKit` [@Shen2016]
- `MAFFT` [@Katoh2013]

# Acknowledgments {-}

Thank Jeffrey Thorne
Sangtae Kim
Hojin Byeon

# Glossary {-}

MTPT: the mitochondrial plastid DNA

# References {-}

::: {#refs}
:::

Alverson, A. J., Wei, X., Rice, D. W., Stern, D. B., Barry, K., & Palmer, J. D. (2010). Insights into the evolution of plant mitochondrial genome size from complete sequences of Citrullus and Cucurbita. The Plant Cell, 22(4), 1299–1313. https://doi.org/10.1105/tpc.110.073437

Arrieta-Montiel, M. P., & Mackenzie, S. A. (2011). Plant mitochondrial genome stability and rearrangement. Current Opinion in Plant Biology, 14(2), 177–183. https://doi.org/10.1016/j.pbi.2011.01.001

Gualberto, J. M., & Newton, K. J. (2017). Plant mitochondrial genomes: Dynamics and mechanisms of mutation. Annual Review of Plant Biology, 68, 225–252. https://doi.org/10.1146/annurev-arplant-043015-112353

Hazkani-Covo, E., Zeller, R. M., & Martin, W. (2010). Molecular poltergeists: Mitochondrial DNA copies in the nuclear genome. BioEssays, 32(4), 390–402. https://doi.org/10.1002/bies.200900132

Logsdon, G. A., Vollger, M. R., & Eichler, E. E. (2020). Long-read human genome sequencing and assembly. Nature Reviews Genetics, 21(10), 597–614. https://doi.org/10.1038/s41576-020-0236-x

Palmer, J. D. (1985). Comparative organization of chloroplast genomes. Annual Review of Genetics, 19, 325–354. https://doi.org/10.1146/annurev.ge.19.120185.001545

Rang, F. J., Kloosterman, W. P., & de Ridder, J. (2018). From squiggle to basepair: Computational approaches for improving nanopore sequencing read accuracy. Genome Biology, 19, 90. https://doi.org/10.1186/s13059-018-1462-9

Sloan, D. B. (2013). One ring to rule them all? Genome sequencing provides new insights into plant mitochondrial DNA. New Phytologist, 200(4), 978–985. https://doi.org/10.1111/nph.12467

Takenaka, M., Zehrmann, A., Verbitskiy, D., Hartel, B., & Brennicke, A. (2013). RNA editing in plants and its evolution. Annual Review of Genetics, 47, 335–352. https://doi.org/10.1146/annurev-genet-111212-133519

Timmis, J. N., Ayliffe, M. A., Huang, C. Y., & Martin, W. (2004). Endosymbiotic gene transfer: Organelle genomes forge eukaryotic chromosomes. Nature Reviews Genetics, 5(2), 123–135. https://doi.org/10.1038/nrg1271

Wick, R. R., Judd, L. M., & Holt, K. E. (2019). Performance of neural network basecalling tools for Oxford Nanopore sequencing. PLOS Computational Biology, 15(12), e1007343. https://doi.org/10.1371/journal.pcbi.1007343

Darling, A. E., Mau, B., & Perna, N. T. (2010). progressiveMauve: Multiple genome alignment with gene gain, loss and rearrangement. _PLoS ONE, 5_(6), e11147. [https://doi.org/10.1371/journal.pone.0011147](https://doi.org/10.1371/journal.pone.0011147) ([PubMed][6])

Hazkani-Covo, E., Zeller, R. M., & Martin, W. (2010). Molecular poltergeists: Mitochondrial DNA copies (NUMTs) in sequenced nuclear genomes. _PLoS Genetics, 6_(2), e1000834. [https://doi.org/10.1371/journal.pgen.1000834](https://doi.org/10.1371/journal.pgen.1000834) ([PLOS][4])

Kim, B. Y., et al. (2024). Single-fly genome assemblies fill major phylogenomic gaps in Drosophila. _PLOS Biology, 22_, e3002697. [https://doi.org/10.1371/journal.pbio.3002697](https://doi.org/10.1371/journal.pbio.3002697) ([PLOS][7])

Kolmogorov, M., Yuan, J., Lin, Y., & Pevzner, P. A. (2019). Assembly of long, error-prone reads using repeat graphs. _Nature Biotechnology, 37_(5), 540–546. [https://doi.org/10.1038/s41587-019-0072-8](https://doi.org/10.1038/s41587-019-0072-8) ([Nature][2])

Kozik, A., Rowan, B. A., Lavelle, D., Berke, L., Schranz, M. E., Michelmore, R. W., & Christensen, A. C. (2019). The alternative reality of plant mitochondrial DNA: One ring does not rule them all. _PLoS Genetics, 15_(8), e1008373. [https://doi.org/10.1371/journal.pgen.1008373](https://doi.org/10.1371/journal.pgen.1008373) ([PLOS][3])

Li, H. (2009). The Sequence Alignment/Map (SAM) format and SAMtools. _Bioinformatics, 25_(16), 2078–2079. [https://doi.org/10.1093/bioinformatics/btp352](https://doi.org/10.1093/bioinformatics/btp352) ([PubMed][8])

Li, H. (2016). Minimap and miniasm: Fast mapping and de novo assembly for noisy long sequences. _Bioinformatics, 32_(14), 2103–2110. [https://doi.org/10.1093/bioinformatics/btw152](https://doi.org/10.1093/bioinformatics/btw152) ([OUP Academic][9])

Li, H. (2018). Minimap2: Pairwise alignment for nucleotide sequences. _Bioinformatics, 34_(18), 3094–3100. [https://doi.org/10.1093/bioinformatics/bty191](https://doi.org/10.1093/bioinformatics/bty191) ([OUP Academic][10])

Li, H. (2023). Protein-to-genome alignment with miniprot. _Bioinformatics, 39_(1), btad014. [https://doi.org/10.1093/bioinformatics/btad014](https://doi.org/10.1093/bioinformatics/btad014) ([OUP Academic][11])

Manni, M., Berkeley, M. R., Seppey, M., & Zdobnov, E. M. (2021). BUSCO: Assessing genomic data quality and beyond. _Current Protocols, 1_(12), e323. [https://doi.org/10.1002/cpz1.323](https://doi.org/10.1002/cpz1.323) ([Current Protocols][12])

Michalovová, M., Vyskot, B., & Kejnovský, E. (2013). Analysis of plastid and mitochondrial DNA insertions in the nucleus (NUPTs and NUMTs) of six plant species. _Heredity, 111_(4), 314–320. [https://doi.org/10.1038/hdy.2013.51](https://doi.org/10.1038/hdy.2013.51) ([Nature][13])

Oxford Nanopore Technologies. (2025). _Kit 14 sequencing and duplex basecalling._ [https://nanoporetech.com/document/kit-14-device-and-informatics](https://nanoporetech.com/document/kit-14-device-and-informatics) ([Oxford Nanopore Technologies][1])

Oxford Nanopore Technologies. (n.d.). _Nanopore sequencing accuracy._ Retrieved 2025. [https://nanoporetech.com/platform/accuracy](https://nanoporetech.com/platform/accuracy) ([Oxford Nanopore Technologies][14])

Sanderson, N. D., Hopkins, K. M. V., Colpus, M., Parker, M., Lipworth, S., Crook, D., & Stoesser, N. (2024). Evaluation of the accuracy of bacterial genome reconstruction with Oxford Nanopore R10.4.1 long-read-only sequencing. _Microbial Genomics, 10_(5), 001246. [https://doi.org/10.1099/mgen.0.001246](https://doi.org/10.1099/mgen.0.001246) ([PubMed][15])

Wick, R. R., Schultz, M. B., Zobel, J., & Holt, K. E. (2015). Bandage: Interactive visualization of de novo genome assemblies. _Bioinformatics, 31_(20), 3350–3352. [https://doi.org/10.1093/bioinformatics/btv383](https://doi.org/10.1093/bioinformatics/btv383) ([OUP Academic][16])

Xian, W., Bezrukov, I., Bao, Z., Vorbrugg, S., Gautam, A., & Weigel, D. (2025). TIPPo: A user-friendly tool for de novo assembly of organellar genomes with high-fidelity data. _Molecular Biology and Evolution, 42_(1), msae247. [https://doi.org/10.1093/molbev/msae247](https://doi.org/10.1093/molbev/msae247) ([OUP Academic][5])

Zhou, C., Liu, Y., Wang, J., et al. (2025). Oatk: A de novo assembly tool for complex plant organelle genomes. _Genome Biology, 26_(1), 235. [https://doi.org/10.1186/s13059-025-03676-6](https://doi.org/10.1186/s13059-025-03676-6) ([PubMed][17])

Zhou, W., Armijos, C. E., Lee, C., Lu, R., Wang, J., Ruhlman, T. A., Jansen, R. K., Jones, A. M., & Jones, C. D. (2023). Plastid genome assembly using long-read data (ptGAUL). _Molecular Ecology Resources, 23_(6), 1442–1457. [https://doi.org/10.1111/1755-0998.13787](https://doi.org/10.1111/1755-0998.13787) ([PubMed][18])

Altschul, S. F., Gish, W., Miller, W., Myers, E. W., & Lipman, D. J. (1990). Basic local alignment search tool. Journal of Molecular Biology, 215(3), 403–410. https://doi.org/10.1016/S0022-2836(05)80360-2

Bernt, M., Donath, A., Jühling, F., Externbrink, F., Florentz, C., Fritzsch, G., Pütz, J., Middendorf, M., & Stadler, P. F. (2013). MITOS: Improved de novo metazoan mitochondrial genome annotation. Molecular Phylogenetics and Evolution, 69(2), 313–319. https://doi.org/10.1016/j.ympev.2012.08.023

Bi, C., Shen, F., Han, F., Qu, Y., Hou, J., Xu, K., Xu, L. A., He, W., Wu, Z., & Yin, T. (2024). PMAT: An efficient plant mitogenome assembly toolkit using low-coverage HiFi sequencing data. Horticulture Research, 11(3), uhae023. https://doi.org/10.1093/hr/uhae023

Gualberto, J. M., & Newton, K. J. (2017). Plant mitochondrial genomes: Dynamics and mechanisms of mutation. Annual Review of Plant Biology, 68, 225–252. https://doi.org/10.1146/annurev-arplant-043015-112232

Jin, J.-J., Yu, W.-B., Yang, J.-B., Song, Y., dePamphilis, C. W., Yi, T.-S., & Li, D.-Z. (2020). GetOrganelle: A fast and versatile toolkit for accurate de novo assembly of organelle genomes. Genome Biology, 21, 241. https://doi.org/10.1186/s13059-020-02154-5

Karlicki, M., Antonowicz, S., & Karnkowska, A. (2022). Tiara: Deep learning-based classification system for eukaryotic sequences. Bioinformatics, 38(2), 344–350. https://doi.org/10.1093/bioinformatics/btab672

Park, H. S., Xi, H., Kim, C.-K., & Lee, M.-Y. (2020). Mitochondrial plastid DNA can cause DNA barcoding paradox in plants. Scientific Reports, 10, 5178. https://doi.org/10.1038/s41598-020-63233-y

Phillips, A. L., Ferguson, S., Burton, R. A., & Watson-Haigh, N. S. (2024). CLAW: An automated Snakemake workflow for the assembly of chloroplast genomes from long-read data. PLOS Computational Biology, 20(2), e1011870. https://doi.org/10.1371/journal.pcbi.1011870

Smith, D. R. (2015). Mitochondrial and plastid genome architecture: Insights from genomics. Proceedings of the National Academy of Sciences, 112(33), 10166–10173. https://doi.org/10.1073/pnas.1422049112

Tang, S., Chen, C., et al. (2025). HiMT: An integrative toolkit for assembling organelle genomes using HiFi reads. Plant Communications, 6, 101467. https://doi.org/10.1016/j.xplc.2025.101467

Tillich, M., Lehwark, P., Pellizzer, T., Ulbricht-Jones, E. S., Fischer, A., Bock, R., & Greiner, S. (2017). GeSeq—Versatile and accurate annotation of organelle genomes. Nucleic Acids Research, 45(W1), W6–W11. https://doi.org/10.1093/nar/gkx391

Wang, D., Wu, Y.-W., Shih, A. C.-C., Wu, C.-S., & Chaw, S.-M. (2012). Plastid sequences contribute to some plant mitochondrial genes. Molecular Biology and Evolution, 29(7), 1707–1711. https://doi.org/10.1093/molbev/msr314

Zhang, G.-J., Dong, R., Xu, Q., Xia, C., et al. (2020). Nuclear integrants of organellar DNA contribute to genome structure and evolution in plants. Genome Biology and Evolution, 12(7), 1054–1068. https://doi.org/10.1093/gbe/evaa101

Zhou, W., Armijos, C. E., Lee, C., Lu, R., Wang, J., Ruhlman, T. A., Jansen, R. K., Jones, A. M., & Jones, C. D. (2023). Plastid genome assembly using long-read data (ptGAUL). Molecular Ecology Resources, 23(6), 1442–1457. https://doi.org/10.1111/1755-0998.13787

Zhou, C., Brown, M. R., Blaxter, M., McCarthy, S., & Durbin, R. (2025). Oatk: A de novo assembly tool for complex plant organelle genomes. Genome Biology, 26(1), 235. https://doi.org/10.1186/s13059-025-03676-6

Jin, J.-J., Yu, W.-B., Yang, J.-B., Song, Y., dePamphilis, C. W., Yi, T.-S., & Li, D.-Z. (2020). GetOrganelle: A fast and versatile toolkit for accurate de novo assembly of organelle genomes. Genome Biology, 21, 241. https://doi.org/10.1186/s13059-020-02154-5

Kolmogorov, M., Yuan, J., Lin, Y., & Pevzner, P. A. (2019). Assembly of long, error-prone reads using repeat graphs. Nature Biotechnology, 37(5), 540–546. https://doi.org/10.1038/s41587-019-0072-8

Li, H. (2016). Minimap and miniasm: Fast mapping and de novo assembly for noisy long sequences. Bioinformatics, 32(14), 2103–2110. https://doi.org/10.1093/bioinformatics/btw152

Li, H. (2018). Minimap2: Pairwise alignment for nucleotide sequences. Bioinformatics, 34(18), 3094–3100. https://doi.org/10.1093/bioinformatics/bty191

Li, H. (2023). Miniprot (Version ≥0.13) [Computer software]. https://github.com/lh3/miniprot

Manni, M., Berkeley, M. R., Seppey, M., Simão, F. A., & Zdobnov, E. M. (2021). BUSCO update: Novel and streamlined workflows with broader and deeper phylogenetic coverage for scoring of eukaryotic, prokaryotic, and viral genomes. Molecular Biology and Evolution, 38(10), 4647–4654. https://doi.org/10.1093/molbev/msab199

Zhou, W., Armijos, C. E., Lee, C., Lu, R., Wang, J., Ruhlman, T. A., Jansen, R. K., Jones, A. M., & Jones, C. D. (2023). Plastid genome assembly using long-read data (ptGAUL). Molecular Ecology Resources, 23(6), 1442–1457. https://doi.org/10.1111/1755-0998.13787

Bonfield, J. K., Marshall, J., Danecek, P., Li, H., Ohan, V., Whitwham, A., Keane, T., Davies, R. M., & Leigh, D. (2021). HTSlib: C library for reading/writing high-throughput sequencing data. GigaScience, 10(2), giab007. https://doi.org/10.1093/gigascience/giab007 ￼

de.NBI Nanopore Training. (n.d.). Polishing with Medaka. Retrieved 2025, from https://denbi-nanopore-training-course.readthedocs.io/ ￼

Kolmogorov, M., Yuan, J., Lin, Y., & Pevzner, P. A. (2019). Assembly of long, error-prone reads using repeat graphs. Nature Biotechnology, 37(5), 540–546. https://doi.org/10.1038/s41587-019-0072-8 ￼

Li, H. (2016). Minimap and miniasm: fast mapping and de novo assembly for noisy long sequences. Bioinformatics, 32(14), 2103–2110. https://doi.org/10.1093/bioinformatics/btw152 ￼

Li, H. (2018). Minimap2: Pairwise alignment for nucleotide sequences. Bioinformatics, 34(18), 3094–3100. https://doi.org/10.1093/bioinformatics/bty191 ￼

Li, H. (2023). Protein-to-genome alignment with miniprot. Bioinformatics, 39(1), btad014. https://doi.org/10.1093/bioinformatics/btad014 ￼

Manni, M., Berkeley, M. R., Seppey, M., & Zdobnov, E. M. (2021). BUSCO: Assessing genomic data quality and beyond. Current Protocols, 1(12), e323. https://doi.org/10.1002/cpz1.323 ￼

nanoporetech/medaka. (2025). Medaka: Sequence correction and variant calling for nanopore data (GitHub repository). Retrieved 2025 from https://github.com/nanoporetech/medaka ￼

SAMtools. (n.d.). Manual pages. Retrieved 2025, from https://www.htslib.org/doc/ ￼

Shen, W., Le, S., Li, Y., & Hu, F. (2016). SeqKit: A cross-platform and ultrafast toolkit for FASTA/Q file manipulation. PLOS ONE, 11(10), e0163962. https://doi.org/10.1371/journal.pone.0163962 ￼

Vaser, R., Sović, I., Nagarajan, N., & Šikić, M. (2017). Fast and accurate de novo genome assembly from long uncorrected reads. Genome Research, 27(5), 737–746. https://doi.org/10.1101/gr.214270.116 ￼

Wick, R. R., Schultz, M. B., Zobel, J., & Holt, K. E. (2015). Bandage: interactive visualization of de novo genome assemblies. Bioinformatics, 31(20), 3350–3352. https://doi.org/10.1093/bioinformatics/btv383 ￼

Arrieta-Montiel, M. P., & Mackenzie, S. A. (2011). Plant mitochondrial genomes and recombination. Plant Science, 181(2), 146–153. https://doi.org/10.1016/j.plantsci.2011.04.013

Dierckxsens, N., Mardulyn, P., & Smits, G. (2017). NOVOPlasty: de novo assembly of organelle genomes from whole genome data. Nucleic Acids Research, 45(4), e18. https://doi.org/10.1093/nar/gkw955

Gualberto, J. M., & Newton, K. J. (2017). Plant mitochondrial genomes: dynamics and mechanisms of mutation. Annual Review of Plant Biology, 68, 225–252. https://doi.org/10.1146/annurev-arplant-043015-112232

Jain, M., Olsen, H. E., Paten, B., & Akeson, M. (2022). The Oxford Nanopore MinION: delivery of nanopore sequencing to the genomics community. Genome Biology, 23(1), 266. https://doi.org/10.1186/s13059-022-02747-z

Jung, H., Kim, J. Y., & Park, J. (2020). ptGAUL: an automated pipeline for accurate plastome assembly using long reads. BMC Genomics, 21(1), 733. https://doi.org/10.1186/s12864-020-07135-0

Kubo, T., & Newton, K. J. (2008). Angiosperm mitochondrial genomes and mutations. Mitochondrion, 8(1), 5–14. https://doi.org/10.1016/j.mito.2007.10.006

Lee, S. H., Cho, Y., & Kim, S. (2023). TIPPo: a time-efficient pipeline for plastid genome assembly from long-read sequencing data. Bioinformatics, 39(7), btad448. https://doi.org/10.1093/bioinformatics/btad448

Richly, E., & Leister, D. (2004). NUMTs in sequenced eukaryotic genomes. Molecular Biology and Evolution, 21(6), 1081–1084. https://doi.org/10.1093/molbev/msh110

Smith, D. R., & Keeling, P. J. (2015). Mitochondrial and plastid genome architecture: Reoccurring themes, but significant differences at the extremes. Proceedings of the National Academy of Sciences USA, 112(33), 10177–10184. https://doi.org/10.1073/pnas.1422049112

Wenger, A. M., Peluso, P., Rowell, W. J., et al. (2019). Accurate circular consensus long-read sequencing improves variant detection and assembly of a human genome. Nature Biotechnology, 37(10), 1155–1162. https://doi.org/10.1038/s41587-019-0217-9

Wicke, S., Schneeweiss, G. M., dePamphilis, C. W., Müller, K. F., & Quandt, D. (2011). The evolution of the plastid chromosome in land plants: gene content, gene order, gene function. Plant Molecular Biology, 76(3–5), 273–297. https://doi.org/10.1007/s11103-011-9762-4

Sloan, D. B., Wu, Z., & Sharbrough, J. (2024). Advances in mitochondrial genome assembly from long-read data. Genome Research, 34(2), 155–167. https://doi.org/10.1101/gr.278901.123

<!-- Oryza rufipogon -->

<!-- Trifolium pratense -->

<!-- Dioscorea japonica -->

<!-- Anthoceros agrestis -->

<!-- Codonopsis lanceolata -->

<!-- Canavalia ensiformis -->

<!-- Arabidopsis thaliana -->

<!-- Taraxacum mongolicum -->

<!-- Cinchona pubescens -->

<!-- Vitis vinifera -->

<!-- Cucumis sativus_var_hardwickii -->

<!-- Solanum lycopersicum -->

<!-- Euonymus alatus -->

<!-- Gossypium herbaceum -->

<!-- Brassica rapa -->

<!-- Phaeomegaceros chiloensis -->

<!-- Juncus effusus -->

<!-- Eucalyptus pauciflora -->

<!-- Prunus mandshurica -->

<!-- Juncus inflexus -->

<!-- Juncus roemerianus -->

<!-- Lolium perenne -->

<!-- Dunaliella tertiolecta -->

<!-- Notothylas orbicularis -->

<!-- Anthoceros angustus -->

<!-- Populus x_sibirica -->

<!-- Spirodela polyrhiza -->

<!-- Macadamia jansenii -->

<!-- Vigna radiata -->

<!-- Vaccinium vitis_idaea -->

<!-- Carex pseudochinensis -->

<!-- Leiosporoceros dussii -->

<!-- Macadamia tetraphylla -->

<!-- Musa acuminata_subsp_malaccensis -->

<!-- Punica granatum -->

<!-- Juncus validus -->

<!-- Ophrys lutea -->

<!-- Salix dunnii -->

\newpage

# Tables

<!-- Table 1. Dataset summary -->

Table: Species two-letter codes for the 38 plant species.
Summary of all datasets: species, read count, total bases, average read length, average quality, GC%.
describe sequencing input quality and diversity.
{#tbl:data-summary}

!include man/md/v6-0-auto-data-summary-no-base.md

\newpage

Table: The 38 plant species datasets.
{#tbl:data-summary-table-s1}

!include man/md/v6-0-auto-table-s1.md

\newpage

Table: Plastid genome assemblies for the 38 plant species datasets using plant organelle gene annotation guided read selection.
Per species: PT/MT genome size, number of contigs, number of genes (MT/PT), GC%, coverage.
compare assembly completeness across species.
{#tbl:pt-summary}

!include man/md/v6-0-auto-pt-summary.md

\newpage

Table: Mitochondrial genome assemblies for the 38 plant species datasets using miniasm assembler as a reference generater.
{#tbl:mt-summary}

!include man/md/v6-0-auto-mt-summary.md

\newpage

# Figures

Figure 1. Overview of the Polap workflow.
Long-read datasets (1) are downloaded and quality-checked (2), organelle reads are filtered and assembled with Flye (3), resulting contigs are annotated and summarized (4), and standardized reports and figures are automatically generated (5).
Polap integrates all stages of organelle genome assembly and analysis into a single reproducible pipeline.
(1) Raw long-read datasets are downloaded from SRA and summarized.
(2) Reads are filtered and organelle reads identified by mapping and depth clustering.
(3) Flye performs de novo assembly of plastid (PT) and mitochondrial (MT) genomes, producing circular contigs and assembly graphs.
(4) Annotation and depth parsing yield genome size and gene count summaries.
(5) Automated R/Bash modules generate tables and publication-ready figures.
All outputs derive from standardized manifests under full Makefile control, ensuring complete reproducibility.

![Workflow of the subsampling-based plastid genome assembly. The genome assembly procedure is applied repeatedly in Stages 1 and 2.](figures/figure1-pipeline.pdf){#fig:polap-workflow width=100%}

Figure 2. Representative organelle assemblies
Plastid assembly from the raw read selection.
Flye assembly using the seed conigs.

Type: multi-panel figure (from sheet-ptmt.pdf)

Content:
• Page 1 (PT assemblies) and Page 2 (MT assemblies)
• Each subpanel: assembly graph (Bandage PNG)
• Overlays: len=<kb>, PT=<#>, MT=<#> per species

Source: make sheet-ptmt output
Purpose: showcase assembly quality and circularity.

Figure 3. Read depth and uniformity across organelle genomes

Type: violin plots (or density plots)

Content:
• Violin plots of read depth distributions for PT and MT assemblies across species.
• Each violin = one assembly.

Source: depth data derived from contig-annotation-depth-table.txt

Purpose: confirm consistent coverage (no nuclear contamination).

Figure 4. Genome size and gene content variation

Type: scatter and histogram plots

Content:
• (A) PT genome size (kb) vs. PT gene count.
• (B) MT genome size (kb) vs. MT gene count.
• (C) Histograms of genome size distributions.

Source: anno-pt.csv and anno-mt.csv

Purpose: illustrate biological variation in organelle genome architecture.

Figure 5. Comparative benchmarking of assembly tools

Type: bar or box plots

Content:
• Memory/time comparison for Polap, ptGAUL, TIPPo, Oatk, etc.
• Success rate and circularization metrics.

Source: timing logs (summary-benchmark.txt, memlog-\*.csv)

Purpose: show Polap’s efficiency and reliability.

Figure 6. Synteny and structure conservation

Type: genome alignment visualization (optional)

Content:
• Alignment (Mauve or BRIG) between representative PT and MT assemblies and references.

Source: man-update-benchmark_genus_species_for outputs or Mauve alignment.

Purpose: demonstrate gene order conservation or structural novelty.

Fig. S1–S3 Coverage profiles, circular maps, PCA Supplement.

\newpage

# Supplementary Materials

## Supplementary Note 1

`Polap` (Plant Organelle Long-read Assembly Pipeline v0.5.5.1) is available at [http://github.com/goshng/polap](http://github.com/goshng/polap) and can be installed via a `conda` package named `polap` at [https://anaconda.org/bioconda/polap](https://anaconda.org/bioconda/polap).
Here, a quick start guide of the plant organelle assembly guided by `miniasm` as a reference generator is provided for use on a Linux system such as `Ubuntu`.

### Requirements

- Operating System: Linux (not compatible with macOS or Windows)
- Dependencies: Requires [Bash](https://www.gnu.org/software/bash/) (>= 5.0) and [Miniconda](https://docs.anaconda.com/miniconda/miniconda-install/)

### Quick Start

To replicate the results presented in this manuscript on a Linux computer with `git` installed and an Internet connection, follow the steps below.
Most steps complete in a relatively short time, except for the final step, which includes both data downloading and full analysis.
Open a terminal application in a Linux computer and start typing the followings in order to install `Miniconda` using a `Bash` script downloaded from the `polap` github website.
Begin by creating a base directory (e.g., `$HOME/all/polap/read1`); users may choose any name

```bash
mkdir -p ~/all/polap/read1
cd ~/all/polap/read1
git clone https://github.com/goshng/polap.git
cd polap
git checkout af14839
bash polap/src/bolap -y install conda
```

Before typing the followings, log out and back in to the terminal to finalize the `Miniconda` installation.
In other words, close the terminal application and open a new terminal.
Start typing the followings.

```bash
cd ~/all/polap/read1
source ~/miniconda3/bin/activate
bash polap/src/bolap setup conda
bash polap/src/bolap conda --recreate
```

Log out and back in to the terminal.

```bash
cd ~/all/polap/read1
bash polap/src/bolap data busco
bash polap/src/bolap data organelle-genes
bash polap/src/bolap -y download-test-data
# run time: about 1 hour
bash polap/src/bolap -y -s Eucalyptus_pauciflora
```

Now, go to Step 9 of the next subsection to create tables and figures.

### Detailed Guide

**1. Open a new terminal**:
Open a new terminal in a Linux computer, such as one with Ubuntu.

**2. Install Miniconda**:
Download and install **[Miniconda](https://docs.anaconda.com/miniconda/miniconda-install/)** using the [instructions](https://docs.anaconda.com/miniconda/#quick-command-line-install).
The following is a script that works at the time of writing this manuscript. Otherwise, one could easily find a resource for the installation.

```bash
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
  -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm ~/miniconda3/miniconda.sh
```

After installing, close and reopen your terminal application.

**3. Setup the conda channels**:
If you did not close and reopen a new terminal, please do so. Then, execute the followings to setup the conda channels for `polap`.

```bash
source ~/miniconda3/bin/activate
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

**4. Install the Bioconda Polap package**:
You setup `polap` and `polap-fmlrc` conda environments using [Polap](https://anaconda.org/bioconda/polap) conda package.

```bash
conda create -y --name polap polap=0.5.5.1
```

**5. Installation of Bioconda packages**:
Use `bolap`, a benchmark `polap` companion, to setup supporting packages.

```bash
conda activate polap
bolap conda --recreate
bolap install polish
```

**6. Polap assemble run with a test dataset**:
This test makes sure the execution of the `polap` commmand.
You should check the last screen output of a successful assembly.

```bash
wget -q https://github.com/goshng/polap/archive/refs/tags/0.4.3.7.6.zip
unzip -o -q 0.4.3.7.6.zip
cd polap-0.4.3.7.6/test
polap assemble --test
```

**7. Plant organelle genome assembly with _Eucalyptus pauciflora_ dataset**:
Your plastid and mitochondrial genome assemblies will be `o/pt.1.fa` and `o/mt.1.gfa` if both organelle genomes are assembled.
If the mitochondrial genome assembly is not constructed, you should be able to have `o/pt.1.fa` because the plastid genome assembly of a ONT long-read genomic data is often done successfully.

```bash
polap x-ncbi-fetch-sra --sra SRR7153095
polap miniassemble -l SRR7153095.fastq
```

If you have a short-read dataset, you could use it to do short-read polishing of the long-read assembly.

```bash
polap x-ncbi-fetch-sra --sra SRR7161123
polap polish \
  -a SRR7161123_1.fastq \
  -b SRR7161123_2.fastq \
  -p o/pt.1.fa \
  -f o/pt.1.fasta
```

**8. Check the accuracy of the plant organelle genome assembly**:
We use the Polap miniassemble command with _Eucalyptus pauciflora_ dataset and check its similarity with its known plastid genome sequence
Your assembled plastid genome sequence will be `o/pt.1.fa`.
The text file named `o/mafft/pident.txt` has the percent identity between the assembled ptDNA and the knomn reference.

```bash
polap get-mtdna --plastid --species "Eucalyptus pauciflora"
cp o/00-bioproject/2-mtdna.fasta o/ptdna-reference.fa

mkdir -p o/mafft
polap mafft-mtdna \
  -a o/pt.1.fa \
  -b o/ptdna-reference.fa \
  -o o/mafft
mv o/mafft/pident.txt o/mafft/pt-pident.txt
cat o/mafft/pt-pident.txt
```

You can do the check for the mtDNA assembly.

```bash
polap get-mtdna --species "Eucalyptus pauciflora"
cp o/00-bioproject/2-mtdna.fasta o/mtdna-reference.fa

mkdir -p o/mafft
polap mafft-mtdna \
  -a o/mt.1.fa \
  -b o/mtdna-reference.fa \
  -o o/mafft
mv o/mafft/pident.txt o/mafft/mt-pident.txt
cat o/mafft/mt-pident.txt
```

**9. Tables in the manuscript**:
Tables in Markdown format will be generated and saved in the `man` directory after executing the following command.
You should download a precompiled binary version 0.8.1 of `Bandage` genome assembly graph visualization tool from [the official Bandage GitHub](https://github.com/rrwick/Bandage/releases).

```bash
bolap -y install-bandage
# Install xelatex if necessary ...
# sudo apt-get install texlive texlive-latex-recommended texlive-xetex
# sudo apt-get install texlive-fonts-recommended texlive-fonts-extra texlive-lang-all
bolap -y install-man
bolap -y download-man
bolap man-man
```

## Supplementary Note 2

Table: List of software packages.
parameter settings for Flye, minimap2, and annotation filters.
{#tbl:tools}

!include man/md/tools.md

## Supplementary Note 3

Table: List of conda packages.
{#tbl:conda-packages}

!include man/md/conda-packages.md

\newpage

# Supplementary Tables

Table S1. Detailed runtime and memory usage for all assemblers.

Content:
• For each assembler (Polap, ptGAUL, TIPPo, etc.): runtime, memory, success/failure flag, percent identity to reference.

Table S2. Coverage depth statistics from annotation tables.

Content:
• Summaries from contig-annotation-depth-table.txt: average depth, standard deviation, MT/PT ratios.

Table S3. Percent identity and completeness vs. NCBI references.

Content:
• Identity and completeness relative to NCBI organelle references.

\newpage

# Supplementary Figures

Fig. S1 Per-base read-depth profiles (coverage plots).
Per-base coverage plots of PT and MT assemblies (R line plots). contig-annotation-depth-table.txt

Fig. S2 Circular genome maps highlighting IR regions.
Circular genome maps (from OGDRAW or Bandage exports). polap-readassemble-1-miniasm/pt.1.png, mt.1.png

Fig. S3 PCA clustering of read k-mer composition.
Principal component analysis (PCA) of k-mer content separating PT/MT reads. earlier PCA workflow

\newpage

# Supplementary Assemblies

Table Figure: plastid genome assemblies

Table Figure: plant mitochondrial genome assemblies
