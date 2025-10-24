# Plant Organelle Genome Assemblies — Review

---

## 1. Data Types and Challenges

- **Oxford Nanopore Technologies (ONT)**

  - ✅ Very long reads — span repeats
  - ❌ Historically lower accuracy [Wick et al. 2019]
  - ❌ Difficult for mtDNA: repeats, NUMTs, NUPTs

- **PacBio HiFi (CCS)**
  - ✅ High accuracy (~99.9%) [Wenger et al. 2019]
  - ✅ Good for repeat resolution
  - ❌ Shorter than ONT
  - ❌ More expensive

::: notes
ONT produces very long reads, helpful for resolving repeats, but accuracy has been lower compared to PacBio.  
HiFi reads are shorter but extremely accurate, and this usually makes assemblies cleaner.  
:::

---

## 2. Organelle Genomes

- **Plastid DNA (ptDNA)**

  - ~150 kb, quadripartite structure (LSC, SSC, IRs) [Palmer 1985; Daniell et al. 2016]
  - Relatively conserved, easier to assemble

- **Mitochondrial DNA (mtDNA)**
  - Hundreds of kb to several Mb [Sloan 2013]
  - Highly complex repeats and dynamic structures
  - Nuclear/ptDNA contamination complicates assemblies [Hazkani-Covo et al. 2010]

::: notes
Plastids have conserved structure and are relatively easy to assemble.  
Plant mitochondria are much larger, dynamic, and often contaminated with nuclear/ptDNA reads.  
:::

---

## 3. Current Strategies in `polap`

### ptDNA assembly

- **ONT:** protein-guided read selection -> Flye v2.9.6 [Kolmogorov et al. 2019]
- **HiFi:** protein-guided read selection -> Flye  
  OR use **Oatk** [Choi et al. 2025, in prep/published]

### mtDNA assembly

- **ONT:** remove ptDNA/nuclear reads -> seed contigs (Miniasm [Li 2016]) -> Flye
- **HiFi:** assemble directly with **Oatk**

::: notes
Flye is robust for noisy long reads.  
Miniasm provides quick seed contigs.  
Oatk, developed for HiFi, handles complex organellar graphs more gracefully.  
:::

---

## 4. Benchmarking with `bolap`

- Provides standardized tests:
  - `bolap run ptgaul`, `mtgaul`, `oatk`, `tippo`
  - `bolap run polap-readassemble-pt` / `nt`
- Enables comparison of:
  - ONT vs HiFi assemblies
  - Alternative assemblers
  - Assembly quality

::: notes
`bolap` makes results reproducible and comparable across datasets and pipelines.  
It is useful for deciding which tool or dataset is giving the most complete assemblies.  
:::

---

## 5. Complementary Tools

- **Read selection**

  - Minimap2 [Li 2018] / Miniprot [Li 2023] (protein-guided alignments)
  - BUSCO [Manni et al. 2021] (nuclear markers)
  - K-mer depth (Meryl [Miller et al. 2008], GenomeScope2 [Ranallo-Benavidez et al. 2020])

- **Assemblers**

  - Flye v2.9.6 — long-read assembler [Kolmogorov et al. 2019]
  - Miniasm — seed contig generation [Li 2016]
  - Oatk — syncmer graph + pathfinder [Choi et al. 2025]
  - TIPPo — HiFi targeted [Xie et al. 2022]

- **Post-assembly**
  - Bandage — graph visualization [Wick et al. 2015]
  - Quast [Gurevich et al. 2013], CheckM [Parks et al. 2015], GeSeq [Tillich et al. 2017]

::: notes
These complementary tools ensure we start with clean organelle-origin reads, choose the right assembler, and finish with quality validation and annotation.  
:::

---

## 6. Future Directions

- **Improved ONT accuracy**

  - Duplex ONT + new basecalling [Oxford Nanopore Tech. 2023] -> Oatk viable for ONT

- **Automated read partitioning**

  - EM/k-mer based depth clustering (e.g., GenomeScope2 extensions)

- **Graph-based approaches**

  - Repeat detection, linear vs circular path classification

- **Annotation integration**

  - Use mt/pt gene annotations to guide graph pathfinding

- **Metagenome-like strategies**
  - Treat plant mtDNA assemblies as "pangenomes" of isoforms [Wu & Sloan 2019]

::: notes
With ONT improving in quality, Oatk could unify HiFi and ONT workflows.  
Future pipelines may integrate annotation directly into assembly graphs.  
Plant mitochondria may be best understood as dynamic isoform collections.  
:::

---

## 7. Recommendations for Practice

- **Plastid DNA**

  - HiFi -> Oatk (preferred)
  - ONT -> protein-guided + Flye

- **Mitochondrial DNA**

  - HiFi -> Oatk (preferred)
  - ONT -> filtering + Miniasm seeds + Flye

- Always benchmark with **bolap** to compare strategies

::: notes
In teaching or practice:

- For plastids, HiFi+Oatk is best.
- For mitochondria, HiFi+Oatk again, but ONT workflows still work with filtering and seeding.  
  :::

---

## 8. Workflow Overview

![](figs/polap_organelle_assembly.svg){fig-alt="Flowchart of polap organelle assembly (ONT vs HiFi, ptDNA vs mtDNA)" width=90%}

::: notes
This figure provides a quick overview of the four cases.  
Note that ONT requires extra filtering and seeding, while HiFi often goes directly into Oatk.  
:::

---

# References

- Wick RR et al. (2019). _Completing bacterial genome assemblies with multiplex MinION sequencing._ Microbial Genomics.
- Wenger AM et al. (2019). _Accurate circular consensus long-read sequencing improves variant detection and assembly._ Nat Biotech.
- Palmer JD (1985). _Comparative organization of chloroplast genomes._ Ann Rev Genet.
- Daniell H et al. (2016). _Chloroplast genomes: diversity and evolution._ F1000Res.
- Sloan DB (2013). _One ring to rule them all? Genome sequencing provides new insights into plant mitochondrial DNA._ Genome Biol.
- Hazkani-Covo E et al. (2010). _Mitochondrial DNA insertion into nuclear genome._ Trends Genet.
- Kolmogorov M et al. (2019). _Assembly of long, error-prone reads using repeat graphs._ Nat Biotech.
- Li H (2016). _Minimap and miniasm: fast mapping and de novo assembly for noisy long sequences._ Bioinformatics.
- Li H (2018). _Minimap2: pairwise alignment for nucleotide sequences._ Bioinformatics.
- Li H (2023). _Miniprot: protein-to-genome alignment with splicing._ Bioinformatics.
- Manni M et al. (2021). _BUSCO update: streamlined workflows and new phylogenetic datasets._ Mol Biol Evol.
- Miller JR et al. (2008). _Aggressive assembly of pyrosequencing reads with mates._ Bioinformatics. (Meryl)
- Ranallo-Benavidez TR et al. (2020). _GenomeScope 2.0 and Smudgeplot for reference-free profiling._ Nat Commun.
- Xie Y et al. (2022). _TIPPo: targeted assembly of plastid and mitochondrial genomes._ Bioinformatics.
- Wick RR et al. (2015). _Bandage: interactive visualization of de novo genome assemblies._ Bioinformatics.
- Gurevich A et al. (2013). _QUAST: quality assessment tool for genome assemblies._ Bioinformatics.
- Parks DH et al. (2015). _CheckM: assessing the quality of microbial genomes._ Genome Res.
- Tillich M et al. (2017). _GeSeq: versatile and accurate annotation of organelle genomes._ Nucleic Acids Res.
- Wu Z, Sloan DB (2019). _Recombination and intraspecific polymorphism in plant mitochondrial genomes._ Mol Biol Evol.
- Oxford Nanopore Technologies (2023). Duplex sequencing documentation.
- Choi SC et al. (2025). _Oatk: syncmer-based assembly of plant organelle genomes._ (preprint/published)

---
