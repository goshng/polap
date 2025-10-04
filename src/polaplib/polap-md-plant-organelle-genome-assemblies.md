Here’s a **short review-style manuscript section** based on your expanded review. I structured it in a style suitable for publication, with inline citations (you can put them in a `.bib` file).

---

## Plant Organelle Genome Assembly: Current Strategies and Future Directions

Long-read sequencing technologies have transformed the assembly of plant organelle genomes, yet the strategies differ depending on the platform and the organelle target. Two major sequencing technologies are in use: Oxford Nanopore Technologies (ONT) and PacBio HiFi circular consensus sequencing. ONT provides very long reads, advantageous for spanning large repeats, but has historically suffered from lower base accuracy, which complicates assembly of highly repetitive and rearranged plant mitochondrial DNA (mtDNA) [Wick et al. 2019]. By contrast, PacBio HiFi reads are shorter but highly accurate (~99.9%), providing a strong foundation for resolving complex repeats and delivering high-quality assemblies [Wenger et al. 2019].

Plastid DNA (ptDNA) is relatively conserved across plants, with a typical size of ~150 kb and a quadripartite architecture consisting of large single-copy (LSC), small single-copy (SSC), and inverted repeat (IR) regions [Palmer 1985; Daniell et al. 2016]. This structural conservation generally makes ptDNA easier to assemble once plastid-origin reads are identified. Plant mtDNA, however, presents greater challenges. Genome sizes range from several hundred kilobases to several megabases, and they exhibit complex repeat-mediated recombination and isoform dynamics [Sloan 2013; Wu & Sloan 2019]. In addition, nuclear mitochondrial DNA (NUMTs) and plastid-derived sequences (NUPTs) introduce contamination that must be filtered prior to assembly [Hazkani-Covo et al. 2010].

In the **polap** pipeline, current strategies are tailored to both sequencing technology and organelle target. For ONT ptDNA assemblies, plastid-origin reads are selected using protein-guided alignments (e.g., Minimap2 or Miniprot [Li 2018; Li 2023]) and assembled with Flye v2.9.6 [Kolmogorov et al. 2019]. ONT mtDNA assemblies require additional preprocessing: plastid- and nuclear-origin reads are filtered out (often using BUSCO [Manni et al. 2021] or k-mer depth profiles [Ranallo-Benavidez et al. 2020]), seed contigs are generated with Miniasm [Li 2016], and Flye is used for final assembly. For HiFi data, ptDNA can be assembled either with the ONT-style protein-guided workflow or directly with the organelle-specific Oatk assembler [Choi et al. 2025], while mtDNA assembly with HiFi is most robustly handled by Oatk, which integrates syncmer-based graphs and the Pathfinder algorithm.

To support reproducibility and comparison, the companion benchmarking framework **bolap** provides standardized assembly runs (e.g., `bolap run ptgaul`, `mtgaul`, `oatk`, `tippo`). This enables evaluation of ONT versus HiFi performance, cross-assembler comparisons (Flye, Miniasm, Oatk, TIPPo [Xie et al. 2022]), and downstream quality assessment using tools such as QUAST [Gurevich et al. 2013], CheckM [Parks et al. 2015], and GeSeq [Tillich et al. 2017]. Visualization of assembly graphs with Bandage [Wick et al. 2015] remains essential for confirming plastid isomers or mitochondrial isoforms.

Looking forward, several trends are likely to shape the field. First, the improving accuracy of ONT duplex sequencing and advanced basecalling suggests that Oatk, originally designed for HiFi, may soon become the default assembler for ONT organelle genomes as well [Oxford Nanopore Tech. 2023]. Second, automated read partitioning based on depth distributions or expectation–maximization clustering could replace protein-guided approaches for organelle read selection. Third, deeper integration of annotation into assembly graphs—for example, guiding pathfinding decisions with known organelle gene loci—could improve the biological relevance of assemblies. Finally, conceptualizing plant mtDNA not as a static genome but as a population of isoforms may be a more accurate framework for analysis [Wu & Sloan 2019].

In practice, the most reliable strategy for plastid genomes is to use HiFi data with Oatk; ONT remains viable with protein-guided selection and Flye. For mitochondria, HiFi with Oatk is strongly recommended, while ONT requires careful filtering and seed contig generation prior to Flye assembly. In all cases, benchmarking with **bolap** ensures that results are reproducible and comparable across datasets and methods.

---

### References

(Would be provided as a `.bib` file for Quarto/Pandoc. Inline placeholders here are consistent with standard citations.)

---

Do you want me to generate a **ready-to-use `.bib` file** for these citations so you can immediately compile this section into Quarto/Pandoc?

Below is a **short review‑style manuscript section** you can drop into a paper or white paper. I expanded the scope beyond our earlier notes, but kept it grounded in what _polap_ does today and what _bolap_ benchmarks. Citations are embedded and point to primary literature, tool papers, and official resources.

---

## Plant organelle genome assembly: strategies, pitfalls, and a practical path forward

Long‑read sequencing has reshaped how plastid (ptDNA) and mitochondrial (mtDNA) genomes are assembled in plants. Oxford Nanopore Technologies (ONT) offers ultra‑long reads that span complex repeats but historically exhibited lower read accuracy; PacBio HiFi (CCS) yields shorter, highly accurate reads that simplify polishing and repeat resolution. These technology traits drive different, organelle‑specific strategies and error‑handling choices in practice. Recent ONT chemistry and basecalling advances (R10.4.1/Kit14, duplex basecalling via Dorado) further narrow the accuracy gap and are beginning to change default choices in organelle pipelines.([PubMed][1])

### Distinct assembly targets: plastids versus mitochondria

Plastid genomes are typically ~120–160 kb and retain a conserved quadripartite architecture (LSC–IR–SSC–IR). Their modest size, conservation, and high cellular copy number make them comparatively tractable if plastid‑origin reads are reliably selected. By contrast, plant mitogenomes span hundreds of kilobases to multiple megabases, often exist as multipartite/isoform populations, and accumulate NUMTs/NUPTs that confound read partitioning and graph resolution—conditions that amplify the benefits of accurate long reads and careful decontamination.([PLOS][2])

### What we do in **polap** today (four-case workflow)

- **ptDNA from ONT:** select plastid‑origin reads using **protein‑to‑genome alignment** (e.g., _miniprot_) and assemble with **Flye v2.9.6**. This leverages ONT read length to bridge repeats while insulating the graph from nuclear/mtDNA reads.([Oxford Academic][3])
- **mtDNA from ONT:** aggressively filter ptDNA/nuclear reads (e.g., protein markers; BUSCO to identify nuclear content), bootstrap **seed contigs with Miniasm**, then complete with **Flye**. This ONT‑specific staging improves graph tractability before the final repeat‑graph assembly.([Oxford Academic][4])
- **ptDNA from HiFi:** either mirror the protein‑guided selection + Flye route or assemble directly with **Oatk**, which couples syncmer‑based assembly with organelle gene–aware pathfinding.([BioMed Central][5])
- **mtDNA from HiFi:** **Oatk** is generally preferred; HiFi accuracy plus Oatk’s organelle‑aware graph traversal tends to reduce tangles and mis‑joins in repeat‑rich plant mitochondria.([BioMed Central][5])

The companion framework **bolap** standardizes benchmarking across these routes (e.g., _polap_ vs. **Oatk** vs. other organelle assemblers) and keeps comparisons reproducible on shared datasets and metrics (assembly continuity, gene content, graph structure). _(Internal framework; methods align with community metrics below.)_

### Read selection & decontamination: why it matters and how to do it

Misassigned reads propagate into graph artifacts. In practice, **protein‑guided alignment** against curated organelle gene sets is powerful for selecting organelle reads from WGS, with **miniprot** providing fast, accurate protein→genome mappings; **minimap2** remains the workhorse for DNA→DNA alignments. Complementary **k‑mer profiling** can estimate copy‑number/coverage enrichments and heterozygosity (useful for catching nuclear contamination and polyploid signals), while **BUSCO** offers an ortholog‑based lens on residual nuclear content.([Oxford Academic][3]) ([Nature][6]) ([Oxford Academic][4])

### Assemblers and organelle‑aware post‑processing

Repeat‑graph assemblers such as **Flye** are robust across ONT/HiFi and remain a practical default, especially when paired with ONT. **Miniasm** is often used to draft “seed” contigs rapidly from noisy data (no consensus), after which Flye can refine the repeat graph. For HiFi, **Oatk** (syncmer‑based assembly + _hmmannot_ + **Pathfinder**) is expressly designed for plant organelles and has now been peer‑reviewed, reporting strong performance across large panels of species. Newly published **TIPPo** targets HiFi organelle assembly with a modern classifier‑guided read selection, and **PMAT** focuses on **low‑coverage HiFi** mitogenomes; both broaden the toolset beyond historical short‑read specialists (e.g., **GetOrganelle** for Illumina).([PubMed][7]) ([BioMed Central][5]) ([PubMed][8]) ([PubMed][9]) ([BioMed Central][10])

> **Context for short‑read toolkits:** _GetOrganelle_ has been the go‑to for Illumina “genome skims,” producing highly reliable plastomes; for HiFi/ONT, its role is typically complementary (read recruitment/validation) rather than the final assembler.([BioMed Central][10])

### Validation: from graphs to genes (and claims you can defend)

After graph construction, **graph‑level inspection** (e.g., **Bandage**) helps confirm repeat resolution, circular candidates, and alternative isoforms. **Sequence‑level QC** with **QUAST** remains standard; **annotation‑aware validation** (**GeSeq**) adds organelle gene content checks and HMM‑based feature calls; **OGDRAW/Chloroplot** produce publication‑quality maps that also reveal structural inconsistencies (e.g., IR size drift, SSC inversions). For mitochondria, interpreting **isoform mixtures** (not just a “master circle”) is increasingly expected—reporting stoichiometry or multiple supported paths when appropriate.([Oxford Academic][11]) ([PubMed][12]) ([Oxford Academic][13]) ([PubMed][14]) ([PMC][15])

### Where the field is heading (and why it affects your defaults)

**ONT accuracy is rising** (duplex Q30 claims; improved R10.4.x chemistry). As ONT simplex/duplex distributions and dorado‑SUP models continue to improve, **HiFi‑centric organelle assemblers (e.g., Oatk)** are likely to become **platform‑agnostic defaults**, even for ONT reads, especially on challenging plant mitogenomes. Meanwhile, **automated read partitioning** (deep models; HMM‑guided filters) is reducing reliance on manual bait sets, and recent tools like **TIPPo** and **PMAT** illustrate how classifier‑guided selection or ultra‑low HiFi coverage can still yield complete organelles. For mitochondria specifically, **reporting isoform networks** rather than a single “master circle” better reflects biological reality and should be standard in results and figures.([PMC][16]) ([BioMed Central][5]) ([PMC][15])

### Practical recommendations (as of Sept 2025)

- **Plastids (ptDNA).**
  **HiFi:** Prefer **Oatk** end‑to‑end; alternatively, protein‑guided selection → Flye remains robust.
  **ONT:** Protein‑guided read selection (**miniprot**) → **Flye 2.9.6**; expect better outcomes with duplex‑rich datasets in newer chemistries.([Oxford Academic][3])
- **Mitochondria (mtDNA).**
  **HiFi:** Prefer **Oatk**; evaluate multiple graph paths and report isoforms.
  **ONT:** Filter pt/nuclear reads (protein markers; BUSCO), draft **Miniasm** seeds, then assemble with **Flye**; use **Bandage** to vet repeats/isoforms and annotate with **GeSeq** before claiming closure.([Oxford Academic][4])
- **Benchmarking & reporting.**
  Use a standardized harness (e.g., _bolap_ in our ecosystem) and publish: (i) graph snapshots (Bandage PNG/SVG), (ii) QUAST summaries, (iii) gene content via GeSeq, and (iv) if applicable, **multiple** supported mtDNA configurations with coverage support.([Oxford Academic][11])

---

### Selected sources (core methods & tools)

**Assemblers & aligners:** Flye long‑read repeat‑graph assembler; Miniasm OLC seeding; Minimap2 DNA alignments; Miniprot protein‑to‑genome mapping.([PubMed][7])

**Organelle‑specific toolkits:** Oatk (syncmer + HMM + pathfinding) for HiFi/ONT; TIPPo (HiFi, classifier‑guided); PMAT (low‑coverage HiFi mitogenomes); GetOrganelle (Illumina skims).([BioMed Central][5])

**QC, annotation, and visualization:** Bandage for assembly graphs; QUAST for assembly metrics; GeSeq (HMM‑assisted annotation); OGDRAW & Chloroplot for figure‑ready maps.([Oxford Academic][11])

**Biological context:** Plastid quadripartite structure; plant mitogenome size/isoform complexity; NUMTs/NUPTs as confounders.([PLOS][2])

**Platform developments:** HiFi accuracy (CCS/HiFi); ONT R10.4.1/Kit14, duplex basecalling/Dorado.([PubMed][1]) ([PMC][17])

---

#### Notes for readers who will reproduce these results

- **Protein‑guided selection.** We prefer **miniprot** for read triage against curated pt/mt protein sets because it’s genuinely protein→genome (unlike DNA‑only mappers). Minimap2 remains our default for DNA‑level mapping and coverage checks.([Oxford Academic][3])
- **“Closure” versus “representation.”** For plant mitochondria, we report **graphs and supported isoforms** (with coverage) rather than forcing a single circular representation. This aligns with recent critiques of the “master circle” model.([PMC][15])

---

If you’d like, I can now tailor this section into journal‑specific formatting (e.g., word count, structured headings, and reference style), or expand the **Methods** box with exact command‑line invocations for the four **polap** cases (including _bolap_ benchmarking stanzas and figure generation).

[1]: https://pubmed.ncbi.nlm.nih.gov/31406327/?utm_source=chatgpt.com "Accurate circular consensus long-read sequencing ..."
[2]: https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0263253&utm_source=chatgpt.com "The complete chloroplast genome sequencing analysis ..."
[3]: https://academic.oup.com/bioinformatics/article/39/1/btad014/6989621?utm_source=chatgpt.com "Protein-to-genome alignment with miniprot - Oxford Academic"
[4]: https://academic.oup.com/mbe/article/38/10/4647/6329644?utm_source=chatgpt.com "BUSCO Update: Novel and Streamlined Workflows along with ..."
[5]: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-025-03676-6?utm_source=chatgpt.com "Oatk: a de novo assembly tool for complex plant organelle ..."
[6]: https://www.nature.com/articles/s41467-020-14998-3?utm_source=chatgpt.com "GenomeScope 2.0 and Smudgeplot for reference-free ..."
[7]: https://pubmed.ncbi.nlm.nih.gov/30936562/?utm_source=chatgpt.com "Assembly of long, error-prone reads using repeat graphs"
[8]: https://pubmed.ncbi.nlm.nih.gov/39800935/ "TIPPo: A User-Friendly Tool for De Novo Assembly of Organellar Genomes with High-Fidelity Data - PubMed"
[9]: https://pubmed.ncbi.nlm.nih.gov/38469379/?utm_source=chatgpt.com "PMAT: an efficient plant mitogenome assembly toolkit ..."
[10]: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02154-5?utm_source=chatgpt.com "GetOrganelle: a fast and versatile toolkit for accurate de novo ..."
[11]: https://academic.oup.com/bioinformatics/article/31/20/3350/196114?utm_source=chatgpt.com "Bandage: interactive visualization of de novo genome ..."
[12]: https://pubmed.ncbi.nlm.nih.gov/23422339/?utm_source=chatgpt.com "QUAST: quality assessment tool for genome assemblies"
[13]: https://academic.oup.com/nar/article/45/W1/W6/3806659?utm_source=chatgpt.com "GeSeq – versatile and accurate annotation of organelle ..."
[14]: https://pubmed.ncbi.nlm.nih.gov/30949694/?utm_source=chatgpt.com "OrganellarGenomeDRAW (OGDRAW) version 1.3.1"
[15]: https://pmc.ncbi.nlm.nih.gov/articles/PMC6742443/?utm_source=chatgpt.com "The alternative reality of plant mitochondrial DNA: One ring ..."
[16]: https://pmc.ncbi.nlm.nih.gov/articles/PMC11077942/?utm_source=chatgpt.com "Closing the gap: Oxford Nanopore Technologies R10 ..."
[17]: https://pmc.ncbi.nlm.nih.gov/articles/PMC9973852/?utm_source=chatgpt.com "Comparison of R9.4.1/Kit10 and R10/Kit12 Oxford ..."
