# Polap: Plant Mitochondrial DNA Assembly Pipeline

Polap is a specialized pipeline designed to assemble plant mitochondrial DNA (mtDNA) using the **Flye** long-read assembler, supplemented by organellar gene annotation to guide contig selection. Inspired by **ptGAUL** (a chloroplast assembly pipeline), Polap integrates Flye's capabilities with gene annotation to identify and assemble mtDNA contigs from whole-genome data.

NOTE: use the version 0.3.7.3 from the conda package.

```bash
conda create -y --name polap polap=0.3.7.3
```

<!-- TOC START -->

## Table of Contents

- [Quick Start](#quick-start)
  - [Requirements](#requirements)
  - [1. Open a new terminal](#1-open-a-new-terminal)
  - [2. Install Miniconda (5 mins)](#2-install-miniconda-5-mins)
  - [3. Setup the conda channels](#3-setup-the-conda-channels)
  - [4. Install the Bioconda Polap package (5 mins)](#4-install-the-bioconda-polap-package-5-mins)
  - [5. Activate the polap conda environment and setup polap-fmlrc environment (2 min)](#5-activate-the-polap-conda-environment-and-setup-polap-fmlrc-environment-2-min)
  - [6. Polap run with a test dataset (10 mins)](#6-polap-run-with-a-test-dataset-10-mins)
  - [7. Polish the test dataset (10 mins)](#7-polish-the-test-dataset-10-mins)
  - [8. Test BioProject SRA data (2 hours)](#8-test-bioproject-sra-data-2-hours)
  - [9. Use your data to assemble a plant mtDNA sequence](#9-use-your-data-to-assemble-a-plant-mtdna-sequence)
- [Description](#description)
- [Key Features](#key-features)
- [How It Works](#how-it-works)
- [Usage with main Polap menus](#usage-with-main-polap-menus)
  - [Help Message](#help-message)
  - [Input Data](#input-data)
  - [Whole-Genome Assembly](#whole-genome-assembly)
  - [Plant Organelle Annotation](#plant-organelle-annotation)
  - [Organelle Contig Selection](#organelle-contig-selection)
  - [Organelle-Genome Assembly](#organelle-genome-assembly)
  - [Extracting a mitochondrial DNA sequence using Bandage](#extracting-a-mitochondrial-dna-sequence-using-bandage)
  - [Polishing](#polishing)
  - [Sample Datasets](#sample-datasets)
  - [MT Contig Name File](#mt-contig-name-file)
  - [Genome assembly numbers](#genome-assembly-numbers)
- [Using BioProject data](#using-bioproject-data)
- [Detailed instructions of installation](#detailed-instructions-of-installation)
  - [Using Bioconda packages](#using-bioconda-packages)
  - [Short-read polishing Requirement: FMLRC](#short-read-polishing-requirement-fmlrc)
  - [Using github source](#using-github-source)
  - [Updating the conda Polap package with the github source](#updating-the-conda-polap-package-with-the-github-source)
- [Uninstalling Polap](#uninstalling-polap)
  - [Step 1: Uninstall Miniconda](#step-1-uninstall-miniconda)
  - [Step 2: Remove Environment Variables](#step-2-remove-environment-variables)
  - [Step 3: Remove Conda Configuration Files](#step-3-remove-conda-configuration-files)
  - [Step 4: Verify Removal](#step-4-verify-removal)
- [More](#more)
  - [Options](#options)
    - [Menu options](#menu-options)
    - [General options](#general-options)
    - [Detailed description of some of the options](#detailed-description-of-some-of-the-options)
  - [Authors](#authors)

<!-- TOC END -->

## Quick Start

#### Requirements

- **Operating System**: Linux (not compatible with macOS or Windows)
- **Dependencies**: Requires [Bash](https://www.gnu.org/software/bash/) (>= 5.0) and **[Miniconda](https://docs.anaconda.com/miniconda/miniconda-install/)**
- **Installation**: Can be installed via [Bioconda](https://bioconda.github.io/)'s [Polap conda package](https://anaconda.org/bioconda/polap) or manually from the Github source available at [Polap](https://github.com/goshng/polap)

#### 1. Open a new terminal

Open a new terminal in a Linux computer, such as one with Ubuntu. Follow the next steps by copying and pasting the scripts to your favorite terminal.

#### 2. Install Miniconda (5 mins)

Download and install **[Miniconda](https://docs.anaconda.com/miniconda/miniconda-install/)** using the [instructions](https://docs.anaconda.com/miniconda/#quick-command-line-install)

```bash
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm ~/miniconda3/miniconda.sh
```

After installing, close and reopen your terminal application.

#### 3. Setup the conda channels

If you did not close and reopen a new terminal, please do so. Then, execute the followings to setup the conda channels for `polap`.

```bash
source ~/miniconda3/bin/activate
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

#### 4. Install the Bioconda Polap package (5 mins)

Setup `polap` and `polap-fmlrc` conda environments using [Polap](https://anaconda.org/bioconda/polap) conda package.

```bash
conda create -y --name polap polap=0.3.7.3
```

#### 5. Activate the polap conda environment and setup polap-fmlrc environment (2 min)

```bash
conda activate polap
base_dir=$(dirname "$(command -v polap)") && conda env create -f $base_dir/polap-conda-environment-fmlrc.yaml
```

#### 6. Polap run with a test dataset (10 mins)

```bash
wget -q https://github.com/goshng/polap/archive/refs/tags/0.3.7.3.zip
unzip -o -q 0.3.7.3.zip
cd polap-0.3.7.3/test
polap assemble --test
```

If you see this output, then you are ready to use the long-read assembly pipeline.

```bash
output: the assembly graph: 2-oga.gfa
```

After completing this step, you can move on to step 7 of [testing the short-read polishing](#7-polish-the-test-dataset) process.

#### 7. Polish the test dataset (10 mins)

Then, polish the provided sequence file, `mt.0.fasta`, which is a test file for the short-read polishing.

```bash
wget -q https://github.com/goshng/polap/archive/refs/tags/0.3.7.3.zip
unzip -o -q 0.3.7.3.zip
cd polap-0.3.7.3/test
polap prepare-polishing
polap polish
diff mt.0.fa mt.1.fa
```

If you see no output from the command, `diff mt.0.fa mt.0.fa`, then you are ready to use the polishing procedure as well.

#### 8. Test BioProject SRA data (2 hours)

Because you have done the test data run, you do not need this step and proceed to the next step where you could try Polap with your own data set.
Note that using Polap's `assemble` menu is experimential or not tested yet! You might have a very simple plant mitochondrial DNA genome structure like this case: mtDNA assembly of _Carex pseudochinensis_ at a linux compuer with Intel(R) Xeon(R) CPU E5-2690 v4 @ 2.60GHz and 256 GB memory.

```bash
mkdir Carex_pseudochinensis
cd Carex_pseudochinensis
polap x-ncbi-fetch-sra --sra SRR30757341
polap x-ncbi-fetch-sra --sra SRR30757340
cp -s SRR30757341.fastq l.fq
cp -s SRR30757340_1.fastq s1.fq
cp -s SRR30757340_2.fastq s2.fq
polap assemble --polap-reads
```

#### 9. Use your data to assemble a plant mtDNA sequence

Now, if you have your Illumina paired-end short- and Oxford Nanopore long-read data files, put those in a folder and rename them as follows. Assume that you have three fastq files:

```bash
SRR30757340_1.fastq <- the first short-read FASTQ file
SRR30757340_2.fastq <- the second short-read FASTQ file
SRR30757341.fastq   <- the long-read FASTQ file
```

Rename the three files as follows:

```bash
mv SRR30757340_1.fastq s1.fq
mv SRR30757340_2.fastq s2.fq
mv SRR30757341.fastq l.fq
```

Then, execute the Polap command as you do in the provided test data:

```bash
polap assemble1
```

Once it completes a whole-genome assembly, then try to find mtDNA seed contigs by displaying a plant organelle gene annotation table for the genes that were roughly annotated by Polap using the NCBI organelle protein sequences that were downloaded May, 2023 using the following command:

```bash
polap annotate view table
```

| Contig   |  Length | Depth | Copy |  MT |  PT | Edge |
| :------- | ------: | ----: | ---: | --: | --: | ---: |
| edge_47  |  286811 |    69 |   10 |  25 |   1 |   47 |
| edge_729 |   50873 |    64 |    9 |   9 |   1 |  729 |
| edge_46  | 1532068 |    14 |    2 |   3 |   0 |   46 |
| edge_75  | 4781975 |    12 |    2 |   3 |   1 |   75 |
| edge_55  | 4374598 |    13 |    2 |   1 |   0 |   55 |
| edge_718 |   27942 |    78 |   11 |   1 |   0 |  718 |
| edge_732 |   17165 |    74 |   11 |   1 |   0 |  732 |

Second, open the whole-genome assembly file named `o/0/30-contigger/graph_final.gfa` using **[Bandage](https://rrwick.github.io/Bandage/)** and search the genome assembly for those that were annotated with more mtDNA than ptDNA genes. In the plant organelle gene annotation table, each row represents a contig sequence annotated with at least one plant organelle gene. Although Polap does not automatically select mtDNA-related contigs, it could help you guess which ones might originate from plant mtDNA. The table has several columns for a row:

1. Contig: the name of a contig sequence
2. Length: the length of the contig sequence
3. Depth: the read coverage of the contig sequence
4. Copy: the depth value of the contig divided by the median of all the depth values
5. MT: the number of rough mitochondrial genes annotated
6. PT: the number of rough plastid genes annotated
7. Edge: the edge number of the contig sequence (the same number of the contig name in the row)

You would choose contig names with the following properties:

- MT > PT
- Depth or Copy values are roughly in between possible depth values from nuclear and those plastid.

Now, go back to the Bandage genome assembly graph to locate one or some of your choice of candidate seed mtDNA contigs.
Watch [YouTube](https://youtu.be/2dViWNEqueU) to see how you could copy contig names using the Bandage software.
Here, you would select contigs with a depth range that you think mtDNA contigs should have.

```
edge_264, edge_266, edge_47, edge_718, edge_729, edge_732
```

In the example of [YouTube](https://youtu.be/2dViWNEqueU), all depth values of the selected contigs are inclusivly in between 61x and 78x.
Come back to the terminal where you have executed `polap assemble1`, and execute the following:

```bash
polap seeds bandage
```

Then, paste or type in your copy of contig names to the question:

```bash
Enter edges using Bandage (e.g., edge_265, edge_520, edge_425):
```

This would create `o/0/mt.contig.name-1` with following content.

```
edge_47
edge_729
edge_718
edge_732
edge_264
edge_266
```

You could see seed contigs that were added to the previous annotation table:

```bash
polap seeds view
```

Under the column `Seed` of the annotation table, `A` denotes that the contig was part of the initial annotation, `G` that the contig was added to the seed set via the genome assembly graph, and `X` the contig was not part of the seed set though it was in the initial annotation.

| Contig    | Length | Depth | Copy |  MT |  PT | Seed |
| :-------- | -----: | ----: | ---: | --: | --: | ---: |
| edge_28   | 385224 |     7 |    2 |  11 |   2 |    A |
| edge_30   | 129617 |     8 |    3 |  11 |   1 |    A |
| edge_33   | 100315 |     7 |    2 |   7 |   0 |    A |
| edge_2864 |  11463 |   284 |   95 |   6 |   4 |    X |
| edge_31   |  60195 |     8 |    3 |   4 |   1 |    A |
| edge_2803 |  79961 |     8 |    3 |   4 |   0 |    A |
| edge_29   | 129920 |     6 |    2 |   3 |   2 |    A |
| edge_2    | 384522 |     3 |    1 |   2 |   0 |    X |
| edge_1282 |   3250 |   517 |  172 |   2 |   1 |    X |
| edge_2817 |  47958 |     8 |    3 |   2 |   1 |    A |
| edge_1888 | 231323 |     3 |    1 |   1 |   0 |    X |
| edge_1510 | 123019 |     3 |    1 |   0 |   0 |    G |
| edge_2802 |  20149 |    13 |    4 |   0 |   0 |    G |
| edge_32   |   8681 |    18 |    6 |   0 |   0 |    G |

You have prepared a list of seed contigs for an organelle-genome assembly. Execute the following:

```bash
polap assemble2
```

If you are very lucky, you might already have a good candidate for your mtDNA from the Flye whole-genome assembly. You might be able to locate a circular mtDNA assembly graph or a graph that you could traverse to generate a circular mtDNA sequence. In that case, you may not need the Polap procedure, and proceed to polish the mtDNA from the whole-genome assembly using your input short-read data. If you are not that lucky to have a good mtDNA candidate in the whole-genome assembly but are lucky enough to have contigs clustered together in the genome assembly graph like the one example of the [YouTube](https://youtu.be/2dViWNEqueU).

After the Polap's organelle-genome assembly step like above, you have two paths to follow. One is to stop there and extract an mtDNA sequence by following a path in the organelle-genome assembly graph using the Bandage software. The other is to consider the organelle-genome assembly as a whole-genome assembly and repeat the same procedure of seed contig selection and one more organelle-genome assembly step. We will consider the other scenario later and assume that you are lucky to have an mtDNA sequence from the organelle-genome assembly.

Flye creates your long-read polished assembly graph as a file named `o/1/assembly_graph.gfa`.
Extract the mitochondrial genome sequence by opening `o/1/assembly_graph.gfa` in [Bandage](https://rrwick.github.io/Bandage/). Save this sequence as `mt.0.fasta` for the long-read mitochondrial genome assembly. You should be able to extract a draft organelle genome sequence from the assembly by watching [YouTube](https://youtu.be/UF-0UIc2ZDY) that I would use in an undergraduate lecture in Korean. The audio does not matter, so just watch the YouTube if you have never used the Bandage software before.

Now, you have prepared a mtDNA draft sequence in the FASTA file named `mt.0.fasta`. From here to the finished mtDNA genome, the credit is due to **ptGAUL**, which uses **FMLRC** for the short-read polishing step. Polap has a handy menu called `polap prepare-polishing` and `polap polish` for the procedure. The former command converts your short-read data to a file in a format **FMLRC** could use to polish your draft mtDNA genome sequence in the latter command. You could have executed the first command `polap prepare-polishing` before the Flye whole-genome assembly. Execute the following:

```bash
polap prepare-polishing
polap polish
```

Your final mtDNA genome assembly is in the file named `mt.1.fa`, that you use [GeSeq](https://chlorobox.mpimp-golm.mpg.de/geseq.html) to actually annotate.
For more detailed installation and usage instructions, please read along with this **README**.

## Description

**Polap** is a specialized meta-assembly pipeline that enhances the functionality of the [Flye](https://github.com/fenderglass/Flye) assembler. Serving as a Flye helper pipeline, Polap incorporates organellar gene annotation, aiding in the educated guess of selection of organelle-origin contigs from Flye's whole-genome assembly. The name "Polap" is inspired by the [ptGAUL](https://github.com/Bean061/ptgaul) chloroplast genome assembly pipeline created by Wenbin Zhou, as Polap similarly focuses on organellar genomes, specifically plant mitochondrial DNA (mtDNA).

The pipeline combines the **Flye** assembler with [FMLRC](https://github.com/holtjma/fmlrc) (a sequence polishing tool by Matt Holt) and applies sequencing depth analysis to identify contigs belonging to the mitochondrial genome. Leveraging both Flye-generated assemblies and organellar gene annotations, Polap helps users identify candidate mtDNA contigs based on sequencing depth and gene presence, enabling targeted selection and assembly of mtDNA.

Polap originated as an extension of ptGAUL, with a specific focus on assembling plant mtDNA from long-read sequencing data. The workflow begins by conducting a Flye whole-genome assembly, generating contigs without a reference genome. These contigs are then analyzed with [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) to confirm the presence of mtDNA or chloroplast (ptDNA) genes. Using depth and copy number metrics, Polap allows users to identify likely mtDNA contigs, collects them, and performs a final targeted mtDNA assembly. The pipeline concludes with polishing using FMLRC, resulting in an accurate, high-quality mtDNA assembly.

## Key Features

- **Rough organellar gene annotation**: Enables Polap to select contigs from mitochondrial and chloroplast origins, focusing specifically on the target mtDNA.
- **Organellar seed contig selection**: (Not tested yet!) Polap adds custom logic to filter and identify organelle-specific contigs.

## How It Works

1. **Initial whole-genome assembly**: Flye performs a whole-genome assembly on the input long reads. We need paired-end short-read data to estimate the genome size using **JellyFish**. The current version of Flye does not require a genome size for the genome assembly, but we have not test it yet. Because Polap still needs short-read data for a final polishing using FMLRC, it uses short-read data to estimate the genome size that might help Flye perform a whole-genome assembly with less memory and computing time.
2. **Organellar seed contig selection**: Using BLAST, Polap identifies mtDNA-related contigs based on the presence of organelle genes. Users need to manually select mitochondrial-origin contigs using organelle gene annotation tables that are made by **Flye** and modified by **Polap**.
3. **Organelle-gename assembly**: Selected mtDNA seed contigs and the input long-read sequencing data are aligned using **[Minimap2](https://github.com/lh3/minimap2)** and those reads that are well mapped on the input long reads are fed into **Flye** to assembly the organellar genome. We borrow the approach of **ptGAUL** in this long-read selection.
4. **Final Polishing**: The selected mtDNA contigs undergo polishing with FMLRC for accurate genome assembly. Short-read sequence polishing uses **FMLRC** for high-quality, error-corrected assemblies of the target mitochondrial genome as suggested by **ptGAUL**.

## Usage with main Polap menus

After successfully running the test dataset with Polap, you have two options for executing Polap, depending on your installation method.

1. **If you installed Polap via [Bioconda](https://anaconda.org/bioconda/polap)**:  
   Simply type `polap` in the command line, provided you have the Polap Conda environment activated:

   ```bash
   (polap) $ polap
   ```

   Be sure the prompt shows `(polap)`, indicating the correct environment is active.

2. **If you installed Polap from the GitHub source**:  
   Clone the Polap repository into your data folder and run the main script. Ensure the Polap Conda environment is still active.

   ```bash
   (polap) $ git clone https://github.com/goshng/polap.git
   (polap) $ cd polap
   (polap) $ src/polap.sh
   ```

For simplicity, the instructions below assume you installed Polap using Bioconda. In cases where you clone the GitHub repository, replace `polap` with `src/polap.sh` as needed. If you're experienced with the command line, this setup will feel straightforward. Otherwise, follow the steps in this README closely.

### Help Message

Polap operates primarily through subcommands, also called "menus." Key menus include `assemble1`, `annotate`, and `assemble2`, each with its own help message.

To view the general help message, run Polap without any options:

```bash
(polap) $ polap
```

To see help for a specific menu, such as `assemble1`, use:

```bash
(polap) $ polap assemble1 help
```

### Input Data

1. **Set Up Your Input Folder**:  
   Create a folder and place your input dataset within it. This dataset should consist of:

   - One Oxford Nanopore long-read FASTQ file
   - Two short-read FASTQ files

   A set of sample FASTQ files is available in the `test` folder if you've cloned the Polap repository. Currently, only configurations with two short-read files are supported.

2. **File Requirements**:  
   If you don't specify input file paths with options, Polap will expect three files named `l.fq`, `s1.fq`, and `s2.fq` in the current working directory where you execute the Polap script. Make sure these files are uncompressed. You can also specify the files using the options `-l`, `-a`, and `-b`.

### Whole-Genome Assembly

With your three FASTQ files prepared, you’re ready for the first step in Polap: running a whole-genome assembly using Flye.

1. **Menu Command Structure**:  
   In Polap, the primary command argument is referred to as the "menu," which specifies the task. You can view help for any menu by adding `help` as a second argument. To see the options for the `assemble1` menu, use:

   ```bash
   (polap) $ polap assemble1 help
   ```

2. **Run Assembly**:  
   Start the whole-genome assembly with:
   ```bash
   (polap) $ polap assemble1
   ```
   Assembly times vary based on the dataset, so please be patient.

### Plant Organelle Annotation

After completing the whole-genome assembly, you can visualize the genome assembly graph in `o/0/30-contigger/graph_final.gfa` using **[Bandage](https://rrwick.github.io/Bandage/)**. This graph shows contigs, some originating from mitochondrial or plastid DNA, alongside nuclear DNA. You could use BLAST your assembly against organelle genes and provide those BLAST results to the Bandage software to locate potential mtDNA-related contigs. Polap can facilitate such process with the command named `polap annotate`, which is part of `polap assemble1`. So, you could display the annotation table that focuses on mtDNA-related contigs.

```bash
(polap) $ polap annotate view table
```

| Contig   |  Length | Depth | Copy |  MT |  PT | Edge |
| :------- | ------: | ----: | ---: | --: | --: | ---: |
| edge_47  |  286811 |    69 |   10 |  25 |   1 |   47 |
| edge_729 |   50873 |    64 |    9 |   9 |   1 |  729 |
| edge_46  | 1532068 |    14 |    2 |   3 |   0 |   46 |
| edge_75  | 4781975 |    12 |    2 |   3 |   1 |   75 |
| edge_55  | 4374598 |    13 |    2 |   1 |   0 |   55 |
| edge_718 |   27942 |    78 |   11 |   1 |   0 |  718 |
| edge_732 |   17165 |    74 |   11 |   1 |   0 |  732 |

If you can identify mitochondrial DNA contigs manually, you can proceed with those to the organelle-genome assembly using `polap assemble2`. But, take your time to choose seed contigs in the following section.

### Organelle Contig Selection

After completing the `polap annotate` step, you'll find organelle gene annotations in text files within the `o/0/` directory:

- `contig-annotation-depth-table.txt`
- `contig-annotation-table.txt`
- `assembly_info_organelle_annotation_count-all.txt`
- `assembly_info_organelle_annotation_count.txt`

You could inspect these files to identify contigs with a higher number of mitochondrial (MT) genes compared to plastid (PT) genes. Execute the following command to display organelle gene annotation with a particular focus on mtDNA genes:

```bash
polap annotate view table
```

It displays a table that you could inspect to determine mtDNA seed contigs that you could try to use later in the Polap's organelle-genome assembly.

This seed contig selection step is subjective, so you may need to carefully review the annotation table files to choose contigs that likely originate from plant mitochondrial DNA. These annotation tables are derived from `o/0/30-contigger/graph_final.gfa` that Flye has generated. We modify `contigs_stats.txt` from Flye to create a similar file and provide information about each contig assembled from the whole-genome assembly.

You’ll use both the graph visualization of the genome assembly `o/0/30-contigger/graph_final.gfa`, and the annotation table files to select candidate contigs. Once identified, prepare a text file (`o/0/mt.contig.name-1`) with each selected edge sequence name on a new line. Note that edge sequence names should start with `edge_`, not `contig_`. An example file might look like this:

```
edge_47
edge_729
edge_718
edge_732
edge_264
edge_266
```

Your could visualize the whole-genome assembly graph using the Bandage software to select contigs. You would choose contig names with the following properties:

- MT > PT
- Depth or Copy values are roughly in between possible depth values from nuclear and those plastid.

Watch [YouTube](https://youtu.be/2dViWNEqueU) to see how you could copy contig names using the Bandage software.
Here, you would select and copy contig names with a depth range that you think mtDNA contigs should have.
Come back to the terminal where you have executed `polap assemble1`, and execute the following:

```bash
polap seeds bandage
```

Then, paste or type in your copy of contig names to the question:

```bash
Enter edges using Bandage (e.g., edge_265, edge_520, edge_425):
```

It generates a text file (`o/0/mt.contig.name-1`) with edge names that you paste. For more information on preparing a mitochondrial contig file, refer to the [MT contig name](#mt-contig-name) section below.

### Organelle-Genome Assembly

With the contig candidates selected, you are ready to proceed with the organelle-genome assembly:

```bash
(polap) $ polap assemble2
```

This command will initiate the second Flye assembly, generating an assembly graph file `o/1/assembly_graph.gfa`. You can view this organelle-genome assembly graph using Bandage.

### Extracting a mitochondrial DNA sequence using Bandage

Extract the mitochondrial genome sequence by opening `o/1/assembly_graph.gfa` in [Bandage](https://rrwick.github.io/Bandage/). Save this sequence as `mt.0.fasta` for the long-read mitochondrial genome assembly. Your long-read polished assembly graph is `o/1/assembly_graph.gfa` You could extract a draft organelle genome sequence from the assembly graph by using [Bandage](https://rrwick.github.io/Bandage/). You could watch [YouTube](https://youtu.be/UF-0UIc2ZDY) (in Korean) to see how you could extract a DNA sequence through a circular path.

### Polishing

Polap follows the [ptGAUL](https://github.com/Bean061/ptgaul) approach, which uses **[FMLRC](https://github.com/holtjma/fmlrc)** to polish long-read assemblies with short-read data. To polish your assembly, activate the FMLRC environment and run the following commands:

```bash
$ conda activate polap-fmlrc
(polap-fmlrc) $ src/polap.sh prepare-polishing
(polap-fmlrc) $ src/polap.sh polish
```

Make sure that you have your draft genome sequence file ready before calling the `polap polish` command.
You could use `-p` option to specify your draft sequence FASTA file.

The final polished mitochondrial genome sequence will be saved as `mt.1.fa`.

### Sample Datasets

#### mtDNA reference-generating method

A sample dataset demonstrating Polap's mtDNA assembly is available on Figshare: [Polap data analysis of 11 datasets](https://figshare.com/s/07305fd4e18c74080fbc).

#### ptDNA subsampling-based method

A sample dataset demonstrating Polap's ptDNA assembly is available on Figshare: [Polap data analysis of 28 datasets](https://figshare.com/s/ec1cb394870c7727a2d4).

```bash
git clone https://github.com/goshng/polap.git
conda activate
bash polap/src/polap-data-v2.sh uninstall
bash polap/src/polap-data-v2.sh install-polap
bash polap/src/polap-data-v2.sh bleeding-edge-polap
conda activate polap
polap-data-v2.sh install-fmlrc
polap-data-v2.sh install-cflye
polap-data-v2.sh delete-polap-github
polap-data-v2.sh sample-csv polap-data-v2.csv test
polap-data-v2.sh download-test-data
polap-data-v2.sh local-batch Eucalyptus_pauciflora t off
polap-data-v2.sh sample-csv polap-data-v2.csv all
polap-data-v2.sh local-batche each
```

### MT Contig Name File

To prepare an MT contig name file, list the contig names you suspect are of mitochondrial origin. Select these contigs based on three key features:

1. **Genome assembly graph**
2. **Presence of organelle genes**
3. **Contig copy numbers**

Contigs are labeled as edges, each suffixed with a number. To examine the whole-genome assembly graph, use Bandage to open `o/0/30-contigger/graph_final.gfa`. You can also view organelle annotations and copy numbers with the following command:

```bash
polap annotate view table
```

For additional guidance on identifying mtDNA contigs, a [YouTube](https://youtu.be/29OaiFCqAzI) tutorial in Korean is available. The following example shows a Flye assembly table with gene counts, where the `Copy` column indicates the copy number, `MT` represents mitochondrial genes, and `PT` represents plastid genes. The `Edge` column lists edge numbers that make up each contig.

Example table:

| Contig   |  Length | Depth | Copy |  MT |  PT | Edge |
| :------- | ------: | ----: | ---: | --: | --: | ---: |
| edge_47  |  286811 |    69 |   10 |  25 |   1 |   47 |
| edge_729 |   50873 |    64 |    9 |   9 |   1 |  729 |
| edge_46  | 1532068 |    14 |    2 |   3 |   0 |   46 |
| edge_75  | 4781975 |    12 |    2 |   3 |   1 |   75 |
| edge_55  | 4374598 |    13 |    2 |   1 |   0 |   55 |
| edge_718 |   27942 |    78 |   11 |   1 |   0 |  718 |
| edge_732 |   17165 |    74 |   11 |   1 |   0 |  732 |

After reviewing the table, edit `o/0/mt.contig.name-1` to include the names of candidate mtDNA contigs. Each line should list an edge, without any empty lines. Here’s an example:

```
edge_47
edge_729
edge_718
edge_732
edge_264
edge_266
```

### Genome assembly numbers

Polap generates all results in an output folder that you can change from the default `o` using option `-o` or `--outdir`.
In the output folder, Polap uses number 0 to represent a whole-genome assembly and positive integer numbers for organelle-genome assemblies.
You always have number 0 folder in your folder whenever you run the whole-genome assembly.
Because a genome assembly number folder is the base folder that Flye creates, it has Flye's output folders including `10-assembly`, `20-consensus`, `30-contigger`, and `40-polishing`. If you have used Flye before, then the folder structure will be familiar to you.
Polap adds more folders by following Flye's output folder structure; e.g., `01-contig`.

Polap's has a few options that are often used or implicitly assumed when you execute. If you execute Polap with whole-genome and organelle-genome assemblies, then you do not need to specify them. But, if you need one more organelle-genome assembly, then you need to specify these options. If you create another mtDNA seed contig file, its name should following accordingly as well.

- `--inum 0 --jnum 1`: default option values, 0 is the whole-genome assembly, and 1 is the organelle-genome assembly. You must have `0/mt.contig.name-1` for the seed contig name file.
- `--inum 0 --jnum 2`: 0 is the whole-genome assembly, and 1 is the organelle-genome assembly. You must have '0/mt.contig.name-2' for the seed file.
- `--inum 1 --jnum 2`: 1 is the whole-genome assembly, and 2 is the organelle-genome assembly. You must have '1/mt.contig.name-2' for the seed file.

The whole-genome assembly number should be less than the organelle-genome assembly although Polap would allow the case of the other way around.

## Using BioProject data

Note that using Polap's `assemble` menu is experimential or not tested yet!
One could try to use BioProject datasets to assemble mtDNA genomes. The following is a template bash script where you could replace `long` with SRA accession for an Oxford Nanopore long-read dataset, and `short` with that of an Illumina short-read dataset.

```bash
polap x-ncbi-fetch-sra --sra long
polap x-ncbi-fetch-sra --sra short
cp -s long.fastq l.fq
cp -s short_1.fastq s1.fq
cp -s short_2.fastq s2.fq
polap assemble
```

The following examples are the cases where one could relativetly easily locate the mitochondrial-origin seed contigs in the Flye whole-genome assembly.

1. mtDNA assembly of _Spirodela polyrhiza_:

```bash
polap x-ncbi-fetch-sra --sra SRR11472010
polap x-ncbi-fetch-sra --sra SRR11472009
cp -s SRR11472010.fastq l.fq
cp -s SRR11472009_1.fastq s1.fq
cp -s SRR11472009_2.fastq s2.fq
polap assemble
```

2. mtDNA assembly of _Taraxacum mongolicum_:

```bash
polap x-ncbi-fetch-sra --sra SRR19182970
polap x-ncbi-fetch-sra --sra SRR19182971
cp -s SRR19182970.fastq l.fq
cp -s SRR19182971_1.fastq s1.fq
cp -s SRR19182971_2.fastq s2.fq
polap assemble
```

3. mtDNA assembly of _Anthoceros agrestis_:

```bash
polap x-ncbi-fetch-sra --sra SRR10190639
polap x-ncbi-fetch-sra --sra SRR10250248
cp -s SRR10190639.fastq l.fq
cp -s SRR10250248_1.fastq s1.fq
cp -s SRR10250248_2.fastq s2.fq
polap assemble
```

## Detailed instructions of installation

### Using Bioconda packages

1. Install Miniconda by following the installation procedure:

   ```bash
   $ curl -OL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
   $ bash Miniconda3-latest-Linux-x86_64.sh
   ```

2. Log out and log in back to activate a conda environment. The following prompt with base indicates that you are ready to use conda environments:
   ```bash
   (base) $
   ```
3. Configure Conda so that it would allow to install [Bioconda](https://bioconda.github.io/) packages:
   ```bash
   (base) $ conda update -y -n base conda
   (base) $ conda config --prepend channels bioconda
   (base) $ conda config --prepend channels conda-forge
   ```
4. Check your default channels and the order of the three channels:
   ```bash
   (base) $ conda config --show channels
   channels:
      - conda-forge
      - bioconda
      - defaults
   ```
5. Install Polap:
   ```bash
   (base) $ conda create --name polap bioconda::polap
   ```
6. Test Polap:
   ```bash
   (base) $ conda activate polap
   (polap) $ git clone https://github.com/goshng/polap.git
   (polap) $ cd polap/test
   (polap) $ polap assemble --test
   ```
7. Check the output:
   ```bash
   output: the assembly graph: 2-oga.gfa
   ```

We need one more step for the installation.
Although it is not actually a part of Polap, we need a short-read polishing step with **FMLRC**, which required a separate conda environment.
It may be complicated, but read along with this **README** where I explain why currently we have to do this one more installation step.

### Short-read polishing Requirement: FMLRC

**FMLRC** required a python version 2, which was too old to be installed along with other Bioconda packages for Polap.
So, we could not create a single Bioconda package for Polap. Instead, users could install one more conda environment for the short-read polishing step.
Here's how to create a separate Conda environment for **FMLRC** while maintaining compatibility with **Polap's** environment. This setup ensures that FMLRC and Polap can function independently without conflicts. Following the previous setup, you'll need to create a new environment specifically for FMLRC. This will prevent any compatibility issues with Polap's main environment.

1. **Create the FMLRC Environment** (note: we need the specific conda environment name of polap-fmlrc for Polap to use it in the short-read polishing step):

   ```bash
   (polap) $ conda deactivate
   (base) $ git clone https://github.com/goshng/polap.git
   (base) $ cd polap
   (base) $ conda env create -f src/polap-conda-environment-fmlrc.yaml
   ```

2. **Activate the FMLRC Environment**:

   ```bash
   (base) $ conda activate polap-fmlrc
   ```

3. **Test-run of FMLRC for Sequence Polishing**:
   Execute the FMLRC polishing commands as needed. Once complete, you can switch back to the Polap environment:
   ```bash
   (polap-fmlrc) $ cd test
   (polap-fmlrc) $ ../src/polap.sh prepare-polishing
   (polap-fmlrc) $ ../src/polap.sh polish
   (polap-fmlrc) $ conda deactivate
   (base) $
   ```

This way, you can seamlessly transition between the two environments, keeping dependencies separate for both Polap and FMLRC.

### Using github source

If you're unable to install Polap via the Bioconda package, you can install it directly from the GitHub source. Ensure that **Miniconda** is installed beforehand, as Polap relies on tools, including **Flye**, **Minimap2**, **Jellyfish**, and others.

0. Make sure that you have **Miniconda** ready before following next steps.

1. Download the source code of Polap available at [Polap](https://github.com/goshng/polap)'s github website:

   ```bash
   (base) $ git clone https://github.com/goshng/polap.git
   (base) $ cd polap
   ```

2. Setup the conda environments for [Polap](https://github.com/goshng/polap).

   ```bash
   (base) $ conda env create -f src/polap-conda-environment.yaml
   ```

3. Setup the conda environment for [FMLRC](https://github.com/holtjma/fmlrc).

   ```bash
   (base) $ conda env create -f src/polap-conda-environment-fmlrc.yaml
   ```

4. Activate polap conda environment.

   ```bash
   (base) $ conda activate polap
   ```

5. Run Polap with a test dataset.

   ```bash
   (polap) $ cd test
   (polap) $ ../src/polap.sh assemble --test
   ```

6. If you see a screen output ending with the something like this:

   ```bash
   output: the assembly graph: 2-oga.gfa
   ```

7. Polishing with a short-read data:
   ```bash
   (polap) $ conda deactivate
   (base) $ conda activate polap-fmlrc
   (polap-fmlrc) $ ../src/polap.sh prepare-polishing
   (polap-fmlrc) $ ../src/polap.sh polish
   ```

Your final mitochondrial genome sequence is `mt.1.fa`.

### Updating the conda Polap package with the github source

Note: We currently have a problem in the Bioconda package. For the time being, follow this script instead.

```bash
git clone https://github.com/goshng/polap.git
polap/src/polap.sh get-revision1
bash polap-revision1.sh patch
bash polap-revision1.sh install-fmlrc
bash polap-revision1.sh test
conda activate polap
```

## Uninstalling Polap

To completely remove Polap and its associated environments, follow these steps:

1. **Deactivate the Polap environment**:

   ```bash
   (polap) $ conda deactivate
   ```

2. **Remove the Polap environment**:

   ```bash
   (base) $ conda remove -n polap --all
   ```

3. **Remove the Polap-FMLRC environment**:
   ```bash
   (base) $ conda remove -n polap-fmlrc --all
   ```

To completely remove Miniconda from an Ubuntu system, follow these steps:

### Step 1: Uninstall Miniconda

1. **Open Terminal** and deactivate any active Conda environments:

   ```bash
   conda deactivate
   ```

2. **Remove Miniconda Installation Directory**:
   Locate your Miniconda installation folder. By default, it is usually in your home directory, e.g., `~/miniconda3`.

   To remove it:

   ```bash
   rm -rf ~/miniconda3
   ```

### Step 2: Remove Environment Variables

1. **Edit Your Shell Profile File** (`~/.bashrc`, `~/.bash_profile`, `~/.zshrc`, etc.), and remove any lines that reference Miniconda. Typically, you’ll see a line like:

   ```bash
   export PATH="/home/yourusername/miniconda3/bin:$PATH"
   ```

   Delete or comment out this line.

2. **Update the Profile File**:
   After editing, reload your shell profile to apply changes:
   ```bash
   source ~/.bashrc
   ```
   Or, if using another shell profile, replace `.bashrc` with the appropriate profile name.

### Step 3: Remove Conda Configuration Files

1. **Remove Conda Configuration Files**:
   Delete the `.conda` and `.condarc` files from your home directory:
   ```bash
   rm -rf ~/.conda ~/.condarc
   ```

### Step 4: Verify Removal

1. **Verify if Conda is Fully Removed**:
   Ensure there is no Conda left in your system path by running:

   ```bash
   which conda
   ```

   It should return nothing.

2. **Check Environment Variables**:
   Ensure the Conda binary is no longer in your `PATH`:
   ```bash
   echo $PATH
   ```
   Confirm that the path to Miniconda is no longer present.

## More

### Options

#### Menu options

- `assemble1`: Flye whole-genome assembly
- `annotate`: Organelle gene annotation. It is part of `assemble1` menu.
- `assemble2`: Flye organelle-genome assembly
- `prepare-polishing`: [FMLRC](https://github.com/holtjma/fmlrc) short-read polishing preparation
- `polish`: [FMLRC](https://github.com/holtjma/fmlrc) short-read polishing

#### General options

```text
POLAP - Plant organelle DNA long-read assembly pipeline.
version v0.3.7-eadc43f

Usage: polap [<menu> [<menu2> [<menu3>]]] [-o|--outdir <arg>]
    [-i|--inum <arg>] [-j|--jnum <arg>] [-w|--single-min <arg>]
    [-l|--long-reads <arg>] [-a|--short-read1 <arg>] [-b|--short-read2 <arg>]
    [-m|--min-read-length <arg>] [-t|--threads <arg>] [--test] [--log <arg>]
    [--random-seed <arg>] [--version] [-h|--help]

Assemble mitochondrial DNA (mtDNA) in a single command (not tested yet):
    polap -l <arg> -a <arg> [-b <arg>]
    polap assemble -l <arg> -a <arg> [-b <arg>]

Perform a polishing of the mtDNA sequence utilizing the FMLRC protocol:
    polap prepare-polishing -a <arg> [-b <arg>]
    polap polish -p <arg> -f <arg>

To assemble mitochondrial DNA (mtDNA), follow a two-step process involving
manual seed contig selection:
    polap assemble1 -l <arg> -a <arg> [-b <arg>] [-m <arg>]
    polap assemble2 -i <arg> -j <arg> [-w <arg>] [-c <arg>]

To assemble mitochondrial DNA (mtDNA), follow a three-step process
that utilizes semi-automatic seed contig selection (not working yet):
    polap assemble1 -l <arg> -a <arg> [-b <arg>] [-m <arg>]
    polap seeds -i <arg> -j <arg>
    polap assemble2 -i <arg> -j <arg> [-w <arg>] [-c <arg>]

To assemble mitochondrial DNA (mtDNA), follow a series of sequential steps:
    polap init -o <arg>
    polap total-length-long -l <arg>
    polap find-genome-size -a <arg> [-b <arg>]
    polap reduce-data -l <arg> [-m <arg>]
    polap summary-reads -a <arg> [-b <arg>]
    polap flye1 [-t <arg>]
    polap edges-stats -i <arg>
    polap annotate -i <arg>
    polap seeds [-i <arg>] -j <arg>
    polap map-reads [-i <arg>] -j <arg>
    polap test-reads [-i <arg>] -j <arg> -w <arg> [-c <arg>]
    polap select-reads [-i <arg>] -j <arg> -w <arg> [-c <arg>]
    polap flye2 [-i <arg>] -j <arg>

Others menus:
    polap blast-genome -i <arg>
    polap count-genes -i <arg>
    polap flye-polishing -j <arg>
    polap make-menus
    polap clean-menus
    polap list

BioProject menus:
    polap get-bioproject --bioproject <arg>
    polap bioproject-prepare -o <arg>
    polap get-bioproject-sra --sra <arg>
    polap get-mtdna --species <arg>

Other options:
    [-p|--unpolished-fasta <arg>] [-f|--final-assembly <arg>]
    [-c|--coverage <arg>] [--flye-asm-coverage <arg>]
    [--bioproject <arg>] [--species <arg>] [--accession <arg>]
    [--query <arg>] [--subject <arg>]
    [--no-reduction-reads] [--no-coverage-check]
    [--plastid]
    [--archive <arg>]
    [--sra <arg>] [-g|--genomesize <arg>]

menu: assemble, assemble1, annotate, assemble2, flye-polishing,
    make-menus, list, clean-menus, cleanup, init,
    summary-reads, total-length-long, find-genome-size, reduce-data, flye1
    blast-genome, count-gene, seeds,
    prepare-seeds, map-reads, test-reads, select-reads, flye2,
    flye-polishing, prepare-polishing, polish,
    version

Options:
  -o, --outdir: output folder name (default: o)
    The option '-o' or '--outdir' specifies the output folder name,
    with a default value of 'o'. The output folder typically contains input
    files that are long-read and short-read data files. Input data files can
    be specified using the options provided by -l, -a, and -b.

  -l, --long-reads: long-reads data file in fastq format (default: l.fq)
    The option '-l' or '--long-reads' specifies the location of a long-reads
    data file in fastq format, with a default filename of 'l.fq'.

  -a, --short-read1: short-read fastq file 1 (default: s1.fq)
    The option '-a' or '--short-read1' specifies the first short-read fastq
    file to be used, with a default value of "s1.fq".

  -b, --short-read2: short-read fastq file 2 (default: s2.fq)
    The option '-b' or '--short-read2' specifies a short-read fastq file 2,
    with a default value of 's2.fq'. The second short-read data file,
    if provided, is considered optional.
    (Note: not tested yet; -a & -b are required.)

  -i, --inum: previous output number of organelle-genome assembly (default: 0)
    The option '-i' or '--inum' specifies the previous output number
    of an organelle-genome assembly, with a default value of '0'.
    The zero for this option specifies the whole-genome assembly.

  -j, --jnum: current output number of organelle-genome assembly (default: 1)
    The option '-j' or '--jnum' allows users to specify the current output number
    for an organelle-genome assembly, with a default value of '1'.

  -w, --single-min: minimum mapped bases or PAF 11th column (default: 3000)
    This parameter ensures that the alignment level between a long-read and a
    seed contig is properly controlled. For plant mitochondrial DNAs, a DNA
    fragment size of approximately 3 kilobases appears to be more effective
    than the smaller 1-kilobase fragment. In the case of plastid DNAs, a
    fragment size of 1 kilobase (kb) might be more suitable, requiring an
    adjustment to the -m option accordingly.

  --polap-reads: use intra- and inter-contig read selection (default: off)
    The default read selection is ptGAUL's approach.
    This option allows long reads that are mapped within a seed contig and
    between two contigs.

  -m, --min-read-length: minimum length of long reads (default: 3000)
    The option '-m' or '--min-read-length' specifies the minimum length of
    long reads, with a default value of 3000.

  -t, --threads: number of CPUs (default: maximum number of cores)
    The option '-t' or '--threads' specifies the number of CPU threads to
    utilize, with a default value equal to the maximum number of available
    cores.

  -c, --coverage: coverage for the organelle-genome assembly (default: 50)
    The option '-c' or '--coverage' specifies the coverage percentage for the
    organelle-genome assembly for controlling the data size.

  -g, --genomesize: expected genome size (default: estimated with a short-read dataset)
    Users can assemble an organelle genome when they have a genome size
    estimate. But, we require a short-read dataset to determine the genome size
    for the whole-genome assembly process. Polishing a long-read assembly
    necessitates the use of a short-read dataset.

  -p, --unpolished-fasta: polishing sequence in fasta format (default: mt.0.fasta)
    The option enables the polishing of sequences in a FASTA format,
    with the default output file being named 'mt.0.fasta'.

  -f, --final-assembly: final assembly in fasta format (default: mt.1.fa)
    The final assembly in FASTA format, with a default file name of 'mt.1.fa'.

  --no-reduction-reads: reduction of long-read data before assemble1
    In the process of whole-genome assembly, we utilize a reduced amount of
    long-read data. By default, we reduce the size of a long-read dataset prior
    to performing a whole-genome assembly.
    Note: The size of coverage is set by --coverage (default: 50)

  --no-coverage-check: coverage check before assemble2 step
    By default, in the process of assembling organelle genomes, we reduce
    the size of selected seed reads.
    Note: The size of coverage is set by --coverage (default: 50)

  --yes: always yes for a question or deletes output completely (off by default)

  --redo: redo a POLAP pipeline (off by default)
    The command specifies that any previously generated intermediate results
    should be disregarded and new calculations performed from scratch.

  --select-read-range: start,end,number for the range of read selection (default: 3000,39000,7)
    It specifies the values for ptGAUL read-selection minimum number of
    bases or ratios. For the start and end values of a ratio, real numbers must
    fall within the range of 0 to 1.
    Note: refer to the menu "test-reads" for help.

  --start-index: used by test-reads

  --random-seed: 5-digit number (default automatically assigned)
    To ensure reproducibility, you can supply a random number seed
    to facilitate sampling of reads.
    seqkit sample random seed; 11 used in seqkit sample.

  --flye-asm-coverage: Flye --asm-coverage (default: 30)
    Flye --asm-coverage is a parameter used with the assembly coverage of Flye.

  --no-flye-asm-coverage: no use of Flye --asm-coverage
    The flag '--no-flye-asm-coverage' indicates that we use Flye option
    neither --asm-coverage nor --genome-size in flye execution.
    This option is the same as --flye-asm-coverage set to 0.
    Note: not tested yet!

  --bridge-same-strand: (default: off)
    When linking two inverted repeats, enabling this feature ensures that
    the strands are equal for the two mapped IR contigs.
    Note: currently only plus strand is used.

  --log: log file (default: <output>/polap.log)
    The log file option allows users to specify a custom log file location,
    with a default setting of '<output>/polap.log'.

  --clock: display the start and ending time (default: off)
    The clock option allows users to display both the start and end times.

  --species: Species scientific name (no default)
  --sra: SRA data (no default)
  -v, --verbose: use multiple times to increase the verbose level
  --version: Prints version
  -h: Prints polap global help
  --help: Prints menu help

Example:
  git clone https://github.com/goshng/polap.git
  cd polap/test
  polap --test

Place your long-read and short-read files at a folder:
  long-read file: l.fq
  short-read file: s1.fq, s2.fq

Execute: polap init
```

#### Detailed description of some of the options

`-o` _DIRECTORY_, `--outdir` _DIRECTORY_ : Write output files to _DIRECTORY_. The default output directory is
`o`.

`-l` _FILE_, `--long-reads` _FILE_ : Specify the long-read FASTQ file. If this option is not specified,
the default long-read file will be used. It will be the `l.fq` at
the current directory. The FASTQ file must be a regular file not
compressed one.

`-a` _FILE_, `--short-read1` _FILE_ : Specify the first short-read FASTQ file. If this option is not
specified, the default long-read file will be used. It will be the
`s1.fq` at the current directory. The FASTQ file must be a regular
file not compressed one. Two short-read data files are required.

`-b` _FILE_, `--short-read2` _FILE_ : Specify the second short-read FASTQ file. If this option is not
specified, the default 1st short-read file will be used. It will be
the `s2.fq` at the current directory. The FASTQ file must be a
regular file not compressed one. Two short-read data files are
required.

`-w` _NUMBER_, `--single-min` _NUMBER_ : Specify the minimum mapped bases or PAF 11th column (default:
`3000`).

`-i` _NUMBER_, `--inum` _NUMBER_ : Specify the previous output number of organelle-genome assembly
(default: `0`).

`-j` _NUMBER_, `--jnum` _NUMBER_ : Specify the current output number of organelle-genome assembly
(default: `1`).

`-m` _NUMBER_, `--min-read-length` _NUMBER_ : Specify the minimum length of long-read sequences. If this option is
not specified, the default minimum length will be used. It will be
the `3000`.

`-t` _NUMBER_, `--threads` _NUMBER_ : Specify the numbre CPU cores to use (default: `0`).

`-c` _NUMBER_, `--coverage` _NUMBER_ : Specify coverage for the 2nd assembly (default: `30`).

`-p` _FILE_, `--unpolished-fasta` _FILE_ : Specify the unpolished FASTA file. If this option is not specified,
the default 2nd short-read file will be used. It will be the
`mt.0.fasta` at the current directory.

`-f` _FILE_, `--final-assembly` _FILE_ : Specify the final genome assembly FASTA file. If this option is not
specified, the default long-read file will be used. It will be the
`mt.1.fa` at the current directory.

`--version` : Print version.

`--help` : Show usage message of a menu.

`-h` : Show usage message.

# Authors

Sang Chul Choi

# License

This project includes components licensed under different licenses: See the [LICENSE](https://github.com/goshng/polap/blob/main/LICENSE) file for details.
It is mainly released under the [GNU General Public License version 3](https://www.gnu.org/licenses/gpl-3.0-standalone.html),
and include other 3rd-party components released under the [BSD-3-clause](https://opensource.org/license/bsd-3-clause):
this software carries no warranty of any kind.
