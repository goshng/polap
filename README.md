# Polap: Plant Mitochondrial DNA Assembly Pipeline

Polap is a specialized pipeline designed to assemble plant mitochondrial DNA (mtDNA) using the **Flye** long-read assembler, supplemented by organellar gene annotation to guide contig selection. Inspired by **ptGAUL** (a chloroplast assembly pipeline), Polap integrates Flye's capabilities with gene annotation to identify and assemble mtDNA contigs from whole-genome data.

## Description

**Polap** is a specialized meta-assembly pipeline that enhances the functionality of the [Flye](https://github.com/fenderglass/Flye) assembler. Serving as a Flye helper pipeline, Polap incorporates organellar gene annotation, aiding in the educated guess of selection of organelle-origin contigs from Flye's whole-genome assembly. The name "Polap" is inspired by the [ptGAUL](https://github.com/Bean061/ptgaul) chloroplast genome assembly pipeline created by Wenbin Zhou, as Polap similarly focuses on organellar genomes, specifically plant mitochondrial DNA (mtDNA).

The pipeline combines the **Flye** assembler with [FMLRC](https://github.com/holtjma/fmlrc) (a sequence polishing tool by Matt Holt) and applies sequencing depth analysis to identify contigs belonging to the mitochondrial genome. Leveraging both Flye-generated assemblies and organellar gene annotations, Polap helps users identify candidate mtDNA contigs based on sequencing depth and gene presence, enabling targeted selection and assembly of mtDNA.

Polap originated as an extension of ptGAUL, with a specific focus on assembling plant mtDNA from long-read sequencing data. The workflow begins by conducting a Flye whole-genome assembly, generating contigs without a reference genome. These contigs are then analyzed with [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) to confirm the presence of mtDNA or chloroplast (ptDNA) genes. Using depth and copy number metrics, Polap allows users to identify likely mtDNA contigs, collects them, and performs a final targeted mtDNA assembly. The pipeline concludes with polishing using FMLRC, resulting in an accurate, high-quality mtDNA assembly.

## Key Features

- **Flexible Assembly**: Polap utilizes Flye as a "black-box" assembler, adding custom logic to filter and identify organelle-specific contigs.
- **Rough Organellar Gene Annotation**: Enables Polap to select contigs from mitochondrial and chloroplast origins, focusing specifically on the target mtDNA.

## How It Works

1. **Initial whole-genome assembly**: Flye performs a whole-genome assembly on the input long reads. We need paired-end short-read data to estimate the genome size using **JellyFish**.
2. **Seed contig selection**: Using BLAST, Polap identifies mtDNA-related contigs based on the presence of organelle genes. Users need to manually select mitochondrial-origin contigs using organelle the gene annotation table that is made by **Flye** and modified by **Polap**.
3. **Organelle-gename assembly**: The selected mtDNA seed contigs and long-read sequencing data are aligned using **Minimap2** and those reads that are well mapped are fed into **Flye** to assembly the organellar genome. We borrow the approach of **ptGAUL** in this long-read selection.
4. **Final Polishing**: The selected mtDNA contigs undergo polishing with FMLRC for accurate genome assembly. **Sequence Polishing** uses **FMLRC** for high-quality, error-corrected assemblies of the target mitochondrial genome as suggested by **ptGAUL**.

## Requirements

- **Operating System**: Linux (not compatible with macOS or Windows)
- **Dependencies**: Requires BASH (>= 3.0) and **Miniconda**
- **Installation**: Can be installed via Bioconda or manually from the GitHub source.

## Quick Start

### Using Bioconda package

1. Install Miniconda:
   ```bash
   $ curl -OL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
   $ bash Miniconda3-latest-Linux-x86_64.sh
   ```
2. Log out and log in again to activate a conda environment:
   ```bash
   (base) $
   ```
3. Configure Conda:
   ```bash
   (base) $ conda update -y -n base conda
   (base) $ conda config --prepend channels bioconda
   (base) $ conda config --prepend channels conda-forge
   ```
4. Check your default channels:
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
   ERROR: no seed contig file of the contig selection types!
   ```

For more detailed installation and usage instructions, please read along with this **README**.

### Short-read polishing Requirement: FMLRC

Here's how to create a separate Conda environment for **FMLRC** while maintaining compatibility with **Polap's** environment. This setup ensures that FMLRC and Polap can function independently without conflicts. Following the previous setup, you'll need to create a new environment specifically for FMLRC. This will prevent any compatibility issues with Polap's main environment.

1. **Create the FMLRC Environment**:

   ```bash
   (polap) $ conda deactivate
   (polap) $ git clone https://github.com/goshng/polap.git
   (polap) $ cd polap
   (base) $ conda env create -f src/polap-conda-environment-fmlrc.yaml
   ```

2. **Activate the FMLRC Environment**:

   ```bash
   (base) $ conda activate polap-fmlrc
   ```

3. **Run FMLRC for Sequence Polishing**:
   When you reach the polishing step in Polap, activate the FMLRC environment:
   ```bash
   (polap) $ conda activate polap-fmlrc
   ```
   Execute the FMLRC polishing commands as needed. Once complete, you can switch back to the Polap environment:
   ```bash
   (polap-fmlrc) $ src/polap.sh prepare-polishing
   (polap-fmlrc) $ conda deactivate
   ```

This way, you can seamlessly transition between environments, keeping dependencies separate for both Polap and FMLRC.

### Using github source

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
   ERROR: no seed contig file of the contig selection types!
   ```

7. Polishing with a short-read data:
   ```bash
   (polap) $ conda deactivate
   (base) $ conda activate polap-fmlrc
   (polap-fmlrc) $ ../src/polap.sh prepare-polishing
   (polap-fmlrc) $ ../src/polap.sh polish
   ```

Your final mitochondrial genome sequence is `mt.1.fa`.

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
   (polap) $ src/polap.sh
   ```

For simplicity, the instructions below assume you installed Polap using Bioconda. In cases where you cloned the GitHub repository, replace `polap` with `src/polap.sh` as needed. If you're experienced with the command line, this setup will feel straightforward. Otherwise, follow the steps in this README closely.

### Input Data

1. **Set Up Your Input Folder**:  
   Create a folder and place your input dataset within it. This dataset should consist of:

   - One Oxford Nanopore long-read FASTQ file
   - Two short-read FASTQ files

   A set of sample FASTQ files is available in the `test` folder if you've cloned the Polap repository. Currently, only configurations with two short-read files are supported.

2. **File Requirements**:  
   If you don't specify input file paths with options, Polap will expect three files named `l.fq`, `s1.fq`, and `s2.fq` in the current working directory where you execute the Polap script. Make sure these files are uncompressed. You can also specify the files using the options `-l`, `-a`, and `-b`.

3. **Help Message**:  
   Run Polap without any options to display a help message:
   ```bash
   (polap) $ polap
   ```

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

After completing the whole-genome assembly, you can visualize the genome assembly graph in `o/0/30-contigger/graph_final.gfa` using **[Bandage](https://rrwick.github.io/Bandage/)**. This graph shows contigs, some originating from mitochondrial or plastid DNA, alongside nuclear DNA.

If you can identify mitochondrial DNA contigs manually, you can proceed with those. Otherwise, continue the Polap pipeline to automatically annotate organellar regions:

```bash
(polap) $ polap annotate
```

This step also requires some time to complete, so be patient.

### Organelle Contig Selection

After completing the `polap annotate` step, you'll find organelle gene annotations in several files within the `o/0/` directory:

- `contig-annotation-depth-table.txt`
- `contig-annotation-table.txt`
- `assembly_info_organelle_annotation_count-all.txt`
- `assembly_info_organelle_annotation_count.txt`

Inspect these files to identify contigs with a higher number of mitochondrial (MT) genes compared to plastid (PT) genes. This selection step can be subjective, so you may need to carefully review the annotation table files to choose contigs that likely originate from plant mitochondrial DNA. These annotation tables are derived from `o/0/30-contigger/graph_final.gfa` and considered modified versions of `contigs_stats.txt` from Flye and provide information about each contig assembled from the whole-genome assembly.

You’ll use both the graph visualization of the genome assembly `o/0/30-contigger/graph_final.gfa`, and the annotation table files to select candidate contigs. Once identified, prepare a text file (`o/0/mt.contig.name-1`) with each selected edge sequence name on a new line. Note that edge sequence names should start with `edge_`, not `contig_`. An example file might look like this:

```
edge_1
edge_2
edge_3
```

For more information on preparing a mitochondrial contig file, refer to the [MT contig name](#mt-contig-name) section below.

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

The final polished mitochondrial genome sequence will be saved as `mt.1.fa`.

### Uninstalling Polap

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

### Sample Datasets

A sample dataset demonstrating Polap’s capabilities is available on Figshare: [Polap data analysis of 11 datasets](https://figshare.com/s/07305fd4e18c74080fbc).

### Creating the MT Contig Name File

To prepare an MT contig name file, list the contig names you suspect are of mitochondrial origin. Select these contigs based on three key features:

1. **Genome assembly graph**
2. **Presence of organelle genes**
3. **Contig copy numbers**

Contigs are labeled as edges, each suffixed with a number. To examine the whole-genome assembly graph, use Bandage to open `o/0/30-contigger/graph_final.gfa`. You can also view organelle annotations and copy numbers with the following command:

```bash
$ column -t o/0/30-contigger/contigs_stats.txt
```

For additional guidance on identifying mtDNA contigs, a [YouTube](https://youtu.be/29OaiFCqAzI) tutorial in Korean is available. The following example shows a Flye assembly table with gene counts, where the `Copy` column indicates the copy number, `MT` represents mitochondrial genes, and `PT` represents plastid genes. The `Edge` column lists edge numbers that make up each contig.

```bash
$ column -t o/assembly_info_organelle_annotation_count.txt
```

Example table:

| Contig   | Length | V3  | V4  | V5  | Copy | V7   | V8  | MT  | PT  | Edge   |
| -------- | ------ | --- | --- | --- | ---- | ---- | --- | --- | --- | ------ |
| contig_3 | 70078  | 6   | N   | N   | 1    | both | \*  | 19  | 41  | 1,3,-1 |
| contig_1 | 39736  | 6   | N   | Y   | 1    | left | \*  | 16  | 16  | 1      |
| contig_2 | 19223  | 3   | N   | N   | 1    | both | \*  | 11  | 13  | 2      |

After reviewing the table, edit `o/0/mt.contig.name-1` to include the names of candidate mtDNA contigs. Each line should list an edge, without any empty lines. Here’s an example:

```
edge_1
edge_2
edge_3
```

# Options

## Menu options

`assemble1`

: Flye whole-genome assembly

`annotate`

: Organelle gene annotation

`assemble2`

: Flye organelle-genome assembly

`prepare-polishing`

: [FMLRC](https://github.com/holtjma/fmlrc) short-read polishing
preparation

`polish`

: [FMLRC](https://github.com/holtjma/fmlrc) short-read polishing

## General options

POLAP - Plant organelle DNA long-read assembly pipeline.
version v0.3.7-eadc43f

Usage: polap [<menu> [<menu2> [<menu3>]]] [-o|--outdir <arg>]
[-l|--long-reads <arg>] [-a|--short-read1 <arg>] [-b|--short-read2 <arg>]
[-i|--inum <arg>] [-j|--jnum <arg>] [-w|--single-min <arg>]
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
polap summary-reads -a <arg> [-b <arg>]
polap total-length-long -l <arg>
polap find-genome-size -a <arg> [-b <arg>]
polap reduce-data -l <arg> [-m <arg>]
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
[--sra <arg>] [-x|--bridge-min <arg>] [-g|--genomesize <arg>]

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

-w, --single-min: minimum mapped bases or PAF 11th column (default: 3000)
This parameter ensures that the alignment level between a long-read and a
seed contig is properly controlled. For plant mitochondrial DNAs, a DNA
fragment size of approximately 3 kilobases appears to be more effective
than the smaller 1-kilobase fragment. In the case of plastid DNAs, a
fragment size of 1 kilobase (kb) might be more suitable, requiring an
adjustment to the -m option accordingly.

-x, --bridge-min: minimum bridging read length or PAF 2nd column (default: 0)
Note: it is not tested yet.

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

--polap-reads: use intra- and inter-contig read selection (default: off)
The default read selection is ptGAUL's approach.
This option allows long reads that are mapped within a seed contig and
between two contigs.

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

`-o` _DIRECTORY_, `--outdir` _DIRECTORY_

: Write output files to _DIRECTORY_. The default output directory is
`o`.

`-l` _FILE_, `--long-reads` _FILE_

: Specify the long-read FASTQ file. If this option is not specified,
the default long-read file will be used. It will be the `l.fq` at
the current directory. The FASTQ file must be a regular file not
compressed one.

`-a` _FILE_, `--short-read1` _FILE_

: Specify the first short-read FASTQ file. If this option is not
specified, the default long-read file will be used. It will be the
`s1.fq` at the current directory. The FASTQ file must be a regular
file not compressed one. Two short-read data files are required.

`-b` _FILE_, `--short-read2` _FILE_

: Specify the second short-read FASTQ file. If this option is not
specified, the default 1st short-read file will be used. It will be
the `s2.fq` at the current directory. The FASTQ file must be a
regular file not compressed one. Two short-read data files are
required.

`-p` _FILE_, `--unpolished-fasta` _FILE_

: Specify the unpolished FASTA file. If this option is not specified,
the default 2nd short-read file will be used. It will be the
`mt.0.fasta` at the current directory.

`-f` _FILE_, `--final-assembly` _FILE_

: Specify the final genome assembly FASTA file. If this option is not
specified, the default long-read file will be used. It will be the
`mt.1.fa` at the current directory.

`-m` _NUMBER_, `--min-read-length` _NUMBER_

: Specify the minimum length of long-read sequences. If this option is
not specified, the default minimum length will be used. It will be
the `3000`.

`-t` _NUMBER_, `--threads` _NUMBER_

: Specify the numbre CPU cores to use (default: `0`).

`-c` _NUMBER_, `--coverage` _NUMBER_

: Specify coverage for the 2nd assembly (default: `30`).

`-r` _NUMBER_, `--pair-min` _NUMBER_

: Specify the minimum mapped bases or PAF 11th column (default:
`3000`).

`-x` _NUMBER_, `--bridge-min` _NUMBER_

: Specify the minimum bridging read length or PAF 7th column (default:
`3000`).

`-w` _NUMBER_, `--single-min` _NUMBER_

: Specify the minimum mapped bases or PAF 11th column (default:
`3000`).

`-i` _NUMBER_, `--inum` _NUMBER_

: Specify the previous output number of organelle-genome assembly
(default: `0`).

`-j` _NUMBER_, `--jnum` _NUMBER_

: Specify the current output number of organelle-genome assembly
(default: `1`).

`-v`, `--version`

: Print version.

`-h`, `--help`

: Show usage message.

# Authors

Copyright 2024 Sungshin Women's University. Released
under the
[GNU General Public License version 3](https://www.gnu.org/licenses/gpl-3.0-standalone.html)
This software carries no warranty of any kind.
