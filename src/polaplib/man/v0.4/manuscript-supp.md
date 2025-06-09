<!-- \blandscape -->

\newpage

# Tables

<!--table1-->
<!-- polap-data-cflye man table-benchmark some 2 -->

Table: Plastid genome assemblies for 23 plant species datasets using subsampled sequencing data with a maximum subsampling rate of 5% (Run Setting A).
All datasets were downsampled to 10x genome coverage, and Stage 1 included 10 subsampling steps (N).
Depending on the dataset, different maximum sampling rates (P) were used in Stage 1, and different replicate sizes (R) were applied in Stages 2 and 3.
NA represents no assemblies in the subsampling-based method and no comparison available.
{#tbl:benchmark-polap-some-2}

!include figures/table-benchmark-polap-some-2.md

<!-- \elandscape -->

<!--table2-->
<!-- polap-data-cflye man table-polap-disassemble Eucalyptus_pauciflora 2 -->

\newpage

Table: Three stages of subsampling-based plastid genome assembly for the _Eucalyptus pauciflora_ dataset with Run Setting A.
The configuration includes an increasing subsample size up to a maximum subsampling rate of 5%, a step size of 10 in Stage 1, 5 replicates in Stages 2 and 3, and a maximum memory limit of 16 GB.
Abbreviations are as follows:
iteration in each Stage (I), subsampling rate (Rate) and read-coverage threshold (Alpha);
assembly metrics including the number of segments in the assembly (N),
the total length of these segments (L), and
the number of circular genome paths detected (C); and
the draft plastid genome assembly length (Length).
Alpha at Stage 3 is the percent identity values between consecutive indices.
{#tbl:polap-disassemble-Eucalyptus_pauciflora-2}

!include figures/table-polap-disassemble-Eucalyptus_pauciflora-2.md

\newpage

# Figures

<!--figure1-->

![Workflow of the subsampling-based plastid genome assembly. The genome assembly procedure is applied repeatedly in Stages 1 and 2.](figures/polap2-figure1.pdf){#fig:mainfigure1 width=100%}

\setcounter{table}{0}
\setcounter{figure}{0}
\renewcommand{\thetable}{S\arabic{table}}
\renewcommand{\thefigure}{S\arabic{figure}}

\blandscape

# Supplementary Materials

<!--tableS1-->
<!-- polap-data-cflye man table-benchmark some 2 -->

Table: Sequencing data for the datasets, including species names and their corresponding taxonomic ranks studied. {#tbl:benchmark-data-some-2}

!include figures/table-benchmark-data-some-2.md

<!-- \begingroup -->
<!-- \tiny -->
<!-- \setlength{\parskip}{2pt} -->
<!---->
<!-- Abbreviations are as follows: -->
<!-- SRA ID for the long-read dataset (L_SRA), -->
<!-- the size of the long-read dataset (L_size), -->
<!-- the long-read coverage (L_cov), -->
<!-- SRA ID for the short-read dataset (S_SRA), -->
<!-- the size of the short-read dataset (S_size), and -->
<!-- the short-read coverage (S_cov). -->
<!---->
<!-- \endgroup -->
<!---->

\elandscape

<!--tableS:computer-->
<!-- polap-data-cflye man table-benchmark some 2 -->

\newpage

Table: Computer setup for the 23 datasets.
{#tbl:benchmark-computer}

!include figures/table-benchmark-computer-some-2.md

<!--tableS2-->
<!-- polap-data-cflye man table-benchmark some 1 -->

\newpage

<!-- \blandscape -->

Table: Replicate of plastid genome assemblies for 23 plant species datasets using subsampled sequencing data with a maximum subsampling rate of 5% (Run Setting A).
All datasets were downsampled to 10x genome coverage, and Stage 1 included 10 subsampling steps (N).
Depending on the dataset, different maximum sampling rates (P) were used in Stage 1, and different replicate sizes (R) were applied in Stages 2 and 3.
{#tbl:benchmark-polap-some-1}

!include figures/table-benchmark-polap-some-1.md

<!-- \elandscape -->

<!--table:benchmark-time-->
<!-- polap-data-cflye man table-benchmark some 2 -->

\newpage

\blandscape

Table: Benchmark of `GetOrganelle`, `ptGAUL`, `PMAT`, `TIPPo`, `Oatk` and the method (Run Setting A) presented here in terms of data processing time.
NA at the column of `NextDenovo` represents no error-corrected long-read results, resulting in no assemblies in the correction-then-assembly pipelines including `PMAT`, `TIPPo`, and `Oatk`.
Abbreviations are as follows:
`GetOrganelle` (GO),
`ptGAUL` (pG),
`NextDenovo` (ND),
`PMAT` with `-fc 0.1` (P0.1),
`PMAT` with `-fc 1.0` (P1.0),
`TIPPo` with `-p onthq` (Thq),
`TIPPo` with `-p ont` (Tont),
`Oatk` with `-c 30` (O30), and
`Oatk` with `-c 20` (O20).
{#tbl:benchmark-time-2}

!include figures/table-benchmark-time-some-2.md

\elandscape

<!--table:benchmark-memory-->
<!-- polap-data-cflye man table-benchmark some 2 -->

\newpage

\blandscape

Table: Benchmark of `GetOrganelle`, `ptGAUL`, `PMAT`, `TIPPo`, `Oatk` and the method (Run Setting A) in terms of peak memory.
NA at the column of `NextDenovo` represents no error-corrected long-read results, resulting in no assemblies in the correction-then-assembly pipelines including `PMAT`, `TIPPo`, and `Oatk`.
Abbreviations are as follows:
`GetOrganelle` (GO),
`ptGAUL` (pG),
`NextDenovo` (ND),
`PMAT` with `-fc 0.1` (P0.1),
`PMAT` with `-fc 1.0` (P1.0),
`TIPPo` with `-p onthq` (Thq),
`TIPPo` with `-p ont` (Tont),
`Oatk` with `-c 30` (O30), and
`Oatk` with `-c 20` (O20).
{#tbl:benchmark-memory-2}

!include figures/table-benchmark-memory-some-2.md

\elandscape

<!--tableS3-->
<!-- polap-data-cflye man table-benchmark some 0 -->

\newpage

<!-- \blandscape -->

Table: Plastid genome assemblies for 23 plant species datasets using subsampled sequencing data with a maximum subsampling rate of 10% (Run Setting B).
All datasets were downsampled to 10x genome coverage, and Stage 1 included 10 subsampling steps (N).
Depending on the dataset, different maximum sampling rates (P) were used in Stage 1, and different replicate sizes (R) were applied in Stages 2 and 3.
NA represents no assemblies in the subsampling-based method and no comparison available.
{#tbl:benchmark-polap-some-0}

!include figures/table-benchmark-polap-some-0.md

<!-- !include manuscript-table1-footnote.md -->

<!-- \elandscape -->

<!--tableS4-->
<!-- polap-data-cflye man table-benchmark some 4 -->

\newpage

Table: Plastid genome assemblies for 23 plant species datasets using subsampled sequencing data with a maximum subsampling rate of 1%.
All datasets were downsampled to 10x genome coverage, and Stage 1 included 10 subsampling steps (N).
Depending on the dataset, different maximum sampling rates (P) were used in Stage 1, and different replicate sizes (R) were applied in Stages 2 and 3.
NA represents no assemblies in the subsampling-based method and no comparison available.
{#tbl:benchmark-polap-some-4}

!include figures/table-benchmark-polap-some-4.md

<!--tableS6-->
<!-- polap-data-cflye man table-polap-disassemble Eucalyptus_pauciflora 0 -->

\newpage

Table: Three stages of subsampling-based plastid genome assembly for the Eucalyptus pauciflora dataset with Run Setting B.
The configuration includes an increasing subsample size up to a maximum subsampling rate of 10%, a step size of 10 in Stage 1, 5 replicates in Stages 2 and 3, and a maximum memory limit of 16 GB.
Abbreviations are as follows:
iteration in each Stage (I), subsampling rate (Rate) and read-coverage threshold (Alpha);
assembly metrics including the number of segments in the assembly (N),
the total length of these segments (L), and
the number of circular genome paths detected (C); and
the draft plastid genome assembly length (Length).
Alpha at Stage 3 is the percent identity values between consecutive indices.
{#tbl:polap-disassemble-Eucalyptus_pauciflora-0}

!include figures/table-polap-disassemble-Eucalyptus_pauciflora-0.md

<!--tableS7-->
<!-- polap-data-cflye man table-polap-disassemble Eucalyptus_pauciflora 3 -->

\newpage

Table: Three stages of subsampling-based plastid genome assembly for the Eucalyptus pauciflora dataset with Run Setting C.
The configuration includes an increasing subsample size up to a maximum subsampling rate of 10%, a step size of 50 in Stage 1, 5 replicates in Stages 2 and 3, and a maximum memory limit of 16 GB.
Abbreviations are as follows:
iteration in each Stage (I), subsampling rate (Rate) and read-coverage threshold (Alpha);
assembly metrics including the number of segments in the assembly (N),
the total length of these segments (L), and
the number of circular genome paths detected (C); and
the draft plastid genome assembly length (Length).
Alpha at Stage 3 is the percent identity values between consecutive indices.
{#tbl:polap-disassemble-Eucalyptus_pauciflora-3}

!include figures/table-polap-disassemble-Eucalyptus_pauciflora-3.md

<!--figureS1-->
<!-- polap-data-cflye man figure-alpha Eucalyptus_pauciflora -->

\newpage

![Line plot of read-coverage thresholds versus subsample size index in Stage 1 of the subsampling-based assemblies for _Eucalyptus pauciflora_.](figures/alpha0.pdf){#fig:alpha0 width=100%}

<!--figureS2-->
<!-- polap-data-cflye man figure-delta Eucalyptus_pauciflora -->

\newpage

![Line plot of increment size versus subsample size index in Stage 1 of the subsampling-based assemblies for _Eucalyptus pauciflora_.](figures/delta.pdf){#fig:delta width=100%}

<!--figureS3-->

<!-- Supporting Material - plastid genome assemblies using the six pipelines -->

\newpage

# Code

`Polap` (Plant Organelle Long-read Assembly Pipeline v0.4.3.7) is available at [http://github.com/goshng/polap](http://github.com/goshng/polap).
A quick start guide of the subsampling-based plastome assembly is provided for use on a Linux system with an Internet connection.
A detailed guide is also available, offering a step-by-step explanation of some of the procedures outlined in the quick start.

### Requirements

- **Operating System**: Linux (not compatible with macOS or Windows)
- **Dependencies**: Requires [Bash](https://www.gnu.org/software/bash/) (>= 5.0) and **[Miniconda](https://docs.anaconda.com/miniconda/miniconda-install/)**

### Quick Start

To replicate the results presented in this manuscript on a Linux computer with `git` installed and an Internet connection, follow the steps below. Most steps complete in a relatively short time, except for the final step, which includes both data downloading and full analysis:

```bash
mkdir -p all/polap/cflye1
cd all/polap/cflye1
git clone https://github.com/goshng/polap.git
bash polap/src/polap-data-cflye -y install conda
```

Log out and back in to the terminal.

```bash
cd all/polap/cflye1
source ~/miniconda3/bin/activate
bash polap/src/polap-data-cflye setup conda
bash polap/src/polap-data-cflye -y install minimal
bash polap/src/polap-data-cflye setup polap
```

Log out and back in to the terminal.

```bash
conda activate polap
polap-data-cflye delete-polap-github
polap-data-cflye sample-csv polap-data-v2.csv test
polap-data-cflye -y download-test-data
# run time: about 1 hour
polap-data-cflye local-batch Eucalyptus_pauciflora t off
polap-data-cflye -y install-getorganelle
# polap-data-cflye -y download-pmat
# polap-data-cflye -y install-pmat
polap-data-cflye sample-csv polap-data-v2.csv all on
# edit the CSV file if necessary
polap-data-cflye local-batch each
```

Now, go to step 10 of the next subsection to create tables and figures.

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
conda create -y --name polap polap=0.4.3.7.6
```

**5. Installation of Flye for disjointig filtering**:
Note that `Flye` with disjointig filtering feature is a slightly modidified version of the original `Flye`.
You activate the `polap` conda environment and setup `polap-fmlrc` environment

```bash
conda activate polap
conda install -y goshng::cflye
base_dir=$(dirname "$(command -v polap)") && \
  conda env create -f $base_dir/polap-conda-environment-fmlrc.yaml
```

**6. Polap assemble run with a test dataset**:
This tests the basic execution of the `polap` commmand.

```bash
wget -q https://github.com/goshng/polap/archive/refs/tags/0.4.3.7.6.zip
unzip -o -q 0.4.3.7.6.zip
cd polap-0.4.3.7.6/test
polap assemble --test
```

**7. Plastid genome assembly with _Eucalyptus pauciflora_ dataset**:
Your assembled plastid genome sequence will be `o/ptdna.0.fa`.

```bash
polap x-ncbi-fetch-sra --sra SRR7153095
polap x-ncbi-fetch-sra --sra SRR7161123
polap disassemble -l SRR7153095.fastq \
  -a SRR7161123_1.fastq \
  -b SRR7161123_2.fastq
```

**8. Check the accuracy of the plastid genome assembly**:
We use the Polap disassemble command with _Eucalyptus pauciflora_ dataset and check its similarity with its known plastid genome sequence
Your assembled plastid genome sequence will be `o/ptdna.ref.0.fa`. The text file named `o/0/mafft/pident.txt` has the percent identity between the assembled ptDNA and the knomn reference.

```bash
polap get-mtdna --plastid --species "Eucalyptus pauciflora"
cp o/00-bioproject/2-mtdna.fasta o/ptdna-reference.fa
polap disassemble \
  --disassemble-i 1 \
  --stages-include 3 \
  -l SRR7153095.fastq \
  -a SRR7161123_1.fastq \
  -b SRR7161123_2.fastq \
  --disassemble-align-reference \
  --disassemble-c o/ptdna-reference.fa

mkdir -p o/0/mafft
polap mafft-mtdna \
  -a o/ptdna-reference.fa \
  -b o/0/disassemble/2/pt.subsample-polishing.reference.aligned.1.fa \
  -o o/0/mafft
cat o/0/mafft/pident.txt
```

**9. Batch script that creates the results in the manuscript**:

```bash
polap-data-cflye -y install-getorganelle
# polap-data-cflye -y download-pmat
# polap-data-cflye -y install-pmat
polap-data-cflye example-data polap-data-v2.csv all on
polap-data-cflye local-batch each
```

**10. Tables in the manuscript**:
Tables in Markdown format will be generated and saved in the `man` directory after executing the following command.
You should download a precompiled binary version 0.8.1 of `Bandage` genome assembly graph visualization tool from [the official Bandage GitHub](https://github.com/rrwick/Bandage/releases).

```bash
polap-data-cflye -y install-bandage
# Install xelatex if necessary ...
# sudo apt-get install texlive texlive-latex-recommended texlive-xetex
# sudo apt-get install texlive-fonts-recommended texlive-fonts-extra texlive-lang-all
polap-data-cflye -y install-man
polap-data-cflye -y download-man
polap-batch-v2.sh
polap-data-cflye -y make-man
```

<!-- benchmark assembly graph figure -->
<!-- polap-data-cflye man figure-sheet some 2 -->
