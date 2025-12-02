<!-- \blandscape -->

\newpage

# Tables

<!--table1-->
<!-- polap-data-cflye man table-benchmark some 2 -->
<!-- polap-data-cflye man table-benchmark some 2 polap view -->

Table: Plastid genome assemblies for 23 plant species datasets using subsampled sequencing data with a maximum subsampling rate of 5% (Run Setting A).
All datasets were downsampled to 10x genome coverage, and Stage 1 included 10 subsampling steps (N).
Depending on the dataset, different maximum sampling rates (P) were used in Stage 1, and different replicate sizes (R) were applied in Stages 2 and 3.
NA represents no assemblies in the subsampling-based method and no comparison available.
{#tbl:benchmark-polap-some-2}

!include figures/table-benchmark-polap-some-2.md

<!-- \elandscape -->

<!--table2-->
<!-- polap-data-cflye man table-polap-disassemble Eucalyptus_pauciflora 2 -->
<!-- polap-data-cflye man table-polap-disassemble Eucalyptus_pauciflora 2 view -->

\newpage

Table: Three stages of subsampling-based plastid genome assembly for the _Eucalyptus pauciflora_ dataset with Run Setting A.
The configuration includes an increasing subsample size up to a maximum subsampling rate of 5%, a step size of 10 in Stage 1, 5 replicates in Stages 2 and 3, and a maximum memory limit of 16 GB.
Abbreviations are as follows:
iteration in each Stage (Index), subsampling rate (Rate) and read-coverage threshold (Alpha);
assembly metrics including the number of segments in the assembly (N),
the total length of these segments (L), and
the number of circular genome paths detected (C);
memory usage in GB (Memory) and
the draft plastid genome assembly length (Length).
Alpha at Stage 3 is the percent identity values between consecutive indices.
NAs at the columns of N, L, and C in Stage 3 are because of no assembly in the stage but just polishing of the draft genome sequences.
{#tbl:polap-disassemble-Eucalyptus-pauciflora-2}

!include figures/table-polap-disassemble-Eucalyptus_pauciflora-2.md

\newpage

# Figures

<!--figure2-->
<!-- polap-data-cflye man figure-benchmark some 2 time -->
<!-- polap-data-cflye man figure-benchmark some 2 time view -->

![Benchmark of `GetOrganelle`, `ptGAUL`, `Oatk`, `TIPPo`, `PMAT` and the method (Run Setting A) presented here in terms of data processing time.
See @fig:time-nextdenovo-some-2 for the same results, including the processing time of NextDenovo.
Abbreviations are as follows:
Genome size estimation by `JellyFish` (G),
Short-read polishing preparation (MSBWT),
Short-read polishing (FMLRC),
`Oatk` with `-c 30` (Oatk-30),
`Oatk` with `-c 20` (Oatk-20),
`TIPPo` with `-p onthq` (TIPPo-hq),
`TIPPo` with `-p ont` (TIPPo-ont),
`PMAT` with `-fc 0.1` (PMAT-0.1), and
`PMAT` with `-fc 1.0` (PMAT-1.0).
](figures/figure-time-some-2.pdf){#fig:time-some-2 width=100%}

<!--figure3-->
<!-- polap-data-cflye man figure-benchmark some 2 memory -->
<!-- polap-data-cflye man figure-benchmark some 2 memory view -->

![Benchmark of `GetOrganelle`, `ptGAUL`, `Oatk`, `TIPPo`, `PMAT` and the method (Run Setting A) presented here in terms of peak memory usage.
Abbreviations are as follows:
Genome size estimation by `JellyFish` (G),
Short-read polishing preparation (MSBWT),
Short-read polishing (FMLRC),
Long-read error-correction (NextDenovo),
`Oatk` with `-c 30` (Oatk-30),
`Oatk` with `-c 20` (Oatk-20),
`TIPPo` with `-p onthq` (TIPPo-hq),
`TIPPo` with `-p ont` (TIPPo-ont),
`PMAT` with `-fc 0.1` (PMAT-0.1), and
`PMAT` with `-fc 1.0` (PMAT-1.0).
](figures/figure-memory-some-2.pdf){#fig:memory-some-2 width=100%}

<!--figure1-->

![Workflow of the subsampling-based plastid genome assembly. The genome assembly procedure is applied repeatedly in Stages 1 and 2.](figures/polap2-figure1.pdf){#fig:mainfigure1 width=100%}

\setcounter{table}{0}
\setcounter{figure}{0}
\renewcommand{\thetable}{S\arabic{table}}
\renewcommand{\thefigure}{S\arabic{figure}}

\blandscape

# Supplementary Materials

<!--tableS3-->
<!-- polap-data-cflye man table-benchmark some 1 -->
<!-- polap-data-cflye man table-benchmark some 1 polap view -->

<!-- \newpage -->

<!-- \blandscape -->

Table: Replicate of plastid genome assemblies for 23 plant species datasets using subsampled sequencing data with a maximum subsampling rate of 5% (Run Setting A).
All datasets were downsampled to 10x genome coverage, and Stage 1 included 10 subsampling steps (N).
Depending on the dataset, different maximum sampling rates (P) were used in Stage 1, and different replicate sizes (R) were applied in Stages 2 and 3.
NA represents no assemblies in the subsampling-based method and no comparison available.
{#tbl:benchmark-polap-some-1}

!include figures/table-benchmark-polap-some-1.md

<!-- \elandscape -->

<!--table:benchmark-time-->
<!-- polap-data-cflye man table-benchmark some 2 -->
<!-- polap-data-cflye man table-benchmark some 2 time view # Table for Figure 2 and S3 -->

<!-- \newpage -->
<!---->
<!-- \blandscape -->
<!---->
<!-- Table: Benchmark of `GetOrganelle`, `ptGAUL`, `PMAT`, `TIPPo`, `Oatk` and the method (Run Setting A) presented here in terms of data processing time. -->
<!-- NA at the column of `NextDenovo` represents no error-corrected long-read results, resulting in no assemblies in the correction-then-assembly pipelines including `PMAT`, `TIPPo`, and `Oatk`. -->
<!-- Abbreviations are as follows: -->
<!-- `GetOrganelle` (GO), -->
<!-- `ptGAUL` (ptG), -->
<!-- `NextDenovo` (ND), -->
<!-- `PMAT` with `-fc 0.1` (P0.1), -->
<!-- `PMAT` with `-fc 1.0` (P1.0), -->
<!-- `TIPPo` with `-p onthq` (Thq), -->
<!-- `TIPPo` with `-p ont` (Tont), -->
<!-- `Oatk` with `-c 30` (O30), and -->
<!-- `Oatk` with `-c 20` (O20). -->
<!-- {#tbl:benchmark-time-2} -->
<!---->
<!-- !include figures/table-benchmark-time-some-2.md -->
<!---->
<!-- \elandscape -->

<!--table:benchmark-memory-->
<!-- polap-data-cflye man table-benchmark some 2 -->
<!-- polap-data-cflye man table-benchmark some 2 memory view # Table for Figure 3 -->

<!-- \newpage -->
<!---->
<!-- \blandscape -->
<!---->
<!-- Table: Benchmark of `GetOrganelle`, `ptGAUL`, `PMAT`, `TIPPo`, `Oatk` and the method (Run Setting A) in terms of peak memory. -->
<!-- NA at the column of `NextDenovo` represents no error-corrected long-read results, resulting in no assemblies in the correction-then-assembly pipelines including `PMAT`, `TIPPo`, and `Oatk`. -->
<!-- Abbreviations are as follows: -->
<!-- `GetOrganelle` (GO), -->
<!-- `ptGAUL` (ptG), -->
<!-- `NextDenovo` (ND), -->
<!-- `PMAT` with `-fc 0.1` (P0.1), -->
<!-- `PMAT` with `-fc 1.0` (P1.0), -->
<!-- `TIPPo` with `-p onthq` (Thq), -->
<!-- `TIPPo` with `-p ont` (Tont), -->
<!-- `Oatk` with `-c 30` (O30), and -->
<!-- `Oatk` with `-c 20` (O20). -->
<!-- {#tbl:benchmark-memory-2} -->
<!---->
<!-- !include figures/table-benchmark-memory-some-2.md -->
<!---->
<!-- \elandscape -->

<!--tableS4-->
<!-- polap-data-cflye man table-benchmark some 0 -->
<!-- polap-data-cflye man table-benchmark some 0 polap view -->

\newpage

Table: Plastid genome assemblies for 23 plant species datasets using subsampled sequencing data with a maximum subsampling rate of 10% (Run Setting B).
All datasets were downsampled to 10x genome coverage, and Stage 1 included 10 subsampling steps (N).
Abbreviations are as follows:
maximum sampling rates (P) in Stage 1, and replicate sizes (R) in Stages 2 and 3.
NA represents no assemblies in the subsampling-based method and no comparison available.
{#tbl:benchmark-polap-some-0}

!include figures/table-benchmark-polap-some-0.md

<!--tableS5-->
<!-- polap-data-cflye man table-benchmark some 4 -->
<!-- polap-data-cflye man table-benchmark some 4 polap view -->

\newpage

Table: Plastid genome assemblies for 23 plant species datasets using subsampled sequencing data with a maximum subsampling rate of 1%.
All datasets were downsampled to 10x genome coverage, and Stage 1 included 10 subsampling steps (N).
Abbreviations are as follows:
maximum sampling rates (P) in Stage 1, and replicate sizes (R) in Stages 2 and 3.
NA represents no assemblies in the subsampling-based method and no comparison available.
{#tbl:benchmark-polap-some-4}

!include figures/table-benchmark-polap-some-4.md

<!--tableS6-->
<!-- polap-data-cflye man table-polap-disassemble Eucalyptus_pauciflora 0 -->
<!-- polap-data-cflye man table-polap-disassemble Eucalyptus_pauciflora 0 view -->

\newpage

Table: Three stages of subsampling-based plastid genome assembly for the _Eucalyptus pauciflora_ dataset with Run Setting B.
The configuration includes an increasing subsample size up to a maximum subsampling rate of 10%, a step size of 10 in Stage 1, 5 replicates in Stages 2 and 3, and a maximum memory limit of 16 GB.
Abbreviations are as follows:
iteration in each Stage (Index), subsampling rate (Rate) and read-coverage threshold (Alpha);
assembly metrics including the number of segments in the assembly (N),
the total length of these segments (L), and
the number of circular genome paths detected (C);
memory usage in GB (Memory) and
the draft plastid genome assembly length (Length).
Alpha at Stage 3 is the percent identity values between consecutive indices.
NAs at the columns of N, L, and C in Stage 3 are because of no assembly in the stage but just polishing of the draft genome sequences.
{#tbl:polap-disassemble-Eucalyptus-pauciflora-0}

!include figures/table-polap-disassemble-Eucalyptus_pauciflora-0.md

<!--tableS7-->
<!-- polap-data-cflye man table-polap-disassemble Eucalyptus_pauciflora 3 -->
<!-- polap-data-cflye man table-polap-disassemble Eucalyptus_pauciflora 3 view -->

\newpage

Table: Three stages of subsampling-based plastid genome assembly for the _Eucalyptus pauciflora_ dataset with Run Setting C.
The configuration includes an increasing subsample size up to a maximum subsampling rate of 10%, a step size of 50 in Stage 1, 5 replicates in Stages 2 and 3, and a maximum memory limit of 16 GB.
Abbreviations are as follows:
iteration in each Stage (Index), subsampling rate (Rate) and read-coverage threshold (Alpha);
assembly metrics including the number of segments in the assembly (N),
the total length of these segments (L), and
the number of circular genome paths detected (C);
memory usage in GB (Memory) and
the draft plastid genome assembly length (Length).
Alpha at Stage 3 is the percent identity values between consecutive indices.
NAs at the columns of N, L, and C in Stage 3 are because of no assembly in the stage but just polishing of the draft genome sequences.
{#tbl:polap-disassemble-Eucalyptus-pauciflora-3}

!include figures/table-polap-disassemble-Eucalyptus_pauciflora-3.md

\newpage

<!--tableS1-->
<!-- polap-data-cflye man table-benchmark some 2 -->
<!-- polap-data-cflye man table-benchmark some 2 data view -->

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

<!--tableS2-->
<!--tableS:computer-->
<!-- polap-data-cflye man table-benchmark some 2 -->
<!-- polap-data-cflye man table-benchmark some 2 computer view -->
<!-- polap-data-cflye man table-benchmark some 2 hostname view -->

\newpage

Table: Computer setup for the 23 datasets.
{#tbl:benchmark-computer}

!include figures/table-benchmark-computer-some-2.md

<!--figureS3-->
<!-- polap-data-cflye man figure-benchmark some 2 time-nextdenovo -->
<!-- polap-data-cflye man figure-benchmark some 2 time-nextdenovo view # Figure S3 -->

![Benchmark of `GetOrganelle`, `ptGAUL`, `Oatk`, `TIPPo`, `PMAT` and the method (Run Setting A) presented here in terms of data processing time.
See @fig:time-some-2 for the same results, excluding the processing time of `NextDenovo`.
Abbreviations are as follows:
Genome size estimation by `JellyFish` (G),
Short-read polishing preparation (MSBWT),
Short-read polishing (FMLRC),
Long-read error-correction (NextDenovo),
`Oatk` with `-c 30` (Oatk-30),
`Oatk` with `-c 20` (Oatk-20),
`TIPPo` with `-p onthq` (TIPPo-hq),
`TIPPo` with `-p ont` (TIPPo-ont),
`PMAT` with `-fc 0.1` (PMAT-0.1), and
`PMAT` with `-fc 1.0` (PMAT-1.0).
](figures/figure-time-nextdenovo-some-2.pdf){#fig:time-nextdenovo-some-2 width=100%}

<!--figureS1-->
<!-- polap-data-cflye man figure-alpha Eucalyptus_pauciflora -->
<!-- polap-data-cflye man figure-alpha Eucalyptus_pauciflora view # Figure S1 -->

\newpage

![Line plot of read-coverage thresholds versus subsample size index in Stage 1 over different initial values of the read-coverage threshold of the subsampling-based assemblies for the _Eucalyptus pauciflora_ dataset.](figures/alpha0.pdf){#fig:alpha0 width=100%}

<!--figureS2-->
<!-- polap-data-cflye man figure-delta Eucalyptus_pauciflora -->
<!-- polap-data-cflye man figure-delta Eucalyptus_pauciflora view # Figure S2 -->

\newpage

![Line plot of read-coverage thresholds versus subsample size index in Stage 1 over different different increment sizes of the subsampling-based assemblies for the _Eucalyptus pauciflora_ dataset.](figures/delta.pdf){#fig:delta width=100%}

<!-- Supporting Material - plastid genome assemblies using the six pipelines -->

\newpage

# Code

`Polap` (Plant Organelle Long-read Assembly Pipeline v0.4.3) is available at [http://github.com/goshng/polap](http://github.com/goshng/polap).
A guide of the subsampling-based plastome assembly is provided for use on a Linux system with an Internet connection.
It offers a step-by-step description of the procedures used for the results presented in the manuscript.

### Requirements

- **Operating System**: Linux with Ubuntu 24.04 LTS (not compatible with macOS or Windows)
- **Dependencies**: Requires [Bash](https://www.gnu.org/software/bash/) (>= 5.0) and **[Miniconda](https://docs.anaconda.com/miniconda/miniconda-install/)**

### Steps

To replicate the results presented in this manuscript on a Linux computer with `git` installed and an Internet connection, follow the steps below.

**1. Open a new terminal**:
Open a new terminal in a Linux computer, such as one with Ubuntu. Make sure that the `Bash` version is at least 5.0.

```bash
bash --version
lsb_release -a
```

The expected screen output includes the following two lines.

```txt
GNU bash, version 5.2.21(1)-release (x86_64-pc-linux-gnu)
Description: Ubuntu 24.04 LTS
```

**2. Install Miniconda**:
Use the following script to download and install **[Miniconda](https://docs.anaconda.com/miniconda/miniconda-install/)**. It uses the [instructions](https://docs.anaconda.com/miniconda/#quick-command-line-install) for the Miniconda installation available [https://docs.anaconda.com/miniconda/#quick-command-line-install](https://docs.anaconda.com/miniconda/#quick-command-line-install).
Use either `git` or `wget` (along with `unzip` and `mv`) to prepare a folder named `polap` at the current directory.
Replace, if necessary, the version number 0.4.3.7.9 with something recommended at [Polap's github website](http://github.com/goshng/polap).

```bash
cd
mkdir -p all/polap/cflye1
cd all/polap/cflye1
rm -rf polap
# option 1: use git command
git clone https://github.com/goshng/polap.git
# option 2: use wget, unzip, and mv
wget -q https://github.com/goshng/polap/archive/refs/tags/0.4.3.7.9.zip
unzip -o -q 0.4.3.7.9.zip
mv polap-0.4.3.7.9 polap
# Now, we have the polap folder.
bash polap/src/polap-data-cflye -y install conda
```

After installing, close and reopen the terminal application to log out and back into the terminal.
Then, execute the followings to setup the conda channels for `polap`.

```bash
cd all/polap/cflye1
source ~/miniconda3/bin/activate
bash polap/src/polap-data-cflye setup conda
```

**3. Install Bioconda packages**:
Install conda packages including `polap` and others in their own conda environments.
We execute `polap-data-cflye -y install polap` if we want the latest version of `polap`.

```bash
# option 1: to install all tools necessary
bash polap/src/polap-data-cflye -y install all=0.4.3
# option 2: to install polap only
bash polap/src/polap-data-cflye -y install polap=0.4.3
bash polap/src/polap-data-cflye setup polap
bash polap/src/polap-data-cflye setup pmat
```

After the Bioconda package installation, log out and back into the terminal.
To use the subsampling-based plastid assembly Pipeline, go to Step 7.
Otherwise, go to the next step.

**4. Test the benchmark script that creates the results in the manuscript**:
Before executing the scripts for the results presented in the manuscript, we test the installation in three parts. Each part ensures that the actual analyses are performed as intended. First, we test the main data script that we used for the results presented in the manuscript. Second, we test Polap's subcommand `assemble` using a simple test dataset. Third, we test Polap's subcommand `disassemble` using the real dataset for _Eucalyptus pauciflora_.
This is the test of the first part.

```bash
cd
cd all/polap/cflye1
source ~/miniconda3/bin/activate
conda activate polap
polap-data-cflye -y download test-data-cflye
polap-data-cflye benchmark Test_species
polap-data-cflye benchmark Taxon_genus
for i in 3 11 12 21 22; do
  polap-data-cflye run polap-disassemble Taxon_genus $i
done
# run time: about 12 hour on a 56-core Ubuntu computer
```

**5. Tables and figures in the manuscript using the test data**:
Tables in Markdown format will be generated and saved in the `man` directory after executing the following command.
Download a precompiled binary version 0.8.1 of `Bandage` genome assembly graph visualization tool from [the official Bandage GitHub](https://github.com/rrwick/Bandage/releases).
To generate a report in PDF format, we need `latex` installed in the Ubuntu computer.
We may skip the installation of `latex`, but we won't be able to generate a PDF-format report.

```bash
polap-data-cflye -y install bandage
polap-data-cflye setup bandage
polap-data-cflye install latex
polap-data-cflye man init
for i in 0 1 2 4; do polap-data-cflye man table-benchmark test $i; done
for i in {0..4}; do
  polap-data-cflye man table-polap-disassemble Taxon_genus $i
done
polap-data-cflye man figure-sheet test 2 bandage # Requires latex
polap-data-cflye man figure-sheet-pmat test # Requires latex
polap-data-cflye man figure-sheet-oatk test # Requires latex
polap-data-cflye man figure-sheet-tippo test # Requires latex
polap-data-cflye man figure-benchmark test 2 time
polap-data-cflye man figure-benchmark test 2 memory
polap-data-cflye man figure-benchmark test 2 time-nextdenovo
polap-data-cflye man figure-alpha Taxon_genus
polap-data-cflye man figure-delta Taxon_genus
polap-data-cflye man pdf test # Requires latex
```

Now, we have `manuscript-test.pdf` for the test report. Use the following to display the result if `latex` is not installed into the Ubuntu computer.

```bash
polap-data-cflye man table-benchmark test 2 computer view
polap-data-cflye man table-benchmark test 2 time view
polap-data-cflye man table-benchmark test 2 memory view
polap-data-cflye man table-polap-disassemble Taxon_genus 2 view
```

**6. Polap assemble run with a test dataset**:
This tests the basic execution of the `polap` commmand.

```bash
cd
cd all/polap/cflye1
source ~/miniconda3/bin/activate
conda activate polap
cd polap/test
polap assemble --test
```

There is a file named `2-oga.gfa` at the curent directory.
One could extract a sequence file from `2-oga.gfa` and name it `mt.0.fasta` using `Bandage` software.
For testing purpose, we polish the extracted sequence.

```bash
polap prepare-polishing
polap polish
diff mt.0.fa mt.1.fa
```

If there is no screen output from the last command: `diff mt.0.fa mt.1.fa`, then the Polap's polishing step works fine. Move on to the testing of subsampling-based assembly.

**7. Plastid genome assembly with _Eucalyptus pauciflora_ dataset**:
We could download the sequencing data files from the NCBI SRA database.
Instead, we use the test dataset we downloaded earlier in Step 4.
Because the test dataset is too small, we downsample the data to the 100x of the genome size estimate and use 0.9 of the maximum subsampling rate.
We use the line-continuation backslash character to prevent the command from exceeding the page width; that is, the three lines together form a single command starting with `polap disassemble` below.

```bash
cd
cd all/polap/cflye1
# option 1: download ONT long- and Illumina short-read datasets
# polap x-ncbi-fetch-sra --sra SRR7153095
# polap x-ncbi-fetch-sra --sra SRR7161123
# polap disassemble -l SRR7153095.fastq \
#  -a SRR7161123_1.fastq \
#  -b SRR7161123_2.fastq
polap disassemble -l l.fastq -a s_1.fastq -b s_2.fastq \
  --downsample 100 --disassemble-p 90 \
  -o a
```

The assembled plastid genome sequence should be `a/ptdna.0.fa`.

**8. Check the accuracy of the plastid genome assembly**:
We use the Polap disassemble command with _Eucalyptus pauciflora_ dataset and check its similarity with its known plastid genome sequence
The assembled plastid genome sequence will be `a/ptdna.ref.0.fa`. The text file named `a/0/mafft/pident.txt` has the percent identity between the assembled ptDNA and the knomn reference.
One could skip this step because it is not the main part of the subsampling-based method but a part for benchmarking.

```bash
polap get-mtdna --plastid --species "Eucalyptus pauciflora" -o a
cp a/00-bioproject/2-mtdna.fasta a/ptdna-reference.fa
polap disassemble \
  --disassemble-i 1 \
  --stages-include 3 \
  -l l.fastq \
  -a s_1.fastq \
  -b s_2.fastq \
  --disassemble-align-reference \
  --disassemble-c a/ptdna-reference.fa \
  --downsample 100 --disassemble-p 90 -o a

mkdir -p a/0/mafft
polap mafft-mtdna \
  -a a/ptdna-reference.fa \
  -b a/0/disassemble/1/pt.subsample-polishing.reference.aligned.1.fa \
  -o a/0/mafft
cat a/0/mafft/pident.txt
```

**9. Batch script that creates the results in the manuscript**:
We are ready to perform the subsampling-based plastid genome assemblies for the datasets.
One could do the analyses for either all datasets or each dataset depending on available computing resources.
One one can download the results available at Figshare: [https://figshare.com/s/ec1cb394870c7727a2d4](https://figshare.com/s/ec1cb394870c7727a2d4) or using the `wget` command to download them as a single zip file as shown below.

```bash
cd
cd all/polap/cflye1
source ~/miniconda3/bin/activate
conda activate polap
# option 1: use some for all datasets: this could take very long time
polap-data-cflye benchmark some
# option 2: use a species folder name for a particular dataset
polap-data-cflye benchmark Eucalyptus_pauciflora
polap-data-cflye benchmark Spirodela_polyrhiza
for i in 3 {11..29}; do
  polap-data-cflye run polap-disassemble Eucalyptus_pauciflora $i
done
# option 3: download the archived results
# https://figshare.com/ndownloader/articles/28740323?private_link=ec1cb394870c7727a2d4
wget -O 28740323.zip https://figshare.com/ndownloader/articles/28740323?private_link=ec1cb394870c7727a2d4
unzip 28740323.zip
for i in *-a.tar.gz; do polap-data-cflye recover ${i%-a.tar.gz}; done
```

**10. Tables and figures in the manuscript**:
We create the report in PDF format.

```bash
polap-data-cflye -y man init
for i in 0 1 2 4; do polap-data-cflye man table-benchmark some $i; done
for i in {0..4}; do
  polap-data-cflye man table-polap-disassemble Eucalyptus_pauciflora $i
done
polap-data-cflye man figure-sheet some 2 bandage # Requires latex
polap-data-cflye man figure-sheet-pmat some # Requires latex
polap-data-cflye man figure-sheet-oatk some # Requires latex
polap-data-cflye man figure-sheet-tippo some # Requires latex
polap-data-cflye man figure-benchmark some 2 time
polap-data-cflye man figure-benchmark some 2 memory
polap-data-cflye man figure-benchmark some 2 time-nextdenovo
polap-data-cflye man figure-alpha Eucalyptus_pauciflora
polap-data-cflye man figure-delta Eucalyptus_pauciflora
polap-data-cflye man pdf # Requires latex
```

Now, we have `manuscript.pdf` for the report.
Use the following to display the result if `latex` is not installed into the Ubuntu computer.

```bash
polap-data-cflye man table-benchmark some 2 computer view # Table S2
polap-data-cflye man table-benchmark some 2 data view # Table S1
polap-data-cflye man table-benchmark some 2 time view # Table for Figure 2 and S3
polap-data-cflye man table-benchmark some 2 memory view # Table for Figure 3
polap-data-cflye man table-benchmark some 2 polap view # Table 1
polap-data-cflye man table-benchmark some 1 polap view # Table S3
polap-data-cflye man table-benchmark some 0 polap view # Table S4
polap-data-cflye man table-benchmark some 4 polap view # Table S5
polap-data-cflye man table-polap-disassemble Eucalyptus_pauciflora 2 view # Table 2
polap-data-cflye man table-polap-disassemble Eucalyptus_pauciflora 0 view # Table S6
polap-data-cflye man table-polap-disassemble Eucalyptus_pauciflora 3 view # Table S7
polap-data-cflye man figure-benchmark some 2 time view # Figure 2
polap-data-cflye man figure-benchmark some 2 memory view # Figure 3
polap-data-cflye man figure-alpha Eucalyptus_pauciflora view # Figure S1
polap-data-cflye man figure-delta Eucalyptus_pauciflora view # Figure S2
polap-data-cflye man figure-benchmark some 2 time-nextdenovo view # Figure S3
polap-data-cflye man figure-sheet some 2 no-bandage # Genome assembly graphs
```

<!-- benchmark assembly graph figure -->
<!-- polap-data-cflye man figure-sheet some 2 -->
