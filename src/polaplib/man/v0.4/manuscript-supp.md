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
NAs at the columns of N, L, and C in Stage 3 are because of no assembly in the stage but just polishing of the draft genome sequences.
{#tbl:polap-disassemble-Eucalyptus_pauciflora-2}

!include figures/table-polap-disassemble-Eucalyptus_pauciflora-2.md

\newpage

# Figures

<!--figure1-->

![Workflow of the subsampling-based plastid genome assembly. The genome assembly procedure is applied repeatedly in Stages 1 and 2.](figures/polap2-figure1.pdf){#fig:mainfigure1 width=100%}

<!--figure2-->

![Benchmark of `GetOrganelle`, `ptGAUL`, `PMAT`, `TIPPo`, `Oatk` and the method (Run Setting A) presented here in terms of data processing time.
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

![Benchmark of `GetOrganelle`, `ptGAUL`, `PMAT`, `TIPPo`, `Oatk` and the method (Run Setting A) presented here in terms of peak memory usage.
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
NA represents no assemblies in the subsampling-based method and no comparison available.
{#tbl:benchmark-polap-some-1}

!include figures/table-benchmark-polap-some-1.md

<!-- \elandscape -->

<!--table:benchmark-time-->
<!-- polap-data-cflye man table-benchmark some 2 -->

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
NAs at the columns of N, L, and C in Stage 3 are because of no assembly in the stage but just polishing of the draft genome sequences.
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
NAs at the columns of N, L, and C in Stage 3 are because of no assembly in the stage but just polishing of the draft genome sequences.
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
Description:    Ubuntu 24.04.1 LTS
```

**2. Install Miniconda**:
Use the following script to download and install **[Miniconda](https://docs.anaconda.com/miniconda/miniconda-install/)**. It uses the [instructions](https://docs.anaconda.com/miniconda/#quick-command-line-install) for the installation.
Use either `git` or `wget` (along with `unzip` and `mv`) to prepare a folder named `polap` at the current directory.
Replace the version number 0.4.3.7.7 with something recommended at [Polap's github website](http://github.com/goshng/polap).

```bash
mkdir -p all/polap/cflye1
cd all/polap/cflye1
# option 1: use git command
git clone https://github.com/goshng/polap.git
# option 2: use wget, unzip, and mv
wget -O polap.zip -q https://github.com/goshng/polap/archive/refs/tags/0.4.3.7.7.zip
unzip -o -q polap.zip
mv polap-0.4.3.7.7 polap
# Now, we have the polap folder.
bash polap/src/polap-data-cflye -y install conda
```

After installing, close and reopen your terminal application to log out and back into the terminal.
Then, execute the followings to setup the conda channels for `polap`.

```bash
cd all/polap/cflye1
source ~/miniconda3/bin/activate
bash polap/src/polap-data-cflye setup conda
```

**3. Install Bioconda packages**:
Install conda packages incluing `polap` and others in their own conda environments.
Remove the version number `0.4.3.7.7` from `polap-data-cflye -y install all` if you want the latest version of `polap`.

```bash
bash polap/src/polap-data-cflye -y install all 0.4.3.7.7
bash polap/src/polap-data-cflye setup polap
bash polap/src/polap-data-cflye setup pmat
```

After the Bioconda package installation, log out and back into the terminal.

**4. Test the benchmark script that creates the results in the manuscript**:
Before executing the scripts for the results presented in the manuscript, we test the installation in three parts. Each part ensures that the actual analyses are performed as intended. First, we test the main data script that we used for the results presented in the manuscript. Second, we test Polap's command `assemble` using a simple test dataset. Third, we test Polap's command `disassemble` using the real dataset for _Eucalyptus pauciflora_.
This is the test of the first part.

```bash
cd all/polap/cflye1
source ~/miniconda3/bin/activate
conda activate polap
polap-data-cflye -y download test-data-cflye
polap-data-cflye benchmark Test_species
polap-data-cflye benchmark Taxon_genus
for i in 3 11 12 21 22; do
  polap-data-cflye run polap-disassemble Taxon_genus $i
done
# run time: about 3 hour on a 56-core Ubuntu computer
```

**5. Tables and figures in the manuscript using the test data**:
Tables in Markdown format will be generated and saved in the `man` directory after executing the following command.
You should download a precompiled binary version 0.8.1 of `Bandage` genome assembly graph visualization tool from [the official Bandage GitHub](https://github.com/rrwick/Bandage/releases).
To generate a report in PDF format, we need `latex` installed in the Ubuntu computer.
You may skip the installation of `latex`, but you won't be able to generate a PDF-format report.

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
polap-data-cflye man figure-benchmark test 2 time
polap-data-cflye man figure-benchmark test 2 memory
polap-data-cflye man figure-alpha Taxon_genus
polap-data-cflye man figure-delta Taxon_genus
polap-data-cflye man pdf test # Requires latex
```

Now, you have `manuscript-test.pdf` for the test report. Use the following to display the result if `latex` is not installed into the Ubuntu computer.

```bash
polap-data-cflye man table-benchmark test 2 computer view
polap-data-cflye man table-benchmark test 2 time view
polap-data-cflye man table-benchmark test 2 memory view
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

You should check if there is a file named `2-oga.gfa` at the curent directory.
One could extract a sequence file from `2-oga.gfa` and name it `mt.0.fasta` using `Bandage` software.
For testing purpose, we polish the extracted sequence.

```bash
polap prepare-polishing
polap polish
diff mt.0.fa mt.1.fa
```

If there is no screen output from the last `diff mt.0.fa mt.1.fa`, then the Polap's polishing step works fine. Move on to the testing of subsampling-based assembly.

**7. Plastid genome assembly with _Eucalyptus pauciflora_ dataset**:
We could download the sequencing data files from the NCBI SRA database.
Instead, we use the test dataset we downloaded earlier in Step 4.
Because the test dataset is too small, we downsample the data to the 100x of the genome size estimate and use 0.9 of the maximum subsampling rate.
We use the line-continuation backslash character to prevent the command from exceeding the page width; that is, the three lines together form a single command starting with `polap disassemble` below.

```bash
# polap x-ncbi-fetch-sra --sra SRR7153095
# polap x-ncbi-fetch-sra --sra SRR7161123
# polap disassemble -l SRR7153095.fastq \
#  -a SRR7161123_1.fastq \
#  -b SRR7161123_2.fastq
polap disassemble -l l.fastq -a s_1.fastq -b s_2.fastq \
  --downsample 100 --disassemble-p 90 \
  -o a
```

Your assembled plastid genome sequence will be `a/ptdna.0.fa`.

**8. Check the accuracy of the plastid genome assembly**:
We use the Polap disassemble command with _Eucalyptus pauciflora_ dataset and check its similarity with its known plastid genome sequence
Your assembled plastid genome sequence will be `o/ptdna.ref.0.fa`. The text file named `o/0/mafft/pident.txt` has the percent identity between the assembled ptDNA and the knomn reference.
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

```bash
cd all/polap/cflye1
source ~/miniconda3/bin/activate
conda activate polap
# option 1: use some for all datasets
polap-data-cflye benchmark some
# option 2: use a species folder name for a particular dataset
polap-data-cflye benchmark Eucalyptus_pauciflora
polap-data-cflye benchmark Spirodela_polyrhiza
for i in 3 {11..29}; do
  polap-data-cflye run polap-disassemble Eucalyptus_pauciflora $i
done
```

**10. Tables and figures in the manuscript**:
We create the report in PDF format.

```bash
rm -rf man
polap-data-cflye man init
for i in 0 1 2 4; do polap-data-cflye man table-benchmark some $i; done
for i in {0..4}; do
  polap-data-cflye man table-polap-disassemble Eucalyptus_pauciflora $i
done
polap-data-cflye man figure-sheet some 2 bandage
polap-data-cflye man figure-benchmark some 2 time
polap-data-cflye man figure-benchmark some 2 memory
polap-data-cflye man figure-alpha Eucalyptus_pauciflora
polap-data-cflye man figure-delta Eucalyptus_pauciflora
polap-data-cflye man pdf
```

Now, you have `manuscript.pdf` for the report.
Use the following to display the result if `latex` is not installed into the Ubuntu computer.

```bash
polap-data-cflye man table-benchmark some 2 computer view
polap-data-cflye man table-benchmark some 2 time view
polap-data-cflye man table-benchmark some 2 memory view
```

<!-- benchmark assembly graph figure -->
<!-- polap-data-cflye man figure-sheet some 2 -->
