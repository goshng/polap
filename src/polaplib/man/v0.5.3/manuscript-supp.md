<!-- \blandscape -->

\newpage

# Tables

<!--table1-->
<!-- bolap man table-benchmark some 2 -->
<!-- bolap man table-benchmark some 2 polap view -->

Table: Plastid genome assemblies for 23 plant species datasets using subsampled sequencing data with a maximum subsampling rate of 5% (Run Setting A).
All datasets were downsampled to 10x genome coverage, and Stage 1 included 10 subsampling steps (N).
Depending on the dataset, different maximum sampling rates (P) were used in Stage 1, and different replicate sizes (R) were applied in Stages 2 and 3.
NA represents no assemblies in the subsampling-based method and no comparison available.
{#tbl:benchmark-polap-some-0}

!include figures/table-benchmark-polap-some-0.md

<!-- \elandscape -->

<!--table2-->
<!-- bolap man table-polap-disassemble Eucalyptus_pauciflora 2 -->
<!-- bolap man table-polap-disassemble Eucalyptus_pauciflora 2 view -->

\newpage

# Figures

<!--figure1-->

![Workflow of the subsampling-based plastid genome assembly. The genome assembly procedure is applied repeatedly in Stages 1 and 2.](figures/polap2-figure1.pdf){#fig:mainfigure1 width=100%}

<!--figure2-->
<!-- bolap man figure-benchmark some 2 time -->
<!-- bolap man figure-benchmark some 2 time view -->

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
](figures/figure-time-some-0.pdf){#fig:time-some-0 width=100%}

<!--figure3-->
<!-- bolap man figure-benchmark some 2 memory -->
<!-- bolap man figure-benchmark some 2 memory view -->

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
](figures/figure-memory-some-0.pdf){#fig:memory-some-0 width=100%}

\setcounter{table}{0}
\setcounter{figure}{0}
\renewcommand{\thetable}{S\arabic{table}}
\renewcommand{\thefigure}{S\arabic{figure}}

\blandscape

# Supplementary Materials

<!--tableS1-->
<!-- bolap man table-benchmark some 2 -->
<!-- bolap man table-benchmark some 2 data view -->

Table: Sequencing data for the datasets, including species names and their corresponding taxonomic ranks studied. {#tbl:benchmark-data-some-0}

!include figures/table-benchmark-data-some-0.md

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
<!-- bolap man table-benchmark some 2 -->
<!-- bolap man table-benchmark some 2 computer view -->
<!-- bolap man table-benchmark some 2 hostname view -->

\newpage

Table: Computer setup for the 23 datasets.
{#tbl:benchmark-computer-some-0}

!include figures/table-benchmark-computer-some-0.md

<!--tableS3-->
<!-- bolap man table-benchmark some 1 -->
<!-- bolap man table-benchmark some 1 polap view -->

<!--figureS1-->
<!-- bolap man figure-alpha Eucalyptus_pauciflora -->
<!-- bolap man figure-alpha Eucalyptus_pauciflora view # Figure S1 -->

# Code

`Polap` (Plant Organelle Long-read Assembly Pipeline v0.5.3) is available at [http://github.com/goshng/polap](http://github.com/goshng/polap).
A guide of the homology-based plastome assembly is provided for use on a Linux system with an Internet connection.
It offers a step-by-step description of the procedures used for the results presented in the manuscript.

### Requirements

- **Operating System**: Linux with Ubuntu 24.04 LTS (not compatible with macOS or Windows)
- **Dependencies**: Requires [Bash](https://www.gnu.org/software/bash/) (>= 5.0) and **[Miniconda](https://docs.anaconda.com/miniconda/miniconda-install/)**

### Steps

To replicate the results presented in this manuscript on a Linux computer with `git` installed and an Internet connection, follow the steps below.

**1. Install Miniconda**:
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

Use the following script to download and install **[Miniconda](https://docs.anaconda.com/miniconda/miniconda-install/)**. It uses the [instructions](https://docs.anaconda.com/miniconda/#quick-command-line-install) for the Miniconda installation available at [https://docs.anaconda.com/miniconda/#quick-command-line-install](https://docs.anaconda.com/miniconda/#quick-command-line-install).
Use either `git` or `wget` (along with `unzip` and `mv`) to prepare a folder named `polap` at the current directory.
Replace, if necessary, the version number {{version}} with something recommended at [Polap's github website](http://github.com/goshng/polap).

```bash
cd
mkdir -p all/polap/read1
cd all/polap/read1
rm -rf polap
# option 1: use git command
git clone https://github.com/goshng/polap.git
# option 2: use wget, unzip, and mv
wget -q https://github.com/goshng/polap/archive/refs/tags/{{version}}.zip
unzip -o -q {{version}}.zip
mv polap-{{version}} polap

# Now, we have the polap folder.
bash polap/src/bolap -y install conda
```

After installing, close and reopen the terminal application to log out and back into the terminal.
Then, execute the followings to setup the conda channels for `polap`.

```bash
cd all/polap/read1
source ~/miniconda3/bin/activate
bash polap/src/bolap setup conda
```

**2. Install Bioconda packages**:
Install conda packages including `polap` and others in their own conda environments.
We execute `bolap -y install polap` if we want the latest version of `polap`.

```bash
# option 1: to install all tools necessary
bash polap/src/bolap -y install all={{version2}
# option 2: to install polap only
bash polap/src/bolap -y install polap={{version2}
bash polap/src/bolap setup polap
```

After the Bioconda package installation, log out and back into the terminal.
To use the subsampling-based plastid assembly Pipeline, go to Step 7.
Otherwise, go to the next step.

**3. Benchamrk to reproduce the tables and figures in the manuscript**:
We create the report in PDF format.

```bash
bolap -f benchmark
bolap -y man init
bolap man table-benchmark
bolap man figure-sheet # Requires latex
bolap man figure-benchmark some 0 time
bolap man figure-benchmark some 0 memory
bolap man pdf # Requires latex
```

Now, we have `manuscript.pdf` for the report.

<!-- benchmark assembly graph figure -->
<!-- bolap man figure-sheet some 2 -->
