# Description

Polap is a meta-assembler pipeline that treats
[Flye](https://github.com/fenderglass/Flye) as a blackbox tool.
Therefore, it can be described as a
[Flye](https://github.com/fenderglass/Flye) helper pipeline with
additional features of organellar gene annotation, which aids in the
selection of organelle-origin contigs obtained from the initial
[Flye](https://github.com/fenderglass/Flye) whole-genome assembly. The
name "polap" was chosen as it was inspired by the chloroplast genome
assembly pipeline known as [ptGAUL](https://github.com/Bean061/ptgaul),
developed by Wenbin Zhou. Polap utilizes the long read genome assembler,
[Flye](https://github.com/fenderglass/Flye), in combination with the
sequence polishing tool, [FMLRC](https://github.com/holtjma/fmlrc),
developed by Matt Holt. It also leverages the sequencing depths obtained
from a flye genome assembly and organelle gene annotation to identify
contigs that belong to the user's target mitochondrial genome.

Polap is a project that originated from the idea of
[ptGAUL](https://github.com/Bean061/ptgaul) and aims to investigate how
to assemble plant mitochondrial DNA (mtDNA) using long reads. The
process begins by performing a
[Flye](https://github.com/fenderglass/Flye) whole-genome assembly
without a reference genome to generate contigs that are then subjected
to [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) to confirm the
presence of mtDNA and ptDNA genes. Depth and Copy Number are also
determined after [Flye](https://github.com/fenderglass/Flye) assembly.
Using this information and the mtDNA and ptDNA genes, the pipeline helps
make educated guess about which contigs originate from mtDNA.

Contigs that are suspected to originate from mtDNA are collected, and
long reads that map well on to these contigs are used for mtDNA
long-read assembly. Once mtDNA assembly is complete, the FMLRC polishing
pipeline is employed to finalize the mtDNA.

# Requirements

Please, note that polap runs on only a Linux OS with
[Miniconda](https://docs.anaconda.com/miniconda/#quick-command-line-install).
I do not know how I format the recipe for
[Bioconda's polap package](https://anaconda.org/bioconda/polap), and
it looks like that polap could be executed in other platforms.
But, it does not run on either macOS or Windows computer.
Although [Bioconda's polap package](https://anaconda.org/bioconda/polap) can be
installed into either of them, they would fail to execute.

It requires [BASH](https://www.gnu.org/software/bash/) (bash\>=3.0).

    bash --version

You must install a conda package manager like
[Miniconda](https://docs.anaconda.com/miniconda/#quick-command-line-install)
on your favorite Linux Operating System in a command-line interface.

    $ curl -OL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    $ bash Miniconda3-latest-Linux-x86_64.sh

Your terminal should look like this after installing
[Miniconda](https://docs.anaconda.com/miniconda/#quick-command-line-install),
exiting, and relogin back to your Linux terminal.
The prompt with (base) indicates your successful
[Miniconda](https://docs.anaconda.com/miniconda/#quick-command-line-install)
installation.

    (base) $

Now, configure your conda with the following:

    (base) $ conda update -y -n base conda
    (base) $ conda config --prepend channels bioconda
    (base) $ conda config --prepend channels conda-forge

# Installing Polap

## Using Biocoda polap package

You could install
[Bioconda's polap package](https://anaconda.org/bioconda/polap)
in a new conda environment, and then
run [Polap](https://github.com/goshng/polap) with a test dataset.

    (base) $ conda create --name polap bioconda::polap
    (base) $ conda activate polap
    (polap) $ git clone https://github.com/goshng/polap.git
    (polap) $ cd polap/test
    (polap) $ polap assemble --test

Because FMLRC's conda environment was incompatible with that of Polap's,
we need to create one for Polap's FMLRC conda environment.
Continued from the previous command lines,

    (polap) $ conda deactivate
    (base) $ conda env create -f ../src/polap-conda-environment-fmlrc.yaml
    (base) $ conda activate polap-fmlrc
    (polap-fmlrc) $ polap prepare-polishing
    (polap-fmlrc) $ polap polish

If you see a screen output ending with the following:

    NEXT: src/polap.sh prepare-polishing [-a s1.fq] [-b s2.fq]

You may be ready to use [Polap](https://github.com/goshng/polap). But,
make sure that you do not have error messages from
[Polap](https://github.com/goshng/polap).

## Using github source

Download the source code of Polap available at
[Polap](https://github.com/goshng/polap)'s github website:

    (base) $ git clone https://github.com/goshng/polap.git

Setup the conda environments for
[Polap](https://github.com/goshng/polap).

    (base) $ conda env create -f src/polap-conda-environment.yaml

Setup the conda environment for
[FMLRC](https://github.com/holtjma/fmlrc).

    (base) $ conda env create -f src/polap-conda-environment-fmlrc.yaml

Activate polap conda environment.

    (base) $ conda activate polap

Run Polap with a test dataset.

    (polap) $ cd test
    (polap) $ ../src/polap.sh assemble --test

# Using Polap

Once you passed the test dataset run with polap,
you could execute polap two different ways depending on how you
installed it.
If you installed polap using
[Bioconda's polap package](https://anaconda.org/bioconda/polap),
you could simply type in `polap` at the prompt as long as you
have the `polap` conda environment activated.

    (polap) $ polap

Note that you have `(polap)` at your prompt.
If you installed polap using the github source, then
it would be easier to `git-clone`
[Polap](https://github.com/goshng/polap) at your data folder.
You still must have the `polap` conda environment activated.

    (polap) $ git clone https://github.com/goshng/polap.git
    (polap) $ src/polap.sh

Let me simply use `polap` command by assuming that you installed
polap using
[Bioconda's polap package](https://anaconda.org/bioconda/polap).
You would replace `polap` with the `git-clone` and `src/polap.sh`
commands. If you are an experienced command line user,
this introduction would be easy. If not, try to follow this
README.

Now, you create a folder and place
an input dataset in it: a single Oxford Nanopore long-read fastq file,
and two short-read fastq files.
The `test` folder from the `git-clone` source folder has a set of
such fastq files.
The case of a single short-read file
is not tested, yet.
If no _input-files_ are specified by options with fastq input file paths,
three input files must be at the current directory where you execute
the polap script:
`l.fq`, `s1.fq`, and `s2.fq`.
The three fastq files must be extracted if they are compressed.
You could use options such as
`-l`, `-a`, and `-b` to set your fastq files.
Please, execute polap without any option to see a help message.

Now that you have the three fastq files, the first polap step is
to run a whole-genome assembly using
[Flye](https://github.com/fenderglass/Flye)
with menu `assemble1`. The first argument of the command `polap`
is called menu. I do not know other better name for the first
argument. A menu has its help message with a second argument
of `help`.
You would see options for the menu `assemble1` with the following
command:

    (polap) $ polap assemble1 help

Assuming that you have the three fastq files:
`l.fq`, `s1.fq`, and `s2.fq`, you could execute a whole-genome
assembly step with the following command:

    (polap) $ polap assemble1

It would take some times depending on your dataset.
Please, be patient.
After the execution of the whole-genome assembly using
[Flye](https://github.com/fenderglass/Flye),
you should be able to open `o/0/30-contigger/graph_final.gfa` using
[Bandage](https://rrwick.github.io/Bandage/) to view the genome assembly
graph.
The graph would display contigs, some of which originate from
mitochondrial or plastid DNAs and most from nuclear DNA molecules.
You could look through the graph to locate contigs from
plant mitochondrial DNAs. If you could do so, then you could use
the contigs for plant mitochondrial DNA candidates.
If not, you continue the `polap` pipeline
by annotating your initial whole-genome assembly.

    (polap) $ polap annotate

This takes some time as well. So, be patient.
After the execution of the `polap` annotation step,
you should be able to see contigs with organelle gene annotations in
`o/0/contig-annotation-table.txt`.
There are two more similar text files:
`o/0/assembly_info_organelle_annotation_count-all.txt` and
`o/0/assembly_info_organelle_annotation_count.txt`.
After inspecting the three files,
choose contigs with more MT genes
than PT genes in the assembly graph.
Admittedly, this choice of contigs is subjective step.
I might need another section to describe what one could do in this
contig selection step. Let's assume that you have contig condidates
that potentially originate from plant mitochondrial DNAs.
The three text file above are a modified version of
`o/0/30-contigger/contigs_stats.txt` output file from
[Flye](https://github.com/fenderglass/Flye).
The file is a table with rows for contigs assembled from the whole-genome
assembly. Each contig may consist of one or more edge sequences.
The [Bandage](https://rrwick.github.io/Bandage/) graph shows
edge sequences rather than contig sequences.
So choose candidate contigs using
the [Bandage](https://rrwick.github.io/Bandage/) graph
and the three annotation table text file.

Now, you need to prepare a text file,
`o/0/mt.contig.name-1`, in your folder. Add one edge sequence name per
line. Note that the edge sequence names should start with `edge_` not
`contig_`. An example file for a MT contig name looks like this:

    edge_1
    edge_2
    edge_3

For more information on how you could prepare a MT contig file, see [MT
contig name](#mt-contig-name) below. Now, you are ready to run an
organelle-genome assembly.

    (polap) $ polap assemble2

Wait for a while.
This second genome assembly using
[Flye](https://github.com/fenderglass/Flye)
produces an assembly graph file called
`o/1/30-contigger/graph_final.gfa`.
Use
[Bandage](https://rrwick.github.io/Bandage/) to view the organelle-genome
assembly graph.
We hope that this organelle-genome assembly graph is simpler than
that of the whole-genome assembly.
If so, you
would want to finish your organelle-genome assembly upto
the polishing stage of
[Flye](https://github.com/fenderglass/Flye).

    (polap) $ polap flye-polishing

Now, extract a mitochondrial genome sequence by opening
`o/1/assembly_graph.gfa` using
[Bandage](https://rrwick.github.io/Bandage/). Save the sequence with a
file name called `mt.0.fasta`.
This is the long-read mitochondrial genome assembly
in fasta format.

We use the approach
[ptGAUL](https://github.com/Bean061/ptgaul), which uses
[FMLRC](https://github.com/holtjma/fmlrc) to polish
long-read assembly sequences using short-read data.
The next step is just a wrapper of the script
recommended by
[ptGAUL](https://github.com/Bean061/ptgaul).

    (polap) $ conda activate polap-fmlrc
    (polap-fmlrc) $ src/polap.sh prepare-polishing
    (polap-fmlrc) $ src/polap.sh polish

Your final mitochondrial genome sequence is `mt.1.fa`.

# Uninstalling Polap

    (polap) $ conda deactivate
    (base) $ conda remove -n polap --all
    (base) $ conda remove -n polap-fmlrc --all

# MT contig name

An MT contig name file should be edited to list contig names.
You select contigs that you suspect might have been from mitochondrial DNAs.
You may use the following for the three features of contigs.

1. genome assembly graph,
2. presence of organelle gene,
3. copy numbers of contigs.

Note that contigs are named by edge suffixed with a number.
You should use [Bandage](https://rrwick.github.io/Bandage/)
to view whole-genome assembly graph: o/0/30-contigger/graph_final.gfa.
You could see organelle annotation and copy numbers in

    $ column -t o/0/30-contigger/contigs_stats.txt

You could watch [YouTube](https://youtu.be/29OaiFCqAzI) (in Korean) to see how you could
select contigs that might have been from mitochondrial DNAs.
The following is an example of the Flye assembly table with organelle gene counts.
In the table, the column copy is the copy number for a contig,
MT is the number of mitochondrial genes, and
PT is the number of plastid genes.
Edge is edeg numbers that constitute the contig at the corresponding line.

    $ column -t o/assembly_info_organelle_annotation_count.txt

| Contig   | Length | V3  | V4  | V5  | Copy | V7   | V8  | MT  | PT  | Edge   |
| -------- | ------ | --- | --- | --- | ---- | ---- | --- | --- | --- | ------ |
| contig_3 | 70078  | 6   | N   | N   | 1    | both | \*  | 19  | 41  | 1,3,-1 |
| contig_1 | 39736  | 6   | N   | Y   | 1    | left | \*  | 16  | 16  | 1      |
| contig_2 | 19223  | 3   | N   | N   | 1    | both | \*  | 11  | 13  | 2      |

Now that you have seen the table, edit `o/0/mt.contig.name-1`
to add lines of mtDNA contig candidates.

    edge_1
    edge_2
    edge_3

You must have no empty lines in `o/0/mt.contig.name-1`.
Once you have your `mt.contig.name-1` file ready, then:
Your next command is:

    (polap) $ polap assemble2

Your draft assembly graph is `o/1/30-contigger/graph_final.gfa`

    (polap) $ polap flye-polishing

Your long-read polished assembly graph is `o/1/assembly_graph.gfa`
You could extract a draft organelle genome sequence from the assembly graph
by using [Bandage](https://rrwick.github.io/Bandage/).
You could watch [YouTube](https://youtu.be/UF-0UIc2ZDY) (in Korean)
to see how you could extract a DNA sequence through a circular path.

# Options

## Menu options

`assemble1`

: Flye whole-genome assembly

`annotate`

: Organelle gene annotation

`assemble2`

: Flye organelle-genome assembly

`flye-polishing`

: Flye assembly long-read polishing

`prepare-polishing`

: [FMLRC](https://github.com/holtjma/fmlrc) short-read polishing
preparation

`polish`

: [FMLRC](https://github.com/holtjma/fmlrc) short-read polishing

## General options

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
