Overview
--------

This tool assembles organelle genomes using seed contig sequences and
long-read read FASTQ data.

Inputs
------

-  **Whole-genome assembly graph**: a Flye’s ``gfa`` format data file.
-  **MT seed edge IDs**: ``mt.contig.name-1``
-  **Minimum length of long-read sequences**: the long-read sequence
   length threshold
-  **Minimum mapping length for read selection**: the sequence length
   threshold
-  **Use Polap’s read selection**: default is ptGAUL’s read selection
-  **Coverage**: maximum coverage of reads for the organelle-genome
   assembly
-  **No coverage check**: no data reduction in an organelle-genome
   assembly

Outputs
-------

-  **organelle-genome assembly**: a Flye’s ``gfa`` format data file.

Usage
-----

To use this tool, upload a Flye’s ``gfa`` assembly file and seed
contigs.

**Example:**

polap assemble2 -w 3000 –polap-reads
