Overview
--------

This tool estimates the genome size by applying a paired-end short-read
dataset to JellyFish.

Inputs
------

-  **short-read fastq 1**: a short-read ``fastq`` data file
-  **short-read fastq 2**: another short-read ``fastq`` data file

Outputs
-------

-  **Output File**: a text file named ``short_expected_genome_size.txt``
   with the genome size in base pairs.

Usage
-----

To use this tool, upload two short-read ``fastq`` files containing the
data you want to process.

**Example:**

polap find-genome-size -a s1.fq -b s2.fq
