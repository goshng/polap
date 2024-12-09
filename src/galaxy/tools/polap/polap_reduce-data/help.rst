Overview
--------

This tool reduces the long-read data for a Flye genome assembly.

Inputs
------

-  **Long-read FASTQ**: a long-read ``fastq`` format data file.
-  **Long-read data size**: ``long_total_length.txt`` (output of polap
   total-length-long)
-  **Genome size**: ``short_expected_genome_size.txt`` (output of polap
   find-genome-size)
-  **Minimum length of long-read sequences**: the long-read sequence
   length threshold
-  **Long-read sampling minimum coverage**: the target coverage for
   sampling the long-read data
-  **No long-read sampling**: **Long-read sampling minimum coverage** is
   ignored.
-  **Random seed for long-read sampling**: seqkit default seed
-  **Use test data**: used only with the test dataset

Outputs
-------

-  **Long-read FASTQ file for organelle-genome assembly**: a compressed
   FASTQ file named ``lk.fq.gz`` with read equences less than **Minimum
   length of long-read sequences**
-  **Long-read FASTQ file for whole-genome assembly**: a compressed
   FASTQ file named ``nk.fq.gz`` with read equences less than **Minimum
   length of long-read sequences** and sampled upto the **Long-read
   sampling minimum coverage** considering **Genome size** and
   **Long-read data size**

Usage
-----

To use this tool, set the **Long-read data size** and **Genome size**
and upload a long-read ``fastq`` file containing the data you want to
process.

**Example:**

polap reduce-data -l l.fq
