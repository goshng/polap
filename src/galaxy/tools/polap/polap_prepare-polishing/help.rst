Overview
--------

This tool prepares polishing data before actual polishing.

Inputs
------

-  **short-read fastq 1**: a short-read ``fastq`` data file
-  **short-read fastq 2**: another short-read ``fastq`` data file

Outputs
-------

-  **polshing data**: ``msbwt/comp_msbwt.npy``

Usage
-----

To use this tool, upload two short-read ``fastq`` files containing the
data you want to process.

**Example:**

polap prepare-polishing -a s1.fq -b s2.fq
