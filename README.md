## TagSeq pipeline
This directory includes snakemake pipeline for processing TagSeq data:
The pipeline includes the following steps.
1. QC using FastQC and multiQC for summarizing and visualizing
2. Adapter trimming by cutadapt and multiQC for summarizing and visualizing
3. Removing Swith oligos sequences in the reads by using UMI-Tools. Here reads are switch oligo removed and separated into 3 parts since there are 3 types of oligos.
4. Mapping with STAR for each part.
5. Deduplicate with UMI-Tools for each part
6. Quantify by featureCounts for each part
7. Add up the counts of the 3 parts by an in-house Python script.

## epiRIL
The scripts include differential expression on >100 epiRILs.
