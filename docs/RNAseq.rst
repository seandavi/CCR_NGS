Overview
########
Our goals with RNA-seq data analysis are manifold.  In general, we want to pipeline where possible and document where pipelining does not allow a complete solution.

We are going to need a "public space" for storing genomes and transcripts.  I typically create a "public" folder just for this type of thing.  

Input Files
-----------
The basic set of input files includes the FASTQ file(s) for each
sample, a tab-separated file with the sample names and file names, and
a YAML-formatted configuration file.

Quality Control
---------------
The first step is to analyze the basic output of the sequencing run - the reads.

FastQC
======
The FastQC program reports metrics about any type of sequencing reads.
It is required to have the FastQC program installed and available on your PATH.

The input files are the tab-separated sample file and the YAML
configuration file.

The output directory must be set in the configuration file.

Picard
======
We will use the Picard sub-program CollectRnaSeqMetrics to get further QC information:

http://picard.sourceforge.net/command-line-overview.shtml#CollectRnaSeqMetrics

The REQUIRED input files for this are:
1. SAM/BAM alignment files for each sample
2. a refFlat formatted file of gene annotations (http://genome.ucsc.edu/goldenPath/gbdDescriptionsOld.html#RefFlat). An example of this file can be downloaded using a script found in the example_data/picard/ directory.

All parameters besides the SAM/BAM input alignment files should be specified in the YAML config file.



Alignment
---------

GSNAP
=====

Set up indexes and transcript models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



