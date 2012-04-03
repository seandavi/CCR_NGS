Overview
########
Our goals with RNA-seq data analysis are manifold.  In general, we want to pipeline where possible and document where pipelining does not allow a complete solution.

We are going to need a "public space" for storing genomes and transcripts.  I typically create a "public" folder just for this type of thing.  

Prerequisites
-------------

There are a number of dependencies that need to be installed to run this pipeline:

Python dependencies
===================

1. Ruffus

Task dependencies
=================
1. Fastqc (>= v0.10.0) must be installed, and the 'fastqc' executable must be on your PATH.
   If running on Biowulf, then this has already been installed.
2. Picard must be installed, and the location of the JAR files known.
   If on Biowulf, this is already installed.
3. RUM (>= version 1.11) must be installed and the EXECUTABLE RUM_runner.pl must be on your PATH.
   If on Biowulf, this will NOT work with the installed version, since they have altered it to only run from a perl interpreter.


Example Run
-----------
First, an example command line run of the pipeline. This must be run on a PBS submit host!

scripts/./rnaseq_qc_qsub.py --config_file=example_data/config/test.yaml --sample_file=example_data/config/test.tsv

Input Files
-----------
The basic set of input files includes the FASTQ file(s) for each sample, a tab-separated file with the sample names and file names, and a YAML-formatted configuration file.
For examples, see the example_data/config/ folder.

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

Alignment
---------

We will then align the sequencing reads to the genome (and transcriptome).
Our current options for alignment include GSNAP (ongoing implementation) and RUM (completed implementation).

RUM
===

The RNAseq Unified Mapper (RUM) has been implemented in the pipeline.

The REQUIRED input files are:

1. A sorted FASTQ file (or two sorted FASTQ files for paired end) of reads to align. In the pipeline, this is determined from the tab-separated sample file. The location of these files should be specified in the YAML config file, and all files need to be in the same directory.
2. A RUM config file (if running on Biowulf, these files are already available for the majority of model organisms), specified in the YAML config file.
3. A sample name, specified in the YAML config file
4. Number of chunks (see RUM documentation for more information on chunking), specified in the YAML config file
5. Amount of RAM to use per chunk (again, see RUM documentation for information), specified in the YAML config file


Picard
======
We will use the Picard sub-program CollectRnaSeqMetrics to get further QC information:

http://picard.sourceforge.net/command-line-overview.shtml#CollectRnaSeqMetrics

The REQUIRED input files for this are:
1. SAM/BAM alignment files for each sample
2. a refFlat formatted file of gene annotations (http://genome.ucsc.edu/goldenPath/gbdDescriptionsOld.html#RefFlat). An example of this file can be downloaded using a script found in the example_data/picard/ directory.
3. a FASTA formatted genome sequence file. If running on Biowulf, these files are already available and do not need to be downloaded.
All parameters besides the SAM/BAM input alignment files should be specified in the YAML config file.



Alignment
---------

GSNAP
=====

Set up indexes and transcript models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



