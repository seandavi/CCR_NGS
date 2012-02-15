Overview
########
Our goals with RNA-seq data analysis are manifold.  In general, we want to pipeline where possible and document where pipelining does not allow a complete solution.

We are going to need a "public space" for storing genomes and transcripts.  I typically create a "public" folder just for this type of thing.  

Input Files
-----------
Need to decide what format to take from an end user and parse this to a YAML file that will serve as entry into the pipeline.

Quality Control
---------------

The first step is to analyze the basic output of the sequencing run - the reads.

FastQC
======
The FastQC program reports metrics about any type of sequencing reads.
It is required to have the FastQC program installed and available on your PATH.

Picard
======


Alignment
---------

GSNAP
=====

Set up indexes and transcript models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



