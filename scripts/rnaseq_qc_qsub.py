#!/usr/bin/env python
"""Run the fastqc pipeline.

This runs FastQC on a number of files, using qsub to submit to the cluster.

It requires a YAML configuration file with parameters for FastQC (output directory, etc.)
It also requires a samples file that has at least a column named 'sample' and 'filename'.

"""

import argparse
import csv
import sys
import os
import argparse
import subprocess

from ruffus import *
import yaml

from ccrngspy.tasks import FastQC
from ccrngspy.tasks import Picard
from ccrngspy.pipeline import fastqc_helpers

parser = argparse.ArgumentParser(description="Run fastqc on files.")

parser.add_argument("--print_only", dest="print_only", action="store_true", default=False,
                    help="Don't run the pipeline, just print what will be run.")

parser.add_argument("--no_log_dir", dest="no_create_log_dir", action="store_true", default=False,
                    help="Don't recreate the output log dir.")

parser.add_argument('--config_file', dest="config_file", type=str,
                    help="A YAML configuration file for pipeline.")

parser.add_argument('--sample_file', dest="sample_file", type=str,
                    help="A YAML configuration file for pipeline.")

# add options for the fastqc task
parser = FastQC.FastQC().argparse(parser)

# Parse the options
opts = parser.parse_args()

# Load the bootstrap config file
with open(opts.config_file, 'r') as configfile:
    config = yaml.load(configfile)

# Load the samples tab-separated file
with open(opts.sample_file, 'r') as samplefile:
    reader = csv.DictReader(samplefile, delimiter="\t")
    samples = list(reader)


test_task_params = fastqc_helpers.make_fastqc_param_list(samples=samples, config=config)

#----------------------------------------------
# begin tasks here
#----------------------------------------------

def run_setup_log_dir(input=None, output=None, params=None):
    if not opts.no_create_log_dir:
        os.mkdir(config['general_params']['log_file_dir'])

@follows(run_setup_log_dir)
@files(test_task_params)
def run_fastqc(input, output, params=None):
    """Set up and run the fastqc program.
    
    """

    fastqc_task = FastQC.FastQC(input_files=[input], output_directory=config['fastqc_params']['output_dir'])
    fastqc_command = fastqc_task.make_command()

    jobid, err = utils.safe_qsub_run(fastqc_command)
    
    # post task, touch output file!
    of = file(output, mode="w")
    of.close()

def run_gsnap(input, output, params=None):
    """Run gsnap.
    
    """

    pass

def run_collect_rnaseq_metrics(input, output, params=None):
    """Set up and run the Picard CollectRnaSeqMetrics program.
    
    """
    
    # Let a parser argument handle setting up arguments and options
    parser = argparse.ArgumentParser()
    
    # Add Picard arguments
    picard = Picard.PicardBase()
    parser = picard.argparse(parser)

    # Update input and output from global config object
    picard_params = config['picard_params']
    picard_params['input'] = input
    picard_params['output'] = output
    
    # Set up using the default arguments, specifying the input and output files since they are required!
    args = parser.parse_args("CollectRnaSeqMetrics --input=%(input) --output=%(output)s --ribosomal_intervals=%(ribosomal_intervals)s --minimum_length=%(minimum_length)s --chart_output=%(chart_output)%s --rrna_fragment_percentage=%(rrna_fragment_percentage)s --metric_accumulation_level=%(metric_accumulation_level)s --stop_after=%(stop_after)s --ref_flat=%(ref_flat) --ref_file=%(ref_file)" % picard_params)
    
    # Run the function for collecting RNASeq metrics
    args.func(args)
    
if opts.print_only:
    pipeline_printout(sys.stdout, [run_setup_log_dir, run_fastqc])
else:
    pipeline_run([run_setup_log_dir, run_fastqc], multiprocess=5)

