#!/usr/bin/env python
"""Run the fastqc pipeline.

This runs FastQC on a number of files.

It requires a YAML configuration file with parameters for FastQC (output directory, etc.)
It also requires a samples file that has at least a column named 'sample' and 'filename'.

"""

import csv
import sys
import os
import argparse
import subprocess

from ruffus import *
import yaml

from ccrngspy.tasks import FastQC
from ccrngspy import utils

from ccrngspy.pipeline import fastqc_helpers

logger = utils.make_local_logger("FastQC logging", level="debug", color="green")

parser = argparse.ArgumentParser(description="Run fastqc on files.")

parser.add_argument("--print_only", dest="print_only", action="store_true", default=False,
                    help="Don't run the pipeline, just print what will be run.")

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


test_task_params = make_fastqc_param_list(samples=samples, config=config)

#----------------------------------------------
# begin tasks here
#----------------------------------------------

@files(test_task_params)
def run_fastqc(input, output, params=None):
    """Set up and run the fastqc program.
    
    """

    fastqc_task = FastQC.FastQC(input_files=[input], output_directory=config['fastqc_params']['output_dir'])
    fastqc_task.run_fastqc()

    # post task, touch output file!
    of = file(output, mode="w")
    of.close()

if opts.print_only:
    pipeline_printout(sys.stdout, [run_fastqc])
else:
    pipeline_run([run_fastqc], multiprocess=5)

