#!/usr/bin/env python

import sys
import os
import argparse
import subprocess

from ruffus import *
import yaml

from ccrngspy.tasks import FastQC

fastqc_task = FastQC.FastQC()

parser = argparse.ArgumentParser(description="Run fastqc on files.")

parser.add_argument("--print_only", dest="print_only", action="store_true", default=False,
                    help="Don't run the pipeline, just print what will be run.")

# parser.add_argument('config',
#                     help="A YAML configuration file for pipeline.")

## add options for the fastqc task
parser = fastqc_task.argparse(parser)

tmp_args = "-o /tmp/"
opts = parser.parse_args()

print opts

# with open(opts.config,'r') as configfile:
#     config = yaml.load(configfile)

#----------------------------------------------
# begin tasks here
#----------------------------------------------

def make_fastqc_param_list(file_list):
    """Helper function to turn a file list into a list for ruffus.

    Needs to be a list of [input, output, params]; for the fastqc,
    the output is none, while the params are taken from the global opts variable
    (and possibly from the YAML config file).
    
    """
    
    return map(lambda x: [x, None, None], file_list)

file_list = ["/var/preserve/git/CCR_NGS/example_data/fastq/s_1_1.fastq.gz",
             "/var/preserve/git/CCR_NGS/example_data/fastq/s_1_2.fastq.gz"]

test_task_params = make_fastqc_param_list(file_list)

@files(test_task_params)
def run_fastqc(input, output, params=None):
    """Set up and run the fastqc program.
    
    """

    fastqc_task = FastQC.FastQC(input_files=[input], output_directory=opts.output_directory)

    cmd = fastqc_task.make_command()

    subprocess.call(cmd, shell=True)
    # fastqc_task.run_fastqc()


if opts.print_only:
    pipeline_printout(sys.stdout, [run_fastqc])
else:
    pipeline_run([run_fastqc])

# def sam2pindel(input,output):
#     print input,output
#     params = {'samtools':config['samtools'],
#               'bamfile':opts.filename,
#               'sam2pindel':config['sam2pindel'],
#               'outfile':str(output),
#               'insertsize':opts.insertsize,
#               'tag':opts.tag}
#     cmd = "%(samtools)s view %(bamfile)s | %(sam2pindel)s - %(outfile)s %(insertsize)d %(tag)s 0" % params
#     subprocess.call(cmd,shell=True)
