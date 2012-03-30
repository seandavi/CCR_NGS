#!/usr/bin/env python
"""Run the RNAseq pipeline.

This runs the RNAseq pipeline on a number of files, using qsub to submit to the cluster.

It requires a YAML configuration file with parameters (output directory, etc.)

It also requires a samples file that has at least columns:
    'samplename', 'sample1', 'sample2', 'filename1', and 'filename2'

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
from ccrngspy.tasks import RUM
from ccrngspy.pipeline import fastqc_helpers
from ccrngspy.pipeline import rum_helpers
from ccrngspy.pipeline import picard_helpers
from ccrngspy import utils

logger = utils.make_local_logger("Ruffus RNASeq QC Logger", level="debug", color=True)

parser = argparse.ArgumentParser(description="Run fastqc on files.")

parser.add_argument("--print_only", dest="print_only", action="store_true", default=False,
                    help="Don't run the pipeline, just print what will be run.")

parser.add_argument("--no_log_dir", dest="no_create_log_dir", action="store_true", default=False,
                    help="Don't recreate the output log dir.")

parser.add_argument("--no_output_dir", dest="no_create_output_dir", action="store_true", default=False,
                    help="Don't recreate the output dirs.")

parser.add_argument('--config_file', dest="config_file", type=str,
                    help="A YAML configuration file for pipeline.")

parser.add_argument('--sample_file', dest="sample_file", type=str,
                    help="A tab separated file about the samples to run.")

# add options for the fastqc task
# parser = FastQC.FastQC().argparse(parser)

# Parse the options
opts = parser.parse_args()

# Load the bootstrap config file
with open(opts.config_file, 'r') as configfile:
    config = yaml.load(configfile)

# Load the samples tab-separated file
with open(opts.sample_file, 'r') as samplefile:
    reader = csv.DictReader(samplefile, delimiter="\t")
    samples = list(reader)

# setup fastqc specific params
fastqc_test_task_params = fastqc_helpers.make_fastqc_param_list(samples=samples, config=config)

# setup rum specific params
rum_test_task_params = rum_helpers.make_rum_param_list(samples=samples, config=config, params=None)

#----------------------------------------------
# begin tasks here
#----------------------------------------------

@follows(mkdir(config['general_params']['log_file_dir']),
         mkdir(config['fastqc_params']['output_dir']),
         mkdir(config['rum_params']['output_dir']),
         mkdir(config['picard_params']['output_dir']))
def run_setup_dir(input=None, output=None, params=None):
    pass

@follows(run_setup_dir)
def run_mk_output_dir(input=None, output=None, params=None):
    if not opts.no_create_output_dir:
        # # Make output directories for each task
        # os.mkdir(config['fastqc_params']['output_dir'])
        # os.mkdir(config['rum_params']['output_dir'])
        # os.mkdir(config['picard_params']['output_dir'])
        
        # Make RUM output directory for each sample
        for sample in samples:
            sample_output_dir = os.path.join(config['rum_params']['output_dir'], sample['samplename'])
            os.mkdir(sample_output_dir)

        # # Make picard output directory for each sample
        # for sample in samples:
        #     sample_output_dir = os.path.join(config['picard_params']['output_dir'], sample['samplename'])
        #     os.mkdir(sample_output_dir)

@follows(run_mk_output_dir)
@files(fastqc_test_task_params)
def run_fastqc(input, output, params=None):
    """Set up and run fastqc.
    
    """

    # Let a parser argument handle setting up arguments and options
    parser = argparse.ArgumentParser()
    
    # Add Picard arguments
    fastqc = FastQC.FastQC()
    parser = fastqc.argparse(parser)
    
    # Update input and output from global config object
    fastqc_params = config['fastqc_params']
    fastqc_params['input'] = input

    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['log_file_dir']
    stderr = config['general_params']['log_file_dir']
    
    cmdline = "--outdir=%(output_dir)s --threads=%(threads)s %(input)s" % fastqc_params

    args = parser.parse_args(cmdline.split())
    fastqc.set_options(args)

    # Final command to run
    fastqc_command = fastqc.make_command()
    
    # if fastqc_params['run_type'] == 'remote':
    #     stdout, stderr = utils.safe_qsub_run(fastqc_command, jobname="run_fastqc")
    # elif fastqc_params['run_type'] == 'local':
    job_stdout, job_stderr = utils.safe_qsub_run(fastqc_command, jobname="fastqc",
                                                 nodes="1:m1:c2",
                                                 stdout=stdout, stderr=stderr)
    
    logger.debug("stdout = %s, stderr = %s" % (job_stdout, job_stderr))

    # post task, touch output file!
    of = file(output, mode="w")
    of.close()

@follows(run_mk_output_dir)
@files(rum_test_task_params)
def run_rum(input, output, params=None):
    """Run RUM on paired reads.
    
    """
    
    # Let a parser argument handle setting up arguments and options
    parser = argparse.ArgumentParser()
    
    # Add Picard arguments
    rum = RUM.RUMrunner()
    parser = rum.argparse(parser)
    
    # Update input and output from global config object
    rum_params = config['rum_params']
    rum_params['input'] = input

    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['log_file_dir']
    stderr = config['general_params']['log_file_dir']

    ## fastq read files
    rum_params['file1'] = input[0]
    rum_params['file2'] = input[1]
    rum_params['sample'] = params['sample']
    
    cmdline = "--rum_config_file=%(config_file)s --rum_run_name=%(sample)s --rum_outdir=%(output_dir)s/%(sample)s --rum_read_files %(file1)s %(file2)s --rum_chunks=%(chunks)s --rum_ram=4" % rum_params
    args = parser.parse_args(cmdline.split())

    rum.set_options(args)
    
    rum_command = rum.make_command()

    # stdout, stderr = utils.safe_run(rum_command, shell=False)
    # logger.debug("stdout = %s, err = %s" % (stdout, stderr))

    job_stdout, job_stderr = utils.safe_qsub_run(rum_command, jobname="rum_%s" % params['sample'],
                                                 nodes="1:m4:c2",
                                                 stdout=stdout, stderr=stderr)
    
    logger.debug("stdout = %s, stderr = %s" % (job_stdout, job_stderr))


# def run_gsnap(input, output, params=None):
#     """Run gsnap.
    
#     """

#     pass

@transform(run_rum, regex(r"(.*).sam"), r"\1.sorted.sam")
def run_sort_sam(input, output, params=None):
    """Set up and run the Picard SortSam program.
    
    """
    
    # # Let a parser argument handle setting up arguments and options
    # parser = argparse.ArgumentParser()
    
    # # Add Picard arguments
    # picard = Picard.PicardBase()
    # parser = picard.argparse(parser)

    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['log_file_dir']
    stderr = config['general_params']['log_file_dir']

    # Update input and output from global config object
    picard_params = config['picard_sortsam_params']

    picard_params['input'] = input
    picard_params['output'] = output

    logger.debug("picard_params = %s" % (picard_params,))
    # Set up using the default arguments, specifying the input and output files since they are required!
    cmdline = "--jar=%(jar_file)s --input=%(input)s --output=%(output)s --sort_order=%(sort_order)s SortSam" % picard_params

    picard_cmd = "python -m ccrngspy.tasks.Picard %s" % cmdline

    # stdout, stderr = utils.safe_run(picard_cmd, shell=False)
    # logger.debug("stdout = %s, err = %s" % (stdout, stderr))
    
    logger.debug("params = %s" % (params, ))
    job_stdout, job_stderr = utils.safe_qsub_run(picard_cmd, jobname="sort_sam",
                                                 nodes="1",
                                                 stdout=stdout, stderr=stderr)
    
    logger.debug("stdout = %s, stderr = %s" % (job_stdout, job_stderr))

@transform(run_sort_sam, regex(r".*/(.*)/RUM.sorted.sam"), r"%s/\1.tsv" % config['picard_params']['output_dir'], r"\1")
def run_collect_rnaseq_metrics(input, output, sample):
    """Set up and run the Picard CollectRnaSeqMetrics program.
    
    """
    
    # # Let a parser argument handle setting up arguments and options
    # parser = argparse.ArgumentParser()
    
    # # Add Picard arguments
    # picard = Picard.PicardBase()
    # parser = picard.argparse(parser)

    # Output dir for qsub stdout and stderr
    stdout = config['general_params']['log_file_dir']
    stderr = config['general_params']['log_file_dir']

    # Update input and output from global config object
    picard_params = config['picard_params']
    picard_params['input'] = input
    picard_params['output'] = output
    
    # Set up using the default arguments, specifying the input and output files since they are required!
    cmdline = "--jar=%(jar_file)s --input=%(input)s --output=%(output)s --ref_flat=%(ref_flat)s --ref_file=%(ref_file)s CollectRnaSeqMetrics --minimum_length=%(minimum_length)s --chart_output=%(chart_output)s --metric_accumulation_level=%(metric_accumulation_level)s --stop_after=%(stop_after)s" % picard_params

    # args = parser.parse_args(cmdline.split())
    
    # # Run the function for collecting RNASeq metrics
    # args.func(args)
    
    picard_cmd = "python -m ccrngspy.tasks.Picard %s" % cmdline

    # stdout, stderr = utils.safe_run(picard_cmd, shell=False)
    # logger.debug("stdout = %s, err = %s" % (stdout, stderr))

    job_stdout, job_stderr = utils.safe_qsub_run(picard_cmd, jobname="rum_%s" % sample,
                                                 nodes="1",
                                                 stdout=stdout, stderr=stderr)
    
    logger.debug("stdout = %s, stderr = %s" % (job_stdout, job_stderr))

@merge(run_collect_rnaseq_metrics, os.path.join(config["picard_params"]["output_dir"], "CollectRNASeqMetrics.tsv"))
def run_merge_rnaseq_metrics(input_files, summary_file):
    """Merge the outputs of collectrnaseqmetrics into one file.

    """

    metrics = []
  
    for fn in input_files:
        metrics.extend(picard_helpers.parse_picard_rnaseq_metrics(fn))

    fieldnames = metrics[0].keys()
    
    with open(summary_file, 'w') as fou:
        dw = csv.DictWriter(fou, delimiter='\t', fieldnames=fieldnames)
        dw.writeheader()
        dw.writerows(metrics)

job_list = [run_setup_dir, run_mk_output_dir, run_fastqc, run_rum, run_sort_sam, run_collect_rnaseq_metrics, run_merge_rnaseq_metrics]

if opts.print_only:
    pipeline_printout(sys.stdout, job_list, verbose=3)
else:
    pipeline_run(job_list, multiprocess=10, logger=logger)
