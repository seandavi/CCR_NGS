"""Helper files for pipeline fastqc steps.

"""

import os

from ccrngspy import utils

logger = utils.make_local_logger("FastQC helper logging", level="debug", color="green")

def make_fastqc_param_list(samples, config, params=None):
    """Helper function to turn the sample file into a list of files.

    Needs to be a list of [input, output, params]; for the fastqc script.
    The output is a log file (a sentinel that can be used to check completion),
    while the params are taken from the global opts variable
    (and possibly from the YAML config file).
    
    """

    final_list = []

    fastq_dir = config['general_params']['fastq_input_dir']
    output_dir = config['fastqc_params']['output_dir']
    
    for sample in samples:

        params = dict(sample=sample['samplename'])

        tmp1 = [os.path.join(fastq_dir, sample['filename1']),
                os.path.join(output_dir, "%s_fastqc.zip" % sample['sample1']),
                params]
        tmp2 = [os.path.join(fastq_dir, sample['filename2']),
                os.path.join(output_dir, "%s_fastqc.zip" % sample['sample2']),
                params]
        
        final_list.extend([tmp1, tmp2])

    return final_list
