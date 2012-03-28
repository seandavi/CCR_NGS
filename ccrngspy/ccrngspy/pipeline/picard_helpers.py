"""Helper files for pipeline fastqc steps.

"""

import os

from ccrngspy import utils

logger = utils.make_local_logger("Picard helper logging", level="debug", color="green")

def make_rum_param_list(samples, config, params=None):
    """Helper function to turn the sample file into a list of files.

    Needs to be a list of [[input1, input2], output, params]; for the fastqc script.
    The output is a log file (a sentinel that can be used to check completion),
    while the params are taken from the global opts variable
    (and possibly from the YAML config file).
    
    """

    final_list = []

    rum_dir = config['rum_params']['output_dir']
    log_dir = config['general_params']['log_file_dir']
    
    for sample in samples:
        tmp = [os.path.join(fastq_dir, sample['samplename', RUM.sam]),
               os.path.join(log_dir, "%s.picard.LOG" % sample['samplename']),
               params]
        
        final_list.append(tmp)

    return final_list
