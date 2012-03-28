"""Wrapper for FastQC.

Kenny Daily, 2012

"""

import subprocess
import argparse
import os
import Task

from ccrngspy import utils

logger = utils.make_local_logger("FastQC logging", level="debug", color="green")
    
    
class FastQC(Task.Task):
    """Container for FastQC tasks.
    
    """
    
    _cmd = "fastqc"
    
    _opt_lookup = dict(output_directory="-o %s",
                       threads="-t %s",
                       casava="--casava")
    
    def __init__(self, input_files=None, output_directory=None, threads=1, casava=False):
        self.input_files = input_files
        self.output_directory = output_directory
        self.threads = threads
        self.casava = casava

    def argparse(self, parser):
        """Add FastQC option group to an OptionParser.
        
        """

        group = parser.add_argument_group("FastQC Options")
        group.add_argument("-o", "--outdir", dest="output_directory", type=str, default=None,
                           help="FastQC output directory.")
        group.add_argument("--casava", dest="casava", action="store_true", default=False,
                           help="Files come from raw casava output.")
        group.add_argument("-t", "--threads", dest="threads", type=int, default=1,
                           help="Specifies the number of files which can be processed simultaneously.")
        group.add_argument("input_files", type=str, nargs="*",
                           help="Input fastq files.")

        return parser
    
    def set_options(self, args):
        """Use args from argparse.parse_args to populate the class.

        """

        logger.debug("input_files = %s" % (args.input_files, ))
        self.__init__(input_files=args.input_files, output_directory=args.output_directory,
                      threads=args.threads, casava=args.casava)

    def make_option_string(self):
        """Use currently set options to create the option string.
        
        """
        
        option_string = []
        
        # add output directory
        output_dir = "-o %s" % self.output_directory

        #threads
        threads = "-t %s" % self.threads

        # using casava files?
        if self.casava:
            casava = "--casava"
        else:
            casava = ""

        return " ".join([output_dir, threads, casava])

    def make_command(self):
        if self.input_files:
            cmd = " ".join([self._cmd, self.make_option_string(), " ".join(self.input_files)])
        else:
            raise ValueError("Did not specify any input files!")
        
        logger.debug("Command to run: %s" % (cmd, ))

        return cmd
    
    def run_fastqc(self):
        """Run the fastqc program from the command line.
        
        Assumes that it is on your PATH.
        
        """

        cmd = self.make_command()
        utils.safe_run(cmd, shell=False)

def main():    
    _test()

def _test():

    fastqc = FastQC()

    usage = "%(prog)s [options] input_files"
    parser = argparse.ArgumentParser(usage=usage)
    parser = fastqc.argparse(parser)
    
    args = parser.parse_args()

    fastqc.set_options(args)
    
    fastqc.run_fastqc()

if __name__ == "__main__":
    main()
