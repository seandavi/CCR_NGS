"""Wrapper for RUM (RNAseq Unified Mapper).

Kenny Daily, 2012

"""

import subprocess
import argparse
import os
import Task

from ccrngspy import utils

logger = utils.make_local_logger("RUM logging", level="debug", color="green")

class RUMrunner(Task.Task):
    """Container for RUM tasks.
    
    """
    
    _cmd = "RUM_runner.pl"
    
    _opt_lookup = dict(output_directory="-o %s",
                       threads="-t %s",
                       casava="--casava")
    
    def __init__(self, config_file=None, name=None, input_file_list=None, output_directory=None, chunks=1, ram=6):

        self.config_file = config_file
        self.name = name

        self.input_file_list = input_file_list
        
        self.output_directory = output_directory
        self.chunks = chunks
        self.ram = ram
        
    def argparse(self, parser):
        """Add RUM option group to an OptionParser.
        
        """

        group = parser.add_argument_group("FastQC Options")
        group.add_argument("--rum_config_file", dest="rum_config_file", type=str, default=None,
                           help="RUM configuration file.")
        group.add_argument("--rum_run_name", dest="rum_run_name", type=str, default=None,
                           help="RUM run name.")
        group.add_argument("-o", "--rum_outdir", dest="rum_output_directory", type=str, default=None,
                           help="RUM output directory.")
        group.add_argument("--rum_chunks", dest="rum_chunks", type=int, default=1,
                           help="Specifies the number of chunks for RUM.")
        group.add_argument("--rum_ram", dest="rum_ram", type=int, default=1,
                           help="Specifies the amount of ram per chunk.")
        group.add_argument("--rum_read_files", type=str, nargs="+",
                           help="Input fastq files to RUM.")
        
        return parser
    
    def set_options(self, args):
        """Use args from argparse.parse_args to populate the class.

        """
        
        self.__init__(config_file=args.rum_config_file,
                      name=args.rum_run_name,
                      input_file_list=args.rum_read_files,
                      output_directory=args.rum_output_directory,
                      chunks=args.rum_chunks,
                      ram=args.rum_ram)


    def make_command(self):
        _cmd_string = "%(prog)s %(config)s %(files)s %(outputdir)s %(chunks)s %(name)s -ram %(ram)s"

        if self.input_file_list and self.name and self.config_file:
            cmd = _cmd_string % dict(prog=self._cmd, config=self.config_file,
                                     files=",,,".join(self.input_file_list),
                                     outputdir=self.output_directory,
                                     chunks=self.chunks,
                                     ram=self.ram,
                                     name=self.name)
        else:
            raise ValueError("Did not specify any input files!")

        
        logger.debug("Command to run: %s" % (cmd, ))

        return cmd
    
    def run_rum(self):
        """Run the fastqc program from the command line.
        
        Assumes that it is on your PATH.
        
        """

        cmd = self.make_command()
        utils.safe_run(cmd, shell=False)

def main():    
    _test()

def _test():

    rumrunner = RUMrunner()

    usage = "%(prog)s [options] input_files"
    parser = argparse.ArgumentParser(usage=usage)
    parser = rumrunner.argparse(parser)
    
    args = parser.parse_args()

    rumrunner.set_options(args)
    
    rumrunner.run_rum()

if __name__ == "__main__":
    main()
