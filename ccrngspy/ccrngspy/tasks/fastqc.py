"""Wrapper for FastQC.

Kenny Daily, 2012

"""
import argparse
import os

from ccrngspy import utils

logger = utils.make_local_logger("FastQC logging", level="debug", color="green")

class Task(object):
    def argparse(self, parser):
        """Add option group for the particular task to an OptionParser object.
        
        """
        
        raise NotImplementedError
    
    def set_options(self, opts, args):
        """Use opts and args from argparse.parse_args to set options for the class.
        
        """
        
        raise NotImplementedError
    
    
class FastQC(Task):
    """Container for FastQC tasks.
    
    """
    
    _cmd = "fastqc"
    
    _opt_lookup = dict(output_directory="-o %s",
                       threads="-t %s",
                       casava="--casava")
    
    def __init__(self, input_files=None, output_directory=None):
        self.input_files = input_files
        self.output_directory = output_directory

    def argparse(self, parser):
        """Add FastQC option group to an OptionParser.
        
        """

        group = parser.add_argument_group("FastQC Options")
        group.add_argument("-o", "--outdir", dest="output_directory", type=str,
                           help="FastQC output directory.")
        group.add_argument("--casava", dest="casava", action="store_true", default=False,
                           help="Files come from raw casava output.")
        group.add_argument("-t", "--threads", dest="threads", type=int, default=1,
                           help="Specifies the number of files which can be processed simultaneously.")
        group.add_argument("input_files", type=str, nargs="+",
                           help="Input fastq files.")

        return parser
    
    def set_options(self, args):
        """Use opts and args from argparse.parse_args to populate the class.

        """

        self.__init__(input_files=args.input_files, output_directory=args.output_directory)
        self.threads = args.threads
        self.casava = args.casava

    def make_option_string(self):
        """Use currently set options to create the option string.
        
        """
        
        option_string = []
        for (optname, optstring) in self._opt_lookup.items():
            try:
                option_string.append(optstring % self.__dict__[optname])
            except TypeError: # if the option doesn't have a parameter
                option_string.append(optstring)

        return " ".join(option_string)

    def run_fastqc(self):
        """Run the fastqc program from the command line.

        Assumes that it is on your PATH.
        
        """

        if self.input_files:
            cmd = " ".join([self._cmd, self.make_option_string(), " ".join(self.input_files)])
        else:
            raise ValueError("Did not specify any input files!")
        
        logger.debug("Command to run: %s" % (cmd, ))
        utils.safe_run(cmd)

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
