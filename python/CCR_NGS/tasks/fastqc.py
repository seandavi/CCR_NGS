"""Wrapper for FastQC.

Kenny Daily, 2012

"""
import optparse
import logging
import os

logger = logging.getLogger("FastQC Logging")
# logger.setLevel(logging.DEBUG)

from CCR_NGS import utils

class Task(object):
    def optparse(self, parser):
        """Add option group for the particular task to an OptionParser object.
        
        """
        
        raise NotImplementedError

    def set_options(self, opts, args):
        """Use opts and args from optparse.parse_args to set options for the class.
        
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

    def optparse(self, parser):
        """Add FastQC option group to an OptionParser.
        
        """

        group = optparse.OptionGroup(parser, "FastQC Options")
        group.add_option("-o", "--outdir", dest="output_directory", type="string", default="",
                          help="FastQC output directory.")
        group.add_option("--casava", dest="casava", action="store_true", default=False,
                         help="Files come from raw casava output.")
        group.add_option("-t", "--threads", dest="threads", type="int", default=1,
                         help="Specifies the number of files which can be processed simultaneously.")
        parser.add_option_group(group)
        return parser
    
    def set_options(self, opts, args):
        """Use opts and args from optparse.parse_args to populate the class.

        """

        self.input_files = args
        self.output_directory = opts.output_directory
        self.threads = opts.threads
        self.casava = opts.casava

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

        cmd = " ".join([self._cmd, self.make_option_string(), " ".join(self.input_files)])
        
        logger.debug("Command to run: %s" % (cmd, ))
        utils.safe_run(cmd)

def main():    
    _test()

def _test():

    fastqc = FastQC()

    usage = "%prog [options] input_files"
    parser = optparse.OptionParser(usage=usage)
    parser = fastqc.optparse(parser)
    
    (opts, args) = parser.parse_args()

    fastqc.set_options(opts, args)
    fastqc.run_fastqc()
    
if __name__ == "__main__":
    main()
