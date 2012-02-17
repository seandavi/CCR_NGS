"""Generic Task class.

"""

class Task(object):
    def argparse(self, parser):
        """Add option group for the particular task to an OptionParser object.
        
        """
        
        raise NotImplementedError
    
    def set_options(self, opts, args):
        """Use opts and args from argparse.parse_args to set options for the class.
        
        """
        
        raise NotImplementedError
