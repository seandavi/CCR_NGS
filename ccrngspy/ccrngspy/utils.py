import subprocess
import shlex
import logging

def safe_run(cmd, shell=False):
    res = None

    if (isinstance(cmd, list)):
        res = subprocess.call(cmd)
    else:
        if shell:
            res = subprocess.call(cmd, shell=True)
        else:
            res = subprocess.call(shlex.split(cmd))
    
    return(res)

_LOGGING_LEVEL = {'debug': logging.DEBUG,
                  'info': logging.INFO,
                  'warning': logging.WARNING,
                  'error': logging.ERROR,
                  'critical': logging.CRITICAL}

_LOGGING_COLOR = {'yellow' : '\033[93m',
                  'green' : '\033[92m',
                  'blue' : '\033[94m',
                  'red' : '\033[91m',
                  'grey' : '\033[90m',
                  True : '\033[90m',
                  False : '\033[90m',}

def make_local_logger(logger_name, level="info", color=False):
    """Helper function to make local loggers with color.
    
    """
    
    logger = logging.getLogger(logger_name)

    try:
        logger.setLevel(_LOGGING_LEVEL[level])
    except KeyError:
        logger.setLevel(level)
    
    format = "%(asctime)s - %(name)s - " + _LOGGING_COLOR[color] + "%(levelname)s:%(module)s.%(funcName)s\033[0m - %(message)s"
        
    formatter = logging.Formatter(format)

    chandler = logging.StreamHandler()

    try:
        chandler.setLevel(_LOGGING_LEVEL[level])
    except KeyError:
        chandler.setLevel(level)
        
    chandler.setFormatter(formatter)
    logger.addHandler(chandler)
    
    return logger


