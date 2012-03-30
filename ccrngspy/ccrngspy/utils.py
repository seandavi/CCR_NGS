import subprocess
import shlex
import logging
import tempfile
import time

def safe_run(cmd, shell=False):
    proc = None

    if (isinstance(cmd, list)):
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    else:
        if shell:
            proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        else:
            proc = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    so, se = proc.communicate()

    return (so, se)

# Could have these in the header too
# #PBS -N %(jobname)s
# #PBS -k oe

_script_header = """
#!/bin/bash
"""

def safe_qsub_run(cmd, jobname, script_header=_script_header, nodes=1, params="", stdout=None, stderr=None, shell=False):
    """Run a command via qsub in blocking mode so that the command waits to exit.

    Requires a header string and a job name.
   
    """
    
    scriptfile = tempfile.NamedTemporaryFile()
    scriptfile.write("%(header)s\n%(command)s\n" % dict(header=script_header, command=cmd))
    scriptfile.file.flush()

    qsub_cmd = "qsub -N %(jobname)s -l nodes=%(nodes)s -W block=true %(params)s" % dict(jobname=jobname,
                                                                                        nodes=nodes,
                                                                                        params=params)

    if stdout:
        qsub_cmd += " -o %s" % stdout

    if stderr:
        qsub_cmd += " -e %s" % stderr
    
    qsub_cmd = "%(cmd)s %(script)s" % dict(cmd=qsub_cmd, script=scriptfile.name)

    proc = subprocess.Popen(qsub_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    so, se = proc.communicate()
    jobid = so.rstrip()
    scriptfile.close()

    return jobid, se

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


