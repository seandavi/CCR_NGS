import subprocess
import shlex
import logging
import tempfile
import time

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

_script_header = """
#!/bin/bash
#PBS -N %(jobname)s
#PBS -k oe
"""

def check_job(jobid):
    """Check if a PBS job with 'jobid' is finished or not.

    Does not currently implement any error checking, only if the
    job is finished (it will not show up in qstat).

    This only works with single jobs as well!
    
    """

    # get the job from qstat -f, get the job state line, assume its the first one, cut out the status
    cmd = 'qstat -f %(jobid)s | grep job_state | head -n 1 | sed "s/.*job_state = \(.*\)/\1/"' % dict(jobid=jobid)

    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    status, stderr = proc.communicate()

    status = status.rstrip()

    # This is important - this will appear when the job no longer exists
    # It could be that it crashed or otherwise failed!
    test = "qstat: Unknown Job Id %s\n" % jobid
    done = False

    if stderr == test:
        done = True

    return status, done

def safe_qsub_run(cmd, shell=False):
    """Run a command via qsub.

    Requires a header string, a job name.

    Uses a hacky testing command using qstat to determine if the script is finished.
    
    """
    
    foo = tempfile.NamedTemporaryFile()
    foo.write(_script_header % {'jobname': 'foobarbaz'})

    foo.write("%s\n" % cmd)
    foo.file.flush()

    qsub_cmd = "qsub -l nodes=1 %s" % foo.name

    print foo.name
    proc = subprocess.Popen(qsub_cmd, shell=True, stdout=subprocess.PIPE)
    jobid = proc.communicate()[0].rstrip()

    done = False
    while not done:
        job_status, done = check_job(jobid)
        time.sleep(5)

    foo.close()

    return proc

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


