import os
from subprocess import run


def collect_logs(path):
    # collect pepsirf logs into a directory at `path`
    # check for there is not a directory for logs already
    if not os.path.isdir(log_dir):
        os.mkdir(log_dir)
    # collect log files
    run((f"mv *.log ./{log_dir}").split(), shell=True, check=True)
