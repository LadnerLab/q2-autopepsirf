import os
import subprocess


def collect_logs(path):
    # collect pepsirf logs into a directory at `path`
    # check for there is not a directory for logs already
    if not os.path.isdir(path):
        os.mkdir(path)
    # collect log files
    subprocess.run((f"mv *.log {path}").split(), shell=True, check=True)
