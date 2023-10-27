import os
import glob
import shutil


def collect_logs(path):
    # collect pepsirf logs into a directory at `path`
    # check for there is not a directory for logs already
    if not os.path.isdir(path):
        os.mkdir(path)
    # collect log files
    files = glob.glob("*.log", recursive=True)
    for file in files:
        shutil.move(file, path)
