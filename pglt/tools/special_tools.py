#! /bin/usr/env python
# D.J. Bennett
# 06/11/2014
"""
pglt special tools
"""

# PACKAGES
import os
import time
import sys
import subprocess
import shutil
import pickle


# FUNCTIONS
def timeit(func, **kwargs):
    """Time a function (platform independent)"""
    # http://stackoverflow.com/questions/10884751/python-timeit-and-program-output
    if sys.platform == "win32":
        # On Windows, the best timer is time.clock()
        default_timer = time.clock
    else:
        # On most other platforms the best timer is time.time()
        default_timer = time.time
    t0 = default_timer()
    output = func(**kwargs)
    t1 = default_timer()
    return output, t1-t0


def getThreads(wd=None):
    """Find number of cores on machine (platform independent). If local, \
search cwd for pickled .threads.p."""
    if wd:
        if os.path.isfile(os.path.join('.threads.p')):
            with open('.threads.p', "rb") as file:
                nthreads = pickle.load(file)
        else:
            nthreads = 1
        return nthreads
    if sys.platform == 'darwin':
        # find threads for a mac
        cmd = "system_profiler SPHardwareDataType | grep Cores | awk \
'{print $5}'"
    elif sys.platform == 'linux2':
        # find threads for a unix
        cmd = "cat /proc/cpuinfo | grep processor | tail -1 | awk '{print $3}'"
    else:
        # TODO: work out n cores for windows
        # TODO: add a threads parameter to parameters.csv
        return 1
    s = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         shell=True)
    nthreads = s.stdout.readline()
    try:
        nthreads = int(nthreads)
    except:
        nthreads = 1
    return nthreads


def clean():
    """Remove all pglt files and folders"""
    # remove log.txt in parent folder
    try:
        os.remove('log.txt')
    except OSError:
        pass

    # go through all subdirs and remove files in list
    files_to_remove = ['info.txt', 'log.txt', '.paradict.p', '.allrankids.p',
                       '.namesdict.p', '.genedict.p', '.terms.p',
                       '.constraint.tre', '.partitions.txt',
                       '.phylogeny_in.phylip', '.alignment_out.fasta',
                       '.threads.p']
    folders = os.listdir(os.getcwd())
    while folders:
        temp_folder = folders.pop()
        temp_files_to_remove = files_to_remove[:]
        while temp_files_to_remove:
            temp_file_to_remove = temp_files_to_remove.pop()
            try:
                os.remove(os.path.join(temp_folder, temp_file_to_remove))
            except OSError:
                pass

    # go through all subdirs and remove folders in list
    folders_to_remove = ['1_names', '2_download', '3_alignment', '4_phylogeny',
                         'tempfiles']
    folders = os.listdir(os.getcwd())
    while folders:
        temp_folder = folders.pop()
        temp_folders_to_remove = folders_to_remove[:]
        while temp_folders_to_remove:
            temp_folder_to_remove = temp_folders_to_remove.pop()
            try:
                shutil.rmtree(os.path.join(temp_folder, temp_folder_to_remove))
            except OSError:
                pass
