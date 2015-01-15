#! /bin/usr/env python
# 25/11/2014
# G. Gorman's script for running pG-lt stages 3-4 in parallel with MPI

# Warning start
print 'WARNING: This script is for running pG-lt in parallel. It must only be \
run from a .sh script and for folders for which pG-lt stages 1-2 have already \
been executed. It requires MPI and mpi4py.'

import os
import sys
import pickle
from mpi4py import MPI
from pglt.tools.system_tools import Stager


# GLOBALS
comm = MPI.COMM_WORLD
rank = comm.Get_rank()


# FUNCTIONS
def generate_trees(folder):
    '''Takes a setup folder containing a list of names and sequence data,
        generates a distribution of trees.'''
    print "running job ", folder
    # Run stages 3 (alignment) and 4 (phylogeny)
    Stager.run_all(wd=folder, stages=['3', '4'])

    # Check with Dom that the data has also been written to disk at this point.
    return


def dummy_job(folder):
    '''This only exists for testing the task farm.'''
    print "running job ", folder


def master(worklist):
    '''The task master takes a work list and farms it out to the workers.'''

    # Check we are not using more ranks than what we need.
    nranks = comm.Get_size()
    if len(worklist) < (nranks-1):
        print "WARNING: There are fewer tasks than ranks. You are wasting \
resources!"
        for rank in range(len(worklist)+1, nranks):
            # Catch the send
            buff = comm.recv(source=rank)

            # Send completion signal
            comm.send("", dest=rank)

        # Reset number of ranks.
        nranks = len(worklist)+1

    # Farm out tasks
    for task in worklist:
        # Probe for any idle worker.
        status = MPI.Status()
        buff = comm.recv(source=MPI.ANY_SOURCE, status=status)
        idle_worker = status.Get_source()

        # Give worker task
        comm.send(task, dest=idle_worker)

    # Send null task to signal completion.
    for rank in range(1, nranks):
        comm.send("", dest=rank)

    return


def worker(task_operation):
    '''The worker takes as an input the operation it should apply to the tasks
       that it receives from the task master.'''

    while True:
        # Alert master that we are waiting.
        comm.send("")

        # Wait for a task
        task = comm.recv()

        if task == "":
            break
        else:
            task_operation(task)
    return


def getWorklist(folder):
    '''Return list of stage folders'''
    worklist = os.listdir(folder)
    res = []
    for each in worklist:
        directory = os.path.join(folder, each, 'tempfiles')
        try:
            with open(os.path.join(directory, 'progress.p'), "rb") as file:
                progress = pickle.load(file)
            if progress['2'] != 'not run':
                res.append(each)
        except:
            pass
    return res


if __name__ == '__main__':
    # read parent folder here (use abs path)!
    parent_folder = sys.argv[1]
    print 'Using: ', parent_folder
    # obtain list of stage folders within parent_folder
    worklist = getWorklist(parent_folder)
    worklist = [os.path.join(parent_folder, e) for e in worklist]
    if rank == 0:
        master(worklist)
    else:
        worker(generate_trees)
