#! /bin/usr/env python
# G. Gorman

import os
from pglt.tools.system_tools import Stager


def generate_trees(folder):
    '''Takes a setup folder containing a list of names and sequence data,
  generates a distribution of trees.

    '''
    print "running job ", folder
    # Run stages 3 (alignment) and 4 (phylogeny)
    Stager.run_all(wd=folder, stages=['3', '4'])

    # Check with Dom that the data has also been written to disk at this point.
    return


def dummy_job(folder):
    '''This only exists for testing the task farm.'''
    print "running job ", folder

from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()


def master(worklist):
    '''The task master takes a work list and farms it out to the workers.'''

    # Check we are not using more ranks than what we need.
    nranks = comm.Get_size()
    if len(worklist) < (nranks-1):
        print "WARNING: There are fewer tasks than ranks. You are wasting resources!"
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
    '''The worker takes as an input the operation it should apply to the tasks that it receives from the task master.'''

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

if __name__ == '__main__':
    # set parent folder here (use abs path)!
    parent_folder = '/work/djb208/predicts_data1'
    # obtain list of folders within parent_folder
    worklist = os.listdir(parent_folder)
    worklist = [os.path.join(parent_folder, e) for e in worklist]
    worklist = [e for e in worklist if os.path.isdir(e)]
    if rank == 0:
        master(worklist)
    else:
        worker(generate_trees)
