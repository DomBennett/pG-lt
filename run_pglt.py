#! /bin/usr/env python
# D.J. Bennett
# 26/05/2014

# PACKAGES
import os
import sys
import pickle
from pglt.tools.setup_tools import setUpLogging
from pglt.tools.setup_tools import printHeader
from pglt.tools.setup_tools import getFolders
from pglt.tools.setup_tools import calcWorkers
from pglt.tools.setup_tools import parseArguments
from pglt.tools.setup_tools import logMessage
from pglt.tools.system_tools import Runner


def main(restart, retry, email, threads, verbose, debug, stages, base_logger):
    '''Run pG-lt'''
    # setup folder to hold pickled args
    temp_dir = os.path.join(os.getcwd(), 'tempfiles')
    if not os.path.isdir(temp_dir):
        os.mkdir(temp_dir)
    argspath = os.path.join(temp_dir, 'arguments.p')
    if restart:
        if not os.path.isfile(argspath):
            sys.exit('Cannot restart, are you sure you have already run \
pG-lt?')
        logMessage('program-restart', logger=base_logger)
        with open(argspath, 'r') as file:
            nworkers, threads_per_worker, folders, email, threads, verbose,\
                debug, stages = pickle.load(file)
        runner = Runner(folders=folders, nworkers=nworkers, stages=stages,
                        threads_per_worker=threads_per_worker, wd=os.getcwd(),
                        email=email, verbose=verbose, debug=debug,
                        logger=base_logger, retry=retry)
    else:
        if verbose:
            printHeader()
        # search cwd for folders that contain names and parameter files
        folders = getFolders()
        # calculate nworkers
        nworkers, threads_per_worker, spare_threads =\
            calcWorkers(threads=threads, nfolders=len(folders))
        # start message
        logMessage('program-start', logger=base_logger, folders=folders,
                   threads=threads_per_worker*nworkers, stages=stages,
                   spare_threads=spare_threads, email=email)
        # setup runner
        base_logger.info('Setting up files and folders ....')
        runner = Runner(folders=folders, nworkers=nworkers, stages=stages,
                        threads_per_worker=threads_per_worker, wd=os.getcwd(),
                        email=email, verbose=verbose, debug=debug,
                        logger=base_logger, retry=False)
        runner.setup()
        base_logger.info('Done.')
        # pickle runner for restart
        with open(argspath, 'w') as file:
            pickle.dump((nworkers, threads_per_worker, folders, email,
                         threads, verbose, debug, stages), file)
    # run stages
    runner.run()


if __name__ == '__main__':
    # parse args
    restart, retry, email, threads, verbose, debug, stages = parseArguments()
    # create base logger -- messages in parent folder log.txt
    base_logger = setUpLogging(verbose, debug, 'base')
    try:
        main(restart, retry, email, threads, verbose, debug, stages,
             base_logger)
    except KeyboardInterrupt:
        base_logger.info('Execution halted by user')
        if verbose:
            sys.exit()
        else:
            sys.exit('Execution halted by user')
