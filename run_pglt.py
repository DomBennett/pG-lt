#! /bin/usr/env python
# D.J. Bennett
# 26/05/2014

# PACKAGES
import os
import sys
import pickle
from pglt import __version__ as pglt_version
from pglt import __doc__ as pglt_doc
from pglt import __year__ as pglt_year
import pglt.tools.setup_tools as stools
from pglt.tools.system_tools import Runner


# GLOBALS
stools.pglt_version = pglt_version
stools.pglt_doc = pglt_doc
stools.pglt_year = pglt_year


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
        stools.logMessage('program-restart', logger=base_logger, retry=retry)
        with open(argspath, 'r') as file:
            nworkers, threads_per_worker, folders, email, threads, verbose,\
                debug, stages = pickle.load(file)
        runner = Runner(folders=folders, nworkers=nworkers, stages=stages,
                        threads_per_worker=threads_per_worker, wd=os.getcwd(),
                        email=email, verbose=verbose, debug=debug,
                        logger=base_logger, retry=retry)
    else:
        if verbose:
            stools.printHeader()
        # search cwd for folders that contain names and parameter files
        folders = stools.getFolders()
        # calculate nworkers
        nworkers, threads_per_worker, spare_threads =\
            stools.calcWorkers(threads=threads, nfolders=len(folders))
        # start message
        stools.logMessage('program-start', logger=base_logger, folders=folders,
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
    restart, retry, email, threads, verbose, debug, stages = \
        stools.parseArguments()
    # create base logger -- messages in parent folder log.txt
    base_logger = stools.setUpLogging(verbose, debug, 'base')
    try:
        main(restart, retry, email, threads, verbose, debug, stages,
             base_logger)
    except KeyboardInterrupt:
        base_logger.info('Execution halted by user')
        if verbose:
            sys.exit()
        else:
            sys.exit('Execution halted by user')
