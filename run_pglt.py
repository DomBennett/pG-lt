#! /bin/usr/env python
# D.J. Bennett
# 26/05/2014
"""
pG-lt is a pipeline for the automated generation of phylogenies through
'Mass Phylogeny Estimation'. This program is built on top of
phyloGenerator (C) 2013 and was written by D.J. Bennett with
additional help from W.D. Pearse and L. Hudson. It makes
use of external programs for phylogeny generation and bioinformatics
these are: RAxML (Copyright (C) Stamatakis 2013) , MAFFT (Copyright
(C) 2013 Kazutaka Katoh) the NCBI's standalone BLAST suite 2.2.29+ and
online API services (Copyright NCBI (C) 2009). It also uses the
following python packages: Biopython (Copyright Cook (C) 2009),
Dendropy (Copyright Sukumaran and Holder (C) 2010) and Taxon Names
Resovler (Copyright (C) Bennett 2014).

Copyright (C) 2014  Dominic John Bennett

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
"""

# PACAKGES
import os
import sys
import argparse
from pglt.tools.setup_tools import setUpLogging
from pglt.tools.setup_tools import printHeader
from pglt.tools.setup_tools import getFolders
from pglt.tools.setup_tools import logMessage
from pglt.tools.system_tools import Runner
from pglt.tools.special_tools import getThreads
from pglt.tools.special_tools import clean


# FUNCTIONS
def parseArguments():
    """Read command-line arguments"""
    # Add new arg for threads
    parser = argparse.ArgumentParser()
    parser.add_argument("-email", "-e", help="please provide email \
for NCBI")
    parser.add_argument("-threads", "-t", help="number of threads, default -1 \
will use all available on machine.", default=-1, type=int)
    parser.add_argument("-restart", "-r", help="restart from \
specified stage")
    parser.add_argument("--verbose", help="increase output verbosity",
                        action="store_true")
    parser.add_argument("--debug", help="log warnings (developer only)",
                        action="store_true")
    parser.add_argument("--clean", help="remove all pG-lt files and \
folders (developer only)", action="store_true")
    # get args
    args = parser.parse_args()
    # check them
    if args.clean:
        clean()
        sys.exit('Files and folders deleted')
    if not args.email:
        # stop if no email
        print 'An email address must be provided. Use \'-e\'.'
        sys.exit()
    # check threads is a valid argument
    if args.threads == 0 or args.threads < -1:
        sys.exit('Invalid threads argument.')
    return args.email, args.threads, args.verbose, args.debug


def calcWorkers(threads, nfolders, min_threads_per_worker=2,
                max_threads_per_worker=100):
    """Calculate the number of workers for parallel running of folders"""
    # get available threads on machine
    available_threads = getThreads()
    # make sure threads arg is not greater than those available
    if threads > available_threads:
        sys.exit('More threads specified than avaiable on machine')
    if threads == -1:
        threads = available_threads
    # calc min_threads_per_worker if it is greater than threads
    if min_threads_per_worker > threads:
        min_threads_per_worker = threads
    # calc max_threads_per_worker if it is greater than threads
    if max_threads_per_worker > threads:
        max_threads_per_worker = threads
    # calc nworkers and threads_per_worker
    # increase workers before threads_per_worker
    threads_per_worker = min_threads_per_worker
    for i in range(nfolders):
        if (float(i)*threads_per_worker) > threads:
            nworkers = i-1
            break
    else:
        nworkers = nfolders
        for i in range(min_threads_per_worker, max_threads_per_worker):
            if (float(nworkers)*i) > threads:
                threads_per_worker = i-1
                break
        else:
            threads_per_worker = max_threads_per_worker
    spare_threads = int(threads - (float(nworkers)*threads_per_worker))
    return nworkers, threads_per_worker, spare_threads


def start():
    """Start pG-lt start!"""
    email, threads, verbose, debug = parseArguments()
    if verbose:
        printHeader()
    # create base logger -- messages in parent folder log.txt
    base_logger = setUpLogging(not verbose, debug, 'base')
    # search cwd for folders that contain names and parameter files
    folders = getFolders()
    return email, threads, verbose, debug, base_logger, folders


def run(email, threads, verbose, debug, base_logger, folders):
    """Run pG-lt run!"""
    # calculate nworkers
    nworkers, threads_per_worker, spare_threads =\
        calcWorkers(threads=threads, nfolders=len(folders))
    # start message
    logMessage('program-start', logger=base_logger, folders=folders,
               threads=threads_per_worker*nworkers,
               spare_threads=spare_threads, email=email)
    # setup runner
    runner = Runner(folders=folders, nworkers=nworkers,
                    threads_per_worker=threads_per_worker, wd=os.getcwd(),
                    email=email, verbose=verbose, debug=debug)
    # logMessage('start', logger=base_logger, directory=dirs[i])
    runner.setup(folders)
    # RUN STAGES 1 & 2
    # don't run these in parallel
    logMessage('stage-start', logger=base_logger, stage='1')
    runner.run(folders=runner.folders, stage='1', parallel=False)
    logMessage('stage-end', logger=base_logger)
    logMessage('stage-start', logger=base_logger, stage='2')
    runner.run(folders=runner.folders, stage='2', parallel=False)
    logMessage('stage-end', logger=base_logger)
    # RUN STAGES 3 & 4
    # do run these in parallel
    logMessage('stage-start', logger=base_logger, stage='3')
    runner.run(folders=runner.folders, stage='3', parallel=True)
    logMessage('stage-end', logger=base_logger)
    logMessage('stage-start', logger=base_logger, stage='4')
    runner.run(folders=runner.folders, stage='4', parallel=True)
    logMessage('stage-end', logger=base_logger)
    # finish message
    logMessage('program-end', logger=base_logger)
    sys.exit()

# MAIN
if __name__ == '__main__':
    email, threads, verbose, debug, base_logger, folders = start()
    try:
        run(email, threads, verbose, debug, base_logger, folders)
    except KeyboardInterrupt:
        base_logger.info('Execution halted by user')
        if not verbose:
            sys.exit()
        else:
            sys.exit('Execution halted by user')
