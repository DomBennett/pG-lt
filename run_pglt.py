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


def main(restart, email, threads, verbose, debug, stages, base_logger):
    '''Run pG-lt'''
    # setup folder to hold pickled runner
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
                        logger=base_logger)
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
                        logger=base_logger)
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
    restart, email, threads, verbose, debug, stages = parseArguments()
    # create base logger -- messages in parent folder log.txt
    base_logger = setUpLogging(verbose, debug, 'base')
    try:
        main(restart, email, threads, verbose, debug, stages, base_logger)
    except KeyboardInterrupt:
        base_logger.info('Execution halted by user')
        if verbose:
            sys.exit()
        else:
            sys.exit('Execution halted by user')
