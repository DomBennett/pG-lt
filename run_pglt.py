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
from pglt.tools.setup_tools import setUpLogging
from pglt.tools.setup_tools import printHeader
from pglt.tools.setup_tools import getFolders
from pglt.tools.setup_tools import calcWorkers
from pglt.tools.setup_tools import parseArguments
from pglt.tools.setup_tools import logMessage
from pglt.tools.system_tools import Runner


def start():
    """Start pG-lt start!"""
    email, threads, verbose, debug, stages = parseArguments()
    if verbose:
        printHeader()
    # create base logger -- messages in parent folder log.txt
    base_logger = setUpLogging(verbose, debug, 'base')
    # search cwd for folders that contain names and parameter files
    folders = getFolders()
    return email, threads, verbose, debug, base_logger, folders, stages


def run(email, threads, verbose, debug, base_logger, folders, stages):
    """Run pG-lt run!"""
    # calculate nworkers
    nworkers, threads_per_worker, spare_threads =\
        calcWorkers(threads=threads, nfolders=len(folders))
    # start message
    logMessage('program-start', logger=base_logger, folders=folders,
               threads=threads_per_worker*nworkers, stages=stages,
               spare_threads=spare_threads, email=email)
    # logMessage('start', logger=base_logger, directory=dirs[i])
    base_logger.info('Setting up files and folders ....')
    # setup runner
    runner = Runner(folders=folders, nworkers=nworkers,
                    threads_per_worker=threads_per_worker, wd=os.getcwd(),
                    email=email, verbose=verbose, debug=debug,
                    logger=base_logger)
    runner.setup(folders=folders)
    base_logger.info('Done.')
    # run stages
    for stage in stages:
        if stage in ['3', '4']:
            parallel = True
        else:
            parallel = False
        logMessage('stage-start', logger=base_logger, stage=stage)
        runner.run(folders=runner.folders, stage=stage, parallel=parallel)
        logMessage('stage-end', logger=base_logger, stage=stage)
    # finish
    logMessage('program-end', logger=base_logger)
    sys.exit('Exiting ....')

if __name__ == '__main__':
    email, threads, verbose, debug, base_logger, folders, stages = start()
    try:
        run(email=email, threads=threads, verbose=verbose, debug=debug,
            base_logger=base_logger, folders=folders, stages=stages)
    except KeyboardInterrupt:
        base_logger.info('Execution halted by user')
        if verbose:
            sys.exit()
        else:
            sys.exit('Execution halted by user')
