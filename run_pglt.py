#! /bin/usr/env python
## D.J. Bennett
## 26/05/2014
"""
pG-lt is a pipeline for the automated generation of phylogenies through
'Mass Phylogeny Estimation'. This program is built on top of 
phyloGenerator (C) 2013 and was written by D.J. Bennett with
additional help from W.D. Pearse and L. Hudson. This program makes
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

## Messages
priming_msg = '\nERROR: The program was unable to start mpe due to a \
problem with the files and folders in the study directory. Check the \
parameters and gene parameters .csv for any potential conflicts.'
toofewspecies_msg = '\nERROR: The program halted as there are too few \
species left of phylogeny building -- five is the minimum. You may \
have started with too few names, or names given could not be \
taxonomically resolved or there may be too little sequence data \
available.'
taxonomicrank_msg =  '\nERROR: It is likely that one or more names \
have been resolved incorrectly, as such the parent taxonomic group \
has been set to Eukaryotes which is too high a taxonomic rank for \
phylogenetic analysis. Consider adding a parent ID to the \
parameters.csv to prevent incorrect names resolution or reducing the \
taxonomic diversity of the analysis names.'
outgroup_msg = '\nERROR: The outgroup has been dropped. This may be \
due to too few sequence data available for outgroup or a failure to \
align sequences that are available. If outgroup has been \
automatically selected, consider manually choosing an outgroup.'
raxml_msg = '\nERROR: Generated maxtrys poor phylogenies \
consecutively, consider reducing maxrttsd.'
unexpected_msg = '\nERROR: The following unexpected error occurred:\n\
\"{0}\" \n\
Please email details to the program maintainer for help.'

## Packages
import sys
from pglt import _PARS as default_pars_file
from pglt import _GPARS as default_gpars_file
from pglt.tools.system_tools import Stager
from pglt.tools.system_tools import TooFewSpeciesError
from pglt.tools.system_tools import PrimingError
from pglt.tools.system_tools import TaxonomicRankError
from pglt.tools.system_tools import OutgroupError
from pglt.tools.system_tools import RAxMLError
from pglt.tools.setup_tools import printHeader
from pglt.tools.setup_tools import parseArgs
from pglt.tools.setup_tools import getDirs
from pglt.tools.setup_tools import setUpLogging
from pglt.tools.setup_tools import tearDownLogging
from pglt.tools.setup_tools import logMessage
from pglt.tools.setup_tools import logError
from pglt.tools.init_tools import prime
from pglt.tools.init_tools import sortArgs
from pglt.tools.special_tools import clean

def main():
	"""Run pG-lt with user defined arguments from command-line"""
	# read arguments
	parser = parseArgs()
	args = parser.parse_args()
	if args.clean:
		clean()
		sys.exit('Files and folders deleted')
	if not args.email:
		# stop if no email
		print 'An email address must be provided. Use \'-e\'.'
		sys.exit()
	# print program header
	printHeader()
	# create base logger -- messages in parent folder log.txt
	base_logger = setUpLogging(args.verbose, args.debug)
	# search cwd for folders that contain names and parameter files
	dirs = getDirs(base_logger)
	logMessage('begin', logger = base_logger, directory = dirs)
	# loop through each folder
	for i in range(len(dirs)):
		if not args.verbose:
			print 'Woking on [{0}]'.format(dirs[i])
		logMessage('start', logger = base_logger, directory = dirs[i])
		error_raised = False
		# set up a root logger, so now default logging refers to this
		folder_logger = setUpLogging(args.verbose, args.debug,\
			logname = '', directory = dirs[i])
		try:
			# get list of arguments
			arguments = sortArgs(dirs[i], args.email, folder_logger,\
				default_pars_file, default_gpars_file)
			# initialise hidden files
			stage = prime(dirs[i], arguments)
			# run Stager
			Stager.run_all(dirs[i], stage = stage, verbose = \
				args.verbose)
		# if error raised handle it accordingly ...
		except PrimingError:
			error_raised = logError(priming_msg, folder_logger)
		except TooFewSpeciesError:
			error_raised = logError(toofewspecies_msg, folder_logger)
		except TaxonomicRankError:
			error_raised = logError(taxonomicrank_msg, folder_logger)
		except OutgroupError:
			error_raised = logError(outgroup_msg, folder_logger)
		except RAxMLError:
			error_raised = logError(raxml_msg, folder_logger)
		except KeyboardInterrupt:
			folder_logger.info('Execution halted by user')
			sys.exit('Execution halted by user')
		except Exception as unexpected_error:
			error_raised = logError(unexpected_msg.format(\
				unexpected_error), folder_logger)
		# disable folder_logger
		tearDownLogging('')
		if error_raised:
			logMessage('folder-error', logger = base_logger,\
				directory = dirs[i])
		else:
			logMessage('finish', logger = base_logger, directory =\
				dirs[i])
	logMessage('end', logger = base_logger)

if __name__ == '__main__':
	main()