#! /bin/usr/env python
## D.J. Bennett
## 26/05/2014
##TODO: bring back restart
"""
mpe is a pipeline for the automated generation of phylogenies through 'Mass
Phylogeny Estimation'. This program is built on top of phyloGenerator (C) 2013
and was written by D.J. Bennett with additional help from W.D. Pearse and L. Hudson.
This program makes use of external programs for phylogeny generation and bioinformatics
these are: RAxML (Copyright (C) Stamatakis 2013) , MAFFT (Copyright (C) 2013 Kazutaka
Katoh) the NCBI's standalone BLAST suite 2.2.29+ and online API services
 (Copyright NCBI (C) 2009). It also uses the following python packages:
 Biopython (Copyright Cook (C) 2009), Dendropy (Copyright Sukumaran and Holder (C)
 2010) and Taxon Names Resovler (Copyright (C) Bennett 2014).

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

## Packages
import argparse,sys,os,re,logging,platform
from datetime import datetime
from mpe.tools.system import clean
from mpe.tools.system import Stager
from mpe.tools.system import sortArgs
from mpe.tools.system import prime
from mpe.tools.system import TooFewSpeciesError
from mpe.tools.system import PrimingError
from mpe.tools.system import TaxonomicRankError
from mpe.tools.system import OutgroupError

## Description
description = """MPE D.J. Bennett (C) 2014

A pipeline for the automated generation of phylogenies from taxonomic
names through Mass Phylogeny Estimation."""

## Error messages
nonamestxt_msg = 'No folders containing \'names.txt\' files found! \
All taxonomic names should be placed in subdirectories and called: \
\'names.txt\''
priming_msg = 'The program was unable to start mpe due to a problem \
with the files and folders in the study directory. Check the parameters \
and gene parameters .csv for any potential conflicts.'
toofewspecies_msg = 'The program halted as there are too few species left of \
phylogeny building -- five is the minimum. You may have started with too few names, \
or names given could not be taxonomically \
resolved or there may be too little sequence data available.'
taxonomicrank_msg =  'It is likely that one or more names have \
been resolved incorrectly, as such the parent taxonomic \
group has been set to Eukaryotes which is too high a \
taxonomic rank for phylogenetic analysis. Consider \
adding a parent ID to the parameters.csv to prevent \
incorrect names resolution or reducing the taxonomic diversity \
of the analysis names.'
outgroup_msg = 'The outgroup has been dropped. This may be due to too few \
sequence data available for outgroup or a failure to align sequences that are \
available. If outgroup has been automatically selected, consider manually choosing \
and outgroup.'
unexpected_msg = 'An unexpected error occurred. Unfortunately, this is not a professional \
program and as such these errors are likely. Please email details to the program \
maintainer for help.'

## Functions
def printHeader():
	"""Print a nice program description header"""
	# use 70 cols as I think this is standard
	print '\n' + '#' * 70
	print description
	print '#' * 70 + '\n'

def parseArgs():
	"""Read command-line arguments"""
	parser = argparse.ArgumentParser()
	parser.add_argument("-email", "-e", help="please provide email for NCBI")
	parser.add_argument("-restart", "-r", help="restart from specified stage")
	parser.add_argument("--verbose", help="increase output verbosity",
					action="store_true")
	parser.add_argument("--debug", help="log warnings (developer only)",
					action="store_true")
	parser.add_argument("--clean", help="remove all mpe files and folders (developer only)",
					action="store_true")
	return parser

def getDirs(logger):
	"""Return folders in directory with names.txt files"""
	# list all folders
	unchecked_dirs = [f for f in os.listdir('.') if not\
	os.path.isfile(f)]
	# remove hidden folders
	unchecked_dirs = [d for d in unchecked_dirs if not\
	re.match('^\.', d)]
	# loop through each and check they contain a names.txt
	checked_dirs = []
	for each in unchecked_dirs:
		path = os.path.join (os.getcwd(), each)
		files = os.listdir(path)
		if 'names.txt' in files:
			checked_dirs.append(each)
	if len(checked_dirs) > 0:
		return checked_dirs
	else:
		logger.error(nonamestxt_msg)
		sys.exit()

def setUpLogging(verbose, debug, logname = 'base',\
	directory = os.getcwd()):
	"""Set up logging : direct and control log statements"""
	# get logger
	logger = logging.getLogger(logname)
	# if this logger already has handlers, remove them -- prevents propogation
	if logger.handlers:
		handlers = logger.handlers[:]
		for h in handlers:
			logger.removeHandler(h)
	if debug:
		# log all statements above DEBUG level
		logger.setLevel(logging.DEBUG)
	else:
		# log all statements above INFO level
		# (which is higher than DEBUG)
		logger.setLevel(logging.INFO)
	# add file hander to root
	logfile = os.path.join(directory, 'log.txt')
	loghandler = logging.FileHandler(logfile, 'a')
	# set statement format -- I only want the message
	loghandler.setFormatter(logging.Formatter('%(message)s'))
	logger.addHandler(loghandler)
	if verbose:
		# if verbose, copy all info statements to console
		console = logging.StreamHandler()
		console.setFormatter(logging.Formatter('%(message)s'))
		logger.addHandler(console)
	return logger

def logMessage(phase, logger, directory = None):
	if phase == 'begin':
		# begin running the program
		# directory is list
		logger.info('\n' + '#' * 70)
		logger.info(description)
		logger.info('#' * 70 + '\n')
		logger.info('Running on [{0}] [{1}]'.format(platform.node(),
					platform.platform()))
		logger.info('Python [{0}]\n'.format(sys.version))
		logger.info('Working with the following directories:')
		# convert dirs to string
		dir_string = ''
		chars_counter = 0
		for each in directory[:-1]:
			chars_counter += len(each)
			if chars_counter > 70:
				# stop at 70 columns
				dir_string += '[' + each + '],\n'
				chars_counter = 0
			else:
				dir_string += '[' + each + '], '
		dir_string += '[' + directory[-1] + ']'
		logger.info('[{0}]'.format(dir_string))
	elif phase == 'start':
		# start for one folder
		# directory is a string
		logger.info('\n' + '#' * 70 + '\n')
		logger.info('Folder [{0}] started at [{1}]'.format(directory,
			datetime.today().strftime("%A, %d %B %Y %I:%M%p")))
	elif phase == 'finish':
		# when a folder is finished running
		logger.info('Folder [{0}] finished at [{1}]'.format(directory,\
			datetime.today().strftime("%A, %d %B %Y %I:%M%p")))
	elif phase == 'folder-error':
		# when a folder is unable to run
		logger.info('Unfinished for folder [{0}] at [{1}]'.format(directory,\
			datetime.today().strftime("%A, %d %B %Y %I:%M%p")))
		logger.info('Check [{0}] for more details'.format(os.path.join(\
			directory, 'log.txt')))
	elif phase == 'end':
		logger.info('\nFIN')
	else:
		raise(ValueError('Unrecognised phase'))

def logError(msg, logger):
	"""Return true when error raised, log informative message"""
	logger.error(msg)
	logger.info('\n\n Moving to the next folder')
	return True

def main():
	"""Run mpe with user defined arguments from command-line"""
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
		# set up a root logger, so now default logging refers to this
		folder_logger = setUpLogging(args.verbose, args.debug, logname = '',\
			directory = dirs[i])
		if not args.verbose:
			print 'Woking on [{0}]'.format(dirs[i])
		logMessage('start', logger = base_logger, directory = dirs[i])
		error_raised = False
		try:
			# get list of arguments
			arguments = sortArgs(dirs[i], args.email, folder_logger)
			# initialise hidden files
			prime(dirs[i], arguments)
			# run Stager
			Stager.run_all(dirs[i], stage = 0, verbose = args.verbose)
		# if error raised handle it accordingly ...
		except PrimingError:
			error_raised = logError(priming_msg, folder_logger)
		except TooFewSpeciesError:
			error_raised = logError(toofewspecies_msg, folder_logger)
		except TaxonomicRankError:
			error_raised = logError(taxonomicrank_msg, folder_logger)
		except OutgroupError:
			error_raised = logError(outgroup_msg, folder_logger)
		except KeyboardInterrupt:
			folder_logger.info('Execution halted by user')
			sys.exit('Execution halted by user')
		except:
			error_raised = logError(unexpected_msg, folder_logger)
		if error_raised:
			logMessage('folder-error', logger = base_logger, directory = dirs[i])
		else:
			logMessage('finish', logger = base_logger, directory = dirs[i])
	logMessage('end', logger = base_logger)

if __name__ == '__main__':
	main()