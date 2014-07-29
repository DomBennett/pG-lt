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
import argparse,pickle,sys,csv,os,re,logging,platform
from datetime import datetime
from mpe.tools.system import Stager
from mpe.tools.system import TooFewSpeciesError
from mpe import _PARS as default_pars_file
from mpe import _GPARS as default_gpars_file

## Informative messages
ioerror_msg = "[{0}] file could not be opened in [{1}]. Check that \
it is not opened by another program"
nonamestxt_msg = 'No folders containing \'names.txt\' files found! \
All taxonomic names should be placed in subdirectories and called: \
\'names.txt\''
description = """MPE D.J. Bennett (C) 2014

A pipeline for the automated generation of phylogenies from taxonomic
names through Mass Phylogeny Estimation."""

## Classes
class PrimingError(Exception):
	pass

## Root functions
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

## Folder functions
def recordPars(paradict):
	"""Return mpe parameters string"""
	record = '\nUsing the following parameters:\n'
	for key in paradict.keys():
		record += '    [{0}] = [{1}]\n'.format(key, paradict[key])
	return record

def recordGpars(genedict):
	"""Return gene parameters string"""
	record = '\nUsing the following genes and gene parameters:\n'
	for gene in genedict.keys():
		record += '  Gene: [{0}]\n'.format(gene)
		for par in genedict[gene]:
			record += '    [{0}] = [{1}]\n'.format(par,\
				genedict[gene][par])
	return record

def prime(directory, arguments):
	"""Write pickle files, print arguments"""
	# Write pickle files
	with open(os.path.join(directory, ".genedict.p"),\
		"wb") as file:
		pickle.dump(arguments['genedict'], file)
	with open(os.path.join(directory,".paradict.p"),\
		"wb") as file:
		pickle.dump(arguments['paradict'], file)
	with open(os.path.join(directory,".terms.p"),\
		"wb") as file:
		pickle.dump(arguments['terms'], file)
	# Print arguments and parameters to file
	record = 'Working with [{0}] names\n'.format(len(arguments['terms']))
	record += recordPars(arguments['paradict'])
	record += recordGpars(arguments['genedict'])
	with open(os.path.join(directory, 'info.txt'),'wd') as file:
		file.write(record)

def readInNames(directory):
	"""Read names from text file in dir"""
	terms = []
	with open(os.path.join(directory, 'names.txt')) as names:
		for name in names:
			terms.append(name.strip())
	terms = [term for term in terms if not term == '']
	return terms

def readInGenePars(gpars_file):
	"""Read gene_parameters.csv. Return list of dictionaries."""
	def split(element):
		# split elements by ':'
		if ':' in element:
			return element.split(':')
		return element
	def _read(gpars_file, template, genes = None):
		# open csv file and replace parameters in template
		#  if they are None. If genes specified, only read
		#  rows for those genes.
		with open(gpars_file, 'rb') as csvfile:
			reader = csv.DictReader(csvfile)
			for row in reader:
				if genes:
					if not row['gene'] in genes:
						continue
				temp = template.copy()
				for key in temp.keys():
					if row[key]:
						if temp[key] is None:
							temp[key] = split(row[key])
				genedict[row['gene']] = temp
		return genedict
	# check if file exists, else use default
	if not os.path.isfile(gpars_file):
		return readInGenePars(default_gpars_file)
	# genedicts
	genedict = {}
	# template of dict in genedict
	template = {'names' : None, 'taxid' : None,\
	'mingaps' : None, 'minoverlap' : None, 'minfails'\
	: None, 'maxtrys' : None, 'minseedsize' : None,\
	'maxseedtrys': None}
	# open file, read each row and fill in template
	genedict = _read(gpars_file, template)
	# if Nones, use defaults
	nones = False
	for gene in genedict.keys():
		for par in genedict[gene].keys():
			if genedict[gene][par] is None:
				nones = True
				break
	if nones:
		# run _read for defaults and limit to genes in genedict
		genedict = _read(default_gpars_file, template, genedict.keys())
	return genedict

def readInPars(pars_file):
	"""Read gene_parameters.csv. Return dictionary."""
	def _read(pars_file, paradict):
		# open csv, and replace all Nones
		with open(pars_file, 'rb') as csvfile:
			reader = csv.DictReader(csvfile)
			for row in reader:
				if paradict[row["Parameter"]] is None:
					paradict[row["Parameter"]] = row["Value"]
		return paradict
	# check if file exists, else use default
	if not os.path.isfile(pars_file):
		return readInPars(default_pars_file)
	# template
	paradict = {'nseqs' : None, 'naligns' : None,\
	'ntrees' : None, 'download_thoroughness' : None,\
	'maxlen' : None, 'minlen' : None, 'maxtrys' : None,\
	'maxrttsd': None, 'parentid' : None}
	# open file, read each row, extract value
	paradict = _read(pars_file, paradict)
	# if Nones remain, use default
	nones = False
	for key in paradict.keys():
		if paradict[key] is None:
			nones = True
			break
	if nones:
		paradict = _read(default_pars_file, paradict)
	return paradict

def sortArgs(directory, email, logger):
	"""Search for relevant files in dir, return list of arguments"""
	# find text file and read, raise error if fail
	try:
		terms = readInNames(directory)
	except IOError:
		logger.error(ioerror_msg.format('names.txt', directory))
		raise PrimingError()
	# find gene parameter file and read, raise error if fail
	try:
		genedict = readInGenePars(os.path.join(directory, 'gene_parameters.csv'))
	except IOError:
		logger.error(ioerror_msg.format('gene_parameters.csv', directory))
		raise PrimingError()
	# find parameter file and read, raise error if fail
	try:
		paradict = readInPars(os.path.join(directory, 'parameters.csv'))
	except IOError:
		logger.error(ioerror_msg.format('parameters.csv', directory))
		raise PrimingError()
	# add email to paradict
	paradict['email'] = email
	return {'terms' : terms, 'genedict' : genedict, 'paradict' : paradict}

def main():
	"""Run mpe with user defined arguments from command-line"""
	# print program header
	printHeader()
	# read arguments
	parser = parseArgs()
	args = parser.parse_args()
	if not args.email:
		# stop if no email
		print 'An email address must be provided. Use \'-e\'.'
		sys.exit()
	# create base logger -- messages in parent folder log.txt
	base_logger = setUpLogging(args.verbose, args.debug)
	# search cwd for folders that contain names and parameter files
	dirs = getDirs(base_logger)
	logMessage('begin', logger = base_logger, directory = dirs)
	# loop through each folder
	for i in range(len(dirs)):
		# set up a root logger, so now default logging refers to this
		logger = setUpLogging(args.verbose, args.debug, logname = '',\
			directory = dirs[i])
		if not args.verbose:
			print 'Woking on [{0}]'.format(dirs[i])
		logMessage('start', logger = base_logger, directory = dirs[i])
		error_raised = False
		try:
			# get list of arguments
			arguments = sortArgs(dirs[i], args.email, logger)
			# initialise hidden files
			prime(dirs[i], arguments)
			# run Stager
			Stager.run_all(dirs[i], stage = 0, verbose = args.verbose)
		except PrimingError:
			error_raised = True
			logger.error('A priming error occurred for [{0}]! Check log \
file for more details. Moving to next folder'.format(dirs[i]))
		except TooFewSpeciesError:
			error_raised = True
			logger.error('Too few species left for [{0}]. Moving to \
next folder'.format(dirs[i]))
		except KeyboardInterrupt:
			error_raised = True
			logger.info('Execution halted by user')
			sys.exit('Execution halted by user')
		except:
			error_raised = True
			logger.exception('Unexpected error occurred for [{0}]. \
Moving to next folder'.format(dirs[i]))
		if error_raised:
			logMessage('folder-error', logger = base_logger, directory = dirs[i])
		else:
			logMessage('finish', logger = base_logger, directory = dirs[i])
	logMessage('end', logger = base_logger)

if __name__ == '__main__':
	main()