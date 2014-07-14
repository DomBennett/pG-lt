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

## Functions
def printHeader():
	"""Print a nice program description header"""
	# use 70 cols as I think this is standard
	print '\n' + '#' * 70
	print description
	print '#' * 70 + '\n'

def logPars(paradict):
	"""Print all mpe parameters"""
	logging.info('\nUsing the following parameters:')
	for key in paradict.keys():
		logging.info('    [{0}] = [{1}]'.format(key, paradict[key]))

def logGpars(genedict):
	"""Print gene parameters"""
	logging.info('\nUsing the following genes and gene parameters:')
	for gene in genedict.keys():
		logging.info('  Gene: [{0}]'.format(gene))
		for par in genedict[gene]:
			logging.info('    [{0}] = [{1}]'.format(par,\
				genedict[gene][par]))

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
	# Print arguments
	logging.info('Working with [{0}] names'.format(len(arguments['terms'])))
	logPars(arguments['paradict'])
	logGpars(arguments['genedict'])

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
	'maxpedge': None, 'parentid' : None}
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

def sortArgs(directory, email):
	"""Search for relevant files in dir, return list of arguments"""
	# find text file and read, raise error if fail
	try:
		terms = readInNames(directory)
	except IOError:
		logging.error(ioerror_msg.format('names.txt', directory))
		raise PrimingError()
	# find gene parameter file and read, raise error if fail
	try:
		genedict = readInGenePars(os.path.join(directory, 'gene_parameters.csv'))
	except IOError:
		logging.error(ioerror_msg.format('gene_parameters.csv', directory))
		raise PrimingError()
	# find parameter file and read, raise error if fail
	try:
		paradict = readInPars(os.path.join(directory, 'parameters.csv'))
	except IOError:
		logging.error(ioerror_msg.format('parameters.csv', directory))
		raise PrimingError()
	# add email to paradict
	paradict['email'] = email
	return {'terms' : terms, 'genedict' : genedict, 'paradict' : paradict}

def getDirs():
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
		logging.error(nonamestxt_msg)
		sys.exit()

def setUpLogging(verbose, debug):
	"""Set up logging : direct and control log statements"""
	# root logger
	root_log = logging.getLogger()
	if debug:
		# log all statements above DEBUG level
		root_log.setLevel(logging.DEBUG)
	else:
		# log all statements above INFO level
		# (which is higher than DEBUG)
		root_log.setLevel(logging.INFO)
	# add file hander to root
	logfile = os.path.join(os.getcwd(), 'log.txt')
	loghandler = logging.FileHandler(logfile, 'a')
	# set statement format -- I only want the message
	loghandler.setFormatter(logging.Formatter('%(message)s'))
	root_log.addHandler(loghandler)
	if verbose:
		# if verbose, copy all info statements to console
		console = logging.StreamHandler()
		console.setFormatter(logging.Formatter('%(message)s'))
		root_log.addHandler(console)

def logStartFolder(directory):
	"""Log start of new folder"""
	logging.info('\n' + '#' * 70 + '\n')
	logging.info('Folder [{0}] started at [{1}]'.format(directory,
		datetime.today().strftime("%A, %d %B %Y %I:%M%p")))

def logEndFolder(directory):
	"""Log end of folder"""
	logging.info('Folder [{0}] finished at [{1}]'.format(directory,
		datetime.today().strftime("%A, %d %B %Y %I:%M%p")))

def logStart(directories):
	"""Log start"""
	logging.info('\n' + '#' * 70)
	logging.info(description)
	logging.info('#' * 70 + '\n')
	logging.info('Running on [{0}] [{1}]'.format(platform.node(),
				platform.platform()))
	logging.info('Python [{0}]\n'.format(sys.version))
	logging.info('Working with the following directories:')
	# convert dirs to string
	print_dirs = ''
	chars_counter = 0
	for each in directories[:-1]:
		chars_counter += len(each)
		if chars_counter > 70:
			# stop at 70 columns
			print_dirs += each + ',\n'
			chars_counter = 0
		else:
			print_dirs += each + ', '
	print_dirs += directories[-1]
	logging.info('[{0}]'.format(print_dirs))

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
	# set up logging
	setUpLogging(args.verbose, args.debug)
	# log exceptions if they occur:
	# search cwd for folders that contain names and parameter files
	dirs = getDirs()
	logStart(dirs)
	# loop through each folder
	for each in dirs:
		if not args.verbose:
			print 'Woking on [{0}]'.format(each)
		logStartFolder(each)
		try:
			# get list of arguments
			arguments = sortArgs(each, args.email)
			# initialise hidden files
			prime(each, arguments)
			# run Stager
			Stager.run_all(each, stage = 0, verbose = args.verbose)
		except PrimingError:
			logging.error('A priming error occurred for [{0}]! Check log \
file for more details. Moving to next folder'.format(each))
		except TooFewSpeciesError:
			logging.error('Too few species left for [{0}]. Moving to \
next folder'.format(each))
		except KeyboardInterrupt:
			logging.info('Execution halted by user')
			sys.exit('Execution halted by user')
		except:
			logging.exception('Unexpected error occurred for [{0}]. \
Moving to next folder'.format(each))
		logEndFolder(each)

if __name__ == '__main__':
	main()