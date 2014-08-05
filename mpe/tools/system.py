#! /bin/usr/env python
## D.J. Bennett
## 24/03/2014
"""
mpe system tools
"""

## Packages
import subprocess,threading,os,logging,shutil,pickle,csv
from mpe import _PARS as default_pars_file
from mpe import _GPARS as default_gpars_file
from datetime import datetime

## Error messages
ioerror_msg = "[{0}] file could not be opened in [{1}]. Check that \
it is not opened by another program"

## Error classes
class StageError(Exception):
	pass

class TooFewSpeciesError(Exception):
	pass

class TaxonomicRankError(Exception):
	pass

class PrimingError(Exception):
	pass

class OutgroupError(Exception):
	pass

## Other classes
class Stager(object):
	"""Stager class : runs each file in stage folder. Adapted from\
 code written by L. Hudson."""

	def __init__(self, wd, stage):
		if stage not in self.STAGES:
			raise StageError('Stage [{0}] not recognised'.format(stage))
		else:
			self.wd = wd
			self.stage = stage
			# dir is second element of tuple
			self.output_dir = os.path.join(wd, self.STAGES[stage][1])

	def _start(self):
		logging.info('-' * 70)
		logging.info('Stage [{0}] started at [{1}]'.format(self.stage,
			self._time_string()))
		logging.info('-' * 70)

	def _end(self):
		logging.info('-' * 70)
		logging.info('Stage [{0}] finished at [{1}]'.format(self.stage,
			self._time_string()))
		logging.info('-' * 70 + '\n\n')

	def _time_string(self):
		return datetime.today().strftime("%A, %d %B %Y %I:%M%p")

	def _cmd(self):
		# function is first element of tuple
		# pass working directory to it
		self.STAGES[self.stage][0](self.wd)

	def run(self):
		if not os.path.isdir(self.output_dir):
			os.mkdir(self.output_dir)
		self._start() # log system info
		self._cmd() # run stage
		self._end() # log end time

	@classmethod
	def run_all(klass, wd, stage, verbose):
		for s in sorted(Stager.STAGES.keys()[stage:]):
			if not verbose:
				print '  Stage [{0}]'.format(Stager.STAGES[s][1])
			Stager(wd, s).run()

class TerminationPipe(object):
	"""Adapted pG object : exectute background programs."""
	def __init__(self, cmd, timeout = 99999, silent = True):
		self.cmd = cmd
		self.timeout = timeout
		self.process = None
		self.output = None
		self.failure = False
		self.stderr = 'EMPTY'
		self.stdout = 'EMPTY'
		self.silent = silent
	
	def run(self):
		def silentTarget():
			self.process = subprocess.Popen(self.cmd, stdout = subprocess.PIPE,\
				shell = True, stderr = subprocess.PIPE)
			self.output = self.process.communicate()
		def loudTarget():
			self.process = subprocess.Popen(self.cmd, shell = False)
			self.output = self.process.communicate()
		if self.silent:
			thread = threading.Thread(target = silentTarget)
		else:
			thread = threading.Thread(target = loudTarget)
		thread.start()
		thread.join(self.timeout)
		if thread.is_alive():
			self.process.terminate()
			thread.join()
			self.failure = True

## Priming functions
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
	'maxseedtrys': None, 'type': None}
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
	'maxrttsd': None, 'parentid' : None, 'outgroupid' : None}
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

# Special functions
def clean():
	"""Remove all files and folders created by mpe"""
	## Remove log.txt in parent folder
	os.remove('log.txt')

	## Go through all subdirs and remove files in list
	files_to_remove = ['info.txt', 'log.txt','.paradict.p','.allrankids.p',\
	'.namesdict.p', '.genedict.p', '.terms.p', '.constraint.tre', '.partitions.txt',\
	'.phylogeny_in.phylip', '.alignment_out.fasta']
	folders = os.listdir(os.getcwd())
	while folders:
		temp_folder = folders.pop()
		temp_files_to_remove = files_to_remove[:]
		while temp_files_to_remove:
			temp_file_to_remove = temp_files_to_remove.pop()
			try:
				os.remove(os.path.join(temp_folder, temp_file_to_remove))
			except OSError:
				pass

	## Go through all subdirs and remove folders in list
	folders_to_remove = ['1_names', '2_download', '3_alignment', '4_phylogeny']
	folders = os.listdir(os.getcwd())
	while folders:
		temp_folder = folders.pop()
		temp_folders_to_remove = folders_to_remove[:]
		while temp_folders_to_remove:
			temp_folder_to_remove = temp_folders_to_remove.pop()
			try:
				shutil.rmtree(os.path.join(temp_folder,temp_folder_to_remove))
			except OSError:
				pass