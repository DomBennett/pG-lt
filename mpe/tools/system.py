#! /bin/usr/env python
## D.J. Bennett
## 24/03/2014
"""
MPE system tools
"""

## Packages
import subprocess,threading,sys,os,platform,csv,pickle,logging
from datetime import datetime
from mpe import _PARS as default_pars
from mpe import _GPARS as default_gpars

## Classes
class LogFilter(object):
	#http://stackoverflow.com/questions/8162419/python-logging-specific-level-only
    def __init__(self, level):
        self.__level = level

    def filter(self, logRecord):
        return logRecord.levelno <= self.__level

class StageError(Exception):
	pass

class Stager(object):
	"""Stager class : runs each file in stage folder. Adapted from\
 code written by L. Hudson."""

	def __init__(self, stage, verbose = True, dev = False):
		self.verbose = verbose
		self.dev = dev
		if stage not in self.STAGES:
			raise StageError('Stage [{0}] not recognised'.format(stage))
		else:
			self.stage = stage
			# dir is second element of tuple
			self.output_dir = self.STAGES[stage][1]

	def _prime(self):
		if not os.path.isdir(self.output_dir):
			os.mkdir(self.output_dir)
		# root logger
		root_log = logging.getLogger()
		root_log.setLevel(logging.INFO)
		# add file hander to root
		logfile = os.path.join(self.output_dir, 'log.txt')
		loghandler = logging.FileHandler(logfile, 'a')
		loghandler.setFormatter(logging.Formatter('%(message)s'))
		root_log.addHandler(loghandler)
		if self.verbose:
			# if verbose
			console = logging.StreamHandler()
			console.setFormatter(logging.Formatter('%(message)s'))
			root_log.addHandler(console)
		# if not dev, filter out everything but info
		if not self.dev:
			root_log.addFilter(LogFilter(logging.INFO))

	def _second(self):
		root_log = logging.getLogger()
		root_handlers = root_log.handlers[:]
		# remove all handlers
		for each in root_handlers:
			root_log.removeHandler(each)

	def _start(self):
		logging.info('\n' + '#' * 70)
		logging.info('Stage [{0}] started at [{1}]'.format(self.stage,
			self._time_string()))
		logging.info('Running on [{0}] [{1}]'.format(platform.node(),
					platform.platform()))
		logging.info('Python [{0}]'.format(sys.version))
		logging.info('#' * 70 + '\n')

	def _end(self):
		logging.info('\n' + '#' * 70)
		logging.info('Stage [{0}] finished at [{1}]'.format(self.stage,
			self._time_string()))
		logging.info('#' * 70 + '\n')

	def _time_string(self):
		return datetime.today().strftime("%A, %d %B %Y %I:%M%p")

	def _cmd(self):
		# function is first element of tuple
		self.STAGES[self.stage][0]()

	def run(self):
		self._prime() # add handlers to root logger
		self._start() # log system info
		self._cmd() # run stage
		self._end() # log end time
		self._second() # remove handlers from root logger

	@classmethod
	def run_all(klass, stage, verbose = True, dev = False):
		for s in sorted(Stager.STAGES.keys()[stage:]):
			if not verbose:
				print 'Stage [{0}]'.format(Stager.STAGES[s][1])
			Stager(s, verbose, dev).run()

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

## Functions
def readArgs(args):
	"""Safely read arguments write pickle files, return stage and verbosity"""
	if args.restart:
		try:
			stage = int(args.restart) - 1
		except TypeError:
			print "Restart must be numeric!"
			sys.exit()
		print "Restarting from stage: [{0}]\n".format(args.restart)
		return stage,args.verbose,args.development
	if args.names:
		try:
			terms = []
			with open(args.names) as names:
				for name in names:
					terms.append(name.strip())
			terms = [term for term in terms if not term == '']
		except IOError:
			print "Names file could not be opened. File: \
[{0}]".format(args.names)
			sys.exit()
	else:
		print "No names file given! Type 'MPE.py -h' for help."
		sys.exit()
	if args.parameters:
		parameters = args.parameters
	else:
		parameters = default_pars
	if args.genes:
		genes = args.genes
	else:
		genes = default_gpars
	genedict = {}
	try:
		with open(genes, 'rb') as csvfile:
			reader = csv.DictReader(csvfile)
			for row in reader:
				temp = {}
				temp["names"] = row["names"].split(":")
				temp["taxid"] = row["taxid"]
				temp["mingaps"] = row["mingaps"]
				temp["minoverlap"] = row["minoverlap"]
				temp["minfails"] = row["minfails"]
				temp["maxtrys"] = row["maxtrys"]
				temp["minseedsize"] = row["minseedsize"]
				temp["maxseedtrys"] = row["maxseedtrys"]
				genedict[row['gene']] = temp
	except IOError:
		print "Gene parameters file could not be opened. File: [{0}]".genes
		sys.exit()
	except KeyError:
		print "Unexpected headers in gene parameters."
		sys.exit()
	paradict = {}
	try:
		with open(parameters, 'rb') as csvfile:
			reader = csv.DictReader(csvfile)
			for row in reader:
				paradict[row["Parameter"]] = row["Value"]
	except IOError:
		print "Parameters file could not be opened. File: [{0}]".parameters
		sys.exit()
	# Output
	with open(".genedict.p", "wb") as file:
		pickle.dump(genedict, file)
	with open(".paradict.p", "wb") as file:
		pickle.dump(paradict, file)
	with open(".terms.p", "wb") as file:
		pickle.dump(terms, file)
	# Print
	if args.parameters:
		parameters = args.parameters
	else:
		parameters = 'DEFAULT'
	if args.genes:
		genes = args.genes
	else:
		genes = 'DEFAULT'
	print 'Running with names [{0}], parameters [{1}] and gene \
parameters [{2}]'.format(args.names, parameters, genes)
	return 0,args.verbose,args.development