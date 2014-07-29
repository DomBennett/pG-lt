#! /bin/usr/env python
## D.J. Bennett
## 24/03/2014
"""
mpe system tools
"""

## Packages
import subprocess,threading,os,logging
from datetime import datetime

## Classes
class StageError(Exception):
	pass

class TooFewSpeciesError(Exception):
	pass

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
		logging.info('\n' + '#'*70)
		logging.info('Stage [{0}] started at [{1}]\n'.format(self.stage,
			self._time_string()))

	def _end(self):
		logging.info('\nStage [{0}] finished at [{1}]'.format(self.stage,
			self._time_string()))
		logging.info('#'*70)

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