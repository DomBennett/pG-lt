#! /bin/usr/env python
## D.J. Bennett
## 24/03/2014
"""
MPE system tools
"""

## Packages
import subprocess,threading,sys,os,platform
from datetime import datetime
from mpe import stages

class StageError(Exception):
	pass

class Tee(object):
	#http://stackoverflow.com/questions/616645/how-do-i-duplicate-sys-stdout-to-a-log-file-in-python/616686#616686
	def __init__(self, path):
		self.file = open(path, 'w')
		self.stdout,sys.stdout = sys.stdout,self
		self.stderr,sys.stderr = sys.stderr,self

	def __enter__(self):
		pass

	def __exit__(self, type, value, traceback):
		sys.stdout,sys.stderr = self.stdout,self.stderr
		self.file.close()

	def write(self, data):
		self.file.write(data)
		self.file.flush()
		self.stdout.write(data)
		self.stdout.flush()

class Stager(object):
	"""Stager class : runs each file in stage folder.\
Adapted from code written by L. Hudson."""
	STAGES = stages.STAGES

	def __init__(self, stage):
		if stage not in self.STAGES:
			raise StageError('Stage [{0}] not recognised'.format(stage))
		else:
			self.stage = stage
			output_dir = os.path.split(self.STAGES[stage])[1]
			self.output_dir = os.path.splitext(output_dir)[0]

	def _prime(self):
		# Creates self.output_dir. An error is raised if it already exists.
		if not os.path.isdir(self.output_dir):
			os.mkdir(self.output_dir)
		return Tee(os.path.join(self.output_dir, 'log.txt'))

	def _time_string(self):
		return datetime.today().strftime("%A, %d %B %Y %I:%M%p")

	def _cmd(self, arg):
		s = subprocess.Popen(['python', '-u', arg], stdout=subprocess.PIPE,
			stderr=subprocess.STDOUT)
		while s.poll() is None:
			line = s.stdout.readline()
			if len(line) > 1:
				print line[:-1]
			else:
				print line
		if s.poll():
			sys.exit(s.poll())

	def run(self):
		with self._prime():
			print 'Stage [{0}] started at [{1}]'.format(self.stage,
			self._time_string())
			print 'Running on [{0}] [{1}]'.format(platform.node(),
				platform.platform())
			pyfile = self.STAGES[self.stage]
			self._cmd('--version')
			self._cmd(pyfile)
			print 'Stage finished at [{0}]'.format(self._time_string())

	@classmethod
	def run_all(klass, stage):
		for s in sorted(Stager.STAGES.keys()[stage:]):
			Stager(s).run()


class TerminationPipe(object):
	"""Adapted pG object : exectute programs."""
	#Background process class
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