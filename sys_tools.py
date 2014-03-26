#!/usr/bin/python
## MPE System tools
## D.J. Bennet
## 24/03/2014

## Pacakges
import subprocess,threading,os,platform,socket,stat,sys
from datetime import datetime
#import stages

## Objects
#class StageError(Exception):
#	pass

#class Tee(object):
	# http://stackoverflow.com/questions/616645/how-do-i-duplicate-sys-stdout-to-a-log-file-in-python/616686#616686
#	def __init__(self, path):
#		self.file = open(path, 'w')
#		self.stdout,sys.stdout = sys.stdout,self
#		self.stderr,sys.stderr = sys.stderr,self
#
#	def __enter__(self):
#		pass
#
#	def __exit__(self, type, value, traceback):
#		sys.stdout,sys.stderr = self.stdout,self.stderr
#		self.file.close()
#
#	def write(self, data):
#		self.file.write(data)
#		self.file.flush()
#		self.stdout.write(data)
#		self.stdout.flush()

# class Stage(object):
	# STAGES = stages.STAGES

	# def __init__(self, stage):
	# 	if stage not in self.STAGES:
	# 		raise StageError('Stage [{0}] not recognised'.format(stage))
	# 	else:
	# 		self.stage = stage
	# 		self.output_dir = os.path.splitext(self.STAGES[stage][-1])[0]

	# def _prime(self):
	# 	# Creates self.output_dir. An error is raised if it already exists.
	# 	if os.path.isdir(self.output_dir):
	# 		raise StageError('Output directory [' + self.output_dir + '] ' + 
	# 						 'already exists. Has this stage already been run?')
	# 	else:
	# 		os.mkdir(self.output_dir)
	# 		return Tee(os.path.join(self.output_dir, 'log.txt'))

	# def _time_string(self):
	# 	return datetime.today().strftime("%A, %d %B %Y %I:%M%p")

	# def _finished(self):
	# 	# Make output files read only
	# 	for root,dirs,files in os.walk(self.output_dir):
	# 		for f in files:
	# 			os.chmod(os.path.join(root, f), stat.S_IREAD)

	# def _cmd(self, args):
	# 	s = subprocess.Popen(args, stdout=subprocess.PIPE, 
	# 						 stderr=subprocess.STDOUT)
	# 	while True:
	# 		line = s.stdout.readline()
	# 		exitcode = s.poll()
	# 		line = line[:-1]
	# 		if (not line) and (exitcode is not None):
	# 			break
	# 		elif line:
	# 			print line

	# 	if exitcode:
	# 		sys.exit(exitcode)

	# def run(self):
	# 	with self._prime():
	# 		print 'Stage [{0}] started at [{1}]'.format(self.stage, 
	# 													self._time_string())
	# 		print 'Running on [{0}] [{1}]'.format(platform.node(), 
	# 											  platform.platform())
	# 		args = self.STAGES[self.stage]
	# 		self._cmd( [args[0], '--version'] )
	# 		self._cmd(args)
	# 		print 'Stage finished at [{0}]'.format(self._time_string())

	# 	self._finished()

	# @classmethod
	# def run_all(klass):
	# 	print 'Running all stages'
	# 	for s in sorted(Stage.STAGES.keys()):
	# 		Stage(s).run()

class TerminationPipe(object):
	"""Adapted pG object : exectute programs."""
	#Background process class
	def __init__(self, cmd, timeout, silent=True):
		self.cmd = cmd
		self.timeout = timeout
		self.process = None
		self.output = None
		self.failure = False
		self.stderr = 'EMPTY'
		self.stdout = 'EMPTY'
		self.silent = silent
	
	def run(self, silent=None, changeDir=False):
		def silentTarget():
			if sys.platform == 'win32':
				if changeDir:
					self.process = subprocess.Popen("requires\\" + self.cmd, stdout=subprocess.PIPE, shell=True, stderr=subprocess.PIPE)
				else:
					self.process = subprocess.Popen(self.cmd, stdout=subprocess.PIPE, shell=True, stderr=subprocess.PIPE)
			else:
				if changeDir:
					self.process = subprocess.Popen("./requires/" + self.cmd, stdout=subprocess.PIPE,shell=True, stderr=subprocess.PIPE)
				else:
					self.process = subprocess.Popen(self.cmd, stdout=subprocess.PIPE,shell=True, stderr=subprocess.PIPE)
			self.output = self.process.communicate()
		
		def loudTarget():
			if sys.platform == 'win32':
				if changeDir:
					self.process = subprocess.Popen("requires\\" + self.cmd, shell=False)
				else:
					self.process = subprocess.Popen(self.cmd, shell=False)
			else:
				if changeDir:
					self.process = subprocess.Popen("./requires/" + self.cmd, shell=False)
				else:
					self.process = subprocess.Popen(self.cmd, shell=False)
			self.output=self.process.communicate()
		
		if silent: self.silent = silent
		if self.silent:
			thread = threading.Thread(target=silentTarget)
		else:
			thread = threading.Thread(target=loudTarget)
		thread.start()
		thread.join(self.timeout)
		if thread.is_alive():
			self.process.terminate()
			thread.join()
			self.failure = True