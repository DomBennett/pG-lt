#!/usr/bin/python
## MPE System tools
## D.J. Bennett
## 24/03/2014

## Pacakges
import subprocess,threading,sys

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