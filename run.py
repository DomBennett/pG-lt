#!/usr/bin/python
## MPE: Run stage
## D.J. Bennett
## 25/03/2014
#TODO: Run for multiple names files, check parameters

import argparse,pickle,csv,os,platform,socket,stat,sys,subprocess
from datetime import datetime
import stages

## Objects
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

class Stage(object):
	STAGES = stages.STAGES

	def __init__(self, stage):
		if stage not in self.STAGES:
			raise StageError('Stage [{0}] not recognised'.format(stage))
		else:
			self.stage = stage
			self.output_dir = os.path.splitext(self.STAGES[stage][-1])[0]

	def _prime(self):
		# Creates self.output_dir. An error is raised if it already exists.
		if os.path.isdir(self.output_dir):
			raise StageError('Output directory [' + self.output_dir + '] ' + 
							 'already exists. Has this stage already been run?')
		else:
			os.mkdir(self.output_dir)
			return Tee(os.path.join(self.output_dir, 'log.txt'))

	def _time_string(self):
		return datetime.today().strftime("%A, %d %B %Y %I:%M%p")

	def _finished(self):
		# Make output files read only
		for root,dirs,files in os.walk(self.output_dir):
			for f in files:
				os.chmod(os.path.join(root, f), stat.S_IREAD)

	def _cmd(self, args):
		s = subprocess.Popen(args, stdout=subprocess.PIPE, 
							 stderr=subprocess.STDOUT)
		while True:
			line = s.stdout.readline()
			exitcode = s.poll()
			line = line[:-1]
			if (not line) and (exitcode is not None):
				break
			elif line:
				print line

		if exitcode:
			sys.exit(exitcode)

	def run(self):
		with self._prime():
			print 'Stage [{0}] started at [{1}]'.format(self.stage, 
														self._time_string())
			print 'Running on [{0}] [{1}]'.format(platform.node(), 
												  platform.platform())
			args = self.STAGES[self.stage]
			self._cmd( [args[0], '--version'] )
			self._cmd(args)
			print 'Stage finished at [{0}]'.format(self._time_string())

		self._finished()

	@classmethod
	def run_all(klass):
		print 'Running all stages'
		for s in sorted(Stage.STAGES.keys()):
			Stage(s).run()

def main():
	# Check args
	args = parser.parse_args()
	if args.parameters:
		parameters = args.parameters
	else:
		parameters = "parameters.csv"
	if args.genes:
		genes = args.genes
	else:
		genes = "gene_parameters.csv"
	# Input
	try:
		terms = []
		with open(args.names) as names:
			for name in names:
				terms.append(name.strip())
		terms = [term for term in terms if not term == '']
	except IOError:
		print "Names file could not be opened. File: [{0}]".args.names
		sys.exit()
	try:
	 	genedict = {}
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
	try:
	 	paradict = {}
	 	with open(parameters, 'rb') as csvfile:
	 		reader = csv.DictReader(csvfile)
	 		for row in reader:
	 			paradict[row["Parameter"]] = row["Value"]
	except IOError:
		print "Parameters file could not be opened. File: [{0}]".parameters
	 	sys.exit()
	# Build programdict
	programdict = {}
	programdict['raxml'] = os.path.join('requires', 'raxml')
	programdict['mafft'] = 'mafft'
	# Output
 	with open("genedict.p", "wb") as file:
 		pickle.dump(genedict, file)
 	with open("paradict.p", "wb") as file:
 		pickle.dump(paradict, file)
 	with open("terms.p", "wb") as file:
 		pickle.dump(terms, file)
 	with open("programdict.p", "wb") as file:
 		pickle.dump(programdict, file)
	Stage.run_all()

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Mass Phylogeny Estimation (MPE) - an automated pipeline\
		for the generation of phylogenies from taxonomic names (Author: D.J. Bennett).")
	parser.add_argument("-names", "-n", help=".txt file of taxonomic names.")
	parser.add_argument("-parameters", "-p", help=".csv file of parameters.")
	parser.add_argument("-genes", "-g", help=".csv file of gene parameters")
	main()