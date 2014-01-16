import sys, os, re, random
import numpy as np
sys.path.append(os.path.join(os.getcwd(), 'functions'))
import phyloGenerator_adapted as pG


class OutgroupError(Exception):
	"""Raised whenever an outgroup is dropped"""
	pass

class SeqObj(dict):
	"""Store species' gene sequences with functions for pulling sequences for alignments and adding penalties for sequences that did not align"""
	def __init__(self, genedir, seqfiles, outgroup, minfails):
		self.minfails = minfails # minimum number of fails in a row
		self.attempts = 0 # number of attempts
		self.dspp = [] # species dropped
		self.nseqs = 0 # counter for seqs
		for i, seqfile in enumerate(seqfiles):
			name = re.sub('\.fasta$', '', seqfile)
			if name != outgroup:
				name = "tx" + name
			else:
				name = "outgroup"
			seqdir = os.path.join(genedir, seqfile)
			seqs = []
			lengths = []
			with open(seqdir, "rU") as infile:
				for record in pG.SeqIO.parse(infile, "fasta"):
					if record.seq.count('N') == 0 and record.seq.count('n') == 0: # removing sequences w/ Ns TODO: move this to download
						if 200 < len(record) < 1000:
							record.id = name
							lengths.append(len(record))
							seqs.append([record, 0]) # seqrecord + nfails
							self.nseqs += 1
			if len(seqs) > 0:
				self[name] = (seqs, np.median(lengths))
			
	def start(self):
		"""Return starting random sp sequences, update spp_pool"""
		self.attempts = 1
		self.spp_pool = self.keys()
		self.spp_pool = [e for e in self.spp_pool if e != "outgroup"]
		self.align_obj = []
		for i in range(3):
			self._add()
		self.spp_pool.append("outgroup")
		return self.align_obj
		
	#def restart(self, index):
	#	""""""
	#	self.align_obj[0][1] += 1
	#	self._check()
	#	return self.start()
		
	def back(self, align):
		"""Add to nfails for random species, return new random species"""
		self.attempts += 1
		self.align_obj[-1][1] += 1
		del self.align_obj[-1]
		self.spp_pool.append(self.next_sp)
		self._check()
		if len(self.spp_pool) != 0:
			self._add()
		return self.align_obj
			
	def next(self, align):
		"""Set nfails to 0, add additional random species"""
		for i in range(len(self.align_obj)):
			self.align_obj[i][1] = 0
		self._add()
		self.attempts = 1
		return self.align_obj
		
	def _check(self):
		"""Check nfails, drop sequences and species"""
		spp = self.keys()
		for sp in spp:
			to_drop = []
			for i in range(len(self[sp][0])):
				if self[sp][0][i][1] > self.minfails:
					print "Dropping [{0}] for [{1}] as nfails is [{2}]".format(i,sp, self[sp][0][i][1])
					to_drop.append(i)
			for i in to_drop:
				del self[sp][0][i]
			if len(self[sp][0]) < 1:
				if sp == "outgroup":
					raise OutgroupError
				else:
					print "Dropped [{0}]".format(sp)
					self.dspp.append(sp)
					self.spp_pool = [e for e in self.spp_pool if e != sp]
					del self[sp]
					
	def _add(self):
		"""Add new random species' sequence"""
		rand_int = random.randint(0, (len(self.spp_pool)-1))
		self.next_sp = self.spp_pool.pop(rand_int)
		next_seq = random.sample(self[self.next_sp][0], 1)[0]
		self.align_obj.append(next_seq)

def cleanAlignment(align, method='trimAl-automated', tempStem='temp', timeout=None):
	if 'trimAl' in method:
		options = ""
		if 'automated' in method:
			options += " -automated1"
		output = []
		for i,gene in enumerate(align):
			geneOutput = []
			for j,method in enumerate(gene):
				pG.AlignIO.write(method, tempStem+"Input.fasta", "fasta")
				fileLine = " -in " + tempStem + "Input.fasta -out " + tempStem + "Output.fasta -fasta"
				trimalVersion = "trimal"
				commandLine = trimalVersion + fileLine + options
				if timeout:
                                        pipe = pG.TerminationPipe(commandLine, timeout)
                                        if 'darwin' == sys.platform: # so that I can run it on mac
                                                pipe.run(changeDir = True)
                                        else:
                                                pipe.run()
					if not pipe.failure:
						geneOutput.append(pG.AlignIO.read(tempStem + "Output.fasta", "fasta"))
						os.remove(tempStem + "Output.fasta")
						os.remove(tempStem + "Input.fasta")
					else:
						raise RuntimeError("Either trimAl failed, or it ran out of time")
			output.append(geneOutput)
		return output
	else:
		raise RuntimeError("Only automated trimAl methods supported at this time.")

					
def incrAlign(seqobj, max_pgap):
	"""Incrementally build an alignment from outgroup sequence"""
	# parameters
	max_trys = 100 # prevents loops continuing forever
	def runAlignment(align_obj):
		align_struct = [[e[0]] for e in align_obj]
		align = pG.alignSequences(align_struct, method= 'mafft', nGenes = 1)
		align = cleanAlignment(align, timeout = 99999)[0][0] #trim with trimal
		al = align.get_alignment_length()
		if al == 0:
			return align,[False],0
		else:
			ngaps = [e.seq.count('-') for e in align]
			pgaps = [float(e)/al for e in ngaps]
			pgaps_bool = [e < max_pgap for e in pgaps]
			return align,pgaps_bool,al
	def returnBestAlign(alignments):
		# test if there are any alignments
		if len(alignments) < 1:
			return None
		# keep only alignments with lots of records
		nrecords = [len(e._records) for e in alignments]
		alignments = [alignments[i] for i,e in enumerate(nrecords)\
				      if e == max(nrecords)]
		# return longest alignment
		alignments_lens = [len(e) for e in alignments]
		max_i = alignments_lens.index(max(alignments_lens))
		print "Returning alignment of length: [{0}]".\
		    format(alignments[max_i].get_alignment_length())
		return alignments[max_i]
	# run alignments until alignment for all species made
	trys = 0
	min_align_len = min([seqobj[e][1] for e in seqobj.keys()])
	# voting in numbers, the more sequences in the starting alignment the better
	while True:
		min_align_len = min([seqobj[e][1] for e in seqobj.keys()])
		align_obj = seqobj.start()
		align,pgaps_bool,al = runAlignment(align_obj)
		if all(pgaps_bool) and al > min_align_len:
			break
		else:
			for i in range(len(align_obj)):
				setaside = align_obj.pop(i)
				# drop each sequence, align again....
				align,pgaps_bool,al = runAlignment(align_obj)
				if all(pgaps_bool) and al > min_align_len:
					setaside[1] += 1
					seqobj._check()
					if len(seqobj.keys()) < 4:
						# TODO: research what this section does (16/01/2014)
						# 4 for now.... random issue if 3 otherwise
						print "Fewer than 4 species"
						return None
					break
				# add the setside where it was before
				align_obj.insert(i,setaside)
			trys += 1
		if max_trys < trys:
			print "Max trys hit"
			return None
	trys = 0
	align_store = [] # stores successful aligns generated
	while True:
		min_align_len = min([seqobj[e][1] for e in seqobj.keys()])
		align,pgaps_bool,al = runAlignment(align_obj)
		if all(pgaps_bool) and al > min_align_len:
			counter = 0
			if len(seqobj.spp_pool) == 0:
				return align
			else:
				align_store.append(align)
				align_obj = seqobj.next(align)
		elif trys < max_trys:
			try:
				align_obj = seqobj.back(align)
			except OutgroupError:
				print "Outgroup error"
				return returnBestAlign(align_store)
			if len(seqobj.spp_pool) == 0:
				# here a species has been dropped and now all species are present
				return returnBestAlign(align_store)
			trys += 1
		else:
			# when the maximum number of species is not reached...
			# ... return the best alignment in the alignment store
			return returnBestAlign(align_store)
