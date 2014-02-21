import sys, os, re, random
import numpy as np
sys.path.append(os.path.join(os.getcwd(), 'functions'))
import phyloGenerator_adapted as pG


class OutgroupError(Exception):
	"""Raised whenever an outgroup is dropped"""
	pass

class MinSpeciesError(Exception):
	"""Raised whenever there are too few species for phylogeny generation"""
	pass

class SeqObj(dict):
	"""Store species' gene sequences with functions for pulling sequences for alignments and adding penalties for sequences that did not align"""
	def __init__(self, genedir, seqfiles, outgroup, minfails):
		self.minfails = minfails # minimum number of fails in a row
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
				self[name] = [seqs, np.median(lengths)]

	def start(self, n):
		"""Return n starting random sp sequences, update spp_pool"""
		self.spp_pool = self.keys()
		self.align_obj = []
		for i in range(n):
			self._add()
		return self.align_obj

	def _add(self):
		"""Add new random species' sequence"""
		rand_int = random.randint(0, (len(self.spp_pool)-1))
		self.next_sp = self.spp_pool.pop(rand_int)
		next_seq = random.sample(self[self.next_sp][0], 1)[0]
		self.align_obj.append(next_seq)
		
	def back(self, align):
		"""Add to nfails for random species, return new random species"""
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
		return self.align_obj
		
	def _check(self):
		"""Check nfails, drop sequences and species"""
		spp = self.keys()
		for sp in spp:
			to_drop = []
			for i in range(len(self[sp][0])):
				if self[sp][0][i][1] > self.minfails:
					print "Dropping [{0}] for [{1}] as nfails is [{2}]".\
					    format(i,sp, self[sp][0][i][1])
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
		if len(self.keys()) < 5: # Raise error if fewer than 5 species left in spp pool
			raise MinSpeciesError

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
					
def incrAlign(seqobj, pintgapmax, pextgapmax, method, nstart):
	"""Incrementally build an alignment from outgroup sequence"""
	# parameters
	# nstart the number of sequences to use in start
	max_nstart_trys = 20 #  the number of trys before nstart -= 1
	min_nstart = 5 # the minimum nstart
	max_trys = 100 # prevents loops continuing forever
	def runAlignment(align_obj):
		align_struct = [[e[0]] for e in align_obj]
		seq_descriptions = [e[0].description for e in align_obj]
		align = pG.alignSequences(align_struct, method = method, nGenes = 1)
		align = cleanAlignment(align, timeout = 99999)[0][0] #trim with trimal
		# ensure sequence names used are in description
		for i,e in enumerate(seq_descriptions):
			align._records[i].description = e
		return align
	def alignCheck(align):
		if align.get_alignment_length() < min_align_len:
			return False
		for each in align:
			sequence = each.seq.tostring()
			totgaps = sequence.count('-')
			totnucs = len(sequence) - totgaps
			extgaps = 0
			start_extgaps = re.search("^-+", sequence)
			end_extgaps = re.search("-+$", sequence)
			if not start_extgaps is None:
				extgaps += start_extgaps.end()
			if not end_extgaps is None:
				extgaps += len(sequence) - end_extgaps.start()
			intgaps = totgaps - extgaps
			overlap = intgaps + totnucs
			if float(intgaps)/overlap > pintgapmax:
				return False
			if float(extgaps)/len(sequence) > pextgapmax:
				return False
		return True
	def returnBestAlign(alignments):
		# keep alignments with outgroups
		# keep alignments with more than 5 species
		alignments = [al for al in alignments if "outgroup" in\
				      [e.id for e in al._records]]
		alignments = [al for al in alignments if len(al) >= 5]
		if len(alignments) == 0:
			return None
		# keep only alignments with lots of records
		nrecords = [len(e._records) for e in alignments]
		alignments = [alignments[i] for i,e in enumerate(nrecords)\
				      if e == max(nrecords)]
		# return longest alignment
		alignments_lens = [len(e) for e in alignments]
		max_i = alignments_lens.index(max(alignments_lens))
		return alignments[max_i]
	## run alignments until alignment for all species made
	trys = 0
	nstart_trys = 0
	align_store = [] # stores successful aligns generated
	# start with as many sequences at start as possible
	while True:
		min_align_len = min([seqobj[e][1] for e in seqobj.keys()])
		align_obj = seqobj.start(nstart)
		align = runAlignment(align_obj)
		if alignCheck(align):
			if len(seqobj.spp_pool) == 0:
				return align, nstart
			else:
				align_store.append(align)
				align_obj = seqobj.next(align)
			break
		trys += 1
		# drop nstart by 1 if max_nstart_trys
		if nstart > min_nstart:
			nstart_trys += 1
			if max_nstart_trys < nstart_trys:
				nstart -= 1
		if max_trys < trys:
			return None, nstart
	trys = 0
	while True:
		min_align_len = min([seqobj[e][1] for e in seqobj.keys()])
		align = runAlignment(align_obj)
		if alignCheck(align):
			trys = 0
			if len(seqobj.spp_pool) == 0:
				return align, nstart
			else:
				align_store.append(align)
				align_obj = seqobj.next(align)
		elif trys < max_trys:
			align_obj = seqobj.back(align)
			if len(seqobj.spp_pool) == 0:
				# here a species has been dropped and now all species are present
				return returnBestAlign(align_store), nstart
			trys += 1
		else:
			# when the maximum number of species is not reached...
			# ... return the best alignment in the alignment store
			return returnBestAlign(align_store), nstart
