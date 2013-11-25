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

					
def incrAlign(seqobj, max_pgap):
	"""Incrementally build an alignment from outgroup sequence"""
	def runAlignment(align_obj):
		align_struct = [[e[0]] for e in align_obj]
		align = pG.alignSequences(align_struct, method= 'mafft', nGenes = 1)
		align = pG.cleanAlignment(align, timeout = 99999)[0][0] #trim with trimal
		al = align.get_alignment_length()
		if al == 0:
			return align,[False],0
		else:
			ngaps = [e.seq.count('-') for e in align]
			pgaps = [float(e)/al for e in ngaps]
			pgaps_bool = [e < max_pgap for e in pgaps]
			return align,pgaps_bool,al
	# run alignments until alignment for all species made
	counter = 0 # prevents loops continuing forever
	# print [seqobj[e][1] for e in seqobj.keys()]
	min_align_len = min([seqobj[e][1] for e in seqobj.keys()])
	print "Min length: [{0}]".format(min_align_len)
	while True: # voting in numbers, the more sequences in the starting alignment the better
		min_align_len = min([seqobj[e][1] for e in seqobj.keys()])
		align_obj = seqobj.start()
		print len(align_obj)
		align,pgaps_bool,al = runAlignment(align_obj)
		if all(pgaps_bool) and al > min_align_len:
			break
		else:
			for i in range(len(align_obj)):
				setaside = align_obj.pop(i) # drop each sequence, align again....
				print len(align_obj)
				align,pgaps_bool,al = runAlignment(align_obj)
				if all(pgaps_bool) and al > min_align_len:
					setaside[1] += 1
					seqobj._check()
					if len(seqobj.keys()) < 4: # 4 for now.... random issue if 3 otherwise
						raise RuntimeError("Failed to align at start")
					break
				align_obj.insert(i,setaside) # add the setside where it was before
			counter += 1
		if counter > 100:
			raise RuntimeError("Failed to align at start")
	counter = 0
	while True:
		min_align_len = min([seqobj[e][1] for e in seqobj.keys()])
		align,pgaps_bool,al = runAlignment(align_obj)
		print len(align)
		print align_obj[-1][0].id
		if all(pgaps_bool) and al > min_align_len:
			counter = 0
			if len(seqobj.spp_pool) == 0:
				return align
			else:
				align_obj = seqobj.next(align)
		else:
			align_obj = seqobj.back(align)
			if len(seqobj.spp_pool) == 0: # here a species has been dropped and now all species are
				return align
			counter += 1
		if counter > 100:
			raise RuntimeError("Failed to align")
				#alignment_store.append(align)
			#if seqobj.attempts > max_ats:
			#	align_obj = seqobj.restart()
			#	alignment_store.append(previous_align)
			#	counter += 1
			#	print "Counter is now [{0}]".format(counter)
		#previous_align = align
	# when the maximum number of species is not reached...
	# ... return the best alignment in the alignment store
	#alignment_store_lens = [len(e) for e in alignment_store]
	#max_i = alignment_store_lens.index(max(alignment_store_lens))
	#print "#"
	#print "Returning alignment of length: [{0}]".format(alignment_store[max_i].get_alignment_length())
	#return alignment_store[max_i]
