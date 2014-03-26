#!/usr/bin/python
## MPE Alignment tools
## D.J. Bennet
## 24/03/2014

## Packages
import sys, os, re, random
import numpy as np
from Bio import SeqIO
from Bio import AlignIO
from sys_tools import *

## Objects
class OutgroupError(Exception):
	"""Raised whenever an outgroup is dropped"""
	pass

class MinSpeciesError(Exception):
	"""Raised whenever there are too few species for phylogeny generation"""
	pass

class SeqObj(dict):
	"""Store species' gene sequences with functions for pulling sequences for \
	alignments and adding penalties for sequences that did not align"""
	def __init__(self, genedir, seqfiles, minfails):
		self.minfails = minfails # minimum number of fails in a row
		self.dspp = [] # species dropped
		self.nseqs = 0 # counter for seqs
		for i, seqfile in enumerate(seqfiles):
			name = re.sub('\.fasta$', '', seqfile)
			seqdir = os.path.join(genedir, seqfile)
			seqs = []
			lengths = []
			with open(seqdir, "rU") as infile:
				for record in SeqIO.parse(infile, "fasta"):
					record.id = name
					lengths.append(len(record))
					seqs.append([record, 0]) # seqrecord + nfails
					self.nseqs += 1
			if len(seqs) > 0:
				self[name] = [seqs, np.min(lengths)]

	def start(self, n):
		"""Return n starting random sp sequences, update sppool"""
		self.sppool = self.keys()
		self.alignobj = []
		for i in range(n):
			self._add()
		return [e[0] for e in self.alignobj]

	def _add(self):
		"""Add new random species' sequence"""
		rand_int = random.randint(0, (len(self.sppool)-1))
		self.next_sp = self.sppool.pop(rand_int)
		next_seq = random.sample(self[self.next_sp][0], 1)[0]
		self.alignobj.append(next_seq)
		
	def back(self):
		"""Add to nfails for random species, return new random species"""
		self.alignobj[-1][1] += 1
		print [e[1] for e in self.alignobj]
		del self.alignobj[-1]
		self.sppool.append(self.next_sp)
		self._check()
		if len(self.sppool) != 0:
			self._add()
		return [e[0] for e in self.alignobj]
			
	def next(self, align):
		"""Set nfails to 0, add additional random species"""
		for i in range(len(self.alignobj)):
			self.alignobj[i][1] = 0
		self._add()
		return [e[0] for e in self.alignobj]
		
	def _check(self):
		"""Check nfails, drop sequences and species"""
		spp = self.keys()
		for sp in spp:
			to_drop = []
			for i in range(len(self[sp][0])):
				if self[sp][0][i][1] > self.minfails:
					print "Dropping [{0}] for [{1}] as nfails is [{2}]".\
						format(self[sp][0][i][0].description,sp, self[sp][0][i][1])
					to_drop.append(i)
			for i in to_drop:
				del self[sp][0][i]
			if len(self[sp][0]) < 1:
				if sp == "outgroup":
					raise OutgroupError
				else:
					print "Dropped [{0}]".format(sp)
					self.dspp.append(sp)
					self.sppool = [e for e in self.sppool if e != sp]
					del self[sp]
		if len(self.keys()) < 5: # Raise error if fewer than 5 species left in spp pool
			raise MinSpeciesError

## Functions
def alignCheck(align, mingaps, minoverlap, minlen):
	def calcOverlap(columns):
		pcolgaps = []
		for i in columns:
			ith =  float(align[:,i].count("-"))/(len(align) - 1)
			pcolgaps.append(ith)
		overlap = len(columns) - (len(columns) * np.mean(pcolgaps))
		return overlap
	align_len = align.get_alignment_length()
	if align_len < minlen:
		return False
	pintgaps = []
	for each in align:
		sequence = each.seq.tostring()
		columns = [ei for ei,e in enumerate(sequence) if e != "-"]
		totnucs = len(columns)
		totgaps = align_len - totnucs
		overlap = calcOverlap(columns)
		if overlap < minoverlap:
			print "Not enough overlap"
			return False
		extgaps = 0
		start_extgaps = re.search("^-+", sequence)
		end_extgaps = re.search("-+$", sequence)
		if not start_extgaps is None:
			extgaps += start_extgaps.end()
		if not end_extgaps is None:
			extgaps += len(sequence) - end_extgaps.start()
		intgaps = totgaps - extgaps
		overlap = intgaps + totnucs
		pintgap = float(intgaps)/overlap
		if each.id == "outgroup":
			pintgap_outgroup = pintgap
		else:
			pintgaps.append(pintgap)
		if pintgap > mingaps:
			print "Too many gaps ..."
			return False
	#try:
	#	if any([e > pintgap_outgroup for e in pintgaps]):
	#		print "Outgroup has fewer gaps than rest ..."
	#		return False # if any seqs have more gaps than outgroup, false
	#except UnboundLocalError:
	#	pass
	return True

def alignSequences(sequences, timeout=99999999, silent=False, verbose=True):
	"""Adapted pG function : aligns sequences using mafft (external program)"""
	input_file = "alignment_in" + '.fasta'
	output_file = "alignment_out" + '.fasta'
	command_line = 'mafft --auto ' + input_file + " > " + output_file
	with open(input_file, "w") as file:
		count = SeqIO.write(sequences, file, "fasta")
	pipe = TerminationPipe(command_line, timeout)
	pipe.run(silent=silent)
	os.remove(input_file)
	if not pipe.failure:
		try:
			res = AlignIO.read(output_file, 'fasta')
		except:
			raise RuntimeError("No MAFFT output.")
		os.remove(output_file)
	else:
		raise RuntimeError("Mafft alignment not complete in time allowed")
	return res

def incrAlign(seqobj, mingaps, minoverlap, seedsize, minseedsize, maxseedtrys):
	"""Incrementally build an alignment"""
	def returnBestAlignment(alignments):
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
	maxtrys = maxseedtrys*2 # twice as many attempts with minseedsize
	seedtrys = 0
	alignstore = [] # stores successful aligns generated
	# start with as many sequences at start as possible
	while True:
		print "Seed phase: [{0}] seedsize".format(seedsize)
		print seedtrys
		minlen = min([seqobj[e][1] for e in seqobj.keys()])
		sequences = seqobj.start(seedsize)
		align = alignSequences(sequences)
		if alignCheck(align, mingaps, minoverlap, minlen):
			if len(seqobj.sppool) == 0:
				return align, seedsize
			else:
				alignstore.append(align)
				sequences = seqobj.next(align)
			break
		trys += 1
		# drop seedsize by 1 if maxseedtrys
		if seedsize > minseedsize:
			seedtrys += 1
			if maxseedtrys < seedtrys:
				seedtrys = 0
				seedsize -= 1
		if maxtrys < trys:
			return None, seedsize
	trys = 0
	while True:
		minlen = min([seqobj[e][1] for e in seqobj.keys()])
		align = alignSequences(sequences)
		if alignCheck(align, mingaps, minoverlap, minlen):
			print "success"
			trys = 0
			if len(seqobj.sppool) == 0:
				return align, seedsize
			else:
				alignstore.append(align)
				alignobj = seqobj.next(align)
		elif trys < maxtrys:
			print "fail"
			sequences = seqobj.back()
			if len(seqobj.sppool) == 0:
				# here a species has been dropped and now all species are present
				return returnBestAlignment(alignstore), seedsize
			trys += 1
		else:
			print "fail"
			# when the maximum number of species is not reached...
			# ... return the best alignment in the alignment store
			return returnBestAlignment(alignstore), seedsize

if __name__ == '__main__':
	pass