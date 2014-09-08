#! /bin/usr/env python
## D.J. Bennett
## 24/03/2014
"""
mpe alignment tools
"""

## Packages
import os,re,random,logging
import numpy as np
from Bio import SeqIO
from Bio import AlignIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Application import ApplicationError
from Bio.Blast import NCBIXML
from StringIO import StringIO
from system_tools import TerminationPipe
from system_tools import OutgroupError
from system_tools import TooFewSpeciesError
from system_tools import MafftError

## Objects
class SeqStore(dict):
	"""Store species' gene sequences with functions for pulling \
sequences for alignments and adding penalties for sequences that did \
not align"""
	def __init__(self, genedir, seqfiles, minfails, mingaps, \
minoverlap):
		self.minfails = minfails # minimum number of fails in a row
		self.dspp = [] # species dropped
		self.nseqs = 0 # counter for seqs
		self.blast_prop = 0.99 # the proportion of sequences in the alignment that a sequence must overlap with
		self.mingaps = mingaps
		self.minoverlap = minoverlap
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

	def _add(self, sequences = None):
		"""Return a random sequence for alignment"""
		if len(self.sequences_in_alignment) > 5:
			# if there are more than 5 sequences, use blast alignment
			rand_ints = range(len(self.sppool))
			random.shuffle(rand_ints)
			for i in rand_ints:
				sp = self.sppool[i]
				next_seqs = [e[0] for e in self[sp][0]]
				# blast next_seqs against sequences in alignment
				next_seq_is,next_seqs = self._alignment_blast(\
					next_seqs,sequences)
				# if any overlap, break
				if next_seq_is:
					break
			else:
				return None
			self.next_sp = self.sppool.pop(i)
			# choose next_seq at random
			if (len(next_seq_is) > 1):
				next_seq_i = random.sample(next_seq_is, 1)[0]
			else:
				next_seq_i = next_seq_is[0]
			# get blasted next_seq
			next_seq = next_seqs[next_seq_is.index(next_seq_i)]
			# record sequence + nfails in sequence_in_alignment
			self.sequences_in_alignment.append(self[self.next_sp][0]\
				[next_seq_i])
		elif len(self.sequences_in_alignment) > 0:
			# if there are fewer than minseedsize sequences, use
			#  normal blast -- not alignment blast. This avoids single
			#  sequence alignments which are likely to fail.
			rand_ints = range(len(self.sppool))
			random.shuffle(rand_ints)
			for i in rand_ints:
				sp = self.sppool[i]
				next_seqs = [e[0] for e in self[sp][0]]
				# blast sequences in alignment against next_seqs
				next_seq_bools,_ = blast(query = sequences, subj = \
					next_seqs, minoverlap = self.minoverlap, mingaps\
					= self.mingaps)
				# if all overlap, break -- all next seqs overlap
				if all(next_seq_bools):
					break
			else:
				return None
			self.next_sp = self.sppool.pop(i)
			# choose next_seq at random
			next_seq = random.sample(next_seqs, 1)[0]
			next_seq_i = next_seqs.index(next_seq)
			# record sequence + nfails in sequence_in_alignment
			self.sequences_in_alignment.append(self[self.next_sp][0]\
				[next_seq_i])
		else:
			# without any sequences in alignment
			#  simply return a random seq
			rand_int = random.randint(0, (len(self.sppool)-1))
			self.next_sp = self.sppool.pop(rand_int)
			result = random.sample(self[self.next_sp][0], 1)[0]
			# record sequence + nfails in sequence_in_alignment
			self.sequences_in_alignment.append(result)
			next_seq = result[0]
		return next_seq

	def _alignment_blast(self, query, sequences_in_alignment):
		"""Return indexes and overlapping sequences for each sequence 
in query that overlaps with more than prop sequences in subj given 
set parameters using NCBI's BLAST"""
		# return query is and overlapping sequences if blast successful
		query_is = []
		sequences = []
		# loop through each sequence in query
		for i in range(len(query)):
			bools,positions = blast(query[i], sequences_in_alignment,\
				self.minoverlap, self.mingaps)
			# if more than prop overlap
			overlap = (float(sum(bools))/len(sequences_in_alignment))\
			> self.blast_prop
			if overlap:
				query_is.append(i)
				# determine overlap as the max and min for all
				#  successful blast results
				sequences.append(query[i][min(positions):\
					max(positions)])
		return query_is,sequences

	def start(self, n):
		"""Return n starting random sp sequences, update sppool"""
		# pool starts with all species
		self.sppool = self.keys()
		# sequences + nfails used in alignment
		self.sequences_in_alignment = []
		# actual sequences for alignment
		sequences = []
		# add a random seq to sequences_in_alignment and sequences
		sequences.append(self._add())
		while True:
			sequence = self._add(sequences)
			if sequence:
				sequences.append(sequence)
			# break from loop if n is hit or sequence is None
			if not sequence or len(self.sequences_in_alignment) == n:
				break
		return sequences
		
	def back(self, alignment):
		"""Add to nfails for random sequence, return new random \
species"""
		self.sequences_in_alignment[-1][1] += 1
		del self.sequences_in_alignment[-1]
		self.sppool.append(self.next_sp)
		self._check()
		if len(self.sppool) != 0:
			sequences = [e for e in alignment]
			return self._add(sequences)
		else:
			return None
			
	def next(self, alignment):
		"""Set nfails to 0, add additional random species"""
		for i in range(len(self.sequences_in_alignment)):
			self.sequences_in_alignment[i][1] = 0
		sequences = [e for e in alignment]
		return self._add(sequences)
		
	def _check(self):
		"""Check nfails, drop sequences and species"""
		spp = self.keys()
		for sp in spp:
			to_drop = [ei for ei,e in enumerate(self[sp][0])\
				if e[1] > self.minfails]
			for i in to_drop:
				logging.info("Dropping [{0}] for [{1}] as nfails is \
[{2}]".format(self[sp][0][i][0].description,sp,self[sp][0][i][1]))
			self[sp][0] = [e for ei,e in enumerate(self[sp][0])\
				if ei not in to_drop]
			if len(self[sp][0]) < 1:
				if sp == "outgroup":
					raise OutgroupError
				else:
					logging.info("Dropped [{0}]".format(sp))
					self.dspp.append(sp)
					self.sppool = [e for e in self.sppool if e != sp]
					del self[sp]
		if len(self.keys()) < 5:
			raise TooFewSpeciesError

class Aligner(object):
	"""Build alignments from seqstore"""
	def __init__(self, seqstore, mingaps, minoverlap, minseedsize,\
		maxtrys, maxseedtrys, gene_type, outgroup = True):
		self.seqstore = seqstore
		self.mingaps = mingaps
		self.minoverlap = minoverlap
		self.minseedsize = minseedsize
		self.maxtrys = maxtrys
		self.buffer = maxseedtrys
		self.buffer_counter = 0 # seedsize buffer counter
		self.seedsize = len(seqstore)
		self.timeout = 99999999
		self.silent = False
		self.type = gene_type
		self.outgroup = outgroup # bool for whether outgroup is necessary

	def _check(self, alignment):
		return checkAlignment(alignment, self.mingaps, \
			self.minoverlap, self.minlen)

	def _return(self, store):
		"""Return best alignment from a list of alignments based on:\
presence of outgroup, number of species and length of alignment"""
		if self.outgroup:
			# keep alignments with outgroups
			# keep alignments with more than 5 species
			store = [a for a in store if "outgroup" in\
						  [e.id for e in a._records]]
			if len(store) == 0:
				logging.info("........ no outgroup")
				return None
		store = [e for e in store if len(e) >= 5]
		if len(store) == 0:
			logging.info("........ too few species")
			return None
		# keep only alignments with lots of records
		nrecords = [len(e._records) for e in store]
		store = [store[i] for i,e in enumerate(nrecords)\
					  if e == max(nrecords)]
		# return longest alignment
		lens = [len(e) for e in store]
		max_i = lens.index(max(lens))
		return store[max_i]

	def _calcSeedsize(self, success = True):
		"""Calculate seedsize based on buffer and success of current \
seedsize. Return 1 if trys must increase, else 0."""
		# increase seedsize if successful buffer times in a row
		# decrease seedsize if unsuccessful buffer times in a row
		if success:
			self.buffer_counter += 1
		else:
			self.buffer_counter -= 1
		if self.buffer_counter >= self.buffer:
			if self.seedsize < len(self.seqstore):
				self.seedsize += 1
				self.buffer_counter = 0# reset every time seedsize changes
		if self.buffer_counter <= (self.buffer * -1):
			if self.seedsize > self.minseedsize:
				self.seedsize -= 1
				self.buffer_counter = 0
			else:
				# here seedsize must be reduced, but has hit
				#  midseedsize, add 1 to trys
				return 1
		return 0

	def run(self):
		"""Incrementally build an alignment by adding sequences to a \
seed alignment"""
		trys = 0
		store = []
		logging.info("........ seed phase: [{0}] seed size".format(\
			self.seedsize))
		while True:
			self.minlen = min([self.seqstore[e][1] for e in self.\
				seqstore.keys()])
			sequences = self.seqstore.start(self.seedsize)
			# make sure there are enough seqs for alignment
			if len(sequences) >= self.minseedsize:
				command = version(sequences, self.type)
				try:
					alignment = align(command, sequences)
				except MafftError:
					continue
				else:
					success = self._check(alignment)
			else:
				success = False
			# add to trys if unsuccessful or sequences are fewer than
			#  seedsize
			trys += self._calcSeedsize(success and len(sequences) == \
				self.seedsize)
			if self.maxtrys < trys:
				return None
			# if alignment is successful, return alignment or move to
			#  next phase
			if success:
				if len(self.seqstore.sppool) == 0:
					return alignment
				else:
					# if a sequence can be added move to add phase
					sequence = self.seqstore.next(alignment)
					if sequence:
						store.append(alignment)
						break
					# if an alignment can be returned, return it
					returnable_alignment = self._return([alignment])
					if returnable_alignment:
						return returnable_alignment
		trys = 0
		logging.info("........ add phase : [{0}] species".format(len(\
			alignment)))
		while True:
			self.minlen = min([self.seqstore[e][1] for e in self.\
				seqstore.keys()])
			try:
				alignment = add(alignment, sequence)
			except MafftError:
				pass
			if self._check(alignment):
				trys = 0
				if len(self.seqstore.sppool) == 0:
					return alignment
				else:
					store.append(alignment)
					sequence = self.seqstore.next(alignment)
					if not sequence:
						return self._return(store)
			elif trys < self.maxtrys:
				alignment = store[-1]
				sequence = self.seqstore.back(alignment)
				if not sequence:
					# here a species has been dropped and now all
					#  species are present
					return self._return(store)
				trys += 1
			else:
				logging.info("............ maxtrys hit")
				# when the maximum number of species is not reached
				#  return the best alignment in the alignment store
				return self._return(store)

## Functions
def version(sequences, gene_type):
	"""Return command for best MAFFT version given sequences"""
	# determine auto, qinsi or xinsi based on:
	# http://mafft.cbrc.jp/alignment/software/source66.html
	# always using default algorithms
	if gene_type != 'deep' or len(sequences) > 250:
		return 'mafft --auto'
	seqlens = [len(s) for s in sequences]
	if max(seqlens) > 1500:
		return 'mafft --auto'
	if len(sequences) > 60:
		return 'mafft-qinsi'
	return 'mafft-xinsi'

def align(command, sequences):
	"""Align sequences using mafft (external program)"""
	input_file = ".sequences_in.fasta"
	output_file = ".alignment_out.fasta"
	command_line = '{0} {1} > {2}'.format(command, input_file, \
		output_file)
	with open(input_file, "w") as file:
		SeqIO.write(sequences, file, "fasta")
	pipe = TerminationPipe(command_line)
	pipe.run()
	os.remove(input_file)
	if not pipe.failure:
		try:
			res = AlignIO.read(output_file, 'fasta')
		except:
			logging.info(pipe.output)
			raise MafftError()
		else:
			os.remove(output_file)
	else:
		raise RuntimeError("MAFFT alignment not complete in time \
allowed")
	return res

def add(alignment, sequence):
	"""Align sequence(s) to an alignment using mafft (external program)"""
	alignment_file = ".alignment_in.fasta"
	sequence_file = ".sequence_in.fasta"
	output_file = "alignment_out.fasta" + '.fasta'
	command_line = 'mafft --auto --add {0} {1} > {2}'.format(\
		sequence_file, alignment_file, output_file)
	with open(sequence_file, "w") as file:
		SeqIO.write(sequence, file, "fasta")
	with open(alignment_file, "w") as file:
		AlignIO.write(alignment, file, "fasta")
	pipe = TerminationPipe(command_line)
	pipe.run()
	os.remove(alignment_file)
	os.remove(sequence_file)
	if not pipe.failure:
		try:
			res = AlignIO.read(output_file, 'fasta')
		except:
			logging.info(pipe.output)
			raise MafftError()
		else:
			os.remove(output_file)
	else:
		raise RuntimeError("MAFFT alignment not complete in time allowed")
	return res

def blast(query, subj, minoverlap, mingaps):
	"""Return bool and positions of query sequences that overlapped
with subject given parameters."""
	SeqIO.write(query, ".query.fasta", "fasta")
	SeqIO.write(subj, ".subj.fasta", "fasta")
	try:
		cline = NcbiblastnCommandline(query = ".query.fasta",\
			subject = ".subj.fasta", outfmt = 5)
		output = cline()[0]
	except ApplicationError:
		logging.warn("---- BLAST Error ----")
		return [],[]
	finally:
		os.remove(".query.fasta")
		os.remove(".subj.fasta")
	# list of T or F for success of alignment between queries and
	#  subject
	bools = []
	# record start and end position to avoid composite sequence
	#  problems
	positions = []
	# BLAST records for each query sequence matched against subj
	bresults = NCBIXML.parse(StringIO(output))
	for record in bresults:
		# if there is a hit....
		if record.alignments:
			res = record.alignments[0].hsps[0]
			# work out how long it is ignoring Ns
			ns = res.query.count("N")
			length = res.align_length - ns
			if length > minoverlap:
				# if it's above the minoverlap
				identities = res.identities
				pgaps = 1 - float(length)/identities
				if pgaps < mingaps:
					# and it has more identities than mingaps give
					#  it a True
					bools.append(True)
					positions.append(res.sbjct_start)
					positions.append(res.sbjct_end)
					continue
		bools.append(False)
	return bools,positions

def checkAlignment(alignment, mingaps, minoverlap, minlen):
	"""Determine if an alignment is good or not based on given \
parameters. Return bool"""
	def calcOverlap(columns):
		pcolgaps = []
		for i in columns:
			ith = float(alignment[:,i].count("-"))/(len(alignment)-1)
			pcolgaps.append(ith)
		overlap = len(columns) - (len(columns) * np.mean(pcolgaps))
		return overlap
	alen = alignment.get_alignment_length()
	if alen < minlen:
		return False
	pintgaps = []
	for each in alignment:
		sequence = str(each.seq)
		columns = [ei for ei,e in enumerate(sequence) if e != "-"]
		totnucs = len(columns)
		totgaps = alen - totnucs
		overlap = calcOverlap(columns)
		if overlap < minoverlap:
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
		#if each.id == "outgroup":
		#	pintgap_outgroup = pintgap
		#else:
		pintgaps.append(pintgap)
		if pintgap > mingaps:
			return False
	#try:
	#	if any([e > pintgap_outgroup for e in pintgaps]):
	#		print "Outgroup has fewer gaps than rest ..."
	#		return False # if any seqs have more gaps than outgroup, false
	#except UnboundLocalError:
	#	pass
	return True