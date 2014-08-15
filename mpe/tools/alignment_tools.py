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
from Bio.Align import AlignInfo
from Bio.SeqRecord import SeqRecord
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
minoverlap, runtime = 10000):
		self.minfails = minfails # minimum number of fails in a row
		self.dspp = [] # species dropped
		self.nseqs = 0 # counter for seqs
		self.mingaps = mingaps
		self.minoverlap = minoverlap
		self.runtime = runtime # limit to the number of loops a while loop can make
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

	def _add(self):
		"""Add a random sequence to sequneces_in_alignment"""
		rand_int = random.randint(0, (len(self.sppool)-1))
		self.next_sp = self.sppool.pop(rand_int)
		next_seq = random.sample(self[self.next_sp][0], 1)[0]
		self.sequences_in_alignment.append(next_seq)

	def start(self, n):
		"""Return n starting random sp sequences, update sppool"""
		def countSeqsLeft(overlap, available):
			# counts number of sequences that could be used in alignment
			if overlap:
				# note, this only calcs an esimtate, the seq may be sampled multiple times
				available = sum([len(self[e][0]) for e in self.sppool])
			else:
				# if blast overlap fails, minus 1 from current available
				available -= 1
			return available
		# pool starts with all species
		self.sppool = self.keys()
		# list of seqs in alignments is empty at first
		self.sequences_in_alignment = []
		# add a random seq to sequences_in_alignment
		self._add()
		# calc availble
		available = countSeqsLeft(True, 0)
		while True:
			# add
			self._add()
			subjseqs = [e[0] for e in self.sequences_in_alignment[:-1]]
			queryseq = self.sequences_in_alignment[-1][0]
			# blast all against last
			blast_bool = blast(subj = subjseqs, query = queryseq,\
				minoverlap = self.minoverlap, mingaps = self.mingaps)
			# did majority overlap?
			overlap = (float(sum(blast_bool))/len(blast_bool)) > 0.5
			# if the majority didn't overlap
			if not overlap:
				del self.sequences_in_alignment[-1]
				self.sppool.append(self.next_sp)
			available = countSeqsLeft(overlap, available)
			# break from loop if all seqs have been tried or n is hit
			if available == 0 or len(self.sequences_in_alignment) == n:
				break
		return [e[0] for e in self.sequences_in_alignment]

	def _blastAdd(self, alignment):
		"""Add new random species' sequence ensuring overlap with \
BLAST"""
		rand_ints = range(len(self.sppool))
		random.shuffle(rand_ints)
		for i in rand_ints:
			sp = self.sppool[i]
			next_seqs = [e[0] for e in self[sp][0]]
			next_seq_is = self._blastAlignment(alignment, next_seqs)
			if next_seq_is:
				next_seq_i = random.sample(next_seq_is, 1)[0]
				break
		else:
			return False
		self.next_sp = self.sppool.pop(i)
		#print self.next_sp
		self.sequences_in_alignment.append(self[self.next_sp][0][\
next_seq_i])
		return True
		
	def back(self):
		"""Add to nfails for random species, return new random \
species"""
		self.sequences_in_alignment[-1][1] += 1
		del self.sequences_in_alignment[-1]
		self.sppool.append(self.next_sp)
		self._check()
		if len(self.sppool) != 0:
			self._add()
			return self.sequences_in_alignment[-1][0]
		else:
			return False
			
	def next(self, alignment):
		"""Set nfails to 0, add additional random species"""
		for i in range(len(self.sequences_in_alignment)):
			self.sequences_in_alignment[i][1] = 0
		if self._blastAdd(alignment):
			return self.sequences_in_alignment[-1][0]
		else:
			return False

	def _blastAlignment(self, alignment, sequences):
		"""Match a sequence to an alignment with stand-alone blast \
to determine if sequences overlap. Return indexes of overlapping \
sequences."""
		# convert alignment into a consensus
		summary_align = AlignInfo.SummaryInfo(alignment)
		consensus = SeqRecord(summary_align.gap_consensus(ambiguous \
			= 'N', threshold = 0.5), id = "con", name = \
		"Alignment consensus", description = \
		"ambiguous = N, threshold = 0.5")
		results = blast (subj = sequences, query = consensus, \
			minoverlap = self.minoverlap, mingaps = self.mingaps)
		return [i for i,e in enumerate(results) if e]
		
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
		maxtrys, maxseedtrys, gene_type):
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

	def _check(self, alignment):
		return checkAlignment(alignment, self.mingaps, \
			self.minoverlap, self.minlen)

	def _return(self, store):
		"""Return best alignment from a list of alignments based on:\
presence of outgroup, number of species and length of alignment"""
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
			if len(sequences) >= self.minseedsize: # make sure there are enough seqs for alignment
				command = version(sequences, self.type)
				try:
					alignment = align(command, sequences)
				except MafftError:
					pass
			else:
				success = False
			success = self._check(alignment)
			trys += self._calcSeedsize(success) # add to trys if seedsize is too small
			if self.maxtrys < trys:
				return None
			if success: # if alignment is successful, return alignment or move to next phase
				if len(self.seqstore.sppool) == 0:
					return alignment
				else:
					sequence = self.seqstore.next(alignment)
					if sequence:
						store.append(alignment)
						break
		trys = 0
		logging.info("........ add phase : [{0}] species".format(len(\
			alignment)))
		while True:
			self.minlen = min([self.seqstore[e][1] for e in self.\
				seqstore.keys()])
			#print "Number of species: {0}".format(len(alignment))
			#TODO: I assume I can't use qinsi or xinsi with --add
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
				sequence = self.seqstore.back()
				if not sequence:
					# here a species has been dropped and now all species are present
					return self._return(store)
				alignment = store[-1]
				trys += 1
			else:
				logging.info("............ maxtrys hit")
				# when the maximum number of species is not reached...
				# ... return the best alignment in the alignment store
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

def blast(subj, query, minoverlap, mingaps):
	"""Return True or False for each sequence in subj that overlaps \
with sequences in query given set parameters using NCBI's BLAST"""
	SeqIO.write(query, ".query.fasta", "fasta")
	SeqIO.write(subj, ".subj.fasta", "fasta")
	output = NcbiblastnCommandline(query = ".query.fasta",\
	subject = ".subj.fasta", outfmt = 5)()[0]
	try:
		output = NcbiblastnCommandline(query = ".query.fasta",\
			subject = ".subj.fasta", outfmt = 5)()[0]
	except ApplicationError:
		logging.warn("---- BLAST Error ----")
		return False
	finally:
		os.remove(".query.fasta")
		os.remove(".subj.fasta")
	bools = []
	bresults = NCBIXML.parse(StringIO(output)) # BLAST records for each sequence
	for record in bresults:
		# if there is a hit....
		if record.alignments:
			# work out how long it is ignoring Ns
			ns = record.alignments[0].hsps[0].query.count("N")
			length = record.alignments[0].hsps[0].align_length - ns
			if length > minoverlap:
				# if it's above the minoverlap
				identities = record.alignments[0].hsps[0].identities
				pgaps = 1 - float(length)/identities
				if pgaps < mingaps:
					# and it has more identities than mingaps give it a True
					bools.append(True)
					continue
		bools.append(False)
	# return True if all in query overlap with subj
	return bools

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
			#print "not enough overlap"
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
			#print "too many gaps"
			return False
	#try:
	#	if any([e > pintgap_outgroup for e in pintgaps]):
	#		print "Outgroup has fewer gaps than rest ..."
	#		return False # if any seqs have more gaps than outgroup, false
	#except UnboundLocalError:
	#	pass
	return True