#!/usr/bin/python
## MPE Alignment tools
## D.J. Bennett
## 24/03/2014

## Packages
import sys, os, re, random
import numpy as np
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from StringIO import StringIO
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
		self.sequences_in_alignment = []
		for i in range(n):
			rand_int = random.randint(0, (len(self.sppool)-1))
			self.next_sp = self.sppool.pop(rand_int)
			next_seq = random.sample(self[self.next_sp][0], 1)[0]
			self.sequences_in_alignment.append(next_seq)
		return [e[0] for e in self.sequences_in_alignment]

	def _add(self, alignment):
		"""Add new random species' sequence"""
		rand_ints = range(len(self.sppool))
		random.shuffle(rand_ints)
		for i in rand_ints:
			sp = self.sppool[i]
			next_seqs = [e[0] for e in self[sp][0]]
			next_seq_is = findOverlappingSeqs(alignment, next_seqs)
			if next_seq_is:
				next_seq_i = random.sample(next_seq_is, 1)[0]
				break
		else:
			return False
		self.next_sp = self.sppool.pop(i)
		print self.next_sp
		self.sequences_in_alignment.append(self[self.next_sp][0][next_seq_i])
		return True
		
	def back(self):
		"""Add to nfails for random species, return new random species"""
		self.sequences_in_alignment[-1][1] += 1
		del self.sequences_in_alignment[-1]
		self.sppool.append(self.next_sp)
		self._check()
		if len(self.sppool) != 0:
			rand_int = random.randint(0, (len(self.sppool)-1))
			self.next_sp = self.sppool.pop(rand_int)
			next_seq = random.sample(self[self.next_sp][0], 1)[0]
			self.sequences_in_alignment.append(next_seq)
		return self.sequences_in_alignment[-1][0]
			
	def next(self, alignment):
		"""Set nfails to 0, add additional random species"""
		for i in range(len(self.sequences_in_alignment)):
			self.sequences_in_alignment[i][1] = 0
		if self._add(alignment):
			return self.sequences_in_alignment[-1][0]
		else:
			return False
		
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
		if len(self.keys()) < 5:
			raise MinSpeciesError

## Functions
def alignmentCheck(align, mingaps, minoverlap, minlen):
	"""Determine if an alignment is good or not based on given parameters"""
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
		if each.id == "outgroup":
			pintgap_outgroup = pintgap
		else:
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

def findOverlappingSeqs (alignment, sequences):
	"""Match a sequence to an alignment with stand-alone blast to determine if sequences overlap"""
	summary_align = AlignInfo.SummaryInfo(alignment)
	consensus = SeqRecord(summary_align.gap_consensus(ambiguous = 'N', threshold = 0.5),\
		id = "con", name = "Alignment consensus", description = "ambiguous = N, threshold = 0.5")
	#print consensus.format("fasta")
	#AlignIO.write(alignment, "subj.fasta", "fasta")
	SeqIO.write(consensus, "query.fasta", "fasta")
	SeqIO.write(sequences, "subj.fasta", "fasta")
	sresults = []
	try:
		output = NcbiblastnCommandline(query = "query.fasta", subject = "subj.fasta", outfmt = 5)()[0]
	except Bio.Application.ApplicationError:
		print "BLAST error"
	else:
		bresults = NCBIXML.parse(StringIO(output)) # BLAST records for each sequence
		i = 0
		for record in bresults:
			if record.alignments:
				ns = record.alignments[0].hsps[0].query.count("N")
				length = record.alignments[0].hsps[0].align_length - ns
				if length > 200:
					identities = record.alignments[0].hsps[0].identities
					pgaps = 1 - float(length)/identities
					if pgaps < 0.1:
						sresults.append(i)
			i += 1
	os.remove("query.fasta")
	os.remove("subj.fasta")
	return sresults

def sequenceCheck (alignment, sequence, minoverlap, mingaps):
	"""Match a sequence to a alignment with stand-alone blast"""
	summary_align = AlignInfo.SummaryInfo(alignment)
	consensus = SeqRecord(summary_align.gap_consensus(ambiguous = 'N', threshold = 0.5),\
		id = "con", name = "Alignment consensus", description = "ambiguous = N, threshold = 0.5")
	#print consensus.format("fasta")
	#AlignIO.write(alignment, "subj.fasta", "fasta")
	SeqIO.write(consensus, "subj.fasta", "fasta")
	SeqIO.write(sequence, "query.fasta", "fasta")
	success = False
	try:
		output = NcbiblastnCommandline(query = "query.fasta", subject = "subj.fasta", outfmt = 5)()[0]
	except Bio.Application.ApplicationError:
		print "BLAST error"
	else:
		results = NCBIXML.parse(StringIO(output))
		for record in results:
			if record.alignments:
				ns = record.alignments[0].hsps[0].sbjct.count("N")
				length = record.alignments[0].hsps[0].align_length - ns
				if length > minoverlap:
					identities = record.alignments[0].hsps[0].identities
					pgaps = 1 - float(length)/identities
					if pgaps < mingaps:
						success = True
				break
	os.remove("query.fasta")
	os.remove("subj.fasta")
	if success:
		return True
	else:
		return False

def alignSequences(sequences, timeout=99999999, silent=False, verbose=True):
	"""Align sequences using mafft (external program)"""
	input_file = "sequences_in.fasta"
	output_file = "alignment_out.fasta"
	command_line = 'mafft --auto {0} > {1}'.format(input_file,output_file)
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

def addToAlignment(alignment, sequence, timeout=99999999, silent=False, verbose=True):
	"""Align sequence(s) to an alignment using mafft (external program)"""
	alignment_file = "alignment_in.fasta"
	sequence_file = "sequence_in.fasta"
	output_file = "alignment_out.fasta" + '.fasta'
	command_line = 'mafft --auto --add {0} {1} > {2}'.format(sequence_file,alignment_file,output_file)
	with open(sequence_file, "w") as file:
		count = SeqIO.write(sequence, file, "fasta")
	with open(alignment_file, "w") as file:
		count = AlignIO.write(alignment, file, "fasta")
	pipe = TerminationPipe(command_line, timeout)
	pipe.run(silent=silent)
	os.remove(alignment_file)
	os.remove(sequence_file)
	if not pipe.failure:
		try:
			res = AlignIO.read(output_file, 'fasta')
		except:
			raise RuntimeError("No MAFFT output.")
		os.remove(output_file)
	else:
		raise RuntimeError("Mafft alignment not complete in time allowed")
	return res

def returnBestAlignment(alignments):
	"""Return best alignment from a list of alignments based on: presence of outgroup,
	number of species and length of alignment"""
	# keep alignments with outgroups
	# keep alignments with more than 5 species
	alignments = [al for al in alignments if "outgroup" in\
				  [e.id for e in al._records]]
	if len(alignments) == 0:
		print "No outgroup!"
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

def incrAlign(seqobj, mingaps, minoverlap, seedsize, minseedsize, maxtrys, maxseedtrys):
	"""Incrementally build an alignment by adding sequences to a seed alignment"""
	trys = seedtrys = 0
	print maxtrys
	alignment_store = []
	print " ........ seed phase: [{0}] seed size".format(seedsize)
	while True:
		minlen = min([seqobj[e][1] for e in seqobj.keys()])
		sequences = seqobj.start(seedsize)
		alignment = alignSequences(sequences)
		if alignmentCheck(alignment, mingaps, minoverlap, minlen):
			if len(seqobj.sppool) == 0:
				return alignment, seedsize
			else:
				alignment_store.append(alignment)
				sequence = seqobj.next(alignment)
			break
		if seedsize > minseedsize:
			seedtrys += 1
			if maxseedtrys < seedtrys:
				seedtrys = 0
				seedsize -= 1
		else:
			trys += 1
			if maxtrys < trys:
				return None, seedsize
	trys = 0
	print  " ........ add phase : [{0}] species".format(len(alignment))
	while True:
		minlen = min([seqobj[e][1] for e in seqobj.keys()])
		#print "Number of species: {0}".format(len(alignment))
		alignment = addToAlignment(alignment, sequence)
		if alignmentCheck(alignment, mingaps, minoverlap, minlen):
			trys = 0
			if len(seqobj.sppool) == 0:
				return alignment, seedsize
			else:
				alignment_store.append(alignment)
				sequence = seqobj.next(alignment)
				if not sequence:
					return returnBestAlignment(alignment_store), seedsize
		elif trys < maxtrys:
			sequence = seqobj.back()
			if len(seqobj.sppool) == 0:
				# here a species has been dropped and now all species are present
				return returnBestAlignment(alignment_store), seedsize
			alignment = alignment_store[-1]
			trys += 1
		else:
			print "maxtrys hit!"
			# when the maximum number of species is not reached...
			# ... return the best alignment in the alignment store
			return returnBestAlignment(alignment_store), seedsize

if __name__ == '__main__':
	pass