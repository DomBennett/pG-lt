#!/usr/bin/env python
## D. J. Bennett
## 16/06/2013
##TODO: Jun's idea of counting Ns and ?s in downloaded seqs and taking them away from seq length.

import random
import sys
sys.path.append('/home/dominic/programs/phyloGenerator-master')
import phyloGenerator_adapted as pG

class SequenceParser(object):
	"""
Parse sequencer class for selecting sequences downloaded with findGenes()
in seq_obj format.
	
Doctest:
>>> mock_seq_obj = ([[['1111111111', '111'], ['2222222222', '2222222222'],\
['3333333333']], [['1111111111', '111'], ['2222222222', '2222222222'],\
['3333333333']], [['1111111111', '111'], ['2222222222', '2222222222'],\
['3333333333']], [['1111111111', '111'], ['22', '22'], ['3333333333']],\
[['1111111111', '111'], ['22', '22'], ['3333333333']], [[()], [()], [()]],\
[['1111111111', '111'], ['2222222222', '2222222222'], ['3333333333']]],\
['gene1', 'gene2', 'gene3'])
>>> gene_lens = [10, 10, 10]
>>> seq_parser = SequenceParser(mock_seq_obj, gene_lens)
>>> seq_parser.parse()
[['1111111111', '2222222222', '3333333333'], ['1111111111', '2222222222',\
'3333333333'], ['1111111111', '2222222222', '3333333333'], ['1111111111',\
'2222222222', '3333333333']]"""
	def __init__(self, seq_obj, gene_lens, gene_overlap = 0.1,
	gene_coverage = 0.1, min_ngene = 1,	min_nspp = 2):
		self.seqs = seq_obj[0]
		self.genes = seq_obj[1]
		self.ids = range(len(seq_obj[0]))
		self.min_bps = [each - (each * gene_overlap) for each in gene_lens]
		self.max_bps = [each + (each * gene_overlap) for each in gene_lens]
		self.gene_coverage = gene_coverage
		self.min_ngene = min_ngene
		self.min_nspp = min_nspp
	
	def parse(self, thorough = False):
		"""
Parse seq_obj to produce a list of sequences for alignment.
Returns: list of sequences or None (not enough sequence data)"""
		# write out taxids used
		if thorough:
			bools_spp, bool_len = self._seqLenBool(greater = True)
			if bool_len:
				self.seqs = self._findGeneInSeq(bools_spp)
		bools_spp, bool_len = self._seqLenBool()
		self.seqs = self._dropSeqs(bools_spp) # note list of list to list
		enough = self._checkEnough()
		if not enough:
			return None
		self.seqs, self.genes = self._dropGenes()
		self.seqs, self.ids = self._dropSpp(len(self.genes))
		
		if self._checkFinished():
			return self.seqs
		else:
			print "Parse failed."
			return None
		
		#N_max = len(self.genes)
		#N_list = range(1, N_max + 1)
		#N_list.reverse()
		#for N in N_list:
		#	seqs, ids = self._dropSpp(N)
		#	print seqs
		#	if len(seqs) < 5:
		#		continue
		#	else:
		#		self.seqs, self.ids = (seqs, ids)
		#	if N >= 1:
		#		self.seqs, self.genes = self._dropGenes()
		#		self.seqs, self.ids = self._dropSpp(N)
		#	if self._checkFinished():
		#		return self.seqs
		#	else:
		#		continue
		#return None

	def _seqLenBool(self, greater = False):
		"""
Determine length of sequences list of list of sequences.
Return 1 for each sequence with length greater than min_bps and less than
max_bps. Or if greater is True, return 1 for all sequences with lengths greater
than max_bp.
Returns: list of 0s and 1s"""
		output = []
		bool_len = False
		seqs_spp = self.seqs[:] # seqs for each spp: [[g1, g2], [g1, g2] ...]
		# iterate through spp
		for i in range(len(seqs_spp)):
			seqs_gen = seqs_spp[i] # seqs for each gene: [[seq1, seq2], ...]
			# iterate through genes
			gen_bool = []
			for j in range(len(seqs_gen)):
				# iterate through seq returns
				seqs_lens = [0 if seq == () else len(seq) for seq in seqs_gen[j]]
				# or change 'if isinstance(seq, pG.SeqRecord)'
				if greater:
					seqs_bool = [1 if l > self.max_bps[j] else 0 for l in seqs_lens]
				else:
					seqs_bool = [1 if self.max_bps[0] > l > self.min_bps[j] \
					else 0 for l in seqs_lens]
				if sum(seqs_bool) > 0:
					bool_len = True
				gen_bool.append(seqs_bool)
			output.append(gen_bool)
		return output, bool_len
		
	def _findGeneInSeq(self, bools_spp):
		"""
For large sequences find genes in sequence
Returns: list of list of seqs"""
		seqs_spp = self.seqs[:]
		for i in range(len(seqs_spp)):
			for j in range(len(seqs_spp[i])):
				for k in range(len(seqs_spp[i][j])):
					if bools_spp[i][j][k] == 1:
						seqs_spp[i][j][k] = pG.findGeneInSeq(seqs_spp[i][j][k],\
						self.genes[j])
		return seqs_spp
	
	def _dropSeqs(self, bools_spp):
		"""
Drop sequences from list of list of sequences based on bools_spp.
Returns: list of seqs"""
		output = []
		seqs_spp = self.seqs[:]
		for i in range(len(seqs_spp)):
			seqs_gen = seqs_spp[i]
			bools_gen = bools_spp[i]
			gen_seq = []
			for j in range(len(seqs_gen)):
				enum_bool = enumerate(bools_gen[j])
				keep_seqs = [seqs_gen[j][k] for k, bl in enum_bool if bl == 1]
				if len(keep_seqs) > 1:
					keep_seqs = random.sample(keep_seqs, 1)[0]
				elif len(keep_seqs) == 0:
					keep_seqs = 0
				else:
					keep_seqs = keep_seqs[0]
				gen_seq.append(keep_seqs)
			output.append(gen_seq)
		return output
	
	def _checkEnough(self):
		"""
Check if there is enough sequence data based on min_nspp and min_ngene.
Return: bool"""
		min_ngene = self.min_ngene
		ngenes = len(self.genes)
		for i in range(ngenes):
			presents = [1 for s in self.seqs if s[i] != 0]
			if sum(presents) >= self.min_nspp:
				min_ngene -= 1
				if min_ngene == 0:
					return True
		return False

	def _dropSpp(self, N):
		"""
Drop species in seqs list that have fewer than N gene sequences.
Returns: tuple of list of sequences and ids kept"""
		seqs = self.seqs[:]
		ids = self.ids[:]
		keep = []
		for sp in seqs:
			keep_bool = [1 for s in sp if s != 0]
			if len(keep_bool) < N:
				keep.append(0)
			else:
				keep.append(1)
		seqs = [seqs[i] for i, kb in enumerate(keep) if kb == 1]
		ids = [ids[i] for i, kb in enumerate(keep) if kb == 1]
		return seqs, ids
	
	def _genCoverage(self, seqs, ngenes, nspp):
		"""
Calculate coverages of genes in seqs list.
Returns: list of coverages"""
		gen_coverages = []
		for i in range(ngenes):
			present = sum([1 for s in seqs if s[i] != 0])
			gen_coverages.append(float(present)/float(nspp))
		return gen_coverages
	
	def _dropGenes(self):
		"""
Drop genes that are not represented within gene_coverage (default is 10%) the
coverage of the best represented gene.
Returns: list of sequences"""
		output = []
		seqs = self.seqs[:]
		ngenes = len(self.genes)
		nspp = len(self.seqs)
		gcovers = self._genCoverage(seqs, ngenes, nspp)
		cutoff = max(gcovers) - max(gcovers)*self.gene_coverage
		gene_bool = [1 if gc > cutoff else 0 for gc in gcovers]
		for i in range(nspp):
			newseqs = [seqs[i][j] for j, bl in enumerate(gene_bool) if bl == 1]
			output.append(newseqs)
		genes = [self.genes[j] for j, bl in enumerate(gene_bool) if bl == 1]
		return (output, genes)

	def _checkFinished(self):
		"""
Check all genes and species are fully represented.
Returns: bool"""
		seqs = self.seqs[:]
		for i in range(len(seqs)):
			for j in range(len(seqs[i])):
				if seqs[i][j] == 0:
					return False
		return True
		
if __name__ == "__main__":
    import doctest
    doctest.testmod()
