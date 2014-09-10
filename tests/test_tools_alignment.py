#! /bin/usr/env python
## D.J. Bennett
## 26/05/2014
"""
Tests for alignment tools.
"""

import unittest,os,re,copy,pickle,random
import mpe.tools.alignment_tools as atools
from Bio import SeqIO
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

## Dirs
working_dir = os.path.dirname(__file__)

## Test data
test_seqs = []
with open(os.path.join(working_dir,'data',\
	'test_sequences.faa'),"rU") as infile:
	for record in SeqIO.parse(infile, "fasta"):
		test_seqs.append(record)

with open(os.path.join(working_dir,'data',\
	'test_sequences.faa'),"rU") as infile:
	test_alignment = AlignIO.read(infile, "fasta")

with open(os.path.join(working_dir,"data",\
	"test_alignment.p"),"rb") as file:
	real_alignment = pickle.load(file)

## Mock
def dummy_blast(query, subj, minoverlap, mingaps):
	# should return bools and positions
	bools = [True for e in subj]
	positions = [0 for e in subj]
	max_positions = [len(e) for e in subj]
	positions.extend(max_positions)
	return bools, positions

class AlignmentTestSuite(unittest.TestCase):

	def setUp(self):
		self.true_blast = atools.blast
		atools.blast = dummy_blast
		genedir = os.path.join(working_dir,'data',\
			'test_sequences')
		seqfiles = sorted(os.listdir(genedir))
		seqfiles = [e for e in seqfiles if not re.search("^\.|^log\.txt$",\
		 e)]
		self.store = atools.SeqStore(genedir=genedir,\
			seqfiles = seqfiles, minfails = 10, mingaps = 0.5,\
			minoverlap = 50)
		self.aligner = atools.Aligner(self.store, mingaps = 0.5,\
			minoverlap = 50, minseedsize = 3, maxtrys = 10,\
			maxseedtrys = 10, gene_type = 'shallow')

	def tearDown(self):
		atools.blast = self.true_blast

	def test_version(self):
		# try different combinations of sequences
		def genSequences(n, length):
			s = 'A' * length
			return [s for i in range(n)]
		self.assertEqual(atools.version(genSequences(30, 800),'deep'),\
		 'mafft-xinsi')
		self.assertEqual(atools.version(genSequences(90, 800),'deep'),\
		 'mafft-qinsi')
		self.assertEqual(atools.version(genSequences(90, 2000),'deep'),\
		 'mafft --auto')

	def test_align(self):
		# align and check if results exist
		res = atools.align('mafft --auto', test_seqs)
		self.assertTrue(res.__class__.__name__ ==\
			'MultipleSeqAlignment')

	def test_add(self):
		# add to test_alignment and check if result
		#  exists
		res = atools.add(test_alignment, test_seqs[0])
		self.assertTrue(res.__class__.__name__ ==\
			'MultipleSeqAlignment')

	def test_blast(self):
		# real_alignment and real_sequences are identical
		real_seqs = [e for e in real_alignment]
		# choose sample of real seqs for speed
		subjseqs = random.sample(real_seqs, 5)
		queryseqs = random.sample(real_seqs, 1)
		res,_ = self.true_blast(query = queryseqs, subj = subjseqs,\
			minoverlap = 50, mingaps = 0.5)
		# all true
		self.assertTrue(all(res))

	def test_checkalignment_arg_mingaps(self):
		# check mingaps argument (proportion of internal gaps)
		# check with good alignment
		res = atools.checkAlignment(test_alignment,\
			mingaps = 0.5, minoverlap = 1, minlen = 1)
		self.assertTrue(res)
		# check with bad alignment
		# this sequence has a large internal gap covering
		#  50% of the alignment
		bad_seq = 'A' + '-' * 51 + 'A' * 48
		bad_seq = SeqRecord(Seq(bad_seq), id = 'bad')
		bad_alignment = test_alignment[:]
		bad_alignment.append(bad_seq)
		res = atools.checkAlignment(bad_alignment,\
			mingaps = 0.5, minoverlap = 1, minlen = 1)
		self.assertFalse(res)

	def test_checkalignment_arg_minoverlap(self):
		# check minoverlap argument (minimum number of 
		#	overlapping nucleotides)
		# bad seq has a large external gap -- the first 51
		#  nucleotides are missing
		bad_seq = '-' * 51 + 'A' * 49
		bad_seq = SeqRecord(Seq(bad_seq), id = 'bad')
		bad_alignment = test_alignment[:]
		bad_alignment.append(bad_seq)
		res = atools.checkAlignment(bad_alignment,\
			mingaps = 0.5, minoverlap = 50, minlen = 1)
		self.assertFalse(res)

	def test_checkalignment_arg_minlen(self):
		# check minlen argument
		res = atools.checkAlignment(test_alignment,\
			mingaps = 0.5, minoverlap = 1, minlen = 101)
		self.assertFalse(res)

	def test_seqstore_private_alignmentblast(self):
		# use real blast
		atools.blast = self.true_blast
		real_seqs = [e for e in real_alignment]
		next_seqs = random.sample(real_seqs, 5)
		seqs_in_alignment = random.sample(real_seqs, 1)
		# same as blast test but with a dodgy sequence
		bad_seq = 'T' * 100
		bad_seq = SeqRecord(Seq(bad_seq), id = 'bad')
		next_seqs_w_bad = next_seqs[:]
		next_seqs_w_bad.append(bad_seq)
		res = self.store._alignmentBlast(query = \
			next_seqs_w_bad, sequences_in_alignment =\
			seqs_in_alignment)
		# bad seq should not be returned
		self.assertFalse(5 in res[0])
		corrected_seq = bad_seq + random.sample(real_seqs, 1)[0]
		next_seqs_w_corrected = next_seqs[:]
		next_seqs_w_corrected.append(corrected_seq)
		res = self.store._alignmentBlast(query = \
			next_seqs_w_corrected, sequences_in_alignment =\
			seqs_in_alignment)
		# corrected seq should now be returned
		self.assertTrue(5 in res[0])
		# and its returned sequence shouldn't have the bad sequence
		self.assertTrue(str(res[1][-1].seq) != str(corrected_seq.seq))
		# switch back to dummy blast
		atools.blast = dummy_blast

	def test_seqstore_private_add(self):
		store = copy.deepcopy(self.store)
		# add lists to obj for add to work
		store.sppool = store.keys()
		store.sequences_in_alignment = []
		store._add()
		res = store.sequences_in_alignment[0]
		# the species should no longer be in the pool
		self.assertFalse(res[0].id in store.sppool)

	def test_seqstore_start(self):
		store = copy.deepcopy(self.store)
		seqs = store.start(3)
		# none of the species should be in the pool
		res = [e.id in store.sppool for e in seqs]
		self.assertFalse(all(res))
		# each sequence returned should be from different species
		res = list(set([e.id for e in seqs]))
		self.assertEqual(len(res), 3)

	def test_seqstore_back(self):
		store = copy.deepcopy(self.store)
		before_seqs = store.start(5)
		after_seq = store.back(before_seqs[:-1])
		# after_seq sp should not be in before_seqs
		res = [after_seq.id == e.id for e in before_seqs]
		self.assertFalse(all(res))
		# last before_seq should have penalty
		res = store[before_seqs[-1].id][0]
		res = [e[1] for e in res]
		self.assertEqual(sum(res), 1)

	def test_seqstore_next(self):
		store = copy.deepcopy(self.store)
		before_seqs = store.start(5)
		# add a mock penalty to seqs in alignment
		seqs = store.sequences_in_alignment
		for each in seqs:
			each[1] += 1
		after_seq = store.next(test_alignment)
		# after_seq sp should not be in before_seqs
		res = [after_seq.id == e.id for e in before_seqs]
		self.assertFalse(all(res))
		# after_seq sp should not be in sppool
		res = [after_seq.id == e for e in\
		store.sppool]
		self.assertFalse(all(res))
		# no penalties in seqs in alignment
		seqs = store.sequences_in_alignment
		res = [e[1] == 0 for e in seqs]
		self.assertTrue(all(res))

	def test_seqstore_private_check_dropsequence(self):
		# add maxfails to a sequence
		store = copy.deepcopy(self.store)
		seq_to_drop = store['sp1'][0][0][0].name
		store['sp1'][0][0][1] = 11
		store._check()
		res = [e[0].name for e in store['sp1'][0]]
		# ._check() should have removed it
		self.assertFalse(seq_to_drop in res)

	def test_seqstore_private_check_dropspecies(self):
		# add maxfails to all seqs of one species
		store = copy.deepcopy(self.store)
		# add a sppool, so that dropped species can
		#  removed
		store.sppool = store.keys()
		sp_to_drop = 'sp1'
		for seq in store[sp_to_drop][0]:
			seq[1] = 11
		store._check()
		# ._check() should have removed it
		self.assertFalse(sp_to_drop in store.keys())

	def test_aligner_return(self):
		# create different types of alignments and see what
		#  is returned from _return
		outgroup_seq = SeqRecord(Seq('T' * 3141), id = 'outgroup')
		real_alignment_w_out = real_alignment[:]
		real_alignment_w_out.append(outgroup_seq)
		outgroup_seq = SeqRecord(Seq('T' * 100), id = 'outgroup')
		test_alignment_w_out = test_alignment[:]
		test_alignment_w_out.append(outgroup_seq)
		mock_multi_alignment = [real_alignment,test_alignment,\
		real_alignment_w_out,real_alignment[:5], test_alignment_w_out]
		self.aligner.store = mock_multi_alignment
		alignment = self.aligner._return()
		# should return real_alignment_w_out (has 12 records)
		# it has an outgroup and it is the longest
		self.assertEqual(len(alignment), 12)

	def test_aligner_run(self):
		# all the outgroup seqs should not align, alignment
		#  cannot have outgroup species
		res = self.aligner.run()
		spp = [e.id for e in res]
		self.assertFalse('outgroup' in spp)

if __name__ == '__main__':
    unittest.main()