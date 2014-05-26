#! /bin/usr/env python
## D.J. Bennett
## 26/05/2014
"""
Tests for alignment tools.
"""

import unittest,os,re,copy,pickle,random
import mpe.tools.alignment as atools
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
class Dummy_SeqStore(atools.SeqStore):
	pass

def dummy_blast(self, alignment, sequences):
	# would normally return the seq indexes for sequences
	#  that overlap with alignment according to BLAST.
	#  Instead return all indexess.
	return range(len(sequences))

Dummy_SeqStore._blast = dummy_blast

class AlignmentTestSuite(unittest.TestCase):

	def setUp(self):
		genedir = os.path.join(working_dir,'data',\
			'test_sequences')
		seqfiles = sorted(os.listdir(genedir))
		seqfiles = [e for e in seqfiles if not re.search("^\.|^log\.txt$",\
		 e)]
		self.store = atools.SeqStore(genedir=genedir,\
			seqfiles=seqfiles, minfails = 10, mingaps = 0.5,\
			minoverlap = 50, verbose = False)
		self.dummy_store = Dummy_SeqStore(genedir=genedir,\
			seqfiles=seqfiles, minfails = 10, mingaps = 0.5,\
			minoverlap = 50, verbose = False)
		self.aligner = atools.Aligner(self.dummy_store, mingaps = 0.5,\
			minoverlap = 50, minseedsize = 3, maxtrys = 10,\
			maxseedtrys = 10, verbose = False)

	def tearDown(self):
		pass

	def test_align(self):
		# align and check if results exist
		res = atools.align(test_seqs)
		self.assertTrue(res.__class__.__name__ ==\
			'MultipleSeqAlignment')

	def test_add(self):
		# add to test_alignment and check if result
		#  exists
		res = atools.add(test_alignment, test_seqs[0])
		self.assertTrue(res.__class__.__name__ ==\
			'MultipleSeqAlignment')

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

	def test_seqstore_private_add(self):
		store = copy.deepcopy(self.store)
		# add lists to obj for add to work
		store.sppool = store.keys()
		store.sequences_in_alignment = []
		store._add()
		res = store.sequences_in_alignment[0]
		# the species should no longer be in the pool
		self.assertFalse(res[0].id in store.sppool)
		## Additional if I implement popping 26/05/2014
		# the seq should no longer be in the species' seqs
		#spseqs = store[res[0].id][0]
		#spseqnames = [e[0].name for e in spseqs]
		#self.assertFalse(res[0].name in spseqnames)

	def test_seqstore_start(self):
		store = copy.deepcopy(self.store)
		seqs = store.start(3)
		# none of the species should be in the pool
		res = [e.id in store.sppool for e in seqs]
		self.assertFalse(all(res))
		# each sequence returned should be from different species
		res = list(set([e.id for e in seqs]))
		self.assertEqual(len(res), 3)

	def test_seqstore_private_blastadd(self):
		# use dummy store to avoid blasting
		store = copy.deepcopy(self.dummy_store)
		before_seqs = store.start(5)
		del before_seqs
		res = store._blastAdd(test_alignment)
		# should be true, all sequence indexes are
		#  returned with dummy
		self.assertTrue(res)

	def test_seqstore_back(self):
		store = copy.deepcopy(self.store)
		before_seqs = store.start(5)
		after_seq = store.back()
		# after_seq sp should not be in before_seqs
		res = [after_seq.id == e.id for e in before_seqs]
		self.assertFalse(all(res))
		# last before_seq should have penalty
		res = store[before_seqs[-1].id][0]
		res = [e[1] for e in res]
		self.assertEqual(sum(res), 1)

	def test_seqstore_next(self):
		# use dummy store to avoid blasting
		store = copy.deepcopy(self.dummy_store)
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

	def test_seqstore_private_blast(self):
		# real_alignment and real_sequences are identical
		real_seqs = []
		for i in range(len(real_alignment)):
			real_seqs.append(real_alignment[i])
		# choose sample of real seqs for speed
		real_seqs = random.sample(real_seqs, 4)
		res = self.store._blast(real_alignment, real_seqs)
		# all indexes should be returned
		self.assertEqual(res, range(len(real_seqs)))
		bad_seq = 'T' * 100
		bad_seq = SeqRecord(Seq(bad_seq), id = 'bad')
		real_seqs_w_bad = real_seqs[:]
		real_seqs_w_bad.append(bad_seq)
		res = self.store._blast(real_alignment,\
			real_seqs_w_bad)
		# bad seq should not be returned
		self.assertFalse(len(real_seqs_w_bad) in res)

	def test_seqstore_private_check_dropseqeunce(self):
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
		alignment = self.aligner._return(mock_multi_alignment)
		# should return real_alignment_w_out (has 12 records)
		# it has an outgroup and it is the longest
		self.assertEqual(len(alignment), 12)

	def test_aligner_run(self):
		# all the outgroup seqs should not align, should return None
		# I used a small seedsize to test both phases: seed and add
		res = self.aligner.run()
		self.assertFalse(res)

if __name__ == '__main__':
    unittest.main()