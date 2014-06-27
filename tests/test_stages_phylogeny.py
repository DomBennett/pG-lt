#! /bin/usr/env python
## D.J. Bennett
## 26/05/2014
"""
Test phylogeny stage.
"""

import unittest,os,shutil,pickle
from mpe.stages import phylogeny_stage
from Bio import AlignIO
from Bio import Phylo
from cStringIO import StringIO

## Dirs
working_dir = os.path.dirname(__file__)

## Dummies
def dummy_concatenateAlignments(alignments):
	return None,None

def dummy_genConstraintTree(alignment, path):
	pass

def dummy_RAxML(alignment,constraint,outgroup,partitions):
	treedata = "(outgroup, (B, C), (D, E))"
	handle = StringIO(treedata)
	tree = Phylo.read(handle, "newick")
	return tree
	
def dummy_goodPhylogenyTest(phylogeny, maxpedge):
	return True

## Test data
with open(os.path.join(working_dir, 'data','test_alignment_ref.faa'), 'r') \
as file:
	alignment = AlignIO.read(file, 'fasta')
paradict = {'ntrees':1,'maxtrys':1,'maxpedge':0.5}


class PhylogenyStageTestSuite(unittest.TestCase):

	def setUp(self):
		# stub out
		self.true_concatenateAlignments = phylogeny_stage.ptools.concatenateAlignments
		self.true_genConstraintTree = phylogeny_stage.ptools.genConstraintTree
		self.true_RAxML = phylogeny_stage.ptools.RAxML
		self.true_goodPhylogenyTest = phylogeny_stage.ptools.goodPhylogenyTest
		phylogeny_stage.ptools.concatenateAlignments = dummy_concatenateAlignments
		phylogeny_stage.ptools.genConstraintTree = dummy_genConstraintTree
		phylogeny_stage.ptools.RAxML = dummy_RAxML
		phylogeny_stage.ptools.goodPhylogenyTest = dummy_genConstraintTree
		# create input data
		with open(".paradict.p", "wb") as file:
			pickle.dump(paradict, file)
		os.mkdir('3_alignment')
		os.mkdir('4_phylogeny')
		os.mkdir(os.path.join('3_alignment','COI'))
		os.mkdir(os.path.join('3_alignment','rbcl'))
		with open(os.path.join('3_alignment', 'rbcl',\
			'test_alignment_rbl.faa'), 'w') as file:
			count = AlignIO.write(alignment, file, "fasta")
			del count
		with open(os.path.join('3_alignment', 'COI',\
			'test_alignment_COI.faa'), 'w') as file:
			count = AlignIO.write(alignment, file, "fasta")
			del count

	def tearDown(self):
		os.remove('.paradict.p')
		shutil.rmtree('3_alignment')
		#stub in
		phylogeny_stage.ptools.concatenateAlignments = self.true_concatenateAlignments
		phylogeny_stage.ptools.genConstraintTree = self.true_genConstraintTree
		phylogeny_stage.ptools.RAxML = self.true_RAxML
		phylogeny_stage.ptools.goodPhylogenyTest = self.true_goodPhylogenyTest

	def test_phylogeny_stage(self):
		# run
		res = phylogeny_stage.run()
		# clean dir
		os.remove(os.path.join('4_phylogeny', 'distribution.tre'))
		os.remove(os.path.join('4_phylogeny', 'consensus.tre'))
		os.rmdir('4_phylogeny')
		# assert
		self.assertIsNone(res)

if __name__ == '__main__':
    unittest.main()