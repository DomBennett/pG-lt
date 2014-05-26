#! /bin/usr/env python
## D.J. Bennett
## 26/05/2014
"""
Tests for phylogeny tools.
"""

import unittest,pickle,os
from Bio import Phylo
import mpe.tools.phylogeny as ptools

## Dirs
working_dir = os.path.dirname(__file__)

## Mock data
with open(os.path.join(working_dir,"data",\
	"test_alignment.p"),"rb") as file:
	alignment = pickle.load(file)

with open(os.path.join(working_dir,"data",\
	"test_alignments.p"),"rb") as file:
	alignments = pickle.load(file)

with open(os.path.join(working_dir,"data",\
	"test_phylo.p"),"rb") as file:
	phylo = pickle.load(file)

partitions = [0, 1761, 3141]

class PhylogenyTestSuite(unittest.TestCase):

	def setUp(self):
		self.partitions = partitions
		self.alignment = alignment
		self.phylo = phylo
		self.alignment = alignment
		self.alignments = alignments
		self.constraint_opt = ' -g constraint.tre'

	def test_getbranchlengths(self):
		# 15 tips, one in-group, one root
		# lengths are all scaled to 1
		# total branch length should equal 17
		res = ptools.getBranchLengths(self.phylo)
		self.assertEqual(sum(res), 17)

	def test_goodphylogenytest(self):
		res = ptools.goodPhylogenyTest(self.phylo, 0.5)
		self.assertTrue(res)

	def test_concatenatealignments(self):
		res_alignment,res_partitions = \
			ptools.concatenateAlignments(self.alignments)
		self.assertEqual(res_partitions, self.partitions)
		seq = self.alignment[0]
		for each in res_alignment:
			if each.id == seq.id:
				res_seq = each
		self.assertEqual(str(res_seq.seq), str(seq.seq))
		
	def test_genconstrainttree(self):
		# 4 tips not present in alignment
		res_opt = ptools.genConstraintTree(self.alignment,\
			os.path.join(working_dir, 'data',\
				'test_phylo.tre'))
		with open("constraint.tre", "r") as file:
			res_tree = Phylo.read(file, "newick")
		self.assertEqual(len(res_tree.get_terminals()), 11)
		self.assertEqual(res_opt, self.constraint_opt)

	def test_raxml(self):
		phylo = ptools.RAxML(self.alignment,\
			partitions = self.partitions, outgroup = \
			"Ignatius_tetrasporus", constraint = self.constraint_opt)
		self.assertTrue(phylo)

if __name__ == '__main__':
	unittest.main()