"""
Tests for Phylogeny tools.
"""

import unittest,pickle,os
import mpe.tools.phylogeny as ptools

with open(os.path.join("data","test_alignment.p"),\
	"rb") as file:
	alignment = pickle.load(file)

# TODO: create dummy phylogeny
# TODO: create dummy alignments
# TODO: create taxonomic tree

partitions = [0, 1760, 3141]

class DownloadTestSuite(unittest.TestCase):

	def setUp(self):
		self.partitions = partitions
		self.alignment = alignment

	def test_getbranchlengths(self):
		pass

	def test_goodphylogenytest(self):
		pass

	def test_genconstrainttree(self):
		pass

	def test_concatenatealignments(self):
		pass

	def test_raxml(self):
		# TODO: add outgroup + constraint
		phylo = ptools.RAxML(self.alignment,\
			partitions = self.partitions)
		self.assertTrue(phylo)

if __name__ == '__main__':
    unittest.main()