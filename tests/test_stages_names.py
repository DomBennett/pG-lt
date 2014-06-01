#! /bin/usr/env python
## D.J. Bennett
## 26/05/2014
"""
Test names stage.
"""

import unittest,pickle,os
from mpe.stages import names_stage
from Bio import Phylo
from cStringIO import StringIO

## Dummies
class Dummy_Resolver(object):
	def __init__(self, terms, datasource, taxon_id):
		pass
	def main(self):
		pass

def dummy_genNamesDict(resolver):
	namesdict = {}
	namesdict['query_name'] = {"txids" : [1,2],\
	"unique_name" : 'returned_name', "rank" : 'species'}
	return namesdict, None

def dummy_genTaxTree(resolver, namesdict):
	treedata = "(A, (B, C), (D, E))"
	handle = StringIO(treedata)
	tree = Phylo.read(handle, "newick")
	return tree, None

class NamesStageTestSuite(unittest.TestCase):

	def setUp(self):
		# stub functions and class
		self.True_Resolver = names_stage.Resolver
		self.true_genNamesDict = names_stage.ntools.genNamesDict
		self.true_genTaxTree = names_stage.ntools.genTaxTree
		names_stage.Resolver = Dummy_Resolver
		names_stage.ntools.genNamesDict = dummy_genNamesDict
		names_stage.ntools.genTaxTree = dummy_genTaxTree
		# write out necessary files to run
		paradict = {'email' : '', 'parentid' : ''}
		with open(".paradict.p", "wb") as file:
			pickle.dump(paradict, file)
		with open(".terms.p", "wb") as file:
			pickle.dump({}, file)
		os.mkdir('resolved_names')

	def tearDown(self):
		names_stage.Resolver = self.True_Resolver
		names_stage.ntools.genNamesDict = self.true_genNamesDict
		names_stage.ntools.genTaxTree = self.true_genTaxTree

	def test_names_stage(self):
		# run
		res = names_stage.run()
		# remove files and folders
		os.remove('.paradict.p')
		os.remove('.terms.p')
		os.remove(os.path.join('1_names','resolved_names.csv'))
		os.rmdir('1_names')
		os.remove(os.path.join('4_phylogeny','taxontree.tre'))
		os.rmdir('4_phylogeny')
		# nothing is returned
		self.assertIsNone(res)

if __name__ == '__main__':
    unittest.main()