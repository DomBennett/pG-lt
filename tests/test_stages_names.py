#! /bin/usr/env python
## D.J. Bennett
## 26/05/2014
"""
Test names stage.
"""

import unittest,pickle,os,shutil
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
	return namesdict, [], None

def dummy_getOutgroup(namesdict, parentid, outgroupid):
	namesdict['outgroup'] = {"txids" : [3],\
	"unique_name" : 'outgroup', "rank" : 'genus'}
	return namesdict

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
		self.true_getOutgroup = names_stage.ntools.getOutgroup
		names_stage.Resolver = Dummy_Resolver
		names_stage.ntools.genNamesDict = dummy_genNamesDict
		names_stage.ntools.getOutgroup = dummy_getOutgroup
		names_stage.ntools.genTaxTree = dummy_genTaxTree
		# write out necessary files to run
		paradict = {'email' : '', 'parentid' : '', 'outgroupid': ''}
		with open(".paradict.p", "wb") as file:
			pickle.dump(paradict, file)
		with open(".terms.p", "wb") as file:
			pickle.dump(['t1', 't2', 't3', 't4', 't5'], file)
		os.mkdir('resolved_names')

	def tearDown(self):
		# stub in
		names_stage.Resolver = self.True_Resolver
		names_stage.ntools.genNamesDict = self.true_genNamesDict
		names_stage.ntools.genTaxTree = self.true_genTaxTree
		names_stage.ntools.getOutgroup = self.true_getOutgroup
		# remove all files potentially generated by names stage
		names_files = ['.paradict.p', '.terms.p', '.allrankids.p',\
		'.namesdict.p']
		while names_files:
			try:
				names_file = names_files.pop()
				os.remove(names_file)
			except OSError:
				pass
		# remove all folders potentially generated by names stage
		names_folders = ['1_names', '4_phylogeny']
		while names_folders:
			try:
				names_folder = names_folders.pop()
				shutil.rmtree(names_folder)
			except OSError:
				pass

	def test_names_stage(self):
		# run
		res = names_stage.run()
		# nothing is returned
		self.assertIsNone(res)

if __name__ == '__main__':
    unittest.main()