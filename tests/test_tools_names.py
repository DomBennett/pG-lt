#! /bin/usr/env python
## D.J. Bennett
## 26/05/2014
"""
Tests for names tools.
"""

import unittest,pickle,os,json
import mpe.tools.names as ntools
import taxon_names_resolver as tnr

ntools.etools.Entrez.email = "python.unittests@mpe.program"

## Dirs
working_dir = os.path.dirname(__file__)

## Test data
with open(os.path.join(working_dir,'data',\
	'test_search.json'),'r') as file:
	res = json.load(file)

with open(os.path.join(working_dir,\
	'data','test_namesdict.p'),'r') as file:
	exp_namesdict = pickle.load(file)

terms = ['GenusA speciesA', 'GenusA speciesB', 'GenusA speciesC',\
'GenusB speciesD','GenusB speciesE','GenusC speciesF','GenusD \
speciesG','GenusE speciesH','GenusF speciesI','GenusG speciesJ']
exp_allrankids = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,\
15, 16, 17, 21, 22, 23, 24, 25, 26, 31, 32, 33, 34, 35, 41, 42, 43,\
44, 51]

## Stubs
def dummy_eFetch(taxid, db):
	return [{'ParentTaxId':1}]


def dummy_findChildren(taxid, next):
	return [60]


class Dummy_GnrDataSources(object):
	def __init__(self):
		pass
	def byName(self, names, invert = False):
		if invert:
			return [1, 2, 3]
		else:
			return [4]

class NamesTestSuite(unittest.TestCase):

	def setUp(self):
		self.True_GnrDataSources = tnr.gnr_tools.GnrDataSources
		self.true_eFetch = ntools.etools.eFetch
		self.true_findChildren = ntools.etools.findChildren
		tnr.gnr_tools.GnrDataSources = Dummy_GnrDataSources
		ntools.etools.eFetch = dummy_eFetch
		ntools.etools.findChildren = dummy_findChildren
		self.resolver = tnr.resolver.Resolver(terms = terms, taxon_id\
			= 51)
		test_store = tnr.gnr_tools.GnrStore(terms)
		test_store.add(res)
		self.resolver._store = test_store

	def tearDown(self):
		tnr.gnr_tools.GnrDataSources = self.True_GnrDataSources
		ntools.etools.eFetch = self.true_eFetch
		ntools.etools.findChildren = self.true_findChildren

	def test_gennamesdict(self):
		namesdict,allrankids = ntools.genNamesDict(self.resolver)
		self.assertEqual(allrankids, exp_allrankids)
		self.assertEqual(namesdict, exp_namesdict)

	def test_gentaxtree(self):
		tree,line = ntools.genTaxTree(self.resolver, exp_namesdict,\
			draw = False)
		self.assertTrue(tree)
		self.assertEqual(line, [u'51'])

if __name__ == '__main__':
    unittest.main()