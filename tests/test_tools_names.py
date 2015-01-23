#! /bin/usr/env python
# D.J. Bennett
# 26/05/2014
"""
Tests for names tools.
"""

import unittest
import pickle
import os
import json
import pglt.tools.names_tools as ntools
import taxon_names_resolver as tnr

ntools.etools.Entrez.email = "python.unittests@pglt.program"

# DIRS
working_dir = os.path.dirname(__file__)

# TEST DATA
with open(os.path.join(working_dir, 'data', 'test_search.json'), 'r') as file:
    res = json.load(file)

with open(os.path.join(working_dir, 'data', 'test_namesdict.p'), 'r') as file:
    exp_namesdict = pickle.load(file)

# create one with and one without an outgroup
exp_namesdict_wo = exp_namesdict.copy()
del exp_namesdict['outgroup']
terms = ['GenusA speciesA', 'GenusA speciesB', 'GenusA speciesC',
         'GenusB speciesD', 'GenusB speciesE', 'GenusC speciesF',
         'GenusD speciesG', 'GenusE speciesH', 'GenusF speciesI',
         'GenusG speciesJ']
exp_allrankids = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
                  21, 22, 23, 24, 25, 26, 31, 32, 33, 34, 35, 41, 42, 43, 44,
                  51]
exp_parentid = 51


# STUBS
class dummy_Logger(object):

    def __init__(self):
        pass

    def info(self, msg):
        pass

    def debug(self, msg):
        pass

    def warn(self, msg):
        pass


def dummy_eFetch(taxid, db, logger):
    # rank and sci name for the outgroup when getmetadata in findbestgenes
    return [{'ParentTaxId': '1', 'Rank': 'outgroup',
             'ScientificName': 'outgroup'}]


def dummy_eSearch(term, logger):
    return {'Count': '10000'}


def dummy_findChildren(taxid, next, logger):
    return [60]


class Dummy_GnrDataSources(object):

    def __init__(self, logger):
        pass

    def byName(self, names, invert=False):
        if invert:
            return [1, 2, 3]
        else:
            return [4]


class NamesTestSuite(unittest.TestCase):

    def setUp(self):
        self.logger = dummy_Logger()
        self.True_GnrDataSources = tnr.gnr_tools.GnrDataSources
        self.true_eFetch = ntools.etools.eFetch
        self.true_eSearch = ntools.etools.eSearch
        self.true_findChildren = ntools.etools.findChildren
        tnr.gnr_tools.GnrDataSources = Dummy_GnrDataSources
        ntools.etools.eFetch = dummy_eFetch
        ntools.etools.eSearch = dummy_eSearch
        ntools.etools.findChildren = dummy_findChildren
        self.resolver = tnr.resolver.Resolver(terms=terms, taxon_id=51,
                                              logger=self.logger)
        # add results to resolver
        test_store = tnr.gnr_tools.GnrStore(terms, logger=self.logger)
        test_store.add(res)
        self.resolver._store = test_store

    def tearDown(self):
        if os.path.isdir('resolved_names'):
            os.rmdir('resolved_names')
        if os.path.isfile('resolved_names.csv'):
            os.remove('resolved_names.csv')
        tnr.gnr_tools.GnrDataSources = self.True_GnrDataSources
        ntools.etools.eFetch = self.true_eFetch
        ntools.etools.findChildren = self.true_findChildren
        ntools.etools.eSearch = self.true_eSearch

    def test_gennamesdict(self):
        namesdict, allrankids, parentid = \
            ntools.genNamesDict(self.resolver, logger=self.logger)
        self.assertEqual(allrankids, exp_allrankids)
        self.assertEqual(namesdict, exp_namesdict)
        self.assertEqual(parentid, exp_parentid)

    def test_getoutgroup(self):
        namesdict = ntools.getOutgroup(exp_namesdict.copy(), exp_parentid,
                                       logger=self.logger)
        self.assertEqual(namesdict, exp_namesdict_wo)

    def test_gentaxtree(self):
        tree = ntools.genTaxTree(self.resolver, exp_namesdict,
                                 logger=self.logger, taxonomy=None, draw=False)
        self.assertTrue(tree)

    def test_write_namesdict(self):
        ntools.writeNamesDict(directory='.', namesdict=exp_namesdict_wo)
        self.assertTrue(os.path.isfile('resolved_names.csv'))

if __name__ == '__main__':
    unittest.main()
