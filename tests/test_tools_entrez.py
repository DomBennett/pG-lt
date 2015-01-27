#! /bin/usr/env python
# D.J. Bennett
# 26/05/2014
"""
Tests for Entrez tools.
"""

# PACKAGES
import unittest
import pglt.tools.entrez_tools as etools


# GLOBALS
etools.Entrez.email = "entrez.unittests@pglt.program"


# FUNCTIONS
def foo(arg1):
    '''Simple function to test safeConnect'''
    raise IOError()
    return arg1


# DUMMIES
class dummy_Logger(object):

    def __init__(self):
        pass

    def info(self, msg):
        pass

    def debug(self, msg):
        pass


class EntrezTestSuite(unittest.TestCase):

    def setUp(self):
        self.logger = dummy_Logger()

    def test_safeconnect(self):
        # test safeconnect by passing it a failing function
        res = etools.safeConnect(efunc=foo, logger=self.logger, arg1='arg1',
                                 waittime=0.1, max_check=1)
        # if fails to connect, returns ()
        self.assertEqual(res, ())

    def test_efetch_taxonomy(self):
        # 9606 is humans
        res = etools.eFetch(ncbi_id='9606', db='taxonomy', logger=self.logger)
        self.assertEqual(res[0]['ScientificName'], 'Homo sapiens')

    def test_efetch_nucleotide(self):
        # should return a list of SeqRecords
        res = etools.eFetch(ncbi_id='DQ373091', logger=self.logger)
        self.assertEqual(res[0].__class__.__name__, 'SeqRecord')

    def test_esearch_nucleotide(self):
        # there are over 500 caulimovirus sequences
        term = 'caulimovirus [PORGN]'
        res = etools.eSearch(term, logger=self.logger)
        self.assertGreater(int(res['Count']), 500)

    def test_findchildren_next_is_true(self):
        # Catarrhines have 2 clades ranked below it: apes and old world monkeys
        res = etools.findChildren(9526, next=True, logger=self.logger)
        self.assertEqual(len(res), 2)

    def test_findchildren_next_is_false(self):
        # there are 8 ape genera
        res = etools.findChildren(314295, logger=self.logger)
        self.assertEqual(len(res), 8)

if __name__ == '__main__':
    unittest.main()
