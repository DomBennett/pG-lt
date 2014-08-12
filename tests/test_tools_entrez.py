#! /bin/usr/env python
## D.J. Bennett
## 26/05/2014
"""
Tests for Entrez tools.
"""

import unittest
import mpe.tools.entrez_tools as etools

etools.Entrez.email = "python.unittests@mpe.program"

class EntrezTestSuite(unittest.TestCase):

	def test_efetch_taxonomy(self):
		# 9606 is humans
		res = etools.eFetch(ncbi_id = '9606', db = 'taxonomy')
		self.assertEqual(res[0]['ScientificName'], 'Homo sapiens')

	def test_efetch_nucleotide(self):
		# should return a list of SeqRecords
		res = etools.eFetch(ncbi_id = 'DQ373091')
		self.assertEqual(res[0].__class__.__name__, 'SeqRecord')

	def test_esearch_nucleotide(self):
		# there are over 500 caulimovirus sequences
		term = 'caulimovirus [PORGN]'
		res = etools.eSearch(term)
		self.assertGreater(int(res['Count']), 500)

	def test_findchildren_next_is_true(self):
		# eutheria has 3 clades ranked below it
		res = etools.findChildren(9347, next = True)
		self.assertEqual(len(res), 3)

	def test_findchildren_next_is_false(self):
		# there are 8 ape genera
		res = etools.findChildren(314295)
		self.assertEqual(len(res), 8)

if __name__ == '__main__':
    unittest.main()