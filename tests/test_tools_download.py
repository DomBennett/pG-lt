"""
Tests for Download tools.
"""

import unittest
import mpe.tools.download as dtools

# Mock data
taxids = ['1', '2']
gene_names = ['name1', 'name2']
nseqs = 100
thoroughness = 3
maxpn = 0.1
seedsize = 3
mingaps = 0.01
minoverlap = 200
maxtrys = 100
maxlen = 2000
minlen = 300
downloader = dtools.Downloader (gene_names = gene_names,\
			nseqs = nseqs, thoroughness = thoroughness,\
			maxpn = maxpn, seedsize = seedsize, maxtrys \
			= maxtrys, mingaps = mingaps, minoverlap = minoverlap,\
			maxlen = maxlen, minlen = minlen)
t1_term = '(txid1[PORGN]) OR txid2[PORGN] AND (("name1"[GENE]) \
OR "name2"[GENE]) NOT predicted[TI] NOT genome[TI] NOT unverified[TI]'
t1_search_res = {'Count':0}
t2_term = '(txid1[PORGN]) OR txid2[PORGN] AND (("name1"[TI]) \
OR "name2"[TI]) NOT predicted[TI] NOT genome[TI] NOT unverified[TI]'
t2_search_res = {'Count':2, 'IdList':['seq1', 'seq2']}
t3_term = '(txid1[PORGN]) OR txid2[PORGN] AND (("name1") \
OR "name2") NOT predicted[TI] NOT shotgun[TI] NOT scaffold[TI] \
NOT assembly[TI] NOT unverified[TI]'
t3_search_res = {'Count':3, 'IdList':['seq1', 'seq2', 'seq3']}

#Stubs
def dummry_eSearch (term, retStart=0, retMax=1, usehistory="n",\
	db = "nucleotide"):
	if term == t1_term:
		return t1_search_res
	if term == t2_term:
		return t2_search_res
	if term == t3_term:
		return t3_search_res
def dummy_align (sequences):
	pass
def dummy_checkAlignment(alignment, mingaps, minoverlap, minlen):
	return True
dtools.etools.eSearch = dummry_eSearch
dtools.atools.align = dummy_align
dtools.atools.checkAlignment = dummy_checkAlignment

class DownloadTestSuite(unittest.TestCase):

	def setUp(self):
		# mock Downloader instance
		self.downloader = downloader
		# expected search terms at different thoroughnesses
		self.t1_term = t1_term
		self.t2_term = t2_term
		self.t3_term = t3_term
		self.taxids = taxids

	def test_downloader_private_buildsearchterm_thoroughness1(self):
		res = self.downloader._buildSearchTerm(self.taxids, 1)
		self.assertEqual(res, self.t1_term)

	def test_downloader_private_buildsearchterm_thoroughness2(self):
		res = self.downloader._buildSearchTerm(self.taxids, 2)
		self.assertEqual(res, self.t2_term)

	def test_downloader_private_buildsearchterm_thoroughness3(self):
		res = self.downloader._buildSearchTerm(self.taxids, 3)
		self.assertEqual(res, self.t3_term)

	def test_downloader_private_search(self):
		# expect to only find 2, 3 and 0 sequences
		res1 = self.downloader._search(self.taxids)
		res2 = self.downloader._search(self.taxids)
		res3 = self.downloader._search(self.taxids)
		self.assertEqual([len(res1),len(res2),len(res3)], [2, 1, 0])

	def test_downloader_private_filter(self):
		pass

if __name__ == '__main__':
    unittest.main()