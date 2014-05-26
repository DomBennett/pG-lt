#! /bin/usr/env python
## D.J. Bennett
## 26/05/2014
"""
Tests for download tools.
"""

import unittest,pickle,os
import mpe.tools.download as dtools

## Dirs
working_dir = os.path.dirname(__file__)

## Dummy data and stubs

# expected terms and term search results
t1_term = '(txid1[PORGN]) OR txid2[PORGN] AND (("name1"[GENE]) \
OR "name2"[GENE]) NOT predicted[TI] NOT genome[TI] NOT unverified[TI]'
t2_term = '(txid1[PORGN]) OR txid2[PORGN] AND (("name1"[TI]) \
OR "name2"[TI]) NOT predicted[TI] NOT genome[TI] NOT unverified[TI]'
t3_term = '(txid1[PORGN]) OR txid2[PORGN] AND (("name1") \
OR "name2") NOT predicted[TI] NOT shotgun[TI] NOT scaffold[TI] \
NOT assembly[TI] NOT unverified[TI]'
t1_search_res = {'Count':0}
t2_search_res = {'Count':2, 'IdList':['seq1', 'seq2']}
t3_search_res = {'Count':3, 'IdList':['seq1', 'seq2', 'seq3']}

# Example seqrecord for findgeneinseq
with open(os.path.join(os.path.dirname(__file__),'data',\
	"test_findgeneinseq_examplesequence.p"),"rb") as file:
	sequence = pickle.load(file)

# Dummy seq records for download
class dummy_Seq(object):
	def __init__(self):
		pass
	def tostring(self):
		# Just for parsing
		return "A" * 500

class dummy_SeqRecord(object):
	def __init__(self, description, length = 500):
		self.description = description
		self.length = length
		self.seq = dummy_Seq()
	def __len__(self):
		return self.length

seq1 = dummy_SeqRecord(description = "A sequence of NAME1")
seq2 = dummy_SeqRecord(description = "A sequence of NAME2")
seq3 = [dummy_SeqRecord(description = "A sequence of NAME3"),\
	dummy_SeqRecord(description = "A sequence of NAME4"),\
	dummy_SeqRecord(description = "A sequence of NAME5"),\
	dummy_SeqRecord(description = "A sequence of NAME1")]


# Sequences -- just Ts and Fs -- for testing filter
sequences = [True for i in range(80)]
sequences.extend([False for i in range(20)])

# Dependent stubs
def dummy_eSearch(term, retStart=0, retMax=1, usehistory="n",\
	db = "nucleotide"):
	if term == t1_term:
		return t1_search_res
	if term == t2_term:
		return t2_search_res
	if term == t3_term:
		return t3_search_res

def dummy_eFetch(ncbi_id, db = "nucleotide"):
	if ncbi_id == 'seq1':
		return seq1
	elif ncbi_id == 'seq2':
		return seq2
	elif ncbi_id == 'seq3':
		return seq3
	else:
		# return all as list
		return [seq1, seq2, seq3]

def dummy_align(sequences):
	return all (sequences)

def dummy_checkAlignment(alignment, mingaps, minoverlap, minlen):
	return alignment

# downloader init variables
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

# dictionary variables
namesdict = {"species1":{'txids':['1','2']}}
allrankids = [1, 2, 3]
genedict = {'gene1':{'taxid':'3','names':['name1', 'name2']}}

class DownloadTestSuite(unittest.TestCase):

	def setUp(self):
		self.true_eSearch = dtools.etools.eSearch
		self.true_eFetch = dtools.etools.eFetch
		self.true_align = dtools.atools.align
		self.true_checkAlignment = dtools.atools.checkAlignment
		dtools.etools.eSearch = dummy_eSearch
		dtools.etools.eFetch = dummy_eFetch
		dtools.atools.align = dummy_align
		dtools.atools.checkAlignment = dummy_checkAlignment
		# mock Downloader instance
		self.downloader = dtools.Downloader (gene_names = gene_names,\
			nseqs = nseqs, thoroughness = thoroughness,\
			maxpn = maxpn, seedsize = seedsize, maxtrys \
			= maxtrys, mingaps = mingaps, minoverlap = minoverlap,\
			maxlen = maxlen, minlen = minlen)
		# expected search terms at different thoroughnesses
		self.t1_term = t1_term
		self.t2_term = t2_term
		self.t3_term = t3_term
		self.seqids = ['seq1', 'seq2', 'seq3']
		self.seq1 = seq1
		self.seq2 = seq2
		self.seq3 = seq3
		self.sequences = sequences
		self.taxids = taxids
		self.record = sequence
		self.namesdict = namesdict
		self.allrankids = allrankids
		self.genedict = genedict

	def tearDown(self):
		#repatch
		dtools.etools.eSearch = self.true_eSearch
		dtools.etools.eFetch = self.true_eFetch
		dtools.atools.align = self.true_align
		dtools.atools.checkAlignment = self.true_checkAlignment

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
		# weeds out Falses from sequences, should not be any Falses
		sequences = self.sequences[:]
		res_filtered,res_downloaded = self.downloader._filter(sequences)
		self.assertTrue(all(res_filtered))

	def test_downloader_private_findgeneinseq(self):
		# change gene names for test
		gene_names = self.downloader.gene_names
		self.downloader.gene_names = ['COI']
		res = self.downloader._findGeneInSeq(self.record)
		self.downloader.gene_names = gene_names
		# I know that the COI sequence is 1545bp (5350..6894)
		# (http://www.ncbi.nlm.nih.gov/nuccore/AM711897.1)
		self.assertEqual(len(res), 1545)

	def test_downloader_private_parse(self):
		# seq3, a list of SeqRecords of which the third is the right
		#  size and contains no Nste
		res = self.downloader._parse(self.seq3)
		self.assertIsNotNone(res)

	def test_downloader_private_download(self):
		res = self.downloader._download(self.seqids)
		self.assertEqual(len(res), 3)

	def test_downloader_run(self):
		# reset thoroughness and deja_vues
		self.downloader.thoroughness = 1
		self.downloader.deja_vues = []
		res = self.downloader.run(self.taxids)
		self.assertEqual(len(res), 3)

	def test_findbestgenes(self):
		res = dtools.findBestGenes(self.namesdict, self.genedict, 3,\
		 self.allrankids, minnseq = 1, minpwithseq = 0.5, verbose = False)
		self.assertEqual(res[0], 'gene1')

if __name__ == '__main__':
    unittest.main()