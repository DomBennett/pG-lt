#! /bin/usr/env python
## D.J. Bennett
## 26/05/2014
"""
Test download stage.
"""

import unittest,pickle,os
from mpe.stages import download_stage
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

## Dummies
def dummy_findBestGenes(namesdict, genedict, thoroughness,\
	allrankids, minnseq, target, minnspp):
	# return a list of genes
	return ['rbcl', 'COI']

class Dummy_Downloader(object):
	def __init__(self,gene_names, nseqs, thoroughness,\
				maxpn, seedsize, maxtrys, mingaps,\
				minoverlap, maxlen, minlen):
		pass

	def run(self, taxids):
		seq = 'A' * 500
		seq = SeqRecord(Seq(seq), id = 'testseq')
		return [seq]

## Test data
genedict = {'rbcl':{'names':['rbcl']},'COI':{'names':['COI']}}
paradict = {'email':'', 'nseqs':100,'download_thoroughness':\
3,'maxlen':2000,'minlen':300}
namesdict = {}
namesdict['query_name'] = {"txids" : [1,2],\
	"unique_name" : 'returned_name', "rank" : 'species'}
allrankids = []

class DownloadStageTestSuite(unittest.TestCase):

	def setUp(self):
		# stub functions and class
		self.True_Downloader = download_stage.dtools.Downloader
		self.true_findBestGenes = download_stage.dtools.findBestGenes
		download_stage.dtools.Downloader = Dummy_Downloader
		download_stage.dtools.findBestGenes = dummy_findBestGenes
		# write out necessary files to run
		with open(".genedict.p", "wb") as file:
			pickle.dump(genedict, file)
		with open(".paradict.p", "wb") as file:
			pickle.dump(paradict, file)
		with open(".namesdict.p", "wb") as file:
			pickle.dump(namesdict, file)
		with open(".allrankids.p", "wb") as file:
			pickle.dump(allrankids, file)

	def tearDown(self):
		download_stage.dtools.Downloader = self.True_Downloader
		download_stage.dtools.findBestGenes = self.true_findBestGenes

	def test_download_stage(self):
		# run
		res = download_stage.run()
		# remove files and folders
		os.remove('.paradict.p')
		os.remove('.namesdict.p')
		os.remove('.genedict.p')
		os.remove('.allrankids.p')
		os.remove(os.path.join('2_download', 'COI', 'query_name.fasta'))
		os.remove(os.path.join('2_download', 'rbcl', 'query_name.fasta'))
		os.rmdir(os.path.join('2_download', 'COI'))
		os.rmdir(os.path.join('2_download', 'rbcl'))
		os.rmdir('2_download')
		self.assertIsNone(res)

if __name__ == '__main__':
    unittest.main()
