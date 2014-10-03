#! /bin/usr/env python
## D.J. Bennett
## 26/05/2014
"""
Test alignment stage.
"""

import unittest,pickle,os,shutil
from mpe.stages import alignment_stage
from Bio import AlignIO

## Dirs
working_dir = os.path.dirname(__file__)

## Dummies
class Dummy_SeqStore(object):
	def __init__(self,gene_dir,seq_files,minfails,\
		mingaps,minoverlap):
		pass
	def __len__(self):
		return 0

class Dummy_Aligner(object):
	def __init__(self, seqstore, mingaps, minoverlap, minseedsize,\
		maxseedsize, maxtrys, maxseedtrys, gene_type):
		pass
	def run(self):
		return alignment

## Test data
with open(os.path.join(working_dir, 'data','test_alignment_ref.faa'),\
	'r') as file:
	alignment = AlignIO.read(file, 'fasta')
genedict = {'rbcl':{'mingaps':0.1,'minoverlap':0.1,'maxtrys':100,\
'minseedsize':5,'maxseedsize':20,'maxseedtrys':10,'minfails':10,\
'type':'both'},'COI':{'mingaps':0.1,'minoverlap':0.1,'maxtrys':\
100,'minseedsize':5,'maxseedsize':20,'maxseedtrys':10,'minfails'\
:10,'type':'shallow'}}
paradict = {'naligns':1} # don't let it run more than once
# all the names in reference alignment
namesdict = {}
namesdict['Ignatius_tetrasporus'] = {"txids" : [1,2],\
	"unique_name" : 'Ignatius_tetrasporus',\
	"rank" : 'species','genes':2}
namesdict['Oltmannsiellopsis_viridis'] = {"txids" : [1,2],\
	"unique_name" : 'Oltmannsiellopsis_viridis',\
	"rank" : 'species','genes':2}
namesdict['Scotinosphaera_austriaca_H5304'] = {"txids" : [1,2],\
	"unique_name" : 'Scotinosphaera_austriaca_H5304',\
	"rank" : 'species','genes':2}
namesdict['Scotinosphaera_facciolae_H5309'] = {"txids" : [1,2],\
	"unique_name" : 'Scotinosphaera_facciolae_H5309',\
	"rank" : 'species','genes':2}
namesdict['Scotinosphaera_gibberosa_H5301'] = {"txids" : [1,2],\
	"unique_name" : 'Scotinosphaera_gibberosa_H5301',\
	"rank" : 'species','genes':2}
namesdict['Scotinosphaera_gibberosa_H5302'] = {"txids" : [1,2],\
	"unique_name" : 'Scotinosphaera_gibberosa_H5302',\
	"rank" : 'species','genes':2}
namesdict['Scotinosphaera_lemnae_H5303a'] = {"txids" : [1,2],\
	"unique_name" : 'Scotinosphaera_lemnae_H5303a',\
	"rank" : 'species','genes':2}
namesdict['Scotinosphaera_lemnae_H5303b'] = {"txids" : [1,2],\
	"unique_name" : 'Scotinosphaera_lemnae_H5303b',\
	"rank" : 'species','genes':2}
namesdict['Scotinosphaera_sp_H5305'] = {"txids" : [1,2],\
	"unique_name" : 'Scotinosphaera_sp_H5305',\
	"rank" : 'species','genes':2}
namesdict['Scotinosphaera_sp_H5306'] = {"txids" : [1,2],\
	"unique_name" : 'Scotinosphaera_sp_H5306',\
	"rank" : 'species','genes':2}
namesdict['Ulothrix_zonata'] = {"txids" : [1,2],\
	"unique_name" : 'Ignatius_tetrasporus',\
	"rank" : 'species','genes':2}

class AlignmentStageTestSuite(unittest.TestCase):

	def setUp(self):
		# stub functions and class
		self.True_SeqStore = alignment_stage.atools.SeqStore
		self.True_Aligner = alignment_stage.atools.Aligner
		alignment_stage.atools.SeqStore = Dummy_SeqStore
		alignment_stage.atools.Aligner = Dummy_Aligner
		# write out necessary files to run
		with open(".genedict.p", "wb") as file:
			pickle.dump(genedict, file)
		with open(".paradict.p", "wb") as file:
			pickle.dump(paradict, file)
		with open(".namesdict.p", "wb") as file:
			pickle.dump(namesdict, file)
		# download files so it can 'read' in sequences
		os.mkdir('2_download')
		os.mkdir(os.path.join('2_download','COI'))
		os.mkdir(os.path.join('2_download','rbcl'))

	def tearDown(self):
		# remove all files and folders potentially generated by alignment stage
		alignment_files = ['.paradict.p','.namesdict.p','.genedict.p']
		while alignment_files:
			try:
				alignment_file = alignment_files.pop()
				os.remove(alignment_file)
			except OSError:
				pass
		alignment_folders = ['2_download','3_alignment']
		while alignment_folders:
			try:
				alignment_folder = alignment_folders.pop()
				shutil.rmtree(alignment_folder)
			except OSError:
				pass
		alignment_stage.atools.SeqStore = self.True_SeqStore
		alignment_stage.atools.Aligner = self.True_Aligner

	def test_alignment_stage(self):
		# run
		res = alignment_stage.run()
		self.assertIsNone(res)

if __name__ == '__main__':
    unittest.main()