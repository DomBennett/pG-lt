#! /bin/usr/env python
# D.J. Bennett
# 26/05/2014
"""
Tests for download tools.
"""

import unittest
import pickle
import logging
import os
import pglt.tools.download_tools as dtools

# DIRS
working_dir = os.path.dirname(__file__)

# DUMMIES

# expected terms and term search results
t1_term = '(txid1[PORGN]) OR txid2[PORGN] AND (("name1"[GENE]) \
OR "name2"[GENE]) NOT predicted[TI] NOT genome[TI] NOT unverified[TI]'
t2_term = '(txid1[PORGN]) OR txid2[PORGN] AND (("name1"[TI]) \
OR "name2"[TI]) NOT predicted[TI] NOT genome[TI] NOT unverified[TI]'
t3_term = '(txid1[PORGN]) OR txid2[PORGN] AND (("name1") \
OR "name2") NOT predicted[TI] NOT shotgun[TI] NOT scaffold[TI] \
NOT assembly[TI] NOT unverified[TI]'
t1_search_res = {'Count': 0}
t2_search_res = {'Count': 2, 'IdList': ['seq1', 'seq2']}
t3_search_res = {'Count': 3, 'IdList': ['seq1', 'seq2', 'seq3']}
outgroup_res = {'Count': 3, 'IdList': ['seq4', 'seq5', 'seq6']}

# Example seqrecord for findgeneinseq
with open(os.path.join(working_dir, 'data', "test_findgeneinseq_examplesequence\
.p"), "rb") as file:
    sequence = pickle.load(file)


# Dummy seq records for download
class dummy_Seq(object):

    def __init__(self):
        pass

    def __str__(self):
        # Just for parsing
        return "A" * 500


class dummy_SeqRecord(object):

    def __init__(self, description, length=500):
        self.description = description
        self.length = length
        self.seq = dummy_Seq()
        self.features = None

    def __len__(self):
        return self.length

seq1 = dummy_SeqRecord(description="A sequence of NAME1")
seq2 = dummy_SeqRecord(description="A sequence of NAME2")
seq3 = [dummy_SeqRecord(description="A sequence of NAME3"),
        dummy_SeqRecord(description="A sequence of NAME4"),
        dummy_SeqRecord(description="A sequence of NAME5"),
        dummy_SeqRecord(description="A sequence of NAME1")]


# Sequences -- just Ts and Fs -- for testing filter
sequences = [True for i in range(80)]
sequences.extend([False for i in range(20)])


# Dependent stubs
def dummy_eSearch(term, logger, retStart=0, retMax=1, usehistory="n",
                  db="nucleotide"):
    if term == t1_term:
        return t1_search_res
    if term == t2_term:
        return t2_search_res
    if term == t3_term:
        return t3_search_res
    else:
        return outgroup_res


def dummy_eFetch(ncbi_id, logger, db="nucleotide"):
    if ncbi_id == 'seq1':
        return seq1
    elif ncbi_id == 'seq2':
        return seq2
    elif ncbi_id == 'seq3':
        return seq3
    else:
        # return all as list
        return [seq1, seq2, seq3]


def dummy_blast(query, subj, minoverlap, logger, wd, threads):
    # should return bools and positions
    # pretend they've matched from 0-100 base positions
    return query, [0, 100]


def dummy_checkAlignment(alignment, mingaps, minoverlap, minlen, logger):
    return alignment

# downloader init variables
taxids = ['1', '2']
gene_names = ['name1', 'name2']
nseqs = 2
thoroughness = 3
maxpn = 0.1
votesize = 3
mingaps = 0.01
minoverlap = 200
maxtrys = 100
maxlen = 2000
minlen = 300

# dictionary variables
namesdict = {"species1": {'txids': ['1', '2']}, 'outgroup': {'txids': ['4']}}
allrankids = [1, 2, 3]
genedict = {'gene1': {'taxid': '3', 'names': ['name1', 'name2'],
                      'type': 'deep'}}


class DownloadTestSuite(unittest.TestCase):

    def setUp(self):
        self.logger = logging.getLogger()
        self.wd = os.getcwd()
        self.true_eSearch = dtools.etools.eSearch
        self.true_eFetch = dtools.etools.eFetch
        self.true_blast = dtools.atools.blast
        self.true_checkAlignment = dtools.atools.checkAlignment
        dtools.etools.eSearch = dummy_eSearch
        dtools.etools.eFetch = dummy_eFetch
        dtools.atools.blast = dummy_blast
        dtools.atools.checkAlignment = dummy_checkAlignment
        # mock Downloader instance
        self.downloader = dtools.Downloader(gene_names=gene_names,
                                            nseqs=nseqs,
                                            thoroughness=thoroughness,
                                            maxpn=maxpn, votesize=votesize,
                                            maxtrys=maxtrys,
                                            minoverlap=minoverlap,
                                            maxlen=maxlen, minlen=minlen,
                                            logger=self.logger, wd=self.wd)
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
        # repatch
        dtools.etools.eSearch = self.true_eSearch
        dtools.etools.eFetch = self.true_eFetch
        dtools.atools.blast = self.true_blast
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
        # expect to only find 2, 1 and 0 sequences
        # it should search until it finds two sequences (nseqs = 2),
        #  and then on the next search after raising its thoroughness
        #  it should find the last sequence. Searching again will
        #  find no more.
        res1 = self.downloader._search(self.taxids)
        res2 = self.downloader._search(self.taxids)
        res3 = self.downloader._search(self.taxids)
        self.assertEqual([len(res1), len(res2), len(res3)], [2, 1, 0])

    def test_downloader_private_filter(self):
        # weeds out Falses from sequences, should not be any Falses
        # self.sequences is 80 Ts and 20 Fs
        sequences = self.sequences[:]
        res_filtered, res_downloaded = self.downloader._filter(sequences)
        self.assertEqual(len(res_filtered), 80)
        self.assertEqual(len(res_downloaded), 20)

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
        #  size and contains no Ns
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
        res = dtools.findBestGenes(self.namesdict, self.genedict, 3,
                                   self.allrankids, logger=self.logger,
                                   minnseq=1, target=1, minnspp=0)
        self.assertEqual(res[0], 'gene1')

    def test_get_clusters(self):
        # make a gene_sequences: [(name, sequence), ...]
        # should return 80 sequences
        names = ['sp1', 'sp2', 'sp3', 'sp4', 'sp5', 'sp6', 'sp7', 'sp8', 'sp9',
                 'sp10']*10
        gene_sequences = zip(names, self.sequences)
        res = dtools.getClusters(gene_sequences, 0.5, self.logger, self.wd)
        self.assertEqual(len(res[0]), 80)


if __name__ == '__main__':
    unittest.main()
