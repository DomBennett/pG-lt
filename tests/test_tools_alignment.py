#! /bin/usr/env python
# D.J. Bennett
# 26/05/2014
"""
Tests for alignment tools.
"""

import unittest
import os
import re
import copy
import pickle
import random
import pglt.tools.alignment_tools as atools
from pglt import _MAFFT as mafft
from pglt import _MAFFTQ as mafftq
from pglt import _MAFFTX as mafftx
from Bio import SeqIO
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# DIRS
working_dir = os.path.dirname(__file__)

# TEST DATA
test_seqs = []
with open(os.path.join(working_dir, 'data', 'test_sequences.faa'),
          "rU") as infile:
    for record in SeqIO.parse(infile, "fasta"):
        test_seqs.append(record)

with open(os.path.join(working_dir, 'data', 'test_sequences.faa'),
          "rU") as infile:
        test_alignment = AlignIO.read(infile, "fasta")

with open(os.path.join(working_dir, "data", "test_alignment.p"),
          "rb") as file:
    real_alignment = pickle.load(file)


# MOCK
# dummy logger needed because can't take copy of an object with a real logger
class dummy_Logger(object):

    def __init__(self):
        pass

    def info(self, msg):
        pass

    def debug(self, msg):
        pass


def dummy_blast(query, subj, minoverlap, logger, wd, threads):
    # should return bools and positions
    bools = [True for e in query]
    positions = [0 for e in subj]
    max_positions = [len(e) for e in subj]
    positions.extend(max_positions)
    return bools, positions


def dummy_align(command, sequences, timeout, logger, wd, threads):
    return test_alignment


def dummy_add(alignment, sequence, timeout, logger, wd, threads):
    return test_alignment


def dummy_check_alignment(alignment, mingaps, minoverlap, minlen, logger):
    return True


class AlignmentTestSuite(unittest.TestCase):

    def setUp(self):
        self.wd = os.getcwd()
        self.logger = dummy_Logger()
        self.true_align = atools.align
        self.true_add = atools.add
        self.true_check_alignment = atools.checkAlignment
        self.true_blast = atools.blast
        atools.blast = dummy_blast
        genedir = os.path.join(working_dir, 'data', 'test_sequences')
        seqfiles = sorted(os.listdir(genedir))
        seqfiles = [e for e in seqfiles if not re.search("^\.|^log\.txt$", e)]
        self.store = atools.SeqStore(genedir=genedir, seqfiles=seqfiles,
                                     minfails=10, mingaps=0.5, minoverlap=50,
                                     logger=self.logger)
        self.aligner = atools.Aligner(self.store, mingaps=0.5, minoverlap=50,
                                      minseedsize=3, maxseedsize=20,
                                      maxtrys=10, maxseedtrys=10,
                                      gene_type='shallow', outgroup=False,
                                      logger=self.logger)

    def tearDown(self):
        atools.blast = self.true_blast
        atools.align = self.true_align
        atools.add = self.true_add
        atools.checkAlignment = self.true_check_alignment
        del self.logger

    def test_gennonalignment(self):
        alignment = atools.genNonAlignment(1, 100)
        self.assertEqual(len(alignment[0]), 100)

    @unittest.skipIf(not all([atools.mafft, atools.mafftq, atools.mafftx]),
                     "Requires MAFFT, MAFFT-QINSI and MAFFT-XINSI")
    def test_version(self):
        # try different combinations of sequences
        def genSequences(n, length):
            s = 'A' * length
            return [s for i in range(n)]
        res = atools.version(genSequences(30, 800), 'deep')
        self.assertEqual(res, mafftx)
        res = atools.version(genSequences(90, 800), 'deep')
        self.assertEqual(res, mafftq)
        res = atools.version(genSequences(90, 2000), 'deep')
        self.assertEqual(res, mafft + ' --auto')

    @unittest.skipIf(not atools.mafft, "Requires MAFFT")
    def test_align(self):
        # align and check if results exist
        res = atools.align(command='mafft --auto', sequences=test_seqs,
                           timeout=99999, logger=self.logger, threads=2,
                           wd=self.wd)
        self.assertTrue(res.__class__.__name__ == 'MultipleSeqAlignment')

    @unittest.skipIf(not atools.mafft, "Requires MAFFT")
    def test_add(self):
        # add to test_alignment and check if result
        #  exists
        res = atools.add(alignment=test_alignment, sequence=test_seqs[0],
                         timeout=99999, logger=self.logger, wd=self.wd,
                         threads=2)
        self.assertTrue(res.__class__.__name__ == 'MultipleSeqAlignment')

    @unittest.skipIf(not atools.blastn, "Requires BLASTN")
    def test_blast(self):
        # real_alignment and real_sequences are identical
        real_seqs = [e for e in real_alignment]
        # choose sample of real seqs for speed
        queryseqs = random.sample(real_seqs, 5)
        subjseq = random.sample(real_seqs, 1)
        res, _ = self.true_blast(query=queryseqs, subj=subjseq, minoverlap=50,
                                 logger=self.logger, wd=self.wd, threads=2)
        # all true
        self.assertTrue(all(res))

    def test_checkalignment_arg_mingaps(self):
        # check mingaps argument (proportion of internal gaps)
        # check with good alignment
        res = atools.checkAlignment(test_alignment, mingaps=0.5, minoverlap=1,
                                    minlen=1, logger=self.logger)
        self.assertTrue(res)
        # check with bad alignment
        # this sequence has a large internal gap covering
        #  50% of the alignment
        bad_seq = '-A' * 50
        bad_seq = SeqRecord(Seq(bad_seq), id='bad')
        bad_alignment = test_alignment[:]
        bad_alignment.append(bad_seq)
        res = atools.checkAlignment(bad_alignment, mingaps=0.1, minoverlap=50,
                                    minlen=1, logger=self.logger)
        self.assertFalse(res)

    def test_checkalignment_arg_minoverlap(self):
        # check minoverlap argument (minimum number of
        #    overlapping nucleotides)
        # bad seq has a large external gap -- the first 51
        #  nucleotides are missing
        bad_seq = '-' * 51 + 'A' * 49
        bad_seq = SeqRecord(Seq(bad_seq), id='bad')
        bad_alignment = test_alignment[:]
        bad_alignment.append(bad_seq)
        res = atools.checkAlignment(bad_alignment, mingaps=0.5, minoverlap=50,
                                    minlen=1, logger=self.logger)
        self.assertFalse(res)

    def test_checkalignment_arg_minlen(self):
        # check minlen argument
        res = atools.checkAlignment(test_alignment, mingaps=0.5, minoverlap=1,
                                    minlen=101, logger=self.logger)
        self.assertFalse(res)

    @unittest.skipIf(not atools.blastn, "Requires BLASTN")
    def test_seqstore_private_alignmentblast(self):
        # use real blast
        atools.blast = self.true_blast
        real_seqs = [e for e in real_alignment]
        next_seqs = random.sample(real_seqs, 5)
        alignment = random.sample(real_seqs, 1)
        # same as blast test but with a dodgy sequence
        bad_seq = 'T' * 100
        bad_seq = SeqRecord(Seq(bad_seq), id='bad')
        next_seqs_w_bad = next_seqs[:]
        next_seqs_w_bad.append(bad_seq)
        for i in range(len(next_seqs_w_bad)):
            res = self.store._alignmentBlast(query=[next_seqs_w_bad[i]],
                                             sequences_in_alignment=alignment)
            if i < 5:
                self.assertTrue(res)
            else:
                # bad seq should not be returned
                self.assertFalse(res)
        corrected_seq = bad_seq + random.sample(real_seqs, 1)[0]
        next_seqs_w_corrected = next_seqs[:]
        next_seqs_w_corrected.append(corrected_seq)
        for i in range(len(next_seqs_w_corrected)):
            res = self.store._alignmentBlast(query=next_seqs_w_corrected,
                                             sequences_in_alignment=alignment)
            # corrected seq should now be returned
            self.assertTrue(res)
        # and its returned sequence shouldn't have the bad sequence
        self.assertTrue(str(res[1].seq) != str(corrected_seq.seq))
        # switch back to dummy blast
        atools.blast = dummy_blast

    def test_seqstore_private_add(self):
        store = copy.deepcopy(self.store)
        # add lists to obj for add to work
        store.sppool = store.keys()
        store.sequences_in_alignment = []
        store._add()
        res = store.sequences_in_alignment[0]
        # the species should no longer be in the pool
        self.assertFalse(res[0].id in store.sppool)

    def test_seqstore_start(self):
        store = copy.deepcopy(self.store)
        seqs = store.start(3)
        # none of the species should be in the pool
        res = [e.id in store.sppool for e in seqs]
        self.assertFalse(all(res))
        # each sequence returned should be from different species
        res = list(set([e.id for e in seqs]))
        self.assertEqual(len(res), 3)

    def test_seqstore_back(self):
        store = copy.deepcopy(self.store)
        before_seqs = store.start(5)
        after_seq = store.back(before_seqs[:-1])
        # after_seq sp should not be in before_seqs
        res = [after_seq.id == e.id for e in before_seqs]
        self.assertFalse(all(res))
        # last before_seq should have penalty
        res = store[before_seqs[-1].id][0]
        res = [e[1] for e in res]
        self.assertEqual(sum(res), 1)

    def test_seqstore_next(self):
        store = copy.deepcopy(self.store)
        before_seqs = store.start(5)
        # add a mock penalty to seqs in alignment
        seqs = store.sequences_in_alignment
        for each in seqs:
            each[1] += 1
        after_seq = store.next(test_alignment)
        # after_seq sp should not be in before_seqs
        res = [after_seq.id == e.id for e in before_seqs]
        self.assertFalse(all(res))
        # after_seq sp should not be in sppool
        res = [after_seq.id == e for e in store.sppool]
        self.assertFalse(all(res))
        # no penalties in seqs in alignment
        seqs = store.sequences_in_alignment
        res = [e[1] == 0 for e in seqs]
        self.assertTrue(all(res))

    def test_seqstore_private_check_dropsequence(self):
        # add maxfails to a sequence
        store = copy.deepcopy(self.store)
        seq_to_drop = store['sp1'][0][0][0].name
        store['sp1'][0][0][1] = 11
        store._check()
        res = [e[0].name for e in store['sp1'][0]]
        # ._check() should have removed it
        self.assertFalse(seq_to_drop in res)

    def test_seqstore_private_check_dropspecies(self):
        # add maxfails to all seqs of one species
        store = copy.deepcopy(self.store)
        # add a sppool, so that dropped species can
        #  removed
        store.sppool = store.keys()
        sp_to_drop = 'sp1'
        for seq in store[sp_to_drop][0]:
            seq[1] = 11
        store._check()
        # ._check() should have removed it
        self.assertFalse(sp_to_drop in store.keys())

    def test_aligner_private_return(self):
        # create different types of alignments and see what
        #  is returned from _return
        outgroup_seq = SeqRecord(Seq('T' * 3141), id='outgroup')
        real_alignment_w_out = real_alignment[:]
        real_alignment_w_out.append(outgroup_seq)
        outgroup_seq = SeqRecord(Seq('T' * 100), id='outgroup')
        test_alignment_w_out = test_alignment[:]
        test_alignment_w_out.append(outgroup_seq)
        mock_multi_alignment = [real_alignment, test_alignment,
                                real_alignment_w_out, real_alignment[:5],
                                test_alignment_w_out]
        self.aligner.store = mock_multi_alignment
        alignment = self.aligner._return()
        # should return real_alignment_w_out (has 12 records)
        # it has an outgroup and it is the longest
        self.assertEqual(len(alignment), 12)

    def test_aligner_private_calctimeout(self):
        # dummy of 10 sequences of length 100
        dummy_alignment = [[1 for b in range(100)] for s in range(10)]
        self.aligner._calcTimeout(seconds=10, alignment=dummy_alignment)
        self.aligner._calcTimeout(seconds=10, alignment=dummy_alignment,
                                  align=False)
        # both should be 0.1 seconds -- the average number of seconds per b*10
        # 1000 bases takes 10 seconds. So per base it is 0.01.
        # Multiplying by buffer gives 0.1 seconds
        self.assertEqual(self.aligner.tadd, 0.1)
        self.assertEqual(self.aligner.talign, 0.1)

    def test_aligner_private_gettimeout(self):
        # dummy of 10 sequences of length 100
        dummy_sequences = [[1 for b in range(100)] for s in range(10)]
        dummy_sequence = [1 for b in range(100)]
        self.aligner.tadd = 0.1
        self.aligner.talign = 0.1
        tadd = self.aligner._getTimeout(sequences=dummy_sequences,
                                        sequence=dummy_sequence)
        talign = self.aligner._getTimeout(sequences=dummy_sequences)
        # number of nucleotides by tadd and talign
        self.assertEqual(tadd, 1100*0.1)
        self.assertEqual(talign, 1000*0.1)

    def test_aligner_private_calcseedsize(self):
        seedsize = self.aligner.seedsize
        # fail buffer times in a row, seedsize should drop
        for i in range(self.aligner.buffer):
            # will return 0s
            self.assertFalse(self.aligner._calcSeedsize(False))
        # now seedsize will have dropped
        self.assertEqual(self.aligner.seedsize, seedsize-1)
        # success buffer times in a row, seedsize should increase
        for i in range(self.aligner.buffer):
            # will return 0s
            self.assertFalse(self.aligner._calcSeedsize(True))
        # after buffer times, seedsize should have increased
        self.assertEqual(self.aligner.seedsize, seedsize)

    def test_aligner_private_seed(self):
        # switch to dummies
        atools.align = dummy_align
        atools.checkAlignment = dummy_check_alignment
        self.aligner.store = []
        # run seed which will be a success and add test_alignment to store
        success, trys = self.aligner._seed(0)
        self.assertTrue(success)
        self.assertFalse(trys)
        self.assertEqual(self.aligner.store[0], test_alignment)

    def test_aligner_private_add(self):
        # switch to dummies
        atools.add = dummy_add
        atools.checkAlignment = dummy_check_alignment
        # set up so that SeqStore methods work
        self.aligner.seqstore.sppool = self.aligner.seqstore.keys()
        self.aligner.seqstore.sequences_in_alignment = []
        self.aligner.store = [test_alignment[-1]]
        # run add which will not be finished
        finished, trys = self.aligner._add(0)
        self.assertFalse(finished)
        self.assertFalse(trys)
        # last store entry will be test_alignment
        self.assertEqual(self.aligner.store[-1], test_alignment)

    def test_aligner_run(self):
        # all the outgroup seqs should not align, alignment
        #  cannot have outgroup species
        res = self.aligner.run()
        spp = [e.id for e in res]
        self.assertFalse('outgroup' in spp)

if __name__ == '__main__':
    unittest.main()
