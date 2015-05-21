#! /bin/usr/env python
# D.J. Bennett
# 06/11/2014
"""
Tests for reseter tools.
"""

# PACKAGES
import os
import unittest
import shutil
import pglt.tools.reseter_tools as rtools
from pglt.tools.setup_tools import readInPars
from pglt.tools.setup_tools import readInGenePars

working_dir = os.path.dirname(__file__)
test_folder = os.path.join(working_dir, 'data', 'test_reset')


class ReseterTestSuite(unittest.TestCase):

    def setUp(self):
        self.wd = os.path.join(working_dir, 'test_reset')
        shutil.copytree(test_folder, self.wd)
        self.paradict = readInPars('')
        self.genedict = readInGenePars('')
        self.reseter = rtools.Reseter(paradict=self.paradict,
                                      genedict=self.genedict,
                                      wd=self.wd)
        self.reseter.verbose = False

    def test_reseter_private_deletestagefolder(self):
        study_folder = os.path.join(self.wd, 'names_0')
        self.reseter._deleteStageFolder(study_folder=study_folder,
                                        stage_folder='1_names')
        self.assertFalse(os.path.isdir(os.path.join(study_folder, '1_names')))

    def test_reseter_private_readPickledFile(self):
        study_folder = os.path.join(self.wd, 'names_0')
        progress = self.reseter._readPickledFile(folder=study_folder,
                                                 filename='progress.p')
        self.assertTrue(isinstance(progress, dict))

    def test_reseter_private_writePickledFile(self):
        study_folder = os.path.join(self.wd, 'names_0')
        test = {'this': 'is a test'}
        self.reseter._writePickledFile(folder=study_folder, filename='test.p',
                                       pickled=test)
        del test
        test = self.reseter._readPickledFile(folder=study_folder,
                                             filename='test.p')
        self.assertEqual(test['this'], 'is a test')

    def test_reseter_private_resetstage(self):
        study_folder = os.path.join(self.wd, 'names_0')
        # CAREFUL: make sure real stage!
        self.reseter._resetstage(stage='2')
        self.assertFalse(os.path.isdir(os.path.join(study_folder,
                                                    '4_phylogeny')))
        progress = self.reseter._readPickledFile(folder=study_folder,
                                                 filename='progress.p')
        self.assertEqual(progress['1'], 'success')
        self.assertEqual(progress['2'], 'success')
        self.assertEqual(progress['3'], 'not run')
        self.assertEqual(progress['4'], 'not run')

    def test_reseter_private_resetparameters(self):
        study_folder = os.path.join(self.wd, 'names_0')
        # CAREFUL: make sure real key and value!
        self.reseter._resetparameters(key='nseqs', value=1)
        paradict = self.reseter._readPickledFile(folder=study_folder,
                                                 filename='paradict.p')
        self.assertEqual(paradict['nseqs'], 1)
        # CAREFUL: make sure real key and value!
        self.reseter._resetparameters(key='threads', value=10)
        threads = self.reseter._readPickledFile(folder=study_folder,
                                                filename='threads.p')
        self.assertEqual(threads, 10)

    def test_reseter_private_resetgeneparameters(self):
        study_folder = os.path.join(self.wd, 'names_0')
        # CAREFUL: make sure real gene, key and value!
        self.reseter._resetgeneparameters(gene='COI', key='minlen', value=1)
        genedict = self.reseter._readPickledFile(folder=study_folder,
                                                 filename='genedict.p')
        self.assertEqual(genedict['COI']['minlen'], 1)
        self.assertNotEqual(genedict['COII']['minlen'], 1)

    def test_reseter_private_setfolders(self):
        user_folders = '0_names 1_names'
        self.reseter._setfolders(user_folders)
        res = [os.path.basename(e) for e in self.reseter.folders]
        res = ' '.join(res)
        self.assertEqual(res, user_folders)

    def test_reseter_private_restore(self):
        # repalce paradict and genedict with nonsense
        study_folder = os.path.join(self.wd, 'names_0')
        paradict = {'this': 'is a test'}
        genedict = {'this': 'is a test'}
        self.reseter._writePickledFile(folder=study_folder,
                                       filename='paradict.p',
                                       pickled=paradict)
        self.reseter._writePickledFile(folder=study_folder,
                                       filename='genedict.p',
                                       pickled=genedict)
        self.reseter._restore(for_genedict='y', for_paradict='y', check=False)
        paradict = self.reseter._readPickledFile(folder=study_folder,
                                                 filename='paradict.p')
        genedict = self.reseter._readPickledFile(folder=study_folder,
                                                 filename='genedict.p')
        # test that the nonsense has been replaced
        with self.assertRaises(KeyError):
            paradict['this']
        with self.assertRaises(KeyError):
            genedict['this']
        key = self.paradict.keys()[0]
        self.assertEqual(paradict[key], self.paradict[key])
        gene = self.genedict.keys()[0]
        key = self.genedict[gene].keys()[0]
        self.assertEqual(genedict[gene][key], self.genedict[gene][key])

    def tearDown(self):
        shutil.rmtree(self.wd)

if __name__ == '__main__':
    unittest.main()
