#! /bin/usr/env python
# D.J. Bennett
# 06/11/2014
"""
Tests for special tools.
"""

# PACKAGES
import os
import unittest
import shutil
import time
import pickle
import pglt.tools.special_tools as stools

working_dir = os.path.dirname(__file__)
test_folder = os.path.join(working_dir, 'data', 'test_reset')


# FUNCTIONS
def sleep(seconds):
    time.sleep(seconds)


class SpecialTestSuite(unittest.TestCase):

    def setUp(self):
        self.wd = os.path.join(working_dir, 'test_reset')
        shutil.copytree(test_folder, self.wd)
        self.reseter = stools.Reseter(wd=self.wd)
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
        self.assertEqual(progress['2'], 'not run')
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

    def test_timeit(self):
        # time time.sleep for 0.1 second
        res, s = stools.timeit(sleep, seconds=0.1)
        self.assertAlmostEqual(s, 0.1, places=1)

    def test_get_threads_local_false(self):
        res = stools.getThreads()
        self.assertTrue(isinstance(res, int))
        self.assertTrue(res > 0)

    def test_get_threads_local_true(self):
        # create .threads.p
        with open("threads.p", "wb") as file:
            pickle.dump(2, file)
        res = stools.getThreads(wd=os.getcwd())
        self.assertEqual(res, 2)
        os.remove("threads.p")

    def test_clean(self):
        # create dummy folders and files
        log_file = os.path.join('afolder', 'log.txt')
        names_folder = os.path.join('afolder', '1_names')
        os.mkdir('afolder')
        os.mkdir(names_folder)
        with open(log_file, "wb") as file:
            file.write('words words words ...\n')
        # clean
        stools.clean()
        # make sure they don't exist
        self.assertFalse(os.path.isfile(log_file))
        self.assertFalse(os.path.isdir(names_folder))

    def tearDown(self):
        shutil.rmtree(self.wd)
        if os.path.isdir('afolder'):
            shutil.rmtree('afolder')

if __name__ == '__main__':
    unittest.main()
