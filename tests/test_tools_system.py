#! /bin/usr/env python
# D.J. Bennett
# 06/11/2014
"""
Tests for system tools.
"""

# PACKAGES
import unittest
import re
import os
import shutil
import pickle
import pglt.tools.system_tools as stools


# DUMMIES
class dummy_Logger(object):

    def __init__(self):
        pass

    def info(self, msg):
        pass

    def debug(self, msg):
        pass


# FUNCTIONS
def dummyNamesStageWError(wd, logger):
    raise(stools.TooFewSpeciesError)


def dummySortArgs(directory, email, logger, default_pars_file,
                  default_gpars_file):
    pass


def dummyPrime(directory, arguments, threads):
    progress = {'1': 'not run', '2': 'not run', '3': 'not run', '4': 'not run'}
    temp_dir = os.path.join(directory, 'tempfiles')
    if not os.path.isdir(temp_dir):
        os.mkdir(temp_dir)
    with open(os.path.join(temp_dir, 'progress.p'), "wb") as file:
        pickle.dump(progress, file)


class SystemTestSuite(unittest.TestCase):

    def setUp(self):
        self.logger = dummy_Logger()
        # replace names stage tuple
        self.true_names_stage = stools.Stager.STAGES['1']
        stools.Stager.STAGES['1'] = (dummyNamesStageWError, '1_names')
        # replace prime function
        self.true_prime = stools.prime
        self.true_sortArgs = stools.sortArgs
        stools.prime = dummyPrime
        stools.sortArgs = dummySortArgs

    def tearDown(self):
        stools.Stager.STAGES['1'] = self.true_names_stage
        stools.prime = self.true_prime
        stools.sortArgs = self.true_sortArgs
        system_files = []
        while system_files:
            try:
                system_file = system_files.pop()
                os.remove(system_file)
            except OSError:
                pass
        system_folders = ['1_names', '4_phylogeny', 'folder1', 'folder2',
                          'folder3', 'folder4', 'folder5', 'folder6',
                          'folder7', 'folder8', 'folder9', 'folder10',
                          'tempfiles']
        # system_folders = []
        while system_folders:
            try:
                system_folder = system_folders.pop()
                shutil.rmtree(system_folder)
            except OSError:
                pass

    def test_stager(self):
        dummyPrime(directory='.', arguments=None, threads=None)
        # init should raise error for '6'
        with self.assertRaises(stools.StageError):
            stools.Stager('.', '6')
        # should get TooFewSpeciesError
        stager = stools.Stager('.', '1')
        stager.run()
        # log.txt should exist and have the error
        with open(os.path.join('1_names', 'log.txt')) as file:
            lines = file.readlines()
        logtext = ''
        for line in lines:
            logtext += line
        self.assertTrue(re.search('ERROR', logtext))

    def test_runner(self):
        # run for 10 folders using dummy names stage
        folders = ['folder1', 'folder2', 'folder3', 'folder4', 'folder5',
                   'folder6', 'folder7', 'folder8', 'folder9', 'folder10',
                   'tempfiles']
        for folder in folders:
            os.mkdir(folder)
        stages = ['1']
        runner = stools.Runner(folders=folders, nworkers=10, stages=stages,
                               threads_per_worker=1, wd='.',
                               email='an.email', verbose=False, debug=False,
                               logger=self.logger, retry=False)
        runner.setup()
        runner.run()
        # check each folder has a 1_names
        for folder in folders:
            self.assertTrue(os.path.isdir(os.path.join(folder, '1_names')))

    def test_termination_pipe(self):
        # simple command, make sure it has been executed
        pipe = stools.TerminationPipe(cmd='mkdir folder1')
        pipe.run()
        self.assertTrue(os.path.isdir('folder1'))

if __name__ == '__main__':
    unittest.main()
