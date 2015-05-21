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


# FUNCTIONS
def sleep(seconds):
    time.sleep(seconds)


class SpecialTestSuite(unittest.TestCase):

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
        if os.path.isdir('afolder'):
            shutil.rmtree('afolder')

if __name__ == '__main__':
    unittest.main()
