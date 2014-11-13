#! /bin/usr/env python
# D.J. Bennett
# 06/11/2014
"""
Tests for setup tools.
"""

# PACKAGES
import unittest
import os
import shutil
import logging
import pglt.tools.setup_tools as stools


# DUMMIES
def dummyGetThreads():
    return 100  # pretend the machine has 100 cores


class DummyArgs(object):
    stages = '1-4'
    clean = False
    email = 'an.email.address'
    threads = 10
    verbose = False
    debug = True

    def __init__(self):
        pass


class SetupTestSuite(unittest.TestCase):

    def setUp(self):
        self.true_getThreads = stools.getThreads
        stools.getThreads = dummyGetThreads

    def tearDown(self):
        stools.getThreads = self.true_getThreads
        if os.path.isdir('folder_1'):
            shutil.rmtree('folder_1')
        if os.path.isdir('folder_2'):
            shutil.rmtree('folder_2')
        if os.path.isfile('log.txt'):
            os.remove('log.txt')

    def test_calc_workers(self):
        # for 100 folder:
        #  -- 50 workers
        #  -- 2 threads_per_worker
        #  -- 0 spare
        res = stools.calcWorkers(threads=100, nfolders=100)
        self.assertEqual(res, (50, 2, 0))
        # for 10 folder:
        #  -- 10 workers
        #  -- 10 threads_per_worker
        #  -- 0 spare
        res = stools.calcWorkers(threads=100, nfolders=10)
        self.assertEqual(res, (10, 10, 0))
        # for 9 folder:
        #  -- 9 workers
        #  -- 11 threads_per_worker
        #  -- 1 spare
        res = stools.calcWorkers(threads=100, nfolders=9)
        self.assertEqual(res, (9, 11, 1))

    def test_create_parser(self):
        # http://dustinrcollins.com/testing-python-command-line-apps
        args = stools.createParser().parse_args(['-e', 'an.email.address'])
        self.assertEqual(args.stages, '1-4')
        self.assertFalse(args.verbose)
        self.assertFalse(args.debug)

    def test_parse_arguments(self):
        # test the stages are returned correctly
        args = DummyArgs()
        _, _, _, _, stages = stools.parseArguments(args)
        self.assertEqual(stages, ['1', '2', '3', '4'])
        # test bad args
        with self.assertRaises(SystemExit):
            args = DummyArgs()
            args.email = None
            stools.parseArguments(args)
        with self.assertRaises(SystemExit):
            args = DummyArgs()
            args.threads = -10
            stools.parseArguments(args)
        with self.assertRaises(SystemExit):
            args = DummyArgs()
            args.stages = '1-8'
            stools.parseArguments(args)
        with self.assertRaises(SystemExit):
            args = DummyArgs()
            args.stages = '4-1'
            stools.parseArguments(args)

    def test_get_folders(self):
        # make some folders with names.txts
        os.mkdir('folder_1')
        open(os.path.join('folder_1', 'names.txt'), 'w').close()
        os.mkdir('folder_2')
        open(os.path.join('folder_2', 'names.txt'), 'w').close()
        folders = stools.getFolders()
        self.assertEqual(folders, ['folder_1', 'folder_2'])
        shutil.rmtree('folder_1')
        shutil.rmtree('folder_2')

    def test_setup_logging(self):
        logger = stools.setUpLogging(verbose=True, debug=True,
                                     logname='a_logger')
        # more than two handlers: console + log.txt
        self.assertEqual(len(logger.handlers), 2)
        # debug level
        self.assertTrue(logger.level == 10)
        # a log file has been created
        self.assertTrue(os.path.isfile('log.txt'))
        os.remove('log.txt')
        del logger

    def test_teardown_logging(self):
        # create simple logger
        logger = logging.getLogger('a_logger')
        logger.addHandler(logging.StreamHandler())
        stools.tearDownLogging('a_logger')
        self.assertEqual(len(logger.handlers), 0)
        del logger

    def test_log_message(self):
        # not much to test here ...
        logger = logging.getLogger('a_logger')
        res = stools.logMessage(phase='program-start', logger=logger,
                                folders=['a'])
        self.assertIsNone(res)

    def test_prime(self):
        pass

    def test_sort_args(self):
        pass

    def test_record_pars(self):
        pass

    def test_record_gpars(self):
        pass

    def test_read_in_names(self):
        pass

    def test_read_in_gene_pars(self):
        pass

    def test_read_in_pars(self):
        pass


if __name__ == '__main__':
    unittest.main()
