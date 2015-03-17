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

# GLOBALS
working_dir = os.path.dirname(__file__)


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
    restart = False
    details = False
    retry = False
    reset = False
    stats = False

    def __init__(self):
        pass


class SetupTestSuite(unittest.TestCase):

    def setUp(self):
        self.true_getThreads = stools.getThreads
        stools.getThreads = dummyGetThreads

    def tearDown(self):
        stools.getThreads = self.true_getThreads
        setup_folders = ['folder_1', 'folder_2', 'tempfiles']
        while setup_folders:
            try:
                setup_folder = setup_folders.pop()
                shutil.rmtree(setup_folder)
            except OSError:
                pass

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
        _, _, _, _, _, _, stages = stools.parseArguments(args)
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
        self.assertTrue('folder_1' in folders)
        self.assertTrue('folder_2' in folders)
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
        arguments = {'genedict': {}, 'paradict': {}, 'terms': []}
        directory = '.'
        stools.prime(directory, arguments, 2)
        self.assertTrue(os.path.isfile(os.path.join('tempfiles',
                                                    'genedict.p')))
        self.assertTrue(os.path.isfile(os.path.join('tempfiles',
                                                    'paradict.p')))
        self.assertTrue(os.path.isfile(os.path.join('tempfiles',
                                                    'terms.p')))
        self.assertTrue(os.path.isfile(os.path.join('tempfiles',
                                                    'threads.p')))
        shutil.rmtree('tempfiles')
        os.remove('info.txt')

    def test_record_pars(self):
        paradict = {'a_parameter': 'a_value'}
        self.assertTrue(isinstance(stools.recordPars(paradict), str))

    def test_record_gpars(self):
        genedict = {'a_gene': {'a_parameter': 'a_value'}}
        self.assertTrue(isinstance(stools.recordGpars(genedict), str))

    def test_read_in_names(self):
        # write a list of names and read in
        with open('names.txt', 'w') as file:
            file.write('a_name\n')
        names = stools.readInNames('.')
        self.assertEqual(names[0], 'a_name')
        os.remove('names.txt')

    def test_read_in_gene_pars(self):
        # read in test parameters in /data
        default_file = os.path.join(working_dir, 'data',
                                    'test_default_gene_parameters.csv')
        test_file = os.path.join(working_dir, 'data',
                                 'test_gene_parameters.csv')
        res = stools.readInGenePars(test_file, default_file)
        # COI minlen should be 'a different value'
        self.assertEqual(res['COI']['minlen'], 'a different value')

    def test_read_in_pars(self):
        # read in test parameters in /data
        default_file = os.path.join(working_dir, 'data',
                                    'test_default_parameters.csv')
        test_file = os.path.join(working_dir, 'data', 'test_parameters.csv')
        res = stools.readInPars(test_file, default_file)
        # nseqs should be 'a different value'
        self.assertEqual(res['nseqs'], 'a different value')

    def test_sort_args(self):
        # make sure error is raised if no names.txt
        directory = '.'
        email = 'an.email'
        logger = stools.setUpLogging(verbose=False, debug=False,
                                     logname='testlogger')
        default_pars_file = os.path.join(working_dir, 'data',
                                         'test_default_parameters.csv')
        default_gpars_file = os.path.join(working_dir, 'data',
                                          'test_default_gene_parameters.csv')
        with self.assertRaises(stools.PrimingError):
            stools.sortArgs(directory, email, logger, default_pars_file,
                            default_gpars_file)
        # create a names.txt and test arguments returned
        names = ['name1', 'name2', 'name3', 'name4', 'name5', 'name6']
        with open('names.txt', 'w') as file:
            for name in names:
                file.write(name + '\n')
        res = stools.sortArgs(directory, email, logger, default_pars_file,
                              default_gpars_file)
        self.assertTrue(isinstance(res['terms'], list))
        self.assertTrue(isinstance(res['genedict'], dict))
        self.assertTrue(isinstance(res['paradict'], dict))
        os.remove('names.txt')
        os.remove('log.txt')

if __name__ == '__main__':
    unittest.main()
