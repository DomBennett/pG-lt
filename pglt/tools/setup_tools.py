#! /bin/usr/env python
# D.J. Bennett
# 07/11/2014
"""
pglt setup tools
"""

# PACKAGES
import argparse
import sys
import os
import re
import pickle
import csv
import logging
import platform
from datetime import datetime
from special_tools import clean
from special_tools import getThreads

# MESSAGES
description = """
----------------------------------------------------------------------
pG-lt version 1, Copyright (C) 2014 Bennett
----------------------------------------------------------------------
This program comes with ABSOLUTELY NO WARRANTY. This is free software,
and you are welcome to redistribute it under certain conditions.
For more details, type `run_pglt.py --details`.
----------------------------------------------------------------------
"""
nonamestxt_msg = '\nERROR: No folders containing \'names.txt\' files \
found! All taxonomic names should be placed in subdirectories and \
called: \'names.txt\''
ioerror_msg = "[{0}] file could not be opened in [{1}]. Check that \
it is not opened by another program"
priming_msg = '\nERROR: The program was unable to start due to a \
problem with the files and folders in the study directory. Check the \
parameters and gene parameters .csv for any potential conflicts.'

# PROGRESS DICT
progress = {'1': 'not run', '2': 'not run', '3': 'not run', '4': 'not run'}


# ERROR CLASSES
class PrimingError(Exception):
    pass


# FUNCTIONS
def printHeader():
    """Print a nice program description header"""
    print description


def calcWorkers(threads, nfolders, min_threads_per_worker=2,
                max_threads_per_worker=100):
    """Calculate the number of workers for parallel running of folders"""
    # get available threads on machine
    available_threads = getThreads()
    # make sure threads arg is not greater than those available
    if threads > available_threads:
        sys.exit('More threads specified than avaiable on machine')
    if threads == -1:
        threads = available_threads
    # calc min_threads_per_worker if it is greater than threads
    if min_threads_per_worker > threads:
        min_threads_per_worker = threads
    # calc max_threads_per_worker if it is greater than threads
    if max_threads_per_worker > threads:
        max_threads_per_worker = threads
    # calc nworkers and threads_per_worker
    # increase workers before threads_per_worker
    threads_per_worker = min_threads_per_worker
    for i in range(nfolders):
        if (float(i)*threads_per_worker) > threads:
            nworkers = i-1
            break
    else:
        nworkers = nfolders
        for i in range(min_threads_per_worker, max_threads_per_worker):
            if (float(nworkers)*i) > threads:
                threads_per_worker = i-1
                break
        else:
            threads_per_worker = max_threads_per_worker
    spare_threads = int(threads - (float(nworkers)*threads_per_worker))
    return nworkers, threads_per_worker, spare_threads


def parseArguments(args=None):
    """Read command-line arguments"""
    stages_err_msg = 'Invalid stage argument. Use \'-s [from]-[to]\' for \
numbers 1 through 4.'
    # get args
    if not args:
        args = createParser().parse_args()
    if args.details:
        print '\nThis is pG-lt version: ', open(os.path.join(
                                                os.path.dirname(__file__),
                                                '..', 'VERSION.txt')).read()
        print open(os.path.join(os.path.dirname(__file__), '..',
                                'DESCRIPTION.txt')).read()

        sys.exit()
    # check them
    if args.clean:
        clean()
        sys.exit('Files and folders deleted')
    if args.restart:
        return True, None, None, None, None, None
    if not args.email:
        # stop if no email
        sys.exit('An email address must be provided. Use \'-e\'.')
    # extract stages
    if not re.match('[1-4]-[1-4]', args.stages):
        sys.exit(stages_err_msg)
    startend = [int(e) for e in args.stages.split('-')]
    stages = [str(e) for e in range(startend[0], startend[1]+1)]
    if not stages:
        sys.exit(stages_err_msg)
    # check threads is a valid argument
    if args.threads == 0 or args.threads < -1:
        sys.exit('Invalid threads argument, must be -1 or >0.')
    return False, args.email, args.threads, args.verbose, args.debug, stages


def getFolders():
    """Return folders in directory with names.txt files"""
    # list all folders
    unchecked_dirs = [f for f in os.listdir('.') if not os.path.isfile(f)]
    # remove hidden folders
    unchecked_dirs = [d for d in unchecked_dirs if not re.match('^\.', d)]
    # loop through each and check they contain a names.txt
    checked_dirs = []
    for each in unchecked_dirs:
        path = os.path.join(os.getcwd(), each)
        files = os.listdir(path)
        if 'names.txt' in files:
            checked_dirs.append(each)
        # TODO: change this to have folders with any of the stage folders too
    if len(checked_dirs) > 0:
        return checked_dirs
    else:
        sys.exit(nonamestxt_msg)


def setUpLogging(verbose, debug, logname, directory=os.getcwd()):
    """Set up logging : direct and control log statements"""
    # get logger
    logger = logging.getLogger(logname)
    if debug:
        # log all statements above DEBUG level
        logger.setLevel(logging.DEBUG)
    else:
        # log all statements above INFO level
        # (which is higher than DEBUG)
        logger.setLevel(logging.INFO)
    # add file hander to root
    logfile = os.path.join(directory, 'log.txt')
    loghandler = logging.FileHandler(logfile, 'a')
    # set statement format -- I only want the message
    loghandler.setFormatter(logging.Formatter('%(message)s'))
    logger.addHandler(loghandler)
    if verbose:
        # if verbose, copy all info statements to console
        console = logging.StreamHandler()
        console.setFormatter(logging.Formatter('%(message)s'))
        logger.addHandler(console)
    logger.propagate = False
    return logger


def tearDownLogging(logname):
    """Remove a logger"""
    # get logger
    logger = logging.getLogger(logname)
    # remove handlers
    handlers = logger.handlers[:]
    for h in handlers:
        logger.removeHandler(h)


def createParser():
    """Create parser for command-line"""
    # Add new arg for threads
    parser = argparse.ArgumentParser()
    parser.add_argument("-email", "-e", help="please provide email \
for NCBI")
    parser.add_argument('--restart', help='restart pipeline if stopped',
                        action='store_true')
    parser.add_argument("-threads", "-t", help="number of threads, default\
 \'-1\', will use all available on machine", default=-1, type=int)
    parser.add_argument("-stages", "-s", help="stages to run, default \
\'1-4\'", default='1-4')
    parser.add_argument("--verbose", help="increase output verbosity",
                        action="store_true")
    parser.add_argument('--details', help='display information about the \
program', action='store_true')
    parser.add_argument("--debug", help="log warnings (developer only)",
                        action="store_true")
    parser.add_argument("--clean", help="remove all pG-lt files and \
folders (developer only)", action="store_true")
    return parser


def logMessage(phase, logger, folders=None, stage=None, threads=None,
               spare_threads=None, email=None, stages=None, counter=None):
    if phase == 'program-start':
        logger.info(description)
        logger.info('-' * 28 + ' Run details ' + '-' * 29)
        logger.info('Running on [{0}] [{1}]'.format(platform.node(),
                    platform.platform()))
        logger.info('Python [{0}]'.format(sys.version))
        logger.info('Using [{0}] threads with [{1}] spare'.
                    format(threads, spare_threads))
        logger.info('Using [{0}] as Entrez email'.format(email))
        logger.info('Running stages {0}'.format(stages))
        logger.info('Working with the following [{0}] folders:'.
                    format(len(folders)))
        # convert folders to string
        folder_string = ''
        chars_counter = 0
        for each in folders[:-1]:
            chars_counter += len(each)
            if chars_counter > 70:
                # stop at 70 columns
                folder_string += each + ',\n'
                chars_counter = 0
            else:
                folder_string += each + ', '
        folder_string += folders[-1]
        logger.info('[{0}]'.format(folder_string))
        logger.info('-' * 70 + '\n')
        logger.info('-' * 31 + ' Start ' + '-' * 32)
    elif phase == 'program-end':
        logger.info('-' * 32 + ' End ' + '-' * 33)
    elif phase == 'stage-start':
        logger.info('Stage [{0}] started at [{1}]'.format(stage, timestamp()))
    elif phase == 'stage-end':
        logger.info('Stage [{0}] finished at [{1}] for [{2}] folders'.
                    format(stage, timestamp(), counter))
    elif phase == 'program-restart':
        logger.info('{0}- Restarting [{1}] {0}'.format('-' * 11, timestamp()))
    else:
        raise(ValueError('Unrecognised phase'))


def prime(directory, arguments, threads):
    """Write pickle files, print arguments"""
    # Write pickle files
    temp_dir = os.path.join(directory, 'tempfiles')
    if not os.path.isdir(temp_dir):
        os.mkdir(temp_dir)
    with open(os.path.join(temp_dir, "genedict.p"), "wb") as file:
        pickle.dump(arguments['genedict'], file)
    with open(os.path.join(temp_dir, "paradict.p"), "wb") as file:
        pickle.dump(arguments['paradict'], file)
    with open(os.path.join(temp_dir, "terms.p"), "wb") as file:
        pickle.dump(arguments['terms'], file)
    with open(os.path.join(temp_dir, 'threads.p'), "wb") as file:
        pickle.dump(threads, file)
    with open(os.path.join(temp_dir, 'progress.p'), "wb") as file:
        pickle.dump(progress, file)
    # Print arguments and parameters to file
    record = 'Working with [{0}] names\n'.format(len(arguments['terms']))
    record += recordPars(arguments['paradict'])
    record += recordGpars(arguments['genedict'])
    with open(os.path.join(directory, 'info.txt'), 'wd') as file:
        file.write(record)


def timestamp():
    timestamp = datetime.today().strftime("%A, %d %B %Y %I:%M%p")
    return timestamp


def recordPars(paradict):
    """Return mpe parameters string"""
    record = '\nUsing the following parameters:\n'
    for key in paradict.keys():
        record += '    [{0}] = [{1}]\n'.format(key, paradict[key])
    return record


def recordGpars(genedict):
    """Return gene parameters string"""
    record = '\nUsing the following genes and gene parameters:\n'
    for gene in genedict.keys():
        record += '  Gene: [{0}]\n'.format(gene)
        for par in genedict[gene]:
            record += '    [{0}] = [{1}]\n'.format(par, genedict[gene][par])
    return record


def readInNames(directory):
    """Read names from text file in dir"""
    terms = []
    with open(os.path.join(directory, 'names.txt')) as names:
        for name in names:
            terms.append(name.strip())
    terms = [term for term in terms if not term == '']
    return terms


def readInGenePars(gpars_file, default_gpars_file):
    """Read gene_parameters.csv. Return list of dictionaries."""
    # TODO: too complex, consider breaking up
    def _read(gpars_file, template, genes=None):
        # open csv file and replace parameters in template
        # if they are None. If genes specified, only read
        # rows for those genes.
        with open(gpars_file, 'rb') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                if genes:
                    if not row['gene'] in genes:
                        continue
                temp = template.copy()
                for key in temp.keys():
                    if row[key]:
                        if temp[key] is None:
                            if key == 'names':
                                # for names, split into a list of syns
                                temp[key] = row[key].split(':')
                            else:
                                temp[key] = row[key]
                                genedict[row['gene']] = temp
        return genedict
    # check if file exists, else use default
    if not os.path.isfile(gpars_file):
        return readInGenePars(default_gpars_file, None)
    # genedicts
    genedict = {}
    # template of dict in genedict
    template = {'names': None, 'taxid': None, 'minlen': None, 'maxlen': None,
                'mingaps': None, 'minoverlap': None, 'minfails': None,
                'maxtrys': None, 'minseedsize': None, 'maxseedsize': None,
                'maxseedtrys': None, 'partition': None, 'type': None}
    # open file, read each row and fill in template
    genedict = _read(gpars_file, template)
    # if Nones, use defaults
    nones = False
    for gene in genedict.keys():
        for par in genedict[gene].keys():
            if genedict[gene][par] is None:
                nones = True
                break
    if nones:
        # run _read for defaults and limit to genes in genedict
        genedict = _read(default_gpars_file, template, genedict.keys())
    return genedict


def readInPars(pars_file, default_pars_file):
    """Read gene_parameters.csv. Return dictionary."""
    def _read(pars_file, paradict):
        # open csv, and replace all Nones
        with open(pars_file, 'rb') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                if paradict[row["Parameter"]] is None:
                    paradict[row["Parameter"]] = row["Value"]
        return paradict
    # check if file exists, else use default
    if not os.path.isfile(pars_file):
        return readInPars(default_pars_file, None)
    # template
    paradict = {'nseqs': None, 'naligns': None, 'nphylos': None,
                'thoroughness': None, 'maxtrys': None, 'rttpvalue': None,
                'parentid': None, 'outgroupid': None}
    # open file, read each row, extract value
    paradict = _read(pars_file, paradict)
    # if Nones remain, use default
    nones = False
    for key in paradict.keys():
        if paradict[key] is None:
            nones = True
            break
    if nones:
        paradict = _read(default_pars_file, paradict)
    return paradict


def sortArgs(directory, email, logger, default_pars_file, default_gpars_file):
    """Search for relevant files in dir, return list of arguments"""
    # find text file and read, raise error if fail
    try:
        terms = readInNames(directory)
    except IOError:
        logger.error(ioerror_msg.format('names.txt', directory))
        raise PrimingError()
    # find gene parameter file and read, raise error if fail
    try:
        genedict = readInGenePars(os.path.join(directory,
                                               'gene_parameters.csv'),
                                  default_gpars_file)
    except IOError:
        logger.error(ioerror_msg.format('gene_parameters.csv', directory))
        raise PrimingError()
    # find parameter file and read, raise error if fail
    try:
        paradict = readInPars(os.path.join(directory,
                                           'parameters.csv'),
                              default_pars_file)
    except IOError:
        logger.error(ioerror_msg.format('parameters.csv', directory))
        raise PrimingError()
    # add email to paradict
    paradict['email'] = email
    return {'terms': terms, 'genedict': genedict, 'paradict': paradict}
