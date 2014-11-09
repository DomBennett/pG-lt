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
import logging
import platform
from datetime import datetime

# MESSAGES
description = """pG-lt: A pipeline for the automated generation of \
phylogenies from taxonomic names through Mass Phylogeny Estimation."""
nonamestxt_msg = '\nERROR: No folders containing \'names.txt\' files \
found! All taxonomic names should be placed in subdirectories and \
called: \'names.txt\''


# FUNCTIONS
def printHeader():
    """Print a nice program description header"""
    # use 70 cols as I think this is standard
    print '\n' + '#' * 70
    print description
    print '#' * 70 + '\n'


def parseArgs():
    """Read command-line arguments"""
    # Add new arg for threads
    parser = argparse.ArgumentParser()
    parser.add_argument("-email", "-e", help="please provide email \
for NCBI")
    parser.add_argument("-restart", "-r", help="restart from \
specified stage")
    parser.add_argument("--verbose", help="increase output verbosity",
                        action="store_true")
    parser.add_argument("--debug", help="log warnings (developer only)",
                        action="store_true")
    parser.add_argument("--clean", help="remove all pG-lt files and \
folders (developer only)", action="store_true")
    return parser


def getDirs(logger):
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
        logger.error(nonamestxt_msg)
        sys.exit()


def setUpLogging(verbose, debug, logname='base', directory=os.getcwd()):
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
    return logger


def tearDownLogging(logname):
    """Remove a logger"""
    # get logger
    logger = logging.getLogger(logname)
    # remove handlers
    handlers = logger.handlers[:]
    for h in handlers:
        logger.removeHandler(h)


def logMessage(phase, logger, directory=None):
    if phase == 'begin':
        # begin running the program
        # directory is list
        logger.info('#' * 70)
        logger.info(description)
        logger.info('#' * 70 + '\n')
        logger.info('-' * 28 + ' Run details ' + '-' * 29)
        logger.info('Running on [{0}] [{1}]'.format(platform.node(),
                    platform.platform()))
        logger.info('Python [{0}]'.format(sys.version))
        logger.info('Working with the following directories:')
        # convert dirs to string
        dir_string = ''
        chars_counter = 0
        for each in directory[:-1]:
            chars_counter += len(each)
            if chars_counter > 70:
                # stop at 70 columns
                dir_string += each + ',\n'
                chars_counter = 0
            else:
                dir_string += each + ', '
        dir_string += directory[-1]
        logger.info('[{0}]'.format(dir_string))
        logger.info('-' * 70 + '\n')
        logger.info('-' * 31 + ' Start ' + '-' * 32)
    elif phase == 'start':
        # start for one folder
        # directory is a string
        logger.info('Folder [{0}] started at [{1}]'.
                    format(directory,
                           datetime.today().strftime("%A, %d %B %Y %I:%M%p")))
    elif phase == 'finish':
        # when a folder is finished running
        logger.info('.... finished at [{1}]'.
                    format(directory,
                           datetime.today().strftime("%A, %d %B %Y %I:%M%p")))
    elif phase == 'folder-error':
        # when a folder is unable to run
        logger.info('.... unfinished at [{1}]'.
                    format(directory,
                           datetime.today().strftime("%A, %d %B %Y %I:%M%p")))
        logger.info('.... check [{0}] for details'.
                    format(os.path.join(directory, 'log.txt')))
    elif phase == 'end':
        logger.info('-' * 32 + ' End ' + '-' * 33)
    else:
        raise(ValueError('Unrecognised phase'))


def logError(msg, logger):
    """Return true when error raised, log informative message"""
    logger.error(msg)
    logger.info('.... Moving to next folder')
    return True
