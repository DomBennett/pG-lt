#! /usr/bin/env python
# D.J. Bennett
# 29/04/2015
"""
Set dependencies for pG-lt
"""
import sys
import os
import subprocess
import re
import argparse
import pickle
from tabulate import tabulate
from pglt import _RAXML as raxml
from pglt import _MAFFT as mafft
from pglt import _MAFFTQ as mafftq
from pglt import _MAFFTX as mafftx
from pglt import _BLASTN as blastn
from pglt import _ROOT as root

# GLOBALS
abserr_msg = 'Absolute paths to dependencies only. [{0}] is not an absolute \
path.'


# FUNCTIONS
def createParser():
    """Create parser for command-line"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-mafft", help="path to `mafft`, by default will search \
current working directory", type=str)
    parser.add_argument("-mafftq", help="path to `mafftq`, by default will \
search current working directory", type=str)
    parser.add_argument("-mafftx", help="path to `mafftx`, by default will \
search current working directory", type=str)
    parser.add_argument("-raxml", help="path to `raxml`, by default will search \
current working directory", type=str)
    parser.add_argument("-blastn", help="path to `blastn`, by default will \
search current working directory", type=str)
    parser.add_argument('--overwrite', help='overwrite dependencies already \
added to pG-lt', action='store_true')
    parser.add_argument('--local', help='add deps to local pglt for testing \
purposes', action='store_true')
    return parser


def runCommand(args):
    """Run and return command"""
    # run, read and kill
    process = subprocess.Popen(args, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
    info = process.stdout.read()
    process.kill()
    return info.strip()


def getVersion(args):
    """Return version number of program"""
    # first check its absolute
    if not os.path.isabs(args[0]):
        print abserr_msg.format(args[0])
        return
    # check if program exists
    try:
        info = runCommand(args)
        # find version number, extract and strip of digits
        pattern = '(v|version)?\s?[0-9]\.[0-9]+'
        res = re.search(pattern, info)
        version = info[res.span()[0]:res.span()[1]]
        non_decimal = re.compile(r'[^\d.]+')
        version = non_decimal.sub('', version)
        return float(version)
    except OSError:
        return False


if __name__ == '__main__':
    args = createParser().parse_args()
    print '\nSearching and setting dependencies to pG-lt ....'
    if args.local:
        print '.... adding deps to local copy of pglt, not installed version.'
        folders = os.listdir(os.getcwd())
        if 'pglt' in folders and 'tests' in folders:
            root = os.path.join(os.getcwd(), 'pglt')
        else:
            sys.exit('Error: to run locally, ensure you`re in downloaded pG-lt/')
    with open(os.path.join(root, 'dependencies.p'), "rb") as file:
        depsdict = pickle.load(file)
    changes = False
    # MAFFT
    if args.mafft and not mafft or args.overwrite:
        version = getVersion([args.mafft, '-u'])
        if version and version > 7.0:
            depsdict['mafft'] = args.mafft
            changes = True
        else:
            print 'No MAFFT detected -- requires MAFFT v7+'
    # MAFFTQ
    if args.mafftq and not mafftq or args.overwrite:
        version = getVersion([args.mafftq, '-u'])
        if version and version > 7.0:
            depsdict['mafftq'] = args.mafftq
            changes = True
        else:
            print 'No mafft-qinsi detected'
    # MAFFTX
    if args.mafftx and not mafftx or args.overwrite:
        version = getVersion([args.mafftx, '-u'])
        if version and version > 7.0:
            depsdict['mafftx'] = args.mafftx
            changes = True
        else:
            print 'No mafft-xinsi detected'
    # RAXML
    if args.raxml and not raxml or args.overwrite:
        version = getVersion([args.raxml, '-version'])
        if version and version > 7.0:
            depsdict['raxml'] = args.raxml
            changes = True
        else:
            print 'No RAxML detected -- requires RAxML v7+'
    # BLAST
    if args.blastn and not blastn or args.overwrite:
        version = getVersion([args.blastn, '-h'])
        if version and version > 2.0:
            depsdict['blastn'] = args.blastn
            changes = True
        else:
            print 'No Stand-alone BLAST detected -- requires BLAST suite v2+'
    if changes:
        with open(os.path.join(root, 'dependencies.p'), "wb") as file:
            pickle.dump(depsdict, file)
        print 'Done.\n'
    else:
        print 'No dependencies given.\n'
    # CHECK
    # create table
    print('-'*70)
    print(' ' * 23 + 'Set dependencies' + ' ' * 23)
    print('-'*70)
    deps = [('mafft', '3'), ('mafftq', '-'), ('mafftx', '-'),
            ('raxml', '4'), ('blastn', '2-3')]
    rows = []
    headers = ['Dependency', 'Status', 'Stages required', 'Path']
    for d, e in deps:
        if depsdict[d]:
            rows.append([d, 'Present', e, depsdict[d]])
        else:
            rows.append([d, 'Absent', e, '-'])
    print(tabulate(rows, headers, tablefmt="simple") + '\n')
    if args.local:
        print('\nUse `python setup.py test` to test functionality of pG-lt on \
your machine with these deps.')
