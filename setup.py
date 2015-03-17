#! /usr/bin/env python
# D.J. Bennett
# 26/05/2014
"""
setup.py for pglt
"""
import os
import subprocess
import re
import sys
import shutil
import pglt
from setuptools import setup, find_packages


# FUNCTIONS
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


def moveExternal(name):
    """Find location, and move to pglt folder"""
    try:
        location = runCommand(['which', name])
    except OSError:
        pass
    shutil.copy(location, os.path.join('pglt', name))


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

# CHECK FOR RAxML, MAFFT and BLAST
# '-u' prevents printing to shell
all_present = True
if not getVersion(['mafft', '-u']):
    print 'No MAFFT detected -- requires MAFFT v7+'
    all_present = False
if getVersion(['mafft', '-u']) < 7.0:
    print 'MAFFT detected too old -- requires MAFFT v7+'
if not getVersion(['mafft-xinsi', '-u']):
    print 'No mafft-xinsi detected -- requires installation of\
 mafft with RNA structural alignments'
    all_present = False
if not getVersion(['mafft-qinsi', '-u']):
    print 'No mafft-qinsi detected -- requires installation of\
 mafft with RNA structural alignments'
    all_present = False
if getVersion(['raxml', '-version']) < 7.0:
    print 'No RAxML detected -- requires RAxML v7+'
    all_present = False
if getVersion(['blastn', '-h']) < 2.0:
    print 'No Stand-alone BLAST detected -- requires BLAST suite v2+'
    all_present = False
if not all_present:
    sys.exit('Unable to install/test! Please install missing external programs')

# MOVE EXTERNALS TO PGLT FOLDER
moveExternal('raxml')
moveExternal('mafft')
moveExternal('mafft-xinsi')
moveExternal('mafft-qinsi')
moveExternal('blastn')


# PACAKAGE INFO
PACKAGES = find_packages()
PACKAGE_DIRS = [p.replace(".", os.path.sep) for p in PACKAGES]

# SETUP
setup(
    name="pglt",
    version=pglt.__version__,
    author="Dominic John Bennett",
    author_email="dominic.john.bennett@gmail.com",
    description=("pG-lt: An automated pipeline for phylogeney generation."),
    license="LICENSE.txt",
    keywords="phylogenetics ecology evolution conservation",
    url="https://github.com/DomBennett/MassPhylogenyEstimation",
    packages=PACKAGES,
    package_dir=dict(zip(PACKAGES, PACKAGE_DIRS)),
    package_data={'pglt': ['parameters.csv', 'gene_parameters.csv', 'raxml',
                           'mafft', 'mafft-qinsi', 'mafft-xinsi', 'blastn']},
    scripts=['run_pglt.py', 'pglt_farm.py'],
    test_suite='tests',
    long_description=pglt.__doc__,
    classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 2.7",
        "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
    ],
    install_requires=['setuptools', 'taxon_names_resolver', 'biopython',
                      'dendropy', 'numpy', 'tabulate'],
)
