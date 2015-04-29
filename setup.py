#! /usr/bin/env python
# D.J. Bennett
# 26/05/2014
"""
setup.py for pglt
"""
import os
import pglt
from setuptools import setup, find_packages

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
    package_data={'pglt': ['parameters.csv', 'gene_parameters.csv',
                           'dependencies.p']},
    scripts=['run_pglt.py', 'pglt_set_dependencies.py'],
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
print('''
Congratulations -- you`ve installed pglt!
Now use `pglt_set_dependencies.py` to add external programs in order to run
the pipeline.
''')
