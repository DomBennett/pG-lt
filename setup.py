import os
from setuptools import setup, find_packages

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

PACKAGES = find_packages()
PACKAGE_DIRS = [p.replace(".", os.path.sep) for p in PACKAGES]

setup(
    name = "mpe",
    version = "0.0.1",
    author = "Dominic John Bennett",
    author_email = "dominic.john.bennett@gmail.com",
    description = ("An automated pipeline for phylogeney generation."),
    license = "No license",
    keywords = "ecology evolution phylogenetics",
    url = "https://github.com/DomBennett/MassPhylogenyEstimation",
    packages = PACKAGES,
    package_dir = dict(zip (PACKAGES, PACKAGE_DIRS)),
    package_data = {'mpe':['parameters.csv','gene_parameters.csv']},
    scripts = ['MPE.py'],
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 1 - Planning",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 2.7",
    ],
    install_requires=[
          # -*- Extra requirements: -*-
          'setuptools',
          'taxon_names_resolver',
          'biopython',
          'dendropy',
          'numpy',
      ],
)