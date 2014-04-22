#!/usr/bin/python

from distutils.core import setup

setup(name = 'MPE', version = 1.0,\
	description = 'Automated phylogeny generation',\
	author = "D. J. Bennett",\
	author_email = 'dominic.john.bennett@gmail.com',\
	py_modules = ['run', '1_names', '2_download', '3_alignment',\
	'4_phylogeny', 'sys_tools', 'stages', 'binaries', 'names_tools',\
	'taxon_names_resolver', 'entrez_tools', 'download_tools',\
	'phylogeny_tools'],
	package_data = [('parameters': ['parameters.csv']),\
	('gene_parameters', ['gene_parameters.csv'])])