#! /bin/usr/env python
## D.J. Bennett
## 26/05/2014
"""
MPE is a pipeline for the automated generation of phylogenies through 'Mass
Phylogeny Estimation'. This program is built on top of phyloGenerator (C) 2013
and was written by D.J. Bennett with additional help from W.D. Pearse and L. Hudson.
This program makes use of external programs for phylogeny generation and bioinformatics
these are: RAxML (Copyright (C) Stamatakis 2013) , MAFFT (Copyright (C) 2013 Kazutaka
Katoh) the NCBI's standalone BLAST suite 2.2.29+ and online API services
 (Copyright NCBI (C) 2009). It also uses a variety of python packages including:
 Biopython (Copyright Cook (C) 2009) and Dendropy (Copyright Sukumaran and Holder (C)
 2010).

Copyright (C) 2014  Dominic John Bennett

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
"""

## Packages
import argparse
from mpe.tools.system import Stager
from mpe.tools.system import readArgs

## Variables
description = """
MPE (D.J. Bennett (C) 2014)

A pipeline for the automated generation of phylogenies from taxonomic
names through Mass Phylogeny Estimation.
"""

##TODO: python interaction

## Terminal interaction
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument("-names", "-n", help=".txt file of taxonomic names.")
	parser.add_argument("-parameters", "-p", help=".csv file of parameters.")
	parser.add_argument("-genes", "-g", help=".csv file of gene parameters")
	parser.add_argument("--verbose", help="increase output verbosity",
					action="store_true")
	parser.add_argument("-restart", "-r", action='store_true', default = False,\
		help="restart process from last completed stage")
	# Read args
	args = parser.parse_args()
	stage = readArgs(args)
	# Print + run
	if args.parameters:
		parameters = args.parameters
	else:
		parameters = 'DEFAULT'
	if args.genes:
		genes = args.genes
	else:
		genes = 'DEFAULT'
	print '\n' + description + '\n'
	print 'Using ...'
	print '.... [{0}] names'.format(args.names)
	print '.... [{0}] parameters'.format(parameters)
	print '.... [{0}] genes\n'.format(genes)
	Stager.run_all(stage = stage, verbose = args.verbose)