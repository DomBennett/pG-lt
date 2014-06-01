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
 (Copyright NCBI (C) 2009). It also uses the following python packages:
 Biopython (Copyright Cook (C) 2009), Dendropy (Copyright Sukumaran and Holder (C)
 2010) and Taxon Names Resovler (Copyright (C) Bennett 2014).

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

description = """MPE D.J. Bennett (C) 2014

A pipeline for the automated generation of phylogenies from taxonomic
names through Mass Phylogeny Estimation."""

def parseArgs():
	"""Read arguments"""
	parser = argparse.ArgumentParser()
	parser.add_argument("-names", "-n", help=".txt file of taxonomic names.")
	parser.add_argument("-parameters", "-p", help=".csv file of parameters.")
	parser.add_argument("-genes", "-g", help=".csv file of gene parameters")
	parser.add_argument("-restart", "-r", help="restart from specified stage")
	parser.add_argument("--verbose", help="increase output verbosity",
					action="store_true")
	parser.add_argument("--development", help="log warnings (developer only)",
					action="store_true")
	return parser

if __name__ == '__main__':
	print '\n' + '#' * 70
	print description
	print '#' * 70 + '\n'
	parser = parseArgs()
	stage,verbose,dev = readArgs(parser.parse_args())
	Stager.run_all(stage = stage, verbose = verbose, dev = dev)