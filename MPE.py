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
import argparse,pickle,csv,sys,os
from mpe.tools.system import Stager
from mpe import _PARS as default_pars
from mpe import _GPARS as default_gpars

def main():
	args = parser.parse_args()
	if args.names:
		try:
			terms = []
			with open(args.names) as names:
				for name in names:
					terms.append(name.strip())
			terms = [term for term in terms if not term == '']
		except IOError:
			print "Names file could not be opened. File: \
[{0}]".args.names
			sys.exit()
	else:
		print "No names file given! Type 'MPE.py -h' for help."
		sys.exit()
	if args.parameters:
		parameters = args.parameters
	else:
		parameters = default_pars
	if args.genes:
		genes = args.genes
	else:
		genes = default_gpars
	try:
	 	genedict = {}
	 	with open(genes, 'rb') as csvfile:
	 		reader = csv.DictReader(csvfile)
	 		for row in reader:
	 			temp = {}
	 			temp["names"] = row["names"].split(":")
	 			temp["taxid"] = row["taxid"]
	 			temp["mingaps"] = row["mingaps"]
	 			temp["minoverlap"] = row["minoverlap"]
	 			temp["minfails"] = row["minfails"]
	 			temp["maxtrys"] = row["maxtrys"]
	 			temp["minseedsize"] = row["minseedsize"]
	 			temp["maxseedtrys"] = row["maxseedtrys"]
	 			genedict[row['gene']] = temp
	except IOError:
		print "Gene parameters file could not be opened. File: [{0}]".genes
		sys.exit()
	except KeyError:
		print "Unexpected headers in gene parameters."
		sys.exit()
	try:
	 	paradict = {}
	 	with open(parameters, 'rb') as csvfile:
	 		reader = csv.DictReader(csvfile)
	 		for row in reader:
	 			paradict[row["Parameter"]] = row["Value"]
	except IOError:
		print "Parameters file could not be opened. File: [{0}]".parameters
	 	sys.exit()
	# Output
 	with open(".genedict.p", "wb") as file:
 		pickle.dump(genedict, file)
 	with open(".paradict.p", "wb") as file:
 		pickle.dump(paradict, file)
 	with open(".terms.p", "wb") as file:
 		pickle.dump(terms, file)
 	if args.restart:
	 	# Find current stage
	 	if os.path.isfile('.stage.p'):
	 		with open('.stage.p', 'rb') as file:
	 			pickle.load(file)
	 	else:
	 		stage = 0
	 	if stage < 4:
			Stager.run_all(stage)
		else:
			print "All stages have already been run."
	else:
		Stager.run_all(0)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Mass Phylogeny Estimation (MPE) -\
an automated pipeline for the generation of phylogenies from taxonomic names (D.J.\
 Bennett (C) 2014).")
	parser.add_argument("-names", "-n", help=".txt file of taxonomic names.")
	parser.add_argument("-parameters", "-p", help=".csv file of parameters.")
	parser.add_argument("-genes", "-g", help=".csv file of gene parameters")
	parser.add_argument("-restart", "-r", action='store_true', default = False,\
		help="restart process from last completed stage")
	main()