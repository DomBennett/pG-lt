#!/usr/bin/python
## MPE: Run stage
## D.J. Bennett
## 25/03/2014

import argparse,pickle,csv,sys,os
from mpe.tools.system import Stager
from mpe import _PARS as default_pars
from mpe import _GPARS as default_gpars

def main():
	# Check args
	args = parser.parse_args()
	if args.parameters:
		parameters = args.parameters
	else:
		parameters = default_pars
	if args.genes:
		genes = args.genes
	else:
		genes = default_gpars
	# Input
	try:
		terms = []
		with open(args.names) as names:
			for name in names:
				terms.append(name.strip())
		terms = [term for term in terms if not term == '']
	except IOError:
		print "Names file could not be opened. File: [{0}]".args.names
		sys.exit()
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
	parser = argparse.ArgumentParser(description="Mass Phylogeny Estimation (MPE) - an automated pipeline\
		for the generation of phylogenies from taxonomic names (Author: D.J. Bennett).")
	parser.add_argument("-names", "-n", help=".txt file of taxonomic names.")
	parser.add_argument("-parameters", "-p", help=".csv file of parameters.")
	parser.add_argument("-genes", "-g", help=".csv file of gene parameters")
	parser.add_argument("-restart", "-r", action='store_true', default = False, help="restart process from last completed stage")
	main()