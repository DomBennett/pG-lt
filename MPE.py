#!/usr/bin/python
## MPE: Run stage
## D.J. Bennet
## 25/03/2014
#TODO: Run each stage, run for multiple names files, check parameters

import argparse,pickle,csv,sys

def main():
	# Check args
	args = parser.parse_args()
	if args.parameters:
		parameters = args.parameters
	else:
		parameters = "parameters.csv"
	if args.genes:
		genes = args.genes
	else:
		genes = "gene_parameters.csv"
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
 	with open("genedict.p", "wb") as file:
 		pickle.dump(genedict, file)
 	with open("paradict.p", "wb") as file:
 		pickle.dump(paradict, file)
 	with open("terms.p", "wb") as file:
 		pickle.dump(terms, file)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Mass Phylogeny Estimation (MPE) - an automated pipeline\
		for the generation of phylogenies from taxonomic names (Author: D.J. Bennett).")
	parser.add_argument("-names", "-n", help=".txt file of taxonomic names.")
	parser.add_argument("-parameters", "-p", help=".csv file of parameters.")
	parser.add_argument("-genes", "-g", help=".csv file of gene parameters")
	main()