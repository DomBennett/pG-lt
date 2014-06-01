#! /usr/bin/env python
## D.J. Bennett
## 24/03/2014
"""
MPE Stage 1: Names resolution
"""

## Packages
import os,csv,pickle,shutil,logging
import mpe.tools.names as ntools
from taxon_names_resolver import Resolver

def run():
	## print stage
	logging.info("\nStage 1: names resolution\n")

	## Dirs
	names_dir = '1_names'
	phylogeny_dir = '4_phylogeny'
	if not os.path.isdir(names_dir):
		os.mkdir(names_dir)
	if not os.path.isdir(phylogeny_dir):
		os.mkdir(phylogeny_dir)

	## Input
	with open(".paradict.p", "rb") as file:
		paradict = pickle.load(file)
	with open(".terms.p", "rb") as file:
		terms = pickle.load(file)

	## Parameters
	ntools.etools.Entrez.email = paradict["email"]

	## Process
	logging.info('Searching for taxids')
	try:
		parentid = int(paradict["parentid"])
	except:
		parentid = False
	resolver = Resolver(terms = terms, datasource = "NCBI", taxon_id = parentid)
	resolver.main()
	logging.info("Generating names dictionary")
	namesdict,allrankids = ntools.genNamesDict(resolver)
	logging.info('Generating taxonomic tree')
	taxontree,shared_lineages = ntools.genTaxTree(resolver, namesdict)

	## Output
	shutil.rmtree("resolved_names")
	with open(".namesdict.p", "wb") as file:
		pickle.dump(namesdict, file)
	with open(".allrankids.p", "wb") as file:
		pickle.dump(allrankids, file)
	headers = ["name", "unique_name", "rank", "NCBI_Taxids"]
	with open(os.path.join(names_dir, 'resovled_names.csv'), 'wb') as file:
		writer = csv.writer(file)
		writer.writerow(headers)
		for key in namesdict.keys():
			temp = namesdict[key]
			row = [key, temp["unique_name"], temp["rank"]]
			if len(temp["txids"]) > 1:
				ids = ""
				for each in temp["txids"]:
					ids += str(each) + "|"
			else:
				ids = temp["txids"][0]
			row.append(ids)
			writer.writerow(row)
	ntools.Phylo.write(taxontree, os.path.join(phylogeny_dir, "taxontree.tre"), "newick")
	logging.info('Stage finished. Resolved [{0}] names including outgroup.'.format(len(namesdict.keys())))