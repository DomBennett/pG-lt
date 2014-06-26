#! /usr/bin/env python
## D.J. Bennett
## 24/03/2014
"""
MPE Stage 1: Names resolution
"""

## Packages
import os,pickle,shutil,logging
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
	resolver = Resolver(terms = terms, datasource = "NCBI",\
		taxon_id = parentid)
	resolver.main()
	logging.info("Generating names dictionary")
	namesdict,allrankids = ntools.genNamesDict(resolver)
	logging.info('Generating taxonomic tree')
	taxontree,shared_lineages = ntools.genTaxTree(resolver,\
		namesdict)

	## Output
	# remove temp TNR folder
	shutil.rmtree("resolved_names")
	# write out changes to hidden pickled files
	with open(".namesdict.p", "wb") as file:
		pickle.dump(namesdict, file)
	with open(".allrankids.p", "wb") as file:
		pickle.dump(allrankids, file)
	# write namesdict as csv
	ntools.writeNamesDict(names_dir, namesdict)
	# write taxon tree
	ntools.Phylo.write(taxontree, os.path.join(phylogeny_dir,\
		"taxontree.tre"), "newick")

	## Finish message
	logging.info('Stage finished. Resolved [{0}] names \
including outgroup.'.format(len(namesdict.keys())))