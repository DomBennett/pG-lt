#! /usr/bin/env python
## D.J. Bennett
## 24/03/2014
"""
mpe Stage 1: Names resolution
"""

## Packages
import os,pickle,shutil,logging
import mpe.tools.names_tools as ntools
from mpe.tools.system_tools import TooFewSpeciesError
from taxon_names_resolver import Resolver

def run(wd = os.getcwd()):
	## print stage
	logging.info("Stage 1: Names resolution")

	## Dirs
	names_dir = os.path.join(wd, '1_names')
	phylogeny_dir = os.path.join(wd, '4_phylogeny')
	if not os.path.isdir(names_dir):
		os.mkdir(names_dir)
	if not os.path.isdir(phylogeny_dir):
		os.mkdir(phylogeny_dir)

	## Input
	with open(os.path.join(wd, ".paradict.p"), "rb") as file:
		paradict = pickle.load(file)
	with open(os.path.join(wd, ".terms.p"), "rb") as file:
		terms = pickle.load(file)

	## Parameters
	outgroupid = paradict["outgroupid"]
	ntools.etools.Entrez.email = paradict["email"]
	minspecies = 5

	## Process
	logging.info('Searching for taxids ....')
	logging.info('------TaxonNamesResolver:Start------')
	try:
		parentid = int(paradict["parentid"])
	except:
		parentid = False
	if len(terms) < minspecies:
		raise TooFewSpeciesError
	resolver = Resolver(terms = terms, datasource = "NCBI",\
		taxon_id = parentid) ## TODO: make tnr accept strings
	resolver.main()
	logging.info('------TaxonNamesResolver:End------')
	logging.info("Generating names dictionary ....")
	namesdict,allrankids,parentid = ntools.genNamesDict(resolver)
	logging.info("Finding an outgroup ....")
	namesdict = ntools.getOutgroup(namesdict, parentid, outgroupid)
	# add outgroup ids to allrankids
	allrankids.extend(namesdict['outgroup']['txids'])
	logging.info('Generating taxonomic tree ....')
	taxontree,shared_lineages = ntools.genTaxTree(resolver,\
		namesdict)

	## Output
	# remove temp TNR folder
	shutil.rmtree("resolved_names")
	# write out changes to hidden pickled files
	with open(os.path.join(wd, ".namesdict.p"), "wb") as file:
		pickle.dump(namesdict, file)
	with open(os.path.join(wd, ".allrankids.p"), "wb") as file:
		pickle.dump(allrankids, file)
	# write namesdict as csv
	ntools.writeNamesDict(names_dir, namesdict)
	# write taxon tree
	ntools.Phylo.write(taxontree, os.path.join(phylogeny_dir,\
		"taxontree.tre"), "newick")

	## Finish message
	logging.info('Stage finished. Resolved [{0}] names \
including outgroup.'.format(len(namesdict.keys())))