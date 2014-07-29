#! /usr/bin/env python
## D.J. Bennett
## 24/03/2014
"""
mpe Stage 1: Names resolution
"""

## Packages
import os,pickle,shutil,logging,sys
import mpe.tools.names as ntools
from mpe.tools.system import TooFewSpeciesError
from taxon_names_resolver import Resolver

## Informative error msgs
error_msg =  'It is likely that one or more names have \
been resolved incorrectly, as such the parent taxonomic \
group has been set to Eukaryotes which is too high a \
taxonomic rank for phylogenetic analysis. Consider \
adding a parent ID to the parameters.csv to prevent \
incorrect names resolution.'

def run(wd = os.getcwd()):
	## print stage
	logging.info("Names resolution\n")

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
	ntools.etools.Entrez.email = paradict["email"]
	minspecies = 5

	## Process
	logging.info('Searching for taxids....')
	logging.info('------TaxonNamesResolver:Start------')
	try:
		parentid = int(paradict["parentid"])
	except:
		parentid = False
	if len(terms) < minspecies:
		raise TooFewSpeciesError
	resolver = Resolver(terms = terms, datasource = "NCBI",\
		taxon_id = parentid)
	resolver.main()
	logging.info('------TaxonNamesResolver:End------')
	logging.info("Generating names dictionary....")
	try:
		namesdict,allrankids = ntools.genNamesDict(resolver)
	except ntools.TaxonomicRankError:
		print error_msg
		sys.exit()
	logging.info('Generating taxonomic tree....')
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
	logging.info('Resolved [{0}] names \
including outgroup.'.format(len(namesdict.keys())))