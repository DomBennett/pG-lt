#!/usr/bin/python
## MPE Stage 1: Names resolution
## D.J. Bennett
## 24/03/2014

## Packages
import os, csv, pickle, shutil
import mpe.tools.names as ntools
from taxon_names_resolver import Resolver

## Print stage
print "\n\nStage 1: names resolution\n"

## Dirs
names_dir = '1_names'
phylogeny_dir = '4_phylogeny'
if not os.path.isdir(names_dir):
	os.mkdir(names_dir)
if not os.path.isdir(phylogeny_dir):
	os.mkdir(phylogeny_dir)

## Input
with open(".genedict.p", "rb") as file:
	genedict = pickle.load(file)
with open(".paradict.p", "rb") as file:
	paradict = pickle.load(file)
with open(".terms.p", "rb") as file:
	terms = pickle.load(file)

## Parameters
ntools.etools.Entrez.email = paradict["email"]

## Process
print 'Searching for taxids'
try:
	parentid = int(paradict["parentid"])
except:
	parentid = False
resolver = Resolver(terms = terms, datasource = "NCBI", taxon_id = parentid)
resolver.main()
print "Generating names dictionary"
namesdict,allrankids = ntools.genNamesDict(resolver)
print 'Generating taxonomic tree'
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
print 'Stage finished. Resolved [{0}] names including outgroup.'.format(len(namesdict.keys()))
stage = 1
with open(".stage.p", "wb") as file:
	pickle.dump(stage, file)