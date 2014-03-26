#!/usr/bin/python
## MPE Stage 1: Names resolution
## D.J. Bennett
## 24/03/2014

## Print stage
print "\n\nStage 1: names resolution\n"

## Packages
import os, re, sys, csv, pickle, shutil
from names_tools import *
from taxon_names_resolver import TaxonNamesResolver

## Dirs
names_dir = os.path.join(os.getcwd(),'1_names')
if not os.path.isdir(names_dir):
	os.mkdir(names_dir)

## Input
with open("genedict.p", "rb") as file:
	genedict = pickle.load(file)
with open("paradict.p", "rb") as file:
	paradict = pickle.load(file)
with open("terms.p", "rb") as file:
	terms = pickle.load(file)

## Parameters
Entrez.email = paradict["email"]

## Process
print 'Searching for taxids'
try:
	parentid = int(paradict["parentid"])
except:
	parentid = False
resolver = TaxonNamesResolver(terms = terms, datasource = "NCBI", taxon_id = parentid)
resolver.main()
print "Generating names dictionary"
namesdict,all_ids = genNamesDict(resolver)
print 'Generating taxonomic tree'
taxontree,shared_lineages = genTaxTree(resolver, namesdict)

## Output
shutil.rmtree("resolved_names")
with open("namesdict.p", "wb") as file:
	pickle.dump(namesdict, file)
with open("all_ids.p", "wb") as file:
	pickle.dump(all_ids, file)
headers = ["MPE_id", "name", "unique_name", "rank", "NCBI_Taxids"]
with open(os.path.join(names_dir, 'resovled_names.csv'), 'wb') as file:
	writer = csv.writer(file)
	writer.writerow(headers)
	for key in namesdict.keys():
		temp = namesdict[key]
		row = [key, temp["name"], temp["unique_name"], temp["rank"]]
		if len(temp["ids"]) > 1:
			ids = ""
			for each in temp["ids"]:
				ids += str(each) + "|"
		else:
			ids = temp["ids"][0]
		row.append(ids)
		writer.writerow(row)
Phylo.write(taxontree, "taxontree.tre", "newick")

## Print stageout
print 'Stage finished. Resolved [{0}] names including outgroup.'.format(len(namesdict.keys()))