#!/usr/bin/python
## MPE Stage 1: Names resolution
## D.J. Bennet
## 24/03/2014

## Print stage
print "\n\nStage 1: names resolution\n"

## Packages + Functions
import os, re, sys, csv, pickle, shutil
sys.path.append(os.path.join(os.getcwd(), 'functions'))
from names_tools import *
from taxon_names_resolver import TaxonNamesResolver

## Dirs
info_dir = os.path.join(os.getcwd(), '0_info')
names_dir = os.path.join(os.getcwd(), '1_names')
if not os.path.isdir(names_dir):
	os.mkdir(names_dir)
	
## Input
files = os.listdir(info_dir)
taxnames_file = sorted([e for e in files if re.search("^.*taxnames.txt$", e)])[0]
parentid_file = sorted([e for e in files if re.search("^.*parentID.txt$", e)])[0]
current_study = re.sub("_taxnames.txt$", "", taxnames_file)

## Process
print 'Searching for taxids'
with open(os.path.join(info_dir, parentid_file), "r") as file:
	parent_id = int(file.read())
input_file = os.path.join(info_dir, taxnames_file)
resolver = TaxonNamesResolver(input_file = input_file, datasource = "NCBI", taxon_id = parent_id)
resolver.main()
print "Generating names dictionary"
namesdict = genNamesDict(resolver)
print 'Generating taxonomic tree'
taxontree,shared_lineages = genTaxTree(resolver, namesdict, draw = True)

## Output
pickle.dump(namesdict, open(os.path.join(info_dir, "namesdict.p"), "wb"))
shutil.rmtree("resolved_names")
headers = ["MPE_id", "name", "unique_name", "rank", "Genbank_ids"]
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
Phylo.write(taxontree, os.path.join(info_dir, "taxontree.tre"), "newick")
print 'Stage finished. Resolved [{0}] names.'.format(len(namesdict.keys()))