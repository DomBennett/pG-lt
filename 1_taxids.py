#!/usr/bin/python
## MRes Project 2013
## Stage 1: Generating taxids
## In: 0_names | Out: 1_taxids
## 23/07/2013

## Print stage
print "\n\nThis is stage 1: taxids\n"

## Packages + Functions
import os, re, sys, csv
from Bio import Phylo
sys.path.append(os.path.join(os.getcwd(), 'functions'))
from taxon_names_resolver import TaxonNamesResolver
from taxon_names_resolver_tools import *

## Dirs
input_dir = os.path.join(os.getcwd(), '0_names')
output_dir = os.path.join(os.getcwd(), '1_taxids')
if not os.path.isdir(output_dir):
	os.mkdir(output_dir)
	
## Reading in taxadata
print "Reading in taxadata.csv ..."
taxadict = {}
with open(os.path.join(input_dir,'taxadata.csv'), 'rb') as csvfile:
	taxreader = csv.DictReader(csvfile)
	for row in taxreader:
		taxadict[row['study']] = [row['parentID'], row['sisterID']]
print "Done. Read in taxadata for [{0}] studies.".format(len(taxadict))
	
## Input + get taxids + output
print '\nLooping through studies ...'
files = os.listdir(input_dir)
taxnames_files = sorted([e for e in files if re.search("^.*taxnames.txt$", e)])
counter = 0
for taxnames_file in taxnames_files:
	current_study = re.sub("_taxnames.txt$", "", taxnames_file)
	print '\n\nWorking on: [{0}]'.format(current_study)

	## Get taxids
	print 'Searching for taxids ...'
	datasource = 'NCBI'
	taxon_id = int(taxadict[current_study][0])
	input_file = os.path.join(input_dir, taxnames_file)
	resolver = TaxonNamesResolver(input_file, datasource, taxon_id)
	resolver.main()
	#resolver = lineageMerge(resolver) # prevent duplicates
	qnames,taxids = extractHighestClade(resolver, by_ids = True)
	print '... resolved [{0}] names ...'.format(len(taxids))
	if len(taxids) < 3:
		print 'Too few names resovled. Dropping study.\n'
		continue
	print '... generating taxonomic tree ...'
	taxontree,shared_lineages = genTaxTree(resolver, by = 'taxids', draw = True)        
	print 'Done.'
	print 'Outputting ...'
	taxids.append(taxadict[current_study][1])
	qnames.append('outgroup') # add outgroup
	taxontree_file = re.sub('taxnames', 'taxontree', taxnames_file)
	taxids_file = re.sub('taxnames', 'taxids', taxnames_file)
	shared_file = re.sub('taxnames', 'shared', taxnames_file)
	qnames_file = re.sub('taxnames', 'qnames', taxnames_file)
	with file(os.path.join(output_dir, taxids_file), 'wb') as outfile:
		for taxid in taxids:
		  outfile.write("%s\n" % taxid)
	Phylo.write(taxontree, os.path.join(output_dir, taxontree_file), "newick")
	with file(os.path.join(output_dir, shared_file), 'wb') as outfile:
		for each in shared_lineages:
		  outfile.write("%s\n" % each)
	with file(os.path.join(output_dir, qnames_file), 'wb') as outfile:
		for each in qnames:
		  outfile.write("%s\n" % each)
	print 'Done.'
	counter += 1
print 'Stage finished. [{0}] studies with data.'.format(counter)
