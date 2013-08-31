#!/usr/bin/python
## MRes Project 2013
## Stage 1: Generating taxids
## In: 0_names | Out: 1_taxids
## 23/07/2013

## Print stage
print "\n\nThis is stage 1: taxids\n"

## Packages + Functions
import os, re, sys, csv
sys.path.append(os.path.join(os.getcwd(), 'functions'))
from taxon_names_resolver import TaxonNamesResolver
from gen_tax_tree import genTaxTree

## Dirs
input_dirs = os.path.join(os.getcwd(), '0_names')
output_dir = os.path.join(os.getcwd(), '1_taxids')
if not os.path.isdir(output_dir):
    os.mkdir(output_dir)
    
## Reading in taxadata
print "Reading in taxadata.csv ..."
taxadict = {}
with open(os.path.join(input_dir, 'taxadata.csv'), 'rb') as csvfile:
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
	print '###### TaxonNamesResolver ######'
	resolver = TaxonNamesResolver(taxnames_file, datasource, indir =\
		input_dir, outdir = output_dir, taxon_id = taxon_id)
	resolver.main()
	taxids = resolver.extract('taxonids')
	lineages = resolver.extract('lineages')
	ranks = resolver.extract('rank_paths')
	qnames = resolver.extract('qnames')
	print '################################'
	print '... resolved [{0}] names ...'.format(len(taxids))
	if len(taxids) < 3:
		print 'Too few names resovled. Dropping study.\n'
		continue
	print 'Done.'
	
	## Generate taxon trees
	#print 'Creating taxon tree ...'
	#taxontree, shared_lineages = genTaxTree(taxids, lineages, ranks)
	#print 'Done.'
	
	## Output
	print 'Outputting ...'
	taxids.append(taxadict[current_study][1])
	qnames.append('outgroup') # add outgroup
	#taxontree_file = re.sub('taxnames', 'taxontree', taxnames_file)
	taxids_file = re.sub('taxnames', 'taxids', taxnames_file)
	#shared_file = re.sub('taxnames', 'shared', taxnames_file)
	qnames_file = re.sub('taxnames', 'qnames', taxnames_file)
	with file(os.path.join(output_dir, taxids_file), 'wb') as outfile:
		for taxid in taxids:
		  outfile.write("%s\n" % taxid)
	#with file(os.path.join(output_dir, taxontree_file), 'wb') as outfile:
	#	outfile.write(taxontree)
	#with file(os.path.join(output_dir, shared_file), 'wb') as outfile:
	#	for each in shared_lineages:
	#	  outfile.write("%s\n" % each)
	with file(os.path.join(output_dir, qnames_file), 'wb') as outfile:
		for each in qnames:
		  outfile.write("%s\n" % each)
	print 'Done.'
	counter += 1
print 'Stage finished. [{0}] studies with data.'.format(counter)
