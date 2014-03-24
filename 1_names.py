#!/usr/bin/python
## MPE Stage 1: Names resolution
## D.J. Bennet
## 24/03/2014

## Print stage
print "\n\nStage 1: names resolution\n"

## Packages + Functions
import os, re, sys, csv, pickle
sys.path.append(os.path.join(os.getcwd(), 'functions'))
from Bio import Phylo
from taxon_names_resolver import TaxonNamesResolver
from names_tools import *

## Dirs
info_dir = os.path.join(os.getcwd(), '0_info')
names_dir = os.path.join(os.getcwd(), '1_names')
if not os.path.isdir(names_dir):
	os.mkdir(names_dir)
	
## Input
files = os.listdir(info_dir)
input_file = sorted([e for e in files if re.search("^.*taxnames.txt$", e)])
current_study = re.sub("_taxnames.txt$", "", taxnames_file)
print '\n\nWorking on: [{0}]'.format(current_study)

## Process
print 'Searching for taxids ...'
taxon_id = int(taxadict[current_study][0])
input_file = os.path.join(input_dir, taxnames_file)
resolver = TaxonNamesResolver(input_file = input_file, datasource = "NCBI", taxon_id = taxon_id)
resolver.main()
namesdict = genNamesDict(resolver)
print '... generating taxonomic tree ...'
taxontree,shared_lineages = genTaxTree(resolver, namesdict, draw = True)

## Output

print 'Stage finished. Resolved [{0}] names.'.format(len(namesdict.keys()))