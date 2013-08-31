#!/usr/bin/python
## MRes Project 2013
## Stage 4: Generate phylogenies (YAY!!)
## In: 3_alignments | Out: 4_phylogenies
## 14/08/2013

## Print stage
print "\n\nThis is stage 4: phylogenies\n"

## Packages
import sys, os, re, random
import numpy as np
sys.path.append(os.path.join(os.getcwd(), 'functions'))
import phyloGenerator_adapted as pG

## dirs
input_dir = os.path.join(os.getcwd(), '3_alignments')
output_dir = os.path.join(os.getcwd(), '4_phylogenies')
if not os.path.isdir(output_dir):
    os.mkdir(output_dir)


## Loop through studies
a_files = os.listdir(input_dir)
a_files = [e for e in a_files if not re.search("^log\.txt$", e)]
counter = 0
studies = list(set([re.sub("_gene_.*_[0-9]*\.faa$", "", e) for e in a_files]))
print 'Looping through studies ...'
nphylos_all = 0
for i in range(len(studies)):

	## what study?
	print '\n\nWorking on: [{0}]'.format(studies[i])
	
	## reading in alignments
	print "Reading in alignments ..."
	pattern = "^" + studies[i]
	study_files = [e for e in a_files if re.search(pattern, e)]
	aligns = []
	for study_file in study_files:
		gene = re.split("_gene_|_[0-9]*\.faa", study_file)[-2]
		align = pG.AlignIO.read(os.path.join(input_dir, study_file), "fasta")
		aligns.append((align, gene)) # a list of tuples of the phylogeny AND the gene that made it
	naligns = len(aligns)
	print "Done. Read in [{0}] alignments".format(naligns)
	
	## run phylogenies
	print "Generating [{0}] phylogenies ...".format(naligns)
	nphylos = 0
	for j in range(naligns):
		align_obj = aligns[j][0]
		gene = aligns[j][1]
		try:
			phylo = pG.RAxML(align_obj)
		except:
			print '... Unexpected error: RAxML ...'
			raxml_files = os.listdir(os.getcwd())
			raxml_files = [f for f in raxml_files if \
				re.search("(^RAxML.*$|^.*\.phylip$|^.*\.fasta$)", f)]
			for f in raxml_files:
				os.remove(f)
			continue	
		raxml_files = os.listdir(os.getcwd())
		raxml_files = [f for f in raxml_files if \
			re.search("(^RAxML.*$|^.*\.phylip$|^.*\.fasta$)", f)]
		for f in raxml_files:
			os.remove(f)
		pG.Phylo.write(phylo, os.path.join(output_dir, studies[i] + "_gene_" +\
			gene + "_" + str(j) + '.tre'), 'newick')
		nphylos += 1
	print 'Done. [{0}] phylogenies for [{1}].'.format(nphylos, studies[i])
	counter += 1
	nphylos_all += nphylos
print '\n\nStage finished. Generated [{0}] phylogenies across [{1}] studies.'.format(nphylos_all, counter)
