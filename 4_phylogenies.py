#!/usr/bin/python
## MRes Project 2013
## Stage 4: Generate phylogenies (YAY!!)
## In: 3_alignments | Out: 4_phylogenies
## 14/08/2013

## Parameters
max_ngenes = 4 # the maximum number genes to be used in phylogeny estimation
min_ngenes = 1 # the minimum ...

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
studies = sorted(os.listdir(input_dir))
studies = [st for st in studies if not re.search("^log\.txt$", st)]
counter = 0
print 'Looping through studies ...'
nphylos_all = 0
for i in range(len(studies)):

	## what study?
	print '\n\nWorking on: [{0}]'.format(studies[i])
	study_dir = os.path.join(os.getcwd(), input_dir, studies[i])
	genes = sorted(os.listdir(study_dir))

	## reading in alignments
	print "Reading in alignments ..."
	align_obj = {}
	for gene in genes:
		align_obj[gene] = []
		gene_dir = os.path.join(study_dir, gene)
		a_files = os.listdir(gene_dir)
		for a_file in a_files:
			a_file_path = os.path.join(gene_dir, a_file)
			align = pG.AlignIO.read(a_file_path, "fasta")
			align_obj[gene].append(align)
	naligns = [len(align_obj[e]) for e in align_obj.keys()]
	print "Done. Read in [{0}] alignments for [{1}]".format(naligns,align_obj.keys())
	
	## working out combinations
	#try:
	#	niterations = int(np.prod(naligns)) # N.B. this is not the true number of combinations (depends on min and max genes)
	#except OverflowError:
	#	niterations = 100
	#if niterations > 100: # limit to 100
	#	niterations = 100
	niterations = 100
	
	## run phylogenies
	print "Generating [{0}] phylogenies ...".format(niterations)
	# reset min and max ngenes if too many or too few
	temp_max_ngenes = max_ngenes
	temp_min_ngenes = min_ngenes
	if len(align_obj) < max_ngenes:
		temp_max_ngenes = len(align_obj)

	if len(align_obj) < min_ngenes:
		temp_min_ngenes = len(align_obj)

	ngenes_range = range(temp_min_ngenes, temp_max_ngenes+1)
	nphylos = 0
	for j in range(niterations):
		ngenes = random.sample(ngenes_range, 1)[0]
		genes = random.sample(align_obj.keys(), ngenes)
		alignment = [random.sample(align_obj[e], 1)[0] for e in genes]
		try:
			phylo = pG.RAxML(alignment, method = "raxmlHPC")
			#pG.Phylo.draw_ascii(phylo) # what does it look like?
			# TODO: phylo.depths()
		except:
			print '... Unexpected error: RAxML ...'
			raxml_files = os.listdir(os.getcwd())
			raxml_files = [f for f in raxml_files if \
				re.search("(^RAxML.*$|^.*\.phylip$|^.*\.fasta$)", f)]
			for f in raxml_files:
				os.remove(f)
			continue
		gene_str = "|".join(genes)
		pG.Phylo.write(phylo, os.path.join(output_dir, studies[i] + "_gene_" +\
			gene_str + "_" + str(j) + '.tre'), 'newick')
		nphylos += 1
	print 'Done. [{0}] phylogenies for [{1}].'.format(nphylos, studies[i])
	counter += 1
	nphylos_all += nphylos
print '\n\nStage finished. Generated [{0}] phylogenies across [{1}] studies.'.format(nphylos_all, counter)
