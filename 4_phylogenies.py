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

def getRTTDists(phylo):
    names = []
    for terminal in phylo.get_terminals():
        names.append(terminal.name)
    rtt_dists = []
    for name in names:
        rtt_dists.append(phylo.distance(name))
    return rtt_dists

def getTBP(phylo):
    term_lens = []
    for terminal in phylo.get_terminals():
        term_lens.append(terminal.branch_length)
    total_len = phylo.total_branch_length()
    tbps = []
    for term_len in term_lens:
        tbps.append(term_len/total_len)
    return tbps

def getBranchLengths(phylo):
    lens = []
    depths =  phylo.depths(unit_branch_lengths = True)
    for branch in depths.keys():
        if branch.branch_length:
            lens.append(branch.branch_length)
    return lens
    

def phyloTest(phylo, max_branch):
    lens = getBranchLengths(phylo)
    total_len = sum(lens)
    lens_bool = [e/total_len > max_branch for e in lens]
    if any(lens_bool):
        return False
    else:
        return True

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
	niterations = 10
	
	## run phylogenies
	print "Generating [{0}] phylogenies ...".format(niterations)
	nphylos = 0
        counter = 0
        while nphylos < niterations:
            counter += 1
            print counter
            alignment = [random.sample(align_obj[e], 1)[0] for e in align_obj.keys()]
            phylo = pG.RAxML(alignment, method = "raxmlHPC")
            phylo.root_with_outgroup("outgroup")
            if phyloTest(phylo, 0.5):
                pG.Phylo.draw_ascii(phylo) # what does it look like?
                gene_str = "|".join(genes)
                pG.Phylo.write(phylo, os.path.join(output_dir, studies[i] + "_gene_" +\
                                                       gene_str + "_" + str(counter) + '.tre'), 'newick')
                nphylos += 1
	print 'Done. [{0}] phylogenies for [{1}].'.format(nphylos, studies[i])
	counter += 1
	nphylos_all += nphylos
print '\n\nStage finished. Generated [{0}] phylogenies across [{1}] studies.'.format(nphylos_all, counter)
