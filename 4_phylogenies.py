#!/usr/bin/python
## MRes Project 2013
## Stage 4: Generate phylogenies (YAY!!)
## In: 3_alignments | Out: 4_phylogenies
## 14/08/2013

## Parameters
rand_ngenes = False
min_ngenes = 2
niterations = 100
max_branch = 0.5 # the maximum proportion of a tree a single branch can represent
max_nsplits = False
max_trys = 1000

## Print stage
print "\n\nThis is stage 4: phylogenies\n"

## Packages
import sys, os, re, random
import dendropy as dp
import numpy as np
sys.path.append(os.path.join(os.getcwd(), 'functions'))
import phyloGenerator_adapted as pG
from phylogeny_gen_tools import *

## dirs
names_dir = os.path.join(os.getcwd(), '1_taxids')
names_files = os.listdir(names_dir)
input_dir = os.path.join(os.getcwd(), '3_alignments')
output_dir = os.path.join(os.getcwd(), '4_phylogenies')
if not os.path.isdir(output_dir):
    os.mkdir(output_dir)

## Loop through studies
studies = sorted(os.listdir(input_dir))
studies = [st for st in studies if not re.search("^\.|^log\.txt$", st)]
study_counter = 0
print 'Looping through studies ...'
nphylos_all = 0
for i in range(len(studies)):
	print '\n\nWorking on: [{0}]'.format(studies[i])
        taxontree_file = os.path.join(names_dir, "{0}_taxontree.txt".\
                                          format(studies[i]))
        study_names_files = [e for e in names_files if\
                                 re.search("^{0}".format(studies[i]), e)]
        taxids_file = [e for e in study_names_files if\
                           re.search("taxids", e)][0]
        qnames_file = [e for e in study_names_files if\
                           re.search("qnames", e)][0]
        names = {}
        with file(os.path.join(names_dir, taxids_file), 'rb') as taxfile:
                    with file(os.path.join(names_dir, qnames_file), 'rb')\
                            as qfile:
                        for txid,qname in zip(taxfile, qfile):        
                            txid = "tx" + txid.strip()
                            # spaces need to be repalced with _
                            names[txid] = re.sub("\\s", "_", qname.strip())

	study_dir = os.path.join(input_dir, studies[i])
	genes = sorted(os.listdir(study_dir))
        genes = [e for e in genes if not re.search("^\.", e)]
	print "Reading in alignments ..."
	align_obj = {}
	for gene in genes:
		align_obj[gene] = []
		gene_dir = os.path.join(study_dir, gene)
		a_files = os.listdir(gene_dir)
                a_files = [e for e in a_files if not re.search("^\.", e)]
		for a_file in a_files:
			a_file_path = os.path.join(gene_dir, a_file)
			align = pG.AlignIO.read(a_file_path, "fasta")
			align_obj[gene].append(align)
	naligns = [len(align_obj[e]) for e in align_obj.keys()]
	print "Done. Read in [{0}] alignments for [{1}]".format(naligns,align_obj.keys())
	## run phylogenies
	print "Generating [{0}] phylogenies ...".format(niterations)
	nphylos = 0
        counter = 0
        phylos = []
        while nphylos < niterations:
            if counter > max_trys:
                break
            counter += 1
            print counter
            if rand_ngenes:
                ngenes = random.sample(range(min_ngenes, len(align_obj.keys()) + 1), 1)[0]
                genes = random.sample(align_obj.keys(), ngenes)
            else:
                genes = align_obj.keys()
            print genes
            alignment = [random.sample(align_obj[e], 1)[0] for e in genes]
            if "darwin" == sys.platform:
                phylo = pG.RAxML(alignment, method = 'localVersion')
            else:
                phylo = pG.RAxML(alignment, method = "raxmlHPC")
            pG.Phylo.draw_ascii(phylo)
            try:
                phylo.root_with_outgroup("outgroup")
            except ValueError:
                continue
            phylo.prune("outgroup")
            phylo = renameTips(phylo, names)
            pG.Phylo.draw_ascii(phylo)
            if phyloTest(phylo, taxontree_file, max_nsplits, max_branch):
                print "A good tree!"
                phylos.append(phylo)
                nphylos += 1
            else:
                print "A bad tree..."
        if len(phylos) > 1:
            print 'Done. [{0}] phylogenies for [{1}].'.format(nphylos, studies[i])
            print 'Writing out phylogenies and generating consensus.'
            filepath = os.path.join(output_dir, studies[i] + "_phylos" + '.tre')
            pG.Phylo.write(phylos, filepath, 'newick')
            phylos = dp.TreeList()
            phylos.read(open(filepath, "rU"), "newick")
            consensus = phylos.consensus(min_freq = 0.5)
            consensus.write_to_path(os.path.join(output_dir, studies[i] + "_consensus" + '.tre'),\
                                        "newick", suppress_edge_lengths = True)
            study_counter += 1
            nphylos_all += nphylos
print '\n\nStage finished. Generated [{0}] phylogenies across [{1}] studies.'.format(nphylos_all, study_counter)
