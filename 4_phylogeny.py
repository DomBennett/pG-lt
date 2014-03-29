#!/usr/bin/python
## MPE Stage 4: Phylogeny generation
## D.J. Bennett
## 24/03/2014

## Print stage
print "\n\nThis is stage 4: phylogenies\n"

## Packages
import os, re, random, pickle
from Bio import AlignIO
from Bio import Phylo
import dendropy as dp
from phylogeny_tools import *

## Dirs
alignment_dir = os.path.join(os.getcwd(),'3_alignment')
phylogeny_dir = os.path.join(os.getcwd(),'4_phylogeny')
if not os.path.isdir(phylogeny_dir):
	os.mkdir(phylogeny_dir)

## Input
with open("genedict.p", "rb") as file:
	genedict = pickle.load(file)
with open("paradict.p", "rb") as file:
	paradict = pickle.load(file)
with open("namesdict.p", "rb") as file:
	namesdict = pickle.load(file)

## Parameters
nphylos = int(paradict["nphylos"])
maxtrys = int(paradict["maxtrys"])
maxpedge = float(paradict["maxpedge"])
constraint = False
phylocounter = 0

## Process
genes = sorted(os.listdir(alignment_dir))
genes = [e for e in genes if not re.search("^\.|^log\.txt$", e)]
print "Reading in alignments"
alignobj = {}
for gene in genes:
	alignobj[gene] = []
	gene_dir = os.path.join(alignment_dir, gene)
	alignment_files = os.listdir(gene_dir)
	alignment_files = [e for e in alignment_files if not re.search("^\.", e)]
	for alignment_file in alignment_files:
		with open(os.path.join(gene_dir, alignment_file), "r") as file:
			alignment = AlignIO.read(file, "fasta")
		alignobj[gene].append(alignment)
print "Generating [{0}] phylogenies".format(nphylos)
phylogenies = []
trys = 0
while phylocounter < nphylos:
	print "Iteration [{0}]".format(phylocounter)
	if trys > maxtrys:
		break
	trys += 1
	alignments = [random.sample(alignobj[e], 1)[0] for e in genes]
	alignment,partitions = concatenateAlignments(alignments)
	if constraint:
		constraint = genConstraintTree(alignment, "taxontree.tre")
	phylogeny = RAxML(alignment, constraint = constraint, outgroup = "outgroup",\
		partitions = partitions)
	phylogeny.root_with_outgroup("outgroup")
	phylogeny.prune("outgroup")
	phylogeny = renameTips(phylogeny, namesdict)
	#Phylo.draw_ascii(phylogeny)
	if goodPhylogenyTest(phylogeny, maxpedge):
		phylogenies.append(phylogeny)
		phylocounter += 1
print 'Generating consensus.'
filepath = os.path.join(phylogeny_dir, 'distribution.tre')
with open(filepath, "w") as file:
	Phylo.write(phylogenies, file, 'newick')
phylos = dp.TreeList()
phylos.read_from_path(open(filepath, "rU"), "newick", as_rooted = True)
consensus = phylos.consensus(min_freq = 0.5, suppress_edge_lengths = True, rooted = True)
consensus.write_to_path(os.path.join(phylogeny_dir, "consensus.tre"), "newick")
print '\n\nStage finished. Generated [{0}] phylogenies.'.format(phylocounter)