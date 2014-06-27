#! /usr/bin/env python
## D.J. Bennett
## 24/03/2014
"""
MPE Stage 4: Phylogeny generation
"""

## Packages
import os,re,random,pickle,logging
from Bio import AlignIO
from Bio import Phylo
import dendropy as dp
import mpe.tools.phylogeny as ptools

def run(wd = os.getcwd()):
	## Print stage
	logging.info("\nGenerating phylogenies\n")

	## Dirs
	alignment_dir = os.path.join(wd, '3_alignment')
	phylogeny_dir = os.path.join(wd,'4_phylogeny')

	## Input
	with open(os.path.join(wd, ".paradict.p"), "rb") as file:
		paradict = pickle.load(file)

	## Parameters
	nphylos = int(paradict["ntrees"])
	maxtrys = int(paradict["maxtrys"])
	maxpedge = float(paradict["maxpedge"])
	constraint = True
	phylocounter = 0

	## Process
	genes = sorted(os.listdir(alignment_dir))
	genes = [e for e in genes if not re.search("^\.|^log\.txt$", e)]
	logging.info("Reading in alignments")
	alignobj = {}
	for gene in genes:
		alignobj[gene] = []
		gene_dir = os.path.join(alignment_dir, gene)
		alignment_files = os.listdir(gene_dir)
		alignment_files = [e for e in alignment_files if not\
		re.search("^\.", e)]
		for alignment_file in alignment_files:
			with open(os.path.join(gene_dir, alignment_file), "r")\
			as file:
				alignment = AlignIO.read(file, "fasta")
			alignobj[gene].append((alignment, alignment_file))
	logging.info("Generating [{0}] phylogenies".format(nphylos))
	phylogenies = []
	trys = 0
	while phylocounter < nphylos:
		logging.info("Iteration [{0}]".format(phylocounter))
		if trys > maxtrys:
			break
		trys += 1
		alignments,alignment_files = zip(*[random.sample(alignobj[e],\
			1)[0] for e in genes])
		logging.info("Using alignments.....")
		for gene,alignment_file in zip(genes,alignment_files):
			logging.info("....[{0}]:[{1}]".format(gene, alignment_file))
		alignment,partitions = ptools.concatenateAlignments(list(alignments))
		if constraint:
			constraint = ptools.genConstraintTree(alignment,\
				os.path.join(phylogeny_dir, "taxontree.tre"))
		phylogeny = ptools.RAxML(alignment, constraint = constraint,\
			outgroup = "outgroup",partitions = partitions)
		phylogeny.root_with_outgroup("outgroup")
		phylogeny.prune("outgroup")
		if ptools.goodPhylogenyTest(phylogeny, maxpedge):
			phylogenies.append(phylogeny)
			phylocounter += 1
	logging.info('Generating consensus.')
	filepath = os.path.join(phylogeny_dir, 'distribution.tre')
	with open(filepath, "w") as file:
		Phylo.write(phylogenies, file, 'newick')
	phylos = dp.TreeList()
	phylos.read_from_path(filepath, "newick", as_rooted = True)
	consensus = phylos.consensus(min_freq = 0.5, suppress_edge_lengths = True,\
		rooted = True)
	consensus.write_to_path(os.path.join(phylogeny_dir, "consensus.tre"),\
		"newick")
	logging.info('\n\nStage finished. Generated [{0}] phylogenies.'.\
		format(phylocounter))