#! /usr/bin/env python
# D.J. Bennett
# 24/03/2014
"""
pglt Stage 4: Phylogeny generation
"""

# PACKAGES
import os
import re
import pickle
import logging
from Bio import Phylo
import pglt.tools.phylogeny_tools as ptools


# RUN
def run(wd=os.getcwd(), logger=logging.getLogger('')):
    # PRINT STAGE
    logging.info("Stage 4: Phylogeny generation")

    # DIRS
    alignment_dir = os.path.join(wd, '3_alignment')
    phylogeny_dir = os.path.join(wd, '4_phylogeny')
    outfile = os.path.join(phylogeny_dir, 'distribution.tre')

    # INPUT
    with open(os.path.join(wd, ".paradict.p"), "rb") as file:
        paradict = pickle.load(file)
    with open(os.path.join(wd, ".genedict.p"), "rb") as file:
        genedict = pickle.load(file)
    with open(os.path.join(wd, ".allrankids.p"), "rb") as file:
        allrankids = pickle.load(file)

    # PARAMETERS
    nphylos = int(paradict["ntrees"])
    maxtrys = int(paradict["maxtrys"])
    rttpvalue = float(paradict["rttpvalue"])
    ptools.logger = logger

    # READ ALIGMENTS
    genes = sorted(os.listdir(alignment_dir))
    genes = [e for e in genes if not re.search("^\.|^log\.txt$", e)]
    genekeys = {}
    for gene in genes:
        genekeys[gene] = re.sub('_cluster[0-9]+', '', gene)
    logging.info("Reading in alignments ....")
    alignment_store = ptools.AlignmentStore(genes=genes, genedict=genedict,
                                            genekeys=genekeys,
                                            allrankids=allrankids,
                                            indir=alignment_dir)

    # GENERATE TREE DIST
    logging.info("Generating [{0}] phylogenies ....".format(nphylos))
    generator = ptools.Generator(alignment_store=alignment_store,
                                 rttpvalue=rttpvalue, outdir=phylogeny_dir,
                                 maxtrys=maxtrys)
    for i in range(nphylos):
        logging.info(".... Iteration [{0}]".format(i + 1))
        success = False
        while not success:
            success = generator.run()
    with open(outfile, "w") as file:
        counter = Phylo.write(generator.phylogenies, file, 'newick')

    # GENERATE CONSENSUS
    logging.info('Generating consensus ....')
    ptools.consensus(generator.phylogenies, phylogeny_dir, min_freq=0.5,
                     is_rooted=True, trees_splits_encoded=False)

    # FINISH MESSAGE
    logging.info('Stage finished. Generated [{0}] phylogenies.'.
                 format(counter))
