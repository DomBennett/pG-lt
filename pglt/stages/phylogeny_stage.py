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
from pglt.tools.system_tools import MissingDepError


# RUN
def run(wd=os.getcwd(), logger=logging.getLogger('')):
    # PRINT STAGE
    logging.info("Stage 4: Phylogeny generation")

    # DIRS
    alignment_dir = os.path.join(wd, '3_alignment')
    phylogeny_dir = os.path.join(wd, '4_phylogeny')
    outfile = os.path.join(phylogeny_dir, 'distribution.tre')
    outfile_unconstrained = os.path.join(phylogeny_dir,
                                         'distribution_unconstrained.tre')
    temp_dir = os.path.join(wd, 'tempfiles')

    # CHECK DEPS
    if not ptools.raxml:
        raise MissingDepError('raxml')

    # INPUT
    with open(os.path.join(temp_dir, "paradict.p"), "rb") as file:
        paradict = pickle.load(file)
    with open(os.path.join(temp_dir, "genedict.p"), "rb") as file:
        genedict = pickle.load(file)
    with open(os.path.join(temp_dir, "allrankids.p"), "rb") as file:
        allrankids = pickle.load(file)

    # PARAMETERS
    nphylos = int(paradict["nphylos"])
    maxtrys = int(paradict["maxtrys"])
    rttstat = float(paradict["rttstat"])
    constraint = int(paradict["constraint"])
    ptools.logger = logger

    # READ ALIGMENTS
    clusters = sorted(os.listdir(alignment_dir))
    clusters = [e for e in clusters if not re.search("^\.|^log\.txt$", e)]
    logging.info("Reading in alignments ....")
    alignment_store = ptools.AlignmentStore(clusters=clusters,
                                            genedict=genedict,
                                            allrankids=allrankids,
                                            indir=alignment_dir, logger=logger)

    # GENERATE TREE DIST
    logging.info("Generating [{0}] phylogenies ....".format(nphylos))
    generator = ptools.Generator(alignment_store=alignment_store,
                                 rttstat=rttstat, outdir=phylogeny_dir,
                                 maxtrys=maxtrys, logger=logger, wd=temp_dir)
    if 1 == constraint:
        generator.constraint = False
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

    # RUN UNCONSTRAINED
    if 3 == constraint:
        logging.info('Repeating unconstrained ....')
        generator.phylogenies = []
        generator.constraint = False
        for i in range(nphylos):
            logging.info(".... Iteration [{0}]".format(i + 1))
            success = False
            while not success:
                success = generator.run()
        with open(outfile_unconstrained, "w") as file:
            counter = Phylo.write(generator.phylogenies, file, 'newick')

    # FINISH MESSAGE
    logging.info('Stage finished. Generated [{0}] phylogenies.'.
                 format(counter))
