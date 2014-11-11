#! /usr/bin/env python
# D.J. Bennett
# 24/03/2014
"""
pglt Stage 1: Names resolution
"""

# PACKAGES
import os
import pickle
import shutil
import logging
import pglt.tools.names_tools as ntools
from pglt.tools.system_tools import TooFewSpeciesError
from taxon_names_resolver import Resolver


# RUN
def run(wd=os.getcwd(), logger=logging.getLogger('')):
    # PRINT STAGE
    logger.info("Stage 1: Names resolution")

    # DIRS
    names_dir = os.path.join(wd, '1_names')
    phylogeny_dir = os.path.join(wd, '4_phylogeny')
    if not os.path.isdir(names_dir):
        os.mkdir(names_dir)
    if not os.path.isdir(phylogeny_dir):
        os.mkdir(phylogeny_dir)

    # INPUT
    with open(os.path.join(wd, ".paradict.p"), "rb") as file:
        paradict = pickle.load(file)
    with open(os.path.join(wd, ".terms.p"), "rb") as file:
        terms = pickle.load(file)

    # PARAMETERS
    outgroupid = paradict["outgroupid"]
    ntools.etools.Entrez.email = paradict["email"]
    minspecies = 5
    ntools.logger = logger

    # PROCESS
    logger.info('Searching for taxids ....')
    logger.info('------TaxonNamesResolver:Start------')
    try:
        parentid = paradict["parentid"]
    except:
        parentid = False
    if len(terms) < minspecies:
        raise TooFewSpeciesError
    resolver = Resolver(terms=terms, datasource="NCBI", taxon_id=parentid)
    resolver.main()
    logger.info('------TaxonNamesResolver:End------')
    logger.info("Generating names dictionary ....")
    namesdict, allrankids, parentid = ntools.genNamesDict(resolver, parentid)
    logger.info("Finding an outgroup ....")
    namesdict = ntools.getOutgroup(namesdict, parentid, outgroupid)
    # add outgroup ids to allrankids
    allrankids.extend(namesdict['outgroup']['txids'])
    logger.info('Generating taxonomic tree ....')
    taxontree, shared_lineages = ntools.genTaxTree(resolver, namesdict)

    # OUTPUT
    # remove temp TNR folder
    shutil.rmtree("resolved_names")
    # write out changes to hidden pickled files
    with open(os.path.join(wd, ".namesdict.p"), "wb") as file:
        pickle.dump(namesdict, file)
    with open(os.path.join(wd, ".allrankids.p"), "wb") as file:
        pickle.dump(allrankids, file)
    # write namesdict as csv
    ntools.writeNamesDict(names_dir, namesdict)
    # write taxon tree
    ntools.Phylo.write(taxontree, os.path.join(phylogeny_dir, "taxontree.tre"),
                       "newick")

    # FINISH MESSAGE
    logger.info('Stage finished. Resolved [{0}] names including outgroup.'.
                format(len(namesdict.keys())))