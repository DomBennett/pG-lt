#! /usr/bin/env python
# D.J. Bennett
# 24/03/2014
"""
pglt names tools
"""

# PACKAGES
import re
import csv
import os
import random
from Bio import Phylo
from cStringIO import StringIO
import entrez_tools as etools
from taxon_names_resolver import TaxDict
from taxon_names_resolver import taxTree
from system_tools import TaxonomicRankError


# FUNCTIONS
def genTaxTree(resolver, namesdict, logger, taxonomy=None, draw=False):
    """Return Phylo from TaxonNamesResolver class."""
    ranks = resolver.retrieve('classification_path_ranks')
    qnames = resolver.retrieve('query_name')
    lineages = resolver.retrieve('classification_path')
    # replace ' ' with '_' for taxon tree
    qnames = [re.sub("\s", "_", e) for e in qnames]
    resolved_names_bool = [e in namesdict.keys() for e in qnames]
    ranks = [ranks[ei] for ei, e in enumerate(resolved_names_bool) if e]
    lineages = [lineages[ei] for ei, e in enumerate(resolved_names_bool) if e]
    # identify unresolved names
    unresolved_names = [qnames[ei] for ei, e in enumerate(resolved_names_bool)
                        if not e]
    idents = [qnames[ei] for ei, e in enumerate(resolved_names_bool) if e]
    statement = "Unresolved names: "
    for each in unresolved_names:
        statement += " " + each
    logger.debug(statement)
    # make taxdict
    taxdict = TaxDict(idents=idents, ranks=ranks, lineages=lineages,
                      taxonomy=taxonomy)
    # make treestring
    treestring = taxTree(taxdict)
    if not taxonomy:
        d = 22  # default_taxonomy + 1 in tnr
    else:
        d = len(taxonomy) + 1
    # add outgroup
    treestring = '({0},outgroup:{1});'.format(treestring[:-1], float(d))
    tree = Phylo.read(StringIO(treestring), "newick")
    if draw:
        Phylo.draw_ascii(tree)
    return tree


def genNamesDict(resolver, logger, parentid=None):
    """Return a dictionary containing all names and metadata"""
    # extract lists from resolver
    qnames = resolver.retrieve('query_name')
    qnames = [re.sub("\s", "_", e) for e in qnames]  # no spaces
    ranks = resolver.retrieve('classification_path_ranks')
    # in order to get name, ID and rank for the context let's use a zipped
    #  lineage IDs
    lids = resolver.retrieve('classification_path_ids')
    lineages = zip(resolver.retrieve('classification_path'), lids,
                   resolver.retrieve('classification_path_ranks'))
    lineages = [zip(e1, e2, e3) for e1, e2, e3 in lineages]
    # make taxdict
    taxdict = TaxDict(idents=qnames, ranks=ranks, lineages=lineages)
    # get all ranks
    allrankids = []
    [allrankids.extend(e) for e in lids]
    allrankids = [int(e) for e in allrankids]  # make sure ints
    allrankids = list(set(allrankids))
    # init namesdict
    namesdict = {}
    # loop through taxdict
    for key in taxdict.keys():
        cident = taxdict[key]['cident']  # Contextual data
        if cident:
            # if there is context data ....
            cname, cident, crank = cident  # unpack context name, ID and rank
            namesdict[key] = {"txids": [int(cident)], "unique_name": cname,
                              "rank": crank}
        else:
            # find unclaimed IDs using allrankids and searching children
            rident = taxdict[key]['ident']  # Resolved ID
            rank = taxdict[key]['rank']  # Resolved rank
            # find ids in the next level
            children = etools.findChildren(rident, logger=logger,
                                           next=True)
            if children:
                unclaimed = [int(e) for e in children]
                unclaimed = [e for e in unclaimed if e not in allrankids]
                # choose random subset of unclaimed
                if len(unclaimed) > 5:
                    txids = random.sample(unclaimed, 5)
                # if there are no unclaimed, just use children
                if not unclaimed:
                    txids = children
            else:
                # if no children, just use rident
                txids = [str(rident)]
            namesdict[key] = {"txids": txids,
                              "unique_name": 'Non-unique resolution',
                              "rank": rank}
    # if no parent id given, work one out
    if not parentid:
        shared_bool = []
        for each in lineages[0]:
            shared_bool.append(all([each in e for e in lineages]))
        parentid = lineages[0][shared_bool.index(False) - 1]
        parentid = int(parentid[1])  # second one in tuple
    return namesdict, allrankids, parentid


def getOutgroup(namesdict, parentid, logger, outgroupid=None, minrecords=1000):
    """Return namesdict with suitable outgroup"""
    # TODO: too complex, consider breaking up
    def findParent(parentid):
        return etools.eFetch(parentid, logger=logger,
                             db="taxonomy")[0]['ParentTaxId']

    def getTaxIdMetaData(ncbi_id):
        etal_bool = False
        if len(ncbi_id) > 1:
            ncbi_id = ncbi_id[0]
            etal_bool = True
        record = etools.eFetch(ncbi_id, logger=logger, db="taxonomy")[0]
        metadata = [record['Rank'], record['ScientificName']]
        if etal_bool:
            metadata = [e + ' et al.' for e in metadata]
        return metadata[0], metadata[1]
    # loop until a suitable outgroup is found. Criteria are:
    #  1. ids returned must belong to a sister group of all ids of
    #   names given
    #  2. ids must have nucleotide data (i.e.avoid returning extinct organisms)
    # assumptions:
    #  1. NCBI taxonomy is not paraphyletic
    # make sure parentid is string
    if not outgroupid:
        parentid = str(parentid)
        outgroup_ids = []
        while not outgroup_ids:
            # if parent id are Cellular Orgs, likely name resolution error
            #  or names given are too diverse
            if parentid == '131567':
                raise TaxonomicRankError()
            # get parent of parent
            grandparentid = findParent(parentid)
            # find all children
            candidates = etools.findChildren(grandparentid, logger=logger,
                                             next=True)
            # filter out children that are in ingroup
            candidates = [e for e in candidates if e != parentid]
            # search genbank for nuc records
            for candidate in candidates:
                term = 'txid' + str(candidate) + '[PORGN]'
                nuc_record = etools.eSearch(term, logger=logger)
                # there must be more than 1000 nuc records
                if int(nuc_record['Count']) > minrecords:
                    outgroup_ids.append(candidate)
            # make grandparentid the new parentid
            parentid = grandparentid
    else:
        outgroup_ids = [outgroupid]
    # add outgroup_ids to namesdict
    rank, unique_name = getTaxIdMetaData(outgroup_ids)
    # convert to ints
    outgroup_ids = [int(e) for e in outgroup_ids]
    namesdict["outgroup"] = {"txids": outgroup_ids, "unique_name": unique_name,
                             "rank": rank}
    return namesdict


def writeNamesDict(directory, namesdict):
    headers = ["name", "unique_name", "rank", "NCBI_Taxids"]
    with open(os.path.join(directory, 'resolved_names.csv'), 'wb') as file:
        writer = csv.writer(file)
        writer.writerow(headers)
        for key in namesdict.keys():
            temp = namesdict[key]
            row = [key, temp["unique_name"], temp["rank"]]
            if len(temp["txids"]) > 1:
                ids = ""
                for each in temp["txids"]:
                    ids += str(each) + "|"
            else:
                ids = temp["txids"][0]
            row.append(ids)
            writer.writerow(row)
