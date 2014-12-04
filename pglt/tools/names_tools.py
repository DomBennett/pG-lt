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
from Bio import Phylo
from cStringIO import StringIO
import entrez_tools as etools
from taxon_names_resolver.misc_tools import taxTree
from system_tools import TaxonomicRankError


# FUNCTIONS
def genTaxTree(resolver, namesdict, logger, taxonomy=None, draw=False):
    """Generate Newick tree from TaxonNamesResolver class.

        Arguments:
         resolver = TaxonNamesResolver class
         by = Tip labels, either 'qnames', 'taxids' or 'name_string'
         draw = Draw ascii tree (logical)

         Return:
          Phylo"""
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
    treestring = taxTree(idents=idents, ranks=ranks, lineages=lineages,
                         taxonomy=taxonomy)
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
    """Return a dictionary containtaining all names and metadata"""
    # TODO: too complex, consider breaking up
    q_names = resolver.retrieve('query_name')
    q_names = [re.sub("\s", "_", e) for e in q_names]
    r_names = resolver.retrieve('classification_path')
    ranks = resolver.retrieve('classification_path_ranks')
    lineages = resolver.retrieve('classification_path_ids')
    lineages = [[int(e2) for e2 in e1] for e1 in lineages]
    allrankids = []
    [allrankids.extend(e) for e in lineages]
    allrankids = list(set(allrankids))
    namesdict = {}
    non_unique_lineages = []
    for i in range(len(lineages)):
        query_line = lineages.pop(i)
        best_bool = [True] * len(lineages)
        j = 0
        while j < len(lineages):
            subj_line = lineages[j]
            match_bool = [e not in subj_line for e in query_line]
            if sum(match_bool) < sum(best_bool):
                best_bool = match_bool
            if not any(best_bool):
                non_unique_lineages.append(i)
                break
            j += 1
        lineages.insert(i, query_line)
        if i in non_unique_lineages:
            continue
        txid = query_line[best_bool.index(True)]
        rank = ranks[i][best_bool.index(True)]
        rname = r_names[i][best_bool.index(True)]
        namesdict[q_names[i]] = {"txids": [txid], "unique_name": rname,
                                 "rank": rank}
    if non_unique_lineages:
        nul_ids = [lineages[e][-1] for e in non_unique_lineages]
        nul_qnames = [q_names[e] for e in non_unique_lineages]
        nul_ranks = [ranks[e][-1] for e in non_unique_lineages]
        i = 0
        while len(nul_ids) > 0:
            temp_id = nul_ids.pop(0)
            # find ids in the next level
            temp_children = etools.findChildren(str(temp_id), logger=logger,
                                                next=True)
            temp_children = [int(e) for e in temp_children]
            # if none are in allrankids, must be unique
            temp_children = [e for e in temp_children if e not in allrankids]
            if len(temp_children) > 0:
                namesdict[nul_qnames[i]] = {"txids": temp_children, "unique_na\
me": "Non-unique resolution", "rank": nul_ranks[i]}
            i += 1
    # if no parent id given, work one out
    if not parentid:
        shared_bool = []
        for each in lineages[0]:
            shared_bool.append(all([each in e for e in lineages]))
        parentid = lineages[0][shared_bool.index(False) - 1]
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
