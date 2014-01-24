import sys, os, re, random
import dendropy as dp
from Bio import Phylo

def renameTips(phylo, names):
    for each in phylo.get_terminals():
        try:
            each.name = names[each.name]
        except KeyError:
            pass
    return phylo

def getRTTDists(phylo):
    """Calcualte root to tips distance"""
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
    
def phyloTest(query, ref_path, max_nsplits, max_branch):
    lens = getBranchLengths(query)
    total_len = sum(lens)
    lens_bool = [e/total_len > max_branch for e in lens]
    if any(lens_bool):
        return False
    elif max_nsplits:
        tax_tree = dp.Tree()
        tax_tree.read_from_path(ref_path, "newick")
        temp_tree = dp.Tree()
        Phylo.write(query, "temp_phylo.txt", 'newick')
        temp_tree.read_from_path("temp_phylo.txt", "newick")
        os.remove("temp_phylo.txt")
        # I don't think it's necessary for them to have the same number of taxa
        tax_list = [e for e in tax_tree.taxon_set]
        temp_list = [e for e in temp_tree.taxon_set]
        drop_list = [e1 for e1 in tax_list if e1.label not in [e2.label for e2 in temp_list]]
        tax_tree.prune_taxa(drop_list)
        print tax_tree.as_ascii_plot()
        # extract number of splits in the tax_tree not in temp_tree
        # http://pythonhosted.org/DendroPy/tutorial/treestats.html#frequency-of-a-split-in-a-collection-of-trees
        nsplits = tax_tree.false_positives_and_negatives(temp_tree)
        print "... [{0}] tax splits, [{1}] temp splits...".format(nsplits[0], nsplits[1])
        if nsplits[1] > max_nsplits:
            return False
        else:
            return True
    else:
        return True


