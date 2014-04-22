#!/usr/bin/python
## MPE Phylogeny tools
## D.J. Bennett
## 24/03/2014

## Packages
import sys, os, re, random, pickle
import dendropy as dp
from Bio import Phylo
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio import Phylo
from Bio import AlignIO
import numpy as np
from sys_tools import *

## Globals
with open("programdict.p", "rb") as file:
 	programdict = pickle.load(file)
raxmlpath = programdict['raxml']

## Functions
def renameTips(phylo, names):
	for each in phylo.get_terminals():
		try:
			each.name = names[each.name]["name"]
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
	
def goodPhylogenyTest(query, max_branch):
	lens = getBranchLengths(query)
	total_len = sum(lens)
	lens_bool = [e/total_len > max_branch for e in lens]
	if any(lens_bool):
		return False
	else:
		return True

# def phyloTest(query, ref_path, max_nsplits, max_branch):
# 	lens = getBranchLengths(query)
# 	total_len = sum(lens)
# 	lens_bool = [e/total_len > max_branch for e in lens]
# 	if any(lens_bool):
# 		return False
# 	elif max_nsplits:
# 		tax_tree = dp.Tree()
# 		tax_tree.read_from_path(ref_path, "newick")
# 		temp_tree = dp.Tree()
# 		Phylo.write(query, "temp_phylo.txt", 'newick')
# 		temp_tree.read_from_path("temp_phylo.txt", "newick")
# 		os.remove("temp_phylo.txt")
# 		# I don't think it's necessary for them to have the same number of taxa
# 		tax_list = [e for e in tax_tree.taxon_set]
# 		temp_list = [e for e in temp_tree.taxon_set]
# 		drop_list = [e1 for e1 in tax_list if e1.label not in [e2.label for e2 in temp_list]]
# 		tax_tree.prune_taxa(drop_list)
# 		print tax_tree.as_ascii_plot()
# 		# extract number of splits in the tax_tree not in temp_tree
# 		# http://pythonhosted.org/DendroPy/tutorial/treestats.html#frequency-of-a-split-in-a-collection-of-trees
# 		nsplits = tax_tree.false_positives_and_negatives(temp_tree)
# 		print "... [{0}] tax splits, [{1}] temp splits...".format(nsplits[0], nsplits[1])
# 		if nsplits[1] > max_nsplits:
# 			return False
# 		else:
# 			return True
# 	else:
# 		return True

def genConstraintTree(alignment, taxontree_file):
	tip_names = []
	for record in alignment:
		tip_names.append(record.id)
	with open(taxontree_file, "r") as file:
		constraint = Phylo.read(file, "newick")
	constraint_tips = []
	for terminal in constraint.get_terminals():
		constraint_tips.append(terminal.name)
	tips_to_drop = [e for e in constraint_tips if not e in tip_names]
	for tip in tips_to_drop:
		constraint.prune(tip)
	with open("constraint.tre", "w") as file:
		Phylo.write(constraint, file, "newick")
	if constraint.is_bifurcating():
		return " -r constraint.tre"
	else:
		return " -g constraint.tre"

def concatenateAlignments(alignments):
	if len(alignments) == 1:
		return alignments[0],False
	# Sort IDs
	alignment_ids = []
	for gene in alignments:
		gene_ids = []
		for rec in gene:
			gene_ids.append(rec.id)
		alignment_ids.append(gene_ids)
	all_ids = []
	[all_ids.extend(e) for e in alignment_ids]
	all_ids = list(set(all_ids))
	# Concatenate
	alignment = MultipleSeqAlignment([])
	for txid in all_ids:
		sequence = ""
		for i,gene in enumerate(alignments):
			if txid in alignment_ids[i]:
				sequence += gene[alignment_ids[i].index(txid)].seq
			else:
				sequence += "-" * gene.get_alignment_length()
		sequence = SeqRecord(sequence, id = txid, description = "multigene sequence")
		alignment.append(sequence)
	# Get partitions
	lengths = [e.get_alignment_length() for e in alignments]
	partitions = [0].append(np.cumsum(lengths))
	return alignment,partitions

def RAxML(alignment, outgroup=None, partitions=None, constraint=None, timeout=999999999):
	"""Adapted pG function: Generate phylogeny from alignment using RAxML (external program)."""
	input_file = 'phylogeny_in.phylip'
	output_file = 'phylogeny_out'
	file_line = ' -s ' + input_file + ' -n ' + output_file
	options = ' -p ' + str(random.randint(0,10000000))
	if outgroup:
		options += ' -o ' + outgroup
	with open(input_file, "w") as file:
		AlignIO.write(alignment, file, "phylip-relaxed")
	if len(alignment) > 100:
		dnamodel = ' -m GTRCAT'
	else:
		dnamodel = ' -m GTRGAMMA'
	if partitions:
		with open("partitions.txt", 'w') as f:
			for i in range(0, len(partitions)-1):
				f.write("DNA, position" + str(partitions[i]) + " = " + str(partitions[i]+1) +\
					"-" + str(partitions[i+1]) + "\n")
		options += " -q " + "partitions.txt"
	if constraint:
		options += constraint
	command_line = raxmlpath + file_line + dnamodel + options
	#print command_line
	pipe = TerminationPipe(command_line)
	pipe.run()
	if not pipe.failure:
		with open('RAxML_bestTree.' + output_file, "r") as file:
			tree = Phylo.read(file, "newick")	
		if constraint:
			os.remove('constraint.tre')
		if partitions:
			os.remove("partitions.txt")
		os.remove(input_file)
		all_files = os.listdir(os.getcwd())
		for each in all_files:
			if re.search("(RAxML)", each):
				os.remove(each)
			if "In.phylip.reduced" == each:
				os.remove(each)
		return tree
	else:
		raise RuntimeError("Either phylogeny building program failed, or ran out of time")

if __name__ == '__main__':
	pass