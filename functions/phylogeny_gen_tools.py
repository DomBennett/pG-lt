import sys, os, re, random
import dendropy as dp
from Bio import Phylo
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
import numpy as np
sys.path.append(os.path.join(os.getcwd(), 'functions'))
import phyloGenerator_adapted as pG

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

def genConstraintTree(alignment, taxontree_file):
	tip_names = []
	for record in alignment:
		tip_names.append(record.id)
	constrainttree = Phylo.read(taxontree_file, "newick")
	taxontree_tips = []
	for terminal in constrainttree.get_terminals():
		taxontree_tips.append(terminal.name)
	print taxontree_tips
	print tip_names
	tips_to_drop = [e for e in taxontree_tips if not e in tip_names]
	for tip in tips_to_drop:
		constrainttree.prune(tip)
	print pG.Phylo.draw_ascii(constrainttree)
	Phylo.write(constrainttree, "constraint.txt", "newick")
	if constrainttree.is_bifurcating():
		return " -r constraint.txt"
	else:
		return " -g constraint.txt"

def concatenateAlignments(alignments):
	if len(alignments) == 1:
		return alignments,False
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

def RAxML(alignment, method='localVersion', tempStem='temp', outgroup=None, timeout=999999999,\
		partitions=None, constraint=None, cleanup=True):
	inputFile = tempStem + 'In.phylip'
	outputFile = tempStem + 'Out'
	fileLine = ' -s ' + inputFile + ' -n ' + outputFile
	options = ' -p ' + str(random.randint(0,10000000))
	if outgroup:
		options += ' -o ' + outgroup
	pG.AlignIO.write(alignment, inputFile, "phylip-relaxed")
	if 'localVersion' in method:
		raxmlVersion = 'raxml'
	else:
		raxmlVersion = 'raxmlHPC'
	#DNA model
	if len(alignment) > 100:
		DNAmodel = ' -m GTRCAT'
	else:
		DNAmodel = ' -m GTRGAMMA'
	#Partitions
	if partitions:
		with open(tempStem+"_partitions.txt", 'w') as f:
			for i in range(0, len(partitions)-1):
				f.write("DNA, position" + str(partitions[i]) + " = " + str(partitions[i]+1) + "-" + str(partitions[i+1]) + "\n")
		options += " -q " + tempStem + "_partitions.txt"
	#Constraint
	if constraint:
		options += constraint
	commandLine = raxmlVersion + fileLine + DNAmodel + options
	print commandLine
	pipe = pG.TerminationPipe(commandLine, timeout)
	if 'localVersion' in method:
			pipe.run(changeDir = True)
	else:
			pipe.run()
	if not pipe.failure:
		tree = Phylo.read('RAxML_bestTree.' + outputFile, "newick")
		if cleanup:
			if constraint:
				os.remove('constraint.txt')
			if partitions:
				os.remove(tempStem+"_partitions.txt")
			os.remove(inputFile)
			dirList = os.listdir(os.getcwd())
			for each in dirList:
				if re.search("(RAxML)", each):
					os.remove(each)
				if tempStem+"In.phylip.reduced"==each:
					os.remove(each)
		return tree
	else:
		raise RuntimeError("Either phylogeny building program failed, or ran out of time")


