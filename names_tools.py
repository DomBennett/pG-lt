#!/usr/bin/python
## MPE Names tools
## D.J. Bennett
## 24/03/2014

## Packages
import collections, re, copy
from Bio import Phylo
from cStringIO import StringIO
from entrez_tools import *

## Functions
def genTaxTree(resolver, namesdict, draw = False):
	"""Generate Newick tree from TaxonNamesResolver class.
		
		Arguments:
		 resolver = TaxonNamesResolver class
		 by = Tip labels, either 'qnames', 'taxids' or 'name_string'
		 draw = Draw ascii tree (logical)

		 Return:
		  (Newick Tree Object, [shared lineage])"""
	ranks = resolver.retrieve('classification_path_ranks')
	qnames = resolver.retrieve('query_name')
	idents = []
	namesdicts_ids = [e for e in namesdict.keys() if e != "outgroup"]
	for each in namesdicts_ids:
		idents.append(qnames.index(namesdict[each]["name"]))
	idents = [namesdicts_ids[e] for e in idents]
	lineages = resolver.retrieve('classification_path_ids')
	for i, lineage in enumerate(lineages):
		lineage.reverse()
		lineages[i] = lineage
	for i, rank in enumerate(ranks):
		rank.reverse()
		ranks[i] = rank
	# make lineages of same ranks
	all_ranks = [e2 for e1 in ranks for e2 in e1]
	rank_freq =  collections.Counter(all_ranks).items()
	shared_ranks = [e for e, f in rank_freq if f == len(idents)]
	line_bool = [[1 if e2 in shared_ranks else 0 for e2 in e1] for e1 in ranks]
	lineages = [[lineages[i1][i2] for i2, e2 in enumerate(e1) if e2 == 1] for \
					i1, e1 in enumerate(line_bool)]
	all_lines = [e2 for e1 in lineages for e2 in e1]
	line_freq =  collections.Counter(all_lines).items()
	shared_lineage = [e for e, f in line_freq if f == len(idents)]
	# create line_obj, a tuple of ident and lineage
	line_obj = zip(idents, lineages)
	for i in range(len(lineages[0])):
		for uniq in set([each[1][i] for each in line_obj]):
			# find shared taxonomic groups
			new_node = [each[0] for each in line_obj if each[1][i] == uniq]
			if len(new_node) > 1:
				# extract shared lineage
				lineage = [each[1] for each in line_obj if each[0] == new_node[0]]
				# remove shareds from line_obj
				line_obj = [each for each in line_obj if not each[0] in new_node]
				# convert to strings
				new_node = [str(each) for each in new_node]
				new_node = [re.sub("\\s", "_", each) for each in new_node]
				# add new node to line_obj
				new_node = ('(' + ','.join(new_node) + ')' + str(uniq), lineage[0])
				line_obj.append(new_node)
		if len(line_obj) < 1:
			break
	tree = Phylo.read(StringIO(line_obj[0][0] + ';'), "newick")
	if draw:
		Phylo.draw_ascii(tree)
	return tree, shared_lineage

def genNamesDict(resolver):
	q_names = resolver.retrieve('query_name')
	r_names = resolver.retrieve('classification_path')
	ranks = resolver.retrieve('classification_path_ranks')
	lineages = resolver.retrieve('classification_path_ids')
	lineages = [[int(e2) for e2 in e1] for e1 in lineages]
	all_ids = []
	[all_ids.extend(e) for e in lineages]
	all_ids = list(set(all_ids))
	namesdict = {}
	non_uniques = []
	for i in range(len(lineages)):
		query_line = lineages.pop(i)
		best_bool = [len(lineages)]
		j = 0
		while j < len(lineages):
			subj_line = lineages[j]
			match_bool = [e not in subj_line for e in query_line]
			if sum(match_bool) < sum(best_bool):
				best_bool = match_bool
			if not any(best_bool):
				non_uniques.append(i)
				break
			j += 1
		lineages.insert(i, query_line)
		if i in non_uniques:
			continue
		txid = query_line[best_bool.index(True)]
		rank = ranks[i][best_bool.index(True)]
		rname = r_names[i][best_bool.index(True)]
		namesdict["id{}".format(i)] = {"name" : q_names[i], "ids" : [txid],\
			"unique_name" : rname, "rank" : rank}
	if non_uniques:
		non_unique_ids = [lineages[e][-1] for e in non_uniques]
		non_unique_qnames = [q_names[e] for e in non_uniques]
		non_unique_ranks = [ranks[e][-1] for e in non_uniques]
		i = 0
		while len(non_unique_ids) > 0:
			temp_id = non_unique_ids.pop(i)
			# find ids in the next level
			temp_children = findChildren(str(temp_id), next = True)
			temp_children = [int(e) for e in temp_children]
			# if none are in all_ids, must be unique
			temp_children = [e for e in temp_children if e not in all_ids]
			if len(temp_children) == 0:
				continue # no unique children
			namesdict["id{0}".format(len(namesdict.keys()) + 1)] = {"name" : non_unique_qnames[i],\
				 "ids" : temp_children, "unique_name" : "Non-unique resolution", "rank" : non_unique_ranks[i]}
			i += 1
	# get outgroup
	if resolver.taxon_id:
		parentid = resolver.taxon_id
	else:
		shared_bool = []
		for each in lineages[0]:
			shared_bool.append(all([each in e for e in lineages]))
		parentid = lineages[0][shared_bool.index(False) - 1]
	above_id = eFetch(str(parentid), db = "taxonomy")[0]['ParentTaxId']
	temp_children = findChildren(above_id, next = True)
	temp_children = [int(e) for e in temp_children]
	# if none are in all_ids, must be unique
	temp_children = [e for e in temp_children if e != parentid]
	namesdict["outgroup"] = {"name" : "outgroup", "ids" : temp_children, "unique_name" : "outgroup",\
		"rank" : "outgroup"}
	return namesdict,all_ids

if __name__ == '__main__':
	pass