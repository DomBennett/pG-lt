#! /bin/usr/env python
## D.J. Bennett
## 24/03/2014
"""
mpe names tools
"""

## Packages
import collections,re,logging, csv, os
from Bio import Phylo
from cStringIO import StringIO
import entrez as etools

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
	lineages = resolver.retrieve('classification_path_ids')
	# replace ' ' with '_' for taxon tree
	qnames = [re.sub("\s", "_", e) for e in qnames]
	resolved_names_bool = [e in namesdict.keys() for e in qnames]
	ranks = [ranks[ei] for ei,e in enumerate(resolved_names_bool)\
	if e]
	lineages = [lineages[ei] for ei,e in enumerate(resolved_names_bool)\
	if e]
	# identify unresolved names
	unresolved_names = [qnames[ei] for ei,e in\
	enumerate(resolved_names_bool) if not e]
	idents = [qnames[ei] for ei,e in enumerate(resolved_names_bool) if e]
	statement = "Unresolved names: "
	for each in unresolved_names:
		statement += " " + each
	logging.debug(statement)
	# reverse ranks and lineages
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
			new_node = [each[0] for each in line_obj if each[1][i]\
			== uniq]
			if len(new_node) > 1:
				# extract shared lineage
				lineage = [each[1] for each in line_obj if each[0]\
				== new_node[0]]
				# remove shareds from line_obj
				line_obj = [each for each in line_obj if not each[0]\
				in new_node]
				# convert to strings
				new_node = [str(each) for each in new_node]
				new_node = [re.sub("\\s", "_", each) for each\
				in new_node]
				# add new node to line_obj
				new_node = ('(' + ','.join(new_node) + ')' +\
					str(uniq), lineage[0])
				line_obj.append(new_node)
		if len(line_obj) < 1:
			break
	tree = Phylo.read(StringIO(line_obj[0][0] + ';'), "newick")
	if draw:
		Phylo.draw_ascii(tree)
	return tree, shared_lineage

def genNamesDict(resolver):
	"""Return a dictionary containtaining all names and metadata"""
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
		namesdict[q_names[i]] = {"txids" : [txid], "unique_name" : rname,\
		"rank" : rank}
	if non_unique_lineages:
		nul_ids = [lineages[e][-1] for e in non_unique_lineages]
		nul_qnames = [q_names[e] for e in non_unique_lineages]
		nul_ranks = [ranks[e][-1] for e in non_unique_lineages]
		i = 0
		while len(nul_ids) > 0:
			temp_id = nul_ids.pop(0)
			# find ids in the next level
			temp_children = etools.findChildren(str(temp_id), next = True)
			temp_children = [int(e) for e in temp_children]
			# if none are in allrankids, must be unique
			temp_children = [e for e in temp_children if e not in allrankids]
			if len(temp_children) > 0:
				namesdict[nul_qnames[i]] = {"txids" : temp_children,\
				"unique_name" : "Non-unique resolution",\
					"rank" : nul_ranks[i]}
			i += 1
	# if no parent id given, work one out
	if resolver.taxon_id:
		parentid = resolver.taxon_id
	else:
		shared_bool = []
		for each in lineages[0]:
			shared_bool.append(all([each in e for e in lineages]))
		parentid = lineages[0][shared_bool.index(False) - 1]
	return namesdict,allrankids,parentid

def getOutgroup(namesdict, parentid, minrecords = 1000):
	"""Return namesdict with suitable outgroup"""
	def findParent(parentid):
		return etools.eFetch(parentid, db = "taxonomy")[0]['ParentTaxId']
	# loop until a suitable outgroup is found. Criteria are:
	#  1. ids returned must belong to a sister group of all ids of
	#   names given
	#  2. ids must have nucleotide data (i.e. avoid returning extinct organisms)
	# assumptions:
	#  1. NCBI taxonomy is not paraphyletic
	# make sure parentid is string
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
		candidates = etools.findChildren(grandparentid, next = True)
		# filter out children that are in ingroup
		candidates = [e for e in candidates if e != parentid]
		# search genbank for nuc records
		for candidate in candidates:
			term = 'txid' + str(candidate) + '[PORGN]'
			nuc_record = etools.eSearch(term)
			# there must be more than 1000 nuc records
			if int(nuc_record['Count']) > minrecords:
				outgroup_ids.append(candidate)
		# make grandparentid the new parentid
		parentid = grandparentid
	# add outgroup_ids to namesdict
	namesdict["outgroup"] = {"txids" : outgroup_ids, "unique_name" :\
	"outgroup", "rank" : "outgroup"}
	return namesdict

def writeNamesDict(directory, namesdict):
	headers = ["name", "unique_name", "rank", "NCBI_Taxids"]
	with open(os.path.join(directory, 'resolved_names.csv'), 'wb')\
	as file:
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