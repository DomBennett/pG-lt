#!/usr/bin/python
## MPE Download tools
## D.J. Bennett
## 24/03/2014

## Packages
import re,random
from entrez_tools import *
from alignment_tools import *
from Bio.SeqFeature import SeqFeature, FeatureLocation

## Functions
def filterSequences(sequences, filter_seed, maxtrys, mingaps, minoverlap, minlen):
	if len(sequences) < filter_seed:
		return [], sequences
	def filterByAlignment(sequences):
		align = alignSequences(sequences)
		return alignCheck(align, mingaps, minoverlap, minlen)
	filtered = []
	trys = 0
	while filter_seed > 2 and trys < maxtrys:
		temp = []
		for i in range(filter_seed):
			randn = random.randint(0, len(sequences)-1)
			temp.append(sequences.pop(randn))
		if filterByAlignment(temp):
			filtered.extend(temp)
		else:
			sequences.extend(temp)
			trys += 1
		if len(sequences) < filter_seed*2:
			filter_seed -= 1
		if len(sequences) < filter_seed:
			filter_seed = len(sequences) - 1
	return filtered,sequences

def findGeneInSeq(record, gene_names):
	"""Adapted pG function. Extract gene sequence from larger sequence by searching sequence features.

	Arguments:
	 record = SeqRecord object
	 gene_names = list of synonyms for gene of interest

	Return:
	 SeqRecord object or False (if feature not found)"""
	try:
		if record.features:
			for feature in record.features:
				feature_names = []
				if 'gene' in feature.qualifiers.keys():
					feature_names.extend(feature.qualifiers['gene'])
				if 'gene_synonym' in feature.qualifiers.keys():
					feature_names.extend(feature.qualifiers['gene_synonym'])
				if 'product' in feature.qualifiers.keys():
					feature_names.extend(feature.qualifiers['product'])
				gene_names = [e.lower() for e in gene_names]
				feature_names = [e.lower() for e in feature_names]
				if set(gene_names) & set(feature_names):
					extractor = SeqFeature(feature.location)
					found_seq = extractor.extract(record)
					return found_seq # dropped TrimSeq (as I didn't know what it did...)
			else:
				print "... findGeneInSeq: can't find gene [{0}] in sequence [{1}]".\
					format(gene_names[0], record.id)
				return ()
		else:
			print "... findGeneInSeq: can't find features in in sequence [{1}]".\
				format(record.id)
			return ()
	except ValueError: # catch value errors raised for sequences with "fuzzy" positions
		return ()

def findBestGenes(namesdict, genedict, thoroughness, allrankids, minnseq = 1, minpwithseq = 0.5):
	alltipids = [namesdict[e]["ids"] for e in namesdict.keys()]
	genes = []
	for gene in genedict.keys():
		print "Checking [{0}]".format(gene)
		taxid = genedict[gene]["taxid"]
		if int(taxid) in allrankids:
			gene_bool = []
			for tipids in alltipids:
				res = sequenceSearch(tipids, genedict[gene]["names"], thoroughness)
				gene_bool.append(len(res) >= minnseq)
			pwithseq = float(sum(gene_bool))/len(alltipids)
			if pwithseq > minpwithseq:
				genes.append(gene)
	return genes

def sequenceSearch(taxids, gene_names, thoroughness, deja_vues = []):
	"""Search for sequences given a taxid and gene names.

	Arguments:
	 taxids = taxon ids (list)
	 gene_names = list of synonyms for gene of interest
	 thoroughness = a number between 1 and 3 determining how thorough
	  a search should be

	 Return:
	  List of Sequence IDs"""
	def buildSearchTerm(taxids, gene_names, phase):
		taxids_term = "txid{0}[PORGN]".format(taxids[0])
		for taxid in taxids[1:]:
			taxids_term = ("({0}) OR txid{1}[PORGN]".\
						 format(taxids_term, taxid))
		# for field see: http://www.ncbi.nlm.nih.gov/books/NBK49540/
		gene_names = ["\"{0}\"".format(e) for e in gene_names]
		if phase == 1: # use gene field and ignore: predicted, genome
			gene_term = ("{0}[GENE]".format(gene_names[0]))
			for gene_name in gene_names[1:]:
				gene_term = ("({0}) OR {1}[GENE]".\
							 format(gene_term, gene_name))
			search_term = ("{0} AND ({1}) NOT predicted[TI] NOT genome[TI] NOT unverified[TI]".\
						   format(taxids_term, gene_term))
		elif phase == 2: # use title for gene name and ignore: predicted, genome
			gene_term = ("{0}[TI]".format(gene_names[0]))
			for gene_name in gene_names[1:]:
				gene_term = ("({0}) OR {1}[TI]".\
							 format(gene_term, gene_name))
			search_term = ("{0} AND ({1}) NOT predicted[TI] NOT genome[TI] NOT unverified[TI]".\
						   format(taxids_term, gene_term))
		else: # all fields, ignore whole genome shotguns scaffolds and genome assemblies
			gene_term = ("{0}".format(gene_names[0]))
			for gene_name in gene_names[1:]:
				gene_term = ("({0}) OR {1}".\
							 format(gene_term, gene_name))
			search_term = ("{0} AND ({1}) NOT predicted[TI] NOT shotgun[TI] NOT scaffold[TI] NOT assembly[TI] NOT unverified[TI]".\
						   format(taxids_term, gene_term))
		return search_term
	seq_ids = []
	count = 0
	phase = 1
	while len(seq_ids) == 0:
		if phase > thoroughness:
			return []
		search_term = buildSearchTerm(taxids, gene_names, phase)
		seqcount = eSearch(search_term)['Count']
		if int(seqcount) >= 1:
			seq_ids.extend(eSearch(search_term, retMax = seqcount)['IdList'])
			seq_ids = [e for e in seq_ids if not e in deja_vues]
		phase += 1
	return list(set(seq_ids))

def sequenceDownload(taxids, gene_names, deja_vues, minlen, maxlen,\
				 maxpn, thoroughness, nseqs, seq_ids = None):
	"""Download and parse sequences given a taxid and gene names or list of sequence IDs.

	Arguments:
	 taxids = taxon ids (list)
	 gene_names = list of synonyms for gene of interest
	 minlen = min sequence length
	 maxlen = max sequence length
	 maxpn = max proportion of ambiguous nucleotides
	 thoroughness = a number between 1 and 3 determining how thorough
	  a search should be
	 nseqs = max number of sequences to download (default 100)
	 seq_ids = list of pre-specified sequence IDs

	Return:
	 List of SeqRecord objects"""
	def parse(record):
		if isinstance(record, list):
			# find which sequence in the list has the gene
			for each in record:
				for gene in gene_names:
					if gene in each.description.lower():
						record = each
						break
			else:
				return False
		if len(record) > maxlen:
			record = findGeneInSeq(record, gene_names)
		if maxlen > len(record) > minlen:
			# find proportion of ambiguous bases
			nn = len(pattern.findall(record.seq.tostring()))
			pn = float(len(record)) - float(nn)
			pn = pn/float(len(record))
			if pn < maxpn:
				return record
		return False
	if seq_ids is None:
		seq_store = sequenceSearch(taxids, gene_names, thoroughness, deja_vues)
		deja_vues = seq_store[:]
	else:
		seq_store = seq_ids[:]
	pattern = re.compile("[ACTGactg]")
	records = []
	i = 1
	while i <= nseqs and len(seq_store) > 0:
		if len(seq_store) > 100:
			n = 100 # Download in chunks of 100
		else:
			n = len(seq_store)
		seqs = []
		for _ in range(n):
			randi = random.randint(0, len(seq_store)-1)
			seqs.append(seq_store.pop(randi))
		for record in eFetch(seqs):
			print record.id
			record = parse(record)
			if record:
				records.append(record)
				i += 1
	return records,deja_vues

if __name__ == "__main__":
	pass