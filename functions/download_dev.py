#!/usr/bin/python
## Dominic John Bennett
## 10/12/2013
## Sequence download functions, main credit phyloGenerator (W.D. Pearse)
## TODO: look into using pG 'thorough' search
## TODO: instead of my thorough search, perhaps find the RNA field?
## TODO: Add doctest

## Packages
from Bio import Entrez #Taxonomy lookup
from Bio.Seq import Seq #Sequence manipulation
from Bio.SeqRecord import SeqRecord #Sequence manipulation
from Bio.SeqFeature import SeqFeature, FeatureLocation #Sequence manipulation
from Bio import SeqIO #Sequence manipulation
import random #Randomly select from multiple matches when downloading sequences
import re #Search for files to delete
import time #For waiting between sequence downloads

## Globals
maxCheck = 4

## Functions
def eSearch(term, retStart=0, retMax=1, usehistory="n"):
	"""Use Entrez.esearch to check for results in the nucleotide database given term.

	Arguments:
	 term = string of term used in search
	 retStart = minimum returned ID of matching sequences IDs
	 retMax = maximum returned ID of matching sequences IDs
	 usehistory = record search in NCBI database ("y" or "n")
	 
	Return:
	 dictionary"""
	finished = 0
	while finished <= maxCheck:
		try:
			handle = Entrez.esearch(db="nucleotide",term=term, usehistory=usehistory,\
							retStart=retStart, retMax=retMax, retmode="text")
			results = Entrez.read(handle)
			handle.close()
			finished = maxCheck + 1
		except:
			if finished == 0:
				print "!!!Server error checking", term, " - retrying..."
				time.sleep(10)
			elif finished == maxCheck:
				print "!!!!!!Unreachable. Returning nothing."
				return()
			else:
				finished += 1
				time.sleep(10)
	return results

def eFetchSeqID(seq_id):
	"""Download sequence using sequence ID number.

	Arguments:
	 seqID = sequence identifier

	Return:
	 SeqRecord object"""
	finished = 0
	while finished <= maxCheck:
		try:
			handle = Entrez.efetch(db = "nucleotide", rettype = 'gb',\
						       retmode="text", id = seq_id)
			results_iter = SeqIO.parse(handle, 'gb')
			results = [x for x in results_iter]
			if len(results) == 1:
				results = results[0]
			handle.close()
			finished = maxCheck + 1
		except:
			if finished == 0:
				print "!!!Server error - retrying..."
				finished += 1
				time.sleep(10)
			if finished == maxCheck:
				print "!!!!!!Unreachable. Returning nothing."
				return(tuple())
	return results

def findGeneInSeq(record, gene_names):
	"""Extract gene sequence from larger sequence by searching sequence features.

	Arguments:
	 record = SeqRecord object
	 gene_names = list of synonyms for gene of interest

	Return:
	 SeqRecord object or False (if feature not found)"""
	if record.features:
		for feature in record.features:
			feature_names = []
			if 'gene' in feature.qualifiers.keys():
				feature_names.extend(feature.qualifiers['gene'])
			if 'gene_synonym' in feature.qualifiers.keys():
				feature_names.extend(feature.qualifiers['gene_synonym'])
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

def sequenceDownload(txid, gene_names, nseqs = 100, minlen = 400, maxlen = 2000,\
			     maxpn = 0.1, thorough = True):
	"""Download all sequences for gene of interest for given taxon id.

	Arguments:
	 txid = taxon id
	 gene_names = list of synonyms for gene of interest
	 nseqs = max number of sequences to download (default 100)
	 minlen = min sequence length (default 400)
	 maxlen = max sequence length (default 2000)
	 maxpn = max proportion of ambiguous nucleotides (default 0.1)
	 thorough = logical, wider search if tight search pulls no results
	  (default True)

	Return:
	 List of SeqRecord objects"""
	def buildSearchTerm(gene_names, tight = True):
		# for field see: http://www.ncbi.nlm.nih.gov/books/NBK49540/
		if tight: # use gene field
			gene_term = ("{0}[GENE]".format(gene_names[0]))
			for gene_name in gene_names[1:]:
				gene_term = ("({0}) OR {1}[GENE]".\
						     format(gene_term, gene_name))
			search_term = ("txid{0}[PORGN] AND ({1})".\
					       format(txid, gene_term))
		else: # use all fields for gene name
			gene_term = ("{0}".format(gene_names[0]))
			for gene_name in gene_names[1:]:
				gene_term = ("({0}) OR {1}".\
						     format(gene_term, gene_name))
			search_term = ("txid{0}[PORGN] AND ({1})".\
					       format(txid, gene_term))
		return search_term
		
	def search(gene_names):
		seq_ids = []
		search_term = buildSearchTerm(gene_names, tight = True)
		count = eSearch(search_term)['Count']
		if int(count) < 1:
			if thorough:
				search_term = buildSearchTerm(gene_names, tight = False)
				count = eSearch(search_term)['Count']
				if int(count) < 1:
					return ()
			else:
				return ()
		seq_ids.extend(eSearch(search_term, retMax = count)['IdList'])
		return list(set(seq_ids))
	def parse(record):
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
	pattern = re.compile("[ACTGactg]")
	seq_store = search(gene_names)
	records = []
	i = 1
	while i <= nseqs and len(seq_store) > 0:
		randi = random.randint(0, len(seq_store)-1)
		seq = seq_store.pop(randi)
		record = eFetchSeqID(seq)
		record = parse(record)
		if record:
			records.append(record)
			i += 1
	return records

if __name__ == "__main__":
	#doctest to come
	pass
