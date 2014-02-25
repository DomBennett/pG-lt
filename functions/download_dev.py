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
import phyloGenerator_adapted as pG
import sys, os

## Globals
maxCheck = 4
download_counter = 0 # Limit the number of requests to NCBI for a given time period

## Functions
def cleanAlignment(align, method='trimAl-automated', tempStem='temp', timeout=None):
	if 'trimAl' in method:
		options = ""
		if 'automated' in method:
			options += " -automated1"
		output = []
		for i,gene in enumerate(align):
			geneOutput = []
			for j,method in enumerate(gene):
				pG.AlignIO.write(method, tempStem+"Input.fasta", "fasta")
				fileLine = " -in " + tempStem + "Input.fasta -out " + tempStem + "Output.fasta -fasta"
				trimalVersion = "trimal"
				commandLine = trimalVersion + fileLine + options
				if timeout:
                                        pipe = pG.TerminationPipe(commandLine, timeout)
                                        if 'darwin' == sys.platform: # so that I can run it on mac
                                                pipe.run(changeDir = True)
                                        else:
                                                pipe.run()
					if not pipe.failure:
						geneOutput.append(pG.AlignIO.read(tempStem + "Output.fasta", "fasta"))
						os.remove(tempStem + "Output.fasta")
						os.remove(tempStem + "Input.fasta")
					else:
						raise RuntimeError("Either trimAl failed, or it ran out of time")
			output.append(geneOutput)
		return output
	else:
		raise RuntimeError("Only automated trimAl methods supported at this time.")

def filterSequences(sequences, filter_seed, pintgapmax, pextgapmax, max_trys, minlen):
	if len(sequences) < filter_seed:
		return [], sequences
	def filterByAlignment(sequences):
		sequences = [[e] for e in sequences]
		align = pG.alignSequences(sequences, method = "mafft", nGenes = 1)
		align = cleanAlignment(align, timeout = 99999)[0][0] #trim with trimal
		if align.get_alignment_length() < minlen:
				return False
		for each in align:
			sequence = each.seq.tostring()
			totgaps = sequence.count('-')
			totnucs = len(sequence) - totgaps
			if totnucs < minlen:
				return False
			extgaps = 0
			start_extgaps = re.search("^-+", sequence)
			end_extgaps = re.search("-+$", sequence)
			if not start_extgaps is None:
				extgaps += start_extgaps.end()
			if not end_extgaps is None:
				extgaps += len(sequence) - end_extgaps.start()
			intgaps = totgaps - extgaps
			overlap = intgaps + totnucs
			if float(intgaps)/overlap > pintgapmax:
				return False
			if float(extgaps)/len(sequence) > pextgapmax:
				return False
		return True
	filtered = []
	trys = 0
	while filter_seed > 2 and trys < max_trys:
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
	 seq_id = sequence identifier

	Return:
	 SeqRecord object"""
	finished = 0
	global download_counter
	print seq_id
	while finished <= maxCheck:
		if download_counter > 1000:
			print "Download counter hit. Waiting for 60 seconds ..."
			download_counter = 0
			time.sleep(60)
		try:
			handle = Entrez.efetch(db = "nucleotide", rettype = 'gb',\
						       retmode="text", id = seq_id)
			results_iter = SeqIO.parse(handle, 'gb')
			results = [x for x in results_iter]
			if len(results) == 1:
				results = results[0]
			handle.close()
			download_counter += 1
			finished = maxCheck + 1
		except ValueError: # if parsing fails, value error raised
			handle.close()
			results = ()
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

def sequenceDownload(txid, gene_names, deja_vues, minlen, maxlen,\
			     maxpn, thoroughness, nseqs = 100):
	"""Download all sequences for gene of interest for given taxon id.

	Arguments:
	 txid = taxon id
	 gene_names = list of synonyms for gene of interest
	 nseqs = max number of sequences to download (default 100)
	 minlen = min sequence length
	 maxlen = max sequence length
	 maxpn = max proportion of ambiguous nucleotides
	 thoroughness = a number between 1 and 3 determining how thorough
	  a search should be

	Return:
	 List of SeqRecord objects"""
	def buildSearchTerm(gene_names, phase):
		# for field see: http://www.ncbi.nlm.nih.gov/books/NBK49540/
		gene_names = ["\"{0}\"".format(e) for e in gene_names]
		if phase == 1: # use gene field and ignore: predicted, genome
			gene_term = ("{0}[GENE]".format(gene_names[0]))
			for gene_name in gene_names[1:]:
				gene_term = ("({0}) OR {1}[GENE]".\
						     format(gene_term, gene_name))
			search_term = ("txid{0}[PORGN] AND ({1}) NOT predicted[TI] NOT genome[TI] NOT unverified[TI]".\
					       format(txid, gene_term))
		elif phase == 2: # use title for gene name and ignore: predicted, genome
			gene_term = ("{0}[TI]".format(gene_names[0]))
			for gene_name in gene_names[1:]:
				gene_term = ("({0}) OR {1}[TI]".\
						     format(gene_term, gene_name))
			search_term = ("txid{0}[PORGN] AND ({1}) NOT predicted[TI] NOT genome[TI] NOT unverified[TI]".\
					       format(txid, gene_term))
		else: # all fields, ignore whole genome shotguns scaffolds and genome assemblies
			gene_term = ("{0}".format(gene_names[0]))
			for gene_name in gene_names[1:]:
				gene_term = ("({0}) OR {1}".\
						     format(gene_term, gene_name))
			search_term = ("txid{0}[PORGN] AND ({1}) NOT predicted[TI] NOT shotgun[TI] NOT scaffold[TI] NOT assembly[TI] NOT unverified[TI]".\
					       format(txid, gene_term))
		print search_term
		return search_term
	def search(gene_names, phase):
		seq_ids = []
		count = 0
		while int(count) < 1:
			if phase > thoroughness:
				return (), phase
			search_term = buildSearchTerm(gene_names, phase = phase)
			count = eSearch(search_term)['Count']
			phase += 1
		seq_ids.extend(eSearch(search_term, retMax = count)['IdList'])
		return list(set(seq_ids)),phase
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
	pattern = re.compile("[ACTGactg]")
	phase = 1
	while phase <= thoroughness:
		seq_store,phase = search(gene_names, phase)
		seq_store = [e for e in seq_store if not e in deja_vues]
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
		if len(records) == 0: # Search again at next phase if no suitable records
			phase += 1
		else:
			break
	return records,seq_store

if __name__ == "__main__":
	#doctest to come
	pass
