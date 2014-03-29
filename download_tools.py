#!/usr/bin/python
## MPE Download tools
## D.J. Bennett
## 24/03/2014

## Packages
import re,random
from entrez_tools import *
from alignment_tools import *
from Bio.SeqFeature import SeqFeature, FeatureLocation

## Objects
class SequenceDownloader(object):
	"""Download sequences given taxids and gene_names"""
	def __init__(self, gene_names, nseqs, thoroughness, maxpn, seedsize, maxtrys, \
			mingaps, minoverlap, maxlen, minlen):
		self.gene_names = gene_names
		self.nseqs = nseqs
		self.max_thoroughness = thoroughness
		self.maxpn = maxpn
		self.seedsize = seedsize
		self.maxtrys = maxtrys
		self.mingaps = mingaps
		self.minoverlap = minoverlap
		self.maxlen = maxlen
		self.minlen = minlen
		self.temp_taxids = None
		self.seqs_deja_vues = []
		self.taxids_deja_vues = []
		self.pattern = re.compile("[ACTGactg]")

	def _buildSearchTerm(self, taxids, thoroughness):
		"""Generate NCBI GenBank query given taxids, gene_names and thoroughness"""
		taxids_term = "txid{0}[PORGN]".format(taxids[0])
		for taxid in taxids[1:]:
			taxids_term = ("({0}) OR txid{1}[PORGN]".\
						 format(taxids_term, taxid))
		# for fields see: http://www.ncbi.nlm.nih.gov/books/NBK49540/
		gene_names = ["\"{0}\"".format(e) for e in self.gene_names]
		if thoroughness == 1: # use gene field and ignore: predicted, genome
			gene_term = ("{0}[GENE]".format(gene_names[0]))
			for gene_name in gene_names[1:]:
				gene_term = ("({0}) OR {1}[GENE]".\
							 format(gene_term, gene_name))
			search_term = ("{0} AND ({1}) NOT predicted[TI] NOT genome[TI] NOT unverified[TI]".\
						   format(taxids_term, gene_term))
		elif thoroughness == 2: # use title for gene name and ignore: predicted, genome
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
		print search_term
		return search_term

	def _search(self, taxids):
		"""Search GenBank for matches, increasing thoroughness if no matches"""
		seqids = []
		self.thoroughness = 1
		while len(seqids) == 0:
			if self.thoroughness > self.max_thoroughness:
				return None
			search_term = self._buildSearchTerm(taxids, self.thoroughness)
			seqcount = eSearch(search_term)['Count']
			if int(seqcount) >= 1:
				seqids.extend(eSearch(search_term, retMax = seqcount)['IdList'])
				seqids = [e for e in seqids if not e in self.seqs_deja_vues]
			self.thoroughness += 1
		self._DejaVue(seqids)
		return list(set(seqids))

	def _filterSequences(self, sequences):
		"""Filter sequences by aligning sequences from closely organisms"""
		if len(sequences) < self.seedsize:
			return sequences
		def filterByAlignment(sequences):
			alignment = alignSequences(sequences)
			return alignmentCheck(alignment, self.mingaps, self.minoverlap, self.minlen)
		filtered = []
		trys = 0
		seedsize = self.seedsize
		while seedsize > 2 and trys < self.maxtrys:
			temp = []
			for i in range(seedsize):
				randn = random.randint(0, len(sequences)-1)
				temp.append(sequences.pop(randn))
			if filterByAlignment(temp):
				filtered.extend(temp)
			else:
				sequences.extend(temp)
				trys += 1
			if len(sequences) < seedsize*2:
				seedsize -= 1
			if len(sequences) < seedsize:
				seedsize = len(sequences) - 1
		return filtered

	def _findGeneInSeq(self, record):
		"""Extract gene sequence from larger sequence by searching features."""
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
					gene_names = [e.lower() for e in self.gene_names]
					feature_names = [e.lower() for e in feature_names]
					if set(gene_names) & set(feature_names):
						extractor = SeqFeature(feature.location)
						found_seq = extractor.extract(record)
						return found_seq
				else:
					return None
			else:
				return None
		except ValueError: # catch value errors raised for sequences with "fuzzy" positions
			return None

	def _parse(self, record):
		"""Parse record returned from GenBank"""
		if isinstance(record, list):
			# find which sequence in the list has the gene
			for each in record:
				for gene in gene_names:
					if gene in each.description.lower():
						record = each
						break
			else:
				return None
		if len(record) > self.maxlen:
			record = self._findGeneInSeq(record)
		if self.maxlen > len(record) > self.minlen:
			# find proportion of ambiguous bases
			nn = len(self.pattern.findall(record.seq.tostring()))
			pn = float(len(record)) - float(nn)
			pn = pn/float(len(record))
			if pn < self.maxpn:
				return record
		return None

	def _download(self, seqids):
		"""Download records from genbank given sequence ids"""
		records = []
		i = 1
		while i <= self.nseqs and len(seqids) > 0:
			if len(seqids) > 100:
				n = 100 # Download in chunks of 100
			else:
				n = len(seqids)
			seqs = []
			for _ in range(n):
				randi = random.randint(0, len(seqids)-1)
				seqs.append(seqids.pop(randi))
			for record in eFetch(seqs):
				record = self._parse(record)
				if record:
					records.append(record)
					i += 1
		return records

	def _findGenera(self, taxids):
		"""Randomly return a taxid of decendant genera (or below). None if genera deja-vue."""
		if self.temp_taxids:
			temp_taxid = self.temp_taxids.pop(-1)
		else:
			taxid = random.sample(taxids, 1)[0]
			temp_taxids = findChildren(taxid, target = 1)
			temp_taxid = random.sample(temp_taxids, 1)[0]
		if temp_taxid in self.taxids_deja_vues:
			return None
		else:
			self._DejaVue(temp_taxid, sequences = False)
			return [temp_taxid]

	def _DejaVue(self, vues, sequences = True):
		"""Add new elements to lists of deja_vues"""
		if sequences:
			self.seqs_deja_vues.extend(vues)
			self.seqs_deja_vues = list(set(self.seqs_deja_vues))
		else:
			self.taxids_deja_vues.extend(vues)
			self.taxids_deja_vues = list(set(self.taxids_deja_vues))

	def main(self, taxids):
		"""Dynamic sequence download"""
		sequences = []
		seqids = self._search(taxids)
		if seqids:
			if len(seqids) >= 1000: # Only filter if more than 1000 mathces in genbank
				print "Lots of sequences -- filtering ...."
				while len(sequences) < self.nseqs:
					downloaded = []
					temp_taxid = self._findGenera(taxids)
					if not temp_taxid:
						break
					seqids = self._search(temp_taxid)
					print seqids
					if seqids:
						downloaded = self._download(seqids)
						if not downloaded:
							break
						filtered = self._filterSequences(downloaded)
						if filtered:
							sequences.extend(filtered)
				return sequences
			else:
				sequences = self._download(seqids)
				if sequences:
					return sequences
				for _ in range(self.thoroughness+1, self.max_thoroughness+1):
					seqids = self._search(taxids)
					sequences = self._download(seqids)
				return sequences
		else:
			return None

## Functions
def findBestGenes(namesdict, genedict, thoroughness, allrankids, minnseq = 1, minpwithseq = 0.5):
	"""Return suitable genes for phylogeny by searching for matches in GenBank"""
	alltipids = [namesdict[e]["txids"] for e in namesdict.keys()]
	genes = []
	for gene in genedict.keys():
		print "Checking [{0}]".format(gene)
		taxid = genedict[gene]["taxid"]
		if int(taxid) in allrankids:
			gene_bool = []
			for tipids in alltipids:
				sequence_downloader = SequenceDownloader(gene_names = genedict[gene]["names"],\
					nseqs = 0, thoroughness = thoroughness, maxpn = 0, seedsize = 0, maxtrys = 0,\
					mingaps = 0, minoverlap = 0, maxlen = 0, minlen = 0)
				res = sequence_downloader._search(tipids)
				gene_bool.append(len(res) >= minnseq)
			pwithseq = float(sum(gene_bool))/len(alltipids)
			if pwithseq > minpwithseq:
				genes.append(gene)
	return genes

if __name__ == "__main__":
	pass