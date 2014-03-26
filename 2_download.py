#!/usr/bin/python
## MPE Stage 2: Sequence Download
## D.J. Bennett
## 24/03/2014

## Print stage
print "\n\nThis is stage 2: sequence download\n"

## Packages
import sys, os, re, random, pickle
from download_tools import *
from entrez_tools import *

## Dirs
download_dir = os.path.join(os.getcwd(),'2_download')
if not os.path.isdir(download_dir):
	os.mkdir(download_dir)

## Input
with open("genedict.p", "rb") as file:
	genedict = pickle.load(file)
with open("paradict.p", "rb") as file:
	paradict = pickle.load(file)
with open("namesdict.p", "rb") as file:
	namesdict = pickle.load(file)
with open("all_ids.p", "rb") as file:
	allrankids = pickle.load(file)

## Parameters
Entrez.email = paradict["email"]
nseqs = int(paradict['nseqs'])
thoroughness = int(paradict['download_thoroughness'])
maxlen = int(paradict['maxlen'])
minlen = int(paradict['minlen'])
filter_threshold = int(paradict['filter_threshold'])
filter_seed = 5
mingaps = 0.01
minoverlap = 200
maxtrys = 100
seqcounter = basecounter = spcounter = 0

## Process
print 'Determining best genes'
#genes = findBestGenes(namesdict, genedict, thoroughness, allrankids)
genes = ["COI", "18S"]
statement = 'Using genes:'
for gene in genes:
	statement += " " + gene
print statement
for gene in genes:
	seqcounter_gene = noseqcounter_gene = spcounter_gene = 0
	gene_names = genedict[gene]["names"]
	print 'Downloading and outputting for [{0}]...'.format(gene)
	gene_dir  = os.path.join(download_dir, str(gene))
	if not os.path.isdir(gene_dir):
		os.mkdir(gene_dir)
	for nameid in namesdict.keys():
		print "..... [{0}]: [{1}]".format(nameid,namesdict[nameid]["name"])
		taxids = namesdict[nameid]["ids"]
		sequences = []
		seqids = sequenceSearch(taxids, gene_names, thoroughness)
		if len(seqids) >= nseqs: # Only filter if more than 100 in genbank
			print "Lots of sequences for [{0}]. Downloading and filtering.".format(nameid)
			for taxid in taxids:
				temp_taxids = findChildren(taxid)
				downloaded = []
				deja_vues = []
				while len(sequences) < nseqs:
					downloaded,temp_deja_vues = sequenceDownload(temp_taxids, gene_names, deja_vues, minlen, maxlen,\
						maxpn = 0.01, thoroughness = thoroughness, nseqs = nseqs, seq_ids = seqids)
					if len(downloaded) == 0: # keeps looping until no more new sequences are being downloaded
						break
					deja_vues.extend(temp_deja_vues)
					deja_vues = list(set(deja_vues))
					if len(downloaded) < filter_seed and len(temp_taxids) < 1: # only skip the filtering stage if few temp_taxids
						sequences.extend(downloaded)
						break
					else:
						filtered,downloaded = filterSequences(sequences = downloaded, filter_seed = filter_seed,\
							mingaps = mingaps, minoverlap = minoverlap, maxtrys = maxtrys, minlen = minlen)
					if len(filtered) > 0:
						sequences.extend(filtered)
		else:
			sequences,_ = sequenceDownload(taxids, gene_names, [], minlen, maxlen, maxpn = 0.01,\
				thoroughness = thoroughness, nseqs = nseqs, seq_ids = seqids)
		if len(sequences) < 1:
			noseqcounter_gene += 1
			print "No sequences found for [id{0}].".format(nameid)
			continue
		if len(sequences) > nseqs:
			sequences = random.sample(sequences, nseqs)
		gene_seqs = []
		for seq in sequences:
			gene_seqs.append(seq.format('fasta'))
		with open(os.path.join(gene_dir, str(nameid) + '.fasta'), 'wb') \
		as outfile:
			for gene_seq in gene_seqs:
				outfile.write("%s\n" % gene_seq)
				seqcounter_gene += 1
				basecounter += len(gene_seq)
		spcounter_gene += 1
	if noseqcounter_gene == len(namesdict.keys()):
		print "No sequences were downloaded for gene [{0}]".format(gene)
		os.rmdir(gene_dir)
	else:
		seqcounter += seqcounter_gene
		spcounter += spcounter_gene
		print "Downloaded [{0}] sequences for gene [{1}] representing [{2}] species".\
			format(seqcounter_gene,gene,spcounter_gene)
print '\n\nStage finished. Downloaded [{0}] bases for [{1}] sequences for [{2}] species.'.\
	format(basecounter, seqcounter, spcounter)