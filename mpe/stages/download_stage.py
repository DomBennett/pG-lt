#! /usr/bin/env python
## D.J. Bennett
## 24/03/2014
"""
mpe Stage 2: Sequence Download
"""

## Packages
import os, pickle,logging
import mpe.tools.download_tools as dtools

def run(wd = os.getcwd()):
	## Print stage
	logging.info("Stage 2: Sequence download")

	## Dirs
	download_dir = os.path.join(wd, '2_download')
	if not os.path.isdir(download_dir):
		os.mkdir(download_dir)

	## Input
	with open(os.path.join(wd, ".genedict.p"), "rb") as file:
		genedict = pickle.load(file)
	with open(os.path.join(wd, ".paradict.p"), "rb") as file:
		paradict = pickle.load(file)
	with open(os.path.join(wd, ".namesdict.p"), "rb") as file:
		namesdict = pickle.load(file)
	with open(os.path.join(wd, ".allrankids.p"), "rb") as file:
		allrankids = pickle.load(file)

	## Parameters
	dtools.etools.Entrez.email = paradict["email"]
	nseqs = int(paradict['nseqs'])
	thoroughness = int(paradict['thoroughness'])
	seedsize = 10
	mingaps = 5
	minoverlap = 300
	maxtrys = 100
	minnseq = 1
	minnspp = 1
	target = 'all'
	maxpn = 0.1
	seqcounter = basecounter = 0

	## Process
	logging.info('Determining best genes ....')
	genes = dtools.findBestGenes(namesdict, genedict, thoroughness,\
		allrankids, minnseq, target, minnspp)
	statement = 'Using genes:'
	for gene in genes:
		statement += " [" + gene + "]"
	logging.info(statement)
	# Add genes to namesdict
	for key in namesdict.keys():
		namesdict[key]['genes'] = 0
	for gene in genes:
		seqcounter_gene = noseqcounter_gene = spcounter_gene = 0
		gene_names = genedict[gene]["names"]
		minlen = int(genedict[gene]["minlen"])
		maxlen = int(genedict[gene]["maxlen"])
		logging.info('Downloading and outputting for [{0}] ....'.\
			format(gene))
		gene_dir  = os.path.join(download_dir, str(gene))
		if not os.path.isdir(gene_dir):
			os.mkdir(gene_dir)
		for name in namesdict.keys():
			logging.info("..... [{0}]".format(name))
			taxids = namesdict[name]["txids"]
			downloader = dtools.Downloader(gene_names, nseqs, thoroughness,\
				maxpn, seedsize, maxtrys, mingaps, minoverlap, maxlen, minlen)
			sequences = downloader.run(taxids)
			if not sequences:
				noseqcounter_gene += 1
				logging.info("........ no sequences found")
				continue
			logging.info("........ downloaded [{0}] sequences".\
				format(len(sequences)))
			gene_seqs = []
			for seq in sequences:
				gene_seqs.append(seq.format('fasta'))
			with open(os.path.join(gene_dir, "{0}.fasta".format(name)), 'wb') \
			as outfile:
				for gene_seq in gene_seqs:
					outfile.write("%s\n" % gene_seq)
					seqcounter_gene += 1
					basecounter += len(gene_seq)
			namesdict[name]['genes'] += 1
			spcounter_gene += 1
		if noseqcounter_gene == len(namesdict.keys()):
			logging.info("No sequences were downloaded for gene [{0}]".format(gene))
			os.rmdir(gene_dir)
		else:
			seqcounter += seqcounter_gene
			logging.info("Downloaded [{0}] sequences for gene [{1}] representing \
[{2}] species".format(seqcounter_gene, gene, spcounter_gene))
	with open(os.path.join(wd, ".namesdict.p"), "wb") as file:
		pickle.dump(namesdict, file)
	logging.info('Stage finished. Downloaded [{0}] bases for [{1}] \
sequences for [{2}] species.'.format(basecounter, seqcounter,\
sum([namesdict[e]['genes'] > 0 for e in namesdict.keys()])))