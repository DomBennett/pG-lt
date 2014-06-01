#! /usr/bin/env python
## D.J. Bennett
## 24/03/2014
"""
MPE Stage 3: Aligning sequences
"""

## Packages
import os, re, pickle,logging
import numpy
from Bio import SeqIO
import mpe.tools.alignment as atools

def run():
	## print stage
	logging.info("\nStage 3: alignment\n")

	## Dirs
	download_dir = '2_download'
	alignment_dir = '3_alignment'
	if not os.path.isdir(alignment_dir):
		os.mkdir(alignment_dir)

	## Input
	with open(".genedict.p", "rb") as file:
		genedict = pickle.load(file)
	with open(".paradict.p", "rb") as file:
		paradict = pickle.load(file)
	with open(".namesdict.p", "rb") as file:
		namesdict = pickle.load(file)

	## Parameters
	naligns = int(paradict["naligns"])
	aligncounter = 0

	## Process
	# add alignments to namesdict
	for key in namesdict.keys():
		namesdict[key]['alignments'] = 0
	genes = sorted(os.listdir(download_dir))
	genes = [e for e in genes if not re.search("^\.|^log\.txt$", e)]
	logging.info('Reading in sequences')
	genestore = []
	for gene in genes:
		gene_dir = os.path.join(download_dir, gene)
		seq_files = os.listdir(gene_dir)
		seqstore = atools.SeqStore(gene_dir, seq_files, minfails = int(genedict[gene]["minfails"]),\
			mingaps = float(genedict[gene]["mingaps"]), minoverlap = int(genedict[gene]["minoverlap"]))
		genestore.append((gene, seqstore))
	logging.info("Running alignments")
	for gene,seqstore in genestore:
		gene_dir = os.path.join(alignment_dir, gene)
		if not os.path.isdir(gene_dir):
			os.mkdir(gene_dir)
		logging.info("Aligning gene [{0}] for [{1}] species ...".format(gene, len(seqstore)))
		mingaps = float(genedict[gene]["mingaps"])
		minoverlap = int(genedict[gene]["minoverlap"])
		minfails = int(genedict[gene]["minfails"])
		maxtrys = int(genedict[gene]["maxtrys"])
		minseedsize = int(genedict[gene]["minseedsize"])
		maxseedtrys = int(genedict[gene]["maxseedtrys"])
		aligner = atools.Aligner(seqstore, mingaps, minoverlap, minseedsize, maxtrys,\
			maxseedtrys)
		trys = 0
		i = 1
		try:
			while i <= naligns:
				logging.info("...iteration [{0}]".format(i))
				alignment = aligner.run()
				if alignment is None:
					trys += 1
					if trys > maxtrys:
						logging.info("Max trys with no alignments hit!")
						break
					else:
						continue
				logging.info("... alignment length [{0}] for [{1}] species".\
					format(alignment.get_alignment_length(), len(alignment)))
				for record in alignment:
					namesdict[record.id]['alignments'] += 1
				align_len = alignment.get_alignment_length()
				output_file = "{0}_nspp{1}_len{2}.faa".format(i,len(alignment),align_len)
				output_path = os.path.join(gene_dir, output_file)
				with open(output_path, "w") as file:
					count = SeqIO.write(alignment, file, "fasta")
					del count
				aligncounter += 1
				trys = 0
				i += 1
		except atools.OutgroupError:
			logging.info("... outgroup dropped")
		except atools.MinSpeciesError:
			logging.info("... too few species left in sequence pool")
		if aligncounter < naligns:
			logging.info("... too few alignments generated")
			os.rmdir(gene_dir)
			continue
	naligns_name = [namesdict[e]['alignments'] for e in namesdict.keys() if namesdict[e]['genes'] > 0]
	paligns_name = [len(naligns_name) * float(e)/aligncounter for e in naligns_name]
	with open(".namesdict.p", "wb") as file:
		pickle.dump(namesdict, file)
	logging.info('Stage finished. Generated [{n1}] alignments for mean [{n2:.{d}f}](sd[{n3:.{d}f}]) species.'.\
		format(n1 = aligncounter, n2 = numpy.mean(paligns_name), n3 = numpy.std(paligns_name), d = 2))