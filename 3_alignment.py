#!/usr/bin/python
## MPE Stage 3: Aligning sequences
## D.J. Bennett
## 24/03/2014

## Print stage
print "\n\nStage 3: alignment\n"

## Packages
import os, re, pickle
from Bio import SeqIO
from alignment_tools import *

## Dirs
download_dir = os.path.join(os.getcwd(),'2_download')
alignment_dir = os.path.join(os.getcwd(),'3_alignment')
if not os.path.isdir(alignment_dir):
	os.mkdir(alignment_dir)

## Input
with open("genedict.p", "rb") as file:
	genedict = pickle.load(file)
with open("paradict.p", "rb") as file:
	paradict = pickle.load(file)

## Parameters
naligns = int(paradict["naligns"])
aligncounter = 0

## Process
genes = sorted(os.listdir(download_dir))
genes = [e for e in genes if not re.search("^\.|^log\.txt$", e)]
print 'Reading in sequences'
genestore = []
for gene in genes:
	gene_dir = os.path.join(download_dir, gene)
	seq_files = os.listdir(gene_dir)
	seqstore = SeqStore(gene_dir, seq_files, minfails = int(genedict[gene]["minfails"]),\
		mingaps = float(genedict[gene]["mingaps"]), minoverlap = int(genedict[gene]["minoverlap"]))
	genestore.append((gene, seqstore))
print "Running alignments"
for gene,seqstore in genestore:
	gene_dir = os.path.join(alignment_dir, gene)
	if not os.path.isdir(gene_dir):
		os.mkdir(gene_dir)
	print "Aligning gene [{0}] for [{1}] species ...".format(gene, len(seqstore))
	mingaps = float(genedict[gene]["mingaps"])
	minoverlap = int(genedict[gene]["minoverlap"])
	minfails = int(genedict[gene]["minfails"])
	maxtrys = int(genedict[gene]["maxtrys"])
	minseedsize = int(genedict[gene]["minseedsize"])
	maxseedtrys = int(genedict[gene]["maxseedtrys"])
	aligner = Aligner(seqstore, mingaps, minoverlap, minseedsize, maxtrys,\
		maxseedtrys)
	alignments = []
	trys = 0
	i = 1
	try:
		while i <= naligns:
			print "...iteration [{0}]".format(i)
			alignment = aligner.run()
			if alignment is None:
				trys += 1
				if trys > maxtrys:
					print "Max trys with no alignments hit!"
					break
				else:
					continue
			print "... alignment length [{0}] for [{1}] species".\
				format(alignment.get_alignment_length(), len(alignment))
			align_len = alignment.get_alignment_length()
			output_file = "{0}_nspp{1}_len{2}.faa".format(i,len(alignment),align_len)
			output_path = os.path.join(gene_dir, output_file)
			with open(output_path, "w") as file:
				count = SeqIO.write(alignment, file, "fasta")
			aligncounter += 1
			trys = 0
			i += 1
	except OutgroupError:
		print "... outgroup dropped"
		continue
	except MinSpeciesError:
		print "... too few species left in sequence pool"
		continue
	if aligncounter == 0:
		print "... no alignments generated"
		continue
print 'Stage finished. Generated [{0}] alignments.'.format(aligncounter)