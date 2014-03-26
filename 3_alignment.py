#!/usr/bin/python
## MPE Stage 3: Aligning sequences
## D.J. Bennett
## 24/03/2014

## Print stage
print "\n\nThis is stage 3: alignment\n"

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
#naligns = int(paradict["naligns"])
naligns = 10
aligncounter = 0

## Process
genes = sorted(os.listdir(download_dir))
genes = [e for e in genes if not re.search("^\.|^log\.txt$", e)]
print 'Reading in sequences'
geneobj = []
for gene in genes:
	gene_dir = os.path.join(download_dir, gene)
	seq_files = os.listdir(gene_dir)
	seqobj = SeqObj(gene_dir, seq_files, minfails = int(genedict[gene]["minfails"]))
	geneobj.append((gene, seqobj))
print "Running alignments"
for gene,seqobj in geneobj:
	print "Aligning gene [{0}] for [{1}] species ...".format(gene, len(seqobj))
	mingaps = float(genedict[gene]["mingaps"])
	minoverlap = int(genedict[gene]["minoverlap"])
	minfails = int(genedict[gene]["minfails"])
	maxtrys = int(genedict[gene]["maxtrys"])
	minseedsize = int(genedict[gene]["minseedsize"])
	maxseedtrys = int(genedict[gene]["maxseedtrys"])
	alignments = []
	seedsize = len(seqobj)
	trys = 0
	i = 1
	try:
		while i <= naligns:
			print "iteration [{0}]".format(i)
			alignment, seedsize = incrAlign(seqobj, mingaps, minoverlap, seedsize,\
				minseedsize, maxseedtrys)
			if alignment is None:
				trys += 1
				if trys > maxtrys:
					print "Max trys with no alignments hit!"
					break
				else:
					continue
			print "... alignment length [{0}] for [{1}] species".\
				format(alignment.get_alignment_length(), len(alignment))
			alignments.append(alignment)
			trys = 0
			i += 1
	except OutgroupError:
		print "... outgroup dropped"
		continue
	except MinSpeciesError:
		print "... too few species left in sequence pool"
		continue
	if len(alignments) < 1:
		print "... no alignments generated"
		continue
	print "... writing out alignments for [{0}] alignments".\
		format(len(alignments))
	for i,alignment in enumerate(alignments):
		gene_dir = os.path.join(alignment_dir, gene)
		if not os.path.isdir(gene_dir):
			os.mkdir(gene_dir)
		align_len = alignment.get_alignment_length()
		output_file = "{0}_nspp{1}_len{2}.faa".format(i,align_len,len(alignment))
		output_path = os.path.join(gene_dir, output_file)
		with open(output_path, "w") as file:
			count = SeqIO.write(alignment, file, "fasta")
		aligncounter += 1
print 'Stage finished. Generated [{0}] alignments.'.format(aligncounter)