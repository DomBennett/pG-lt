#!/usr/bin/python
## MRes Project 2013
## Stage 3: Generate alignments
## In: 0_names, 2_download | Out: 3_alignments
## 14/08/2013

## Parameters
minfails = 20 # the minimum sequence quality
max_pgap = 0.5 # the proportion of gaps in a sequence for a good alignment
#min_align_len = 200 # minimum alignment length
iterations = 100 # number of iterations to perform
#max_attempts = 10 # the maximum number of failed in a row alignments



## Print stage
print "\n\nThis is stage 3: alignment\n"

## Packages
import sys, os, re, random, csv
import time # let's see if the iterations get faster
import numpy as np
sys.path.append(os.path.join(os.getcwd(), 'functions'))
import phyloGenerator_adapted as pG
from alignment_dev import *

## Dirs
input_dirs = [os.path.join(os.getcwd(), '2_download'), os.path.join(os.getcwd(), '0_names')]
output_dir = os.path.join(os.getcwd(), '3_alignments')
if not os.path.isdir(output_dir):
    os.mkdir(output_dir)

## Taxadata for identifying outgroup
print "Reading in taxadata.csv ..."
taxadict = {}
with open(os.path.join(input_dirs[1], 'taxadata.csv'), 'rb') as csvfile:
	taxreader = csv.DictReader(csvfile)
	for row in taxreader:
		taxadict[row['study']] = row['sisterID']
print "Done. Read in taxadata for [{0}] studies.".format(len(taxadict))

## Loop through studies
studies = sorted(os.listdir(input_dirs[0]))
studies = [st for st in studies if not re.search("^log\.txt$", st)]
counter = 0
print '\nLooping through studies ...'
naligns_all = 0
for i in range(len(studies)):

	## what study?
	print '\n\nWorking on: [{0}]\n'.format(studies[i])

	## determine what genes to use
	print 'Working out how many genes to use ...'
	study_dir = os.path.join(os.getcwd(), input_dirs[0], studies[i])
	genes = sorted(os.listdir(study_dir))
	print 'Done. Working with [{0}] ...'.format(genes)

	## read in seqs
	print 'Reading in sequences ...'
	outgroup = taxadict[studies[i]]
	seqs_obj = []
	for gene in genes:
		gene_dir = os.path.join(study_dir, gene)
		seq_files = os.listdir(gene_dir)
		seq_obj = SeqObj(gene_dir, seq_files, outgroup, minfails)
		if len(seq_obj) < 5 or seq_obj.nseqs < 5:
			continue
		else:
			seqs_obj.append((gene, seq_obj))
	if len(seqs_obj) == 0:
		print "Too little sequence data -- dropping study"
		continue
	nseqs = [e[1].nseqs for e in seqs_obj]
	ntaxa = [len(e[1]) for e in seqs_obj]
	if max(ntaxa) > 50: # to keep things fast for now
		print "Dropping large studies for now..."
		continue
	print 'Done. Read in [{0}] sequences for [{1}] gene(s) and between [{2}] to [{3}] species'\
		.format(sum(nseqs), len(seqs_obj), min(ntaxa), max(ntaxa))
		
	## alignment
	print "\nRunning alignments ..."
	all_alignments = []
        study_dir = os.path.join(output_dir, studies[i])
	if not os.path.isdir(study_dir):
		os.mkdir(study_dir)
	for gene,seq_obj in seqs_obj:
		print "Aligning gene [{0}] for [{1}] species ...".\
			format(gene, len(seq_obj))
		gene_alignments = []
                nstart = len(seq_obj)
                # Generate alignments
                for j in range(iterations):
                    print "iteration [{0}]".format(j)
                    t0 = time.clock()
                    alignment, nstart = incrAlign(seq_obj, max_pgap, nstart)
                    t1 = time.clock() - t0
                    if alignment is None:
                        continue
                    print "... alignment length [{0}] for [{1}] species in [t{2}]".\
                        format(alignment.get_alignment_length(), len(alignment), t1)
                    gene_alignments.append(alignment)
                # Write out alignments
                if len(gene_alignments) < 1:
                    print "... no alignments generated"
                    continue
                print "... writing out alignments for [{0}] alignments".\
                    format(len(gene_alignments))
                for j,alignment in enumerate(gene_alignments):
                    gene_dir = os.path.join(study_dir, gene)
                    if not os.path.isdir(gene_dir):
                        os.mkdir(gene_dir)
                    alength = alignment.get_alignment_length()
                    ngap = sum([e.seq.count("-") for e in alignment])
                    output_file = "a{0}_ngap{1}_length{2}.faa".format(j,ngap,alength)
                    output_path = os.path.join(gene_dir, output_file)
                    with open(output_path, "w") as outfile:
                        count = pG.SeqIO.write(alignment, outfile, "fasta")
                    naligns_all += 1   
        counter += 1

## Remove mafft files
mafft_files = os.listdir(os.getcwd())
mafft_files = [f for f in mafft_files if re.search("\.fasta$", f)]
for f in mafft_files:
	os.remove(f)
print 'Stage finished. Generated [{0}] alignments across [{1}] studies.'.\
    format(naligns_all, counter)
