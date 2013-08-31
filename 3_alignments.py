#!/usr/bin/python
## MRes Project 2013
## Stage 3: Generate alignments
## In: 0_names, 2_download | Out: 3_alignments
## 14/08/2013

## Parameters
# the sizes of the shrinking window of overlap
gene_overlaps = [0.5, 0.25, 0.1, 0.05]
# the sizes of subsamples of the species
props_nspp = [1.0, 0.90, 0.75, 0.5, 0.25, 0.1]
# the max size of the alignment with respect to the median
align_len_max = 0.5
# the minimum number of species for alignment
min_nspp = 5
# the number of fails in a row before dropping alignment attempt
nfails = 10

## Print stage
print "\n\nThis is stage 3: alignment\n"

## Packages
import sys, os, re, random, csv
import numpy as np
sys.path.append(os.path.join(os.getcwd(), 'functions'))
import phyloGenerator_adapted as pG

## dirs
input_dirs = [os.path.join(os.getcwd(), '2_download'), os.path.join(os.getcwd(), '0_names')]
output_dir = os.path.join(os.getcwd(), '3_alignments')
if not os.path.isdir(output_dir):
    os.mkdir(output_dir)
    
## Blacklist
print "Checking for a blacklist ..."
blacklist_dir = os.path.join(input_dirs[1], 'blacklist.csv')
if os.path.exists(blacklist_dir):
	blacklist = {}
	nbgenes = nbids = 0
	with open(blacklist_dir, 'rb') as csvfile:
		blreader = csv.DictReader(csvfile)
		for row in blreader:
			nbids += 1
			if row['study'] in blacklist.keys():
				if row['genes'] in blacklist[row['study']].keys():
					blacklist[row['study']][row['genes']].append(row['taxids'])
				else:
					blacklist[row['study']][row['genes']] = [row['taxids']]
			else:
				blacklist[row['study']] = {row['genes'] : [row['taxids']]}
				nbgenes += 1
		print "Done. Read in a blacklist for [{0}] genes and [{1}] ids.".\
			format(nbgenes,nbids)
else:
	print "No blacklist found ..."
	blacklist = False

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
if blacklist:
	studies = [st for st in studies if st in blacklist.keys()]
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
	if blacklist:
		blacklist_st = blacklist[studies[i]]
	print 'Done. Working with [{0}] ...'.format(genes)

	## read in seqs
	print 'Reading in sequences ...'
	outgroup = taxadict[studies[i]]
	seqs_obj = []
	nseqs = 0
	ntaxa = []
	for gene in genes:
		gene_dir = os.path.join(study_dir, gene)
		taxid_files = sorted(os.listdir(gene_dir))
		taxids = [re.sub('\.fasta$', '', t) for t in taxid_files]
		if blacklist:
			blackids = [e for e in taxids if str(e) in blacklist_st[gene]]
			taxids = [e for e in taxids if e not in blackids]
			print "Dropped [{0}] ids on blacklist for gene [{1}]:".\
					format(len(blackids),gene)
		ntaxa.append(len(taxids))
		gene_seqs = []
		for j, taxid_file in enumerate(taxid_files):
			tax_seqs = []
			handle = open(os.path.join(gene_dir, taxid_file), "rU")
			for record in pG.SeqIO.parse(handle, "fasta"):
				if taxids[j] == outgroup:
					record.id = "outgroup"
				else:
					record.id = "tx" + taxids[j] # rename seqs with taxids
				tax_seqs.append(record)
				nseqs += 1
			handle.close()
			gene_seqs.append(tax_seqs)
		seqs_obj.append(gene_seqs)
	if nseqs < 5:
		print "Too few sequences.\nDropping study."
		continue
	print 'Done. Read in [{0}] sequences for [{1}] genes and between [{2}] to [{3}] species'\
		.format(nseqs, len(genes), min(ntaxa), max(ntaxa))
	
	## Parse seqs
	prop_ind = overlap_ind = 0
	print 'Parsing sequences ... with gene overlap [{0}]'.\
		format(gene_overlaps[overlap_ind])
	parsed_seqs_obj,genes,nseqs,medians,ntaxa = parseSeqsObj(seqs_obj,\
		gene_overlaps[overlap_ind], genes, min_nspp)
	if len(parsed_seqs_obj) < 1:
		print "Too few species with sequences.\nDropping study."
		continue
	print 'After parsing, working with [{0}] genes for [{1}] to [{2}] species...'\
		.format(genes, min(ntaxa), max(ntaxa))
	print 'Done. Parsed [{0}] sequences.'.format(nseqs)
	
	## run alignments
	print "\nRunning alignments ..."
	aligns_std = []
	genes_used = []
	noaligns = False
	naligns_std = 0
	for j in range(len(parsed_seqs_obj)):
		temp_gene = genes[j]
		gene_obj = parsed_seqs_obj[j]
		temp_median = medians[j]
		print "Aligning gene [{0}] for [{1}] species ...".\
			format(temp_gene, ntaxa[j])
		reparse = False
		while True:
			aligns_gene, nruns = alignSeqsObj(gene_obj, align_len_max, nfails,\
				temp_median, props_nspp[prop_ind])
			print "-- ran [{0}] iterations".format(nruns)
			if len(aligns_gene) > 1:
				break
			if prop_ind < len(props_nspp) - 1:
				prop_ind += 1
				nspp_prop = int(len(gene_obj) * props_nspp[prop_ind])
				print "Too many failed alignments ... re-running with subsample [{0}] for mean [{1}] species".\
					format(props_nspp[prop_ind], nspp_prop)
				if nspp_prop < min_nspp:
					print "... too few species at that subsample ..."
					reparse = True
			else:
				reparse = True
			if reparse:
				reparse = False		
				if overlap_ind < len(gene_overlaps) - 1:
					overlap_ind += 1
					prop_ind = 0
					print "... re-parsing with gene overlap [{0}]".\
						format(gene_overlaps[overlap_ind])
					gene_obj,temp_gene,temp_nseqs,temp_median,temp_ntaxa =\
						parseSeqsObj([gene_obj], gene_overlaps[overlap_ind],\
						[temp_gene], min_nspp)
					if len(gene_obj) < 1:
						print "... no sequences left at that overlap!"
						noaligns = True
						break
					gene_obj = gene_obj[0]
					temp_gene = temp_gene[0]
					temp_median = temp_median[0]
					if len(gene_obj) < min_nspp:
						print "... too few species left!"
						noaligns = True
						break
					print 'Done. Parsed [{0}] sequences for [{1}] species.'.\
						format(temp_nseqs, temp_ntaxa[0])
				else:
					print "... minimum gene overlap hit!"
					noaligns = True
					break
		if noaligns:
			print "No alignments could be performed for gene [{0}]".format(genes[j])
			noaligns = False
		else:
			print "Done. Generated [{0}] alignments.".format(len(aligns_gene))
			naligns_std += len(aligns_gene)
			genes_used.append(temp_gene)
			aligns_std.append(aligns_gene)
		prop_ind = overlap_ind = 0
	if naligns_std < 1:
		print "No alignments could be generated for any genes.\nDropping study."
		continue
	print "Done. Generated [{0}] alignments in total for [{1}].".\
		format(naligns_std,studies[i])
	naligns_all += naligns_std
	counter += 1
	for j in range(len(aligns_std)):
		for k in range(len(aligns_std[j])):
			align = aligns_std[j][k]
			handle = open(os.path.join(output_dir, studies[i] + "_gene_" +\
				genes_used[j] + "_" + str(k) + '.faa'), "w")
			count = pG.SeqIO.write(align, handle, "fasta")
			handle.close()

## Remove mafft files
mafft_files = os.listdir(os.getcwd())
mafft_files = [f for f in mafft_files if re.search("\.fasta$", f)]
for f in mafft_files:
	os.remove(f)
print 'Stage finished. Generated [{0}] alignments across [{1}] studies.'.format(naligns_all, counter)
