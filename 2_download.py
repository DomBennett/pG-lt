#!/usr/bin/python
## MRes Project 2013
## Stage 2: Downloading sequences
## In: 0_names, 1_taxids | Out: 2_download
## 07/08/2013

## Print stage
print "\n\nThis is stage 2: download\n"

## Parameters
dwnload_nseqs = 100

## Packages
import sys, os, re, csv
sys.path.append(os.path.join(os.getcwd(), 'functions'))
import phyloGenerator_adapted as pG

## Entrez email
pG.Entrez.email = 'dominic.john.bennett@gmail.com'
print 'Using email [{0}] for NCBI Entrez'.format(pG.Entrez.email)

## Dirs
input_dirs = [os.path.join(os.getcwd(), '0_names'),\
	os.path.join(os.getcwd(), '1_taxids')]
output_dir = os.path.join(os.getcwd(), '2_download')
if not os.path.isdir(output_dir):
    os.mkdir(output_dir)
files = os.listdir(input_dirs[1])
taxids_files = sorted([e for e in files if re.search("^.*taxids\.txt$", e)])
#shared_files = sorted([e for e in files if re.search("^.*shared\.txt$", e)])

## Taxadata
print "Reading in taxadata.csv ..."
taxadict = {}
with open(os.path.join(input_dirs[0], 'taxadata.csv'), 'rb') as csvfile:
	taxreader = csv.DictReader(csvfile)
	for row in taxreader:
		temp_genes = row['genes'].split("|")
		temp_genes = [e1.split(",") for e1 in temp_genes]
		taxadict[row['study']] = temp_genes
print "Done. Read in taxadata for [{0}] studies.".format(len(taxadict))

## Loop through studies
print '\nLooping through studies ...'
counter = nseqs = nbases = nspp = 0
for i in range(len(taxids_files)):

	## What study?
	current_study = re.sub("_taxids.txt$", "", taxids_files[i])
	print '\n\nWorking on: [{0}]'.format(current_study)
	print 'Reading data and determining genes to search for ...'
	
	## read in taxids
	taxids = []
	with file(os.path.join(input_dirs[1], taxids_files[i]), 'rb') as infile:
		for each in infile:
			taxids.append(each.strip())
	if len(taxids) < 2:
		print 'Too few taxids. Dropping study.'
		continue
	print '... Found [{0}] taxids in input file ...'.format(len(taxids))
	taxids = [int(taxid) for taxid in taxids] # cos' they're imported as strs!
	
	## work out genes
	print 'Extracting genes ...'
	genes = taxadict[current_study]
	print 'Using genes: [{0}] for sequence download.'.format(genes)
	
	## Sequence download
	nseqs_study = 0
	nspp_study = []
	study_dir  = os.path.join(output_dir, str(current_study))
	if not os.path.isdir(study_dir):
			os.mkdir(study_dir)
	for gene in genes:
		nseqs_gene = no_seqs_gene = nspp_gene = 0
		print 'Downloading and outputting for [{0}]...'.format(gene[0])
		gene_dir  = os.path.join(study_dir, str(gene[0]))
		if not os.path.isdir(gene_dir):
			os.mkdir(gene_dir)
		for taxid in taxids:
			# download
			try:
				seqs_obj = pG.findGenes(speciesList = [taxid], geneNames =\
					[gene], download = True, thorough = False, noSeqs =\
					dwnload_nseqs, verbose = False, retMax = dwnload_nseqs)
			except:
				print "pG.findGenes called an unexpected error. Dropping taxid [{0}].".\
					format(taxid)
				continue
			no_seqs = isinstance(seqs_obj[0][0][0], tuple)
			if no_seqs:
				try:
					seqs_obj = pG.findGenes(speciesList = [taxid], geneNames =\
					[gene], download = True, thorough = True, noSeqs =\
					dwnload_nseqs, verbose = False, retMax = dwnload_nseqs)
				except:
					print "pG.findGenes called an unexpected error. Dropping taxid [{0}].".\
					format(taxid)
					continue
				no_seqs = isinstance(seqs_obj[0][0][0], tuple)
				if no_seqs:
					no_seqs_gene += 1
					print "No sequences found for taxid [{0}].".\
					format(taxid)
					continue
			# unpack + write
			seqs = seqs_obj[0][0][0]
			if isinstance(seqs, list):
				gene_seqs = []
				for s in seqs:
					if isinstance(s, list): # sometimes a list of a list is returned
						gene_seqs.append([e.format('fasta') for e in s])
					else:
						gene_seqs.append(s.format('fasta'))
			elif isinstance(seqs, pG.SeqRecord):
				gene_seqs = seqs.format('fasta')
			else:
				print "Error seqs is an unexpected object"
				continue
			with file(os.path.join(gene_dir, str(taxid) + '.fasta'), 'wb') \
				as outfile:
				if isinstance(gene_seqs, list):
					for gene_seq in gene_seqs:
						outfile.write("%s\n" % gene_seq)
						nseqs_gene += 1
						nbases += len(gene_seq)
				else:
					outfile.write(gene_seqs)
					nseqs_gene += 1
					nbases += len(gene_seqs)
			nspp_gene += 1
		if no_seqs_gene == len(taxids):
			print "No sequences were downloaded for gene [{0}]".format(gene[0])
			os.rmdir(gene_dir)
		else:
			nseqs_study += nseqs_gene
			nspp_study.append(nspp_gene)
			print "Downloaded [{0}] sequences for gene [{1}] representing [{2}] species".\
				format(nseqs_gene,gene[0],nspp_gene)
	print 'Done. Downloaded [{0}] sequences for study [{1}] representing [{2}] species.'.\
		format(nseqs_study,current_study,max(nspp_study))
	counter += 1
	nseqs += nseqs_study
	nspp += max(nspp_study)
print '\n\nStage finished. Downloaded [{0}] bases for [{1}] sequences for [{2}] studies for [{3}] species.'.\
	format(nbases, nseqs, counter, nspp)
