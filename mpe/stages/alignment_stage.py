#! /usr/bin/env python
## D.J. Bennett
## 24/03/2014
"""
mpe Stage 3: Aligning sequences
"""

## Packages
import os,re,pickle,logging,shutil
import numpy
from Bio import SeqIO
import mpe.tools.alignment_tools as atools

## Functions
def writeAlignment(alignment, i, namesdict, gene_dir):
	# log alignment details
	logging.info(".... alignment length [{0}] for [{1}] species".\
		format(alignment.get_alignment_length(), len(alignment)))
	# record records in alignment
	for record in alignment:
		namesdict[record.id]['alignments'] += 1
	# write out
	align_len = alignment.get_alignment_length()
	output_file = "{0}_nspp{1}_len{2}.faa".format(i,len(alignment),\
		align_len)
	output_path = os.path.join(gene_dir, output_file)
	with open(output_path, "w") as file:
		count = SeqIO.write(alignment, file, "fasta")
		del count
	return namesdict, None

def run(wd = os.getcwd()):
	## Print stage
	logging.info("Stage 3: Sequence alignment")

	## Dirs
	download_dir = os.path.join(wd, '2_download')
	alignment_dir = os.path.join(wd, '3_alignment')
	if not os.path.isdir(alignment_dir):
		os.mkdir(alignment_dir)

	## Input
	with open(os.path.join(wd, ".genedict.p"), "rb") as file:
		genedict = pickle.load(file)
	with open(os.path.join(wd, ".paradict.p"), "rb") as file:
		paradict = pickle.load(file)
	with open(os.path.join(wd, ".namesdict.p"), "rb") as file:
		namesdict = pickle.load(file)

	## Parameters
	naligns = int(paradict["naligns"])
	all_counter = 0

	## Read in sequences
	# add alignments to namesdict
	for key in namesdict.keys():
		namesdict[key]['alignments'] = 0
	genes = sorted(os.listdir(download_dir))
	genes = [e for e in genes if not re.search("^\.|^log\.txt$", e)]
	genekeys = {}
	for gene in genes:
		genekeys[gene] = re.sub('_cluster[0-9]+', '', gene)
	logging.info('Reading in sequences ....')
	genestore = []
	for gene in genes:
		gene_dir = os.path.join(download_dir, gene)
		seq_files = os.listdir(gene_dir)
		seqstore = atools.SeqStore(gene_dir, seq_files, minfails =\
			int(genedict[genekeys[gene]]["minfails"]), mingaps = \
			float(genedict[genekeys[gene]]["mingaps"]), minoverlap =\
			int(genedict[genekeys[gene]]["minoverlap"]))
		genestore.append((gene, seqstore))

	## Run alignments
	logging.info("Running alignments ....")
	# loop through genes
	for gene,seqstore in genestore:
		# set up dir
		gene_dir = os.path.join(alignment_dir, gene)
		if not os.path.isdir(gene_dir):
			os.mkdir(gene_dir)
		logging.info("Aligning gene [{0}] for [{1}] species ....".\
			format(gene, len(seqstore)))
		# set up aligner obj
		aligner = atools.Aligner(seqstore, mingaps = \
			float(genedict[genekeys[gene]]["mingaps"]), minoverlap = \
			int(genedict[genekeys[gene]]["minoverlap"]), minseedsize = \
			int(genedict[genekeys[gene]]["minseedsize"]), maxseedsize = \
			int(genedict[genekeys[gene]]["maxseedsize"]), maxtrys = \
			int(genedict[genekeys[gene]]["maxtrys"]), maxseedtrys = \
			int(genedict[genekeys[gene]]["maxseedtrys"]), gene_type = \
			genedict[genekeys[gene]]['type'], outgroup = \
			'outgroup' in seqstore.keys())
		# run for naligns
		each_counter = 0
		alignment = None
		try:
			for i in range(1, naligns + 1):
				logging.info(".... iteration [{0}]".format(i))
				while not alignment:
					alignment = aligner.run()
				namesdict, alignment = writeAlignment(alignment, i,\
					namesdict, gene_dir)
				each_counter += 1
		except atools.TrysError:
			logging.info(".... max trys hit")
		except atools.OutgroupError:
			logging.info(".... outgroup dropped")
		except atools.TooFewSpeciesError:
			logging.info(".... too few species left in sequence pool")
		if each_counter < naligns:
			logging.info(".... too few alignments generated")
			shutil.rmtree(gene_dir)
			continue
		all_counter += each_counter

	## Wrap-up
	# the number of alignments per name in namesdict
	naligns_name = [namesdict[e]['alignments'] for e in \
	namesdict.keys() if namesdict[e]['genes'] > 0]
	# the proportion of alignments each name has of all alignments
	paligns_name = [len(naligns_name) * float(e)/all_counter for e\
	in naligns_name]
	with open(os.path.join(wd, ".namesdict.p"), "wb") as file:
		pickle.dump(namesdict, file)
	logging.info('Stage finished. Generated [{n1}] alignments for \
mean [{n2:.{d}f}](sd[{n3:.{d}f}]) species.'.format(n1 = all_counter,\
	n2 = numpy.mean(paligns_name), n3 = numpy.std(paligns_name),\
	d = 2))