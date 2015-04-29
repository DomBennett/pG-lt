#! /usr/bin/env python
# D.J. Bennett
# 24/03/2014
"""
pglt Stage 3: Aligning sequences
"""

# PACAKAGES
import os
import re
import pickle
import logging
import shutil
import numpy
from Bio import SeqIO
import pglt.tools.alignment_tools as atools
from pglt.tools.system_tools import MissingDepError


# FUNCTIONS
def readSequences(download_dir, namesdict, genedict, logger, wd):
    """Read sequences into a genestore"""
    # add alignments key to namesdict
    for key in namesdict.keys():
        namesdict[key]['alignments'] = 0
    genes = sorted(os.listdir(download_dir))
    genes = [e for e in genes if not re.search("^\.|^log\.txt$", e)]
    genekeys = {}
    for gene in genes:
        genekeys[gene] = re.sub('_cluster[0-9]+', '', gene)
    genestore = []
    for gene in genes:
        gene_dir = os.path.join(download_dir, gene)
        seq_files = os.listdir(gene_dir)
        minfails = int(genedict[genekeys[gene]]["minfails"])
        mingaps = float(genedict[genekeys[gene]]["mingaps"])
        minoverlap = int(genedict[genekeys[gene]]["minoverlap"])
        seqstore = atools.SeqStore(gene_dir, seq_files, minfails=minfails,
                                   mingaps=mingaps, minoverlap=minoverlap,
                                   logger=logger, wd=wd)
        genestore.append((gene, seqstore))
    return namesdict, genestore, genekeys


def setUpAligner(gene, genedict, genekeys, seqstore, logger, wd):
    """Set-up Aligner class for gene"""
    # get parameters
    mingaps = float(genedict[genekeys[gene]]["mingaps"])
    minoverlap = int(genedict[genekeys[gene]]["minoverlap"])
    minseedsize = int(genedict[genekeys[gene]]["minseedsize"])
    maxseedsize = int(genedict[genekeys[gene]]["maxseedsize"])
    maxtrys = int(genedict[genekeys[gene]]["maxtrys"])
    maxseedtrys = int(genedict[genekeys[gene]]["maxseedtrys"])
    gene_type = genedict[genekeys[gene]]['type']
    outgroup = 'outgroup' in seqstore.keys()
    # set up aligner obj
    aligner = atools.Aligner(seqstore, mingaps=mingaps, minoverlap=minoverlap,
                             minseedsize=minseedsize, maxseedsize=maxseedsize,
                             maxtrys=maxtrys, maxseedtrys=maxseedtrys,
                             gene_type=gene_type, outgroup=outgroup,
                             logger=logger, wd=wd)
    return aligner


def runAligner(aligner, naligns, namesdict, gene_dir, logger):
    """Run aligner, save alignments as their made"""
    # run for naligns
    each_counter = 0
    alignment = None
    try:
        for i in range(1, naligns + 1):
            logger.info(".... iteration [{0}]".format(i))
            while not alignment:
                alignment = aligner.run()
            # log alignment details
            logger.info(".... alignment length [{0}] for [{1}] species".
                        format(alignment.get_alignment_length(),
                               len(alignment)))
            namesdict, alignment = writeAlignment(alignment, i, namesdict,
                                                  gene_dir)
            each_counter += 1
    except atools.TrysError:
        logger.info(".... max trys hit")
    except atools.OutgroupError:
        logger.info(".... outgroup dropped")
    except atools.TooFewSpeciesError:
        logger.info(".... too few species left in sequence pool")
    if each_counter < naligns:
        logger.info(".... too few alignments generated")
        shutil.rmtree(gene_dir)
        return 0
    return each_counter


def writeAlignment(alignment, i, namesdict, gene_dir):
    """Save alignment"""
    # record records in alignment
    for record in alignment:
        namesdict[record.id]['alignments'] += 1
    # write out
    align_len = alignment.get_alignment_length()
    output_file = "{0}_nspp{1}_len{2}.faa".format(i, len(alignment), align_len)
    output_path = os.path.join(gene_dir, output_file)
    with open(output_path, "w") as file:
        count = SeqIO.write(alignment, file, "fasta")
        del count
    return namesdict, None


# RUN
def run(wd=os.getcwd(), logger=logging.getLogger('')):
    # PRINT STAGE
    logger.info("Stage 3: Sequence alignment")

    # DIRS
    download_dir = os.path.join(wd, '2_download')
    alignment_dir = os.path.join(wd, '3_alignment')
    temp_dir = os.path.join(wd, 'tempfiles')
    if not os.path.isdir(alignment_dir):
        os.mkdir(alignment_dir)

    # CHECK DEPS
    if not atools.blastn:
        raise MissingDepError('blastn')
    if not atools.mafft:
        raise MissingDepError('mafft')

    # INPUT
    with open(os.path.join(temp_dir, "genedict.p"), "rb") as file:
        genedict = pickle.load(file)
    with open(os.path.join(temp_dir, "paradict.p"), "rb") as file:
        paradict = pickle.load(file)
    with open(os.path.join(temp_dir, "namesdict.p"), "rb") as file:
        namesdict = pickle.load(file)

    # PARAMETERS
    naligns = int(paradict["naligns"])
    all_counter = 0

    # READ IN SEQUENCES
    logger.info('Reading in sequences ....')
    namesdict, genestore, genekeys = readSequences(download_dir, namesdict,
                                                   genedict, logger, temp_dir)

    # RUN ALIGNMENTS
    logger.info("Running alignments ....")
    # loop through genes
    for gene, seqstore in genestore:
        logger.info("Aligning gene [{0}] for [{1}] species ....".
                    format(gene, len(seqstore)))
        # set up dir
        gene_dir = os.path.join(alignment_dir, gene)
        if not os.path.isdir(gene_dir):
            os.mkdir(gene_dir)
        aligner = setUpAligner(gene, genedict, genekeys, seqstore, logger,
                               temp_dir)
        all_counter += runAligner(aligner, naligns, namesdict, gene_dir,
                                  logger)

    # CALC STATS
    # the number of alignments per name in namesdict
    naligns_name = [namesdict[e]['alignments'] for e in namesdict.keys() if
                    namesdict[e]['genes'] > 0]
    # the proportion of alignments each name has of all alignments
    paligns_name = [len(naligns_name) * float(e)/all_counter for
                    e in naligns_name]

    # OUTPUT
    with open(os.path.join(temp_dir, "namesdict.p"), "wb") as file:
        pickle.dump(namesdict, file)

    # FINISH MESSAGE
    logger.info('Stage finished. Generated [{n1}] alignments for \
mean [{n2:.{d}f}](sd[{n3:.{d}f}]) species.'.
                format(n1=all_counter, n2=numpy.mean(paligns_name),
                       n3=numpy.std(paligns_name), d=2))
