#! /usr/bin/env python
# D.J. Bennett
# 24/03/2014
"""
pglt Stage 2: Sequence Download
"""

# PACKAGES
import os
import pickle
import logging
import pglt.tools.download_tools as dtools


def run(wd=os.getcwd(), logger=logging.getLogger('')):
    # TODO: too complex, consider breaking up
    # PRINT STAGE
    logger.info("Stage 2: Sequence download")

    # DIRS
    download_dir = os.path.join(wd, '2_download')
    if not os.path.isdir(download_dir):
        os.mkdir(download_dir)

    # INPUT
    with open(os.path.join(wd, ".genedict.p"), "rb") as file:
        genedict = pickle.load(file)
    with open(os.path.join(wd, ".paradict.p"), "rb") as file:
        paradict = pickle.load(file)
    with open(os.path.join(wd, ".namesdict.p"), "rb") as file:
        namesdict = pickle.load(file)
    with open(os.path.join(wd, ".allrankids.p"), "rb") as file:
        allrankids = pickle.load(file)

    # PARAMETERS
    dtools.etools.Entrez.email = paradict["email"]
    nseqs = int(paradict['nseqs'])
    thoroughness = int(paradict['thoroughness'])
    dtools.logger = logger
    dtools.atools.logger = logger
    seedsize = 10
    maxtrys = 100
    minnseq = 1
    minnspp = 5
    target = 6
    maxpn = 0.1
    seqcounter = basecounter = 0

    # PROCESS
    logger.info('Determining best genes ....')
    genes = dtools.findBestGenes(namesdict, genedict, thoroughness, allrankids,
                                 minnseq, target, minnspp)
    statement = 'Using genes:'
    for gene in genes:
        statement += " [" + gene + "]"
    logger.info(statement)
    # Add genes to namesdict
    for key in namesdict.keys():
        namesdict[key]['genes'] = 0
    for gene in genes:
        gene_sequences = []
        seqcounter_gene = noseqcounter_gene = spcounter_gene = 0
        gene_names = genedict[gene]["names"]
        minlen = int(genedict[gene]["minlen"])
        maxlen = int(genedict[gene]["maxlen"])
        minoverlap = int(genedict[gene]['minoverlap'])
        logger.info('Downloading and outputting for [{0}] ....'.format(gene))
        for name in namesdict.keys():
            logger.info("..... [{0}]".format(name))
            taxids = namesdict[name]["txids"]
            downloader = dtools.Downloader(gene_names, nseqs, thoroughness,
                                           maxpn, seedsize, maxtrys,
                                           minoverlap, maxlen, minlen)
            sequences = downloader.run(taxids)
            if not sequences:
                noseqcounter_gene += 1
                logger.info("........ no sequences found")
                continue
            logger.info("........ downloaded [{0}] sequences".
                        format(len(sequences)))
            gene_sequences.extend(zip([name] * len(sequences), sequences))
        if noseqcounter_gene == len(namesdict.keys()):
            logger.info("No sequences were downloaded for gene [{0}]".
                        format(gene))
            continue
        else:
            seqcounter += seqcounter_gene
        logger.info('Checking for distinct clusters ....')
        gene_sequences = dtools.getClusters(gene_sequences, minoverlap)
        if not gene_sequences:
            logger.info('.... could not find any clustering sequences')
            continue
        logger.info('.... found [{0}] clusters'.format(len(gene_sequences)))
        gene_dir = os.path.join(download_dir, str(gene))
        for i in range(len(gene_sequences)):
            outdir = '{0}_cluster{1}'.format(gene_dir, i)
            if not os.path.isdir(outdir):
                os.mkdir(outdir)
            for name in namesdict.keys():
                seqs = [e[1] for e in gene_sequences[i] if e[0] == name]
                if not seqs:
                    continue
                seqs = [e.format('fasta') for e in seqs]
                with open(os.path.join(outdir, "{0}.fasta".format(name)), 'wb')\
                        as outfile:
                    for seq in seqs:
                        outfile.write("{0}\n".format(seq))
                        seqcounter_gene += 1
                        basecounter += len(seq)
                namesdict[name]['genes'] += 1
                spcounter_gene += 1
        logger.info("Downloaded [{0}] sequences for gene [{1}] representing \
[{2}] species".format(seqcounter_gene, gene, spcounter_gene))
    with open(os.path.join(wd, ".namesdict.p"), "wb") as file:
        pickle.dump(namesdict, file)
    logger.info('Stage finished. Downloaded [{0}] bases for [{1}] \
sequences for [{2}] species.'.format(basecounter, seqcounter,
                                     sum([namesdict[e]['genes'] > 0 for e in
                                          namesdict.keys()])))
