#! /bin/usr/env python
# D.J. Bennett
# 24/03/2014
"""
pglt download tools
"""

# PACKAGES
import os
import re
import random
import entrez_tools as etools
import alignment_tools as atools
from special_tools import getThreads
from Bio.SeqFeature import SeqFeature


# CLASSES
class Downloader(object):
    """Download sequences given taxids and gene_names"""
    def __init__(self, gene_names, nseqs, thoroughness, maxpn, votesize,
                 maxtrys, minoverlap, maxlen, minlen, logger, wd=os.getcwd()):
        self.wd = wd
        self.logger = logger
        self.threads = getThreads(wd=wd)
        self.gene_names = gene_names
        self.nseqs = nseqs
        self.max_thoroughness = thoroughness
        self.maxpn = maxpn
        self.votesize = votesize
        self.maxtrys = maxtrys
        self.minoverlap = minoverlap
        self.maxlen = maxlen
        self.minlen = minlen
        self.thoroughness = 1
        self.deja_vues = []
        self.pattern = re.compile("[ACTGactg]")

    def _buildSearchTerm(self, taxids, thoroughness):
        """Generate NCBI GenBank query given taxids, gene_names and
thoroughness"""
        taxids_term = "txid{0}[PORGN]".format(taxids[0])
        for taxid in taxids[1:]:
            taxids_term = ("({0}) OR txid{1}[PORGN]".
                           format(taxids_term, taxid))
        # for fields see: http://www.ncbi.nlm.nih.gov/books/NBK49540/
        gene_names = ["\"{0}\"".format(e) for e in self.gene_names]
        # use gene field and ignore: predicted, genome
        if thoroughness == 1:
            gene_term = ("{0}[GENE]".format(gene_names[0]))
            for gene_name in gene_names[1:]:
                gene_term = ("({0}) OR {1}[GENE]".format(gene_term, gene_name))
            search_term = ("{0} AND ({1}) NOT predicted[TI] NOT \
genome[TI] NOT unverified[TI]".format(taxids_term, gene_term))
        # use title for gene name and ignore: predicted, genome
        elif thoroughness == 2:
            gene_term = ("{0}[TI]".format(gene_names[0]))
            for gene_name in gene_names[1:]:
                gene_term = ("({0}) OR {1}[TI]".format(gene_term, gene_name))
            search_term = ("{0} AND ({1}) NOT predicted[TI] NOT \
genome[TI] NOT unverified[TI]".format(taxids_term, gene_term))
        # all fields, ignore whole genome shotguns scaffolds and
        #  genome assemblies
        else:
            gene_term = ("{0}".format(gene_names[0]))
            for gene_name in gene_names[1:]:
                gene_term = ("({0}) OR {1}".format(gene_term, gene_name))
            search_term = ("{0} AND ({1}) NOT predicted[TI] NOT \
shotgun[TI] NOT scaffold[TI] NOT assembly[TI] NOT unverified[TI]".
                           format(taxids_term, gene_term))
        # self.logger.debug(search_term)
        return search_term

    def _search(self, taxids):
        """Search GenBank for matches, increasing thoroughness if
fewer matches than target nseqs"""
        # record seqids of all sequences that match query
        # move to next thoroughness if too few sequences
        # if _search is run again, will not return seqids of sequences
        #  that have already been found -- self.deja_vues records seen
        #  sequences
        # limit number of taxids to 5, search recursively if more than 5
        if len(taxids) > 5:
            res = self._search(taxids[:5])
            self.thoroughness = 1
            res.extend(self._search(taxids[5:]))
            return res
        seqids = []
        while len(seqids) < self.nseqs:
            if self.thoroughness > self.max_thoroughness:
                # keep searching until enough sequences have been
                #  found or until the max thoroughness has been hit
                break
            search_term = self._buildSearchTerm(taxids, self.thoroughness)
            seqcount = etools.eSearch(search_term, logger=self.logger)['Count']
            if int(seqcount) >= 1:
                # return ALL matching seqids if more than 0
                seqids.extend(etools.eSearch(search_term, logger=self.logger,
                                             retMax=seqcount)['IdList'])
                # filter those that have already been seen
                seqids = [e for e in seqids if e not in self.deja_vues]
            self.thoroughness += 1
        self.deja_vues.extend(seqids)
        self.deja_vues = list(set(self.deja_vues))
        return list(set(seqids))

    def _filter(self, sequences):
        """Filter sequences by BLASTing"""
        # choose random species for query
        randn = random.randint(0, len(sequences)-1)
        query = sequences
        subj = [sequences[randn]]
        # blast rand seq against all other seqs
        blast_bool, _ = atools.blast(query, subj, self.minoverlap, self.logger,
                                     wd=self.wd, threads=self.threads)
        # filtered are all sequences that are true
        filtered = [sequences[i] for i, e in enumerate(blast_bool) if e]
        # sequence pool are all sequences that are false
        seqpool = [sequences[i] for i, e in enumerate(blast_bool) if not e]
        # return filtered if there are more than votesize sequences in
        #  filtered
        if len(filtered) > self.votesize:
            return filtered, seqpool
        # else return empty list of filtered and the sequences
        else:
            return [], sequences

    def _findGeneInSeq(self, record):
        """Extract gene sequence from larger sequence (e.g. genomes)
by searching features."""
        if not record.features:
            # if there aren't any features, just return the record
            return record
        for feature in record.features:
            feature_names = []
            if 'gene' in feature.qualifiers.keys():
                feature_names.extend(feature.qualifiers['gene'])
            if 'gene_synonym' in feature.qualifiers.keys():
                feature_names.extend(feature.qualifiers['gene_synonym'])
            if 'product' in feature.qualifiers.keys():
                feature_names.extend(feature.qualifiers['product'])
            gene_names = [e.lower() for e in self.gene_names]
            feature_names = [e.lower() for e in feature_names]
            if set(gene_names) & set(feature_names):
                try:
                    extractor = SeqFeature(feature.location)
                    found_seq = extractor.extract(record)
                except ValueError:
                    # catch value errors raised for sequences
                    #  with "fuzzy" positions
                    # TODO: what are fuzzy positions and can I use
                    #  them?
                    return record
                else:
                    return found_seq
        return record

    def _parse(self, record):
        """Parse record returned from GenBank"""
        if isinstance(record, list):
            # find which sequence in the list has the gene
            for each in record:
                for gene in self.gene_names:
                    if gene in each.description.lower():
                        record = each
                        break
        if isinstance(record, list):
            return None
        # Always search through features for gene
        record = self._findGeneInSeq(record)
        if self.maxlen > len(record) > self.minlen:
            # find proportion of ambiguous bases
            nn = len(self.pattern.findall(str(record.seq)))
            pn = float(len(record)) - float(nn)
            pn = pn/float(len(record))
            if pn < self.maxpn:
                return record
        return None

    def _download(self, seqids):
        """Download records from GenBank given sequence ids"""
        records = []
        i = 1
        while i <= self.nseqs and len(seqids) > 0:
            if len(seqids) > 100:
                n = 100  # Download in chunks of 100
            else:
                n = len(seqids)
            seqs = []
            for _ in range(n):
                randi = random.randint(0, len(seqids)-1)
                seqs.append(seqids.pop(randi))
            for record in etools.eFetch(seqs, logger=self.logger):
                record = self._parse(record)
                if record:
                    records.append(record)
                    i += 1
        return records

    def run(self, taxids):
        """Dynamic sequence download"""
        # TODO: too complex, consider breaking up
        sequences = []
        while self.thoroughness < self.max_thoroughness:
            seqids = self._search(taxids)
            # filter if there are 10 times target nseqs
            if len(seqids) >= self.nseqs*10:
                self.logger.info("........ filtering")
                downloaded = []
                lower = 0
                while len(sequences) < self.nseqs and lower != len(seqids):
                    upper = min(len(seqids), lower + 1000)
                    downloaded.extend(self._download(seqids[lower:upper]))
                    if not downloaded:
                        break
                    filtered, downloaded = self._filter(downloaded)
                    if filtered:
                        sequences.extend(filtered)
                    lower = upper
                if len(sequences) > 0:
                    return sequences
            if len(seqids) > 0:
                sequences = self._download(seqids)
                if sequences:
                    return sequences
        return None


# FUNCTIONS
def findBestGenes(namesdict, genedict, thoroughness, allrankids, logger,
                  minnseq=1, target='all', minnspp=5):
    """Return suitable genes for phylogeny by searching for \
matches in GenBank"""
    # TODO: too complex, consider breaking up
    def nextBest(searchlist, types):
        # return best genes based on n and types
        ns = [e[1] for e in searchlist]
        # if types is None, ignore types and return gene with most n
        if types is None:
            gene_i = ns.index(max(ns))
            types = None
        # else find all types for each gene
        else:
            search_types = [genedict[e[0]]['type'] for e in searchlist]
            # if no types or search types and types have same number
            #  as unique types, return gene with max n
            if len(types) == 0 or len(list(set(search_types))) == \
                    len(list(set(types))):
                gene_i = ns.index(max(ns))
            # else, loop until the next type with the highest n is
            #  reached
            else:
                while True:
                    gene_i = ns.index(max(ns))
                    gene = searchlist[gene_i]
                    if genedict[gene]['type'] not in types:
                        break
        # return gene, shrunk searchlist and types
        gene = searchlist.pop(gene_i)
        return gene[0], searchlist, types
    alltipids = [namesdict[e]["txids"] for e in namesdict.keys()]
    outgroupids = namesdict['outgroup']['txids']
    searchlist = []
    # list of bools for number of genes w/o outgroup seqs
    outgroup_bool = []
    for gene in genedict.keys():
        logger.info('.... checking [{0}]'.format(gene))
        # first check if its suitable for this taxonomic group
        taxid = genedict[gene]["taxid"]
        gene_type = genedict[gene]["type"]
        if int(taxid) in allrankids:
            # if it is search genbank
            gene_bool = []
            for tipids in alltipids:
                downloader = Downloader(gene_names=genedict[gene]["names"],
                                        nseqs=minnseq + 1,
                                        thoroughness=thoroughness, maxpn=0,
                                        votesize=0, maxtrys=0, minoverlap=0,
                                        maxlen=0, minlen=0, logger=logger)
                res = downloader._search(tipids)
                # if gene is deep or both, then make sure it has
                #  outgroup
                if gene_type != 'shallow':
                    # if outgroupids do not have sequences, move to
                    #  next genes
                    if tipids == outgroupids:
                        if len(res) < minnseq:
                            outgroup_bool.append(False)
                            continue
                        else:
                            outgroup_bool.append(True)
                gene_bool.append(len(res) >= minnseq)
            nspp = sum(gene_bool)
            if nspp > minnspp:
                # if more than minnspp species, add it to searchlist
                searchlist.append((gene, nspp))
    # if all outgroup_bool are false, raise error
    if not any(outgroup_bool):
        # TODO: allow the program to run without an outgroup?
        raise atools.OutgroupError
    if target == 'all' or target > len(searchlist):
        return [e[0] for e in searchlist]
    else:
        if 'type' in genedict[genedict.keys()[0]].keys():
            # only if types have been provided in gene parameters
            types = []
        else:
            types = None
        genes = []
        while len(genes) < target:
            gene, searchlist, types = nextBest(searchlist, types)
            genes.append(gene)
    return genes


def getClusters(gene_sequences, minoverlap, logger, wd):
    """Identify clusters in sequences"""
    def findClusters(gene_sequences):
        # blast all against 1
        sequences = [e[1] for e in gene_sequences]
        randi = random.randint(0, len(sequences)-1)
        bools, _ = atools.blast(sequences, sequences[randi], minoverlap,
                                logger, wd, threads)
        # how many species had sequences in the cluster?
        cluster_sequences = [gene_sequences[i] for i, e in enumerate(bools)
                             if e]
        nspp = len(set([e[0] for e in cluster_sequences]))
        pspp = float(nspp)/tot_nspp
        # if more than 50% and 5 species ...
        if pspp > 0.5 and nspp > 5:
            # return cluster, remove those sequences from gene_sequences
            gene_sequences = [gene_sequences[i] for i, e in enumerate(bools)
                              if not e]
            return cluster_sequences, gene_sequences
        return None, gene_sequences
    threads = getThreads(wd=wd)
    res = []
    tot_nspp = len(set([e[0] for e in gene_sequences]))
    # try max 5 times to get a cluster from randomly selecting a seq
    for i in range(5):
        cluster_sequences, gene_sequences = findClusters(gene_sequences)
        if cluster_sequences:
            res.append(cluster_sequences)
        # if gene sequences has not enough seqs left, break
        pspp = float(len(set([e[0] for e in gene_sequences])))/tot_nspp
        if pspp < 0.5:
            break
    return res
