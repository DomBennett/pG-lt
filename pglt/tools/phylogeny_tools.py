#! /bin/usr/env python
# D.J. Bennett
# 24/03/2014
"""
pglt phylogeny tools
"""

# Packages
import os
import re
import random
import logging
from collections import Counter
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio import Phylo
from Bio import AlignIO
import dendropy as dp
from math import sqrt
from system_tools import TerminationPipe
from system_tools import RAxMLError
from special_tools import getThreads
from pglt import _RAXML as raxml


# CLASSES
class StopCodonRetriever(object):
    """Stop codon retrival class"""
    # ref: http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
    # Multiple lists referring to:
    #  [0] Vertebrates (excl. Ascidia)
    #  [1] Ciliates, green algae and diplomonads
    #  [2] All other Eukaryotes (incl. Ascidia)
    mt_txids = [[7742], [5878, 3135, 5738], [2759, 30275]]
    mt_fpattern = ['(taa|tag|aga|agg)', 'tag', '(taa|tag)']
    mt_rpattern = ['(tta|cta|tct|cct)', 'cta', '(tta|cta)']
    # TODO: add Nuc stop codons (Problem no use though -- low priority)

    def __init__(self):
        pass

    def pattern(self, ids, logger, genome_type='mt'):
        """Return stop pattern for lowest matching id"""
        # assumes lower taxonomic levels are at higher indexes
        if genome_type == 'mt':
            # find all matching ids in mt_txids
            matches = [i for i, e in enumerate(ids) if
                       e in sum(self.mt_txids, [])]
            if not matches:
                logger.debug('No taxonomic IDs matched!')
                return None
            match = max(matches)
            i = [i for i, e in enumerate(self.mt_txids) if ids[match] in e][0]
            # return reobjs
            return re.compile(self.mt_fpattern[i], flags=re.IGNORECASE), \
                re.compile(self.mt_rpattern[i], flags=re.IGNORECASE)
        else:
            return None


class AlignmentStore(dict):
    """Alignment holding class"""
    retriever = StopCodonRetriever()

    def __init__(self, clusters, genedict, allrankids, indir, logger):
        self.logger = logger
        # Read in alignments for each cluster
        # first find corresponding gene name for each cluster
        genes = []
        for cluster in clusters:
            genes.append(re.sub('_cluster[0-9]+', '', cluster))
        for cluster, gene in zip(clusters, genes):
            # add a key to the AlignmentStore dict
            self[cluster] = {'alignments': [], 'files': [], 'counters': []}
            # retrieve its stop codon if it's mt
            self[cluster]['stop'] =\
                self.retriever.pattern(ids=allrankids, logger=self.logger,
                                       genome_type=genedict[gene]
                                       ['partition'].lower())
            # find corresponding input dir and read in alignments
            cluster_dir = os.path.join(indir, cluster)
            alignment_files = os.listdir(cluster_dir)
            alignment_files = [e for e in alignment_files if
                               not re.search("^\.", e)]
            for alignment_file in alignment_files:
                with open(os.path.join(cluster_dir, alignment_file), "r") \
                        as file:
                    alignment = AlignIO.read(file, "fasta")
                self[cluster]['alignments'].append(alignment)
                self[cluster]['files'].append(alignment_file)
                self[cluster]['counters'].append(0)

    def count(self):
        """Add to alignments' counters"""
        for each in self.counters:
            each += 1

    def report(self):
        """Return info on alignments used"""
        # TODO: find a way to record alignments used
        pass

    def pull(self):
        """Randomly select an alignment for each gene. Return
list of alignments and stop pattern if partition"""
        self.counters = []
        alignments = []
        stops = []
        self.logger.info("........ Using alignments:")
        for gene in self.keys():
            genedata = self[gene]
            i = random.randint(0, len(genedata['alignments']) - 1)
            alignments.append(genedata['alignments'][i])
            stops.append(genedata['stop'])
            self.counters.append(genedata['counters'][i])
            afile = genedata['files'][i]
            if genedata['stop']:
                self.logger.info("............ {0}(codon partitioned):\
[{1}]".format(gene, afile))
            else:
                self.logger.info("............ {0}:[{1}]".format(gene, afile))
        return alignments, stops


class Generator(object):
    """Phylogeny generating class"""
    def __init__(self, alignment_store, rttstat, outdir, maxtrys, logger,
                 wd=os.getcwd()):
        self.logger = logger
        self.wd = wd
        self.threads = getThreads(wd=wd)
        self.threads = self.threads + 1  # RAxML runs a small master process
        self.trys = 0
        self.phylogenies = []
        self.maxtrys = maxtrys
        self.alignment_store = alignment_store
        self.genes = alignment_store.keys()
        self.rttstat = rttstat
        self.outdir = outdir
        self.taxontree = os.path.join(outdir, "taxontree.tre")
        self.constraint = os.path.isfile(self.taxontree)

    def _test(self, phylogeny):
        """Return phylogeny if RTT stat is below max RTT stat"""
        if phylogeny:
            # remove outgroup, reduces RTT variance
            try:
                phylogeny.prune('outgroup')
            except:
                pass
            # calc root to tip distance (rtt.dist) for each tip
            rtt_dists = []
            for terminal in phylogeny.get_terminals():
                rtt_dists.append(phylogeny.distance(terminal))
            # calculate CoV
            mean = sum(rtt_dists)/len(rtt_dists)
            sd = sqrt(sum([(e-mean)**2 for e in rtt_dists]) /
                      len(rtt_dists))
            rttstat = sd/mean  # CoV
            # Phylo.draw_ascii(phylogeny)
            self.logger.debug('..... [{0}] RTT stat'.format(rttstat))
            if rttstat < self.rttstat:
                return phylogeny
            else:
                self.logger.info('........ poor phylogeny, retrying')
        else:
            self.logger.debug('.... no phylogeny, retrying')

    def _concatenate(self, alignments):
        """Return single alignment from list of alignments for
multiple genes."""
        if len(alignments) == 1:
            return alignments[0]
        # sort IDs
        alignment_ids = []
        for gene in alignments:
            gene_ids = []
            for rec in gene:
                gene_ids.append(rec.id)
            alignment_ids.append(gene_ids)
        all_ids = []
        [all_ids.extend(e) for e in alignment_ids]
        all_ids = list(set(all_ids))
        # concatenate
        alignment = MultipleSeqAlignment([])
        for txid in all_ids:
            sequence = ""
            for i, gene in enumerate(alignments):
                if txid in alignment_ids[i]:
                    sequence += gene[alignment_ids[i].index(txid)].seq
                else:
                    sequence += "-" * gene.get_alignment_length()
            sequence = SeqRecord(sequence, id=txid,
                                 description="multigene sequence")
            alignment.append(sequence)
        return alignment

    def _constraint(self, alignment):
        """Generate constraint tree using taxontree, return arg"""
        if not self.constraint:
            return False
        # drop tips from taxontree if not in alignment
        tip_names = []
        for record in alignment:
            tip_names.append(record.id)
        with open(self.taxontree, "r") as file:
            constraint = Phylo.read(file, "newick")
        constraint_tips = []
        for terminal in constraint.get_terminals():
            constraint_tips.append(terminal.name)
        tips_to_drop = [e for e in constraint_tips if e not in tip_names]
        for tip in tips_to_drop:
            constraint.prune(tip)
        # write out tree
        with open(os.path.join(self.wd, "constraint.tre"), "w") as file:
            Phylo.write(constraint, file, "newick")
        # return arg
        if constraint.is_bifurcating():
            return " -r constraint.tre"
        else:
            return " -g constraint.tre"

    def _outgroup(self, alignment):
        """Return arg for outgroup"""
        spp = [e.id for e in alignment]
        if 'outgroup' in spp:
            return 'outgroup'
        # otherwise find the species(s) with the fewest shared
        #  taxonomic groups
        if self.constraint:
            with open(os.path.join(self.wd, "constraint.tre"), "r") as file:
                constraint = Phylo.read(file, "newick")
            distances = [constraint.distance(e) for e in spp]
            index = [i for i, e in enumerate(distances) if e ==
                     min(distances)]
            # choose one at random
            ri = random.sample(index, 1)[0]
            return spp[ri]
        else:
            return False

    def _findORF(self, alignment, stop):
        """Return ORF of alignment based on absence of stop codons"""
        # TODO: too complex, consider breaking up
        def reframe(alignment, frame):
            # return alignment from point where frame starts
            alignment = alignment[:, frame:]
            offset = alignment.get_alignment_length() % 3
            if offset > 0:
                alignment = alignment[:, :-offset]
            return alignment
        if not stop:
            return alignment, False
        # Unpack stop patterns
        fstop, rstop = stop
        frame_stops = [0, 0, 0, 0, 0, 0]
        for record in alignment:
            # convert to string
            seq = str(record.seq)
            # search for stop codons ignoring last 50bps; expect
            #  a stop codon at the end of a sequence
            frame_stops[0] += sum([bool(fstop.match(seq[:-50]
                                  [e:e + 3])) for e in range(0, len(seq), 3)])
            frame_stops[1] += sum([bool(fstop.match(seq[:-50]
                                  [e:e + 3])) for e in range(1, len(seq), 3)])
            frame_stops[2] += sum([bool(fstop.match(seq[:-50]
                                  [e:e + 3])) for e in range(2, len(seq), 3)])
            frame_stops[3] += sum([bool(rstop.match(seq[50:]
                                  [e:e + 3])) for e in range(0, len(seq), 3)])
            frame_stops[4] += sum([bool(rstop.match(seq[50:]
                                  [e:e + 3])) for e in range(1, len(seq), 3)])
            frame_stops[5] += sum([bool(rstop.match(seq[50:]
                                  [e:e + 3])) for e in range(2, len(seq), 3)])
        # if more than one frame wo stop codon
        #  return wo codon partitions
        if sum([e == 0 for e in frame_stops]) > 1:
            return alignment, False
        # if no frames wo stop codons, return wo codon partitions
        if sum([e == 0 for e in frame_stops]) == 0:
            return alignment, False
        # else return frame and reframed alignment
        if frame_stops[0] == 0 or frame_stops[3] == 0:
            return reframe(alignment, 0), True
        if frame_stops[1] == 0 or frame_stops[4] == 0:
            return reframe(alignment, 1), True
        if frame_stops[2] == 0 or frame_stops[5] == 0:
            return reframe(alignment, 2), True

    def _partition(self, alignments, stops):
        """Return partition argument, write out partition postitions
to partitions.txt"""
        if len(alignments) == 1:
            if not stops[0]:
                return alignments, None
        begin = 1
        ngene = 1
        text = ''
        reframed = []
        for alignment, stop in zip(alignments, stops):
            # if stop pattern gets ORF, partition by codon ...
            alignment, partitioned = self._findORF(alignment, stop)
            if partitioned:
                end = alignment.get_alignment_length() + begin - 1
                text += 'DNA, gene{0}codon1 = {1}-{2}\\3\n'.\
                    format(ngene, begin, end)
                text += 'DNA, gene{0}codon2 = {1}-{2}\\3\n'.\
                    format(ngene, begin + 1, end)
                text += 'DNA, gene{0}codon3 = {1}-{2}\\3\n'.\
                    format(ngene, begin + 2, end)
            else:
                # ... else just for the whole gene
                end = alignment.get_alignment_length() + begin - 1
                text += 'DNA, gene{0} = {1}-{2}\n'.\
                    format(ngene, begin, end)
            begin = end + 1
            ngene += 1
            reframed.append(alignment)
        logging.debug([e.get_alignment_length() for e in alignments])
        logging.debug(text)
        with open(os.path.join(self.wd, 'partitions.txt'), 'w') as file:
            file.write(text)
        return reframed, ' -q partitions.txt'

    def _setUp(self, alignments, stops):
        """Set up for RAxML"""
        # partition
        alignments, parg = self._partition(alignments, stops)
        # create supermatrix alignment
        alignment = self._concatenate(alignments)
        # create constraint
        carg = self._constraint(alignment)
        # get outgroup arg
        outgroup = self._outgroup(alignment)
        return alignment, carg, outgroup, parg

    def run(self):
        """Generate phylogeny from alignments"""
        if self.trys > self.maxtrys:
            raise RAxMLError()
        # choose random alignment for each gene
        alignments, stops = self.alignment_store.pull()
        # set up
        alignment, carg, outgroup, parg = self._setUp(alignments, stops)
        # run RAxML
        phylogeny = RAxML(alignment, wd=self.wd, logger=self.logger,
                          threads=self.threads, constraint=carg,
                          outgroup=outgroup, partitions=parg)
        phylogeny = self._test(phylogeny)
        if phylogeny:
            self.phylogenies.append(phylogeny)
            self.trys = 0
            return True
        else:
            self.trys += 1
            return False


def RAxML(alignment, wd, logger, threads, outgroup=None, partitions=None,
          constraint=None, timeout=999999999):
    """Adapted pG function: Generate phylogeny from alignment using
RAxML (external program)."""
    # TODO: too complex, consider breaking up
    input_file = 'phylogeny_in.phylip'
    output_file = 'phylogeny_out'
    file_line = ' -s ' + input_file + ' -n ' + output_file
    options = ' -p ' + str(random.randint(0, 10000000)) + ' -T ' + str(threads)
    if outgroup:
        options += ' -o ' + outgroup
    with open(os.path.join(wd, input_file), "w") as file:
        AlignIO.write(alignment, file, "phylip-relaxed")
    # only use GTRCAT for more than 100 taxa (ref RAxML manual)
    if len(alignment) > 100:
        dnamodel = ' -m GTRCAT'
    else:
        dnamodel = ' -m GTRGAMMA'
    if partitions:
        options += partitions
    if constraint:
        options += constraint
    command_line = raxml + file_line + dnamodel + options
    logger.debug(command_line)
    pipe = TerminationPipe(command_line, silent=True, cwd=wd)
    pipe.run()
    if not pipe.failure:
        try:
            with open(os.path.join(wd, 'RAxML_bestTree.' + output_file), "r") \
                    as file:
                tree = Phylo.read(file, "newick")
        except IOError:
            return None
        finally:
            if constraint:
                os.remove(os.path.join(wd, 'constraint.tre'))
            if partitions:
                os.remove(os.path.join(wd, "partitions.txt"))
            os.remove(os.path.join(wd, input_file))
            all_files = os.listdir(wd)
            for each in all_files:
                if re.search("(RAxML)", each):
                    os.remove(os.path.join(wd, each))
                if re.search("\.reduced$", each):
                    os.remove(os.path.join(wd, each))
        return tree
    else:
        raise RuntimeError()


def consensus(phylogenies, outdir, min_freq=0.5, is_rooted=True,
              trees_splits_encoded=False):
    """Generate a rooted consensus tree"""
    # first ensure that all trees in the distribution have same number
    # of taxa, otherwise, make it so by dropping taxa not present in
    # all trees
    all_tip_names = []
    for phylogeny in phylogenies:
        terminals = phylogeny.get_terminals()
        all_tip_names.append([e.name for e in terminals])
    counted = Counter(sum(all_tip_names, []))
    to_drop = [e for e in counted.keys() if counted[e] < len(phylogenies)]
    for tip_names, phylogeny in zip(all_tip_names, phylogenies):
        dropping = [e for e in tip_names if e in to_drop]
        for tip_name in dropping:
            phylogeny.prune(tip_name)
    with open('.for_consensus.tre', "w") as file:
        Phylo.write(phylogenies, file, 'newick')
    # create dendropy list
    trees = dp.TreeList()
    trees.read_from_path('.for_consensus.tre', "newick", as_rooted=True)
    os.remove('.for_consensus.tre')
    # https://groups.google.com/forum/#!topic/dendropy-users/iJ32ibnS5Bc
    sd = dp.treesplit.SplitDistribution(taxon_set=trees.taxon_set)
    sd.is_rooted = is_rooted
    tsum = dp.treesum.TreeSummarizer()
    tsum.count_splits_on_trees(trees, split_distribution=sd,
                               trees_splits_encoded=trees_splits_encoded)
    consensus = tsum.tree_from_splits(sd, min_freq=min_freq)
    consensus.write_to_path(os.path.join(outdir, 'consensus.tre'), "newick")
