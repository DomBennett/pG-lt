#! /bin/usr/env python
# D.J. Bennett
# 24/03/2014
"""
pglt alignment tools
"""

# PACKAGES
import os
import re
import random
import numpy as np
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Application import ApplicationError
from Bio.Blast import NCBIXML
from StringIO import StringIO
from system_tools import TerminationPipe
from system_tools import OutgroupError
from system_tools import TooFewSpeciesError
from system_tools import MafftError
from system_tools import TrysError
from special_tools import timeit
from special_tools import getThreads
from pglt import _MAFFT as mafft
from pglt import _MAFFTQ as mafftq
from pglt import _MAFFTX as mafftx
from pglt import _BLASTN as blastn


# OBEJECTS
class SeqStore(dict):
    """Store species' gene sequences with functions for pulling \
sequences for alignments and adding penalties for sequences that did \
not align"""
    def __init__(self, genedir, seqfiles, minfails, mingaps, minoverlap,
                 logger, wd=os.getcwd()):
        self.wd = wd
        self.logger = logger
        self.threads = getThreads(wd)
        self.minfails = minfails  # minimum number of fails in a row
        self.dspp = []  # species dropped
        self.nseqs = 0  # counter for seqs
        self.blast_prop = 0.5  # the p sequences a sequence must overlap
        self.mingaps = mingaps
        self.minoverlap = minoverlap
        for i, seqfile in enumerate(seqfiles):
            name = re.sub('\.fasta$', '', seqfile)
            seqdir = os.path.join(genedir, seqfile)
            seqs = []
            lengths = []
            with open(seqdir, "rU") as infile:
                for record in SeqIO.parse(infile, "fasta"):
                    record.id = name
                    lengths.append(len(record))
                    seqs.append([record, 0])  # seqrecord + nfails
                    self.nseqs += 1
            if len(seqs) > 0:
                self[name] = [seqs, np.min(lengths)]

    def _add(self, sequences=None, limit=None):
        """Return a random sequence for alignment"""
        # if there are sequences, use blast alignment
        if self.sequences_in_alignment:
            # if add is limited to certain species, use only its i
            if limit and limit in self.sppool:
                rand_ints = [ei for ei, e in enumerate(self.sppool) if
                             e == limit]
            else:
                rand_ints = range(len(self.sppool))
                random.shuffle(rand_ints)
            for i in rand_ints:
                sp = self.sppool[i]
                next_seqs = [e[0] for e in self[sp][0]]
                # blast next_seqs against sequences in alignment
                res = self._alignmentBlast(next_seqs, sequences)
                # if success break
                if res:
                    break
            else:
                return None
            self.next_sp = self.sppool.pop(i)
            next_seq_i = res[0]
            next_seq = res[1]
            # record sequence + nfails in sequence_in_alignment
            self.sequences_in_alignment.append(self[self.next_sp][0]
                                               [next_seq_i])
        else:
            # without any sequences in alignment
            #  simply return a random seq
            rand_int = random.randint(0, (len(self.sppool)-1))
            self.next_sp = self.sppool.pop(rand_int)
            result = random.sample(self[self.next_sp][0], 1)[0]
            # record sequence + nfails in sequence_in_alignment
            self.sequences_in_alignment.append(result)
            next_seq = result[0]
        return next_seq

    def _alignmentBlast(self, query, sequences_in_alignment):
        """Return indexes and overlapping sequences for each sequence
in query that overlaps with more than prop sequences in
sequences_in_alignment given set parameters using NCBI's BLAST"""
        # loop through each sequence in query, if success, return
        #  overlapping sequence and its index
        # make sure indexes are randomised to avoid biased sampling
        indexes = random.sample(range(len(query)), len(query))
        for i in indexes:
            # blast prospective next sequence against all sequences
            #  in alignment
            bools, positions = blast(query[i], sequences_in_alignment,
                                     self.minoverlap, self.logger, self.wd,
                                     self.threads)
            # if more than prop overlap ...
            overlap = (float(sum(bools))/len(sequences_in_alignment)) >\
                self.blast_prop
            if overlap:
                # ... return its index in the seqstore and the
                #  overlapping sequence
                return i, query[i][min(positions):max(positions)]

    def start(self, n):
        """Return n starting random sp sequences, update sppool"""
        # pool starts with all species
        self.sppool = self.keys()
        # sequences + nfails used in alignment
        self.sequences_in_alignment = []
        # actual sequences for alignment
        sequences = []
        # add a random seq to sequences_in_alignment and sequences
        sequences.append(self._add())
        while True:
            sequence = self._add(sequences)
            if sequence:
                sequences.append(sequence)
            # break from loop if n is hit or sequence is None
            if not sequence or len(self.sequences_in_alignment) == n:
                break
        return sequences

    def back(self, alignment):
        """Add to nfails for random sequence, return new random \
species"""
        self.sequences_in_alignment[-1][1] += 1
        del self.sequences_in_alignment[-1]
        self.sppool.append(self.next_sp)
        self._check()
        if len(self.sppool) != 0:
            sequences = [e for e in alignment]
            return self._add(sequences)
        else:
            return None

    def next(self, alignment, limit=None):
        """Set nfails to 0, add additional random or limited species"""
        for i in range(len(self.sequences_in_alignment)):
            self.sequences_in_alignment[i][1] = 0
        sequences = [e for e in alignment]
        if limit:
            return self._add(sequences=sequences, limit=limit)
        else:
            return self._add(sequences=sequences)

    def _check(self):
        """Check nfails, drop sequences and species"""
        spp = self.keys()
        for sp in spp:
            to_drop = [ei for ei, e in enumerate(self[sp][0]) if e[1] >
                       self.minfails]
            for i in to_drop:
                self.logger.info("Dropping [{0}] for [{1}] as nfails is \
[{2}]".format(self[sp][0][i][0].description, sp, self[sp][0][i][1]))
            self[sp][0] = [e for ei, e in enumerate(self[sp][0])
                           if ei not in to_drop]
            if len(self[sp][0]) < 1:
                if sp == "outgroup":
                    raise OutgroupError
                else:
                    self.logger.info("Dropped [{0}]".format(sp))
                    self.dspp.append(sp)
                    self.sppool = [e for e in self.sppool if e != sp]
                    del self[sp]
        if len(self.keys()) < 5:
            raise TooFewSpeciesError


class Aligner(object):
    """Build alignments from seqstore"""
    def __init__(self, seqstore, mingaps, minoverlap, minseedsize,
                 maxseedsize, maxtrys, maxseedtrys, gene_type, outgroup,
                 logger, wd=os.getcwd()):
        self.wd = wd
        self.logger = logger
        self.threads = getThreads(wd=wd)
        self.seqstore = seqstore
        self.mingaps = mingaps
        self.minoverlap = minoverlap
        self.minseedsize = minseedsize
        self.maxtrys = 2  # trys for alignment attempts
        self.buffer = maxseedtrys  # trys for a seedsize
        self.buffer_counter = 0  # seedsize buffer counter
        self.seedsize = len(seqstore)
        self.timeout = 99999999
        self.talign = False
        self.tadd = False
        self.silent = False
        self.total_trys = 0  # counter for total number of trys
        self.type = gene_type
        self.outgroup = outgroup

    def _calcTimeout(self, seconds, alignment, align=True):
        """Calculate the timeout"""
        # sequences that align takes less time to finish
        # make the most of this fact by capping the time MAFFT will
        # run for based on time taken for successful alignments
        # size if mean length * number of sequences
        mean_len = np.mean([len(e) for e in alignment])
        size = mean_len*len(alignment)
        # add a *10 buffer
        timeout = (seconds/size)*10
        if align:
            # only update if new timeout is bigger
            if self.talign < timeout:
                self.talign = timeout
        else:
            if self.tadd < timeout:
                self.tadd = timeout

    def _getTimeout(self, sequences, sequence=None):
        """Return seconds it should take for mafft to run"""
        # Only return modified timeout if _calcTimeout has been run
        # Use size of alignment to calc seconds per nuc
        mean_len = np.mean([len(e) for e in sequences])
        size = mean_len*len(sequences)
        if sequence and self.tadd:
            return self.tadd*(size + len(sequence))
        elif self.talign:
            return self.talign*size
        return self.timeout

    def _check(self, alignment):
        return checkAlignment(alignment, self.mingaps, self.minoverlap,
                              self.minlen, self.logger)

    def _return(self):
        """Return best alignment from a list of alignments based on:\
presence of outgroup, number of species and length of alignment"""
        if not self.store:
            self.total_trys += 1
            return None
        if self.outgroup:
            # keep alignments with outgroups
            # keep alignments with more than 5 species
            self.store = [a for a in self.store if "outgroup" in
                          [e.id for e in a._records]]
            if len(self.store) == 0:
                self.logger.debug("........ no outgroup")
                self.total_trys += 1
                return None
        self.store = [e for e in self.store if len(e) >= 5]
        if len(self.store) == 0:
            self.logger.debug("........ too few species")
            self.total_trys += 1
            return None
        # keep only alignments with lots of records
        nrecords = [len(e._records) for e in self.store]
        self.store = [self.store[i] for i, e in enumerate(nrecords)
                      if e == max(nrecords)]
        # return longest alignment
        lens = [len(e) for e in self.store]
        max_i = lens.index(max(lens))
        self.total_trys = 0
        return self.store[max_i]

    def _calcSeedsize(self, success=True):
        """Calculate seedsize based on buffer and success of current \
seedsize. Return 1 if trys must increase, else 0."""
        # increase seedsize if successful buffer times in a row
        # decrease seedsize if unsuccessful buffer times in a row
        if success:
            self.buffer_counter += 1
        else:
            self.buffer_counter -= 1
        if self.buffer_counter >= self.buffer:
            if self.seedsize < len(self.seqstore):
                self.seedsize += 1
                # reset every time seedsize changes
                self.buffer_counter = 0
        if self.buffer_counter <= (self.buffer * -1):
            if self.seedsize > self.minseedsize:
                self.seedsize -= 1
                self.buffer_counter = 0
            else:
                # here seedsize must be reduced, but has hit
                #  minseedsize, add 1 to trys
                return 1
        return 0

    def _seed(self, trys):
        """BLAST seedsize sequences together, align and return True
if successful"""
        self.minlen = min([self.seqstore[e][1] for e in self.seqstore.keys()])
        sequences = self.seqstore.start(self.seedsize)
        # make sure there are enough seqs for alignment
        if len(sequences) >= self.minseedsize:
            command = version(sequences, self.type)
            try:
                alignment, seconds = timeit(func=align, command=command,
                                            sequences=sequences,
                                            logger=self.logger, wd=self.wd,
                                            threads=self.threads,
                                            timeout=self._getTimeout(sequences)
                                            )
            except MafftError:
                self.logger.debug('MAFTT error raised')
                success = False
            else:
                success = self._check(alignment)
        else:
            success = False
        # add to trys if unsuccessful or sequences are fewer than
        #  seedsize
        trys += self._calcSeedsize(success and len(sequences) == self.seedsize)
        if success:
            self._calcTimeout(seconds, alignment)
            self.store.append(alignment)
        return success, trys

    def _add(self, trys, limit=None):
        """Add sequence to alignment, return True if successful"""
        alignment = self.store[-1]
        self.minlen = min([self.seqstore[e][1] for e in self.seqstore.keys()])
        if len(self.seqstore.sppool) == 0:
            return True, trys
        else:
            sequence = self.seqstore.next(alignment, limit=limit)
            if not sequence:
                # if no sequence is returned, nothing more can be
                #  added
                self.logger.debug('No new sequence added')
                return True, trys
        try:
            new_alignment, seconds = timeit(func=add, alignment=alignment,
                                            sequence=sequence,
                                            logger=self.logger, wd=self.wd,
                                            threads=self.threads,
                                            timeout=self._getTimeout(alignment,
                                                                     sequence))
        except MafftError:
            self.logger.debug('MAFTT error raised')
            success = False
        else:
            success = self._check(new_alignment)
        if success:
            self._calcTimeout(seconds, alignment, align=False)
            self.store.append(new_alignment)
        else:
            trys += 1
            sequence = self.seqstore.back(alignment)
            if not sequence:
                # here a species has been dropped and now all
                #  species are present
                return True, trys
        return False, trys

    def run(self):
        """Incrementally build an alignment by adding sequences to a \
seed alignment"""
        self.store = []
        if self.total_trys > self.maxtrys:
            # if run has been run maxtrys in a row wo success
            #  raise trys error, if outgroup isn't being used
            if self.outgroup:
                self.outgroup = False
            else:
                raise TrysError
        # seed
        self.logger.info("........ seed phase: [{0}] seed size".
                         format(self.seedsize))
        trys = 0
        success = False
        while not success:
            success, trys = self._seed(trys)
            if trys > self.maxtrys:
                self.logger.debug("............ maxtrys hit")
                return self._return()
        # add
        self.logger.info("........ add phase : [{0}] species".
                         format(len(self.store[-1])))
        trys = 0
        finished = False
        # if outgroup, force it for the first add after seed
        if self.outgroup:
            finished, trys = self._add(trys, limit='outgroup')
        while not finished:
            finished, trys = self._add(trys)
            if trys > self.maxtrys:
                self.logger.debug("............ maxtrys hit")
                return self._return()
        return self._return()


# FUNCTIONS
def version(sequences, gene_type):
    """Return command for best MAFFT version given sequences"""
    # determine auto, qinsi or xinsi based on:
    # http://mafft.cbrc.jp/alignment/software/source66.html
    # always using default algorithms
    if not mafftq and mafftx:  # return mafft, if no mafftq or x
        return mafft
    if gene_type != 'deep' or len(sequences) > 250:
        return mafft + ' --auto'
    seqlens = [len(s) for s in sequences]
    if max(seqlens) > 1500:
        return mafft + ' --auto'
    if len(sequences) > 60:
        return mafftq
    return mafftx


def genNonAlignment(nseqs, alen):
    """Return non-alignment, for when align or add timeout"""
    seqs = [SeqRecord(Seq('-' * alen), id='Seq{0}'.format(e),
                      description='timeout non-sequence') for
            e in range(nseqs)]
    return MultipleSeqAlignment(seqs)


def align(command, sequences, timeout, logger, wd, threads):
    """Adapted pG function: Align sequences using mafft (external
program)"""
    input_file = "sequences_in.fasta"
    output_file = "alignment_out.fasta"
    command_line = '{0} --thread {1} {2} > {3}'.format(command, threads,
                                                       input_file, output_file)
    with open(os.path.join(wd, input_file), "w") as file:
        SeqIO.write(sequences, file, "fasta")
    logger.debug(command_line)
    pipe = TerminationPipe(command_line, timeout=timeout, cwd=wd)
    pipe.run()
    os.remove(os.path.join(wd, input_file))
    if not pipe.failure:
        try:
            res = AlignIO.read(os.path.join(wd, output_file), 'fasta')
        except:
            logger.info(pipe.output)
            raise MafftError()
        else:
            os.remove(os.path.join(wd, output_file))
    else:
        # if pipe.failure, runtime error, return non-alignment
        logger.debug('.... align timeout ....')
        return genNonAlignment(len(sequences), len(sequences[0]))
    return res


def add(alignment, sequence, timeout, logger, wd, threads):
    """Align sequence(s) to an alignment using mafft (external
program)"""
    alignment_file = "alignment_in.fasta"
    sequence_file = "sequence_in.fasta"
    output_file = "alignment_out.fasta" + '.fasta'
    command_line = '{0} --auto --thread {1} --add {2} {3} > {4}'.\
                   format(mafft, threads, sequence_file, alignment_file,
                          output_file)
    with open(os.path.join(wd, sequence_file), "w") as file:
        SeqIO.write(sequence, file, "fasta")
    with open(os.path.join(wd, alignment_file), "w") as file:
        AlignIO.write(alignment, file, "fasta")
    pipe = TerminationPipe(command_line, timeout=timeout, cwd=wd)
    pipe.run()
    os.remove(os.path.join(wd, alignment_file))
    os.remove(os.path.join(wd, sequence_file))
    if not pipe.failure:
        try:
            res = AlignIO.read(os.path.join(wd, output_file), 'fasta')
        except:
            logger.info(pipe.output)
            raise MafftError()
        else:
            os.remove(os.path.join(wd, output_file))
    else:
        logger.debug('.... add timeout ....')
        return genNonAlignment(len(alignment) + 1,
                               len(alignment.get_alignment_length()))
    return res


def blast(query, subj, minoverlap, logger, wd, threads):
    """Return bool and positions of query sequences that overlapped
with subject given parameters."""
    query_file = os.path.join(wd, 'query.fasta')
    subj_file = os.path.join(wd, 'subj.fasta')
    SeqIO.write(query, query_file, "fasta")
    SeqIO.write(subj, subj_file, "fasta")
    try:
        # options: http://www.ncbi.nlm.nih.gov/books/NBK1763/
        cline = NcbiblastnCommandline(query=query_file, subject=subj_file,
                                      outfmt=5, cmd=blastn, word_size=8,
                                      num_threads=threads)
        logger.debug(cline)
        output = cline()[0]
    except ApplicationError:  # as error_msg:
        # logger.debug(error_msg)
        # logger.warn("---- BLAST Error ----")
        # TODO: work out why this is happening, doesn't seem to affect
        #  results though, low priority
        return [], []
    finally:
        os.remove(query_file)
        os.remove(subj_file)
    # list of T or F for success of alignment between queries and
    #  subject
    bools = []
    # record start and end position to avoid composite sequence
    #  problems
    positions = []
    # BLAST records for each query sequence matched against subj
    bresults = NCBIXML.parse(StringIO(output))
    for record in bresults:
        if record.alignments:
            res = record.alignments[0].hsps[0]
            # if identities > minoverlap, keep
            if res.identities > minoverlap:
                bools.append(True)
                positions.append(res.query_start)
                positions.append(res.query_end)
                continue
        bools.append(False)
    return bools, positions


def checkAlignment(alignment, mingaps, minoverlap, minlen, logger):
    """Determine if an alignment is good or not based on given \
parameters. Return bool"""
    # TODO: too complex, consider breaking up
    # internals
    def calcOverlap(columns):
        if not columns:
            return 0
        # what proportion of columns have nucs in other seqs
        pcolgaps = []
        for i in columns:
            ith = float(alignment[:, i].count("-"))/(len(alignment)-1)
            pcolgaps.append(ith)
        # overlap is the mean proportion of columns shared
        overlap = len(columns) - (len(columns) * np.mean(pcolgaps))
        return overlap

    def calcNgap(sequence):
        # count the number of gaps
        gaps = re.subn('-+', '', sequence)[1]
        return float(gaps)

    # process
    if alignment is None:
        return False
    alen = alignment.get_alignment_length()
    if alen < minlen:
        logger.debug('........ alignment too small')
        return False
    for each in alignment:
        sequence = str(each.seq)
        columns = [ei for ei, e in enumerate(sequence) if e != "-"]
        overlap = calcOverlap(columns)
        if overlap < minoverlap:
            logger.debug('........ alignment too little overlap')
            return False
        ngap = calcNgap(sequence)
        if ngap > mingaps:
            logger.debug('........ alignment too many gaps')
            return False
    return True
