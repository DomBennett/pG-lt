#! /bin/usr/env python
# 10/02/2015
# Check alignments stats for a fasta file

# PACKAGES
import sys
import re
from numpy import mean
from Bio import AlignIO


# FUNCTIONS
def calcOverlap(columns):
    if not columns:
        return 0
    # what proportion of columns have nucs in other seqs
    pcolgaps = []
    for i in columns:
        ith = float(alignment[:, i].count("-"))/(len(alignment)-1)
        pcolgaps.append(ith)
    # overlap is the mean proportion of columns shared
    overlap = len(columns) - (len(columns) * mean(pcolgaps))
    return overlap


def calcNgap(sequence):
    # count the number of gaps
    gaps = re.subn('-+', '', sequence)[1]
    return float(gaps)


# MAIN
if __name__ == '__main__':
    alignment_file = sys.argv[1]
    with open(alignment_file, 'r') as f:
        alignment = AlignIO.read(f, 'fasta')
    alen = alignment.get_alignment_length()
    print('Alignment length [{0}]'.format(alen))
    for each in alignment:
        sequence = str(each.seq)
        columns = [ei for ei, e in enumerate(sequence) if e != "-"]
        overlap = calcOverlap(columns)
        ngap = calcNgap(sequence)
        print('{0}: [{1}] gaps and [{2}] overlap'.
              format(each.name, ngap, overlap))
    print('Done.')
