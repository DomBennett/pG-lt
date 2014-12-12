#! /bin/usr/env python
# D.J. Bennett
# 24/03/2014
"""
pglt Entrez tools
"""

# PACKAGES
import time
import random
from Bio import Entrez
from Bio import SeqIO

# GLOBALS
max_check = 4
download_counter = 0


# FUNCTIONS
def eSearch(term, logger, retStart=0, retMax=1, usehistory="n",
            db="nucleotide"):
    """Use Entrez.esearch to search a term in an NCBI database.

    Arguments:
     term = string of term used in search
     logger = logging object
     retStart = minimum returned ID of matching sequences IDs
     retMax = maximum returned ID of matching sequences IDs
     usehistory = record search in NCBI database ("y" or "n")
     db = NCBI database

    Return:
     dictionary

    Adapted pG code written by W.D. Pearse."""
    # TODO: too complex, consider breaking up
    finished = 0
    global download_counter
    while finished <= max_check:
        if download_counter > 1000:
            logger.info(" ---- download counter hit: waiting 5 minutes ----")
            download_counter = 0
            time.sleep(300)
        try:
            if db is "nucleotide":
                handle = Entrez.esearch(db="nucleotide", term=term,
                                        usehistory=usehistory,
                                        retStart=retStart, retMax=retMax,
                                        retmode="text")
                results = Entrez.read(handle)
                handle.close()
            elif db is "taxonomy":
                handle = Entrez.esearch(db="taxonomy", retStart=retStart,
                                        retmax=retMax, term=term)
                results = Entrez.read(handle)
                handle.close()
            else:
                raise(ValueError('Invalid db argument!'))
            download_counter += 1
            return results
        except ValueError:  # if parsing fails, value error raised
            handle.close()
            logger.warn('Parsing failed!')
            return ()
        except:  # else server error
            if finished == 0:
                logger.debug(" ---- server error: retrying ----")
            elif finished == max_check:
                logger.debug(" ----- server error: no records retrieved ----")
                return()
        time.sleep(60)
        finished += 1


def eFetch(ncbi_id, logger, db="nucleotide"):
    """Download NCBI record(s) using ID number(s).

    Arguments:
     ncbi_id = sequence identifier (list or string)
     logger = logging object
     db = NCBI database (default is nucleotide)

    Return:
     SeqRecord object

    Adapted pG code written by W.D. Pearse."""
    # TODO: too complex, consider breaking up
    finished = 0
    global download_counter
    while finished <= max_check:
        if download_counter > 1000:
            logger.info(" ---- download counter hit: waiting 5 minutes ----")
            download_counter = 0
            time.sleep(300)
        try:
            if db is "nucleotide":
                handle = Entrez.efetch(db="nucleotide", rettype='gb',
                                       retmode="text", id=ncbi_id)
                results_iter = SeqIO.parse(handle, 'gb')
                results = [x for x in results_iter]
                handle.close()
            elif db is "taxonomy":
                handle = Entrez.efetch(db="taxonomy", id=ncbi_id,
                                       retmode="xml")
                results = Entrez.read(handle)
                handle.close()
            else:
                raise(ValueError('Invalid db argument!'))
            download_counter += len(ncbi_id)
            return results
        except ValueError:  # if parsing fails, value error raised
            handle.close()
            logger.warn('Parsing failed!')
            return ()
        except:  # else server error
            if finished == 0:
                logger.debug(" ----- server error: retrying ----")
            elif finished == max_check:
                logger.debug(" ----- server error: no seqs retrieved ----")
                return ()
        time.sleep(60)
        finished += 1


def findChildren(taxid, logger, target=100, next=False):
    """
    Return all decendant genera (or below) of a taxonmic ID.

    Args:
     taxid = taxonomic ID
     logger = logging object
     target = the target number of children returned (default 100)
     next = stop at all children in the rank below given id's rank

    Returns:
     taxid
    """
    # internals
    def findNext(frecord):
        term = "{0}[Next Level] AND {1}[Division]".\
            format(frecord[0]['ScientificName'], frecord[0]['Division'])
        count = eSearch(term, logger=logger, db="taxonomy")["Count"]
        srecord = eSearch(term, logger=logger, db="taxonomy", retMax=count)
        return srecord['IdList']

    def findTillTarget(taxids):
        res = []
        taxids = random.sample(taxids, len(taxids))
        while len(taxids) > 0:
            if len(res) > target:
                break
            taxid = taxids.pop()
            frecord = eFetch(taxid, logger=logger, db="taxonomy")
            if frecord[0]['Rank'] in target_ranks:
                res.append(taxid)
            else:
                res.extend(findTillTarget(findNext(frecord)))
        return res

    # process
    taxid = str(taxid)
    target_ranks = ['genus', 'subgenus', 'species', 'subspecies']
    if next:
        frecord = eFetch(taxid, logger=logger, db="taxonomy")
        return findNext(frecord)
    else:
        return findTillTarget([taxid])
