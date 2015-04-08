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
Entrez.tool = 'pglt'


# FUNCTIONS
def safeConnect(efunc, logger, max_check=100, waittime=1, power=2,
                **kwargs):
    '''Return Entrez results safely'''
    # waitime should be min 1/3 second
    # no more than 3 URL requests per second
    # http://www.ncbi.nlm.nih.gov/books/NBK25497/
    i = 0
    results = ()
    while i < max_check:
        try:
            # open handle with Entrez function
            handle = efunc(**kwargs)
            # print(handle.url)
            # if rettype is GenBank, read each seq into a list
            if 'rettype' in kwargs.keys() and 'gb' == kwargs['rettype']:
                results_iter = SeqIO.parse(handle, 'gb')
                results = [x for x in results_iter]
            else:
                results = Entrez.read(handle)
            handle.close()
            i = max_check
        # catch IOErrors and RuntimeErrors; servers turns down occasionally
        except (IOError, RuntimeError) as errmsg:
            logger.debug(" ---- server error [{0}]: retrying in [{1}s]----".
                         format(errmsg, waittime))
            if i == max_check:
                logger.debug(" ----- max attempts: no records retrieved ----")
            time.sleep(waittime)
            waittime = waittime * power
            i += 1
    return results


def eSearch(term, logger, retStart=0, retMax=1, db="nucleotide"):
    """Use Entrez.esearch to search a term in an NCBI database.

    Arguments:
     term = string of term used in search
     logger = logging object
     retStart = minimum returned ID of matching sequences IDs
     retMax = maximum returned ID of matching sequences IDs
     db = NCBI database

    Return:
     dictionary

    Adapted pG code written by W.D. Pearse."""
    if db not in ['nucleotide', 'taxonomy']:
        raise(ValueError('Invalid db argument!'))
    results = ()
    results = safeConnect(efunc=Entrez.esearch, logger=logger, db=db,
                          term=term, usehistory='n', retStart=retStart,
                          retMax=retMax, retmode="text")
    return results


def eFetch(ncbi_id, logger, db="nucleotide"):
    """Download NCBI record(s) using ID number(s).

    Arguments:
     ncbi_id = sequence identifier (list or string)
     logger = logging object
     db = NCBI database (default is nucleotide)

    Return:
     List of SeqRecords (db = 'nucleotide')
     List of dictionaries (db = 'taxonomy')

    Adapted pG code written by W.D. Pearse."""
    if db not in ['nucleotide', 'taxonomy']:
        raise(ValueError('Invalid db argument!'))
    results = ()
    if db == 'taxonomy':
        results = safeConnect(efunc=Entrez.efetch, logger=logger, db=db,
                              retmode='xml', id=ncbi_id)
    else:
        results = safeConnect(efunc=Entrez.efetch, logger=logger, db=db,
                              rettype='gb', retmode='text', id=ncbi_id)
    return results


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
